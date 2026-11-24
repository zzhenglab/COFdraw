"""Test the extracted Mac download, including the native window and a CPU job."""
from collections import defaultdict
import json
from math import atan2, degrees, hypot
from pathlib import Path
import platform
import subprocess
import sys
from tempfile import TemporaryDirectory
import time
from urllib.request import Request, urlopen

from rdkit import Chem


def downloaded_molecule(url, drawing):
    with urlopen(url.rstrip("/") + drawing["downloads"]["mol"], timeout=30) as response:
        return Chem.MolFromMolBlock(response.read().decode(), sanitize=False, removeHs=False)


def check_hcb_imines(api, url, presets):
    smiles = {row["id"]: row["smiles"] for row in presets}
    recipe = api("prepare", {"monomer_a": smiles["tfb"], "monomer_b": smiles["pa"],
                            "linkage": "imine", "template_id": "hcb_3_2", "instances": 12})["recipe"]
    drawing = api("render", {"recipe": recipe, "size": 640})
    mol = downloaded_molecule(url, drawing)
    xy = mol.GetConformer().GetPositions()[:, :2]
    owners = drawing["audit"]["atom_origin"]
    # Six TFB nodes each have two joined imines and one outward continuation (=*).
    # The six diamine linkers each have two joined nitrogen ends.
    for pattern, arms in (("[#6]-[#6]=[#7,#0]", 3), ("[#6]-[#7]=[#6]", 2)):
        by_instance = defaultdict(list)
        for anchor, site, tip in mol.GetSubstructMatches(Chem.MolFromSmarts(pattern)):
            incoming = xy[site] - xy[anchor]
            outgoing = xy[tip] - xy[site]
            cross = incoming[0] * outgoing[1] - incoming[1] * outgoing[0]
            turn = degrees(atan2(cross, incoming @ outgoing))
            owner = owners[site]["instance"]
            assert abs(abs(turn) - 60) < 1, f"{owner} imine angle {180 - abs(turn):.1f}°"
            by_instance[owner].append(turn)
        assert len(by_instance) == 6, f"Expected six instances matching {pattern}"
        for owner, turns in by_instance.items():
            assert len(turns) == arms and min(turns) * max(turns) > 0, (owner, turns)

    # A phenyl substituent makes the two para aldehyde ends inequivalent.
    # The four translated pores share nineteen distinct straight linkers.
    recipe = api("prepare", {"monomer_a": "O=Cc1c(-c2ccccc2)cc(C=O)cc1", "monomer_b": smiles["tapb"],
                            "linkage": "imine", "template_id": "hcb_3_2", "instances": 12})["recipe"]
    sheet = api("repeat", {"recipe": recipe, "repeats": 2, "size": 640})
    assert sheet["summary"]["pore_count"] == 4
    mol = downloaded_molecule(url, sheet)
    xy = mol.GetConformer().GetPositions()[:, :2]
    para_linker = Chem.MolFromSmarts("[#7,#0]=[#6]-[#6]1~[#6]~[#6]~[#6](~[#6]~[#6]~1)-[#6]=[#7,#0]")
    matches = mol.GetSubstructMatches(para_linker)
    assert len(matches) == 19, f"Expected 19 shared and boundary linkers, found {len(matches)}"
    for left_tip, left, *ring, right, right_tip in matches:
        axis = xy[right] - xy[left]
        normal = [-axis[1], axis[0]]
        # Round-off can put two nearly collinear bonds on opposite sides of zero.
        sides = []
        for site, tip in ((left, left_tip), (right, right_tip)):
            bond = xy[tip] - xy[site]
            sides.append((bond @ normal) / (hypot(*axis) * hypot(*bond)))
        assert min(sides) < -.05 and max(sides) > .05, f"Unresolved or same-side bends at {left} and {right}"


def check_service(url):
    def api(path, payload=None):
        request = Request(url + "api/drawing/" + path,
                          data=json.dumps(payload).encode() if payload is not None else None,
                          headers={"Content-Type": "application/json"})
        with urlopen(request, timeout=90) as response:
            return json.load(response)

    library = api("library")
    assert library["monomer_presets"] and not library["source_registry_errors"]
    assert any(m["synthetic_accessibility"]["available"] for m in library["monomer_presets"])
    groups = {row["id"]: {g["id"] for g in row["linkage_functional_groups"]}
              for row in library["monomer_presets"]}
    assert groups["dmtp"] == {"aldehyde"} and groups["tapt"] == {"primary_amine"}
    assert groups["bfbapy"] == {"aldehyde", "primary_amine"}
    check_hcb_imines(api, url, library["monomer_presets"])
    drawing = api("render", {"source_recipe_id": "cof-300-projection-v1", "size": 640})
    with urlopen(url.rstrip("/") + drawing["downloads"]["png"], timeout=30) as response:
        assert response.read(8) == bytes.fromhex("89504e470d0a1a0a")
    model = api("science/model", {"template_id": "cof-300-c5"})
    assert model["structure"]["periodic"] and model["viewer"]["bonds"]
    assert "<svg" in model["chemical_circuit_svg"] and "<svg" in model["chemical_drawing_svg"]
    with urlopen(url.rstrip("/") + model["downloads"]["cif"], timeout=30) as response:
        assert b"_cell_length_a" in response.read()
    variant = api("science/model", {"template_id": "design-dia-dmta"})
    assert variant["structure"]["provenance"]["monomer_ids"] == ["tapm", "dmtp"]
    assert variant["structure"]["atoms"] != model["structure"]["atoms"] and variant["viewer"]["periodic"]
    assert api("science/adsorption/engine")["available"]
    structure = api("science/model", {"template_id": "hcb-imine-aa"})["structure"]
    job = api("science/adsorption/jobs", {"structure": structure, "gases": ["Ar"], "charge_mode": "zero",
              "pressures_Pa": [100000], "cycles": 100, "init_cycles": 100, "widom_cycles": 100,
              "calculate_heat": False})
    deadline = time.monotonic() + 240
    while job["status"] in {"queued", "running"} and time.monotonic() < deadline:
        time.sleep(1)
        job = api("science/adsorption/jobs/" + job["job_id"])
    assert job["status"] == "complete", (job["status"], job.get("message"))
    print(f"Verified {len(library['monomer_presets'])} monomers, hcb imine angles and node symmetry, "
          "opposed imines across four neighboring hcb pores, "
          "COF-300 PNG and CIF, and argon Widom/GCMC.")


def main():
    archive = Path(sys.argv[1]).resolve()
    with TemporaryDirectory(prefix="COFdraw Mac test ") as directory:
        extracted = Path(directory)
        subprocess.run(["ditto", "-x", "-k", str(archive), directory], check=True)
        app = extracted / "COFdraw.app"
        exe = app / "Contents/MacOS/COFdraw"
        assert subprocess.check_output(["lipo", "-archs", str(exe)], text=True).strip() == platform.machine()
        subprocess.run(["codesign", "--verify", "--deep", "--strict", str(app)], check=True)
        status = extracted / "native.json"
        launch = subprocess.run([str(exe), "--smoke-window", "--status-file", str(status)], timeout=120)
        proof = json.loads(status.read_text())
        assert launch.returncode == 0 and proof["status"] == "desktop_smoke_passed" and proof["native_save_bridge"], proof
        service = extracted / "service.json"
        proc = subprocess.Popen([str(exe), "--headless", "--status-file", str(service)])
        try:
            deadline = time.monotonic() + 60
            while proc.poll() is None and time.monotonic() < deadline:
                try:
                    url = json.loads(service.read_text())["url"]
                    with urlopen(url + "api/desktop", timeout=2) as response:
                        assert json.load(response)["offline"]
                    break
                except (OSError, ValueError, KeyError):
                    time.sleep(.25)
            else:
                raise RuntimeError(f"Packaged service did not start, exit {proc.poll()}")
            check_service(url)
        finally:
            proc.terminate()
            proc.wait(timeout=15)
        subprocess.run(["codesign", "--verify", "--deep", "--strict", str(app)], check=True)


if __name__ == "__main__":
    main()
