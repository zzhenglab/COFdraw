"""Test the extracted Mac download, including the native window and a CPU job."""
import json
from pathlib import Path
import platform
import subprocess
import sys
from tempfile import TemporaryDirectory
import time
from urllib.request import Request, urlopen


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
    drawing = api("render", {"source_recipe_id": "cof-300-projection-v1", "size": 640})
    with urlopen(url.rstrip("/") + drawing["downloads"]["png"], timeout=30) as response:
        assert response.read(8) == bytes.fromhex("89504e470d0a1a0a")
    model = api("science/model", {"template_id": "cof-300-c5"})
    assert model["structure"]["periodic"] and model["viewer"]["bonds"]
    with urlopen(url.rstrip("/") + model["downloads"]["cif"], timeout=30) as response:
        assert b"_cell_length_a" in response.read()
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
    print(f"Verified {len(library['monomer_presets'])} monomers, COF-300 PNG and CIF, and argon Widom/GCMC.")


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
        subprocess.run([str(exe), "--smoke-window", "--status-file", str(status)], check=True, timeout=120)
        proof = json.loads(status.read_text())
        assert proof["status"] == "desktop_smoke_passed" and proof["native_save_bridge"]
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
