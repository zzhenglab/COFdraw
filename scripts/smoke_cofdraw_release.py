"""Exercise the frozen executable, its native window, and real browser actions."""
import argparse
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sys
import time
from urllib.request import Request, urlopen
import xml.etree.ElementTree as ET
from smoke_cofdraw_geometry import check_geometry_examples

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--exe", type=Path, required=True)
    ap.add_argument("--data-template", type=Path, required=True)
    ap.add_argument("--output-dir", type=Path, required=True)
    args = ap.parse_args()
    output = args.output_dir.resolve(); output.mkdir(parents=True, exist_ok=True)
    runtime = output / "runtime-data"
    shutil.copytree(args.data_template, runtime, dirs_exist_ok=True)
    executable = args.exe.resolve()
    report = {"exe_sha256": hashlib.sha256(executable.read_bytes()).hexdigest(), "checks": [], "failure": None, "training": False, "model_inference": False}
    hidden = getattr(subprocess, "CREATE_NO_WINDOW", 0)
    server = None
    try:
        status_file = output / "server.json"
        server = subprocess.Popen([str(executable), "--headless", "--data-dir", str(runtime), "--status-file", str(status_file)], creationflags=hidden)
        for _ in range(120):
            if server.poll() is not None:
                raise RuntimeError(f"Packaged server exited {server.returncode}")
            try:
                status = json.loads(status_file.read_text(encoding="utf-8"))
                if status["pid"] != server.pid:
                    raise ValueError("stale status file")
                base = status["url"].rstrip("/")
                with urlopen(base + "/api/desktop", timeout=2) as response:
                    assert json.load(response)["offline"]
                break
            except (OSError, ValueError, KeyError):
                time.sleep(.25)
        else:
            raise RuntimeError("Packaged service did not become ready")
        def api(path, data=None):
            request = Request(base + path, data=None if data is None else json.dumps(data).encode(), headers={"Content-Type": "application/json"})
            with urlopen(request, timeout=120) as response:
                return json.load(response)
        library = api("/api/drawing/library")
        assert len(library["source_recipes"]) == 17 and not library["source_registry_errors"]
        assert {'gradient_stack','sphere_cutaway'} <= {row['id'] for row in library["styles"]}
        expected_catalog=json.loads((args.data_template/'framework-drawing/catalog/current.json').read_text(encoding='utf-8'))
        expected_monomers=len(expected_catalog['monomer_presets'])
        assert {r['id'] for r in library['monomer_presets']}=={r['id'] for r in expected_catalog['monomer_presets']}
        assert len(library["templates"]) == 7
        smiles = {row["id"]: row["smiles"] for row in library["monomer_presets"]}
        for a, b, linkage, template, method in (
            ("ta", "pytta", "imine", "sql_4_2", "phenyl_torsions"),
            ("tp", "deth", "keto_hydrazide", "hcb_3_2", "restrained_mmff94"),
        ):
            drawing = api("/api/drawing/prepare", {"monomer_a": smiles[a], "monomer_b": smiles[b],
                          "linkage": linkage, "template_id": template})
            built = api("/api/drawing/science/build-periodic", {"recipe": drawing["recipe"]})
            assert built["available"] and built["structure"]["provenance"]["fragment_geometry"] == method
            assert built["cell_validation"]["short_contacts"] == 0 and built["downloads"]["cif"]
        report["recovered_periodic_cells"] = ["ta+pytta", "tp+deth"]
        report["library"] = {key: len(library[key]) for key in ("styles", "monomer_presets", "templates", "topology_catalog", "linkages")}
        def download(path):
            with urlopen(base + path, timeout=30) as response:
                return response.read()
        report["geometry_checks"] = check_geometry_examples(api, download, output / "geometry")
        for source in library["source_recipes"]:
            start = time.monotonic()
            result = api("/api/drawing/render", {"source_recipe_id": source["recipe_id"], "style_id": "red_blue", "size": 800})
            assert result["status"] == "source_bound_provisional"
            if result["audit"].get("cofex_string"):
                assert " COFex." in result["audit"]["cofex_string"]
            assert set(result["downloads"]) == {"svg", "png", "mol", "sdf", "cdxml", "json"}
            files = {}
            for ext, path in result["downloads"].items():
                with urlopen(base + path, timeout=30) as response:
                    payload = response.read()
                saved = output / f"{source['recipe_id']}.{ext}"; saved.write_bytes(payload)
                files[ext] = hashlib.sha256(payload).hexdigest()
                if ext == "cdxml":
                    xml = ET.fromstring(payload)
                    nodes = list(xml.iter("n"))
                    assert len(nodes) == result["audit"]["atoms"]
                    assert not any(n.get("Element") == "0" for n in nodes)
                    assert any(n.get("NodeType") == "ExternalConnectionPoint" for n in nodes)
                    pages = list(xml.iter("page"))
                    assert len(pages) == 1 and pages[0].get("BoundingBox") == "0 0 540 720"
                    assert pages[0].get("WidthPages") == pages[0].get("HeightPages") == "1"
                    for node in nodes:
                        x, y = map(float, node.get("p").split())
                        assert 0 < x < 540 and 0 < y < 720
                    assert result["audit"]["chemical_exports"]["page_geometry"]["paper_size_points"] == [612, 792]
            report["checks"].append({"source": source["recipe_id"], "atoms": result["audit"]["atoms"], "elapsed_seconds": round(time.monotonic()-start, 2), "downloads_sha256": files})
        subprocess.run([sys.executable, str(Path(__file__).with_name("smoke_cofdraw_browser.py")), "--base-url", base,
                        "--output-dir", str(output / "browser")], check=True, timeout=600, creationflags=hidden)
        report["browser_report"] = json.loads((output / "browser/cofdraw-browser-report.json").read_text(encoding="utf-8"))
        agent_input = output / "solubility-request.json"
        agent_input.write_text(json.dumps({"monomer_id": "tfb", "solvent_ids": ["dioxane"], "solvent_ml": 2}), "utf-8")
        agent_output = output / "agent-solubility"
        subprocess.run([sys.executable, str(executable.parent / "skills/cofdraw/scripts/cofdraw_client.py"),
                        "--exe", str(executable), "--output", str(agent_output),
                        "science", "solubility", "--request", str(agent_input)],
                       check=True, timeout=120, creationflags=hidden)
        agent_result = json.loads((agent_output / "result.json").read_text(encoding="utf-8"))
        assert agent_result["solvents"][0]["status"] == "reported"
        assert (agent_output / "structure.csv").is_file()
        report["agent_solubility"] = {"passed": True, "report_id": agent_result["report_id"]}
        subprocess.run([sys.executable, str(Path(__file__).with_name("smoke_cofdraw_compatibility.py")), "--base-url", base,
                        "--expected-monomers", str(expected_monomers), "--output-dir", str(output / "compatibility")], check=True, timeout=600, creationflags=hidden)
        report["compatibility_report"] = json.loads((output / "compatibility/compatibility-report.json").read_text(encoding="utf-8"))
        assert report["compatibility_report"]["passed"]
        subprocess.run([sys.executable, str(Path(__file__).with_name("smoke_cofdraw_identity.py")), "--base-url", base,
                        "--output-dir", str(output / "identity-auto")], check=True, timeout=600, creationflags=hidden)
        report["identity_auto_report"] = json.loads((output / "identity-auto/identity-browser-report.json").read_text(encoding="utf-8"))
        assert report["identity_auto_report"]["passed"]
        subprocess.run([sys.executable, str(Path(__file__).with_name("smoke_cofdraw_synthesis.py")), "--base-url", base,
                        "--output-dir", str(output / "synthesis")], check=True, timeout=600, creationflags=hidden)
        report["synthesis_report"] = json.loads((output / "synthesis/synthesis-browser-report.json").read_text(encoding="utf-8"))
        assert report["synthesis_report"]["passed"]
        subprocess.run([sys.executable, str(Path(__file__).with_name("smoke_cofdraw_studio.py")), "--base-url", base,
                        "--output-dir", str(output / "studio")], check=True, timeout=600, creationflags=hidden)
        report["studio_report"] = json.loads((output / "studio/studio-smoke-report.json").read_text(encoding="utf-8"))
        subprocess.run([sys.executable, str(Path(__file__).with_name("smoke_cofdraw_structures.py")), "--base-url", base,
                        "--output-dir", str(output / "structures")], check=True, timeout=600, creationflags=hidden)
        report["structures_report"] = json.loads((output / "structures/structures-smoke-report.json").read_text(encoding="utf-8"))
        assert report["structures_report"]["passed"]
        subprocess.run([sys.executable, str(Path(__file__).with_name("smoke_cofdraw_measure.py")), "--base-url", base,
                        "--output-dir", str(output / "measure")], check=True, timeout=600, creationflags=hidden)
        report["measure_report"] = json.loads((output / "measure/measure-smoke-report.json").read_text(encoding="utf-8"))
        assert report["measure_report"]["failure"] is None and not report["measure_report"]["javascript_errors"]
        subprocess.run([sys.executable, str(Path(__file__).with_name("smoke_cofdraw_porosity.py")), "--base-url", base,
                        "--output-dir", str(output / "geometry")], check=True, timeout=600, creationflags=hidden)
        report["geometry_report"] = json.loads((output / "geometry/geometry-smoke-report.json").read_text(encoding="utf-8"))
        assert report["geometry_report"]["failure"] is None and not report["geometry_report"]["javascript_errors"]
        subprocess.run([sys.executable, str(Path(__file__).with_name("smoke_cofdraw_adsorption.py")), "--base-url", base,
                        "--output-dir", str(output / "adsorption")], check=True, timeout=600, creationflags=hidden)
        report["adsorption_report"] = json.loads((output / "adsorption/adsorption-smoke-report.json").read_text(encoding="utf-8"))
        assert report["adsorption_report"]["failure"] is None and not report["adsorption_report"]["javascript_errors"]
        subprocess.run([sys.executable, str(Path(__file__).with_name("smoke_cofdraw_materials.py")), "--base-url", base,
                        "--output-dir", str(output / "materials")], check=True, timeout=600, creationflags=hidden)
        report["materials_report"] = json.loads((output / "materials/materials-smoke-report.json").read_text(encoding="utf-8"))
        subprocess.run([sys.executable, str(Path(__file__).with_name("smoke_cofdraw_describe.py")), "--base-url", base,
                        "--output-dir", str(output / "describe")], check=True, timeout=600, creationflags=hidden)
        report["description_report"] = json.loads((output / "describe/describe-smoke-report.json").read_text(encoding="utf-8"))
        subprocess.run([sys.executable, str(Path(__file__).with_name("smoke_cofdraw_chat.py")), "--base-url", base,
                        "--output-dir", str(output / "chat")], check=True, timeout=600, creationflags=hidden)
        report["chat_report"] = json.loads((output / "chat/chat-smoke-report.json").read_text(encoding="utf-8"))
        assert report["chat_report"]["passed"] and not report["chat_report"]["javascript_errors"]
        native_status = output / "native-window.json"
        with subprocess.Popen([str(executable), "--smoke-window", "--data-dir", str(runtime), "--status-file", str(native_status)], creationflags=hidden) as native:
            try:
                native.wait(timeout=60)
            except subprocess.TimeoutExpired:
                native.terminate(); native.wait(timeout=10); raise
        native_result = json.loads(native_status.read_text(encoding="utf-8"))
        assert native_result["status"] == "desktop_smoke_passed" and native_result["native_save_bridge"]
        report["native_window"] = native_result
    except Exception as exc:
        report["failure"] = f"{type(exc).__name__}: {exc}"
        raise
    finally:
        if server is not None and server.poll() is None:
            server.terminate(); server.wait(timeout=10)
        (output / "release-smoke-report.json").write_text(json.dumps(report, indent=2), encoding="utf-8")
        print(json.dumps({"report": str(output / "release-smoke-report.json"), "failure": report["failure"], "source_checks": len(report["checks"])}))

if __name__ == "__main__":
    main()
