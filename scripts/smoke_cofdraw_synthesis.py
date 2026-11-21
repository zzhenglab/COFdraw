"""Real-browser SA Score checks; also accepts a running frozen loopback app."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import shutil
import socket
import sys
import threading
import time
import traceback
from urllib.parse import urlsplit

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--base-url")
    args = parser.parse_args()
    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=True)
    server = worker = None
    if args.base_url:
        parsed = urlsplit(args.base_url)
        if parsed.scheme != "http" or parsed.hostname not in ("127.0.0.1", "localhost", "::1") or parsed.username or parsed.password or parsed.query or parsed.fragment or parsed.path not in ("", "/"):
            parser.error("--base-url must be an HTTP loopback origin")
        base = args.base_url.rstrip("/")
    else:
        runtime = output / "runtime-data"
        shutil.copytree(ROOT / "data", runtime, ignore=shutil.ignore_patterns("exports", "*.log"))
        import uvicorn
        from cof_vision.cofdraw_desktop import create_desktop_app
        with socket.socket() as sock:
            sock.bind(("127.0.0.1", 0))
            port = sock.getsockname()[1]
        server = uvicorn.Server(uvicorn.Config(create_desktop_app(runtime), host="127.0.0.1", port=port, log_level="warning"))
        worker = threading.Thread(target=server.run, daemon=True)
        worker.start()
        deadline = time.monotonic() + 15
        while not server.started and time.monotonic() < deadline:
            time.sleep(.05)
        base = f"http://127.0.0.1:{port}"
    checks, errors = [], []
    failure = None
    try:
        from playwright.sync_api import sync_playwright
        with sync_playwright() as p:
            executable = next((path for path in (Path(r"C:\Program Files\Google\Chrome\Application\chrome.exe"), Path(r"C:\Program Files (x86)\Microsoft\Edge\Application\msedge.exe")) if path.is_file()), None)
            browser = p.chromium.launch(executable_path=str(executable) if executable else None, headless=True)
            page = browser.new_page(viewport={"width": 1512, "height": 1060}, bypass_csp=True)
            page.set_default_timeout(30000)
            page.on("pageerror", lambda error: errors.append(str(error)))
            page.goto(base)
            page.wait_for_function("window.COFdrawStudio?.state.online && !document.getElementById('build').disabled && window.COFdrawStudio.state.monomers.a?.synthetic_accessibility?.available")
            page.locator("#mode-manual").click()
            assert " / 10" in page.locator("#sa-a").inner_text()
            assert " / 10" in page.locator("#sa-b").inner_text()
            assert page.locator("#precursor-sa-card").is_hidden()
            checks.append("Checked finite monomers have SA Score; no drawing aggregate appears before rendering")

            library = page.request.get(base + "/api/drawing/library").json()
            report = page.request.get(base + "/api/drawing/science/synthesis-library").json()
            assert report["count"] == len(library["monomer_presets"])
            assert report["available_count"] == sum(m["synthetic_accessibility"]["available"] for m in library["monomer_presets"])
            page.locator("#open-monomer-library").click()
            assert page.locator("#library-list .library-row").count() == report["count"]
            for row in library["monomer_presets"]:
                score = row["synthetic_accessibility"]
                if score["available"]:
                    assert score["fragment_table_verified"]
                label=f'SA Score {score["score"]:.2f} / 10' if score["available"] else 'SA Score unavailable'
                assert page.locator(f'#library-list .library-row[data-id="{row["id"]}"] [data-sa-score]').inner_text() == label
            page.locator("#library-guide > summary").click()
            assert "not a synthesis route" in page.locator("#library-guide").inner_text()
            page.locator("#library-guide > summary").click()
            page.screenshot(path=str(output / "monomer-score-library.png"))
            page.locator('[aria-label="Close monomer library"]').click()
            checks.append("Every catalog card shows its score or the unavailable state for out-of-scope elements")

            page.locator("#appearance-options > summary").click()
            page.locator("#image-size").select_option("1000")
            page.locator("#appearance-options > summary").click()
            page.locator("#build").click()
            page.wait_for_function("window.COFdrawStudio.state.result && !window.COFdrawStudio.state.dirty && !window.COFdrawStudio.state.busy")
            original = page.evaluate("window.COFdrawStudio.state.result")
            audit = original["audit"]["precursor_synthetic_accessibility"]
            card = page.locator("#precursor-sa-card")
            assert card.is_visible() and audit["complete"] and audit["total_count"] == 2
            assert f'Mean precursor score {audit["mean"]:.2f}' in card.inner_text()
            for precursor in audit["precursors"]:
                assert f'SA Score {precursor["synthetic_accessibility"]["score"]:.2f} / 10' in card.inner_text()
            assert page.request.get(base + original["downloads"]["json"]).json()["audit"]["precursor_synthetic_accessibility"] == audit
            card.scroll_into_view_if_needed()
            page.screenshot(path=str(output / "drawing-precursor-scores.png"))
            checks.append("Fresh drawing shows distinct declared precursors and matching mean/max, preserved in its exported JSON")

            # Exercise real input events: an unvalidated edit must immediately
            # remove both the old single-monomer value and drawing aggregate.
            new_smiles = "O=Cc1ccc(C=O)cc1"
            page.locator("#smiles-a").fill(new_smiles)
            assert card.is_hidden() and card.inner_text() == ""
            assert "updates after" in page.locator("#sa-a").inner_text()
            page.wait_for_function("smiles => window.COFdrawStudio.state.monomers.a?.input===smiles", arg=new_smiles)
            assert " / 10" in page.locator("#sa-a").inner_text()
            assert card.is_hidden()
            checks.append("Changing a precursor clears stale aggregate and stale monomer score immediately; revalidation alone does not resurrect the old drawing score")
            page.locator("#label-a").fill("TA")
            page.locator("#build").click()
            page.wait_for_function("old => window.COFdrawStudio.state.result?.drawing_id!==old && !window.COFdrawStudio.state.dirty && !window.COFdrawStudio.state.busy", arg=original["drawing_id"])
            fresh = page.evaluate("window.COFdrawStudio.state.result.audit.precursor_synthetic_accessibility")
            assert card.is_visible() and fresh["mean"] != audit["mean"]
            assert new_smiles in {row["smiles"] for row in fresh["precursors"]}
            assert f'Mean precursor score {fresh["mean"]:.2f}' in card.inner_text()
            checks.append("Rebuilding recomputes scores from the changed finite precursor graphs")

            # Editing raw recipe JSON invalidates the visible aggregate even
            # while the last rendered result remains available in memory.
            page.locator("#recipe-json").evaluate("el=>{for(let p=el.parentElement;p;p=p.parentElement)if(p.tagName==='DETAILS')p.open=true;}")
            page.locator("#recipe-json").fill(page.locator("#recipe-json").input_value() + "\n")
            assert page.evaluate("window.COFdrawStudio.state.jsonPending")
            assert card.is_hidden() and card.inner_text() == ""
            checks.append("Unapplied recipe JSON clears the old aggregate instead of displaying a stale assessment")

            page.select_option("#cof-preset","source:ucof-source-layout-v1")
            page.wait_for_function("""COFdrawStudio.state.result?.audit.original_recipe_id==='ucof-source-layout-v1'
                && !COFdrawStudio.state.busy && !COFdrawStudio.state.dirty""")
            source = page.evaluate("window.COFdrawStudio.state.result.audit.precursor_synthetic_accessibility")
            assert card.is_visible() and source["complete"] and source["total_count"] == 2
            assert "not a synthesis route or COF reaction/yield prediction" in card.inner_text()
            checks.append("Explicit source drawings retain a qualified finite-precursor summary")
            assert not errors, errors
            browser.close()
    except BaseException:
        failure = traceback.format_exc()
        raise
    finally:
        if server:
            server.should_exit = True
            worker.join(10)
        result = {"passed": failure is None, "checks": checks, "javascript_errors": errors, "failure": failure}
        (output / "synthesis-browser-report.json").write_text(json.dumps(result, indent=2) + "\n", "utf-8")
        print(json.dumps({"passed": result["passed"], "checks": len(checks), "report": str(output / "synthesis-browser-report.json")}))


if __name__ == "__main__":
    main()
