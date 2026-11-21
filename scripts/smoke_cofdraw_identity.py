"""Real-browser regression for visible COFex records and Auto/manual chemistry."""
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
            executable = next((p for p in (Path(r"C:\Program Files\Google\Chrome\Application\chrome.exe"), Path(r"C:\Program Files (x86)\Microsoft\Edge\Application\msedge.exe")) if p.is_file()), None)
            browser = p.chromium.launch(executable_path=str(executable) if executable else None, headless=True)
            page = browser.new_page(viewport={"width": 1512, "height": 1060}, bypass_csp=True)
            page.set_default_timeout(30000)
            page.on("pageerror", lambda error: errors.append(str(error)))
            page.goto(base)
            page.wait_for_function("window.COFdrawStudio?.state.online && !document.getElementById('build').disabled")
            page.locator("#mode-manual").click()
            assert page.locator("#auto-chemistry").is_checked()
            assert page.locator("#linkage").is_disabled() and page.locator("#template").is_disabled()
            assert "Imine" in page.locator("#chemistry-match").inner_text()
            assert page.locator("#identity-string").is_visible()
            assert page.locator("#identity-string").evaluate("el => !el.closest('details')")
            checks.append("Auto is the real initial mode; resolved choices are disabled and notation is not hidden in details")

            page.locator("#smiles-b").fill("Nc1ccc(N)cc1")
            page.locator("#label-b").fill("Pa")
            page.wait_for_function("window.COFdrawStudio.state.monomers.b?.input==='Nc1ccc(N)cc1' && document.getElementById('template').value==='hcb_3_2'")
            library = page.request.get(base + "/api/drawing/library").json()
            tetratopic = next(r for r in library["monomer_presets"] if r["id"] == "tfppy_f")
            page.locator("#smiles-a").fill(tetratopic["smiles"])
            page.locator("#label-a").fill("F-Tfppy")
            page.wait_for_function("smiles => window.COFdrawStudio.state.monomers.a?.input===smiles && document.getElementById('template').value==='sql_4_2'", arg=tetratopic["smiles"])
            assert page.locator("#linkage").input_value() == "imine"
            assert page.locator('#template option[value="hcb_3_3"]').is_disabled()
            checks.append("Auto reacts to real B and A graph changes: 3+3 to 3+2 to 4+2; incompatible connectivity remains disabled")

            page.locator("#auto-chemistry").uncheck()
            assert page.locator("#linkage").is_enabled() and page.locator("#template").is_enabled()
            reference = next(key for key, row in library["linkages"].items() if row["supported"] is False)
            page.locator("#linkage").select_option(reference)
            assert page.locator("#build").is_disabled()
            checks.append("A reference-only manual linkage stays inspectable but cannot be constructed")
            page.locator("#linkage").select_option("secondary_amine")
            page.locator("#template").select_option("alternating_hexagon")
            page.locator('.validate-monomer[data-block="b"]').click()
            page.wait_for_function("!window.COFdrawStudio.state.busy && document.getElementById('notice').textContent.includes('checked')")
            assert page.locator("#linkage").input_value() == "secondary_amine"
            assert page.locator("#template").input_value() == "alternating_hexagon"
            checks.append("Manual linkage and connectivity survive monomer revalidation")

            page.locator("#auto-chemistry").check()
            assert page.locator("#linkage").input_value() == "imine"
            assert page.locator("#template").input_value() == "sql_4_2"
            assert page.locator("#linkage").is_disabled() and page.locator("#template").is_disabled()
            checks.append("Returning to Auto replaces the previous reduced-product/custom-template choices")

            page.locator("#appearance-options > summary").click()
            page.locator("#image-size").select_option("1000")
            page.locator("#appearance-options > summary").click()
            with page.expect_response(lambda r: r.url.endswith("/api/drawing/prepare") and r.request.method == "POST") as request:
                page.locator("#build").click()
            assert request.value.request.post_data_json["linkage"] == "imine"
            assert request.value.request.post_data_json["template_id"] == "sql_4_2"
            page.wait_for_function("window.COFdrawStudio.state.result && !window.COFdrawStudio.state.dirty && !window.COFdrawStudio.state.busy")
            expected = page.evaluate("window.COFdrawStudio.state.recipe.design_cofex_notation")
            assert page.locator("#identity-string").inner_text() == expected and " COFex.sql.cat0;" in expected
            page.locator('.identity-more > summary').click()
            assert "assumption" in page.locator("#identity-cat").inner_text()
            assert "unknown" in page.locator("#identity-cat").inner_text()
            page.locator("#identity-heading").scroll_into_view_if_needed()
            page.screenshot(path=str(output / "visible-design-record.png"))
            checks.append("Real prepare/render uses Auto choices and displays complete design SMILES/net/cat with explicit assumptions")

            page.evaluate("Object.defineProperty(navigator,'clipboard',{configurable:true,value:{writeText:async text=>{window.__copiedRecord=text;}}})")
            page.locator("#copy-identity").click()
            page.wait_for_function("window.__copiedRecord")
            record = json.loads(page.evaluate("window.__copiedRecord"))
            assert record["design_cofex_notation"] == expected and record["cofex_notation"] is None
            assert record["linkage"] == "imine" and record["selected_template"]["id"] == "sql_4_2"
            assert record["catenation"] == "cat0" and "assumption" in record["catenation_qualification"]
            checks.append("Copy full record includes exact notation, current linkage/template and catenation qualification")

            page.click("#copy-drawing-cofex")
            exchange = page.evaluate("window.__copiedRecord")
            assert "COFex.sql.cat0;COFdraw/1?" in exchange
            page.click("#open-cofex")
            page.fill("#cofex-input",exchange)
            page.click("#cofex-draw")
            page.wait_for_function("!document.querySelector('#cofex-dialog').open && !COFdrawStudio.state.busy")
            assert page.evaluate("COFdrawStudio.state.result.audit.drawing_cofex.text")==exchange
            checks.append("Copy COFex and paste reproduce the editable square design and selected linkage")

            page.select_option("#cof-preset","source:tppa-pore-v1")
            page.wait_for_function("window.COFdrawStudio.state.sourceId==='tppa-pore-v1' && !window.COFdrawStudio.state.busy")
            expected_source = page.evaluate("window.COFdrawStudio.state.sourceMetadata.cofex_notation")
            assert page.locator("#identity-string").inner_text() == expected_source
            assert " COFex.unk.cat0;" in expected_source
            assert "not independently verified" in page.locator("#identity-cat").inner_text()
            assert page.locator("#auto-chemistry").is_hidden()
            page.locator("#identity-heading").scroll_into_view_if_needed()
            page.screenshot(path=str(output / "visible-source-record.png"))
            checks.append("Source record displays its exact full notation, unknown net and unverified catenation without constructor rematching")

            page.select_option("#cof-preset","source:pyurea-source-layout-v1")
            page.wait_for_function("window.COFdrawStudio.state.sourceId==='pyurea-source-layout-v1' && !window.COFdrawStudio.state.busy")
            assert page.locator("#identity-string").inner_text() == "Not assigned"
            assert page.locator("#identity-cat").inner_text().startswith("Unknown")
            assert page.locator("#identity-key").inner_text() == "Not supplied by source"
            page.locator("#copy-identity").click()
            page.wait_for_function("JSON.parse(window.__copiedRecord).source_recipe_id==='pyurea-source-layout-v1'")
            record = json.loads(page.evaluate("window.__copiedRecord"))
            assert record["cofex_notation"] is None and record["cofex_id"] is None and record["catenation"] is None
            checks.append("Unassigned source identity stays unknown in both the visible card and copied JSON")

            page.locator("#smiles-a").fill("O=Cc1ccc(C=O)cc1")
            assert page.locator("#identity-string").inner_text() == "Not assigned"
            assert page.locator("#copy-identity").is_disabled()
            assert page.evaluate("window.COFdrawStudio.state.sourceId") is None
            checks.append("Edited inputs clear active source identity and disable copying a stale structure record")
            page.select_option("#cof-preset","crystal:cof-300-c5")
            page.wait_for_function("Boolean(window.COFdrawCrystalDesign?.report)",timeout=120000)
            page.locator("#crystal-drawing-cofex button").click()
            exchange=page.evaluate("window.__copiedRecord")
            assert "COFex.dia.cat0;COFdraw/1?" in exchange
            page.click("#rotate-crystal-design")
            assert page.locator("#model-cofex .code-value").inner_text()==exchange.split(';')[0]
            page.click("#close-science")
            page.click("#open-cofex")
            page.fill("#cofex-input",exchange)
            page.click("#cofex-draw")
            page.wait_for_function("!document.querySelector('#cofex-dialog').open && COFdrawCrystalDesign.report",timeout=120000)
            assert page.evaluate("COFdrawCrystalDesign.report.structure.provenance.drawing_cofex.text")==exchange
            checks.append("Both crystal views expose the generated pair and network count; pasting rebuilds the periodic model")
            page.click("#open-cofex")
            # The original COF-300 reference used five networks. Opening this
            # record should recover its published phase rather than the editable one-net design.
            page.fill("#cofex-input",exchange.split(" COFex.",1)[0]+" COFex.dia.cat5;COF-300")
            page.click("#cofex-draw")
            page.wait_for_function("!document.querySelector('#cofex-dialog').open && COFdrawCrystalDesign.report",timeout=120000)
            reference=page.evaluate("COFdrawCrystalDesign.report.structure.provenance")
            assert reference["kind"]=="reported_crystal" and reference["drawing_cofex"]["network_count"]==5
            assert page.input_value("#smiles-a") and page.input_value("#smiles-b")
            checks.append("The source fivefold phase reopens with its DOI and editable monomers")
            assert not errors, errors
            browser.close()
    except BaseException:
        failure = traceback.format_exc()
        raise
    finally:
        if server:
            server.should_exit = True
            worker.join(10)
        report = {"passed": failure is None, "checks": checks, "javascript_errors": errors, "failure": failure}
        (output / "identity-browser-report.json").write_text(json.dumps(report, indent=2) + "\n", "utf-8")
        print(json.dumps({"passed": report["passed"], "checks": len(checks), "report": str(output / "identity-browser-report.json")}))


if __name__ == "__main__":
    main()
