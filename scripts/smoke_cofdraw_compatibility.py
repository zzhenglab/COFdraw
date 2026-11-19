"""Exercise real monomer selection against the local chemistry engine in Chromium.

Requires development Playwright and Chrome/Edge. Copies runtime data into the
explicit output directory; never writes to the source catalog or installed app.
"""
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
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--base-url", help="Existing loopback service, including the frozen release; never stopped by this script")
    parser.add_argument("--expected-monomers", type=int)
    args = parser.parse_args()
    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=True)
    from playwright.sync_api import sync_playwright
    server = worker = None
    if args.base_url:
        parsed = urlsplit(args.base_url)
        if parsed.scheme != "http" or parsed.hostname not in ("127.0.0.1", "localhost", "::1") or parsed.username or parsed.password or parsed.query or parsed.fragment or parsed.path not in ("", "/"):
            parser.error("--base-url must be an HTTP loopback origin without credentials, path, query, or fragment")
        base = args.base_url.rstrip("/")
    else:
        runtime = output / "runtime-data"
        if runtime.exists():
            raise ValueError("Use a fresh output directory for this browser regression.")
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
    checks, errors, api_errors = [], [], []
    failure = None
    browser = page = None
    try:
        assert server is None or server.started, "Local test server did not start"
        with sync_playwright() as p:
            executable = next((path for path in [Path(r"C:\Program Files\Google\Chrome\Application\chrome.exe"), Path(r"C:\Program Files (x86)\Microsoft\Edge\Application\msedge.exe")] if path.is_file()), None)
            browser = p.chromium.launch(executable_path=str(executable) if executable else None, headless=True)
            page = browser.new_page(viewport={"width": 1512, "height": 1060}, bypass_csp=True)
            page.set_default_timeout(30000)
            page.on("pageerror", lambda error: errors.append(str(error)))
            page.on("response", lambda response: api_errors.append({"status": response.status, "url": response.url}) if "/api/drawing/" in response.url and response.status >= 400 else None)
            page.goto(base)
            page.wait_for_function("window.COFdrawStudio?.state.online && !document.getElementById('build').disabled")
            page.locator("#mode-manual").click()
            library = page.request.get(f"{base}/api/drawing/library").json()
            monomers = {item["id"]: item for item in library["monomer_presets"]}
            if args.expected_monomers is not None:
                assert len(monomers) == args.expected_monomers, (len(monomers), args.expected_monomers)

            def ready():
                page.wait_for_function("!document.getElementById('build').disabled")

            def inspected(block):
                page.wait_for_function("b => window.COFdrawStudio.state.monomers[b]?.input === document.getElementById('smiles-'+b).value.trim()", arg=block)

            def open_b():
                page.locator('.library-open[data-block="b"]').click()

            def choose(item):
                page.locator(f'.library-row[data-id="{item}"] .library-use').click()

            def close_library():
                page.keyboard.press("Escape")
                page.locator("#monomer-library-dialog").wait_for(state="hidden")

            def ids():
                return page.locator(".library-row").evaluate_all("rows => rows.map(row=>row.dataset.id)")

            search=page.locator("#combination-search")
            starts=page.locator("#cof-preset")
            assert starts.locator('option[value^="pair:"]').count()>len(library["cof_presets"])
            search.fill("lzu1")
            choices=page.locator("#combination-suggestions")
            assert choices.is_visible() and "COF-LZU1" in choices.inner_text()
            search.press("ArrowDown")
            assert search.get_attribute("aria-activedescendant")
            with page.expect_response("**/api/drawing/render",timeout=90000):
                search.press("Enter")
            assert "COF-LZU1" in starts.locator("option:checked").inner_text()
            assert not choices.is_visible()
            search.fill("TpPa")
            assert "Tp + Pa" in starts.locator('option[value="tppa-pore-v1"]').inner_text()
            starts.select_option("tppa-pore-v1")
            ready()
            assert page.locator("#linkage").input_value()=="beta_ketoenamine"
            search.fill("COF300")
            assert starts.locator('option[value="source:cof-300-projection-v1"]').count()==1
            search.fill("Pyurea")
            starts.select_option("source:pyurea-source-layout-v1")
            page.wait_for_function("window.COFdrawStudio.state.sourceId==='pyurea-source-layout-v1' && !window.COFdrawStudio.state.busy")
            assert starts.input_value()=="source:pyurea-source-layout-v1"
            search.fill("nAnOtExTnOmAtCh")
            assert "No matches" in page.locator("#combination-count").inner_text()
            assert starts.input_value()=="source:pyurea-source-layout-v1"
            search.fill("imine")
            candidate=starts.locator('option[value^="pair:"]').first.get_attribute("value")
            starts.select_option(candidate)
            ready()
            assert page.evaluate("window.COFdrawStudio.state.sourceId===null")
            with page.expect_response("**/api/drawing/prepare") as prepared:
                page.locator("#build").click()
            assert prepared.value.ok
            page.wait_for_function("!window.COFdrawStudio.state.busy")
            search.fill("")
            starts.select_option(library["cof_presets"][0]["id"])
            ready()
            checks.append("Search and build matched pairs across linkages, open a named recorded COF, and return to a new design")

            # The old automatic fallback incorrectly allowed TFB to partner TFB.
            page.locator("#open-monomer-library").click()
            assert page.locator("#library-current-a").inner_text() == "TFB"
            assert page.locator("#library-current-b").inner_text() == "TAPB"
            assert not page.locator("#library-guide").evaluate("element => element.open")
            assert str(len(monomers)) in page.locator("#library-count").inner_text()
            page.locator("#library-step-b").focus()
            page.keyboard.press("Enter")
            assert page.locator("#library-step-b").get_attribute("aria-pressed") == "true"
            assert "tfb" not in ids() and "tapb" in ids()
            page.locator("#library-step-a").focus()
            page.keyboard.press("Enter")
            assert page.locator("#library-step-a").get_attribute("aria-pressed") == "true"
            assert len(ids()) == len(monomers)
            checks.append("Current A/B pair and dynamic catalog count are visible; keyboard selection moves between the guided steps")
            choose("tfb")
            page.wait_for_function("document.getElementById('library-target').value==='b' && !document.getElementById('library-target').disabled")
            assert page.locator("#library-step-b").get_attribute("aria-pressed") == "true"
            assert page.locator("#library-current-a").inner_text() == "TFB"
            assert "tfb" not in ids() and "tapb" in ids() and "deth" in ids()
            assert page.locator(".library-use:disabled").count() == 0
            assert page.locator('#preset-b option[value="tfb"]').count() == 0
            page.screenshot(path=str(output / "compatible-partners.png"))
            checks.append("Selecting TFB as A defaults to real complementary partners in both the grid and B dropdown")

            page.locator("#library-show-all").check()
            assert len(ids()) == len(monomers)
            references = [item for item in monomers.values() if item.get("catalog_use") == "reference_only"]
            assert references, "The bundled catalog includes qualified reference-only monomers"
            reference = references[0]
            ref_card = page.locator(f'.library-row[data-id="{reference["id"]}"]')
            assert "Reference only" in ref_card.locator(".library-card-status").inner_text()
            assert ref_card.locator(".library-use").is_disabled()
            ref_card.locator(".library-card-details").click()
            assert page.locator("#library-edit").is_enabled()
            assert page.locator("#library-apply").is_disabled()
            assert not page.locator("#library-evidence-panel").evaluate("element => element.open")
            page.wait_for_function("document.activeElement === document.querySelector('#library-detail-panel > summary')")
            page.locator("#library-evidence-panel > summary").click()
            assert "not ground truth" in page.locator("#library-audit-qualifier").inner_text()
            checks.append("Reference-only entries keep disabled construction, working inspection/editor access, and expandable qualified source evidence")
            bound = next((item for item in monomers.values() if item.get("master_binding")), None)
            if bound:
                page.locator(f'.library-row[data-id="{bound["id"]}"] .library-card-details').click()
                page.locator("#library-evidence-panel > summary").click()
                page.locator("#library-evidence-records details > summary").click()
                audit_text = page.locator("#library-evidence-records").inner_text()
                assert "not a human ground-truth admission" in audit_text
                assert bound["master_binding"]["master_sha256"] in audit_text
                assert all(str(value) in audit_text for value in bound["master_binding"]["record_ids"])
                assert "not ground truth" in page.locator("#library-audit-qualifier").inner_text()
                page.screenshot(path=str(output / "source-review-details.png"))
                checks.append("Master-key evidence displays actual record IDs and hashes with its non-ground-truth qualification")
            assert page.locator('[data-id="tfb"] .library-use').is_disabled()
            assert page.locator('[data-id="tfb"] .library-card-status').bounding_box()["height"] >= 12
            page.locator('[data-id="tfb"] .library-card-details').click()
            assert page.locator("#library-apply").is_disabled()
            assert "No implemented linkage" in page.locator("#library-compatibility").inner_text()
            assert page.locator("#library-edit").is_enabled()
            page.locator("#library-show-all").uncheck()
            assert page.locator("#library-detail-panel").is_hidden()
            page.locator("#library-search").fill("TFB")
            assert page.locator(".library-row").count() == 0
            assert "No compatible partners" in page.locator("#library-list").inner_text()
            page.locator("#library-search").fill("")
            page.set_viewport_size({"width": 390, "height": 844})
            page.wait_for_function("document.querySelector('.library-row')?.dataset.previewState === 'ready'")
            assert page.locator("#monomer-library-dialog").evaluate("element => element.scrollWidth <= element.clientWidth + 1")
            for selector in ("#library-step-a", "#library-step-b", "#library-search", "#library-sites", "#library-show-all"):
                rect = page.locator(selector).bounding_box()
                assert rect and 0 <= rect["x"] and rect["x"] + rect["width"] <= 391
                assert 0 <= rect["y"] and rect["y"] + rect["height"] <= 844
            page.screenshot(path=str(output / "guided-partners-390.png"))
            page.locator("#library-search").fill("DETH")
            assert page.locator('.library-row[data-id="deth"] [data-sa-score]').count() == 1
            page.wait_for_function("document.querySelector('.library-row[data-id=deth]')?.dataset.previewState === 'ready'")
            choose("deth")
            page.locator("#monomer-library-dialog").wait_for(state="hidden")
            ready()
            page.set_viewport_size({"width": 1512, "height": 1060})
            assert page.locator("#linkage").input_value() == "acylhydrazone"
            assert page.locator("#template").input_value() == "hcb_3_2"
            checks.append("Show-all preserves inspection with disabled reasons; choosing hydrazide updates actual linkage and 3+2 template")
            checks.append("At 390px the pair steps, filters and compatibility control fit; a real preview/SA card selects a complementary partner")

            with page.expect_response("**/api/drawing/render") as rendered:
                with page.expect_response("**/api/drawing/prepare") as prepared:
                    page.locator("#build").click()
                assert prepared.value.status == 200, prepared.value.text()
                assert prepared.value.request.post_data_json["linkage"] == "acylhydrazone"
            assert rendered.value.status == 200, rendered.value.text()
            page.wait_for_function("!window.COFdrawStudio.state.busy")
            checks.append("Selected compatible pair completes real prepare and chemical rendering")

            # Explicit UI-only fixture: keep the real rendered chemistry/exports,
            # vary the reported quality to exercise warning and clearing states.
            quality = {"value": "needs_review"}
            def quality_fixture(route):
                response = route.fetch()
                payload = response.json()
                payload["audit"]["layout_refinement"] = {"method": "bond_angle_pose_v2", "quality": quality["value"]}
                route.fulfill(response=response, json=payload)

            page.route("**/api/drawing/render", quality_fixture)
            page.locator("#appearance-options").evaluate("element => element.open=true")
            with page.expect_response("**/api/drawing/render"):
                page.locator('[data-style="red_blue"]').click()
            page.wait_for_function("!window.COFdrawStudio.state.busy")
            assert "linkage junctions still need geometry review" in page.locator("#notice").inner_text()
            assert "Move / adjust" in page.locator("#geometry-message").inner_text()
            quality["value"] = "regular"
            with page.expect_response("**/api/drawing/render"):
                page.locator('[data-style="standard"]').click()
            page.wait_for_function("!window.COFdrawStudio.state.busy")
            assert "need geometry review" not in page.locator("#notice").inner_text()
            assert "need geometry review" not in page.locator("#geometry-message").inner_text()
            page.unroute("**/api/drawing/render", quality_fixture)
            checks.append("Explicit audit-quality fixture shows linkage-review guidance and clears it after a regular result")

            page.locator("#auto-chemistry").uncheck()
            page.locator("#linkage").select_option("imine")
            assert page.locator("#build").is_disabled()
            open_b()
            assert "deth" not in ids() and "pa" in ids()
            page.locator("#library-show-all").check()
            assert page.locator('[data-id="deth"] .library-use').is_disabled()
            choose("pa")
            page.locator("#monomer-library-dialog").wait_for(state="hidden")
            ready()
            assert page.locator("#linkage").input_value() == "imine"
            page.locator("#preset-a").select_option("tfpm")
            inspected("a")
            ready()
            assert page.locator("#template").input_value() == "alternating_hexagon"
            assert page.locator('#template option[value="sql_4_2"]').is_disabled()
            checks.append("Fixed linkage invalidates prior hydrazide B; tetrahedral A invalidates the planar template without assigning a false net")

            page.locator("#auto-chemistry").check()
            page.locator("#preset-a").select_option("tapb")
            inspected("a")
            assert page.locator("#build").is_disabled()
            assert page.locator("#build-from-canvas").is_disabled()
            assert page.locator('#preset-b option[value="pa"]').is_disabled()
            open_b()
            assert "tapb" not in ids() and "pa" not in ids() and "tfb" in ids()
            choose("tfb")
            page.locator("#monomer-library-dialog").wait_for(state="hidden")
            ready()
            assert page.locator("#smiles-a").input_value() == monomers["tapb"]["smiles"]
            assert page.locator("#smiles-b").input_value() == monomers["tfb"]["smiles"]
            assert page.locator("#template").input_value() != "sql_4_2"
            with page.expect_response("**/api/drawing/render") as drawn:
                with page.expect_response("**/api/drawing/prepare") as prepared:
                    page.locator("#build").click()
                assert prepared.value.status == 200
                assert prepared.value.request.post_data_json['monomer_a'] == monomers['tfb']['smiles']
                assert prepared.value.request.post_data_json['monomer_b'] == monomers['tapb']['smiles']
            assert drawn.value.status == 200
            page.wait_for_function("!window.COFdrawStudio.state.busy")
            assert page.locator("#smiles-a").input_value() == monomers["tapb"]["smiles"]
            assert page.locator("#smiles-b").input_value() == monomers["tfb"]["smiles"]
            checks.append("Choosing aldehyde B preserves the user's amine A; real preparation maps reactive roles without moving the visible selections")

            # In-flight and unsupported custom A must not use stale previous metadata.
            page.locator("#smiles-a").fill("c1ccccc1")
            assert page.locator("#build").is_disabled()
            assert page.locator("#build-from-canvas").is_disabled()
            open_b()
            assert page.locator(".library-row").count() == 0
            inspected("a")
            assert page.locator(".library-row").count() == 0
            page.locator("#library-show-all").check()
            assert page.locator(".library-use:enabled").count() == 0
            close_library()
            assert page.locator("#build").is_disabled()
            page.locator("#smiles-a").fill("")
            open_b()
            assert len(ids()) == len(monomers)
            assert page.locator("#library-show-all-field").is_hidden()
            assert page.locator("#build").is_disabled()
            checks.append("Pending/unsupported custom A cannot reuse stale compatibility; an empty A permits full catalog browsing")
            close_library()
            page.locator("#smiles-a").fill(monomers["tapb"]["smiles"])
            inspected("a")
            ready()
            open_b()
            assert not page.locator("#library-show-all").is_checked()
            assert "tapb" not in ids() and "tfb" in ids()
            assert page.locator('.library-thumbnail img').first.evaluate("element => getComputedStyle(element).transform") == 'none'
            checks.append("Restoring valid A refreshes B choices and resets the compatible-only default")
            close_library()
            quality["value"] = "needs_review"
            page.route("**/api/drawing/render", quality_fixture)
            page.locator(".source-library").evaluate("element => element.open=true")
            page.locator("#source-select").select_option("f-cof-1-source-v1")
            with page.expect_response("**/api/drawing/render"):
                page.locator("#load-source").click()
            page.wait_for_function("!window.COFdrawStudio.state.busy && window.COFdrawStudio.state.sourceId==='f-cof-1-source-v1'")
            assert "need geometry review" not in page.locator("#notice").inner_text()
            assert "need geometry review" not in page.locator("#geometry-message").inner_text()
            page.unroute("**/api/drawing/render", quality_fixture)
            checks.append("Preserved source drawing does not acquire an automatic-design alignment warning")
            assert not errors, errors
            assert not api_errors, api_errors
            browser.close()
    except Exception:
        failure = traceback.format_exc()
        if page:
            try:
                page.screenshot(path=str(output / "failure.png"), full_page=True)
            except Exception:
                pass
    finally:
        if server is not None:
            server.should_exit = True
            worker.join(timeout=10)
        report = {"passed": failure is None and not errors and not api_errors, "checks": checks, "browser_errors": errors, "api_errors": api_errors, "failure": failure, "base_url": base, "server_ownership": "existing_service" if args.base_url else "owned_python_thread", "expected_monomers": args.expected_monomers}
        (output / "compatibility-report.json").write_text(json.dumps(report, indent=2), encoding="utf-8")
        print(json.dumps(report, indent=2))
    if not report["passed"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
