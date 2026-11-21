"""Exercise the packaged 2D/3D switch, pore layers and scientific exports."""
import argparse
import json
from pathlib import Path
from urllib.parse import urljoin

import numpy as np
from playwright.sync_api import expect, sync_playwright


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--base-url", required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    browsers = [
        Path(r"C:\Program Files\Google\Chrome\Application\chrome.exe"),
        Path(r"C:\Program Files (x86)\Microsoft\Edge\Application\msedge.exe"),
    ]
    installed = next((path for path in browsers if path.is_file()), None)
    with sync_playwright() as playwright:
        options = {"executable_path": str(installed)} if installed else {}
        browser = playwright.chromium.launch(**options)
        page = browser.new_page(viewport={"width": 1450, "height": 960})
        errors = []
        page.on("pageerror", lambda error: errors.append(str(error)))
        page.goto(args.base_url)
        expect(page.locator("main.studio")).to_be_visible()
        expect(page.locator("#science-dialog")).to_be_hidden()
        expect(page.locator("#framework-image")).to_be_visible(timeout=90000)
        expect(page.locator("#auto-chemistry")).to_be_checked()
        expect(page.locator("#linkage")).to_be_hidden()
        page.locator("#auto-chemistry").uncheck()
        expect(page.locator("#linkage")).to_be_visible()
        with page.expect_response("**/api/drawing/render",timeout=90000):
            page.locator("#auto-chemistry").check()

        # Repeated clicks must change the depiction, including when palettes
        # happen to be different on two members of the same style family.
        families={"red_blue":"lines", "hetero_colors":"lines", "standard":"lines",
                  "gradient_stack":"layers", "half_space_stack":"layers"}
        history=["lines"]
        graph=None
        for i in range(6):
            with page.expect_response("**/api/drawing/render",timeout=90000) as rendered:
                page.locator("#random-style").click()
            result=rendered.value.json()
            audit=result["audit"]
            chosen=audit["drawing_style"]["id"]
            assert graph is None or graph==audit["fragment_serialization_sha256"]
            graph=audit["fragment_serialization_sha256"]
            family=families.get(chosen,chosen)
            assert family not in history[-4:], (history,chosen)
            history.append(family)
            png=page.request.get(urljoin(args.base_url,result["downloads"]["png"]))
            assert png.ok
            (args.output_dir/f"style-{i+1}-{chosen}.png").write_bytes(png.body())
        print("Six successive random treatments:",", ".join(history[1:]),flush=True)

        def draw_pair(first, second):
            page.locator("#preset-a").select_option(first)
            with page.expect_response("**/api/drawing/render", timeout=90000) as rendered:
                page.locator("#preset-b").select_option(second)
            assert rendered.value.ok
            with page.expect_response("**/api/drawing/science/build-periodic", timeout=90000) as model:
                page.locator("#open-science" if page.locator("#science-dialog").is_hidden() else "#close-science").click()
            result = model.value.json()
            assert result.get("available"), result.get("reason") or result.get("detail")
            return result

        result = draw_pair("dmtp", "tapb")
        assert len(result["structure"]["atoms"]) == 150
        expect(page.locator("#map-pores")).to_be_enabled(timeout=30000)
        assert page.evaluate("window.COFdrawScience.structure.periodic")
        assert page.locator(":modal").count() == 0
        canvas = page.locator("#structure-canvas")
        panel = page.locator("#render-panel")
        bounds = canvas.bounding_box()
        assert panel.bounding_box()["x"] >= bounds["x"] + bounds["width"]
        page.locator("#render-toggle").click()
        assert not panel.is_visible()
        page.locator("#render-toggle").click()
        expect(page.locator("#pore-visible")).not_to_be_checked()
        page.click("#map-pores")
        expect(page.locator("#pore-visible")).to_be_checked(timeout=120000)
        page.click('[data-science-tab="model"]')
        assert page.evaluate("window.COFdrawPores.triangleCount > 0")
        assert page.evaluate("window.COFdrawPores.tileCount") == 8
        vertices = page.evaluate("""() => {
            const faces=[];
            COFdrawPores.addFaces(faces, ([x,y,z])=>({x,y,z}), 1);
            return faces.flatMap(face=>face.p.map(p=>[p.x,p.y,p.z]));
        }""")
        display_cell = page.evaluate("COFdrawScience.viewer.model.cell")
        fractional = np.array(vertices) @ np.linalg.inv(display_cell)
        assert np.all(fractional >= -1e-8) and np.all(fractional <= 1+1e-8)
        assert np.all(fractional.min(axis=0) < .3) and np.all(fractional.max(axis=0) > .7)
        signature = page.evaluate("COFdrawPores.structureHash")
        assert isinstance(signature, str) and len(signature) == 64
        assert page.locator('[data-science-tab]').evaluate_all("els=>els.slice(0,2).map(e=>e.dataset.scienceTab)") == ["monomer", "model"]
        expect(page.locator('[data-science-tab="model"]')).to_have_attribute('aria-selected','true')
        expect(page.locator('#science-model #science-inspect summary')).to_be_visible()
        page.locator("#display-settings summary").click()
        page.locator("#display-a").select_option("1")
        expect(page.locator("#cell-view-status")).to_contain_text("600 atoms displayed")
        assert page.evaluate("COFdrawPores.tileCount") == 4
        assert page.evaluate("COFdrawPores.structureHash") == signature
        page.locator("#display-a").select_option("2")
        expect(page.locator("#cell-view-status")).to_contain_text("1200 atoms displayed")
        assert page.evaluate("COFdrawPores.tileCount") == 8
        page.locator("#camera-perspective").focus()
        page.keyboard.press("End")
        assert page.evaluate("COFdrawScience.viewer.projected.some(p=>p.factor>1.01)")
        page.locator("#viewer-reset").click()
        page.locator("#display-settings summary").click()
        visible = page.locator("#pore-visible")
        assert visible.is_checked()
        full = canvas.evaluate("e => e.toDataURL()")
        visible.uncheck()
        assert canvas.evaluate("e => e.toDataURL()") != full
        visible.check()
        page.locator("#pore-layer").select_option("connected")
        connected = canvas.evaluate("e => e.toDataURL()")
        page.locator("#pore-cut").evaluate(
            "e => {e.value = '0.6'; e.dispatchEvent(new Event('input', {bubbles:true}))}"
        )
        assert canvas.evaluate("e => e.toDataURL()") != connected
        page.screenshot(path=str(args.output_dir / "channel-section.png"))

        page.locator("#render-panel summary").filter(has_text="External analysis").click()
        with page.expect_download() as download:
            page.locator("#mof-notebook").click()
        notebook_path = args.output_dir / "analysis.ipynb"
        download.value.save_as(notebook_path)
        notebook = json.loads(notebook_path.read_text("utf8"))
        for cell in notebook["cells"]:
            if cell["cell_type"] == "code":
                compile(cell["source"], "notebook", "exec")
        assert "do_periodic_segmentation_rust" in json.dumps(notebook)
        with page.expect_download() as download:
            page.locator("#builder-cif").click()
        cif_path = args.output_dir / "exchange.cif"
        download.value.save_as(cif_path)
        assert "_atom_site_fract_x" in cif_path.read_text("utf8")

        page.locator('[data-science-tab="pxrd"]').click()
        page.locator("#xrd-max").fill("20")
        page.locator("#compare-pxrd-stacks").click()
        expect(page.locator("#pattern-spacing")).to_be_visible(timeout=120000)
        page.locator("#pattern-separation").fill("1")
        page.locator("#zoom-pattern").click()
        page.locator("#pattern-limit").fill("8")
        page.locator("#pattern-limit").press("Enter")
        expect(page.locator("#pattern-peak-list tr").nth(1)).to_be_visible()
        assert page.locator("#pattern-separation").input_value()=="1"
        rows=page.locator("#pattern-peak-list tr").all()[1:]
        assert all(float(row.locator("td").nth(1).inner_text())<=8 for row in rows)
        page.screenshot(path=str(args.output_dir/"low-angle-stacking.png"))

        page.locator("#open-science" if page.locator("#science-dialog").is_hidden() else "#close-science").click()
        assert not panel.is_visible()
        drawer = page.locator("#appearance-options > summary")
        drawer.click()
        side = page.locator(".drawing-render-body").bounding_box()
        sheet = page.locator("#canvas-viewport").bounding_box()
        assert side["x"] >= sheet["x"] + sheet["width"] - 2, (side, sheet)
        drawer.click()
        page.locator("#open-science" if page.locator("#science-dialog").is_hidden() else "#close-science").click()
        assert page.evaluate("COFdrawPores.structureHash") == signature
        page.locator('[data-science-tab="porosity"]').click()
        page.locator("#probe-radius").fill("2.1")
        assert page.evaluate("COFdrawPores.triangleCount") == 0
        assert not visible.is_checked()
        page.locator("#open-science" if page.locator("#science-dialog").is_hidden() else "#close-science").click()
        result=draw_pair("tfppy","tapa")
        expect(page.locator("#viewer-title")).to_contain_text("COF-76")
        assert result["structure"]["provenance"]["kind"]=="reported_crystal"
        expect(page.locator("#pore-visible")).not_to_be_checked()
        page.click("#map-pores")
        expect(page.locator("#pore-visible")).to_be_checked(timeout=120000)
        page.click('[data-science-tab="model"]')
        assert page.evaluate("COFdrawPores.triangleCount > 0")
        assert page.evaluate("COFdrawPores.structureHash")!=signature
        page.screenshot(path=str(args.output_dir/"cof-76-ribbon.png"))
        page.locator("#open-science" if page.locator("#science-dialog").is_hidden() else "#close-science").click()
        result = draw_pair("cv_18fdac58ad343c87", "hth_s")
        assert len(result["structure"]["atoms"]) == 348
        assert result["cell_validation"]["short_contacts"] == 0
        assert "Hth" in page.locator("#viewer-provenance").inner_text()
        page.screenshot(path=str(args.output_dir / "hth-current-drawing.png"))
        page.locator('[data-science-tab="model"]').click()
        page.locator("#structure-source").select_option("periodic")
        page.locator("#crystal-template").select_option("hcb-imine-aa")
        with page.expect_response("**/api/drawing/science/model"):
            page.locator("#build-crystal").click()
        expect(page.locator("#viewer-provenance")).to_contain_text("Independent periodic model")
        page.locator("#open-science" if page.locator("#science-dialog").is_hidden() else "#close-science").click()
        with page.expect_response("**/api/drawing/science/build-periodic"):
            page.locator("#open-science" if page.locator("#science-dialog").is_hidden() else "#close-science").click()
        expect(page.locator("#viewer-provenance")).to_contain_text("Hth")
        page.set_viewport_size({"width": 390, "height": 844})
        assert page.locator("#science-dialog").evaluate("e => e.scrollWidth <= e.clientWidth + 1")
        assert not errors, errors
        browser.close()
        print("Verified pore layers, section, invalidation, docked workspaces and both exports.")


if __name__ == "__main__":
    main()
