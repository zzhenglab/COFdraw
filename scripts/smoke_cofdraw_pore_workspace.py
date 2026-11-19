"""Exercise the packaged 2D/3D switch, pore layers and scientific exports."""
import argparse
import json
from pathlib import Path

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
        page.goto(args.base_url.rstrip("/") + "/?workspace=crystal")
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
        page.locator("#map-pores").click()
        expect(page.locator("#pore-visible")).to_be_checked(timeout=120000)
        assert page.evaluate("window.COFdrawPores.triangleCount > 0")
        signature = page.evaluate("COFdrawPores.structureHash")
        assert isinstance(signature, str) and len(signature) == 64
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

        page.locator("#render-panel details summary").click()
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

        page.locator("#switch-space").click()
        assert not panel.is_visible()
        drawer = page.locator("#appearance-options > summary")
        drawer.click()
        side = page.locator(".drawing-render-body").bounding_box()
        sheet = page.locator("#canvas-viewport").bounding_box()
        assert side["x"] >= sheet["x"] + sheet["width"] - 2, (side, sheet)
        drawer.click()
        page.locator("#switch-space").click()
        assert page.evaluate("COFdrawPores.structureHash") == signature
        page.locator("#probe-radius").fill("2.1")
        assert page.evaluate("COFdrawPores.triangleCount") == 0
        assert not visible.is_checked()
        page.set_viewport_size({"width": 390, "height": 844})
        assert page.locator("#science-dialog").evaluate("e => e.scrollWidth <= e.clientWidth + 1")
        assert not errors, errors
        browser.close()
        print("Verified pore layers, section, invalidation, docked workspaces and both exports.")


if __name__ == "__main__":
    main()
