"""Exercise the 2D → periodic 3D → linker change path in a running desktop app."""
import argparse
import time
from pathlib import Path
from playwright.sync_api import sync_playwright


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--base-url", required=True)
    parser.add_argument("--screenshots", type=Path)
    args = parser.parse_args()
    with sync_playwright() as p:
        chrome = Path(r"C:\Program Files\Google\Chrome\Application\chrome.exe")
        browser = p.chromium.launch(executable_path=str(chrome) if chrome.is_file() else None)
        page = browser.new_page(viewport={"width": 1440, "height": 930})
        errors = []
        page.on("pageerror", lambda error: errors.append(str(error)))
        page.set_default_timeout(90000)
        def wait_for(expression):
            deadline=time.monotonic()+85
            while not page.evaluate(expression):
                if time.monotonic()>deadline:
                    raise AssertionError((expression,page.locator("#crystal-drawing-progress").inner_text(), errors))
                page.wait_for_timeout(100)
        page.goto(args.base_url)
        wait_for("document.querySelectorAll('#cof-preset option[value^=\u0022crystal:\u0022]').length >= 42")
        wait_for("!document.getElementById('cof-preset').disabled")
        page.select_option("#cof-preset", "crystal:cof-300-c5")

        def ready():
            wait_for("!document.getElementById('rotate-crystal-design').disabled && document.getElementById('crystal-drawing-image').naturalWidth > 0")
            assert page.locator("#science-dialog").is_hidden()

        def screenshot(name):
            if args.screenshots:
                args.screenshots.mkdir(parents=True, exist_ok=True)
                page.screenshot(path=str(args.screenshots / name))

        ready()
        assert page.locator("#crystal-drawing-view").input_value()=="Annotated framework"
        screenshot("COF-300-drawing.png")
        page.click("#rotate-crystal-design")
        wait_for("COFdrawScience.structure?.provenance?.template==='cof-300-c5'")
        assert page.locator("#structure-canvas").is_visible()
        assert page.locator("#crystal-scene").input_value()=="cell"
        before=page.evaluate("COFdrawScience.viewer.theta")
        box=page.locator("#structure-canvas").bounding_box()
        x,y=box["x"]+.5*box["width"],box["y"]+.45*box["height"]
        page.mouse.move(x,y); page.mouse.down(); page.mouse.move(x+60,y+25,steps=5); page.mouse.up()
        assert page.evaluate("COFdrawScience.viewer.theta")!=before
        page.select_option("#crystal-scene","spotlight")
        page.select_option("#crystal-net-index",index=1)
        assert page.evaluate("COFdrawScience.viewer.model.highlightNetwork")==1
        screenshot("COF-300-network.png")
        page.click("#switch-space")
        ready()
        page.select_option("#crystal-linker","design-dia-dmta")
        ready()
        assert page.evaluate("COFdrawCrystalDesign.report.structure.provenance.monomer_ids.includes('dmtp')")
        screenshot("TAPM-DMTA-drawing.png")
        page.click("#rotate-crystal-design")
        assert page.evaluate("COFdrawScience.structure===COFdrawCrystalDesign.report.structure")
        assert page.evaluate("COFdrawScience.viewer.model.periodic")
        screenshot("TAPM-DMTA-cell.png")
        page.click("#close-science")
        page.select_option("#crystal-linker","curated-13040N3")
        page.select_option("#crystal-linker","cof-300-c5")
        ready()
        assert page.evaluate("COFdrawCrystalDesign.report.structure.provenance.template==='cof-300-c5'")
        page.select_option("#crystal-linker","aryl_f")
        ready()
        assert page.evaluate("COFdrawCrystalDesign.report.structure.provenance.kind==='linker_design'")
        page.click("#leave-crystal-design")
        assert page.locator(".canvas-panel").is_visible() and page.locator("#crystal-drawing-panel").is_hidden()
        assert not errors,errors
        print("Passed default 2D selection, rotation, network focus, linker swap, fast reselection, aryl substitution and editor return.")
        browser.close()


if __name__=="__main__":
    main()
