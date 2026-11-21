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
        calls = []
        page.on("pageerror", lambda error: errors.append(str(error)))
        page.on("request",lambda req: calls.append((req.url.rsplit('/',1)[-1],req.post_data_json))
                if req.method=='POST' and '/api/drawing/science/' in req.url else None)
        page.set_default_timeout(90000)
        def wait_for(expression):
            deadline=time.monotonic()+85
            while not page.evaluate(expression):
                if time.monotonic()>deadline:
                    raise AssertionError((expression,page.locator("#crystal-drawing-progress").inner_text(), errors))
                page.wait_for_timeout(100)
        page.goto(args.base_url)
        wait_for("document.querySelector('#cof-preset option[value=\u0022crystal:cof-300-c5\u0022]') !== null")
        assert page.locator('#cof-preset option[value="crystal:curated-07011N3"]').count()==0
        assert page.locator('#crystal-template option[value="curated-07011N3"]').count()==1
        wait_for("!document.getElementById('cof-preset').disabled")

        def ready():
            wait_for("!document.getElementById('rotate-crystal-design').disabled && document.getElementById('crystal-drawing-image').naturalWidth > 0")
            assert page.locator("#science-dialog").is_hidden()

        def screenshot(name):
            if args.screenshots:
                args.screenshots.mkdir(parents=True, exist_ok=True)
                page.screenshot(path=str(args.screenshots / name))

        presets=page.locator('#cof-preset option[value^="crystal:"]').evaluate_all("(os)=>os.map(o=>o.value)")
        presets=[s for s in presets if s not in ('crystal:custom-dia','crystal:cof-300-c5')]+['crystal:cof-300-c5']
        for preset in presets:
            start=len(calls)
            page.select_option("#cof-preset",preset)
            ready()
            builds=[p for name,p in calls[start:] if name=='build-dia']
            assert len(builds)==1 and not any(name=='model' for name,_ in calls[start:]),(preset,calls[start:])
            assert [builds[0]['monomer_'+b] for b in ('a','b')]==[page.locator('#smiles-'+b).input_value() for b in ('a','b')]
            assert page.evaluate("COFdrawCrystalDesign.report.structure.provenance.kind==='monomer_dia_design'")
            page.click('#rotate-crystal-design')
            assert page.evaluate("COFdrawScience.structure===COFdrawCrystalDesign.report.structure")
            page.click('#close-science')
            assert not any(name=='model' for name,_ in calls[start:]),(preset,calls[start:])
        print(f'{len(presets)} named presets generated from editable SMILES; zero catalog loads during selection or 3D switch.',flush=True)
        assert page.locator("#crystal-drawing-view").input_value()=="Framework drawing"
        assert page.locator('#crystal-drawing-view option',has_text='Literature drawing (comparison)').count()==1
        assert page.locator('#card-a').is_visible() and page.locator('#card-b').is_visible()
        assert page.locator('#preset-a').input_value()=='tapm' and page.locator('#preset-b').input_value()=='ta'
        screenshot("COF-300-drawing.png")
        page.click("#rotate-crystal-design")
        wait_for("COFdrawScience.structure?.provenance?.kind==='monomer_dia_design'")
        assert page.locator("#structure-canvas").is_visible()
        assert page.locator("#crystal-scene").input_value()=="cell"
        before=page.evaluate("COFdrawScience.viewer.theta")
        box=page.locator("#structure-canvas").bounding_box()
        x,y=box["x"]+.5*box["width"],box["y"]+.45*box["height"]
        page.mouse.move(x,y); page.mouse.down(); page.mouse.move(x+60,y+25,steps=5); page.mouse.up()
        assert page.evaluate("COFdrawScience.viewer.theta")!=before
        for scene in ['blockA','blockB','blocknet','block_overlay','cage_overlay']:
            page.select_option('#crystal-scene',scene)
            assert page.evaluate("COFdrawScience.viewer.atoms.length>0 && COFdrawScience.structure===COFdrawCrystalDesign.report.structure"),scene
            assert page.evaluate("COFdrawScience.viewer.model.illustration"),scene
        page.select_option('#crystal-scene','cell')
        for mode in ['spacefill','halfspace','ballstick']:
            page.select_option('#viewer-mode',mode)
            assert page.evaluate("COFdrawScience.viewer.atoms.length>0 && COFdrawScience.structure===COFdrawCrystalDesign.report.structure")
        page.select_option('#crystal-scene','pores')
        wait_for("document.getElementById('map-status').textContent.includes('Mapped across')")
        assert page.locator('#pore-visible').is_checked()
        print('Generated-cell scenes verified: A, B, block net, two overlays, CPK, half space filling, and computed pore surface.',flush=True)
        page.click("#close-science")
        page.click("#compare-crystal-reference")
        wait_for("COFdrawScience.structure?.provenance?.template==='cof-300-c5'")
        assert any(name=='model' for name,_ in calls)
        page.select_option("#crystal-scene","spotlight")
        page.select_option("#crystal-net-index",index=1)
        assert page.evaluate("COFdrawScience.viewer.model.highlightNetwork")==1
        screenshot("COF-300-network.png")
        page.click("#switch-space")
        ready()
        page.select_option("#crystal-linker","design-dia-dmta")
        ready()
        assert page.locator('#preset-b').input_value()=='dmtp'
        assert page.evaluate("COFdrawCrystalDesign.report.structure.provenance.kind==='monomer_dia_design'")
        screenshot("TAPM-DMTA-drawing.png")
        page.click("#rotate-crystal-design")
        assert page.evaluate("COFdrawScience.structure===COFdrawCrystalDesign.report.structure")
        assert page.evaluate("COFdrawScience.viewer.model.periodic")
        screenshot("TAPM-DMTA-cell.png")
        page.click("#close-science")
        page.select_option("#preset-b","ta")
        ready()
        assert page.evaluate("COFdrawCrystalDesign.report.structure.provenance.kind==='monomer_dia_design'")
        assert page.locator("#cof-preset").input_value()=="crystal:custom-dia"
        page.select_option('#reticular-net','lon')
        ready()
        assert page.evaluate("COFdrawCrystalDesign.report.structure.provenance.net==='lon' && COFdrawCrystalDesign.report.cell_validation.junction_coordination_sequence[2]===25")
        assert page.locator('#crystal-drawing-view').input_value()=='Framework drawing'
        page.select_option('#dia-net-count','3')
        ready()
        assert page.evaluate("COFdrawCrystalDesign.report.structure.provenance.network_count===3")
        print('Custom lon drawing and three independent periodic networks verified.',flush=True)
        page.click("#rotate-crystal-design")
        assert page.evaluate("COFdrawScience.structure===COFdrawCrystalDesign.report.structure")
        page.click("#close-science")
        page.fill("#smiles-b","CC")
        wait_for("COFdrawCrystalDesign.report===null && document.getElementById('rotate-crystal-design').disabled")
        wait_for("document.getElementById('crystal-drawing-progress').textContent.length>0 && !document.getElementById('crystal-drawing-progress').textContent.includes('Preparing')")
        assert page.locator("#crystal-drawing-image").is_hidden()
        assert page.evaluate("COFdrawScience.structure===null")
        page.select_option("#cof-preset","crystal:cof-300-c5")
        ready()
        page.select_option("#crystal-linker","curated-13040N3")
        page.select_option("#crystal-linker","cof-300-c5")
        ready()
        assert page.locator('#preset-b').input_value()=='ta'
        assert page.evaluate("COFdrawCrystalDesign.report.structure.provenance.kind==='monomer_dia_design'")
        page.select_option("#crystal-linker","aryl_f")
        wait_for("COFdrawScience.structure?.provenance?.kind==='linker_design'")
        page.click('#close-science')
        ready()
        assert page.evaluate("COFdrawCrystalDesign.report.structure.provenance.kind==='monomer_dia_design'")
        page.click("#leave-crystal-design")
        assert page.locator(".canvas-panel").is_visible() and page.locator("#crystal-drawing-panel").is_hidden()
        assert not errors,errors
        print("Passed default 2D selection, rotation, network focus, linker swap, fast reselection, aryl substitution and editor return.")
        browser.close()


if __name__=="__main__":
    main()
