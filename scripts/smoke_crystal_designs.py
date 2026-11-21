"""Exercise the 2D → periodic 3D → linker change path in a running desktop app."""
import argparse
import time
from pathlib import Path
from xml.etree import ElementTree
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
                    raise AssertionError((expression,page.locator("#drawing-summary").inner_text(), errors))
                page.wait_for_timeout(100)
        page.goto(args.base_url)
        wait_for("document.querySelector('#cof-preset option[value=\u0022crystal:cof-300-c5\u0022]') !== null")
        assert page.locator('#crystal-template option[value="curated-07011N3"]').count()==1
        wait_for("!document.getElementById('cof-preset').disabled")

        def ready():
            wait_for("window.COFdrawCrystalDesign?.visible && document.getElementById('framework-image').naturalWidth > 0")
            assert page.locator("#science-dialog").is_hidden()

        def screenshot(name):
            if args.screenshots:
                args.screenshots.mkdir(parents=True, exist_ok=True)
                page.screenshot(path=str(args.screenshots / name))

        catalog=page.request.get(args.base_url.rstrip('/')+'/api/drawing/science/library').json()
        entries={row['id']:row for group in ['templates','linker_models','curated_models','reported_models','starter_designs'] for row in catalog.get(group,[])}
        presets=page.locator('#cof-preset option[value^="crystal:"]').evaluate_all("(os)=>os.map(o=>o.value)")
        presets=[s for s in presets if s not in ('crystal:custom-dia','crystal:cof-300-c5')]+['crystal:cof-300-c5']
        generated=references=0
        for preset in presets:
            start=len(calls)
            page.select_option("#cof-preset",preset)
            ready()
            starter=entries[preset.split(':',1)[1]].get('starter',{})
            if starter.get('kind')=='reference':
                assert [(name,p['template_id']) for name,p in calls[start:] if name=='model']==[('model',preset.split(':',1)[1])]
                assert page.evaluate("COFdrawCrystalDesign.report.structure.periodic")
                reaction=page.evaluate("COFdrawCrystalDesign.referenceReaction()")
                assert len(reaction['feeds'])==len(starter['monomer_ids']) and reaction['linkage']==starter['linkage']
                arrow=page.locator('.canvas-panel .scheme-toggle')
                assert arrow.is_enabled()
                if preset=='crystal:curated-07013N3':
                    page.select_option('#extent-options','3')
                    wait_for("COFdrawCrystalDesign.visible?.count===3")
                    screenshot('COF-108-three-connected-circuits.png')
                    arrow.click()
                    wait_for("COFdrawScheme.composed && document.querySelector('.reaction-figure')?.naturalWidth>0")
                    arrow.click()
                references+=1
            else:
                builds=[p for name,p in calls[start:] if name in ('build-dia','build-reticular')]
                assert len(builds)==1 and not any(name=='model' for name,_ in calls[start:]),(preset,calls[start:])
                assert [builds[0]['monomer_'+b] for b in ('a','b')]==[page.locator('#smiles-'+b).input_value() for b in ('a','b')]
                assert page.evaluate("COFdrawCrystalDesign.report.structure.provenance.kind==='monomer_dia_design'")
                assert page.locator('#crystal-drawing-view option',has_text='Chemical drawing + periodic network').count()==1
                generated+=1
            print(preset, 'ready',flush=True)
            page.click('#open-science')
            assert page.evaluate("COFdrawScience.structure===COFdrawCrystalDesign.report.structure")
            assert page.evaluate("COFdrawScience.viewer.atoms.length>=COFdrawScience.structure.atoms.length")
            if starter.get('kind')=='reference':
                page.check('#crystal-topology')
                assert page.evaluate("COFdrawScience.viewer.model.spatial_guides.length>0")
                if preset=='crystal:curated-07013N3':screenshot('COF-108-periodic-topology.png')
                page.uncheck('#crystal-topology')
            page.click('#close-science')
            if starter.get('kind')!='reference':
                assert not any(name=='model' for name,_ in calls[start:]),(preset,calls[start:])
        print(f'{generated} editable designs and {references} reported CIFs on the shared canvas.',flush=True)
        assert page.locator("#crystal-drawing-view").input_value()=="Framework drawing"
        page.click('#zoom-in')
        assert page.locator('#zoom-label').inner_text()=='125%'
        page.click('#zoom-fit')
        assert page.locator('#zoom-label').inner_text()=='100%'
        page.click('#random-style')
        assert page.evaluate("COFdrawStudio.state.style!=='reticular_depth'")
        page.click('#measure-pore')
        assert page.locator('#science-dialog').is_hidden()
        wait_for("COFdrawCrystalDesign.visible?.text.includes('Periodic cavity (LCD)')")
        page.click('#measure-pore')
        page.locator('#appearance-options summary').first.click()
        page.click('#style-list [data-style="reticular_blocks"]')
        wait_for("COFdrawCrystalDesign.visible?.text.includes('data-drawing-style=\\\"reticular_blocks\\\"')")
        with page.expect_download() as faded:
            page.locator('#exports').get_by_role('button',name='SVG',exact=True).click()
        drawing=ElementTree.parse(faded.value.path()).getroot()
        marks=[float(n.attrib['opacity']) for n in drawing.iter() if 'data-depth' in n.attrib]
        assert min(marks)<.25 and max(marks)>.95
        screenshot('COF-300-depth-drawing.png')
        page.click('.canvas-panel .scheme-toggle')
        for n in [1,2,3]:
            page.select_option('#extent-options',str(n))
            wait_for(f"COFdrawCrystalDesign.visible?.count==={n}")
            page.wait_for_selector('#canvas-sheet .reaction-figure')
            with page.expect_download() as saved:
                page.locator('#exports').get_by_role('button',name='SVG',exact=True).click()
            root=ElementTree.parse(saved.value.path()).getroot()
            framework=root.find(".//*[@data-part='framework']")
            assert framework is not None and framework.get('preserveAspectRatio')=='xMidYMid meet'
            assert float(framework.get('width'))>0 and float(framework.get('height'))>0
        if page.locator('#appearance-options').evaluate('el=>el.open'):page.locator('#appearance-options > summary').click()
        screenshot('COF-300-three-circuit-scheme.png')
        page.click('.canvas-panel .scheme-toggle')
        assert page.locator('#card-a').is_visible() and page.locator('#card-b').is_visible()
        assert page.locator('#preset-a').input_value()=='tapm' and page.locator('#preset-b').input_value()=='ta'
        screenshot("COF-300-drawing.png")
        if not page.locator('#appearance-options').evaluate('el=>el.open'):page.click('#appearance-options > summary')
        page.select_option('#crystal-drawing-view','Chemical drawing + periodic network')
        wait_for("COFdrawCrystalDesign.visible?.label==='Chemical drawing + periodic network'")
        page.click('#appearance-options > summary')
        with page.expect_download() as download:
            page.locator('#exports').get_by_role('button',name='SVG',exact=True).click()
        svg=ElementTree.parse(download.value.path()).getroot()
        network=svg.find(".//*[@data-view='periodic-network']")
        assert int(network.attrib['data-source-sites'])*8==int(network.attrib['data-displayed-sites'])
        assert len(network)>16 and svg.find('{http://www.w3.org/2000/svg}svg') is not None
        screenshot('COF-300-drawing-and-network.png')
        page.click("#open-science")
        wait_for("COFdrawScience.structure?.provenance?.kind==='monomer_dia_design'")
        assert page.locator("#structure-canvas").is_visible()
        assert page.locator("#crystal-scene").input_value()=="cell"
        page.check('#crystal-topology')
        assert page.evaluate("COFdrawScience.viewer.model.spatial_guides.length>20")
        screenshot('COF-300-full-topology.png')
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
        page.click("#crystal-interpenetration")
        assert page.locator('#crystal-scene').input_value()=='spotlight'
        page.select_option("#crystal-net-index",index=1)
        assert page.evaluate("COFdrawScience.viewer.model.highlightNetwork")==1
        screenshot("COF-300-network.png")
        page.click("#close-science")
        ready()
        page.select_option("#crystal-linker","design-dia-dmta")
        ready()
        assert page.locator('#preset-b').input_value()=='dmtp'
        assert page.evaluate("COFdrawCrystalDesign.report.structure.provenance.kind==='monomer_dia_design'")
        screenshot("TAPM-DMTA-drawing.png")
        page.click("#open-science")
        assert page.evaluate("COFdrawScience.structure===COFdrawCrystalDesign.report.structure")
        assert page.evaluate("COFdrawScience.viewer.model.periodic")
        screenshot("TAPM-DMTA-cell.png")
        if not page.locator('#display-settings').evaluate('(panel) => panel.open'):
            page.click('#display-settings > summary')
        page.select_option('#display-a','1')
        wait_for("COFdrawScience.viewer.atoms.length===4*COFdrawScience.structure.atoms.length")
        page.select_option('#crystal-scene','blockA')
        page.click('#close-science')
        page.click('#open-science')
        assert page.locator('#display-a').input_value()=='1' and page.locator('#crystal-scene').input_value()=='blockA'
        page.click("#close-science")
        page.select_option("#preset-b","ta")
        ready()
        assert page.evaluate("COFdrawCrystalDesign.report.structure.provenance.kind==='monomer_dia_design'")
        assert page.locator("#cof-preset").input_value()=="crystal:custom-dia"
        page.select_option('#reticular-net','lon')
        ready()
        assert page.evaluate("COFdrawCrystalDesign.report.structure.provenance.net==='lon' && COFdrawCrystalDesign.report.cell_validation.junction_coordination_sequence[2]===25")
        assert page.locator('#crystal-drawing-view').input_value()=='Chemical drawing + periodic network'
        screenshot('lon-drawing-and-network.png')
        page.select_option('#dia-net-count','3')
        ready()
        assert page.evaluate("COFdrawCrystalDesign.report.structure.provenance.network_count===3")
        print('Custom lon drawing and three independent periodic networks verified.',flush=True)
        page.click("#open-science")
        assert page.evaluate("COFdrawScience.structure===COFdrawCrystalDesign.report.structure")
        page.click('#crystal-interpenetration')
        page.select_option('#crystal-build-nets','1')
        wait_for("COFdrawScience.structure?.provenance?.network_count===1 && COFdrawScience.structure===COFdrawCrystalDesign.report?.structure")
        assert page.locator('#crystal-scene').input_value()=='cell'
        page.click('#crystal-interpenetration')
        page.select_option('#crystal-build-nets','3')
        wait_for("COFdrawScience.structure?.provenance?.network_count===3 && COFdrawScience.structure===COFdrawCrystalDesign.report?.structure")
        assert page.locator('#crystal-scene').input_value()=='spotlight'
        page.click("#close-science")
        page.fill("#smiles-b","CC")
        wait_for("COFdrawCrystalDesign.report===null && !window.COFdrawCrystalDesign?.visible")
        wait_for("document.getElementById('drawing-summary').textContent.length>0 && !document.getElementById('drawing-summary').textContent.includes('Preparing')")
        assert page.locator("#framework-image").is_hidden()
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
        assert page.locator(".canvas-panel").is_visible() and not page.evaluate("COFdrawCrystalDesign.active")
        # The optional repeated-cell view can fail independently of a successful
        # construction. Both the diagram and the primitive model remain usable.
        repeat_url='**/api/drawing/science/reports/*/view'
        page.route(repeat_url,lambda route: route.fulfill(status=503,content_type='application/json',body='{"detail":"repeat service unavailable"}'))
        page.select_option('#cof-preset','crystal:cof-300-c5')
        ready()
        assert 'Showing the source cell' in page.locator('#drawing-summary').inner_text()
        assert page.evaluate("COFdrawCrystalDesign.report.structure.provenance.kind==='monomer_dia_design'")
        page.click('#open-science')
        assert page.evaluate("COFdrawScience.viewer.atoms.length===COFdrawScience.structure.atoms.length")
        page.unroute(repeat_url)
        assert not errors,errors
        print("Passed default 2D selection, rotation, network focus, linker swap, fast reselection, aryl substitution and editor return.")
        browser.close()


if __name__=="__main__":
    main()
