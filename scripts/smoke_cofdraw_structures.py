"""Actual-browser structure provenance, stale-response and arrangement checks."""
import argparse
import json
from pathlib import Path
import shutil
import socket
import sys
import threading
import time
import traceback

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--base-url');parser.add_argument('--output-dir',type=Path,required=True);args=parser.parse_args()
    output=args.output_dir.resolve();output.mkdir(parents=True,exist_ok=True);server=worker=None
    if args.base_url:base=args.base_url.rstrip('/')
    else:
        import uvicorn
        from cof_vision.cofdraw_desktop import create_desktop_app
        runtime=output/'runtime-data';shutil.copytree(ROOT/'data',runtime,dirs_exist_ok=True,ignore=shutil.ignore_patterns('exports','*.log'))
        with socket.socket() as sock:sock.bind(('127.0.0.1',0));port=sock.getsockname()[1]
        server=uvicorn.Server(uvicorn.Config(create_desktop_app(runtime),host='127.0.0.1',port=port,log_level='warning'));worker=threading.Thread(target=server.run,daemon=True);worker.start()
        for _ in range(200):
            if server.started:break
            time.sleep(.05)
        assert server.started;base=f'http://127.0.0.1:{port}'
    checks=[];errors=[];failure=None;page=None
    try:
        from playwright.sync_api import sync_playwright
        with sync_playwright() as p:
            executable=next(x for x in [Path(r'C:\Program Files\Google\Chrome\Application\chrome.exe'),Path(r'C:\Program Files (x86)\Microsoft\Edge\Application\msedge.exe')] if x.is_file())
            browser=p.chromium.launch(executable_path=str(executable),headless=True);page=browser.new_page(viewport={'width':1600,'height':1100},bypass_csp=True);page.set_default_timeout(30000)
            page.on('pageerror',lambda e:errors.append(str(e)))
            page.goto(base);page.wait_for_function("window.COFdrawScience && !document.getElementById('build').disabled")
            page.locator('#mode-manual').click()
            def build():
                with page.expect_response('**/api/drawing/render') as response:page.locator('#build').click()
                assert response.value.status==200,response.value.text()
                page.wait_for_function('!window.COFdrawStudio.state.busy')
            build()
            held=[]
            def delay_first_cell(route):
                if not held:held.append((route,route.fetch()))
                else:route.continue_()
            page.route('**/api/drawing/science/build-periodic',delay_first_cell)
            page.locator('#open-science').click()
            page.wait_for_function("document.getElementById('science-status').textContent==='Calculating locally…'")
            for _ in range(50):
                if held:break
                page.wait_for_timeout(100)
            assert held
            assert page.locator('#structure-source').input_value()=='drawing'
            assert page.locator('#science-knowledge').count()==0
            page.locator('#close-science').click()
            page.locator('#preset-b').select_option('pa')
            page.wait_for_function("!document.getElementById('build').disabled && document.getElementById('preset-b').value==='pa'")
            build()
            expected=66
            with page.expect_response('**/api/drawing/science/build-periodic') as fresh:page.locator('#open-science').click()
            assert fresh.value.status==200
            page.wait_for_function('n=>window.COFdrawScience.structure?.atoms.length===n',arg=expected)
            assert 'Pa' in page.locator('#viewer-provenance').inner_text()
            held[0][0].fulfill(response=held[0][1]);page.wait_for_timeout(150)
            assert page.evaluate('window.COFdrawScience.structure.atoms.length')==expected
            assert 'Pa' in page.locator('#viewer-provenance').inner_text()
            page.unroute('**/api/drawing/science/build-periodic',delay_first_cell)
            checks.append('Current drawing opens by default; delayed A response cannot replace the freshly generated B drawing or its provenance')

            page.locator('#finite-view-options > summary').click()
            page.locator('#compare-stacks').click();page.wait_for_function("document.querySelectorAll('#stack-comparison button').length===5")
            initial=page.evaluate('window.COFdrawScience.viewer.atoms.map(a=>[a.x,a.y,a.z])')
            for mode in ['AB','ABC','slipped','twisted']:
                page.locator(f'[data-stacking="{mode}"]').click()
                positions=page.evaluate('window.COFdrawScience.viewer.atoms.map(a=>[a.x,a.y,a.z])')
                assert positions!=initial
                assert page.locator('#patch-stacking').input_value()==mode
            for mode in ['ballstick','spacefill','wireframe','mixed']:
                page.locator('#viewer-mode').select_option(mode)
                for color in ['element','component','layer']:page.locator('#viewer-color').select_option(color)
            assert 'Layer 3' in page.locator('#viewer-legend').inner_text()
            page.screenshot(path=str(output/'stacking-comparison.png'))
            checks.append('Five actual stacking calculations are selectable; four rendering styles and element/block/layer coloring use supplied groups')
            page.locator('#patch-spacing').fill('4.2')
            assert page.evaluate('window.COFdrawScience.viewer.atoms.length')==0
            assert page.locator('#structure-exports button').count()==0
            assert page.locator('#stack-comparison button').count()==0
            with page.expect_response('**/api/drawing/science/patch') as refreshed:page.locator('#build-patch-3d').click()
            assert refreshed.value.json()['parameters']['spacing_angstrom']==4.2
            checks.append('Changing layer parameters clears stale comparison, view and exports until regeneration')

            page.locator('#structure-source').select_option('periodic')
            page.locator('#crystal-template').select_option('dia-imine')
            with page.expect_response('**/api/drawing/science/model'):page.locator('#build-crystal').click()
            page.wait_for_function("window.COFdrawScience.structure?.periodic")
            assert 'Independent periodic model' in page.locator('#viewer-provenance').inner_text()
            original_count=page.evaluate('window.COFdrawScience.viewer.atoms.length')
            page.locator('#network-hypothesis').evaluate('e=>e.open=true')
            with page.expect_response('**/api/drawing/science/networks') as networks:page.locator('#build-networks').click()
            assert networks.value.status==200,networks.value.text()
            page.wait_for_function("window.COFdrawScience.origin==='networks'")
            assert page.evaluate('window.COFdrawScience.viewer.atoms.length')==2*original_count
            assert page.evaluate('window.COFdrawScience.viewer.bonds.every(([a,b])=>window.COFdrawScience.viewer.atoms[a].network===window.COFdrawScience.viewer.atoms[b].network)')
            assert 'not an experimentally assigned interpenetration' in page.locator('#structure-warnings').text_content()
            assert page.locator('#structure-warnings details').count()==1
            assert 'Network copy 2' in page.locator('#viewer-legend').inner_text()
            page.locator('#network-view').select_option('base')
            assert page.evaluate('window.COFdrawScience.viewer.atoms.length')==original_count
            page.locator('#network-view').select_option('hypothesis')
            assert page.evaluate('window.COFdrawScience.viewer.atoms.length')==2*original_count
            page.screenshot(path=str(output/'network-hypothesis.png'))
            checks.append('Periodic source remains separate; translated-copy comparison retains groups and warnings, with no inferred intercopy bonds')

            page.locator('#network-view').select_option('base')
            delayed=[]
            def delay_pores(route):delayed.append((route,route.fetch()))
            page.route('**/api/drawing/science/geometry',delay_pores)
            page.locator('[data-science-tab="porosity"]').click();page.locator('#pore-samples').select_option('0.9');page.locator('#calculate-pores').click()
            for _ in range(50):
                if delayed:break
                page.wait_for_timeout(100)
            assert delayed
            page.locator('[data-science-tab="model"]').click();page.locator('#crystal-template').select_option('hcb-imine-aa')
            with page.expect_response('**/api/drawing/science/model'):page.locator('#build-crystal').click()
            page.wait_for_function("window.COFdrawScience.structure?.atoms.length===66")
            delayed[0][0].fulfill(response=delayed[0][1]);page.wait_for_timeout(150)
            assert page.locator('#porosity-results .metric-card').count()==0
            page.unroute('**/api/drawing/science/geometry',delay_pores)
            checks.append('Changing periodic model cancels pending analysis; late pore results cannot populate the new structure')
            with page.expect_response('**/api/drawing/science/networks') as contact_model:page.locator('#build-networks').click()
            assert contact_model.value.status==200
            page.wait_for_function("window.COFdrawScience.origin==='networks'")
            assert page.locator('#structure-warnings .metric-note.warning').filter(has_text='Short intercopy contacts').is_visible()
            checks.append('Short intercopy contact warning remains visible while routine model notes are collapsed')
            page.set_viewport_size({'width':390,'height':844})
            assert page.locator('#science-dialog').evaluate('e=>e.scrollWidth<=e.clientWidth+1')
            page.screenshot(path=str(output/'structures-mobile.png'))
            checks.append('Structure controls, warnings and selectable views fit a 390px viewport')
            assert not errors,errors;browser.close()
    except Exception:
        failure=traceback.format_exc()
        if page:
            try:failure+='\nUI status: '+page.locator('#science-status').inner_text()
            except Exception:pass
    finally:
        if server:server.should_exit=True;worker.join(timeout=10)
        report={'passed':failure is None and not errors,'checks':checks,'javascript_errors':errors,'failure':failure}
        (output/'structures-smoke-report.json').write_text(json.dumps(report,indent=2),encoding='utf-8');print(json.dumps(report,indent=2))
    if failure:raise SystemExit(1)


if __name__=='__main__':main()
