"""Exercise the real 0.4 studio, layout controls, 3D viewer, and analysis API.

Requires development Playwright and installed Chrome/Edge. All working data
and screenshots stay under the explicitly supplied --output-dir.
"""
import argparse
import json
from pathlib import Path
import shutil
import socket
import sys
import threading
import time

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--output-dir',type=Path,required=True);parser.add_argument('--base-url');args=parser.parse_args()
    output=args.output_dir.resolve();output.mkdir(parents=True,exist_ok=True)
    from playwright.sync_api import sync_playwright
    server=worker=None
    if args.base_url:base=args.base_url.rstrip('/')
    else:
        from cof_vision.cofdraw_desktop import create_desktop_app
        import uvicorn
        runtime=output/'runtime-data';shutil.copytree(ROOT/'data',runtime,dirs_exist_ok=True,ignore=shutil.ignore_patterns('exports','*.log'))
        with socket.socket() as sock:sock.bind(('127.0.0.1',0));port=sock.getsockname()[1]
        server=uvicorn.Server(uvicorn.Config(create_desktop_app(runtime),host='127.0.0.1',port=port,log_level='warning'))
        worker=threading.Thread(target=server.run,daemon=True);worker.start()
        for _ in range(200):
            if server.started:break
            time.sleep(.05)
        if not server.started:raise RuntimeError('Test service failed to start')
        base=f'http://127.0.0.1:{port}'
    errors=[];checks=[];failure=None
    try:
        with sync_playwright() as p:
            executable=next((path for path in [Path(r'C:\Program Files\Google\Chrome\Application\chrome.exe'),Path(r'C:\Program Files (x86)\Microsoft\Edge\Application\msedge.exe')] if path.is_file()),None)
            browser=p.chromium.launch(executable_path=str(executable) if executable else None,headless=True)
            # Playwright's polling evaluates predicates; bypass CSP only in this
            # test context. The application keeps its production CSP unchanged.
            page=browser.new_page(viewport={'width':1600,'height':1100},bypass_csp=True);page.set_default_timeout(30000)
            page.on('pageerror',lambda e:errors.append(str(e)))
            page.goto(base);page.wait_for_function("window.COFdrawStudio?.state.online && !document.getElementById('build').disabled")
            page.locator('#mode-manual').click()
            page.locator('#appearance-options').evaluate('e=>e.open=true')
            page.locator('#build').click();page.wait_for_function("window.COFdrawStudio.state.result && !window.COFdrawStudio.state.busy")
            page.screenshot(path=str(output/'studio-standard.png'));checks.append('Default drawing rendered with real backend')
            page.locator('#measure-pore').click();page.locator('#pore-measurement').wait_for(state='visible')
            assert 'Drawing-plane opening' in page.locator('#pore-measurement').inner_text()
            checks.append('Bounded 2D pore opening measurement displayed')
            page.locator('#smiles-b').fill('NNC(=O)c1ccc(C(=O)NN)cc1')
            page.wait_for_function("document.getElementById('linkage').value==='acylhydrazone'")
            assert page.locator('#template').input_value()=='hcb_3_2'
            checks.append('Monomer edit automatically changed linkage and connectivity')
            page.locator('.source-library summary').click()
            page.locator('#source-select').select_option('f-cof-1-source-v1');page.locator('#load-source').click()
            page.wait_for_function("window.COFdrawStudio.state.sourceId==='f-cof-1-source-v1' && !window.COFdrawStudio.state.busy")
            page.locator('[data-style="source"]').click();page.wait_for_function("!window.COFdrawStudio.state.busy")
            page.screenshot(path=str(output/'studio-source-halo.png'));checks.append('Recovered source alignment and halo style rendered')
            page.locator('#edit-geometry').click();handles=page.locator('#geometry-overlay [data-atom]');assert handles.count()>50
            before=page.evaluate('JSON.stringify(window.COFdrawStudio.state.recipe.explicit_graph.atoms[0].xy)')
            handles.first.scroll_into_view_if_needed()
            box=handles.first.bounding_box();page.mouse.move(box['x']+box['width']/2,box['y']+box['height']/2);page.mouse.down();page.mouse.move(box['x']+18,box['y']+12,steps=4);page.mouse.up()
            after=page.evaluate('JSON.stringify(window.COFdrawStudio.state.recipe.explicit_graph.atoms[0].xy)');assert before!=after
            assert page.evaluate('window.COFdrawStudio.state.sourceId') is None
            page.locator('#render-geometry').click();page.wait_for_function('!window.COFdrawStudio.state.busy')
            assert page.evaluate('window.COFdrawStudio.state.result.audit.atoms')==90
            checks.append('Direct atom drag preserves graph and clears source identity')
            page.locator('#undo-geometry').click();page.wait_for_function('!window.COFdrawStudio.state.busy')
            assert page.evaluate('window.COFdrawStudio.state.sourceId')=='f-cof-1-source-v1'
            layered=page.evaluate("window.COFdrawStudio.state.library.styles.find(s=>s.category==='layered').id")
            page.locator(f'[data-style="{layered}"]').click();page.wait_for_function('!window.COFdrawStudio.state.busy')
            assert page.locator('#geometry-overlay [data-atom]').count()==0
            page.locator('[data-style="source"]').click();page.wait_for_function('!window.COFdrawStudio.state.busy')
            checks.append('Undo restores source identity; projected styles disable atom handles')
            page.locator('#edit-geometry').click()
            page.locator('#open-science').click();page.locator('#structure-source').select_option('periodic');page.locator('#crystal-template').select_option('dia-imine');page.locator('#build-crystal').click()
            page.wait_for_function("document.getElementById('science-status').textContent.startsWith('Periodic model generated')")
            page.screenshot(path=str(output/'studio-3d-dia.png'));assert page.locator('#structure-exports').inner_text().find('CIF')>=0
            checks.append('Real periodic DIA model generated and shown with CIF/XYZ/PDB exports')
            with page.expect_download() as download:
                page.locator('#structure-exports button').filter(has_text='CIF').click()
            download.value.save_as(output/'viewer-export.cif')
            page.locator('[data-science-tab="porosity"]').click();page.locator('#pore-samples').select_option('0.9');page.locator('#calculate-pores').click()
            page.wait_for_function("document.getElementById('science-status').textContent.startsWith('Geometric properties complete')")
            assert page.locator('#porosity-results .metric-card').count()>=6
            page.screenshot(path=str(output/'studio-pore-analysis.png'));checks.append('Periodic geometric pore/volume/surface analysis completed')
            page.locator('[data-science-tab="model"]').click();page.locator('#crystal-template').select_option('hcb-imine-aa');page.locator('#build-crystal').click()
            page.wait_for_function("document.getElementById('science-status').textContent.startsWith('Periodic model generated')")
            assert page.locator('#porosity-results .metric-card').count()==0
            checks.append('Changing a periodic model clears previous pore results')
            page.locator('[data-science-tab="bet"]').click()
            rows=[(x,5*80*x/((1-x)*(1+79*x))) for x in [.03,.05,.07,.1,.15,.2,.25,.3,.35]]
            page.locator('#bet-data').fill('\n'.join(f'{x},{n}' for x,n in rows));page.locator('#fit-bet').click()
            page.wait_for_function("document.getElementById('science-status').textContent.startsWith('BET fit calculated')")
            assert '487' in page.locator('#bet-results').inner_text();checks.append('BET fit on synthetic known-capacity fixture recovers expected area')
            page.locator('#close-science').click();page.locator('#smiles-a').fill('O=Cc1cc(C=O)cc(C=O)c1')
            page.locator('#open-science').click();page.locator('[data-science-tab="monomer"]').click();page.locator('#calculate-monomer').click()
            page.wait_for_function("document.getElementById('science-status').textContent.startsWith('Monomer descriptors calculated')")
            assert page.locator('#monomer-results .metric-card').count()>=10
            page.screenshot(path=str(output/'studio-monomer-properties.png'));checks.append('RDKit descriptors and actual optimized monomer conformer displayed')
            page.locator('#close-science').click();page.set_viewport_size({'width':390,'height':844});page.screenshot(path=str(output/'studio-mobile.png'),full_page=True)
            assert page.evaluate('document.documentElement.scrollWidth<=window.innerWidth+1')
            checks.append('390px viewport has no horizontal page overflow')
            assert not errors,errors
            browser.close()
    except Exception as exc:failure=f'{type(exc).__name__}: {exc}';raise
    finally:
        if server:server.should_exit=True;worker.join(timeout=10)
        (output/'studio-smoke-report.json').write_text(json.dumps({'checks':checks,'javascript_errors':errors,'failure':failure},indent=2),encoding='utf-8')
        print(json.dumps({'passed':len(checks),'failure':failure,'javascript_errors':errors}))

if __name__=='__main__':main()
