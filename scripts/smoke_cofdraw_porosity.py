"""Real Level-1 porosity UI/API/export regression, including unresolved channels."""
import argparse
import json
from pathlib import Path
import shutil
import socket
import sys
import threading
import time

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'src'))


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--output-dir',type=Path,required=True)
    parser.add_argument('--base-url')
    args=parser.parse_args(); output=args.output_dir.resolve(); output.mkdir(parents=True,exist_ok=True)
    from playwright.sync_api import sync_playwright
    server=worker=None
    if args.base_url: base=args.base_url.rstrip('/')
    else:
        from cof_vision.cofdraw_desktop import create_desktop_app
        import uvicorn
        runtime=output/'runtime-data'
        shutil.copytree(ROOT/'data',runtime,dirs_exist_ok=True,ignore=shutil.ignore_patterns('exports','*.log'))
        with socket.socket() as sock:sock.bind(('127.0.0.1',0));port=sock.getsockname()[1]
        server=uvicorn.Server(uvicorn.Config(create_desktop_app(runtime),host='127.0.0.1',port=port,log_level='warning'))
        worker=threading.Thread(target=server.run,daemon=True);worker.start()
        for _ in range(200):
            if server.started:break
            time.sleep(.05)
        if not server.started:raise RuntimeError('Local test service failed to start')
        base=f'http://127.0.0.1:{port}'
    checks=[];errors=[];failure=None
    try:
        with sync_playwright() as p:
            executable=next((path for path in [Path(r'C:\Program Files\Google\Chrome\Application\chrome.exe'),
                Path(r'C:\Program Files (x86)\Microsoft\Edge\Application\msedge.exe')] if path.is_file()),None)
            browser=p.chromium.launch(executable_path=str(executable) if executable else None,headless=True)
            page=browser.new_page(viewport={'width':1440,'height':1100},bypass_csp=True);page.set_default_timeout(45000)
            page.on('pageerror',lambda error:errors.append(str(error)))
            requests=[]
            page.on('request',lambda request:requests.append(request.post_data_json) if request.url.endswith('/science/geometry') else None)
            page.goto(base);page.wait_for_function('window.COFdrawStudio?.state.online && !document.getElementById("build").disabled')
            page.evaluate("window.addEventListener('cofdraw:geometry',event=>window.geometryResult=event.detail)")
            page.locator('#mode-manual').click();page.locator('#build').click()
            page.wait_for_function('COFdrawStudio.state.result && !COFdrawStudio.state.busy')
            page.locator('#open-science').click()
            page.locator('#structure-source').select_option('periodic')
            page.locator('#crystal-template').select_option('hcb-imine-aa')
            page.locator('#build-crystal').click()
            page.wait_for_function('COFdrawScience.structure?.periodic===true')
            page.locator('[data-science-tab="porosity"]').click()
            assert page.locator('#pore-samples').input_value()=='0.6'
            page.locator('#calculate-pores').click()
            page.wait_for_function('window.geometryResult?.operation==="geometry"')
            result=page.evaluate('window.geometryResult')
            assert requests[-1]['grid_spacing']==.6 and requests[-1]['surface_samples']==64
            assert result['engine']=='cofdraw-grid-v1'
            text=page.locator('#porosity-results').inner_text()
            for required in ['PLD','LCD','geometric bounds','Accessible area','Accessible probe-center volume','Density','void fraction','Local exposed area','Periodic grid']:
                assert required in text, required
            assert result['connectivity']['status']=='certified_channel'
            checks.append('Default Level 1 uses deterministic geometry endpoint and displays bounds, accessibility, density and effective grid')
            page.screenshot(path=str(output/'geometry-level1.png'))
            with page.expect_download() as download:
                page.locator('#porosity-results button').filter(has_text='Save geometric-properties report').click()
            report_path=output/'geometry-report.json';download.value.save_as(report_path)
            exported=json.loads(report_path.read_text('utf8'))
            assert exported==result and exported['inputs']['grid_spacing']==.6
            reread=page.request.get(base+'/api/drawing/science/reports/'+result['report_id']).json()
            assert reread==exported
            checks.append('Saved geometry JSON and report readback retain exact inputs, structure hash and method')
            page.locator('#probe-radius').fill('1.0')
            assert page.locator('#porosity-results .metric-card').count()==0
            page.locator('#pore-refinement').check()
            page.evaluate('window.geometryResult=null');page.locator('#calculate-pores').click()
            page.wait_for_function('window.geometryResult!==null')
            assert page.evaluate('geometryResult.convergence.requested') is True
            assert 'Refined PLD' in page.locator('#porosity-results').inner_text()
            checks.append('Probe/refinement edits clear old results and the finer-grid comparison is shown')
            # A periodic FCC fixture has an open diagonal channel which the
            # conservative coarse axis grid cannot certify at radius 1 A.
            cif=output/'fcc-geometry-fixture.cif'
            cif.write_text('''data_fcc_numerical_fixture
_cell_length_a 8
_cell_length_b 8
_cell_length_c 8
_cell_angle_alpha 90
_cell_angle_beta 90
_cell_angle_gamma 90
_space_group_name_H-M_alt 'P 1'
loop_
_space_group_symop_operation_xyz
'x,y,z'
loop_
_atom_site_label
_atom_site_type_symbol
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
_atom_site_occupancy
C1 C 0 0 0 1
C2 C 0 0.5 0.5 1
C3 C 0.5 0 0.5 1
C4 C 0.5 0.5 0 1
''',encoding='utf8')
            page.locator('[data-science-tab="model"]').click();page.locator('#cif-file').set_input_files(cif)
            page.wait_for_function('COFdrawScience.structure?.atoms.length===4')
            page.locator('[data-science-tab="porosity"]').click();page.locator('#pore-refinement').uncheck()
            page.evaluate('window.geometryResult=null');page.locator('#calculate-pores').click()
            page.wait_for_function('window.geometryResult!==null')
            assert page.evaluate('geometryResult.connectivity.status')=='unresolved_at_this_resolution'
            values=page.locator('#porosity-results .metric-card strong').all_text_contents()
            assert values.count('Unresolved')>=2
            assert 'not evidence of a closed pore' in page.locator('#porosity-results').inner_text()
            page.screenshot(path=str(output/'geometry-unresolved-channel.png'))
            checks.append('An actual unresolved diagonal channel is labeled Unresolved, not falsely closed or zero-accessible')
            def delayed(route):
                response=route.fetch()
                page.locator('#probe-radius').fill('1.1')
                route.fulfill(response=response)
            page.route('**/api/drawing/science/geometry',delayed)
            page.evaluate('window.geometryResult=null');page.locator('#calculate-pores').click()
            for _ in range(100):
                if page.locator('#probe-radius').input_value()=='1.1':break
                page.wait_for_timeout(100)
            page.wait_for_timeout(250)
            assert page.locator('#porosity-results .metric-card').count()==0
            assert page.evaluate('window.geometryResult===null')
            page.unroute('**/api/drawing/science/geometry',delayed)
            checks.append('A late real response cannot restore properties calculated with superseded settings')
            count=len(requests)
            page.locator('[data-science-tab="model"]').click();page.locator('#structure-source').select_option('drawing')
            page.locator('[data-science-tab="porosity"]').click();page.locator('#calculate-pores').click()
            page.wait_for_function('document.getElementById("science-status").classList.contains("error")')
            assert len(requests)==count
            assert 'periodic' in page.locator('#science-status').inner_text()
            checks.append('Finite drawing layers are rejected before periodic geometry requests')
            assert not errors,errors
            browser.close()
    except Exception as error:
        failure=f'{type(error).__name__}: {error}'
    finally:
        if server:server.should_exit=True;worker.join(10)
        report={'checks':checks,'javascript_errors':errors,'failure':failure,
                'scope':'Real local periodic geometry; synthetic FCC fixture tests an unresolved channel, not an experimental material.'}
        (output/'geometry-smoke-report.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf8')
        print(json.dumps(report,indent=2))
    return 1 if failure is not None else 0


if __name__=='__main__':raise SystemExit(main())
