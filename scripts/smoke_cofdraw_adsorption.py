"""Real CPU RASPA browser workflow plus explicitly separate UI race fixtures."""
import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import shutil
import socket
import sys
import threading
import time
import zipfile

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--output-dir',type=Path,required=True);parser.add_argument('--base-url')
    args=parser.parse_args();output=args.output_dir.resolve();output.mkdir(parents=True,exist_ok=True)
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
        if not server.started:raise RuntimeError('Local service did not start')
        base=f'http://127.0.0.1:{port}'
    checks=[];errors=[];failure=None;real_jobs=[];fixtures=[]
    try:
        with sync_playwright() as p:
            executable=next((path for path in [Path(r'C:\Program Files\Google\Chrome\Application\chrome.exe'),
                Path(r'C:\Program Files (x86)\Microsoft\Edge\Application\msedge.exe')] if path.is_file()),None)
            browser=p.chromium.launch(executable_path=str(executable) if executable else None,headless=True)
            page=browser.new_page(viewport={'width':1440,'height':1100},bypass_csp=True);page.set_default_timeout(45000)
            page.on('pageerror',lambda error:errors.append(str(error)))
            started=[]
            def observe(response):
                if response.url.endswith('/science/adsorption/jobs') and response.request.method=='POST' and response.status==202:
                    started.append(response.json())
            page.on('response',observe)
            page.goto(base);page.wait_for_function('window.COFdrawStudio?.state.online && !document.getElementById("build").disabled')
            engine=page.request.get(base+'/api/drawing/science/adsorption/engine').json()
            assert engine.get('available'),engine
            page.locator('#mode-manual').click();page.locator('#build').click()
            page.wait_for_function('COFdrawStudio.state.result && !COFdrawStudio.state.busy')
            page.locator('#open-science').click();page.locator('#structure-source').select_option('periodic')
            page.locator('#crystal-template').select_option('hcb-imine-aa');page.locator('#build-crystal').click()
            page.wait_for_function('COFdrawScience.structure?.periodic===true')
            page.locator('[data-science-tab="adsorption"]').click()
            page.locator('#ads-gas').select_option('Ar');page.locator('#ads-charges').select_option('zero')
            page.locator('#ads-heat').uncheck();page.locator('#ads-pressures').fill('0.1, 1')
            page.locator('#ads-temperature').fill('298.15')
            page.locator('#science-adsorption details').first.evaluate('element=>element.open=true')
            for field in ('ads-cycles','ads-initial','ads-widom'):page.locator('#'+field).fill('100')
            page.locator('#ads-seed').fill('42137')

            def run(label):
                before=len(started);page.locator('#ads-start').click()
                deadline=time.monotonic()+240;last=None
                while time.monotonic()<deadline:
                    page.wait_for_timeout(300)
                    if len(started)<=before:continue
                    identifier=started[-1]['job_id']
                    status=page.request.get(base+'/api/drawing/science/adsorption/jobs/'+identifier).json()
                    stage=(status['status'],status['completed_stages'],status['stage'])
                    if stage!=last:print(json.dumps({'run':label,'status':stage}),flush=True);last=stage
                    if status['status'] in ('complete','failed','cancelled','interrupted'):
                        assert status['status']=='complete',status.get('message')
                        page.wait_for_function('document.getElementById("ads-progress").textContent.startsWith("complete")')
                        real_jobs.append({'label':label,'job_id':identifier,'input_sha256':status['input_sha256'],
                                          'completed_stages':status['completed_stages'],'engine':engine})
                        (output/(label+'-job.json')).write_text(json.dumps(status,indent=2)+'\n',encoding='utf8')
                        return status
                raise TimeoutError(f'{label} did not complete within 240 seconds: '+page.locator('#ads-progress').inner_text())

            def assert_visible_conditions(composition, reservoir):
                # Check the visible paragraphs preceding the plot, not the raw
                # diagnostics disclosure where all inputs are also available.
                notes=page.locator('#ads-results').evaluate('element => {const result=[];for(const child of element.children){if(child.tagName.toLowerCase()==="svg")break;if(child.tagName==="P")result.push(child.textContent);}return result.join("\\n");}')
                assert '298.15 K' in notes and composition in notes,notes
                assert 'Uncharged framework approximation' in notes,notes
                assert 'Host\u2013guest electrostatics from framework partial charges are omitted' in notes,notes
                assert engine['profile_name'] in notes and reservoir in notes,notes
                assert 'nuclear quantum effects' not in notes,notes

            pure=run('argon')
            result=pure['result'];assert len(result['isotherm'])==2 and pure['completed_stages']==3
            assert result['inputs']['gases']==['Ar'] and result['inputs']['charge_mode']=='zero'
            assert result['zero_loading_heat'] is None
            assert result['input_sha256']==hashlib.sha256(json.dumps(result['inputs'],sort_keys=True).encode()).hexdigest()
            for stage in [result['widom'],*result['isotherm']]:
                assert stage['manifest']['engine']=='RASPA3' and stage['manifest']['engine_version']=='3.1.0'
                assert stage['manifest']['cycles']==100 and stage['manifest']['charge_mode']=='zero'
                assert stage['manifest']['structure_sha256']==result['structure_sha256']
            for row in result['isotherm']:
                component=row['components'][0]
                assert component['gas']=='Ar' and math.isfinite(component['loading_mol_kg']) and component['loading_mol_kg']>=0
            assert page.locator('#ads-results svg[aria-label="Simulated absolute adsorption isotherm"]').count()==1
            assert page.locator('#ads-results table').count()>=2
            assert_visible_conditions('Ar 100%', 'Ideal gas reservoir')
            checks.append('Real RASPA3 Ar Widom plus two 100-cycle GCMC pressure points show the plot, tables, temperature, composition, force field and uncharged/ideal-reservoir assumptions')
            page.locator('#science-adsorption details').first.evaluate('element=>element.open=false')
            page.locator('#ads-results').scroll_into_view_if_needed()
            page.screenshot(path=str(output/'adsorption-argon-results.png'))
            downloads={}
            for ext,label in [('json','Save JSON'),('csv','Save CSV'),('zip','Save complete run')]:
                with page.expect_download() as download:page.locator('#ads-results button').filter(has_text=label).click()
                destination=output/f'argon-run.{ext}';download.value.save_as(destination);downloads[ext]=destination
            assert json.loads(downloads['json'].read_text('utf8'))==result
            rows=list(csv.DictReader(downloads['csv'].open(encoding='utf8',newline='')))
            assert len(rows)==2 and {float(row['pressure_Pa']) for row in rows}=={10000.,100000.}
            with zipfile.ZipFile(downloads['zip']) as archive:
                assert archive.testzip() is None
                hashes=json.loads(archive.read('sha256.json'))
                assert all(hashlib.sha256(archive.read(name)).hexdigest()==digest for name,digest in hashes.items())
                assert all(f'{stage}/engine-input.json' in archive.namelist() for stage in ['widom','pressure-01','pressure-02'])
                request=json.loads(archive.read('request.json'))
                assert request==result['inputs']
            checks.append('Actual JSON/CSV/ZIP downloads retain settings, engine manifests, full stage inputs and verified archive checksums')

            page.locator('#ads-gas').select_option('CO2,N2');page.locator('#ads-fraction').fill('0.15')
            page.locator('#ads-pressures').fill('1');page.locator('#science-adsorption details').first.evaluate('element=>element.open=true')
            page.locator('#ads-fugacity').fill('1, 1')
            mixture=run('co2-n2')
            mixed=mixture['result'];assert mixed['inputs']['mole_fractions']==[.15,.85]
            assert mixed['inputs']['fugacity_coefficients']==[1.,1.]
            assert len(mixed['isotherm'])==1 and len(mixed['isotherm'][0]['components'])==2
            for component in mixed['isotherm'][0]['components']:
                assert math.isfinite(component['loading_mol_kg']) and component['loading_mol_kg']>=0
            assert 'Henry selectivity' in page.locator('#ads-results').inner_text()
            assert 'Selectivity' in page.locator('#ads-results').inner_text()
            manifest=mixed['isotherm'][0]['manifest']
            assert manifest['partial_fugacities_Pa']==[15000.,85000.]
            assert_visible_conditions('CO2 15% / N2 85%', 'Supplied fugacity coefficients')
            checks.append('A real CO2/N2 mixture retains physical composition/fugacity and visibly shows the temperature, 15/85 composition, supplied coefficients, uncharged approximation and qualified selectivity')
            page.locator('#science-adsorption details').first.evaluate('element=>element.open=false')
            page.locator('#ads-results').scroll_into_view_if_needed();page.screenshot(path=str(output/'adsorption-mixture-results.png'))
            page.set_viewport_size({'width':390,'height':844})
            assert page.evaluate('document.documentElement.scrollWidth<=innerWidth+1')
            assert page.evaluate('document.getElementById("science-dialog").scrollWidth<=document.getElementById("science-dialog").clientWidth+1')
            page.locator('#ads-results svg').evaluate('element=>element.scrollIntoView({block:"center"})')
            page.wait_for_timeout(150)
            plot_box=page.locator('#ads-results svg').bounding_box()
            if not (plot_box and plot_box['y']>=0 and plot_box['y']+plot_box['height']<=844):
                print(json.dumps({'mobile_plot':plot_box,'scroll_containers':page.evaluate('Array.from(document.querySelectorAll(".science-workspace,.science-controls,#ads-results")).map(el=>({class:el.className,y:el.getBoundingClientRect().y,height:el.clientHeight,scrollHeight:el.scrollHeight,scrollTop:el.scrollTop}))')}),flush=True)
            assert plot_box and plot_box['y']>=0 and plot_box['y']+plot_box['height']<=844
            page.screenshot(path=str(output/'adsorption-phone.png'))
            checks.append('The actual adsorption result fits a 390px viewport without page or dialog horizontal overflow')
            page.set_viewport_size({'width':1440,'height':1100})

            # These two transport fixtures do not stand in for scientific runs.
            # They test cancellation/epoch handling only, after real runs above.
            fake_id='f'*32;state={'cancelled':False,'stale':False,'cancel_before_reply':True,'calls':0}
            def transport(route):
                path=route.request.url.split('/science/adsorption/',1)[1]
                if path=='jobs' and route.request.method=='POST':
                    if state['cancel_before_reply']:page.locator('#ads-cancel').click()
                    route.fulfill(status=202,json={'job_id':fake_id,'status':'queued'})
                elif path==f'jobs/{fake_id}/cancel':
                    state['cancelled']=True;state['calls']+=1;route.fulfill(json={'status':'cancelled'})
                elif path==f'jobs/{fake_id}':
                    terminal='complete' if state['stale'] else 'cancelled' if state['cancelled'] else 'running'
                    route.fulfill(json={'status':terminal,'completed_stages':0,'total_stages':1,'stage':'UI-only delayed transport fixture',
                        'message':'Not a scientific simulation','result':mixed if state['stale'] else None})
                else:route.continue_()
            page.route('**/api/drawing/science/adsorption/**',transport)
            page.locator('#ads-start').click()
            page.wait_for_function('document.getElementById("ads-progress").textContent.startsWith("cancelled")')
            assert state['cancelled'] and state['calls']==1
            assert page.locator('#ads-results table').count()==0
            fixtures.append('Pending-POST cancel transport fixture; no engine job launched')
            checks.append('Cancel during the pending start request is retained and sent once the job identifier arrives')
            state.update(cancelled=False,stale=False,cancel_before_reply=False,calls=0)
            page.locator('#ads-start').click()
            page.wait_for_function('document.getElementById("ads-progress").textContent.startsWith("running")')
            page.locator('[data-science-tab="model"]').click()
            page.locator('#crystal-template').select_option('hcb-imine-ab');page.locator('#build-crystal').click()
            page.wait_for_function('COFdrawScience.structure?.atoms.length===132')
            assert state['cancelled']
            state['stale']=True
            page.locator('[data-science-tab="adsorption"]').click()
            page.wait_for_function('document.getElementById("ads-progress").textContent.includes("Previous simulation stopped")')
            assert page.locator('#ads-results table').count()==0
            fixtures.append('Late completion fixture reuses already verified data only to check stale-result suppression')
            checks.append('Changing the periodic structure cancels pending work and prevents a late result from attaching to the new model')
            page.unroute('**/api/drawing/science/adsorption/**',transport)
            assert not errors,errors
            browser.close()
    except Exception as error:failure=f'{type(error).__name__}: {error}'
    finally:
        if server:server.should_exit=True;worker.join(10)
        report={'checks':checks,'javascript_errors':errors,'failure':failure,'real_jobs':real_jobs,'ui_only_fixtures':fixtures,
                'scope':'Real 100-cycle RASPA CPU workflows validate execution and reporting, not convergence or experimental accuracy. UI race fixtures are separately identified.'}
        (output/'adsorption-smoke-report.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf8')
        print(json.dumps(report,indent=2))
    return 1 if failure is not None else 0


if __name__=='__main__':raise SystemExit(main())
