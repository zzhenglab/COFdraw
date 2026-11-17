"""Exercise COF-level inspection, PXRD, catalog queries and report retrieval."""
import argparse,csv,io,json,shutil,socket,sys,threading,time
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--base-url');ap.add_argument('--output-dir',type=Path,required=True);args=ap.parse_args()
 out=args.output_dir.resolve();out.mkdir(parents=True,exist_ok=True);server=worker=None
 if args.base_url:base=args.base_url.rstrip('/')
 else:
  from cof_vision.cofdraw_desktop import create_desktop_app
  import uvicorn
  runtime=out/'runtime-data';shutil.copytree(ROOT/'data',runtime,dirs_exist_ok=True,ignore=shutil.ignore_patterns('exports','*.log'))
  with socket.socket() as sock:sock.bind(('127.0.0.1',0));port=sock.getsockname()[1]
  server=uvicorn.Server(uvicorn.Config(create_desktop_app(runtime),host='127.0.0.1',port=port,log_level='warning'))
  worker=threading.Thread(target=server.run,daemon=True);worker.start()
  for _ in range(200):
   if server.started:break
   time.sleep(.05)
  assert server.started;base=f'http://127.0.0.1:{port}'
 errors=[];checks=[];failure=None
 try:
  from playwright.sync_api import sync_playwright
  with sync_playwright() as p:
   exe=next(x for x in [Path(r'C:\Program Files\Google\Chrome\Application\chrome.exe'),Path(r'C:\Program Files (x86)\Microsoft\Edge\Application\msedge.exe')] if x.is_file())
   browser=p.chromium.launch(executable_path=str(exe),headless=True);page=browser.new_page(viewport={'width':1600,'height':1100},bypass_csp=True);page.set_default_timeout(30000)
   page.on('pageerror',lambda e:errors.append(str(e)))
   page.goto(base);page.wait_for_function('window.COFdrawScience && window.COFdrawStudio.state.online')
   page.locator('#open-science').click();page.locator('#crystal-template').select_option('hcb-imine-aa');page.locator('#build-crystal').click()
   page.wait_for_function("document.getElementById('science-status').textContent.startsWith('Periodic model generated')")
   page.locator('#viewer-axis').select_option('c');page.locator('#viewer-hide-h').check();page.locator('#show-atom-table').click()
   assert page.locator('#atom-table-panel tr').count()==67
   page.locator('#show-atom-table').click();checks.append('Lattice-axis view, hydrogen visibility and fractional atom table')
   page.locator('[data-science-tab="inspect"]').click();page.locator('#screen-pores').check();page.locator('#screen-current').click()
   page.wait_for_function("document.getElementById('science-status').textContent.startsWith('Structure check complete')")
   assert '2D' in page.locator('#inspection-results').inner_text()
   assert page.locator('#inspection-results .metric-card').count()>=9
   page.screenshot(path=str(out/'cofdraw-structure-check.png'));checks.append('Periodic composition,2D network and geometric pore screening displayed')
   with page.expect_download() as download:page.locator('#comparison-results button').filter(has_text='Save CSV').click()
   download.value.save_as(out/'screen.csv');assert 'density_g_cm3' in (out/'screen.csv').read_text('utf-8')
   checks.append('Comparison CSV downloads through actual UI')
   with page.expect_response('**/api/drawing/science/screen') as batch:
    page.locator('#screen-cif-files').set_input_files([str(ROOT/'data/crystal-models/hcb-imine-aa.cif'),str(ROOT/'data/crystal-models/dia-imine.cif')])
   assert batch.value.status==200
   page.wait_for_function("document.getElementById('science-status').textContent.includes('2 of 2')")
   assert 'dia-imine.cif' in page.locator('#comparison-results').inner_text()
   checks.append('Multiple CIF upload compares two actual structures under common settings')
   page.locator('[data-science-tab="pxrd"]').click();page.locator('#calculate-pxrd').click()
   page.wait_for_function("document.getElementById('science-status').textContent.startsWith('Model PXRD calculated')")
   assert page.locator('#pxrd-plot svg path').count()>5
   assert page.locator('#pxrd-results tr').count()>3
   page.screenshot(path=str(out/'cofdraw-pxrd.png'))
   with page.expect_download() as download:page.locator('#pxrd-results button').filter(has_text='Save CSV').click()
   download.value.save_as(out/'pxrd.csv');rows=list(csv.DictReader((out/'pxrd.csv').read_text('utf-8').splitlines()));assert len(rows)>1000
   checks.append('Physical model PXRD preview, peak table and numeric CSV export')
   page.locator('#xrd-fwhm').fill('0.02');page.locator('#calculate-pxrd').click()
   page.wait_for_function("document.getElementById('science-status').textContent.startsWith('Model PXRD calculated')")
   assert page.locator('#pxrd-plot svg').count()==1
   checks.append('Narrow PXRD width automatically uses adequate profile sampling')
   page.locator('[data-science-tab="knowledge"]').click();page.locator('#knowledge-query').fill('Pyurea');page.locator('#knowledge-type').select_option('source');page.locator('#query-knowledge').click()
   page.wait_for_function("document.querySelectorAll('#knowledge-results .knowledge-card').length>0")
   assert 'source:pyurea-source-layout-v1' in page.locator('#knowledge-results').inner_text()
   page.screenshot(path=str(out/'cofdraw-knowledge.png'));checks.append('Source-information query returns exact record, qualification and evidence')
   page.locator('#query-reports').click();page.wait_for_function("document.querySelectorAll('#saved-reports > button').length>=3")
   page.locator('#saved-reports > button').first.click();page.wait_for_function("document.querySelectorAll('#saved-reports pre').length>0")
   report=json.loads(page.locator('#saved-reports pre').first.text_content());assert report['report_id'] and report['input_sha256'] and report['inputs']
   retrieved=page.request.get(base+'/api/drawing/science/reports/'+report['report_id']).json();assert retrieved==report
   caps=page.request.get(base+'/api/drawing/science/capabilities').json();assert 'pxrd' in caps['operations'] and caps['workflows']
   checks.append('Agents discover schemas and retrieve identical saved reports with input hashes')
   page.locator('#knowledge-results button').filter(has_text='Open source drawing').click();page.wait_for_function("window.COFdrawStudio.state.sourceId==='pyurea-source-layout-v1' && !window.COFdrawStudio.state.busy")
   assert not page.locator('#science-dialog').is_visible();checks.append('Information result opens its actual saved source drawing')
   page.set_viewport_size({'width':390,'height':844});page.locator('#open-science').click();page.locator('[data-science-tab="knowledge"]').click()
   assert page.locator('#science-dialog').evaluate('e=>e.scrollWidth<=e.clientWidth+1');checks.append('Materials tools fit390px viewport without horizontal overflow')
   assert not errors,errors;browser.close()
 except Exception as exc:failure=f'{type(exc).__name__}: {exc}';raise
 finally:
  if server:server.should_exit=True;worker.join(timeout=10)
  (out/'materials-smoke-report.json').write_text(json.dumps({'checks':checks,'javascript_errors':errors,'failure':failure},indent=2),encoding='utf-8')
  print(json.dumps({'passed':len(checks),'failure':failure,'javascript_errors':errors}))
if __name__=='__main__':main()
