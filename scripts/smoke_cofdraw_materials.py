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
   def quick_check():
    page.locator('[data-science-tab="model"]').click()
    box=page.locator('#science-inspect')
    if not box.evaluate('el=>el.open'):box.locator('summary').click()
   page.goto(base);page.wait_for_function('window.COFdrawScience && window.COFdrawStudio.state.online')
   page.locator('#mode-manual').click();page.locator('#build').click()
   page.wait_for_function('window.COFdrawStudio.state.result && !window.COFdrawStudio.state.busy')
   page.locator('#open-science').click()
   page.wait_for_function("window.COFdrawScience.origin==='current-periodic'")
   assert page.evaluate('window.COFdrawScience.structure.atoms.length')==57
   assert page.evaluate('window.COFdrawScience.viewer.bonds.length')>0
   assert 'unrelaxed periodic cell' in page.locator('#viewer-provenance').inner_text()
   with page.expect_download() as download:page.locator('#structure-exports button').filter(has_text='CIF').click()
   download.value.save_as(out/'current-drawing.cif')
   assert 'COFdraw unrelaxed AA starting model' in (out/'current-drawing.cif').read_text('utf8')
   checks.append('Opening 3D & properties produces a downloadable matched CIF for the current drawing')
   page.locator('#finite-view-options > summary').click()
   page.locator('#build-patch-3d').click()
   page.wait_for_function('window.COFdrawScience.finiteModel && !document.querySelector(".science-controls").classList.contains("science-busy")')
   quick_check();page.locator('#screen-current').click()
   page.wait_for_function("document.getElementById('science-status').textContent.startsWith('Checked the displayed finite')")
   assert 'Displayed atom inventory' in page.locator('#inspection-results').inner_text()
   assert page.locator('#screen-pores').is_disabled()
   checks.append('Current finite drawing checks composition, bonds and contacts without requiring CIF')
   page.locator('[data-science-tab="porosity"]').click();page.locator('#measure-current-finite-pore').click()
   page.wait_for_function("document.getElementById('science-status').textContent.startsWith('Measured the current')")
   assert page.locator('#porosity-results .metric-card').count()==2
   assert 'lower bound for the projected opening' in page.locator('#porosity-results').inner_text()
   assert page.locator('#calculate-pores').is_disabled()
   checks.append('Finite drawing opening can be measured directly from pore tools with explicit periodic-property limits')
   page.locator('[data-science-tab="pxrd"]').click();page.locator('#calculate-pxrd').click()
   page.wait_for_function("document.getElementById('science-status').textContent.startsWith('Finite-stack scattering calculated')")
   assert 'Debye scattering' in page.locator('#pxrd-results').inner_text()
   assert page.locator('#pxrd-results table').count()==0 and page.locator('#xrd-fwhm').is_disabled()
   with page.expect_download() as download:page.locator('#pxrd-results button').filter(has_text='Save CSV').click()
   download.value.save_as(out/'finite-powder.csv')
   assert 'intensity_electron_squared' in (out/'finite-powder.csv').read_text()
   checks.append('Displayed finite atoms yield actual Debye scattering and CSV without fictitious hkl')
   page.locator('#compare-pxrd-stacks').click()
   page.wait_for_function("document.getElementById('science-status').textContent.startsWith('5 stacking patterns compared')",timeout=120000)
   assert page.locator('#pxrd-plot [data-series-label]').count()==5
   slider=page.locator('#pattern-separation')
   slider.focus();page.keyboard.press('End')
   assert slider.input_value()=='1' and page.locator('#pxrd-plot [data-trace]').count()==5
   page.screenshot(path=str(out/'cofdraw-finite-stacking-patterns.png'))
   checks.append('Five finite stacking patterns compare the same drawing under chosen transforms')
   quick_check();page.locator('#inspect-model-state .build-cif-from-drawing').click()
   page.wait_for_function("document.getElementById('science-status').textContent.startsWith('Periodic model and CIF ready')")
   assert page.evaluate('window.COFdrawScience.structure.atoms.length')==57
   assert page.evaluate('window.COFdrawScience.viewer.bonds.length')>0
   assert 'unrelaxed periodic cell' in page.locator('#viewer-provenance').inner_text()
   quick_check()
   page.locator('#screen-current').click()
   page.wait_for_function("document.getElementById('science-status').textContent.startsWith('Structure check complete')")
   assert '2D' in page.locator('#inspection-results').inner_text()
   checks.append('Exact current TFB/TAPB design opens its matched periodic cell, bonds and checks without importing CIF')
   page.locator('[data-science-tab="model"]').click();page.locator('#structure-source').select_option('periodic')
   assert page.locator('#pxrd-plot svg').count()==0
   page.locator('#crystal-template').select_option('hcb-imine-aa');page.locator('#build-crystal').click()
   page.wait_for_function("COFdrawScience.structure?.periodic && !document.querySelector('.science-busy')")
   page.locator('#viewer-axis').select_option('c');page.locator('#viewer-hide-h').check();page.locator('#show-atom-table').click()
   assert page.locator('#atom-table-panel tr').count()==67
   page.locator('#show-atom-table').click();checks.append('Lattice-axis view, hydrogen visibility and fractional atom table')
   canvas=page.locator('#structure-canvas')
   canvas.scroll_into_view_if_needed()
   point=page.evaluate("""() => {
     const v=window.COFdrawScience.viewer, r=document.getElementById('structure-canvas').getBoundingClientRect();
     const i=v.atoms.findIndex((a,i)=>a.element!=='H' && v.projected[i].x>20 && v.projected[i].x<r.width-20
       && v.projected[i].y>20 && v.projected[i].y<r.height-20);
     return {x:v.projected[i].x,y:v.projected[i].y};
   }""")
   page.locator('#viewer-measure').check()
   canvas.click(position=point);canvas.click(position=point)
   assert 'different atom' in page.locator('#viewer-distance').inner_text()
   page.locator('#viewer-measure').uncheck()
   assert page.locator('#viewer-distance').inner_text()==''
   canvas.scroll_into_view_if_needed();box=canvas.bounding_box()
   page.mouse.move(box['x']+point['x'],box['y']+point['y']);page.mouse.down()
   page.mouse.move(box['x']+point['x']+15,box['y']+point['y']+12,steps=3);page.mouse.up()
   assert page.locator('#viewer-axis').input_value()=='free'
   checks.append('Dragging clears the lattice axis and distance measurement requires two distinct atoms')
   quick_check();page.locator('#screen-pores').check();page.locator('#screen-current').click()
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
   angle=page.locator('#pxrd-results tr').nth(1).locator('td').first.inner_text()
   page.locator('#peak-check summary').click()
   page.locator('#observed-peaks').fill(angle+', 0.9')
   with page.expect_response('**/api/drawing/science/pxrd-peaks') as compared:page.locator('#check-peaks').click()
   assert compared.value.ok and any(r['candidate_indices'] for r in compared.value.json()['observed'])
   page.get_by_text('outside range',exact=True).wait_for()
   assert page.locator('#pxrd-plot svg line').count()>=1
   page.screenshot(path=str(out/'observed-peaks.png'))
   assert page.locator('#pattern-view').is_visible() and page.locator('#model-view').is_hidden()
   with page.expect_download() as exported:page.locator('#save-pattern-svg').click()
   svg_path=out/'observed-peaks.svg';exported.value.save_as(svg_path)
   import xml.etree.ElementTree as ET
   metadata=json.loads(ET.parse(svg_path).getroot().find('{http://www.w3.org/2000/svg}metadata').text)
   assert metadata['report_id']==compared.value.json()['report_id'] and len(metadata['observed'])==2
   page.locator('#return-to-model').click();assert page.locator('#structure-canvas').is_visible()
   page.locator('#return-to-pattern').click();assert page.locator('#pattern-view').is_visible()
   checks.append('Measured peak ticks, candidates and out-of-range observations are visible for the active periodic model')
   page.locator('#xrd-fwhm').fill('0.02');page.locator('#calculate-pxrd').click()
   page.wait_for_function("document.getElementById('science-status').textContent.startsWith('Model PXRD calculated')")
   assert page.locator('#pxrd-plot svg').count()==1
   checks.append('Narrow PXRD width automatically uses adequate profile sampling')
   page.locator('#xrd-fwhm').fill('0.15');page.locator('#compare-pxrd-stacks').click()
   page.wait_for_function("document.getElementById('science-status').textContent.startsWith('4 stacking patterns compared')",timeout=120000)
   assert page.locator('#pxrd-plot [data-series-label]').count()==4
   assert 'explicit periodic AA' in page.locator('#pxrd-results').inner_text()
   assert slider.input_value()=='0'
   paths=page.locator('#pxrd-plot [data-trace]')
   bounds="paths => paths.map(p=>{const r=p.getBoundingClientRect();return [r.top,r.bottom]})"
   overlay=paths.evaluate_all(bounds)
   assert all(abs(bottom-overlay[0][1])<1 for _,bottom in overlay)
   originals=paths.evaluate_all("paths=>paths.map(p=>p.getAttribute('d'))")
   slider.fill('0.5')
   middle=paths.evaluate_all(bounds)
   assert middle[0][0]<middle[1][0]<middle[0][1]
   slider.focus();page.keyboard.press('End')
   split=paths.evaluate_all(bounds)
   assert all(right<split[i+1][0] for i,(_,right) in enumerate(split[:-1]))
   assert originals==paths.evaluate_all("paths=>paths.map(p=>p.getAttribute('d'))")
   with page.expect_download() as separated:page.locator('#save-pattern-svg').click()
   split_svg=out/'stacking-vertical.svg';separated.value.save_as(split_svg)
   saved=json.loads(ET.parse(split_svg).getroot().find('{http://www.w3.org/2000/svg}metadata').text)
   assert saved['vertical_separation']==1 and len(saved['series'])==4
   page.screenshot(path=str(out/'cofdraw-periodic-stacking-pxrd.png'))
   slider.focus();page.keyboard.press('Home');assert slider.input_value()=='0'
   assert paths.evaluate_all(bounds)==overlay
   page.locator('#pxrd-results button').filter(has_text='View this arrangement').nth(1).click()
   assert slider.is_hidden() and page.locator('#viewer-title').inner_text().startswith('AB')
   page.locator('#compare-pxrd-stacks').click()
   page.wait_for_function("document.getElementById('science-status').textContent.startsWith('4 stacking patterns compared')",timeout=120000)
   checks.append('Vertical separation preserves angles and exports the stacked traces; comparing again from AB uses the parent layer')
   checks.append('Four actual periodic stacking cells produce overlaid Bragg PXRD and individual report downloads')
   assert page.locator('[data-science-tab="knowledge"]').count()==0
   queried=page.request.post(base+'/api/drawing/science/query',data={'query':'Pyurea','record_types':['source'],'limit':20}).json()
   assert any(r['id']=='source:pyurea-source-layout-v1' for r in queried['records'])
   checks.append('Information panel removed; catalog query API still returns qualified source records')
   reports=page.request.get(base+'/api/drawing/science/reports').json()['reports'];assert len(reports)>=3
   report=page.request.get(base+'/api/drawing/science/reports/'+reports[0]['report_id']).json();assert report['report_id'] and report['input_sha256'] and report['inputs']
   caps=page.request.get(base+'/api/drawing/science/capabilities').json();assert 'pxrd' in caps['operations'] and caps['workflows']
   checks.append('Agents discover schemas and retrieve saved reports with input hashes through the API')
   page.locator('[data-science-tab="model"]').click()
   page.locator('#cif-file').set_input_files({'name':'oversize.cif','mimeType':'text/plain','buffer':b' ' * 2000001})
   page.wait_for_function("document.getElementById('science-status').textContent.includes('2 MB')")
   page.locator('#network-hypothesis > summary').click()
   page.locator('#build-networks').click()
   page.wait_for_function("document.getElementById('science-status').textContent.startsWith('Generate a periodic model')")
   assert page.evaluate('window.COFdrawScience.structure===null')
   checks.append('A failed replacement CIF import clears the previous frame and network-copy parent')
   page.set_viewport_size({'width':390,'height':844});quick_check()
   assert page.locator('#science-dialog').evaluate('e=>e.scrollWidth<=e.clientWidth+1');checks.append('Materials tools fit390px viewport without horizontal overflow')
   assert not errors,errors;browser.close()
 except Exception as exc:failure=f'{type(exc).__name__}: {exc}';raise
 finally:
  if server:server.should_exit=True;worker.join(timeout=10)
  (out/'materials-smoke-report.json').write_text(json.dumps({'checks':checks,'javascript_errors':errors,'failure':failure},indent=2),encoding='utf-8')
  print(json.dumps({'passed':len(checks),'failure':failure,'javascript_errors':errors}))
if __name__=='__main__':main()
