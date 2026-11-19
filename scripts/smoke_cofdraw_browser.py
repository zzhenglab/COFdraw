"""Actual local COFdraw browser flow using only bundled resources."""
from __future__ import annotations
import argparse, json, socket, sys, threading, time, traceback
from pathlib import Path
from urllib.parse import urlsplit
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--output-dir",type=Path,required=True)
    parser.add_argument("--data-dir",type=Path,help="Data directory for an owned Python server; required unless --base-url is given")
    parser.add_argument("--base-url",help="Existing local server, e.g. http://127.0.0.1:8765; it is never stopped by this script")
    parser.add_argument("--browser",type=Path)
    args=parser.parse_args()
    if not args.base_url and not args.data_dir:parser.error("--data-dir is required without --base-url")
    if args.base_url:
      parsed=urlsplit(args.base_url)
      if parsed.scheme!="http" or parsed.hostname not in ("127.0.0.1","localhost","::1") or parsed.username or parsed.password or parsed.query or parsed.fragment or parsed.path not in ("","/"):
        parser.error("--base-url must be an HTTP loopback server origin without credentials, path, query, or fragment")
    out=args.output_dir.resolve();out.mkdir(parents=True,exist_ok=True)
    candidates=[args.browser] if args.browser else [Path(r"C:\Program Files\Google\Chrome\Application\chrome.exe"),Path(r"C:\Program Files (x86)\Microsoft\Edge\Application\msedge.exe")]
    executable=next((p for p in candidates if p and p.is_file()),None)
    if executable is None: raise RuntimeError("An installed Chrome or Edge browser is required.")
    from playwright.sync_api import sync_playwright
    server=worker=None
    if args.base_url:
      base_url=args.base_url.rstrip("/")
    else:
      import uvicorn
      from cof_vision.cofdraw_desktop import create_desktop_app as create_app
      with socket.socket() as sock:sock.bind(("127.0.0.1",0));port=sock.getsockname()[1]
      server=uvicorn.Server(uvicorn.Config(create_app(args.data_dir),host="127.0.0.1",port=port,log_level="warning"))
      worker=threading.Thread(target=server.run,daemon=True);worker.start()
      deadline=time.monotonic()+15
      while not server.started and time.monotonic()<deadline:time.sleep(.05)
      if not server.started:
        server.should_exit=True;worker.join(timeout=10)
        raise RuntimeError("Local test server did not start.")
      base_url=f"http://127.0.0.1:{port}"
    checks=[];errors=[];api_errors=[];failure=None;page=None;monomer_requests=[]
    try:
      with sync_playwright() as p:
        browser=p.chromium.launch(executable_path=str(executable),headless=True)
        page=browser.new_page(viewport={"width":1512,"height":1060},device_scale_factor=1)
        page.set_default_timeout(120000)
        page.on("pageerror",lambda error:errors.append(str(error)))
        def record_response(response):
          if "/api/drawing/" in response.url and response.status>=400:
            try:api_errors.append({"url":response.url,"status":response.status,"body":response.text()[:3000]})
            except Exception:api_errors.append({"url":response.url,"status":response.status})
        page.on("response",record_response)
        def record_request(request):
          if request.url.endswith("/api/drawing/monomer") and request.method=="POST":
            try:monomer_requests.append(request.post_data_json.get("smiles"))
            except Exception:pass
        page.on("request",record_request)
        def build():
          with page.expect_response("**/api/drawing/render") as rendered:
            with page.expect_response("**/api/drawing/prepare") as prepared:page.locator("#build").click()
            assert prepared.value.status==200,prepared.value.text()
          return rendered
        page.goto(f"{base_url}/cofdraw.html")
        page.get_by_text("Local chemistry engine · ready",exact=True).wait_for()
        page.locator('#mode-manual').click()
        page.locator('#appearance-options').evaluate('e=>e.open=true')
        assert page.locator('.identity-summary').is_visible()
        library=page.request.get(f"{base_url}/api/drawing/library").json()
        page.wait_for_function("!document.getElementById('build').disabled")
        page.locator('#tool-jump').select_option('balance')
        page.locator('#calculate-balance').click()
        page.get_by_text('Weigh block B',exact=True).wait_for()
        assert page.locator('#balance-results tbody tr,#balance-results tr').count()>=3
        page.screenshot(path=str(out/'reagent-balance.png'))
        page.locator('#balance-mass').fill('20')
        assert not page.locator('#balance-results').inner_text()
        page.locator('#balance-dialog button[aria-label]').click()
        checks.append('Tools opens the current reagent pair before rendering; editing the mass clears the previous weigh-out')
        page.wait_for_function("() => document.querySelector('#preview-a').complete && !document.querySelector('#preview-a').hidden")
        page.wait_for_function("() => document.querySelector('.brand-logo').complete && document.querySelector('.brand-logo').naturalWidth>0")
        if library.get("topology_catalog"):
          original_a=page.locator("#smiles-a").input_value()
          page.locator("#open-monomer-library").click()
          page.wait_for_function("() => document.querySelectorAll('.library-row[data-preview-state=ready]').length>=5")
          assert page.locator(".library-thumbnail img:not([hidden])").count()>=5
          assert page.locator(".library-row:not([data-preview-state])").count()>0,"Off-screen cards should remain lazy"
          assert not page.locator("#library-detail-panel").is_visible(),"The structure grid is primary; details are optional"
          grid_columns=lambda:page.locator("#library-list").evaluate("e=>getComputedStyle(e).gridTemplateColumns.split(' ').length")
          assert grid_columns()==5
          page.screenshot(path=str(out/"cofdraw-library-grid.png"),full_page=True)
          for width,columns in [(1920,6),(1100,3),(460,2),(390,1)]:
            page.set_viewport_size({"width":width,"height":1060})
            assert grid_columns()==columns
            assert page.locator("#monomer-library-dialog").evaluate("e=>e.scrollWidth<=e.clientWidth+1")
            if width in (1920,390):page.screenshot(path=str(out/f"cofdraw-library-grid-{width}.png"),full_page=True)
          page.set_viewport_size({"width":1512,"height":1060})
          last=page.locator(".library-row").last
          last.scroll_into_view_if_needed()
          page.wait_for_function("() => document.querySelector('.library-row:last-child').dataset.previewState==='ready'")
          assert last.locator(".library-thumbnail img").evaluate("e=>e.complete&&e.naturalWidth>0")
          page.locator("#library-list").evaluate("e=>e.scrollTop=0")
          first=page.locator(".library-row[data-preview-state=ready]").first
          cached_id=first.get_attribute("data-id")
          cached_item=next(item for item in library["monomer_presets"] if item["id"]==cached_id)
          cached_url=first.locator("img").get_attribute("src")
          prior_calls=monomer_requests.count(cached_item["smiles"])
          page.locator("#library-search").fill(cached_item.get("abbreviation") or cached_id)
          page.wait_for_function("() => document.querySelector('.library-row').dataset.previewState==='ready'")
          assert page.locator(f'[data-id="{cached_id}"] img').get_attribute("src")==cached_url
          assert monomer_requests.count(cached_item["smiles"])==prior_calls
          page.locator("#library-search").fill("")
          checks.append("Real SVG card grid: lazy off-screen loading, cached filter reuse, scrolling, and 1/2/3/5/6 responsive columns")
          # This legacy case checks fixed-linkage compatibility; automatic
          # switching is covered by smoke_cofdraw_studio.py.
          page.keyboard.press("Escape")
          page.locator("#auto-chemistry").uncheck()
          page.locator("#open-monomer-library").click()
          page.locator("#library-sites").select_option("4")
          assert page.locator(".library-row").count()>0
          assert all(text=="4c" for text in page.locator(".library-site-badge").all_text_contents())
          page.locator("#library-step-b").click()
          page.locator("#library-show-all").check()
          assert page.locator(".library-row.unavailable .library-use:disabled").count()>0
          page.locator(".library-row.unavailable .inspect-monomer").first.click()
          assert page.locator("#library-apply").is_disabled()
          assert "Requires" in page.locator("#library-compatibility").inner_text()
          page.wait_for_function("() => !document.querySelector('#library-preview').hidden && document.querySelector('#library-preview').complete")
          page.screenshot(path=str(out/"cofdraw-library.png"),full_page=True)
          page.locator("#library-editor-mode").click()
          page.locator("#library-edit").click()
          page.wait_for_function("() => document.querySelector('#sketch-formula').textContent.length>0")
          with page.expect_download() as single_mol:page.locator("#sketch-export").click()
          single_mol.value.save_as(out/"library-single-monomer.mol")
          page.locator("#sketch-cancel").click()
          assert page.locator("#smiles-a").input_value()==original_a
          checks.append("Broad library4c filter, visible disabled partner reasons, and single-monomer MOL export without changing A")
          page.locator("#open-monomer-library").click()
          page.locator("#library-search").fill("Trip-CHO")
          assert page.locator('[data-id="hfptp_h"]').count()==1
          page.locator("#library-search").fill("TFB")
          page.locator('[data-id="tfb"] .inspect-monomer').click()
          with page.expect_response("**/api/drawing/science/solubility") as first_screen:
            page.locator("#library-solubility").click()
          first = first_screen.value.json()
          measured = next(row for row in first["solvents"] if row["solvent_id"]=="dioxane")
          assert measured["status"]=="reported" and measured["source_doi"]
          transfer = next(row for row in first["solvents"] if row["solvent_id"]=="dmso")
          assert transfer["estimate_method"]=="shared_descriptor_regression" and len(first["solvents"])==46
          page.get_by_text("Estimated · descriptors", exact=True).first.wait_for()
          page.locator("#solubility-output .solubility-reported").first.wait_for()
          page.locator("#solubility-volume").fill("2")
          with page.expect_response("**/api/drawing/science/solubility") as next_screen:
            page.locator("#solubility-run").click()
          scaled = next(row for row in next_screen.value.json()["solvents"] if row["solvent_id"]=="dioxane")
          assert abs(scaled["mass_mg"]-2*measured["mass_mg"])<1e-8
          with page.expect_download() as solvent_csv:
            page.locator("#solubility-exports button",has_text="Save CSV").click()
          solvent_csv.value.save_as(out/"tfb-solubility.csv")
          page.screenshot(path=str(out/"cofdraw-solubility.png"),full_page=True)
          page.get_by_role("button",name="Close solubility screening").click()
          page.locator("#library-search").fill("")
          checks.append("Library aliases and reported solubility, solvent-volume scaling, and CSV export")
          page.locator('[data-id="tfb"] .library-use').click()
          page.wait_for_function("() => document.querySelector('#library-target').value==='b'")
          assert page.locator('[data-id="tfb"]').count()==0
          page.locator("#library-show-all").check()
          assert page.locator('[data-id="tfb"] .library-use').is_disabled()
          page.locator('[data-id="tapb"] .library-use').click()
          page.locator("#monomer-library-dialog").wait_for(state="hidden")
          checks.append("Choose A then B: incompatible partner actions remain visible and disabled")
          page.locator("#open-topology-library").click()
          assert page.locator(".topology-card").count()>0
          page.locator("#topology-search").fill("3D")
          assert page.locator(".topology-card").count()>0
          assert all("Reference only" in text for text in page.locator(".reference-badge").all_text_contents())
          page.screenshot(path=str(out/"cofdraw-topology-reference.png"),full_page=True)
          page.locator('[aria-label="Close topology references"]').click()
          checks.append("3D topology references remain browsable and explicitly nonconstructing")
        page.screenshot(path=str(out/"cofdraw-start.png"),full_page=True)
        checks.append("Local catalog loaded and real RDKit monomer SVG displayed")
        response=build()
        result=response.value.json();assert response.value.status==200,result
        assert result["audit"]["building_block_instances"]==6 and result["status"]=="designed_variant"
        page.wait_for_function("() => !document.querySelector('#render-busy').classList.contains('hidden') === false")
        page.locator("#framework-image").evaluate("i=>i.decode()")
        page.screenshot(path=str(out/"cofdraw-imine.png"),full_page=True)
        checks.append("Monomers prepared and six-block imine patch rendered")
        formats=["png","svg","mol","sdf","json"]
        if result["downloads"].get("cdxml"):formats.insert(4,"cdxml")
        else:
          assert result["audit"]["chemical_exports"]["cdxml"]
          assert page.locator("#export-note").is_visible()
        for ext in formats:
          with page.expect_download() as dl:page.locator("#exports button",has_text=ext.upper()).filter(has_text=f"{ext.upper()}").get_by_text(ext.upper(),exact=True).click()
          dl.value.save_as(out/f"export-imine.{ext}")
        checks.append("/".join(e.upper() for e in formats)+" downloads saved")
        with page.expect_response("**/api/drawing/render") as styled:page.locator('[data-style="red_blue"]').click()
        assert styled.value.status==200,styled.value.text()
        assert styled.value.json()["audit"]["drawing_style"]["id"]=="red_blue"
        page.locator("#framework-image").evaluate("i=>i.decode()")
        page.screenshot(path=str(out/"cofdraw-colors.png"),full_page=True)
        checks.append("Actual colored drawing style rendered")
        if library.get("style_controls"):
          page.locator("#style-controls").evaluate("e=>e.open=true")
          page.locator("#style-palette").select_option("teal_orange")
          page.locator("#style-font_family").select_option("Cambria")
          page.locator("#style-edge_style").select_option("dashed")
          page.locator("#style-transparent_background").check()
          with page.expect_response("**/api/drawing/render") as tuned:page.locator("#apply-style-options").click()
          assert tuned.value.status==200,tuned.value.text()
          opts=tuned.value.json()["audit"]["drawing_style"]["options"]
          assert opts["palette"]=="teal_orange" and opts["font_family"]=="Cambria" and opts["edge_style"]=="dashed" and opts["transparent_background"] is True
          page.locator("#framework-image").evaluate("i=>i.decode()")
          page.screenshot(path=str(out/"cofdraw-style-controls.png"),full_page=True)
          page.locator("#style-controls").evaluate("e=>e.open=false")
          checks.append("Actual palette/font/continuation/transparency overrides applied")
        with page.expect_response("**/api/drawing/render") as random:page.locator("#random-style").click()
        assert random.value.status==200,random.value.text()
        assert random.value.json()["audit"]["drawing_style"]["seed"] is not None
        random_options=random.value.json()["audit"]["drawing_style"].get("options")
        checks.append("Random style returns reproducible seed")
        with page.expect_download() as saved:page.locator("#save-project").click()
        project=out/"COFdraw-project.json";saved.value.save_as(project)
        data=json.loads(project.read_text(encoding="utf-8"));assert data["format"]=="cofdraw-project" and len(data["monomers"])==2
        if library.get("style_controls"):assert data["style_options"]==random_options
        page.locator("#project-file").set_input_files(project)
        page.get_by_text("Project opened as an editable design.",exact=False).wait_for()
        assert page.locator("#scope-badge").inner_text()=="DESIGNED PATCH"
        assert "Tp plus" not in page.locator("#chemistry-note").inner_text()
        checks.append("Project saves and reopens as an editable design without inheriting source identity")
        page.locator('.card-heading .sketch-open[data-block="a"]').click()
        page.wait_for_function("() => document.querySelector('#sketch-formula').textContent.length>0")
        page.screenshot(path=str(out/"cofdraw-monomer-editor.png"),full_page=True)
        page.locator('#sketch-clear').click()
        page.locator('#classic-templates').click()
        page.get_by_role('button', name='Phenyl', exact=True).click()
        page.locator('#classic-shelf').wait_for(state='hidden')
        sheet=page.locator('#monomer-canvas')
        cursor=sheet.evaluate("e=>{const p=new DOMPoint(220,130).matrixTransform(e.getScreenCTM());return [p.x,p.y]}")
        page.mouse.click(*cursor)
        page.locator('#sketch-formula').filter(has_text='C6H6').wait_for()
        page.locator('#classic-templates').click()
        page.locator('#classic-category').select_option('Functional groups')
        page.get_by_role('button', name='Formyl', exact=True).click()
        page.locator('#classic-shelf').wait_for(state='hidden')
        page.mouse.click(*cursor)
        page.locator('#sketch-formula').filter(has_text='C7H6O').wait_for()
        checks.append('Phenyl stays at the placement point and accepts a formyl template at that atom')
        page.locator('#sketch-smiles').fill('Nc1ccc(O)cc1')
        with page.expect_response('**/api/drawing/monomer'):
            page.locator('#sketch-load-smiles').click()
        page.locator('#monomer-canvas text').filter(has_text='NH2').wait_for()
        assert 'OH' in page.locator('#monomer-canvas').text_content()
        checks.append('Validated aminophenol shows the terminal amine and hydroxyl hydrogens in the sketch')
        page.locator("#sketch-clear").click()
        page.locator('[data-tool="ring-aromatic"]').click()
        canvas=page.locator("#monomer-canvas").bounding_box();page.mouse.click(canvas["x"]+canvas["width"]*.5,canvas["y"]+canvas["height"]*.5)
        midpoint=page.evaluate("""() => {const b=sketcher.graph.bonds[0],a=sketcher.graph.atoms[b.a],c=sketcher.graph.atoms[b.b];
            const p=sketcher.screen({x:(a.x+c.x)/2,y:(a.y+c.y)/2});const v=new DOMPoint(...p).matrixTransform(sketcher.svg.getScreenCTM());return [v.x,v.y];}""")
        page.mouse.move(*midpoint)
        assert page.locator('#monomer-canvas polyline').count()==1
        page.mouse.click(*midpoint)
        graph=page.evaluate("sketcher.exportGraph()")
        fused=page.request.post(base_url+'/api/drawing/monomer',data={'graph':graph})
        assert fused.ok and fused.json()['formula']=='C10H8' and '.' not in fused.json()['smiles']
        page.locator('#monomer-canvas').press('Control+z')
        checks.append("Fusing benzene onto a bond previews the outside ring, produces connected naphthalene and undoes as one edit")
        with page.expect_response("**/api/drawing/monomer") as checked:page.locator("#sketch-use").click()
        assert checked.value.status==200,checked.value.text()
        assert checked.value.json()["formula"]=="C6H6"
        assert "Monomers changed" in page.locator("#dirty-status").inner_text()
        checks.append("Interactive benzene ring drawing validates and replaces one monomer")
        # Assigned stereo remains byte-equivalent in graph fields when a vertex moves.
        page.locator('.card-heading .sketch-open[data-block="a"]').click()
        page.locator("#sketch-smiles").fill("N[C@@H](C)C(=O)O")
        with page.expect_response("**/api/drawing/monomer") as stereo:page.locator("#sketch-load-smiles").click()
        assert stereo.value.status==200
        original=stereo.value.json()["smiles"]
        page.locator('[data-tool="bond"][data-value="1"]').click()
        canvas=page.locator("#monomer-canvas").bounding_box();page.mouse.click(canvas["x"]+20,canvas["y"]+20)
        assert "assigned stereochemistry" in page.locator("#sketch-message").inner_text()
        page.locator("#sketch-use").click()
        assert page.locator("#smiles-a").input_value()==original
        checks.append("Stereo graph retained; unsupported stereo connectivity edit is blocked explicitly")
        page.locator('.card-heading .sketch-open[data-block="a"]').click()
        page.wait_for_function("() => document.querySelector('#sketch-formula').textContent.length>0")
        page.locator("#sketch-clear").click()
        page.locator("#sketch-elements").click()
        page.locator("#element-search").fill("Zn")
        page.get_by_role("button",name="Zinc (Zn)",exact=True).click()
        canvas = page.locator("#monomer-canvas").bounding_box()
        cx, cy = canvas["x"]+canvas["width"]/2, canvas["y"]+canvas["height"]/2
        page.mouse.click(cx+65,cy+15)
        page.locator("#atom-charge").fill("2")
        page.locator("#atom-charge").press("Tab")
        page.locator('[data-tool="atom"][data-value="N"]').click()
        page.mouse.click(cx-65,cy)
        page.get_by_role("button",name="Coordination bond",exact=True).click()
        page.mouse.move(cx-65,cy);page.mouse.down();page.mouse.move(cx+65,cy+15,steps=6);page.mouse.up()
        arrow = page.locator("#monomer-canvas .sketch-dative")
        original_span = float(arrow.get_attribute("x2"))-float(arrow.get_attribute("x1"))
        with page.expect_response("**/api/drawing/monomer") as zinc:
          page.keyboard.press("Control+Shift+K")
        assert zinc.value.status==200 and zinc.value.json()["smiles"]=="[NH3]->[Zn+2]"
        page.wait_for_function("() => document.querySelector('#sketch-message').textContent.startsWith('Layout cleaned')")
        clean_span = float(arrow.get_attribute("x2"))-float(arrow.get_attribute("x1"))
        assert abs(clean_span-original_span)>10
        draft_graph = zinc.value.request.post_data_json["graph"]
        page.locator("#sketch-undo").click()
        with page.expect_response("**/api/drawing/science/solubility") as metal_screen:
          with page.expect_response("**/api/drawing/monomer") as undo_checked:
            page.locator("#sketch-solubility").click()
        assert undo_checked.value.request.post_data_json["graph"]==draft_graph
        assert all(row["status"]=="unavailable" for row in metal_screen.value.json()["solvents"])
        page.get_by_role("button",name="Close solubility screening").click()
        page.locator("#sketch-cancel").click()
        checks.append("Periodic picker, formal charge, metal coordination, Ctrl+Shift+K with undo, and graph-bound screening")
        page.locator("#linkage").select_option("beta_ketoenamine")
        page.locator("#cof-preset").select_option("tppa-pore-v1")
        keto=build()
        assert keto.value.status==200,keto.value.text();assert keto.value.json()["audit"]["building_block_instances"]==12
        checks.append("Actual twelve-block Tp ketoenamine patch")
        # Catalog contracts, exact-site template filtering, and new styles are
        # exercised only when exposed by the engine under test.
        if library.get("templates"):
          page.locator("#linkage").select_option("imine")
          page.wait_for_function("() => document.querySelector('#groups-a').textContent.length>0")
          assert page.locator("#template-field").is_visible()
          assert page.locator("#role-b").inner_text()
          options=page.locator("#preset-b option").all_text_contents()
          assert any("—" in text for text in options),options
          checks.append("Expanded catalog shows full monomer names, reactive groups and chosen template")
          unsupported=next((key for key,value in library["linkages"].items() if value.get("supported") is False),None)
          if unsupported:
            page.locator("#linkage-search").fill(unsupported.replace("_"," "))
            page.locator("#linkage").select_option(unsupported)
            assert page.locator("#build").is_disabled()
            assert page.locator("#chemistry-note").inner_text()
            checks.append("Unsupported linkage is searchable and explicitly prevents generation")
            page.locator("#linkage-search").fill("")
          # Select a catalog 4-connected template using independently cataloged
          # monomers; expected site counts come from the actual local inspector.
          square=next((item for item in library.get("cof_presets",[]) if item.get("template_id") in ("sql_4_2","sql_4_4")),None)
          if square:
            page.locator("#linkage").select_option(square["linkage"])
            page.locator("#cof-preset").select_option(square["id"])
            page.wait_for_function("() => !document.querySelector('#build').disabled")
            square_render=build()
            assert square_render.value.status==200,square_render.value.text()
            square_result=square_render.value.json()
            assert square_result["audit"]["building_block_instances"]==square["instances"]
            if not page.locator('.identity-more').evaluate('element=>element.open'):
              page.locator('.identity-more > summary').click()
            assert page.locator("#identity-key").inner_text()=="Unassigned design"
            page.locator("#framework-image").evaluate("i=>i.decode()")
            page.screenshot(path=str(out/"cofdraw-four-connected.png"),full_page=True)
            checks.append("Actual four-connected catalog monomers render their explicitly chosen template")
          else:
            initial=next(item for item in library["cof_presets"] if item["linkage"]=="imine")
            page.locator("#linkage").select_option("imine")
            page.locator("#cof-preset").select_option(initial["id"])
            restored=build()
            assert restored.value.status==200,restored.value.text()
          expanded=next((item for item in library.get("cof_presets",[]) if item.get("linkage") in ("acylhydrazone","keto_hydrazide","azine")),None)
          if expanded:
            page.locator("#linkage").select_option(expanded["linkage"])
            page.locator("#cof-preset").select_option(expanded["id"])
            page.wait_for_function("() => !document.querySelector('#build').disabled")
            expanded_render=build()
            assert expanded_render.value.status==200,expanded_render.value.text()
            page.locator("#framework-image").evaluate("i=>i.decode()")
            page.screenshot(path=str(out/"cofdraw-expanded-linkage.png"),full_page=True)
            checks.append(f"Actual {expanded['linkage']} monomer pair rendered")
          for category in ("layered","schematic"):
            style=next((item for item in library.get("styles",[]) if item.get("category")==category),None)
            if style:
              with page.expect_response("**/api/drawing/render") as presentation:page.locator(f'[data-style="{style["id"]}"]').click()
              assert presentation.value.status==200,presentation.value.text()
              page.locator("#representation-note").wait_for(state="visible")
              page.locator("#framework-image").evaluate("i=>i.decode()")
              page.screenshot(path=str(out/f"cofdraw-{category}.png"),full_page=True)
              checks.append(f"Actual {category} presentation rendered with explicit scope notice")
        page.locator(".source-library").evaluate("e=>e.open=true")
        page.locator("#source-select").select_option("tfb-tab-pore-v1")
        with page.expect_response("**/api/drawing/render") as source:page.locator("#load-source").click()
        assert source.value.status==200,source.value.text();assert source.value.json()["audit"]["cofex_string"]
        if not page.locator('.identity-more').evaluate('element=>element.open'):
          page.locator('.identity-more > summary').click()
        assert "COFEXM-" in page.locator("#identity-key").inner_text()
        assert page.locator("#linkage").input_value()=="imine"
        assert "Tp plus" not in page.locator("#chemistry-note").inner_text()
        assert page.locator("#export-note").is_visible()==("cdxml" not in source.value.json()["downloads"])
        page.locator("#framework-image").evaluate("i=>i.decode()")
        page.screenshot(path=str(out/"cofdraw-source.png"),full_page=True)
        checks.append("Source recipe displays qualified COFex ID and full notation")
        page.locator("#geometry-options > summary").click()
        page.locator(".layout-details > summary").click()
        coordinate=page.locator("#layout-rows input").first;coordinate.fill(str(float(coordinate.input_value())+.1));coordinate.press("Tab")
        assert page.locator("#identity-key").inner_text()=="Unassigned design"
        with page.expect_response("**/api/drawing/render") as edited:page.locator("#render-changes").click()
        assert edited.value.status==200,edited.value.text();assert edited.value.json()["audit"]["cofex_id"] is None
        checks.append("Source coordinate edits strip active reference identity before re-render")
        page.set_viewport_size({"width":390,"height":844});page.locator("#geometry-options").evaluate("e=>e.open=false")
        page.screenshot(path=str(out/"cofdraw-mobile.png"),full_page=True)
        assert page.evaluate("() => document.documentElement.scrollWidth<=innerWidth+1")
        checks.append("390 px layout has no horizontal overflow")
        # Exercise the omission UI with a declared fixture after actual renders.
        # This is not an actual failed chemical export or another model result.
        omitted=edited.value.json()
        omitted["downloads"].pop("cdxml",None)
        omitted["audit"]["chemical_exports"]={"cdxml":"omitted_validation_or_writer_failure","reason":"Fixture: stereo tables did not roundtrip"}
        page.route("**/api/drawing/render",lambda route:route.fulfill(status=200,json=omitted))
        with page.expect_response("**/api/drawing/render"):page.locator('[data-style="standard"]').click()
        page.locator("#export-note").wait_for(state="visible")
        assert "Fixture: stereo tables did not roundtrip" in page.locator("#export-note").inner_text()
        assert page.locator("#exports button").get_by_text("CDXML",exact=True).count()==0
        checks.append("Explicit fixture: omitted CDXML reason visible; unsupported download absent")
        assert not errors,errors
        browser.close()
    except Exception as exc:
      failure=str(exc);traceback.print_exc()
      if page:
        try:page.screenshot(path=str(out/"cofdraw-failure.png"),full_page=True)
        except Exception:pass
      raise
    finally:
      if server is not None:
        server.should_exit=True;worker.join(timeout=10)
      report={"checks":checks,"javascript_errors":errors,"api_errors":api_errors,"failure":failure,"data_dir":str(args.data_dir) if args.data_dir else None,"base_url":base_url,"server_ownership":"external_untouched" if args.base_url else "owned_python_thread"}
      (out/"cofdraw-browser-report.json").write_text(json.dumps(report,indent=2),encoding="utf-8")
      print(json.dumps(report,indent=2))

if __name__=="__main__":main()
