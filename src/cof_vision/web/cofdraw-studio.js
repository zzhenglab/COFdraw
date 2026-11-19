/* Local structure laboratory and direct geometry editing. No remote assets. */
"use strict";
(()=>{
  const S=window.COFdrawStudio, $=id=>document.getElementById(id), NS='http://www.w3.org/2000/svg';
  const copy=value=>JSON.parse(JSON.stringify(value));
  let scienceBusy=false, currentStructure=null, geometryMode=false, geometryDrag=null, geometryHistory=[], lastRenderedRecipe=null;
  let revision=0, activeJob=null, origin=null, basePeriodic=null, hypothesis=null, drawingKey='', currentLayer=null, comparisonParent=null;
  const fmt=(x,d=2)=>Number.isFinite(x)?Number(x).toLocaleString(undefined,{maximumFractionDigits:d}):'—';
  const stackLabels={AA:'AA · aligned',AB:'AB · alternating',ABC:'ABC · three-step',slipped:'Slipped · progressive',twisted:'Twisted · rotated'};
  function status(message,error=false){$('science-status').textContent=message;$('science-status').classList.toggle('error',error);}
  function signature(){return JSON.stringify([S.state.sourceId,S.state.recipe,S.state.inputsDirty,S.state.jsonPending,$('smiles-a').value,$('smiles-b').value]);}
  function drawingName(){if(S.state.recipe?.assembly)return "Current design · "+Object.values(S.state.recipe.components).map(c=>c.label).join(" + ");return S.state.sourceId?`Source drawing · ${S.state.sourceMetadata?.label||S.state.recipe?.label||S.state.sourceId}`:`Current design · ${$('label-a').value||'A'} + ${$('label-b').value||'B'}`;}
  function drawingReady(){return Boolean(S.state.recipe&&!S.state.inputsDirty&&!S.state.jsonPending);}
  function syncSource(){
    const valid=drawingReady();
    $('current-drawing-label').textContent=valid?drawingName():'Create or apply a drawing to generate its model.';
    for(const id of ['build-current-cell','build-patch-3d','compare-stacks']) $(id).disabled=!valid;
  }
  function invalidate(message,clearViewer=true){revision++;activeJob=null;scienceBusy=false;document.querySelector('.science-controls').classList.remove('science-busy');if(clearViewer){syncDisplayCells(null);currentStructure=null;origin=null;viewer.set({atoms:[],bonds:[]});$('structure-exports').replaceChildren();$('structure-warnings').replaceChildren();$('structure-scope').textContent=message||'';$('viewer-title').textContent='Ready for a new structure';$('viewer-provenance').textContent=message||'';}for(const id of ['porosity-results','monomer-results'])$(id).replaceChildren();window.dispatchEvent(new CustomEvent('cofdraw:structure',{detail:{invalidated:true}}));if(message)status(message);syncSource();}
  async function work(task){
    if(scienceBusy)return;const job={revision};activeJob=job;scienceBusy=true;document.querySelector('.science-controls').classList.add('science-busy');status('Calculating locally…');
    const check=()=>{if(activeJob!==job||revision!==job.revision){const e=new Error('An outdated result was discarded.');e.name='ObsoleteScienceResult';throw e;}};
    const api=async(path,payload)=>{check();const result=await S.api(path,payload);check();return result;};
    try{await task(api,check);}catch(e){if(e.name!=='ObsoleteScienceResult'&&activeJob===job)status(e.message,true);}finally{if(activeJob===job){activeJob=null;scienceBusy=false;document.querySelector('.science-controls').classList.remove('science-busy');syncSource();}}
  }
  function tab(name){if(name==='inspect'){name='model';$('science-inspect').open=true;}document.querySelectorAll('[data-science-tab]').forEach(b=>{const active=b.dataset.scienceTab===name;b.setAttribute('aria-selected',String(active));b.setAttribute('aria-controls','science-'+b.dataset.scienceTab);b.tabIndex=active?0:-1;});document.querySelectorAll('.science-pane').forEach(p=>p.hidden=p.id!==`science-${name}`);viewer.draw();window.dispatchEvent(new Event("cofdraw:science-tab"));}
  document.querySelectorAll('[data-science-tab]').forEach(b=>b.addEventListener('click',()=>tab(b.dataset.scienceTab)));
  document.querySelector('.science-tabs').addEventListener('keydown',event=>{
    const buttons=[...event.currentTarget.querySelectorAll('[role=tab]')], i=buttons.indexOf(event.target);
    if(i<0||!['ArrowLeft','ArrowRight','Home','End'].includes(event.key))return;
    const idx=event.key==='Home'?0:event.key==='End'?buttons.length-1:(i+(event.key==='ArrowRight'?1:-1)+buttons.length)%buttons.length;
    event.preventDefault();buttons[idx].click();buttons[idx].focus();
  });
  function workspace(visible){
    $('space-workspace').hidden=true;
    $('space-open').setAttribute('aria-pressed','false');
    $('science-dialog').hidden=!visible;
    document.querySelector('main.studio').hidden=visible;
    $('switch-space').textContent=visible?'2D drawing':'3D workspace';
    $('switch-space').setAttribute('aria-pressed',String(visible));
    window.scrollTo(0,0);
    if(visible)requestAnimationFrame(()=>viewer.draw());
  }
  $('open-science').addEventListener('click',()=>{
    workspace(true);tab('model');
    const activeDrawing=drawingReady();
    const otherModel=activeDrawing&&($('structure-source').value!=='drawing'||origin==='monomer');
    if(otherModel||(!currentStructure&&!origin&&!scienceBusy)){
      tab('model');
      if(!activeDrawing){
        $('structure-source').value='drawing';sourceMode();
        return;
      }
      $('structure-source').value='drawing';
      sourceMode();
      generateCurrentCell(true);
    }
    requestAnimationFrame(()=>viewer.draw());
  });
  $('close-science').addEventListener('click',()=>{workspace(false);$('switch-space').focus();});
  $('switch-space').addEventListener('click',()=>$('science-dialog').hidden?$('open-science').click():$('close-science').click());
  window.addEventListener('cofdraw:rendered',()=>{
    if(!$('science-dialog').hidden&&$('structure-source').value==='drawing'&&!origin&&!scienceBusy&&drawingReady())generateCurrentCell();
  });
  function metric(container,label,value){const row=document.createElement('div');row.className='metric-card';const title=document.createElement('span'),number=document.createElement('strong');title.textContent=label;number.textContent=value;row.append(title,number);container.append(row);}
  function note(container,text,warning=false){if(!text)return;const p=document.createElement('p');p.className=`metric-note${warning?' warning':''}`;p.textContent=text;container.append(p);}
  function details(container,data){const el=document.createElement('details'),summary=document.createElement('summary'),pre=document.createElement('pre');summary.textContent='Method, units, and diagnostics';pre.textContent=JSON.stringify(data,null,2);el.append(summary,pre);container.append(el);}
  function report(container,data,labels={}){container.replaceChildren();for(const [key,value] of Object.entries(data||{}))if(typeof value==='number'&&Number.isFinite(value))metric(container,labels[key]||key.replaceAll('_',' '),fmt(value,3));}
  function exportsFor(result){$('structure-exports').replaceChildren();for(const [ext,path] of Object.entries(result.downloads||{})){const b=document.createElement('button');b.className='button outline';b.textContent=ext.toUpperCase();b.addEventListener('click',()=>S.saveExport(path,ext));$('structure-exports').append(b);}if(result.viewer){const b=document.createElement('button');b.className='button outline';b.textContent='View PNG';b.addEventListener('click',async()=>{try{const data=$('structure-canvas').toDataURL('image/png');if(window.pywebview?.api?.save_image)await window.pywebview.api.save_image('COFdraw-3D-view.png',data);else{const a=document.createElement('a');a.href=data;a.download='COFdraw-3D-view.png';a.click();}}catch(e){status(e.message,true);}});$('structure-exports').append(b);}}
  let automaticNetworkColor=false;
  $('viewer-color').addEventListener('change',()=>{automaticNetworkColor=false;});
  function display(result,title,kind='periodic',parentLayer=null){
    if(result.structure?.provenance?.network_dimension===3&&$('viewer-color').value==='element'){
      $('viewer-color').value='network';automaticNetworkColor=true;
    }else if(automaticNetworkColor&&result.structure?.provenance?.network_dimension!==3){
      $('viewer-color').value='element';automaticNetworkColor=false;
    }
    comparisonParent=parentLayer;
    syncDisplayCells(result);
    currentStructure=result.structure||null;origin=kind;viewer.set(result.viewer||{atoms:[],bonds:[]});$('viewer-title').textContent=title;
    const provenance = currentStructure?.provenance;
    $('viewer-provenance').textContent=kind==='drawing'?`${drawingName()} · finite layers from the active recipe`:kind==='current-periodic'?`${drawingName()} · unrelaxed periodic cell`:kind==='monomer'?'Selected monomer · force-field conformer':kind==='networks'?'User-defined translated copies · interpenetration hypothesis':'Independent periodic model / imported CIF';
    $('structure-scope').textContent=result.scope||result.qualification||result.viewer?.scope||(kind==='drawing'?'Finite layers scaled from the 2D drawing; no periodic crystal or optimized stacking is inferred.':'Model coordinates in Å; inspect geometry and provenance before simulation.');
    const warnings=$('structure-warnings'),routine=[];warnings.replaceChildren();for(const warning of result.warnings||[]){if(/^(Finite illustrative layers|No interlayer covalent|Do not use this preview|Translated copies are|The requested number|Inspect contacts and relaxation|Imported\/proposed coordinates|CIF\/extended XYZ retain|PDB is a coordinate)/.test(warning))routine.push(warning);else note(warnings,warning,true);}if(routine.length){const box=document.createElement('details'),summary=document.createElement('summary');box.className='model-notes';summary.textContent=`Model notes (${routine.length})`;box.append(summary);for(const warning of routine)note(box,warning);warnings.append(box);}
    if(Number.isFinite(result.minimum_intercopy_distance_angstrom))metric($('structure-warnings'),'Closest intercopy atoms / Å',fmt(result.minimum_intercopy_distance_angstrom,3));
    if(kind==='current-periodic'&&provenance?.spacing_angstrom) {
      note(warnings,`${provenance.stacking || 'AA'} layers · ${fmt(provenance.spacing_angstrom)} Å spacing`);
    }
    if(result.viewer?.periodic===true&&!(result.viewer.bonds?.length))note($('structure-warnings'),'Viewer sticks are inferred from distances; they do not assign bond orders or validate the periodic network.');
    $('porosity-results').replaceChildren();exportsFor(result);window.dispatchEvent(new CustomEvent('cofdraw:structure',{detail:result}));
  }
  let displayedReport=null, displayedCounts=[1,1,1], repeatRevision=0;
  const cellSelects=['a','b','c'].map(axis=>$('display-'+axis));
  function showCounts(){
    cellSelects.forEach((select,i)=>select.value=String(displayedCounts[i]));
  }
  function syncDisplayCells(result){
    repeatRevision++;
    displayedReport=result?.viewer?.periodic?result.report_id:null;
    displayedCounts=result?.repeats||[1,1,1];
    showCounts();
    $('display-cells').disabled=!displayedReport;
    $('cell-view-status').textContent=displayedReport?'Display only · calculations use the source cell.':'Available for periodic structures.';
  }
  async function changeCells(){
    const job=++repeatRevision, counts=cellSelects.map(select=>Number(select.value));
    if(!displayedReport)return;
    $('cell-view-status').textContent='Updating displayed cells…';
    try{
      const result=await S.api('science/reports/'+displayedReport+'/view',{repeats:counts});
      if(job!==repeatRevision)return;
      displayedCounts=result.repeats;
      viewer.set(result.viewer,true);
      window.dispatchEvent(new Event('cofdraw:display-cells'));
      $('cell-view-status').textContent=counts.join(' × ')+' cells · '+viewer.atoms.length+' atoms displayed.';
    }catch(error){
      if(job!==repeatRevision)return;
      showCounts();$('cell-view-status').textContent=error.message;
    }
  }
  cellSelects.forEach(select=>select.addEventListener('change',changeCells));
  function drawingRequest(){if(!S.state.recipe)throw new Error('Create or open a drawing first.');if(S.state.inputsDirty||S.state.jsonPending)throw new Error('Apply your input changes first.');return S.state.sourceId?{source_recipe_id:S.state.sourceId}:{recipe:copy(S.state.recipe)};}
  const repeats=()=>['a','b','c'].map(x=>Number($(`repeat-${x}`).value));
  function patchRequest(stacking=$('patch-stacking').value){return {...drawingRequest(),layers:Number($('patch-layers').value),layer_spacing:Number($('patch-spacing').value),stacking,slip:[Number($('patch-slip-x').value),Number($('patch-slip-y').value)],twist_degrees:Number($('patch-twist').value)};}
  async function showFinite(api){
    const request=patchRequest(), result=await api('science/patch',request);
    display(result,`Current drawing · ${stackLabels[request.stacking]}`,'drawing');
  }
  function generatePatch(){return work(async api=>{await showFinite(api);status('Finite layers ready.');});}
  function cellNote(message){
    $('current-cell-result').textContent=message;
    $('current-cell-result').hidden=!message;
  }
  function generateCurrentCell(useFiniteFallback=true){
    return work(async api=>{
      const value=$('current-cell-spacing').value;
      const result=await api('science/build-periodic',{
        ...drawingRequest(), layer_spacing:value===''?null:Number(value), repeats:[2,2,2]
      });
      if(!result.available){
        cellNote(result.reason);
        if(useFiniteFallback){
          $('finite-view-options').open=true;
          await showFinite(api);
        }
        status(useFiniteFallback?'Showing finite layers. See the model panel for CIF options.':result.reason);
        return;
      }
      currentLayer=result.structure;
      const provenance=result.structure.provenance;
      const reported=provenance.kind==='reported_crystal';
      display(result,reported?provenance.model_name+(provenance.network_dimension===1?' · published ribbon packing':' · published crystal model'):'Current drawing · periodic AA','current-periodic');
      if(reported)$('viewer-provenance').textContent=drawingName()+' · published coordinates · DOI '+provenance.doi;
      cellNote(reported?provenance.qualification:`CIF ready · ${result.structure.atoms.length} atoms per unit cell · ${fmt(provenance.spacing_angstrom)} Å layer spacing`);
      $('finite-view-options').open=false;
      status('Periodic model and CIF ready. PXRD and pore geometry are available.');
    });
  }
  function clearComparison(){$('stack-comparison').replaceChildren();}
  function sourceMode(){
    const periodic=$('structure-source').value==='periodic';
    $('drawing-model-controls').hidden=periodic;
    $('periodic-model-controls').hidden=!periodic;
    invalidate(periodic?'Choose a model template or import a CIF.':'Build a periodic cell or view finite layers for this drawing.');
    cellNote('');
    clearComparison();
    currentLayer=null;
    basePeriodic=null;
    hypothesis=null;
    $('network-comparison').hidden=true;
  }
  function openCurrentCell(){
    tab('model');
    if($('structure-source').value!=='drawing'){
      $('structure-source').value='drawing';
      sourceMode();
    }
    return generateCurrentCell();
  }
  $('structure-source').addEventListener('change',sourceMode);
  $('build-current-cell').addEventListener('click',()=>generateCurrentCell());
  $('current-cell-spacing').addEventListener('input',()=>{
    currentLayer=null;
    invalidate('Layer spacing changed. Build the CIF again.',origin==='current-periodic');
    cellNote('');
  });
  $('build-patch-3d').addEventListener('click',generatePatch);
  for(const id of ['patch-layers','patch-spacing','patch-stacking','patch-slip-x','patch-slip-y','patch-twist'])$(id).addEventListener('input',()=>{
    invalidate('Layer settings changed. Generate a preview or compare stackings.',origin==='drawing');
    clearComparison();
  });
  function thumbnail(canvas,model){const c=canvas.getContext('2d'),w=canvas.width,h=canvas.height,atoms=model.atoms||[],pts=atoms.map(a=>[a.x+.45*a.z,a.y+.65*a.z]);c.fillStyle='#f3f7f8';c.fillRect(0,0,w,h);if(!pts.length)return;const lo=[0,1].map(k=>Math.min(...pts.map(p=>p[k]))),hi=[0,1].map(k=>Math.max(...pts.map(p=>p[k]))),scale=Math.min((w-18)/(hi[0]-lo[0]+1),(h-18)/(hi[1]-lo[1]+1)),point=i=>[w/2+(pts[i][0]-(lo[0]+hi[0])/2)*scale,h/2-(pts[i][1]-(lo[1]+hi[1])/2)*scale],colors=['#2e8a91','#de995e','#8279b0','#72a881'];c.lineWidth=.8;for(const bond of model.bonds||[]){const a=bond.a??bond[0],b=bond.b??bond[1],p=point(a),q=point(b);c.strokeStyle=colors[(atoms[a].layer||0)%colors.length];c.beginPath();c.moveTo(...p);c.lineTo(...q);c.stroke();}}
  $('compare-stacks').addEventListener('click',()=>work(async(api,check)=>{const requested=patchRequest(),rows=[];for(const mode of Object.keys(stackLabels))rows.push({mode,result:await api('science/patch',{...requested,stacking:mode})});check();$('stack-comparison').replaceChildren();for(const {mode,result} of rows){const button=document.createElement('button'),canvas=document.createElement('canvas'),label=document.createElement('span');button.className='stack-option';button.dataset.stacking=mode;canvas.width=180;canvas.height=120;label.textContent=stackLabels[mode];button.append(canvas,label);thumbnail(canvas,result.viewer);button.addEventListener('click',()=>{$('patch-stacking').value=mode;display(result,`Current drawing · ${stackLabels[mode]}`,'drawing');for(const b of $('stack-comparison').children)b.setAttribute('aria-pressed',String(b===button));status(`${stackLabels[mode]} comparison selected. All five models use the same drawing and layer parameters.`);});$('stack-comparison').append(button);}$('stack-comparison').firstElementChild.click();}));
  let crystalTemplates=[];
  function templateNote(){const item=crystalTemplates.find(t=>t.id===$('crystal-template').value);$('crystal-template-note').textContent=item?.description||'';}
  $('crystal-template').addEventListener('change',()=>{templateNote();invalidate('Periodic template changed. Generate it to replace the displayed model.');basePeriodic=null;hypothesis=null;$('network-comparison').hidden=true;});
  async function loadScience(){try{
    const data=await S.api('science/library'),select=$('crystal-template');
    const groups=[['Published models',data.reported_models],['3D structures · CURATED-COFs',data.curated_models],['Reticular designs',data.templates]];
    crystalTemplates=[];select.replaceChildren();
    for(const [label,entries] of groups){
      if(!entries?.length)continue;
      const group=document.createElement('optgroup');group.label=label;
      for(const t of entries){group.append(new Option(t.label,t.id));crystalTemplates.push(t);}
      select.append(group);
    }
    select.value=crystalTemplates.find(t=>t.id==='cof-300-c5')?.id||crystalTemplates[0]?.id;templateNote();
    window.dispatchEvent(new CustomEvent('cofdraw:crystal-catalog',{detail:data.curated_models||[]}));
    if(new URLSearchParams(location.search).get('workspace')==='crystal'&&!currentStructure&&!origin){
      workspace(true);tab('model');$('structure-source').value='periodic';sourceMode();
      $('crystal-template').value='hcb-imine-aa';templateNote();$('build-crystal').click();
    }
  }catch(e){status(e.message,true);}}
  window.addEventListener('cofdraw:open-crystal',e=>{
    if(!crystalTemplates.some(t=>t.id===e.detail))return;
    workspace(true);tab('model');$('structure-source').value='periodic';sourceMode();
    $('crystal-template').value=e.detail;templateNote();$('build-crystal').click();
  });
  $('build-crystal').addEventListener('click',()=>work(async api=>{const r=await api('science/model',{template_id:$('crystal-template').value,repeats:repeats()});basePeriodic={result:r,title:crystalTemplates.find(t=>t.id===$('crystal-template').value)?.label||'Periodic model'};hypothesis=null;$('network-comparison').hidden=true;display(r,basePeriodic.title);if(['reported_crystal','curated_crystal'].includes(r.structure.provenance.kind))$('viewer-provenance').textContent=(r.structure.provenance.dataset_id||'Published coordinates')+' · DOI '+r.structure.provenance.doi;status('Periodic model ready.');}));
  $('import-cif').addEventListener('click',()=>$('cif-file').click());
  $('cif-file').addEventListener('change',e=>{const file=e.target.files[0];if(file){basePeriodic=null;hypothesis=null;$('network-comparison').hidden=true;invalidate('Importing a periodic CIF…');work(async api=>{if(file.size>2000000)throw new Error('CIF limit is 2 MB.');const r=await api('science/model',{cif_text:await file.text(),repeats:repeats()});basePeriodic={result:r,title:file.name};hypothesis=null;$('network-comparison').hidden=true;display(r,file.name);status('CIF imported. Cell parameters, symmetry, and occupancies were checked.');});}e.target.value='';});
  $('build-networks').addEventListener('click',()=>work(async api=>{if(!basePeriodic?.result.structure)throw new Error('Generate a periodic model or import a CIF first.');const offsets=$('network-offsets').value.trim().split(/\r?\n/).filter(Boolean).map(row=>row.trim().split(/[\s,;]+/).map(Number));if(offsets.some(row=>row.length!==3||row.some(x=>!Number.isFinite(x))))throw new Error('Enter one fractional x, y, z offset per line.');const r=await api('science/networks',{structure:basePeriodic.result.structure,offsets,repeats:repeats()});hypothesis=r;$('network-comparison').hidden=false;$('network-view').value='hypothesis';$('viewer-color').value='network';display(r,'Translated-network hypothesis','networks');status('Network-copy hypothesis generated. Review contacts; the number of copies does not establish interpenetration.');}));
  $('network-view').addEventListener('change',()=>{invalidate('',false);if($('network-view').value==='base'&&basePeriodic)display(basePeriodic.result,basePeriodic.title);else if(hypothesis)display(hypothesis,'Translated-network hypothesis','networks');});
  $('network-offsets').addEventListener('input',()=>{invalidate('Network offsets changed. Regenerate the hypothesis.');hypothesis=null;$('network-comparison').hidden=true;});
  for(const axis of ['a','b','c'])$(`repeat-${axis}`).addEventListener('input',()=>{invalidate('Cell repeats changed. Generate the model again.');basePeriodic=null;hypothesis=null;$('network-comparison').hidden=true;});
  $('calculate-monomer').addEventListener('click',()=>work(async api=>{const monomer=S.selectedDescriptor();const r=await api('science/monomer',{smiles:monomer.smiles,conformers:Number($('descriptor-conformers').value)});display(r,`Monomer ${monomer.label}`,'monomer');const el=$('monomer-results');report(el,r.descriptors,{molecular_weight_g_mol:'Molecular weight / g mol⁻¹',exact_mass_da:'Exact mass / Da',tpsa_angstrom2:'Topological polar surface / Å²',crippen_logp:'Crippen logP',hydrogen_bond_donors:'H-bond donors',hydrogen_bond_acceptors:'H-bond acceptors',rotatable_bonds:'Rotatable bonds',formal_charge:'Formal charge'});if(r.selected_conformer?.energy_kcal_mol!==null)metric(el,'Conformer energy / kcal mol⁻¹',fmt(r.selected_conformer?.energy_kcal_mol,3));note(el,`${r.method?.embedding||'ETKDG'} · ${r.method?.force_field||'No available force field'} · ${r.selected_conformer?.converged?'converged':'not converged'}`,!r.selected_conformer?.converged);for(const warning of r.warnings||[])note(el,warning);details(el,{method:r.method,conformers:r.conformers,partial_charges:r.partial_charges});status('Monomer descriptors calculated. Conformer energies are force-field values, not quantum orbital energies.');}));
  $('calculate-pores').addEventListener('click',()=>work(async(api,check)=>{
    if(!currentStructure||currentStructure.periodic!==true)throw new Error('Generate a periodic model or import a CIF first. Finite drawing layers are not a periodic pore model.');
    const r=await window.COFdrawPores.calculate();
    check();
    if(!r)return;
    const el=$('porosity-results');el.replaceChildren();el.dataset.engine=r.engine;el.dataset.reportId=r.report_id||'';
    const range=pair=>Array.isArray(pair)&&pair.length===2?`${fmt(pair[0],2)} – ${fmt(pair[1],2)}`:'Unavailable';
    metric(el,'PLD · geometric bounds / Å',range(r.pld_bounds_angstrom));
    metric(el,'LCD · geometric bounds / Å',range(r.lcd_bounds_angstrom));
    note(el,'PLD is the largest passing-sphere diameter; LCD is the largest included-sphere diameter. Ranges bound the chosen hard-sphere geometry at this grid resolution, not experimental uncertainty.');
    const unresolved=r.connectivity?.status==='unresolved_at_this_resolution',surface=r.method?.surface||{};
    const unresolvedArea=unresolved||(surface.exposed_samples>0&&surface.connected_samples===0&&r.connectivity?.certified_winding);
    metric(el,'Accessible area / m² g⁻¹',unresolvedArea?'Unresolved':Number.isFinite(r.transport_accessible_surface_m2_g)?fmt(r.transport_accessible_surface_m2_g,1):'Unavailable');
    metric(el,'Accessible probe-center volume / cm³ g⁻¹',unresolved?'Unresolved':fmt(r.transport_accessible_probe_center_volume_cm3_g,3));
    note(el,`Accessible-volume bounds: ${range(r.transport_accessible_probe_center_volume_bounds_cm3_g)} cm³ g⁻¹. Accessibility requires a periodic winding path for the selected probe.`,unresolved);
    metric(el,'Density / g cm⁻³',fmt(r.density_g_cm3,3));
    metric(el,'Point-probe void fraction',`${fmt(r.void_fraction*100,1)}%`);
    metric(el,'Geometric void volume / cm³ g⁻¹',fmt(r.geometric_void_volume_cm3_g,3));
    metric(el,'Local exposed area / m² g⁻¹',Number.isFinite(r.local_probe_exposed_surface_m2_g)?fmt(r.local_probe_exposed_surface_m2_g,1):'Unavailable');
    note(el,'Local exposed area includes closed cavities. Accessible area excludes samples without a certified channel connection; neither area is BET.');
    metric(el,'Certified channel dimensionality',unresolved?'Unresolved':String(r.connectivity?.maximum_periodic_rank??'Unavailable'));
    note(el,`Periodic grid ${r.method.grid_shape.join(' × ')} (${r.method.grid_points.toLocaleString()} points); effective axis steps ${r.method.effective_axis_steps_angstrom.map(x=>fmt(x,3)).join(', ')} Å. Surface: ${surface.samples??0} deterministic Halton samples.`,Boolean(r.method.grid_capped));
    if(surface.available===false)note(el,surface.reason,true);
    if(r.convergence?.requested){if(r.convergence.performed){metric(el,'Refined PLD bounds / Å',range(r.convergence.refined.pld_bounds_angstrom));metric(el,'Refined LCD bounds / Å',range(r.convergence.refined.lcd_bounds_angstrom));note(el,`Second grid ${r.convergence.refined_grid_shape.join(' × ')}. Compare both results; a second grid alone does not establish convergence. Surface sampling is unchanged.`);}else note(el,r.convergence.reason||'A finer grid was unavailable.',true);}
    for(const warning of r.warnings||[])if(/coarsened|budget|No transport|no exposed surface|No explicit hydrogen/i.test(warning))note(el,warning,true);
    const {pore_map,...diagnostics}=r;
    details(el,diagnostics);
    if(r.downloads?.json){const b=document.createElement('button');b.className='button outline';b.textContent='Save geometric-properties report';b.addEventListener('click',()=>S.saveExport(r.downloads.json,'json'));el.append(b);}
    window.dispatchEvent(new CustomEvent('cofdraw:geometry',{detail:r}));
    status('Geometric properties complete. Review diameter bounds, effective grid, and accessibility qualifications.');
  }));
  $('fit-bet').addEventListener('click',()=>work(async api=>{const rows=[];
    $('bet-results').replaceChildren();
    for(const [i,raw] of $('bet-data').value.split(/\r?\n/).entries()){
      const line=raw.trim();if(!line||line.startsWith('#'))continue;
      if(!rows.length && /^(p\s*\/\s*p[0₀o]|relative[ _]pressure)/i.test(line))continue;
      const fields=line.split(line.includes(',')?',':line.includes(';')?';':/\s+/).map(field=>field.trim());
      if(fields.length!==2||fields.some(field=>!field||!Number.isFinite(Number(field))))throw new Error('BET data line '+(i+1)+' needs two numbers: relative pressure and uptake in mmol/g.');
      rows.push(fields.map(Number));
    }
    const r=await api('science/bet',{relative_pressures:rows.map(r=>r[0]),uptake_mmol_g:rows.map(r=>r[1]),fit_min:Number($('bet-min').value),fit_max:Number($('bet-max').value),cross_section_nm2:Number($('bet-cross-section').value)});const el=$('bet-results');report(el,r,{surface_area_m2_g:'BET area / m² g⁻¹',bet_area_m2_g:'BET area / m² g⁻¹',monolayer_capacity_mmol_g:'Monolayer capacity / mmol g⁻¹',C:'BET C',r_squared:'R²',points_used:'Fit points'});note(el,r.all_basic_checks_pass?'Basic consistency checks pass.':'One or more basic consistency checks fail. Do not report this fit as a validated BET area.',!r.all_basic_checks_pass);for(const warning of r.warnings||[])note(el,warning,true);details(el,r);if(r.downloads?.json){const b=document.createElement('button');b.className='button outline';b.textContent='Save BET report';b.addEventListener('click',()=>S.saveExport(r.downloads.json,'json'));el.append(b);}status('BET fit calculated from your supplied isotherm. Check the consistency criteria before reporting the area.');}));

  // The 3D canvas displays actual supplied coordinates. It never infers a cell.
  const elementColors={C:'#73818c',H:'#e8edf0',N:'#5597dc',O:'#e76769',S:'#e7c356',P:'#d99a55',B:'#d5ac85',F:'#72cfa4',Cl:'#66bb8c',Br:'#b88b65',I:'#a77ac0',Si:'#bcae9a',Cu:'#d5a477',Zn:'#939eb7',Fe:'#cb895b','*':'#b8bbc3'};
  const radii={H:.31,C:.76,N:.71,O:.66,F:.57,P:1.07,S:1.05,Cl:1.02,Br:1.2,I:1.39,B:.84,Si:1.11,Cu:1.32,Zn:1.22,Fe:1.32};
  const vdw={H:1.2,C:1.7,N:1.55,O:1.52,F:1.47,P:1.8,S:1.8,Cl:1.75,Br:1.85,I:1.98,B:1.92,Si:2.1,Cu:1.4,Zn:1.39,Fe:1.4};
  const groupPalette=['#28949c','#d89755','#8475b4','#6c9b68','#c97589','#5d8abe','#ae9666','#967b9f'];
  function viewerColors(atoms){const mode=$('viewer-color').value,root=$('viewer-legend'),colors=new Map(),groups=new Map();root.replaceChildren();if(mode==='element'){for(const a of atoms)colors.set(a.index,elementColors[a.element]||'#a59bc0');
      const keys=new Map(atoms.filter(a=>a.element!=='*' && !(viewer.hideHydrogen && a.element==='H')).map(a=>[(a.isotope||'')+a.element,colors.get(a.index)]));
      for(const [label,color] of keys){const chip=document.createElement('span'),swatch=document.createElement('i');swatch.style.background=color;chip.append(swatch,label);root.append(chip);}return colors;}const value=a=>mode==='network'?(a.network??a.inferred_component):a[mode];for(const a of atoms){const key=value(a);if(key!==undefined&&!groups.has(key))groups.set(key,groupPalette[groups.size%groupPalette.length]);colors.set(a.index,groups.get(key)||elementColors[a.element]||'#a59bc0');}for(const [key,color]of groups){const chip=document.createElement('span'),swatch=document.createElement('i');swatch.style.background=color;chip.append(swatch,document.createTextNode(mode==='component'?`Block ${key}`:mode==='network'&&atoms.every(a=>a.network===undefined)?`Inferred component ${Number(key)+1}`:`${mode==='layer'?'Layer':'Network copy'} ${Number(key)+1}`));root.append(chip);}if(!groups.size&&atoms.length)root.textContent=`No ${mode==='component'?'building-block':mode} groups supplied; showing element colors.`;return colors;}
  function tint(hex,target,fraction){
    const n=parseInt(hex.slice(1),16);
    return 'rgb('+[n>>16,(n>>8)&255,n&255].map(v=>Math.round(v+(target-v)*fraction)).join(',')+')';
  }
  const viewer={model:null,atoms:[],bonds:[],theta:-.3,phi:.6,zoom:1,selected:[],projected:[],set(model,retainZoom=false){this.model=model;this.projected=[];this.atoms=(model.atoms||[]).map((a,i)=>({...a,index:i,x:a.x??a.position?.[0]??a.cartesian?.[0]??0,y:a.y??a.position?.[1]??a.cartesian?.[1]??0,z:a.z??a.position?.[2]??a.cartesian?.[2]??0}));this.selected=[];$('viewer-distance').textContent='';$('viewer-legend').replaceChildren();if(!retainZoom)this.zoom=1;this.bonds=Array.isArray(model.bonds)&&!(model.periodic===true&&model.bonds.length===0)?model.bonds.map(b=>Array.isArray(b)?b:[b.a,b.b]):this.inferBonds();this.draw();},inferBonds(){const bins=new Map(),bond=[],width=3;for(let i=0;i<this.atoms.length;i++){const a=this.atoms[i],base=[a.x,a.y,a.z].map(x=>Math.floor(x/width));for(let x=-1;x<=1;x++)for(let y=-1;y<=1;y++)for(let z=-1;z<=1;z++){const key=[base[0]+x,base[1]+y,base[2]+z].join(',');for(const j of bins.get(key)||[]){const b=this.atoms[j];if(a.network!==undefined&&b.network!==undefined&&a.network!==b.network)continue;if(a.layer!==undefined&&b.layer!==undefined&&a.layer!==b.layer)continue;const distance=Math.hypot(a.x-b.x,a.y-b.y,a.z-b.z),cut=1.18*((radii[a.element]||.8)+(radii[b.element]||.8));if(distance>.35&&distance<cut)bond.push([i,j]);}}const key=base.join(',');if(!bins.has(key))bins.set(key,[]);bins.get(key).push(i);}return bond;},draw(){const canvas=$('structure-canvas'),rect=canvas.getBoundingClientRect();if(rect.width<1)return;const dpr=Math.min(devicePixelRatio||1,2),w=rect.width,h=rect.height;canvas.width=Math.round(w*dpr);canvas.height=Math.round(h*dpr);const c=canvas.getContext('2d');c.scale(dpr,dpr);c.fillStyle='#f8fbfa';c.fillRect(0,0,w,h);if(!this.atoms.length){c.fillStyle='#91a59b';c.font='14px Segoe UI';c.textAlign='center';c.fillText('Your 3D structure will appear here',w/2,h/2);return;}const coords=this.atoms.map(a=>[a.x,a.y,a.z]),lo=[0,1,2].map(k=>Math.min(...coords.map(p=>p[k]))),hi=[0,1,2].map(k=>Math.max(...coords.map(p=>p[k]))),center=lo.map((v,k)=>(v+hi[k])/2);const rotate=p=>{const x=p[0]-center[0],y=p[1]-center[1],z=p[2]-center[2],x1=x*Math.cos(this.theta)+z*Math.sin(this.theta),z1=-x*Math.sin(this.theta)+z*Math.cos(this.theta);return {x:x1,y:y*Math.cos(this.phi)-z1*Math.sin(this.phi),z:y*Math.sin(this.phi)+z1*Math.cos(this.phi)};};
      const rotated=coords.map(rotate);
      const extent=Math.max(1,Math.hypot(...hi.map((v,k)=>v-lo[k]))/2);
      const perspective=Number($('camera-perspective').value)*.006/extent;
      const fading=Number($('camera-fade').value)*.008;
      function camera(r){
        const factor=1/Math.max(.25,1-perspective*r.z);
        return {x:r.x*factor,y:r.y*factor,z:r.z,factor};
      }
      const plane=rotated.map(camera),xs=plane.map(p=>p.x),ys=plane.map(p=>p.y);
      const xlo=Math.min(...xs),xhi=Math.max(...xs),ylo=Math.min(...ys),yhi=Math.max(...ys);
      const scale=Math.min(w*.82/(xhi-xlo+4),h*.82/(yhi-ylo+4))*this.zoom;
      const screen=r=>({...r,x:w/2+(r.x-(xlo+xhi)/2)*scale,y:h/2-(r.y-(ylo+yhi)/2)*scale});
      const project=p=>screen(camera(rotate(p)));
      this.projected=plane.map(screen);
      const mode=$('viewer-mode').value,colors=viewerColors(this.atoms),color=a=>colors.get(a.index)||elementColors[a.element]||'#a59bc0';
      const layers=new Map();
      if(mode==='mixed')this.atoms.forEach((a,i)=>{
        if(a.layer===undefined)return;
        const stats=layers.get(a.layer)||[0,0];stats[0]+=this.projected[i].z;stats[1]++;layers.set(a.layer,stats);
      });
      const frontLayer=[...layers].sort(([,a],[,b])=>b[0]/b[1]-a[0]/a[1])[0]?.[0];
      const featured=a=>a.layer!==undefined?a.layer===frontLayer:a.network!==undefined?a.network===0:a.component!==undefined?a.component===this.atoms[0]?.component:!['C','H'].includes(a.element);
      if(this.model?.cell){c.strokeStyle='#9cb4aa';c.lineWidth=.8;c.setLineDash([4,4]);const cell=this.model.cell,corner=[];for(let i=0;i<8;i++)corner.push(project([0,1,2].map(k=>((i&1)?cell[0][k]:0)+((i&2)?cell[1][k]:0)+((i&4)?cell[2][k]:0))));for(let i=0;i<8;i++)for(const bit of [1,2,4])if(!(i&bit)){const a=corner[i],b=corner[i|bit];c.beginPath();c.moveTo(a.x,a.y);c.lineTo(b.x,b.y);c.stroke();}c.setLineDash([]);}
      const objects=this.atoms.map((a,i)=>({kind:'atom',index:i,z:this.projected[i].z})).filter(o=>!this.hideHydrogen||this.atoms[o.index].element!=='H');
      const zmin=Math.min(...this.projected.map(p=>p.z)),zmax=Math.max(...this.projected.map(p=>p.z));
      if(mode!=='spacefill'){
        for(const [a,b] of this.bonds){
          const p=this.projected[a],q=this.projected[b];
          if(!p||!q||this.hideHydrogen&&(this.atoms[a].element==='H'||this.atoms[b].element==='H'))continue;
          const mid={x:(p.x+q.x)/2,y:(p.y+q.y)/2,z:(p.z+q.z)/2};
          for(const [index,tip] of [[a,p],[b,q]])objects.push({kind:'bond',index,tip,mid,z:(tip.z+mid.z)/2});
        }
      }
      window.COFdrawPores?.addFaces(objects,project,scale);
      objects.sort((a,b)=>a.z-b.z);
      for(const o of objects){
        const depth=Math.max(0,Math.min(1,(zmax-o.z)/Math.max(.01,zmax-zmin)));
        const visibility=1-fading*depth;
        if(o.kind==='pore'){window.COFdrawPores.paint(c,o,visibility);continue;}
        const a=this.atoms[o.index],p=this.projected[o.index],base=color(a);
        const localScale=scale*p.factor;
        if(o.kind==='bond'){
          const width=mode==='ballstick'?Math.max(1.15,localScale*.17):Math.max(.85,localScale*.075);
          c.lineWidth=width;c.lineCap='round';c.strokeStyle=base;c.globalAlpha=(1-.30*depth)*visibility;
          c.beginPath();c.moveTo(o.tip.x,o.tip.y);c.lineTo(o.mid.x,o.mid.y);c.stroke();
        }else{
          const skeletal=mode==='wireframe'||mode==='mixed'&&!featured(a);
          const r=skeletal?Math.max(.85,localScale*.06):Math.max(1.2,localScale*(vdw[a.element]||1.7)*(mode==='ballstick'?.25:mode==='halfspace'?.5:1));
          c.globalAlpha=visibility;c.beginPath();c.arc(p.x,p.y,r,0,Math.PI*2);
          if(skeletal){c.fillStyle=base;}
          else{
            const g=c.createRadialGradient(p.x-r*.34,p.y-r*.38,r*.03,p.x,p.y,r*1.07);
            g.addColorStop(0,tint(base,255,.61));g.addColorStop(.34,tint(base,255,.20+depth*.12));
            g.addColorStop(.75,base);g.addColorStop(1,tint(base,20,.43));c.fillStyle=g;
          }
          c.fill();if(this.selected.includes(o.index)){c.strokeStyle='#ec9e30';c.lineWidth=2;c.stroke();}
        }
      }
      c.globalAlpha=1;
      if(this.selected.length===2){const [a,b]=this.selected.map(i=>this.projected[i]);c.strokeStyle='#d98c22';c.lineWidth=2;c.setLineDash([4,3]);c.beginPath();c.moveTo(a.x,a.y);c.lineTo(b.x,b.y);c.stroke();c.setLineDash([]);}c.fillStyle='#82988d';c.font='10px Segoe UI';c.fillText(`${this.atoms.length} atoms · Å`,13,h-14);
    }};
  $('viewer-preset').addEventListener('change',()=>{
    const key=$('viewer-preset').value;
    const presets={pore:[0,0,'spacefill'],layers:[-.06,1.38,'ballstick'],oblique:[-.36,.72,'mixed']};
    if(key==='projection'&&currentStructure){
      $('viewer-mode').value='wireframe';viewer.hideHydrogen=true;$('viewer-hide-h').checked=true;
      $('viewer-axis').value='c';$('viewer-axis').dispatchEvent(new Event('change'));
      viewer.zoom=1;viewer.selected=[];$('viewer-distance').textContent='';
      $('viewer-preset').value='';viewer.draw();return;
    }
    if(!presets[key])return;
    [viewer.theta,viewer.phi,$('viewer-mode').value]=presets[key];
    viewer.hideHydrogen=true;$('viewer-hide-h').checked=true;
    $('viewer-axis').value='free';viewer.zoom=1;viewer.selected=[];$('viewer-distance').textContent='';
    viewer.draw();$('viewer-preset').value='';
  });
  let rotation=null;
  const canvas=$('structure-canvas');
  function clearDistance(){
    viewer.selected=[];$('viewer-distance').textContent='';viewer.draw();
  }
  $('viewer-measure').addEventListener('change',clearDistance);
  canvas.addEventListener('pointerdown',event=>{
    if(rotation||event.button!==0)return;
    rotation={pointer:event.pointerId,x:event.clientX,y:event.clientY,moved:0};
    canvas.setPointerCapture(event.pointerId);
  });
  canvas.addEventListener('pointermove',event=>{
    if(rotation?.pointer!==event.pointerId)return;
    const dx=event.clientX-rotation.x,dy=event.clientY-rotation.y;
    if(dx||dy)$('viewer-axis').value='free';
    viewer.theta+=dx*.008;viewer.phi+=dy*.008;rotation.moved+=Math.abs(dx)+Math.abs(dy);
    rotation.x=event.clientX;rotation.y=event.clientY;viewer.draw();
  });
  canvas.addEventListener('pointerup',event=>{
    if(rotation?.pointer!==event.pointerId)return;
    const click=rotation.moved<6;rotation=null;
    if(!click||!$('viewer-measure').checked)return;
    const rect=canvas.getBoundingClientRect(), x=event.clientX-rect.left,y=event.clientY-rect.top;
    let best=-1,dist=15;
    viewer.projected.forEach((p,i)=>{
      if(viewer.hideHydrogen&&viewer.atoms[i].element==='H')return;
      const d=Math.hypot(p.x-x,p.y-y);
      if(d<dist-.05||(best>=0&&Math.abs(d-dist)<=.05&&p.z>viewer.projected[best].z)){best=i;dist=d;}
    });
    if(best<0)return;
    if(viewer.selected.length===2)viewer.selected=[];
    if(viewer.selected[0]===best){$('viewer-distance').textContent='Choose a different atom.';return;}
    viewer.selected.push(best);
    if(viewer.selected.length===2){
      const [a,b]=viewer.selected.map(i=>viewer.atoms[i]);
      const distance=Math.hypot(a.x-b.x,a.y-b.y,a.z-b.z);
      $('viewer-distance').textContent=(a.isotope||'')+a.element+(a.index+1)+' ↔ '+(b.isotope||'')+b.element+(b.index+1)+': '+fmt(distance,3)+' Å · distance between displayed coordinates';
    }else $('viewer-distance').textContent='Select a second atom.';
    viewer.draw();
  });
  canvas.addEventListener('lostpointercapture',e=>{if(rotation?.pointer===e.pointerId)rotation=null;});
  canvas.addEventListener('wheel',e=>{e.preventDefault();viewer.zoom=Math.max(.3,Math.min(4,viewer.zoom*Math.exp(-e.deltaY*.001)));viewer.draw();},{passive:false});
  function updateCamera(){
    for(const id of ['camera-fade','camera-perspective']){
      document.querySelector('output[for="'+id+'"]').value=$(id).value+'%';
    }
    viewer.draw();
  }
  for(const id of ['camera-fade','camera-perspective'])$(id).addEventListener('input',updateCamera);
  $('viewer-mode').addEventListener('change',()=>viewer.draw());
  $('viewer-color').addEventListener('change',()=>viewer.draw());
  $('viewer-reset').addEventListener('click',()=>{$('camera-perspective').value=0;$('camera-fade').value=0;updateCamera();viewer.theta=-.3;viewer.phi=.6;viewer.zoom=1;$('viewer-axis').value='free';clearDistance();});
  new ResizeObserver(()=>viewer.draw()).observe(canvas);

  // Edit the exact displayed coordinates; graph connectivity is not rewritten.
  function remember(){if(!S.state.recipe)return;geometryHistory.push({recipe:copy(S.state.recipe),sourceId:S.state.sourceId,metadata:copy(S.state.sourceMetadata)});if(geometryHistory.length>20)geometryHistory.shift();$('undo-geometry').disabled=false;}
  function geometryMessage(message){$('geometry-message').textContent=message;$('geometry-message').hidden=!message;}
  function drawHandles(){const overlay=$('geometry-overlay');overlay.replaceChildren();overlay.toggleAttribute('hidden',!geometryMode);$('render-geometry').hidden=!geometryMode;if(!geometryMode)return;const layout=S.state.result?.audit?.drawing_style?.layout;if(!layout||layout.display_projection||lastRenderedRecipe!==S.state.recipe){geometryMessage('Use a flat chemical style to adjust atom positions.');return;}overlay.setAttribute('viewBox',`0 0 ${layout.canvas_size.join(' ')}`);overlay.setAttribute('preserveAspectRatio','xMidYMid meet');for(const [i,p] of (layout.atom_pixels||[]).entries()){const circle=document.createElementNS(NS,'circle');circle.setAttribute('cx',p[0]);circle.setAttribute('cy',p[1]);circle.setAttribute('r',Math.max(4,layout.canvas_size[0]*.004));circle.classList.add('geometry-handle');circle.dataset.atom=i;overlay.append(circle);}geometryMessage('Drag an atom to adjust its position, then Render changes. Connectivity stays fixed. Undo restores the previous layout.');}
  $('edit-geometry').addEventListener('click',()=>{if(!S.state.result||S.state.inputsDirty||S.state.jsonPending||lastRenderedRecipe!==S.state.recipe)return S.notice('Render the current drawing before adjusting its atoms.',true);geometryMode=!geometryMode;$('edit-geometry').setAttribute('aria-pressed',String(geometryMode));drawHandles();});
  $('geometry-overlay').addEventListener('pointerdown',e=>{const target=e.target.closest('[data-atom]');if(!target||S.state.busy||S.state.inputsDirty||S.state.jsonPending||lastRenderedRecipe!==S.state.recipe)return;remember();geometryDrag={index:Number(target.dataset.atom),layout:copy(S.state.result.audit.drawing_style.layout),circle:target};e.currentTarget.setPointerCapture(e.pointerId);e.preventDefault();});
  $('geometry-overlay').addEventListener('pointermove',e=>{if(!geometryDrag)return;const overlay=e.currentTarget,p=new DOMPoint(e.clientX,e.clientY).matrixTransform(overlay.getScreenCTM().inverse()),{index,layout,circle}=geometryDrag,origin=layout.atom_pixels[index],world=layout.atom_coordinates[index],ratio=layout.pixels_per_unit;const xy=[world[0]+(p.x-origin[0])/ratio,world[1]-(p.y-origin[1])/ratio];if(xy.some(v=>!Number.isFinite(v)||Math.abs(v)>500))return;const recipe=S.state.recipe;S.removeIdentity();recipe.layout_refinement='manual';if(recipe.explicit_graph){recipe.explicit_graph.atoms[index].xy=xy;}else{if(!recipe.atom_coordinates)recipe.atom_coordinates=layout.atom_keys.map((key,i)=>{const at=S.state.result.audit.atom_origin[i];return {instance:at.instance,fragment_atom_index:at.fragment_atom_index,xy:layout.atom_coordinates[i]};});const key=S.state.result.audit.atom_origin[index],row=recipe.atom_coordinates.find(a=>a.instance===key.instance&&a.fragment_atom_index===key.fragment_atom_index);if(!row)return;row.xy=xy;}circle.setAttribute('cx',p.x);circle.setAttribute('cy',p.y);circle.classList.add('active');S.state.dirty=true;S.syncJSON();S.sync();});
  for(const event of ['pointerup','pointercancel'])$('geometry-overlay').addEventListener(event,()=>{if(geometryDrag){geometryDrag=null;geometryMessage('Position adjusted. Render changes to check bonds and save the revised drawing.');}});
  $('undo-geometry').addEventListener('click',()=>{const previous=geometryHistory.pop();if(!previous)return;S.acceptRecipe(previous.recipe,previous.sourceId,previous.metadata);$('undo-geometry').disabled=!geometryHistory.length;S.render();});
  $('render-geometry').addEventListener('click',()=>S.render());
  $('clean-geometry').addEventListener('click',async()=>{if(S.state.busy)return;try{remember();const r=await S.api('layout',{...drawingRequest(),operation:'clean'});S.acceptRecipe(r.recipe);await S.render();const review=S.geometryReviewMessage(r.audit?.layout_refinement);geometryMessage(review||'Geometry cleaned as a new design. Source alignment is preserved in the original source recipe.');if(review)S.notice(`Drawing created. ${review}`);}catch(e){S.notice(e.message,true);}});
  $('measure-pore').addEventListener('click',async()=>{
    const recipe=S.state.recipe,rendered=S.state.result,key=signature();
    try{
      if(!rendered||S.state.dirty||S.state.inputsDirty||S.state.jsonPending)throw new Error('Render the current drawing before measuring its pore.');
      const result=await S.api('science/measure',drawingRequest());
      if(S.state.recipe!==recipe||S.state.result!==rendered||signature()!==key||S.state.dirty||S.state.inputsDirty||S.state.jsonPending)return;
      const el=$('pore-measurement'),title=document.createElement('strong');
      title.textContent=`Pore opening · ${fmt(result.diameter_nm,2)} nm`;
      const method=document.createElement('details'),summary=document.createElement('summary'),body=document.createElement('p');
      method.className='measurement-method';
      summary.textContent='Measurement method';
      body.textContent=result.scope;
      method.append(summary,body);
      el.replaceChildren(title,method);
      el.hidden=false;
      window.dispatchEvent(new CustomEvent('cofdraw:pore-measured',{detail:{result,audit:rendered.audit,recipe}}));
    }catch(e){if(signature()===key)S.notice(e.message,true);}
  });
  function drawSourceGallery(library){const root=$('source-gallery');root.replaceChildren();for(const source of library.source_recipes||[]){const b=document.createElement('button'),icon=document.createElement('span'),text=document.createElement('div'),small=document.createElement('small');icon.textContent='⬡';text.textContent=source.label||source.recipe_id;small.textContent=source.layout_origin==='source_coordinates'?'Source-aligned coordinates':'Qualified source drawing';text.append(small);b.append(icon,text);b.addEventListener('click',async()=>{if(S.state.busy)return;try{const r=await S.api(`sources/${encodeURIComponent(source.recipe_id)}`);S.acceptRecipe(r.recipe,source.recipe_id,r.metadata);S.showSourceMonomers(r.recipe,r.metadata);S.state.style=r.recipe.source_style?'source':'standard';S.state.styleOptions={};await S.render();}catch(e){S.notice(e.message,true);}});root.append(b);}}
  window.addEventListener('cofdraw:layoutchanged',()=>{lastRenderedRecipe=null;drawHandles();$('pore-measurement').hidden=true;});
  window.addEventListener('cofdraw:library',e=>{drawSourceGallery(e.detail);loadScience();});
  window.addEventListener('cofdraw:rendered',()=>{lastRenderedRecipe=S.state.recipe;$('pore-measurement').hidden=true;drawHandles();$('layout-details').hidden=Boolean(S.state.recipe?.explicit_graph);const audit=S.state.result?.audit,review=S.state.sourceId?'':S.geometryReviewMessage(audit?.layout_refinement),messages=[review,...(audit?.layout_diagnostics?.warnings||[])].filter(Boolean);if(messages.length||!geometryMode)geometryMessage(messages.join(' '));});
  for(const [ids,target] of [[['probe-radius','pore-samples','pore-refinement'],'porosity-results'],[['bet-data','bet-min','bet-max','bet-cross-section'],'bet-results'],[['descriptor-block','descriptor-conformers','smiles-a','smiles-b'],'monomer-results']])for(const id of ids)$(id).addEventListener('input',()=>{invalidate('',false);$(target).replaceChildren();});
  function drawingChanged(){
    const key=signature();
    syncSource();
    if(key===drawingKey)return;
    drawingKey=key;
    clearComparison();
    cellNote('');
    currentLayer=null;
    invalidate('The current drawing changed. Generate its model again.', $('structure-source').value==='drawing'||origin==='monomer'||origin==='current-periodic');
  }
  window.addEventListener('cofdraw:state',drawingChanged);
  window.addEventListener('cofdraw:layoutchanged',drawingChanged);
  drawingKey=signature();syncSource();
  window.COFdrawScience={get structure(){return currentStructure;},get finiteModel(){return !currentStructure&&origin&&viewer.model?.periodic===false?viewer.model:null;},get revision(){return revision;},get origin(){return origin;},get stackingLayer(){return comparisonParent||(origin==='current-periodic'?currentLayer||currentStructure:currentStructure);},patchRequest,drawingRequest,drawingReady,openCurrentCell,viewer,work,status,metric,note,details,display,tab,invalidate};
  if(S.state.library){drawSourceGallery(S.state.library);loadScience();}
})();
