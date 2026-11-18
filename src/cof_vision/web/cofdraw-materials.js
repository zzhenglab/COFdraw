/* Structure inspection and evidence retrieval, using the same local API as agents. */
"use strict";
(() => {
  const S=window.COFdrawStudio, C=window.COFdrawScience;
  const $=id=>document.getElementById(id);
  const number=(value,digits=3)=>Number.isFinite(value)?Number(value).toLocaleString(undefined,{maximumFractionDigits:digits}):'—';
  const text=(tag,value,className)=>{const el=document.createElement(tag);el.textContent=value;if(className)el.className=className;return el;};
  function structure(){if(!C.structure)throw new Error('Generate a periodic model or import a CIF first.');return C.structure;}
  function finite(){if(!C.finiteModel)throw new Error('Generate or select the structure to inspect first.');return C.finiteModel;}
  function syncModelState(){
    const periodic=Boolean(C.structure),available=periodic||Boolean(C.finiteModel);
    for(const id of ['inspect-model-state','pxrd-model-state','porosity-model-state']){
      const root=$(id);root.replaceChildren();root.append(text('span',periodic?'Current view: periodic atomic model with a unit cell.':available?'Current view: finite atoms from the drawing. Checks and finite-cluster scattering are available without a CIF. Density, PLD and periodic Bragg peaks require a periodic cell.':'Generate your drawing in 3D, or select a periodic model.'));
      if(!periodic){const b=text('button','Choose a periodic model','button outline');b.addEventListener('click',()=>{C.tab('model');$('structure-source').value='periodic';$('structure-source').dispatchEvent(new Event('change'));});root.append(b);}
      if(!periodic&&available&&C.origin==='drawing'){
        const find=text('button','Find a periodic design for these monomers','button outline');find.classList.add('find-matching-periodic');
        find.addEventListener('click',()=>C.work(async api=>{
          const request=C.drawingRequest(),result=await api('science/periodic-options',request),choices=text('div','','matching-periodic-options');
          root.querySelector('.matching-periodic-options')?.remove();choices.append(text('p',result.reason));
          for(const option of result.options){const button=text('button',`Use ${option.label}`,'button outline');button.dataset.templateId=option.id;button.addEventListener('click',()=>C.buildMatchedPeriodic(option,request));choices.append(button);}
          if(result.options.length)choices.append(text('small','Same precursors and chosen connectivity; new unrelaxed coordinates and unit cell.'));
          root.append(choices);C.status(result.reason);
        }));root.append(find);
      }
      if(id==='porosity-model-state'&&!periodic&&available&&C.origin==='drawing'){const b=text('button','Measure this drawing’s pore','button outline');b.id='measure-current-finite-pore';b.addEventListener('click',()=>C.work(async api=>{const result=await api('science/measure',C.drawingRequest()),target=$('porosity-results');target.replaceChildren();C.metric(target,'Drawing-plane opening / Å',number(result.diameter_angstrom,2));C.metric(target,'Drawing-plane opening / nm',number(result.diameter_angstrom/10,2));C.note(target,result.scope);C.note(target,'This estimates one opening in the drawing plane. PLD, LCD, density and periodic pore volume require an explicit periodic structure.');downloads(target,result);C.details(target,result);C.status('Measured the current drawing’s opening without a CIF.');}));root.append(b);}
    }
    $('screen-current').disabled=!available;$('calculate-pxrd').disabled=!available;
    $('compare-pxrd-stacks').disabled=!available||(!periodic&&C.origin!=='drawing');
    $('screen-pores').disabled=!periodic;
    $('calculate-pores').disabled=!periodic;
    $('calculate-pxrd').textContent=periodic?'Calculate periodic PXRD':'Calculate finite scattering';
    $('xrd-fwhm').disabled=!periodic;
  }
  function table(parent,headers,rows){
    const wrap=text('div','','materials-table'),grid=document.createElement('table'),head=document.createElement('tr');
    for(const label of headers)head.append(text('th',label));grid.append(head);
    for(const row of rows){const tr=document.createElement('tr');for(const value of row)tr.append(text('td',String(value??'—')));grid.append(tr);}
    wrap.append(grid);parent.append(wrap);
  }
  function downloads(parent,result){
    const row=text('div','','materials-downloads');
    for(const [format,path] of Object.entries(result.downloads||{})){
      const button=text('button',`Save ${format.toUpperCase()}`,'button outline');
      button.addEventListener('click',()=>S.saveExport(path,format));row.append(button);
    }
    if(result.report_id)row.append(text('small',`Report ${result.report_id}`));parent.append(row);
  }
  function warnings(parent,items){for(const warning of items||[])C.note(parent,warning,true);}
  function showInspection(parent,entry){
    const summary=entry.summary,network=entry.network;
    C.metric(parent,'Cell formula',summary.formula);
    C.metric(parent,'Density / g cm⁻³',number(summary.density_g_cm3));
    C.metric(parent,'Cell volume / Å³',number(summary.cell_volume_angstrom3,1));
    C.metric(parent,'Atoms / cell',summary.atom_count);
    C.metric(parent,'Distance-inferred network',`${network.maximum_dimension}D`);
    C.metric(parent,'Components in cell graph',network.quotient_component_count);
    const cell=summary.cell_parameters;
    if(cell)C.note(parent,`a = ${number(cell.a)} Å · b = ${number(cell.b)} Å · c = ${number(cell.c)} Å; α = ${number(cell.alpha,2)}° · β = ${number(cell.beta,2)}° · γ = ${number(cell.gamma,2)}°`);
    if(summary.elements)table(parent,['Element','Count','Mass %','Element sites / mmol g⁻¹'],summary.elements.map(e=>[e.element,e.count,number(e.mass_fraction*100,2),number(e.site_density_mmol_g)]));
    if(entry.porosity){
      for(const [key,label] of Object.entries({largest_sampled_included_sphere_diameter_angstrom:'Sampled included sphere / Å',geometric_void_volume_cm3_g:'Geometric void volume / cm³ g⁻¹',local_probe_exposed_surface_m2_g:'Local probe-exposed surface / m² g⁻¹'}))C.metric(parent,label,number(entry.porosity[key]));
      C.note(parent,'Geometric screening of the supplied model; the area is not BET.');
    }
    warnings(parent,summary.warnings);warnings(parent,network.warnings);
    C.details(parent,{summary,network,porosity:entry.porosity});
  }
  function showComparison(result){
    const single=$('inspection-results'),comparison=$('comparison-results');single.replaceChildren();comparison.replaceChildren();
    if(result.entries.length===1&&result.entries[0].status==='ok')showInspection(single,result.entries[0]);
    else table(comparison,['Model','Formula','Density / g cm⁻³','Network','Void / cm³ g⁻¹','Status'],result.entries.map(e=>[
      e.name,e.summary?.formula||'—',number(e.summary?.density_g_cm3),e.network?`${e.network.maximum_dimension}D`:'—',number(e.porosity?.geometric_void_volume_cm3_g),e.status==='ok'?'Checked':e.error]));
    downloads(comparison,result);C.details(comparison,result);
    C.status(`Structure check complete: ${result.successful} of ${result.entries.length} models processed. Review inferred bonds and quality diagnostics.`);
  }
  $('screen-current').addEventListener('click',()=>C.work(async api=>{
    if(!C.structure){
      const result=await api('science/finite-inspect',{model:finite()}),root=$('inspection-results');root.replaceChildren();$('comparison-results').replaceChildren();
      C.metric(root,'Displayed atom inventory',result.displayed_atom_inventory);C.metric(root,'Real atoms',result.atom_count);
      C.metric(root,'Continuation sites',result.continuation_site_count);C.metric(root,'Supplied real-atom bonds',result.known_real_bond_count);
      C.metric(root,'Finite graph components',result.known_graph_components);C.metric(root,'Extent x / y / z · Å',result.extent_angstrom.map(x=>number(x,1)).join(' / '));
      if(result.contact_checks.available){C.metric(root,'Short contacts',result.contact_checks.short_contact_count);C.metric(root,'Closest atoms / Å',number(result.contact_checks.minimum_distance_angstrom));}
      C.note(root,result.scope);warnings(root,result.warnings);downloads(root,result);C.details(root,result);
      C.status('Checked the displayed finite structure. No CIF was required.');return;
    }
    const model=structure(),result=await api('science/screen',{entries:[{name:$('viewer-title').textContent,structure:model}],include_porosity:$('screen-pores').checked});showComparison(result);
    if(result.entries[0]?.network&&C.viewer.atoms.every(a=>a.network===undefined)){const groups=new Map();result.entries[0].network.components.forEach((group,index)=>group.atom_indices.forEach(i=>groups.set(i,index)));for(const atom of C.viewer.atoms)atom.inferred_component=groups.get(atom.index%model.atoms.length);C.viewer.draw();C.note($('inspection-results'),'Network coloring can show distance-inferred quotient-graph components; this does not establish an interpenetration degree.');}
  }));
  $('screen-files').addEventListener('click',()=>$('screen-cif-files').click());
  $('screen-cif-files').addEventListener('change',e=>{
    const files=[...e.target.files];e.target.value='';if(!files.length)return;
    C.work(async api=>{
      if(files.length>8)throw new Error('Compare at most eight CIF files at once.');
      const entries=[],failures=[];
      for(const file of files){
        try{if(file.size>2_000_000)throw new Error('CIF limit is 2 MB.');const result=await api('science/model',{cif_text:await file.text(),repeats:[1,1,1]});entries.push({name:file.name,structure:result.structure});}
        catch(error){failures.push(`${file.name}: ${error.message}`);}
      }
      if(!entries.length)throw new Error(failures.join('\n'));
      const result=await api('science/screen',{entries,include_porosity:$('screen-pores').checked});showComparison(result);warnings($('comparison-results'),failures);
    });
  });
  function drawPattern(profile,series=[]){
    const root=$('pxrd-plot');root.replaceChildren();const ns='http://www.w3.org/2000/svg';
    const svg=document.createElementNS(ns,'svg');svg.setAttribute('viewBox','0 0 600 285');svg.setAttribute('role','img');svg.setAttribute('aria-label','Calculated model powder diffraction pattern');
    const xs=profile.two_theta_degrees,ys=profile.intensity,min=xs[0],max=xs.at(-1),peak=Math.max(...ys,1);
    const point=(x,y)=>[52+(x-min)/(max-min)*525,235-y/peak*205];
    function element(tag,attrs,label){const el=document.createElementNS(ns,tag);for(const [key,value] of Object.entries(attrs))el.setAttribute(key,value);if(label!==undefined)el.textContent=label;svg.append(el);return el;}
    element('path',{d:'M52 20V235H577',fill:'none',stroke:'#9aafa6','stroke-width':1});
    for(let i=0;i<=5;i++){const value=min+(max-min)*i/5,x=point(value,0)[0];element('path',{d:`M${x} 235v5`,stroke:'#9aafa6'});element('text',{x,y:255,'text-anchor':'middle',fill:'#657a70','font-size':11},value.toFixed(1));}
    element('text',{x:315,y:278,'text-anchor':'middle',fill:'#375e52','font-size':12},'2θ / degrees');
    element('text',{x:15,y:130,transform:'rotate(-90 15 130)','text-anchor':'middle',fill:'#375e52','font-size':12},'Relative intensity');
    const colors=['#217d70','#c66a46','#7165b5','#2d80b5','#ad5382'];
    const lines=series.length?series:[{profile}];
    lines.forEach((row,index)=>element('path',{d:row.profile.two_theta_degrees.map((x,i)=>(i?'L':'M')+point(x,row.profile.intensity[i]).join(',')).join(' '),fill:'none',stroke:colors[index%colors.length],'stroke-width':1.5,'stroke-linejoin':'round'}));
    root.append(svg);
    if(series.length){const legend=text('div','','diffraction-legend');series.forEach((row,index)=>{const label=text('span',row.label);label.style.borderColor=colors[index%colors.length];legend.append(label);});root.append(legend);}
  }
  function diffractionParameters(){return {wavelength:Number($('xrd-wavelength').value),two_theta_min:Number($('xrd-min').value),two_theta_max:Number($('xrd-max').value)};}
  async function pattern(api,model,periodic){return api(periodic?'science/pxrd':'science/finite-powder',{...diffractionParameters(),...(periodic?{structure:model,fwhm:Number($('xrd-fwhm').value),step:Math.min(.02,Number($('xrd-fwhm').value)/2)}:{model,step:.1})});}
  $('calculate-pxrd').addEventListener('click',()=>C.work(async api=>{
    const periodic=Boolean(C.structure),result=await pattern(api,periodic?structure():finite(),periodic);
    drawPattern(result.profile);const root=$('pxrd-results');root.replaceChildren();
    if(result.peaks)table(root,['2θ / °','d / Å','Relative I','h k l'],[...result.peaks].sort((a,b)=>b.intensity-a.intensity).slice(0,12).map(p=>[number(p.two_theta_degrees),number(p.d_spacing_angstrom),number(p.intensity,1),p.hkls.map(h=>Array.isArray(h)?h.join(' '):JSON.stringify(h)).join('; ')]));
    else C.note(root,result.scope);
    warnings(root,result.warnings);downloads(root,result);C.details(root,{method:result.method,parameters:result.parameters});
    C.status(periodic?'Model PXRD calculated from the periodic cell. CSV and JSON are ready.':'Finite-stack scattering calculated from the displayed atoms. CSV and JSON are ready; no periodic crystal was inferred.');
  }));
  $('compare-pxrd-stacks').addEventListener('click',()=>C.work(async(api,check)=>{
    const periodic=Boolean(C.structure),parent=periodic?structure():C.patchRequest(),modes=periodic?['AA','AB','ABC','slipped']:['AA','AB','ABC','slipped','twisted'],rows=[];
    for(const mode of modes){
      const model=await api(periodic?'science/stacking':'science/patch',periodic?{structure:parent,stacking:mode,spacing:Number($('patch-spacing').value),slip:[Number($('patch-slip-x').value),Number($('patch-slip-y').value)]}:{...parent,stacking:mode});
      const result=await pattern(api,periodic?model.structure:model.viewer,periodic);rows.push({label:mode,profile:result.profile,result,model});
    }
    check();drawPattern(rows[0].profile,rows);const root=$('pxrd-results');root.replaceChildren();
    C.note(root,periodic?'Same supplied planar layer, explicit periodic AA / AB / ABC / progressive-slip cells. Geometry is unrelaxed. Arbitrary twist needs a commensurate cell.':'Same finite drawing and layer count under five chosen transforms. These are Debye scattering curves, without reciprocal-lattice peaks or hkl assignments.');
    C.note(root,'Each curve is normalized to its own maximum. This compares shape and positions, not absolute scattering power or phase fractions.');
    for(const row of rows){const box=text('div','','metric-note');box.append(text('strong',row.label));downloads(box,row.result);const view=text('button','View this arrangement','button outline');view.addEventListener('click',()=>{C.invalidate('',false);C.display(row.model,`${row.label} · ${periodic?'periodic stacking hypothesis':'finite drawing stack'}`,periodic?'periodic':'drawing');$('patch-stacking').value=row.label;C.status(`${row.label} comparison selected. Recalculate to inspect this arrangement.`);});box.append(view);root.append(box);}
    warnings(root,rows[0].result.warnings);C.status(`${rows.length} stacking patterns compared under identical calculation settings.`);
  }));
  $('viewer-axis').addEventListener('change',()=>{
    const axis=$('viewer-axis').value;if(axis==='free')return;const index={a:0,b:1,c:2}[axis],direction=C.viewer.model?.cell?.[index];
    if(!direction)return C.status('Axis views require a periodic cell.',true);
    const norm=Math.hypot(...direction),n=direction.map(v=>v/norm);C.viewer.phi=Math.asin(n[1]);C.viewer.theta=Math.atan2(-n[0],n[2]);C.viewer.draw();
  });
  $('viewer-hide-h').addEventListener('change',()=>{C.viewer.hideHydrogen=$('viewer-hide-h').checked;C.viewer.selected=[];$('viewer-distance').textContent='';C.viewer.draw();});
  $('show-atom-table').addEventListener('click',()=>{
    const root=$('atom-table-panel');root.hidden=!root.hidden;if(root.hidden)return;root.replaceChildren();
    if(C.structure)table(root,['Site','Element','Fractional x','y','z'],C.structure.atoms.map((a,i)=>[i+1,a.element,...a.fractional.map(x=>number(x,5))]));
    else table(root,['Atom','Element','x / Å','y / Å','z / Å'],C.viewer.atoms.map((a,i)=>[i+1,a.element,number(a.x,5),number(a.y,5),number(a.z,5)]));
  });
  function clearAnalysis(){for(const id of ['inspection-results','comparison-results','pxrd-results','pxrd-plot','atom-table-panel'])$(id).replaceChildren();$('atom-table-panel').hidden=true;$('viewer-axis').value='free';}
  window.addEventListener('cofdraw:structure',()=>{clearAnalysis();syncModelState();});
  for(const id of ['xrd-wavelength','xrd-fwhm','xrd-min','xrd-max'])$(id).addEventListener('input',()=>{C.invalidate('',false);$('pxrd-results').replaceChildren();$('pxrd-plot').replaceChildren();});
  $('screen-pores').addEventListener('change',()=>{C.invalidate('',false);$('inspection-results').replaceChildren();$('comparison-results').replaceChildren();});
  syncModelState();
})();
