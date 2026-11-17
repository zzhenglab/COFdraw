/* Structure inspection and evidence retrieval, using the same local API as agents. */
"use strict";
(() => {
  const S=window.COFdrawStudio, C=window.COFdrawScience;
  const $=id=>document.getElementById(id);
  const number=(value,digits=3)=>Number.isFinite(value)?Number(value).toLocaleString(undefined,{maximumFractionDigits:digits}):'—';
  const text=(tag,value,className)=>{const el=document.createElement(tag);el.textContent=value;if(className)el.className=className;return el;};
  function structure(){if(!C.structure)throw new Error('Generate a periodic model or import a CIF first.');return C.structure;}
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
  $('screen-current').addEventListener('click',()=>C.work(async()=>{
    const result=await S.api('science/screen',{entries:[{name:$('viewer-title').textContent,structure:structure()}],include_porosity:$('screen-pores').checked});showComparison(result);
  }));
  $('screen-files').addEventListener('click',()=>$('screen-cif-files').click());
  $('screen-cif-files').addEventListener('change',e=>{
    const files=[...e.target.files];e.target.value='';if(!files.length)return;
    C.work(async()=>{
      if(files.length>8)throw new Error('Compare at most eight CIF files at once.');
      const entries=[],failures=[];
      for(const file of files){
        try{if(file.size>2_000_000)throw new Error('CIF limit is 2 MB.');const result=await S.api('science/model',{cif_text:await file.text(),repeats:[1,1,1]});entries.push({name:file.name,structure:result.structure});}
        catch(error){failures.push(`${file.name}: ${error.message}`);}
      }
      if(!entries.length)throw new Error(failures.join('\n'));
      const result=await S.api('science/screen',{entries,include_porosity:$('screen-pores').checked});showComparison(result);warnings($('comparison-results'),failures);
    });
  });
  function drawPattern(profile){
    const root=$('pxrd-plot');root.replaceChildren();const ns='http://www.w3.org/2000/svg';
    const svg=document.createElementNS(ns,'svg');svg.setAttribute('viewBox','0 0 600 285');svg.setAttribute('role','img');svg.setAttribute('aria-label','Calculated model powder diffraction pattern');
    const xs=profile.two_theta_degrees,ys=profile.intensity,min=xs[0],max=xs.at(-1),peak=Math.max(...ys,1);
    const point=(x,y)=>[52+(x-min)/(max-min)*525,235-y/peak*205];
    function element(tag,attrs,label){const el=document.createElementNS(ns,tag);for(const [key,value] of Object.entries(attrs))el.setAttribute(key,value);if(label!==undefined)el.textContent=label;svg.append(el);return el;}
    element('path',{d:'M52 20V235H577',fill:'none',stroke:'#9aafa6','stroke-width':1});
    for(let i=0;i<=5;i++){const value=min+(max-min)*i/5,x=point(value,0)[0];element('path',{d:`M${x} 235v5`,stroke:'#9aafa6'});element('text',{x,y:255,'text-anchor':'middle',fill:'#657a70','font-size':11},value.toFixed(1));}
    element('text',{x:315,y:278,'text-anchor':'middle',fill:'#375e52','font-size':12},'2θ / degrees');
    element('text',{x:15,y:130,transform:'rotate(-90 15 130)','text-anchor':'middle',fill:'#375e52','font-size':12},'Relative intensity');
    element('path',{d:xs.map((x,i)=>(i?'L':'M')+point(x,ys[i]).join(',')).join(' '),fill:'none',stroke:'#217d70','stroke-width':1.7,'stroke-linejoin':'round'});
    root.append(svg);
  }
  $('calculate-pxrd').addEventListener('click',()=>C.work(async()=>{
    const result=await S.api('science/pxrd',{structure:structure(),wavelength:Number($('xrd-wavelength').value),two_theta_min:Number($('xrd-min').value),two_theta_max:Number($('xrd-max').value),fwhm:Number($('xrd-fwhm').value),step:Math.min(.02,Number($('xrd-fwhm').value)/2)});
    drawPattern(result.profile);const root=$('pxrd-results');root.replaceChildren();
    table(root,['2θ / °','d / Å','Relative I','h k l'],[...result.peaks].sort((a,b)=>b.intensity-a.intensity).slice(0,12).map(p=>[number(p.two_theta_degrees),number(p.d_spacing_angstrom),number(p.intensity,1),p.hkls.map(h=>Array.isArray(h)?h.join(' '):JSON.stringify(h)).join('; ')]));
    warnings(root,result.warnings);downloads(root,result);C.details(root,{method:result.method,parameters:result.parameters});
    C.status('Model PXRD calculated. Export the numeric CSV for comparison or publication plotting.');
  }));
  async function query(){
    const type=$('knowledge-type').value,result=await S.api('science/query',{query:$('knowledge-query').value,record_types:type?[type]:null,limit:20});
    const root=$('knowledge-results');root.replaceChildren(text('p',`${result.matched_count} matching records`,'small-note'));
    const save=text('button','Save search results','button outline');save.addEventListener('click',()=>S.saveText('cofdraw-catalog-query.json',JSON.stringify(result,null,2)));root.append(save);
    for(const record of result.records){
      const card=text('article','','knowledge-card');card.append(text('strong',record.label),text('small',`${record.record_type} · ${record.id}`));
      if(record.support?.status)card.append(text('p',`Availability: ${record.support.status.replaceAll('_',' ')}`));
      const qualifiers=Array.isArray(record.qualifications)?record.qualifications:[record.qualifications];for(const qualifier of qualifiers.filter(Boolean))card.append(text('p',typeof qualifier==='string'?qualifier:JSON.stringify(qualifier)));
      if(record.record_type==='source'){
        const b=text('button','Open source drawing','button outline');b.addEventListener('click',async()=>{try{const r=await S.api(`sources/${encodeURIComponent(record.record_id)}`);S.acceptRecipe(r.recipe,record.record_id,r.metadata);S.showSourceMonomers(r.recipe,r.metadata);S.state.style=r.recipe.source_style?'source':'standard';S.state.styleOptions={};$('science-dialog').close();await S.render();}catch(error){C.status(error.message,true);}});card.append(b);
      }else if(record.record_type==='crystal'){
        const b=text('button','Open periodic model','button outline');b.addEventListener('click',()=>{$('crystal-template').value=record.record_id;$('crystal-template').dispatchEvent(new Event('change'));C.tab('model');$('build-crystal').click();});card.append(b);
      }
      C.details(card,record);root.append(card);
    }
    for(const error of result.errors||[])C.note(root,error,true);
  }
  $('query-knowledge').addEventListener('click',()=>C.work(query));
  $('knowledge-query').addEventListener('keydown',e=>{if(e.key==='Enter'){e.preventDefault();C.work(query);}});
  $('query-reports').addEventListener('click',()=>C.work(async()=>{
    const result=await S.api('science/reports'),root=$('saved-reports');root.replaceChildren();
    if(!result.reports.length)C.note(root,'No scientific reports have been saved in this data directory yet.');
    for(const report of result.reports){const b=text('button',`${report.operation} · ${new Date(report.created_at_utc).toLocaleString()}`,'button outline');b.addEventListener('click',async()=>{try{const r=await S.api(`science/reports/${report.report_id}`);C.details(root,r);downloads(root,r);}catch(error){C.status(error.message,true);}});root.append(b);}
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
  window.addEventListener('cofdraw:structure',clearAnalysis);
  for(const id of ['xrd-wavelength','xrd-fwhm','xrd-min','xrd-max'])$(id).addEventListener('input',()=>{$('pxrd-results').replaceChildren();$('pxrd-plot').replaceChildren();});
  $('screen-pores').addEventListener('change',()=>{$('inspection-results').replaceChildren();$('comparison-results').replaceChildren();});
})();
