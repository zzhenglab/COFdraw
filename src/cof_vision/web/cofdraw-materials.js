/* Structure inspection and evidence retrieval, using the same local API as agents. */
"use strict";
(() => {
  const S=window.COFdrawStudio, C=window.COFdrawScience;
  const $=id=>document.getElementById(id);
  const number=(value,digits=3)=>Number.isFinite(value)?Number(value).toLocaleString(undefined,{maximumFractionDigits:digits}):'—';
  const text=(tag,value,className)=>{const el=document.createElement(tag);el.textContent=value;if(className)el.className=className;return el;};
  function structure(){if(!C.structure)throw new Error('Generate a periodic model or import a CIF first.');return C.structure;}
  function finite(){if(!C.finiteModel)throw new Error('Generate or select the structure to inspect first.');return C.finiteModel;}
  let customDiffractionRange=false;
  for(const id of ['xrd-min','xrd-max'])$(id).addEventListener('input',()=>{customDiffractionRange=true;});
  function syncModelState(){
    const periodic=Boolean(C.structure),available=periodic||Boolean(C.finiteModel);
    if(!customDiffractionRange&&available){
      // The complete curated collection has been screened at this angular range.
      $('xrd-max').value=C.structure?.provenance?.dataset_id?'20':'50';
    }
    $('viewer-preset').querySelector('[value=projection]').disabled=!periodic;
    for(const id of ['inspect-model-state','pxrd-model-state','porosity-model-state']){
      const root=$(id);root.replaceChildren();root.append(text('span',periodic?'Current view: periodic atomic model with a unit cell.':available?'Current view: finite atoms from the drawing. Checks and finite-cluster scattering are available without a CIF. Density, PLD and periodic Bragg peaks require a periodic cell.':'Generate your drawing in 3D, or select a periodic model.'));
      if(!periodic){const b=text('button','Choose a periodic model','button outline');b.addEventListener('click',()=>{C.tab('model');$('structure-source').value='periodic';$('structure-source').dispatchEvent(new Event('change'));});root.append(b);}
      if(!periodic&&C.drawingReady()){
        const generate=text('button','Build CIF from current drawing','button outline');
        generate.classList.add('build-cif-from-drawing');
        generate.addEventListener('click',C.openCurrentCell);
        root.append(generate);
      }
      if(id==='porosity-model-state'&&!periodic&&available&&C.origin==='drawing'){const b=text('button','Measure this drawing’s pore','button outline');b.id='measure-current-finite-pore';b.addEventListener('click',()=>C.work(async api=>{const result=await api('science/measure',C.drawingRequest()),target=$('porosity-results');target.replaceChildren();C.metric(target,'Drawing-plane opening / Å',number(result.diameter_angstrom,2));C.metric(target,'Drawing-plane opening / nm',number(result.diameter_angstrom/10,2));C.note(target,result.scope);C.note(target,'This estimates one opening in the drawing plane. PLD, LCD, density and periodic pore volume require an explicit periodic structure.');downloads(target,result);C.details(target,result);C.status('Measured the current drawing’s opening without a CIF.');}));root.append(b);}
    }
    $('screen-current').disabled=!available;$('calculate-pxrd').disabled=!available;
    $('viewer-axis').disabled=!periodic;$('check-peaks').disabled=!periodic;
    $('check-peaks').title=periodic?'':'Build a periodic cell to check measured reflections.';
    const provenance=C.structure?.provenance;
    const nonLayer=provenance?.network_dimension===1||provenance?.network_dimension===3||provenance?.template==='dia-imine';
    $('compare-pxrd-stacks').disabled=!available||(!periodic&&C.origin!=='drawing')||nonLayer;
    $('compare-pxrd-stacks').title=nonLayer?'Layer stacking requires a two-dimensional sheet.':'';
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
  let patternVisible=false;
  let positionTraces=null;
  const separation=$('pattern-separation');
  separation.addEventListener('input',()=>positionTraces?.(Number(separation.value)));
  function updatePatternView(){
    const hasPattern=Boolean($('pxrd-plot').querySelector('svg'));
    const active=document.querySelector('[data-science-tab="pxrd"]').getAttribute('aria-selected')==='true';
    const show=active&&hasPattern&&patternVisible;
    $('pattern-view').hidden=!show;$('model-view').hidden=show;
    $('return-to-pattern').hidden=!active||!hasPattern||show;
    if(!show)requestAnimationFrame(()=>C.viewer.draw());
  }
  $('return-to-model').addEventListener('click',()=>{patternVisible=false;updatePatternView();$('return-to-pattern').focus();});
  $('return-to-pattern').addEventListener('click',()=>{patternVisible=true;updatePatternView();$('return-to-model').focus();});
  window.addEventListener('cofdraw:science-tab',updatePatternView);
  $('save-pattern-svg').addEventListener('click',async()=>{
    const svg=$('pxrd-plot').querySelector('svg');if(!svg)return;
    try{await S.saveText('COFdraw-powder-pattern.svg',new XMLSerializer().serializeToString(svg),'image/svg+xml');}
    catch(error){C.status(error.message,true);}
  });
  function drawPattern(profile,series=[],observations=[],report={},detail=false,selected=0,offset=0){
    positionTraces=null;
    const comparison=series.length>1;
    $('pattern-spacing').hidden=!comparison;
    separation.value=String(offset);
    const root=$('pxrd-plot');root.replaceChildren();const ns='http://www.w3.org/2000/svg';
    const svg=document.createElementNS(ns,'svg');svg.setAttribute('role','img');
    const xs=profile.two_theta_degrees,min=xs[0],end=xs.at(-1);
    const limit=$('pattern-limit');
    limit.min=String(Number((min+.1).toFixed(2)));limit.max=String(end);
    if(!Number.isFinite(limit.valueAsNumber)||limit.valueAsNumber<=min)limit.value=String(Math.min(end,min+8));
    const max=detail?Math.min(end,limit.valueAsNumber):end;
    const redraw=(zoom,which=selected)=>drawPattern(profile,series,observations,report,zoom,which,Number(separation.value));
    $('zoom-pattern').textContent=detail?'Full range':'Low-angle detail';
    $('zoom-pattern').setAttribute('aria-pressed',String(detail));
    $('zoom-pattern').onclick=()=>redraw(!detail);
    limit.onchange=()=>{if(limit.reportValidity())redraw(true);};
    const point=(x,y)=>[(52+(x-min)/(max-min)*525).toFixed(2),(235-y/100*205).toFixed(2)];
    function element(tag,attrs,label){const el=document.createElementNS(ns,tag);for(const [key,value] of Object.entries(attrs))el.setAttribute(key,value);if(label!==undefined)el.textContent=label;svg.append(el);return el;}
    const title=C.structure?'Simulated powder diffraction':'Finite-cluster scattering';
    $('pattern-title').textContent=title;
    svg.setAttribute('aria-label',title);
    svg.setAttribute('viewBox','0 0 600 322');
    svg.setAttribute('font-family','Arial, sans-serif');
    element('title',{},title);
    const provenance=$('viewer-provenance').textContent+' · '+$('viewer-title').textContent;
    $('pattern-source').textContent=provenance;
    const metadata={source:provenance,operation:report.operation,parameters:report.parameters,
      model_sha256:report.structure_sha256,report_id:report.report_id,observed:observations.map(r=>r.observed_two_theta),
      series:series.map(r=>({label:r.label,report_id:r.result.report_id,structure_sha256:r.result.structure_sha256})),normalization:'Each curve normalized to 100'};
    metadata.display_two_theta_degrees=[min,max];
    const metadataNode=element('metadata',{},JSON.stringify(metadata));
    element('text',{x:52,y:14,fill:'#254b43','font-size':12},title+' · λ '+$('xrd-wavelength').value+' Å');
    element('path',{d:'M52 20V235',fill:'none',stroke:'#9aafa6','stroke-width':1});
    for(const intensity of [0,50,100]){
      const y=Number(point(min,intensity)[1]);
      element('line',{x1:48,x2:52,y1:y,y2:y,stroke:'#9aafa6'});
      element('text',{x:44,y:y+4,'text-anchor':'end',fill:'#657a70','font-size':10},String(intensity));
    }
    const sharedAxis=element('g',{});
    sharedAxis.append(element('path',{d:'M52 235H577',stroke:'#9aafa6'}));
    for(let i=0;i<=5;i++){const value=min+(max-min)*i/5,x=point(value,0)[0];sharedAxis.append(element('path',{d:`M${x} 235v5`,stroke:'#9aafa6'}),element('text',{x,y:255,'text-anchor':'middle',fill:'#657a70','font-size':11},value.toFixed(1)));}
    const axisTitle=element('text',{x:315,y:278,'text-anchor':'middle',fill:'#375e52','font-size':12},'2θ / degrees');
    element('text',{x:15,y:130,transform:'rotate(-90 15 130)','text-anchor':'middle',fill:'#375e52','font-size':12},'Relative intensity');
    const colors=['#217d70','#c66a46','#7165b5','#2d80b5','#ad5382'];
    const lines=series.length?series:[{profile}];
    const traces=lines.map((row,index)=>element('path',{'data-trace':row.label||'current',d:row.profile.two_theta_degrees.map((x,i)=>x<=max?(i?'L':'M')+point(x,row.profile.intensity[i]).join(','):'').join(' ').trim(),fill:'none',stroke:colors[index%colors.length],'stroke-width':1.5,'stroke-linejoin':'round','vector-effect':'non-scaling-stroke'}));
    for(const observation of observations){
      if(observation.status==='outside range'||observation.observed_two_theta>max)continue;
      const x=point(observation.observed_two_theta,0)[0];
      const tick=element('line',{x1:x,x2:x,y1:223,y2:235,stroke:observation.status==='unmatched'?'#ae4833':'#4f548b','stroke-width':2});
      const title=document.createElementNS(ns,'title');title.textContent=observation.observed_two_theta+'° · '+observation.status;tick.append(title);
    }
    const legend=element('g',{});
    if(series.length)series.forEach((row,i)=>{
      const x=52+(i+0.24)*525/series.length;legend.append(element('line',{x1:x,x2:x+16,y1:305,y2:305,stroke:colors[i%colors.length],'stroke-width':2}),
        element('text',{'data-series-label':row.label,x:x+20,y:309,'font-size':10,fill:'#375e52'},row.label));
    });
    if(observations.length)element('text',{x:52,y:309,'font-size':11,fill:'#375e52'},'Measured peak ticks · violet: candidate center within tolerance · red: unmatched');
    const peakBox=$('pattern-peak-list');peakBox.hidden=!detail;
    const markers=[];
    if(detail){
      const row=lines[selected],src=row.profile,data=row.result||report;
      const entries=(data.peaks||src.two_theta_degrees.flatMap((x,i)=>{
        const y=src.intensity[i];
        return i>0&&i+1<src.intensity.length&&y>src.intensity[i-1]&&y>=src.intensity[i+1]?[{two_theta_degrees:x,intensity:y}]:[];
      })).filter(p=>p.two_theta_degrees>=min&&p.two_theta_degrees<=max).sort((a,b)=>a.two_theta_degrees-b.two_theta_degrees);
      const body=peakBox.querySelector('div');body.replaceChildren();
      peakBox.open=true;
      peakBox.querySelector('summary').textContent=(row.label?row.label+' · ':'')+entries.length+(data.peaks?' reflections':' profile maxima')+' in view';
      if(comparison){
        const select=document.createElement('select');select.setAttribute('aria-label','Inspect peaks for stacking');
        lines.forEach((v,i)=>select.add(new Option(v.label,String(i))));select.value=String(selected);
        select.onchange=()=>redraw(true,Number(select.value));body.append(select);
      }
      table(body,['Peak','2θ / °','Relative I',data.peaks?'h k l':'Type'],entries.map((p,i)=>[i+1,p.two_theta_degrees.toFixed(3),number(p.intensity,1),data.peaks?p.hkls.map(h=>h.join(' ')).join('; '):'Profile maximum']));
      const labels=[];
      for(const {p,num} of entries.map((p,i)=>({p,num:i+1})).sort((a,b)=>b.p.intensity-a.p.intensity)){
        const x=Number(point(p.two_theta_degrees,0)[0]);
        if(labels.some(v=>Math.abs(v-x)<32)||labels.length>=8)continue;
        labels.push(x);
        const idx=Math.max(0,Math.min(src.intensity.length-1,Math.round((p.two_theta_degrees-xs[0])/(xs[1]-xs[0]))));
        const intensity=src.intensity[idx];
        const mark=element('text',{x,y:Math.max(26,230-2.05*intensity),'text-anchor':'middle',fill:colors[selected],'font-size':10},String(num));
        markers.push({mark,intensity});
      }
      metadata.peak_series=row.label||'current';metadataNode.textContent=JSON.stringify(metadata);
    }
    if(comparison){
      const ticks=svg.querySelectorAll('line[x1="48"], text[x="44"]');
      positionTraces=value=>{
        const n=traces.length, height=205, gap=6;
        const amplitude=1-value*(1-(height-gap*(n-1))/(n*height));
        const stride=(height-amplitude*height)/(n-1);
        const label=value===0?'0 · Overlay':value===1?'1 · Stacked':value.toFixed(2);
        $('pattern-separation-value').textContent=label;
        separation.setAttribute('aria-valuetext',label);
        metadata.vertical_separation=value;
        metadataNode.textContent=JSON.stringify(metadata);
        ticks.forEach(tick=>tick.setAttribute('visibility',value?'hidden':'visible'));
        traces.forEach((trace,i)=>{
          const baseline=235-(n-1-i)*stride;
          trace.setAttribute('transform',`translate(0 ${baseline}) scale(1 ${amplitude}) translate(0 -235)`);
        });
        markers.forEach(({mark,intensity})=>mark.setAttribute('y',Math.max(25,230-(n-1-selected)*stride-intensity*2.05*amplitude)));
      };
      positionTraces(offset);
    }
    root.append(svg);
    patternVisible=true;updatePatternView();
  }
  function diffractionParameters(){return {wavelength:Number($('xrd-wavelength').value),two_theta_min:Number($('xrd-min').value),two_theta_max:Number($('xrd-max').value)};}
  async function pattern(api,model,periodic){return api(periodic?'science/pxrd':'science/finite-powder',{...diffractionParameters(),...(periodic?{structure:model,fwhm:Number($('xrd-fwhm').value),step:Math.min(.02,Number($('xrd-fwhm').value)/2)}:{model,step:.1})});}
  $('calculate-pxrd').addEventListener('click',()=>C.work(async api=>{
    const periodic=Boolean(C.structure),result=await pattern(api,periodic?structure():finite(),periodic);
    drawPattern(result.profile,[],[],result);const root=$('pxrd-results');root.replaceChildren();
    if(result.peaks)table(root,['2θ / °','d / Å','Relative I','h k l'],[...result.peaks].sort((a,b)=>b.intensity-a.intensity).slice(0,12).map(p=>[number(p.two_theta_degrees),number(p.d_spacing_angstrom),number(p.intensity,1),p.hkls.map(h=>Array.isArray(h)?h.join(' '):JSON.stringify(h)).join('; ')]));
    else C.note(root,result.scope);
    warnings(root,result.warnings);downloads(root,result);C.details(root,{method:result.method,parameters:result.parameters});
    C.status(periodic?'Model PXRD calculated from the periodic cell. CSV and JSON are ready.':'Finite-stack scattering calculated from the displayed atoms. CSV and JSON are ready; no periodic crystal was inferred.');
  }));
  $('check-peaks').addEventListener('click',()=>C.work(async api=>{
    const input=$('observed-peaks').value.trim();
    const observed=input?input.split(/[\s,;]+/).map(Number):[];
    if(!observed.length||observed.length>128||observed.some(v=>!Number.isFinite(v)||v<=0||v>=180))throw new Error('Enter 1–128 peak centers in degrees, separated by spaces or commas.');
    const fwhm=Number($('xrd-fwhm').value);
    const result=await api('science/pxrd-peaks',{structure:structure(),...diffractionParameters(),fwhm,step:Math.min(.02,fwhm/2),
      observed_two_theta:observed,tolerance_degrees:Number($('peak-tolerance').value),min_relative_intensity:Number($('peak-floor').value)});
    drawPattern(result.profile,[],result.observed,result);
    const root=$('pxrd-results');root.replaceChildren();
    C.note(root,'Ticks mark measured peak centers. Red ticks have no calculated center within tolerance.');
    table(root,['Observed 2θ / °','Candidate centers / °','Result'],result.observed.map(row=>[number(row.observed_two_theta),
      row.candidate_indices.map(i=>number(result.peaks[i].two_theta_degrees)).join(', ')||'—',row.status]));
    C.note(root,result.calculated_only_indices.length+' calculated centers above the intensity cutoff have no measured counterpart.');
    if(result.calculated_only_indices.length){
      const details=text('details','');details.append(text('summary','Unmatched calculated centers'));
      const missing=result.calculated_only_indices.map(i=>result.peaks[i]).sort((a,b)=>b.intensity-a.intensity);
      table(details,['2θ / °','Relative I / %','Representative h k l'],missing.slice(0,20).map(p=>[number(p.two_theta_degrees),number(p.intensity,1),p.hkls[0].join(' ')]));
      if(missing.length>20)C.note(details,'Showing the strongest 20. The CSV contains every center.');
      root.append(details);
    }
    C.note(root,result.scope);downloads(root,result);C.status('Peak comparison complete. CSV includes candidate Miller indices.');
  }));
  $('compare-pxrd-stacks').addEventListener('click',()=>C.work(async(api,check)=>{
    const periodic=Boolean(C.structure), origin=C.origin;
    const parent=periodic?C.stackingLayer:C.patchRequest();
    const modes=periodic?['AA','AB','ABC','slipped']:['AA','AB','ABC','slipped','twisted'];
    const spacing=origin==='current-periodic'&&parent.provenance?.spacing_angstrom
      ?parent.provenance.spacing_angstrom:Number($('patch-spacing').value);
    const slip=[Number($('patch-slip-x').value),Number($('patch-slip-y').value)];
    const corrugated_layer=origin==='current-periodic'&&Boolean(parent.provenance?.corrugated_layer || parent.provenance?.tilted_phenyl_ring_count);
    const rows=[], rejected=[];
    for(const mode of modes){
      try{
        const input=periodic?{structure:parent,stacking:mode,spacing,slip,corrugated_layer}:{...parent,stacking:mode};
        const model=await api(periodic?'science/stacking':'science/patch',input);
        const result=await pattern(api,periodic?model.structure:model.viewer,periodic);
        rows.push({label:mode,profile:result.profile,result,model});
      }catch(error){
        if(error.name==='ObsoleteScienceResult')throw error;
        rejected.push(`${mode} · ${error.message}`);
      }
    }
    if(!rows.length)throw new Error(rejected[0]||'No stacking pattern could be calculated.');
    check();drawPattern(rows[0].profile,rows,[],rows[0].result);const root=$('pxrd-results');root.replaceChildren();
    C.note(root,periodic?'Same supplied planar layer, explicit periodic AA / AB / ABC / progressive-slip cells. Geometry is unrelaxed. Arbitrary twist needs a commensurate cell.':'Same finite drawing and layer count under five chosen transforms. These are Debye scattering curves, without reciprocal-lattice peaks or hkl assignments.');
    C.note(root,'Each curve is normalized to its own maximum. This compares shape and positions, not absolute scattering power or phase fractions.');
    for(const row of rows){const box=text('div','','metric-note');box.append(text('strong',row.label));downloads(box,row.result);const view=text('button','View this arrangement','button outline');view.addEventListener('click',()=>{C.invalidate('',false);C.display(row.model,`${row.label} · ${periodic?'periodic stacking hypothesis':'finite drawing stack'}`,origin==='current-periodic'?origin:periodic?'periodic':'drawing',periodic?parent:null);if(!periodic)$('patch-stacking').value=row.label;C.status(`${row.label} comparison selected. Recalculate to inspect this arrangement.`);});box.append(view);root.append(box);}
    warnings(root,rejected);
    warnings(root,rows[0].result.warnings);
    C.status(`${rows.length} stacking patterns compared under identical calculation settings.`);
  }));
  $('viewer-axis').addEventListener('change',()=>{
    const axis=$('viewer-axis').value;if(axis==='free')return;const index={a:0,b:1,c:2}[axis],direction=C.viewer.model?.cell?.[index];
    if(!direction)return C.status('Axis views require a periodic cell.',true);
    const norm=Math.hypot(...direction),n=direction.map(v=>v/norm);C.viewer.phi=Math.asin(n[1]);C.viewer.theta=Math.atan2(-n[0],n[2]);C.viewer.draw();
  });
  $('viewer-hide-h').addEventListener('change',()=>{C.viewer.hideHydrogen=$('viewer-hide-h').checked;C.viewer.selected=[];$('viewer-distance').textContent='';C.viewer.draw();});
  $('show-atom-table').addEventListener('click',()=>{
    const root=$('atom-table-panel');root.hidden=!root.hidden;if(root.hidden)return;root.replaceChildren();
    if(C.structure)table(root,['Site','Element','Fractional x','y','z'],C.structure.atoms.map((a,i)=>[i+1,(a.isotope||'')+a.element,...a.fractional.map(x=>number(x,5))]));
    else table(root,['Atom','Element','x / Å','y / Å','z / Å'],C.viewer.atoms.map((a,i)=>[i+1,(a.isotope||'')+a.element,number(a.x,5),number(a.y,5),number(a.z,5)]));
  });
  function clearAnalysis(){positionTraces=null;$('pattern-spacing').hidden=true;for(const id of ['inspection-results','comparison-results','pxrd-results','pxrd-plot','atom-table-panel'])$(id).replaceChildren();$('atom-table-panel').hidden=true;$('viewer-axis').value='free';}
  window.addEventListener('cofdraw:structure',()=>{clearAnalysis();patternVisible=false;updatePatternView();syncModelState();});
  for(const id of ['xrd-wavelength','xrd-fwhm','xrd-min','xrd-max','observed-peaks','peak-tolerance','peak-floor'])$(id).addEventListener('input',()=>{C.invalidate('',false);$('pxrd-results').replaceChildren();$('pxrd-plot').replaceChildren();});
  $('screen-pores').addEventListener('change',()=>{C.invalidate('',false);$('inspection-results').replaceChildren();$('comparison-results').replaceChildren();});
  const balanceDialog=$('balance-dialog');
  let balanceRevision=0;
  function balanceAvailable(){return !multi()&&!sourceMode() && Boolean(checkedMonomer('a') && checkedMonomer('b'));}
  function clearBalance(){
    balanceRevision++;$('balance-results').replaceChildren();$('balance-status').textContent='';
    $('calculate-balance').disabled=!balanceAvailable();
  }
  function refreshBalancePair(){
    clearBalance();
    const design=balanceAvailable();
    $('balance-pair').textContent=design?['a','b'].map(slot=>slot.toUpperCase()+' · '+($('label-'+slot).value||'custom monomer')).join('   +   '):'Choose a pair in the manual builder to calculate a feed.';
  }
  function openBalance(){
    refreshBalancePair();
    balanceDialog.showModal();$('balance-mass').focus();
  }
  window.addEventListener('cofdraw:state',()=>{if(balanceDialog.open)refreshBalancePair();});
  $('open-balance').addEventListener('click',openBalance);
  $('balance-form').addEventListener('input',clearBalance);
  balanceDialog.addEventListener('close',clearBalance);
  $('balance-form').addEventListener('submit',async event=>{
    event.preventDefault();if(!balanceAvailable())return;clearBalance();const revision=balanceRevision;
    $('calculate-balance').disabled=true;
    try{
      const result=await S.api('science/reagent-balance',{monomer_a:$('smiles-a').value,monomer_b:$('smiles-b').value,linkage:$('linkage').value,
        mass_a_mg:Number($('balance-mass').value),purity_a_percent:Number($('purity-a').value),purity_b_percent:Number($('purity-b').value),b_site_equivalents:Number($('balance-eq').value)});
      if(revision!==balanceRevision||!balanceDialog.open)return;
      const root=$('balance-results'),[a,b]=result.reagents;
      C.metric(root,'Weigh block B',Number(b.weigh_mg.toPrecision(5))+' mg');
      C.note(root,'Balanced monomer ratio A : B = '+result.balanced_molar_ratio_a_b.join(' : '));
      table(root,['Block','Formula','Sites','M / g mol⁻¹','Weigh / mg','Molecule / mmol','Sites / mmol'],[a,b].map(r=>[r.block,r.formula,r.sites_per_molecule,number(r.molar_mass_g_mol),number(r.weigh_mg,5),number(r.monomer_mmol,6),number(r.site_mmol,6)]));
      C.note(root,result.scope);downloads(root,result);
    }catch(error){if(revision===balanceRevision)$('balance-status').textContent=error.message;}
    finally{if(revision===balanceRevision)$('calculate-balance').disabled=false;}
  });
  $('tool-jump').addEventListener('change',event=>{
    const name=event.target.value;event.target.value='';
    if(name==='balance'){openBalance();return;}
    if(name==='library'||name.startsWith('sketch-')){
      $('mode-manual').click();
      if(name==='library')$('open-monomer-library').click();
      else document.querySelector('.sketch-open[data-block="'+name.slice(-1)+'"]').click();
      return;
    }
    if(!name)return;
    if($('science-dialog').hidden)$('open-science').click();
    C.tab(name);
    (name==='inspect'?$('science-inspect').querySelector('summary'):document.querySelector('[data-science-tab="'+name+'"]')).focus();
  });
  syncModelState();
})();
