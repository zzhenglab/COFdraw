/* Local CPU RASPA jobs. Results belong to the captured model and settings. */
"use strict";
(() => {
  const S=window.COFdrawStudio,C=window.COFdrawScience,$=id=>document.getElementById(id);
  const pane=document.createElement('div');pane.id='science-adsorption';pane.className='science-pane';pane.hidden=true;
  pane.innerHTML=`<h3>Adsorption with RASPA</h3>
    <p class="science-intro">Run a rigid periodic CIF locally. Built-in structures are unrelaxed model hypotheses; inspect your geometry and charges first.</p>
    <p id="ads-engine" class="small-note">Checking the bundled CPU engine…</p>
    <label class="field"><span>Gas or mixture</span><select id="ads-gas"><option>CO2</option><option>N2</option><option>CH4</option><option>H2</option><option>Ar</option><option value="CO2,N2">CO₂ / N₂</option><option value="CO2,CH4">CO₂ / CH₄</option></select></label>
    <label id="ads-fraction-field" class="field" hidden><span>CO₂ gas-phase mole fraction</span><input id="ads-fraction" type="number" value="0.15" min="0.001" max="0.999" step="0.01"></label>
    <div class="input-row"><label class="field"><span>Temperature / K</span><input id="ads-temperature" type="number" value="298.15" min="20" max="600" step="0.1"></label><label class="field"><span>Pressures / bar</span><input id="ads-pressures" value="0.1, 1" aria-label="Pressure points in bar, separated by commas"></label></div>
    <label class="field"><span>Framework partial charges</span><select id="ads-charges"><option value="zero">Explicit uncharged-framework approximation</option><option value="cif">Use complete supplied CIF charges</option></select></label>
    <label class="auto-chemistry"><input id="ads-heat" type="checkbox" checked> Estimate zero-loading heat (three Widom runs)</label>
    <details><summary>Sampling and reservoir</summary>
      <p class="small-note">Default runs are short screening runs. Increase sampling and compare independent seeds before interpreting material trends.</p>
      <div class="input-row"><label class="field"><span>Production cycles</span><input id="ads-cycles" type="number" value="2000" min="100" max="100000" step="100"></label><label class="field"><span>Equilibration cycles</span><input id="ads-initial" type="number" value="1000" min="100" max="100000" step="100"></label></div>
      <div class="input-row"><label class="field"><span>Widom cycles</span><input id="ads-widom" type="number" value="2000" min="100" max="100000" step="100"></label><label class="field"><span>Random seed</span><input id="ads-seed" type="number" value="42137" min="1" max="2147483000"></label></div>
      <label class="field"><span>Fugacity coefficients, in gas order</span><input id="ads-fugacity" placeholder="1 for each gas (ideal reservoir)"></label>
      <p class="small-note">An empty field uses φ = 1. Supply validated coefficients for nonideal gas conditions. Classical H₂ omits nuclear quantum effects, especially important at low temperature.</p>
    </details>
    <button id="ads-start" class="button primary full">Run adsorption simulation</button><button id="ads-cancel" class="button outline full" hidden>Cancel running simulation</button>
    <p id="ads-progress" class="small-note" role="status" aria-live="polite"></p><div id="ads-results" class="science-results"></div>`;
  $('science-bet').before(pane);
  let active=null,epoch=0,engine=null;
  const num=(v,d=3)=>Number.isFinite(v)?Number(v).toLocaleString(undefined,{maximumSignificantDigits:d}):'Unavailable';
  const node=(tag,text,cls)=>{const e=document.createElement(tag);e.textContent=text;if(cls)e.className=cls;return e;};
  const values=value=>value.split(/[\s,;]+/).filter(Boolean).map(Number);
  function sync(){$('ads-start').disabled=Boolean(active)||!C.structure||!engine?.available;$('ads-cancel').hidden=!active;$('ads-fraction-field').hidden=!$('ads-gas').value.includes(',');}
  function table(root,headers,rows){const wrap=node('div','','materials-table'),table=document.createElement('table'),head=document.createElement('tr');headers.forEach(v=>head.append(node('th',v)));table.append(head);for(const row of rows){const tr=document.createElement('tr');row.forEach(v=>tr.append(node('td',String(v))));table.append(tr);}wrap.append(table);root.append(wrap);}
  function plot(root,result){
    const entries=result.isotherm||[],gases=result.inputs.gases,points=entries.flatMap(r=>(r.components||[]).map(c=>({p:r.pressure_Pa/1e5,q:c.loading_mol_kg,gas:c.gas})));
    if(!points.length)return;const maxX=Math.max(...points.map(p=>p.p)),maxY=Math.max(...points.map(p=>Number.isFinite(p.q)?p.q:0));if(maxY<=0)return;
    const ns='http://www.w3.org/2000/svg',svg=document.createElementNS(ns,'svg');svg.setAttribute('viewBox','0 0 380 240');svg.setAttribute('role','img');svg.setAttribute('aria-label','Simulated absolute adsorption isotherm');svg.style.width='100%';
    const add=(tag,attrs,text)=>{const el=document.createElementNS(ns,tag);Object.entries(attrs).forEach(([k,v])=>el.setAttribute(k,v));if(text!==undefined)el.textContent=text;svg.append(el);return el;};
    add('path',{d:'M52 12V192H368',stroke:'#9aafa6',fill:'none'});
    for(let i=0;i<=4;i++){const x=52+316*i/4,y=192-180*i/4;add('text',{x,y:210,'text-anchor':'middle','font-size':10,fill:'#526860'},num(maxX*i/4));add('text',{x:46,y:y+3,'text-anchor':'end','font-size':10,fill:'#526860'},num(maxY*i/4));}
    add('text',{x:210,y:232,'text-anchor':'middle','font-size':11},'Pressure / bar');add('text',{x:12,y:110,transform:'rotate(-90 12 110)','text-anchor':'middle','font-size':11},'Absolute loading / mol kg⁻¹');
    gases.forEach((gas,i)=>{const rows=points.filter(p=>p.gas===gas&&Number.isFinite(p.q)).sort((a,b)=>a.p-b.p),color=['#087f8c','#b45828','#6b56a8'][i];add('path',{d:rows.map((p,j)=>`${j?'L':'M'}${52+316*p.p/maxX},${192-180*p.q/maxY}`).join(' '),stroke:color,fill:'none','stroke-width':1.8});rows.forEach(p=>add('circle',{cx:52+316*p.p/maxX,cy:192-180*p.q/maxY,r:3,fill:color}));add('text',{x:65+i*90,y:23,'font-size':11,fill:color},gas);});root.append(svg);
  }
  function display(job){
    const root=$('ads-results');root.replaceChildren();const result=job.result;if(!result)return;
    C.note(root,result.scope);
    const settings=result.inputs,manifest=result.widom?.manifest;
    C.note(root,`${num(settings.temperature_K,5)} K · ${settings.gases.map((gas,i)=>`${gas} ${num(settings.mole_fractions[i]*100)}%`).join(' / ')} · ${settings.charge_mode==='zero'?'Uncharged framework approximation':'Supplied CIF charges'}.`);
    C.note(root,`${engine?.profile_name||manifest?.profile_id||'Documented force-field profile'} · ${settings.fugacity_coefficients?'Supplied fugacity coefficients':'Ideal gas reservoir (φ = 1)'}.`);
    if(settings.charge_mode==='zero')C.note(root,'Host–guest electrostatics from framework partial charges are omitted. Guest charge models remain active.',true);
    if(settings.gases.includes('H2'))C.note(root,'Classical H₂ model: nuclear quantum effects are omitted, especially limiting at low temperature.',true);
    plot(root,result);
    table(root,['P / bar','Gas','mol kg⁻¹','95% error'],(result.isotherm||[]).flatMap(r=>(r.components||[]).map(c=>[num(r.pressure_Pa/1e5),c.gas,num(c.loading_mol_kg),num(c.loading_confidence95_mol_kg)])));
    if(result.widom)table(root,['Gas','K_H / mol kg⁻¹ Pa⁻¹','95% error'],(result.widom.components||[]).map(c=>[c.gas,num(c.henry_mol_kg_Pa),num(c.henry_confidence95_mol_kg_Pa)]));
    if(result.zero_loading_heat){table(root,['Gas','Qst⁰ / kJ mol⁻¹','Input-interval envelope'],result.zero_loading_heat.components.map(c=>[c.gas,num(c.qst0_kJ_mol),c.input_interval_envelope_kJ_mol?.map(v=>num(v)).join(' – ')||'Unresolved']));C.note(root,result.zero_loading_heat.scope);}
    if(result.inputs.gases.length>1){table(root,['P / bar','Pair','Selectivity'],(result.isotherm||[]).flatMap(r=>(r.selectivity||[]).map(c=>[num(r.pressure_Pa/1e5),`${c.numerator}/${c.denominator}`,num(c.selectivity)])));table(root,['Pair','Henry selectivity'],(result.henry_selectivity||[]).map(c=>[`${c.numerator}/${c.denominator}`,num(c.selectivity)]));}
    for(const warning of result.warnings||[])C.note(root,warning,true);
    const buttons=node('div','','materials-downloads');for(const [ext,path] of Object.entries(job.downloads||{})){const b=node('button',ext==='zip'?'Save complete run':`Save ${ext.toUpperCase()}`,'button outline');b.addEventListener('click',()=>S.saveExport(path,ext));buttons.append(b);}root.append(buttons);C.details(root,result);
  }
  async function poll(job){
    if(active!==job)return;
    try{const result=await S.api(`science/adsorption/jobs/${job.id}`);if(active!==job)return;
      $('ads-progress').textContent=`${result.status} · ${result.completed_stages}/${result.total_stages} stages · ${result.stage}. ${result.message||''}`;
      if(['complete','failed','cancelled','interrupted'].includes(result.status)){active=null;if(job.epoch===epoch)display(result);else $('ads-progress').textContent='Previous simulation stopped after the model or settings changed. Start a new run for the current structure.';sync();return;}
      setTimeout(()=>poll(job),900);
    }catch(error){$('ads-progress').textContent=`Cannot read simulation status: ${error.message}. Retrying…`;setTimeout(()=>poll(job),2000);}
  }
  async function cancel(){if(!active)return;active.cancelRequested=true;if(!active.id){$('ads-progress').textContent='Cancellation requested; waiting for the queued job identifier…';return;}try{await S.api(`science/adsorption/jobs/${active.id}/cancel`,{});}catch(error){$('ads-progress').textContent=error.message;}}
  function invalidate(){epoch++;$('ads-results').replaceChildren();if(active)cancel();else $('ads-progress').textContent='';sync();}
  $('ads-cancel').addEventListener('click',cancel);
  $('ads-start').addEventListener('click',async()=>{
    if(active||!C.structure)return;const captured=epoch;active={id:null,epoch:captured};sync();$('ads-results').replaceChildren();
    try{const gases=$('ads-gas').value.split(','),fraction=Number($('ads-fraction').value),phis=values($('ads-fugacity').value);
      const payload={structure:C.structure,gases,mole_fractions:gases.length===1?[1]:[fraction,1-fraction],temperature_K:Number($('ads-temperature').value),pressures_Pa:values($('ads-pressures').value).map(p=>p*1e5),cycles:Number($('ads-cycles').value),init_cycles:Number($('ads-initial').value),widom_cycles:Number($('ads-widom').value),calculate_heat:$('ads-heat').checked,seed:Number($('ads-seed').value),charge_mode:$('ads-charges').value,fugacity_coefficients:phis.length?phis:null};
      const result=await S.api('science/adsorption/jobs',payload),cancelRequested=active?.cancelRequested;active={id:result.job_id,epoch:captured};if(captured!==epoch||cancelRequested)await cancel();poll(active);
    }catch(error){active=null;$('ads-progress').textContent=error.message;sync();}
  });
  for(const input of pane.querySelectorAll('input,select'))input.addEventListener('input',invalidate);
  window.addEventListener('cofdraw:structure',invalidate);
  S.api('science/adsorption/engine').then(result=>{engine=result;$('ads-engine').textContent=result.available?`RASPA ${result.version||result.engine_version||'3'} · local CPU · model details are retained in every run`:`Engine unavailable: ${result.reason||'Check the installation'}`;sync();}).catch(error=>{$('ads-engine').textContent=error.message;sync();});
  sync();
})();
