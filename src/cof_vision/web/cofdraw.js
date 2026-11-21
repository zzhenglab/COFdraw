"use strict";
const $ = id => document.getElementById(id);
const NS = "http://www.w3.org/2000/svg";
const L = COFdrawLibrary;
const state = {online: false, busy: false, recipe: null, sourceId: null, sourceMetadata: null, result: null, dirty: false, inputsDirty: false, jsonPending: false, style: "standard", styleCategory:"all", styleOptions:{}, resolvedStyleOptions:{}, seed: null, zoom: 1, library: null, monomers: {}, layoutScale: 1, layoutCenter: [0, 0], autoChemistry:true};
const fallbackPresets = [{id:"imine",label:"Imine · TFB + TAPB",linkage:"imine",monomer_a:"O=Cc1cc(C=O)cc(C=O)c1",monomer_b:"Nc1ccc(-c2cc(-c3ccc(N)cc3)cc(-c3ccc(N)cc3)c2)cc1",label_a:"TFB",label_b:"TAPB",instances:6},{id:"keto",label:"Ketoenamine · Tp + Pa",linkage:"beta_ketoenamine",monomer_a:"O=Cc1c(O)c(C=O)c(O)c(C=O)c1O",monomer_b:"Nc1ccc(N)cc1",label_a:"Tp",label_b:"Pa",instances:12}];
const previewURLs = {};
let combinations=[], combinationFilter=null;
function clone(value) { return JSON.parse(JSON.stringify(value)); }
function sourceMode(){return Boolean((state.sourceId||state.recipe?.explicit_graph)&&!state.inputsDirty);}
function linkageInfo(){return state.library?.linkages?.[$("linkage").value]||{};}
function templates(){const value=state.library?.templates;return Array.isArray(value)?value:value?Object.entries(value).map(([id,t])=>({id,...t})):[];}
function selectedTemplate(){return templates().find(t=>t.id===$("template").value);}
function checkedMonomer(block){
  const value=state.monomers[block];
  if(!value||value.input!==$('smiles-'+block).value.trim())return null;
  const named=state.library?.monomer_presets?.find(row=>row.canonical_smiles===value.smiles);
  return named?{...named,...value,geometry_hint:named.geometry_hint||value.geometry_hint,binary_roles:named.binary_roles,binary_eligible:named.binary_eligible,canonical_smiles:value.smiles}:value;
}
function checkedMonomers(){return {a:checkedMonomer("a"),b:checkedMonomer("b")};}
function matchedCrystalStarter(a,b){
  if(!a?.smiles||!b?.smiles||!state.library)return null;
  const pair=[a.smiles,b.smiles].sort().join('.');
  for(const row of crystalRows){
    const p=row.starter;
    if(!p?.common_name||p.monomer_ids.length!==2)continue;
    const smiles=p.monomer_ids.map(id=>state.library.monomer_presets.find(m=>m.id===id)?.canonical_smiles);
    if(smiles.every(Boolean)&&smiles.sort().join('.')===pair)return row;
  }
  return null;
}
function engineOrder(){const pair=checkedMonomers(),match=sourceMode()?null:L.compatibleLinkages(state.library||{},pair.a,pair.b).find(row=>row.id===$("linkage").value);return match?.swapped?['b','a']:['a','b'];}
function templateCheck(template=selectedTemplate()){if(!template)return {compatible:true};const order=engineOrder(),check=L.templateCompatibility(template,$("linkage").value,{a:checkedMonomer(order[0]),b:checkedMonomer(order[1])},state.library||{});if(order[0]==='b'&&check.reason)check.reason=check.reason.replace(/Block ([AB])/g,(_,block)=>`Block ${block==='A'?'B':'A'}`);return check;}
function option(value,label){const item=document.createElement("option");item.value=value;item.textContent=label;return item;}
function roleLabel(role,index){const info=linkageInfo();return info.role_labels?.[index]||role?.replaceAll("_"," ")||"Reactive building block";}
function periodicLinkage(a,b){
  if(!a||!b)return null;
  const supported=['imine','secondary_amine','acylhydrazone'];
  const selected=$('linkage').value;
  for(const linkage of [...new Set([selected,...(state.autoChemistry?supported:[])])]){
    if(!supported.includes(linkage))continue;
    const roles=state.library?.linkages?.[linkage]?.roles;
    if(roles?.length!==2)continue;
    for(const [node,link] of [[a,b],[b,a]]){
      if(node.geometry_hint!=='tetrahedral')continue;
      if(roles.some((r,i)=>L.reactiveCount(node,linkage,r)===4&&
        (L.reactiveCount(link,linkage,roles[1-i])===2||link.geometry_hint==='tetrahedral'&&L.reactiveCount(link,linkage,roles[1-i])===4)))return linkage;
    }
  }
  return null;
}
function monomerChoice(item,block,other){
  
  if(block!=='editor'&&periodicLinkage(item,other))return {compatible:true,reason:'This pair supports a tetrahedral periodic framework.'};
  if(block!=='editor'&&item.binary_eligible===false)return {compatible:false,reason:item.qualification||'Outside the supported pair constructor.'};
  if(block==='b'&&!other&&$("smiles-a").value.trim())return {compatible:false,reason:'Building block A is being checked or has an invalid structure. Update its preview before selecting B.'};
  if(!state.autoChemistry&&block!=='editor'&&other){const pair=block==='a'?[item,other]:[other,item],match=L.compatibleLinkages(state.library||{},...pair).find(row=>row.id===$("linkage").value);if(match)return {compatible:true,reason:`Compatible with this partner for ${match.label}. Your A/B choices stay in place; construction assigns the reactive roles.`};}
  return state.autoChemistry?L.automaticChoice(state.library||{},item,block,other):L.monomerChoice(state.library||{},item,$("linkage").value,block,other);
}
function assemblyMode(){return $('assembly-count').value;}
const customAssemblyRows=new Map();
function multi(){return assemblyMode()!=='2';}
function assemblyExamples(){
  return (state.library?.assembly_examples||[]).filter(e=>e.monomer_ids.length===Number(assemblyMode()));
}
function assemblyStructures(items){
  const root=$('assembly-structures'), key=JSON.stringify(items.map(r=>r?.smiles));
  if(root.dataset.selection===key)return;
  root.dataset.selection=key;root.replaceChildren();
  for(const [i,item] of items.entries()){
    if(!item)continue;
    const fig=document.createElement('figure'),img=document.createElement('img'),label=document.createElement('figcaption');
    img.alt=L.displayName(item);label.textContent=item.abbreviation||item.label;fig.title=img.alt;
    const edit=document.createElement('button');edit.type='button';edit.className='text-button';edit.textContent='Draw / edit';
    const block='abc'[i];edit.setAttribute('aria-label','Edit precursor '+block.toUpperCase());
    edit.onclick=()=>openSketch(block,item,true);
    fig.append(img,label,edit);root.append(fig);
    libraryPreviews.request(item.smiles).then(preview=>{if(img.isConnected)img.src=preview.url;}).catch(()=>img.hidden=true);
  }
}
function fillAssembly(saved=null){
  const count=assemblyMode(), mixed=count==='3', rows=[...(state.library?.monomer_presets||[]),...customAssemblyRows.values()];
  $('binary-fields').hidden=multi();$('assembly-fields').hidden=!multi();
  $('mixed-fields').hidden=!mixed;
  if(!multi())return;
  const one=rows.filter(r=>r.assembly_roles?.self_kind||r.assembly_roles?.self_imine);
  const countFor=(r,role)=>r.assembly_roles?.[role]||0;
  const nodes=rows.filter(r=>!r.assembly_roles?.self_imine&&!r.geometry_hint?.startsWith('tetrahedral')&&['aryl_aldehyde','primary_arylamine'].some(k=>[3,4].includes(countFor(r,k))));
  function choose(key,values,fallback){
    const select=$('assembly-'+key), previous=saved?.monomer_ids?.['abc'.indexOf(key)]||select.value;
    select.replaceChildren(...values.map(r=>option(r.id,L.displayName(r))));
    select.value=[previous,fallback,values[0]?.id].find(id=>values.some(r=>r.id===id))||'';
    return values.find(r=>r.id===select.value);
  }
  $('assembly-a-label').textContent=mixed?'Node A':'Precursor';
  const node=choose('a',mixed?nodes:one,mixed?'tape':'bdba');
  if(!node)return;
  let chosen=[node];
  const faces=mixed?[]:node.assembly_roles.self_pores||[];
  const face=$('self-face'), desiredFace=saved?.pore_view||face.value||'hexagon';
  face.replaceChildren(...faces.map(f=>option(f.id,f.label)));
  face.value=faces.some(f=>f.id===desiredFace)?desiredFace:faces.some(f=>f.id==='hexagon')?'hexagon':faces[0]?.id||'';
  $('self-face-field').hidden=faces.length<2;
  if(mixed){
    const nodeAmine=countFor(node,'primary_arylamine')>=3;
    const other=nodeAmine?'aryl_aldehyde':'primary_arylamine';
    const layout=$('mixed-layout'),hasKagome=assemblyExamples().some(e=>e.pore_layout==='kgm'&&e.monomer_ids[0]===node.id);
    $('mixed-layout-field').hidden=!hasKagome;
    if(saved||!hasKagome)layout.value=hasKagome&&saved?.pore_layout==='kgm'?'kgm':'local';
    const kgm=layout.value==='kgm', n=countFor(node,nodeAmine?'primary_arylamine':'aryl_aldehyde')===3?6:4;
    const links=rows.filter(r=>r.id!==node.id&&!r.assembly_roles?.self_imine&&countFor(r,other)===2&&(!kgm||r.assembly_roles.kgm_linear_dialdehyde));
    const b=choose('b',links,nodeAmine?'ta':'pa');
    const c=choose('c',links.filter(r=>r.canonical_smiles!==b?.canonical_smiles),nodeAmine?'bpda_cho':'benzidine');
    chosen=[node,b,c];
    const pattern=$('assembly-pattern'),arr=$('assembly-arrangement');
    const choices=kgm?[{sequence:'BCBCBC',label:'Hexagon · alternating B and C'},{sequence:'BBB',label:'Triangle · three B linkers'},{sequence:'CCC',label:'Triangle · three C linkers'}]:state.library?.mixed_pore_patterns?.[String(n)]||[];
    const current=saved?.linker_sequence??pattern.value;
    const free=!kgm&&arr.value==='custom'&&pattern.dataset.length===String(n)&&!saved;
    arr.replaceChildren(...choices.map(p=>option(p.sequence,p.label)));
    if(!kgm)arr.add(option('custom','Enter a clockwise sequence…'));
    const valid=kgm?choices.some(p=>p.sequence===current):current.length===n&&/^[BC]+$/.test(current)&&current.includes('B')&&current.includes('C');
    pattern.value=free||valid?current:kgm?'BCBCBC':'BC'.repeat(n/2);
    pattern.dataset.length=String(n);
    arr.value=free?'custom':choices.some(p=>p.sequence===pattern.value)?pattern.value:'custom';
    $('sequence-field').hidden=arr.value!=='custom';
    $('sequence-note').textContent=kgm?'The two triangles and the hexagon share the same network.': 'Letters give the clockwise order of the linkers. Rotations and mirror images are grouped.';
  }
  assemblyStructures(chosen);
  const examples=assemblyExamples(), selected=examples.find(e=>e.monomer_ids.every((id,i)=>id===chosen[i]?.id)&&(!mixed||(e.pore_layout==='kgm'?$('mixed-layout').value==='kgm':$('mixed-layout').value!=='kgm'&&e.linker_sequence===$('assembly-pattern').value)));
  const menu=$('assembly-example');menu.replaceChildren(option('','Custom selection'),...examples.map(e=>option(e.id,e.label)));menu.value=selected?.id||'';
  const evidence=$('assembly-evidence');evidence.replaceChildren();evidence.hidden=!selected?.doi;
  if(selected?.doi){const a=document.createElement('a');a.href='https://doi.org/'+selected.doi;a.target='_blank';a.rel='noopener noreferrer';a.textContent='Source paper';evidence.append(a);}
  $('assembly-note').textContent=selected?.description||(mixed?'A selected linker sequence defines a local design. Bulk order depends on how the COF forms.':'A finite drawing from the selected precursor.');
}
function assemblyReady(){
  if(!$('assembly-a').value)return false;
  if(assemblyMode()==='1')return true;
  const p=$('assembly-pattern');
  const valid=$('mixed-layout').value==='kgm'?['BBB','CCC','BCBCBC'].includes(p.value):p.value.length===Number(p.dataset.length)&&/^[BC]+$/.test(p.value)&&p.value.includes('B')&&p.value.includes('C');
  return Boolean($('assembly-b').value&&$('assembly-c').value&&valid);
}
function assemblyFeeds(recipe){
  if(!recipe?.assembly)return null;
  if(recipe.assembly.precursors)return recipe.assembly.precursors;
  return Object.values(recipe.components).map(c=>({label:c.label,smiles:c.monomer_smiles}));
}
function restoreAssembly(recipe){
  const saved=recipe?.assembly;
  $('assembly-count').value=saved?String(saved.monomer_ids.length):'2';
  fillAssembly(saved);
  const descriptor=$('descriptor-block'),previous=descriptor.value;
  const feeds=assemblyFeeds(recipe),keys=feeds?feeds.map((_,i)=>'ABC'[i]):['A','B'];
  descriptor.replaceChildren(...keys.map((k,i)=>option(k.toLowerCase(),feeds?feeds[i].label:'Building block '+k)));
  if(keys.includes(previous.toUpperCase()))descriptor.value=previous;
}
function selectedDescriptor(){
  if(window.COFdrawCrystalDesign?.active){
    const item=window.COFdrawCrystalDesign.precursor($('descriptor-block').value);
    return {label:item?.label||'Unassigned precursor',smiles:item?.smiles||''};
  }
  const block=$('descriptor-block').value, feed=assemblyFeeds(state.recipe)?.['abc'.indexOf(block)];
  return feed||{label:$('label-'+block).value,smiles:$('smiles-'+block).value.trim()};
}
function assemblyChanged(saved=null){
  clearTimeout(liveDrawingTimer);selectionRevision++;
  $('cof-preset').value='custom';
  if(multi())window.COFdrawCrystalDesign?.leave();
  state.inputsDirty=true;state.dirty=true;state.sourceId=null;state.sourceMetadata=null;
  fillAssembly(saved);sync();updateIdentity();queueDrawing();
}
$('assembly-count').addEventListener('change',()=>assemblyChanged(assemblyExamples()[0]||null));
for(const id of ['assembly-a','assembly-b','assembly-c','self-face','mixed-layout'])$(id).addEventListener('change',()=>assemblyChanged());
$('assembly-example').addEventListener('change',()=>assemblyChanged(assemblyExamples().find(e=>e.id===$('assembly-example').value)||null));
$('assembly-arrangement').addEventListener('change',()=>{const value=$('assembly-arrangement').value;if(value!=='custom')$('assembly-pattern').value=value;assemblyChanged();if(value==='custom')$('assembly-pattern').focus();});
$('assembly-pattern').addEventListener('input',()=>{$('assembly-pattern').value=$('assembly-pattern').value.toUpperCase();assemblyChanged();});
$('assembly-build').onclick=()=>$('build').click();
function constructionCheck(){
  if(!checkedMonomer('a')||!checkedMonomer('b'))return {compatible:false,reason:'Choose both building blocks and wait for their structure checks before creating a drawing.'};
  const info=linkageInfo();
  if(info.supported!==true)return {compatible:false,reason:'This linkage is reference-only; no automatic construction is implemented.'};
  const order=engineOrder();for(const [index,block] of order.entries()){
    const evidence=L.roleEvidence(checkedMonomer(block),$("linkage").value,L.roleFor(state.library||{},$("linkage").value,index===0?'a':'b'));
    if(evidence?.compatible!==true)return {compatible:false,reason:`Building block ${block.toUpperCase()} is incompatible with the selected linkage. Choose a compatible partner or change the linkage.`};
  }
  return templateCheck();
}
const recentTreatments=[];
let liveDrawingTimer, selectionRevision=0;
let standardPair='';
function currentPairKey(){return JSON.stringify(['a','b'].map(k=>$('smiles-'+k).value.trim()));}
window.addEventListener('cofdraw:prefer-standard',()=>{standardPair=currentPairKey();});
function queueDrawing(){
  clearTimeout(liveDrawingTimer);
  if(window.COFdrawCrystalDesign?.active){
    if(window.COFdrawCrystalDesign.active==='custom-dia')window.COFdrawCrystalDesign.inputsChanged();
    return;
  }
  const pair=JSON.stringify([$("smiles-a").value,$("smiles-b").value]);
  liveDrawingTimer=setTimeout(()=>{
    if(pair!==JSON.stringify([$("smiles-a").value,$("smiles-b").value])||state.jsonPending||sourceMode())return;
    if(state.busy){queueDrawing();return;}
    if((state.inputsDirty||!state.recipe)&&!$("build").disabled)$("build").click();
  },250);
}
function autoMatchChemistry(){
  if(window.COFdrawCrystalDesign?.active)return;
  if(!state.autoChemistry||sourceMode())return;
  const pair=checkedMonomers();if(!pair.a||!pair.b){$("chemistry-match").textContent='Choose both building blocks and wait for their structure checks.';return;}
  const selected=combinations.find(p=>p.id===$("cof-preset").value && p.monomer_a===pair.a.input && p.monomer_b===pair.b.input);
  const molecularIDs=[pair.a,pair.b].map(m=>state.library.monomer_presets.find(p=>p.canonical_smiles===m.smiles)?.id);
  const preferred=selected||state.library.named_pairs?.find(p=>molecularIDs.includes(p.a)&&molecularIDs.includes(p.b));
  const match=L.automaticPlan(state.library,pair.a,pair.b,preferred);
  if(!match){$("chemistry-match").textContent="No supported transform matches this pair. Change one monomer or load a source drawing.";return;}
  populateLinkages(match.id);
  const template=templates().find(item=>item.id===match.template_id);
  $("chemistry-match").textContent=`${match.label} · ${template?.label||'finite patch'}`;
  return match.template_id;
}
function populateLinkages(desired=$("linkage").value){
  const search=$("linkage-search").value.trim().toLowerCase().replaceAll("_"," "),rows=Object.entries(state.library?.linkages||{});
  const select=$("linkage");select.replaceChildren();const groups=new Map();
  for(const[id,item]of rows){if(search&&id!==desired&&!JSON.stringify([id,item.label,item.category,item.roles,item.description,item.keywords,item.search_terms,item.conceptual_motif]).toLowerCase().replaceAll("_"," ").includes(search))continue;
    const category=item.category||"Linkage chemistry";if(!groups.has(category)){const group=document.createElement("optgroup");group.label=category;groups.set(category,group);select.append(group);}
    const o=option(id,`${item.label||id.replaceAll("_"," ")}${item.supported===false?" · reference only":""}`);groups.get(category).append(o);
  }
  if(rows.some(([id])=>id===desired))select.value=desired;
  if(!select.value&&select.options.length)select.selectedIndex=0;
}
function refreshPresetLists(){
  if(!state.library)return;
  for(const block of["a","b"]){const select=$(`preset-${block}`);select.replaceChildren(option("","Custom / pasted structure"));
    for(const item of state.library.monomer_presets||[]){
      if(item.binary_eligible===false)continue;
      const choice=monomerChoice(item,block,checkedMonomer(block==="a"?"b":"a")),selected=item.smiles===$(`smiles-${block}`).value.trim();
      if(block==='b'&&$("smiles-a").value.trim()&&!choice.compatible&&!selected)continue;
      const o=option(item.id,`${L.displayName(item)}${choice.compatible?"":" — unavailable"}`);o.disabled=!choice.compatible;o.title=choice.reason;select.append(o);
    }
    const matched=(state.library.monomer_presets||[]).find(item=>item.smiles===$(`smiles-${block}`).value.trim());
    select.value=matched?.id||"";
  }
  refreshCombinations();
}
let matchingCombinations=[], crystalRows=[], suggestionIndex=-1;
const comboInput=$('combination-search'),suggestions=$('combination-suggestions');
window.addEventListener('cofdraw:crystal-catalog',e=>{
  crystalRows=e.detail;
  refreshCombinations(true);
  if(document.activeElement===comboInput)showSuggestions();
});
function crystalStarts(){
  const blocks=new Map((state.library?.monomer_presets||[]).map(r=>[r.id,r]));
  const names={tbpm_boronic:'TBPM',tbps_boronic:'TBPS'};
  const rows=[];
  for(const row of crystalRows){
    const p=row.starter;
    if(!p)continue;
    const labels=p.monomer_ids.map(id=>names[id]||blocks.get(id)?.abbreviation||id);
    const commonName=p.common_name,feed=labels.join(' + '),net=p.net,reference=p.kind==='reference';
    const label=feed+(commonName?' ('+commonName+')':'');
    rows.push({id:'crystal:'+row.id,crystal:true,net,reference,commonName,label,precursors:p.monomer_ids,
      search:[label,row.label,row.model_name,net,p.doi,...p.monomer_ids,reference?'reference CIF':'generated design'].filter(Boolean).join(' ')});
  }
  return rows;
}
function closeSuggestions(){
  suggestions.hidden=true;suggestionIndex=-1;
  comboInput.setAttribute('aria-expanded','false');comboInput.removeAttribute('aria-activedescendant');
}
function showSuggestions(){
  closeSuggestions();suggestions.replaceChildren();
  for(const row of matchingCombinations.slice(0,12)){
    const item=document.createElement('div');
    item.id='suggestion-'+suggestions.childElementCount;item.dataset.value=row.id;item.setAttribute('role','option');
    item.setAttribute('aria-selected','false');
    const name=document.createElement('strong'),group=document.createElement('small');
    name.textContent=row.label;
    group.textContent=row.crystal?'3D · '+row.net+(row.reference?' · reference CIF':' · editable design'):row.source?'Recorded drawing':row.assembly?'Assembly':state.library?.linkages?.[row.linkage]?.label||'Drawing';
    item.append(name,group);suggestions.append(item);
  }
  if(suggestions.childElementCount){suggestions.hidden=false;comboInput.setAttribute('aria-expanded','true');}
}
function selectSuggestion(item){
  if(!item)return;
  const id=item.dataset.value;closeSuggestions();comboInput.value='';
  refreshCombinations(true,id);
  $('cof-preset').value=id;$('cof-preset').dispatchEvent(new Event('change'));
}
function refreshCombinations(force=false,desired=null){
  const query=comboInput.value.trim(), linkage=state.autoChemistry?'':$('linkage').value;
  const select=$('cof-preset'),current=desired||(window.COFdrawCrystalDesign?.active?'crystal:'+window.COFdrawCrystalDesign.active:state.sourceId?'source:'+state.sourceId:select.value);
  const custom=Boolean(periodicLinkage(checkedMonomer('a'),checkedMonomer('b'))),key=JSON.stringify([query,linkage,current,custom]);
  if(!force&&combinationFilter===key)return;
  combinationFilter=key;
  const binary=combinations.filter(row=>row.source||row.assembly||(linkage?row.linkage===linkage:row.linkage!=='secondary_amine'));
  const {visible,matches}=L.starterSections([...binary,...crystalStarts()],state.library||{},query,current);
  const fragment=document.createDocumentFragment();
  fragment.append(option(current==='crystal:custom-dia'?current:'custom','Custom monomers'));
  for(const [name,rows] of visible){
    const group=document.createElement('optgroup');group.label=name;
    for(const row of rows)group.append(option(row.id,row.label));
    fragment.append(group);
  }
  select.replaceChildren(fragment);select.value=[...select.options].some(o=>o.value===current)?current:'custom';
  matchingCombinations=matches;
  $('combination-count').textContent=!matches.length?'No matches. Try a COF, monomer or linkage name.':query?matches.length+' matched choices':'Search all '+matches.length.toLocaleString()+' choices';
}
function refreshTemplates(desired=null){
  if(sourceMode()){$("template-field").hidden=true;$("pore-field").hidden=true;$("template-description").textContent="This drawing uses its recorded connections and layout. Edit a monomer to begin a new construction.";return;}
  const rows=templates(),available=rows.filter(t=>!t.linkages||t.linkages.includes($("linkage").value));
  $("template-field").hidden=!rows.length;$("pore-field").hidden=Boolean(rows.length);
  if(!rows.length){$("template-description").textContent="Alternating monomers with unconnected outward continuation ports.";return;}
  const select=$("template"),current=desired||select.value;select.replaceChildren();
  for(const t of available){const check=templateCheck(t),o=option(t.id,t.label||t.id);o.disabled=!check.compatible;o.title=check.reason||t.description||"";select.append(o);}
  const validCurrent=Array.from(select.options).find(o=>o.value===current&&!o.disabled);const first=Array.from(select.options).find(o=>!o.disabled);
  select.value=validCurrent?.value||first?.value||available[0]?.id||"";
  if(!select.options.length)select.append(option("","No compatible template available"));
  const item=selectedTemplate(),check=templateCheck();$("template-description").textContent=item?[item.description,check.compatible?item.limits:check.reason].filter(Boolean).join(" "):"Choose supported linkage chemistry to view its templates.";
  if(item?.instances)$("pore").value=String(item.instances);
}
function refreshInspectors(){
  const source=sourceMode(),info=linkageInfo(),roles=info.roles||[];$("linkage-roles").textContent=source?(state.sourceId?"Source-associated precursors":"Recorded graph — associated precursors"):roles.length?`${roleLabel(roles[0],0)} + ${roleLabel(roles[1],1)}`:(info.description||info.conceptual_motif||"Select a linkage category.");
  for(const[block,index]of[["a",0],["b",1]]){
    const roleIndex=engineOrder().indexOf(block),value=checkedMonomer(block),evidence=source?null:L.roleEvidence(value,$("linkage").value,roles[roleIndex]);$(`role-${block}`).textContent=source?`${state.sourceId?'Source':'Associated'} precursor ${block.toUpperCase()}`:`Required: ${roleLabel(roles[roleIndex],roleIndex)}`;
    const named=state.library?.monomer_presets?.find(item=>item.smiles===$(`smiles-${block}`).value.trim()||value?.smiles&&item.smiles===value.smiles);$(`preset-${block}`).title=named?L.displayName(named):"Choose a compatible monomer, or paste a custom structure.";
    let fullName=$(`name-${block}`);if(!fullName){fullName=document.createElement("p");fullName.id=`name-${block}`;fullName.className="monomer-full-name";$(`inspector-${block}`).prepend(fullName);}fullName.textContent=named?L.displayName(named):"";fullName.hidden=!named;
    const chips=$(`groups-${block}`);chips.replaceChildren();for(const group of L.groupRows(value)){const chip=document.createElement("span");chip.className="group-chip";chip.textContent=`${group.count} × ${group.label||group.kind?.replaceAll("_"," ")||"reactive group"}`;chips.append(chip);}
    $(`inspector-${block}`).classList.toggle("incompatible",evidence?.compatible===false);
    $(`compatibility-${block}`).textContent=source?($(`smiles-${block}`).value.trim()?"Associated precursor for the recorded graph. An automatic construction transform is not implied.":"No associated precursor structure is supplied for this part of the graph."):(!value?"Update preview to inspect reactive functionality.":evidence?.compatible===false?evidence.reason||"This structure cannot supply the selected role.":evidence?`${evidence.reactive_sites??evidence.count??"Supported"} reactive sites for this role. Attachment positions are checked during construction.`:"Structure valid. Reactive compatibility will be checked during construction.")+(value?.geometry_note?` ${value.geometry_note}`:"");
    if(window.COFdrawCrystalDesign?.active){
      const choice=value&&window.COFdrawCrystalDesign.choice(value,checkedMonomer(block==='a'?'b':'a'));
      $(`role-${block}`).textContent=$('reticular-net').value+' framework · '+$('dia-linkage').selectedOptions[0].textContent;
      $(`compatibility-${block}`).textContent=choice?.reason||'Choose or draw a precursor.';
      $(`inspector-${block}`).classList.toggle('incompatible',Boolean(choice&&!choice.compatible));
    }
  }
}
function refreshChemistry(desiredTemplate=null){const automatic=autoMatchChemistry();refreshPresetLists();refreshInspectors();refreshTemplates(state.autoChemistry?automatic:desiredTemplate);updateIdentity();sync();if($("monomer-library-dialog").open)drawLibrary();}
function notice(text, error=false) { $("notice").textContent=text; $("notice").className=`notice${error ? " error" : ""}`; }
async function api(path, payload) {
  const response=await fetch(`./api/drawing/${path}`,payload===undefined?{}:{method:"POST",headers:{"Content-Type":"application/json"},body:JSON.stringify(payload)});
  let value; try { value=await response.json(); } catch { throw new Error("The local drawing engine returned an unreadable response."); }
  if(!response.ok) throw new Error(typeof value.detail==="string"?value.detail:Array.isArray(value.detail)?value.detail.map(d=>d.msg).join("; "):"The drawing could not be processed.");
  return value;
}
function sync() {
  const source=sourceMode(),linkage=linkageInfo();
  $("chemistry-note").textContent=source?`${state.sourceId?'Source':'Recorded'} graph loaded. Change presentation or edit the recorded layout; change a precursor to start a new construction.`:linkage.limits||"The pore sequence is an explicit design choice. It does not assign a periodic topology.";
  const hideChoices=source||state.autoChemistry;
  document.querySelector('.linkage-first').hidden=hideChoices;
  document.querySelector('.chemistry-options').hidden=hideChoices;
  for(const id of ['linkage-roles','template-description'])$(id).hidden=hideChoices;
  document.querySelector('.template-scope').hidden=hideChoices;
  $('chemistry-mode-field').hidden=source;$("combination-search").disabled=state.busy;
  if(source)$("chemistry-match").textContent="The recorded recipe supplies this graph directly; constructor compatibility does not apply.";
  else if($("chemistry-match").textContent.startsWith("The recorded recipe supplies"))$("chemistry-match").textContent="Linkages and connectivity update as you choose building blocks.";
  const blocked=!state.online||state.busy;
  const construction=constructionCheck();
  $("build").disabled=blocked||source||!construction.compatible||(templates().length>0&&!selectedTemplate());$("build").title=source?'The recorded graph is already loaded.':construction.compatible?'Create the selected finite chemical patch.':construction.reason; $("random-style").disabled=blocked;
  if(multi()){
    $("build").disabled=blocked||!assemblyReady();
    $("build").title=assemblyReady()?'Draw the selected arrangement.':'Choose valid monomers and an explicit linker order.';
  }
  const crystal=window.COFdrawCrystalDesign;
  const eligible=!multi()&&!crystal?.active&&periodicLinkage(checkedMonomer('a'),checkedMonomer('b'));
  const periodic=$('build-periodic-from-pair');
  const matched=!multi()&&!crystal?.active&&matchedCrystalStarter(checkedMonomer('a'),checkedMonomer('b'));
  const ref=!eligible&&matched?.starter?.kind==='reference'&&(state.autoChemistry||$('linkage').value===matched.starter.linkage);
  periodic.hidden=!eligible&&!ref;periodic.disabled=blocked;periodic.dataset.linkage=eligible||'';
  periodic.dataset.reference=ref?matched.id:'';
  periodic.dataset.net=matched?.starter?.kind==='generated'?matched.starter.net:'dia';
  periodic.textContent=ref?'View '+matched.starter.net+' reference · '+matched.starter.common_name:
    matched?.starter?.kind==='generated'?'Build '+matched.starter.net+' design · '+matched.starter.common_name+' precursors':'Build a 3D framework from this pair';
  if((eligible||ref)&&(!state.recipe||state.inputsDirty&&state.autoChemistry))$('build').disabled=blocked;
  if(crystal?.active){
    $('build').disabled=blocked||!$('smiles-a').value.trim()||!$('smiles-b').value.trim();
    $('build').textContent='Create drawing →';$('build').title='Draw the selected framework.';
    if(crystal.isReference&&crystal.report)$('build').disabled=blocked;
  }else if($('build').textContent==='Build 3D framework')$('build').textContent='Create drawing →';
  $("assembly-build").disabled=$("build").disabled;
  for(const id of ['assembly-count','assembly-example','assembly-a','assembly-b','assembly-c','assembly-pattern','self-face','mixed-layout','assembly-arrangement'])$(id).disabled=blocked;
  $("build-from-canvas").disabled=$("build").disabled;
  $("build-from-canvas").title=$("build").title;
  $("save-project").disabled=window.COFdrawCrystalDesign?.active?!window.COFdrawCrystalDesign.report:!state.recipe||state.busy||state.jsonPending||(multi()&&state.inputsDirty);
  $("render-changes").disabled=blocked||!state.recipe||state.inputsDirty||state.jsonPending;
  $("apply-recipe").disabled=blocked;
  $("render-busy").classList.toggle("hidden",!state.busy);
  $("scope-badge").textContent=state.sourceId?"SOURCE · PROVISIONAL":"DESIGNED PATCH";
  $("scope-badge").classList.toggle("source",Boolean(state.sourceId));
  document.querySelector(".canvas-panel").classList.toggle("is-dirty",state.dirty);
  $("dirty-status").textContent=state.inputsDirty?"Monomers changed · create a new drawing":state.jsonPending?"Apply the edited recipe":state.dirty?"Changes awaiting render":"";
  document.querySelectorAll("#exports button").forEach(b=>b.disabled=crystal?.active?!crystal.visible:!state.result||state.dirty||state.busy);
  document.querySelectorAll(".style-choice").forEach(b=>{b.disabled=blocked;b.classList.toggle("selected",b.dataset.style===state.style);b.setAttribute("aria-pressed",String(b.dataset.style===state.style));});
  document.querySelectorAll(".sketch-open,.validate-monomer").forEach(b=>b.disabled=blocked);
  document.querySelectorAll(".library-open,#open-monomer-library,#new-monomer").forEach(b=>b.disabled=blocked);
  $("open-topology-library").disabled=blocked||!(state.library?.topology_catalog?.length);
  document.querySelectorAll("#style-controls input,#style-controls select,#style-controls button").forEach(b=>b.disabled=blocked);
  $("apply-style-options").disabled=blocked||!state.recipe||state.inputsDirty||state.jsonPending;
  for(const id of ["cof-preset","linkage","linkage-search","pore","template","image-size","label-a","label-b","smiles-a","smiles-b","preset-a","preset-b","recipe-json","open-project"]) $(id).disabled=state.busy;
  $("auto-chemistry").checked=state.autoChemistry;$("auto-chemistry").disabled=state.busy;
  for(const id of ['linkage','linkage-search','template','pore'])$(id).disabled=state.busy||state.autoChemistry;
  $("linkage-mode-label").textContent=state.autoChemistry?'Linkage · Auto':'Linkage · Manual';
  $("template-mode-label").textContent=state.autoChemistry?'Connectivity · Auto':'Connectivity · Manual';
  if(!source&&!state.autoChemistry)$("chemistry-match").textContent='Choose a linkage and a connectivity template.';
  if(crystal?.active){
    $('scope-badge').textContent='3D FRAMEWORK';
    $('dirty-status').textContent='';
    $('measure-pore').disabled=!crystal.report;$('measure-pore').title='Calculate pore dimensions from the current periodic cell';
    document.querySelector('.canvas-panel').classList.remove('is-dirty');
    $('random-style').disabled=!crystal.visible;
    $('apply-style-options').disabled=!crystal.visible;
    for(const button of $('exports').children)button.disabled=!crystal.visible;
    $('canvas-empty').classList.add('hidden');
  }else $('measure-pore').removeAttribute('title');
  window.dispatchEvent(new CustomEvent('cofdraw:state'));
}
function removeIdentity() {
  if(!state.recipe)return;
  const wasSource=Boolean(state.sourceId);
  state.recipe.derived_from_recipe_id=state.recipe.derived_from_recipe_id||state.recipe.recipe_id;
  const r=state.recipe;
  if(!r.source_registration_history&&(r.explicit_graph?.source_coordinates_preserved||r.source_coordinate_provenance||r.layout_evidence||r.source_style?.source_registered))r.source_registration_history={scope:'historical_parent_only',explicit_graph_coordinate_provenance:clone(r.explicit_graph?.coordinate_provenance||null),source_coordinate_provenance:clone(r.source_coordinate_provenance||null),layout_evidence:clone(r.layout_evidence||null),source_style_registration:clone(r.source_style||null)};
  r.source_registered=false;if('source_coordinates_preserved' in r)r.source_coordinates_preserved=false;if('source_representation_kind' in r)r.source_representation_kind='user_edited_finite_patch';if('source_scope' in r)r.source_scope='User-edited finite chemical patch; original source scope is historical only.';
  if(r.explicit_graph){r.explicit_graph.source_coordinates_preserved=false;r.explicit_graph.coordinate_provenance={kind:'user_manual_coordinates',source_registered:false,scope:'Edited design; parent registration is historical evidence only.'};}
  if(r.source_style){r.source_style.source_registered=false;r.source_style.exact_style_equivalence=false;r.source_style.qualification='Inherited presentation on an edited drawing; no source registration asserted.';delete r.source_style.coordinate_evidence_sha256;delete r.source_style.review_evidence_sha256;}
  for(const key of ['source_coordinate_provenance','coordinate_provenance','layout_evidence'])if(r[key])r[key]={kind:'user_manual_coordinates',source_registered:false,scope:'Current design coordinates; see source_registration_history for original evidence.'};
  for(const k of ["cofex_id","cofex_string","cofex_notation","parent_cofex_id","source","reference","topology","topology_evidence","design_cofex_notation","design_topology"])state.recipe[k]=null;
  Object.assign(state.recipe,{designed_variant:true,ground_truth:false,review_status:"user_design"});state.sourceId=null;state.sourceMetadata=null;
  updateIdentity();if(wasSource)refreshInspectors();
}
function markInputsChanged() { selectionRevision++;const wasSource=sourceMode();if(state.recipe){removeIdentity();state.dirty=true;}state.inputsDirty=true;$("cof-preset").value="custom";window.COFdrawCrystalDesign?.inputsChanged();if(wasSource)refreshChemistry();else{refreshInspectors();sync();} }
function syncJSON(){ $("recipe-json").value=state.recipe?JSON.stringify(state.recipe,null,2):"";state.jsonPending=false; }
function validateRecipeShape(recipe) {
  if(recipe?.explicit_graph){if(!Array.isArray(recipe.explicit_graph.atoms)||!Array.isArray(recipe.explicit_graph.bonds))throw new Error("This source drawing needs explicit atoms and bonds.");return;}
  if(!recipe||!recipe.components||!Array.isArray(recipe.instances)||!Array.isArray(recipe.connections)||recipe.instances.length<3||recipe.instances.length>24)throw new Error("A recipe needs components, 3–24 instances, and explicit connections.");
  for(const item of recipe.instances)if(!recipe.components[item.component]||!Array.isArray(item.layout_position)||item.layout_position.length!==2||item.layout_position.some(x=>typeof x!=="number"||!Number.isFinite(x)||Math.abs(x)>200))throw new Error("Each block needs a component and finite x/y coordinates between −200 and 200.");
}
function acceptRecipe(recipe,sourceId=null,metadata=null,engineDesign=false){validateRecipeShape(recipe);restoreAssembly(recipe);state.recipe=clone(recipe);state.sourceId=sourceId;state.sourceMetadata=metadata;state.inputsDirty=false;state.dirty=true;if(!sourceId){const designNotation=state.recipe.design_cofex_notation,designTopology=state.recipe.design_topology;removeIdentity();if(engineDesign){state.recipe.design_cofex_notation=designNotation;state.recipe.design_topology=designTopology;}}syncJSON();drawLayout();drawLayoutTable();updateIdentity();sync();window.dispatchEvent(new CustomEvent('cofdraw:layoutchanged'));}
function updateIdentity(){
  if(window.COFdrawCrystalDesign?.active)return;
  const audit=state.result&&!state.dirty?state.result.audit:null,metadata=state.sourceMetadata||audit?.literature_metadata||{};
  const key=state.sourceId?(metadata.cofex_id||audit?.cofex_id||state.recipe?.cofex_id):null;
  const notation=state.inputsDirty||state.jsonPending?null:state.sourceId?(metadata.cofex_notation||metadata.cofex_string||audit?.cofex_string||state.recipe?.cofex_notation):(audit?.design_cofex_notation||state.recipe?.design_cofex_notation||audit?.drawing_cofex?.compact);
  const pending=!state.recipe||state.inputsDirty;
  $("identity-heading").textContent=state.sourceId?'Source COFex record':'Design COFex record';
  $("identity-linkage").textContent=pending?($("linkage").selectedOptions[0]?.textContent||"Not selected"):(state.library?.linkages?.[state.recipe?.linkage]?.label||state.recipe?.linkage||(sourceMode()?(state.sourceId?"Source-recorded bonds":"Explicit recorded bonds"):$("linkage").selectedOptions[0]?.textContent)||"Unknown");
  $("identity-key").textContent=key||(state.sourceId?"Not supplied by source":"Unassigned design");$("identity-string").textContent=notation||"Not assigned";
  const topology=audit?audit.design_topology:state.recipe?.design_topology,template=pending?selectedTemplate():audit?audit.selected_template:state.recipe?.selected_template;$("identity-template").textContent=state.sourceId?"Source-visible finite patch":template?.label||(typeof topology==="object"?(topology?.label||topology?.template_id||topology?.name):topology)||"Explicit editable connections";
  const cat=notation?.match(/\sCOFex\.[^.\s]+\.(cat\d+)(?:;|$)/)?.[1]||metadata.catenation||null;
  $("identity-cat").textContent=cat?(state.sourceId?`${cat} · source convention; not independently verified`:`${cat} · design assumption; observed interpenetration unknown`):'Unknown · no catenation assignment';
  $("identity-note").textContent=state.recipe?.assembly&&!state.inputsDirty?state.recipe.assembly.scope:state.inputsDirty?'Inputs changed. Create a drawing to update its notation and explicit connectivity.':state.sourceId?(notation?'Full source precursor SMILES and recorded COFex suffix. Provisional finite patch; unknown net codes and unverified catenation remain qualified.':'This source has no complete COFex assignment. Its recorded graph and source provenance are retained; missing identities are not inferred.'):notation?'Full design precursor SMILES with the chosen template and assumed cat0. These are design choices, not an observed material identity or a periodic crystal.':'Choose monomers and create a drawing to obtain a qualified design notation.';
  $("copy-identity").disabled=!state.recipe||state.inputsDirty||state.jsonPending;
}
async function saveText(name,text,type="application/json"){
  const bridge=window.pywebview?.api;
  if(bridge?.save_text){const value=await bridge.save_text(name,text,type);if(value?.saved)notice(`Saved ${name}.`);return;}
  const url=URL.createObjectURL(new Blob([text],{type}));const a=document.createElement("a");a.href=url;a.download=name;a.click();setTimeout(()=>URL.revokeObjectURL(url),3000);
}
function showExports(downloads){
  const bar=$("exports");bar.replaceChildren();
  for(const ext of ["png","svg","mol","sdf","cdxml","json"]){
    if(!downloads?.[ext])continue;
    const button=document.createElement("button");
    button.className="button outline";button.textContent=ext.toUpperCase();button.title="Save "+ext.toUpperCase();
    button.onclick=()=>saveExport(downloads[ext],ext);bar.append(button);
  }
}
async function saveExport(path,ext){
  try{const name=`cofdraw-${ext==="json"?"drawing":"patch"}.${ext}`;const bridge=window.pywebview?.api;if(bridge?.save_export){const r=await bridge.save_export(path,name);if(r?.saved)notice(`Saved ${name}.`);return;}const response=await fetch(`.${path}`);if(!response.ok)throw new Error("This export is unavailable; render the drawing again.");const blob=await response.blob(),url=URL.createObjectURL(blob),a=document.createElement("a");a.href=url;a.download=name;a.click();setTimeout(()=>URL.revokeObjectURL(url),3000);}catch(e){notice(e.message,true);}
}
function setPreview(block,result){
  if(previewURLs[block])URL.revokeObjectURL(previewURLs[block]);previewURLs[block]=URL.createObjectURL(new Blob([result.svg],{type:"image/svg+xml"}));
  $(`preview-${block}`).src=previewURLs[block];$(`preview-${block}`).hidden=false;$(`preview-label-${block}`).hidden=true;$(`formula-${block}`).textContent=result.formula||"Structure checked";state.monomers[block]={...result,input:$(`smiles-${block}`).value.trim()};
  refreshChemistry();
}
async function validateMonomer(block,quiet=false){
  const input=$(`smiles-${block}`).value.trim();if(!input)throw new Error(`Enter a SMILES for building block ${block.toUpperCase()}.`);
  const result=await api("monomer",{smiles:input});if($(`smiles-${block}`).value.trim()===input)setPreview(block,result);if(!quiet)notice(`Building block ${block.toUpperCase()} checked · ${result.formula}.`);return result;
}
async function refreshMonomers(){const results=await Promise.allSettled([validateMonomer("a",true),validateMonomer("b",true)]);for(const r of results)if(r.status==="rejected")notice(r.reason.message,true);}
function setCrystalMonomers(items){
  clearTimeout(liveDrawingTimer);selectionRevision++;
  $('assembly-count').value='2';fillAssembly();
  Object.assign(state,{recipe:null,result:null,sourceId:null,sourceMetadata:null,inputsDirty:true,dirty:false,jsonPending:false});
  $('framework-image').hidden=true;$('framework-image').removeAttribute('src');$('canvas-empty').classList.remove('hidden');$('exports').replaceChildren();
  $('drawing-summary').textContent='';syncJSON();drawLayout();drawLayoutTable();updateIdentity();
  if(items){
    for(const [i,b] of ['a','b'].entries()){
      const item=items[i];$('smiles-'+b).value=item?.smiles||'';$('label-'+b).value=item?.abbreviation||item?.label||'';
      delete state.monomers[b];$('preview-'+b).hidden=true;$('preview-label-'+b).hidden=false;
      $('formula-'+b).textContent=item?'Checking structure…':'No assigned precursor';
      if(item?.smiles)validateMonomer(b,true).catch(e=>notice(e.message,true));
    }
  }
  $('linkage-search').value='';populateLinkages('imine');refreshChemistry();
}
function setCombination(preset){
  if(multi()){$('assembly-count').value='2';fillAssembly();}
  for(const block of ["a","b"]){$(`smiles-${block}`).value=preset[`monomer_${block}`];$(`label-${block}`).value=preset[`label_${block}`]||block.toUpperCase();$(`preset-${block}`).value="";}
  $("linkage-search").value="";populateLinkages(preset.linkage);$("pore").value=preset.instances;markInputsChanged();const wanted=preset.template_id||templates().find(t=>t.instances===preset.instances&&(!t.linkages||t.linkages.includes(preset.linkage)))?.id;refreshChemistry(wanted);$("cof-preset").value=preset.id;updateIdentity();refreshMonomers().then(queueDrawing);
}
function showSourceMonomers(recipe,metadata){
  const layer=metadata?.periodic_design_inputs;
  let components=Object.entries(recipe.components||{}).filter(([,c])=>c.monomer_smiles);
  if(!components.length)components=(metadata?.monomers||[]).filter(m=>m.smiles||m.monomer_smiles).map((m,i)=>[m.label||`precursor-${i+1}`,{...m,monomer_smiles:m.smiles||m.monomer_smiles}]);
  const a=components.find(([key])=>key===layer?.precursor_components[0])||components.find(([,c])=>["aryl_aldehyde","tp_trialdehyde"].includes(c.role))||components[0];
  const b=components.find(([key])=>key===layer?.precursor_components[1])||components.find(([,c])=>c.role==="primary_arylamine")||components.find(x=>x!==a);
  for(const [block,entry] of [["a",a],["b",b]]){
    const smiles=entry?.[1]?.monomer_smiles||"";
    const listed=(metadata?.monomers||[]).find(item=>smiles&&(item.smiles===smiles||item.canonical_smiles===smiles));
    $(`smiles-${block}`).value=smiles;$(`label-${block}`).value=listed?.label||entry?.[1]?.label||entry?.[0]||block.toUpperCase();
    $(`preset-${block}`).value="";
  }
  for(const block of ['a','b'])if(!$(`smiles-${block}`).value){delete state.monomers[block];$(`preview-${block}`).hidden=true;$(`preview-label-${block}`).hidden=false;$(`formula-${block}`).textContent="No precursor supplied";}
  const link=layer?.linkage||recipe.linkage||a?.[1]?.linkage;if(link&&state.library?.linkages?.[link]){$("linkage-search").value="";populateLinkages(link);}if([6,12].includes(recipe.instances?.length))$("pore").value=recipe.instances.length;refreshChemistry(layer?.template||recipe.template_id);if(state.sourceId)refreshCombinations(true);else $("cof-preset").value="custom";
  for(const block of ['a','b'])if($(`smiles-${block}`).value.trim())validateMonomer(block,true).catch(e=>notice(e.message,true));
}
async function render(){
  if(!state.recipe||state.busy||state.inputsDirty||state.jsonPending)return;
  state.busy=true;sync();
  try{
    const revision=selectionRevision,recipe=state.recipe;
    const wasRandom=state.style==="random",payload=state.sourceId?{source_recipe_id:state.sourceId}:{recipe:state.recipe};Object.assign(payload,{size:Number($("image-size").value),style_id:state.style});if(state.library?.style_controls)payload.style_options=state.styleOptions;if(state.seed!==null)payload.style_seed=state.seed;
    if(wasRandom)payload.recent_styles=[...recentTreatments];
    const result=await api("render",payload);if(revision!==selectionRevision||recipe!==state.recipe)return;state.result=result;state.dirty=false;
    if(!state.sourceId){for(const key of["design_cofex_notation","design_topology","selected_template"])state.recipe[key]=clone(result.audit[key]??null);syncJSON();}
    $("framework-image").src=`.${result.downloads.svg}`;$("framework-image").hidden=false;$("canvas-empty").classList.add("hidden");
    const audit=result.audit;$("drawing-summary").textContent=`${audit.building_block_instances} blocks · ${audit.atoms} atoms · ${audit.continuation_ports.length} outward ports · ${state.sourceId?"provisional source patch":"designed finite patch"}`;
    showExports(result.downloads);
    const exportNote=$("export-note"),exportAudit=audit.chemical_exports||{};
    exportNote.hidden=Boolean(result.downloads.cdxml);
    exportNote.textContent=result.downloads.cdxml?"":exportAudit.cdxml==="unavailable_in_this_rdkit_build"?"CDXML is unavailable in this local chemistry engine. Use MOL or SDF for chemical editing.":`CDXML was omitted${exportAudit.reason?`: ${exportAudit.reason}`:" because its export could not be validated"}. MOL and SDF remain available.`;
    const chosen=audit.drawing_style;if(chosen){state.style=chosen.id||state.style;state.seed=chosen.seed??state.seed;state.resolvedStyleOptions=clone(chosen.options||{});if(wasRandom)state.styleOptions=clone(chosen.options||{});if(recentTreatments.at(-1)!==state.style){recentTreatments.push(state.style);if(recentTreatments.length>4)recentTreatments.shift();}$("canvas-sheet").classList.toggle("transparent-preview",Boolean(chosen.options?.transparent_background));setStyleDescription();drawStyleControls();}
    updateIdentity();window.dispatchEvent(new CustomEvent('cofdraw:rendered',{detail:result}));const review=state.sourceId?'':geometryReviewMessage(audit.layout_refinement);notice(state.sourceId?"Source drawing rendered with its qualified reference record. Style changes leave the graph intact.":review?`Drawing created. ${review}`:"Drawing ready. Inspect the connections and continuation ports before saving.");return result;
  }catch(error){notice(error.message,true);return null;}finally{state.busy=false;sync();}
}
$("build").addEventListener("click",async()=>{
  clearTimeout(liveDrawingTimer);
  if(window.COFdrawCrystalDesign?.active){window.COFdrawCrystalDesign.build();return;}
  if(currentPairKey()!==standardPair&&(!state.recipe||state.inputsDirty&&state.autoChemistry)&&!multi()&&!$('build-periodic-from-pair').hidden){$('build-periodic-from-pair').click();return;}
  if(state.busy)return;state.busy=true;sync();
  const revision=selectionRevision;
  try{
    if(multi()){
      const mixed=assemblyMode()==='3';
      const payload={monomer_ids:(mixed?['a','b','c']:['a']).map(k=>$('assembly-'+k).value),linker_sequence:mixed?$('assembly-pattern').value:''};
      payload.custom_precursors=payload.monomer_ids.filter(k=>customAssemblyRows.has(k)).map(k=>{const r=customAssemblyRows.get(k);return {label:r.label,smiles:r.smiles};});
      if(!mixed&&!$('self-face-field').hidden)payload.pore_view=$('self-face').value;
      if(mixed&&$('mixed-layout').value==='kgm')payload.pore_layout='kgm';
      const result=await api('prepare/assembly',payload);
      if(revision!==selectionRevision)return;
      acceptRecipe(result.recipe,null,null,true);state.busy=false;await render();return;
    }
    await Promise.all(["a","b"].map(b=>checkedMonomer(b)||validateMonomer(b,true)));const t=selectedTemplate(),check=constructionCheck();if(!check.compatible)throw new Error(check.reason);const order=engineOrder(),payload={linkage:$("linkage").value,monomer_a:$(`smiles-${order[0]}`).value.trim(),monomer_b:$(`smiles-${order[1]}`).value.trim(),label_a:$(`label-${order[0]}`).value.trim(),label_b:$(`label-${order[1]}`).value.trim(),instances:t?.instances||Number($("pore").value)};if(t)payload.template_id=t.id;const result=await api("prepare",payload);if(revision!==selectionRevision)return;acceptRecipe(result.recipe,null,null,true);$("chemistry-note").textContent=result.limits;state.busy=false;await render();}catch(error){notice(error.message,true);}finally{state.busy=false;sync();}
});
$("build-from-canvas").addEventListener("click",()=>$("build").click());
function geometryReviewMessage(refinement){return refinement?.quality==='needs_review'?'Some linkage junctions still need geometry review. Inspect the connections and use Move / adjust or Clean bond geometry before saving.':'';}
for(const block of ["a","b"]){
  let validationTimer; $(`smiles-${block}`).addEventListener('input',()=>{clearTimeout(validationTimer);validationTimer=setTimeout(()=>{if($(`smiles-${block}`).value.trim()&&!state.busy)validateMonomer(block,true).then(queueDrawing).catch(()=>{});},650);});
  $(`smiles-${block}`).addEventListener("input",()=>{$(`formula-${block}`).textContent="Preview needs update";$(`preview-${block}`).hidden=true;$(`preview-label-${block}`).hidden=false;delete state.monomers[block];$(`preset-${block}`).value="";markInputsChanged();refreshChemistry();});
  $(`label-${block}`).addEventListener("input",markInputsChanged);
  $(`preset-${block}`).addEventListener("change",()=>{const item=state.library?.monomer_presets.find(p=>p.id===$(`preset-${block}`).value);if(!item)return;const choice=monomerChoice(item,block,checkedMonomer(block==='a'?'b':'a'));if(!choice.compatible){refreshPresetLists();notice(choice.reason,true);return;}$(`smiles-${block}`).value=item.smiles;$(`label-${block}`).value=item.abbreviation||item.label?.split(" · ")[0]||item.id;delete state.monomers[block];markInputsChanged();refreshChemistry();validateMonomer(block).then(queueDrawing).catch(e=>notice(e.message,true));});
}
document.querySelectorAll(".validate-monomer").forEach(b=>b.addEventListener("click",()=>validateMonomer(b.dataset.block).catch(e=>notice(e.message,true))));
comboInput.addEventListener('focus',()=>{refreshCombinations();showSuggestions();});
comboInput.addEventListener('input',()=>{refreshCombinations();showSuggestions();});
comboInput.addEventListener('blur',closeSuggestions);
suggestions.addEventListener('pointerdown',e=>e.preventDefault());
suggestions.addEventListener('click',e=>selectSuggestion(e.target.closest('[role=option]')));
comboInput.addEventListener('keydown',e=>{
  if(e.key==='Escape'){e.preventDefault();closeSuggestions();return;}
  if(e.key==='Enter'&&!suggestions.hidden){
    e.preventDefault();selectSuggestion(suggestions.children[Math.max(0,suggestionIndex)]);return;
  }
  if(!['ArrowDown','ArrowUp'].includes(e.key))return;
  e.preventDefault();if(suggestions.hidden)showSuggestions();
  const items=[...suggestions.children];if(!items.length)return;
  suggestionIndex=suggestionIndex<0?(e.key==='ArrowDown'?0:items.length-1):
    (suggestionIndex+(e.key==='ArrowDown'?1:-1)+items.length)%items.length;
  items.forEach((item,i)=>item.setAttribute('aria-selected',String(i===suggestionIndex)));
  const active=items[suggestionIndex];comboInput.setAttribute('aria-activedescendant',active.id);active.scrollIntoView({block:'nearest'});
});
$("cof-preset").addEventListener("change",()=>{
  const selection=$('cof-preset').value;
  comboInput.value='';closeSuggestions();
  if(selection.startsWith('crystal:')){
    const id=selection.slice(8);
    window.dispatchEvent(new CustomEvent('cofdraw:open-crystal',{detail:id}));return;
  }
  window.COFdrawCrystalDesign?.leave();
  const row=combinations.find(p=>p.id===selection);
  if(row?.assembly){
    $('assembly-count').value=String(row.assembly.monomer_ids.length);
    assemblyChanged(row.assembly);$('cof-preset').value=selection;return;
  }
  if(row?.source){loadSource(row.source.recipe_id);return;}
  if(row){setCombination(row);notice("Checking the selected pair…");}
  else if(sourceMode()){markInputsChanged();refreshMonomers();}
});
$("auto-chemistry").addEventListener("change",()=>{state.autoChemistry=$("auto-chemistry").checked;if(state.autoChemistry)markInputsChanged();refreshChemistry();if(state.autoChemistry)queueDrawing();});
$("linkage").addEventListener("change",()=>{markInputsChanged();refreshChemistry();updateIdentity();});
$("linkage-search").addEventListener("input",()=>populateLinkages());
$("pore").addEventListener("change",markInputsChanged);
$("template").addEventListener("change",()=>{markInputsChanged();refreshTemplates();updateIdentity();});
async function loadSource(id){
  window.COFdrawCrystalDesign?.leave();
  if(multi()){$('assembly-count').value='2';fillAssembly();}
  clearTimeout(liveDrawingTimer);
  const revision=++selectionRevision;
  state.busy=true;sync();
  try{
    const result=await api('sources/'+encodeURIComponent(id));
    if(revision!==selectionRevision)return;
    acceptRecipe(result.recipe,id,result.metadata);
    showSourceMonomers(result.recipe,result.metadata);
    state.style=result.recipe.source_style?'source':'standard';
    state.styleOptions={};state.seed=null;state.busy=false;
    await render();
    return state.result;
  }catch(e){notice(e.message,true);}
  finally{state.busy=false;sync();}
}
function styleIcon(style){
  const span=document.createElement("span");span.className="style-icon";span.setAttribute("aria-hidden","true");const svg=node("svg",{viewBox:"0 0 56 48",fill:"none"});
  const kind=style.preview_kind||style.id,layered=style.category==="layered"||/stack|layer|3d/.test(kind),schematic=style.category==="schematic"||/cartoon|schematic/.test(kind);
  const colored=/color|red|blue|highlight|filled/.test(kind),color=colored?"#527ba1":"#647c68",points="28,7 44,16 44,34 28,43 12,34 12,16";
  if(kind==='pore_perspective'){
    for(const factor of [.50,.66,.82])svg.append(node('polygon',{points,transform:`translate(28 25) scale(${factor}) translate(-28 -25)`,stroke:'#7791b5','stroke-width':2}));
  }else if(kind==='gradient_stack'||kind==='half_space_stack'){
    for(const [depth,shade] of [[3,'#74a7a9'],[2,'#9ababa'],[1,'#b6d0c5']]){
      svg.append(node('polygon',{points,transform:`translate(${depth},${depth*2-2}) scale(.96,.88)`,stroke:shade,'stroke-width':2}));
    }
  }else if(kind==='sphere_cutaway'){
    for(const [x,y,r] of [[12,16,5],[12,34,5],[20,39,4],[20,11,4],[28,7,3]]){
      svg.append(node('circle',{cx:x,cy:y,r,fill:'#bc8fba','fill-opacity':.72}));
    }
  }else if(layered){svg.append(node("polygon",{points,transform:"translate(-4,-5)",fill:"#83a6c0","fill-opacity":.15,stroke:"#6f98b2","stroke-opacity":.4}),node("polygon",{points,transform:"translate(3,0)",fill:"#92b495","fill-opacity":.2,stroke:"#789b7b","stroke-opacity":.6}));}
  svg.append(node("polygon",{points,stroke:color,"stroke-width":schematic?2.6:1.3,fill:/fill|highlight/.test(kind)?"#c9dce9":"none","fill-opacity":.6}));
  if(!schematic){svg.append(node("path",{d:"M17 18v13 M28 37l11-6 M28 13l11 6",stroke:color,"stroke-width":1}));if(colored)svg.append(node("path",{d:"M28 7 44 16v18",stroke:"#bd7461","stroke-width":2}));
    for(const[x,y,txt]of[[9,18,"N"],[43,37,"N"]]){svg.append(node("rect",{x:x-2,y:y-7,width:9,height:10,fill:"white"}),node("text",{x,y,fill:color,"font-size":9,"font-family":/serif/.test(kind)?"Georgia,serif":"Arial,sans-serif"},txt));}
  }else for(const[x,y]of[[28,7],[44,34],[12,34]])svg.append(node("circle",{cx:x,cy:y,r:4,fill:x===28?"#bf8069":"#658caa",stroke:"white","stroke-width":1}));
  if(/wave/.test(kind))svg.append(node("path",{d:"M27 2q-3 3 0 5t0 5",stroke:"#719178","stroke-width":1}));
  span.append(svg);return span;
}
function drawStyles(){
  const crystal=window.COFdrawCrystalDesign;
  const styles=(crystal?.active?crystal.styles:state.library?.styles)||[{id:"standard",label:"Standard",description:"Conventional chemical bonds"}],categories=["all",...new Set(styles.map(s=>s.category||"chemical"))],names={all:"All",chemical:"Chemical",layered:"Layered / depth",schematic:"Schematic"};
  $("style-categories").replaceChildren();for(const category of categories){const b=document.createElement("button");b.textContent=names[category]||category;b.setAttribute("aria-pressed",String(state.styleCategory===category));b.addEventListener("click",()=>{state.styleCategory=category;drawStyles();sync();});$("style-categories").append(b);}
  $("style-list").replaceChildren();for(const style of L.stylesInCategory(styles,state.styleCategory)){
    const button=document.createElement("button");button.className="style-choice";button.dataset.style=style.id;const label=document.createElement("span"),title=document.createElement("strong"),check=document.createElement("span");title.textContent=style.label;button.title=style.description||style.label;label.append(title);
    check.className="style-check";check.textContent="✓";button.append(styleIcon(style),label,check);button.addEventListener("click",()=>chooseStyle(style.id));$("style-list").append(button);
  }
}
const styleOptionDefaults={hydrogen_labels:"show",element_scheme:"style_default",palette:"classic",font_family:"Arial",font_scale:1,bond_scale:1,layer_count:3,layer_opacity:.24,fill_opacity:.18,edge_style:"wavy_cross",transparent_background:false,abbreviation_mode:"explicit",highlight_width:1,glow_strength:.55,layer_spacing:.04,projection_tilt:35,sphere_scale:1};
function drawStyleControls(){
  const spec=state.library?.style_controls,root=$("style-control-fields");$("style-controls").hidden=!spec;if(!spec)return;root.replaceChildren();
  const styleDefaults=state.style==='element_compact'?{element_scheme:'jewel',hydrogen_labels:'compact'}:
    state.style==='pore_perspective'?{layer_count:5,layer_spacing:.04,layer_opacity:.64,projection_tilt:12,palette:'cobalt_citrus',transparent_background:true}:
    state.style==='half_space_stack'||state.style==='gradient_stack'?{layer_count:5,layer_spacing:.044,layer_opacity:.56,projection_tilt:43,palette:'lagoon'}:{};
  const current={...styleOptionDefaults,...styleDefaults,...state.resolvedStyleOptions,...state.styleOptions};
  const hasLayers=['stacked_eclipsed','stacked_slipped','gradient_stack','half_space_stack','pore_perspective'].includes(state.style);
  const hasSpheres=['molecular_spheres','sphere_cutaway','half_space_stack','pore_perspective'].includes(state.style);
  const lists=[['palette','Color palette',spec.palettes||[]],['element_scheme','Atom colors',spec.element_schemes||['style_default']],['hydrogen_labels','Functional-group H',spec.hydrogen_labels||['show','compact']],['font_family','Atom-label font',spec.fonts||[]],['edge_style','Continuation marks',spec.edge_styles||[]],['abbreviation_mode','Substituent labels',spec.abbreviation_modes||['explicit','condensed','OMe']]];
  const names={show:'Show H labels',compact:'Hide inferable H labels',style_default:"Style default",conventional:"Conventional elements",deep_ink:"Deep ink",jewel:"Jewel tones",classic:"Red / blue",teal_orange:"Teal / orange",purple_gold:"Purple / gold",pastel:"Pastel",ink:"Black",wavy_cross:"Transverse wave",dashed:"Dashed",crossbar:"Crossbar",bond_stub:"Plain bond end",explicit:"Explicit atoms",condensed:"OCH3 / SO3H",OMe:"OMe / OEt / SO3H"};
  const crystal=window.COFdrawCrystalDesign?.active;
  const usesBlockColors=!crystal||['depth_ink','reticular_depth','ring_fill'].includes(state.style)||['reticular_blocks','red_blue','bold_color'].includes(state.style)&&!window.COFdrawCrystalDesign.isReference;
  function changed(key,value){state.styleOptions[key]=value;if(state.recipe&&!crystal)state.dirty=true;sync();}
  for(const[key,label,values]of lists){if(crystal&&(['edge_style','abbreviation_mode'].includes(key)||(key==='palette'&&!usesBlockColors)))continue;const field=document.createElement("label");field.className="field";const span=document.createElement("span");span.textContent=label;const select=document.createElement("select");select.id=`style-${key}`;select.dataset.styleOption=key;for(const value of values)select.append(option(value,names[value]||value));select.value=current[key];select.addEventListener("change",()=>changed(key,select.value));field.append(span,select);root.append(field);}
  const swatches=document.createElement('div');swatches.className='palette-swatches';
  for(const [name,colors] of Object.entries(spec.palette_colors||{})){const button=document.createElement('button');button.className='palette-swatch'+(name===current.palette?' selected':'');button.title=name.replaceAll('_',' ');button.setAttribute('aria-label',button.title+' palette');for(const color of colors){const i=document.createElement('i');i.style.background=color;button.append(i);}button.addEventListener('click',()=>{delete state.styleOptions.custom_palette;changed('palette',name);if(crystal){drawStyleControls();window.COFdrawCrystalDesign.redraw();}else if(state.recipe&&!state.inputsDirty&&!state.jsonPending)render();else drawStyleControls();});swatches.append(button);}if(usesBlockColors)root.prepend(swatches);
  const custom=document.createElement('div');custom.className='custom-colors';const title=document.createElement('span');title.textContent='Your colors';custom.append(title);const selected=current.custom_palette||spec.palette_colors?.[current.palette]||['#167d9a','#ee796b','#6796b5','#cfb86d'];
  for(let i=0;i<(crystal?2:4);i++){const input=document.createElement('input');input.type='color';input.value=selected[i%selected.length];input.setAttribute('aria-label',`Custom building-block color ${i+1}`);input.addEventListener('input',()=>changed('custom_palette',[...custom.querySelectorAll('input')].map(x=>x.value)));custom.append(input);}if(usesBlockColors)root.append(custom);
  for(const[key,label,min,max,step]of[['font_scale','Label size',.6,1.6,.05],['bond_scale','Bond weight',.5,2.5,.1],['fill_opacity','Fill opacity',.05,.6,.05],['highlight_width','Highlight width',.5,2.5,.1],...(state.style==='glow_outline'?[['glow_strength','Glow softness',0,1,.05]]:[]),...(hasSpheres?[['sphere_scale','Sphere size',.5,1.6,.05]]:[]),...(hasLayers?[['layer_count','Illustrative layers',2,5,1],['layer_opacity','Back-layer opacity',.05,.7,.05],['layer_spacing','Layer separation',.015,.09,.005]]:[]),...(hasLayers||state.style==='perspective'?[['projection_tilt','Perspective tilt',0,65,1]]:[])]){
    if(crystal&&(key==='highlight_width'||key==='fill_opacity'&&!['ring_fill','reticular_blocks','reticular_depth'].includes(state.style)))continue;
    const field=document.createElement("label");field.className="field";const span=document.createElement("span"),output=document.createElement("output"),input=document.createElement("input");input.id=`style-${key}`;input.dataset.styleOption=key;input.type="range";input.min=min;input.max=max;input.step=step;input.value=current[key];span.textContent=label;output.textContent=String(current[key]);span.append(output);input.addEventListener("input",()=>{output.textContent=input.value;changed(key,Number(input.value));});field.append(span,input);root.append(field);
  }
  const label=document.createElement("label");label.className="checkbox-field";const check=document.createElement("input");check.type="checkbox";check.id="style-transparent_background";check.checked=Boolean(current.transparent_background);check.addEventListener("change",()=>changed("transparent_background",check.checked));label.append(check,document.createTextNode("Transparent image background"));root.append(label);
}
function setStyleDescription(){
  const crystal=window.COFdrawCrystalDesign,item=(crystal?.active?crystal.styles:state.library?.styles)?.find(s=>s.id===state.style);$("style-description").textContent=[item?.description||"A randomly selected visual style.",item?.limits,state.seed!==null?`Style seed: ${state.seed}.`:null].filter(Boolean).join(" ");
  const special=item&&item.category&&item.category!=="chemical",note=$("representation-note");note.hidden=!special;
  note.textContent=special?(item.representation_scope||item.limits||(item.category==="layered"?"Layered depiction: illustrative offset copies of one finite chemical patch, not a crystal packing prediction.":"Schematic depiction: simplified building blocks. Chemical graph exports retain the explicit atom/bond record.")):"";
}
function chooseStyle(id,seed=null){if(window.COFdrawCrystalDesign?.active){window.COFdrawCrystalDesign.chooseStyle(id);return;}state.style=id;state.seed=seed;state.styleOptions={};state.resolvedStyleOptions={};setStyleDescription();drawStyleControls();if(state.recipe){state.dirty=true;sync();render();}else sync();}
$("apply-style-options").addEventListener("click",()=>window.COFdrawCrystalDesign?.active?window.COFdrawCrystalDesign.redraw():render());
$("reset-style-options").addEventListener("click",()=>chooseStyle(state.style));
$("random-style").addEventListener("click",()=>chooseStyle("random",crypto.getRandomValues(new Uint32Array(1))[0] % 2147483647));
$("image-size").addEventListener("change",()=>{if(state.recipe&&!window.COFdrawCrystalDesign?.active){state.dirty=true;sync();render();}});
$("render-changes").addEventListener("click",render);
function setZoom(value){state.zoom=Math.min(3,Math.max(.5,value));const sheet=$("canvas-sheet");sheet.classList.toggle("zoomed",state.zoom!==1);sheet.style.width=state.zoom===1?"100%":`${state.zoom*100}%`;sheet.style.height=state.zoom===1?"100%":`${state.zoom*100}%`;$("zoom-label").textContent=state.zoom===1?"100%":`${Math.round(state.zoom*100)}%`;}
$("zoom-in").addEventListener("click",()=>setZoom(state.zoom+.25));$("zoom-out").addEventListener("click",()=>setZoom(state.zoom-.25));$("zoom-fit").addEventListener("click",()=>setZoom(1));
function node(name,attrs,text){const e=document.createElementNS(NS,name);Object.entries(attrs).forEach(([k,v])=>e.setAttribute(k,v));if(text!==undefined)e.textContent=text;return e;}
function layoutPoint(xy){return[320+(xy[0]-state.layoutCenter[0])*state.layoutScale,150-(xy[1]-state.layoutCenter[1])*state.layoutScale];}
function blockLayoutPosition(item){
  const rows=state.recipe?.atom_coordinates?.filter(row=>row.instance===item.id);
  if(!rows?.length)return item.layout_position;
  return [0,1].map(axis=>rows.reduce((sum,row)=>sum+row.xy[axis],0)/rows.length);
}
function moveLayoutBlock(item,position){
  const previous=blockLayoutPosition(item),delta=position.map((value,axis)=>value-previous[axis]);
  if(Array.isArray(state.recipe.atom_coordinates)){
    for(const row of state.recipe.atom_coordinates)if(row.instance===item.id)row.xy=row.xy.map((value,axis)=>value+delta[axis]);
  }
  item.layout_position=[...position];
}
function drawLayout(recalculate=true){
  const canvas=$("layout-canvas");canvas.replaceChildren();if(!state.recipe||state.recipe.explicit_graph)return;const items=state.recipe.instances;
  if(recalculate){const positions=items.map(blockLayoutPosition),xs=positions.map(p=>p[0]),ys=positions.map(p=>p[1]);state.layoutCenter=[(Math.min(...xs)+Math.max(...xs))/2,(Math.min(...ys)+Math.max(...ys))/2];state.layoutScale=Math.min(500/Math.max(1,Math.max(...xs)-Math.min(...xs)),210/Math.max(1,Math.max(...ys)-Math.min(...ys)));}
  const map=new Map(items.map(i=>[i.id,i]));for(const edge of state.recipe.connections){const a=map.get(edge.a?.instance),b=map.get(edge.b?.instance);if(!a||!b)continue;const[x1,y1]=layoutPoint(blockLayoutPosition(a)),[x2,y2]=layoutPoint(blockLayoutPosition(b));canvas.append(node("line",{x1,y1,x2,y2,stroke:"#b4c5ad","stroke-width":1.5}));}
  items.forEach((item,i)=>{const[x,y]=layoutPoint(blockLayoutPosition(item)),g=node("g",{class:"block-node","data-index":i});g.append(node("circle",{cx:x,cy:y,r:17,fill:i%2?"#f6edda":"#e3eee0",stroke:i%2?"#b19b69":"#759d79"}),node("text",{x,y:y+4,"text-anchor":"middle","font-size":11,"pointer-events":"none",fill:"#355844"},item.id));canvas.append(g);});
}
function changedLayout(rebuildCoordinates=false){
  if(rebuildCoordinates&&Array.isArray(state.recipe.atom_coordinates)){
    // Attachment or reflection edits need a new assembly. Retain the latest
    // block positions before discarding coordinates whose atom keys may change.
    for(const item of state.recipe.instances)item.layout_position=blockLayoutPosition(item);
    delete state.recipe.atom_coordinates;
  }
  removeIdentity();state.recipe.layout_refinement='manual';state.dirty=true;syncJSON();sync();window.dispatchEvent(new CustomEvent('cofdraw:layoutchanged'));
}
function drawLayoutTable(){
  const body=$("layout-rows");body.replaceChildren();if(!state.recipe||state.recipe.explicit_graph)return;
  state.recipe.instances.forEach(item=>{const row=document.createElement("tr"),label=document.createElement("td");label.textContent=item.id;row.append(label);["x","y","incoming_port","outgoing_port","flip_across_attachment_axis"].forEach((key,axis)=>{const cell=document.createElement("td"),input=document.createElement("input");input.setAttribute("aria-label",`${item.id} ${key}`);input.type=axis===4?"checkbox":"number";if(axis===4)input.checked=Boolean(item[key]);else{input.value=axis<2?Number(blockLayoutPosition(item)[axis].toFixed(3)):item[key];input.step=axis<2?"any":"1";input.min=axis<2?"-200":"1";input.max=axis<2?"200":"100";}
    input.addEventListener("change",()=>{if(state.busy)return;const value=axis===4?input.checked:Number(input.value);if(axis!==4&&(!input.checkValidity()||!Number.isFinite(value))){notice("Enter a valid coordinate or attachment port.",true);drawLayoutTable();return;}if(axis<2){const position=[...blockLayoutPosition(item)];position[axis]=value;moveLayoutBlock(item,position);}else if(axis===4)item[key]=value;else{const previous=item[key];for(const edge of state.recipe.connections)for(const end of[edge.a,edge.b])if(end.instance===item.id&&end.port===previous)end.port=value;item[key]=value;}changedLayout(axis>=2);drawLayout();});cell.append(input);row.append(cell);});body.append(row);});
}
let dragged=null;
$("layout-canvas").addEventListener("pointerdown",e=>{const target=e.target.closest(".block-node");if(!target||state.busy||state.jsonPending)return;dragged=Number(target.dataset.index);e.currentTarget.setPointerCapture(e.pointerId);e.preventDefault();});
$("layout-canvas").addEventListener("pointermove",e=>{if(dragged===null)return;const p=new DOMPoint(e.clientX,e.clientY).matrixTransform(e.currentTarget.getScreenCTM().inverse());const position=[(p.x-320)/state.layoutScale+state.layoutCenter[0],(150-p.y)/state.layoutScale+state.layoutCenter[1]].map(x=>Math.max(-200,Math.min(200,Number(x.toFixed(3)))));moveLayoutBlock(state.recipe.instances[dragged],position);changedLayout();drawLayout(false);});
for(const type of["pointerup","pointercancel"])$("layout-canvas").addEventListener(type,()=>{if(dragged!==null){dragged=null;syncJSON();drawLayoutTable();}});
$("recipe-json").addEventListener("input",()=>{state.jsonPending=true;state.dirty=true;sync();});
$("apply-recipe").addEventListener("click",()=>{try{const input=$("recipe-json").value;if(input.length>200000)throw new Error("Recipe limit is 200 kB.");const value=JSON.parse(input);acceptRecipe(value.recipe||value);notice("Recipe applied as a designed variant. Render to check its chemistry and connections.");}catch(e){notice(e.message,true);}});
$("save-project").addEventListener("click",()=>{if(window.COFdrawCrystalDesign?.active){window.COFdrawCrystalDesign.save();return;}const project={format:"cofdraw-project",version:1,chemistry_mode:state.autoChemistry?"auto":"manual",saved_at:new Date().toISOString(),ground_truth:false,reaction_scheme:window.COFdrawScheme?.snapshot(),drawing_extent:window.COFdrawExtent?.value||1,color_pore_types:$('color-pore-types').checked,monomers:assemblyFeeds(state.recipe)||["a","b"].map(block=>({label:$(`label-${block}`).value,smiles:$(`smiles-${block}`).value,editor_graph:state.monomers[block]?.graph||null})),linkage:$("linkage").value,instances:selectedTemplate()?.instances||Number($("pore").value),template_id:selectedTemplate()?.id||null,recipe:state.recipe,style_id:state.style,style_seed:state.seed,style_options:clone(state.styleOptions),image_size:Number($("image-size").value),derived_from_source_recipe_id:state.sourceId,catalog_version:state.library?.catalog_version||null,inputs_changed_since_recipe:state.inputsDirty};saveText("COFdraw-project.json",JSON.stringify(project,null,2)).catch(e=>notice(e.message,true));});
$("open-project").addEventListener("click",()=>$("project-file").click());
$("project-file").addEventListener("change",async e=>{try{const file=e.target.files[0];if(!file)return;if(file.size>1000000)throw new Error("Project files must be smaller than 1 MB.");const value=JSON.parse(await file.text());window.COFdrawScheme?.restore(value.reaction_scheme);if(value.format==="cofdraw-crystal-project"){await window.COFdrawCrystalDesign.restore(value);return;}window.COFdrawCrystalDesign?.leave();state.autoChemistry=value.chemistry_mode==="auto";
  const imported=value.recipe||value;
  for(const input of (imported.assembly?.custom_precursors||[]).slice(0,3)){
    const row=await api("assembly/precursor",{label:input.label,smiles:input.smiles});customAssemblyRows.set(row.id,row);
  }
  acceptRecipe(imported);if(value.format==="cofdraw-project"&&Array.isArray(value.monomers)){for(const[block,entry]of[["a",value.monomers[0]],["b",value.monomers[1]]])if(entry){$(`label-${block}`).value=String(entry.label||block.toUpperCase()).slice(0,100);$(`smiles-${block}`).value=String(entry.smiles||"").slice(0,4000);}$("linkage-search").value="";populateLinkages(state.library?.linkages?.[value.linkage]?value.linkage:"imine");$("pore").value=[6,12].includes(value.instances)?value.instances:6;refreshChemistry(value.template_id||value.recipe?.template_id);state.style=state.library?.styles.some(s=>s.id===value.style_id)?value.style_id:"standard";state.seed=Number.isInteger(value.style_seed)?value.style_seed:null;state.styleOptions=value.style_options&&typeof value.style_options==="object"&&!Array.isArray(value.style_options)?clone(value.style_options):{};state.resolvedStyleOptions={};$("image-size").value=[1000,1600,2400].includes(value.image_size)?value.image_size:1600;state.inputsDirty=Boolean(value.inputs_changed_since_recipe);}else showSourceMonomers(state.recipe);window.COFdrawExtent?.restore(value.drawing_extent,value.color_pore_types);$("cof-preset").value="custom";setStyleDescription();drawStyleControls();if(!multi())await refreshMonomers();sync();if(!state.inputsDirty)await render();notice("Project opened as an editable design. Source records are reloaded separately from the audited library.");}catch(error){notice(error.message,true);}finally{e.target.value="";}});
$("copy-identity").addEventListener("click",async()=>{
  if(!state.recipe||state.inputsDirty||state.jsonPending)return;
  const source=Boolean(state.sourceId),metadata=state.sourceMetadata||state.result?.audit?.literature_metadata||{};
  const notation=source?(metadata.cofex_notation||metadata.cofex_string||null):state.recipe.design_cofex_notation||state.result?.audit?.drawing_cofex?.text||null;
  const record={cofex_id:source?(metadata.cofex_id||state.recipe.cofex_id||null):null,
    cofex_notation:source?notation:null,source_recipe_id:state.sourceId,
    design_cofex_notation:source?null:notation,design_topology:source?null:state.recipe.design_topology||null,
    linkage:state.recipe.linkage||null,selected_template:state.recipe.selected_template||null,
    catenation:notation?.match(/\sCOFex\.[^.\s]+\.(cat[^;\s]+);/)?.[1]||null,
    catenation_qualification:$("identity-cat").textContent,qualification:$("identity-note").textContent,
    drawing_cofex:state.result?.audit?.drawing_cofex?.text||null,
    ground_truth:false,graph_scope:"finite_patch_with_continuation_ports",
    monomers:assemblyFeeds(state.recipe)||["a","b"].map(b=>({label:$(`label-${b}`).value,smiles:$(`smiles-${b}`).value}))};
  try{await navigator.clipboard.writeText(JSON.stringify(record,null,2));notice("Full structure record copied.");}
  catch{await saveText("COFdraw-structure-record.json",JSON.stringify(record,null,2));}
});
$("help-button").addEventListener("click",()=>$("help-dialog").showModal());

// The partner view starts with compatible candidates; show-all keeps inspection available.
const libraryView={selected:null,observer:null,busy:false};
const libraryPreviews=L.createPreviewCache(async smiles=>{
  const result=await api("monomer",{smiles});
  if(!result.valid||!result.svg)throw new Error("A chemical preview is unavailable for this entry.");
  return {result,url:URL.createObjectURL(new Blob([result.svg],{type:"image/svg+xml"}))};
},3);
function libraryItem(){return state.library?.monomer_presets?.find(item=>item.id===libraryView.selected);}
function libraryChoice(item){const target=$("library-target").value;return monomerChoice(item,target,target==="editor"?null:checkedMonomer(target==="a"?"b":"a"));}
function libraryGroups(item){const groups=[...L.groupRows(item),...(item.catalog_reactive_groups||[])],seen=new Set();return groups.filter(g=>{const key=`${g.kind}:${g.count}`;if(seen.has(key))return false;seen.add(key);return true;});}
function libraryReferenceOnly(item){return item.catalog_use==="reference_only"||item.automatic_prepare_supported===false||item.annotation_only===true;}
function libraryCurrent(block){const smiles=$(`smiles-${block}`).value.trim(),item=state.library?.monomer_presets?.find(row=>row.smiles===smiles||row.canonical_smiles===smiles);return {smiles,item,name:smiles?(item?.abbreviation||$(`label-${block}`).value.trim()||"Custom monomer"):`Choose ${block.toUpperCase()}`};}
function updateLibraryPair(target){
  for(const block of ["a","b"]){const current=libraryCurrent(block),button=$(`library-step-${block}`),active=target===block;button.setAttribute("aria-pressed",String(active));button.disabled=libraryView.busy;button.classList.toggle("has-value",Boolean(current.smiles));$(`library-current-${block}`).textContent=current.name;button.title=current.item?L.displayName(current.item):current.name;button.setAttribute("aria-label",`Choose building block ${block.toUpperCase()}. Current: ${current.name}`);button.querySelector(".library-step-action").textContent=active?"Choosing":current.smiles?"Change":"Choose";}
  $("library-editor-mode").setAttribute("aria-pressed",String(target==="editor"));$("library-editor-mode").disabled=libraryView.busy;$("monomer-library-dialog").dataset.target=target;
}
function drawLibraryEvidence(item){
  const evidence=Array.isArray(item.source_evidence)?item.source_evidence:[],groundTruth=item.ground_truth===true,explicitProvisional=item.ground_truth===false||evidence.some(row=>row.ground_truth===false);
  $("library-audit-qualifier").textContent=groundTruth?"This catalog record is explicitly marked as ground truth; its source records and scope are shown below.":explicitProvisional?"Source-catalog candidate — not ground truth. Review evidence does not establish synthetic feasibility.":"No monomer-level ground-truth status is supplied. Keep the source record's scope and qualifiers.";
  $("library-name-scope").textContent=item.name_scope||"Abbreviations can refer to different molecules in different papers.";
  const records=$("library-evidence-records");records.replaceChildren();
  for(const row of evidence){const entry=document.createElement("li");entry.textContent=[row.annotation_id||row.source_label||"Source record",row.reported_framework,row.source_locator,row.source_page!==undefined?`page ${row.source_page}`:null,row.review_status?row.review_status.replaceAll("_"," "):null,row.ground_truth===false?"provisional":row.ground_truth===true?"marked ground truth":null].filter(Boolean).join(" · ");records.append(entry);}
  if(!evidence.length){const entry=document.createElement("li");entry.textContent=item.review_status?`Catalog review: ${item.review_status.replaceAll("_"," ")}`:"No per-source review record is supplied for this entry.";records.append(entry);}
  if(item.master_binding&&typeof item.master_binding==="object"){const binding=item.master_binding,entry=document.createElement("li"),details=document.createElement("details"),summary=document.createElement("summary"),qualification=document.createElement("p");summary.textContent=binding.status==="source_doi_and_finite_precursor_graph_match"?"Master-key match: source DOI and finite precursor graph":"Master-key audit record";qualification.textContent=[binding.qualification,binding.ground_truth===false?"This match is not a human ground-truth admission.":null].filter(Boolean).join(" ");details.append(summary,qualification);for(const key of ["record_ids","component_ids","reference_index_sha256","master_sha256","mapped_sha256","installation_receipt_sha256"]){const value=binding[key];if(!value)continue;const line=document.createElement("p");line.className="library-audit-record";line.textContent=`${key.replaceAll("_"," ")}: ${Array.isArray(value)?value.join(", "):String(value)}`;details.append(line);}entry.append(details);records.append(entry);}
  const links=$("library-source-links");links.replaceChildren();const dois=[...new Set([...(item.source_dois||[]),...evidence.map(row=>row.source_doi)].filter(value=>typeof value==="string"&&/^10\.\d{4,9}\/\S+$/i.test(value)))];
  for(const doi of dois){const link=document.createElement("a");link.href=`https://doi.org/${doi}`;link.target="_blank";link.rel="noopener noreferrer";link.textContent=doi;links.append(link);}
  $("library-limitations").textContent=[...(Array.isArray(item.limitations)?item.limitations:[]),item.compatibility_scope,item.qualification].filter(value=>typeof value==="string").join(" ");
}
async function loadLibraryThumbnail(card,item){
  if(card.dataset.previewState)return;
  card.dataset.previewState="loading";
  const img=card.querySelector(".library-thumbnail img"),label=card.querySelector(".thumbnail-state");
  label.textContent="Preparing structure…";
  try{
    const preview=await libraryPreviews.request(item.smiles);
    if(!card.isConnected)return;
    img.src=preview.url;await img.decode();
    img.hidden=false;label.hidden=true;card.dataset.previewState="ready";
  }catch(error){if(card.isConnected){card.dataset.previewState="error";label.textContent="Preview unavailable · open details";label.title=error.message;}}
}
async function selectLibraryItem(id,expand=true){
  libraryView.selected=id;const item=libraryItem();if(!item)return;document.querySelectorAll(".library-row").forEach(row=>row.classList.toggle("selected",row.dataset.id===id));
  $("library-detail-panel").hidden=false;$("library-detail-label").textContent=item.abbreviation||item.label||item.id;if(expand)$("library-detail-panel").open=true;
  $("library-name").textContent=item.abbreviation||item.label||item.id;$("library-full-name").textContent=item.full_name||L.displayName(item);
  $("library-site-count").textContent=L.connectionCount(item)===null?"Not cataloged":`${L.connectionCount(item)} reactive connection sites`;
  const symmetry=item.point_group||item.symmetry;$("library-geometry").textContent=[item.geometry_hint?.replaceAll("_"," "),symmetry?`Point group: ${Array.isArray(symmetry)?symmetry.join(", "):symmetry}`:null].filter(Boolean).join("; ")||"Not assigned";
  const chips=$("library-groups");chips.replaceChildren();for(const group of libraryGroups(item)){const chip=document.createElement("span");chip.className="group-chip";chip.textContent=`${group.count} × ${group.label||group.kind.replaceAll("_"," ")}`;chips.append(chip);}
  const choice=libraryChoice(item),target=$("library-target").value;$("library-compatibility").textContent=choice.reason;$("library-compatibility").classList.toggle("unavailable",!choice.compatible);$("library-apply").disabled=!choice.compatible||libraryView.busy;$("library-edit").disabled=libraryView.busy;$("library-apply").textContent=target==="editor"?"Open single monomer":`Use as ${target.toUpperCase()}`;
  $("library-solubility").disabled=libraryView.busy;
  drawLibraryEvidence(item);if(expand){$("library-evidence-panel").open=false;requestAnimationFrame(()=>{$("library-detail-panel").scrollIntoView({block:"nearest"});$("library-detail-panel").querySelector("summary").focus({preventScroll:true});});}
  if(!$("library-detail-panel").open)return;
  $("library-preview").hidden=true;$("library-preview-state").hidden=false;$("library-preview-state").textContent="Preparing local structure preview…";
  try{const preview=await libraryPreviews.request(item.smiles);if(libraryView.selected!==id)return;$("library-preview").src=preview.url;$("library-preview").hidden=false;$("library-preview-state").hidden=true;}catch(e){if(libraryView.selected===id)$("library-preview-state").textContent=e.message;}
}
function drawLibrary(){
  const target=$("library-target").value,partnerView=target==='b'&&Boolean($("smiles-a").value.trim()),showAll=$("library-show-all").checked;
  const visibleCatalog=(state.library?.monomer_presets||[]).filter(r=>r.catalog_use!=="synthetic_intermediate");
  const searched=L.searchMonomers(visibleCatalog,$("library-search").value,$("library-sites").value),rows=partnerView&&!showAll?searched.filter(item=>libraryChoice(item).compatible):searched;
  $("library-show-all-field").hidden=!partnerView;$("library-target").disabled=libraryView.busy;
  updateLibraryPair(target);
  $("library-context").textContent=target==="editor"?"Browse any monomer to edit or export on its own.":partnerView?(!checkedMonomer('a')?"Check A's structure before choosing a partner.":showAll?`All monomers shown. Unavailable partners for ${libraryCurrent('a').name} are disabled.`:`Choose a compatible B for ${libraryCurrent('a').name}. ${state.autoChemistry?'Chemistry matches automatically.':`Linkage: ${linkageInfo().label||$("linkage").value}.`}`):target==="b"?"Choose A first to narrow this list to compatible partners.":"Choose A, then select a compatible partner B.";
  $("library-count").textContent=`${rows.length} ${partnerView&&!showAll?'compatible ':''}of ${visibleCatalog.length} monomers`;
  const list=$("library-list");libraryView.observer?.disconnect();list.replaceChildren();list.scrollTop=0;
  const byId=new Map(rows.map(item=>[item.id,item]));
  libraryView.observer="IntersectionObserver" in window?new IntersectionObserver(entries=>{
    for(const entry of entries)if(entry.isIntersecting){libraryView.observer?.unobserve(entry.target);loadLibraryThumbnail(entry.target,byId.get(entry.target.dataset.id));}
  },{root:list,rootMargin:"120px 0px"}):null;
  for(const item of rows){
    const choice=libraryChoice(item),reference=libraryReferenceOnly(item),row=document.createElement("article");row.className=`library-row library-structure-card${choice.compatible?"":" unavailable"}`;row.dataset.id=item.id;row.dataset.availability=reference?"reference":choice.compatible?"available":"incompatible";
    const heading=document.createElement("div");heading.className="library-card-heading";
    const title=document.createElement("strong"),count=document.createElement("span");title.textContent=item.abbreviation||item.label||item.id;title.title=L.displayName(item);count.className="library-site-badge";count.textContent=L.connectionCount(item)===null?"?c":`${L.connectionCount(item)}c`;heading.append(title,count);
    const inspect=document.createElement("button");inspect.className="inspect-monomer library-thumbnail";inspect.type="button";inspect.setAttribute("aria-label",`Inspect ${L.displayName(item)}`);
    const img=document.createElement("img"),pending=document.createElement("span");img.alt=`${item.abbreviation||item.id} chemical structure`;img.hidden=true;img.decoding="async";pending.className="thumbnail-state";pending.textContent="Structure loads when visible";inspect.append(img,pending);inspect.addEventListener("click",()=>selectLibraryItem(item.id));
    const name=document.createElement("p"),status=document.createElement("p");name.className="library-card-name";name.textContent=item.full_name||item.label||item.id;name.title=name.textContent;status.className="library-card-status";status.id=`library-choice-${item.id}`;status.textContent=reference?"Reference only · inspect or edit":choice.compatible?(target==="editor"?"Available in editor":partnerView?"Compatible partner":"Available for drawing"):choice.reason;status.title=choice.reason;
    const actions=document.createElement("div"),details=document.createElement("button"),apply=document.createElement("button");actions.className="library-card-actions";details.type="button";details.className="library-card-details";details.textContent="Details";details.addEventListener("click",()=>selectLibraryItem(item.id));apply.className="library-use";apply.type="button";apply.textContent=target==="editor"?"Open editor":`Use as ${target.toUpperCase()}`;apply.disabled=!choice.compatible||libraryView.busy;apply.title=choice.reason;apply.setAttribute("aria-label",`Use ${item.abbreviation||item.id} ${target==="editor"?"in editor":`as ${target.toUpperCase()}`}`);apply.addEventListener("click",()=>applyLibraryItem(item));actions.append(details,apply);
    apply.setAttribute("aria-describedby",status.id);const sa=document.createElement('p');sa.className='library-sa';sa.dataset.saScore='true';const score=item.synthetic_accessibility;sa.textContent=score?.available?`SA Score ${score.score.toFixed(2)} / 10`:'SA Score unavailable';sa.title='1 easier–10 harder; finite-molecule heuristic, not a synthesis route or COF yield.';
    row.append(heading,inspect,name,status,sa,actions);list.append(row);if(libraryView.observer)libraryView.observer.observe(row);else loadLibraryThumbnail(row,item);
  }
  if(!rows.length){const p=document.createElement("p");p.className="style-empty";p.textContent=partnerView&&!showAll?"No compatible partners match these filters. Clear the search or site filter, check A, or show all monomers to inspect why a choice is unavailable.":"No monomers match these filters. Clear the name or reactive-site filter to browse all entries.";list.append(p);$("library-detail-panel").hidden=true;$("library-apply").disabled=true;$("library-edit").disabled=true;return;}
  const selected=rows.find(item=>item.id===libraryView.selected);if(selected)selectLibraryItem(selected.id,false);else{$("library-detail-panel").hidden=true;$("library-detail-panel").open=false;libraryView.selected=null;}
}
function openLibrary(target="a"){if(!state.online)return;$("library-target").value=target;$("library-search").value="";$("library-sites").value="all";$("library-show-all").checked=false;libraryView.selected=null;$("library-detail-panel").open=false;$("monomer-library-dialog").showModal();drawLibrary();}
function setLibraryTarget(target){if(libraryView.busy)return;$("library-target").value=target;$("library-show-all").checked=false;drawLibrary();}
async function applyLibraryItem(item=libraryItem()){
  if(!item||libraryView.busy||!libraryChoice(item).compatible)return;const target=$("library-target").value;
  if(target==="editor"){$("monomer-library-dialog").close();openSketch("a",item);return;}
  libraryView.busy=true;$("library-apply").disabled=true;$("library-edit").disabled=true;
  try{$(`smiles-${target}`).value=item.smiles;$(`label-${target}`).value=item.abbreviation||item.label||item.id;delete state.monomers[target];markInputsChanged();await validateMonomer(target,true);refreshChemistry();queueDrawing();notice(`Building block ${target.toUpperCase()} loaded from the monomer library.`);
    if(target==="a"){$("library-target").value="b";$("library-search").value="";$("library-sites").value="all";$("library-show-all").checked=false;libraryView.selected=null;}else $("monomer-library-dialog").close();
  }catch(e){$("library-compatibility").textContent=e.message;$("library-compatibility").classList.add("unavailable");}finally{libraryView.busy=false;drawLibrary();}
}
$("open-monomer-library").addEventListener("click",()=>openLibrary("a"));document.querySelectorAll(".library-open").forEach(b=>b.addEventListener("click",()=>openLibrary(b.dataset.block)));
$("monomer-library-dialog").addEventListener("close",()=>libraryView.observer?.disconnect());
$("library-detail-panel").addEventListener("toggle",()=>{if($("library-detail-panel").open&&libraryView.selected)selectLibraryItem(libraryView.selected,false);});
$("library-solubility").addEventListener("click",()=>{
  const item=libraryItem();
  if(item)window.COFdrawSolubility.open(item.smiles,item.abbreviation||item.label);
});
for(const id of["library-search","library-sites","library-show-all"])$(id).addEventListener(id==="library-search"?"input":"change",drawLibrary);
$("library-target").addEventListener('change',()=>{$("library-show-all").checked=false;drawLibrary();});
document.querySelectorAll("[data-library-target]").forEach(button=>button.addEventListener("click",()=>setLibraryTarget(button.dataset.libraryTarget)));
$("library-editor-mode").addEventListener("click",()=>setLibraryTarget($("library-target").value==="editor"?"a":"editor"));
$("library-apply").addEventListener("click",()=>applyLibraryItem());$("library-edit").addEventListener("click",()=>{const item=libraryItem();if(item){const block=$("library-target").value==="b"?"b":"a";$("monomer-library-dialog").close();openSketch(block,item);}});
$("new-monomer").addEventListener("click",()=>openSketch("a",{smiles:"",full_name:"New monomer"}));
function drawTopologies(){
  const search=$("topology-search").value.toLowerCase().trim(),rows=(state.library?.topology_catalog||[]).filter(row=>!search||JSON.stringify([row.id,row.label,row.description,row.dimension,row.net_code,row.point_group,row.kind]).toLowerCase().includes(search));
  const list=$("topology-list");list.replaceChildren();for(const row of rows){const card=document.createElement("article");card.className="topology-card";const title=document.createElement("h3"),description=document.createElement("p"),semantics=document.createElement("small"),badge=document.createElement("span");title.textContent=row.label||row.id;description.textContent=row.description||row.qualification||"";const pairs=row.connectivity||row.reactive_site_counts||[],symmetry=row.point_group||row.node_symmetry||[];semantics.textContent=[row.dimension?`Dimension: ${row.dimension}`:null,pairs.length?`Reactive-site examples: ${pairs.map(pair=>Array.isArray(pair)?pair.join(" + "):pair).join("; ")}`:null,symmetry.length?`Node symmetry / geometry examples: ${Array.isArray(symmetry)?symmetry.join(", "):symmetry}`:null].filter(Boolean).join(". ");badge.className="reference-badge";badge.textContent=row.construction_supported?"Related explicit template available":"Reference only · no automatic construction";card.append(title,description,semantics,badge);
    if(row.graph){
      const {nodes,edges,dimension}=row.graph,graph=document.createElementNS(NS,'svg');
      graph.setAttribute('viewBox','0 0 200 116');graph.setAttribute('role','img');
      graph.setAttribute('aria-label',row.id+' periodic connectivity graph');
      graph.style.cssText='width:100%;height:112px;display:block';
      const mark=(name,attrs)=>{const element=document.createElementNS(NS,name);for(const [k,v] of Object.entries(attrs))element.setAttribute(k,v);graph.append(element);};
      for(const [a,b] of edges){
        const [x1,y1,d1]=nodes[a],[x2,y2,d2]=nodes[b];
        mark('line',{x1,y1,x2,y2,stroke:'#31545c','stroke-width':1.4,opacity:dimension===2?1:.25+.375*(d1+d2)});
      }
      for(const [x,y,d,degree] of [...nodes].sort((a,b)=>a[2]-b[2])){
        mark('circle',{cx:x,cy:y,r:2.4,fill:degree===3?'#b74e40':degree===4?'#176b76':'#78639a',opacity:dimension===2?1:.4+.6*d});
      }
      title.after(graph);
    }
    const limits=document.createElement("small");limits.textContent=Array.isArray(row.limitations)?row.limitations.join(" "):row.limitations||row.qualification||"";card.append(limits);
    for(const ref of row.source_refs||[]){try{const url=new URL(ref.url);if(!["http:","https:"].includes(url.protocol))continue;const a=document.createElement("a");a.href=url.href;a.target="_blank";a.rel="noopener noreferrer";a.textContent=ref.title||"Primary source";card.append(a);}catch{}}
    list.append(card);
  }
  if(!rows.length){const p=document.createElement("p");p.className="style-empty";p.textContent="No reference entries match this search.";list.append(p);}
}
$("open-topology-library").addEventListener("click",()=>{$("topology-search").value="";drawTopologies();$("topology-dialog").showModal();});$("topology-search").addEventListener("input",drawTopologies);

// Monomer sketch dialog: no optical models or remote editor dependencies.
let activeBlock="a",sketchBusy=false,sketchResult=null,sketchHydrating=false,sketchTextPending=false;
function sketchMessage(text,error=false){$("sketch-message").textContent=text;$("sketch-message").className=`sketch-message${error?" error":""}`;}
function sketchControls(){
  document.querySelectorAll("#sketch-dialog button,#sketch-dialog textarea").forEach(el=>el.disabled=sketchBusy);
  const canvas=$("monomer-canvas");
  canvas.style.pointerEvents=sketchBusy?"none":"";
  canvas.setAttribute("aria-disabled",String(sketchBusy));
  $("atom-charge").disabled=sketchBusy||sketcher.selected?.kind!=="atom";
  $("sketch-dialog").dispatchEvent(new Event("cofdraw:sketch-busy"));
}
const sketcher=new MonomerSketcher($("monomer-canvas"),()=>{
  if(!sketchHydrating){
    sketchResult=null;sketchTextPending=false;$("sketch-smiles").value="";
    $("sketch-formula").textContent="Graph edited";
    sketchMessage("Use this building block to check and apply the sketch.");
  }
},text=>sketchMessage(text,true),atom=>{
  $("atom-charge").disabled=sketchBusy||!atom;
  $("atom-charge").value=atom?.charge??0;
});
function hydrateSketch(result){
  sketchHydrating=true;try{sketcher.load(result.graph);}finally{sketchHydrating=false;}
  sketchResult=result;sketchTextPending=false;$("sketch-smiles").value=result.smiles;
  $("sketch-formula").textContent=result.formula;
  sketchMessage(sketcher.hasStereo()?"Assigned stereo is retained. Move atoms here, or edit the stereospecific SMILES to change connectivity.":"Click an atom to extend a bond, or drag to connect atoms.");
}
async function sketchRequest(payload){sketchBusy=true;sketchControls();try{const result=await api("monomer",payload);hydrateSketch(result);return result;}finally{sketchBusy=false;sketchControls();}}
let sketchAssembly=false;
async function openSketch(block,libraryMonomer=null,assembly=false){
  sketchAssembly=assembly;
  activeBlock=block;
  $("sketch-dialog").dispatchEvent(new Event("cofdraw:open-sketch"));
  $("sketch-title").textContent=libraryMonomer
    ? libraryMonomer.abbreviation||libraryMonomer.full_name||libraryMonomer.id
    : `Building block ${block.toUpperCase()} · ${$(`label-${block}`).value}`;
  $("sketch-use").textContent=libraryMonomer?`Use as building block ${block.toUpperCase()}`:"Use this building block";
  $("sketch-dialog").showModal();sketchResult=null;sketchTextPending=false;sketchHydrating=true;
  sketcher.clear();sketchHydrating=false;sketcher.setTool("bond","1");
  sketchMessage("Draw a bond, import a MOL, or load SMILES.");
  $("sketch-smiles").value=libraryMonomer?libraryMonomer.smiles:$(`smiles-${block}`).value;
  $("sketch-formula").textContent="";
  try{if($("sketch-smiles").value.trim())await sketchRequest({smiles:$("sketch-smiles").value.trim()});}
  catch(e){sketchMessage(e.message,true);}
  sketcher.history=[];sketcher.future=[];$("monomer-canvas").focus({preventScroll:true});
}
document.querySelectorAll(".sketch-open").forEach(b=>b.addEventListener("click",()=>openSketch(b.dataset.block)));
const sketchTools=[...document.querySelectorAll(".sketch-tools button[data-tool]")];
sketchTools.forEach(button=>button.addEventListener("click",()=>{
  sketcher.setTool(button.dataset.tool,button.dataset.value);
  $("monomer-canvas").focus({preventScroll:true});
}));
$("monomer-canvas").addEventListener("sketchtool",({detail:{tool,value}})=>{
  let matched=false;
  for(const button of sketchTools){
    const active=button.dataset.tool===tool&&(button.dataset.value||"")===value;
    button.classList.toggle("active",active);button.setAttribute("aria-pressed",String(active));matched||=active;
  }
  const other=tool==="atom"&&!matched;
  $("sketch-elements").classList.toggle("active",other);
  $("sketch-elements").textContent=other?`Elements · ${value}`:"Elements";
});
$("sketch-elements").addEventListener("click",mountElementPicker($("element-dialog"),symbol=>{
  sketcher.setTool("atom",symbol);$("monomer-canvas").focus({preventScroll:true});
}));
$("sketch-undo").addEventListener("click",()=>sketcher.undo());$("sketch-redo").addEventListener("click",()=>sketcher.redo());$("sketch-clear").addEventListener("click",()=>sketcher.clear());
$("atom-charge").addEventListener("change",event=>{
  const field=event.target;
  if(!field.value||!field.checkValidity()){field.value=sketcher.graph.atoms[sketcher.selected.index].charge;sketchMessage("Enter a whole atom charge from −4 to +8.",true);return;}
  sketcher.charge(field.value);
});
$("sketch-smiles").addEventListener("input",()=>{sketchTextPending=true;sketchMessage("SMILES changed. Load it to replace the sketch, or apply it with Use this building block.");});
$("sketch-load-smiles").addEventListener("click",()=>sketchRequest({smiles:$("sketch-smiles").value.trim()}).catch(e=>sketchMessage(e.message,true)));
async function cleanSketch(){
  if(sketchBusy||sketcher.drag)return;
  try{
    await sketchRequest({...sketchInput(),clean_layout:true});
    sketchMessage("Layout cleaned. Undo restores the previous drawing.");
  }catch(e){sketchMessage(e.message,true);}
}
$("sketch-clean").addEventListener("click",cleanSketch);
$("sketch-dialog").addEventListener("keydown",event=>{
  if((event.ctrlKey||event.metaKey)&&event.shiftKey&&event.key.toLowerCase()==="k"){
    event.preventDefault();cleanSketch();
  }
});
$("sketch-import").addEventListener("click",()=>$("mol-file").click());
$("mol-file").addEventListener("change",async e=>{try{const file=e.target.files[0];if(!file)return;if(file.size>500000)throw new Error("MOL import limit is 500 kB.");const text=await file.text();if(text.includes("$$$$"))throw new Error("Import a single-molecule MOL file; split multi-record SDF files first.");await sketchRequest({molblock:text});}catch(error){sketchMessage(error.message,true);}finally{e.target.value="";}});
function sketchInput(){return sketchTextPending?{smiles:$("sketch-smiles").value.trim()}:{graph:sketcher.exportGraph()};}
async function checkedSketch(){return !sketchTextPending&&sketchResult?sketchResult:sketchRequest(sketchInput());}
$("sketch-solubility").addEventListener("click",async()=>{
  try{
    const result=await checkedSketch();
    window.COFdrawSolubility.open(result.smiles,`${$("sketch-title").textContent} · ${result.formula}`);
  }catch(e){sketchMessage(e.message,true);}
});
$("sketch-export").addEventListener("click",async()=>{try{const r=await checkedSketch();await saveText(`COFdraw-monomer-${activeBlock.toUpperCase()}.mol`,r.molblock,"chemical/x-mdl-molfile");}catch(e){sketchMessage(e.message,true);}});
async function applySketch(){
  const result=await checkedSketch(), block=activeBlock;
  const known=(state.library?.monomer_presets||[]).find(item=>item.canonical_smiles===result.smiles);
  if(sketchAssembly){
    const item=known||await api('assembly/precursor',{smiles:result.smiles,label:'Custom '+block.toUpperCase()});
    const roles=item.assembly_roles||{},ids=(multi()&&assemblyMode()==='3'?['a','b','c']:['a']).map(k=>$('assembly-'+k).value);
    const index='abc'.indexOf(block);
    if(ids.length===1&&!roles.self_kind&&!roles.self_imine)throw new Error('Choose a precursor with a supported self-condensation reaction.');
    if(ids.length===3){
      if(index===0&&(![roles.primary_arylamine,roles.aryl_aldehyde].some(n=>[3,4].includes(n))||item.geometry_hint?.startsWith('tetrahedral')||roles.self_imine))throw new Error('Node A needs three or four planar reactive sites for imine formation.');
      if(index>0){
        const rows=[...state.library.monomer_presets,...customAssemblyRows.values()],node=rows.find(r=>r.id===ids[0]);
        const required=node.assembly_roles.primary_arylamine>=3?'aryl_aldehyde':'primary_arylamine';
        if(roles[required]!==2||roles.self_imine)throw new Error('This linker needs two complementary '+required.replaceAll('_',' ')+' sites.');
        if($('mixed-layout').value==='kgm'&&!roles.kgm_linear_dialdehyde)throw new Error('The Kagome arrangement needs a linear dialdehyde. Use Local pore design for a bent linker.');
        if(rows.find(r=>r.id===ids[3-index])?.canonical_smiles===item.canonical_smiles)throw new Error('Choose two chemically distinct linkers.');
      }
    }
    if(!known)customAssemblyRows.set(item.id,item);
    ids[index]=item.id;
    const saved={monomer_ids:ids,linker_sequence:$('assembly-pattern').value,pore_layout:$('mixed-layout').value,pore_view:$('self-face').value};
    $('sketch-dialog').close();assemblyChanged(saved);return;
  }
  const previous=checkedMonomer(block);
  if(known)$(`label-${block}`).value=known.abbreviation||known.full_name||known.id;
  else if(previous?.smiles!==result.smiles)$(`label-${block}`).value="Custom "+block.toUpperCase();
  $(`smiles-${block}`).value=result.smiles;
  $(`preset-${block}`).value=known?.id||"";
  setPreview(block,result);markInputsChanged();$("sketch-dialog").close();
  notice(`Building block ${block.toUpperCase()} updated · ${result.formula}. Updating drawing.`);queueDrawing();
}
$("sketch-use").addEventListener("click",()=>applySketch().catch(e=>sketchMessage(e.message,true)));
$("sketch-cancel").addEventListener("click",()=>$("sketch-dialog").close());

async function applyDescriptionPlan(plan){
  if(!state.online||state.busy)throw new Error('Wait for the local drawing engine to be ready.');
  if(plan?.status!=='ready'||Boolean(plan.source_recipe_id)===Boolean(plan.prepare_request))throw new Error('Analyze and review a complete drawing plan first.');
  window.COFdrawCrystalDesign?.leave();
  state.busy=true;sync();
  try{
    if(plan.source_recipe_id){
      const source=await api(`sources/${encodeURIComponent(plan.source_recipe_id)}`);
      acceptRecipe(source.recipe,plan.source_recipe_id,source.metadata);showSourceMonomers(source.recipe,source.metadata);
    }else{
      const request=clone(plan.prepare_request),prepared=await api('prepare',request);
      // A reviewed description plan contains explicit choices. Keep those
      // choices visible and editable rather than silently rematching them.
      state.autoChemistry=false;
      for(const block of ['a','b']){$(`smiles-${block}`).value=request[`monomer_${block}`];$(`label-${block}`).value=request[`label_${block}`]||block.toUpperCase();$(`preset-${block}`).value='';}
      $('linkage-search').value='';populateLinkages(request.linkage);$('cof-preset').value='custom';
      acceptRecipe(prepared.recipe,null,null,true);await refreshMonomers();refreshChemistry(request.template_id);state.inputsDirty=false;
    }
    state.style=plan.style_id||'standard';state.seed=null;state.styleOptions=clone(plan.style_options||{});state.resolvedStyleOptions={};setStyleDescription();drawStyleControls();
    state.busy=false;sync();
    const result=await render();
    if(!result)throw new Error($('notice').textContent||'The drawing could not be rendered. Review the plan and try again.');
    return result;
  }finally{state.busy=false;sync();}
}

async function applyCofexAssembly(request){
  if(!state.online||state.busy)throw new Error('Wait for the current drawing to finish.');
  clearTimeout(liveDrawingTimer);
  const revision=++selectionRevision;
  state.busy=true;sync();
  try{
    const prepared=await api('prepare/assembly',request);
    const rows=await Promise.all((request.custom_precursors||[]).map(item=>api('assembly/precursor',item)));
    if(revision!==selectionRevision)throw new Error('The selected precursors changed during import.');
    for(const row of rows)customAssemblyRows.set(row.id,row);
    window.COFdrawCrystalDesign?.leave();
    $('cof-preset').value='custom';
    acceptRecipe(prepared.recipe,null,null,true);
    state.busy=false;sync();
    if(!await render())throw new Error($('notice').textContent||'The assembly could not be drawn.');
  }finally{state.busy=false;sync();}
}

async function initialize(){
  try{
    const library=await api("library");state.library=library;state.online=true;$("engine-dot").classList.add("live");$("engine-status").textContent="Local chemistry engine · ready";$("catalog-version").textContent=library.catalog_version||"Built-in drawing library";
    const presets=library.cof_presets?.length?library.cof_presets:fallbackPresets;library.cof_presets=presets;combinations=L.startingCombinations(library);populateLinkages("imine");refreshChemistry();
    drawStyles();drawStyleControls();
    if(library.source_registry_errors?.length)notice("Some source recipes failed integrity checks and are omitted from the library.",true);
    setCombination(presets[0]);state.inputsDirty=false;setStyleDescription();window.dispatchEvent(new CustomEvent('cofdraw:library',{detail:library}));
  }catch(e){$("engine-dot").classList.add("offline");$("engine-status").textContent="Local engine unavailable";notice("Open COFdraw through its Windows launcher or the local COF Vision server. The drawing engine is not available on this static page.",true);}
  sync();
}
function selectedFigureName(){
  const key=$('cof-preset').value;
  if(window.COFdrawCrystalDesign?.active){
    if(key!=='crystal:'+window.COFdrawCrystalDesign.active)return null;
    const row=crystalRows.find(r=>'crystal:'+r.id===key);
    const name=row?.starter?.common_name;
    return name?(window.COFdrawCrystalDesign.isReference?name:name+' precursor design'):null;
  }
  if(state.inputsDirty||state.jsonPending)return null;
  const row=combinations.find(r=>r.id===key);
  return state.sourceMetadata?.literature_name||row?.commonName||null;
}
function selectedCatalogPair(){
  return ['a','b'].map(block=>{const input=$('smiles-'+block).value.trim();return state.library?.monomer_presets?.find(row=>row.smiles===input||row.canonical_smiles===input)||checkedMonomer(block);});
}
window.COFdrawStudio={state,periodicPair:()=>periodicLinkage(...selectedCatalogPair()),planarPair:()=>L.compatibleLinkages(state.library||{},...selectedCatalogPair()).length>0,selectedFigureName,refreshStyles(){drawStyles();setStyleDescription();drawStyleControls();},loadSource,precursorSVG:async smiles=>(await libraryPreviews.request(smiles)).result.svg,showExports,selectedDescriptor,api,notice,render,acceptRecipe,removeIdentity,syncJSON,sync,drawLayout,drawLayoutTable,changedLayout,saveText,saveExport,chooseStyle,refreshChemistry,showSourceMonomers,drawStyleControls,applyDescriptionPlan,geometryReviewMessage,setCrystalMonomers,queueDrawing};
window.COFdrawStudio.applyCofexAssembly=applyCofexAssembly;
initialize();
