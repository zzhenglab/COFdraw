/* The 2D drawing and the periodic viewer share one selected crystal report. */
(()=>{
  'use strict';
  const $=id=>document.getElementById(id), studio=window.COFdrawStudio, lab=window.COFdrawScience;
  const main=document.querySelector('main.studio'), picker=$('crystal-linker'), picture=$('crystal-drawing-image');
  const entries=new Map(), previews=[];
  const matched={
    'cof-300-c5':['tapm','ta'], 'curated-13040N3':['tapm','bpda_cho'], 'dia-imine':['tapm','ta']
  };
  const custom='custom-dia';
  let id=null, edit=null, revision=0, referenceRevision=0, report=null, svgText='', svgDownload='', views=[], timer, failure='';
  let restoredReference=false, comparing=false;
  function blob(text){
    const url=URL.createObjectURL(new Blob([text],{type:'image/svg+xml'}));previews.push(url);return url;
  }
  function clearPicture(){
    for(const url of previews)URL.revokeObjectURL(url);
    previews.length=0;report=null;svgText='';svgDownload='';views=[];comparing=false;$('crystal-drawing-view').replaceChildren();
    picture.hidden=true;picture.removeAttribute('src');
    $('save-crystal-drawing').disabled=true;$('rotate-crystal-design').disabled=true;
    $('crystal-drawing-caption').textContent='';
  }
  function leave(){
    if(!id)return;
    clearTimeout(timer);revision++;id=null;edit=null;failure='';clearPicture();main.classList.remove('crystal-design-active');
    $('construction-target').value='patch';
    $('crystal-design-controls').hidden=true;$('crystal-drawing-panel').hidden=true;
    if($('cof-preset').value.startsWith('crystal:'))$('cof-preset').value='custom';
    $('crystal-template').querySelector('option[value="custom-dia"]')?.remove();
    if(lab.origin==='catalog-design')lab.invalidateCrystalDesign('');
    studio.refreshChemistry();
  }
  function linkerChoices(row){
    const options=Object.entries(row.linker_series?.choices||{[row.id]:row.model_name});
    picker.replaceChildren();
    for(const [key,label] of options){
      picker.add(new Option(label,key));
    }
    if(row.aryl_f_sites)picker.add(new Option('Compare aryl-F reference derivative · '+row.aryl_f_sites+' sites/cell','aryl_f'));
    picker.value=edit||row.id;
    $('crystal-linker-note').textContent=edit?'Replace matching aromatic C–H sites with C–F in the selected cell.':
      options.length>1?'Shared node · '+row.linker_series.shared_node:
      row.aryl_f_sites?'Aryl substitution uses explicit C–H sites and checks contacts in this crystal.':'No verified alternative linker model for this crystal.';
  }
  function showPeriodic(){
    if(!id)return;
    comparing=false;referenceRevision++;
    if(!report){lab.status(failure||$('crystal-drawing-progress').textContent||'Waiting for the selected crystal coordinates.');return;}
    if(lab.structure!==report.structure||lab.origin!=='catalog-design')
      lab.showCatalogDesign(report,restoredReference?id:custom);
    requestAnimationFrame(()=>lab.viewer.draw());
  }
  function precursor(block){
    return ['a','b'].includes(block)?{smiles:$('smiles-'+block).value.trim(),label:$('label-'+block).value.trim()}:null;
  }
  function choice(item,partner){
    const linkage=$('dia-linkage').value;
    function roles(m){
      if(!m)return [];
      return (studio.state.library?.linkages?.[linkage]?.roles||[]).flatMap(role=>{
        const n=window.COFdrawLibrary.reactiveCount(m,linkage,role);
        if(![2,4].includes(n)||(n===4&&['square_planar','planar_square'].includes(m.geometry_hint)))return [];
        return [{role,n}];
      });
    }
    const mine=roles(item),other=roles(partner);
    const compatible=mine.length>0&&(!partner||mine.some(a=>other.some(b=>a.role!==b.role&&a.n+b.n===6)));
    return {compatible,reason:compatible?'Reactive site counts match. The cell builder checks tetrahedral and linker geometry.':
      'Choose complementary groups for this linkage on a four-arm tetrahedral node and a two-ended linker.'};
  }
  function drawView(){
    const view=views[$('crystal-drawing-view').selectedIndex];
    if(!view)return;
    svgText=view.text||'';svgDownload=view.download||'';
    picture.src=view.url;picture.hidden=false;$('crystal-drawing-caption').textContent=view.caption;
  }
  function enterDrawing(){
    if(!$('science-dialog').hidden)$('close-science').click();
    if(!$('space-workspace').hidden)window.dispatchEvent(new Event('cofdraw:leave-explorer'));
    if(!$('science-dialog').hidden||main.hidden)$('close-science').click();
    $('mode-manual').click();
    main.classList.add('crystal-design-active');$('crystal-design-controls').hidden=false;$('crystal-drawing-panel').hidden=false;
    $('construction-target').value='dia';
  }
  function startCustom(){
    id=custom;edit=null;failure='';restoredReference=false;clearTimeout(timer);revision++;clearPicture();
    $('cof-preset').value='crystal:'+custom;$('compare-crystal-reference').hidden=true;
    $('crystal-drawing-title').textContent=[precursor('a').label||'A',precursor('b').label||'B'].join(' + ');
    const nets=+$('dia-net-count').value;
    $('crystal-drawing-kind').textContent=$('reticular-net').selectedOptions[0].textContent+' · '+nets+(nets===1?' net':' interpenetrating nets');
    picker.replaceChildren(new Option('Current A + B · '+$('reticular-net').value,custom));
    $('crystal-linker-note').textContent='Changing either monomer, linkage or network count rebuilds both views.';
    $('crystal-drawing-progress').textContent='Preparing the selected precursor pair…';
    lab.invalidateCrystalDesign('The monomer pair changed.');studio.sync();
  }
  function inputsChanged(){
    if(!id)return;
    startCustom();
    if(!precursor('a').smiles||!precursor('b').smiles){
      $('crystal-drawing-progress').textContent='Choose or draw both building blocks.';return;
    }
    timer=setTimeout(build,600);
  }
  async function build(){
    if(!id)return;
    startCustom();
    await generate(null,revision);
  }
  async function generate(row,token){
    const a=precursor('a'),b=precursor('b');
    try{
      if(!a.smiles||!b.smiles)throw Error('Choose or draw both building blocks.');
      const nets=+$('dia-net-count').value;
      const topology=$('reticular-net').value,path=topology==='dia'?'science/build-dia':'science/build-reticular';
      const result=await studio.api(path,{monomer_a:a.smiles,monomer_b:b.smiles,label_a:a.label||'A',label_b:b.label||'B',linkage:$('dia-linkage').value,...(topology==='dia'?{}:{topology}),...(nets>1?{network_count:nets}:{})});
      if(token===revision)await present(result,row,token);
    }catch(error){if(token===revision)showError(error);}
  }
  function showError(error){
    failure=error.message;$('crystal-drawing-progress').textContent=failure;
    if(!$('science-dialog').hidden)lab.status(failure);
    studio.sync();
  }
  async function select(key,modification=null,restoreReference=false){
    if(key===custom){enterDrawing();startCustom();studio.setCrystalMonomers();await build();return;}
    const row=entries.get(key);if(!row)return;
    clearTimeout(timer);const token=++revision;
    id=key;edit=modification;failure='';restoredReference=restoreReference;clearPicture();enterDrawing();$('dia-linkage').value='imine';$('dia-net-count').value='1';$('reticular-net').value='dia';
    const ids=edit?[]:(row.monomer_ids||matched[key]||[]);
    const items=ids.map(k=>studio.state.library?.monomer_presets?.find(m=>m.id===k));
    studio.setCrystalMonomers(items);
    $('cof-preset').value='crystal:'+id;linkerChoices(row);
    $('compare-crystal-reference').hidden=restoreReference;
    $('crystal-drawing-title').textContent=restoreReference?row.model_name:[precursor('a').label,precursor('b').label].filter(Boolean).join(' + ');
    $('crystal-drawing-kind').textContent=restoreReference?'Saved reference model':'Generated diamond design · single net · '+row.model_name+' precursors';
    $('crystal-drawing-progress').textContent='Preparing the drawing and unit cell…';
    lab.invalidateCrystalDesign('Loading the selected crystal.');
    if(!restoreReference){await generate(row,token);return;}
    try{
      const result=await studio.api('science/model',{template_id:key,...(edit?{linker_edit:edit}:{})});
      if(token!==revision)return;
      await present(result,row,token);
    }catch(error){if(token===revision)showError(error);}
  }
  async function present(result,row,token){
      const choices=[];
      if(result.chemical_framework_svg)choices.push({label:'Framework drawing',text:result.chemical_framework_svg,url:blob(result.chemical_framework_svg),caption:'Bond-line scheme of a closed 3D framework circuit. Ring geometry is drawn in the plane.'});
      for(const [key,label,caption] of [
        ['chemical_circuit_svg','Projected circuit','Cartesian projection of a closed covalent circuit and its junction branches.'],
        ['chemical_drawing_svg','Linker detail','Local linker connectivity.']
      ]){
        const text=result[key];
        if(text)choices.push({text,label,url:blob(text),caption:caption+' Stars mark continuation. Bond markings are inferred from the CIF.'});
      }
      const generatedChoices=choices.length;
      if(row?.source_recipe_id&&!edit){
        try{
          const rendered=await studio.api('render',{source_recipe_id:row.source_recipe_id,style_id:'standard',size:1000});
          if(token!==revision)return;
          choices.push({label:'Literature drawing (comparison)',download:rendered.downloads.svg,url:rendered.downloads.svg,caption:'Source-annotated drawing for comparison · '+row.model_name});
        }catch(error){if(token!==revision)return;studio.notice('Source annotation unavailable. '+error.message,true);}
      }
      report=result;views=choices;studio.sync();
      for(const choice of choices)$('crystal-drawing-view').add(new Option(choice.label));
      if(!generatedChoices)$('crystal-drawing-view').selectedIndex=-1;
      drawView();$('crystal-drawing-progress').textContent=generatedChoices?'':'Periodic cell ready. A chemical projection could not be resolved.';
      $('save-crystal-drawing').disabled=!generatedChoices;$('rotate-crystal-design').disabled=false;
      if(!$('science-dialog').hidden&&!comparing)showPeriodic();
  }
  window.addEventListener('cofdraw:crystal-catalog',e=>{for(const r of e.detail)entries.set(r.id,r);});
  window.addEventListener('cofdraw:open-crystal',e=>select(e.detail));
  picker.addEventListener('change',()=>{if(picker.value==='aryl_f'){picker.value=id;selectInViewer(id,'aryl_f');}else select(picker.value);});
  $('compare-crystal-reference').onclick=()=>selectInViewer(id);
  $('crystal-drawing-view').addEventListener('change',drawView);
  $('rotate-crystal-design').addEventListener('click',()=>$('open-science').click());
  function returnToPatch(){leave();studio.queueDrawing();}
  $('leave-crystal-design').addEventListener('click',returnToPatch);
  $('construction-target').onchange=()=>{if($('construction-target').value==='dia')select(custom);else returnToPatch();};
  $('dia-linkage').onchange=()=>{inputsChanged();studio.refreshChemistry();};
  $('dia-net-count').onchange=inputsChanged;
  $('reticular-net').onchange=()=>{inputsChanged();studio.refreshChemistry();};
  for(const element of ['mode-describe','load-source','new-project'])$(element)?.addEventListener('click',leave,true);
  $('save-crystal-drawing').addEventListener('click',()=>{
    const name=$('crystal-drawing-title').textContent;
    (svgDownload?studio.saveExport(svgDownload,'svg'):studio.saveText(name+'.svg',svgText,'image/svg+xml')).catch(e=>studio.notice(e.message,true));
  });
  async function restore(project){
    if(project.template_id===custom){
      if(!Array.isArray(project.monomers)||project.monomers.length!==2||project.monomers.some(m=>typeof m?.smiles!=='string'||m.smiles.length>4000||typeof m.label!=='string'||m.label.length>100))
        throw new Error('The diamond project needs two editable precursors.');
      const linkage=project.linkage||'imine';
      if(![...$('dia-linkage').options].some(o=>o.value===linkage))throw Error('Unknown diamond linkage.');
      const count=project.network_count??1;
      if(!Number.isInteger(count)||count<1||count>7)throw Error('Diamond network count must be between 1 and 7.');
      const topology=project.topology||'dia';
      if(!['dia','lon'].includes(topology))throw Error('Unknown reticular net.');
      enterDrawing();$('reticular-net').value=topology;$('dia-linkage').value=linkage;$('dia-net-count').value=String(count);startCustom();studio.setCrystalMonomers(project.monomers);await build();
    }else if(!entries.has(project.template_id)||![null,'aryl_f'].includes(project.linker_edit??null))
      throw new Error('Unknown crystal model or linker modification in this project.');
    else await select(project.template_id,project.linker_edit??null,project.version===1||project.reference_only===true);
    if(report&&project.view&&[...$('crystal-drawing-view').options].some(o=>o.value===project.view)){
      $('crystal-drawing-view').value=project.view;drawView();
    }
  }
  function save(){
    if(report)studio.saveText('COFdraw-crystal-project.json',JSON.stringify({format:'cofdraw-crystal-project',
      version:2,template_id:id,linker_edit:edit,reference_only:restoredReference,...(id===custom?{monomers:['a','b'].map(precursor),linkage:$('dia-linkage').value,network_count:+$('dia-net-count').value,topology:$('reticular-net').value}:{}),view:$('crystal-drawing-view').value},null,2)).catch(e=>studio.notice(e.message,true));
  }
  async function selectInViewer(key,modification=null){
    if(!entries.has(key))return;
    if($('science-dialog').hidden)$('open-science').click();
    const token=revision,ref=++referenceRevision;comparing=true;
    lab.invalidateCrystalDesign('Loading the comparison CIF…');
    try{
      const result=await studio.api('science/model',{template_id:key,...(modification?{linker_edit:modification}:{})});
      if(token!==revision||ref!==referenceRevision||$('science-dialog').hidden)return;
      lab.showCatalogDesign(result,key);
      lab.status('Comparison CIF · '+entries.get(key).model_name+(modification?' · aryl-F design':''));
    }catch(e){if(token===revision&&ref===referenceRevision)lab.status(e.message,true);}
  }
  window.COFdrawCrystalDesign={restore,save,precursor,choice,build,inputsChanged,selectInViewer,has:key=>entries.has(key),get active(){return id;},get report(){return report;},leave,showPeriodic};
})();
