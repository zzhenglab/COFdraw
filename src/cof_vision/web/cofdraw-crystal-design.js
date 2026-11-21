/* The 2D drawing and the periodic viewer share one selected crystal report. */
(()=>{
  'use strict';
  const $=id=>document.getElementById(id), studio=window.COFdrawStudio, lab=window.COFdrawScience;
  const main=document.querySelector('main.studio'), picker=$('crystal-linker'), picture=$('crystal-drawing-image');
  const entries=new Map(), previews=[];
  const matched={
    'cof-300-c5':['tapm','ta'], 'curated-13040N3':['tapm','bpda_cho'], 'dia-imine':['tapm','ta']
  };
  let id=null, edit=null, revision=0, report=null, svgText='', svgDownload='', views=[];
  function blob(text){
    const url=URL.createObjectURL(new Blob([text],{type:'image/svg+xml'}));previews.push(url);return url;
  }
  function clearPicture(){
    for(const url of previews)URL.revokeObjectURL(url);
    previews.length=0;report=null;svgText='';svgDownload='';views=[];$('crystal-drawing-view').replaceChildren();
    picture.hidden=true;picture.removeAttribute('src');
    $('save-crystal-drawing').disabled=true;$('rotate-crystal-design').disabled=true;
    $('crystal-drawing-caption').textContent='';$('crystal-precursors').replaceChildren();
  }
  function leave(){
    if(!id)return;
    revision++;id=null;edit=null;clearPicture();main.classList.remove('crystal-design-active');studio.sync();
    $('crystal-design-controls').hidden=true;$('crystal-drawing-panel').hidden=true;
    if($('cof-preset').value.startsWith('crystal:'))$('cof-preset').value='custom';
    if(lab.origin==='catalog-design')lab.invalidate('');
  }
  function linkerChoices(row){
    const options=Object.entries(row.linker_series?.choices||{[row.id]:row.model_name});
    picker.replaceChildren();
    for(const [key,label] of options){
      picker.add(new Option(label,key));
    }
    if(row.aryl_f_sites)picker.add(new Option('Aryl-F derivative · '+row.aryl_f_sites+' sites/cell','aryl_f'));
    picker.value=edit||row.id;
    $('crystal-linker-note').textContent=edit?'Replace matching aromatic C–H sites with C–F in the selected cell.':
      options.length>1?'Shared node · '+row.linker_series.shared_node:
      row.aryl_f_sites?'Aryl substitution uses explicit C–H sites and checks contacts in this crystal.':'No verified alternative linker model for this crystal.';
  }
  async function showPrecursors(row,token){
    const ids=row.monomer_ids||matched[row.id];
    if(!ids||edit)return;
    for(const key of ids){
      const item=studio.state.library?.monomer_presets?.find(m=>m.id===key);
      if(!item)continue;
      const molecule=await studio.api('monomer',{smiles:item.smiles});
      if(revision!==token)return;
      const fig=document.createElement('figure'), img=document.createElement('img'),caption=document.createElement('figcaption');
      img.alt=item.label;img.src=blob(molecule.svg);caption.textContent=item.label;
      fig.append(img,caption);$('crystal-precursors').append(fig);
    }
  }
  function showPeriodic(){
    if(!id)return;
    if(!report){lab.status('Waiting for the selected crystal coordinates.');return;}
    if(lab.structure!==report.structure||lab.origin!=='catalog-design')lab.showCatalogDesign(report,id);
    requestAnimationFrame(()=>lab.viewer.draw());
  }
  function precursor(block){
    const row=entries.get(id), ids=!edit&&(row?.monomer_ids||matched[id]);
    const key=ids?.[block==='a'?0:1];
    return key?studio.state.library?.monomer_presets?.find(m=>m.id===key):null;
  }
  function drawView(){
    const view=views[$('crystal-drawing-view').selectedIndex];
    if(!view)return;
    svgText=view.text||'';svgDownload=view.download||'';
    picture.src=view.url;picture.hidden=false;$('crystal-drawing-caption').textContent=view.caption;
  }
  async function select(key,modification=null){
    const row=entries.get(key);if(!row)return;
    const token=++revision;
    id=key;edit=modification;clearPicture();
    if(!$('science-dialog').hidden)$('close-science').click();
    if(!$('space-workspace').hidden)window.dispatchEvent(new Event('cofdraw:leave-explorer'));
    if(!$('science-dialog').hidden||main.hidden)$('close-science').click();
    $('mode-manual').click();
    main.classList.add('crystal-design-active');$('crystal-design-controls').hidden=false;$('crystal-drawing-panel').hidden=false;
    $('cof-preset').value='crystal:'+id;linkerChoices(row);
    const design=Boolean(edit||row.linker_label||key==='dia-imine');
    $('crystal-drawing-title').textContent=(key==='dia-imine'?'TAPM + TA':row.model_name)+(edit?' · aryl-F design':'');
    $('crystal-drawing-kind').textContent=design?'Unrelaxed periodic design':'Literature crystal';
    $('crystal-drawing-progress').textContent='Preparing the drawing and unit cell…';
    lab.invalidate('Loading the selected crystal.');
    showPrecursors(row,token).catch(()=>{});
    try{
      const result=await studio.api('science/model',{template_id:key,...(edit?{linker_edit:edit}:{})});
      if(token!==revision)return;
      const choices=[];
      if(row.source_recipe_id&&!edit){
        const rendered=await studio.api('render',{source_recipe_id:row.source_recipe_id,style_id:'standard',size:1000});
        if(token!==revision)return;
        choices.push({label:'Annotated framework',download:rendered.downloads.svg,url:rendered.downloads.svg,caption:'Manually aligned source drawing · COF-300'});
      }
      for(const [key,label,caption] of [
        ['chemical_circuit_svg','Framework circuit','One closed covalent circuit and its junction branches.'],
        ['chemical_drawing_svg','Linker detail','Local linker connectivity.']
      ]){
        const text=result[key];
        if(text)choices.push({text,label,url:blob(text),caption:caption+' Stars mark continuation. Bond markings are inferred from the CIF.'});
      }
      if(!choices.length)throw new Error('A connected chemical drawing could not be resolved for these coordinates.');
      report=result;views=choices;studio.sync();
      for(const choice of choices)$('crystal-drawing-view').add(new Option(choice.label));
      drawView();$('crystal-drawing-progress').textContent='';
      $('save-crystal-drawing').disabled=false;$('rotate-crystal-design').disabled=false;
      if(!$('science-dialog').hidden)showPeriodic();
    }catch(error){
      if(token===revision)$('crystal-drawing-progress').textContent=error.message;
    }
  }
  window.addEventListener('cofdraw:crystal-catalog',e=>{for(const r of e.detail)entries.set(r.id,r);});
  window.addEventListener('cofdraw:open-crystal',e=>select(e.detail));
  picker.addEventListener('change',()=>picker.value==='aryl_f'?select(id,'aryl_f'):select(picker.value));
  $('crystal-drawing-view').addEventListener('change',drawView);
  $('rotate-crystal-design').addEventListener('click',()=>$('open-science').click());
  $('leave-crystal-design').addEventListener('click',leave);
  for(const element of ['mode-describe','load-source','new-project'])$(element)?.addEventListener('click',leave,true);
  $('save-crystal-drawing').addEventListener('click',()=>{
    const name=$('crystal-drawing-title').textContent;
    (svgDownload?studio.saveExport(svgDownload,'svg'):studio.saveText(name+'.svg',svgText,'image/svg+xml')).catch(e=>studio.notice(e.message,true));
  });
  async function restore(project){
    if(!entries.has(project.template_id)||![null,'aryl_f'].includes(project.linker_edit??null))
      throw new Error('Unknown crystal model or linker modification in this project.');
    await select(project.template_id,project.linker_edit??null);
    if(report&&project.view&&[...$('crystal-drawing-view').options].some(o=>o.value===project.view)){
      $('crystal-drawing-view').value=project.view;drawView();
    }
  }
  function save(){
    if(report)studio.saveText('COFdraw-crystal-project.json',JSON.stringify({format:'cofdraw-crystal-project',
      version:1,template_id:id,linker_edit:edit,view:$('crystal-drawing-view').value},null,2)).catch(e=>studio.notice(e.message,true));
  }
  async function selectInViewer(key){
    if(key!==id){
      await select(key);
      if(!report||key!==id)return;
    }
    if($('science-dialog').hidden)$('open-science').click();else showPeriodic();
  }
  window.COFdrawCrystalDesign={restore,save,precursor,selectInViewer,has:key=>entries.has(key),get active(){return id;},get report(){return report;},leave,showPeriodic};
})();
