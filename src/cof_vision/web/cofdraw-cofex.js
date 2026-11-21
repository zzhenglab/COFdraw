'use strict';
document.addEventListener('DOMContentLoaded',()=>{
  const $=id=>document.getElementById(id), studio=window.COFdrawStudio;
  const dialog=$('cofex-dialog'),input=$('cofex-input'),menu=$('cofex-choices'),message=$('cofex-message'),draw=$('cofex-draw');
  const open=document.createElement('button');
  open.type='button';open.id='open-cofex';open.className='button quiet';open.textContent='Paste COFex';
  $('open-project').before(open);
  function footer(id,after){
    const section=document.createElement('section'),head=document.createElement('div'),heading=document.createElement('span');
    const copy=document.createElement('button'),value=document.createElement('div');
    section.id=id;section.hidden=true;section.className='cofex-footer';section.setAttribute('aria-label','COFex structure record');
    head.className='identity-heading-row';heading.className='identity-notation-label';
    copy.type='button';copy.textContent='Copy COFex';copy.className='text-button';value.className='code-value';value.tabIndex=0;
    head.append(heading,copy);section.append(head,value);after.after(section);
    return {section,heading,copy,value};
  }
  const crystal=footer('crystal-drawing-cofex',document.querySelector('.crystal-drawing-bottom'));
  const model=footer('model-cofex',$('viewer-distance'));
  async function copyText(text,element){
    if(!text)return;
    try{await navigator.clipboard.writeText(text);}
    catch{
      const box=document.createElement('textarea');box.className='cofex-clipboard';box.value=text;
      document.body.append(box);box.select();
      let copied=false;
      try{copied=document.execCommand('copy');}finally{box.remove();}
      if(!copied){
        element.textContent=text;element.focus();
        const range=document.createRange();range.selectNodeContents(element);
        const selection=getSelection();selection.removeAllRanges();selection.addRange(range);
        studio.notice('The COFex record is selected. Use Copy to copy it.');return;
      }
    }
    studio.notice('COFex copied.');
  }
  function setFooter(panel,record,structureShown=false){
    panel.section.hidden=!record&&!structureShown;
    panel.copy.hidden=!record;
    if(!record){
      panel.heading.textContent='COFex';panel.value.textContent='No precursor COFex is recorded for this structure.';
      panel.copy.onclick=null;return;
    }
    panel.heading.textContent=record.scope==='source'?'Source COFex'+(record.id?' · '+record.id:record.doi?' · DOI '+record.doi:''):'Design COFex · '+record.linkage.replaceAll('_',' ')+' · '+record.net;
    panel.value.textContent=record.compact;panel.copy.onclick=()=>copyText(record.text,panel.value);
  }
  function drawingRecord(){
    const state=studio.state;
    if(state.inputsDirty||state.jsonPending||state.dirty)return null;
    return state.result?.audit?.drawing_cofex;
  }
  function sync(){
    const design=window.COFdrawCrystalDesign, science=window.COFdrawScience;
    const record=drawingRecord(),button=$('copy-drawing-cofex');
    button.disabled=!record;
    button.onclick=()=>copyText(record?.text,$('identity-string'));
    setFooter(crystal,design?.report?.structure?.provenance?.drawing_cofex,Boolean(design?.report));
    const current=Boolean(science?.structure||science?.finiteModel);
    const fromDrawing=current&&['drawing','current-periodic'].includes(science.origin);
    setFooter(model,science?.structure?.provenance?.drawing_cofex||(fromDrawing?record:null),current&&science.origin!=='monomer');
    open.disabled=!studio.state.online;
  }
  for(const type of ['cofdraw:state','cofdraw:rendered','cofdraw:structure','cofdraw:crystalviewclear'])window.addEventListener(type,sync);
  let plans=[],revision=0;
  function reset(){revision++;plans=[];menu.hidden=true;menu.replaceChildren();message.textContent='';draw.textContent='Draw structure';draw.disabled=false;}
  input.addEventListener('input',reset);
  open.onclick=()=>{reset();dialog.showModal();input.focus();input.select();};
  dialog.addEventListener('close',reset);
  function showDrawing(){
    if(!$('science-dialog').hidden)$('close-science').click();
    if(!$('space-workspace').hidden)window.dispatchEvent(new Event('cofdraw:leave-explorer'));
    $('mode-manual').click();
  }
  draw.onclick=async()=>{
    const token=revision;
    try{
      if(!input.value.trim())throw Error('Paste a COFex record first.');
      draw.disabled=true;
      if(!plans.length){
        message.textContent='Checking precursors and connectivity…';
        const result=await studio.api('cofex/resolve',{text:input.value.trim()});
        if(token!==revision)return;
        plans=result.choices;
        if(plans.length>1){
          menu.replaceChildren(...plans.map((p,i)=>new Option(p.label,String(i))));menu.hidden=false;menu.focus();
          message.textContent='This record allows more than one drawing. Choose the product or arrangement.';
          draw.textContent='Draw selected structure';return;
        }
      }
      const plan=plans[Number(menu.value)||0];
      message.textContent='Building the structure…';showDrawing();
      if(plan.crystal_project){
        await window.COFdrawCrystalDesign.restore(plan.crystal_project);
        if(!window.COFdrawCrystalDesign.report)throw Error($('crystal-drawing-progress').textContent||'The periodic model could not be built.');
      }else if(plan.assembly_request){
        await studio.applyCofexAssembly(plan.assembly_request);
      }else{
        await studio.applyDescriptionPlan({...plan,status:'ready',style_id:studio.state.style==='source'?'standard':studio.state.style,style_options:studio.state.styleOptions});
      }
      if(token===revision){dialog.close();studio.notice('Drawing restored from COFex.');}
    }catch(error){if(token===revision)message.textContent=error.message;}
    finally{if(token===revision)draw.disabled=false;}
  };
  sync();
});
