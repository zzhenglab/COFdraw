/* Show each drawing extent on the same canvas, with exports from the visible graph. */
document.addEventListener('DOMContentLoaded',()=>{
  const studio=window.COFdrawStudio, state=studio.state, $=id=>document.getElementById(id);
  const bar=$('extent-bar'), choices=$('extent-options'), feedback=$('extent-feedback');
  let count=1, revision=0, base=null, busy=false, summary='', singleExportNote='', cache=new Map();
  function clean(){return state.recipe&&state.result===base&&!state.inputsDirty&&!state.dirty&&!state.jsonPending;}
  function options(){return base?.repeat_options?.extents||[{value:1,label:'Original drawing'}];}
  function sync(){
    const ready=clean();
    $('pore-color-control').hidden=count===1||state.recipe?.layout_kind==='open_chain';
    $('color-pore-types').disabled=!ready||state.busy;
    for(const button of choices.children){
      button.setAttribute('aria-pressed',String(Number(button.value)===count));
      button.disabled=!ready||state.busy;
    }
    if(busy||!ready)for(const button of $('exports').children)button.disabled=true;
    if(state.recipe?.layout_kind==='open_chain')$('measure-pore').disabled=true;
    else $('measure-pore').disabled=!ready||state.busy||busy;
  }
  function baseView(){
    $('framework-image').src='.'+base.downloads.svg;
    $('drawing-summary').textContent=summary;studio.showExports(base.downloads);
    $('export-note').textContent=singleExportNote;$('export-note').hidden=!singleExportNote;
    feedback.textContent=base.repeat_options?.available?'':base.repeat_options?.reason||'';
    $('pore-type-legend').replaceChildren();
  }
  function showRepeat(result){
    const s=result.summary, downloads=result.downloads||{svg:result.svg};
    $('framework-image').src='.'+downloads.svg;
    $('drawing-summary').textContent=[s.extent,s.blocks+' blocks',s.atoms+' atoms',s.continuation_sites+' outward sites'].join(' · ');
    studio.showExports(downloads);$('export-note').textContent=s.cdxml_note||'';$('export-note').hidden=!s.cdxml_note;
    const legend=$('pore-type-legend');legend.replaceChildren();
    const familyNames={triangle:'Triangle',hexagon:'Hexagon',square:'Square'};
    for(const [type,n] of Object.entries(s.pore_sequences||s.pore_counts||{})){
      const face=s.template==='mixed_kgm'?$('assembly-arrangement'):s.template==='self_kagome'?$('self-face'):null;
      const canFocus=face&&Array.from(face.options).some(o=>o.value===type);
      const chip=document.createElement(canFocus?'button':'span'),color=s.pore_colors?.[type];
      if(color&&/^#[0-9a-f]{6}$/i.test(color)){
        const swatch=document.createElement('i');swatch.style.backgroundColor=color;chip.append(swatch);
      }
      const name=s.template==='mixed_kgm'?(type.length===3?'Triangle ':'Hexagon ')+type:familyNames[type]||type;
      if(canFocus){
        chip.type='button';chip.title='Show one '+name.toLowerCase();chip.setAttribute('aria-label',chip.title);
        chip.onclick=()=>{select(1);if(face.value!==type){face.value=type;face.dispatchEvent(new Event('change',{bubbles:true}));}};
      }
      chip.append(name+' × '+n);legend.append(chip);
    }
    feedback.textContent=s.source_drawing_basis?'Neighboring pores use the precursor-derived layout.':'';
  }
  async function select(value){
    if(!clean())return;
    const token=++revision, drawing=base;
    count=options().some(o=>Number(o.value)===value)?value:1;
    busy=false;bar.removeAttribute('aria-busy');
    window.dispatchEvent(new Event('cofdraw:extentchanged'));
    if(count===1){baseView();sync();return;}
    if($('edit-geometry').getAttribute('aria-pressed')==='true')$('edit-geometry').click();
    const colored=$('color-pore-types').checked&&state.recipe?.layout_kind!=='open_chain',key=count+':'+colored;
    if(cache.has(key)){showRepeat(cache.get(key));sync();return;}
    busy=true;bar.setAttribute('aria-busy','true');feedback.textContent='Drawing neighboring pores…';sync();
    try{
      const payload=state.sourceId?{source_recipe_id:state.sourceId}:{recipe:state.recipe};
      Object.assign(payload,{repeats:count,color_pores:colored,size:Number($('image-size').value),style_id:state.style,style_options:state.styleOptions});
      if(state.seed!==null)payload.style_seed=state.seed;
      const result=await studio.api('repeat',payload);
      if(token!==revision||drawing!==base||!clean())return;
      cache.set(key,result);showRepeat(result);
    }catch(error){
      if(token!==revision||drawing!==base||!clean())return;
      count=1;baseView();feedback.textContent=error.message;
    }finally{
      if(token===revision){busy=false;bar.removeAttribute('aria-busy');sync();}
    }
  }
  $('color-pore-types').addEventListener('change',()=>select(count));
  function render(event){
    revision++;busy=false;bar.removeAttribute('aria-busy');cache=new Map();base=event.detail;
    summary=$('drawing-summary').textContent;singleExportNote=$('export-note').hidden?'':$('export-note').textContent;
    choices.replaceChildren();bar.hidden=false;
    for(const o of options()){
      const button=document.createElement('button');button.type='button';button.value=o.value;button.textContent=o.label;
      button.addEventListener('click',()=>select(Number(o.value)));choices.append(button);
    }
    select(count);
  }
  // Geometry edits and a diameter belong to the chosen local pore.
  for(const id of ['edit-geometry','clean-geometry','measure-pore','undo-geometry']){
    $(id).addEventListener('click',e=>{if(e.isTrusted&&count>1)select(1);},{capture:true});
  }
  window.addEventListener('cofdraw:rendered',render);
  window.addEventListener('cofdraw:layoutchanged',()=>{revision++;busy=false;cache.clear();base=null;bar.hidden=true;});
  window.addEventListener('cofdraw:state',sync);
  window.COFdrawExtent={get value(){return count;},restore(value,color=true){count=[1,2,3,4].includes(value)?value:1;$('color-pore-types').checked=color!==false;}};
});
