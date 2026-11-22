/* Show each drawing extent on the same canvas, with exports from the visible graph. */
document.addEventListener('DOMContentLoaded',()=>{
  const studio=window.COFdrawStudio, state=studio.state, $=id=>document.getElementById(id);
  const bar=$('extent-bar'), choices=$('extent-options'), feedback=$('extent-feedback');
  let count=1, revision=0, base=null, visible=null, busy=false, summary='', singleExportNote='', cache=new Map();
  let crystalPick=null;
  const inCrystal=()=>Boolean(window.COFdrawCrystalDesign?.active);
  function clean(){return !inCrystal()&&state.recipe&&state.result===base&&!state.inputsDirty&&!state.dirty&&!state.jsonPending;}
  function options(){return base?.repeat_options?.extents||[{value:1,label:'Original drawing'}];}
  function populate(rows,selected){
    choices.replaceChildren(...rows.map(row=>new Option(row.label,String(row.value))));
    choices.value=String(selected);
    bar.hidden=!rows.length;
  }
  function setCrystalOptions(rows,selected,onChange){
    if(!inCrystal())return;
    revision++;busy=false;crystalPick=onChange;visible=null;bar.removeAttribute('aria-busy');
    populate(rows,selected);
    choices.disabled=!onChange||rows.length<2;
    $('pore-color-control').hidden=true;$('pore-type-legend').replaceChildren();feedback.textContent='';
  }
  function sync(){
    if(inCrystal())return;
    if(crystalPick){crystalPick=null;populate(base?options():[],count);}
    const ready=clean();
    $('pore-color-control').hidden=!ready||count===1||state.recipe?.layout_kind==='open_chain';
    $('color-pore-types').disabled=!ready||state.busy;
    choices.value=String(count);choices.disabled=!ready||state.busy||choices.options.length<2;
    if(busy||!ready)for(const button of $('exports').children)button.disabled=true;
    if(state.recipe?.layout_kind==='open_chain')$('measure-pore').disabled=true;
    else $('measure-pore').disabled=!ready||state.busy||busy||(count!==1&&!visible?.summary?.measurement_face);
    $('measure-pore').title=count===1?'Measure the displayed pore':visible?.summary?.measurement_face?
      'Measure focal pore · '+visible.summary.measurement_face.label:visible?.summary?.measurement_reason||'Waiting for the displayed pore';
  }
  function baseView(){
    $('framework-image').src='.'+base.downloads.svg;
    $('drawing-summary').textContent=summary;studio.showExports(base.downloads);
    $('export-note').textContent=singleExportNote;$('export-note').hidden=!singleExportNote;
    feedback.textContent=base.repeat_options?.available?'':base.repeat_options?.reason||'';
    $('pore-type-legend').replaceChildren();
    visible=base;window.dispatchEvent(new CustomEvent('cofdraw:viewready',{detail:visible}));
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
    visible=result;window.dispatchEvent(new CustomEvent('cofdraw:viewready',{detail:visible}));
  }
  async function select(value){
    if(!clean())return;
    const token=++revision, drawing=base;
    const chosen=options().find(o=>o.value===value)||options()[0];
    count=chosen.value;visible=null;
    busy=false;bar.removeAttribute('aria-busy');
    window.dispatchEvent(new Event('cofdraw:extentchanged'));
    if(count===1){baseView();sync();return;}
    if($('edit-geometry').getAttribute('aria-pressed')==='true')$('edit-geometry').click();
    const colored=$('color-pore-types').checked&&state.recipe?.layout_kind!=='open_chain',key=count+':'+colored;
    if(cache.has(key)){showRepeat(cache.get(key));sync();return;}
    busy=true;bar.setAttribute('aria-busy','true');feedback.textContent='Drawing neighboring pores…';sync();
    try{
      const payload=state.sourceId?{source_recipe_id:state.sourceId}:{recipe:state.recipe};
      Object.assign(payload,{repeats:chosen.repeats||count,extent_mode:chosen.extent_mode||'pores',color_pores:colored,size:Number($('image-size').value),style_id:state.style,style_options:state.styleOptions});
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
  choices.addEventListener('change',()=>{
    if(inCrystal()){crystalPick?.(choices.value);return;}
    const chosen=options().find(row=>String(row.value)===choices.value);
    if(chosen)select(chosen.value);
  });
  function render(event){
    if(inCrystal())return;
    crystalPick=null;
    revision++;busy=false;bar.removeAttribute('aria-busy');cache=new Map();base=event.detail;
    summary=$('drawing-summary').textContent;singleExportNote=$('export-note').hidden?'':$('export-note').textContent;
    populate(options(),count);
    select(count);
  }
  // Geometry edits apply to the editable local pore. Measurements use the displayed extent.
  for(const id of ['edit-geometry','clean-geometry','undo-geometry']){
    $(id).addEventListener('click',()=>{
      const closingEditor=id==='edit-geometry'&&$(id).getAttribute('aria-pressed')==='true';
      if(!inCrystal()&&count!==1&&!closingEditor)select(1);
    },{capture:true});
  }
  window.addEventListener('cofdraw:rendered',render);
  window.addEventListener('cofdraw:layoutchanged',()=>{
    if(inCrystal())return;
    revision++;busy=false;cache.clear();base=null;visible=null;crystalPick=null;bar.hidden=true;
    $('pore-color-control').hidden=true;$('pore-type-legend').replaceChildren();feedback.textContent='';
  });
  window.addEventListener('cofdraw:crystalviewclear',()=>{if(inCrystal())setCrystalOptions([],null,null);});
  window.addEventListener('cofdraw:state',sync);
  window.COFdrawExtent={get value(){return count;},get expanded(){return !inCrystal()&&count!==1;},get visible(){return clean()?visible:null;},
    select,setCrystalOptions,restore(value,color=true){count=[1,'three-pores','kagome-star','triangle',2,3,4].includes(value)?value:1;$('color-pore-types').checked=color!==false;}};
});
