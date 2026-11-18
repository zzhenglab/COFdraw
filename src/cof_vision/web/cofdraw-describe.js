/* Description -> reviewed catalog choices -> the existing chemical renderer. */
"use strict";
(() => {
  const S=window.COFdrawStudio, $=id=>document.getElementById(id);
  const state={catalog:null,plan:null,selection:null,selectionEdited:false,busy:false,revision:0,startMode:'describe',context:null,history:[],activeText:'',activeContext:null,applied:false};
  let currentSnapshot=null;
  const fields={source_recipe_id:'describe-source',monomer_a_id:'describe-monomer-a',monomer_b_id:'describe-monomer-b',linkage:'describe-linkage',template_id:'describe-template',style_id:'describe-style',palette:'describe-palette'};
  const clone=value=>JSON.parse(JSON.stringify(value));
  const previews=new Map();
  function message(value,error=false){const el=$('describe-status');el.textContent=value;el.classList.toggle('error',error);}
  function readable(value){if(typeof value==='string')return value;if(value===null||value===undefined)return '';if(typeof value==='object')return value.message||value.reason||value.summary||value.label||value.name||JSON.stringify(value);return String(value);}
  function currentValid(){return Boolean(currentSnapshot&&S.state.recipe&&!S.state.inputsDirty&&!S.state.jsonPending&&JSON.stringify(S.state.recipe)===currentSnapshot.recipe&&S.state.result===currentSnapshot.result&&S.state.style===currentSnapshot.style&&JSON.stringify(S.state.styleOptions)===currentSnapshot.options);}
  function captureCurrent(selection={}){
    const token=crypto.randomUUID().replaceAll('-','');currentSnapshot={recipe:JSON.stringify(S.state.recipe),result:S.state.result,style:S.state.style,options:JSON.stringify(S.state.styleOptions)};
    const selected={};for(const key of Object.keys(fields))if(typeof selection[key]==='string'&&selection[key])selected[key]=selection[key];selected.style_id=S.state.style;if(S.state.styleOptions.palette)selected.palette=S.state.styleOptions.palette;
    state.context={origin:'current_drawing',token,selection:selected,pending:state.context?.pending||[],requirements:state.context?.requirements||[]};
  }
  function contextLabel(){
    if(state.context?.origin==='current_drawing')return 'Context: the current drawing. Appearance changes preserve its graph and coordinates; chemistry changes require a new reviewed design.';
    const selection=state.context?.selection||{},source=state.catalog?.sources.find(row=>row.id===selection.source_recipe_id);
    if(source)return 'Context: '+source.label+' (qualified source drawing).';
    const names=['monomer_a_id','monomer_b_id'].map(key=>state.catalog?.monomers.find(row=>row.id===selection[key])?.abbreviation).filter(Boolean);
    return names.length?'Context: '+names.join(' + ')+'. Changes are checked before applying.':'Choose a source or two monomers. You can give one choice at a time.';
  }
  function openAction(id){
    if(id==='examples'){showInputs(true);$('describe-demos').scrollIntoView({block:'nearest'});return;}
    if(id==='monomer_library'){startMode('manual');$('open-monomer-library').click();return;}
    const tabs={science_model:'model',science_porosity:'porosity',science_adsorption:'adsorption',science_bet:'bet'};
    if(tabs[id]){$('open-science').click();document.querySelector(`[data-science-tab="${tabs[id]}"]`).click();}
  }
  function drawConversation(plan=state.plan){
    const history=$('chat-history');history.replaceChildren();
    for(const turn of state.history){const card=document.createElement('article'),who=document.createElement('strong'),text=document.createElement('p'),summary=document.createElement('p'),outcome=document.createElement('span');card.className='chat-turn';who.textContent='You';text.textContent=turn.text;summary.textContent=turn.summary;outcome.className='chat-outcome';outcome.textContent=turn.applied?'Rendered successfully':turn.status==='ready'?'Reviewed plan · not applied':turn.status==='unsupported'?'Unsupported request retained':'Needs clarification';card.append(who,text,summary,outcome);history.append(card);}
    $('chat-composer').hidden=!state.context;$('chat-context').textContent=contextLabel();
    const actions=$('chat-actions');actions.replaceChildren();for(const action of plan?.actions||[]){const button=document.createElement('button');button.type='button';button.textContent=action.label;button.addEventListener('click',()=>openAction(action.id));actions.append(button);}
    const suggestions=$('chat-suggestions');suggestions.replaceChildren();for(const text of (plan?.conversation?.suggestions||['Make it pastel','Use wavy bonds','How do I check SA scores?']).slice(0,3)){const button=document.createElement('button');button.type='button';button.textContent=text;button.addEventListener('click',()=>{$('chat-text').value=text;$('chat-text').focus();controls();});suggestions.append(button);}
  }
  function newConversation(){if(state.busy||S.state.busy)return;state.history=[];state.context=null;state.activeText='';state.activeContext=null;state.applied=false;currentSnapshot=null;$('describe-text').value='';$('chat-text').value='';invalidate({reset:true});drawConversation();message('New conversation. Choose an example or describe your drawing.');$('describe-text').focus();}
  $('chat-new').addEventListener('click',newConversation);
  $('chat-use-current').addEventListener('click',()=>{
    if(!S.state.recipe||S.state.inputsDirty||S.state.jsonPending||state.busy||S.state.busy)return;
    if(state.plan&&state.activeText)state.history.push({text:state.activeText,summary:state.plan.summary,status:state.plan.status,applied:state.applied});state.history=state.history.slice(-8);
    state.revision++;state.plan=null;state.selectionEdited=false;captureCurrent(S.state.sourceId?{source_recipe_id:S.state.sourceId}:{});state.activeContext=clone(state.context);state.activeText='Use current drawing';$('describe-review').hidden=true;showInputs(false);$('describe-current-text').textContent=state.activeText;drawConversation();message('Current drawing selected. Ask for an appearance change, then review and apply it.');$('chat-text').focus();controls();
  });
  $('chat-text').addEventListener('input',controls);
  $('chat-send').addEventListener('click',()=>analyze({followup:true}));
  $('chat-text').addEventListener('keydown',event=>{if((event.ctrlKey||event.metaKey)&&event.key==='Enter'){event.preventDefault();analyze({followup:true});}});
  function controls(){
    const blocked=state.busy||S.state.busy;
    if(!blocked&&state.plan?.operation==='style_only'&&!currentValid()){state.plan=null;state.revision++;$('describe-plan-state').textContent='Drawing changed';message('The current drawing changed. Choose Use current drawing again before applying a style.',true);}
    $('describe-panel').setAttribute('aria-busy',String(state.busy));
    $('describe-analyze').disabled=blocked||(!state.catalog?false:!S.state.online||!$('describe-text').value.trim());
    $('describe-create').disabled=blocked||!S.state.online||state.plan?.status!=='ready'||state.applied;
    $('chat-new').disabled=blocked;$('chat-use-current').disabled=blocked||!S.state.recipe||S.state.inputsDirty||S.state.jsonPending;
    $('chat-text').disabled=blocked;$('chat-send').disabled=blocked||!state.context||!$('chat-text').value.trim();
    $('chat-send').textContent=$('describe-mode').value==='llm'?'Review with model API':'Review follow-up';
    document.querySelectorAll('#chat-suggestions button,#chat-actions button').forEach(button=>button.disabled=blocked);
    for(const id of ['describe-text','describe-mode','describe-base-url','describe-model','describe-api-key','describe-edit-text',...Object.values(fields),'mode-describe','mode-manual'])$(id).disabled=blocked;
    $('describe-demos').querySelectorAll('button').forEach(button=>button.disabled=blocked);
    for(const id of ['edit-geometry','clean-geometry','measure-pore'])$(id).disabled=blocked||!S.state.recipe;
  }
  function startMode(mode,{focus=false}={}){
    if(state.busy||S.state.busy)return;
    state.startMode=mode==='manual'?'manual':'describe';document.body.dataset.startMode=state.startMode;
    for(const name of ['describe','manual']){const active=name===state.startMode;$(name+'-panel').hidden=!active;const button=$('mode-'+name);button.setAttribute('aria-selected',String(active));button.tabIndex=active?0:-1;}
    if(focus)$(state.startMode==='describe'?'describe-text':'open-monomer-library').focus();
  }
  for(const name of ['describe','manual'])$('mode-'+name).addEventListener('click',()=>startMode(name));
  document.querySelector('.start-modes').addEventListener('keydown',event=>{
    if(!['ArrowLeft','ArrowRight','Home','End'].includes(event.key))return;
    event.preventDefault();const mode=event.key==='Home'?'describe':event.key==='End'?'manual':state.startMode==='describe'?'manual':'describe';startMode(mode);$('mode-'+mode).focus();
  });
  function showInputs(show){$('describe-inputs').hidden=!show;$('describe-edit-text').hidden=show;$('describe-current-text').textContent=$('describe-text').value.trim();provider();}
  $('describe-edit-text').addEventListener('click',()=>{showInputs(true);$('describe-text').focus();});
  function invalidate({reset=false}={}){
    state.revision++;state.plan=null;$('describe-create').disabled=true;
    $('describe-monomer-previews').replaceChildren();$('describe-monomer-previews').hidden=true;
    $('describe-plan-state').textContent='Review needed';$('describe-plan-state').className='plan-state';
    if(reset){state.selection=null;state.selectionEdited=false;state.context=null;currentSnapshot=null;delete document.body.dataset.chatStage;$('chat-composer').hidden=true;$('describe-review').hidden=true;showInputs(true);}
    else{$('describe-summary').textContent='Your choices changed. Analyze the description again to check the complete plan.';$('describe-feedback').replaceChildren();$('describe-evidence').hidden=true;}
    if(state.catalog)message('Analyze the current description and choices before creating a drawing.');
    provider();controls();
  }
  $('describe-text').addEventListener('input',()=>invalidate({reset:true}));
  function provider(){
    const remote=$('describe-mode').value==='llm';$('describe-llm-settings').hidden=!remote;$('describe-local-note').hidden=remote;$('describe-mode-label').textContent=remote?'Model API':'Local rules';
    $('describe-analyze').textContent=state.selectionEdited?'Check revised choices':remote?'Analyze with model API':'Analyze description';
    $('describe-analyze').hidden=Boolean(!state.selectionEdited&&$('describe-inputs').hidden);
  }
  $('describe-mode').addEventListener('change',()=>{provider();invalidate({reset:true});});
  for(const id of ['describe-base-url','describe-model','describe-api-key'])$(id).addEventListener('input',()=>invalidate({reset:true}));
  function fill(id,rows,empty){const select=$(id);select.replaceChildren();const placeholder=document.createElement('option');placeholder.value='';placeholder.textContent=empty;select.append(placeholder);for(const row of rows||[]){const option=document.createElement('option');option.value=row.id||row.recipe_id;option.textContent=row.label||row.full_name||row.abbreviation||option.value;if(row.supported===false)option.textContent+=' · reference only';select.append(option);}}
  function selectValue(id,value){const select=$(id),text=value||'';select.value=Array.from(select.options).some(option=>option.value===text)?text:'';}
  function sourceContext(){
    const id=$('describe-source').value,source=state.catalog?.sources?.find(row=>(row.id||row.recipe_id)===id);
    $('describe-design-fields').hidden=Boolean(id);$('describe-source-note').hidden=!id;
    if(source){const names=(source.monomer_ids||[]).map(key=>state.catalog.monomers.find(row=>row.id===key)?.label||key);$('describe-source-note').textContent=[source.label,names.length?'Associated precursors: '+names.join(' + '):'',source.scope||'This source uses its recorded connections and layout.'].filter(Boolean).join(' · ');}
    else $('describe-source-note').textContent=id?'This source uses its recorded connections and layout.':'';
  }
  function selectionFromFields(){const selection={};for(const[key,id]of Object.entries(fields))selection[key]=$(id).value||null;return selection;}
  function choiceSummary(){
    const parent=$('describe-choice-summary');parent.replaceChildren();
    const source=$('describe-source').value;
    const choices=state.plan?.operation==='style_only'?[]:source?[['Source drawing','describe-source']]:[['Monomer A','describe-monomer-a'],['Monomer B','describe-monomer-b'],['Linkage','describe-linkage'],['Template','describe-template']];
    choices.push(['Style','describe-style'],['Palette','describe-palette']);
    for(const[label,id]of choices){const dt=document.createElement('dt'),dd=document.createElement('dd');dt.textContent=label;dd.textContent=$(id).selectedOptions[0]?.textContent||'Not selected';parent.append(dt,dd);}
  }
  function monomerPreviews(plan){
    const parent=$('describe-monomer-previews');parent.replaceChildren();
    const preset=id=>S.state.library?.monomer_presets?.find(row=>row.id===id);
    let rows=[];
    if(plan.prepare_request)rows=['a','b'].map((block,index)=>({smiles:plan.prepare_request[`monomer_${block}`],label:plan.prepare_request[`label_${block}`]||String.fromCharCode(65+index),caption:`Monomer ${block.toUpperCase()}`}));
    else{const source=state.catalog.sources.find(row=>(row.id||row.recipe_id)===plan.source_recipe_id);rows=(source?.monomer_ids||[]).map(id=>preset(id)).filter(Boolean).map(row=>({smiles:row.smiles,label:row.abbreviation||row.label||row.id,caption:'Associated precursor'}));}
    rows=rows.filter(row=>row.smiles);parent.hidden=!rows.length;const current=state.plan;
    for(const row of rows){const card=document.createElement('div'),caption=document.createElement('span'),waiting=document.createElement('span');card.className='describe-monomer-card';caption.className='describe-monomer-caption';caption.textContent=row.label;caption.title=row.caption;waiting.className='describe-preview-state';waiting.textContent='Loading structure…';card.append(caption,waiting);parent.append(card);
      if(!previews.has(row.smiles))previews.set(row.smiles,S.api('monomer',{smiles:row.smiles}).then(result=>URL.createObjectURL(new Blob([result.svg],{type:'image/svg+xml'}))).catch(error=>{previews.delete(row.smiles);throw error;}));
      previews.get(row.smiles).then(url=>{if(state.plan!==current||!card.isConnected)return;const image=document.createElement('img');image.alt=`${row.caption}: ${row.label}`;image.src=url;waiting.replaceWith(image);}).catch(()=>{if(state.plan===current&&card.isConnected)waiting.textContent='Preview unavailable';});
    }
  }
  for(const[key,id]of Object.entries(fields))$(id).addEventListener('change',()=>{
    if(key==='source_recipe_id'&&$(id).value)for(const name of ['monomer_a_id','monomer_b_id','linkage','template_id'])selectValue(fields[name],null);
    if(['monomer_a_id','monomer_b_id','linkage','template_id'].includes(key))selectValue('describe-source',null);
    sourceContext();choiceSummary();state.selection=selectionFromFields();state.selectionEdited=true;
    invalidate();
  });
  function list(parent,title,items,warning=false){
    const rows=Array.isArray(items)?items:items?[items]:[];if(!rows.length)return;
    const box=document.createElement('div');box.className=warning?'describe-notes warning':'describe-notes';const heading=document.createElement('strong');heading.textContent=title;box.append(heading);const ul=document.createElement('ul');for(const value of rows){const li=document.createElement('li');li.textContent=readable(value);ul.append(li);}box.append(ul);parent.append(box);
  }
  function showPlan(plan){
    if(!['ready','needs_input','unsupported'].includes(plan?.status))throw new Error('The description service returned an incomplete plan. Please try again.');
    if(plan.status==='ready'&&plan.operation!=='style_only'&&Boolean(plan.source_recipe_id)===Boolean(plan.prepare_request))throw new Error('The plan does not identify one drawing workflow. Please revise the description.');
    state.plan=clone(plan);state.selection=clone(plan.selection||{});state.selectionEdited=false;document.body.dataset.chatStage='review';provider();
    state.context=clone(plan.conversation?.context||{selection:plan.selection||{},pending:[],origin:'plan'});state.applied=false;
    for(const[key,id]of Object.entries(fields))selectValue(id,state.selection[key]??(key==='style_id'?plan.style_id:null));
    sourceContext();choiceSummary();showInputs(false);monomerPreviews(plan);$('describe-adjust').open=plan.status!=='ready';$('describe-review').hidden=false;$('describe-summary').textContent=plan.summary||'Review the catalog choices below.';
    $('describe-current-text').textContent=state.activeText;$('describe-create').textContent=plan.operation==='style_only'?'Apply appearance':'Create drawing';
    $('describe-create').hidden=['help','reset'].includes(plan.operation);$('describe-review-title').textContent=plan.operation==='help'?'Local guidance':'Review the drawing';
    $('describe-adjust').hidden=['help','reset','style_only'].includes(plan.operation);
    $('describe-choice-summary').hidden=['help','reset'].includes(plan.operation);
    const badge=$('describe-plan-state');badge.textContent=plan.operation==='help'?'Tool guide':plan.status==='ready'?'Ready to draw':plan.status==='unsupported'?'Not supported':'Needs your input';badge.className='plan-state '+(plan.status==='ready'?'ready':'');
    $('describe-scope').textContent=plan.scope||'This creates a finite chemical drawing. It does not assign a crystal structure or validate synthesis.';
    const feedback=$('describe-feedback');feedback.replaceChildren();list(feedback,'Resolve before drawing',plan.unresolved,true);list(feedback,'Please review',plan.warnings,true);list(feedback,'Available alternatives',plan.alternatives);
    const evidence=$('describe-evidence-content');evidence.replaceChildren();list(evidence,'Matched catalog records',Array.isArray(plan.matches)?plan.matches:Object.entries(plan.matches||{}).map(([key,value])=>`${key.replaceAll('_',' ')}: ${readable(value)}`));list(evidence,'Model interpretation',plan.interpretation?.explanation);list(evidence,'Assumptions',plan.assumptions);$('describe-evidence').hidden=!evidence.childElementCount;
    message(plan.status==='ready'?'The plan is ready. Review the exact choices, then create your drawing.':plan.status==='unsupported'?'This request is outside the supported drawing workflows. Review the alternatives.':'Some choices need clarification. Edit the description or selections, then analyze again.');
    if(plan.operation==='help')message('Local guidance only. Open a tool to supply inputs and run a calculation.');
    if(plan.operation==='reset'){state.context=null;currentSnapshot=null;state.history=[];$('describe-text').value='';showInputs(true);$('describe-review').hidden=true;}
    drawConversation(plan);
  }
  async function analyze(options={}){
    if(state.busy||S.state.busy)return;if(!state.catalog){await loadCatalog();return;}
    const followup=options.followup===true,text=(state.selectionEdited?state.activeText:followup?$('chat-text').value:$('describe-text').value).trim();if(!text)return;
    const context=state.selectionEdited?state.activeContext:followup?state.context:null;
    const secret=$('describe-api-key').value;if(secret&&text.includes(secret)){message('Keep API credentials in the password field, separate from chat messages.',true);return;}
    if(context?.origin==='current_drawing'&&!currentValid()){message('The drawing changed. Choose Use current drawing again.',true);state.plan=null;controls();return;}
    const revision=state.revision,mode=state.selectionEdited?'rules':$('describe-mode').value;
    const payload={text,mode,selection:state.selectionEdited&&state.selection?clone(state.selection):null,context:context?clone(context):null,llm:null};
    if(mode==='llm'){
      const base=$('describe-base-url').value.trim(),model=$('describe-model').value.trim();
      if(!base||!model){message('Enter your provider API base URL and model name first.',true);$('describe-provider').open=true;(!base?$('describe-base-url'):$('describe-model')).focus();return;}
      payload.llm={base_url:base,model,api_key:$('describe-api-key').value};
    }
    const previous=state.plan?{text:state.activeText,summary:state.plan.summary,status:state.plan.status,applied:state.applied}:null;
    state.plan=null;state.busy=true;$('describe-review').hidden=true;message(mode==='llm'?'Sending your new message and catalog IDs to the selected provider…':'Checking the message against the local catalog…');controls();
    try{const plan=await S.api('describe/chat',payload);if(revision!==state.revision)return;if(context?.origin==='current_drawing'&&!currentValid())throw new Error('The drawing changed while the request was running. Choose Use current drawing again.');if(previous&&previous.text)state.history.push(previous);state.history=state.history.slice(-8);state.activeText=text;state.activeContext=context?clone(context):null;$('chat-text').value='';showPlan(plan);}
    catch(error){state.plan=null;$('describe-review').hidden=true;showInputs(true);message(error.message||'The description could not be analyzed.',true);}
    finally{payload.llm=null;state.busy=false;controls();}
  }
  $('describe-analyze').addEventListener('click',analyze);
  $('describe-text').addEventListener('keydown',event=>{if((event.ctrlKey||event.metaKey)&&event.key==='Enter'){event.preventDefault();analyze();}});
  async function create(){
    if(state.busy||S.state.busy||state.plan?.status!=='ready')return;
    const plan=clone(state.plan);state.busy=true;message('Building the checked chemical graph and preparing exports…');controls();
    try{let result;if(plan.operation==='style_only'){if(!currentValid()||plan.context_token!==state.context?.token)throw new Error('The drawing changed. Review a fresh current-drawing request.');S.state.style=plan.style_id;S.state.styleOptions={...S.state.styleOptions,...clone(plan.style_options||{})};S.state.resolvedStyleOptions={};S.drawStyleControls();result=await S.render();}else result=await S.applyDescriptionPlan(plan);if(!result?.downloads?.svg)throw new Error('No completed drawing was returned.');state.applied=true;captureCurrent(state.context?.selection||plan.selection);message('Drawing rendered successfully. Inspect the canvas, then save the image and editable files.');$('describe-plan-state').textContent='Applied';drawConversation(plan);if(window.innerWidth<760)$('canvas-title').scrollIntoView({behavior:matchMedia('(prefers-reduced-motion: reduce)').matches?'auto':'smooth',block:'start'});}
    catch(error){state.plan=null;$('describe-plan-state').textContent='Check the plan';$('describe-plan-state').className='plan-state';message(error.message||'The drawing could not be created.',true);}
    finally{state.busy=false;controls();}
  }
  $('describe-create').addEventListener('click',create);
  async function loadCatalog(){
    state.busy=true;controls();
    try{
      const catalog=await S.api('describe/catalog');if(!Array.isArray(catalog.demos)||!Array.isArray(catalog.monomers))throw new Error('Description choices are unavailable. Retry when the local engine is ready.');state.catalog=catalog;
      fill('describe-source',catalog.sources,'New design from monomers');fill('describe-monomer-a',catalog.monomers,'Choose monomer A');fill('describe-monomer-b',catalog.monomers,'Choose monomer B');fill('describe-linkage',catalog.linkages,'Auto — from monomers and description');fill('describe-template',catalog.templates,'Auto — from reactive sites and description');fill('describe-style',catalog.styles,'Default style');fill('describe-palette',catalog.palettes,'Style default');
      const demos=$('describe-demos');demos.replaceChildren();for(const demo of catalog.demos){const button=document.createElement('button');button.type='button';button.textContent=demo.label;button.dataset.demo=demo.id;button.addEventListener('click',()=>{$('describe-text').value=demo.text;invalidate({reset:true});$('describe-text').focus();});demos.append(button);}
      provider();message('Local rules are ready. Describe your idea or choose an example.');
    }catch(error){state.catalog=null;message(error.message,true);$('describe-analyze').textContent='Retry description tools';$('describe-demos').replaceChildren();}
    finally{state.busy=false;controls();}
  }
  window.addEventListener('cofdraw:state',controls);window.addEventListener('cofdraw:library',controls);
  window.COFdrawDescribe={state,analyze,create,startMode,newConversation,reloadCatalog:loadCatalog};
  provider();loadCatalog();
})();
