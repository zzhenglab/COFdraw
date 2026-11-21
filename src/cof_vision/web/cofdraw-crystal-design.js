/* The 2D drawing and the periodic viewer share one selected crystal report. */
(()=>{
  'use strict';
  const $=id=>document.getElementById(id), studio=window.COFdrawStudio, lab=window.COFdrawScience;
  const main=document.querySelector('main.studio'), picker=$('crystal-linker'), picture=$('framework-image');
  const entries=new Map();
  const matched={
    'cof-300-c5':['tapm','ta'], 'curated-13040N3':['tapm','bpda_cho'], 'dia-imine':['tapm','ta']
  };
  const custom='custom-dia';
  let id=null, edit=null, revision=0, referenceRevision=0, report=null, visible=null, timer, failure='';
  let preparedCells=null, preferredView='Framework drawing',count=1,job=0,patchStyle=null,scope='';
  const cache=new Map();
  const styles=[
    ['reticular_depth','Depth fade · ink and ring wash','layered','Farther bonds fade with their depth in the crystal.'],
    ['reticular_blocks','Depth fade · two building blocks','layered','Blue and amber distinguish the reacted precursors.'],
    ['standard','Standard chemical drawing','chemical','Conventional chemical bonds and continuation marks.'],
    ['red_blue','Colored building blocks','chemical','Trace each assigned precursor with its own ink.'],
    ['ring_fill','Ring highlights','chemical','Shade molecular rings beneath the bond-line drawing.'],
    ['bold_color','Bold chemical drawing','chemical','Emphasize bonds and heteroatom labels.']
  ].map(([id,label,category,description])=>({id,label,category,description,preview_kind:id}));
  let restoredReference=false, comparing=false;
  function availableStyles(){
    return restoredReference?styles.filter(s=>!['reticular_blocks','red_blue'].includes(s.id)):styles;
  }
  function clearImage(){
    if(visible?.url)URL.revokeObjectURL(visible.url);
    visible=null;picture.hidden=true;picture.removeAttribute('src');$('exports').replaceChildren();
    window.dispatchEvent(new Event('cofdraw:crystalviewclear'));
  }
  function clearPicture(){
    job++;report=null;preparedCells=null;cache.clear();comparing=false;clearImage();
    $('crystal-drawing-view').replaceChildren();$('representation-note').hidden=true;
  }
  function leave(){
    if(!id)return;
    clearTimeout(timer);revision++;id=null;edit=null;failure='';clearPicture();main.classList.remove('crystal-design-active');
    $('crystal-design-controls').hidden=true;$('crystal-view-field').hidden=true;
    if(patchStyle){Object.assign(studio.state,patchStyle);patchStyle=null;}
    $('canvas-title').textContent='Framework canvas';$('canvas-sheet').classList.remove('transparent-preview');
    $('extent-bar').hidden=true;$('canvas-empty').classList.remove('hidden');studio.refreshStyles();
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
    if(!report){lab.status(failure||$('drawing-summary').textContent||'Waiting for the selected crystal coordinates.');return;}
    if(lab.structure!==report.structure||lab.origin!=='catalog-design')
      lab.showCatalogDesign(report,restoredReference?id:custom,preparedCells);
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
    const compatible=mine.length>0&&(!partner||mine.some(a=>other.some(b=>a.role!==b.role&&
      ((a.n===4&&item.geometry_hint==='tetrahedral'&&b.n===2)||(b.n===4&&partner.geometry_hint==='tetrahedral'&&a.n===2)
      ||a.n===4&&b.n===4&&[item,partner].every(m=>m.geometry_hint==='tetrahedral')))));
    return {compatible,reason:compatible?'Reactive site counts match. The cell builder checks tetrahedral and linker geometry.':
      'Pair a four-arm tetrahedral node with a complementary two-ended spacer or another tetrahedral node.'};
  }
  function styleSVG(text){
    const svg=new DOMParser().parseFromString(text,'image/svg+xml').documentElement;
    const {style,styleOptions:opt,library}=studio.state;
    const colors=opt.custom_palette||library?.style_controls?.palette_colors?.[opt.palette]||['#2768b0','#bd741e'];
    const depthFade=style.startsWith('reticular_'), wash=depthFade||style==='ring_fill';
    const blockInk=['reticular_blocks','red_blue','bold_color'].includes(style);
    const elementInk={jewel:{N:'#5944a6',O:'#c64059',F:'#10857d',S:'#ab7213'},conventional:{N:'#245bd1',O:'#d82c2c',F:'#268a4b',S:'#b18d0c'},deep_ink:{N:'#244487',O:'#ad3540',F:'#207458',S:'#946c09'}}[opt.element_scheme]||{};
    if(wash)svg.querySelectorAll('[data-ring-washes]').forEach(g=>g.removeAttribute('display'));
    const parentWidth=Number(svg.querySelector('g[stroke-width]')?.getAttribute('stroke-width'))||1.65;
    for(const mark of svg.querySelectorAll('[data-depth]')){
      const block=mark.dataset.block, ink=block==='A'?colors[0]:block==='B'?colors[1]||colors[0]:null;
      const visible=depthFade?.14+.86*Number(mark.dataset.depth)**.85:1;
      if(mark.localName==='polygon'){
        mark.setAttribute('fill',style==='reticular_depth'?'#abbad8':ink||'#94aaca');
        mark.setAttribute('opacity',((opt.fill_opacity??.35)*visible).toFixed(3));
      }else if(mark.localName==='g'){
        mark.setAttribute('opacity',visible.toFixed(3));
        if(blockInk&&ink)mark.setAttribute('stroke',ink);
        mark.setAttribute('stroke-width',parentWidth*(opt.bond_scale||1)*(style==='bold_color'?1.65:1));
      }else if(mark.localName==='path'){
        mark.setAttribute('stroke-opacity',visible.toFixed(3));
        mark.setAttribute('stroke-width',Number(mark.getAttribute('stroke-width'))*(opt.bond_scale||1)*(style==='bold_color'?1.65:1));
        if(blockInk&&ink)mark.setAttribute('stroke',ink);
        else if(elementInk[mark.dataset.element])mark.setAttribute('stroke',elementInk[mark.dataset.element]);
      }else if(mark.localName==='text'){
        mark.setAttribute('opacity',visible.toFixed(3));
        if(opt.font_family)mark.setAttribute('font-family',opt.font_family);
        if(opt.font_scale)mark.setAttribute('font-size',opt.font_scale*Number(mark.getAttribute('font-size')||mark.parentNode.getAttribute('font-size')||14));
        if(opt.hydrogen_labels==='compact')mark.textContent=mark.textContent.replace(/^(N|O|S)H$/,'$1');
        const color=elementInk[mark.textContent[0]];
        if(color)mark.setAttribute('fill',color);
      }
    }
    svg.setAttribute('data-circuit-count',count);svg.setAttribute('data-drawing-style',style);
    if(opt.transparent_background){
      svg.querySelectorAll('rect[data-paper]').forEach(r=>r.remove());
    }else if(!svg.querySelector('rect[data-paper]')){
      const bg=document.createElementNS(svg.namespaceURI,'rect');bg.setAttribute('width','100%');bg.setAttribute('height','100%');bg.setAttribute('fill','white');svg.prepend(bg);
    }
    return new XMLSerializer().serializeToString(svg);
  }
  function circuitOptions(){
    const type=$('crystal-drawing-view').value;
    const traced=['Framework drawing','Projected circuit'].includes(type),counts=traced?(report.circuit_extents||[1]):[1];
    window.COFdrawExtent?.setCrystalOptions(counts.map(n=>({value:String(n),label:traced?n+(n===1?' circuit':' connected circuits'):type==='Linker detail'?'One fragment':'Displayed cells'})),String(count),n=>drawView(Number(n)));
  }
  async function drawView(n=count){
    if(!report)return;
    const token=++job, model=report, type=$('crystal-drawing-view').value;
    count=['Framework drawing','Projected circuit'].includes(type)&&model.circuit_extents?.includes(n)?n:1;
    preferredView=type;clearImage();$('canvas-empty').classList.add('hidden');circuitOptions();
    let source=cache.get(count),text='',caption='';
    try{
      if(!source){
        $('drawing-summary').textContent='Tracing neighboring circuits…';
        source=studio.api('science/reports/'+model.report_id+'/circuits/'+count,{});
        cache.set(count,source);
      }
      source=await source;
      if(token!==job||model!==report)return;
      if(type==='Framework drawing'){
        text=source.scheme_svg||source.projection_svg;caption=source.scheme_svg?'Bond-line geometry · crystal depth':'Cartesian projection';
      }else if(type==='Projected circuit'){text=source.projection_svg;caption='Cartesian projection · actual atomic coordinates';}
      else if(type==='Chemical drawing + periodic network'){text=window.COFdrawChemicalNetworkFigure(model,preparedCells);caption='Bond-line drawing and the corresponding repeated crystal';}
      else if(type==='Periodic packing'){text=window.COFdrawChemicalNetworkFigure(model,preparedCells,true);caption='Periodic bonds in the selected crystal';}
      else{text=model.chemical_drawing_svg;caption='Local linker connectivity';}
      if(!text)throw Error('No chemical drawing is available for this view.');
      text=styleSVG(text);
      const url=URL.createObjectURL(new Blob([text],{type:'image/svg+xml'}));
      visible={text,url,label:type,caption,count};picture.src=url;picture.hidden=false;
      $('drawing-summary').textContent=[scope,caption,preparedCells?.warning].filter(Boolean).join(' · ');
      $('canvas-sheet').classList.toggle('transparent-preview',Boolean(studio.state.styleOptions.transparent_background));
      $('representation-note').hidden=true;
      const bar=$('exports');bar.replaceChildren();
      for(const ext of ['png','svg',...['cif','json'].filter(ext=>model.downloads?.[ext])]){
        const b=document.createElement('button');b.type='button';b.className='button outline';b.textContent=ext.toUpperCase();
        b.onclick=()=>saveFigure(ext).catch(error=>studio.notice(error.message,true));bar.append(b);
      }
      studio.sync();window.dispatchEvent(new CustomEvent('cofdraw:crystalviewready',{detail:visible}));
    }catch(error){
      if(token!==job||report!==model)return;
      cache.delete(count);$('drawing-summary').textContent=error.message;studio.sync();
    }
  }
  async function saveFigure(ext){
    if(!visible)return;
    if(ext==='cif'||ext==='json')return studio.saveExport(report.downloads[ext],ext);
    const name=$('canvas-title').textContent;
    if(ext==='svg')return studio.saveText(name+'.svg',visible.text,'image/svg+xml');
    const img=new Image();img.src=visible.url;await img.decode();
    const box=new DOMParser().parseFromString(visible.text,'image/svg+xml').documentElement.viewBox.baseVal;
    const size=Number($('image-size').value),factor=size/Math.max(box.width,box.height);
    const canvas=document.createElement('canvas');canvas.width=Math.max(1,Math.round(box.width*factor));canvas.height=Math.max(1,Math.round(box.height*factor));
    canvas.getContext('2d').drawImage(img,0,0,canvas.width,canvas.height);
    if(window.pywebview?.api?.save_image)return window.pywebview.api.save_image(name+'.png',canvas.toDataURL());
    const file=await new Promise(resolve=>canvas.toBlob(resolve));const url=URL.createObjectURL(file);
    const a=document.createElement('a');a.href=url;a.download=name+'.png';a.click();setTimeout(()=>URL.revokeObjectURL(url),3000);
  }
  function chooseStyle(style){
    const state=studio.state;
    if(style==='random'){
      const family=id=>['reticular_depth','reticular_blocks'].includes(id)?'depth':id;
      const options=availableStyles().filter(s=>family(s.id)!==family(state.style));style=options[Math.floor(Math.random()*options.length)].id;
    }
    if(!availableStyles().some(s=>s.id===style))return;
    state.style=style;state.styleOptions={};state.resolvedStyleOptions={};studio.refreshStyles();drawView();studio.sync();
  }
  function enterDrawing(){
    if(!$('science-dialog').hidden)$('close-science').click();
    if(!$('space-workspace').hidden)window.dispatchEvent(new Event('cofdraw:leave-explorer'));
    if(!$('science-dialog').hidden||main.hidden)$('close-science').click();
    $('mode-manual').click();
    if(!patchStyle){
      const state=studio.state;
      patchStyle={style:state.style,styleOptions:state.styleOptions,resolvedStyleOptions:state.resolvedStyleOptions,styleCategory:state.styleCategory};
      state.style='reticular_depth';state.styleCategory='all';state.styleOptions={};state.resolvedStyleOptions={};count=1;
    }
    main.classList.add('crystal-design-active');$('crystal-design-controls').hidden=false;$('crystal-view-field').hidden=false;
    studio.refreshStyles();
  }
  function startCustom(){
    id=custom;edit=null;failure='';restoredReference=false;clearTimeout(timer);revision++;clearPicture();
    for(const control of ['dia-linkage','dia-net-count','reticular-net'])$(control).closest('label').hidden=false;
    $('cof-preset').value='crystal:'+custom;$('compare-crystal-reference').hidden=true;
    $('canvas-title').textContent=[precursor('a').label||'A',precursor('b').label||'B'].join(' + ');
    const nets=+$('dia-net-count').value;
    scope=$('reticular-net').selectedOptions[0].textContent+' · '+nets+(nets===1?' net':' interpenetrating nets');
    picker.replaceChildren(new Option('Current A + B · '+$('reticular-net').value,custom));
    $('crystal-linker-note').textContent='Changing either monomer, linkage or network count rebuilds both views.';
    $('drawing-summary').textContent='Preparing the selected precursor pair…';
    lab.invalidateCrystalDesign('The monomer pair changed.');studio.sync();
  }
  function inputsChanged(){
    if(!id)return;
    if(restoredReference&&entries.get(id)?.starter?.kind==='reference'){
      leave();studio.queueDrawing();
      studio.notice('The reference model has been closed. The new precursors are editable in the drawing workspace.');
      return;
    }
    startCustom();
    if(!precursor('a').smiles||!precursor('b').smiles){
      $('drawing-summary').textContent='Choose or draw both building blocks.';return;
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
    failure=error.message;$('drawing-summary').textContent=failure;
    if(!$('science-dialog').hidden)lab.status(failure);
    studio.sync();
  }
  async function select(key,modification=null,restoreReference=false){
    if(key===custom){id=custom;enterDrawing();startCustom();studio.setCrystalMonomers();await build();return;}
    const row=entries.get(key);if(!row)return;
    const preset=row.starter;
    restoreReference=restoreReference||preset?.kind==='reference';
    clearTimeout(timer);const token=++revision;
    id=key;edit=modification;failure='';restoredReference=restoreReference;clearPicture();enterDrawing();$('dia-linkage').value='imine';$('dia-net-count').value='1';$('reticular-net').value='dia';
    for(const control of ['dia-linkage','dia-net-count','reticular-net'])$(control).closest('label').hidden=restoreReference;
    if(preset?.net==='lon')$('reticular-net').value='lon';
    const ids=edit?[]:(preset?.monomer_ids||row.monomer_ids||matched[key]||[]);
    const items=ids.map(k=>studio.state.library?.monomer_presets?.find(m=>m.id===k));
    studio.setCrystalMonomers(items);
    $('cof-preset').value='crystal:'+id;linkerChoices(row);
    $('compare-crystal-reference').hidden=restoreReference||!row.path;
    $('canvas-title').textContent=restoreReference?row.model_name:[precursor('a').label,precursor('b').label].filter(Boolean).join(' + ');
    scope=restoreReference?'Published '+(preset?.net||'crystal')+' model · '+row.model_name:
      'Generated '+(preset?.net||'dia')+' design · single net'+(preset?.common_name?' · '+preset.common_name+' precursors':'');
    $('drawing-summary').textContent='Preparing the drawing and unit cell…';
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
      let cells={viewer:result.viewer,repeats:result.repeats||[1,1,1]};
      if(!restoredReference){
        const n=result.structure.atoms.length,counts=[2,2,2];
        const lengths=result.structure.cell.map(axis=>Math.hypot(...axis));
        // Prefer repeats along shorter lattice vectors when the full set would exceed the viewer limit.
        for(const axis of [0,1,2].sort((a,b)=>lengths[b]-lengths[a])){
          if(n*counts.reduce((p,k)=>p*k,1)<=12000)break;
          counts[axis]=1;
        }
        cells={viewer:result.viewer,repeats:[1,1,1]};
        if(counts.some(k=>k>1)){
          try{cells=await studio.api('science/reports/'+result.report_id+'/view',{repeats:counts});}
          catch(error){cells.warning='Showing the source cell. The repeated view could not be prepared. '+error.message;}
        }
        if(token!==revision)return;
      }
      if(result.chemical_framework_svg)choices.push('Framework drawing');
      if(result.chemical_circuit_svg)choices.push('Projected circuit');
      if(cells&&result.chemical_framework_svg)choices.push('Chemical drawing + periodic network');
      if(cells?.viewer?.bonds?.length)choices.push('Periodic packing');
      if(result.chemical_drawing_svg)choices.push('Linker detail');
      report=result;preparedCells=cells;cache.set(1,{scheme_svg:result.chemical_framework_svg,projection_svg:result.chemical_circuit_svg});
      if(!availableStyles().some(s=>s.id===studio.state.style))studio.state.style='reticular_depth';
      studio.refreshStyles();
      const selector=$('crystal-drawing-view');selector.replaceChildren(...choices.map(label=>new Option(label)));
      if(choices.includes(preferredView))selector.value=preferredView;
      studio.sync();
      if(choices.length)await drawView();
      else $('drawing-summary').textContent='Periodic cell ready. A chemical circuit could not be resolved.';
      if(token!==revision)return;
      if(!$('science-dialog').hidden&&!comparing)showPeriodic();
  }
  window.addEventListener('cofdraw:crystal-catalog',e=>{for(const r of e.detail)entries.set(r.id,r);});
  window.addEventListener('cofdraw:open-crystal',e=>select(e.detail));
  picker.addEventListener('change',()=>{if(picker.value==='aryl_f'){picker.value=id;selectInViewer(id,'aryl_f');}else select(picker.value);});
  $('compare-crystal-reference').onclick=()=>selectInViewer(id);
  $('crystal-drawing-view').addEventListener('change',()=>drawView());
  $('measure-pore').addEventListener('click',e=>{
    if(!id)return;
    e.preventDefault();e.stopImmediatePropagation();
    if(report){$('open-science').click();lab.tab('porosity');}
  },true);
  function returnToPatch(){leave();studio.queueDrawing();}
  $('leave-crystal-design').addEventListener('click',returnToPatch);
  $('build-periodic-from-pair').onclick=()=>{
    const button=$('build-periodic-from-pair');
    if(!button.hidden&&button.dataset.reference){select(button.dataset.reference);return;}
    if(button.hidden||!button.dataset.linkage)return;
    $('dia-linkage').value=button.dataset.linkage;
    $('reticular-net').value=button.dataset.net||'dia';$('dia-net-count').value='1';
    select(custom);
  };
  $('dia-linkage').onchange=()=>{inputsChanged();studio.refreshChemistry();};
  $('dia-net-count').onchange=inputsChanged;
  $('reticular-net').onchange=()=>{inputsChanged();studio.refreshChemistry();};
  for(const element of ['mode-describe','new-project'])$(element)?.addEventListener('click',leave,true);
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
      id=custom;enterDrawing();$('reticular-net').value=topology;$('dia-linkage').value=linkage;$('dia-net-count').value=String(count);startCustom();studio.setCrystalMonomers(project.monomers);await build();
    }else if(!entries.has(project.template_id)||![null,'aryl_f'].includes(project.linker_edit??null))
      throw new Error('Unknown crystal model or linker modification in this project.');
    else await select(project.template_id,project.linker_edit??null,project.version===1||project.reference_only===true);
    if(report){
      const legacy={'Depth fade · ink and ring wash':'reticular_depth','Depth fade · two building blocks':'reticular_blocks'};
      studio.state.style=availableStyles().some(s=>s.id===project.style)?project.style:legacy[project.view]||'reticular_depth';
      studio.state.styleOptions=project.style_options&&typeof project.style_options==='object'?project.style_options:{};
      const view=legacy[project.view]||project.view==='Three connected circuits'?'Framework drawing':project.view;
      const selector=$('crystal-drawing-view');
      if([...selector.options].some(o=>o.value===view))selector.value=view;
      studio.refreshStyles();await drawView(project.circuits||(project.view==='Three connected circuits'?3:1));
    }
  }
  function save(){
    if(report)studio.saveText('COFdraw-crystal-project.json',JSON.stringify({format:'cofdraw-crystal-project',reaction_scheme:window.COFdrawScheme?.snapshot(),
      version:2,circuits:count,style:studio.state.style,style_options:studio.state.styleOptions,template_id:id,linker_edit:edit,reference_only:restoredReference,...(id===custom?{monomers:['a','b'].map(precursor),linkage:$('dia-linkage').value,network_count:+$('dia-net-count').value,topology:$('reticular-net').value}:{}),view:$('crystal-drawing-view').value},null,2)).catch(e=>studio.notice(e.message,true));
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
  window.COFdrawCrystalDesign={restore,save,get styles(){return availableStyles();},chooseStyle,redraw:()=>drawView(),get visible(){return visible;},get scope(){return scope;},get isReference(){return restoredReference;},precursor,choice,build,inputsChanged,selectInViewer,has:key=>entries.has(key),get active(){return id;},get report(){return report;},leave,showPeriodic};
})();
