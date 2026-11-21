/* Coordinate-based publication figures. All calculations retain the parent CIF. */
(()=>{
  'use strict';
  const $=id=>document.getElementById(id), lab=window.COFdrawScience, viewer=lab.viewer;
  const scene=$('crystal-scene'), tools=$('crystal-figure-tools'), sheet=$('crystal-sheet');
  const toggle=$('crystal-figure-toggle'), save=$('crystal-figure-save'), note=$('crystal-scene-note');
  for(const slot of ['a','b']){
    const card=$('card-'+slot),details=document.createElement('button');
    details.type='button';details.className='text-button crystal-card-details';
    details.textContent='Reactivity and properties';details.setAttribute('aria-expanded','false');
    details.setAttribute('aria-controls','inspector-'+slot);
    details.onclick=()=>details.setAttribute('aria-expanded',String(card.classList.toggle('properties-open')));
    card.insertBefore(details,$('inspector-'+slot));
  }
  const colors={C:'#465568',N:'#3576b6',O:'#c95550',S:'#b78d37',B:'#a67e64',Si:'#808597',F:'#398f77'};
  const palette=['#166d83','#bd654a','#8071ac','#a4882e','#39816b','#b45d8a','#5077ac','#998766','#727b4e'];
  const labels={ring:'Molecular circuit',junction:'Tetrahedral junction',net:'Junction net',overlay:'Structure + net',single:'One network',cell:'Crystal packing'};
  labels.spotlight='Highlight a network';
  Object.assign(labels,{blockA:'Building block A',blockB:'Building block B',cage_overlay:'Circuit + topology',layer:'Single layer',stack:'Layer stack',pores:'Pore surface'});
  labels.blocknet='Building-block net';labels.block_overlay='Building blocks + net';
  const descriptions={ring:'Closed circuit traced through a single covalent network. All heavy atoms along its linkers are shown.',
    junction:'Four bonded shells around one four-connected center, unwrapped across cell faces.',
    net:'Each segment follows a traced linker between two four-connected centers. Two junction shells are shown.',
    overlay:'Dashed guides connect junction centers along four traced linker paths.',
    single:'One distance-connected periodic network within the displayed cell boundaries.',
    cell:'Networks are colored by connectivity in the selected coordinate model.'};
  descriptions.spotlight='Selected network in color with the surrounding crystal faded.';
  Object.assign(descriptions,{blockA:'One reacted building block A, unwrapped across periodic boundaries.',
    blockB:'One reacted building block B, unwrapped across periodic boundaries.',
    cage_overlay:'Covalent circuit with its junction-to-junction paths superimposed.',
    layer:'One explicitly assigned layer of the displayed structure.',stack:'All displayed layers colored by layer index.',
    pores:'Probe-center surface calculated on the current periodic cell.'});
  descriptions.blocknet='Centers of reacted building blocks, joined through the actual interblock bonds.';
  descriptions.block_overlay='Two coordination shells of the framework with its building-block graph overlaid.';
  let source=null,full=null,isFigure=false;
  const sceneChoices=new WeakMap();

  function part(model,accept){
    const ids=new Map(), atoms=[], bonds=[];
    model.atoms.forEach((a,i)=>{if(accept(a)){ids.set(i,atoms.length);atoms.push(a);}});
    for(const [a,b] of model.bonds)if(ids.has(a)&&ids.has(b))bonds.push([ids.get(a),ids.get(b)]);
    return {atoms,bonds};
  }
  function choose(){
    const key=scene.value;
    let model=full;
    const selected=Number($('crystal-net-index').value);
    $('crystal-net-field').hidden=!['single','spotlight'].includes(key);
    if(key==='layer'){
      const index=full.atoms.find(a=>a.layer!==undefined)?.layer;
      model=part(full,a=>a.layer===index);
    }else if(key==='single'){
      model={...part(full,a=>a.network===selected),cell:full.cell};
    }else if(key==='spotlight'){
      model={...full,highlightNetwork:selected};
    }else if(!['cell','stack','pores'].includes(key))model=source.crystal_details[key];
    viewer.set(['cell','spotlight','stack','pores'].includes(key)?model:{...model,illustration:true});
    const direction=model.direction;
    if(direction){
      viewer.theta=Math.atan2(-direction[0],direction[2]);
      viewer.phi=Math.asin(Math.min(1,Math.max(-1,direction[1])));
    }
    $('viewer-axis').value='free';
    $('viewer-color').value=['layer','stack'].includes(key)?'layer':['single','spotlight'].includes(key)?'network':key==='blocknet'?'component':'element';
    viewer.hideHydrogen=true;$('viewer-hide-h').checked=true;
    note.textContent=descriptions[key]+' Calculations use the full crystal.';
    window.dispatchEvent(new Event('cofdraw:crystal-scene'));
    if(key==='pores'){
      if(isFigure)mode(false);
      $('viewer-mode').value='wireframe';
      $('pore-visible').checked=true;
      window.COFdrawPores?.calculate().catch(error=>lab.status(error.message,true));
    }
    toggle.disabled=key==='pores';
    viewer.draw();paint();
  }

  function node(tag,attrs,text){
    const e=document.createElementNS('http://www.w3.org/2000/svg',tag);
    for(const [key,val] of Object.entries(attrs))e.setAttribute(key,val);
    if(text!==undefined)e.textContent=text;
    return e;
  }
  function chemicalNetworkFigure(result,displayCells){
    const model=displayCells.viewer,chem=new DOMParser().parseFromString(result.chemical_framework_svg,'image/svg+xml').documentElement;
    if(chem.localName!=='svg'||!model?.bonds?.length)return null;
    const net=result.structure.provenance.net||'';
    const svg=node('svg',{xmlns:'http://www.w3.org/2000/svg',viewBox:'0 0 1160 580',width:1160,height:580,role:'img'});
    svg.append(node('title',{},result.structure.provenance.model_name+' · chemical drawing and periodic network'),
      node('desc',{},'Left: bond-line diagram. Right: oblique projection of the repeated cell using its periodic bond translations. Depth is encoded by line weight and opacity.'),
      node('rect',{width:1160,height:580,fill:'white'}));
    const drawing=document.importNode(chem,true);
    for(const [key,value] of Object.entries({x:20,y:48,width:535,height:500,preserveAspectRatio:'xMidYMid meet'}))drawing.setAttribute(key,value);
    const title=(x,text)=>node('text',{x,y:31,fill:'#203e43','font-family':'sans-serif','font-size':17},text);
    svg.append(title(28,'Chemical drawing'),drawing,title(594,'Periodic network · '+net+' · '+displayCells.repeats.join(' × ')));
    const modelGroup=node('g',{'data-view':'periodic-network','data-displayed-sites':model.atoms.length,
      'data-source-sites':result.structure.atoms.length,'stroke-linecap':'round','stroke-linejoin':'round',fill:'none'});
    svg.append(modelGroup);
    const ct=Math.cos(-.64),st=Math.sin(-.64),cp=Math.cos(.39),sp=Math.sin(.39);
    const projected=model.atoms.map(a=>{
      const depth=-a.x*st+a.z*ct;
      return [a.x*ct+a.z*st,a.y*cp-depth*sp,a.y*sp+depth*cp];
    });
    const visible=projected.filter((pt,i)=>model.atoms[i].element!=='H');
    if(!visible.length)return null;
    const bounds=[0,1,2].map(k=>[Math.min(...visible.map(p=>p[k])),Math.max(...visible.map(p=>p[k]))]);
    const range=bounds.map(([lo,hi])=>Math.max(.01,hi-lo)),mid=bounds.map(([lo,hi])=>(lo+hi)/2);
    const scale=Math.min(525/range[0],475/range[1]);
    const xy=projected.map(([x,y])=>[865+(x-mid[0])*scale,295-(y-mid[1])*scale]);
    const steps=16,paths=Array.from({length:steps},()=>new Map());
    const color=a=>a.element==='C'?(a.component==='A'?'#397280':a.component==='B'?'#88633e':'#526268'):colors[a.element]||'#765888';
    for(const [a,b] of model.bonds){
      if(model.atoms[a].element==='H'||model.atoms[b].element==='H')continue;
      const half=[(xy[a][0]+xy[b][0])/2,(xy[a][1]+xy[b][1])/2];
      for(const [i,j] of [[a,b],[b,a]]){
        const z=(3*projected[i][2]+projected[j][2])/4;
        const k=Math.min(steps-1,Math.max(0,Math.floor(steps*(z-bounds[2][0])/range[2])));
        const shade=color(model.atoms[i]),bucket=paths[k];
        if(!bucket.has(shade))bucket.set(shade,[]);
        bucket.get(shade).push('M'+xy[i][0].toFixed(2)+' '+xy[i][1].toFixed(2)+'L'+half[0].toFixed(2)+' '+half[1].toFixed(2));
      }
    }
    paths.forEach((bucket,i)=>{
      const depth=i/(steps-1);
      for(const [shade,lines] of bucket)modelGroup.append(node('path',{d:lines.join(''),stroke:shade,
        'stroke-opacity':(.26+.74*depth).toFixed(3),'stroke-width':(.70+1.15*depth).toFixed(2)}));
    });
    return new XMLSerializer().serializeToString(svg);
  }
  window.COFdrawChemicalNetworkFigure=chemicalNetworkFigure;
  function paint(){
    if(!isFigure||!source||!viewer.atoms.length)return;
    const atoms=viewer.atoms, key=scene.value, net=['net','blocknet'].includes(key), packing=key==='cell';
    const name=source.structure.provenance?.model_name||'Current structure';
    const svg=node('svg',{xmlns:'http://www.w3.org/2000/svg',viewBox:'0 0 760 570',role:'img','aria-label':name+' · '+labels[key]});
    svg.append(node('title',{},name+' · '+labels[key]),node('desc',{},note.textContent+' '+(source.structure.provenance?.doi||'')),
      node('rect',{width:760,height:570,fill:'white'}));
    const {theta,phi}=viewer, ct=Math.cos(theta),st=Math.sin(theta),cp=Math.cos(phi),sp=Math.sin(phi);
    const pts=atoms.map(({x,y,z})=>{const depth=-x*st+z*ct;return [x*ct+z*st, y*cp-depth*sp,y*sp+depth*cp];});
    const visible=pts.filter((p,i)=>atoms[i].element!=='H');
    if(!visible.length)return;
    const bounds=[0,1,2].map(k=>[Math.min(...visible.map(p=>p[k])),Math.max(...visible.map(p=>p[k]))]);
    const [width,height,depth]=bounds.map(([lo,hi])=>Math.max(.01,hi-lo));
    const scale=Math.min(638/width,412/height);
    const xy=pts.map(([x,y])=>[380+(x-bounds[0][0]-width/2)*scale,272-(y-bounds[1][0]-height/2)*scale]);
    const color=i=>packing&&Number.isInteger(atoms[i].network)?palette[atoms[i].network%palette.length]:colors[atoms[i].element]||'#795f90';
    for(const [a,b] of viewer.model.guides||[])svg.append(node('line',{x1:xy[a][0],y1:xy[a][1],x2:xy[b][0],y2:xy[b][1],
      stroke:'#bc6a36','stroke-width':4,'stroke-dasharray':'7 5',opacity:.65}));
    const guideXY=([x,y,z])=>[380+(x*ct+z*st-bounds[0][0]-width/2)*scale,
      272-(y*cp-(-x*st+z*ct)*sp-bounds[1][0]-height/2)*scale];
    for(const guide of viewer.model.spatial_guides||[]){
      const [a,b]=guide.map(guideXY);
      svg.append(node('line',{x1:a[0],y1:a[1],x2:b[0],y2:b[1],stroke:'#bc6a36','stroke-width':3,'stroke-dasharray':'7 5',opacity:.75}));
    }
    const commands=[];
    for(const [a,b] of viewer.bonds){
      if(atoms[a].element==='H'||atoms[b].element==='H')continue;
      const [x1,y1]=xy[a],[x2,y2]=xy[b],mx=(x1+x2)/2,my=(y1+y2)/2;
      for(const [i,x,y] of [[a,x1,y1],[b,x2,y2]]){
        const z=(2*pts[i][2]+pts[a][2]+pts[b][2])/4;
        commands.push([z,node('line',{x1:x,y1:y,x2:mx,y2:my,stroke:net?'#388293':color(i),
          'stroke-width':net?3.2:Math.max(1,Math.min(2.6,scale*.12)),'stroke-linecap':'round'})]);
      }
    }
    atoms.forEach((a,i)=>{
      if(a.element==='H')return;
      const radius=net?5.2:Math.min(4.8,Math.max(.9,scale*(a.element==='C'?.08:.19)));
      commands.push([pts[i][2],node('circle',{cx:xy[i][0],cy:xy[i][1],r:radius,fill:net?'#266575':color(i)})]);
    });
    for(const [z,el] of commands.sort((a,b)=>a[0]-b[0])){
      el.setAttribute('opacity',(.43+.57*(z-bounds[2][0])/depth).toFixed(2));svg.append(el);
    }
    svg.append(node('text',{x:34,y:518,fill:'#193b46','font-family':'Georgia, serif','font-size':24},name+' · '+labels[key]),
      node('text',{x:35,y:543,fill:'#556573','font-family':'sans-serif','font-size':13},'Orthographic projection · '+(source.structure.provenance?.doi||'Selected coordinate model')));
    sheet.replaceChildren(svg);
  }
  function mode(active){
    isFigure=active;sheet.hidden=!active;save.hidden=!active;
    $('model-view').classList.toggle('figure-active',active);
    toggle.textContent=active?'Rotate in 3D':'2D figure';toggle.setAttribute('aria-pressed',String(active));
    if(!$('science-dialog').hidden)$('switch-space').textContent='2D drawing';
    if(active)paint();else viewer.draw();
  }
  function selectScene(){
    if(source?.structure)sceneChoices.set(source.structure,{key:scene.value,network:$('crystal-net-index').value});
    choose();
  }
  scene.addEventListener('change',selectScene);
  $('crystal-net-index').addEventListener('change',selectScene);
  $('crystal-annotation').addEventListener('click',()=>{
    window.dispatchEvent(new CustomEvent('cofdraw:open-crystal',{detail:'cof-300-c5'}));
  });
  toggle.addEventListener('click',()=>mode(!isFigure));
  save.addEventListener('click',()=>window.COFdrawStudio.saveText('COFdraw-'+scene.value+'.svg',
    new XMLSerializer().serializeToString(sheet.firstElementChild),'image/svg+xml').catch(e=>lab.status(e.message,true)));
  function networks(){
    const ids=[...new Set(full.atoms.map(a=>a.network).filter(Number.isInteger))].sort((a,b)=>a-b);
    const select=$('crystal-net-index');select.replaceChildren();
    for(const id of ids)select.add(new Option(String(id+1),id));
    $('crystal-net-field').hidden=true;
    return ids;
  }
  window.addEventListener('cofdraw:structure',({detail:r})=>{
    if(r.invalidated&&lab.structure===source?.structure)return;
    source=null;full=null;tools.hidden=true;note.hidden=true;sheet.replaceChildren();mode(false);
    if(!r.structure||!r.viewer?.atoms?.length)return;
    source=r;full=r.viewer;
    scene.replaceChildren();
    const ids=networks();
    toggle.disabled=false;
    const layers=new Set(full.atoms.map(a=>a.layer).filter(x=>x!==undefined));
    for(const key of ['cell','layer','stack','blockA','blockB','blocknet','block_overlay','ring','junction','net','overlay','cage_overlay','single','spotlight','pores']){
      const supported=key==='cell'||(key==='single'&&ids.length)||(key==='spotlight'&&ids.length>1)||
        (key==='layer'&&layers.size)||(key==='stack'&&layers.size>1)||(key==='pores'&&r.structure.periodic)||r.crystal_details?.[key];
      if(supported)scene.add(new Option(key==='cell'&&!r.structure.periodic?'Full structure':labels[key],key));
    }
    tools.hidden=false;note.hidden=false;
    $('crystal-annotation').hidden=r.structure.provenance?.template!=='cof-300-c5';
    note.textContent=descriptions.cell;
    if(lab.origin==='catalog-design'){
      viewer.hideHydrogen=true;$('viewer-hide-h').checked=true;
      $('viewer-color').value=ids.length>1?'network':'element';viewer.draw();
      const previous=sceneChoices.get(source.structure);
      if(previous&&[...scene.options].some(o=>o.value===previous.key)){
        scene.value=previous.key;
        if([...$('crystal-net-index').options].some(o=>o.value===previous.network))$('crystal-net-index').value=previous.network;
        choose();
      }
    }
  });
  window.addEventListener('cofdraw:display-cells',()=>{
    if(!source)return;
    const key=scene.value, previous=$('crystal-net-index').value;
    full=viewer.model;networks();
    if([...$('crystal-net-index').options].some(o=>o.value===previous))$('crystal-net-index').value=previous;
    if(key!=='cell')choose();else{note.textContent=descriptions.cell;paint();}
  });
  for(const id of ['viewer-axis','viewer-reset','viewer-preset'])$(id).addEventListener(id==='viewer-reset'?'click':'change',()=>requestAnimationFrame(paint));
})();
