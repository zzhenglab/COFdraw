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
  const labels={ring:'Molecular circuit',junction:'Tetrahedral junction',net:'Junction net',overlay:'Structure + net',single:'One network',cell:'Crystal packing'};
  labels.spotlight='Highlight a network';
  Object.assign(labels,{blockA:'Building block A',blockB:'Building block B',cage_overlay:'Circuit + topology',layer:'Single layer',stack:'Layer stack',pores:'Pore surface'});
  labels.blocknet='Building-block net';labels.block_overlay='Building blocks + net';
  const descriptions={ring:'Closed circuit traced through a single covalent network. All heavy atoms along its linkers are shown.',
    junction:'Four bonded shells around one four-connected center, unwrapped across cell faces.',
    net:'Each segment follows a traced linker between two four-connected centers. Two junction shells are shown.',
    overlay:'Dashed guides connect junction centers along four traced linker paths.',
    single:'One distance-connected periodic network within the displayed cell boundaries.',
    cell:'Atomic packing in the displayed cells.'};
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
  const topology=$('crystal-topology'), interpenetration=$('crystal-interpenetration'), buildNets=$('crystal-build-nets');
  const packingScenes=['cell','single','spotlight','stack','pores'];
  let repeats=[1,1,1], networkIds=[], preferredOverlay=false, preferSpotlight=false;
  for(let count=1;count<=7;count++)buildNets.add(new Option(String(count),String(count)));

  // Clip each translated crystallographic edge at the planes of the displayed
  // supercell. All coordinates here are fractional until converted for display.
  function periodicGuides(net,model,counts,network){
    const cell=model.cell.map((row,k)=>row.map(n=>n/counts[k]));
    const xyz=f=>[0,1,2].map(k=>f.reduce((v,n,i)=>v+n*cell[i][k],0));
    const n=source.structure.atoms.length, guides=[], networks=[], nodes=new Map();
    function networkAt(vertex,shift){
      if(!net.site_ids)return undefined;
      const anchor=shift.map((x,k)=>x+(net.site_shifts?.[vertex][k]||0));
      if(anchor.some((x,k)=>x<0||x>=counts[k]))return undefined;
      const image=(anchor[0]*counts[1]+anchor[1])*counts[2]+anchor[2];
      return model.atoms[image*n+net.site_ids[vertex]]?.network;
    }
    for(const [i,j,...image] of net.edges){
      const p=net.points[i], q=net.points[j].map((x,k)=>x+image[k]), d=q.map((x,k)=>x-p[k]);
      const lo=p.map((x,k)=>Math.ceil(-Math.max(x,q[k]))),hi=p.map((x,k)=>Math.floor(counts[k]-Math.min(x,q[k])));
      for(let a=lo[0];a<=hi[0];a++)for(let b=lo[1];b<=hi[1];b++)for(let c=lo[2];c<=hi[2];c++){
        const shift=[a,b,c], start=p.map((x,k)=>x+shift[k]);
        const endShift=shift.map((x,k)=>x+image[k]), id=networkAt(i,shift)??networkAt(j,endShift);
        if(network!==null&&id!==network)continue;
        let enter=0,exit=1;
        for(let k=0;k<3;k++){
          if(Math.abs(d[k])<1e-10){
            if(start[k]<-1e-8||start[k]>counts[k]+1e-8)exit=-1;
          }else{
            const first=-start[k]/d[k],last=(counts[k]-start[k])/d[k];
            enter=Math.max(enter,Math.min(first,last));exit=Math.min(exit,Math.max(first,last));
          }
        }
        if(exit-enter<=1e-7)continue;
        const [from,to]=[enter,exit].map(t=>xyz(start.map((x,k)=>x+t*d[k])));
        guides.push([from,to]);networks.push(id);
        if(enter<1e-8)nodes.set([i,...shift].join(','),{position:from,network:id});
        if(exit>1-1e-8)nodes.set([j,...endShift].join(','),{position:to,network:id});
      }
    }
    return {spatial_guides:guides,guide_networks:networks,topology_nodes:[...nodes.values()]};
  }
  function syncPackingControls(key){
    const net=source.crystal_details?.periodic_topology;
    $('crystal-topology-field').hidden=!net||!packingScenes.includes(key)||key==='pores';
    const canBuild=source.structure===window.COFdrawCrystalDesign?.report?.structure&&
      source.structure.provenance?.kind==='monomer_dia_design';
    interpenetration.hidden=!(networkIds.length>1||canBuild);
    interpenetration.textContent='Interpenetration · '+networkIds.length+(networkIds.length===1?' net':' nets');
    $('crystal-build-nets-field').hidden=!canBuild||!preferSpotlight;
    if(canBuild)buildNets.value=String(source.structure.provenance.network_count||networkIds.length);
  }

  function part(model,accept){
    const ids=new Map(), atoms=[], bonds=[];
    model.atoms.forEach((a,i)=>{if(accept(a)){ids.set(i,atoms.length);atoms.push(a);}});
    for(const [a,b] of model.bonds)if(ids.has(a)&&ids.has(b))bonds.push([ids.get(a),ids.get(b)]);
    return {atoms,bonds};
  }
  function choose(keepAppearance=false){
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
    syncPackingControls(key);
    const drawOverlay=topology.checked&&!$('crystal-topology-field').hidden;
    if(drawOverlay)model={...model,topologyOverlay:true,...periodicGuides(source.crystal_details.periodic_topology,full,repeats,
      ['single','spotlight'].includes(key)?selected:null)};
    viewer.set(['cell','spotlight','stack','pores'].includes(key)?model:{...model,illustration:true},keepAppearance);
    const direction=model.direction;
    if(direction&&!keepAppearance){
      viewer.theta=Math.atan2(-direction[0],direction[2]);
      viewer.phi=Math.asin(Math.min(1,Math.max(-1,direction[1])));
    }
    if(!keepAppearance){
      $('viewer-axis').value='free';
      $('viewer-color').value=['layer','stack'].includes(key)?'layer':['single','spotlight'].includes(key)||(key==='cell'&&networkIds.length>1)?'network':key==='blocknet'?'component':'element';
    }
    note.textContent=descriptions[key]+(drawOverlay?' Overlay: '+source.crystal_details.periodic_topology.kind.toLowerCase()+'.':'');
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
  function chemicalNetworkFigure(result,displayCells,packingOnly=false){
    const model=displayCells.viewer;
    if(!model?.bonds?.length)return null;
    const net=result.structure.provenance.net||'';
    const svg=node('svg',{xmlns:'http://www.w3.org/2000/svg',viewBox:'0 0 1160 580',width:1160,height:580,role:'img'});
    svg.append(node('title',{},result.structure.provenance.model_name+' · chemical drawing and periodic network'),
      node('desc',{},'Left: bond-line diagram. Right: oblique projection of the repeated cell using its periodic bond translations. Depth is encoded by line weight and opacity.'),
      node('rect',{width:1160,height:580,fill:'white','data-paper':''}));
    const title=(x,text)=>node('text',{x,y:31,fill:'#203e43','font-family':'sans-serif','font-size':17},text);
    if(packingOnly){
      svg.setAttribute('viewBox','580 0 580 580');svg.setAttribute('width','580');
      svg.querySelector('title').textContent=result.structure.provenance.model_name+' · periodic packing';
      svg.querySelector('desc').textContent='Heavy atoms and periodic bonds projected from the selected crystal with depth shading.';
    }else{
      const chem=new DOMParser().parseFromString(result.chemical_framework_svg,'image/svg+xml').documentElement;
      if(chem.localName!=='svg')return null;
      const drawing=document.importNode(chem,true);
      for(const [key,value] of Object.entries({x:20,y:48,width:535,height:500,preserveAspectRatio:'xMidYMid meet'}))drawing.setAttribute(key,value);
      svg.append(title(28,'Chemical drawing'),drawing);
    }
    svg.append(title(594,'Periodic network · '+net+' · '+displayCells.repeats.join(' × ')));
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
    const ink={N:'#24679f',O:'#c14848',F:'#2a8c70',S:'#b68c26',Si:'#808597'};
    const color=a=>a.element==='C'?(a.component==='A'?'#397280':a.component==='B'?'#88633e':'#526268'):ink[a.element]||'#765888';
    for(const [a,b] of model.bonds){
      if(model.atoms[a].element==='H'||model.atoms[b].element==='H')continue;
      const half=[(xy[a][0]+xy[b][0])/2,(xy[a][1]+xy[b][1])/2];
      for(const [i,j] of [[a,b],[b,a]]){
        const z=(3*projected[i][2]+projected[j][2])/4;
        const k=Math.min(steps-1,Math.max(0,Math.floor(steps*(z-bounds[2][0])/range[2])));
        const atom=model.atoms[i],key=atom.element+':'+(atom.component||''),bucket=paths[k];
        if(!bucket.has(key))bucket.set(key,{atom,lines:[]});
        bucket.get(key).lines.push('M'+xy[i][0].toFixed(2)+' '+xy[i][1].toFixed(2)+'L'+half[0].toFixed(2)+' '+half[1].toFixed(2));
      }
    }
    paths.forEach((bucket,i)=>{
      const depth=i/(steps-1);
      for(const {atom,lines} of bucket.values())modelGroup.append(node('path',{d:lines.join(''),stroke:color(atom),
        'data-depth':depth.toFixed(3),'data-block':atom.component||'','data-element':atom.element,
        'stroke-opacity':(.26+.74*depth).toFixed(3),'stroke-width':(.70+1.15*depth).toFixed(2)}));
    });
    return new XMLSerializer().serializeToString(svg);
  }
  window.COFdrawChemicalNetworkFigure=chemicalNetworkFigure;
  function paint(){
    if(!isFigure||!source||!viewer.atoms.length)return;
    const atoms=viewer.atoms, key=scene.value, net=['net','blocknet'].includes(key);
    const name=source.structure.provenance?.model_name||'Current structure';
    const svg=node('svg',{xmlns:'http://www.w3.org/2000/svg',viewBox:'0 0 760 570',role:'img','aria-label':name+' · '+labels[key]});
    svg.append(node('title',{},name+' · '+labels[key]),node('desc',{},note.textContent+' '+(source.structure.provenance?.doi||'')),
      node('rect',{width:760,height:570,fill:'white'}));
    const {theta,phi}=viewer, ct=Math.cos(theta),st=Math.sin(theta),cp=Math.cos(phi),sp=Math.sin(phi);
    const pts=atoms.map(({x,y,z})=>{const depth=-x*st+z*ct;return [x*ct+z*st, y*cp-depth*sp,y*sp+depth*cp];});
    const hidden=i=>viewer.hideHydrogen&&atoms[i].element==='H';
    const visible=pts.filter((p,i)=>!hidden(i));
    if(!visible.length)return;
    const bounds=[0,1,2].map(k=>[Math.min(...visible.map(p=>p[k])),Math.max(...visible.map(p=>p[k]))]);
    const [width,height,depth]=bounds.map(([lo,hi])=>Math.max(.01,hi-lo));
    const scale=Math.min(638/width,412/height);
    const xy=pts.map(([x,y])=>[380+(x-bounds[0][0]-width/2)*scale,272-(y-bounds[1][0]-height/2)*scale]);
    const ink=viewer.colorsFor(atoms), color=i=>ink.get(i)||'#795f90', guides=[];
    for(const [a,b] of viewer.model.guides||[])guides.push(node('line',{x1:xy[a][0],y1:xy[a][1],x2:xy[b][0],y2:xy[b][1],
      stroke:'#bc6a36','stroke-width':3,'stroke-dasharray':'7 5',opacity:.8}));
    const guideXY=([x,y,z])=>[380+(x*ct+z*st-bounds[0][0]-width/2)*scale,
      272-(y*cp-(-x*st+z*ct)*sp-bounds[1][0]-height/2)*scale];
    for(const [i,guide] of (viewer.model.spatial_guides||[]).entries()){
      const [a,b]=guide.map(guideXY);
      guides.push(node('line',{x1:a[0],y1:a[1],x2:b[0],y2:b[1],stroke:viewer.topologyInk(viewer.model.guide_networks?.[i]),'stroke-width':2.6,
        'stroke-dasharray':viewer.model.topologyOverlay?'none':'7 5',opacity:.8}));
    }
    for(const vertex of viewer.model.topology_nodes||[]){
      const [x,y]=guideXY(vertex.position);
      guides.push(node('circle',{cx:x,cy:y,r:2.7,fill:viewer.topologyInk(vertex.network)}));
    }
    const commands=[];
    for(const [a,b] of viewer.bonds){
      if(hidden(a)||hidden(b))continue;
      const [x1,y1]=xy[a],[x2,y2]=xy[b],mx=(x1+x2)/2,my=(y1+y2)/2;
      for(const [i,x,y] of [[a,x1,y1],[b,x2,y2]]){
        const z=(2*pts[i][2]+pts[a][2]+pts[b][2])/4;
        commands.push([z,i,node('line',{x1:x,y1:y,x2:mx,y2:my,stroke:color(i),
          'stroke-width':net?3.2:Math.max(1,Math.min(2.6,scale*.12)),'stroke-linecap':'round'})]);
      }
    }
    atoms.forEach((a,i)=>{
      if(hidden(i))return;
      const radius=net?5.2:Math.min(4.8,Math.max(.9,scale*(a.element==='C'?.08:.19)));
      commands.push([pts[i][2],i,node('circle',{cx:xy[i][0],cy:xy[i][1],r:radius,fill:color(i)})]);
    });
    for(const [z,i,el] of commands.sort((a,b)=>a[0]-b[0])){
      const selected=viewer.model.highlightNetwork;
      const ghost=Number.isInteger(selected)&&atoms[i].network!==selected;
      const fade=ghost ? .1 : viewer.model.topologyOverlay ? .5 : 1;
      el.setAttribute('opacity',(fade*(.43+.57*(z-bounds[2][0])/depth)).toFixed(2));svg.append(el);
    }
    svg.append(...guides);
    svg.append(node('text',{x:34,y:518,fill:'#193b46','font-family':'Georgia, serif','font-size':24},name+' · '+labels[key]),
      node('text',{x:35,y:543,fill:'#556573','font-family':'sans-serif','font-size':13},'Orthographic projection · '+(source.structure.provenance?.doi||'Selected coordinate model')));
    sheet.replaceChildren(svg);
  }
  function mode(active){
    isFigure=active;sheet.hidden=!active;save.hidden=!active;
    $('model-view').classList.toggle('figure-active',active);
    toggle.textContent=active?'Rotate in 3D':'2D figure';toggle.setAttribute('aria-pressed',String(active));
    if(active)paint();else viewer.draw();
  }
  function selectScene(){
    if(source?.structure)sceneChoices.set(source.structure,{key:scene.value,network:$('crystal-net-index').value});
    choose();
  }
  topology.onchange=()=>{preferredOverlay=topology.checked;choose(true);};
  interpenetration.onclick=()=>{
    preferSpotlight=true;
    if(networkIds.length>1)scene.value='spotlight';
    selectScene();
    if(!$('crystal-build-nets-field').hidden&&networkIds.length===1)buildNets.focus();
  };
  buildNets.onchange=()=>{
    const recipe=$('dia-net-count');
    recipe.value=buildNets.value;
    recipe.dispatchEvent(new Event('change'));
  };
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
    networkIds=ids;return ids;
  }
  window.addEventListener('cofdraw:structure',({detail:r})=>{
    if(r.invalidated&&lab.structure===source?.structure)return;
    source=null;full=null;tools.hidden=true;note.hidden=true;sheet.replaceChildren();mode(false);
    if(!r.structure||!r.viewer?.atoms?.length)return;
    source=r;full=r.viewer;repeats=r.repeats||[1,1,1];topology.checked=preferredOverlay;
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
    note.textContent=descriptions.cell;syncPackingControls('cell');
    if(lab.origin==='catalog-design'){
      viewer.hideHydrogen=true;$('viewer-hide-h').checked=true;
      $('viewer-color').value=ids.length>1?'network':'element';viewer.draw();
      const previous=sceneChoices.get(source.structure);
      if(previous&&[...scene.options].some(o=>o.value===previous.key)){
        scene.value=previous.key;
        if([...$('crystal-net-index').options].some(o=>o.value===previous.network))$('crystal-net-index').value=previous.network;
      }else if(preferSpotlight&&ids.length>1)scene.value='spotlight';
    }
    if(scene.value!=='cell'||preferredOverlay)choose();
  });
  window.addEventListener('cofdraw:display-cells',()=>{
    if(!source)return;
    const key=scene.value, previous=$('crystal-net-index').value;
    repeats=['a','b','c'].map(axis=>Number($('display-'+axis).value));
    full=viewer.model;networks();
    if([...$('crystal-net-index').options].some(o=>o.value===previous))$('crystal-net-index').value=previous;
    choose(true);
  });
  for(const id of ['viewer-axis','viewer-reset','viewer-preset','viewer-color','viewer-palette','viewer-hide-h'])$(id).addEventListener(id==='viewer-reset'?'click':'change',()=>requestAnimationFrame(paint));
})();
