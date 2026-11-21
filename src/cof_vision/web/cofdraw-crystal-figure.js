/* Coordinate-based publication figures. All calculations retain the parent CIF. */
(()=>{
  'use strict';
  const $=id=>document.getElementById(id), lab=window.COFdrawScience, viewer=lab.viewer;
  const scene=$('crystal-scene'), tools=$('crystal-figure-tools'), sheet=$('crystal-sheet');
  const toggle=$('crystal-figure-toggle'), save=$('crystal-figure-save'), note=$('crystal-scene-note');
  const colors={C:'#465568',N:'#3576b6',O:'#c95550',S:'#b78d37',B:'#a67e64',Si:'#808597',F:'#398f77'};
  const palette=['#166d83','#bd654a','#8071ac','#a4882e','#39816b','#b45d8a','#5077ac','#998766','#727b4e'];
  const labels={ring:'Molecular circuit',junction:'Tetrahedral junction',net:'Junction net',overlay:'Structure + net',single:'One network',cell:'Crystal packing'};
  labels.spotlight='Highlight a network';
  const descriptions={ring:'Closed circuit traced through a single covalent network. All heavy atoms along its linkers are shown.',
    junction:'Four bonded shells around one four-connected center, unwrapped across cell faces.',
    net:'Each segment follows a traced linker between two four-connected centers. Two junction shells are shown.',
    overlay:'Dashed guides connect junction centers along four traced linker paths.',
    single:'One distance-connected periodic network within the displayed cell boundaries.',
    cell:'Networks are colored by connectivity in the selected coordinate model.'};
  descriptions.spotlight='Selected network in color with the surrounding crystal faded.';
  let source=null,full=null,isFigure=false;

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
    if(key==='single'){
      model={...part(full,a=>a.network===selected),cell:full.cell};
    }else if(key==='spotlight'){
      model={...full,highlightNetwork:selected};
    }else if(key!=='cell')model=source.crystal_details[key];
    viewer.set(key==='cell'||key==='spotlight'?model:{...model,illustration:true});
    const direction=model.direction;
    if(direction){
      viewer.theta=Math.atan2(-direction[0],direction[2]);
      viewer.phi=Math.asin(Math.min(1,Math.max(-1,direction[1])));
    }
    $('viewer-axis').value='free';
    $('viewer-color').value=['cell','single','spotlight'].includes(key)?'network':'element';
    viewer.hideHydrogen=true;$('viewer-hide-h').checked=true;
    note.textContent=descriptions[key]+' Calculations use the full crystal.';
    window.dispatchEvent(new Event('cofdraw:crystal-scene'));
    viewer.draw();paint();
  }

  function node(tag,attrs,text){
    const e=document.createElementNS('http://www.w3.org/2000/svg',tag);
    for(const [key,val] of Object.entries(attrs))e.setAttribute(key,val);
    if(text!==undefined)e.textContent=text;
    return e;
  }
  function paint(){
    if(!isFigure||!source||!viewer.atoms.length)return;
    const atoms=viewer.atoms, key=scene.value, net=key==='net', packing=key==='cell';
    const name=source.structure.provenance.model_name||'Periodic crystal';
    const svg=node('svg',{xmlns:'http://www.w3.org/2000/svg',viewBox:'0 0 760 570',role:'img','aria-label':name+' · '+labels[key]});
    svg.append(node('title',{},name+' · '+labels[key]),node('desc',{},note.textContent+' DOI '+(source.structure.provenance.doi||'')),
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
      node('text',{x:35,y:543,fill:'#556573','font-family':'sans-serif','font-size':13},'Orthographic projection · '+(source.structure.provenance.doi||'Selected coordinate model')));
    sheet.replaceChildren(svg);
  }
  function mode(active){
    isFigure=active;sheet.hidden=!active;save.hidden=!active;
    $('model-view').classList.toggle('figure-active',active);
    toggle.textContent=active?'Rotate in 3D':'2D figure';toggle.setAttribute('aria-pressed',String(active));
    if(!$('science-dialog').hidden)$('switch-space').textContent='2D drawing';
    if(active)paint();else viewer.draw();
  }
  scene.addEventListener('change',choose);
  $('crystal-net-index').addEventListener('change',choose);
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
    if(!r.structure||!(r.color_periodic_networks||r.structure.provenance?.network_dimension===3))return;
    source=r;full=r.viewer;
    scene.replaceChildren();
    const ids=networks();
    for(const key of ['cell','single','spotlight','ring','junction','net','overlay']){
      if(key==='cell'||(key==='single'?ids.length:key==='spotlight'?ids.length>1:r.crystal_details?.[key]))
        scene.add(new Option(labels[key],key));
    }
    tools.hidden=false;note.hidden=false;
    $('crystal-annotation').hidden=r.structure.provenance.template!=='cof-300-c5';
    note.textContent=descriptions.cell;
    if(lab.origin==='catalog-design'){
      viewer.hideHydrogen=true;$('viewer-hide-h').checked=true;
      $('viewer-color').value=ids.length>1?'network':'element';viewer.draw();
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
