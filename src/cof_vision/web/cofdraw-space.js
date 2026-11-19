'use strict';
document.addEventListener('DOMContentLoaded', () => {
  const $ = id => document.getElementById(id);
  const studio=window.COFdrawStudio, repeatDialog=$('repeat-view');
  let repeatRequest=0, repeatFile;
  async function updateRepeat() {
    const generation=++repeatRequest, state=studio.state, count=+$('repeat-count').value;
    const image=$('repeat-image'), message=$('repeat-message'), save=$('repeat-save');
    repeatFile=null;save.disabled=true;image.hidden=true;image.removeAttribute('src');
    if(!state.recipe||!state.result||state.inputsDirty||state.jsonPending||state.busy||state.dirty){
      message.textContent='Render the current drawing to inspect its repeat units.';return;
    }
    const template=state.recipe.selected_template?.id;
    const supports=!state.sourceId&&!state.recipe.assembly&&['hcb_3_3','hcb_3_2','sql_4_2','sql_4_4','chain_2_2'].includes(template);
    for(const option of $('repeat-count').options)if(+option.value>1)option.disabled=!supports;
    if(count===1){
      repeatFile=state.result.downloads.svg;image.src='.'+repeatFile;image.hidden=false;save.disabled=false;
      message.textContent=supports?'Current drawing. Choose a larger extent to join neighboring blocks.':'Original drawing. To expand it, choose a design with declared chain, square or hexagonal connectivity.';return;
    }
    message.textContent='Fitting connected repeat units…';
    try {
      const r=await studio.api('repeat',{recipe:state.recipe,repeats:count,size:+$('image-size').value});
      if(generation!==repeatRequest||!repeatDialog.open)return;
      const s=r.summary;repeatFile=r.svg;image.src='.'+r.svg;image.hidden=false;save.disabled=false;
      message.textContent=s.extent+' · '+s.blocks+' building blocks · '+s.continuation_sites+' boundary sites. Standard skeletal drawing of the selected design.';
    }catch(error){if(generation===repeatRequest)message.textContent=error.message;}
  }
  $('repeat-open').onclick=()=>{$('repeat-count').value='1';repeatDialog.showModal();updateRepeat();};
  $('repeat-count').onchange=updateRepeat;
  $('repeat-save').onclick=()=>{if(repeatFile)studio.saveExport(repeatFile,'svg');};
  $('repeat-close').onclick=()=>repeatDialog.close();
  repeatDialog.addEventListener('close',()=>{repeatRequest++;});
  for(const event of ['cofdraw:layoutchanged','cofdraw:rendered'])window.addEventListener(event,()=>{if(repeatDialog.open)repeatDialog.close();});
  const ns = 'http://www.w3.org/2000/svg';
  const element = (name, text) => { const el = document.createElement(name); if (text !== undefined) el.textContent = text; return el; };
  function svg(name, attrs, text) {
    const node = document.createElementNS(ns, name);
    for (const [k,v] of Object.entries(attrs)) node.setAttribute(k,v);
    if (text !== undefined) node.textContent=text;
    return node;
  }
  const colorList=['#087b72','#c16b3d','#5d6fa9','#9567a5','#948031','#5e89a4','#974f65','#648946'];
  let data, loading, family, ranked, index, colors, zoom=1, pan=[0,0], base, selectedPoint;
  let linking, views, selectors, mapFamily=null, shownPoints=[];
  const tip=element('div'); tip.className='space-tip'; tip.hidden=true; document.body.append(tip);
  function hover(node, message) {
    const move=e=>{ tip.textContent=message; tip.hidden=false; tip.style.left=Math.max(4,Math.min(innerWidth-tip.offsetWidth-8,e.clientX+12))+'px'; tip.style.top=Math.max(4,Math.min(innerHeight-tip.offsetHeight-8,e.clientY+12))+'px'; };
    node.addEventListener('pointerenter',move);node.addEventListener('pointermove',move);
    node.addEventListener('pointerleave',()=>tip.hidden=true);
  }
  function refs(panel, point) {
    for(const [name,doi] of point.refs) {
      const line=element('div',name);
      const link=element('a',doi); link.href='https://doi.org/'+doi;link.target='_blank';link.rel='noopener noreferrer';
      line.append(link);panel.append(line);
    }
  }
  const label=id=>data.blocks[id]||('Block '+id);
  const pairCode=(a,b)=>Math.min(a,b)*data.stride+Math.max(a,b);
  function openPoint(point) {
    selectedPoint=point;
    tip.hidden=true;
    const panel=$('space-point-detail'); panel.hidden=false;
    panel.replaceChildren(element('strong',data.families[point.f].name+' · '+point.ids.map(label).join(' + ')));
    refs(panel,point);$('space-hit').value=String(point.i);
    drawPoints();
  }
  const map=$('space-map');
  const bounds=map.viewBox.baseVal;
  function drawPoints() {
    map.replaceChildren();
    const group=svg('g',{}); map.append(group);
    for(const p of shownPoints) {
      const x=bounds.width/2+(p.xy[0]-base.x)*base.scale*zoom+pan[0], y=bounds.height/2-(p.xy[1]-base.y)*base.scale*zoom+pan[1];
      if(x < -5 || x>bounds.width+5 || y < -5 || y>bounds.height+5)continue;
      const dot=svg('circle',{cx:x,cy:y,r:p===selectedPoint?5:2.6,fill:colors[p.f],opacity:p===selectedPoint?1:.78,class:'space-point',stroke:p===selectedPoint?'#13292f':'none'});
      hover(dot,p.refs[0][0]+' · '+data.families[p.f].name+' · '+p.refs[0][1]);
      dot.addEventListener('click',()=>{if(!moved)openPoint(p);});group.append(dot);
    }
  }
  function filterPoints() {
    const query=$('space-search').value.trim().toLowerCase();
    shownPoints=data.points.filter(p=>(!mapFamily||p.f===mapFamily)&&(!query||p.search.includes(query)));
    $('space-map-summary').textContent=shownPoints.length.toLocaleString()+' indexed compositions · '+Object.keys(data.blocks).length.toLocaleString()+' building blocks in the merged inventory';
    const select=$('space-hit');select.replaceChildren(new Option('Inspect a composition…',''));
    for(const p of shownPoints)select.add(new Option(p.refs[0][0]+' · '+data.families[p.f].name,p.i));
    if(!shownPoints.includes(selectedPoint)){$('space-point-detail').hidden=true;selectedPoint=null;}
    if(selectedPoint)select.value=String(selectedPoint.i);
    drawPoints();
  }
  const pointer=e=>{const p=map.createSVGPoint();p.x=e.clientX;p.y=e.clientY;return p.matrixTransform(map.getScreenCTM().inverse());};
  let start, moved=false;
  map.addEventListener('pointerdown',e=>{moved=false;start=null;if(e.pointerType!=='mouse'||e.button!==0)return;const p=pointer(e);start=[p.x,p.y,pan[0],pan[1]];});
  map.addEventListener('pointermove',e=>{if(!start||!e.buttons)return;const p=pointer(e),dx=p.x-start[0],dy=p.y-start[1];if(Math.hypot(dx,dy)<3&&!moved)return;moved=true;map.setPointerCapture(e.pointerId);pan=[start[2]+dx,start[3]+dy];tip.hidden=true;drawPoints();});
  window.addEventListener('pointerup',()=>start=null);
  map.addEventListener('pointercancel',()=>start=null);
  function magnify(factor){const next=Math.max(.6,Math.min(12,zoom*factor)),ratio=next/zoom;pan=pan.map(v=>v*ratio);zoom=next;tip.hidden=true;drawPoints();}
  $('space-plus').onclick=()=>magnify(1.5);$('space-minus').onclick=()=>magnify(1/1.5);$('space-fit').onclick=()=>{zoom=1;pan=[0,0];drawPoints();};
  $('space-search').oninput=()=>{mapFamily=null;filterPoints();};
  $('space-hit').onchange=e=>{if(e.target.value!=='')openPoint(data.points[+e.target.value]);else{selectedPoint=null;$('space-point-detail').hidden=true;drawPoints();}};
  let galleryScroll=0, galleryFocus=null, reportLimit=6, returnToMap=false, enteredFromGallery=false;
  function figure(key){
    const img=new DOMParser().parseFromString(data.families[key].svg,'image/svg+xml').documentElement;
    img.setAttribute('aria-hidden','true');
    return document.importNode(img,true);
  }
  function renderReports(){
    const query=$('space-reports-search').value.trim().toLowerCase();
    const matches=index[family].filter(p=>!query||p.search.includes(query));
    const root=$('space-reports-list');root.replaceChildren();
    for(const point of matches.slice(0,reportLimit)){
      const item=element('li'),matched=point.refs.filter(ref=>query&&ref.join(' ').toLowerCase().includes(query));
      const ordered=matched.length?[...matched,...point.refs.filter(r=>!matched.includes(r))]:point.refs;
      item.append(element('strong',ordered[0][0]),element('p',point.ids.map(label).join(' + ')));
      const extra=element('details'),summary=element('summary',ordered.length===1?'Source record':'Source records ('+ordered.length+')');
      extra.append(summary);refs(extra,{refs:ordered});item.append(extra);root.append(item);
    }
    $('space-reports-count').textContent=matches.length?Math.min(reportLimit,matches.length)+' of '+matches.length.toLocaleString()+' compositions':'No matching compositions in the index.';
    $('space-more-reports').hidden=matches.length<=reportLimit;
  }
  $('space-more-reports').onclick=()=>{reportLimit+=12;renderReports();};
  $('space-reports-search').oninput=()=>{reportLimit=6;renderReports();};
  function detail(key){
    family=key;tip.hidden=true;
    const f=data.families[key];$('space-family-heading').textContent=f.name;
    $('space-family-note').textContent=f.illustration_note;
    const points=index[key], papers=new Set(points.flatMap(p=>p.refs.map(r=>r[1])));
    $('space-family-coverage').textContent=points.length.toLocaleString()+' indexed compositions · '+papers.size.toLocaleString()+' DOIs';
    $('space-family-figure').replaceChildren(figure(key));
    $('space-overview').hidden=true;$('space-family-page').hidden=false;
    $('space-reports-search').value='';reportLimit=6;
    renderMatrix();renderReports();$('space-family-heading').focus({preventScroll:true});window.scrollTo(0,0);
  }
  function overview(){
    tip.hidden=true;$('space-overview').hidden=false;$('space-family-page').hidden=true;$('space-family-menu').value='';
    if(returnToMap){
      mapFamily=family;$('space-search').value=data.families[family].name;zoom=1;pan=[0,0];filterPoints();
      window.scrollTo(0,0);$('space-search').focus({preventScroll:true});returnToMap=false;
    }else{window.scrollTo(0,galleryScroll);galleryFocus?.focus({preventScroll:true});}
  }
  function openFamily(key,origin){
    galleryScroll=window.scrollY;galleryFocus=origin;
    enteredFromGallery=true;history.pushState({cofLinkage:key},'','#linkage/'+encodeURIComponent(key));detail(key);
  }
  function leaveFamily(){
    if(enteredFromGallery&&history.state?.cofLinkage)history.back();
    else {history.replaceState(null,'','#chemical-space');overview();}
  }
  $('space-family-back').onclick=leaveFamily;
  $('space-map-family').onclick=()=>{returnToMap=true;leaveFamily();};
  function familyCards() {
    const root=$('space-families');root.replaceChildren();
    for(const key of ranked) {
      const f=data.families[key], b=element('button');b.className='space-family';
      const meta=element('small',index[key].length+' indexed');meta.title='Indexed compositions';
      const caption=element('span');caption.className='space-card-caption';caption.append(element('strong',f.name),meta);
      b.append(figure(key),caption);
      $('space-family-menu').add(new Option(f.name,key));
      b.onclick=()=>openFamily(key,b);root.append(b);
    }
    $('space-family-menu').onchange=e=>{if(e.target.value)openFamily(e.target.value,e.target);};
  }
  function orient(a,b) {
    // The polar partner occupies the same axis throughout a functional-class matrix.
    const priority=['imine','beta-ketoenamine'].includes(family)?'amine':['boronate-ester'].includes(family)?'boron':'aldehyde';
    const score=id=>data.roles[id].includes(priority)?1:0;
    return score(a)>score(b)||score(a)===score(b)&&a<b?[a,b]:[b,a];
  }
  function cellState(a,b) {
    const key=pairCode(a,b),p=linking.reports.get(key),allowed=linking.eligible.has(key);
    return [p?allowed?2:3:allowed?1:0,p];
  }
  const messages=['No match in this index','Eligible drawing · no indexed report','Reported and eligible','Reported · outside screened drawing rules'];
  function showPair(a,b) {
    const [level,p]=cellState(a,b);
    const panel=$('space-pair-detail');panel.replaceChildren(element('strong',label(a)+' × '+label(b)),element('div',messages[level]));
    if(p)refs(panel,p);
  }
  function renderMatrix() {
    const f=data.families[family], points=index[family], pairs=points.filter(p=>p.ids.length===2);
    linking={reports:new Map(pairs.map(p=>[pairCode(...p.ids),p])),eligible:new Set(f.eligible)};
    const counts=[new Map(),new Map()];
    function track(ids,weight) {orient(...ids).forEach((id,i)=>counts[i].set(id,(counts[i].get(id)||0)+weight));}
    for(const p of pairs)track(p.ids,1);
    for(const key of f.eligible)track([Math.floor(key/data.stride),key%data.stride],0);
    views=counts.map(m=>Array.from(m.keys()).sort((a,b)=>m.get(b)-m.get(a)||label(a).localeCompare(label(b))||a-b));
    $('space-heat-title').textContent='Building-block pairs';
    const metrics=$('space-metrics');metrics.replaceChildren();
    const screened=f.eligible.length>0;
    for(const [n,title] of [[screened?f.eligible.length:'—',screened?'eligible pairs':'eligibility unscreened'],[pairs.length,'reported pairs'],[screened?pairs.filter(p=>linking.eligible.has(pairCode(...p.ids))).length:'—','reported and eligible']]){
      const v=element('div');v.append(element('b',n.toLocaleString()),element('small',title));metrics.append(v);
    }
    const matrix=$('space-matrix'),pick=$('space-pickers');matrix.replaceChildren();pick.replaceChildren();
    if(!pairs.length&&!f.eligible.length){$('space-pair-detail').textContent='No binary compositions indexed for this linkage. Single- and multi-precursor compositions appear in the map.';return;}
    const table=element('table'),head=table.createTHead().insertRow(),body=table.createTBody();head.append(element('th','A × B'));
    for(const b of views[1].slice(0,10)){const h=element('th',label(b));h.title=label(b);head.append(h);}
    for(const a of views[0].slice(0,10)){
      const row=body.insertRow(),h=element('th',label(a));h.title=label(a);row.append(h);
      for(const b of views[1].slice(0,10)){
        const [level,p]=cellState(a,b),button=element('button'),info=label(a)+' × '+label(b)+' · '+messages[level]+(p?' · '+p.refs[0].join(' · '):'');
        button.className='space-cell';button.dataset.level=level;button.setAttribute('aria-label',info);hover(button,info);
        button.onclick=()=>{selectors[0].value=a;selectors[1].value=b;showPair(a,b);};row.insertCell().append(button);
      }
    }
    matrix.append(table);
    selectors=views.map((ids,i)=>{const field=element('label','Linker '+(i?'B':'A')),s=element('select');s.setAttribute('aria-label',field.textContent);for(const id of ids)s.add(new Option(label(id),id));field.append(s);pick.append(field);s.onchange=()=>showPair(...selectors.map(v=>+v.value));return s;});
    showPair(...selectors.map(v=>+v.value));
  }
  function scatter() {
    const plot=$('space-scatter'), key=$('space-property').value, rows=data.properties;
    plot.replaceChildren();
    const xmax=Math.max(...rows.map(r=>r.density))*1.05,ymax=key==='void'?100:Math.ceil(Math.max(...rows.map(r=>r[key]))/5)*5;
    const x=v=>48+v/xmax*360,y=v=>224-v/ymax*197;
    for(let i=0;i<=4;i++){
      const yy=y(i*ymax/4),xx=x(i*xmax/4);
      plot.append(svg('line',{x1:48,y1:yy,x2:408,y2:yy,stroke:'#dce6df'}),svg('text',{x:40,y:yy+4,'text-anchor':'end','font-size':11},Math.round(i*ymax/4)),svg('text',{x:xx,y:241,'text-anchor':'middle','font-size':11},(i*xmax/4).toFixed(2)));
    }
    plot.append(svg('text',{x:225,y:263,'text-anchor':'middle','font-size':12},'Model density (g/cm³)'),svg('text',{x:48,y:16,'font-size':12},key==='void'?'Void fraction (%)':'Largest sampled cavity (Å)'));
    const chosen=$('space-model').value;
    for(const r of rows){
      const dot=svg('circle',{cx:x(r.density),cy:y(r[key]),r:r.id===chosen?6:4.5,fill:r.id===chosen?'#bd643b':'#24756c',stroke:'white','stroke-width':1,class:'space-point'});
      hover(dot,r.name+' · '+r.density+' g/cm³ · '+r[key]+(key==='void'?'%':' Å'));dot.onclick=()=>{$('space-model').value=r.id;scatter();};plot.append(dot);
    }
    const r=rows.find(r=>r.id===chosen);
    if(r){const panel=$('space-model-detail');panel.replaceChildren(element('strong',r.name),element('div',r.density+' g/cm³ · '+r.void+'% void · '+r.diameter+' Å sampled cavity · '+r.samples+' samples'));refs(panel,{refs:[[r.id,r.doi]]});}
  }
  async function load() {
    const response=await fetch('./chemical-space.json');if(!response.ok)throw Error('Chemical-space inventory could not be loaded.');
    data=await response.json(); index=Object.fromEntries(Object.keys(data.families).map(k=>[k,[]]));
    data.points.forEach((p,i)=>{p.i=i;index[p.f].push(p);p.search=[data.families[p.f].name,...p.ids.map(label),...p.refs.flat()].join(' ').toLowerCase();});
    ranked=Object.keys(index).sort((a,b)=>index[b].length-index[a].length||a.localeCompare(b));family=ranked[0];
    colors=Object.fromEntries(ranked.map((k,i)=>[k,colorList[Math.min(i,colorList.length-1)]]));
    const colorsBox=$('space-colors');ranked.slice(0,7).concat(['other']).forEach(k=>{const spot=element('span'),chip=element('i');chip.style.background=colors[k]||colorList.at(-1);spot.append(chip,document.createTextNode(data.families[k]?.name||'Other'));colorsBox.append(spot);});
    const xs=data.points.map(p=>p.xy[0]),ys=data.points.map(p=>p.xy[1]),xmin=Math.min(...xs),xmax=Math.max(...xs),ymin=Math.min(...ys),ymax=Math.max(...ys);
    base={x:(xmin+xmax)/2,y:(ymin+ymax)/2,scale:Math.min((bounds.width-32)/Math.max(1,xmax-xmin),(bounds.height-20)/Math.max(1,ymax-ymin))};
    $('space-content').hidden=false;$('space-status').textContent='';filterPoints();familyCards();
    for(const r of data.properties)$('space-model').add(new Option(r.name+' · '+r.id,r.id));
    $('space-model-count').textContent=data.properties.length+' bundled CIFs · precomputed geometric screening';scatter();
  }
  $('space-property').onchange=scatter;$('space-model').onchange=scatter;
  function clearRoute(){
    tip.hidden=true;returnToMap=false;
    if(location.hash==='#chemical-space'||location.hash.startsWith('#linkage/'))history.replaceState(null,'',location.pathname+location.search);
  }
  function back(){clearRoute();$('space-workspace').hidden=true;$('space-open').setAttribute('aria-pressed','false');document.querySelector('main.studio').hidden=false; }
  window.addEventListener('cofdraw:leave-explorer',clearRoute);
  $('space-back').onclick=()=>{back();$('space-open').focus();};
  async function openExplorer(hash){
    if(!$('science-dialog').hidden)$('close-science').click();
    document.querySelector('main.studio').hidden=true;$('space-workspace').hidden=false;$('space-open').setAttribute('aria-pressed','true');window.scrollTo(0,0);
    if(!loading){$('space-status').textContent='Loading the chemical-space inventory…'; loading=load().catch(e=>{$('space-status').textContent=e.message;loading=null;});}await loading;
    if(!data||$('space-workspace').hidden)return;
    const key=ranked.find(k=>hash==='#linkage/'+encodeURIComponent(k));
    if(key){detail(key);return;}
    history.replaceState(null,'','#chemical-space');overview();
  }
  $('space-open').onclick=()=>{
    if(!$('space-workspace').hidden){back();return;}
    galleryFocus=null;galleryScroll=0;enteredFromGallery=false;openExplorer('#chemical-space');
  };
  function route(){
    const hash=location.hash;
    if(hash==='#chemical-space'||hash.startsWith('#linkage/'))openExplorer(hash);
    else if(!$('space-workspace').hidden)back();
  }
  window.addEventListener('popstate',route);route();
});
