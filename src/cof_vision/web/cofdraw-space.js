'use strict';
document.addEventListener('DOMContentLoaded', () => {
  const $ = id => document.getElementById(id);
  const studio=window.COFdrawStudio;
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
  let linking, views, selectors, mapFamily=null, shownPoints=[], screen, candidates, modelRows, screenPage=0;
  let campaign;
  const chemistryFields={repeat_mass:'Dry repeat mass (g/mol)',nitrogen_wt:'Nitrogen (wt%)',oxygen_wt:'Oxygen (wt%)',junction_density:'Linkage density (mmol/g)',aromatic_c_pct:'Aromatic share of carbon (%)'};
  const crystalFields={density:'Density (g/cm³)',void:'Void fraction (%)',lcd:'Largest cavity, LCD (Å)',pld:'Limiting diameter, PLD (Å)',pore_volume:'Void volume (cm³/g)',free_volume:'Probe-free volume (cm³/g)',connected_volume:'Connected pore volume (cm³/g)',area:'Connected surface (m²/g)',local_area:'Local exposed surface (m²/g)',cell_volume:'Cell volume (Å³)'};
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
  function picture(id){
    const img=element('img');img.src='./api/drawing/space/blocks/'+id+'.svg';img.alt=label(id)+' structure';return img;
  }
  function structureHint(control,id){
    const show=()=>{if(tip.hidden||tip.querySelector('img')?.dataset.block!==String(id)){const img=picture(id);img.dataset.block=id;tip.replaceChildren(element('strong',label(id)),img);}tip.hidden=false;const r=control.getBoundingClientRect();tip.style.left=Math.max(4,Math.min(innerWidth-tip.offsetWidth-4,r.left))+'px';tip.style.top=Math.max(4,Math.min(innerHeight-tip.offsetHeight-4,r.bottom+4))+'px';};
    control.addEventListener('pointerenter',show);control.addEventListener('pointermove',show);control.addEventListener('focus',show);
    for(const event of ['pointerleave','blur','click'])control.addEventListener(event,()=>tip.hidden=true);
  }
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
    candidates=new Map((screen.families[key]||[]).map(row=>[row[0],row]));
    $('screen-set').value='hypothetical';$('screen-query').value='';$('screen-metric').value='junction_density';screenPage=0;
    renderMatrix();renderReports();screenTable();$('space-family-heading').focus({preventScroll:true});window.scrollTo(0,0);
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
    const panel=$('space-pair-detail'),images=element('div');images.className='pair-molecules';
    for(const id of [a,b]){const fig=element('figure');fig.append(picture(id),element('figcaption',label(id)));images.append(fig);}
    panel.replaceChildren(images,element('strong',messages[level]));
    const code=pairCode(a,b),row=candidates.get(code),model=screen.models[family]?.[code];
    if(row){
      panel.append(element('p',row[2]+' · '+row[3]+' g/mol dry repeat · '+row[4]+' wt% N · '+row[6]+' mmol/g junctions'));
      const draw=element('button','Open pair in 2D');draw.className='button outline';
      draw.onclick=async()=>{draw.disabled=true;const ids=[Math.floor(code/data.stride),code%data.stride];if(!row[1])ids.reverse();
        try{await studio.applyDescriptionPlan({status:'ready',prepare_request:{monomer_a:screen.blocks[ids[0]][0],monomer_b:screen.blocks[ids[1]][0],label_a:label(ids[0]).slice(0,100),label_b:label(ids[1]).slice(0,100),linkage:screen.rules[family],template_id:row[2],instances:6}});back();window.scrollTo(0,0);}
        catch(e){panel.append(element('p',e.message));}finally{draw.disabled=false;}};
      panel.append(draw);
      if(model){
        const save=element('button','Save AA model CIF');save.className='button outline';save.onclick=()=>studio.saveText('COFdraw-'+family+'-'+code+'.cif',model.cif,'chemical/x-cif');panel.append(save);
        panel.append(element('p','AA model · '+model.density+' g/cm³ · LCD '+model.lcd+' Å · PLD '+model.pld+' Å · '+model.channel.replaceAll('_',' ')));
      }else{panel.append(element('p',screen.failed[family]?.[code]?'AA fit did not pass structural checks.':row[2]==='chain_2_2'?'Open-chain design.':'Periodic model pending.'));}
    }
    if(p)refs(panel,p);
    campaign?.showPair(panel, family, code);
  }
  function choosePair(a,b){
    const [left,right]=orient(a,b);
    for(const [i,id] of [[0,left],[1,right]]){if(!Array.from(selectors[i].options).some(o=>+o.value===id))selectors[i].add(new Option(label(id),id));selectors[i].value=id;}
    showPair(left,right);
  }
  function axisHeading(id,axis){
    const h=element('th'),button=element('button',label(id));button.className='axis-link';button.setAttribute('aria-label','Browse partners for '+label(id));structureHint(button,id);
    button.onclick=()=>{
      const opposite=1-axis;
      const ids=views[opposite].filter(other=>cellState(id,other)[0]>0);
      selectors[axis].replaceChildren(...views[axis].map(v=>new Option(label(v),v)));
      selectors[axis].value=id;selectors[opposite].replaceChildren(...ids.map(v=>new Option(label(v),v)));
      if(ids.length)showPair(+selectors[0].value,+selectors[1].value);
      selectors[opposite].focus();selectors[opposite].scrollIntoView({block:'nearest'});
    };h.append(button);return h;
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
    for(const b of views[1].slice(0,10))head.append(axisHeading(b,1));
    for(const a of views[0].slice(0,10)){
      const row=body.insertRow();row.append(axisHeading(a,0));
      for(const b of views[1].slice(0,10)){
        const [level,p]=cellState(a,b),button=element('button'),info=label(a)+' × '+label(b)+' · '+messages[level]+(p?' · '+p.refs[0].join(' · '):'');
        button.className='space-cell';button.dataset.level=level;button.setAttribute('aria-label',info);hover(button,info);
        button.onclick=()=>choosePair(a,b);row.insertCell().append(button);
      }
    }
    matrix.append(table);
    selectors=views.map((ids,i)=>{const field=element('label','Linker '+(i?'B':'A')),s=element('select');s.setAttribute('aria-label',field.textContent);for(const id of ids)s.add(new Option(label(id),id));field.append(s);pick.append(field);s.onchange=()=>showPair(...selectors.map(v=>+v.value));return s;});
    showPair(...selectors.map(v=>+v.value));
  }
  function screenTable(){
    const rows=screen.families[family]||[],models=screen.models[family]||{};
    const result=$('screen-results'), hist=$('screen-hist');
    hist.replaceChildren();result.replaceChildren();
    for(const id of ['screen-controls','screen-pages','screen-export'])$(id).hidden=!rows.length;
    if(!rows.length){$('screen-coverage').textContent='This family has no audited pair-construction rule in the current screen.';return;}
    const metric=$('screen-metric').value, col=screen.columns.indexOf(metric), filter=$('screen-set').value, term=$('screen-query').value.trim().toLowerCase();
    const order=$('screen-order').value==='desc'?-1:1, matches=[];
    const blockMatches=term?new Set(Object.keys(data.blocks).filter(id=>label(id).toLowerCase().includes(term)).map(Number)):null;
    let inSet=0;
    for(const row of rows){
      const code=row[0],model=models[code];
      if((filter==='hypothetical'&&linking.reports.has(code))||(filter==='modeled'&&!model)||(filter==='chains'&&row[2]!=='chain_2_2'))continue;
      const a=Math.floor(code/data.stride), b=code%data.stride;
      if(blockMatches&&!blockMatches.has(a)&&!blockMatches.has(b))continue;
      inSet++;
      const v=col>=0?row[col]:model?.[metric];
      if(Number.isFinite(v))matches.push({code,a,b,v});
    }
    matches.sort((a,b)=>order*(a.v-b.v)||a.code-b.code);
    $('screen-coverage').textContent=rows.length.toLocaleString()+' eligible pairs · '+Object.keys(models).length+' screened AA cells · '+matches.length.toLocaleString()+' of '+inSet.toLocaleString()+' filtered pairs have this property';
    const width=6,last=Math.max(0,Math.ceil(matches.length/width)-1);screenPage=Math.min(screenPage,last);
    $('screen-prev').disabled=screenPage===0;$('screen-next').disabled=screenPage===last;
    $('screen-page').textContent=matches.length?'Page '+(screenPage+1)+' of '+(last+1):'';
    if(!matches.length){result.textContent='No computed values for this selection.';return;}
    const lo=Math.min(matches[0].v,matches.at(-1).v),hi=Math.max(matches[0].v,matches.at(-1).v);
    const bins=Array(hi===lo?1:28).fill(0);
    for(const {v} of matches)bins[Math.min(bins.length-1,Math.floor((v-lo)/(hi-lo||1)*bins.length))]++;
    const chart=svg('svg',{viewBox:'0 0 600 53',preserveAspectRatio:'none',role:'img','aria-label':'Distribution of computed values'}),peak=Math.max(...bins);
    bins.forEach((n,i)=>{const h=46*n/peak,w=600/bins.length,r=svg('rect',{x:i*w,y:52-h,width:w-1,height:h,fill:'#368475'});hover(r,'Range '+(lo+i*(hi-lo)/bins.length).toFixed(2)+'–'+(lo+(i+1)*(hi-lo)/bins.length).toFixed(2)+' · '+n+' pairs');chart.append(r);});
    hist.append(chart,element('small',lo.toLocaleString()+' – '+hi.toLocaleString()+' · '+(chemistryFields[metric]||crystalFields[metric])));
    for(const item of matches.slice(screenPage*width,(screenPage+1)*width)){
      const btn=element('button');btn.className='screen-pair';
      btn.append(element('span',label(item.a)+' × '+label(item.b)),element('b',item.v.toLocaleString(undefined,{maximumFractionDigits:4})));
      btn.onclick=()=>{choosePair(item.a,item.b);$('space-pair-detail').scrollIntoView({block:'center'});};result.append(btn);
    }
  }
  for(const id of ['screen-set','screen-metric','screen-order','screen-query'])$(id).addEventListener(id==='screen-query'?'input':'change',()=>{screenPage=0;screenTable();});
  $('screen-prev').onclick=()=>{screenPage--;screenTable();};$('screen-next').onclick=()=>{screenPage++;screenTable();};
  $('screen-export').onclick=async()=>{const b=$('screen-export');b.disabled=true;
    try{const r=await studio.api('space/export/'+encodeURIComponent(family),{});await studio.saveExport(r.url,'zip');}
    catch(e){$('screen-coverage').textContent=e.message;}finally{b.disabled=false;}};
  function scatter(){
    const plot=$('space-scatter'), key=$('space-property').value, horizontal=$('space-property-x').value;
    const group=$('space-model-set').value;
    const rows=modelRows.filter(r=>(group==='all'||r.source===group)&&Number.isFinite(r[key])&&Number.isFinite(r[horizontal]));
    const selected=$('space-model'),prev=selected.value;selected.replaceChildren();
    for(const r of rows)selected.add(new Option(r.name,r.id));if(rows.some(r=>r.id===prev))selected.value=prev;
    $('space-model-count').textContent=rows.length+' periodic models with both properties';
    plot.replaceChildren();
    if(!rows.length){$('space-model-detail').textContent='No models with both values.';return;}
    const max=k=>k==='void'?100:Math.max(.01,...rows.map(r=>r[k]))*1.05;
    const xmax=max(horizontal),ymax=max(key);
    const x=v=>48+v/xmax*360,y=v=>224-v/ymax*197;
    for(let i=0;i<=4;i++){
      const yy=y(i*ymax/4),xx=x(i*xmax/4);
      plot.append(svg('line',{x1:48,y1:yy,x2:408,y2:yy,stroke:'#dce6df'}),svg('text',{x:40,y:yy+4,'text-anchor':'end','font-size':11},Number((i*ymax/4).toPrecision(3))),svg('text',{x:xx,y:241,'text-anchor':'middle','font-size':11},Number((i*xmax/4).toPrecision(3))));
    }
    plot.append(svg('text',{x:225,y:263,'text-anchor':'middle','font-size':12},crystalFields[horizontal]),svg('text',{x:48,y:16,'font-size':12},crystalFields[key]));
    const chosen=$('space-model').value;
    for(const r of rows){
      const dot=svg('circle',{cx:x(r[horizontal]),cy:y(r[key]),r:r.id===chosen?6:4.5,fill:r.id===chosen?'#bd643b':r.source==='hypothesis'?'#7464a0':'#24756c',stroke:'white','stroke-width':1,class:'space-point'});
      hover(dot,r.name+' · '+r[horizontal]+' '+crystalFields[horizontal]+' · '+r[key]+' '+crystalFields[key]);dot.onclick=()=>{$('space-model').value=r.id;tip.hidden=true;scatter();};plot.append(dot);
    }
    const r=rows.find(r=>r.id===chosen);
    if(r){const panel=$('space-model-detail');panel.replaceChildren(element('strong',r.name),element('p',r.source+' · '+r.channel.replaceAll('_',' ')+' · grid '+r.grid.join(' × ')),element('p','PLD bounds '+r.pld_bounds.map(v=>v.toFixed(2)).join('–')+' Å · LCD bounds '+r.lcd_bounds.map(v=>v.toFixed(2)).join('–')+' Å'));
      if(r.doi)refs(panel,{refs:[[r.id,r.doi]]});
      if(r.family){const b=element('button','Inspect this pair');b.className='text-button';b.onclick=()=>{openFamily(r.family,selected);choosePair(Math.floor(r.code/data.stride),r.code%data.stride);$('space-pair-detail').scrollIntoView({block:'center'});};panel.append(b);}
    }
  }
  async function load() {
    [data,screen]=await Promise.all(['chemical-space.json','chemical-screen.json'].map(async file=>{const r=await fetch('./'+file);if(!r.ok)throw Error('Could not load '+file);return r.json();}));
    index=Object.fromEntries(Object.keys(data.families).map(k=>[k,[]]));
    data.points.forEach((p,i)=>{p.i=i;index[p.f].push(p);p.search=[data.families[p.f].name,...p.ids.map(label),...p.refs.flat()].join(' ').toLowerCase();});
    ranked=Object.keys(index).sort((a,b)=>index[b].length-index[a].length||a.localeCompare(b));family=ranked[0];
    colors=Object.fromEntries(ranked.map((k,i)=>[k,colorList[Math.min(i,colorList.length-1)]]));
    const colorsBox=$('space-colors');ranked.slice(0,7).concat(['other']).forEach(k=>{const spot=element('span'),chip=element('i');chip.style.background=colors[k]||colorList.at(-1);spot.append(chip,document.createTextNode(data.families[k]?.name||'Other'));colorsBox.append(spot);});
    const xs=data.points.map(p=>p.xy[0]),ys=data.points.map(p=>p.xy[1]),xmin=Math.min(...xs),xmax=Math.max(...xs),ymin=Math.min(...ys),ymax=Math.max(...ys);
    base={x:(xmin+xmax)/2,y:(ymin+ymax)/2,scale:Math.min((bounds.width-32)/Math.max(1,xmax-xmin),(bounds.height-20)/Math.max(1,ymax-ymin))};
    $('space-content').hidden=false;$('space-status').textContent='';filterPoints();familyCards();
    modelRows=screen.references.slice();
    for(const [f,models] of Object.entries(screen.models))for(const [c,m] of Object.entries(models)){
      const code=Number(c);modelRows.push({...m,id:f+'-'+c,source:'hypothesis',family:f,code,name:data.families[f].name+' · '+label(Math.floor(code/data.stride))+' + '+label(code%data.stride)});
    }
    for(const [k,v] of Object.entries({...chemistryFields,...crystalFields}))$('screen-metric').add(new Option(v,k));
    for(const id of ['space-property','space-property-x'])for(const [k,v] of Object.entries(crystalFields))$(id).add(new Option(v,k));
    $('space-property').value='void';$('space-property-x').value='density';scatter();
    campaign=window.createCampaignExplorer({
      atlas:data, screen, studio, label, resultsChanged:()=>{
        if(!$('space-family-page').hidden&&selectors?.length===2)showPair(...selectors.map(s=>+s.value));
      }, openPair:(f,code)=>{
        openFamily(f,$('campaign-candidate'));choosePair(Math.floor(code/data.stride),code%data.stride);
        $('space-pair-detail').scrollIntoView({block:'center'});
      },
    });
  }
  for(const id of ['space-property','space-property-x','space-model','space-model-set'])$(id).onchange=scatter;
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
