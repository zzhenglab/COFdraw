'use strict';
window.createCampaignExplorer = function ({atlas, screen, studio, label, openPair, resultsChanged}) {
  const $ = id => document.getElementById('campaign-'+id);
  const el = (tag, text) => {const node=document.createElement(tag);node.textContent=text;return node;};
  const format = n => n.toLocaleString(undefined,{maximumSignificantDigits:5});
  const hex = b => Array.from(new Uint8Array(b),v=>v.toString(16).padStart(2,'0')).join('');
  const pairs = new Map(Object.entries(screen.families).map(([family,rows])=>[family,new Set(rows.map(r=>r[0]))]));
  const recordName = row => label(Math.floor(row.pair/atlas.stride))+' + '+label(row.pair%atlas.stride);
  const pageSize=64;
  let snapshot, byId=new Map(), filtered=[], plotted=[], selected, page=0, request=0;
  let pointPositions=[];

  function check(value) {
    if(value?.format!=='cofdraw-screen-results'||value.schema_version!==1||!Array.isArray(value.records)||
       !value.methods?.screen||!value.metrics||!/^[a-f0-9]{64}$/.test(value.protocol)||
       !/^\d{4}-\d\d-\d\dT/.test(value.snapshot_utc)||!Number.isFinite(Date.parse(value.snapshot_utc))||
       value.atlas_sha256!==screen.atlas_sha256)
      throw Error('Choose a COFdraw campaign snapshot built for this precursor atlas.');
    const keys=Object.keys(value.metrics);
    if(!keys.length||keys.length>100||keys.some(k=>!/^\w+$/.test(k)||typeof value.metrics[k]!=='string'))
      throw Error('The result file has invalid property definitions.');
    const records=new Map();
    for(const r of value.records){
      if(!pairs.get(r.family)?.has(r.pair)||r.id!==r.family+'-'+r.pair||records.has(r.id)||!r.values||
         ![r.status,r.template,...(r.cohorts||[]),r.channel??'unknown'].every(v=>typeof v==='string'&&/^[a-zA-Z0-9_]+$/.test(v))||
         typeof r.indexed!=='boolean'||
         Object.entries(r.values).some(([k,v])=>!Object.hasOwn(value.metrics,k)||(v!==null&&!Number.isFinite(v)))||
         (r.sha256!=null&&!/^[a-f0-9]{64}$/.test(r.sha256))||(r.cif_gzip&&!r.sha256))
        throw Error('Invalid or duplicated candidate in the result file.');
      records.set(r.id,r);
    }
    return records;
  }

  function install(value) {
    const records=check(value);
    const family=$('family').value, x=$('x').value, y=$('y').value;
    byId=records;snapshot=value;
    $('family').replaceChildren(new Option('All linkages',''));
    for(const key of [...new Set(value.records.map(r=>r.family))].sort())
      $('family').add(new Option(atlas.families[key].name,key));
    $('family').value=family;if($('family').selectedIndex<0)$('family').value='';
    const metrics=Object.entries(value.metrics);
    for(const [id,previous,preferred] of [['x',x,'delivery_cm3_STP_cm3'],['y',y,'dilute_CH4_H2']]){
      $(id).replaceChildren(...metrics.map(([key,title])=>new Option(title,key)));
      $(id).value=Object.hasOwn(value.metrics,previous)?previous:Object.hasOwn(value.metrics,preferred)?preferred:metrics[0][0];
    }
    $('body').hidden=false;
    const metadata={...value};delete metadata.records;
    $('methods').textContent=JSON.stringify(metadata,null,2);
    $('status').textContent=value.records.length.toLocaleString()+' candidate outcomes · '+new Date(value.snapshot_utc).toLocaleString();
    update();resultsChanged();
  }

  async function load(file) {
    const current=++request;
    try{
      let payload;
      if(file){
        const stream=file.name.toLowerCase().endsWith('.gz')?file.stream().pipeThrough(new DecompressionStream('gzip')):file.stream();
        payload=await new Response(stream).json();
      }else{
        const response=await fetch('./api/drawing/space/screen-results',{cache:'no-store'});
        if(current!==request)return;
        if(response.status===204){
          $('status').textContent=snapshot?'No locally installed snapshot. The open file is still displayed.':'Open a screening snapshot (.json.gz) to explore calculated properties and download its periodic structures.';
          return;
        }
        if(!response.ok)throw Error('The locally installed screening file could not be loaded.');
        payload=await response.json();
      }
      if(current===request)install(payload);
    }catch(e){if(current===request)$('status').textContent='Could not open results. '+e.message;}
  }

  async function saveCif(r, status) {
    try{
      const bytes=Uint8Array.from(atob(r.cif_gzip),c=>c.charCodeAt(0));
      const plain=await new Response(new Blob([bytes]).stream().pipeThrough(new DecompressionStream('gzip'))).arrayBuffer();
      if(hex(await crypto.subtle.digest('SHA-256',plain))!==r.sha256)throw Error('CIF checksum mismatch.');
      await studio.saveText(r.id+'.cif',new TextDecoder().decode(plain),'chemical/x-cif');
    }catch(e){status.textContent=e.message;}
  }
  function details(panel,r,link) {
    panel.append(el('strong',r.id+' · '+recordName(r)),el('p',atlas.families[r.family].name+' · '+r.template+' · '+r.status.replaceAll('_',' ')));
    if(r.channel)panel.append(el('p',r.channel.replaceAll('_',' ')+' · '+(r.indexed?'Same-family precursor precedent':'No indexed same-family precursor precedent')));
    if(r.reason)panel.append(el('p',r.reason));
    if(r.cohorts)panel.append(el('p','Cohort · '+r.cohorts.join(', ')));
    for(const [stage,title] of [['confirmation','Confirmation'],['deep','Extended run']]){
      const seeds=r[stage+'_seeds'];
      if(seeds)panel.append(el('p',title+' seeds · '+seeds.delivery+' paired CH₄ · '+seeds.mixture+' mixture'));
    }
    if(r.pld_bounds&&r.lcd_bounds)panel.append(el('p','Grid bounds · PLD '+r.pld_bounds.map(format).join('–')+' Å · LCD '+r.lcd_bounds.map(format).join('–')+' Å'));
    const dl=document.createElement('dl');
    for(const [k,title] of Object.entries(snapshot.metrics))
      if(Number.isFinite(r.values[k]))dl.append(el('dt',title),el('dd',format(r.values[k])));
    if(dl.childElementCount)panel.append(dl);
    const note=el('p','');panel.append(note);
    if(r.cif_gzip){
      const b=el('button','Save simulated CIF');b.className='button outline';b.onclick=()=>saveCif(r,note);panel.append(b);
    }
    if(link){
      const b=el('button','Inspect precursor pair');b.className='text-button';b.onclick=()=>openPair(r.family,r.pair);panel.append(b);
    }
  }

  function draw() {
    const canvas=$('plot'),ctx=canvas.getContext('2d'),w=canvas.width,h=canvas.height;
    const xKey=$('x').value,yKey=$('y').value;
    ctx.clearRect(0,0,w,h);pointPositions=[];
    if(!plotted.length){ctx.fillStyle='#44625e';ctx.font='16px sans-serif';ctx.fillText('No candidates with both properties.',58,52);return;}
    function bounds(key){
      let lo=Infinity,hi=-Infinity;
      for(const r of plotted){lo=Math.min(lo,r.values[key]);hi=Math.max(hi,r.values[key]);}
      const margin=(hi-lo||Math.abs(lo)||1)*.04;
      return [lo-margin,hi+margin];
    }
    const [xmin,xmax]=bounds(xKey),[ymin,ymax]=bounds(yKey),left=68,right=w-14,top=12,bottom=h-36;
    ctx.font='14px sans-serif';ctx.fillStyle='#36544f';ctx.strokeStyle='#dce7e2';ctx.lineWidth=1;
    for(let i=0;i<=4;i++){
      const y=bottom-i*(bottom-top)/4,x=left+i*(right-left)/4;
      ctx.beginPath();ctx.moveTo(left,y);ctx.lineTo(right,y);ctx.stroke();ctx.textAlign='right';
      ctx.fillText(Number((ymin+(ymax-ymin)*i/4).toPrecision(3)),left-6,y+4);
      ctx.textAlign='center';ctx.fillText(Number((xmin+(xmax-xmin)*i/4).toPrecision(3)),x,bottom+19);
    }
    ctx.fillStyle='#247c7080';
    for(const r of plotted){
      const x=left+(r.values[xKey]-xmin)/(xmax-xmin)*(right-left), y=bottom-(r.values[yKey]-ymin)/(ymax-ymin)*(bottom-top);
      pointPositions.push([x,y,r]);ctx.beginPath();ctx.arc(x,y,3,0,Math.PI*2);ctx.fill();
    }
    const pt=pointPositions.find(p=>p[2].id===selected);
    if(pt){ctx.strokeStyle='#984425';ctx.lineWidth=2;ctx.beginPath();ctx.arc(pt[0],pt[1],6,0,Math.PI*2);ctx.stroke();}
  }

  function select(id) {
    selected=id;
    const panel=$('detail');panel.replaceChildren();
    const r=byId.get(id);if(r)details(panel,r,true);
    const list=$('candidate'),offset=page*pageSize;
    list.replaceChildren(...filtered.slice(offset,offset+pageSize).map(row=>new Option(row.id+' · '+recordName(row),row.id)));
    list.value=id||'';
    $('prev').disabled=page===0;$('next').disabled=offset+pageSize>=filtered.length;
    draw();
  }
  function update() {
    const family=$('family').value,term=$('query').value.trim().toLowerCase(),x=$('x').value,y=$('y').value;
    filtered=snapshot.records.filter(r=>(!family||r.family===family)&&(!term||r.id.includes(term)||recordName(r).toLowerCase().includes(term)));
    filtered.sort((a,b)=>(b.values[y]??-Infinity)-(a.values[y]??-Infinity)||a.id.localeCompare(b.id));
    plotted=filtered.filter(r=>Number.isFinite(r.values[x])&&Number.isFinite(r.values[y]));
    $('count').textContent=plotted.length.toLocaleString()+' of '+filtered.length.toLocaleString()+' filtered candidates have both properties. Candidate list ordered by '+snapshot.metrics[y]+'.';
    $('export').disabled=!filtered.length;
    const index=filtered.findIndex(r=>r.id===selected);page=Math.floor(Math.max(0,index)/pageSize);
    select(index>=0?selected:filtered[0]?.id);
  }
  $('plot').onclick=e=>{
    const box=$('plot').getBoundingClientRect(),x=(e.clientX-box.left)*$('plot').width/box.width,y=(e.clientY-box.top)*$('plot').height/box.height;
    let best,dist=144;for(const [px,py,r] of pointPositions){const d=(x-px)**2+(y-py)**2;if(d<dist){dist=d;best=r;}}
    if(best){page=Math.floor(filtered.indexOf(best)/pageSize);select(best.id);}
  };
  $('candidate').onchange=()=>select($('candidate').value);
  for(const [id,delta] of [['prev',-1],['next',1]])$(id).onclick=()=>{page+=delta;select(filtered[page*pageSize].id);};
  for(const id of ['family','x','y'])$(id).onchange=update;
  $('query').oninput=update;
  $('refresh').onclick=()=>load();
  $('import').onclick=()=>$('file').click();
  $('file').onchange=()=>{const file=$('file').files[0];$('file').value='';if(file)load(file);};
  $('methods-save').onclick=()=>studio.saveText('cofdraw-screen-methods.json',$('methods').textContent).catch(e=>$('status').textContent=e.message);
  $('export').onclick=async()=>{
    const fields=['id','family','pair','template','status','indexed','sha256','channel'];
    const keys=Object.keys(snapshot.metrics);
    const cell=v=>'"'+String(v??'').replaceAll('"','""')+'"';
    const lines=[[...fields,'cohorts','delivery_confirmation_seeds','mixture_confirmation_seeds','delivery_deep_seeds','mixture_deep_seeds','protocol','snapshot_utc',...keys].join(',')];
    for(const r of filtered)lines.push([...fields.map(k=>r[k]),r.cohorts?.join('|'),r.confirmation_seeds?.delivery,r.confirmation_seeds?.mixture,r.deep_seeds?.delivery,r.deep_seeds?.mixture,snapshot.protocol,snapshot.snapshot_utc,...keys.map(k=>r.values[k])].map(cell).join(','));
    try{await studio.saveText('cofdraw-campaign-results.csv',lines.join('\r\n')+'\r\n','text/csv');}
    catch(e){$('status').textContent=e.message;}
  };
  load();
  return {showPair(panel,family,code){
    const r=byId.get(family+'-'+code);
    if(r){const more=document.createElement('details');more.append(el('summary','Simulation results · '+r.status.replaceAll('_',' ')));details(more,r,false);panel.append(more);}
  }};
};
