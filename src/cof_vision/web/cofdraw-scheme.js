/* Assemble a publication scheme from the exact visible product and its precursor graphs. */
(() => {
  'use strict';
  const studio=window.COFdrawStudio, ns='http://www.w3.org/2000/svg';
  const byId=id=>document.getElementById(id);
  const figure=document.createElement('img');
  figure.className='reaction-figure';
  figure.alt='Precursor structures and reaction arrow leading to the current framework';
  let enabled=false, direction='horizontal', crystalView=null, revision=0, graphic=null;
  const toggles=[...document.querySelectorAll('.scheme-toggle')];
  const directions=[...document.querySelectorAll('.scheme-direction')];

  let nameEnabled=false, customName=null;
  const dialog=byId('figure-name-dialog'), field=byId('figure-name-text'), suffix=byId('figure-name-subscript');
  const subs='₀₁₂₃₄₅₆₇₈₉';
  function names(model){
    const linkage=model.linkage||'', label=studio.state.library?.linkages?.[linkage]?.label?.split(' · ')[0]||linkage.replaceAll('_',' ');
    const counts=model.sites||[], feeds=model.feeds;
    let coefficients=[];
    if(feeds.length===2&&counts.length===2&&counts.every(n=>Number.isInteger(n)&&n>0)){
      const gcd=(a,b)=>b?gcd(b,a%b):a;
      const g=gcd(...counts);coefficients=[counts[1]/g,counts[0]/g];
    }
    const compact=feeds.map((f,i)=>{
      const name=f.label||('ABC'[i]||'Monomer'),n=coefficients[i];
      if(!n)return name;
      if(n===1)return i?'('+name+')':name;
      return (i?'('+name+')':name)+String(n).replace(/\d/g,d=>subs[d]);
    });
    const formula=compact.join(coefficients.length?'':' + ').slice(0,100);
    return {text:formula||model.common||'COF',sub:label.slice(0,45)};
  }
  function displayName(model){
    if(customName?.key===model.key)return customName;
    return model.common?{text:model.common,sub:''}:names(model);
  }
  function openName(formula=false){
    const m=activeModel();if(!m)return;
    const value=formula?names(m):displayName(m);
    field.setCustomValidity('');field.value=value.text;suffix.value=value.sub;
    byId('figure-name-remove').hidden=!nameEnabled;
    if(!dialog.open){dialog.returnValue='';dialog.showModal();}
    field.focus();field.select();
  }
  byId('figure-name').onclick=()=>openName();
  byId('figure-name-formula').onclick=()=>openName(true);
  byId('figure-name-remove').onclick=()=>{nameEnabled=false;customName=null;dialog.close('remove');refresh();};
  for(const edit of [field,suffix])edit.onkeydown=e=>{if(e.key==='Enter'){e.preventDefault();dialog.querySelector('form').requestSubmit(byId('figure-name-apply'));}};
  dialog.querySelector('form').onsubmit=e=>{
    if(e.submitter?.value==='cancel')return;
    e.preventDefault();const m=activeModel();if(!m)return;
    if(!field.value.trim()){field.setCustomValidity('Enter a name or remove the label.');field.reportValidity();return;}
    nameEnabled=true;customName={key:m.key,text:field.value.trim(),sub:suffix.value.trim()};
    dialog.close('apply');refresh();
  };
  field.oninput=()=>field.setCustomValidity('');
  const nameHit=document.createElement('button');
  nameHit.type='button';nameHit.className='figure-name-hit';nameHit.setAttribute('aria-label','Edit the name printed on the figure');
  nameHit.title='Double-click to edit the COF name';
  nameHit.ondblclick=()=>openName();
  nameHit.onclick=e=>{if(e.detail===0)openName();};
  function positionName(){
    if(!graphic?.nameBox||!figure.isConnected)return;
    const {width,height,nameBox:b}=graphic;
    const scale=Math.min(figure.clientWidth/width,figure.clientHeight/height);
    const x=figure.offsetLeft+(figure.clientWidth-width*scale)/2,y=figure.offsetTop+(figure.clientHeight-height*scale)/2;
    Object.assign(nameHit.style,{left:(x+b.x*scale)+'px',top:(y+b.y*scale)+'px',width:b.w*scale+'px',height:b.h*scale+'px'});
  }
  new ResizeObserver(positionName).observe(figure);
  function editing(){return !window.COFdrawCrystalDesign?.active&&byId('edit-geometry').getAttribute('aria-pressed')==='true';}
  function activeModel(){
    const crystal=window.COFdrawCrystalDesign;
    if(crystal?.active){
      const provenance=crystal.report?.structure?.provenance;
      if(!crystalView||!provenance)return null;
      const feeds=(provenance.precursors||[]).map((smiles,i)=>({smiles,label:crystal.precursor('ab'[i])?.label}));
      return {feeds,text:crystalView.text,url:crystalView.download,parent:byId('canvas-sheet'),key:JSON.stringify([provenance.kind,provenance.model_name,provenance.net,provenance.network_count,provenance.precursors,provenance.linkage]),linkage:provenance.linkage,sites:provenance.precursor_connectivities||[],common:studio.selectedFigureName()||(['reported_crystal','curated_crystal'].includes(provenance.kind)?provenance.model_name:null)};
    }
    const view=window.COFdrawExtent?.visible, s=studio.state;
    if(!view||!s.recipe||!view.downloads?.svg)return null;
    let feeds=s.recipe.assembly?.precursors;
    if(!feeds?.length){
      feeds=Object.values(s.recipe.components||{}).filter(c=>c.monomer_smiles)
        .map(c=>({label:c.label,smiles:c.monomer_smiles}));
      if(!feeds.length)feeds=(s.sourceMetadata?.monomers||[]).map(m=>({label:m.label,smiles:m.smiles||m.monomer_smiles}));
    }
    const components=Object.values(s.recipe.components||{});
    return {feeds:feeds||[],url:view.downloads.svg,parent:byId('canvas-sheet'),linkage:s.recipe.linkage,
      sites:!s.recipe.assembly&&components.length===2?components.map(c=>c.connectivity):[],
      common:studio.selectedFigureName(),key:JSON.stringify([s.sourceId,s.recipe.linkage,s.recipe.template_id,s.recipe.assembly?.monomer_ids,s.recipe.assembly?.linker_sequence,(feeds||[]).map(f=>[f.smiles,f.label])])};
  }
  function updateControls(){
    const model=activeModel(),canDraw=Boolean(model?.feeds.length&&model.feeds.length<=3&&model.feeds.every(f=>f.smiles));
    const nameButton=byId('figure-name');
    nameButton.disabled=!model;nameButton.textContent=nameEnabled&&model?'Edit COF name':'Add COF name';
    for(const button of toggles){
      button.disabled=!canDraw;
      button.setAttribute('aria-pressed',String(enabled&&canDraw));
      button.textContent=enabled&&canDraw?'Hide monomers + arrow':'Add monomers + arrow';
      button.title=canDraw?'Include the starting materials in the figure':'Requires a rendered framework with identified precursors';
    }
    for(const control of directions){control.hidden=!enabled||!canDraw;control.value=direction;}
  }
  function clear(){
    revision++;
    figure.parentElement?.classList.remove('showing-scheme');
    byId('canvas-viewport').style.removeProperty('--figure-aspect');
    nameHit.remove();figure.remove();figure.removeAttribute('src');
    if(graphic)URL.revokeObjectURL(graphic.url);
    graphic=null;
  }
  function element(name,attrs={}){
    const node=document.createElementNS(ns,name);
    for(const [key,value] of Object.entries(attrs))node.setAttribute(key,String(value));
    return node;
  }
  function embed(parent,text,x,y,width,height,role){
    const doc=new DOMParser().parseFromString(text,'image/svg+xml');
    if(doc.documentElement.localName!=='svg'||doc.querySelector('parsererror'))throw Error('The chemical drawing could not be read.');
    const svg=document.importNode(doc.documentElement,true);
    if(!svg.getAttribute('viewBox')){
      const w=parseFloat(svg.getAttribute('width')),h=parseFloat(svg.getAttribute('height'));
      if(!(w>0&&h>0))throw Error('The structure drawing has no dimensions.');
      svg.setAttribute('viewBox',[0,0,w,h].join(' '));
    }
    // Each embedded SVG has its own identifier namespace, including marker, clip and gradient references.
    const ids=new Map([...svg.querySelectorAll('[id]')].map(node=>[node.id,role+'-'+node.id]));
    for(const node of [svg,...svg.querySelectorAll('*')]){
      if(node.id&&ids.has(node.id))node.id=ids.get(node.id);
      for(const attr of [...node.attributes]){
        if(attr.value.includes('url(#'))attr.value=attr.value.replace(/url\(#([^)]*)\)/g,(all,id)=>ids.has(id)?'url(#'+ids.get(id)+')':all);
        if(attr.localName==='href'&&attr.value.startsWith('#')&&ids.has(attr.value.slice(1)))attr.value='#'+ids.get(attr.value.slice(1));
      }
    }
    for(const [key,value] of Object.entries({x,y,width,height,'data-part':role,preserveAspectRatio:'xMidYMid meet'}))svg.setAttribute(key,value);
    parent.append(svg);
  }
  function label(root,x,y,value,size,full){
    const t=element('text',{x,y,'text-anchor':'middle','font-family':'Arial, sans-serif','font-size':size,fill:'#253633'});
    t.textContent=value;
    if(full){const title=element('title');title.textContent=full;t.append(title);}
    root.append(t);
  }
  function compose(product,feeds,images,labeled){
    const vertical=direction==='vertical',n=feeds.length,unit=166,gap=22,feedWidth=n*unit+(n-1)*gap;
    const productBox=new DOMParser().parseFromString(product,'image/svg+xml').documentElement.getAttribute('viewBox')?.trim().split(/[ ,]+/).map(Number);
    const ratio=productBox?.[2]>0&&productBox?.[3]>0?productBox[2]/productBox[3]:1;
    const productHeight=560,productWidth=Math.round(productHeight*Math.min(2.5,Math.max(.55,ratio)));
    const width=vertical?Math.max(productWidth+28,feedWidth+28):feedWidth+productWidth+102;
    const height=vertical?productHeight+235:productHeight+20;
    const root=element('svg',{viewBox:'0 0 '+width+' '+height,width,height,role:'img','aria-label':'COF precursor reaction scheme'});
    const startX=vertical?(width-feedWidth)/2:10,top=vertical?8:(height-165)/2;
    feeds.forEach((f,i)=>{
      const x=startX+i*(unit+gap);
      embed(root,images[i],x,top,unit,143,'precursor-'+(i+1));
      const full=f.label||('Monomer '+'ABC'[i]),caption=full.length>23?full.slice(0,22)+'…':full;
      label(root,x+unit/2,top+161,caption,13,full);
      if(i<n-1)label(root,x+unit+gap/2,top+78,'+',22);
    });
    const ax=vertical?width/2:startX+feedWidth+8,ay=vertical?183:height/2;
    const path=vertical?'M '+ax+' '+ay+' v 37 m -5 -8 5 8 5 -8':'M '+ax+' '+ay+' h 59 m -9 -5 9 5 -9 5';
    root.append(element('path',{d:path,stroke:'#263633','stroke-width':2,fill:'none','stroke-linejoin':'round','data-part':'reaction-arrow'}));
    const px=vertical?(width-productWidth)/2:width-productWidth-10,py=vertical?224:10;
    embed(root,product,px,py,productWidth,productHeight,'framework');
    let nameBox=null;
    if(labeled){
      const k=Math.min(productWidth/labeled.width,productHeight/labeled.height),b=labeled.nameBox;
      nameBox={x:px+(productWidth-k*labeled.width)/2+k*b.x,y:py+(productHeight-k*labeled.height)/2+k*b.y,w:k*b.w,h:k*b.h};
    }
    return {text:new XMLSerializer().serializeToString(root),width,height,nameBox};
  }
  function namedProduct(product,model){
    const source=new DOMParser().parseFromString(product,'image/svg+xml').documentElement;
    let [, ,width,height]=(source.getAttribute('viewBox')||'').trim().split(/[ ,]+/).map(Number);
    if(!(width>0&&height>0)){width=parseFloat(source.getAttribute('width'));height=parseFloat(source.getAttribute('height'));}
    if(!(width>0&&height>0))throw Error('The drawing has no usable dimensions.');
    const result=displayName(model);
    const full=result.sub?'{['+result.text+']'+result.sub+'}':result.text;
    const font=Math.min(Math.min(width,height)*.052,width*1.72/Math.max(1,full.length));
    const space=Math.max(Math.min(width,height)*.067,font*1.9);
    const root=element('svg',{viewBox:'0 0 '+width+' '+(height+space),width,height:height+space});
    embed(root,product,0,0,width,height,'chemical-product');
    const caption=element('text',{'data-cof-name':'1',x:width/2,y:height+space*.59,'text-anchor':'middle','font-family':'Arial, sans-serif','font-size':font,'font-weight':600,fill:'#17231f'});
    caption.append(document.createTextNode(result.sub?'{['+result.text+']':result.text));
    if(result.sub){
      const sub=element('tspan',{dy:font*.26,'font-size':font*.65});sub.textContent=result.sub;caption.append(sub);
      const end=element('tspan',{dy:-font*.26,'font-size':font});end.textContent='}';caption.append(end);
    }
    root.append(caption);
    const hitWidth=Math.min(width,Math.max(font*3,font*.73*full.length));
    return {text:new XMLSerializer().serializeToString(root),width,height:height+space,nameBox:{x:(width-hitWidth)/2,y:height,w:hitWidth,h:space}};
  }
  async function refresh(){
    clear();updateControls();
    const model=activeModel();
    const withFeeds=Boolean(enabled&&model?.feeds.length&&model.feeds.length<=3&&model.feeds.every(f=>f.smiles));
    if(!model||editing()||(!withFeeds&&!nameEnabled))return;
    const token=revision;
    try{
      const getProduct=async()=>{
        if(model.text)return model.text;
        const response=await fetch(new URL('.'+model.url,document.baseURI));
        if(!response.ok)throw Error('The displayed framework export is unavailable.');
        return response.text();
      };
      const [product,...images]=await Promise.all([getProduct(),...(withFeeds?model.feeds.map(f=>studio.precursorSVG(f.smiles)):[])]);
      if(token!==revision)return;
      const measured=window.COFdrawCrystalDesign?.active?product:window.COFdrawMeasure?.decorate(product)||product;
      const labeled=nameEnabled?namedProduct(measured,model):null;
      const result=withFeeds?compose(labeled?.text||measured,model.feeds,images,labeled):labeled;
      figure.alt=[withFeeds?'Precursor structures leading to the framework':'Framework drawing',nameEnabled?displayName(model).text:''].filter(Boolean).join(' · ');

      const url=URL.createObjectURL(new Blob([result.text],{type:'image/svg+xml'}));
      graphic={...result,url};
      byId('canvas-viewport').style.setProperty('--figure-aspect',String(result.width/result.height));
      figure.src=url;model.parent.append(figure);model.parent.classList.add('showing-scheme');
      if(result.nameBox){model.parent.append(nameHit);positionName();}
    }catch(e){if(token===revision){enabled=false;nameEnabled=false;clear();updateControls();studio.notice(e.message,true);}}
  }
  function setEnabled(value){
    if(value&&!window.COFdrawCrystalDesign?.active){
      if(byId('edit-geometry').getAttribute('aria-pressed')==='true')byId('edit-geometry').click();
    }
    enabled=value;refresh();
  }
  toggles.forEach(b=>b.onclick=()=>setEnabled(!enabled));
  directions.forEach(s=>s.onchange=()=>{direction=s.value;refresh();});
  byId('edit-geometry').addEventListener('click',()=>{if(enabled||nameEnabled)refresh();});
  window.addEventListener('cofdraw:arrowchanged',()=>{if(enabled||nameEnabled)refresh();});
  window.addEventListener('cofdraw:crystalviewready',e=>{crystalView=e.detail;refresh();});
  window.addEventListener('cofdraw:crystalviewclear',()=>{crystalView=null;clear();updateControls();});
  window.addEventListener('cofdraw:viewready',()=>{if(!window.COFdrawCrystalDesign?.active)refresh();});
  window.addEventListener('cofdraw:extentchanged',()=>{clear();updateControls();});
  window.addEventListener('cofdraw:state',()=>{if(!activeModel()&&graphic)clear();updateControls();});
  async function save(ext){
    if(!graphic)throw Error('Wait for the reaction figure to finish drawing.');
    const current=graphic;
    if(ext==='svg')return studio.saveText('COFdraw-reaction.svg',current.text,'image/svg+xml');
    const img=new Image();img.src=current.url;await img.decode();
    const canvas=document.createElement('canvas'),size=Number(byId('image-size').value)||1600;
    const factor=size/Math.max(current.width,current.height);
    canvas.width=Math.max(1,Math.round(current.width*factor));canvas.height=Math.max(1,Math.round(current.height*factor));
    canvas.getContext('2d').drawImage(img,0,0,canvas.width,canvas.height);
    const bridge=window.pywebview?.api;
    if(bridge?.save_image)return bridge.save_image('COFdraw-reaction.png',canvas.toDataURL('image/png'));
    const blob=await new Promise(resolve=>canvas.toBlob(resolve)),a=document.createElement('a');
    a.href=URL.createObjectURL(blob);a.download='COFdraw-reaction.png';a.click();setTimeout(()=>URL.revokeObjectURL(a.href),3000);
  }
  byId('exports').addEventListener('click',e=>{
    const ext=e.target.closest('button')?.textContent.toLowerCase();
    if(!['svg','png'].includes(ext)||!activeModel()||editing())return;
    if(enabled||nameEnabled){e.stopImmediatePropagation();e.preventDefault();save(ext).catch(err=>studio.notice(err.message,true));}
    else if(window.COFdrawMeasure?.active){e.stopImmediatePropagation();window.COFdrawMeasure.save(ext);}
  },true);
  window.COFdrawScheme={get composed(){return Boolean(graphic);},save,snapshot:()=>({enabled,direction,nameEnabled,customName}),restore(v){
    enabled=v?.enabled===true;direction=v?.direction==='vertical'?'vertical':'horizontal';nameEnabled=v?.nameEnabled===true;
    const n=v?.customName;
    customName=n&&typeof n.key==='string'&&typeof n.text==='string'&&typeof n.sub==='string'?{key:n.key,text:n.text.slice(0,100),sub:n.sub.slice(0,45)}:null;refresh();
  }};
})();
