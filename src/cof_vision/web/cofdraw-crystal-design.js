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
  let poreText='', poreRequest=0;
  const sharedStyles=new Set(['standard','mono_dashed','mono_wave','serif','red_blue','hetero_colors',
    'element_compact','node_highlight','ring_fill','pastel_blocks','bold_color','glow_outline']);
  const depthStyles=[
    ['reticular_depth','Depth fade · ink and ring wash','Fade rear bonds and shade aromatic rings.'],
    ['reticular_blocks','Depth fade · building blocks','Distinguish assigned precursors with the selected palette.'],
    ['depth_ink','Front and back ink','Color near and rear bonds by projected depth.']
  ].map(([id,label,description])=>({id,label,description,category:'layered',preview_kind:id}));
  let restoredReference=false, comparing=false;
  function availableStyles(){
    return [...(studio.state.library?.styles||[]).filter(s=>sharedStyles.has(s.id)),...depthStyles];
  }
  function clearImage(){
    if(visible?.url)URL.revokeObjectURL(visible.url);
    visible=null;picture.hidden=true;picture.removeAttribute('src');$('exports').replaceChildren();
    window.dispatchEvent(new Event('cofdraw:crystalviewclear'));
  }
  function clearPicture(){
    poreText='';poreRequest++;$('crystal-layout').replaceChildren();$('crystal-projection').hidden=true;
    job++;report=null;preparedCells=null;cache.clear();comparing=false;clearImage();
    $('crystal-drawing-view').replaceChildren();$('representation-note').hidden=true;
  }
  function leave(){
    if(!id)return;
    clearTimeout(timer);revision++;id=null;edit=null;failure='';clearPicture();main.classList.remove('crystal-design-active');
    $('crystal-design-controls').hidden=true;$('crystal-view-field').hidden=true;
    if(patchStyle){
      const {style,styleOptions}=studio.state;
      Object.assign(studio.state,patchStyle);
      if(sharedStyles.has(style))Object.assign(studio.state,{style,styleCategory:'all',styleOptions,resolvedStyleOptions:{}});
      else for(const k of ['palette','custom_palette','element_scheme'])if(k in styleOptions)studio.state.styleOptions[k]=styleOptions[k];
      patchStyle=null;
    }
    $('canvas-title').textContent='Framework canvas';$('canvas-sheet').classList.remove('transparent-preview');
    $('extent-bar').hidden=true;$('canvas-empty').classList.remove('hidden');studio.refreshStyles();
    const customEntry=$('cof-preset').querySelector('option[value="crystal:custom-dia"]');
    if(customEntry)customEntry.value='custom';
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
  function referenceReaction(){
    const starter=entries.get(id)?.starter;
    if(!restoredReference||edit||report?.structure?.provenance?.template!==id||!starter?.monomer_ids?.length)return null;
    const catalog=studio.state.library?.monomer_presets||[];
    const feeds=starter.monomer_ids.map(key=>catalog.find(m=>m.id===key));
    if(!feeds.every(m=>m?.smiles))return null;
    return {linkage:starter.linkage,feeds:feeds.map(m=>({smiles:m.smiles,label:m.abbreviation||m.label}))};
  }
  function choice(item,partner){
    const starter=restoredReference&&entries.get(id)?.starter;
    if(starter){
      const roleLabel=[starter.net,starter.linkage?.replaceAll('_',' ')].filter(Boolean).join(' · ');
      return {compatible:true,roleLabel,reason:'Reported precursor for '+(starter.common_name||'this crystal model')+'.'};
    }
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
    const {style,styleOptions,library}=studio.state;
    const spec=library?.style_controls||{};
    const opt={...(spec.defaults?.[style]||{}),...styleOptions};
    const colors=opt.custom_palette||spec.palette_colors?.[opt.palette||'classic']||['#d82c2c','#264bd1'];
    const depthFade=style.startsWith('reticular_')||style==='depth_ink', wash=depthFade||style==='ring_fill';
    const junction=$('crystal-drawing-view').value==='Framework drawing'&&count===1&&Boolean(report?.chemical_junction_svg);
    const blockInk=['reticular_blocks','red_blue','bold_color'].includes(style);
    const scheme=opt.element_scheme&&opt.element_scheme!=='style_default'?opt.element_scheme:
      style==='element_compact'?'jewel':style==='hetero_colors'?'conventional':'';
    const elementInk=spec.element_colors?.[scheme]||{};
    const assigned=Boolean(svg.querySelector('[data-block="A"],[data-block="B"]'));
    const monochrome=['mono_wave','mono_dashed','serif'].includes(style);
    const highlight=['node_highlight','pastel_blocks','glow_outline'].includes(style);
    const halos=document.createElementNS(svg.namespaceURI,'g');
    halos.setAttribute('data-highlights',style);
    halos.setAttribute('stroke-linecap','round');
    halos.setAttribute('stroke-linejoin','round');
    halos.setAttribute('fill','none');
    if(style==='glow_outline'&&opt.glow_strength>0)halos.style.filter='blur('+(1+2*opt.glow_strength)+'px)';
    if(wash)svg.querySelectorAll('[data-ring-washes]').forEach(g=>g.removeAttribute('display'));
    const parentWidth=Number(svg.querySelector('g[stroke-width]')?.getAttribute('stroke-width'))||1.65;
    for(const mark of svg.querySelectorAll('[data-depth]')){
      const block=mark.dataset.block, ink=style==='depth_ink'?colors[Number(mark.dataset.depth)<.48?1:0]:block==='A'?colors[0]:block==='B'?colors[1]||colors[0]:assigned?null:colors[0];
      const floor=junction?.60:.19;
      const visible=depthFade?floor+(1-floor)/(1+Math.exp(-9*(Number(mark.dataset.depth)-.48))):1;
      if(mark.localName==='polygon'){
        mark.setAttribute('fill',ink||'#99a3ad');
        mark.setAttribute('opacity',((opt.fill_opacity??.18)*visible).toFixed(3));
      }else if(mark.localName==='g'){
        if(mark.dataset.cut){
          const cut=opt.edge_style||'wavy_cross',wave=mark.querySelector('polyline');
          if(wave&&cut!=='wavy_cross'){
            wave.remove();
            if(cut==='crossbar'){
              const [x,y]=mark.dataset.cut.split(',').map(Number),[nx,ny]=mark.dataset.cutNormal.split(',').map(Number);
              const bar=document.createElementNS(svg.namespaceURI,'path');
              bar.setAttribute('d','M'+(x-4*nx)+','+(y-4*ny)+'l'+(8*nx)+','+(8*ny));
              mark.append(bar);
            }else if(cut==='dashed')mark.setAttribute('stroke-dasharray','3 2');
          }
        }
        if(monochrome)mark.setAttribute('stroke','#25282b');
        mark.setAttribute('opacity',visible.toFixed(3));
        if((blockInk||style==='depth_ink')&&ink){mark.setAttribute('stroke',ink);mark.setAttribute('color',ink);}
        mark.setAttribute('stroke-width',parentWidth*(opt.bond_scale||1)*(style==='bold_color'?1.65:1));
      }else if(mark.localName==='path'){
        mark.setAttribute('stroke-opacity',visible.toFixed(3));
        mark.setAttribute('stroke-width',Number(mark.getAttribute('stroke-width'))*(opt.bond_scale||1)*(style==='bold_color'?1.65:1));
        if((blockInk||style==='depth_ink')&&ink)mark.setAttribute('stroke',ink);
        else if(elementInk[mark.dataset.element])mark.setAttribute('stroke',elementInk[mark.dataset.element]);
      }else if(mark.localName==='text'){
        mark.setAttribute('opacity',visible.toFixed(3));
        if(opt.font_family)mark.setAttribute('font-family',opt.font_family);
        if(opt.font_scale)mark.setAttribute('font-size',opt.font_scale*Number(mark.getAttribute('font-size')||mark.parentNode.getAttribute('font-size')||14));
        if((opt.hydrogen_labels||(style==='element_compact'?'compact':'show'))==='compact')mark.textContent=mark.textContent.replace(/^(N|O|S)H$/,'$1');
        const symbol=mark.dataset.element||mark.textContent.match(/[A-Z][a-z]?/)?.[0];
        const color=elementInk[symbol]||(monochrome?'#25282b':null);
        if(style==='bold_color')mark.setAttribute('font-weight','bold');
        if(color||blockInk&&ink)mark.setAttribute('fill',color||ink);
      }
      if(highlight&&['g','path'].includes(mark.localName)){
        const halo=mark.cloneNode(true);
        halo.removeAttribute('data-depth');
        const width=parentWidth*(opt.highlight_width||1)*(style==='pastel_blocks'?9:6);
        halo.setAttribute('stroke',ink||'#98a4a7');
        halo.setAttribute('color',ink||'#98a4a7');
        halo.setAttribute('stroke-width',width);
        halo.setAttribute('opacity',opt.fill_opacity??.25);
        halo.removeAttribute('stroke-opacity');
        halos.append(halo);
      }
    }
    if(halos.childElementCount){const paper=svg.querySelector(':scope > rect[data-paper]');if(paper)paper.after(halos);else svg.prepend(halos);}
    if(count>1&&$('crystal-drawing-view').value==='Framework drawing'){
      // At a crossing, let the front phenyl silhouette and its bonds occlude
      // the rear chair. A luminance mask also works on transparent exports.
      const ns=svg.namespaceURI,box=svg.viewBox.baseVal;
      const defs=document.createElementNS(ns,'defs'),mask=document.createElementNS(ns,'mask');
      mask.id='cofdraw-chair-occlusion';mask.setAttribute('maskUnits','userSpaceOnUse');
      for(const [key,value] of Object.entries({x:box.x,y:box.y,width:box.width,height:box.height}))mask.setAttribute(key,value);
      const field=document.createElementNS(ns,'rect');
      for(const key of ['x','y','width','height'])field.setAttribute(key,mask.getAttribute(key));
      field.setAttribute('fill','white');mask.append(field);
      for(const mark of svg.querySelectorAll('[data-depth]')){
        const depth=Number(mark.dataset.depth);
        if(depth<.33){
          mark.setAttribute('mask','url(#cofdraw-chair-occlusion)');
        }else if(depth>.67&&['polygon','g'].includes(mark.localName)){
          const silhouette=mark.cloneNode(true);
          silhouette.removeAttribute('data-depth');
          silhouette.setAttribute('opacity','1');silhouette.setAttribute('color','black');
          silhouette.setAttribute('fill',mark.localName==='polygon'?'black':'none');
          silhouette.setAttribute('stroke','black');silhouette.setAttribute('stroke-width',mark.localName==='polygon'?'4':'5');
          silhouette.setAttribute('stroke-linejoin','round');silhouette.setAttribute('stroke-linecap','round');
          mask.append(silhouette);
        }
      }
      defs.append(mask);svg.prepend(defs);
    }
    svg.setAttribute('data-circuit-count',count);svg.setAttribute('data-drawing-style',style);
    if(opt.transparent_background){
      svg.querySelectorAll('rect[data-paper]').forEach(r=>r.remove());
    }else if(!svg.querySelector('rect[data-paper]')){
      const bg=document.createElementNS(svg.namespaceURI,'rect');bg.setAttribute('width','100%');bg.setAttribute('height','100%');bg.setAttribute('fill','white');svg.prepend(bg);
    }
    if(poreText){
      const box=svg.viewBox.baseVal,h=box.height,w=box.width,strip=Math.max(30,w*.057);
      const label=document.createElementNS(svg.namespaceURI,'text');
      label.setAttribute('x',w/2);label.setAttribute('y',h+strip*.68);label.setAttribute('text-anchor','middle');
      label.setAttribute('font-size',strip*.39);label.setAttribute('font-family','Arial,sans-serif');label.setAttribute('fill','#205558');label.textContent=poreText;
      svg.append(label);svg.setAttribute('viewBox','0 0 '+w+' '+(h+strip));svg.setAttribute('height',h+strip);
    }
    const box=svg.viewBox.baseVal;
    if(box.width&&box.height){svg.setAttribute('width',box.width);svg.setAttribute('height',box.height);}
    return new XMLSerializer().serializeToString(svg);
  }
  function drawGeometry(source){
    const svg=$('crystal-layout');svg.replaceChildren();
    const junction=usesJunction(), knobs=$('crystal-projection');
    if(!['Framework drawing','Projected circuit'].includes($('crystal-drawing-view').value)){knobs.hidden=true;return;}
    knobs.hidden=!junction&&!source.junction_loops?.length;
    const adjust=$('crystal-angle');
    $('crystal-angle-label').textContent=junction?'Tetrahedral arm spread':'Circuit tilt';
    adjust.min=junction?-24:-45;adjust.max=junction?24:45;
    adjust.value=junction?armSpread:tilt;
    $('crystal-angle-value').textContent=adjust.value+'°';
    if(junction){
      for(let k=0;k<4;k++){
        const angle=([30,150,230,310][k]+(-1)**k*.3*armSpread)*Math.PI/180;
        const x=160+68*Math.cos(angle),y=90-68*Math.sin(angle);
        add('line',{x1:160,y1:90,x2:x,y2:y,stroke:'#407077','stroke-width':k===2?5:3,opacity:k===3?.55:1,'stroke-dasharray':k===3?'3 3':'none'});
        add('circle',{cx:x,cy:y,r:4,fill:'#295861'});
      }
      add('circle',{cx:160,cy:90,r:6,fill:'#295861'});
      return;
    }
    if(!source.junction_loops?.length)return;
    const n=source.direction,axis=[0,0,0];
    axis[n.map(Math.abs).indexOf(Math.min(...n.map(Math.abs)))]=1;
    const cross=(a,b)=>[a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]];
    const dot=(a,b)=>a.reduce((s,x,i)=>s+x*b[i],0);
    const u=cross(axis,n),v=cross(n,u),cell=report.structure.cell;
    const loops=source.junction_loops.map(loop=>loop.map(f=>{
      const xyz=n.map((_,i)=>f.reduce((s,t,j)=>s+t*cell[j][i],0));
      return [dot(xyz,u),-dot(xyz,v),dot(xyz,n)];
    }));
    const points=loops.flat(),bounds=[0,1,2].map(i=>[Math.min(...points.map(p=>p[i])),Math.max(...points.map(p=>p[i]))]);
    const scale=Math.min(270/(bounds[0][1]-bounds[0][0]||1),150/(bounds[1][1]-bounds[1][0]||1));
    const pos=p=>[160+scale*(p[0]-(bounds[0][0]+bounds[0][1])/2),90+scale*(p[1]-(bounds[1][0]+bounds[1][1])/2)];
    function add(tag,attrs){const el=document.createElementNS(svg.namespaceURI,tag);for(const [k,v] of Object.entries(attrs))el.setAttribute(k,v);svg.append(el);}
    const edges=loops.flatMap(loop=>loop.map((p,i)=>[p,loop[(i+1)%loop.length]])).sort((a,b)=>a[0][2]+a[1][2]-b[0][2]-b[1][2]);
    for(const [p,q] of edges){const [x1,y1]=pos(p),[x2,y2]=pos(q);add('line',{x1,y1,x2,y2,stroke:'#407077','stroke-width':2,opacity:.25+.75*((p[2]+q[2])/2-bounds[2][0])/(bounds[2][1]-bounds[2][0]||1)});}
    for(const p of points){const [cx,cy]=pos(p);add('circle',{cx,cy,r:3.5,fill:'#295861'});}
  }
  let armSpread=0,tilt=0;
  function usesJunction(){return count===1&&$('crystal-drawing-view').value==='Framework drawing'&&Boolean(report?.chemical_junction_svg);}
  $('crystal-angle').oninput=e=>{$('crystal-angle-value').textContent=e.target.value+'°';};
  $('crystal-angle').onchange=e=>{if(usesJunction())armSpread=Number(e.target.value);else tilt=Number(e.target.value);drawView();};
  $('crystal-angle-reset').onclick=()=>{if(usesJunction())armSpread=0;else tilt=0;drawView();};
  function circuitOptions(){
    const type=$('crystal-drawing-view').value;
    const traced=['Framework drawing','Projected circuit'].includes(type),counts=traced?(report.circuit_extents||[1]):[1];
    window.COFdrawExtent?.setCrystalOptions(counts.map(n=>({value:String(n),label:n===1&&type==='Framework drawing'&&report.chemical_junction_svg?'Tetrahedral junction':traced?n+(n===1?' circuit':' connected circuits'):type==='Linker detail'?'One fragment':'Displayed cells'})),String(count),n=>drawView(Number(n)));
  }
  async function drawView(n=count){
    if(!report)return;
    const token=++job, model=report, type=$('crystal-drawing-view').value;
    count=['Framework drawing','Projected circuit'].includes(type)&&model.circuit_extents?.includes(n)?n:1;
    clearImage();$('canvas-empty').classList.add('hidden');circuitOptions();
    const key=[count,usesJunction()?armSpread:0,usesJunction()?0:tilt].join('/');
    let source=cache.get(key),text='',caption='';
    try{
      if(!source&&!model.circuit_extents?.length&&!model.chemical_junction_svg)source={};
      if(!source){
        $('drawing-summary').textContent='Tracing neighboring circuits…';
        source=studio.api('science/reports/'+model.report_id+'/circuits/'+count,{arm_spread:usesJunction()?armSpread:0,tilt:usesJunction()?0:tilt});
        if(cache.size>=12)cache.delete(cache.keys().next().value);
        cache.set(key,source);
      }
      source=await source;
      if(token!==job||model!==report)return;
      drawGeometry(source);
      if(type==='Framework drawing'){
        text=source.junction_svg||source.scheme_svg||source.projection_svg;
        caption=source.junction_svg?'Tetrahedral junction · four connected arms':source.scheme_svg?'Bond-line geometry · crystal depth':'Cartesian projection';
      }else if(type==='Projected circuit'){text=source.projection_svg;caption='Cartesian projection · actual atomic coordinates';}
      else if(type==='Chemical drawing + periodic network'){text=window.COFdrawChemicalNetworkFigure(model,preparedCells);caption='Bond-line drawing and the corresponding repeated crystal';}
      else if(type==='Periodic packing'){text=window.COFdrawChemicalNetworkFigure(model,preparedCells,true);caption='Periodic bonds in the selected crystal';}
      else{text=model.chemical_drawing_svg;caption='Local linker connectivity';}
      if(!text)throw Error('No chemical drawing is available for this view.');
      text=styleSVG(text);
      const url=URL.createObjectURL(new Blob([text],{type:'image/svg+xml'}));
      visible={text,url,label:type,caption,count};picture.src=url;picture.hidden=false;
      if(Boolean(document.getElementById('style-edge_style'))!==text.includes('data-cut='))studio.drawStyleControls();
      if(['red_blue','reticular_blocks','bold_color','ring_fill','reticular_depth','node_highlight','pastel_blocks','glow_outline'].includes(studio.state.style)){
        const assigned=['A','B'].filter(block=>text.includes('data-block="'+block+'"'));
        const unresolved=Number(text.match(/data-unassigned-sites="(\d+)"/)?.[1]||0);
        if(!assigned.length)caption+=' · Single framework color (precursor assignments unavailable)';
        else if(unresolved)caption+=' · '+unresolved+' heavy-atom sites have unresolved precursor identity';
      }
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
      cache.delete(key);$('drawing-summary').textContent=error.message;studio.sync();
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
    const options=portableInkOptions(style);
    if(style==='random'){
      const family=id=>sharedStyles.has(id)&&['standard','mono_dashed','mono_wave','serif'].includes(id)?'ink':
        ['node_highlight','pastel_blocks'].includes(id)?'highlight':['reticular_depth','reticular_blocks','depth_ink'].includes(id)?'depth':id;
      const options=availableStyles().filter(s=>family(s.id)!==family(state.style));style=options[Math.floor(Math.random()*options.length)].id;
    }
    if(!availableStyles().some(s=>s.id===style))return;
    state.style=style;state.styleOptions=options;state.resolvedStyleOptions={};studio.refreshStyles();drawView();studio.sync();
  }
  function enterDrawing(){
    if(!$('science-dialog').hidden)$('close-science').click();
    if(!$('space-workspace').hidden)window.dispatchEvent(new Event('cofdraw:leave-explorer'));
    if(!$('science-dialog').hidden||main.hidden)$('close-science').click();
    $('mode-manual').click();
    if(!patchStyle){
      const state=studio.state;
      patchStyle={style:state.style,styleOptions:state.styleOptions,resolvedStyleOptions:state.resolvedStyleOptions,styleCategory:state.styleCategory};
      state.style=sharedStyles.has(state.style)&&state.style!=='standard'?state.style:'reticular_depth';
      state.styleCategory='all';state.styleOptions={...state.styleOptions};state.resolvedStyleOptions={};count=1;
    }
    main.classList.add('crystal-design-active');$('crystal-design-controls').hidden=false;$('crystal-view-field').hidden=false;
    studio.refreshStyles();
  }
  function startCustom(){
    id=custom;edit=null;failure='';restoredReference=false;clearTimeout(timer);revision++;clearPicture();
    for(const control of ['dia-linkage','dia-net-count','reticular-net'])$(control).closest('label').hidden=false;
    const customEntry=$('cof-preset').querySelector('option[value="custom"]');
    if(customEntry)customEntry.value='crystal:'+custom;
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
    const match=studio.periodicPair();
    if(!match&&(studio.planarPair()||restoredReference)){
      leave();studio.queueDrawing();
      studio.notice('The reference model has been closed. The new precursors are editable in the drawing workspace.');
      return;
    }
    if(match&&studio.state.autoChemistry&&$('dia-linkage').value!==match&&$('dia-linkage').value!=='secondary_amine')$('dia-linkage').value=match;
    startCustom();
    if(!precursor('a').smiles||!precursor('b').smiles){
      $('drawing-summary').textContent='Choose or draw both building blocks.';return;
    }
    timer=setTimeout(build,600);
  }
  async function build(){
    if(!id)return;
    if(restoredReference)return drawView();
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
  async function select(key,modification=null,restoreReference=false,keepView=false){
    if(key===custom){id=custom;enterDrawing();startCustom();studio.setCrystalMonomers();await build();return;}
    const row=entries.get(key);if(!row)return;
    if(!keepView){count=1;armSpread=0;tilt=0;preferredView='Framework drawing';}
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
      if(result.chemical_junction_svg||result.chemical_circuit_svg)choices.push('Framework drawing');
      if(result.chemical_circuit_svg)choices.push('Projected circuit');
      if(cells&&result.chemical_framework_svg)choices.push('Chemical drawing + periodic network');
      if(cells?.viewer?.bonds?.length)choices.push('Periodic packing');
      if(result.chemical_drawing_svg)choices.push('Linker detail');
      report=result;preparedCells=cells;
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
  picker.addEventListener('change',()=>{if(picker.value==='aryl_f'){picker.value=id;selectInViewer(id,'aryl_f');}else select(picker.value,null,false,true);});
  $('compare-crystal-reference').onclick=()=>selectInViewer(id);
  $('crystal-drawing-view').addEventListener('change',e=>{preferredView=e.target.value;drawView();});
  $('measure-pore').addEventListener('click',async e=>{
    if(!id)return;
    e.preventDefault();e.stopImmediatePropagation();
    if(!report)return;
    const model=report,token=++poreRequest;
    if(poreText){poreText='';await drawView();return;}
    $('drawing-summary').textContent='Measuring the periodic cavity…';
    try{
      const r=await studio.api('science/geometry',{structure:model.structure,grid_spacing:.8,surface_samples:32,refinement:false,include_map:false,probe_radius:1.82});
      if(token!==poreRequest||model!==report)return;
      const range=r.lcd_bounds_angstrom;
      if(!range?.every(Number.isFinite))throw Error('Cavity size was unresolved for this model.');
      poreText='Periodic cavity (LCD)  '+range.map(v=>(v/10).toFixed(2)).join('–')+' nm';
      await drawView();
    }catch(err){if(token===poreRequest)$('drawing-summary').textContent=err.message;}
  },true);
  function returnToPatch(){window.dispatchEvent(new Event('cofdraw:prefer-standard'));leave();studio.queueDrawing();}
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
      if(!['dia','lon','4h','6h'].includes(topology))throw Error('Unknown reticular net.');
      id=custom;enterDrawing();$('reticular-net').value=topology;$('dia-linkage').value=linkage;$('dia-net-count').value=String(count);startCustom();studio.setCrystalMonomers(project.monomers);await build();
    }else if(!entries.has(project.template_id)||![null,'aryl_f'].includes(project.linker_edit??null))
      throw new Error('Unknown crystal model or linker modification in this project.');
    else await select(project.template_id,project.linker_edit??null,project.version===1||project.reference_only===true);
    if(report){
      armSpread=Number.isInteger(project.arm_spread)&&Math.abs(project.arm_spread)<=24?project.arm_spread:0;
      tilt=Number.isInteger(project.tilt)&&Math.abs(project.tilt)<=45?project.tilt:0;
      const legacy={'Depth fade · ink and ring wash':'reticular_depth','Depth fade · two building blocks':'reticular_blocks'};
      studio.state.style=availableStyles().some(s=>s.id===project.style)?project.style:legacy[project.view]||'reticular_depth';
      studio.state.styleOptions=project.style_options&&typeof project.style_options==='object'?project.style_options:{};
      const view=legacy[project.view]||project.view==='Three connected circuits'?'Framework drawing':project.view;
      const selector=$('crystal-drawing-view');
      if([...selector.options].some(o=>o.value===view))selector.value=preferredView=view;
      studio.refreshStyles();await drawView(project.circuits||(project.view==='Three connected circuits'?3:1));
    }
  }
  function save(){
    if(report)studio.saveText('COFdraw-crystal-project.json',JSON.stringify({format:'cofdraw-crystal-project',arm_spread:armSpread,tilt,reaction_scheme:window.COFdrawScheme?.snapshot(),
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
  window.COFdrawCrystalDesign={restore,save,get styles(){return availableStyles();},chooseStyle,redraw:()=>drawView(),get visible(){return visible;},get scope(){return scope;},get isReference(){return restoredReference;},referenceReaction,precursor,choice,build,inputsChanged,selectInViewer,has:key=>entries.has(key),get active(){return id;},get report(){return report;},leave,showPeriodic};
})();
