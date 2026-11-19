/* Optional monomer workspace, using the same editor, history and RDKit checks. */
"use strict";
(() => {
  const dialog=$("sketch-dialog"), canvas=$("monomer-canvas");
  const shelf=$("classic-shelf"), list=$("classic-list"), anchor=$("classic-anchor");
  const preview=$("classic-preview"), category=$("classic-category"), search=$("classic-search");
  const groups=[
    ["Amino","N"],["Hydroxy","O"],["Formyl","C=O"],["Cyano","C#N"],["Carboxyl","C(=O)O"],
    ["Acetyl","C(=O)C"],["Methoxy","OC"],["Ethynyl","C#C"],["Vinyl","C=C"],["Methyl","C"],
    ["Trifluoromethyl","C(F)(F)F"],["Thiol","S"],["Boronic acid","B(O)O"],
    ["Hydrazino","NN"],["Nitro","[NH+](=O)[O-]"],["Fluoro","F"],["Bromo","Br"]
  ];
  const rings=[
    ["Phenyl","c1ccccc1"],["Cyclopentyl","C1CCCC1"],["Cyclohexyl","C1CCCCC1"],
    ["Pyridyl","c1ccncc1"],["Pyrimidinyl","c1cncnc1"],["Triazinyl","c1ncncn1"],
    ["Thienyl","c1cccs1"],["Furyl","c1ccco1"],["Biphenyl","c1ccc(-c2ccccc2)cc1"],
    ["Naphthyl","c1ccc2ccccc2c1"],["Anthracenyl","c1ccc2cc3ccccc3cc2c1"],
    ["Pyrenyl","c1cc2ccc3cccc4ccc(c1)c2c34"]
  ];
  const motifs=rows=>rows.map(([name,smiles])=>({name,smiles}));
  const collections=new Map([["Functional groups",motifs(groups)],["Ring systems",motifs(rings)],["Imported",[]]]);
  collections.set("Bicyclics",motifs([["Norbornane","C1CC2CCC1C2"],["Bicyclo[2.2.2]octane","C1CC2CCC1CC2"],["Decalin","C1CCC2CCCCC2C1"],["Indane","c1ccc2c(c1)CCC2"]]));
  collections.set("Linkers and spacers",motifs([
    ["Ethylene","CC"],["Propylene","CCC"],["Ethylene glycol","OCCO"],["Diethylene glycol","OCCOCCO"],
    ["Diphenyl ether","c1ccc(Oc2ccccc2)cc1"],["Diphenylacetylene","c1ccc(C#Cc2ccccc2)cc1"],
    ["Benzophenone","c1ccc(C(=O)c2ccccc2)cc1"],["Triphenylamine","c1ccc(N(c2ccccc2)c2ccccc2)cc1"]
  ]));
  collections.set("Protecting groups",motifs([
    ["Benzyl","Cc1ccccc1"],["Acetyl","C(=O)C"],["Benzoyl","C(=O)c1ccccc1"],["tert-Butyl","C(C)(C)C"],
    ["Boc","C(=O)OC(C)(C)C"],["Cbz","C(=O)OCc1ccccc1"],["Trimethylsilyl","[SiH](C)(C)C"]
  ]));
  let rows=[], fragment=null, armed=false, revision=0;
  let discovery=null;
  let pageStart=0;
  const previews=new Map(), pageSize=12;
  function templateData(row){
    const key=row.molblock||row.smiles;
    if(!previews.has(key))previews.set(key,api("monomer",row.molblock?{molblock:row.molblock}:{smiles:row.smiles}).catch(e=>{previews.delete(key);throw e;}));
    return previews.get(key);
  }
  function openShelf(open){
    shelf.hidden=!open;$("classic-templates").setAttribute("aria-expanded",String(open));
    if(open){arm(false);drawTiles();search.focus();}
  }
  $("classic-templates").onclick=()=>openShelf(shelf.hidden);
  function drawTiles(){
    const grid=$("classic-tiles");grid.replaceChildren();
    if(shelf.hidden)return;
    const stop=Math.min(rows.length,pageStart+pageSize);
    $("tiles-page").textContent=rows.length?(pageStart+1)+"–"+stop+" of "+rows.length:"No matching templates";
    $("tiles-prev").disabled=!pageStart;$("tiles-next").disabled=stop===rows.length;
    rows.slice(pageStart,stop).forEach((row,index)=>{
      const tile=document.createElement("button"),img=document.createElement("img"),name=document.createElement("span");
      tile.type="button";tile.title=row.name;tile.disabled=true;img.alt="";name.textContent=row.name;tile.append(img,name);grid.append(tile);
      templateData(row).then(r=>{
        if(!tile.isConnected)return;
        const thumbnail=r.svg.replace(/stroke-width:([\d.]+)px/g,(_,v)=>'stroke-width:'+Number(v)*2.5+'px');
        img.src="data:image/svg+xml;charset=utf-8,"+encodeURIComponent(thumbnail);tile.disabled=false;
      }).catch(()=>{tile.title=row.name+" · preview unavailable";});
      tile.onclick=async()=>{
        list.value=String(pageStart+index);await showTemplate();
        if(fragment&&rows[Number(list.value)]===row){openShelf(false);arm(true);canvas.focus();sketchMessage("Click to place. Point at an atom to attach. Esc cancels.");}
      };
    });
  }
  $("tiles-prev").onclick=()=>{pageStart=Math.max(0,pageStart-pageSize);drawTiles();};
  $("tiles-next").onclick=()=>{if(pageStart+pageSize<rows.length)pageStart+=pageSize;drawTiles();};
  const shelfRows=result=>result.templates.map(t=>({...t,name:t.name+" · "+t.formula}));
  function discoverInstalled(){
    if(discovery)return;
    discovery=api("templates/installed").then(result=>{
      if(!result.collections.length)return;
      const section=document.createElement("optgroup");section.label="Installed ChemDraw templates";
      for(const name of result.collections)section.append(new Option("ChemDraw · "+name,"installed:"+name));
      category.append(section);
    }).catch(()=>{discovery=null;});
  }
  const view=new MonomerSketcher(preview);
  preview.setAttribute("aria-disabled","true");
  function arm(value){
    armed=value;
    $("classic-place").setAttribute("aria-pressed",String(value));
    canvas.classList.toggle("placing",value);
    canvas.querySelector(".placement-preview")?.remove();
  }
  function syncControls(){
    const usable=!!fragment && !sketchBusy;
    anchor.disabled=!usable;
    $("classic-place").disabled=!usable;
    for(const el of [category,list,search])el.disabled=sketchBusy;
  }
  dialog.addEventListener("cofdraw:sketch-busy",syncControls);
  function setMode(enabled){
    arm(false);openShelf(false);sketcher.cancelDrag();
    sketcher.classicMode=enabled;
    if(!enabled && sketcher.selected?.kind==="group"){sketcher.selected=null;sketcher.draw();}
    if(!enabled&&["pan","chain"].includes(sketcher.tool))sketcher.setTool("bond","1");
    dialog.classList.toggle("classic",enabled);
    $("classic-switch").textContent=enabled?"Standard layout":"Try classic layout";
    $("classic-switch").setAttribute("aria-pressed",String(enabled));
    if(enabled){discoverInstalled();populate();}
    else {revision++;fragment=null;}
  }
  const glyphs={
    'select':'M4 2L17 13L11 14L9 20Z',
    'bond:1':'M3 18L19 3',
    'bond:2':'M2 16L17 2 M5 20L20 6',
    'bond:3':'M1 14L14 1 M4 17L17 4 M7 20L20 7',
    'erase':'M2 14L13 3L20 10L10 20H8Z M6 11L13 17',
    'chain':'M1 15L6 7L11 15L16 7L21 15',
    'pan':'M11 1V21 M1 11H21 M8 4L11 1L14 4 M8 18L11 21L14 18 M4 8L1 11L4 14 M18 8L21 11L18 14'
  };
  for(const button of dialog.querySelectorAll('.sketch-tools [data-tool]')){
    const path=glyphs[button.dataset.tool==='bond'?'bond:'+button.dataset.value:button.dataset.tool];
    if(!path)continue;
    const icon=sketcher.node('svg',{viewBox:'0 0 22 22',class:'classic-only chemical-glyph','aria-hidden':'true'});
    icon.append(sketcher.node('path',{d:path,fill:'none',stroke:'currentColor','stroke-width':1.7}));
    const caption=document.createElement('span');caption.className='standard-tool-text';caption.textContent=button.textContent;
    if(!button.hasAttribute('aria-label'))button.setAttribute('aria-label',button.title);
    button.replaceChildren(caption,icon);
  }
  $("classic-switch").addEventListener("click",()=>setMode(!dialog.classList.contains("classic")));
  dialog.addEventListener("cofdraw:open-sketch",()=>{
    if(!dialog.dataset.layoutChosen){setMode(true);dialog.dataset.layoutChosen="true";}
  });
  $("launch-classic").addEventListener("click",()=>{setMode(true);openSketch("a");});
  if(new URLSearchParams(location.search).get("sketch")==="classic"){
    if(state.online)$("launch-classic").click();
    else window.addEventListener("cofdraw:library",()=>$("launch-classic").click(),{once:true});
  }
  $("classic-partner").addEventListener("click",async()=>{
    try{const next=activeBlock==="a"?"b":"a";await applySketch();await openSketch(next);}
    catch(error){sketchMessage(error.message,true);}
  });
  $("classic-build").addEventListener("click",async()=>{
    try{
      await applySketch();
      if($("build").disabled)notice($("build").title||"Choose a compatible partner and connectivity to draw this COF.",true);
      else $("build").click();
    }catch(error){sketchMessage(error.message,true);}
  });
  dialog.addEventListener("close",()=>{arm(false);openShelf(false);});
  dialog.addEventListener("cancel",e=>{if(sketchBusy)e.preventDefault();});
  canvas.addEventListener("sketchtool",()=>arm(false));
  function chooseAtom(index){
    anchor.value=String(index);
    view.selected={kind:"atom",index};view.draw();
  }
  async function showTemplate(){
    const token=++revision,row=rows[Number(list.value)];
    arm(false);fragment=null;anchor.replaceChildren();view.clear();syncControls();
    if(!row)return;
    try{
      const result=await templateData(row);
      if(token!==revision)return;
      view.load(result.graph);view.history=[];view.future=[];
      // Keep a small group legible at the scale of the shelf.
      const atoms=view.graph.atoms;
      const span=key=>Math.max(...atoms.map(a=>a[key]))-Math.min(...atoms.map(a=>a[key]));
      view.scale=Math.min(76,520/Math.max(1,span("x")),250/Math.max(1,span("y")));
      fragment=result.graph;
      atoms.forEach((atom,i)=>anchor.add(new Option(atom.element+" "+(i+1),String(i))));
      chooseAtom(0);syncControls();
    }catch(e){if(token===revision)sketchMessage(e.message,true);}
  }
  function populate(){
    const previous=list.selectedOptions[0]?.textContent;
    const all=category.value==="COF monomers"?(state.library?.monomer_presets||[]).map(m=>({
      name:L.displayName(m),smiles:m.smiles,aliases:(m.aliases||[]).join(" ")})):collections.get(category.value)||[];
    const query=search.value.trim().toLowerCase();
    rows=all.filter(r=>(r.name+" "+(r.aliases||"")).toLowerCase().includes(query));
    list.replaceChildren(...rows.map((r,i)=>new Option(r.name,String(i))));
    if(rows.length)list.selectedIndex=Math.max(0,rows.findIndex(r=>r.name===previous));
    pageStart=0;drawTiles();
    showTemplate();
  }
  list.addEventListener("change",showTemplate);
  category.addEventListener("change",async()=>{
    search.value="";const key=category.value;populate();
    if(!key.startsWith("installed:")||collections.has(key))return;
    try{
      const result=await api("templates/installed/"+encodeURIComponent(key.slice(10)));
      collections.set(key,shelfRows(result));
      if(category.value===key){populate();sketchMessage(result.templates.length+" molecular entries from your local ChemDraw installation.");}
    }catch(error){if(category.value===key)sketchMessage(error.message,true);}
  });
  search.addEventListener("input",populate);
  window.addEventListener("cofdraw:library",()=>{if(dialog.classList.contains("classic")&&category.value==="COF monomers")populate();});
  anchor.addEventListener("change",()=>chooseAtom(Number(anchor.value)));
  preview.addEventListener("pointerdown",e=>{
    e.preventDefault();e.stopImmediatePropagation();
    if(fragment && e.button===0 && !sketchBusy){
      const i=view.atomAt(view.position(e));if(i>=0)chooseAtom(i);
      preview.setPointerCapture(e.pointerId);
    }
  },true);
  preview.addEventListener("pointerup",e=>{
    if(!preview.hasPointerCapture(e.pointerId))return;
    preview.releasePointerCapture(e.pointerId);
    const rect=canvas.getBoundingClientRect();
    if(fragment&&e.clientX>=rect.left&&e.clientX<=rect.right&&e.clientY>=rect.top&&e.clientY<=rect.bottom)placeOnSheet(e);
  });
  // The preview is a site picker. Only the main sheet edits molecules.
  $("classic-place").addEventListener("click",()=>{
    if(!fragment||sketchBusy)return;
    if(sketchTextPending){sketchMessage("Load your edited SMILES before placing a template.",true);return;}
    arm(!armed);
    if(armed){openShelf(false);sketchMessage("Click to place. Point at an atom to attach. Esc cancels.");}
    canvas.focus();
  });
  function combine(source, site, target, point){
    const graph=sketcher.exportGraph(),offset=graph.atoms.length;
    if(offset+source.atoms.length>300)throw new Error("The editor limit is 300 atoms.");
    let angle=0,base=point;
    if(target>=0){
      sketcher.connectivityEditable();
      view.connectivityEditable();
      const direction=sketcher.openAngle(target);
      base=sketcher.endAt(graph.atoms[target],direction);
      angle=direction+Math.PI-view.openAngle(site);
    }
    const pivot=source.atoms[site],c=Math.cos(angle),s=Math.sin(angle);
    for(const atom of source.atoms){
      const dx=atom.x-pivot.x,dy=atom.y-pivot.y;
      graph.atoms.push({...atom,x:base.x+dx*c-dy*s,y:base.y+dx*s+dy*c});
    }
    for(const b of source.bonds)graph.bonds.push({...b,a:b.a+offset,b:b.b+offset,
      stereo_atoms:(b.stereo_atoms||[]).map(i=>i+offset)});
    if(target>=0){
      // Implicit H adjusts at sanitization. Replace one bracket-count H at a
      // chosen attachment site as well, e.g. the Si–H cap on trimethylsilane.
      for(const index of [target,site+offset]){
        if(graph.atoms[index].explicit_hs>0)graph.atoms[index].explicit_hs--;
      }
      graph.bonds.push({a:target,b:site+offset,order:1});
    }
    return graph;
  }
  canvas.addEventListener("pointermove",e=>{
    if(!armed||!fragment)return;
    const point=sketcher.position(e),target=sketcher.atomAt(point);
    canvas.querySelector(".placement-preview")?.remove();
    let graph;
    try{graph=combine(fragment,Number(anchor.value),target,point);}catch{return;}
    const ghost=sketcher.node("g",{class:"placement-preview",opacity:.65,"pointer-events":"none"});
    for(const b of graph.bonds.slice(sketcher.graph.bonds.length)){
      const [x1,y1]=sketcher.screen(graph.atoms[b.a]),[x2,y2]=sketcher.screen(graph.atoms[b.b]);
      ghost.append(sketcher.node("line",{x1,y1,x2,y2,stroke:"#2679ba","stroke-width":2}));
    }
    graph.atoms.slice(sketcher.graph.atoms.length).forEach(a=>{
      if(a.element==="C"&&fragment.atoms.length>1)return;
      const [x,y]=sketcher.screen(a);ghost.append(sketcher.node("text",{x,y,fill:"#1764a4","text-anchor":"middle"},a.element));
    });
    canvas.append(ghost);
  });
  async function placeOnSheet(event){
    if(sketchBusy)return;
    if(sketchTextPending){sketchMessage("Load your edited SMILES before placing a template.",true);return;}
    try{
      const point=sketcher.position(event),target=sketcher.atomAt(point);
      const graph=combine(fragment,Number(anchor.value),target,point);
      const center=[...sketcher.center],scale=sketcher.scale;
      // A worksheet can hold several fragments while the chemist joins them.
      // Check valence on attached products before committing the edit.
      if(target<0&&sketcher.graph.atoms.length){
        sketcher.load(graph);
        sketchMessage("Template placed. Join the fragments before applying this monomer.");
      }else{
        try{await sketchRequest({graph});}
        catch(error){
          if(!error.message.includes("one connected structure"))throw error;
          sketcher.load(graph);
        }
        sketchMessage("Template placed. Undo removes this edit.");
      }
      sketcher.center=center;sketcher.scale=scale;sketcher.draw();arm(false);
    }catch(e){sketchMessage(e.message,true);}
    finally{syncControls();}
  }
  canvas.addEventListener("pointerdown",event=>{
    if(!armed||event.button!==0)return;
    event.preventDefault();event.stopImmediatePropagation();placeOnSheet(event);
  },true);
  dialog.addEventListener("keydown",e=>{
    if(e.key==="Escape"&&armed){e.preventDefault();e.stopPropagation();arm(false);sketchMessage("Placement cancelled.");}
  },true);
  $("classic-import").addEventListener("click",()=>$("classic-file").click());
  $("classic-file").addEventListener("change",async e=>{
    const file=e.target.files[0];e.target.value="";if(!file)return;
    try{
      if(file.size>500000)throw new Error("Import molecular template sheets up to 500 kB.");
      let result;
      if(/\.(ctp|cdx)$/i.test(file.name)){
        const bytes=new Uint8Array(await file.arrayBuffer());
        let binary="";for(let i=0;i<bytes.length;i+=8192)binary+=String.fromCharCode(...bytes.subarray(i,i+8192));
        result=await api("templates/binary",{content_base64:btoa(binary)});
      }else result=await api("templates/cdxml",{cdxml:await file.text()});
      const added=shelfRows(result);
      const imported=collections.get("Imported");imported.push(...added);
      category.value="Imported";search.value="";populate();
      list.selectedIndex=imported.length-added.length;await showTemplate();
      sketchMessage(added.length+" structures imported from "+file.name+(result.skipped?" · "+result.skipped+" unsupported fragments skipped.":"."));
    }catch(error){sketchMessage(error.message,true);}
  });
  function zoom(factor){sketcher.scale=Math.max(5,Math.min(120,sketcher.scale*factor));sketcher.draw();}
  const zoomLevel=$('classic-zoom-level');
  zoomLevel.onchange=()=>{if(zoomLevel.value!=='current'){sketcher.scale=35*parseInt(zoomLevel.value,10)/100;sketcher.draw();}};
  canvas.addEventListener('sketchzoom',()=>{zoomLevel.options[0].textContent=Math.round(sketcher.scale/35*100)+'%';zoomLevel.value='current';});
  $("classic-minus").onclick=()=>zoom(1/1.25);$("classic-plus").onclick=()=>zoom(1.25);
  $("classic-fit").onclick=()=>{sketcher.fit();sketcher.draw();};
  shelf.addEventListener("keydown",e=>{if(e.key==="Enter"&&e.target===list)$("classic-place").click();});
  const menus=dialog.querySelector(".classic-menus");
  function closeMenus(keep){
    menus.querySelectorAll("details[open]").forEach(menu=>{if(menu!==keep)menu.open=false;});
  }
  dialog.addEventListener("pointerdown",e=>closeMenus(e.target.closest(".classic-menus details")));
  menus.addEventListener("click",e=>{
    const command=e.target.closest("[data-command]")?.dataset.command;
    if(!command)return;
    closeMenus();arm(false);
    try{
      if(command==="all"){sketcher.selectAll();canvas.focus();return;}
      if(command==="delete"){sketcher.removeSelected();canvas.focus();return;}
      if(command==="keep"||command==="copy-smiles"){
        (async()=>{
          const mol=await checkedSketch();
          if(command==="copy-smiles"){
            try{await navigator.clipboard.writeText(mol.smiles);sketchMessage("SMILES copied.");}
            catch{const field=$("sketch-smiles");field.value=mol.smiles;field.focus();field.select();sketchMessage("Press Ctrl+C or Command+C to copy the selected SMILES.");}
          }else{
            const imported=collections.get("Imported");imported.push({name:"Drawing "+(imported.length+1)+" · "+mol.formula,molblock:mol.molblock});
            category.value="Imported";search.value="";populate();list.selectedIndex=imported.length-1;await showTemplate();
            sketchMessage("Drawing added to this session’s shelf. Save MOL to keep a copy.");
          }
        })().catch(error=>sketchMessage(error.message,true));
        return;
      }
      if(command==="rotate"){
        const sel=sketcher.selected;
        const ids=sel?.kind==="group"?sel.indices:sel?.kind==="bond"?[sketcher.graph.bonds[sel.index].a,sketcher.graph.bonds[sel.index].b]:[];
        if(ids.length<2)throw new Error("Drag a selection box around the atoms to rotate, or choose Edit → Select all atoms.");
        const atoms=ids.map(i=>sketcher.graph.atoms[i]);
        const cx=atoms.reduce((n,a)=>n+a.x,0)/ids.length,cy=atoms.reduce((n,a)=>n+a.y,0)/ids.length;
        const c=Math.cos(Math.PI/12),s=Math.sin(Math.PI/12);
        sketcher.remember();
        for(const a of atoms){const dx=a.x-cx,dy=a.y-cy;a.x=cx+dx*c-dy*s;a.y=cy+dx*s+dy*c;}
        sketcher.changed();canvas.focus();return;
      }
      $(command).click();
    }catch(error){sketchMessage(error.message,true);}
  });
})();
