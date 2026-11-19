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
  }
  function syncControls(){
    const usable=!!fragment && !sketchBusy;
    anchor.disabled=!usable;
    $("classic-place").disabled=!usable;
    for(const el of [category,list,search])el.disabled=sketchBusy;
  }
  function setMode(enabled){
    arm(false);sketcher.cancelDrag();
    sketcher.classicMode=enabled;
    if(!enabled && sketcher.selected?.kind==="group"){sketcher.selected=null;sketcher.draw();}
    if(!enabled&&["pan","chain"].includes(sketcher.tool))sketcher.setTool("bond","1");
    dialog.classList.toggle("classic",enabled);
    $("classic-switch").textContent=enabled?"Standard layout":"Try classic layout";
    $("classic-switch").setAttribute("aria-pressed",String(enabled));
    if(enabled){discoverInstalled();populate();}
    else {revision++;fragment=null;}
  }
  $("classic-switch").addEventListener("click",()=>setMode(!dialog.classList.contains("classic")));
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
  dialog.addEventListener("close",()=>arm(false));
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
      const result=await api("monomer",row.molblock?{molblock:row.molblock}:{smiles:row.smiles});
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
    if(armed)sketchMessage(sketcher.graph.atoms.length?"Click the drawing atom to attach with a single bond. Esc cancels.":"Click the sheet to place the structure.");
    canvas.focus();
  });
  function combine(source, site, target, point){
    const graph=sketcher.exportGraph(),offset=graph.atoms.length;
    if(offset+source.atoms.length>300)throw new Error("The editor limit is 300 atoms.");
    let angle=0,base=point;
    if(offset){
      if(target<0)throw new Error("Click an existing atom to attach the template.");
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
    if(offset){
      // Implicit H adjusts at sanitization. Replace one bracket-count H at a
      // chosen attachment site as well, e.g. the Si–H cap on trimethylsilane.
      for(const index of [target,site+offset]){
        if(graph.atoms[index].explicit_hs>0)graph.atoms[index].explicit_hs--;
      }
      graph.bonds.push({a:target,b:site+offset,order:1});
    }
    return graph;
  }
  async function placeOnSheet(event){
    if(sketchBusy)return;
    if(sketchTextPending){sketchMessage("Load your edited SMILES before placing a template.",true);return;}
    try{
      const point=sketcher.position(event),target=sketcher.atomAt(point);
      const graph=combine(fragment,Number(anchor.value),target,point);
      arm(false);
      await sketchRequest({graph});
      sketchMessage("Template placed. Undo removes this edit.");
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
