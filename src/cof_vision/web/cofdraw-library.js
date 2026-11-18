/* Metadata-driven controls. Chemical recognition belongs to the local engine. */
"use strict";
(function (scope) {
  const list = value => Array.isArray(value) ? value : [];
  function displayName(item) {
    const abbreviation = item.abbreviation || "";
    const name = item.full_name || item.name || "";
    if (abbreviation && name && abbreviation !== name) return `${abbreviation} — ${name}`;
    return name || item.label || abbreviation || item.id || "Unnamed monomer";
  }
  function roleFor(library, linkage, block) {
    return library.linkages?.[linkage]?.roles?.[block === "a" ? 0 : 1] || null;
  }
  function roleEvidence(item, linkage, role) {
    if (!item || !role) return null;
    const rows = list(item.compatible_roles);
    const match = rows.find(row => row.linkage === linkage && row.role === role);
    if (match) {
      const count = match.reactive_sites ?? match.count;
      if (match.compatible !== false && Number.isInteger(count) && count < 2) return {...match,compatible:false,reason:`This role has ${count} reactive site${count===1?"":"s"}; the available patch templates require at least two.`};
      return {...match, compatible: match.compatible !== false};
    }
    // An explicitly supplied (even empty) list is authoritative.
    if (Array.isArray(item.compatible_roles)) return {compatible: false, reason: "No supported transformation for this role."};
    return null;
  }
  function compatiblePresets(library, linkage, block) {
    const role = roleFor(library, linkage, block);
    return list(library.monomer_presets).filter(item => {
      const evidence = roleEvidence(item, linkage, role);
      return evidence ? evidence.compatible : item.annotation_only !== true;
    });
  }
  function groupRows(item) {
    return list(item?.functionality || item?.reactive_groups).filter(row => Number(row.count) > 0);
  }
  function reactiveCount(item, linkage, role) {
    const evidence = roleEvidence(item, linkage, role);
    const count = evidence?.reactive_sites ?? evidence?.count;
    return Number.isInteger(count) && count >= 0 ? count : null;
  }
  function templateCompatibility(template, linkage, monomers, library) {
    if (template.supported === false) return {compatible: false, reason: template.limits || "Template unavailable."};
    if (Array.isArray(template.linkages) && !template.linkages.includes(linkage)) return {compatible: false, reason: "Not available for this linkage."};
    if (template.topology_label && [monomers.a, monomers.b].some(item => item?.geometry_hint?.startsWith("tetrahedral"))) return {compatible:false, reason:"A tetrahedral node needs explicit 3D incidence. Use an unassigned projection rather than a planar net template."};
    for (const block of ["a", "b"]) {
      const role = roleFor(library, linkage, block);
      const evidence = roleEvidence(monomers[block], linkage, role);
      if (evidence?.compatible === false) return {compatible: false, reason: `Block ${block.toUpperCase()}: ${evidence.reason || "incompatible reactive group"}`};
      const allowed = template.required_functionalities?.[block];
      const count = reactiveCount(monomers[block], linkage, role);
      if (count !== null && Array.isArray(allowed) && !allowed.includes(count)) return {compatible: false, reason: `Block ${block.toUpperCase()} has ${count} reactive sites; this template requires ${allowed.join(" or ")}.`};
    }
    const pairs = template.allow_pairs || template.required_functionalities?.allow_pairs;
    const a = reactiveCount(monomers.a, linkage, roleFor(library, linkage, "a"));
    const b = reactiveCount(monomers.b, linkage, roleFor(library, linkage, "b"));
    if (Array.isArray(pairs) && pairs.length && a !== null && b !== null && !pairs.some(pair => pair[0] === a && pair[1] === b)) return {compatible: false, reason: `This template accepts ${pairs.map(pair => pair.join(" + ")).join(" or ")} reactive sites; these monomers have ${a} + ${b}.`};
    return {compatible: true, reason: "The engine checks attachment incidence when the drawing is created."};
  }
  function stylesInCategory(styles, category) {
    return list(styles).filter(style => category === "all" || (style.category || "chemical") === category);
  }
  function connectionCount(item) {
    if (Number.isInteger(item?.connection_count)) return item.connection_count;
    const counts=groupRows(item).map(row=>row.count).filter(Number.isInteger);
    return counts.length?Math.max(...counts):null;
  }
  function searchMonomers(rows,query="",count="all") {
    const tokens=query.toLowerCase().trim().split(/\s+/).filter(Boolean);
    return list(rows).filter(item=>{
      if(count!=="all"&&connectionCount(item)!==Number(count))return false;
      const text=[item.abbreviation,item.full_name,item.label,item.id,...list(item.literature_names),...list(item.aliases)].filter(Boolean).join(" ").toLowerCase();
      return tokens.every(token=>text.includes(token));
    });
  }
  function monomerChoice(library,item,linkage,block,other=null) {
    if(block==="editor")return {compatible:true,reason:"Open this single monomer for editing or chemical export; no partner is required."};
    const rule=library.linkages?.[linkage]||{},role=roleFor(library,linkage,block),label=rule.role_labels?.[block==="a"?0:1]||role?.replaceAll("_"," ")||"supported reactive group";
    if(rule.supported===false||!role)return {compatible:false,reason:"This linkage is reference-only. Open the monomer in the editor to inspect or save it."};
    const evidence=roleEvidence(item,linkage,role);
    if(!evidence)return {compatible:false,reason:"Inspect this monomer before using it in a construction."};
    if(evidence?.compatible===false)return {compatible:false,reason:evidence.reactive_sites===undefined?`Requires ${label}; this monomer has no supported ${label.toLowerCase()} role.`:evidence.reason};
    const raw=library.templates,templates=Array.isArray(raw)?raw:raw?Object.entries(raw).map(([id,t])=>({id,...t})):[];
    const otherRole=roleFor(library,linkage,block==="a"?"b":"a");
    if(block==='b'&&other&&roleEvidence(other,linkage,otherRole)?.compatible!==true)return {compatible:false,reason:"Building block A does not supply the selected linkage role. Change A or enable automatic chemistry matching."};
    if(other&&roleEvidence(other,linkage,otherRole)?.compatible===true&&templates.length){
      const pair=block==="a"?{a:item,b:other}:{a:other,b:item};
      const eligible=templates.filter(t=>templateCompatibility(t,linkage,pair,library).compatible);
      if(!eligible.length)return {compatible:false,reason:"No available connectivity template supports this pair of reactive-site counts and geometry."};
      if(eligible.every(t=>!t.topology_label))return {compatible:true,reason:`Supplies ${label}. Only an unassigned finite projection is available for this pair.`};
    }
    return {compatible:true,reason:`Supplies ${label}. Attachment incidence and the selected template are checked during construction.`};
  }
  function compatibleLinkages(library, a, b) {
    if (!a || !b) return [];
    const raw=library.templates, templates=Array.isArray(raw)?raw:Object.values(raw||{}), results=[];
    for (const [id, rule] of Object.entries(library.linkages||{})) {
      if (rule.supported!==true || !Array.isArray(rule.roles) || rule.roles.length!==2) continue;
      for (const swapped of [false,true]) {
        const pair=swapped?{a:b,b:a}:{a,b};
        if (!['a','b'].every((block,index)=>roleEvidence(pair[block],id,rule.roles[index])?.compatible===true)) continue;
        const eligible=templates.filter(t=>templateCompatibility(t,id,pair,library).compatible);
        if (templates.length && !eligible.length) continue;
        results.push({id,swapped,templates:eligible.map(t=>t.id),label:rule.label||id});
        break;
      }
    }
    // Product-state reduction is never an automatic replacement for condensation.
    const order=['beta_ketoenamine','keto_hydrazide','imine','acylhydrazone','azine','nitrone','secondary_amine'];
    return results.sort((a,b)=>(order.indexOf(a.id)<0?99:order.indexOf(a.id))-(order.indexOf(b.id)<0?99:order.indexOf(b.id)));
  }
  function automaticChoice(library,item,block,other) {
    if(block==='editor') return {compatible:true,reason:'Open, inspect, or edit this monomer.'};
    const pair=block==='a'?[item,other]:[other,item], choices=compatibleLinkages(library,...pair);
    if(choices.length) return {compatible:true,reason:`Compatible with this partner: ${choices.map(c=>c.label).join(', ')}. Chemistry and connectivity update when selected.`};
    if(block==='b'&&other)return {compatible:false,reason:'No implemented linkage and connectivity template connects this monomer with building block A. Choose a different partner or inspect it in the single-monomer editor.'};
    const supported=list(item.compatible_roles).some(r=>library.linkages?.[r.linkage]?.supported===true&&r.compatible!==false&&(r.reactive_sites??r.count)>=2);
    return {compatible:supported,reason:supported?'Choose this monomer, then a compatible partner. Available chemistry will update automatically.':'No implemented transform for this monomer yet; use the single-monomer editor or import a reviewed source drawing.'};
  }
  function automaticPlan(library,a,b,preferred=null) {
    const choices=compatibleLinkages(library,a,b);
    // A reduced product state always requires an explicit manual choice.
    const allowed=choices.filter(choice=>choice.id!=='secondary_amine');
    const chosen=allowed.find(choice=>choice.id===preferred?.linkage)||allowed[0];
    if(!chosen)return null;
    const raw=library.templates,templates=Array.isArray(raw)?raw:Object.values(raw||{});
    const eligible=templates.filter(template=>chosen.templates.includes(template.id));
    const template=eligible.find(item=>item.id===preferred?.template_id)||eligible.find(item=>item.allow_pairs?.length)||eligible[0];
    return {...chosen,template_id:template?.id||null,choice_count:choices.length};
  }
  function startingCombinations(library) {
    const blocks=library.monomer_presets||[], byId=new Map(blocks.map(b=>[b.id,b]));
    const sources=new Map((library.source_recipes||[]).map(s=>[s.recipe_id,s]));
    const searchText=(...parts)=>{
      const terms=parts.flat(Infinity).filter(Boolean).join(' ').toLowerCase();
      return terms+' '+terms.replace(/[-‐‑–·]/g,'');
    };
    const names=b=>b?[b.id,b.abbreviation,b.full_name,b.aliases||[],b.literature_names||[]]:[];
    const pairKey=(a,b,link)=>JSON.stringify([...[a,b].sort(),link]);
    const occupied=new Set(), rows=[];
    for(const preset of library.cof_presets||[]){
      const a=byId.get(preset.monomer_a_id), b=byId.get(preset.monomer_b_id);
      const name=sources.get(preset.source_recipe_id)?.literature_name;
      const pair=preset.label_a+' + '+preset.label_b;
      const label=name?name+' · '+pair:preset.label;
      rows.push({...preset,label,group:name?'Named starting pairs':'Catalog designs',
        search:searchText(label,preset.linkage,library.linkages[preset.linkage]?.label,names(a),names(b))});
      if(a&&b)occupied.add(pairKey(a.canonical_smiles||a.smiles,b.canonical_smiles||b.smiles,preset.linkage));
    }
    for(const source of sources.values()){
      const label=source.literature_name||source.label;
      rows.push({id:'source:'+source.recipe_id,source,group:'Recorded drawings',label,
        search:searchText(label,source.label,source.recipe_id,(source.monomers||[]).map(names),(source.monomers||[]).map(m=>m.label))});
    }
    const templates=new Map(Object.values(library.templates||{}).map(t=>[t.id,t]));
    const reactive=blocks.filter(b=>b.compatible_roles?.length&&!b.catalog_use?.startsWith('reference_only'));
    const generated=[];
    for(let i=0;i<reactive.length;i++)for(let j=i+1;j<reactive.length;j++){
      const a=reactive[i], b=reactive[j], sa=a.canonical_smiles||a.smiles, sb=b.canonical_smiles||b.smiles;
      if(sa===sb)continue;
      const plan=automaticPlan(library,a,b), template=templates.get(plan?.template_id);
      if(!template?.allow_pairs?.length)continue;
      const key=pairKey(sa,sb,plan.id);
      if(occupied.has(key))continue;
      occupied.add(key);
      const [left,right]=plan.swapped?[b,a]:[a,b];
      const labelA=left.abbreviation||left.id, labelB=right.abbreviation||right.id;
      const label=labelA+' + '+labelB;
      generated.push({id:'pair:'+left.id+':'+right.id+':'+plan.id,label,label_a:labelA,label_b:labelB,
        monomer_a:left.smiles,monomer_b:right.smiles,linkage:plan.id,template_id:template.id,instances:template.instances,
        group:(template.id==='chain_2_2'?'Open-chain designs · ':'Matched pairs · ')+plan.label,
        search:searchText(label,plan.id,plan.label,names(left),names(right),template.label)});
    }
    return rows.concat(generated.sort((a,b)=>a.group.localeCompare(b.group)||a.label.localeCompare(b.label)));
  }
  // Share in-flight work as well as completed previews. Only the caller decides
  // which visible cards to request; a long catalog never starts all at once.
  function createPreviewCache(load, concurrency=3) {
    if (!Number.isInteger(concurrency) || concurrency < 1) throw new Error("Preview concurrency must be positive.");
    const entries=new Map(),queue=[];let active=0;
    function pump(){
      while(active<concurrency&&queue.length){
        const job=queue.shift();active++;
        Promise.resolve().then(()=>load(job.key)).then(job.resolve,error=>{entries.delete(job.key);job.reject(error);}).finally(()=>{active--;pump();});
      }
    }
    function request(key){
      if(entries.has(key))return entries.get(key);
      const promise=new Promise((resolve,reject)=>queue.push({key,resolve,reject}));
      entries.set(key,promise);pump();return promise;
    }
    return {request};
  }
  const api = {displayName, roleFor, roleEvidence, compatiblePresets, groupRows, reactiveCount, templateCompatibility, stylesInCategory,connectionCount,searchMonomers,monomerChoice,compatibleLinkages,automaticChoice,automaticPlan,startingCombinations,createPreviewCache};
  scope.COFdrawLibrary = api;
  if (typeof module !== "undefined" && module.exports) module.exports = api;
})(typeof window === "undefined" ? globalThis : window);
