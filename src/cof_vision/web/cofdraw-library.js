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
      const text=[item.abbreviation,item.full_name,item.label,item.id,...list(item.literature_names)].filter(Boolean).join(" ").toLowerCase();
      return tokens.every(token=>text.includes(token));
    });
  }
  function monomerChoice(library,item,linkage,block,other=null) {
    if(block==="editor")return {compatible:true,reason:"Open this single monomer for editing or chemical export; no partner is required."};
    const rule=library.linkages?.[linkage]||{},role=roleFor(library,linkage,block),label=rule.role_labels?.[block==="a"?0:1]||role?.replaceAll("_"," ")||"supported reactive group";
    if(rule.supported===false||!role)return {compatible:false,reason:"This linkage is reference-only. Open the monomer in the editor to inspect or save it."};
    const evidence=roleEvidence(item,linkage,role);
    if(evidence?.compatible===false)return {compatible:false,reason:evidence.reactive_sites===undefined?`Requires ${label}; this monomer has no supported ${label.toLowerCase()} role.`:evidence.reason};
    const raw=library.templates,templates=Array.isArray(raw)?raw:raw?Object.entries(raw).map(([id,t])=>({id,...t})):[];
    const otherRole=roleFor(library,linkage,block==="a"?"b":"a");
    if(other&&roleEvidence(other,linkage,otherRole)?.compatible===true&&templates.length){
      const pair=block==="a"?{a:item,b:other}:{a:other,b:item};
      const eligible=templates.filter(t=>templateCompatibility(t,linkage,pair,library).compatible);
      if(!eligible.length)return {compatible:false,reason:"No available connectivity template supports this pair of reactive-site counts and geometry."};
      if(eligible.every(t=>!t.topology_label))return {compatible:true,reason:`Supplies ${label}. Only an unassigned finite projection is available for this pair.`};
    }
    return {compatible:true,reason:`Supplies ${label}. Attachment incidence and the selected template are checked during construction.`};
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
  const api = {displayName, roleFor, roleEvidence, compatiblePresets, groupRows, reactiveCount, templateCompatibility, stylesInCategory,connectionCount,searchMonomers,monomerChoice,createPreviewCache};
  scope.COFdrawLibrary = api;
  if (typeof module !== "undefined" && module.exports) module.exports = api;
})(typeof window === "undefined" ? globalThis : window);
