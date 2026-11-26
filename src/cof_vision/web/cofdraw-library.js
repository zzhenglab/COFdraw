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
    if(item.binary_eligible===false)return {compatible:false,reason:"Outside the supported two-monomer construction rules."};
    const rows = list(item.binary_roles || item.compatible_roles);
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
    const blocks=['a','b'];
    const roles=blocks.map(block=>roleFor(library,linkage,block));
    const shapes=blocks.map((block,i)=>monomers[block]?.site_shapes?.[roles[i]]);
    const directions=blocks.map((block,i)=>monomers[block]?.attachment_directions?.[roles[i]]);
    const tetrahedral=blocks.map((block,i)=>directions[i]==='tetrahedral'||shapes[i]==='tetrahedral'||(!directions[i]&&!shapes[i]&&monomers[block]?.geometry_hint?.startsWith('tetrahedral')));
    const molecularMatch=(template.id==='pts_4_4'&&directions.includes('square_planar')&&directions.includes('tetrahedral'))||
      (template.id==='dia_4_2'&&directions.some((dir,i)=>dir==='linear'&&tetrahedral[1-i]));
    if(template.topology_label&&shapes.includes('unresolved_alkyl'))return {compatible:false,reason:'The attachment directions of this aliphatic building block are unassigned. Choose a local projection.'};
    if(['pts_4_4','dia_4_2'].includes(template.id)&&!template.shape_pair&&!molecularMatch)
      return {compatible:false,reason:'The attachment directions do not match this spatial net.'};
    if(template.shape_pair&&!molecularMatch&&!template.shape_pair.every((value,i)=>value===shapes[i])&&
       !template.shape_pair.every((value,i)=>value===shapes[1-i]))return {compatible:false,reason:'The reactive arms do not match this spatial arrangement.'};
    if(template.topology_label&&(shapes.includes('spiro')||shapes.includes('spatial_unspecified')))return {compatible:false,reason:'This node has out-of-plane arms. Its geometry needs a separate connectivity assignment.'};
    if(template.topology_label&&shapes.includes('prismatic'))return {compatible:false,reason:'The six arms on this triptycene core form a prism. Choose a spatial net.'};
    if(template.id==='branched_4_3'&&![monomers.a,monomers.b].some(m=>m?.geometry_hint==='tetrahedral'))return {compatible:false,reason:'Requires a tetrahedral four-connected monomer.'};
    if (template.supported === false) return {compatible: false, reason: template.limits || "Template unavailable."};
    if (Array.isArray(template.linkages) && !template.linkages.includes(linkage)) return {compatible: false, reason: "Not available for this linkage."};
    if(template.id==='kgm_4_2'){
      const kagomeRoles=blocks.map((block,i)=>monomers[block]?.kagome_roles?.[roles[i]]);
      if(!directions.includes('linear')&&!(kagomeRoles.includes('node')&&kagomeRoles.includes('spacer')))
        return {compatible:false,reason:'Choose a four-connected node and a linear two-connected spacer.'};
    }
    if(template.topology_label&&template.dimension!==3&&!molecularMatch&&tetrahedral.some(Boolean))
      return {compatible:false,reason:'A tetrahedral node requires a spatial template.'};
    for (const block of ["a", "b"]) {
      const role = roleFor(library, linkage, block);
      const evidence = roleEvidence(monomers[block], linkage, role);
      if (evidence?.compatible === false) return {compatible: false, reason: `Block ${block.toUpperCase()}: ${evidence.reason || "incompatible reactive group"}`};
      const allowed = template.required_functionalities?.[block];
      const count = reactiveCount(monomers[block], linkage, role);
      if (count !== null && Array.isArray(allowed) && !allowed.includes(count)) return {compatible: false, reason: `Block ${block.toUpperCase()} has ${count} reactive sites; this template requires ${allowed.join(" or ")}.`};
    }
    const sizes=['a','b'].map(block=>roleEvidence(monomers[block],linkage,roleFor(library,linkage,block))?.fragment_atoms);
    // A large mtf pair opens on a smaller square when the full neighborhood
    // is too large. Require enough room for its smallest drawable face.
    const minimumBlocks=template.id==='mtf_4_3'?4:template.instances;
    if(sizes.every(Number.isInteger)&&sizes.reduce((a,b)=>a+b)*minimumBlocks/2>1000)return {compatible:false,reason:'This patch exceeds the 1,000-atom drawing limit.'};
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
  function linkageGroupOptions(rows, groupsFor = row=>list(row.linkage_functional_groups)) {
    const groups = new Map();
    for (const row of rows) {
      for (const group of groupsFor(row)) groups.set(group.id, group);
    }
    return [...groups.values()].sort((a,b) => a.label.localeCompare(b.label));
  }
  function searchMonomers(rows,query="",count="all",group="all",groupsFor = row=>list(row.linkage_functional_groups)) {
    const tokens=query.toLowerCase().trim().split(/\s+/).filter(Boolean);
    return list(rows).filter(item=>{
      if(group!=="all" && !groupsFor(item).some(entry=>entry.id===group))return false;
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
      if(eligible.every(t=>!t.topology_label&&!t.parent_net))return {compatible:true,reason:`Supplies ${label}. Only an unassigned finite projection is available for this pair.`};
    }
    return {compatible:true,reason:`Supplies ${label}. Attachment incidence and the selected template are checked during construction.`};
  }
  function compatibleLinkages(library, a, b) {
    if (!a || !b || a.binary_eligible===false || b.binary_eligible===false) return [];
    if(a.canonical_smiles && a.canonical_smiles===b.canonical_smiles)return [];
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
    const order=['beta_ketoenamine','alkyl_ketoenamine','keto_hydrazide','imine','alkyl_imine','acylhydrazone','azine','nitrone','secondary_amine'];
    return results.sort((a,b)=>(order.indexOf(a.id)<0?99:order.indexOf(a.id))-(order.indexOf(b.id)<0?99:order.indexOf(b.id)));
  }
  function periodicLinkageMatches(library,a,b,selected,automatic) {
    const matches = [], supported = ['imine','secondary_amine','acylhydrazone'];
    if (!a || !b) return matches;
    for (const id of new Set([selected,...(automatic ? supported : [])])) {
      const roles = library?.linkages?.[id]?.roles;
      if (!supported.includes(id) || roles?.length !== 2) continue;
      for (const swapped of [false,true]) {
        const pair = swapped ? [b,a] : [a,b];
        const sizes = pair.map((m,i)=>reactiveCount(m,id,roles[i]));
        const hub = pair.map((m,i)=>m.geometry_hint === 'tetrahedral' && sizes[i] === 4);
        if ((hub[0] && (sizes[1]===2 || hub[1])) || (hub[1] && sizes[0]===2)) {
          matches.push({id,swapped});
          break;
        }
      }
    }
    return matches;
  }
  function partnerLinkageGroups(library,a,b,selected,automatic) {
    const groups = list(b?.linkage_functional_groups);
    if (!a || !groups.length) return [];
    const matches = [...compatibleLinkages(library,a,b), ...periodicLinkageMatches(library,a,b,selected,automatic)];
    return groups.filter(group => matches.some(({id,swapped}) => {
      if (!automatic && id !== selected) return false;
      const partnerRole = library.linkages[id].roles[swapped ? 0 : 1];
      return list(group.binary_roles).some(match=>match.linkage===id && match.role===partnerRole);
    }));
  }
  function automaticChoice(library,item,block,other) {
    if(block==='editor') return {compatible:true,reason:'Open, inspect, or edit this monomer.'};
    const pair=block==='a'?[item,other]:[other,item], choices=compatibleLinkages(library,...pair);
    if(choices.length) return {compatible:true,reason:`Compatible with this partner: ${choices.map(c=>c.label).join(', ')}. Chemistry and connectivity update when selected.`};
    if(block==='b'&&other)return {compatible:false,reason:'No implemented linkage and connectivity template connects this monomer with building block A. Choose a different partner or inspect it in the single-monomer editor.'};
    const supported=item.binary_eligible!==false&&list(item.binary_roles||item.compatible_roles).some(r=>library.linkages?.[r.linkage]?.supported===true&&r.compatible!==false&&(r.reactive_sites??r.count)>=2);
    return {compatible:supported,reason:supported?'Choose this monomer, then a compatible partner. Available chemistry will update automatically.':'No implemented transform for this monomer yet; use the single-monomer editor or import a reviewed source drawing.'};
  }
  function automaticPlan(library,a,b,preferred=null) {
    const choices=compatibleLinkages(library,a,b);
    // A reduced product state always requires an explicit manual choice.
    const allowed=choices.filter(choice=>choice.id!=='secondary_amine'&&!library.linkages[choice.id]?.explicit_only);
    const chosen=allowed.find(choice=>choice.id===preferred?.linkage)||allowed[0];
    if(!chosen)return null;
    const raw=library.templates,templates=Array.isArray(raw)?raw:Object.values(raw||{});
    const eligible=templates.filter(template=>chosen.templates.includes(template.id));
      const template=eligible.find(item=>item.id===preferred?.template_id&&item.allow_pairs?.length)||eligible.find(item=>item.id==='kgm_4_2')||eligible.find(item=>item.allow_pairs?.length)||eligible[0];
    return {...chosen,template_id:template?.id||null,choice_count:choices.length};
  }
  function startingCombinations(library) {
    const blocks=library.monomer_presets||[], byId=new Map(blocks.map(b=>[b.id,b]));
    const abbreviationCounts=new Map();
    for(const block of blocks){
      const key=(block.abbreviation||'').toLowerCase();
      abbreviationCounts.set(key,(abbreviationCounts.get(key)||0)+1);
    }
    function pairName(block){
      return abbreviationCounts.get((block.abbreviation||'').toLowerCase())>1
        ? displayName(block) : block.abbreviation||block.id;
    }
    const sources=new Map((library.source_recipes||[]).map(s=>[s.recipe_id,s]));
    const searchText=(...parts)=>{
      const terms=parts.flat(Infinity).filter(Boolean).join(' ').toLowerCase();
      return terms+' '+terms.replace(/[-‐‑–·]/g,'');
    };
    const names=b=>b?[b.id,b.abbreviation,b.full_name,b.aliases||[],b.literature_names||[]]:[];
    const pairKey=(a,b,link)=>JSON.stringify([...[a,b].sort(),link]);
    const occupied=new Set(), rows=[];
    const templates=new Map(Object.values(library.templates||{}).map(t=>[t.id,t]));
    const nameFor=new Map((library.named_pairs||[]).map(n=>[pairKey(n.a,n.b,n.linkage),n]));
    for(const preset of library.cof_presets||[]){
      const a=byId.get(preset.monomer_a_id), b=byId.get(preset.monomer_b_id);
      const named=nameFor.get(pairKey(preset.monomer_a_id,preset.monomer_b_id,preset.linkage));
      const name=sources.get(preset.source_recipe_id)?.literature_name||named?.name;
      if(a&&b&&!compatibleLinkages(library,a,b).some(c=>c.id===preset.linkage&&(!preset.template_id||c.templates.includes(preset.template_id))))continue;
      const pair=preset.label_a+' + '+preset.label_b;
      const label=name?pair+' ('+name+')':pair;
      rows.push({...preset,label,commonName:name,doi:named?.doi,group:name?'Named starting pairs':'Catalog designs',
        search:searchText(label,preset.linkage,library.linkages[preset.linkage]?.label,names(a),names(b))});
      if(a&&b)occupied.add(pairKey(a.canonical_smiles||a.smiles,b.canonical_smiles||b.smiles,preset.linkage));
    }
    for(const source of sources.values()){
      const name=source.literature_name||source.label;
      const pair=(source.monomers||[]).map(m=>m.label).filter(Boolean).join(' + ');
      const label=pair?pair+' ('+name+')':name;
      rows.push({id:'source:'+source.recipe_id,source,commonName:name,group:'Recorded drawings',label,
        search:searchText(label,source.label,source.recipe_id,(source.monomers||[]).map(names),(source.monomers||[]).map(m=>m.label))});
    }
    for(const named of nameFor.values()){
      const a=byId.get(named.a),b=byId.get(named.b),template=templates.get(named.template_id);
      if(!a||!b||!template||!compatibleLinkages(library,a,b).some(c=>c.id===named.linkage&&c.templates.includes(template.id)))continue;
      const key=pairKey(a.canonical_smiles||a.smiles,b.canonical_smiles||b.smiles,named.linkage);
      if(occupied.has(key))continue;
      occupied.add(key);
      const label=pairName(a)+' + '+pairName(b)+' ('+named.name+')';
      rows.push({id:'named:'+named.name,label,commonName:named.name,doi:named.doi,group:'Named starting pairs',label_a:a.abbreviation,label_b:b.abbreviation,
        monomer_a:a.smiles,monomer_b:b.smiles,linkage:named.linkage,template_id:template.id,instances:template.instances,
        search:searchText(label,named.linkage,names(a),names(b),named.doi)});
    }
    const reactive=blocks.filter(b=>b.binary_eligible!==false&&!b.assembly_roles?.self_imine&&b.compatible_roles?.length&&!b.catalog_use?.startsWith('reference_only'));
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
      const label=pairName(left)+' + '+pairName(right);
      generated.push({id:'pair:'+left.id+':'+right.id+':'+plan.id,label,label_a:labelA,label_b:labelB,
        monomer_a:left.smiles,monomer_b:right.smiles,linkage:plan.id,template_id:template.id,instances:template.instances,
        group:(template.id==='chain_2_2'?'Open-chain designs · ':'Matched pairs · ')+plan.label,
        rank:[left,right].reduce((sum,r)=>sum+(r.id.startsWith('cv_')?4:0)+(r.abbreviation?.length||0)/30,0),
        search:searchText(label,plan.id,plan.label,names(left),names(right),template.label)});
    }
    for(const example of library.assembly_examples||[]){
      const parts=example.monomer_ids.map(id=>byId.get(id));
      if(parts.some(r=>!r))continue;
      const name=example.doi?example.label.split(' · ')[0]:null;
      const feed=parts.map(pairName).join(' + '),label=name?feed+' ('+name+')':feed;
      const kind=parts[0].assembly_roles?.self_kind;
      const linkage=kind==='self_boroxine'?'boroxine':kind==='self_triazine'?'triazine':'imine';
      rows.push({id:'assembly:'+example.id,assembly:example,commonName:name,label,linkage,
        dimension:example.id==='aba-chain'||example.id==='latent-imine'?1:2,
        group:'Assembly examples',search:searchText(name,example.label,feed,linkage,parts.map(names),example.doi)});
    }
    const starters=rows.concat(generated.sort((a,b)=>a.group.localeCompare(b.group)||a.rank-b.rank||a.label.localeCompare(b.label)));
    for(const row of starters){
      const template=templates.get(row.template_id);
      if(template?.parent_net)Object.assign(row,{dimension:template.dimension,net:template.parent_net});
    }
    return starters;
  }

  const commonStarts=['named:COF-LZU1','tppa-pore-v1','crystal:cof-300-c5','assembly:cof-1',
    'design-tfb_deth','assembly:ctf-1','tpb-dmtp-cof-pore-v1','named:4PE-1P','crystal:curated-13040N3','named:Py-Azine'];
  function starterSections(rows,library,query,current){
    const words=(query||'').toLowerCase().trim().split(/\s+/).filter(Boolean);
    const compact=s=>s.toLowerCase().replace(/[\s‐‑–-]/g,'');
    const order=id=>{const i=commonStarts.indexOf(id);return i<0?commonStarts.length:i;};
    const groups=new Map(), byId=new Map(), matches=[];
    for(const row of rows){
      const dimension=row.dimension||(row.crystal?3:row.template_id==='chain_2_2'||row.commonName==='COF-76'?1:2);
      const chemistry=library.linkages?.[row.linkage]?.label?.split(' · ')[0]||row.linkage;
      let section;
      if(row.crystal)section=row.id==='crystal:custom-dia'?'Current pair':'3D · '+row.net+(row.reference?' · reference CIF':'');
      else if(row.source)section='Recorded source drawings';
      else if(dimension===3)section='3D · '+(row.net||chemistry||'Frameworks')+' · molecular circuits';
      else if(dimension===1)section='1D · Chains and ribbons';
      else section='2D · '+(chemistry||'Frameworks');
      const text=[row.id,row.search,row.label,section].filter(Boolean).join(' ');
      if(words.length&&!words.every(w=>text.toLowerCase().includes(w)||compact(text).includes(compact(w))))continue;
      matches.push(row);
      if(!groups.has(section))groups.set(section,[]);
      groups.get(section).push(row);byId.set(row.id,row);
    }
    const common=commonStarts.map(id=>byId.get(id)).filter(Boolean);
    const rank=([k])=>k.startsWith('2D')?0:k.startsWith('3D')?1:k.startsWith('1D')?2:k==='Current pair'?4:3;
    const visible=new Map();
    if(common.length)visible.set('Common COFs',common);
    // A cited reference and a generated model may share a feed. Two generated
    // entries for the same feed and net describe a single choice.
    const modelKey=r=>r.crystal&&r.precursors?.length?JSON.stringify([[...r.precursors].sort(),r.net,r.reference]):r.id;
    const seen=new Set(common.map(modelKey));
    for(const [name,items] of [...groups].sort((a,b)=>rank(a)-rank(b)||a[0].localeCompare(b[0]))){
      const unique=items.filter(row=>{
        const key=modelKey(row);
        if(seen.has(key))return false;
        seen.add(key);return true;
      });
      const top=unique.slice(0,words.length?24:8);
      if(top.length)visible.set(name,top);
    }
    const active=rows.find(r=>r.id===current);
    if(active&&![...visible.values()].some(group=>group.includes(active)))visible.set('Current selection',[active]);
    matches.sort((a,b)=>order(a.id)-order(b.id)||Number(Boolean(b.commonName))-Number(Boolean(a.commonName)));
    return {visible,matches};
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
  const api = {displayName, roleFor, roleEvidence, compatiblePresets, groupRows, reactiveCount, templateCompatibility, stylesInCategory,connectionCount,linkageGroupOptions,partnerLinkageGroups,periodicLinkageMatches,searchMonomers,monomerChoice,compatibleLinkages,automaticChoice,automaticPlan,startingCombinations,starterSections,createPreviewCache};
  scope.COFdrawLibrary = api;
  if (typeof module !== "undefined" && module.exports) module.exports = api;
})(typeof window === "undefined" ? globalThis : window);
