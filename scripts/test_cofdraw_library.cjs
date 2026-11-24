"use strict";
const assert = require("node:assert/strict");
const L = require("../src/cof_vision/web/cofdraw-library.js");
const library = {linkages:{imine:{roles:["aldehyde","amine"]}},monomer_presets:[
  {id:"a",abbreviation:"A4",full_name:"Four-site aldehyde",compatible_roles:[{linkage:"imine",role:"aldehyde",reactive_sites:4}]},
  {id:"b",compatible_roles:[{linkage:"imine",role:"amine",reactive_sites:2}]},
  {id:"unsupported",annotation_only:true,compatible_roles:[]}
]};
assert.equal(L.displayName(library.monomer_presets[0]),"A4 — Four-site aldehyde");
assert.deepEqual(L.compatiblePresets(library,"imine","a").map(x=>x.id),["a"]);
assert.deepEqual(L.compatiblePresets(library,"imine","b").map(x=>x.id),["b"]);
assert.equal(L.roleEvidence({compatible_roles:[]},"imine","aldehyde").compatible,false);
assert.equal(L.roleEvidence({smiles:"not chemistry"},"imine","aldehyde"),null);
assert.equal(L.roleEvidence({compatible_roles:[{linkage:"imine",role:"aldehyde",reactive_sites:1}]},"imine","aldehyde").compatible,false);
const monomers = {a:library.monomer_presets[0],b:library.monomer_presets[1]};
const automaticLibrary={...library,linkages:{imine:{supported:true,roles:['aldehyde','amine']},reference:{supported:false,roles:[]}},templates:[{id:'sql',allow_pairs:[[4,2]]}]};
assert.deepEqual(L.compatibleLinkages(automaticLibrary,monomers.a,monomers.b).map(x=>[x.id,x.swapped,x.templates]),[['imine',false,['sql']]]);
assert.equal(L.compatibleLinkages(automaticLibrary,monomers.b,monomers.a)[0].swapped,true);
assert.deepEqual(L.compatibleLinkages(automaticLibrary,monomers.a,{compatible_roles:[]}),[]);
assert.equal(L.automaticChoice(automaticLibrary,monomers.a,'a',monomers.b).compatible,true);
assert.equal(L.automaticChoice(automaticLibrary,{compatible_roles:[]},'a',monomers.b).compatible,false);
assert.equal(L.automaticChoice(automaticLibrary,monomers.a,'b',monomers.a).compatible,false,"An individually supported aldehyde is not an aldehyde partner");
assert.equal(L.automaticChoice(automaticLibrary,monomers.b,'b',monomers.b).compatible,false,"Two amines do not become compatible through the fallback");
assert.equal(L.automaticChoice(automaticLibrary,monomers.a,'a',monomers.a).compatible,true,"A can be replaced without a partner-selection deadlock");
assert.equal(L.automaticChoice(automaticLibrary,monomers.a,'b',null).compatible,true,"A missing first choice keeps single-role browsing available");
assert.equal(L.automaticChoice(automaticLibrary,monomers.a,'b',monomers.b).compatible,true,"Automatic chemistry accepts reversed roles");
assert.equal(L.automaticChoice({...automaticLibrary,templates:[{allow_pairs:[[3,3]]}]},monomers.b,'b',monomers.a).compatible,false,"A matching role is insufficient when every explicit template rejects the pair");
assert.equal(L.automaticChoice(automaticLibrary,{compatible_roles:[{linkage:'reference',role:'x',reactive_sites:3}]},'a',null).compatible,false,"Reference-only chemistry cannot authorize construction");
const nitroneLibrary={linkages:{nitrone:{supported:true,roles:['aryl_aldehyde','aryl_hydroxylamine']}},templates:[{id:'chain_2_2',allow_pairs:[[2,2]]}]};
const nitroneA={compatible_roles:[{linkage:'nitrone',role:'aryl_aldehyde',reactive_sites:2}]};
const nitroneB={compatible_roles:[{linkage:'nitrone',role:'aryl_hydroxylamine',reactive_sites:2}]};
assert.equal(L.automaticPlan(nitroneLibrary,nitroneA,nitroneB).id,'nitrone','The supported hydroxylamine pair has a real automatic construction');
assert.equal(L.automaticPlan(nitroneLibrary,nitroneB,nitroneA).swapped,true,'Visible A/B order may remain stable while payload roles are mapped');
assert.equal(L.automaticChoice(nitroneLibrary,monomers.b,'b',nitroneA).compatible,false,'An ordinary amine is not a hydroxylamine partner');
const planarA = {...monomers.a, attachment_directions:{aldehyde:"square_planar"}};
const tetraB = {compatible_roles:[{linkage:"imine", role:"amine", reactive_sites:4}], attachment_directions:{amine:"tetrahedral"}};
const planarB = {...tetraB, attachment_directions:{amine:"square_planar"}};
const spatialLibrary = {...automaticLibrary, templates:[
  {id:"sql_4_4", topology_label:"sql", allow_pairs:[[4,4]]},
  {id:"pts_4_4", topology_label:"pts", allow_pairs:[[4,4]]}
]};
for (const [first, second] of [[planarA,tetraB], [tetraB,planarA]]) {
  assert.deepEqual(L.compatibleLinkages(spatialLibrary, first, second)[0].templates, ["pts_4_4"]);
}
assert.deepEqual(L.compatibleLinkages(spatialLibrary, planarA, planarB)[0].templates, ["sql_4_4"]);
assert.equal(L.templateCompatibility(spatialLibrary.templates[1], "imine", {a:planarA, b:monomers.b}, library).compatible, false);
const square = {id:"square",linkages:["imine"],required_functionalities:{a:[4],b:[2]}};
assert.equal(L.templateCompatibility(square,"imine",monomers,library).compatible,true);
assert.equal(L.templateCompatibility({...square,required_functionalities:{a:[3],b:[2]}},"imine",monomers,library).compatible,false);
assert.equal(L.templateCompatibility(square,"boronate",monomers,library).compatible,false);
assert.equal(L.templateCompatibility(square,"imine",{a:{compatible_roles:[]}},library).compatible,false);
assert.equal(L.templateCompatibility(square,"imine",{},library).compatible,true); // unknown is checked by engine, not guessed
assert.equal(L.templateCompatibility({allow_pairs:[[3,2],[2,3]]},"imine",monomers,library).compatible,false);
assert.equal(L.templateCompatibility({allow_pairs:[[2,4],[4,2]]},"imine",monomers,library).compatible,true);
assert.equal(L.templateCompatibility({allow_pairs:[]},"imine",monomers,library).compatible,true);
assert.equal(L.templateCompatibility({topology_label:"sql",allow_pairs:[[4,2]]},"imine",{...monomers,a:{...monomers.a,geometry_hint:"tetrahedral"}},library).compatible,false);
assert.deepEqual(L.stylesInCategory([{id:"plain"},{id:"stack",category:"layered"}],"chemical").map(x=>x.id),["plain"]);
const cross={id:"cross",abbreviation:"Cross",full_name:"Four-arm cross",point_group:"C2",connection_count:4,literature_names:["Alternate name"]};
assert.equal(L.connectionCount(cross),4);
const aldehydeGroup = {id:"aldehyde", label:"Aldehyde"};
const amineGroup = {id:"primary_amine", label:"Primary amine"};
const decorated = {...cross, functional_groups:["Aldehyde","Phenol"], linkage_functional_groups:[aldehydeGroup]};
const amine = {id:"diamine", connection_count:2, linkage_functional_groups:[amineGroup]};
const noLinkageRole = {id:"reference", functional_groups:["Phenol"], linkage_functional_groups:[]};
const catalog = [decorated, cross, amine, noLinkageRole];
assert.deepEqual(L.linkageGroupOptions([amine,...catalog]),[aldehydeGroup,amineGroup]);
assert.deepEqual(L.searchMonomers(catalog,"","all","aldehyde"),[decorated]);
assert.deepEqual(L.searchMonomers(catalog,"",4,"all"),[decorated,cross]);
assert.deepEqual(L.searchMonomers(catalog,"alternate",4,"aldehyde"),[decorated]);
assert.deepEqual(L.searchMonomers(catalog,"alternate",3,"aldehyde"),[]);
assert.deepEqual(L.searchMonomers(catalog,"",4,"primary_amine"),[]);
assert.deepEqual(L.searchMonomers(catalog,"","all","phenol"),[]);
// Only the group's role in a viable pair qualifies in the compatible-B view.
const carbonyl = {linkage:"imine",role:"aldehyde"}, amino = {linkage:"imine",role:"amine"};
const boronic = {linkage:"boronate",role:"boronic_acid"};
const groupRows = [
  {...aldehydeGroup,binary_roles:[carbonyl]}, {...amineGroup,binary_roles:[amino]},
  {id:"boronic_acid",label:"Boronic acid",binary_roles:[boronic]},
  {id:"nitrile",label:"Nitrile",binary_roles:[]}
];
const multifunctional = {linkage_functional_groups:groupRows,
  binary_roles:[carbonyl,amino,boronic].map(r=>({...r,reactive_sites:2}))};
const partnerLibrary = {linkages:{
  imine:{supported:true,roles:["aldehyde","amine"]}, boronate:{supported:true,roles:["diol","boronic_acid"]}
}};
const feed = {binary_roles:[{...carbonyl,reactive_sites:3},{linkage:"boronate",role:"diol",reactive_sites:2}]};
const partnerGroupIds = (a,selected,automatic,lib=partnerLibrary) =>
  L.partnerLinkageGroups(lib,a,multifunctional,selected,automatic).map(g=>g.id);
assert.deepEqual(partnerGroupIds(feed,"imine",true),["primary_amine","boronic_acid"]);
assert.deepEqual(partnerGroupIds(feed,"imine",false),["primary_amine"]);
assert.deepEqual(partnerGroupIds({binary_roles:[{...amino,reactive_sites:3}]},"imine",false),["aldehyde"]);
const dualFeed = {binary_roles:multifunctional.binary_roles};
assert.equal(L.compatibleLinkages(partnerLibrary,dualFeed,multifunctional).find(m=>m.id==="imine").swapped,false);
assert.deepEqual(partnerGroupIds(dualFeed,"imine",false),["primary_amine"],
  "With two viable imine assignments, filter B by the end used by the fixed-linkage drawing");
assert.deepEqual(partnerGroupIds({...library.monomer_presets[0],geometry_hint:"tetrahedral"},"imine",false,
  {...partnerLibrary,templates:[{allow_pairs:[[3,3]]}]}),["primary_amine"]);
assert.deepEqual(L.searchMonomers([cross],"alternate",4),[cross]);
assert.deepEqual(L.searchMonomers([cross],"",2),[]); // C2 symmetry is not two reactive sites.
assert.equal(L.monomerChoice(library,library.monomer_presets[0],"imine","b").compatible,false);
assert.match(L.monomerChoice(library,library.monomer_presets[0],"imine","b").reason,/Requires amine/);
assert.equal(L.monomerChoice(library,library.monomer_presets[0],"imine","editor").compatible,true);
assert.equal(L.monomerChoice({...library,templates:[square]},library.monomer_presets[0],"imine","a",{compatible_roles:[]}).compatible,true); // replacing one side must not deadlock after a linkage switch
assert.equal(L.monomerChoice(library,monomers.b,"imine","b",monomers.b).compatible,false,"Fixed chemistry requires A to supply its specified role");
assert.equal(L.monomerChoice(library,{},"imine","b",monomers.a).compatible,false,"Unknown metadata cannot authorize a partner");
async function checkPreviewQueue(){
  const releases=new Map(),calls=[],cache=L.createPreviewCache(key=>{
    calls.push(key);return new Promise(resolve=>releases.set(key,resolve));
  },2);
  const a=cache.request("graph-a"),duplicate=cache.request("graph-a"),b=cache.request("graph-b"),c=cache.request("graph-c");
  assert.strictEqual(a,duplicate,"Repeated visible cards must share in-flight rendering");
  await new Promise(setImmediate);
  assert.deepEqual(calls,["graph-a","graph-b"],"Requests beyond the concurrency limit stay queued");
  releases.get("graph-a")("svg-a");await a;await new Promise(setImmediate);
  assert.deepEqual(calls,["graph-a","graph-b","graph-c"]);
  releases.get("graph-b")("svg-b");releases.get("graph-c")("svg-c");
  assert.deepEqual(await Promise.all([a,b,c]),["svg-a","svg-b","svg-c"]);
  assert.equal(await cache.request("graph-a"),"svg-a");assert.equal(calls.length,3,"Changing filters reuses completed previews");
  let attempts=0;
  const retry=L.createPreviewCache(async()=>{if(++attempts===1)throw new Error("temporary failure");return "recovered";});
  await assert.rejects(retry.request("graph"),/temporary failure/);
  assert.equal(await retry.request("graph"),"recovered","A failed preview can be retried from details");
  console.log("COFdraw metadata filters and bounded preview cache: passed");
}
checkPreviewQueue().catch(error=>{console.error(error);process.exitCode=1;});
