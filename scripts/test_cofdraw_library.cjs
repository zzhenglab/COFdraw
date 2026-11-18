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
