"use strict";
const assert=require('node:assert/strict');
const L=require('../src/cof_vision/web/cofdraw-library.js');
const aldehyde={compatible_roles:[{linkage:'imine',role:'cho',reactive_sites:3},{linkage:'secondary_amine',role:'cho',reactive_sites:3}]};
const amine={compatible_roles:[{linkage:'imine',role:'nh2',reactive_sites:2},{linkage:'secondary_amine',role:'nh2',reactive_sites:2}]};
const library={linkages:{secondary_amine:{supported:true,roles:['cho','nh2'],label:'Reduced'},imine:{supported:true,roles:['cho','nh2'],label:'Imine'}},
  templates:[{id:'unassigned',allow_pairs:[],topology_label:null},{id:'hcb',allow_pairs:[[3,2],[2,3]],topology_label:'hcb'}]};
const normal=L.automaticPlan(library,aldehyde,amine);
assert.equal(normal.id,'imine');
assert.equal(normal.template_id,'hcb','Auto prefers an explicit site-count template over the general projection');
assert.equal(normal.swapped,false);
assert.equal(L.automaticPlan(library,amine,aldehyde).swapped,true);
assert.equal(L.automaticPlan(library,aldehyde,aldehyde),null);
assert.equal(L.automaticPlan({...library,linkages:{secondary_amine:library.linkages.secondary_amine}},aldehyde,amine),null,'Auto never chooses a product-state reduction');
assert.equal(L.automaticPlan(library,{...aldehyde,geometry_hint:'tetrahedral'},amine).template_id,'unassigned','Auto keeps tetrahedral projection topology unknown');
assert.equal(L.automaticPlan(library,null,amine),null);
console.log('COFdraw Auto selection: 8 checks passed');
