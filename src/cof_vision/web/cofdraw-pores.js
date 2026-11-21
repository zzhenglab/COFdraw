/* Render the periodic probe-center field on the same camera as the atoms. */
"use strict";
(() => {
  const C=window.COFdrawScience,S=window.COFdrawStudio,$=id=>document.getElementById(id);
  let field=null,mesh=null,tiles=[],report=null;
  function repetitions(){
    return field.cell.map((row,i)=>Math.max(1,Math.round(Math.hypot(...C.viewer.model.cell[i])/Math.hypot(...row))));
  }
  const bytes=text=>Uint8Array.from(atob(text),ch=>ch.charCodeAt(0));
  function axis(n){
    const cells=repetitions().reduce((a,b)=>a*b,1);
    const count=Math.min(22,n,Math.max(6,Math.floor(44/Math.cbrt(cells)))),items=[[n-1,0,0]];
    for(let i=0;i<count;i++){
      const j=Math.round(i*(n-1)/(count-1));items.push([j,j,(j+.5)/n]);
    }
    items.push([n-1,0,1]);return items;
  }
  function extract(field,connected){
    const [nx,ny,nz]=field.shape,axes=[nx,ny,nz].map(axis);
    const samples=new DataView(bytes(field.clearance_f32le).buffer),bits=bytes(field.connected_bits);
    function scalar(i,j,k){
      const idx=(i*ny+j)*nz+k, v=samples.getFloat32(idx*4,true)-field.probe_radius;
      return connected && !(bits[idx>>3]&(1<<(idx%8))) ? Math.min(-.001,-Math.abs(v)) : v;
    }
    const values=[],points=[],vertices=[],faces=[],edgeCache=new Map();
    // Extend the cell-centered field to the six cell faces by periodic interpolation.
    for(const x of axes[0])for(const y of axes[1])for(const z of axes[2]){
      let v=0;
      for(const i of x.slice(0,2))for(const j of y.slice(0,2))for(const k of z.slice(0,2))v+=scalar(i,j,k);
      values.push(v/8);points.push([x[2],y[2],z[2]]);
    }
    function cross(a,b){
      if(a>b)[a,b]=[b,a];
      const key=a*points.length+b;
      if(edgeCache.has(key))return edgeCache.get(key);
      const t=values[a]/(values[a]-values[b]),p=points[a].map((v,k)=>v+(points[b][k]-v)*t);
      const idx=vertices.push(p)-1;edgeCache.set(key,idx);return idx;
    }
    const dy=axes[2].length,dx=dy*axes[1].length;
    const corners=[0,dx,dx+dy,dy,1,dx+1,dx+dy+1,dy+1];
    // All cubes use the same body diagonal, including adjacent periodic faces.
    const tetra=[[0,1,2,6],[0,2,3,6],[0,3,7,6],[0,7,4,6],[0,4,5,6],[0,5,1,6]];
    for(let i=0;i<axes[0].length-1;i++)for(let j=0;j<axes[1].length-1;j++)for(let k=0;k<dy-1;k++){
      const base=i*dx+j*dy+k, ids=corners.map(o=>o+base);
      if(ids.every(n=>values[n]>=0)||ids.every(n=>values[n]<0))continue;
      for(const tet of tetra){
        const p=[],m=[];
        for(const k of tet){const n=ids[k];(values[n]>=0?p:m).push(n);}
        if(!p.length||!m.length)continue;
        if(p.length===2){
          const a=cross(p[0],m[0]),b=cross(p[0],m[1]),c=cross(p[1],m[0]),d=cross(p[1],m[1]);
          faces.push([a,b,c],[b,d,c]);
        }else{
          const [one,many]=p.length===1?[p,m]:[m,p];
          faces.push(many.map(v=>cross(one[0],v)));
        }
      }
    }
    // The viewer uses CIF's conventional Cartesian orientation. The sampled
    // field is fractional and can come from a differently oriented input cell.
    const repeats=repetitions();
    const cell=C.viewer.model.cell.map((row,i)=>row.map(v=>v/repeats[i]));
    return {vertices,faces,cell,values,axes};
  }
  // Fill the open area in the moving section plane, including disconnected islands.
  function cutPlane(cut,xyz,polygons){
    const {axes,values}=mesh,[x,y,z]=axes;
    const k=z.findIndex(sample=>sample[2]>=cut)-1;
    if(k<0)return;
    const t=(cut-z[k][2])/(z[k+1][2]-z[k][2]);
    const nodes=x.map((a,i)=>y.map((b,j)=>{
      const n=(i*y.length+j)*z.length+k;
      return {p:[a[2],b[2],cut],v:values[n]+t*(values[n+1]-values[n])};
    }));
    for(let i=1;i<x.length;i++)for(let j=1;j<y.length;j++){
      const a=nodes[i-1][j-1],b=nodes[i][j-1],c=nodes[i][j],d=nodes[i-1][j];
      for(const tri of [[a,b,c],[a,c,d]]){
        const poly=[];
        for(let e=0;e<3;e++){
          const u=tri[e],v=tri[(e+1)%3];
          if(u.v>=0)poly.push(xyz(u.p));
          if((u.v>0&&v.v<0)||(u.v<0&&v.v>0)){
            const f=u.v/(u.v-v.v);
            poly.push(xyz(u.p.map((n,axis)=>n+f*(v.p[axis]-n))));
          }
        }
        if(poly.length>=3)polygons.push(poly);
      }
    }
  }
  function section(){
    tiles=[];if(!mesh)return;
    const cut=Number($('pore-cut').value),{vertices,faces,cell}=mesh;
    const xyz=p=>cell[0].map((_,i)=>p[0]*cell[0][i]+p[1]*cell[1][i]+p[2]*cell[2][i]);
    const cart=vertices.map(xyz);
    const [na,nb,nc]=repetitions(),height=cut*nc;
    for(let k=0;k<Math.ceil(height);k++){
      const top=Math.min(1,height-k),polygons=[];
      if(k===Math.ceil(height)-1)cutPlane(top,xyz,polygons);
      for(const face of faces){
        const poly=[];
        for(let n=0;n<3;n++){
          const a=face[n],b=face[(n+1)%3],za=vertices[a][2],zb=vertices[b][2];
          if(za<=top)poly.push(cart[a]);
          if((za<top&&zb>top)||(za>top&&zb<top)){
            const t=(top-za)/(zb-za);poly.push(cart[a].map((v,k)=>v+t*(cart[b][k]-v)));
          }
        }
        if(poly.length>=3)polygons.push(poly);
      }
      for(let i=0;i<na;i++)for(let j=0;j<nb;j++)tiles.push({polygons,offset:xyz([i,j,k])});
    }
  }
  function update(){
    if(C.viewer.model?.illustration){
      mesh=null;tiles=[];
      $('map-status').textContent='Choose Crystal packing to view the pore surface of the full crystal.';
      C.viewer.draw();return;
    }
    const isChannel=$('pore-layer').value==='connected';
    mesh=field?extract(field,isChannel):null;section();
    if(field){
      const text=isChannel?'Certified periodic channels':'All probe-center voids';
      $('map-status').textContent=text+' · '+field.probe_radius+' Å probe · '+field.shape.join(' × ')+' grid. '+
        (mesh.faces.length?'Mapped across the '+repetitions().join(' × ')+' displayed cells.':'No interface resolved for this selection.');
    }
    C.viewer.draw();
  }
  let lastStructure, lastSettings='', pending=null;
  function settings(){return [$('probe-radius').value,$('pore-samples').value,$('pore-refinement').checked].join('|');}
  async function calculate(){
    const structure=C.structure,key=settings();
    if(!structure?.periodic)return null;
    if(report&&structure===lastStructure&&key===lastSettings)return report;
    if(!pending){
      const spacing=Number($('pore-samples').value);
      $('map-status').textContent='Calculating pore surface…';
      pending=S.api('science/geometry',{structure,probe_radius:Number($('probe-radius').value),grid_spacing:spacing,
        surface_samples:spacing>=.9?32:spacing>=.6?64:128,refinement:$('pore-refinement').checked,include_map:true});
    }
    const request=pending;
    try{
      const result=await request;
      if(structure!==C.structure||key!==settings()||request!==pending&&result!==report)return null;
      if(result!==report){field=result.pore_map;report=result;$('pore-visible').disabled=false;$('pore-visible').checked=true;update();}
      return result;
    }catch(e){
      if(request!==pending)return null;
      $('map-status').textContent='Pore map unavailable. '+e.message;throw e;
    }finally{if(request===pending)pending=null;}
  }
  function reset(event){
    const changed=C.structure!==lastStructure;
    if(!changed&&settings()===lastSettings)return;
    lastStructure=C.structure;lastSettings=settings();pending=null;
    field=null;mesh=null;tiles=[];report=null;
    const periodic=C.structure?.periodic===true;
    for(const id of ['map-pores','mof-notebook'])$(id).disabled=!periodic;
    $('builder-cif').disabled=!periodic||![...$('structure-exports').querySelectorAll('button')].some(b=>b.textContent==='CIF');
    $('pore-visible').checked=false;$('pore-visible').disabled=true;
    $('map-status').textContent=periodic?'Calculate the surface for this cell.':'Choose or generate a periodic model to map its pores.';
    C.viewer.draw();
    if(periodic&&changed)calculate().catch(()=>{});
  }
  window.addEventListener('cofdraw:structure',reset);
  window.addEventListener('cofdraw:display-cells',update);
  window.addEventListener('cofdraw:crystal-scene',update);
  $('map-pores').onclick=()=>{C.tab('porosity');$('calculate-pores').click();};
  document.querySelector('[data-science-tab="porosity"]').addEventListener('click',()=>{
    if(C.structure?.periodic&&!$('porosity-results').childElementCount)$('calculate-pores').click();
  });
  $('render-toggle').onclick=()=>{
    const panel=$('render-panel');panel.hidden=!panel.hidden;$('render-toggle').setAttribute('aria-expanded',String(!panel.hidden));
  };
  $('builder-cif').onclick=()=>[...$('structure-exports').querySelectorAll('button')].find(b=>b.textContent==='CIF')?.click();
  $('mof-notebook').onclick=async()=>{
    try{
      const s=C.structure;if(!s?.periodic)return;
      const input={symbols:s.atoms.map(a=>a.element),cell:s.cell,scaled_positions:s.atoms.map(a=>a.fractional)};
      const probe=Number($('probe-radius').value),spacing=Number($('pore-samples').value);
      if(!Number.isFinite(probe)||probe<0||probe>5)throw new Error('Enter a probe radius from 0 to 5 Å.');
      const code=source=>({cell_type:'code',metadata:{},execution_count:null,outputs:[],source});
      const cells=[{cell_type:'markdown',metadata:{},source:'# Mofography pore analysis\nActive COFdraw cell. Run in a Python environment with [Mofography 2026.10.09.1](https://github.com/mofography/mofography/tree/2026.10.09.1) and its compiled Rust backend. Structural hypotheses retain their original qualification. Mofography uses its own radii and segmentation defaults.\n\nModel: '+$('viewer-title').textContent+' · '+$('viewer-provenance').textContent},
        code('import json\nimport mofography as mgr\nfrom ase import Atoms\nimport numpy as np\nmodel = json.loads('+JSON.stringify(JSON.stringify(input))+')\natoms = Atoms(**model, pbc=True)\nspacing = '+spacing+'\nprobe = '+probe),
        code('dgrid = mgr.dgrid_from_atoms_rust(atoms, spacing=spacing).astype(np.float32)\nseg = mgr.do_periodic_segmentation_rust(dgrid, atoms, spacing=spacing, compute_connections=True, return_dgrid=True, disable_tqdm=True)\nprint("Pockets:", len(seg["maxima_pbc"]))'),
        code('windows = [row for row in seg["connections"] if np.any(np.isin(row[0], seg["unique_regions"]))]\nprint("Windows in cell:", len(windows))\nif windows:\n    radii = [float(np.ravel(row[4])[0]) for row in windows]\n    print("Window radius range (Å):", min(radii), max(radii))'),
        code('lcd, pld, dfsd, _ = mgr.get_pore_dia_mgr_winding(atoms, dgrid=dgrid, results=seg)\nprint("LCD, PLD, DFSD (Å):", lcd, pld, dfsd)\nprint(mgr.compute_void_volume_from_results_winding(seg, atoms, probe_radius=probe))\n# ASE uses natural element masses unless isotope masses are specified.'),
        code('typed = None\nif len(seg["maxima_pbc"]):\n    seg = mgr.enrich_segmentation_with_pocket_psds(seg, atoms)\n    typed = mgr.assign_pocket_types(seg, atoms, method="env_cluster")\nfig = mgr.show_mof_plotly(seg, atoms, typed_result=typed, framework_atoms=atoms, stylize=True)\nfig.show()')];
      await S.saveText('COFdraw-Mofography.ipynb',JSON.stringify({nbformat:4,nbformat_minor:4,metadata:{kernelspec:{name:'python3',display_name:'Python 3',language:'python'}},cells},null,2));
    }catch(e){C.status(e.message,true);}
  };
  $('pore-layer').onchange=update;
  $('pore-cut').oninput=()=>{section();C.viewer.draw();};
  for(const id of ['pore-alpha','pore-visible'])$(id).oninput=()=>C.viewer.draw();
  window.COFdrawPores={
    calculate,
    get triangleCount(){return mesh?.faces.length||0;},get tileCount(){return tiles.length;},get structureHash(){return report?.structure_sha256;},
    addFaces(objects,project,scale){
      if(!$('pore-visible').checked||!mesh)return;
      const channel=$('pore-layer').value==='connected',alpha=Number($('pore-alpha').value);
      for(const {polygons,offset} of tiles){
        const projected=new Map();
        for(const polygon of polygons){
          const p=polygon.map(v=>{if(!projected.has(v))projected.set(v,project(v.map((x,i)=>x+offset[i])));return projected.get(v);});
          const a=p[0],b=p[1],c=p[2];
          const ux=b.x-a.x,uy=b.y-a.y,uz=(b.z-a.z)*scale,vx=c.x-a.x,vy=c.y-a.y,vz=(c.z-a.z)*scale;
          const nx=uy*vz-uz*vy,ny=uz*vx-ux*vz,nz=ux*vy-uy*vx;
          if(Math.abs(nz)<.02)continue;
          const light=Math.abs((-.35*nx-.52*ny+.78*nz)/(Math.hypot(nx,ny,nz)||1));
          objects.push({kind:'pore',p,z:p.reduce((s,v)=>s+v.z,0)/p.length,alpha,
            color:'hsl('+(channel?32:194)+' 57% '+Math.round(32+31*light)+'%)'});
        }
      }
    },
    paint(ctx,o,visibility=1){
      const p=o.p;
      ctx.fillStyle=o.color;ctx.globalAlpha=o.alpha*visibility;
      ctx.beginPath();ctx.moveTo(p[0].x,p[0].y);for(const v of p.slice(1))ctx.lineTo(v.x,v.y);ctx.closePath();ctx.fill();
    }
  };
  reset();
})();
