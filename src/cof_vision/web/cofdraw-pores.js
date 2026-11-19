/* Render the periodic probe-center field on the same camera as the atoms. */
"use strict";
(() => {
  const C=window.COFdrawScience,S=window.COFdrawStudio,$=id=>document.getElementById(id);
  let field=null,mesh=null,polygons=[],report=null;
  const bytes=text=>Uint8Array.from(atob(text),ch=>ch.charCodeAt(0));
  function axis(n){
    const count=Math.min(22,n),items=[[n-1,0,0]];
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
    return {vertices,faces,cell:field.cell,values,axes};
  }
  // Fill the open area in the moving section plane, including disconnected islands.
  function cutPlane(cut,xyz){
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
    polygons=[];if(!mesh)return;
    const cut=Number($('pore-cut').value),{vertices,faces,cell}=mesh;
    const xyz=p=>cell[0].map((_,i)=>p[0]*cell[0][i]+p[1]*cell[1][i]+p[2]*cell[2][i]);
    const cart=vertices.map(xyz);
    if(cut>0&&cut<1)cutPlane(cut,xyz);
    for(const face of faces){
      const poly=[];
      for(let n=0;n<3;n++){
        const a=face[n],b=face[(n+1)%3],za=vertices[a][2],zb=vertices[b][2];
        if(za<=cut)poly.push(cart[a]);
        if((za<cut&&zb>cut)||(za>cut&&zb<cut)){
          const t=(cut-za)/(zb-za);poly.push(cart[a].map((v,k)=>v+t*(cart[b][k]-v)));
        }
      }
      if(poly.length>=3)polygons.push(poly);
    }
  }
  function update(){
    const isChannel=$('pore-layer').value==='connected';
    mesh=field?extract(field,isChannel):null;section();
    if(field){
      const text=isChannel?'Certified periodic channels':'All probe-center voids';
      $('map-status').textContent=text+' · '+field.probe_radius+' Å probe · '+field.shape.join(' × ')+' grid. '+
        (mesh.faces.length?'Interpolated display surface in one cell.':'No interface resolved for this selection.');
    }
    C.viewer.draw();
  }
  function reset(){
    field=null;mesh=null;polygons=[];report=null;
    const periodic=C.structure?.periodic===true;
    for(const id of ['map-pores','mof-notebook'])$(id).disabled=!periodic;
    $('builder-cif').disabled=!periodic||![...$('structure-exports').querySelectorAll('button')].some(b=>b.textContent==='CIF');
    $('pore-visible').checked=false;$('pore-visible').disabled=true;
    $('map-status').textContent=periodic?'Calculate the surface for this cell.':'Choose or generate a periodic model to map its pores.';
    C.viewer.draw();
  }
  window.addEventListener('cofdraw:structure',reset);
  window.addEventListener('cofdraw:geometry',e=>{
    if(!e.detail.pore_map||!C.structure?.periodic)return;
    field=e.detail.pore_map;report=e.detail;
    $('pore-visible').disabled=false;$('pore-visible').checked=true;update();
  });
  $('map-pores').onclick=()=>{C.tab('porosity');$('calculate-pores').click();};
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
    get triangleCount(){return mesh?.faces.length||0;},get structureHash(){return report?.structure_sha256;},
    addFaces(objects,project,scale){
      if(!$('pore-visible').checked||!mesh)return;
      const projected=new Map(),channel=$('pore-layer').value==='connected',alpha=Number($('pore-alpha').value);
      for(const polygon of polygons){
        const p=polygon.map(v=>{if(!projected.has(v))projected.set(v,project(v));return projected.get(v);});
        const a=p[0],b=p[1],c=p[2];
        const ux=b.x-a.x,uy=b.y-a.y,uz=(b.z-a.z)*scale,vx=c.x-a.x,vy=c.y-a.y,vz=(c.z-a.z)*scale;
        const nx=uy*vz-uz*vy,ny=uz*vx-ux*vz,nz=ux*vy-uy*vx;
        if(Math.abs(nz)<.02)continue;
        const light=Math.abs((-.35*nx-.52*ny+.78*nz)/(Math.hypot(nx,ny,nz)||1));
        objects.push({kind:'pore',p,z:p.reduce((s,v)=>s+v.z,0)/p.length,alpha,
          color:'hsl('+(channel?32:194)+' 57% '+Math.round(32+31*light)+'%)'});
      }
    },
    paint(ctx,o){
      const p=o.p;
      ctx.fillStyle=o.color;ctx.globalAlpha=o.alpha;
      ctx.beginPath();ctx.moveTo(p[0].x,p[0].y);for(const v of p.slice(1))ctx.lineTo(v.x,v.y);ctx.closePath();ctx.fill();
    }
  };
  reset();
})();
