"""Trace finite pore faces through the heavy-atom graph of a reference CIF."""
from collections import deque

import numpy as np


def bridge_overlap_cost(points,links,cameras):
    """Score nonincident projected bridges, particularly parallel crossings."""
    costs=np.zeros(len(cameras))
    lengths=np.array([np.linalg.norm(points[b]-points[a]) for a,b in links])
    scale=np.median(lengths)
    for k,(a,b) in enumerate(links):
        d=points[b]-points[a]
        dz=cameras@d
        ld=np.maximum(d@d-dz**2,1.e-9)
        for i,j in links[:k]:
            if a in (i,j) or b in (i,j):continue
            e=points[j]-points[i]
            ez=cameras@e
            le=np.maximum(e@e-ez**2,1.e-9)
            cosine=(e@d-ez*dz)/np.sqrt(ld*le)
            # Prefer a view that separates complete bridges before flattening
            # their aryl chemistry. Oblique chairs can often keep all bridges
            # distinct even when a slightly more face-on view crosses them.
            baseline=points[i]-points[a]
            signed=cameras@np.cross(d,e)
            safe=np.where(np.abs(signed)>1.e-8,signed,np.inf)
            along=cameras@np.cross(baseline,e)/safe
            across=cameras@np.cross(baseline,d)/safe
            crossing=(along>.01)&(along<.99)&(across>.01)&(across<.99)
            costs+=2.5*crossing
            # A transverse crossover can be depth-faded. Two aryls running
            # along the same projected edge can appear to share bonds.
            proximity=np.zeros(len(cameras))
            for fraction in (.2,.5,.8):
                p=points[a]-points[i]+fraction*d
                pz=cameras@p
                along=p@e-pz*ez
                t=np.clip(along/le,0,1)
                dist=np.sqrt(np.maximum(0,p@p-pz**2-2*t*along+t**2*le))/scale
                proximity+=np.maximum(.42-dist,0)**2
            costs+=proximity*(.3+1.7*np.clip(cosine**2,0,1))
    return costs


def reference_framework(structure, edges, limit=4):
    """Contract local molecular rings, then follow bonds between branching residues.

    Every vertex includes its lattice image. An edge touching the opposite face
    of the cell cannot close a molecular or framework ring by wrapping alone.
    """
    atoms=structure["atoms"]
    graph=[[] for _ in atoms]
    for e in edges:
        i,j,t=e["a"],e["b"],tuple(e["translation"])
        if "H" in (atoms[i]["element"],atoms[j]["element"]):continue
        graph[i].append((j,*t));graph[j].append((i,*(-n for n in t)))
    def translated(v,origin):
        return (v[0],v[1]+origin[1],v[2]+origin[2],v[3]+origin[3])
    def adjacent(v,table):
        return [translated(w,v) for w in table[v[0]]]
    # Merge local ring systems, including the internal macrocycle of a
    # porphyrin, before tracing cycles of the framework itself.
    local=[[] for _ in atoms]
    for i,rows in enumerate(graph):
        start=(i,0,0,0)
        for target in rows:
            j,x,y,z=target
            if start>=target:continue
            seen={start};front=set(rows)-{target}
            for _ in range(16):
                seen.update(front)
                front={q for p in front for q in adjacent(p,graph) if q not in seen}
                if target in front:
                    local[i].append(target);local[j].append((i,-x,-y,-z));break
                if not front:break
    groups=[];owner={};offset={}
    for i in range(len(atoms)):
        if i in owner or not graph[i]:continue
        group=len(groups);queue=[(i,0,0,0)];members=set()
        while queue:
            site=queue.pop();n=site[0]
            if n in owner:
                if offset[n]!=site[1:]:return {}, None
                continue
            owner[n]=group;offset[n]=site[1:];members.add(site)
            queue.extend(adjacent(site,local))
        groups.append(members)
    residue=[set() for _ in groups]
    for i,g in owner.items():
        for j,x,y,z in graph[i]:
            h=owner[j]
            shift=tuple(a+b-c for a,b,c in zip(offset[i],(x,y,z),offset[j]))
            if h!=g or any(shift):residue[g].add((h,*shift))
    # Pendant substituents do not create framework junctions.
    keep=set(range(len(groups)))
    while True:
        ends={i for i in keep if sum(v[0] in keep for v in residue[i])<2}
        if not ends:break
        keep-=ends
    net=[sorted(v for v in row if v[0] in keep) for row in residue]
    junctions={i for i in keep if len(net[i])>=3}
    if not junctions:return {}, None
    arms={i:{} for i in junctions}
    for i in junctions:
        root=(i,0,0,0)
        for first in net[i]:
            path=[];previous,site=root,first
            while site[0] not in junctions and len(path)<100:
                upcoming=set(adjacent(site,net))-{previous}
                if len(upcoming)!=1:break
                path.append(site);previous,site=site,upcoming.pop()
            if site[0] not in junctions or site==root or site in arms[i]:return {}, None
            arms[i][site]=path
    def neighbors(v):
        return [translated(w,v) for w in arms[v[0]]]
    def other_face(a,b,excluded,max_length):
        # Shortest closed path through this real edge using another arm of b.
        queue=deque([a,b,n] for n in neighbors(b) if n!=a and n not in excluded)
        seen={a,b}|set(excluded)
        while queue and len(seen)<3000:
            path=queue.popleft();v=path[-1]
            if v in seen:continue
            seen.add(v);near=neighbors(v)
            if a in near:return path
            if len(path)<max_length:queue.extend(path+[n] for n in near if n not in seen)
        return None
    start=(min(junctions),0,0,0)
    first=next((p for b in sorted(neighbors(start)) if (p:=other_face(start,b,set(),10))),None)
    if not first:return {}, None
    cell=np.asarray(structure["cell"])
    centers=[np.mean([np.asarray(atoms[i]["fractional"])+[x,y,z] for i,x,y,z in group],axis=0) for group in groups]
    # Wrap each centroid into the primitive cell. Its incident edges must carry
    # the same shift, even when the contracted ring crosses a cell boundary.
    ordered=sorted(junctions)
    indices={g:i for i,g in enumerate(ordered)}
    images={g:np.floor(centers[g]+1e-10).astype(int) for g in ordered}
    spans=set()
    for g in ordered:
        for h,x,y,z in arms[g]:
            shift=np.array([x,y,z])+images[h]-images[g]
            forward=(indices[g],indices[h],*map(int,shift))
            reverse=(indices[h],indices[g],*map(int,-shift))
            spans.add(min(forward,reverse))
    anchors=[min(groups[g],key=lambda site:np.linalg.norm(
        (np.array(atoms[site[0]]['fractional'])+site[1:]-centers[g])@cell)) for g in ordered]
    from .cofdraw_network import _components
    contracted=[dict(a=a,b=b,translation=t) for a,b,*t in spans]
    components,_=_components(len(ordered),contracted,cell)
    if any(group['periodic_rank']!=3 for group in components):
        return {}, None
    topology=dict(kind='Molecular junctions from distance-inferred CIF bonds',
                  points=[np.maximum(0.,centers[g]-images[g]).tolist() for g in ordered],
                  site_ids=[site[0] for site in anchors],
                  site_shifts=[(np.array(site[1:])-images[g]).tolist() for g,site in zip(ordered,anchors)],
                  edges=[list(edge) for edge in sorted(spans)])
    def centroid(loop):
        return np.mean([(centers[i]+[x,y,z])@cell for i,x,y,z in loop],axis=0)
    def ring_edges(loop):
        return {tuple(sorted((a,b))) for a,b in zip(loop,loop[1:]+loop[:1])}
    def sites_for(loop):
        cluster=set(loop)
        for a,b in zip(loop,loop[1:]+loop[:1]):
            target=(b[0],*(v-u for u,v in zip(a[1:],b[1:])))
            cluster.update(translated(p,a) for p in arms[a[0]][target])
        # A complete tetraaryl node includes the two aryl rings that lead out
        # of a circuit. Carry each adjacent cyclic residue to the cut so that
        # continuation marks appear beyond the ring, not on a bare Si or C.
        for g in loop:
            if len(groups[g[0]])!=1 or len(net[g[0]])!=4:
                continue
            atom=next(iter(groups[g[0]]))[0]
            if atoms[atom]["element"] in ("C","Si") and all(atoms[j]["element"]=="C" for j,*_ in graph[atom]):
                cluster.update(s for s in adjacent(g,residue) if len(groups[s[0]])>=6)
        # Complete substituents pendant to these residues, e.g. fluorine or
        # methoxy, without crossing into another ring of the extended net.
        pending=list(cluster)
        while pending:
            newly_found={s for s in adjacent(pending.pop(),residue) if s[0] not in keep}-cluster
            cluster.update(newly_found)
            pending.extend(newly_found)
        return {translated(p,g) for g in cluster for p in groups[g[0]]}
    points=np.array([(centers[i]+[x,y,z])@cell for i,x,y,z in first])
    direction=np.linalg.svd(points-points.mean(0),full_matrices=False)[2][-1].tolist()
    patches={};loops=[first];used=ring_edges(first);site_images=sites_for(first)
    rejected=set()
    while len(loops)<=limit:
        n=len(loops)
        if n>1:
            faces=[np.array([(centers[i]+[x,y,z])@cell for i,x,y,z in loop]) for loop in loops]
            normals=np.array([np.linalg.svd(p-p.mean(0))[2][-1] for p in faces])
            edges_xyz=np.concatenate([np.roll(p,-1,axis=0)-p for p in faces])
            unit=edges_xyz/np.maximum(np.linalg.norm(edges_xyz,axis=1)[:,None],1.e-8)
            k=np.arange(256)+.5
            z=1-2*k/256
            r=np.sqrt(1-z*z)
            cameras=np.column_stack((r*np.cos(k*2.399963),r*np.sin(k*2.399963),z))
            exposure=np.abs(normals@cameras.T)
            shortening=np.sqrt(np.maximum(0,1-(unit@cameras.T)**2))
            vertices=sorted({p for face in loops for p in face})
            coords=np.array([(centers[i]+[x,y,z])@cell for i,x,y,z in vertices])
            links=[(vertices.index(a),vertices.index(b)) for a,b in sorted(used)]
            overlap=bridge_overlap_cost(coords,links,cameras)
            direction=cameras[np.argmax(exposure.min(0)+.65*shortening.min(0)+.18*exposure.mean(0)-2*overlap)].tolist()
        patches[n]=dict(circuits=n,site_images=[[i,[x,y,z]] for i,x,y,z in sorted(site_images)],direction=direction)
        patches[n]["junction_loops"]=[[(centers[i]+[x,y,z]).tolist() for i,x,y,z in loop] for loop in loops]
        if n==limit:break
        options={}
        for loop in loops:
            for a,b in zip(loop,loop[1:]+loop[:1]):
                tag=(a,b,frozenset(loop))
                if tag in rejected:continue
                candidate=other_face(a,b,set(loop)-{a,b},10)
                if not candidate:rejected.add(tag);continue
                links=ring_edges(candidate)
                if links<=used:rejected.add(tag);continue
                options[frozenset(links)]=candidate
        if not options:break
        midpoint=np.mean([centroid(loop) for loop in loops],axis=0)
        selected=min(options.values(),key=lambda v:(-len(ring_edges(v)&used),len(v),np.linalg.norm(centroid(v)-midpoint),v))
        loops.append(selected);used.update(ring_edges(selected));site_images.update(sites_for(selected))
    return patches, topology
