"""Finite molecular circuits traced through a (3,4)-connected bor design graph."""
from collections import Counter
from copy import deepcopy
from functools import lru_cache
import json

import numpy as np
from rdkit import Chem

from .cofdraw_mixed import _block, _orient, join_fragments
from .cofdraw_net_graphs import neighbors, quotient


def loop_edges(loop):
    return {tuple(sorted(pair)) for pair in zip(loop,loop[1:]+loop[:1])}


@lru_cache(maxsize=1)
def circuit_layout():
    adjacency,_,_=quotient('bor')
    def cycles(path,size):
        if len(path)==size:
            if path[0] in neighbors(adjacency,path[-1]):
                yield path
            return
        for site in neighbors(adjacency,path[-1]):
            if site not in path:
                yield from cycles(path+[site],size)

    def unfold(face,start,edge,side):
        angle=side*2*np.pi/len(face)
        c,s=np.cos(angle),np.sin(angle)
        rotate=np.array([[c,-s],[s,c]])
        coords={}
        for site in face:
            coords[site]=start
            start=start+edge
            edge=rotate@edge
        return coords

    first=next(cycles([(0,0,0,0)],6))
    faces=[first]
    points=unfold(first,np.zeros(2),np.array([1.,0]),1)
    edges=loop_edges(first)
    while len(faces)<3:
        # Trace a neighboring six- or eight-cycle across one actual shared
        # edge. Unfold only the selected circuits for this chemical scheme.
        best=None
        for a,b in sorted(edges):
            for size in (6,8):
                for face in cycles([a,b],size):
                    if set(face)&points.keys()!={a,b}:
                        continue
                    for sign in (-1,1):
                        trial=unfold(face,points[a],points[b]-points[a],sign)
                        new=np.array([trial[v] for v in face[2:]])
                        old=np.array([p for k,p in points.items() if k not in (a,b)])
                        gap=np.linalg.norm(new[:,None]-old,axis=2).min()
                        score=(round(gap,5),-size)
                        if best is None or score>best[0]:
                            best=(score,face,trial)
        if best is None:
            raise ValueError('No neighboring circuit with a shared edge was found.')
        _,face,trial=best
        points.update({v:trial[v] for v in face[2:]})
        faces.append(face)
        edges.update(loop_edges(face))
    sites=sorted(points)
    joins=sorted((sites.index(a),sites.index(b)) for a,b in edges)
    return ([[sites.index(v) for v in f] for f in faces],joins,np.array([points[s] for s in sites]),
            [len(adjacency[p[0]]) for p in sites])


def graph_for(components,count,key):
    from .cofdraw_spatial import SPATIAL,patch_graph
    if key in SPATIAL:
        from .cofdraw_chemistry import inspect_monomer
        a=components['A']
        shape=inspect_monomer(a['monomer_smiles'])['site_shapes'].get(a['role'])
        return patch_graph(key,shape,count)
    faces,edges,xy,degrees=circuit_layout()
    by_degree={c['connectivity']:k for k,c in components.items()}
    kinds=[by_degree[n] for n in degrees]
    return faces[:count],edges,xy,kinds


def branch_patch(components,count,key='branched_4_3'):
    fields=('fragment_smiles','monomer_smiles','role','connectivity')
    signature=json.dumps({k:{f:c[f] for f in fields} for k,c in components.items()},sort_keys=True)
    mol,details=_cached_patch(signature,count,key)
    if mol is None:raise ValueError(details)
    return Chem.Mol(mol),deepcopy(details)


@lru_cache(maxsize=48)
def _cached_patch(signature,count,key):
    try:return _draw_patch(json.loads(signature),count,key)
    except ValueError as e:return None,str(e)


def _draw_patch(components,count,key):
    from .cofdraw_spatial import SPATIAL
    faces,edges,xy,kinds=graph_for(components,count,key)
    geo={k:_block(c) for k,c in components.items()}
    sites=sorted({v for face in faces for v in face})
    index={site:i for i,site in enumerate(sites)}
    kept=sorted(set.union(*(loop_edges(face) for face in faces)))
    joins=[(index[a],index[b]) for a,b in kept]
    def place(points,geometry):
        ports={}
        orientations=[]
        for site in sites:
            adjacent=[b if a==site else a for a,b in edges if site in (a,b)]
            axes=points[adjacent]-points[site]
            _,handles,local,_=geometry[kinds[site]]
            slots,theta=_orient(axes,handles,local)
            used=[b if a==site else a for a,b in kept if site in (a,b)]
            orientations.append(([slots[adjacent.index(other)] for other in used],theta))
            ports.update({(site,other):i for i,other in enumerate(used)})
        slot_pairs=[(ports[a,b],ports[b,a]) for a,b in kept]
        return join_fragments(components,[kinds[s] for s in sites],points[sites],joins,slot_pairs,2,orientations,geometry)
    mol,origins=place(xy,geo)
    if key in SPATIAL:
        from .cofdraw_spatial import check_projection,clear_junctions
        try:
            check_projection(mol)
        except ValueError as issue:
            original=Chem.MolToSmiles(mol)
            alternate={k:_block(c,alternate=True) for k,c in components.items()}
            for use,flip,shear in [('',1,0),('A',1,0),('B',1,0),('AB',1,0),('',-1,0),('AB',-1,0),
                                   ('',1,.65),('',1,-.65),('AB',1,.65),('AB',1,-.65)]:
                retry=bool(use) or flip<0 or shear!=0
                if retry:
                    geometry={k:alternate[k] if k in use else g for k,g in geo.items()}
                    mol,origins=place(xy@np.array([[1.,0],[shear,flip]]),geometry)
                    if Chem.MolToSmiles(mol)!=original:continue
                try:
                    if retry:check_projection(mol)
                    else:clear_junctions(mol)
                    break
                except ValueError:
                    if retry:
                        try:clear_junctions(mol);break
                        except ValueError:pass
            else:raise issue
    rims=[[f'{kinds[v]}{index[v]+1}' for v in face] for face in faces]
    labels=[f'{len(f)}-node circuit' for f in faces]
    return mol,dict(pores=count,pore_faces=rims,face_labels=labels,pore_counts=dict(Counter(labels)),
                    atom_origin=origins,order_assumption=f'Finite circuits of a proposed {SPATIAL.get(key,{}).get("code","bor")} net. The crystal phase requires independent validation.')
