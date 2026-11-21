"""Finite molecular circuits traced through a (3,4)-connected bor design graph."""
from collections import Counter
from functools import lru_cache

import numpy as np

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


def branch_patch(components,count):
    faces,edges,xy,degrees=circuit_layout()
    by_degree={c['connectivity']:k for k,c in components.items()}
    kinds=[by_degree[n] for n in degrees]
    geo={k:_block(c) for k,c in components.items()}
    faces=faces[:count]
    sites=sorted({v for face in faces for v in face})
    index={site:i for i,site in enumerate(sites)}
    kept=sorted(set.union(*(loop_edges(face) for face in faces)))
    joins=[(index[a],index[b]) for a,b in kept]
    ports={}
    orientations=[]
    for site in sites:
        adjacent=[b if a==site else a for a,b in edges if site in (a,b)]
        axes=xy[adjacent]-xy[site]
        _,handles,points,_=geo[kinds[site]]
        slots,theta=_orient(axes,handles,points)
        used=[b if a==site else a for a,b in kept if site in (a,b)]
        orientations.append(([slots[adjacent.index(other)] for other in used],theta))
        ports.update({(site,other):i for i,other in enumerate(used)})
    slot_pairs=[(ports[a,b],ports[b,a]) for a,b in kept]
    mol,origins=join_fragments(components,[kinds[s] for s in sites],xy[sites],joins,slot_pairs,2,orientations)
    rims=[[f'{kinds[v]}{index[v]+1}' for v in face] for face in faces]
    labels=[f'{len(f)}-node circuit' for f in faces]
    return mol,dict(pores=count,pore_faces=rims,face_labels=labels,pore_counts=dict(Counter(labels)),
                    atom_origin=origins,order_assumption='Finite circuits of a proposed bor incidence graph. The substrate geometry and crystal phase require independent validation.')
