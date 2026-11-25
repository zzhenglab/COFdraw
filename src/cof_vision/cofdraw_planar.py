"""Finite chemical faces traced in the planar embedding of a periodic net."""
from functools import lru_cache

import numpy as np

from .cofdraw_net_graphs import neighbors, quotient


PLANAR = {
    'mtf_4_3': ('mtf', (4, 3), ('quadrangular', 'triangular')),
    'kgd_6_3': ('kgd', (6, 3), ('hexagonal', 'triangular')),
}

# Explicit finite layouts for schematics whose projected boundary does not
# follow a shortest circuit of a periodic net. No lattice is assigned.
LOCAL_FACES={'local_square':4,'local_hexagon':6,'local_octagon':8}


@lru_cache(maxsize=3)
def local_faces(size):
    angle=np.arange(size)*2*np.pi/size
    vertices=np.column_stack((np.cos(angle),np.sin(angle)))
    points=list(vertices)
    faces=[list(range(size))]
    kinds=['AB'[i%2] for i in range(size)]
    for start in (0,size//2):
        end=(start+1)%size
        axis=vertices[end]-vertices[start]
        axis/=np.linalg.norm(axis)
        reflected=vertices[start]+(vertices-vertices[start])@(2*np.outer(axis,axis)-np.eye(2))
        face=[]
        for i,p in enumerate(reflected):
            if i in (start,end):
                face.append(i)
            else:
                face.append(len(points))
                points.append(p)
                kinds.append('AB'[i%2])
        faces.append(face)
    edges=sorted({tuple(sorted((v,face[(i+1)%size]))) for face in faces for i,v in enumerate(face)})
    return faces,edges,np.asarray(points),kinds


def templates():
    return [dict(id=key, label=f'{code} sheet · {degrees[0]}c + {degrees[1]}c',
                 allow_pairs=[list(degrees), list(reversed(degrees))], instances=32 if code == 'mtf' else 4,
                 shape_pair=list(shapes), dimension=2, topology_label=code, linkages=['imine'],
                 description=f'Connected faces of a proposed {code} sheet. Junction angles can depart from ideal monomer symmetry.')
            for key, (code, degrees, shapes) in PLANAR.items()]


@lru_cache(maxsize=len(PLANAR))
def face_layout(code):
    adjacency, fractional, cell = quotient(code)

    def point(site):
        return (fractional[site[0]] + site[1:]) @ cell

    def ordered(site):
        def angle(other):
            dx, dy = point(other) - point(site)
            return np.arctan2(dy, dx)
        return sorted(neighbors(adjacency, site), key=angle)

    def boundary(a, b):
        start = a, b
        face = []
        while a not in face:
            face.append(a)
            adjacent = ordered(b)
            a, b = b, adjacent[adjacent.index(a) - 1]
        if (a, b) != start:
            raise ValueError('Invalid planar face boundary.')
        return face

    root = (0, 0, 0)
    faces = [boundary(root, ordered(root)[0])]
    for _ in range(2):
        candidates = {}
        seen = {frozenset(face) for face in faces}
        for face in faces:
            for a, b in zip(face, face[1:] + face[:1]):
                adjacent = boundary(b, a)
                key = frozenset(adjacent)
                if key not in seen:
                    candidates[key] = adjacent
        # Grow toward the patch center so an expanded view stays compact.
        center = np.mean([point(v) for face in faces for v in face], axis=0)
        sizes = {len(face) for face in faces}
        def rank(face):
            distance = np.linalg.norm(np.mean([point(v) for v in face], axis=0) - center)
            return len(face) in sizes, round(float(distance), 8), tuple(sorted(face))
        faces.append(min(candidates.values(), key=rank))
    core = {site for face in faces for site in face}
    sites = sorted(core | {v for site in core for v in ordered(site)})
    index = {site: i for i, site in enumerate(sites)}
    edges = {tuple(sorted((index[a], index[b]))) for a in core for b in ordered(a)}
    return ([[index[v] for v in face] for face in faces], sorted(edges),
            np.array([point(v) for v in sites]), [len(adjacency[v[0]]) for v in sites])
