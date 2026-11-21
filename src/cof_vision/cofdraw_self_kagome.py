"""Finite Kagome imine drawings of the two-in-one monomers in Li et al., Fig. 6.

Each vertex has two opposite aldehyde ports and two opposite amine ports.
Three rigidly rotated copies close the triangle, with the opposite arm of
every junction continuing into the neighboring hexagon and triangle.
"""
from functools import lru_cache
from itertools import product
import math

import numpy as np
from rdkit import Chem

from .cofdraw_depiction import compute_2d_coords
from .cofdraw_self import LI_DOI, feed_metadata
from .framework_drawing import _ports


HEX = ((1,0),(0,1),(-1,1),(-1,0),(0,-1),(1,-1))
TRI = ((0,0),(1,0),(0,1))
# Directed edges for the three junction orientations in an elementary triangle.
DIRECTIONS = ((1,0),(-1,1),(0,-1))


def _rotation(theta):
    s,c = math.sin(theta),math.cos(theta)
    return np.array([[c,s],[-s,c]])


@lru_cache(maxsize=8)
def _cell(fragment_smiles, roles):
    mol = Chem.MolFromSmiles(fragment_smiles)
    ports = _ports(mol)
    compute_2d_coords(mol,canonOrient=False,forceRDKit=True)
    xy = mol.GetConformer().GetPositions()[:, :2].copy()
    xy -= xy[[a for _,a in ports.values()]].mean(axis=0)
    aldehydes = [p for p,r in roles if r=='aldehyde']
    amines = [p for p,r in roles if r=='amine']
    if len(aldehydes)!=2 or len(amines)!=2:
        raise ValueError('This Kagome arrangement requires two aldehyde arms and two amine arms.')
    # The opposed ends must overlay by inversion, including the terminal
    # bond direction. This is a molecular symmetry check on these drawings.
    for i,j in (aldehydes,amines):
        if np.max(np.abs(xy[list(ports[i])]+xy[list(ports[j])])) > .05:
            raise ValueError('The selected fragment does not have opposed equivalent reactive arms.')
    rot = _rotation(math.tau/3)
    choices = []
    for aldehyde,amine in product(aldehydes,amines):
        d1,a1 = ports[aldehyde]
        d2,a2 = ports[amine]
        tangent = xy[d1]-xy[a1]
        receiving = (xy[a2]-xy[d2]) @ rot
        cosine = np.dot(tangent,receiving) / (np.linalg.norm(tangent)*np.linalg.norm(receiving))
        if cosine < math.cos(math.radians(12)):
            continue
        direction = tangent+receiving
        step = xy[a1]-xy[a2]@rot + 1.5*direction/np.linalg.norm(direction)
        choices.append((float(np.linalg.norm(step)),aldehyde,amine,step))
    if not choices:
        raise ValueError('No head-to-tail imine placement fits this Kagome drawing.')
    _,ald,ami,step = max(choices,key=lambda choice:choice[:3])
    # Align the row of hexagons with the page while preserving every bond.
    align = _rotation(-math.atan2(step[1],step[0]))
    lattice = np.stack([step, -step@rot@rot])@align
    shapes,sites = [],[]
    for k in range(3):
        shapes.append(xy@_rotation(k*math.tau/3)@align)
        positive = DIRECTIONS[k]
        received = tuple(-x for x in DIRECTIONS[k-1])
        sites.append({positive:ald, tuple(-x for x in positive):next(p for p in aldehydes if p!=ald),
                      received:ami, tuple(-x for x in received):next(p for p in amines if p!=ami)})
    return mol,ports,lattice,shapes,sites


def _hexagon(u,v):
    return [(2*u+1+x,2*v+1+y) for x,y in HEX]


def _kind(point):
    x,y = point
    return x%2 + 2*(y%2)


def _incidence(component, points):
    roles = {row['port']:row['site_role'] for row in component['reactive_sites']}
    mol,ports,lattice,shapes,sites = _cell(component['fragment_smiles'],tuple(sorted(roles.items())))
    centers = {p:np.asarray(p)@lattice for p in points}
    coords = {p:centers[p]+shapes[_kind(p)] for p in points}
    edges = []
    for x,y in points:
        p = (x,y)
        for dx,dy in HEX[:3]:
            q = (x+dx,y+dy)
            if q not in centers:
                continue
            a,b = sites[_kind(p)][dx,dy],sites[_kind(q)][-dx,-dy]
            distance = np.linalg.norm(coords[p][ports[a][1]]-coords[q][ports[b][1]])
            if roles[a]==roles[b] or abs(distance-1.5)>.08:
                raise ValueError('A Kagome edge failed its complementary imine or drawing bond-length check.')
            edges.append((p,a,q,b))
    return mol,ports,centers,coords,edges


def kagome_recipe(row,component,pore):
    points = _hexagon(0,0) if pore=='hexagon' else list(TRI)
    mol,ports,centers,coords,edges = _incidence(component,points)
    by_pair = {(p,q):a for p,a,q,_ in edges}
    by_pair.update({(q,p):b for p,_,q,b in edges})
    ids = {p:f'A{i+1}' for i,p in enumerate(points)}
    instances, coordinates = [],[]
    for i,p in enumerate(points):
        a,b = by_pair[p,points[i-1]],by_pair[p,points[(i+1)%len(points)]]
        instances.append(dict(id=ids[p],component='A',incoming_port=a,outgoing_port=b,layout_position=centers[p].tolist()))
        for atom in mol.GetAtoms():
            j = atom.GetIdx()
            if j not in (ports[a][0],ports[b][0]):
                coordinates.append(dict(instance=ids[p],fragment_atom_index=j,xy=coords[p][j].tolist()))
    joins = [dict(a=dict(instance=ids[p],port=a),b=dict(instance=ids[q],port=b),bond_order=2)
             for p,a,q,b in edges]
    return dict(components={'A':component},instances=instances,connections=joins,atom_coordinates=coordinates,
                coordinate_provenance='generated Kagome layer from rotated precursor fragments',
                layout_kind='closed_pore',linkage='imine',selected_template=None,design_topology=None,
                assembly=dict(monomer_ids=[row['id']],kind='self_kagome',precursors=[feed_metadata(row)],linker_sequence='',
                              pore_view=pore,source_doi=LI_DOI,
                              scope='A '+('hexagonal' if pore=='hexagon' else 'triangular')
                              +' pore in the reported Kagome arrangement. Rigid molecular blocks with head-to-tail imine bonds.'))


def patch_points(repeats,triad=False,star=False):
    if triad or star:
        # Each side of the central hexagon bounds one outward triangular face.
        points=set(_hexagon(0,0))
        for i,(ax,ay) in enumerate(HEX if star else HEX[:2]):
            bx,by=HEX[(i+1)%6]
            points.add((1+ax+bx,1+ay+by))
        return points
    return {p for u,v in product(range(repeats),repeat=2) for p in _hexagon(u,v)}


def kagome_patch(recipe,repeats,triad=False,star=False):
    """Draw the selected connected neighborhood, joining shared vertices only once."""
    points = sorted(patch_points(repeats,triad,star))
    mol,ports,centers,coords,edges = _incidence(recipe['components']['A'],points)
    if len(points)*mol.GetNumAtoms()-2*len(edges)>2400:
        raise ValueError('Choose a smaller multi-pore extent for this molecular size.')
    tagged = {}
    for tag,(p,a,q,b) in enumerate(edges,1):
        tagged[p,a] = tagged[q,b] = tag
    joined = Chem.RWMol()
    for i,p in enumerate(points):
        copy = Chem.Mol(mol)
        for a in copy.GetAtoms():
            a.SetIntProp('fragment_atom_index',a.GetIdx())
            a.SetProp('instance',f'A{i+1}')
            a.SetProp('component','A')
            if a.GetAtomicNum()==0:
                a.SetAtomMapNum(tagged.get((p,a.GetAtomMapNum()),0))
        for j,(x,y) in enumerate(coords[p]):
            copy.GetConformer().SetAtomPosition(j,(float(x),float(y),0.))
        if i==0:
            joined = copy
        else:
            joined = Chem.CombineMols(joined,copy)
    result = Chem.molzip(joined)
    if len(Chem.GetMolFrags(result))!=1 or result.GetNumAtoms()!=len(points)*mol.GetNumAtoms()-2*len(edges):
        raise ValueError('The Kagome patch failed its connectivity check.')
    result.UpdatePropertyCache(strict=True)
    Chem.GetSSSR(result)
    origins = [dict(instance=a.GetProp('instance'),component='A',fragment_atom_index=a.GetIntProp('fragment_atom_index'))
               for a in result.GetAtoms()]
    # E - V + 1 counts all bounded faces in this connected planar junction graph.
    face_count = len(edges)-len(points)+1
    hexagons = 1 if (triad or star) else repeats*repeats
    return result,dict(blocks=len(points),pore_count=face_count,
                       pore_counts={'triangle':face_count-hexagons,'hexagon':hexagons}, atom_origin=origins)
