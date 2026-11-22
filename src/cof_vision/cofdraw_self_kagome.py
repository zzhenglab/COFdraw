"""Finite Kagome imine drawings from the geometry of a two-in-one precursor.

Each vertex has two opposite aldehyde ports and two opposite amine ports.
Three rigidly rotated copies close the triangle, with the opposite arm of
every junction continuing into the neighboring hexagon and triangle.
"""
from functools import lru_cache
from itertools import product
import math

import numpy as np
from rdkit import Chem

from .cofdraw_depiction import compute_2d_coords, coordgen_coords
from .cofdraw_self import feed_metadata
from .depiction_rings import symmetric_depiction_rings
from .framework_drawing import _ports


HEX = ((1,0),(0,1),(-1,1),(-1,0),(0,-1),(1,-1))
TRI = ((0,0),(1,0),(0,1))
# Directed edges for the three junction orientations in an elementary triangle.
DIRECTIONS = ((1,0),(-1,1),(0,-1))


def _rotation(theta):
    s,c = math.sin(theta),math.cos(theta)
    return np.array([[c,s],[-s,c]])


def _opposed_aryl_hub(mol, ports, classes):
    """Four independent aryl arms on a 1,2,4,5 hub; each role occupies opposite positions."""
    for ring in mol.GetRingInfo().AtomRings():
        if len(ring)!=6 or not all(mol.GetAtomWithIdx(i).GetIsAromatic() and
                                    mol.GetRingInfo().NumAtomRings(i)==1 for i in ring):
            continue
        ends={}
        for group in classes:
            attachments=[]
            for port in group:
                target=ports[port][1]
                paths=[Chem.GetShortestPath(mol,i,target) for i in ring]
                route=min(paths,key=len)
                if (len(route)<5 or any(not (mol.GetAtomWithIdx(i).GetIsAromatic()
                        or mol.GetAtomWithIdx(i).GetHybridization() in (Chem.HybridizationType.SP,Chem.HybridizationType.SP2))
                        for i in route[1:-1])):
                    break
                # Require a single direct exit from the hub, and the actual
                # shortest path to this arm must originate at that exit.
                if any(i in ring for i in route[1:]) or sum(len(p)==len(route) for p in paths)!=1:
                    break
                attachments.append(ring.index(route[0]))
                ends[port]=(route[0],route[1])
            if len(attachments)!=2 or abs(attachments[0]-attachments[1])!=3:
                break
        else:
            if len({hub for hub, arm in ends.values()})==4:
                return ends
    return {}


@lru_cache(maxsize=128)
def _hub_branches(smiles, roles):
    mol=Chem.MolFromSmiles(smiles)
    ports=_ports(mol)
    classes=tuple([n for n,r in roles if r==role] for role in ('aldehyde','amine'))
    if any(len(group)!=2 for group in classes):
        return {}
    branches={}
    for n,(hub,first) in _opposed_aryl_hub(mol,ports,classes).items():
        # Cutting this one bond must isolate exactly one arm. This prevents
        # removal through a fused ring or a bridge to another attachment.
        pending=[first]
        branch=set()
        while pending:
            i=pending.pop()
            if i in branch:
                continue
            branch.add(i)
            pending.extend(a.GetIdx() for a in mol.GetAtomWithIdx(i).GetNeighbors()
                           if not (i==first and a.GetIdx()==hub) and a.GetIdx() not in branch)
        dummy=ports[n][0]
        if hub in branch or {i for i,_ in ports.values() if i in branch}!={dummy}:
            return {}
        branches[dummy]=(hub,first,branch-{dummy})
    return branches


def _crop_at_hubs(mol,recipe,origins):
    if (recipe.get('assembly') or {}).get('kind')!='self_kagome':
        return None
    part=recipe.get('components',{}).get('A',{})
    if len(origins)!=mol.GetNumAtoms() or not part.get('reactive_sites'):
        return None
    roles=tuple(sorted((r['port'],r['site_role']) for r in part['reactive_sites'] if 'site_role' in r))
    branches=_hub_branches(part['fragment_smiles'],roles)
    if not branches:
        return None
    indices={(o['instance'],o['fragment_atom_index']):i for i,o in enumerate(origins)}
    clipped=Chem.RWMol(mol)
    removed=set()
    count=0
    for port in mol.GetAtoms():
        if port.GetAtomicNum() or port.GetDegree()!=1:
            continue
        i=port.GetIdx()
        origin=origins[i]
        branch=branches.get(origin['fragment_atom_index'])
        if not branch:
            continue
        instance=origin['instance']
        hub,first,arm=branch
        if (instance,hub) not in indices or any((instance,j) not in indices for j in arm):
            return None
        cut={indices[(instance,j)] for j in arm}
        clipped.AddBond(indices[(instance,hub)],i,Chem.BondType.SINGLE)
        for conf in clipped.GetConformers():
            conf.SetAtomPosition(i,conf.GetAtomPosition(indices[(instance,first)]))
        removed.update(cut)
        count+=1
    if not count:
        return None
    clipped.BeginBatchEdit()
    for i in removed:
        clipped.RemoveAtom(i)
    clipped.CommitBatchEdit()
    clipped.UpdatePropertyCache(strict=True)
    symmetric_depiction_rings(clipped)
    return clipped.GetMol(),removed,count


def hub_boundary(mol,recipe,origins=None):
    """Figure-only copy, retaining the port and visible hub on each open arm."""
    if origins is None:
        if not all(a.HasProp('instance') and a.HasProp('fragment_atom_index') for a in mol.GetAtoms()):
            return None
        origins=[dict(instance=a.GetProp('instance'),fragment_atom_index=a.GetIntProp('fragment_atom_index'))
                 for a in mol.GetAtoms()]
    result=_crop_at_hubs(mol,recipe,origins)
    return result[0] if result else None


def hub_figure(mol,audit,recipe):
    from .cofdraw_boundary import remap_figure
    result=_crop_at_hubs(mol,recipe,audit['atom_origin'])
    if result:
        cropped,removed,count=result
        return (*remap_figure(mol,cropped,removed,audit,recipe),count)
    return None


def _clear_neighborhood(heavy_xy,step,rot):
    # All three Kagome vertex types are related by rotation. Include nearby
    # unconnected vertices and translated copies, particularly for long side chains.
    lattice=np.stack([step,-step@rot@rot])
    radius=np.max(np.linalg.norm(heavy_xy,axis=1))
    reach=2*radius+1.2
    limit=math.ceil(reach/np.linalg.svd(lattice,compute_uv=False).min())
    shapes=(heavy_xy,heavy_xy@rot,heavy_xy@rot@rot)
    for x,y in product(range(-limit,limit+1),repeat=2):
        kind=_kind((x,y))
        if kind==3 or (x==0 and y==0):
            continue
        offset=np.array([x,y])@lattice
        if np.linalg.norm(offset)>reach:
            continue
        delta=heavy_xy[:,None,:]-(shapes[kind]+offset)
        if np.any(np.sum(delta*delta,axis=2)<1.2**2):
            return False
    return True


def _placements(xy, ports, classes, rot, heavy):
    if any(np.max(np.abs(xy[list(ports[i])]+xy[list(ports[j])])) > .05 for i,j in classes):
        return []
    choices=[]
    for aldehyde,amine in product(*classes):
        d1,a1=ports[aldehyde]
        d2,a2=ports[amine]
        tangent=xy[d1]-xy[a1]
        receiving=(xy[a2]-xy[d2]) @ rot
        cosine=np.dot(tangent,receiving)/(np.linalg.norm(tangent)*np.linalg.norm(receiving))
        if cosine>=math.cos(math.radians(12)):
            direction=tangent+receiving
            step=xy[a1]-xy[a2]@rot + 1.5*direction/np.linalg.norm(direction)
            if np.linalg.norm(step)<1.5 or not _clear_neighborhood(xy[heavy],step,rot):
                continue
            choices.append((float(np.linalg.norm(step)),aldehyde,amine,step,xy))
    return choices


def _alternating_imine_ends(mol, xy, ports, classes):
    # Keep the entire chemical scaffold rigid. Each terminal trigonal atom
    # can draw its outward double bond 60 degrees either side of its aryl axis.
    for turns in product((-math.pi/3, math.pi/3), repeat=2):
        trial=xy.copy()
        for group,angle in zip(classes,turns):
            for port in group:
                dummy,anchor=ports[port]
                aryl=next(a.GetIdx() for a in mol.GetAtomWithIdx(anchor).GetNeighbors() if a.GetIdx()!=dummy)
                axis=xy[anchor]-xy[aryl]
                trial[dummy]=xy[anchor]+1.5*(axis@_rotation(angle))/np.linalg.norm(axis)
        yield trial


@lru_cache(maxsize=128)
def _cell(fragment_smiles, roles):
    mol=Chem.MolFromSmiles(fragment_smiles)
    ports=_ports(mol)
    aldehydes=[p for p,r in roles if r=='aldehyde']
    amines=[p for p,r in roles if r=='amine']
    classes=(aldehydes,amines)
    if len(aldehydes)!=2 or len(amines)!=2:
        raise ValueError('This Kagome arrangement requires two aldehyde arms and two amine arms.')
    # A flat drawing of an sp3 junction cannot establish planar four-arm
    # connectivity. Leave substituents outside the conjugated arms unrestricted.
    anchors=[a for _,a in ports.values()]
    for i,a in enumerate(anchors):
        for b in anchors[i+1:]:
            if any(mol.GetAtomWithIdx(j).GetHybridization() not in
                   (Chem.HybridizationType.SP,Chem.HybridizationType.SP2)
                   for j in Chem.GetShortestPath(mol,a,b)):
                raise ValueError('The reactive arms must share a conjugated planar backbone.')
    rot=_rotation(math.tau/3)
    compute_2d_coords(mol,canonOrient=False,forceRDKit=True)
    heavy=[a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum()>1]
    xy=mol.GetConformer().GetPositions()[:,:2].copy()
    xy-=xy[anchors].mean(axis=0)
    choices=_placements(xy,ports,classes,rot,heavy)
    if not choices and _opposed_aryl_hub(mol,ports,classes):
        for use_coordgen in (False,True):
            if use_coordgen:
                coordgen_coords(mol)
                xy=mol.GetConformer().GetPositions()[:,:2].copy()
                xy-=xy[anchors].mean(axis=0)
            for projection in (xy,xy*np.array([-1.,1.])):
                for trial in _alternating_imine_ends(mol,projection,ports,classes):
                    choices.extend(_placements(trial,ports,classes,rot,heavy))
            if choices:
                break
    if not choices:
        raise ValueError('No opposed, head-to-tail imine placement fits this Kagome drawing.')
    _,ald,ami,step,xy=max(choices,key=lambda choice:choice[:3])
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
                              pore_view=pore,
                              scope='A '+('hexagonal' if pore=='hexagon' else 'triangular')
                              +' pore of a proposed Kagome sheet. Rigid molecular blocks with head-to-tail imine bonds.'))


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
