"""Explicit single-precursor and mixed-linker finite imine designs."""
import math

import numpy as np
from rdkit import Chem

from .cofdraw_chemistry import best_port_pair
from .cofdraw_depiction import compute_2d_coords
from .framework_drawing import _ports, linkage_fragment

ALDEHYDE, AMINE = 'aryl_aldehyde', 'primary_arylamine'


def _fragment(row, role):
    smiles = row['smiles']
    if role == 'self_imine' and row.get('self_activation'):
        mol = Chem.MolFromSmiles(smiles)
        rw = Chem.RWMol(mol)
        acetal = row['self_activation'] == 'dimethyl_acetal'
        query = '[c][CH]([O][CH3])[O][CH3]' if acetal else '[c][NH][C](=[O])[O][C]([CH3])([CH3])[CH3]'
        sites = mol.GetSubstructMatches(Chem.MolFromSmarts(query))
        if len(sites) != 2:
            raise ValueError('The precursor does not match the two documented protecting groups.')
        remove = []
        for hit in sites:
            if acetal:
                rw.GetBondBetweenAtoms(hit[1],hit[2]).SetBondType(Chem.BondType.DOUBLE)
                remove.extend(hit[3:])
            else:
                remove.extend(hit[2:])
        for idx in sorted(set(remove),reverse=True):
            rw.RemoveAtom(idx)
        activated = rw.GetMol()
        Chem.SanitizeMol(activated)
        smiles = Chem.MolToSmiles(activated)
    component = linkage_fragment(smiles, 'imine', role)
    if smiles != row['smiles']:
        component['activated_smiles'] = component['monomer_smiles']
        component['monomer_smiles'] = row['canonical_smiles']
        component['activation'] = row['self_activation']
    component['label'] = (row.get('abbreviation') or row['label'])[:100]
    return component


def _port_span(component, angle):
    incoming, outgoing, _ = best_port_pair(component, angle)
    mol = Chem.MolFromSmiles(component['fragment_smiles'])
    compute_2d_coords(mol)
    xy = mol.GetConformer().GetPositions()
    ports = _ports(mol)
    return float(np.linalg.norm(xy[ports[incoming][0]] - xy[ports[outgoing][0]]))


def _polygon(lengths):
    # A cyclic polygon gives different chord lengths a common radius and closes exactly.
    low = max(lengths)/2
    def angles(radius):
        return [2*math.asin(min(1, length/(2*radius))) for length in lengths]
    if sum(angles(low)) < math.tau:
        raise ValueError('These linker spans do not fit a convex pore in this layout.')
    high = sum(lengths)
    for _ in range(48):
        mid = (low+high)/2
        if sum(angles(mid)) > math.tau:
            low = mid
        else:
            high = mid
    turns = np.r_[0., np.cumsum(angles(high)[:-1])] + math.pi/2
    return high*np.column_stack([np.cos(turns), np.sin(turns)])


def assemble(rows, pattern, pore_view=None, pore_layout=None):
    if len({r['canonical_smiles'] for r in rows}) != len(rows):
        raise ValueError('Select distinct precursors for the node and both linkers.')
    components, ordered = {}, []
    if len(rows) == 1:
        if pore_layout:
            raise ValueError('A mixed-linker arrangement requires three precursors.')
        from .cofdraw_self import feed_metadata, trimer_assembly
        row = rows[0]
        kind = row['assembly_roles'].get('self_kind')
        if pore_view is not None and not (kind == 'self_kagome' and pore_view in ('hexagon','triangle')):
            raise ValueError('This pore view is unavailable for the selected self-condensation.')
        if kind in ('self_boroxine','self_triazine'):
            result = trimer_assembly(row,kind)
            result['assembly']['precursors'] = [feed_metadata(row)]
            return result
        if not row['assembly_roles']['self_imine']:
            raise ValueError('Choose a self-condensing aminoaldehyde, linear aryl diboronic acid or linear aryl dinitrile.')
        comp = components['A'] = _fragment(row, 'self_imine')
        if kind == 'self_kagome':
            from .cofdraw_self_kagome import kagome_recipe
            return kagome_recipe(row,comp,pore_view or 'hexagon')
        chain = comp['connectivity'] == 2
        # Coordinates give directions to the rigid-block fitter. Compact unit
        # positions let it solve the length scale for the actual precursor.
        if chain:
            points = [[i, 0] for i in range(6)]
            angles = [180]*len(points)
        else:
            # Two neighboring arms choose the acute and obtuse turns of the
            # same rhombus. Mapping every corner to 90° straightens pyrene.
            _,_,fit = best_port_pair(comp,65)
            acute = min(fit['depicted_site_angle'],180-fit['depicted_site_angle'])
            angle = math.radians(max(50.,acute))
            u,v = np.array([1.,0.]),np.array([math.cos(angle),math.sin(angle)])
            points = [np.array([0.,0.]), u,u+v,v]
            angles = [math.degrees(angle),180-math.degrees(angle)]*2
        ordered = [('A',p,a) for p,a in zip(points,angles)]
        layout = 'open_chain' if chain else 'closed_pore'
        scope = 'Head-to-tail imine chain.' if chain else 'Four corners of one proposed self-condensed imine pore.'
        if comp.get('activation'):
            scope += ' Includes removal of two '+comp['activation']+' protecting groups before imine formation.'
    elif len(rows) == 3:
        if pore_view is not None or pore_layout not in (None,'kgm'):
            raise ValueError('Choose a supported mixed-linker pore layout and B/C sequence.')
        roles = [r['assembly_roles'] for r in rows]
        node_role = next((role for role, other in [(AMINE, ALDEHYDE), (ALDEHYDE, AMINE)]
                          if roles[0][role] in (3,4) and all(r[other] == 2 for r in roles[1:])), None)
        if not node_role:
            raise ValueError('Choose a three- or four-connected node and two complementary ditopic linkers.')
        edge_role = ALDEHYDE if node_role == AMINE else AMINE
        components = {k:_fragment(row, role) for k,row,role in zip('ABC',rows,[node_role,edge_role,edge_role])}
        corners = 6 if components['A']['connectivity'] == 3 else 4
        if any(r.get('geometry_hint','').startswith('tetrahedral') for r in rows):
            raise ValueError('This planar pore design needs a node with planar connectivity.')
        if any(r['assembly_roles']['self_imine'] for r in rows):
            raise ValueError('Competing aminoaldehyde roles require an explicit selective reaction design.')
        if pore_layout=='kgm':
            from .cofdraw_mixed import ETTA, linear_dialdehyde
            if rows[0]['canonical_smiles']!=Chem.MolToSmiles(Chem.MolFromSmiles(ETTA)):
                raise ValueError('This kagome arrangement requires the four planar ethylenetetraaniline (ETTA) arms.')
            if not all(linear_dialdehyde(r['canonical_smiles']) for r in rows[1:]):
                raise ValueError('Choose two linear dialdehydes with outward-facing formyl groups for this kagome design.')
            if pattern not in ('BBB','CCC','BCBCBC','CBCBCB'):
                raise ValueError('Choose the B triangle, C triangle, or alternating hexagon.')
            corners=len(pattern)
        elif len(pattern) != corners or set(pattern) != {'B','C'}:
            raise ValueError(f'Enter {corners} linker positions using both B and C.')
        interior = 180-360/corners
        # Two node radii along each edge, plus the actual linker attachment span.
        node_span = _port_span(components['A'],interior)/math.sin(math.radians(interior/2))
        spans = {k:node_span+_port_span(components[k],180) for k in 'BC'}
        vertices = _polygon([spans[k] for k in pattern])
        # Span ratios determine the polygon's shape. Its absolute scale must
        # remain free for build_patch to fit the actual imine junctions.
        vertices /= min(spans[k] for k in pattern)
        for i,k in enumerate(pattern):
            p,q = vertices[i],vertices[(i+1)%corners]
            u,v = vertices[i-1]-p,q-p
            angle = math.degrees(math.acos(np.clip(u@v/(np.linalg.norm(u)*np.linalg.norm(v)),-1,1)))
            ordered.extend([('A',p,angle),(k,(p+q)/2,180)])
        scope = (f'Kagome design · {pattern}. B and C form distinct triangular faces around an alternating hexagon.'
                 if pore_layout=='kgm' else f'Local imine pore · {pattern}. The ordered neighborhood is one B/C placement hypothesis.')
        layout = 'closed_pore'
    else:
        raise ValueError('Choose one self-condensing precursor or one node and two linkers.')
    instances = []
    for i,(key,point,angle) in enumerate(ordered):
        component = components[key]
        incoming,outgoing,_ = best_port_pair(component,angle)
        instances.append(dict(id=f'{key}{i+1}',component=key,incoming_port=incoming,outgoing_port=outgoing,
                              layout_position=list(map(float,point))))
    edges=[]
    for i in range(len(instances)-(layout=='open_chain')):
        left,right=instances[i],instances[(i+1)%len(instances)]
        edges.append({'a':{'instance':left['id'],'port':left['outgoing_port']},
                      'b':{'instance':right['id'],'port':right['incoming_port']},'bond_order':2})
    return dict(components=components,instances=instances,connections=edges,linkage='imine',layout_kind=layout,
                assembly={'monomer_ids':[r['id'] for r in rows], 'linker_sequence':pattern if len(rows)==3 else '',
                          'scope':scope, **({'kind':'self_imine','precursors':[feed_metadata(rows[0])]} if len(rows)==1 else
                             {'kind':'mixed_imine', **({'pore_layout':pore_layout} if pore_layout else {})})},
                selected_template=None,design_topology=None)
