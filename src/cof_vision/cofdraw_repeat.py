"""Finite drawings with shared molecular blocks and explicitly reacted boundaries."""
from collections import Counter
from copy import deepcopy
from functools import lru_cache
from itertools import product
import hashlib
import json
import math

import numpy as np
from rdkit import Chem
from rdkit.Chem.Draw import rdMolDraw2D

from .cofdraw_chemistry import design_signature
from .cofdraw_depiction import compute_2d_coords
from .cofdraw_periodic_build import SUPPORTED_TEMPLATES, _fit_blocks, _net, _rotate
from .cofdraw_periodic_options import validated_design
from .depiction_rings import symmetric_depiction_rings
from .framework_drawing import _double_bond_geometry_audit, _ports, build_patch

from .cofdraw_spatial import SPATIAL

BRANCHED=('branched_4_3',*SPATIAL)
MAX_ATOMS = 2400
SHEET_TEMPLATES = (*SUPPORTED_TEMPLATES, 'hxl_6_2')


def _planar_block(component):
    mol = Chem.MolFromSmiles(component['fragment_smiles'])
    compute_2d_coords(mol, canonOrient=False, forceRDKit=True)
    ports = list(_ports(mol).values())
    xyz = mol.GetConformer().GetPositions()
    sites = np.array([site for _, site in ports])
    xyz -= xyz[sites].mean(axis=0)
    directions = xyz[[dummy for dummy, _ in ports], :2] - xyz[sites, :2]
    directions /= np.linalg.norm(directions, axis=1)[:, None]
    return dict(mol=mol, xyz=xyz, sites=sites, ports=ports, directions=directions)


def _classification(recipe, source_rule=None, source_sha256=None):
    if source_sha256:
        if not source_rule or source_rule.get('source_sha256') != source_sha256:
            raise ValueError('This source drawing has no reviewed rule connecting its boundary sites to neighboring pores.')
        return source_rule['template']
    if recipe.get('explicit_graph') or recipe.get('design_template_signature') != design_signature(recipe):
        raise ValueError('Choose a connectivity template to describe how neighboring drawing units join.')
    assembly = recipe.get('assembly')
    if not assembly:
        template = (recipe.get('selected_template') or {}).get('id')
        if template in (*SHEET_TEMPLATES, 'chain_2_2','kgm_4_2',*BRANCHED):
            return template
        raise ValueError('This local ring has no declared rule for sharing edges with adjacent pores.')
    kind = assembly.get('kind')
    if kind in ('self_boroxine','self_triazine'):
        parts = recipe['components']
        if (set(parts) != {'A','B'} or [parts[k].get('role') for k in 'AB'] != ['product_junction','product_spacer']
                or recipe.get('linkage') != kind.removeprefix('self_') or len(assembly.get('monomer_ids',[])) != 1):
            raise ValueError('The displayed building blocks disagree with the selected self-condensation.')
        return kind
    if len(assembly.get('monomer_ids', [])) == 3:
        pattern = assembly.get('linker_sequence')
        if pattern != ''.join(i['component'] for i in recipe['instances'] if i['component']!='A'):
            raise ValueError('The declared linker sequence differs from the displayed local pore.')
        if assembly.get('pore_layout') == 'kgm':
            from .cofdraw_mixed import ETTA
            if Chem.MolToSmiles(Chem.MolFromSmiles(ETTA)) != recipe['components']['A']['monomer_smiles']:
                raise ValueError('The selected Kagome rule requires a planar ETTA node.')
            return 'mixed_kgm'
        return 'mixed_linkers'
    if set(recipe['components']) == {'A'} and recipe['components']['A'].get('role') == 'self_imine':
        if kind=='self_kagome':
            from .cofdraw_self import KAGOME_SMILES
            if recipe['components']['A']['monomer_smiles'] not in KAGOME_SMILES:
                raise ValueError('This aminoaldehyde has no registered triangular/hexagonal assembly rule.')
            return kind
        return 'self_chain' if recipe['layout_kind'] == 'open_chain' else 'self_square'
    raise ValueError('This assembly has no defined neighboring-pore rule.')


def _quotient(template, components):
    if template == 'self_square':
        _, cell, centers, edges = _net('sql_4_4', {})
        return ['A', 'A'], cell, centers, edges
    counts = {k:c['connectivity'] for k,c in components.items()}
    if template != 'hxl_6_2':
        return _net(template, counts)
    cell = np.array([[1., 0.], [.5, math.sqrt(3)/2]])
    node = max(counts, key=counts.get)
    linker = 'B' if node == 'A' else 'A'
    directions = [(1,0), (0,1), (-1,1)]
    centers = [np.zeros(2), *[np.array(step) @ cell / 2 for step in directions]]
    edges = [(0, i, (0,0)) for i in (1,2,3)]
    edges += [(i, 0, step) for i,step in enumerate(directions,1)]
    return [node, linker, linker, linker], cell, centers, edges


def _faces(template, n, cell, centers, triad=False):
    """Look up face boundaries in the infinite ideal net, then deduplicate their blocks."""
    def key(xy):
        return tuple(np.round(xy, 6))
    lookup = {key(np.array([u,v]) @ cell + pos):(u,v,k)
              for u,v,k in product(range(-n-2,n+3),range(-n-2,n+3),range(len(centers)))
              for pos in [centers[k]]}
    closed = []
    if template == 'hxl_6_2':
        for i in range(n):
            for j in range(n-i):
                p = [i,j] @ cell
                closed.append([p,p+cell[0],p+cell[1]])
                if i+j < n-1:
                    closed.append([p+cell[0],p+cell[0]+cell[1],p+cell[1]])
    elif template.startswith('hcb'):
        vertices = np.array([[0.,0.],[1.,0.],[1.5,math.sqrt(3)/2],
                             [1.,math.sqrt(3)],[0.,math.sqrt(3)],[-.5,math.sqrt(3)/2]])
        closed = [vertices + [u,v] @ cell for u,v in product(range(n),repeat=2)]
    else:
        vertices = np.array([[0.,0.],[1.,0.],[1.,1.],[0.,1.]])
        closed = [vertices + [u,v] for u,v in product(range(n),repeat=2)]
    # In a 2 × 2 sheet the first three faces share edges around a junction.
    # This compact crop is distinct from the legacy 3 × 3 extent.
    if triad:
        closed = closed[:3]
    insert_linkers = template in ('sql_4_2','hcb_3_2','hxl_6_2')
    result = []
    for vertices in closed:
        face = []
        for i,p in enumerate(vertices):
            face.append(lookup[key(p)])
            if insert_linkers:
                face.append(lookup[key((p + vertices[(i+1)%len(vertices)]) / 2)])
        result.append(face)
    return result


def _id(key):
    return '.'.join(map(str,key))


def _bounded_faces(mol, origins, links):
    """Follow clockwise successors from each directed edge in the block drawing."""
    points = {}
    xy = mol.GetConformer().GetPositions()[:,:2]
    for atom,origin in zip(mol.GetAtoms(),origins):
        if atom.GetAtomicNum():
            points.setdefault(origin['instance'],[]).append(xy[atom.GetIdx()])
    centers = {k:np.mean(v,axis=0) for k,v in points.items()}
    neighbors = {k:[] for k in centers}
    for bond in links:
        a,b = (origins[i]['instance'] for i in (bond.GetBeginAtomIdx(),bond.GetEndAtomIdx()))
        neighbors[a].append(b)
        neighbors[b].append(a)
    for node,adjacent in neighbors.items():
        adjacent.sort(key=lambda v:math.atan2(*(centers[v]-centers[node])[::-1]))
    visited, faces = set(), []
    for start,adjacent in neighbors.items():
        for end in adjacent:
            a,b = start,end
            if (a,b) in visited:
                continue
            cycle = []
            while (a,b) not in visited:
                visited.add((a,b))
                cycle.append(a)
                successors = neighbors[b]
                a,b = b,successors[successors.index(a)-1]
            path = np.array([centers[k] for k in cycle])
            area = np.sum(path[:,0]*np.roll(path[:,1],-1) - path[:,1]*np.roll(path[:,0],-1))
            if area>1e-6 and len(cycle)==len(set(cycle)):
                faces.append(cycle)
    return faces


def _sheet_layout(design, template, repeats, extent_mode):
    kinds, cell, centers, edges = _quotient(template, design['components'])
    if extent_mode == 'triad':
        faces = _faces(template,2,cell,centers,triad=True)
    else:
        faces = _faces(template,repeats,cell,centers) if extent_mode == 'pores' else []
    nodes = set(k for face in faces for k in face) if faces else set(product(range(repeats),range(repeats),range(len(kinds))))
    return kinds, cell, centers, edges, faces, nodes


def _self_states(block):
    # Antiparallel or quarter-turned sublattices cover unlike ports on both
    # alternating and adjacent 2+2 arrangements. Keep each fragment rigid.
    for angle in (0., math.pi/2, math.pi, -math.pi/2):
        other = dict(block, xyz=block['xyz'].copy())
        other['xyz'][:,:2] = _rotate(block['xyz'][:,:2],angle)
        other['directions'] = _rotate(block['directions'],angle)
        yield [block, other]


def _chain(recipe, repeats):
    extended = deepcopy(recipe)
    sequence = recipe['instances']
    count = len(sequence) * repeats
    order = recipe['connections'][0]['bond_order']
    items = [{**sequence[i % len(sequence)], 'id': str(i), 'layout_position': [12*i, 0]} for i in range(count)]
    extended.update(instances=items, atom_coordinates=None, connections=[
        {'a': {'instance': a['id'], 'port': a['outgoing_port']},
         'b': {'instance': b['id'], 'port': b['incoming_port']}, 'bond_order': order}
        for a,b in zip(items,items[1:])])
    mol, audit = build_patch(extended)
    return mol, dict(atom_origin=audit['atom_origin'],blocks=count,pore_count=0)


@lru_cache(maxsize=64)
def _fit_sheet(template,fragments,order):
    # The placement is independent of drawing size and style. Reuse it in
    # the extent picker and the actual drawing, including unsuccessful fits.
    components = {key:dict(fragment_smiles=smiles,connectivity=count) for key,smiles,count in fragments}
    kinds,ref,centers,edges = _quotient(template,components)
    blocks = {key:_planar_block(c) for key,c in components.items()}
    poses = []
    cell_blocks = [blocks[k] for k in kinds]
    states = _self_states(cell_blocks[0]) if template == 'self_square' else [cell_blocks]
    for state in states:
        for pose in _fit_blocks(state,ref,centers,edges,1.35 if order == 2 else 1.5):
            if template == 'self_square':
                # Check the actual endpoint elements at every translated C=N join.
                if any({state[k]['mol'].GetAtomWithIdx(state[k]['ports'][s][1]).GetAtomicNum()
                        for k,s in zip((a,b),slots)} != {6,7}
                       for (a,b,_),slots in zip(edges,pose[5])):
                    continue
            poses.append((pose[0],pose,state))
    return min(poses,key=lambda item:item[0])[1:] if poses else None


def _sheet_placement(design,template):
    fragments = tuple((key,c['fragment_smiles'],c['connectivity']) for key,c in sorted(design['components'].items()))
    placement = _fit_sheet(template,fragments,design['connections'][0]['bond_order'])
    if placement is None:
        raise ValueError('The reactive sites do not align with this planar net. Choose a different connectivity to extend the drawing.')
    return placement


def _sheet(design,template,repeats,extent_mode):
    kinds, _, _, edges, faces, nodes = _sheet_layout(design,template,repeats,extent_mode)
    pose,cell_blocks = _sheet_placement(design,template)
    counts = sum(cell_blocks[k]['mol'].GetNumAtoms() for _,_,k in nodes)
    face_edges = {frozenset((a,b)) for face in faces for a,b in zip(face,face[1:]+face[:1])}
    if (counts - 2*len(face_edges) if faces else counts) > MAX_ATOMS:
        raise ValueError('This repeat exceeds the 2,400-atom drawing limit. Choose a smaller extent.')
    _, lattice, fitted, angles, mirror, slots, *_ = pose
    fragments = {}
    for key in sorted(nodes):
        u,v,k = key
        block = cell_blocks[k]
        mol = Chem.Mol(block['mol'])
        xy = _rotate(block['xyz'][:,:2]*mirror[k],angles[k]) + fitted[k] + [u,v] @ lattice
        for atom,(x,y) in zip(mol.GetAtoms(),xy):
            atom.SetProp('instance',_id(key))
            atom.SetProp('component',kinds[k])
            atom.SetIntProp('fragment_atom_index',atom.GetIdx())
            mol.GetConformer().SetAtomPosition(atom.GetIdx(),(float(x),float(y),0.))
            if atom.GetAtomicNum() == 0:
                atom.SetAtomMapNum(0)
        fragments[key] = mol
    joined_edges = set()
    for u,v,k in fragments:
        for (a,b,(du,dv)),(sa,sb) in zip(edges,slots):
            left,right = (u,v,k),(u+du,v+dv,b)
            pair = frozenset((left,right))
            if k != a or right not in fragments or (faces and pair not in face_edges):
                continue
            joined_edges.add(pair)
            for key,slot in ((left,sa),(right,sb)):
                atom = fragments[key].GetAtomWithIdx(cell_blocks[key[2]]['ports'][slot][0])
                if atom.GetAtomMapNum():
                    raise ValueError('A reactive site would be consumed twice.')
                atom.SetAtomMapNum(len(joined_edges))
    if faces and joined_edges != face_edges:
        raise ValueError('Declared pore boundaries do not match the net edges.')
    parts = iter(fragments.values())
    joined = next(parts)
    for mol in parts:
        joined = Chem.CombineMols(joined,mol)
    mol = Chem.molzip(joined)
    if len(Chem.GetMolFrags(mol)) != 1 or mol.GetNumAtoms() != counts - 2*len(joined_edges):
        raise ValueError('Repeat connectivity could not be verified.')
    mol.UpdatePropertyCache(strict=True)
    symmetric_depiction_rings(mol)
    nfaces = len(joined_edges) - len(fragments) + 1
    if faces and nfaces != len(faces):
        raise ValueError('Shared edges do not reproduce the declared pores.')
    return mol, dict(blocks=len(fragments),pore_count=nfaces,
                     pore_faces=[[_id(k) for k in face] for face in faces],
                     interblock_bonds=len(joined_edges))


def _estimate(recipe, template, n, extent_mode='pores'):
    """Count actual fragment atoms and bounded faces before allocating a molecular graph."""
    components = recipe['components']
    counts = {k:Chem.MolFromSmiles(c['fragment_smiles']).GetNumAtoms() for k,c in components.items()}
    if template in BRANCHED:
        from .cofdraw_branched import graph_for,loop_edges
        faces,_,_,kinds=graph_for(components,n,template)
        sites={v for face in faces for v in face}
        edges=set.union(*(loop_edges(face) for face in faces))
        return sum(counts[kinds[v]] for v in sites)-2*len(edges),len(faces)
    if template in (*SHEET_TEMPLATES,'self_square'):
        kinds,_,_,_,faces,nodes = _sheet_layout(recipe,template,n,extent_mode)
        shared = {frozenset((a,b)) for face in faces for a,b in zip(face,face[1:]+face[:1])}
        return sum(counts[kinds[k]] for _,_,k in nodes) - 2*len(shared),len(faces)
    if template in ('self_chain','chain_2_2'):
        return sum(counts[i['component']] for i in recipe['instances']) * n - 2*(n*len(recipe['instances'])-1),0
    if template=='self_kagome':
        from .cofdraw_self_kagome import HEX,patch_points
        nodes = patch_points(n, triad=extent_mode=='triad',star=extent_mode=='star')
        joins = sum((x+dx,y+dy) in nodes for x,y in nodes for dx,dy in HEX[:3])
        return len(nodes)*counts['A'] - 2*joins, joins-len(nodes)+1
    from .cofdraw_mixed import _faces as mixed_faces
    binary_kgm = template=='kgm_4_2'
    pattern = ('BBB' if extent_mode=='triangle' else 'BCBCBC') if binary_kgm else (recipe['assembly'].get('linker_sequence') or 'BBBBBB')
    faces = mixed_faces(pattern,n,template=='mixed_kgm' or binary_kgm,triad=extent_mode=='triad',
                        single_linker='A' if binary_kgm else None)
    edge_types = {tuple(sorted((p,vertices[(i+1)%len(vertices)]))):sequence[i]
                  for vertices,sequence in faces for i,p in enumerate(vertices)}
    nodes = {p for edge in edge_types for p in edge}
    return len(nodes)*counts['B' if binary_kgm else 'A'] + sum(counts[k] for k in edge_types.values()) - 4*len(edge_types),len(faces)


def repeat_options(recipe, source_rule=None, source_sha256=None):
    """Offer connected extents that fit the current fragment geometry and atom budget."""
    label = 'Source drawing' if source_sha256 else ('One segment' if recipe.get('layout_kind')=='open_chain' else 'One pore')
    options = {'available':False,'extents':[{'value':1,'label':label}],'reason':None,
               'scope':'Extent 1 retains the current drawing and coordinates.'}
    try:
        template = _classification(recipe,source_rule,source_sha256)
        model = recipe
        if source_sha256:
            model,_ = validated_design(recipe,source_rule=source_rule,source_sha256=source_sha256)
        if template in (*SHEET_TEMPLATES,'self_square'):
            _sheet_placement(model,template)
        chain = template in ('self_chain','chain_2_2')
        if template in BRANCHED:
            options['extents'][0]['label']='One circuit'
            for n in (2,3):
                atoms,_=_estimate(model,template,n)
                if atoms<=MAX_ATOMS:
                    if template in SPATIAL:
                        from .cofdraw_branched import branch_patch
                        try:branch_patch(model['components'],n,template)
                        except ValueError as error:
                            options['reason']=str(error)
                            continue
                    options['extents'].append(dict(value=n,label=f'{n} connected circuits',pores=n,atoms=atoms))
            options.update(template=template,available=len(options['extents'])>1,scope='Unfolded circuits from the chosen 3D parent net')
            return options
        if template=='kgm_4_2':
            options['extents'][0]['label']='One hexagon'
            atoms,_=_estimate(model,template,1,'triangle')
            if atoms<=MAX_ATOMS:
                options['extents'].append(dict(value='triangle',repeats=1,extent_mode='triangle',label='One triangle',pores=1,atoms=atoms))
        if not chain:
            atoms,pores = _estimate(model,template,3,'triad')
            if pores==3 and atoms<=MAX_ATOMS:
                options['extents'].append({'value':'three-pores','repeats':3,'extent_mode':'triad',
                                           'label':'3 pores','pores':3,'atoms':atoms})
        if template=='self_kagome':
            atoms,pores=_estimate(model,template,2,'star')
            if pores==7 and atoms<=MAX_ATOMS:
                options['extents'].append(dict(value='kagome-star',repeats=2,extent_mode='star',
                                               label='Hexagon + 6 triangles',pores=pores,atoms=atoms))
        for n in (2,3,4):
            count,pores = _estimate(model,template,n)
            if count > MAX_ATOMS:
                continue
            name = ('Hexagon + 6 triangles' if template in ('mixed_kgm','kgm_4_2') and n==2 else
                    f'{n} chain segments' if chain else f'{pores} pores')
            options['extents'].append({'value':n,'label':name,'pores':pores,'atoms':count})
        options.update(available=len(options['extents'])>1,template=template,
                       scope='Finite chemical design with shared building blocks. Larger extents use a new layout.')
        if not options['available']:
            options['reason']='Even the smallest expanded patch exceeds the 2,400-atom drawing limit.'
    except (ValueError,KeyError,TypeError,RuntimeError) as exc:
        options['reason']=str(exc)
    return options


def repeat_patch(recipe,repeats,*,source_rule=None,source_sha256=None,extent_mode='pores'):
    sizes={'pores':(2,3,4),'cells':(2,3,4),'triad':(3,),'star':(2,),'triangle':(1,)}
    if repeats not in sizes.get(extent_mode,()):
        raise ValueError('Choose an extent available for this drawing.')
    template = _classification(recipe,source_rule,source_sha256)
    if ((extent_mode=='star' and template!='self_kagome') or (extent_mode=='triangle' and template!='kgm_4_2')):
        raise ValueError('Choose a Kagome design for this pore arrangement.')
    assembly = recipe.get('assembly') or {}
    design = recipe
    if not assembly:
        design,_ = validated_design(recipe,source_rule=source_rule,source_sha256=source_sha256)
        if not source_sha256:
            for key,c in design['components'].items():
                c['label'] = recipe['components'][key].get('label',key)
    if extent_mode == 'cells' and assembly:
        raise ValueError('Assembly repeats use complete pore or chain neighborhoods.')
    if template in ('self_chain','chain_2_2'):
        if extent_mode=='triad':
            raise ValueError('Open chains can be extended by segments. This drawing has no closed pores.')
        if _estimate(recipe,template,repeats)[0] > MAX_ATOMS:
            raise ValueError('This chain exceeds the 2,400-atom drawing limit.')
        mol,details = _chain(recipe,repeats)
    elif template in BRANCHED:
        from .cofdraw_branched import branch_patch
        if extent_mode!='pores' or repeats not in (2,3) or _estimate(design,template,repeats)[0]>MAX_ATOMS:
            raise ValueError('Choose an available two- or three-circuit extent.')
        mol,details=branch_patch(design['components'],repeats,template)
    elif template == 'self_kagome':
        from .cofdraw_self_kagome import kagome_patch
        mol,details = kagome_patch(recipe,repeats,triad=extent_mode=='triad',star=extent_mode=='star')
    elif template in ('mixed_linkers','mixed_kgm','kgm_4_2','self_boroxine','self_triazine'):
        from .cofdraw_mixed import pore_patch
        pattern = ('BBB' if extent_mode=='triangle' else 'BCBCBC') if template=='kgm_4_2' else (
            'BBBBBB' if template in ('self_boroxine','self_triazine') else assembly['linker_sequence'])
        if not pattern:
            raise ValueError('The sequence around the selected aperture must be specified.')
        mol,details = pore_patch(design['components'],pattern,repeats,design['connections'][0]['bond_order'],
                                 layout='kgm' if template in ('mixed_kgm','kgm_4_2') else '',triad=extent_mode=='triad')
        if template=='kgm_4_2':
            details['face_labels']=['triangle' if len(s)==3 else 'hexagon' for s in details['face_labels']]
            details.pop('pore_sequences')
            details.update(pore_counts=dict(Counter(details['face_labels'])),
                           order_assumption='A single imine linker family joins the ETTA nodes of the selected Kagome net.')
    else:
        mol,details = _sheet(design,template,repeats,extent_mode)
    if mol.GetNumAtoms() > MAX_ATOMS or len(Chem.GetMolFrags(mol)) != 1:
        raise ValueError('The finite patch must be connected and fit within the 2,400-atom drawing limit.')
    atoms = mol.GetAtoms()
    origins = details.pop('atom_origin',None)
    if origins is None:
        origins = [dict(instance=a.GetProp('instance'),component=a.GetProp('component'),
                        fragment_atom_index=a.GetIntProp('fragment_atom_index')) for a in atoms]
    if len(origins) != mol.GetNumAtoms():
        raise ValueError('Building-block identity is missing for a repeated atom.')
    boundaries = [dict(atom_index=a.GetIdx(),**origins[a.GetIdx()]) for a in atoms if a.GetAtomicNum()==0]
    links = [b for b in mol.GetBonds() if origins[b.GetBeginAtomIdx()]['instance'] != origins[b.GetEndAtomIdx()]['instance']]
    fragcount = len(set(o['instance'] for o in origins))
    nfaces = len(links) - fragcount + 1
    expected={'triad':3,'star':7,'triangle':1}.get(extent_mode)
    if expected and nfaces!=expected:
        raise ValueError('The generated pores do not match the selected neighborhood.')
    if nfaces < 0 or details.get('pore_count',details.get('pores',nfaces)) != nfaces:
        raise ValueError('The counted pores disagree with the independent block-graph cycle count.')
    if template in ('self_boroxine','self_triazine'):
        details.pop('pore_sequences',None)
        details.update(pore_counts={'hexagon':nfaces},
                       order_assumption='One precursor supplies both junction rings and arylene spacers.')
    if nfaces and not details.get('pore_faces'):
        faces = _bounded_faces(mol,origins,links)
        if len(faces)==nfaces:
            details['pore_faces'] = faces
    if assembly and template in ('self_square','self_chain','self_kagome'):
        if any({b.GetBeginAtom().GetAtomicNum(),b.GetEndAtom().GetAtomicNum()} != {6,7}
               or b.GetBondTypeAsDouble() != 2 for b in links):
            raise ValueError('Every aminoaldehyde junction must join an aldehyde carbon to an amine nitrogen.')
    extent = (f'{repeats} chain segments' if template in ('self_chain','chain_2_2') else
              f'{repeats} × {repeats} cells' if extent_mode=='cells' else 'One triangle' if extent_mode=='triangle' else f'{nfaces} connected pores')
    basis = design.get('source_drawing_basis')
    if basis:
        basis = {k:basis[k] for k in ('recipe_id','recipe_sha256','drawing_doi','template','orientation_assumption') if k in basis}
        basis['assumption'] = 'The reviewed precursors and local pore determine the chemistry. Neighboring pores follow the selected 2D net.'
    summary = {**details,'extent':extent,'repeats':repeats,'template':template,'extent_mode':extent_mode,
               'pore_count':nfaces,'blocks':fragcount,'interblock_bonds':len(links),
               'atoms':mol.GetNumAtoms(),'continuation_sites':len(boundaries),'component_counts':dict(Counter(
                   {o['instance']:o['component'] for o in origins}.values())),
               'source_drawing_basis':basis,
               'scope':'Finite connected drawing with explicit covalent joins and external continuation ports. The neighboring arrangement is a design hypothesis.'}
    return mol, origins, boundaries, summary, design['components']


def _focal_aperture(mol,origins,faces,labels,summary,recipe):
    """Pick the central occurrence of the selected pore family in the displayed patch."""
    from .cofdraw_measure import prepare_aperture
    xy=mol.GetConformer().GetPositions()[:,:2]
    owners=[o['instance'] for o in origins]
    centers={name:xy[[i for i,owner in enumerate(owners) if owner==name]].mean(axis=0) for name in set(owners)}
    positions=np.array([np.mean([centers[name] for name in rim],axis=0) for rim in faces])
    middle=positions.mean(axis=0)
    assembly=recipe.get('assembly') or {}
    preferred=assembly.get('pore_view') or assembly.get('linker_sequence')
    if summary['extent_mode']=='star' or (summary['pore_count']==7 and summary['template'] in ('kgm_4_2','mixed_kgm')):
        preferred='BCBCBC' if summary['template']=='mixed_kgm' else 'hexagon'
    index=min(range(len(faces)),key=lambda i:(bool(preferred and labels[i]!=preferred),float(np.linalg.norm(positions[i]-middle)),i))
    kind=labels[index]
    name=({'BBB':'B triangle','CCC':'C triangle','BCBCBC':'Alternating B/C hexagon'}.get(kind)
          or (kind.capitalize() if not set(kind)<=set('BC') else f'{len(kind)}-sided pore · {kind}'))
    aperture=prepare_aperture(mol,owners,faces[index])
    return {**aperture,'label':name,'type':kind,'instances':faces[index],'index':index}


def write_repeat(recipe, directory, key, request, *, source_rule=None,source_sha256=None):
    from .cofdraw_core import export_chemical_formats, render_styled
    mol,origins,boundaries,summary,components = repeat_patch(recipe,request.repeats,source_rule=source_rule,
                                    source_sha256=source_sha256,extent_mode=request.extent_mode)
    faces = summary.pop('pore_faces',[])
    labels = summary.pop('face_labels',[])
    if faces and not labels:
        template = summary['template']
        if template in ('mixed_linkers','mixed_kgm'):
            from .cofdraw_mixed import canonical_sequence
            owner = {o['instance']:o['component'] for o in origins}
            labels = [canonical_sequence(''.join(owner[k] for k in face if owner[k] != 'A')) for face in faces]
        elif template=='self_kagome':
            labels = ['triangle' if len(face)==3 else 'hexagon' for face in faces]
        else:
            shape = ('triangle' if template=='hxl_6_2' else 'quadrilateral' if template in
                     ('sql_4_2','sql_4_4','self_square') else 'hexagon')
            labels = [shape]*len(faces)
    if summary['template'] in ('self_triazine','self_boroxine'):
        labels = ['hexagon']*len(faces)
    if len(labels)!=len(faces):
        raise ValueError('A pore environment could not be matched to its drawing boundary.')
    # Keep colors stable when changing extent, or selecting another example
    # with the same pore environments. For additional mixed sequences assign
    # distinct colors within the current drawing.
    preferred = {'BBB':'#009E73','CCC':'#E69F00','BCBCBC':'#8064B0',
                 'triangle':'#009E73','hexagon':'#8064B0','quadrilateral':'#0072B2'}
    colors = {k:preferred[k] for k in sorted(set(labels)) if k in preferred}
    extras = (c for c in ('#0072B2','#D55E00','#CC79A7','#56B4E9','#117733','#AA4499',
              '#882255','#999933','#44AA99','#DDCC77','#332288') if c not in colors.values())
    for label in sorted(set(labels)-colors.keys()):
        colors[label]=next(extras)
    summary['pore_colors'] = colors if request.color_pores else {}
    drawing = dict(recipe_id=key,components=components,layout_kind=recipe.get('layout_kind'),repeat_faces=faces,
                   repeat_face_colors=[colors[k] for k in labels] if request.color_pores else [],
                   designed_variant=True,ground_truth=False)
    if source_sha256:
        style = recipe.get('source_style') or {}
        # Atom indices and instance names from the witnessed drawing do not
        # address this new molecule. Carry only portable, component-wide rules.
        drawing['source_style'] = {k:style[k] for k in ('halo','edge_style','abbreviation_mode') if k in style}
        for field in ('component_colors','component_fills'):
            if field in style:
                drawing['source_style'][field] = {new:style[field][old] for new,old in
                    zip(('A','B'),source_rule['precursor_components']) if old in style[field]}
    smiles = Chem.MolToSmiles(mol,canonical=False,isomericSmiles=True)
    checksum = hashlib.sha256(json.dumps(drawing,sort_keys=True,separators=(',',':')).encode()).hexdigest()
    audit = dict(recipe_sha256=checksum,atom_origin=origins,continuation_ports=boundaries,
                 fragment_smiles=smiles,fragment_smiles_canonical=False,whole_cof_smiles=None,
                 fragment_serialization_sha256=hashlib.sha256(smiles.encode()).hexdigest(),
                 graph_scope='finite_designed_connected_repeat',ground_truth=False,source_registered=False,
                 atoms=mol.GetNumAtoms(),bonds=mol.GetNumBonds(),defined_double_bond_geometry=_double_bond_geometry_audit(mol))
    audit = render_styled(drawing,directory,request.size,request.style_id,request.style_seed,request.style_options,
                          prepared_patch=(mol,audit),recent_styles=request.recent_styles)
    if faces:
        try:
            aperture=_focal_aperture(mol,origins,faces,labels,summary,recipe)
            audit['drawing_aperture']=aperture
            summary['measurement_face']={k:aperture[k] for k in ('label','type','index','instances')}
        except ValueError as exc:
            summary['measurement_reason']=str(exc)
    extra = export_chemical_formats(directory,key,source_mol=mol,repeated_patch=True)
    chemical = json.loads((directory/f'{key}.chemical-export.json').read_text(encoding='utf8'))
    audit['chemical_exports'] = chemical
    if 'cdxml' not in extra:
        summary['cdxml_note'] = chemical.get('reason') or 'CDXML is unavailable in this RDKit build.'
    urls = {ext:f'/api/drawing/exports/{key}/{ext}' for ext in ('svg','png','mol','json',*extra)}
    (directory/f'{key}.json').write_text(json.dumps({'summary':summary,'audit':audit,'input_recipe':recipe,
        'drawing_id':key,'downloads':urls},ensure_ascii=False,indent=2),encoding='utf8')
    return {'summary':summary,'svg':urls['svg'],'drawing_id':key,'downloads':urls}


def repeated_drawing(recipe, repeats, size):
    """Compatibility for existing Python callers that requested a crop of primitive cells."""
    mol,_,_,summary,_ = repeat_patch(recipe,repeats,extent_mode='cells')
    drawer = rdMolDraw2D.MolDraw2DSVG(size,size)
    drawer.drawOptions().padding = .03
    rdMolDraw2D.PrepareAndDrawMolecule(drawer,mol)
    drawer.FinishDrawing()
    return drawer.GetDrawingText(),summary
