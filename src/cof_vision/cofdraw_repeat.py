"""Finite 2D supercells from declared template edges and reacted molecular ports."""
from copy import deepcopy
from itertools import product

import numpy as np
from rdkit import Chem
from rdkit.Chem.Draw import rdMolDraw2D

from .cofdraw_depiction import compute_2d_coords
from .cofdraw_periodic_build import SUPPORTED_TEMPLATES, _fit_blocks, _net, _rotate
from .cofdraw_periodic_options import validated_design
from .framework_drawing import _ports, build_patch


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


def _chain(recipe, repeats):
    extended = deepcopy(recipe)
    sequence = recipe['instances']
    count = len(sequence) * repeats
    if count * max(Chem.MolFromSmiles(c['fragment_smiles']).GetNumAtoms() for c in recipe['components'].values()) > 2400:
        raise ValueError('This repeat exceeds the 2,400-atom drawing limit.')
    items = [{**sequence[i % len(sequence)], 'id': str(i), 'layout_position': [12*i, 0]} for i in range(count)]
    order = recipe['connections'][0]['bond_order']
    extended.update(instances=items, atom_coordinates=None, connections=[
        {'a': {'instance': a['id'], 'port': a['outgoing_port']},
         'b': {'instance': b['id'], 'port': b['incoming_port']}, 'bond_order': order}
        for a, b in zip(items, items[1:])])
    return build_patch(extended)[0], count


def _sheet(design, template, repeats):
    components = design['components']
    kinds, reference, centers, edges = _net(template, {k:c['connectivity'] for k,c in components.items()})
    blocks = {key:_planar_block(c) for key,c in components.items()}
    cell_blocks = [blocks[key] for key in kinds]
    if sum(b['mol'].GetNumAtoms() for b in cell_blocks)*repeats**2 > 2400:
        raise ValueError('This repeat exceeds the 2,400-atom drawing limit. Choose a smaller extent.')
    order = design['connections'][0]['bond_order']
    poses = _fit_blocks(cell_blocks, reference, centers, edges, 1.35 if order == 2 else 1.5)
    if not poses:
        raise ValueError('These fragments cannot fit the chosen net with consistent bond lengths. The original patch is available.')
    _, lattice, centers, angles, mirror, slots, *_ = poses[0]
    fragments = {}
    for u,v,k in product(range(repeats),range(repeats),range(len(kinds))):
        block = cell_blocks[k]
        mol = Chem.Mol(block['mol'])
        xy = _rotate(block['xyz'][:, :2]*mirror[k],angles[k]) + centers[k] + [u,v] @ lattice
        for i,(x,y) in enumerate(xy):
            mol.GetConformer().SetAtomPosition(i,(float(x),float(y),0.))
        for dummy,_ in block['ports']:
            mol.GetAtomWithIdx(dummy).SetAtomMapNum(0)
        fragments[u,v,k] = mol
    tag = 1
    for u,v in product(range(repeats),repeat=2):
        for (a,b,(du,dv)), (sa,sb) in zip(edges,slots):
            neighbor = (u+du,v+dv,b)
            if neighbor not in fragments:
                continue
            for key,site in [((u,v,a),sa),(neighbor,sb)]:
                dummy = cell_blocks[key[2]]['ports'][site][0]
                fragments[key].GetAtomWithIdx(dummy).SetAtomMapNum(tag)
            tag += 1
    parts = iter(fragments.values())
    joined = next(parts)
    for mol in parts:
        joined = Chem.CombineMols(joined,mol)
    result = Chem.molzip(joined)
    if len(Chem.GetMolFrags(result)) != 1 or result.GetNumAtoms() != joined.GetNumAtoms()-2*(tag-1):
        raise ValueError('Repeat connectivity could not be verified.')
    result.UpdatePropertyCache(strict=True)
    Chem.GetSSSR(result)
    return result, len(fragments)


def repeated_drawing(recipe, repeats, size):
    design,_ = validated_design(recipe)
    template = design['selected_template']['id']
    if template == 'chain_2_2':
        mol,blocks = _chain(recipe,repeats)
        extent = str(repeats)+' chain segments'
    elif template in SUPPORTED_TEMPLATES:
        mol,blocks = _sheet(design,template,repeats)
        extent = f'{repeats} × {repeats} cells'
    else:
        raise ValueError('Repeat drawings require a chain, hexagonal, or square template with declared translation edges.')
    drawer = rdMolDraw2D.MolDraw2DSVG(size,size)
    drawer.drawOptions().padding = .03
    rdMolDraw2D.PrepareAndDrawMolecule(drawer,mol)
    drawer.FinishDrawing()
    return drawer.GetDrawingText(), dict(extent=extent, template=template, blocks=blocks,
        atoms=mol.GetNumAtoms(), continuation_sites=sum(a.GetAtomicNum()==0 for a in mol.GetAtoms()))
