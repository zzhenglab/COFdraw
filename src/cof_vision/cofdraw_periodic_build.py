"""Assemble a periodic starting cell from reacted planar building blocks."""
from __future__ import annotations

from itertools import product
import math

import numpy as np
from rdkit import Chem

from .cofdraw_depiction import compute_2d_coords, free_base_porphyrin_core
from .cofdraw_network import inspect_periodic_network
from .cofdraw_science import normalize_structure
from .framework_drawing import _ports


SUPPORTED_TEMPLATES = ('hcb_3_3', 'hcb_3_2', 'sql_4_2', 'sql_4_4')
MAX_CELL_SITES = 360
ORGANIC_ELEMENTS = {1, 5, 6, 7, 8, 9}


def _net(template, counts):
    """Return block types, an ideal net for port directions, and quotient edges."""
    if template in ('hcb_3_3', 'hcb_3_2'):
        cell = np.array([[1.5, math.sqrt(3)/2], [0., math.sqrt(3)]])
        directions = [(0, 0), (-1, 0), (-1, 1)]
        partner = np.array([1., 0.])
        if template == 'hcb_3_3':
            return ['A', 'B'], cell, [np.zeros(2), partner], [(0, 1, t) for t in directions]
        node = max(counts, key=counts.get)
        linker = 'B' if node == 'A' else 'A'
        midpoints = [(partner + np.array(t) @ cell)/2 for t in directions]
        edges = [(0, i, (0, 0)) for i in range(2, 5)]
        edges += [(i, 1, t) for i, t in enumerate(directions, 2)]
        return [node, node, linker, linker, linker], cell, [np.zeros(2), partner, *midpoints], edges
    if template == 'sql_4_2':
        node = max(counts, key=counts.get)
        linker = 'B' if node == 'A' else 'A'
        edges = [(0, 1, (0, 0)), (1, 0, (1, 0)), (0, 2, (0, 0)), (2, 0, (0, 1))]
        return [node, linker, linker], np.eye(2), [[0., 0.], [.5, 0.], [0., .5]], edges
    if template == 'sql_4_4':
        cell = np.array([[1., 1.], [-1., 1.]])
        edges = [(0, 1, t) for t in ((0, 0), (-1, 0), (-1, -1), (0, -1))]
        return ['A', 'B'], cell, [[0., 0.], [0., 1.]], edges
    raise ValueError('Periodic construction supports selected hexagonal and square planar templates.')


def _tilt_phenyl_arms(mol, core):
    """Tilt para-substituted meso phenyl rings about their unchanged axes."""
    if not core:
        return 0
    conformer = mol.GetConformer()
    xyz = conformer.GetPositions()
    angle = math.radians(60)
    count = 0
    for ring in mol.GetRingInfo().AtomRings():
        if len(ring) != 6 or any(not mol.GetAtomWithIdx(i).GetIsAromatic()
                                 or mol.GetAtomWithIdx(i).GetAtomicNum() != 6 for i in ring):
            continue
        ring_set = set(ring)
        external_heavy = [(i, a.GetIdx()) for i in ring for a in mol.GetAtomWithIdx(i).GetNeighbors()
                          if a.GetIdx() not in ring_set and a.GetAtomicNum() != 1]
        meso = [i for i, a in external_heavy if a in core]
        if len(meso) != 1:
            continue
        first, opposite = meso[0], ring[(ring.index(meso[0])+3) % 6]
        if any(i not in (first, opposite) for i, _ in external_heavy):
            continue
        axis = xyz[opposite] - xyz[first]
        axis /= np.linalg.norm(axis)
        move = [i for i in ring if i not in (first, opposite)]
        move += [a.GetIdx() for i in move for a in mol.GetAtomWithIdx(i).GetNeighbors() if a.GetAtomicNum() == 1]
        for i in move:
            vector = xyz[i]-xyz[first]
            rotated = (vector*math.cos(angle) + np.cross(axis, vector)*math.sin(angle)
                       + axis*(axis @ vector)*(1-math.cos(angle)))
            conformer.SetAtomPosition(i, tuple(xyz[first]+rotated))
        count += 1
    return count


def _block(component):
    mol = Chem.MolFromSmiles(component['fragment_smiles'])
    if mol is None or len(Chem.GetMolFrags(mol)) != 1:
        raise ValueError('Each reacted building block must be one connected molecule.')
    if (any(a.GetAtomicNum() not in ORGANIC_ELEMENTS | {0} or a.GetIsotope()
            or a.GetNumRadicalElectrons() for a in mol.GetAtoms())
            or Chem.GetFormalCharge(mol)):
        raise ValueError('Automatic cells require neutral organic fragments containing C, H, B, N, O or F, without isotope labels or radicals.')
    if any(a.GetChiralTag() != Chem.ChiralType.CHI_UNSPECIFIED for a in mol.GetAtoms()) or any(
            b.GetStereo() not in (Chem.BondStereo.STEREONONE, Chem.BondStereo.STEREOANY) for b in mol.GetBonds()):
        raise ValueError('Building blocks with specified stereochemistry need explicit 3D coordinates.')
    core = free_base_porphyrin_core(mol)
    ports = _ports(mol)
    compute_2d_coords(mol, canonOrient=False, forceRDKit=True)
    conformer = mol.GetConformer()
    xy = conformer.GetPositions()[:, :2]
    xy = (xy - xy[[a for _, a in ports.values()]].mean(axis=0)) * (1.4/1.5)
    for i, (x, y) in enumerate(xy):
        conformer.SetAtomPosition(i, (float(x), float(y), 0.))
    # Port bonds remain in place during hydrogen addition, so a C=N port has
    # the same valence and H count as the actual inter-block C=N bond.
    mol = Chem.AddHs(mol, addCoords=True)
    tilted_rings = _tilt_phenyl_arms(mol, set(core))
    indices = [a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum()]
    index = {i: j for j, i in enumerate(indices)}
    atoms = [mol.GetAtomWithIdx(i) for i in indices]
    xyz = mol.GetConformer().GetPositions()[indices]
    sites, directions = [], []
    for dummy, endpoint in ports.values():
        vector = xy[dummy] - xy[endpoint]
        sites.append(index[endpoint])
        directions.append(vector/np.linalg.norm(vector))
    bonds = [(index[a], index[b], bond.GetBondTypeAsDouble()) for bond in mol.GetBonds()
             if (a := bond.GetBeginAtomIdx()) in index and (b := bond.GetEndAtomIdx()) in index]
    return {'xyz': xyz, 'sites': np.array(sites), 'directions': np.array(directions),
            'elements': [a.GetSymbol() for a in atoms], 'bonds': bonds,
            'charges': [a.GetFormalCharge() for a in atoms], 'tilted_rings': tilted_rings}


def _rotate(vectors, angles):
    cos, sin = np.cos(angles), np.sin(angles)
    x, y = np.moveaxis(vectors, -1, 0)
    return np.stack((cos*x-sin*y, sin*x+cos*y), axis=-1)


def _alignment(radial, directions):
    """Match cyclically ordered ports to the incident net directions."""
    angles = np.arctan2(radial[:, 1], radial[:, 0])
    targets = np.arctan2(directions[:, 1], directions[:, 0])
    order = np.argsort(angles)
    desired = np.argsort(targets)
    candidates = []
    for offset in range(len(order)):
        mapping = np.empty_like(order)
        mapping[desired] = np.roll(order, offset)
        gap = targets-angles[mapping]
        angle = float(np.arctan2(np.sin(gap).sum(), np.cos(gap).sum()))
        error = float(np.max(np.abs(np.angle(np.exp(1j*(gap-angle))))))
        candidates.append((round(error, 7), abs(angle), mapping, angle))
    error, _, mapping, angle = min(candidates, key=lambda item: item[:2])
    if error > math.radians(25):
        raise ValueError('The reactive sites do not point along the chosen planar net.')
    return mapping, angle


def _fit_blocks(blocks, net_cell, positions, edges, join_length):
    """Fit rotations and solve the linear periodic translations for each pose."""
    n = len(blocks)
    incidence = np.zeros((len(edges), n+1))
    incident = [[] for _ in blocks]
    edge_slots = np.empty((len(edges), 2), dtype=int)
    for i, (a, b, shift) in enumerate(edges):
        if a:
            incidence[i, a-1] -= 1
        if b:
            incidence[i, b-1] += 1
        incidence[i, -2:] = shift
        vector = np.asarray(positions[b]) + np.array(shift) @ net_cell - positions[a]
        for j, owner, direction in ((0, a, vector), (1, b, -vector)):
            edge_slots[i, j] = len(incident[owner])
            incident[owner].append(direction/np.linalg.norm(direction))
    solve_translations = np.linalg.pinv(incidence)
    owners = np.array([(a, b) for a, b, _ in edges])
    axes = [np.asarray(values) for values in incident]
    candidates = []
    for flips in product((1, -1), repeat=n):
        mirror = [np.array([1, sign]) for sign in flips]
        try:
            placements = [_alignment(block['xyz'][block['sites'], :2]*sign, target)
                          for block, sign, target in zip(blocks, mirror, axes)]
        except ValueError:
            continue
        initial = np.array([angle for _, angle in placements])
        site_ids = np.array([[placements[owner][0][slot] for owner, slot in zip(pair, slots)]
                             for pair, slots in zip(owners, edge_slots)])
        endpoints = np.array([[blocks[k]['xyz'][blocks[k]['sites'][port], :2]*mirror[k]
                               for k, port in zip(pair, ids)] for pair, ids in zip(owners, site_ids)])
        tangents = np.array([[blocks[k]['directions'][port]*mirror[k]
                              for k, port in zip(pair, ids)] for pair, ids in zip(owners, site_ids)])

        def evaluate(angles):
            anchors = _rotate(endpoints, angles[owners])
            vectors = _rotate(tangents, angles[owners])
            target = anchors[:, 0] - anchors[:, 1] + .5*join_length*(vectors[:, 0]-vectors[:, 1])
            fitted = solve_translations @ target
            residual = np.r_[((incidence @ fitted) - target).ravel(),
                             (.65*join_length*vectors.sum(axis=1)).ravel(), .12*(angles-initial)]
            return residual, fitted, anchors, vectors

        angles = initial.copy()
        for _ in range(25):
            residual, *_ = evaluate(angles)
            step = np.eye(n)*1e-5
            jacobian = np.column_stack([(evaluate(angles+delta)[0]-residual)/1e-5 for delta in step])
            delta = np.linalg.lstsq(jacobian, -residual, rcond=None)[0]
            delta *= min(1., .25/max(float(np.max(abs(delta))), 1e-12))
            score = float(residual @ residual)
            for rate in (1., .5, .25, .1):
                trial = angles + delta*rate
                new_residual, *_ = evaluate(trial)
                if np.max(abs(trial-initial)) <= .75 and new_residual @ new_residual < score-1e-9:
                    angles = trial
                    break
            else:
                break
        residual, fitted, anchors, vectors = evaluate(angles)
        cell = fitted[-2:]
        displacements = incidence @ fitted + anchors[:, 1]-anchors[:, 0]
        lengths = np.linalg.norm(displacements, axis=1)
        if (np.linalg.det(cell) <= 0 or np.linalg.cond(cell) > 3.5 or np.min(np.linalg.norm(cell, axis=1)) < 4
                or np.max(abs(lengths-join_length)) > .14 or np.min(lengths) < 1.15):
            continue
        cosine = np.einsum('ij,ij->i', displacements/lengths[:, None], vectors[:, 0])
        other_cosine = -np.einsum('ij,ij->i', displacements/lengths[:, None], vectors[:, 1])
        minimum_cosine = float(min(cosine.min(), other_cosine.min()))
        if minimum_cosine < math.cos(math.radians(24)):
            continue
        candidate = (float(residual @ residual), cell, np.vstack([np.zeros(2), fitted[:-2]]),
                     angles, mirror, site_ids, lengths, math.degrees(math.acos(min(1., minimum_cosine))))
        candidates.append(candidate)
    return sorted(candidates, key=lambda item: item[0])


def edge_key(a, b, shift):
    """Canonical atom pair and image displacement."""
    t = tuple(map(int, shift))
    return (a, b, t) if a < b else (b, a, tuple(-v for v in t))


def _cell_atoms(blocks, edges, pose, spacing, order):
    _, cell, centers, angles, mirror, ports, lengths, angle_error = pose
    basis = np.zeros((3, 3))
    basis[:2, :2] = cell
    basis[2, 2] = spacing
    inverse = np.linalg.inv(basis)
    atoms, bonds, wraps, offsets, formal_charges = [], [], [], [], []
    for i, block in enumerate(blocks):
        offsets.append(len(atoms))
        xyz = block['xyz'].copy()
        xyz[:, :2] = _rotate(xyz[:, :2]*mirror[i], angles[i]) + centers[i]
        xyz[:, 2] += spacing/2
        fraction = xyz @ inverse
        wrap = np.floor(fraction + 1e-10).astype(int)
        wraps.append(wrap)
        atoms.extend({'element': element, 'fractional': p.tolist()}
                     for element, p in zip(block['elements'], fraction-wrap))
        charges = [{'atom': offsets[i]+j, 'formal_charge': charge}
                   for j, charge in enumerate(block['charges']) if charge]
        formal_charges.extend(charges)
        bonds.extend((offsets[i]+a, offsets[i]+b, wrap[b]-wrap[a], kind)
                     for a, b, kind in block['bonds'])
    for (i, j, shift), (p, q) in zip(edges, ports):
        a, b = blocks[i]['sites'][p], blocks[j]['sites'][q]
        translation = np.array([*shift, 0])+wraps[j][b]-wraps[i][a]
        bonds.append((offsets[i]+int(a), offsets[j]+int(b), translation, order))
    return {'periodic': True, 'cell': basis.tolist(), 'atoms': atoms}, bonds, formal_charges, lengths, angle_error


def build_periodic_layer(design, *, spacing=None):
    """Build and check a single AA layer cell using the validated design."""
    template = design['selected_template']['id']
    if template not in SUPPORTED_TEMPLATES:
        raise ValueError('The chosen connectivity has no automatic periodic cell builder.')
    components = design['components']
    block_types, reference, positions, edges = _net(template, {k: c['connectivity'] for k, c in components.items()})
    unique = {k: _block(c) for k, c in components.items()}
    blocks = [unique[k] for k in block_types]
    if sum(len(b['elements']) for b in blocks) > MAX_CELL_SITES:
        raise ValueError(f'Automatic cells are limited to {MAX_CELL_SITES} explicit sites.')
    order = design['connections'][0]['bond_order']
    poses = _fit_blocks(blocks, reference, positions, edges, 1.29 if order == 2 else 1.42)
    if not poses:
        raise ValueError('The selected building blocks could not close this cell with consistent bond lengths and directions.')
    tilted_rings = sum(block['tilted_rings'] for block in blocks)
    spacings = [spacing] if spacing is not None else [5., 5.8] if tilted_rings else [3.4]
    for height in spacings:
        for pose in poses[:8]:
            structure, bonds, charges, lengths, angle_error = _cell_atoms(blocks, edges, pose, height, order)
            try:
                structure = normalize_structure(structure)
                network = inspect_periodic_network(structure)
            except ValueError:
                continue
            expected = {edge_key(a, b, t) for a, b, t, _ in bonds}
            observed = {edge_key(e['a'], e['b'], e['translation']) for e in network['edges']}
            if (network['short_contact_count'] or network['dimensions'] != [2]
                    or network['quotient_component_count'] != 1 or network['threshold_sensitivity']['classification_changed']
                    or expected != observed or len(expected) != len(bonds)):
                continue
            structure['provenance'] = {
                'kind': 'cofdraw_periodic_starting_model', 'template': template, 'linkage': design['linkage'],
                'precursors': [c['monomer_smiles'] for c in components.values()],
                'geometry_method': 'Rigid reacted fragments fitted to the chosen periodic connection graph',
                'preserves_drawing_coordinates': False, 'stacking': 'AA', 'spacing_angstrom': height,
                'phenyl_tilt_degrees': 60 if tilted_rings else None, 'tilted_phenyl_ring_count': tilted_rings,
                'qualification': 'Unrelaxed periodic starting model with explicit hydrogens and chosen layer spacing.',
            }
            periodic_bonds = [{'a': int(a), 'b': int(b), 'translation': list(map(int, t)), 'order': kind}
                              for a, b, t, kind in bonds]
            return structure, {
                'periodic_bonds': periodic_bonds, 'formal_charges': charges,
                'periodic_rank': 2, 'expected_atom_edges_verified': len(bonds), 'short_contacts': 0,
                'junction_length_range_angstrom': [float(lengths.min()), float(lengths.max())],
                'maximum_port_direction_error_degrees': angle_error,
            }
    raise ValueError('The fitted cell has steric contacts or does not reproduce the intended periodic bonds. Try a validated CIF.')
