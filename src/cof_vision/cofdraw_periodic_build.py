"""Assemble a periodic starting cell from reacted planar building blocks."""
from __future__ import annotations

from itertools import product
import math

import numpy as np
from rdkit import Chem
from rdkit.Chem import rdDistGeom, rdForceFieldHelpers
from rdkit.Geometry import Point3D

from .cofdraw_depiction import compute_2d_coords, free_base_porphyrin_core, normalize_imine_ports
from .cofdraw_network import inspect_periodic_network
from .cofdraw_science import normalize_structure
from .framework_drawing import _ports


SUPPORTED_TEMPLATES = ('hcb_3_3', 'hcb_3_2', 'sql_4_2', 'sql_4_4')
MAX_CELL_SITES = 360
ORGANIC_ELEMENTS = {1, 5, 6, 7, 8, 9, 14, 15, 16, 17, 35, 53}


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


def _tilt_phenyl_arms(mol, core, general_tilt):
    """Rotate unfused para-phenylene rings about their two attachment atoms."""
    conformer = mol.GetConformer()
    xyz = conformer.GetPositions()
    count = 0
    for ring in mol.GetRingInfo().AtomRings():
        if len(ring) != 6 or any(not mol.GetAtomWithIdx(i).GetIsAromatic()
                                 or mol.GetAtomWithIdx(i).GetAtomicNum() != 6 for i in ring):
            continue
        ring_set = set(ring)
        external_heavy = [(i, a.GetIdx()) for i in ring for a in mol.GetAtomWithIdx(i).GetNeighbors()
                          if a.GetIdx() not in ring_set and a.GetAtomicNum() != 1]
        if len(external_heavy) != 2 or any(mol.GetRingInfo().NumAtomRings(i) != 1 for i in ring):
            continue
        meso = [i for i, a in external_heavy if a in core]
        first = meso[0] if meso else external_heavy[0][0]
        opposite = ring[(ring.index(first)+3) % 6]
        degrees = 60 if meso else general_tilt
        if not degrees or {i for i, _ in external_heavy} != {first, opposite}:
            continue
        angle = math.radians(degrees)
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


def _relieve_fragment(mol, ports, cap_element):
    """Minimize a capped fragment while fixing each junction and its tangent."""
    capped = Chem.RWMol(mol)
    fixed = set()
    for dummy, endpoint in ports.values():
        capped.ReplaceAtom(dummy, Chem.Atom(cap_element))
        fixed.update((dummy, endpoint))
    capped = capped.GetMol()
    Chem.SanitizeMol(capped)
    conf = capped.GetConformer()
    conf.Set3D(True)
    # Break the exact reflection symmetry of the drawing so a torsion can
    # move away from its planar saddle point during minimization.
    for i, (x, y, z) in enumerate(conf.GetPositions()):
        if i not in fixed:
            conf.SetAtomPosition(i, (x, y, z + .15*math.sin(2.4*i + .5)))
    capped = Chem.AddHs(capped, addCoords=True)
    props = rdForceFieldHelpers.MMFFGetMoleculeProperties(capped)
    if props is None:
        return False
    ff = rdForceFieldHelpers.MMFFGetMoleculeForceField(capped, props)
    if ff is None:
        return False
    for i in fixed:
        ff.AddFixedPoint(i)
    for atom in capped.GetAtoms():
        i = atom.GetIdx()
        if atom.GetAtomicNum() > 1 and i not in fixed:
            ff.MMFFAddPositionConstraint(i, .4, 3.)
    ff.Initialize()
    if ff.Minimize(maxIts=600) != 0:
        return False
    result = capped.GetConformer().GetPositions()[:mol.GetNumAtoms()]
    if not np.isfinite(result).all():
        return False
    conf = mol.GetConformer()
    for i, point in enumerate(result):
        conf.SetAtomPosition(i, tuple(point))
    return True


def _embedded_block(mol, ports, xy, seed, cap_element):
    """Embed a capped fragment at its junctions, preserving specified stereocentres."""
    capped = Chem.RWMol(mol)
    fixed = {}
    for dummy, endpoint in ports.values():
        capped.ReplaceAtom(dummy, Chem.Atom(cap_element))
        for i in (dummy, endpoint):
            fixed[i] = Point3D(*xy[i], 0.)
    capped = capped.GetMol()
    Chem.SanitizeMol(capped)
    capped = Chem.AddHs(capped)
    settings = rdDistGeom.ETKDGv3()
    settings.randomSeed = seed
    settings.useRandomCoords = True
    settings.maxIterations = 45
    settings.SetCoordMap(fixed)
    if rdDistGeom.EmbedMolecule(capped, settings) < 0:
        raise ValueError('Could not embed the specified stereochemistry at the planar junctions.')
    props = rdForceFieldHelpers.MMFFGetMoleculeProperties(capped)
    ff = rdForceFieldHelpers.MMFFGetMoleculeForceField(capped, props) if props else None
    if ff is None:
        raise ValueError('The constrained fragment needs MMFF parameters.')
    for i in fixed:
        ff.AddFixedPoint(i)
    ff.Initialize()
    if ff.Minimize(maxIts=1200):
        raise ValueError('Constrained 3D fragment relaxation did not converge.')
    observed = Chem.Mol(capped)
    Chem.AssignStereochemistryFrom3D(observed)
    for atom in capped.GetAtoms():
        tag = atom.GetChiralTag()
        if tag != Chem.ChiralType.CHI_UNSPECIFIED and observed.GetAtomWithIdx(atom.GetIdx()).GetChiralTag() != tag:
            raise ValueError('A generated fragment did not retain its assigned stereochemistry.')
    # AddHs appends cap hydrogens as well. Map only hydrogens of real atoms.
    product = Chem.AddHs(mol)
    conf = Chem.Conformer(product.GetNumAtoms())
    embedded = capped.GetConformer()
    def added_hydrogens(graph, i):
        return [a.GetIdx() for a in graph.GetAtomWithIdx(i).GetNeighbors()
                    if a.GetIdx() >= mol.GetNumAtoms() and a.GetAtomicNum() == 1]
    caps = {dummy for dummy, _ in ports.values()}
    for i in range(mol.GetNumAtoms()):
        conf.SetAtomPosition(i, embedded.GetAtomPosition(i))
        if i in caps:
            continue
        for dest, src in zip(added_hydrogens(product, i), added_hydrogens(capped, i), strict=True):
            conf.SetAtomPosition(dest, embedded.GetAtomPosition(src))
    product.RemoveAllConformers()
    product.AddConformer(conf)
    return product


def _block(component, general_tilt=0, cap_element=None, embed_seed=None, imine_mode=None):
    mol = Chem.MolFromSmiles(component['fragment_smiles'])
    if mol is None or len(Chem.GetMolFrags(mol)) != 1:
        raise ValueError('Each reacted building block must be one connected molecule.')
    if (any(a.GetAtomicNum() not in ORGANIC_ELEMENTS | {0}
            or a.GetNumRadicalElectrons() for a in mol.GetAtoms())
            or Chem.GetFormalCharge(mol)):
        raise ValueError('Automatic cells require neutral organic fragments containing H, B, C, N, O, F, Si, P, S, Cl, Br or I and no radicals.')
    if any(b.GetStereo() not in (Chem.BondStereo.STEREONONE, Chem.BondStereo.STEREOANY) for b in mol.GetBonds()):
        raise ValueError('Building blocks with specified bond stereochemistry need explicit 3D coordinates.')
    core = free_base_porphyrin_core(mol)
    ports = _ports(mol)
    compute_2d_coords(mol, canonOrient=False, forceRDKit=True)
    conformer = mol.GetConformer()
    xy = conformer.GetPositions()[:, :2]
    if imine_mode:
        normalize_imine_ports(mol, xy, preserve_uniform=True, straight_only=imine_mode == 'straight')
    xy = (xy - xy[[a for _, a in ports.values()]].mean(axis=0)) * (1.4/1.5)
    for i, (x, y) in enumerate(xy):
        conformer.SetAtomPosition(i, (float(x), float(y), 0.))
    # Port bonds remain in place during hydrogen addition, so a C=N port has
    # the same valence and H count as the actual inter-block C=N bond.
    mol = _embedded_block(mol, ports, xy, embed_seed, cap_element) if embed_seed is not None else Chem.AddHs(mol, addCoords=True)
    tilted_rings = 0 if embed_seed is not None else _tilt_phenyl_arms(mol, set(core), general_tilt)
    if cap_element is not None and embed_seed is None and not _relieve_fragment(mol, ports, cap_element):
        raise ValueError('The restrained fragment geometry did not converge.')
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
            'elements': [a.GetSymbol() for a in atoms], 'isotopes': [a.GetIsotope() for a in atoms], 'bonds': bonds,
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
        # In-plane reversal must be a proper 3D rotation to preserve handedness.
        xyz[:, 2] = xyz[:, 2]*mirror[i][1] + spacing/2
        fraction = xyz @ inverse
        wrap = np.floor(fraction + 1e-10).astype(int)
        wraps.append(wrap)
        atoms.extend({'element': element, 'fractional': p.tolist(),
                      **({'isotope': isotope} if isotope else {}), 'component': block['component']}
                     for element, isotope, p in zip(block['elements'], block['isotopes'], fraction-wrap))
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


def _interlayer_contact_height(blocks, pose):
    """Minimum AA translation keeping adjacent sheets at 90% of summed vdW radii."""
    _, cell, centers, angles, mirror, *_ = pose
    xy, z, radii = [], [], []
    table = Chem.GetPeriodicTable()
    for i, block in enumerate(blocks):
        xy.extend(_rotate(block['xyz'][:, :2]*mirror[i], angles[i]) + centers[i])
        z.extend(block['xyz'][:, 2] * mirror[i][1])
        radii.extend(table.GetRvdw(table.GetAtomicNumber(e)) for e in block['elements'])
    fractional = (np.asarray(xy) @ np.linalg.inv(cell)) % 1
    delta = fractional[:, None, :] - fractional[None, :, :]
    cutoff = .90 * np.add.outer(radii, radii)
    # All images within the largest possible contact distance lie inside this
    # range, even in an oblique cell (bound via the smallest singular value).
    span = math.ceil(1 + cutoff.max()/np.linalg.svd(cell, compute_uv=False)[-1])
    nearest = np.full(cutoff.shape, np.inf)
    for translation in product(range(-span, span+1), repeat=2):
        plane = (delta+translation) @ cell
        nearest = np.minimum(nearest, np.einsum('ijk,ijk->ij', plane, plane))
    projected = nearest < cutoff**2
    vertical = np.subtract.outer(z, z)[projected] + np.sqrt((cutoff**2-nearest)[projected])
    return math.ceil(float(vertical.max())*100)/100


def build_periodic_layer(design, *, spacing=None):
    """Build and check a single AA layer cell using the validated design."""
    template = design['selected_template']['id']
    if template not in SUPPORTED_TEMPLATES:
        raise ValueError('The chosen connectivity has no automatic periodic cell builder.')
    components = design['components']
    block_types, reference, positions, edges = _net(template, {k: c['connectivity'] for k, c in components.items()})
    chiral = any('@' in c['fragment_smiles'] for c in components.values())
    initial_caps = {}
    if chiral:
        for key, component in components.items():
            fragment = Chem.MolFromSmiles(component['fragment_smiles'])
            ends = {fragment.GetAtomWithIdx(i).GetSymbol() for _, i in _ports(fragment).values()}
            if len(ends) != 1:
                raise ValueError('Stereospecific construction requires equivalent attachment elements on each monomer.')
            initial_caps['B' if key == 'A' else 'A'] = ends.pop()
    imine_mode = None
    if design['linkage'] in ('imine', 'alkyl_imine'):
        # Source-derived cells are new constructions from precursors. Keep the
        # existing node pose and use opposed ends on their generated spacers.
        imine_mode = 'straight' if design.get('source_drawing_basis') else 'all'
    unique = {k: {**_block(c, embed_seed=2941 if chiral else None, cap_element=initial_caps.get(k),
                           imine_mode=imine_mode), 'component': k} for k, c in components.items()}
    blocks = [unique[k] for k in block_types]
    if sum(len(b['elements']) for b in blocks) > MAX_CELL_SITES:
        raise ValueError(f'Automatic cells are limited to {MAX_CELL_SITES} explicit sites.')
    order = design['connections'][0]['bond_order']
    poses = _fit_blocks(blocks, reference, positions, edges, 1.29 if order == 2 else 1.42)
    if not poses:
        raise ValueError('The selected building blocks could not close this cell with consistent bond lengths and directions.')
    # The attachment atoms and port tangents are fixed in all recovery poses.
    # Their fitted translations and junction lengths can therefore be reused.
    cap_elements = {}
    for key, block in unique.items():
        partner = unique['B' if key == 'A' else 'A']
        elements = {partner['elements'][i] for i in partner['sites']}
        cap_elements[key] = elements.pop() if len(elements) == 1 else None
    best = None
    reference_height = design.get('source_drawing_basis', {}).get('reference_layer_spacing_angstrom')
    minimum_height = None
    long_bonds = any(e in {'Si', 'P', 'S', 'Cl', 'Br', 'I'} for block in blocks for e in block['elements'])
    methods = ('restrained_mmff94',) if long_bonds else ('rigid', 'phenyl_torsions', 'restrained_mmff94')
    if chiral:
        methods = ('junction_constrained_etkdg_mmff',)
    for method in methods:
        if method in ('restrained_mmff94', 'phenyl_torsions'):
            if method == 'restrained_mmff94' and not all(cap_elements.values()):
                continue
            try:
                adjusted = {k: _block(c, general_tilt=45, imine_mode=imine_mode,
                                      cap_element=cap_elements[k] if method == 'restrained_mmff94' else None)
                            for k, c in components.items()}
            except ValueError:
                continue
            if all(np.allclose(adjusted[k]['xyz'], unique[k]['xyz']) for k in unique):
                continue
            blocks = [{**adjusted[k], 'component': k} for k in block_types]
        tilted_rings = sum(block['tilted_rings'] for block in blocks)
        z = np.concatenate([block['xyz'][:, 2] for block in blocks])
        if spacing is not None:
            spacings = [spacing]
        elif method == 'rigid':
            spacings = [5., 5.8] if tilted_rings else [3.4]
        else:
            spacings = [round(3.4 + float(np.ptp(z)), 2)]
        for height, pose in product(spacings, poses):
            contact_height = None
            if reference_height:
                contact_height = _interlayer_contact_height(blocks, pose)
                minimum_height = min(minimum_height, contact_height) if minimum_height is not None else contact_height
                if spacing is None:
                    height = max(reference_height, contact_height)
                elif height < contact_height:
                    continue
            structure, bonds, charges, lengths, angle_error = _cell_atoms(blocks, edges, pose, height, order)
            try:
                structure = normalize_structure(structure)
                network = inspect_periodic_network(structure)
            except ValueError:
                continue
            expected = {edge_key(a, b, t) for a, b, t, _ in bonds}
            observed = {edge_key(e['a'], e['b'], e['translation']) for e in network['edges']}
            extra, missing = len(observed - expected), len(expected - observed)
            short = network['short_contact_count']
            unstable = network['threshold_sensitivity']['classification_changed']
            if (short or network['dimensions'] != [2] or network['quotient_component_count'] != 1
                    or unstable or extra or missing or len(expected) != len(bonds)):
                score = short + missing + extra + int(unstable)
                if best is None or score < best[0]:
                    best = (score, f'{short} short contacts, {extra} unintended bonds and {missing} missing bonds'
                            f' at the distance cutoff (network ranks {network["dimensions"]},'
                            f' {network["quotient_component_count"]} component(s),'
                            f' cutoff-sensitive: {"yes" if unstable else "no"}).')
                continue
            structure['provenance'] = {
                'kind': 'cofdraw_periodic_starting_model', 'template': template, 'linkage': design['linkage'],
                'precursors': [c['monomer_smiles'] for c in components.values()],
                'geometry_method': 'Reacted fragments fitted to the chosen periodic connection graph',
                'fragment_geometry': method, 'corrugated_layer': bool(np.ptp(z) > .35),
                'preserves_drawing_coordinates': False, 'stacking': 'AA', 'spacing_angstrom': height,
                'phenyl_tilt_degrees': 60 if tilted_rings and method == 'rigid' else None, 'tilted_phenyl_ring_count': tilted_rings,
                'qualification': 'Unrelaxed periodic starting model with explicit hydrogens and chosen layer spacing.',
            }
            periodic_bonds = [{'a': int(a), 'b': int(b), 'translation': list(map(int, t)), 'order': kind}
                              for a, b, t, kind in bonds]
            if contact_height is not None:
                structure['provenance']['interlayer_contact_height_angstrom'] = contact_height
                structure['provenance']['interlayer_vdw_fraction'] = .90
            return structure, {
                'periodic_bonds': periodic_bonds, 'formal_charges': charges,
                'periodic_rank': 2, 'expected_atom_edges_verified': len(bonds), 'short_contacts': 0,
                'junction_length_range_angstrom': [float(lengths.min()), float(lengths.max())],
                'maximum_port_direction_error_degrees': angle_error,
            }
    if spacing is not None and minimum_height is not None and spacing < minimum_height:
        raise ValueError(f'This AA orientation needs at least {minimum_height:.2f} Å to clear the interlayer '
                         'contact screen. Use automatic spacing or supply another packing model.')
    detail = best[1] if best else 'No converged fragment geometry passed the structural checks.'
    raise ValueError(f'The {template} cell failed geometric validation. Best fit: {detail} '
                     'Inspect the precursor geometry and template, or import a corresponding CIF.')
