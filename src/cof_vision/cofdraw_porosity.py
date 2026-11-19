"""Deterministic, bounded grid geometry; not Voronoi, adsorption or BET.

Clearance is 1-Lipschitz in Cartesian distance. This gives conservative
continuous edge capacities and enclosing voxel capacities on a periodic grid.
"""
from __future__ import annotations

import hashlib
import itertools
import json
import math
import numpy as np
from rdkit import rdBase

from .cofdraw_science import (AVOGADRO, ELEMENTS, PT, site_mass, _PeriodicDistances,
                             _finite, _integer, normalize_structure)


class _WindingUnion:
    """Integer translation potentials, including cycles across periodic faces."""
    def __init__(self, count):
        self.parent = np.arange(count)
        self.size = np.ones(count, dtype=np.int64)
        self.shift = np.zeros((count, 3), dtype=np.int64)
        self.basis = {}

    def find(self, vertex):
        vertex = int(vertex)
        parent = int(self.parent[vertex])
        if parent != vertex:
            root, offset = self.find(parent)
            self.shift[vertex] += offset
            self.parent[vertex] = root
        return int(self.parent[vertex]), self.shift[vertex].copy()

    def add_vector(self, root, vector):
        if not np.any(vector):
            return
        basis = self.basis.setdefault(root, [])
        if not basis or (len(basis) == 1 and np.any(np.cross(basis[0], vector))) or (
                len(basis) == 2 and int(np.dot(np.cross(basis[0], basis[1]), vector)) != 0):
            basis.append(vector.copy())

    def union(self, a, b, translation):
        ra, pa = self.find(a)
        rb, pb = self.find(b)
        shift = pa + translation - pb
        if ra == rb:
            self.add_vector(ra, shift)
            return bool(self.basis.get(ra))
        if self.size[ra] < self.size[rb]:
            ra, rb, shift = rb, ra, -shift
        self.parent[rb] = ra
        self.shift[rb] = shift
        self.size[ra] += self.size[rb]
        for vector in self.basis.pop(rb, []):
            self.add_vector(ra, vector)
        return bool(self.basis.get(ra))


def _edges(shape, diagonal=False):
    indices = np.indices(shape).reshape(3, -1).T
    start = np.arange(len(indices))
    steps = [(1, 0, 0), (0, 1, 0), (0, 0, 1)]
    if diagonal:
        steps = [step for step in itertools.product((-1, 0, 1), repeat=3)
                 if next((value for value in step if value), 0) > 0]
    for step in steps:
        moved = indices + step
        yield start, np.ravel_multi_index((moved % shape).T, shape), moved // shape, np.asarray(step)


def _threshold_components(clearance, edges, threshold):
    free = clearance >= threshold
    union = _WindingUnion(len(clearance))
    for starts, ends, translations, capacity in edges:
        for index in np.flatnonzero((capacity >= threshold) & free[starts] & free[ends]):
            union.union(starts[index], ends[index], translations[index])
    labels = np.full(len(clearance), -1, dtype=int)
    records = {}
    for vertex in np.flatnonzero(free):
        root, _ = union.find(vertex)
        labels[vertex] = root
        records[root] = records.get(root, 0) + 1
    winding = np.zeros(len(clearance), dtype=bool)
    components = []
    for root, count in records.items():
        basis = union.basis.get(root, [])
        if basis:
            winding[labels == root] = True
        components.append({'grid_points': count, 'periodic_rank': len(basis),
                           'winding_vectors': [vector.tolist() for vector in basis],
                           'percolating_axes': ['abc'[axis] for axis in range(3)
                                                if any(vector[axis] for vector in basis)]})
    components.sort(key=lambda row: (-row['periodic_rank'], -row['grid_points']))
    return winding, components


def _maximum_winding_capacity(count, edges):
    starts = np.concatenate([row[0] for row in edges])
    ends = np.concatenate([row[1] for row in edges])
    translations = np.concatenate([row[2] for row in edges])
    capacities = np.concatenate([row[3] for row in edges])
    union = _WindingUnion(count)
    for index in np.argsort(-capacities, kind='stable'):
        if capacities[index] <= 0:
            break
        if union.union(starts[index], ends[index], translations[index]):
            return float(capacities[index])
    return 0.


def _grid(distance, spacing):
    cell = distance.cell
    lengths = np.linalg.norm(cell, axis=1)
    requested = np.maximum(4, np.ceil(lengths / spacing).astype(int))
    image_cost = sum(len(offsets) for _, _, offsets in distance.components)
    limit = max(64, min(65536, 80_000_000 // (len(distance.radii) * image_cost)))
    shape = np.minimum(requested, 96)
    while math.prod(shape) > limit:
        axis = int(np.argmax(np.where(shape > 4, shape, 0)))
        shape[axis] -= 1
    points = (np.indices(shape).reshape(3, -1).T + .5) / shape
    clearance = distance.clearance(points)
    step_vectors = cell / shape[:, None]
    cover = float(max(np.linalg.norm(np.asarray(signs) @ step_vectors / 2)
                      for signs in itertools.product((-1, 1), repeat=3)))
    lower_edges = []
    for a, b, translation, step in _edges(shape):
        length = float(np.linalg.norm(step @ step_vectors))
        # The lower envelopes from both endpoint clearances meet halfway in
        # clearance space; their minimum is (ca+cb-length)/2.
        capacity = np.minimum(np.minimum(clearance[a], clearance[b]), (clearance[a] + clearance[b] - length) / 2)
        lower_edges.append((a, b, translation, capacity))
    upper_edges = [(a, b, translation, np.minimum(clearance[a], clearance[b]) + cover)
                   for a, b, translation, _ in _edges(shape, diagonal=True)]
    return {'shape': shape, 'requested_shape': requested, 'points': points, 'clearance': clearance,
            'cover': cover, 'lower_edges': lower_edges, 'upper_edges': upper_edges,
            'step_lengths': np.linalg.norm(step_vectors, axis=1), 'image_cost': image_cost,
            'point_limit': limit}


def _radical_inverse(count, base):
    integer = np.arange(1, count + 1, dtype=np.int64)
    result = np.zeros(count)
    fraction = 1. / base
    while integer.any():
        result += (integer % base) * fraction
        integer //= base
        fraction /= base
    return result


def _segment_image_shifts(distance, probe, maximum_length):
    radius = float(distance.radii.max()) + probe + maximum_length / 2
    fractional_bounds = radius * np.linalg.norm(np.linalg.inv(distance.cell), axis=0)
    bounds = np.floor(fractional_bounds + .5 + 1e-10).astype(int)
    if math.prod(2 * bounds + 1) > 2197:
        raise ValueError('Surface connectivity requires too many periodic images at this grid resolution')
    return np.asarray(list(itertools.product(*(range(-int(b), int(b)+1) for b in bounds)))), fractional_bounds


def _segments_clear(distance, starts, ends, probe, shifts, fractional_bounds):
    """Exact minimum atom-to-segment distances in a finite, proven image box."""
    starts, ends = np.asarray(starts), np.asarray(ends)
    safe = np.ones(len(starts), dtype=bool)
    for begin in range(0, len(starts), 64):
        a, b = starts[begin:begin+64], ends[begin:begin+64]
        midpoint = (a+b)/2
        vectors = (b-a) @ distance.cell
        squared = np.einsum('ij,ij->i', vectors, vectors)
        delta = midpoint[:, None, :] - distance.fractional[None]
        delta -= np.rint(delta)
        good = np.ones(len(a), dtype=bool)
        for shift in shifts:
            relative = delta - shift
            rows, atoms = np.where(good[:, None] & np.all(abs(relative) <= fractional_bounds + 1e-10, axis=2))
            if not len(rows):
                continue
            cartesian = relative[rows, atoms] @ distance.cell
            along = np.clip(-np.einsum('ij,ij->i', cartesian, vectors[rows]) / np.maximum(squared[rows], 1e-30), -.5, .5)
            minimum = np.linalg.norm(cartesian + along[:, None]*vectors[rows], axis=1)
            good[rows[minimum < distance.radii[atoms] + probe - 1e-7]] = False
            if not good.any():
                break
        safe[begin:begin+len(a)] = good
    return safe


def _surface(distance, grid, accessible, probe, requested_per_atom):
    inflated = distance.radii + probe
    areas = 4 * math.pi * inflated**2
    total = float(areas.sum())
    requested = max(512, len(inflated) * requested_per_atom)
    shifts, fractional_bounds = _segment_image_shifts(distance, probe, 2*grid['cover'])
    # Eight possible connections per surface point, each with bounded images.
    limit = min(32768, 160_000_000 // (len(inflated) * max(grid['image_cost'], 8*len(shifts))))
    if limit < 64:
        raise ValueError('The bounded surface-connectivity budget permits fewer than 64 samples; reduce the cell skew or use a dedicated porosity engine')
    count = min(requested, limit)
    chosen = np.searchsorted(np.cumsum(areas)/total, _radical_inverse(count, 2), side='right')
    z = 1 - 2*_radical_inverse(count, 3)
    angle = 2*math.pi*_radical_inverse(count, 5)
    directions = np.column_stack((np.sqrt(1-z*z)*np.cos(angle), np.sqrt(1-z*z)*np.sin(angle), z))
    points = distance.fractional[chosen] + (directions*inflated[chosen, None]) @ np.linalg.inv(distance.cell)
    exposed = distance.clearance(points) >= probe-1e-7
    connected = np.zeros(count, dtype=bool)
    base = np.floor(points*grid['shape']-.5).astype(int)
    for offset in itertools.product((0, 1), repeat=3):
        candidates = base + offset
        ids = np.ravel_multi_index((candidates % grid['shape']).T, grid['shape'])
        use = np.flatnonzero(exposed & ~connected & accessible[ids])
        if len(use):
            ends = (candidates[use]+.5)/grid['shape']
            connected[use] = _segments_clear(distance, points[use], ends, probe, shifts, fractional_bounds)
    return {'available': True, 'samples': count, 'requested_samples': requested, 'sample_limit': limit,
            'exposed_samples': int(exposed.sum()),
            'connected_samples': int(connected.sum()), 'unresolved_exposed_samples': int((exposed & ~connected).sum()),
            'local_area': total*float(exposed.mean()), 'accessible_area': total*float(connected.mean()),
            'method': 'Deterministic Halton(2,3,5) area-weighted sphere quadrature; exact collision-free segments to winding grid components',
            'periodic_segment_images': len(shifts)}


def _evaluate(distance, spacing, probe):
    grid = _grid(distance, spacing)
    c, cover = grid['clearance'], grid['cover']
    accessible, components = _threshold_components(c, grid['lower_edges'], probe)
    lower_edges = [(a, b, t, np.minimum(c[a], c[b])-cover) for a, b, t, _ in grid['lower_edges']]
    lower, _ = _threshold_components(c-cover, lower_edges, probe)
    upper, _ = _threshold_components(c+cover, grid['upper_edges'], probe)
    pld_low = 2*_maximum_winding_capacity(len(c), grid['lower_edges'])
    pld_high = 2*_maximum_winding_capacity(len(c), grid['upper_edges'])
    lcd_low = 2*max(0., float(c.max()))
    lcd_high = 2*max(0., float(c.max())+cover)
    return grid, accessible, {'lcd_angstrom': lcd_low, 'lcd_bounds_angstrom': [lcd_low, lcd_high],
            'pld_angstrom': pld_low, 'pld_bounds_angstrom': [pld_low, min(pld_high, lcd_high)],
            'void_fraction': float((c >= 0).mean()),
            'void_fraction_bounds': [float((c-cover >= 0).mean()), float((c+cover >= 0).mean())],
            'probe_center_free_fraction': float((c >= probe).mean()),
            'transport_accessible_probe_center_fraction': float(accessible.mean()),
            'transport_accessible_probe_center_fraction_bounds': [float(lower.mean()), float(upper.mean())],
            'connectivity': {'status': 'certified_channel' if accessible.any() else
                             'unresolved_at_this_resolution' if upper.any() else 'no_channel_at_this_probe',
                             'certified_winding': bool(accessible.any()), 'possible_winding_in_upper_graph': bool(upper.any()),
                             'components': components[:100], 'component_count': len(components),
                             'components_truncated': len(components)>100,
                             'maximum_periodic_rank': max((row['periodic_rank'] for row in components), default=0)}}


def geometric_porosity(structure, *, probe_radius=1.82, grid_spacing=.6, surface_samples=64, refinement=False, include_map=False):
    """Level-1 deterministic periodic hard-sphere geometry with explicit bounds."""
    normalized = normalize_structure(structure)
    probe = _finite(probe_radius, 'probe_radius', 0, 5)
    spacing = _finite(grid_spacing, 'grid_spacing', .1, 3)
    per_atom = _integer(surface_samples, 'surface_samples', 16, 256)
    if type(refinement) is not bool:
        raise ValueError('refinement must be a boolean')
    atoms = normalized['atoms']
    elements = [atom['element'] for atom in atoms]
    radii = np.asarray([PT.GetRvdw(ELEMENTS[element]) for element in elements])
    if not np.isfinite(radii).all() or (radii <= 0).any():
        raise ValueError('Every element needs a positive van der Waals radius')
    mass = float(sum(site_mass(a) for a in atoms))
    volume = normalized['cell_volume_angstrom3']
    density = mass/AVOGADRO/(volume*1e-24)
    distance = _PeriodicDistances(np.asarray(normalized['cell']), np.asarray([atom['fractional'] for atom in atoms]), radii)
    grid, accessible, result = _evaluate(distance, spacing, probe)
    if include_map:
        import base64
        # X, Y, Z fractional voxel centers, Z varying fastest. The mask is
        # evaluated at full precision before serializing the display field.
        result['pore_map'] = {
            'shape': grid['shape'].tolist(), 'cell': distance.cell.tolist(),
            'probe_radius': probe,
            'clearance_f32le': base64.b64encode(grid['clearance'].astype('<f4').tobytes()).decode('ascii'),
            'connected_bits': base64.b64encode(np.packbits(accessible, bitorder='little').tobytes()).decode('ascii'),
        }
    warnings = ['Fixed, supplied periodic geometry and elemental hard spheres; not adsorption, BET, helium porosimetry or a relaxed-crystal prediction.',
                'PLD is a conservative continuous-path grid estimate; LCD is a sampled lower estimate. Inspect their bounds before comparing structures.',
                'Accessible volume means transport-connected probe-center volume, not probe-occupiable or experimental adsorption volume.',
                'Surface is the boundary of probe-inflated spheres. Samples without a certified connection to a winding grid component are excluded; inaccessible and unresolved samples are not distinguished.',
                'Grid bounds cover discretization of the stated hard-sphere model only; they omit structural and radius uncertainty. Surface quadrature has no statistical confidence interval.']
    try:
        surface = _surface(distance, grid, accessible, probe, per_atom)
        area_factor = 1e-20*AVOGADRO/mass
        result.update(transport_accessible_surface_angstrom2=surface['accessible_area'],
                      transport_accessible_surface_m2_g=surface['accessible_area']*area_factor,
                      local_probe_exposed_surface_angstrom2=surface['local_area'],
                      local_probe_exposed_surface_m2_g=surface['local_area']*area_factor)
        if surface['samples'] < surface['requested_samples']:
            warnings.append('Surface sampling was reduced by the bounded-work budget; inspect the effective count and refine independently.')
    except ValueError as exc:
        surface = {'available': False, 'reason': str(exc)}
        result.update(transport_accessible_surface_angstrom2=None, transport_accessible_surface_m2_g=None,
                      local_probe_exposed_surface_angstrom2=None, local_probe_exposed_surface_m2_g=None)
        warnings.append(str(exc))
    if not np.array_equal(grid['shape'], grid['requested_shape']):
        warnings.append('Grid resolution was coarsened by the per-axis or work limit. Effective shape and rigorous geometric bounds are reported.')
    if result['connectivity']['status'] == 'unresolved_at_this_resolution':
        warnings.append('No transport path was certified at this resolution, but the enclosing graph permits one. Zero accessible estimates are unresolved, not evidence of a closed pore; refine the grid or use a dedicated engine.')
    elif surface.get('exposed_samples', 0) and not surface.get('connected_samples', 0) and result['connectivity']['certified_winding']:
        warnings.append('A channel was certified but no exposed surface sample obtained a certified connection. The zero accessible-area estimate remains unresolved.')
    if 'H' not in elements:
        warnings.append('No explicit hydrogen atoms are listed; missing atoms alter both geometry and mass-normalized descriptors.')
    result.update(kind='deterministic_periodic_grid_geometry', engine='cofdraw-grid-v1',
        structure_sha256=hashlib.sha256(json.dumps(normalized, sort_keys=True, separators=(',', ':')).encode()).hexdigest(),
        atom_count=len(atoms), density_g_cm3=density, cell_volume_angstrom3=volume, cell_molar_mass_g_mol=mass,
        geometric_void_volume_cm3_g=result['void_fraction']/density,
        probe_center_free_volume_cm3_g=result['probe_center_free_fraction']/density,
        transport_accessible_probe_center_volume_angstrom3=result['transport_accessible_probe_center_fraction']*volume,
        transport_accessible_probe_center_volume_cm3_g=result['transport_accessible_probe_center_fraction']/density,
        transport_accessible_probe_center_volume_bounds_cm3_g=[value/density for value in result['transport_accessible_probe_center_fraction_bounds']],
        method={'backend': 'bounded periodic grid, not Zeo++/Voronoi', 'probe_radius_angstrom': probe,
                'requested_grid_spacing_angstrom': spacing, 'grid_shape': grid['shape'].tolist(),
                'requested_grid_shape': grid['requested_shape'].tolist(),
                'grid_points': len(grid['points']), 'grid_point_limit': grid['point_limit'],
                'grid_capped': not np.array_equal(grid['shape'], grid['requested_shape']),
                'periodic_image_searches_per_pair': grid['image_cost'],
                'effective_axis_steps_angstrom': grid['step_lengths'].tolist(),
                'voxel_covering_radius_angstrom': grid['cover'],
                'distance': 'Exact bounded triclinic minimum-image clearance with RDKit elemental van der Waals radii',
                'radii_angstrom': {element: PT.GetRvdw(ELEMENTS[element]) for element in sorted(set(elements))},
                'connectivity': 'Integer translation cycles; 6-neighbor certified line segments for lower paths and 26-neighbor enclosing cells for upper bounds',
                'surface': surface, 'rdkit_version': rdBase.rdkitVersion,
                'references': ['https://www.zeoplusplus.org/about.html', 'https://github.com/SarkisovGitHub/PoreBlazer']},
        convergence={'requested': refinement, 'performed': False}, bet_surface_area_m2_g=None,
        scope='Deterministic geometric descriptors of the supplied periodic hard-sphere model; transport-connected grid approximation, not experimental porosity or adsorption.',
        warnings=warnings)
    if refinement:
        finer, _, refined = _evaluate(distance, max(.05, spacing/2), probe)
        different = not np.array_equal(finer['shape'], grid['shape'])
        result['convergence'] = {'requested': True, 'performed': different,
            'refined_grid_shape': finer['shape'].tolist(), 'refined_grid_points': len(finer['points']),
            'reason': None if different else 'The work cap prevents a finer grid',
            'refined': {key: refined[key] for key in ('pld_angstrom', 'pld_bounds_angstrom', 'lcd_angstrom', 'lcd_bounds_angstrom',
                'void_fraction', 'transport_accessible_probe_center_fraction', 'transport_accessible_probe_center_fraction_bounds')},
            'interpretation': 'A second grid is a resolution diagnostic, not proof of convergence; surface quadrature is not refined by this option.'}
    return result
