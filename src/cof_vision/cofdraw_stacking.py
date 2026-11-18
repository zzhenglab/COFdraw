"""Explicit stacking hypotheses for a supplied single planar periodic layer."""
from __future__ import annotations

import numpy as np

from .cofdraw_current_model import model_hash
from .cofdraw_network import inspect_periodic_network
from .cofdraw_science import normalize_structure, structure_exports, _finite, MAX_ATOMS


def _connected_layer_lifts(network, atom_count):
    neighbors = [[] for _ in range(atom_count)]
    for edge in network['edges']:
        a, b, tz = edge['a'], edge['b'], edge['translation'][2]
        neighbors[a].append((b, tz))
        neighbors[b].append((a, -tz))
    lifts, todo = {0: 0}, [0]
    while todo:
        a = todo.pop()
        for b, tz in neighbors[a]:
            height = lifts[a]+tz
            if b in lifts:
                if height != lifts[b]:
                    raise ValueError('The covalent layer connects periodically along c.')
            else:
                lifts[b] = height
                todo.append(b)
    if len(lifts) != atom_count:
        raise ValueError('Corrugated stacking requires one connected periodic layer in the cell.')
    return np.array([lifts[i] for i in range(atom_count)])


def periodic_stacking(structure, *, stacking='AA', spacing=3.4, slip=(1.5, 0), corrugated_layer=False):
    parent = normalize_structure(structure)
    if stacking not in ('AA', 'AB', 'ABC', 'slipped'):
        raise ValueError('Periodic stacking supports AA, AB, ABC and slipped; arbitrary twist is not a periodic cell construction')
    spacing = _finite(spacing, 'spacing', 2.5, 12)
    shift = np.asarray(slip, dtype=float)
    if shift.shape != (2,) or not np.isfinite(shift).all() or np.max(abs(shift)) > 20:
        raise ValueError('Slip requires two finite in-plane angstrom offsets between -20 and 20')
    cell = np.asarray(parent['cell'])
    normal = np.cross(cell[0], cell[1]); normal /= np.linalg.norm(normal)
    axis_x = cell[0] / np.linalg.norm(cell[0]); axis_y = np.cross(normal, axis_x)
    height = float(np.dot(cell[2], normal))
    fraction = np.asarray([a['fractional'] for a in parent['atoms']])
    # A wrapped layer may straddle z=0. Cut at the largest empty fractional gap.
    heights = np.sort(fraction[:, 2] % 1)
    gaps = np.diff(np.r_[heights, heights[0] + 1])
    start = heights[(int(np.argmax(gaps)) + 1) % len(heights)]
    unwrapped = (fraction[:, 2] - start) % 1
    thickness = float(np.ptp(unwrapped) * height)
    if thickness > .35 and not corrugated_layer:
        raise ValueError('Stack comparison requires one near-planar layer parallel to a,b (thickness ≤0.35 Å). Choose a single-layer model; 3D networks and multilayer cells cannot be restacked automatically.')
    # Lift the boundary-wrapped atoms into one physical layer before checking
    # its connectivity. A raw edge translation along c may merely join atoms
    # straddling the unit-cell boundary; its translation in this lifted gauge
    # is t_c + lift[a] - lift[b].
    lifts = np.rint(start + unwrapped - fraction[:, 2]).astype(int)
    try:
        network = inspect_periodic_network(parent)
    except ValueError as exc:
        raise ValueError('Cannot verify a single periodic layer within the bounded contact search: ' + str(exc)) from exc
    if network['short_contact_count']:
        raise ValueError('Cannot verify a single periodic layer: implausibly short contacts require inspection before restacking')
    if corrugated_layer:
        lifts = _connected_layer_lifts(network, len(fraction))
        lifted_heights = fraction[:, 2]+lifts
        start, thickness = float(lifted_heights.min()), float(np.ptp(lifted_heights)*height)
        if thickness >= spacing-.05:
            raise ValueError(f'Layer thickness is {thickness:.2f} Å. Increase the stacking spacing.')
    contacts = network['edges'] + network['borderline_contacts']
    if any(edge['translation'][2] + int(lifts[edge['a']]) - int(lifts[edge['b']]) != 0 for edge in contacts):
        raise ValueError('Cannot restack: distance-inferred or borderline covalent contacts connect different layers along c')
    if any(component['periodic_rank'] != 2 for component in network['components']):
        raise ValueError('Cannot verify a 2D periodic covalent layer: the distance-inferred graph contains molecular, chain, isolated, or 3D components')
    if network['threshold_sensitivity']['classification_changed']:
        raise ValueError('Cannot verify a 2D periodic covalent layer: connectivity is sensitive to the distance cutoff')
    count = {'AA': 1, 'AB': 2, 'ABC': 3, 'slipped': 1}[stacking]
    if len(parent['atoms']) * count > MAX_ATOMS:
        raise ValueError('Stacked unit cell exceeds the atom limit')
    # Preserve the full Cartesian geometry, including c's in-plane component
    # in oblique cells. Dropping that component would shear a corrugated layer.
    # Only a common translation centers the physical layer on the new origin.
    lifted_fraction = fraction.copy()
    lifted_fraction[:, 2] += lifts
    base = lifted_fraction @ cell
    base -= normal * (start * height + .5 * thickness)
    delta = shift[0] * axis_x + shift[1] * axis_y
    c = count * spacing * normal + (delta if stacking == 'slipped' else 0)
    new_cell = np.vstack([cell[:2], c])
    inverse = np.linalg.inv(new_cell)
    atoms, wraps = [], []
    for layer in range(count):
        translation = layer * spacing * normal + (layer * delta if stacking in ('AB', 'ABC') else 0)
        raw = (base+translation) @ inverse
        wrap = np.floor(raw+1e-10).astype(int)
        wraps.append(wrap)
        for atom, point in zip(parent['atoms'], raw-wrap):
            atoms.append({**atom, 'fractional': point.tolist(), 'layer': layer})
    model = {'periodic': True, 'cell': new_cell.tolist(), 'atoms': atoms,
             'provenance': {'kind': 'explicit_periodic_stacking_hypothesis', 'stacking': stacking,
                            'spacing_angstrom': spacing, 'slip_angstrom_in_ab_frame': shift.tolist(),
                            'input_structure_sha256': model_hash(parent), 'parent': parent.get('provenance', {}),
                            'layer_validation': {'method': 'distance-inferred periodic covalent graph plus geometric thickness',
                                                 'component_count': len(network['components']), 'periodic_rank': 2,
                                                 'cross_layer_contacts': 0,
                                                 'bond_scale': network['method']['bond_scale'],
                                                 'tolerance_angstrom': network['method']['tolerance_angstrom'],
                                                 'sensitivity_margin_angstrom': network['threshold_sensitivity']['distance_margin_angstrom'],
                                                 'experimentally_assigned_connectivity': False},
                            'qualification': 'Chosen unrelaxed stacking of the supplied planar layer; no experimental stacking assignment.'}}
    if corrugated_layer:
        output_network = inspect_periodic_network(model)
        edges = set()
        for layer, wrap in enumerate(wraps):
            offset = layer*len(fraction)
            for edge in network['edges']:
                a, b = edge['a'], edge['b']
                image = np.r_[edge['translation'][:2], 0]+wrap[b]-wrap[a]
                edges.add((offset+a, offset+b, tuple(image)))
        observed = {(e['a'],e['b'],tuple(e['translation'])) for e in output_network['edges']}
        if output_network['short_contact_count'] or edges != observed:
            raise ValueError('This stacking creates interlayer contacts or changes the inferred covalent bonds. Adjust the spacing or slip.')
        model['provenance']['layer_validation']['output_atom_edges_verified'] = len(edges)
    result = structure_exports(model)
    result.update(kind='periodic_stacking_hypothesis',
                  parameters={'stacking': stacking, 'spacing_angstrom': spacing, 'slip_angstrom': shift.tolist(),
                              'layer_thickness_angstrom': thickness, 'layers_per_cell': count},
                  scope='Periodic stacking hypothesis from the same supplied planar layer. AA repeats one layer; AB repeats offsets 0,s; ABC repeats 0,s,2s; slipped has an oblique c vector.')
    result['warnings'].append('Chosen layer spacing and shifts are unrelaxed assumptions. Inspect contacts before interpreting diffraction or pore properties; arbitrary twist needs a commensurate cell and is unavailable here.')
    result['warnings'].append('Input-layer connectivity was screened using covalent radii and distances, including periodic images and cutoff sensitivity; it is not experimentally assigned bonding or a named topology.')
    return result
