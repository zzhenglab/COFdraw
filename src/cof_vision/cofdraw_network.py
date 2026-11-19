"""Bounded distance-inferred periodic connectivity, never a net assignment.

Lattice vectors are rows. An edge (a, b, t) joins atom a in cell zero to atom b
in cell t. All images within the cutoff are considered, including self images;
minimum-image-only or finite-supercell connectivity would lose this information.
"""
from __future__ import annotations

from collections import Counter, deque
import hashlib
import itertools
import json
import math

import numpy as np
from rdkit import Chem, rdBase

from .cofdraw_science import normalize_structure

MAX_IMAGE_PAIRS = 120_000_000
MAX_IMAGES = 1331
MAX_CONTACTS = 30_000
SENSITIVITY_MARGIN = .10


def _number(value, name, low, high):
    if isinstance(value, bool):
        raise ValueError(f"{name} must be a finite number between {low} and {high}")
    try:
        value = float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{name} must be a finite number") from exc
    if not math.isfinite(value) or not low <= value <= high:
        raise ValueError(f"{name} must be between {low} and {high}")
    return value


def _cross(a, b):
    return (a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0])


def _append_independent(basis, vector):
    """Exact rank in Z^3 using Python integers, not a floating SVD tolerance."""
    vector = tuple(int(v) for v in vector)
    if not any(vector) or len(basis) == 3:
        return
    if not basis or (len(basis) == 1 and any(_cross(basis[0], vector))):
        basis.append(vector)
    elif len(basis) == 2 and sum(a*b for a, b in zip(_cross(*basis), vector)) != 0:
        basis.append(vector)


def _translation_basis(vectors):
    """Integer row elimination preserves the subgroup, including nonprimitive periods."""
    rows = sorted(set(vectors) - {(0, 0, 0)})
    basis = []
    for axis in range(3):
        while any(row[axis] for row in rows):
            pivot = min((r for r in rows if r[axis]), key=lambda r: abs(r[axis]))
            rows.remove(pivot)
            if pivot[axis] < 0:
                pivot = tuple(-x for x in pivot)
            rows = [tuple(x - (r[axis] // pivot[axis]) * y for x, y in zip(r, pivot)) for r in rows]
            if not any(r[axis] for r in rows):
                basis.append(pivot)
                break
            rows.append(pivot)
    return basis


def label_periodic_networks(viewer, structure, repeats, edges):
    """Color translated copies by cosets of the full-rank cycle lattice."""
    groups, _ = _components(len(structure['atoms']), edges, np.asarray(structure['cell']), track_translations=True)
    labels = {}
    for group in groups:
        if group['periodic_rank'] != 3:
            continue
        basis = _translation_basis(group['cycles'])
        for image, shift in enumerate(itertools.product(*(range(n) for n in repeats))):
            for site, path in group['paths'].items():
                v = [x-y for x,y in zip(shift,path)]
                for axis,row in enumerate(basis):
                    factor = v[axis] // row[axis]
                    v = [x-factor*y for x,y in zip(v,row)]
                key = (group['id'], *v)
                viewer['atoms'][image*len(structure['atoms'])+site]['network'] = labels.setdefault(key,len(labels))


def _components(atom_count, edges, cell, *, track_translations=False):
    adjacency = [[] for _ in range(atom_count)]
    for edge in edges:
        a, b, shift = edge['a'], edge['b'], tuple(edge['translation'])
        adjacency[a].append((b, shift))
        adjacency[b].append((a, tuple(-v for v in shift)))
    visited = set()
    components = []
    for root in range(atom_count):
        if root in visited:
            continue
        potential = {root: (0, 0, 0)}
        queue = deque([root]); visited.add(root)
        basis = []
        cycles = set()
        while queue:
            a = queue.popleft()
            for b, shift in adjacency[a]:
                expected = tuple(x+y for x, y in zip(potential[a], shift))
                if b not in potential:
                    potential[b] = expected; visited.add(b); queue.append(b)
                else:
                    cycle = tuple(x-y for x, y in zip(expected, potential[b]))
                    _append_independent(basis, cycle)
                    if track_translations:
                        cycles.add(cycle)
        indices = sorted(potential)
        components.append({
            'id': len(components), 'atom_indices': indices, 'atom_count': len(indices),
            'edge_count': sum(len(adjacency[a]) for a in indices)//2,
            'periodic_rank': len(basis), 'inferred_dimension': len(basis),
            'independent_cycle_translations': [list(v) for v in basis],
            'independent_cycle_vectors_angstrom': [(np.asarray(v) @ cell).tolist() for v in basis],
            **({'cycles': cycles, 'paths': potential} if track_translations else {}),
        })
    return components, [len(rows) for rows in adjacency]


def _dimension_summary(components):
    counts = Counter(str(row['inferred_dimension']) for row in components)
    return {'dimensions': sorted(int(x) for x in counts),
            'maximum_dimension': max(int(x) for x in counts),
            'quotient_component_count': len(components),
            'quotient_components_by_dimension': dict(sorted(counts.items()))}


def _candidate_contacts(cell, fractional, radii, scale, tolerance):
    maximum_cutoff = float(2*radii.max()*scale+tolerance+SENSITIVITY_MARGIN)
    # If |d| <= R and d = (delta + t) @ cell, then
    # |delta_k+t_k| <= R * ||inverse(cell)[:,k]||. Fractional positions are
    # wrapped; their coordinate spread provides a complete integer-image bound.
    reciprocal_norm = np.linalg.norm(np.linalg.inv(cell), axis=0)
    bounds = np.floor(maximum_cutoff*reciprocal_norm+np.ptp(fractional, axis=0)+1e-10).astype(int)
    image_count = math.prod(int(2*b+1) for b in bounds)
    n = len(fractional)
    work = image_count*n*n
    if image_count > MAX_IMAGES or work > MAX_IMAGE_PAIRS:
        raise ValueError("Periodic contact search exceeds its work bound; use a reduced/primitive cell or fewer explicit sites")
    contacts = []
    closest = None
    for shift in itertools.product(*(range(-int(b), int(b)+1) for b in bounds)):
        # a < b retains every image; a == b retains one of t and -t.
        positive_self = next((v > 0 for v in shift if v), False)
        for start in range(0, n, 64):
            stop = min(n, start+64)
            a_indices = np.arange(start, stop)[:, None]
            b_indices = np.arange(n)[None, :]
            delta = fractional[None, :, :]+shift-fractional[start:stop, None, :]
            distances = np.linalg.norm(delta @ cell, axis=2)
            sums = radii[start:stop, None]+radii[None, :]
            cutoffs = sums*scale+tolerance
            canonical = a_indices < b_indices
            if positive_self:
                canonical |= a_indices == b_indices
            # A common envelope also captures very short contacts for pairs
            # with smaller element-specific cutoffs without inventing edges.
            candidates = canonical & (distances <= maximum_cutoff+1e-10)
            if np.any(candidates):
                value = float(distances[candidates].min())
                closest = value if closest is None else min(closest, value)
            rows, columns = np.nonzero(canonical & (distances <= cutoffs+SENSITIVITY_MARGIN+1e-10))
            for local, b in zip(rows, columns):
                a = start+int(local); b = int(b)
                length = float(distances[local, b]); radius_sum = float(sums[local, b])
                contacts.append({'a': a, 'b': b, 'translation': list(shift),
                                 'distance_angstrom': length, 'radius_sum_angstrom': radius_sum,
                                 'cutoff_angstrom': float(cutoffs[local, b]),
                                 'short_contact': length < max(.35, .65*radius_sum)})
                if len(contacts) > MAX_CONTACTS:
                    raise ValueError("Too many short-range periodic contacts; inspect cell units/overlaps before network inference")
    return contacts, closest, {'translation_bounds': bounds.tolist(), 'image_count': image_count,
                               'directed_image_pair_upper_bound': work,
                               'search_envelope_angstrom': maximum_cutoff}


def inspect_periodic_network(structure, *, bond_scale=1.15, tolerance_angstrom=.15):
    """Inspect connectivity implied by a radius cutoff in a full periodic cell.

    No bond orders, charges, metal coordination chemistry, RCSR name, or
    interpenetration multiplicity are assigned. Short implausible contacts are
    reported separately and excluded; their presence makes inference unreliable.
    """
    scale = _number(bond_scale, 'bond_scale', .8, 1.4)
    tolerance = _number(tolerance_angstrom, 'tolerance_angstrom', 0, .4)
    normalized = normalize_structure(structure)
    cell = np.asarray(normalized['cell'])
    fractional = np.asarray([a['fractional'] for a in normalized['atoms']])
    elements = [a['element'] for a in normalized['atoms']]
    table = Chem.GetPeriodicTable()
    radii = np.array([table.GetRcovalent(table.GetAtomicNumber(e)) for e in elements])
    if not np.isfinite(radii).all() or np.any(radii <= 0):
        raise ValueError("Positive covalent radii are required for every listed element")
    contacts, closest, search = _candidate_contacts(cell, fractional, radii, scale, tolerance)
    short = [row for row in contacts if row['short_contact']]
    usable = [row for row in contacts if not row['short_contact']]
    edges = [row for row in usable if row['distance_angstrom'] <= row['cutoff_angstrom']+1e-10]
    components, degree = _components(len(elements), edges, cell)
    summary = _dimension_summary(components)
    for row in components:
        row['element_counts'] = dict(sorted(Counter(elements[i] for i in row['atom_indices']).items()))
    sensitivity = {}
    for name, margin in (('minus_margin', -SENSITIVITY_MARGIN), ('plus_margin', SENSITIVITY_MARGIN)):
        selected = [row for row in usable if row['distance_angstrom'] <= row['cutoff_angstrom']+margin+1e-10]
        sensitivity[name] = _dimension_summary(_components(len(elements), selected, cell)[0])
        sensitivity[name]['edge_count'] = len(selected)
    dimension_change = any(row['quotient_components_by_dimension'] != summary['quotient_components_by_dimension']
                           for row in sensitivity.values())
    isolates = [i for i, count in enumerate(degree) if count == 0]
    # Neighbor counts are not valence. These deliberately loose checks identify
    # common organic-atom anomalies; they do not exclude unusual valid chemistry.
    upper = {'H': 1, 'B': 4, 'C': 4, 'N': 4, 'O': 3, 'F': 1, 'Cl': 1, 'Br': 1, 'I': 1}
    unusual = [{'atom_index': i, 'element': e, 'coordination': degree[i],
                'reason': 'above_common_organic_neighbor_count'}
               for i, e in enumerate(elements) if e in upper and degree[i] > upper[e]]
    histogram = {}
    for element in sorted(set(elements)):
        histogram[element] = dict(sorted(Counter(str(degree[i]) for i, e in enumerate(elements) if e == element).items()))
    borderline = [row for row in usable if abs(row['distance_angstrom']-row['cutoff_angstrom']) <= SENSITIVITY_MARGIN+1e-10]
    warnings = [
        'Bonds are inferred from covalent radii and distances, not read as experimentally assigned chemical connectivity. No bond orders or formal valences are assigned.',
        'Dimension describes the inferred periodic atomic graph. It does not assign an RCSR net, COF identity, stacking, synthesis feasibility, or crystal stability.',
        'Component counts refer to the finite quotient graph. Equivalent disconnected translated networks and interpenetration multiplicity are not determined.',
    ]
    if short:
        warnings.append(f'{len(short)} implausibly short periodic contacts were excluded; inspect units, occupancy, and geometry before using inferred dimensions.')
    if isolates:
        warnings.append(f'{len(isolates)} sites have no inferred neighbors; these may be guests, missing connections, unsuitable radii, or incorrect coordinates.')
    if len(components) > 1:
        warnings.append('Several quotient components are present; separate layers, molecular guests, or disconnected fragments may be expected or may indicate incomplete geometry.')
    if unusual:
        warnings.append(f'{len(unusual)} atoms exceed a loose common-organic coordination bound; inspect contacts and the cutoff. Neighbor count is not chemical valence.')
    if borderline:
        warnings.append(f'{len(borderline)} contacts lie within 0.10 angstrom of the cutoff; inspect threshold sensitivity before interpreting connectivity.')
    if dimension_change:
        warnings.append('Component/dimension classification changes when the distance cutoff shifts by 0.10 angstrom; the inferred network is threshold-sensitive.')
    if 'H' not in elements:
        warnings.append('No explicit hydrogen sites are supplied. Missing H atoms can change coordination diagnostics; no hydrogens are added.')
    organic = {'H', 'B', 'C', 'N', 'O', 'F', 'Si', 'P', 'S', 'Cl', 'Br', 'I'}
    if any(e not in organic for e in elements):
        warnings.append('Elements outside the common organic set are present. Covalent-radius cutoffs may miss or invent metal/ionic contacts; coordination chemistry needs specialist review.')
    distances = [row['distance_angstrom'] for row in edges]
    return {
        'kind': 'distance_inferred_periodic_network',
        'structure_sha256': hashlib.sha256(json.dumps(normalized, sort_keys=True, allow_nan=False).encode()).hexdigest(),
        'scope': 'Geometric connectivity screening of an explicit periodic structure; no experimentally assigned bonding or named topology is implied.',
        'units': {'distance': 'angstrom', 'translation': 'integer_cell_vectors', 'dimension': 'rank'},
        'method': {'bond_scale': scale, 'tolerance_angstrom': tolerance,
                   'cutoff': 'distance <= bond_scale * (r_cov[a] + r_cov[b]) + tolerance_angstrom',
                   'radii': 'RDKit periodic-table covalent radii',
                   'radii_angstrom': {e: float(table.GetRcovalent(table.GetAtomicNumber(e))) for e in sorted(set(elements))},
                   'short_contact_rule': 'distance < max(0.35 angstrom, 0.65 * covalent_radius_sum); excluded from inferred edges',
                   'periodic_search': 'All integer images inside a reciprocal-cell bound; self-image edges retained',
                   'dimensionality': 'Exact integer rank of quotient-graph cycle translation vectors',
                   'edge_convention': '(fractional[b] + translation - fractional[a]) @ cell; reverse edge is equivalent',
                   'component_scope': 'Quotient connected components, not the count of connected objects in the infinite crystal',
                   'rdkit_version': rdBase.rdkitVersion, **search},
        'atom_count': len(elements), 'explicit_hydrogen_count': elements.count('H'),
        'edge_count': len(edges), 'edges': edges,
        **summary, 'components': components,
        'coordination': {'per_atom': [{'atom_index': i, 'element': e, 'coordination': degree[i]} for i, e in enumerate(elements)],
                         'histogram_by_element': histogram, 'mean': float(np.mean(degree)),
                         'maximum': max(degree), 'unusual_atoms': unusual},
        'isolated_atom_indices': isolates,
        'minimum_candidate_distance_angstrom': closest,
        'edge_distance_range_angstrom': [min(distances), max(distances)] if distances else None,
        'short_contact_count': len(short), 'short_contacts': short,
        'borderline_contact_count': len(borderline), 'borderline_contacts': borderline,
        'threshold_sensitivity': {'distance_margin_angstrom': SENSITIVITY_MARGIN,
                                  'classification_changed': dimension_change, **sensitivity},
        'warnings': warnings,
    }
