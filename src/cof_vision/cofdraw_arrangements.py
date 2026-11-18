"""Explicit translated-network hypotheses; no automatic catenation assignment."""
from __future__ import annotations

import hashlib
import json
import numpy as np

from .cofdraw_science import MAX_ATOMS, _PeriodicDistances, normalize_structure, structure_exports


def network_hypothesis(structure, offsets, repeats=(1, 1, 1)):
    """Copy an entire supplied periodic geometry with explicit fractional shifts.

    Copies retain their atom coordinates modulo the stated translation. They
    have no invented inter-network bonds. This is a user-chosen candidate,
    not an interpenetration classifier or an atomistic relaxation.
    """
    parent = normalize_structure(structure)
    shifts = np.asarray(offsets, dtype=float)
    if shifts.ndim != 2 or shifts.shape[1] != 3 or not 1 <= len(shifts) <= 4 or not np.isfinite(shifts).all() or np.max(abs(shifts)) > 1:
        raise ValueError('Use one to four explicit fractional network offsets, each with three values in [-1, 1]')
    shifts %= 1.
    if len({tuple(np.round(s, 8)) for s in shifts}) != len(shifts):
        raise ValueError('Network offsets must be distinct modulo the unit cell')
    if len(parent['atoms'])*len(shifts) > MAX_ATOMS:
        raise ValueError(f'Combined networks exceed {MAX_ATOMS} atoms per unit cell')
    if len({a.get('network', 0) for a in parent['atoms']}) > 1:
        raise ValueError('Choose the original single-network structure before making another hypothesis')
    coordinates = np.array([a['fractional'] for a in parent['atoms']])
    cell = np.array(parent['cell'])
    atoms, groups = [], []
    for index, shift in enumerate(shifts):
        moved = (coordinates+shift) % 1.
        groups.append(moved)
        atoms.extend({**a, 'fractional': xyz.tolist(), 'network': index} for a, xyz in zip(parent['atoms'], moved))
    closest = None
    for index, group in enumerate(groups):
        for other in groups[:index]:
            distance = float(_PeriodicDistances(cell, other, np.zeros(len(other))).clearance(group).min())
            closest = distance if closest is None else min(closest, distance)
    if closest is not None and closest < .8:
        raise ValueError(f'Network copies collide: closest intercopy atom distance is {closest:.3f} angstrom; change the offsets')
    warnings = ['Translated copies are a geometric hypothesis. Independent periodic-net connectivity and entanglement have not been established.',
                'The requested number of copies is not an experimentally assigned interpenetration degree or a COFex cat value.',
                'Inspect contacts and relaxation before using this candidate for material-property calculations.']
    if closest is not None and closest < 2.:
        warnings.append('Short intercopy contacts below 2 angstrom require review; no intercopy bonds are assigned.')
    provenance = {'kind': 'user_translated_network_hypothesis', 'source_geometry_sha256': hashlib.sha256(
                     json.dumps(parent, sort_keys=True, separators=(',', ':')).encode()).hexdigest(),
                  'copy_count': len(shifts), 'offsets_fractional': shifts.tolist(), 'interpenetration_degree': None,
                  'parent_provenance': parent.get('provenance', {}), 'qualification': warnings[0]}
    combined = {'cell': parent['cell'], 'periodic': True, 'atoms': atoms, 'provenance': provenance}
    result = structure_exports(combined, repeats=repeats)
    result.update(kind='translated_network_hypothesis', parameters={'copy_count': len(shifts), 'offsets_fractional': shifts.tolist()},
                  minimum_intercopy_distance_angstrom=closest, interpenetration_degree=None,
                  scope='Explicit translated periodic-network hypothesis; copy count does not establish interpenetration or chemical feasibility.')
    result['warnings'].extend(warnings)
    result['exports']['cif'] = '# Network-copy grouping and chosen offsets are recorded in the companion JSON.\n' + result['exports']['cif']
    return result
