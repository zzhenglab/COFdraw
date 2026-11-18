"""Checks and powder scattering for the actual finite model in the viewer.

Finite coordinates remain finite. No periodic cell, missing atom, or crystal
identity is inferred. Continuation placeholders are omitted from atom physics.
"""
from __future__ import annotations

from collections import Counter
import hashlib
import json
import math

import numpy as np

from .cofdraw_characterization import atomic_form_factor, _formula
from .cofdraw_science import PT, _finite
from .cofdraw_xray_factors import DATA_SHA256, SOURCE_VERSION

ELEMENTS = frozenset(PT.GetElementSymbol(n) for n in range(1, 119))


def finite_model(model):
    if not isinstance(model, dict) or model.get('periodic') is not False:
        raise ValueError('Provide the explicitly nonperiodic viewer model; periodic structures use the CIF/model workflow')
    rows = model.get('atoms')
    if not isinstance(rows, list) or not 1 <= len(rows) <= 12000:
        raise ValueError('Finite model requires 1–12000 displayed sites')
    atoms, positions, indices, omitted = [], [], [], []
    for i, row in enumerate(rows):
        if not isinstance(row, dict):
            raise ValueError('Each atom must be an object')
        xyz = [_finite(row.get(k), k, -10000, 10000) for k in ('x', 'y', 'z')]
        element = row.get('element')
        if row.get('continuation_site') is True and element not in ('*', 'X', 'R'):
            raise ValueError('A real element cannot be marked as a continuation placeholder')
        if element in ('*', 'X', 'R'):
            omitted.append(i)
            continue
        if not isinstance(element, str) or element not in ELEMENTS:
            raise ValueError('Unknown finite-model element')
        atoms.append(element); positions.append(xyz); indices.append(i)
    if not atoms:
        raise ValueError('No real atoms are present')
    bonds = []
    supplied_bonds=model.get('bonds', [])
    if not isinstance(supplied_bonds,list) or len(supplied_bonds)>36000:
        raise ValueError('Finite bonds must be a list of at most 36000 edges')
    for bond in supplied_bonds:
        if not isinstance(bond, dict):
            raise ValueError('Finite bonds must contain a and b atom indices')
        a, b = bond.get('a'), bond.get('b')
        if type(a) is not int or type(b) is not int or a == b or not (0 <= a < len(rows) and 0 <= b < len(rows)):
            raise ValueError('Invalid finite-model bond indices')
        bonds.append((min(a, b), max(a, b)))
    if len(bonds) > 36000 or len(set(bonds)) != len(bonds):
        raise ValueError('Duplicate bonds or finite-model bond limit exceeded')
    return atoms, np.asarray(positions), indices, omitted, bonds


def inspect_finite(model):
    atoms, positions, indices, omitted, bonds = finite_model(model)
    counts = Counter(atoms)
    index_map = {old: i for i, old in enumerate(indices)}
    real_bonds = [(index_map[a], index_map[b]) for a, b in bonds if a in index_map and b in index_map]
    parents = list(range(len(atoms)))
    def root(i):
        while parents[i] != i:
            parents[i] = parents[parents[i]]; i = parents[i]
        return i
    for a, b in real_bonds:
        parents[root(a)] = root(b)
    bond_lengths = [float(np.linalg.norm(positions[a] - positions[b])) for a, b in real_bonds]
    contacts = {'available': False, 'reason': 'Contact check is limited to 2000 real atoms; use fewer layers.'}
    if len(atoms) <= 2000:
        a, b = np.triu_indices(len(atoms), 1)
        distance = np.linalg.norm(positions[a] - positions[b], axis=1)
        radii = np.asarray([PT.GetRcovalent(e) for e in atoms])
        short = np.flatnonzero(distance < .65 * (radii[a] + radii[b]))
        contacts = {'available': True, 'minimum_distance_angstrom': float(distance.min()) if len(distance) else None,
                    'short_contact_count': len(short), 'short_contacts': [
                        {'atoms': [indices[a[k]], indices[b[k]]], 'distance_angstrom': float(distance[k])} for k in short[:50]],
                    'rule': 'distance < 0.65 × sum of RDKit covalent radii, including supplied bonds; first 50 listed'}
    return {'kind': 'finite_viewer_check', 'periodic': False, 'atom_count': len(atoms),
            'displayed_atom_inventory': _formula(counts), 'elements': dict(sorted(counts.items())),
            'continuation_site_count': len(omitted), 'known_real_bond_count': len(real_bonds),
            'known_graph_components': len({root(i) for i in range(len(atoms))}),
            'extent_angstrom': np.ptp(positions, axis=0).tolist(),
            'bond_length_range_angstrom': [min(bond_lengths), max(bond_lengths)] if bond_lengths else None,
            'contact_checks': contacts, 'model_sha256': model_hash(model),
            'scope': 'Checks the displayed finite coordinates and supplied bonds. Atom inventory excludes continuation placeholders and does not add missing hydrogens. Components are finite graph components, not periodic dimensionality.',
            'unavailable': {'cell': 'No periodic cell supplied', 'density': 'No periodic cell supplied',
                            'periodic_pore_properties': 'Require an explicit periodic atomic model'},
            'warnings': ['Scaled drawing layers are illustrative coordinates, not a relaxed 3D structure.']}


def model_hash(model):
    return hashlib.sha256(json.dumps(model, sort_keys=True, allow_nan=False).encode()).hexdigest()


def finite_powder(model, *, wavelength=1.5406, two_theta_min=2., two_theta_max=50., step=.1):
    """Exact pair summation of the independent-atom Debye scattering equation."""
    atoms, positions, _, omitted, _ = finite_model(model)
    wavelength = _finite(wavelength, 'wavelength', .3, 3)
    lo = _finite(two_theta_min, 'two_theta_min', .5, 100)
    hi = _finite(two_theta_max, 'two_theta_max', 1, 150)
    step = _finite(step, 'step', .01, .2)
    if hi <= lo:
        raise ValueError('Maximum angle must exceed minimum angle')
    pairs = len(atoms) * (len(atoms) - 1) // 2
    points = math.ceil((hi - lo) / step) + 1
    if len(atoms) > 2000 or points > 15001 or pairs * points > 300_000_000:
        raise ValueError('Finite scattering pair-sum limit exceeded. Choose fewer layers, a narrower angle range, or a larger angular step.')
    angles = np.linspace(lo, hi, points)
    stol = np.sin(np.radians(angles / 2)) / wavelength
    q = 4 * np.pi * stol
    unique = sorted(set(atoms)); codes = np.asarray([unique.index(e) for e in atoms])
    factors = np.asarray([atomic_form_factor(e, stol**2) for e in unique])
    intensities = (factors[codes]**2).sum(axis=0)
    a, b = np.triu_indices(len(atoms), 1)
    distances = np.linalg.norm(positions[a] - positions[b], axis=1)
    # Small pair chunks bound temporary memory independently of model size.
    for start in range(0, pairs, 512):
        stop = min(start + 512, pairs)
        weights = factors[codes[a[start:stop]]] * factors[codes[b[start:stop]]]
        intensities += 2 * (weights * np.sinc(distances[start:stop, None] * q[None, :] / np.pi)).sum(axis=0)
    if np.min(intensities) < -1e-7 * max(1., float(np.max(intensities))):
        raise ValueError('Numerically invalid negative scattering intensity')
    intensities = np.maximum(intensities, 0)
    peak = float(intensities.max())
    relative = intensities / peak * 100 if peak else intensities
    csv = 'two_theta_degrees,Q_inverse_angstrom,intensity_electron_squared,relative_intensity\n'
    csv += ''.join(f'{x:.10g},{v:.10g},{y:.10g},{r:.10g}\n' for x,v,y,r in zip(angles,q,intensities,relative))
    return {'kind': 'finite_debye_powder_scattering', 'periodic': False,
            'profile': {'two_theta_degrees': angles.tolist(), 'intensity': relative.tolist(),
                        'absolute_intensity_electron_squared': intensities.tolist(), 'Q_inverse_angstrom': q.tolist()},
            'parameters': {'wavelength_angstrom': wavelength, 'two_theta_min': lo, 'two_theta_max': hi,
                           'requested_step_degrees': step, 'effective_step_degrees': float(angles[1]-angles[0]),
                           'atom_count': len(atoms), 'omitted_continuation_sites': len(omitted)},
            'method': {'equation': 'I(Q) = sum_i f_i(Q)^2 + 2 sum_(i<j) f_i(Q) f_j(Q) sin(Q r_ij)/(Q r_ij)',
                       'pair_sum': 'All real-atom pairs; no distance histogram or periodic images',
                       'scattering': SOURCE_VERSION, 'coefficient_sha256': DATA_SHA256,
                       'reference': 'https://doi.org/10.1107/S1600576723002121'},
            'model_sha256': model_hash(model), 'exports': {'csv': csv},
            'scope': 'Orientationally averaged scattering of the displayed finite atoms. No unit cell, reciprocal-lattice peaks, or hkl assignment is inferred.',
            'warnings': ['Finite drawing coordinates are illustrative and may omit hydrogens; continuation placeholders scatter no intensity.',
                         'This is finite-cluster Debye scattering, not a periodic crystal PXRD prediction.',
                         'Neutral independent-atom form factors; no thermal motion, instrumental convolution, polarization, absorption, anomalous or incoherent scattering.']}
