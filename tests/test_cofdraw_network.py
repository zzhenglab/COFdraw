"""Periodic graph invariants and real template geometry, not mocked chemistry."""
import itertools
import json

import numpy as np
import pytest
from rdkit import Chem

from cof_vision import cofdraw_network as network
from cof_vision.cofdraw_science import import_cif, structure_exports


def structure(cell, positions, elements=None):
    elements = elements or ['C']*len(positions)
    return {'periodic': True, 'cell': cell,
            'atoms': [{'element': e, 'fractional': list(p)} for e, p in zip(elements, positions)]}


@pytest.mark.parametrize('dimension', [0, 1, 2, 3])
def test_rank_counts_periodic_cycle_translations_and_self_image_neighbors(dimension):
    cell = np.diag([1.54 if axis < dimension else 8 for axis in range(3)])
    result = network.inspect_periodic_network(structure(cell.tolist(), [[.4, .3, .2]]))
    assert result['dimensions'] == [dimension]
    assert result['components'][0]['periodic_rank'] == dimension
    assert result['edge_count'] == dimension
    assert result['coordination']['per_atom'][0]['coordination'] == 2*dimension
    assert all(row['a'] == row['b'] and any(row['translation']) for row in result['edges'])
    json.dumps(result, allow_nan=False)


def test_boundary_crossing_molecule_is_finite_not_one_dimensional():
    result = network.inspect_periodic_network(structure(np.diag([15, 8, 8]).tolist(), [[.95, .3, .2], [.05, .3, .2]]))
    assert result['dimensions'] == [0]
    assert result['quotient_component_count'] == 1
    assert result['edge_count'] == 1
    assert result['edges'][0]['translation'] == [1, 0, 0]
    assert result['edges'][0]['distance_angstrom'] == pytest.approx(1.5)


def test_skew_cell_retains_connections_beyond_neighboring_27_images():
    # (-2,1,0) is the short direction: checking only t in [-1,1]^3 misses it.
    cell = [[6, 0, 0], [12, 1.5, 0], [0, 0, 10]]
    result = network.inspect_periodic_network(structure(cell, [[.2, .3, .4]]))
    assert result['dimensions'] == [1]
    assert result['coordination']['mean'] == 2
    assert result['edges'][0]['translation'] == [2, -1, 0]
    assert result['edges'][0]['distance_angstrom'] == pytest.approx(1.5)
    assert result['method']['translation_bounds'][0] >= 2


def test_skew_contacts_match_wide_brute_force_enumeration():
    rng = np.random.default_rng(423)
    cell = np.array([[4.2, 0, 0], [3.9, 2.4, 0], [1.3, .6, 5]])
    positions = rng.random((8, 3))
    elements = ['C', 'N', 'O', 'H']*2
    result = network.inspect_periodic_network(structure(cell.tolist(), positions, elements))
    table = Chem.GetPeriodicTable()
    radii = [table.GetRcovalent(table.GetAtomicNumber(e)) for e in elements]
    expected = {}
    for a in range(len(positions)):
        for b in range(a, len(positions)):
            for shift in itertools.product(range(-4, 5), repeat=3):
                if a == b and not next((v > 0 for v in shift if v), False):
                    continue
                distance = np.linalg.norm((positions[b]+shift-positions[a]) @ cell)
                radius_sum = radii[a]+radii[b]
                if max(.35, .65*radius_sum) <= distance <= 1.15*radius_sum+.15+1e-10:
                    expected[(a, b, shift)] = distance
    actual = {(r['a'], r['b'], tuple(r['translation'])): r['distance_angstrom'] for r in result['edges']}
    assert actual.keys() == expected.keys()
    for key, value in expected.items():
        assert actual[key] == pytest.approx(value, abs=1e-10)


def test_integer_rank_does_not_confuse_long_nearly_parallel_cycle_vectors():
    basis = []
    network._append_independent(basis, [1_000_000_000, 1, 0])
    network._append_independent(basis, [999_999_999, 1, 0])
    network._append_independent(basis, [0, 0, 1])
    assert len(basis) == 3


def test_dimension_survives_unimodular_cell_basis_change_and_supercell():
    cell = np.diag([1.54, 1.54, 8])
    change = np.array([[1, 0, 0], [2, 1, 0], [0, 0, 1]])
    base = structure(cell.tolist(), [[.25, .35, .5]])
    transformed = structure((change @ cell).tolist(), np.array([[.25, .35, .5]]) @ np.linalg.inv(change))
    a = network.inspect_periodic_network(base)
    b = network.inspect_periodic_network(transformed)
    expanded = import_cif(structure_exports(base, repeats=(3, 2, 1))['exports']['cif'])
    c = network.inspect_periodic_network(expanded)
    assert a['dimensions'] == b['dimensions'] == c['dimensions'] == [2]
    assert a['coordination']['mean'] == b['coordination']['mean'] == c['coordination']['mean'] == 4
    assert a['edge_count']*6 == c['edge_count']


def test_short_contacts_are_reported_excluded_and_not_silently_bonded():
    result = network.inspect_periodic_network(structure(np.diag([10]*3).tolist(), [[.1, .1, .1], [.14, .1, .1]]))
    assert result['short_contact_count'] == 1
    assert result['short_contacts'][0]['distance_angstrom'] == pytest.approx(.4)
    assert result['edge_count'] == 0
    assert result['isolated_atom_indices'] == [0, 1]
    assert any('implausibly short' in warning for warning in result['warnings'])


def test_cutoff_sensitivity_and_explicit_hydrogen_coordination_are_visible():
    sensitive = network.inspect_periodic_network(structure(np.diag([10]*3).tolist(), [[.1, .1, .1], [.29, .1, .1]]))
    assert sensitive['threshold_sensitivity']['classification_changed']
    assert sensitive['borderline_contact_count'] == 1
    result = network.inspect_periodic_network(structure(np.diag([10]*3).tolist(), [[.5, .5, .5], [.4, .5, .5], [.6, .5, .5]], ['H', 'C', 'C']))
    assert result['explicit_hydrogen_count'] == 1
    assert result['coordination']['unusual_atoms'] == [{'atom_index': 0, 'element': 'H', 'coordination': 2, 'reason': 'above_common_organic_neighbor_count'}]


@pytest.mark.parametrize('kwargs', [{'bond_scale': float('nan')}, {'bond_scale': True}, {'bond_scale': 1.6}, {'tolerance_angstrom': -.1}, {'tolerance_angstrom': float('inf')}])
def test_invalid_cutoff_parameters_fail_explicitly(kwargs):
    with pytest.raises(ValueError):
        network.inspect_periodic_network(structure(np.diag([8]*3).tolist(), [[.5]*3]), **kwargs)


def test_work_limit_rejects_instead_of_truncating_periodic_search(monkeypatch):
    monkeypatch.setattr(network, 'MAX_IMAGE_PAIRS', 1)
    with pytest.raises(ValueError, match='work bound'):
        network.inspect_periodic_network(structure(np.diag([1.54]*3).tolist(), [[.5]*3]))
