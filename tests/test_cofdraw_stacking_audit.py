"""Independent geometric controls for stacking skewed and wrapped layers."""
from copy import deepcopy
import numpy as np
import pytest

from cof_vision.cofdraw_stacking import periodic_stacking


def wrapped_graphene():
    # A graphene honeycomb primitive cell, slightly corrugated across z=0.
    # The second atom belongs at z=1.01 relative to the first at z=.99.
    return {'periodic': True,
            'cell': [[2.46, 0, 0], [1.23, 2.46*np.sqrt(3)/2, 0], [1, .4, 4]],
            'atoms': [{'element': 'C', 'fractional': [.2, .2, .99], 'charge': -.1},
                      {'element': 'C', 'fractional': [.2+1/3, .2+1/3, .01], 'charge': .1}]}


@pytest.mark.parametrize('mode', ['AA', 'AB', 'ABC', 'slipped'])
def test_skewed_boundary_layer_retains_cartesian_geometry_and_charges(mode):
    model = wrapped_graphene(); before = deepcopy(model)
    result = periodic_stacking(model, stacking=mode, spacing=3.6, slip=[.7, -.2])
    out = result['structure']; cell = np.array(out['cell'])
    # Independent Cartesian oracle: use the explicitly identified neighboring
    # image of atom 1, then translate both atoms by the same 4 A normal offset.
    original_cell = np.array(model['cell'])
    physical = np.array([[.2,.2,.99], [.2+1/3,.2+1/3,1.01]]) @ original_cell
    physical -= [0,0,4]
    expected_fraction = (physical @ np.linalg.inv(cell)) % 1
    actual_fraction = np.array([a['fractional'] for a in out['atoms'][:2]])
    assert actual_fraction == pytest.approx(expected_fraction, abs=1e-12)
    assert [a['charge'] for a in out['atoms'][:2]] == [-.1,.1]
    assert model == before
    evidence = out['provenance']['layer_validation']
    assert evidence['periodic_rank'] == 2 and evidence['cross_layer_contacts'] == 0
    assert evidence['experimentally_assigned_connectivity'] is False


def test_global_cartesian_rotation_rotates_stack_and_slip_consistently():
    model = wrapped_graphene(); rotated = deepcopy(model)
    t = np.radians(39)
    rotation = np.array([[1,0,0], [0,np.cos(t),-np.sin(t)], [0,np.sin(t),np.cos(t)]])
    rotated['cell'] = (np.array(model['cell']) @ rotation).tolist()
    a = periodic_stacking(model, stacking='slipped', slip=[.7,-.2])['structure']
    b = periodic_stacking(rotated, stacking='slipped', slip=[.7,-.2])['structure']
    assert np.array(b['cell']) == pytest.approx(np.array(a['cell']) @ rotation, abs=1e-12)
    assert np.array([r['fractional'] for r in b['atoms']]) == pytest.approx(np.array([r['fractional'] for r in a['atoms']]), abs=1e-12)


def test_flat_looking_3d_carbon_cell_cannot_be_restacked():
    model = {'cell': [[1.5,0,0], [0,1.5,0], [0,0,1.5]],
             'atoms': [{'element':'C', 'fractional':[0,0,0]}]}
    # Thickness is zero but self-image bonds run in all three directions.
    with pytest.raises(ValueError, match='different layers'):
        periodic_stacking(model)


@pytest.mark.parametrize('cell', [[[20,0,0],[0,20,0],[0,0,4]],
                                [[1.5,0,0],[0,20,0],[0,0,4]]])
def test_molecular_or_chain_cell_is_not_certified_as_a_2d_layer(cell):
    with pytest.raises(ValueError, match='2D periodic covalent layer'):
        periodic_stacking({'cell':cell, 'atoms':[{'element':'C','fractional':[0,0,0]}]})


def test_too_short_vertical_contacts_are_rejected_even_when_not_inferred_as_bonds():
    model = {'cell': [[1.5,0,0], [0,1.5,0], [0,0,.6]],
             'atoms': [{'element':'C','fractional':[0,0,0]}]}
    with pytest.raises(ValueError, match='implausibly short'):
        periodic_stacking(model)
