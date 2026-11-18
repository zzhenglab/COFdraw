"""Analytic geometry and periodic-transport tests for deterministic Level 1."""
import itertools
import json
import math
from copy import deepcopy

import numpy as np
import pytest

from cof_vision.cofdraw_porosity import (geometric_porosity, _WindingUnion, _edges,
    _threshold_components, _segments_clear, _segment_image_shifts)
from cof_vision.cofdraw_science import AVOGADRO, PT, _PeriodicDistances


def sphere_cell(length=10):
    return {'cell': np.diag([length]*3).tolist(), 'atoms': [{'element': 'C', 'fractional': [.5]*3}]}


def partitioned_cell(axes=(0, 1, 2), skew=False):
    # Dense periodically repeated walls close every listed cell-normal axis.
    points = sorted({p for p in itertools.product(range(4), repeat=3) if any(p[a] == 0 for a in axes)})
    return {'cell': [[6, 0, 0], [1 if skew else 0, 6, 0], [0, 0, 6]],
            'atoms': [{'element': 'C', 'fractional': (np.asarray(p)/4).tolist()} for p in points]}


def test_simple_cubic_sphere_has_analytic_lcd_pld_inside_grid_bounds_and_area_units():
    result = geometric_porosity(sphere_cell(), probe_radius=1, grid_spacing=1, surface_samples=16, refinement=True)
    radius = PT.GetRvdw(6)
    exact_lcd = 2*(math.sqrt(3)*5-radius)
    exact_pld = 2*(math.sqrt(2)*5-radius)
    for name, exact in [('lcd', exact_lcd), ('pld', exact_pld)]:
        low, high = result[f'{name}_bounds_angstrom']
        assert low <= exact <= high+1e-10
        refined = result['convergence']['refined'][f'{name}_bounds_angstrom']
        assert refined[0] <= exact <= refined[1]+1e-10
        assert refined[1]-refined[0] < high-low
    area = 4*math.pi*(radius+1)**2
    assert result['transport_accessible_surface_angstrom2'] == pytest.approx(area)
    assert result['local_probe_exposed_surface_angstrom2'] == pytest.approx(area)
    mass = PT.GetAtomicWeight(6)
    assert result['density_g_cm3'] == pytest.approx(mass/AVOGADRO/1e-21)
    assert result['transport_accessible_surface_m2_g'] == pytest.approx(area*1e-20*AVOGADRO/mass)
    assert result['transport_accessible_probe_center_volume_cm3_g'] == pytest.approx(
        result['transport_accessible_probe_center_fraction']/result['density_g_cm3'])
    for bounds, exact in [('void_fraction_bounds', 1-4*math.pi*radius**3/3000),
                          ('transport_accessible_probe_center_fraction_bounds', 1-4*math.pi*(radius+1)**3/3000)]:
        assert result[bounds][0] <= exact <= result[bounds][1]
    assert result['connectivity']['maximum_periodic_rank'] == 3
    assert result['bet_surface_area_m2_g'] is None
    json.dumps(result, allow_nan=False)


def test_closed_cavity_has_local_area_but_no_transport_accessible_area_or_volume():
    result = geometric_porosity(partitioned_cell(), probe_radius=.4, grid_spacing=.5, surface_samples=32)
    assert result['lcd_angstrom'] > 1
    assert result['pld_angstrom'] == 0
    assert result['pld_bounds_angstrom'][1] == 0
    assert result['probe_center_free_fraction'] > 0
    assert result['local_probe_exposed_surface_angstrom2'] > 0
    assert result['transport_accessible_probe_center_fraction'] == 0
    assert result['transport_accessible_surface_angstrom2'] == 0
    assert result['connectivity']['maximum_periodic_rank'] == 0


@pytest.mark.parametrize('skew', [False, True])
def test_periodic_channel_retains_one_dimensional_transport_with_skew_cell(skew):
    structure = partitioned_cell((0, 1), skew=skew)
    result = geometric_porosity(structure, probe_radius=.4, grid_spacing=.6, surface_samples=16)
    assert result['pld_angstrom'] > .8
    assert result['transport_accessible_surface_angstrom2'] > 0
    assert result['transport_accessible_probe_center_fraction'] > 0
    assert result['connectivity']['maximum_periodic_rank'] == 1
    assert all(component['percolating_axes'] == ['c'] for component in result['connectivity']['components'] if component['periodic_rank'])
    # No stochastic state, mutation, or implicit hydrogens in the calculation.
    before = deepcopy(structure)
    assert geometric_porosity(structure, probe_radius=.4, grid_spacing=.6, surface_samples=16) == result
    assert structure == before


def test_periodic_face_touching_is_not_winding_and_diagonal_winding_is_rank_one():
    shape = np.array([6, 6, 6])
    clearance = np.full(math.prod(shape), -1.)
    # A small pocket straddles a boundary but has no noncontractible cycle.
    for xyz in [(0, 2, 2), (5, 2, 2)]:
        clearance[np.ravel_multi_index(xyz, shape)] = 1
    edges = [(a,b,t,np.minimum(clearance[a],clearance[b])) for a,b,t,_ in _edges(shape)]
    accessible, components = _threshold_components(clearance, edges, .5)
    assert not accessible.any()
    assert len(components) == 1 and components[0]['periodic_rank'] == 0
    union = _WindingUnion(3)
    assert not union.union(0, 1, np.array([0, 0, 0]))
    assert not union.union(1, 2, np.array([0, 0, 0]))
    assert union.union(2, 0, np.array([1, 1, 0]))
    root, _ = union.find(0)
    assert len(union.basis[root]) == 1
    union.union(2, 0, np.array([2, 2, 0]))
    assert len(union.basis[root]) == 1


def test_exact_segment_check_catches_midsegment_collision_and_periodic_image():
    distance = _PeriodicDistances(np.diag([10.,10.,10.]), np.array([[.5,.5,.5]]), np.array([1.7]))
    starts = np.array([[.2,.5,.5], [.2,.8,.5], [1.2,.5,.5]])
    ends = np.array([[.8,.5,.5], [.8,.8,.5], [1.8,.5,.5]])
    assert np.all(distance.clearance(starts) > 0)
    assert np.all(distance.clearance(ends) > 0)
    shifts, bounds = _segment_image_shifts(distance, 0, 6)
    assert _segments_clear(distance, starts, ends, 0, shifts, bounds).tolist() == [False, True, False]


def test_triclinic_segment_image_bound_matches_independent_wide_enumeration():
    rng = np.random.default_rng(34)
    cell = np.array([[8.,0,0], [3,6,0], [1,1,5]])
    atoms = np.array([[.1,.2,.3], [.8,.7,.4], [.4,.8,.9]])
    radii = np.array([1.7,1.55,1.2])
    distance = _PeriodicDistances(cell, atoms, radii)
    starts = rng.random((32,3)) + [2,-1,3]
    ends = starts + rng.uniform(-.3,.3,(32,3))
    maximum = float(np.linalg.norm((ends-starts)@cell,axis=1).max())
    shifts, bounds = _segment_image_shifts(distance, .6, maximum)
    actual = _segments_clear(distance, starts, ends, .6, shifts, bounds)
    expected = []
    images = np.asarray(list(itertools.product(range(-4,5),repeat=3)))
    for a,b in zip(starts,ends):
        origin = np.floor((a+b)/2)
        sphere_centers = (atoms[:,None] + images[None] + origin) @ cell
        vector = (b-a) @ cell
        delta = a@cell - sphere_centers
        t = np.clip(-np.einsum('...i,i->...',delta,vector)/(vector@vector),0,1)
        minimum = np.linalg.norm(delta+t[...,None]*vector,axis=2)-radii[:,None]
        expected.append(bool(minimum.min() >= .6-1e-7))
    assert actual.tolist() == expected


def test_work_caps_are_explicit_and_invalid_parameters_do_not_silently_coerce():
    result = geometric_porosity(sphere_cell(40), grid_spacing=.1, surface_samples=16)
    assert result['method']['grid_points'] <= result['method']['grid_point_limit'] <= 65536
    assert max(result['method']['grid_shape']) <= 96
    assert any('coarsened' in warning for warning in result['warnings'])
    assert max(result['method']['effective_axis_steps_angstrom']) > .1
    for keywords in [{'grid_spacing': 0}, {'probe_radius': -1}, {'surface_samples': True}, {'refinement': 'yes'}]:
        with pytest.raises(ValueError):
            geometric_porosity(sphere_cell(), **keywords)
