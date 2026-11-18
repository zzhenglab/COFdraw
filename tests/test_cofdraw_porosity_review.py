"""Independent analytic lattice checks beyond the implementation's cubic fixture."""
import math

import numpy as np
import pytest

from cof_vision.cofdraw_porosity import geometric_porosity
from cof_vision.cofdraw_science import PT


def fcc_spheres(shear=False):
    cell = np.diag([8.0] * 3)
    positions = np.asarray([[0, 0, 0], [0, .5, .5], [.5, 0, .5], [.5, .5, 0]])
    if shear:
        # Integer determinant-one basis change: exactly the same infinite sites.
        transform = np.array([[1, 0, 0], [1, 1, 0], [0, 0, 1]])
        cell = transform @ cell
        positions = positions @ np.linalg.inv(transform)
    return {"cell": cell.tolist(), "atoms": [{"element": "C", "fractional": p.tolist()} for p in positions]}


@pytest.mark.parametrize("shear", [False, True])
def test_fcc_analytic_windows_and_cavities_remain_enclosed_after_lattice_rebasis(shear):
    result = geometric_porosity(fcc_spheres(shear), probe_radius=1, grid_spacing=.6, surface_samples=16)
    radius = PT.GetRvdw(6)
    # FCC void network has octahedral cavities (nearest center a/2) joined
    # through equilateral triangle windows with circumradius a/sqrt(6).
    # Subtract the hard-sphere radius, then double for diameters.
    for name, exact in (("lcd", 2 * (4 - radius)), ("pld", 2 * (8 / math.sqrt(6) - radius))):
        low, high = result[name + "_bounds_angstrom"]
        assert low <= exact <= high
    # At this probe all enlarged spheres are disjoint (2*(r+1)<a/sqrt(2)),
    # and their complement is transport-connected in three dimensions.
    free = 1 - 4 * (4 * math.pi / 3 * (radius + 1) ** 3) / 8 ** 3
    assert 2 * (radius + 1) < 8 / math.sqrt(2)
    low, high = result["transport_accessible_probe_center_fraction_bounds"]
    assert low <= free <= high
    # The axis-aligned coarse grid can miss a real window; a lower estimate
    # of zero therefore cannot certify that the physical void is closed.
    assert result["pld_bounds_angstrom"][1] > 2
    assert result["connectivity"]["status"] == ("certified_channel" if shear else "unresolved_at_this_resolution")
    if not shear:
        assert result["transport_accessible_probe_center_fraction"] == 0
        assert any("unresolved" in message.lower() for message in result["warnings"])


def test_orthorhombic_path_uses_easiest_transport_axis_not_smallest_cavity():
    structure = {"cell": np.diag([6., 8., 10.]).tolist(),
                 "atoms": [{"element": "C", "fractional": [0, 0, 0]}]}
    result = geometric_porosity(structure, probe_radius=.5, grid_spacing=.7, surface_samples=16)
    radius = PT.GetRvdw(6)
    # Translation along the short a vector leaves b/2,c/2 clearances at its
    # bottleneck; every a-crossing must encounter such a face. Paths along
    # other axes have smaller bottlenecks. LCD occurs at the cell center.
    exact_pld = 2 * (math.hypot(4, 5) - radius)
    exact_lcd = 2 * (math.sqrt(3 ** 2 + 4 ** 2 + 5 ** 2) - radius)
    assert result["pld_bounds_angstrom"][0] <= exact_pld <= result["pld_bounds_angstrom"][1]
    assert result["lcd_bounds_angstrom"][0] <= exact_lcd <= result["lcd_bounds_angstrom"][1]
    assert result["connectivity"]["maximum_periodic_rank"] == 3
