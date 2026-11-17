"""Analytic geometry, chemical invariants, and honest-method regression checks."""
import itertools
import json
import math

import numpy as np
import pytest
from rdkit import Chem

from cof_vision.cofdraw_science import (
    AVOGADRO, _PeriodicDistances, _percolation, cell_from_parameters,
    fit_bet, import_cif, monomer_analysis, normalize_structure,
    patch_3d, pore_analysis, structure_exports,
)


def carbon_cell(length=10):
    return {"cell": cell_from_parameters(length, length, length),
            "atoms": [{"element": "C", "fractional": [.5, .5, .5]}]}


def test_monomer_real_descriptors_3d_stereo_and_reproducibility():
    smiles = "N[C@@H](C)C(=O)O"
    first = monomer_analysis(smiles, conformers=2, seed=17)
    second = monomer_analysis(smiles, conformers=2, seed=17)
    assert first["descriptors"]["formula"] == "C3H7NO2"
    assert first["descriptors"]["molecular_weight_g_mol"] == pytest.approx(89.094, abs=.01)
    assert first["method"]["force_field"] == "MMFF94s"
    assert first["partial_charges"]["available"]
    assert first["partial_charges"]["sum"] == pytest.approx(0, abs=1e-7)
    assert len(first["partial_charges"]["values"]) == 13
    assert first["selected_conformer"]["converged"]
    assert first["viewer"]["atoms"] == second["viewer"]["atoms"]
    assert first["exports"]["xyz"].splitlines()[0] == "13"
    restored = Chem.MolFromMolBlock(first["exports"]["sdf"], removeHs=True)
    assert restored.GetConformer().Is3D()
    assert Chem.MolToSmiles(restored, isomericSmiles=True) == Chem.MolToSmiles(Chem.MolFromSmiles(smiles), isomericSmiles=True)
    json.dumps(first, allow_nan=False)


@pytest.mark.parametrize("smiles", ["CC.CC", "*CC", "not chemistry", ""])
def test_monomer_rejects_disconnected_or_open_periodic_graph(smiles):
    with pytest.raises(ValueError):
        monomer_analysis(smiles)


def test_monomer_optimizer_limit_remains_visible():
    result = monomer_analysis("CCCCCCO", conformers=1, max_iterations=1)
    assert result["selected_conformer"]["optimizer_status"] == 1
    assert not result["selected_conformer"]["converged"]
    assert any("not force-field converged" in note for note in result["warnings"])


def test_p1_cif_supercell_roundtrip_preserves_atoms_metric_and_occupancy():
    structure = {"cell": cell_from_parameters(12, 12, 3.4, 90, 90, 120),
                 "atoms": [{"element": "C", "fractional": [.1, .2, .3]},
                           {"element": "N", "fractional": [.6, .7, .8]}]}
    exported = structure_exports(structure, repeats=(2, 1, 2))
    restored = import_cif(exported["exports"]["cif"])
    assert len(restored["atoms"]) == 8
    assert restored["cell_volume_angstrom3"] == pytest.approx(4 * np.linalg.det(structure["cell"]), rel=1e-9)
    assert all(a["occupancy"] == 1 for a in restored["atoms"])
    assert exported["exports"]["xyz"].splitlines()[0] == "8"
    assert 'pbc="T T T"' in exported["exports"]["xyz"]
    pdb_rows = [line for line in exported["exports"]["pdb"].splitlines() if line.startswith("HETATM")]
    assert len(pdb_rows) == 8
    for row, atom in zip(pdb_rows, exported["viewer"]["atoms"]):
        assert [float(row[a:b]) for a, b in ((30, 38), (38, 46), (46, 54))] == pytest.approx([atom[k] for k in ("x", "y", "z")], abs=.001)


def symmetry_cif(operations="'x,y,z'\n'-x,-y,-z'", occupancy="1"):
    return f"""data_test
_cell_length_a 10.000(5)
_cell_length_b 11
_cell_length_c 12
_cell_angle_alpha 90
_cell_angle_beta 90
_cell_angle_gamma 90
_symmetry_space_group_name_H-M 'P -1'
loop_
_symmetry_equiv_pos_as_xyz
{operations}
loop_
_atom_site_label
_atom_site_type_symbol
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
_atom_site_occupancy
C1 C .1 .2 .3 {occupancy}
O1 O 0 0 0 1
"""


def test_cif_symmetry_expands_general_sites_and_deduplicates_special_sites():
    parsed = import_cif(symmetry_cif())
    assert len(parsed["atoms"]) == 3
    assert parsed["provenance"]["asymmetric_sites"] == 2
    assert parsed["provenance"]["supplied_symmetry_operations"] == 2
    assert sorted(a["element"] for a in parsed["atoms"]) == ["C", "C", "O"]


@pytest.mark.parametrize("text,reason", [
    (symmetry_cif(occupancy="0.5"), "partial occupancy"),
    (symmetry_cif(operations="'x+y,y,z'"), "inconsistent"),
    (symmetry_cif(operations="'-x,-y,-z'"), "identity"),
    (symmetry_cif(operations="'x,y,z'\n'x+1/3,y,z'"), "closed"),
    (symmetry_cif(operations="'__import__(x),y,z'"), "CIF"),
    (symmetry_cif().replace("loop_\n_symmetry_equiv_pos_as_xyz\n'x,y,z'\n'-x,-y,-z'\n", ""), "explicit symmetry"),
])
def test_cif_rejects_unresolved_or_invalid_structure(text, reason):
    with pytest.raises(ValueError, match=reason):
        import_cif(text)


def test_cif_does_not_guess_missing_spacegroup_or_duplicate_data_blocks():
    exported = structure_exports(carbon_cell())["exports"]["cif"]
    with pytest.raises(ValueError, match="one CIF"):
        import_cif(exported + exported)
    no_symmetry = exported.replace("_space_group_name_H-M_alt 'P 1'\n_space_group_IT_number 1\nloop_\n_space_group_symop_operation_xyz\n'x,y,z'\n", "")
    with pytest.raises(ValueError, match="declare P1"):
        import_cif(no_symmetry)


def test_structure_rejects_dummy_sites_partial_occupancy_and_invalid_cell():
    for update in ({"element": "*"}, {"occupancy": .8}):
        value = carbon_cell(); value["atoms"][0].update(update)
        with pytest.raises(ValueError):
            normalize_structure(value)
    value = carbon_cell(); value["cell"][2] = value["cell"][0]
    with pytest.raises(ValueError):
        normalize_structure(value)
    with pytest.raises(ValueError, match="periodic"):
        normalize_structure({"periodic": False, **{k: v for k, v in carbon_cell().items()}})


def test_triclinic_minimum_images_match_explicit_wide_image_enumeration():
    cell = np.array(cell_from_parameters(12, 11, 5, 85, 90, 120))
    atoms = np.array([[.05, .2, .4], [.95, .82, .6]])
    radii = np.array([1.7, 1.55])
    points = np.random.default_rng(11).random((35, 3))
    result = _PeriodicDistances(cell, atoms, radii).clearance(points)
    images = np.array(list(itertools.product(range(-4, 5), repeat=3)))
    reference = []
    for point in points:
        distances = np.linalg.norm((point[None, None, :] - atoms[:, None, :] - images[None, :, :]) @ cell, axis=2)
        reference.append(np.min(distances - radii[:, None]))
    assert result == pytest.approx(reference, abs=1e-10)
    translated = _PeriodicDistances(cell, atoms, radii).clearance(points + [2, -3, 5])
    assert translated == pytest.approx(result, abs=1e-10)


def test_hard_sphere_analytic_volume_surface_density_units_and_repeatability():
    structure = carbon_cell()
    result = pore_analysis(structure, probe_radius=1, samples=24000, surface_samples=96, grid_size=12, seed=9)
    radius = Chem.GetPeriodicTable().GetRvdw(6)
    expected_void = 1 - 4 * math.pi * radius**3 / 3 / 1000
    expected_free = 1 - 4 * math.pi * (radius + 1)**3 / 3 / 1000
    assert result["void_fraction"] == pytest.approx(expected_void, abs=.006)
    assert result["probe_center_free_fraction"] == pytest.approx(expected_free, abs=.008)
    expected_area = 4 * math.pi * (radius + 1)**2
    assert result["local_probe_exposed_surface_angstrom2"] == pytest.approx(expected_area)
    mass = Chem.GetPeriodicTable().GetAtomicWeight(6)
    assert result["density_g_cm3"] == pytest.approx(mass / AVOGADRO / 1e-21)
    assert result["local_probe_exposed_surface_m2_g"] == pytest.approx(expected_area * 1e-20 * AVOGADRO / mass)
    assert result["geometric_void_volume_cm3_g"] == pytest.approx(result["void_fraction"] / result["density_g_cm3"])
    assert result["bet_surface_area_m2_g"] is None
    assert result["connectivity"]["components"][0]["percolating_axes"] == ["a", "b", "c"]
    again = pore_analysis(structure, probe_radius=1, samples=24000, surface_samples=96, grid_size=12, seed=9)
    assert result == again
    json.dumps(result, allow_nan=False)


def test_periodic_self_images_occlude_surface_in_small_cells():
    result = pore_analysis(carbon_cell(3), probe_radius=2, samples=1000, grid_size=8)
    assert result["probe_center_free_fraction"] == 0
    assert result["local_probe_exposed_surface_angstrom2"] == 0
    assert result["connectivity"]["component_count"] == 0


def test_percolation_distinguishes_closed_pockets_from_periodic_channel():
    pocket = np.zeros((8, 8, 8), dtype=bool); pocket[2:5, 2:5, 2:5] = True
    result = _percolation(pocket)
    assert result["component_count"] == 1
    assert result["percolating_probe_center_fraction"] == 0
    channel = np.zeros((8, 8, 8), dtype=bool); channel[:, 3, 3] = True
    result = _percolation(channel)
    assert result["components"][0]["percolating_axes"] == ["a"]
    assert result["percolating_probe_center_fraction"] == pytest.approx(8 / 512)


def test_bet_recovers_known_isotherm_and_area_conversion():
    x = np.array([.01, .03, .05, .08, .1, .15, .2, .25, .3])
    capacity, constant = 3.0, 80.0
    uptake = capacity * constant * x / ((1 - x) * (1 + (constant - 1) * x))
    fit = fit_bet(x, uptake, fit_min=.01)
    assert fit["monolayer_capacity_mmol_g"] == pytest.approx(capacity)
    assert fit["C"] == pytest.approx(constant)
    assert fit["surface_area_m2_g"] == pytest.approx(292.676040936)
    assert fit["r_squared"] == pytest.approx(1)
    assert fit["all_basic_checks_pass"]
    invalid_range = fit_bet(x, uptake, fit_min=.2)
    assert not invalid_range["all_basic_checks_pass"]
    assert not invalid_range["consistency_checks"]["monolayer_relative_pressure_in_fit_range"]


@pytest.mark.parametrize("x,n", [([.1, .1, .2], [1, 2, 3]), ([0, .1, .2], [1, 2, 3]), ([.1, .2, .3], [0, 2, 3]), ([.1, .2], [1, 2])])
def test_bet_rejects_invalid_isotherm(x, n):
    with pytest.raises(ValueError):
        fit_bet(x, n)


def test_finite_patch_layers_preserve_boundary_sites_and_never_claim_periodic():
    from cof_vision.cofdraw_chemistry import prepare_template
    from cof_vision.framework_drawing import VERSION, build_patch
    prepared = prepare_template("O=Cc1cc(C=O)cc(C=O)c1", "Nc1ccc(N)cc1", "imine", "TFB", "Pa", 6, "hcb_3_2")
    recipe = {"schema_version": VERSION, "recipe_id": "science-test", **prepared}
    mol, _ = build_patch(recipe)
    value = patch_3d(recipe, layers=2, stacking="AB", slip=[1.2, .4])
    assert not value["periodic"] and not value["viewer"]["periodic"]
    assert len(value["viewer"]["atoms"]) == 2 * mol.GetNumAtoms()
    assert sum(atom["continuation_site"] for atom in value["viewer"]["atoms"]) == 2 * sum(atom.GetAtomicNum() == 0 for atom in mol.GetAtoms())
    first, second = value["viewer"]["atoms"][0], value["viewer"]["atoms"][mol.GetNumAtoms()]
    assert [second[k] - first[k] for k in ("x", "y", "z")] == pytest.approx([1.2, .4, 3.4])
    assert set(value["exports"]) == {"xyz"}
    with pytest.raises(ValueError, match="explicit"):
        patch_3d(recipe, stacking="AB")
