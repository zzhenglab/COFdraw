"""Analytic diffraction, invariant composition, and actual bundled-model checks."""
from copy import deepcopy
import hashlib
import json
import math
from pathlib import Path
import re

import numpy as np
import pytest

from cof_vision.cofdraw_characterization import atomic_form_factor, powder_xrd, structure_summary
from cof_vision.cofdraw_science import AVOGADRO, PT, cell_from_parameters, import_cif, structure_exports
from cof_vision.cofdraw_xray_factors import DATA_SHA256, IT92


def crystal(positions=((0, 0, 0),), elements=None, cell=None):
    elements = elements or ["C"] * len(positions)
    return {"cell": cell or [[4, 0, 0], [0, 4, 0], [0, 0, 4]],
            "atoms": [{"element": e, "fractional": list(p)} for e, p in zip(elements, positions)]}


def test_form_factor_snapshot_matches_independent_gemmi_reference():
    # Gemmi 0.7.5 calculate_sf(float stol2), captured from the dev environment.
    # Runtime/CI do not import Gemmi. These exercise indexing, units, and formula.
    expected = {
        "H": [.9999529719, .8542699218, .1397940814, .01449676231],
        "C": [5.999199867, 5.107644081, 1.801606536, 1.223814249],
        "O": [7.999400139, 7.245515347, 2.646724701, 1.468104720],
        "Fe": [25.99040031, 23.68308067, 12.66082001, 7.052098274],
        "Cu": [28.98590088, 27.09569740, 15.10260010, 7.833151340],
    }
    for element, values in expected.items():
        assert atomic_form_factor(element, [0, .01, .2, .8]) == pytest.approx(values, rel=3e-7)
    canonical = json.dumps(IT92, sort_keys=True, separators=(",", ":")).encode()
    assert hashlib.sha256(canonical).hexdigest() == DATA_SHA256


def test_summary_composition_density_and_atomic_site_inventory():
    s = crystal([(0, 0, 0), (.2, .2, .2), (.7, .7, .7)], ["O", "H", "H"],
                cell=cell_from_parameters(6, 7, 8, 80, 75, 65))
    result = structure_summary(s)
    assert result["formula"] == result["reduced_formula"] == "H2O"
    mass = PT.GetAtomicWeight("O") + 2 * PT.GetAtomicWeight("H")
    volume = np.linalg.det(s["cell"])
    assert result["cell_molar_mass_g_mol"] == pytest.approx(mass)
    assert result["cell_volume_angstrom3"] == pytest.approx(volume)
    assert result["density_g_cm3"] == pytest.approx(mass / AVOGADRO / (volume * 1e-24))
    assert list(result["cell_parameters"].values()) == pytest.approx([6, 7, 8, 80, 75, 65])
    assert sum(row["mass_fraction"] for row in result["elements"]) == pytest.approx(1)
    assert sum(row["atomic_fraction"] for row in result["elements"]) == pytest.approx(1)
    for row in result["elements"]:
        assert row["site_density_mmol_g"] == pytest.approx(row["count"] * 1000 / mass)
    assert result["quality"] == "model_geometry"
    assert result["contact_checks"]["available"]
    json.dumps(result, allow_nan=False)


def test_summary_hill_formula_reduction_and_supercell_invariance():
    source = crystal([(0, 0, 0), (.2, .2, .2), (.4, .4, .4), (.6, .6, .6)], ["O", "H", "C", "H"])
    original = structure_summary(source)
    repeated = import_cif(structure_exports(source, repeats=(2, 1, 1))["exports"]["cif"])
    enlarged = structure_summary(repeated)
    assert original["formula"] == "CH2O"
    assert enlarged["formula"] == "C2H4O2"
    assert enlarged["reduced_formula"] == original["reduced_formula"] == "CH2O"
    assert enlarged["density_g_cm3"] == pytest.approx(original["density_g_cm3"])
    for first, second in zip(original["elements"], enlarged["elements"]):
        assert second["site_density_mmol_g"] == pytest.approx(first["site_density_mmol_g"])
        assert second["mass_fraction"] == pytest.approx(first["mass_fraction"])


def test_summary_flags_periodic_clashes_and_does_not_call_candidate_distance_global():
    report = structure_summary(crystal([(0, 0, 0), (.99, 0, 0)]))
    assert report["contact_checks"]["short_contact_count"] == 1
    assert report["contact_checks"]["minimum_candidate_distance_angstrom"] == pytest.approx(.04)
    assert any("short periodic contacts" in message for message in report["warnings"])
    dilute = structure_summary(crystal(cell=[[50, 0, 0], [0, 50, 0], [0, 0, 50]]))
    assert dilute["contact_checks"]["minimum_candidate_distance_angstrom"] is None


def test_simple_cubic_bragg_positions_multiplicity_and_intensity():
    result = powder_xrd(crystal(), two_theta_max=70)
    first = result["peaks"][:4]
    for peak, norm2, multiplicity in zip(first, [1, 2, 3, 4], [6, 12, 8, 6]):
        d = 4 / math.sqrt(norm2)
        angle = math.degrees(2 * math.asin(1.5406 / (2 * d)))
        assert peak["d_spacing_angstrom"] == pytest.approx(d)
        assert peak["two_theta_degrees"] == pytest.approx(angle)
        assert peak["multiplicity"] == multiplicity == len(peak["hkls"])
        assert all(sum(x*x for x in hkl) == norm2 for hkl in peak["hkls"])
        theta = math.radians(angle / 2)
        lp = (1 + math.cos(2*theta)**2) / (math.sin(theta)**2 * math.cos(theta))
        expected = multiplicity * atomic_form_factor("C", 1 / (4*d*d))**2 * lp
        assert peak["integrated_intensity_arbitrary"] == pytest.approx(float(expected))
    assert max(result["profile"]["intensity"]) == pytest.approx(100)
    assert max(p["intensity"] for p in result["peaks"]) == pytest.approx(100)


@pytest.mark.parametrize("positions,first_norm2,allowed", [
    ([(0, 0, 0), (.5, .5, .5)], 2, lambda h: sum(h) % 2 == 0),
    ([(0, 0, 0), (0, .5, .5), (.5, 0, .5), (.5, .5, 0)], 3,
     lambda h: len({x % 2 for x in h}) == 1),
])
def test_body_and_face_centered_extinction_rules(positions, first_norm2, allowed):
    result = powder_xrd(crystal(positions), two_theta_max=100)
    assert all(allowed(hkl) for peak in result["peaks"] for hkl in peak["hkls"])
    assert result["peaks"][0]["d_spacing_angstrom"] == pytest.approx(4 / math.sqrt(first_norm2))


def test_two_element_structure_factor_interference_is_not_count_squared():
    # CsCl-like two-site basis: sum for even h+k+l, difference for odd.
    source = crystal([(0, 0, 0), (.5, .5, .5)], ["Cs", "Cl"])
    peaks = powder_xrd(source, two_theta_max=40)["peaks"]
    for peak in peaks:
        d = peak["d_spacing_angstrom"]
        f_cs, f_cl = [atomic_form_factor(e, 1/(4*d*d)) for e in ["Cs", "Cl"]]
        sign = 1 if sum(peak["hkls"][0]) % 2 == 0 else -1
        theta = math.radians(peak["two_theta_degrees"] / 2)
        lp = (1 + math.cos(2*theta)**2)/(math.sin(theta)**2 * math.cos(theta))
        expected = peak["multiplicity"] * float(f_cs + sign*f_cl)**2 * lp
        assert peak["integrated_intensity_arbitrary"] == pytest.approx(expected)


def test_triclinic_reciprocal_metric_rotation_and_origin_invariance():
    source = crystal([(0.12, .25, .33), (.61, .39, .77)], ["C", "O"], cell=cell_from_parameters(5, 6, 7, 75, 82, 68))
    original = powder_xrd(source, two_theta_max=65)
    metric_inverse = np.linalg.inv(np.asarray(source["cell"]) @ np.asarray(source["cell"]).T)
    for peak in original["peaks"]:
        for h in peak["hkls"]:
            assert peak["d_spacing_angstrom"] == pytest.approx(1 / math.sqrt(np.asarray(h) @ metric_inverse @ h))
    changed = deepcopy(source)
    rotation = np.array([[0, -1, 0], [1, 0, 0], [0, 0, 1]])
    changed["cell"] = (np.asarray(source["cell"]) @ rotation).tolist()
    for atom in changed["atoms"]:
        atom["fractional"] = (np.asarray(atom["fractional"]) + [.19, -.27, .41]).tolist()
    other = powder_xrd(changed, two_theta_max=65)
    assert other["profile"]["intensity"] == pytest.approx(original["profile"]["intensity"], abs=2e-10)
    assert [p["intensity"] for p in other["peaks"]] == pytest.approx([p["intensity"] for p in original["peaks"]])


def test_xrd_supercell_extinctions_preserve_normalized_pattern():
    source = crystal([(0, 0, 0), (.21, .34, .45)], ["C", "N"])
    repeated = import_cif(structure_exports(source, repeats=(2, 2, 2))["exports"]["cif"])
    original, enlarged = [powder_xrd(s, two_theta_max=70) for s in (source, repeated)]
    assert enlarged["profile"]["intensity"] == pytest.approx(original["profile"]["intensity"], abs=2e-10)
    assert [p["two_theta_degrees"] for p in enlarged["peaks"]] == pytest.approx([p["two_theta_degrees"] for p in original["peaks"]])
    for p, q in zip(original["peaks"], enlarged["peaks"]):
        assert q["integrated_intensity_arbitrary"] == pytest.approx(64 * p["integrated_intensity_arbitrary"])


def test_isotropic_b_factor_changes_intensity_by_correct_amplitude_squared():
    zero, damped = [powder_xrd(crystal(), two_theta_max=70, b_iso=b) for b in (0, 3)]
    for original, changed in zip(zero["peaks"], damped["peaks"]):
        stol2 = 1 / (4 * original["d_spacing_angstrom"]**2)
        assert changed["integrated_intensity_arbitrary"] / original["integrated_intensity_arbitrary"] == pytest.approx(math.exp(-2 * 3 * stol2))


def test_gaussian_fwhm_grid_csv_and_reproducibility():
    result = powder_xrd(crystal(), two_theta_min=21, two_theta_max=24, fwhm=.2, step=.005)
    angles, signal = [np.asarray(result["profile"][key]) for key in ("two_theta_degrees", "intensity")]
    half = angles[signal >= 50]
    assert half[-1] - half[0] == pytest.approx(.2, abs=.01)
    rows = result["exports"]["csv"].splitlines()
    assert rows[0] == "two_theta_degrees,relative_intensity"
    assert len(rows) == len(signal) + 1
    assert np.loadtxt(rows[1:], delimiter=",")[:, 1] == pytest.approx(signal, abs=1e-6)
    assert result == powder_xrd(crystal(), two_theta_min=21, two_theta_max=24, fwhm=.2, step=.005)
    json.dumps(result, allow_nan=False)


def test_empty_angular_window_and_boundary_tails_are_explicit():
    empty = powder_xrd(crystal(), two_theta_min=2, two_theta_max=10)
    assert empty["peaks"] == []
    assert not any(empty["profile"]["intensity"])
    tail = powder_xrd(crystal(), two_theta_min=22.22, two_theta_max=22.4)
    assert tail["peaks"] == []
    assert max(tail["profile"]["intensity"]) == pytest.approx(100)
    assert any("only of neighboring peak tails" in text for text in tail["warnings"])


@pytest.mark.parametrize("kwargs", [
    {"wavelength": 0}, {"wavelength": float("nan")}, {"two_theta_min": 50, "two_theta_max": 20},
    {"step": .2, "fwhm": .15}, {"b_iso": -1}, {"fwhm": True},
    {"wavelength": .1, "two_theta_max": 100}, {"step": .005, "two_theta_max": 170},
])
def test_invalid_or_excessive_requests_fail_explicitly(kwargs):
    with pytest.raises(ValueError):
        powder_xrd(crystal(), **kwargs)


def test_oversized_reciprocal_search_and_missing_factors_rejected():
    with pytest.raises(ValueError, match="bounded workload"):
        powder_xrd(crystal(cell=[[500, 0, 0], [0, 500, 0], [0, 0, 499]]))
    with pytest.raises(ValueError, match="IT92"):
        powder_xrd(crystal(elements=["Es"]))
    with pytest.raises(ValueError, match="periodic"):
        structure_summary({**crystal(), "periodic": False})


MODEL_ROOT = Path(__file__).resolve().parents[1] / "data/crystal-models"
MODELS = json.loads((MODEL_ROOT / "manifest.json").read_text("utf-8"))["templates"]


@pytest.mark.parametrize("record", MODELS, ids=lambda row: row["id"])
def test_actual_cof_templates_have_finite_reproducible_patterns_and_consistent_mass(record):
    source = import_cif((MODEL_ROOT / record["path"]).read_text("utf-8"))
    summary = structure_summary(source)
    composition = lambda text: {e: int(n or 1) for e, n in re.findall(r"([A-Z][a-z]?)(\d*)", text)}
    assert composition(summary["formula"]) == composition(record["formula"])
    assert summary["atom_count"] == record["atoms"]
    assert summary["contact_checks"]["short_contact_count"] == 0
    assert summary["density_g_cm3"] > 0
    result = powder_xrd(source)
    assert result["peak_count"] > 10
    assert np.isfinite(result["profile"]["intensity"]).all()
    assert all(0 < p["intensity"] <= 100 for p in result["peaks"])
    assert result["structure_sha256"] == summary["structure_sha256"]
    assert result == powder_xrd(source)
