"""Drawing-projection geometry and scientific endpoint provenance regressions."""
import hashlib
import json
import math
from pathlib import Path
import shutil

import numpy as np
import pytest
from fastapi.testclient import TestClient
from rdkit import Chem

from cof_vision.cofdraw_desktop import create_desktop_app
from cof_vision.cofdraw_measure import measure_drawing
from cof_vision.framework_drawing import VERSION


def explicit_ring(points, *, hydrogen=None, closed=True):
    atoms = [{"index": i, "element": "C", "xy": list(map(float, xy))} for i, xy in enumerate(points)]
    bonds = [{"a": i, "b": i + 1, "order": 1} for i in range(len(points) - 1)]
    if closed:
        bonds.append({"a": len(points) - 1, "b": 0, "order": 1})
    if hydrogen is not None:
        bonds.append({"a": 0, "b": len(atoms), "order": 1})
        atoms.append({"index": len(atoms), "element": "H", "xy": list(hydrogen)})
    return {"schema_version": VERSION, "recipe_id": "measure-fixture", "explicit_graph": {
        "coordinate_system": "cartesian_y_up", "atoms": atoms, "bonds": bonds}}


def regular_points(count=20, radius=5):
    return np.array([[radius * math.cos(2 * math.pi * i / count), radius * math.sin(2 * math.pi * i / count)]
                     for i in range(count)])


def test_regular_ring_calibration_matches_analytic_circle_and_is_scale_invariant():
    points = regular_points()
    result = measure_drawing(explicit_ring(points))
    calibrated_radius = 1.5 / (2 * math.sin(math.pi / len(points)))
    expected_diameter = 2 * (calibrated_radius - Chem.GetPeriodicTable().GetRvdw(6))
    assert result["diameter_angstrom"] <= expected_diameter + 1e-8
    assert expected_diameter - result["diameter_angstrom"] <= 2 * result["grid_diagonal_angstrom"]
    transformed = measure_drawing(explicit_ring(points * 7.5 + [71.3, -48.2]))
    assert transformed["diameter_angstrom"] == pytest.approx(result["diameter_angstrom"], abs=1e-9)
    assert transformed["projected_ring_area_angstrom2"] == pytest.approx(result["projected_ring_area_angstrom2"])
    assert result["sampling_lower_bound"]
    json.dumps(result, allow_nan=False)


def test_inscribed_circle_cannot_escape_between_long_projected_bonds():
    # Deliberately elongated bond gaps: atom exclusion alone permits a circle
    # wider than the thin ring. Polygon-edge containment must also constrain it.
    points = [(0, 0), (1, 0), (2, 0), (18, 0), (19, 0), (20, 0),
              (20, 2), (19, 2), (18, 2), (2, 2), (1, 2), (0, 2)]
    result = measure_drawing(explicit_ring(points))
    assert 2.9 < result["diameter_angstrom"] <= 3.0
    assert result["projected_ring_area_angstrom2"] == pytest.approx(90)


def test_annotation_endpoints_bind_exact_reported_diameter_to_drawing_coordinates():
    points = regular_points() * 7.5 + [71.3, -48.2]
    recipe = explicit_ring(points)
    before = json.dumps(recipe, sort_keys=True)
    result = measure_drawing(recipe)
    annotation = result['annotation']
    endpoints = np.asarray(annotation['endpoints'])
    scale = annotation['angstrom_per_drawing_unit']
    assert annotation['coordinate_system'] == 'drawing_xy_cartesian_y_up'
    assert endpoints.mean(axis=0) == pytest.approx(annotation['center'])
    assert np.asarray(annotation['center']) * scale == pytest.approx(result['center_angstrom'])
    assert np.linalg.norm(endpoints[1] - endpoints[0]) * scale == pytest.approx(result['diameter_angstrom'])
    assert result['diameter_nm'] * 10 == pytest.approx(result['diameter_angstrom'])
    assert endpoints[0, 1] == endpoints[1, 1]
    assert 'inscribed circle' in annotation['meaning']
    assert annotation['label_qualification'] == 'Pore opening'
    for anchor in annotation['coordinate_anchors']:
        assert anchor['xy'] == pytest.approx(points[anchor['atom_index']])
    # Every point on the diameter stays clear of explicit-atom exclusion disks.
    samples = np.linspace(endpoints[0], endpoints[1], 101) * scale
    clearance = np.linalg.norm(samples[:, None] - points[None] * scale, axis=2) - Chem.GetPeriodicTable().GetRvdw(6)
    assert clearance.min() >= -1e-8
    assert json.dumps(recipe, sort_keys=True) == before


def test_measure_api_retains_annotation_units_and_scope(tmp_path):
    with TestClient(create_desktop_app(tmp_path)) as client:
        response = client.post('/api/drawing/science/measure', json={'recipe': explicit_ring(regular_points())})
    assert response.status_code == 200, response.text
    result = response.json()
    assert len(result['annotation']['endpoints']) == 2
    assert result['diameter_nm'] == pytest.approx(result['diameter_angstrom'] / 10)
    assert result['sampling_lower_bound'] is True
    assert 'lower bound for the projected opening' in result['scope']


def test_explicit_hydrogen_is_an_obstacle_in_the_supplied_projection():
    points = regular_points()
    empty = measure_drawing(explicit_ring(points))
    blocked = measure_drawing(explicit_ring(points, hydrogen=[0, 0]))
    assert blocked["explicit_hydrogens_included"] == 1
    assert blocked["diameter_angstrom"] < empty["diameter_angstrom"] / 2


def test_crossed_macrocycle_and_open_chain_do_not_produce_pore_estimates():
    points = regular_points(12)
    crossed = points[[(i * 5) % 12 for i in range(12)]]
    with pytest.raises(ValueError, match="crossed"):
        measure_drawing(explicit_ring(crossed))
    with pytest.raises(ValueError, match="closed macrocyclic"):
        measure_drawing(explicit_ring(points, closed=False))


def test_edited_recipe_cannot_inherit_reference_identity_in_scientific_patch(tmp_path):
    recipe = explicit_ring(regular_points())
    recipe.update(cofex_id="UNVERIFIED-INHERITED-ID", source={"doi": "10.1000/example"},
                  reference={"name": "Parent material"}, topology="hcb")
    with TestClient(create_desktop_app(tmp_path)) as client:
        response = client.post("/api/drawing/science/patch", json={"recipe": recipe, "layers": 1})
    assert response.status_code == 200, response.text
    audit = response.json()["source_audit"]
    assert audit.get("cofex_id") is None
    assert audit.get("source") is None
    assert audit.get("reference") is None
    assert audit.get("topology") is None


def test_periodic_template_api_preserves_input_cif_hash_with_attribution(tmp_path):
    folder = Path(__file__).resolve().parents[1] / "data/crystal-models"
    shutil.copytree(folder, tmp_path / "crystal-models")
    record = json.loads((folder / "manifest.json").read_text("utf-8"))["templates"][0]
    with TestClient(create_desktop_app(tmp_path)) as client:
        response = client.post("/api/drawing/science/model", json={"template_id": record["id"]})
    assert response.status_code == 200, response.text
    provenance = response.json()["structure"]["provenance"]
    assert provenance["input_sha256"] == hashlib.sha256((folder / record["path"]).read_text("utf-8").encode("utf-8")).hexdigest()
    assert provenance["template"] == record["id"]
    assert "generator" in provenance


@pytest.mark.parametrize("atom", [None, 0, "C"])
def test_porosity_api_rejects_nonobject_atom_entries_with_422(tmp_path, atom):
    structure = {"cell": [[10, 0, 0], [0, 10, 0], [0, 0, 10]], "atoms": [atom]}
    with TestClient(create_desktop_app(tmp_path), raise_server_exceptions=False) as client:
        response = client.post("/api/drawing/science/porosity", json={"structure": structure})
    assert response.status_code == 422
    assert "must be an object" in response.json()["detail"]
