"""Reviewed layout revisions preserve chemistry, source qualifiers and history."""
from copy import deepcopy
import importlib.util
import json
from pathlib import Path

import numpy as np
import pytest

from cof_vision.framework_drawing import build_patch

REPO = Path(__file__).resolve().parents[1]
DATA = REPO / "data/framework-drawing"
spec = importlib.util.spec_from_file_location("source_layout_importer", REPO / "scripts/import_source_layouts.py")
importer = importlib.util.module_from_spec(spec)
spec.loader.exec_module(importer)
REPORT = importer.read(DATA / "catalog/source-layout-revisions.audit.json")
CATALOG = importer.read(DATA / "catalog/current.json")
REGISTRY = importer.read(DATA / "registry.json")


def pair(row):
    entry = next(r for r in REGISTRY["recipes"] if r["recipe_id"] == row["recipe_id"])
    old_path = DATA / "history/recipes" / (row["previous_recipe_sha256"] + ".recipe.json")
    new_path = DATA / entry["path"]
    assert importer.digest(old_path) == row["previous_recipe_sha256"]
    assert importer.digest(new_path) == entry["sha256"] == row["recipe_sha256"]
    return importer.read(old_path), importer.read(new_path)


@pytest.mark.parametrize("row", REPORT["revisions"], ids=lambda r: r["recipe_id"])
def test_revised_coordinates_are_exactly_retained_with_the_same_chemical_graph(row):
    old, new = pair(row)
    before, old_audit = build_patch(old)
    after, new_audit = build_patch(new)
    assert importer.graph_signature(before, old_audit) == importer.graph_signature(after, new_audit) == row["graph_signature_sha256"]
    for key in ("components", "connections", "instances", "source", "cofex_source_records", "cofex_id", "topology"):
        assert old.get(key) == new.get(key)
    assert new["source_registered"] is False
    assert new["coordinate_provenance"]["source_registered"] is False
    keyed = {(c["instance"], c["fragment_atom_index"]): c["xy"] for c in new["atom_coordinates"]}
    expected = np.array([keyed[(o["instance"], o["fragment_atom_index"])] for o in new_audit["atom_origin"]])
    assert np.array_equal(after.GetConformer().GetPositions()[:, :2], expected)
    assert new_audit["layout_refinement"]["applied"] is False
    assert len(keyed) == after.GetNumAtoms() == row["source_overlay_evidence"]["coordinate_count"]
    for name in {o["instance"] for o in old_audit["atom_origin"]}:
        indices = [i for i, o in enumerate(old_audit["atom_origin"]) if o["instance"] == name and before.GetAtomWithIdx(i).GetAtomicNum()]
        a, b = before.GetConformer().GetPositions()[indices], after.GetConformer().GetPositions()[indices]
        assert np.allclose(np.linalg.norm(a[:, None] - a[None, :], axis=-1),
                           np.linalg.norm(b[:, None] - b[None, :], axis=-1), atol=1e-12, rtol=0)
    diagnostics = new_audit["layout_diagnostics"]
    assert not diagnostics["nonbonded_overlap_pairs"] and not diagnostics["nonincident_bond_crossings"]
    assert diagnostics["bond_length_outliers_35percent"] == 0
    assert max(a["deviation_degrees"] for a in diagnostics["junction_angles"]) < .01
    assert max(abs(x - 1) for x in diagnostics["junction_bond_ratios"]) < .0001


def test_layout_admission_preserves_all_other_sources_and_the_monomer_step():
    old_catalog_path = DATA / "catalog/history" / (REPORT["previous_catalog_sha256"] + ".json")
    old_registry_path = DATA / "history/registries" / (REPORT["previous_registry_sha256"] + ".json")
    assert importer.digest(old_catalog_path) == REPORT["previous_catalog_sha256"]
    assert importer.digest(old_registry_path) == REPORT["previous_registry_sha256"]
    old_catalog, old_registry = importer.read(old_catalog_path), importer.read(old_registry_path)
    # The later master-record admission adds entries after this layout revision.
    # Check the layout event against its immutable output snapshot, then require
    # every existing row to remain byte-for-byte equivalent in the current list.
    layout_catalog_path = DATA / "catalog/history" / (REPORT["catalog_sha256"] + ".json")
    layout_catalog = importer.read(layout_catalog_path)
    assert old_catalog["monomer_presets"] == layout_catalog["monomer_presets"]
    current_rows = {row['id']: row for row in CATALOG['monomer_presets']}
    assert all(current_rows[row['id']] == row for row in old_catalog['monomer_presets'])
    revised = {r["recipe_id"] for r in REPORT["revisions"]}
    assert revised == {"tppa-pore-v1", "tpb-dmtp-cof-pore-v1"}
    for row in old_registry["recipes"]:
        if row["recipe_id"] not in revised:
            assert row == next(r for r in REGISTRY["recipes"] if r["recipe_id"] == row["recipe_id"])
            assert importer.digest(DATA / row["path"]) == row["sha256"]
            assert CATALOG["source_metadata"][row["recipe_id"]] == old_catalog["source_metadata"][row["recipe_id"]]
    assert REPORT["catalog_sha256"] == importer.digest(layout_catalog_path)
    assert REPORT["registry_sha256"] == importer.digest(DATA / "registry.json")
    audit = importer.read(DATA / "catalog/current.audit.json")
    assert audit["source_layout_revisions_sha256"] == importer.digest(DATA / "catalog/source-layout-revisions.audit.json")
    assert CATALOG["source_registry_sha256"] == REPORT["registry_sha256"]


@pytest.mark.parametrize("change", ["chemistry", "registration", "coordinates", "designed"])
def test_coordinate_admission_rejects_changed_chemistry_or_unsupported_claims(change):
    row = REPORT["revisions"][0]
    original, candidate = pair(row)
    candidate = deepcopy(candidate)
    if change == "chemistry":
        candidate["linkage"] = "imine"
    elif change == "registration":
        candidate["coordinate_provenance"]["source_registered"] = True
    elif change == "coordinates":
        candidate["atom_coordinates"] = candidate["atom_coordinates"][:-1]
    else:
        candidate["designed_variant"] = True
    with pytest.raises(ValueError):
        importer.revised_recipe(original, candidate, row["source_overlay_evidence"])
