"""Recovered source graphs retain their actual coordinates, chemistry and scope."""
from copy import deepcopy
import hashlib
import json
from pathlib import Path

import numpy as np
import pytest
from rdkit import Chem

from cof_vision.cofdraw_explicit import build_explicit_patch, clear_source_registration, validate_explicit_recipe
from cof_vision.framework_drawing import build_patch


DATA = Path(__file__).resolve().parents[1] / "data/framework-drawing"


def recipe(name):
    return json.loads((DATA / "recovered" / f"{name}.recipe.json").read_text(encoding="utf-8"))


@pytest.mark.parametrize("path", sorted((DATA / "recovered").glob("*.recipe.json")), ids=lambda p: p.stem)
def test_recovered_graph_preserves_coordinates_and_explicit_tables(path):
    source = json.loads(path.read_text(encoding="utf-8"))
    original = deepcopy(source)
    mol, audit = build_explicit_patch(source)
    graph = source["explicit_graph"]
    assert source == original
    assert np.array_equal(mol.GetConformer().GetPositions()[:, :2], [a["xy"] for a in graph["atoms"]])
    assert mol.GetNumAtoms() == len(graph["atoms"])
    assert mol.GetNumBonds() == len(graph["bonds"])
    for row, atom in zip(graph["atoms"], mol.GetAtoms()):
        assert atom.GetSymbol() == row["element"]
        assert atom.GetIsAromatic() == row["aromatic"]
        assert atom.GetFormalCharge() == row["formal_charge"]
        if "explicit_hs" in row:
            assert atom.GetTotalNumHs() == row["explicit_hs"]
    assert audit["defined_double_bond_geometry"]["valid"]
    assert audit["whole_cof_smiles"] is None
    assert audit["ground_truth"] is False
    assert audit["layout_refinement"]["applied"] is False


def test_fluorescent_vinylene_source_retains_six_E_bonds_and_six_cuts():
    source = recipe("f-cof-1-source-v1")
    mol, audit = build_explicit_patch(source)
    assert mol.GetNumAtoms() == 90
    assert audit["defined_double_bond_geometry"]["checked"] == 6
    assert len(audit["continuation_ports"]) == 6
    assert sum(b.GetStereo() == Chem.BondStereo.STEREOE for b in mol.GetBonds()) == 6
    assert source["explicit_graph"]["coordinate_provenance"]["source_registered"]
    assert "halo" in source["source_style"]


def test_projected_cof300_does_not_merge_crossings_or_claim_crystal_geometry():
    source = recipe("cof-300-projection-v1")
    mol, audit = build_explicit_patch(source)
    assert mol.GetNumAtoms() == 402
    assert len(audit["continuation_ports"]) == 16
    assert source["topology"] is None
    assert source["source_representation_kind"] == "source_informed_projected_network"
    assert not source["explicit_graph"]["coordinate_provenance"]["source_registered"]
    assert source["source_style"]["faded_atoms"]


@pytest.mark.parametrize("name", ["cof-bp-ni1-network-v1", "cof-bp-ni2-network-v1"])
def test_coordination_graphs_retain_dative_bonds_and_charge_balance(name):
    mol, audit = build_explicit_patch(recipe(name))
    assert sum(a.GetAtomicNum() == 28 for a in mol.GetAtoms()) == 4
    assert sum(a.GetFormalCharge() for a in mol.GetAtoms()) == 0
    assert sum(b.GetBondType() == Chem.BondType.DATIVE for b in mol.GetBonds()) == 16
    assert audit["cofex_id"] is None


def test_registry_catalog_and_bundle_bind_the_exact_recipe_bytes():
    registry = json.loads((DATA / "registry.json").read_text(encoding="utf-8"))
    catalog = json.loads((DATA / "catalog/current.json").read_text(encoding="utf-8"))
    for entry in registry["recipes"]:
        digest = hashlib.sha256((DATA / entry["path"]).read_bytes()).hexdigest()
        assert digest == entry["sha256"] == catalog["source_metadata"][entry["recipe_id"]]["source_recipe_sha256"]
    manifest = json.loads((DATA / "current-bundle-manifest.json").read_text(encoding="utf-8"))
    for entry in manifest["files"]:
        assert hashlib.sha256((DATA / entry["path"]).read_bytes()).hexdigest() == entry["sha256"]


def test_explicit_invalid_bond_is_rejected_before_rdkit_mutation():
    source = recipe("f-cof-1-source-v1")
    source["explicit_graph"]["bonds"][0]["a"] = 10000
    with pytest.raises(ValueError, match="endpoints"):
        validate_explicit_recipe(source)


@pytest.mark.parametrize("name,atoms,ports", [
    ("pyurea-source-layout-v1", 232, 8),
    ("ucof-source-layout-v1", 144, 6),
])
def test_source_informed_revisions_keep_finite_scope_and_qualified_layout(name, atoms, ports):
    source = recipe(name)
    mol, audit = build_patch(source)
    assert mol.GetNumAtoms() == atoms
    assert len(audit["continuation_ports"]) == ports
    assert source["cofex_id"] is None
    assert source["topology"] is None
    assert source["layout_evidence"]["exact_atom_registration"] is False
    assert source["layout_evidence"]["source_image_sha256"] == source["source"]["crop_sha256"]
    assert source["explicit_graph"]["coordinate_provenance"]["upstream_generic_coordinates_preserved"] is False
    assert source["layout_evidence"]["maximum_within_aromatic_block_distance_change"] < 1e-10
    assert audit["layout_diagnostics"]["nonbonded_overlap_pairs"] == []
    assert audit["layout_diagnostics"]["bond_length_outliers_35percent"] == 0
    assert len(source["components"]) == 2
    assert "private_alignment_input" not in source


def test_edited_registered_graph_audit_uses_manual_coordinates_without_mutating_caller():
    original = recipe("f-cof-1-source-v1")
    _, source_audit = build_explicit_patch(original)
    assert source_audit["source_registered"] is True
    edited = deepcopy(original)
    edited["designed_variant"] = True
    # A boundary atom can be dragged without changing the defined vinylene stereo.
    boundary = next(a for a in edited["explicit_graph"]["atoms"] if a["element"] == "*")
    boundary["xy"][0] += .15
    before = deepcopy(edited)
    mol, audit = build_explicit_patch(edited)
    assert edited == before
    assert np.array_equal(mol.GetConformer().GetPositions()[:, :2], [a["xy"] for a in edited["explicit_graph"]["atoms"]])
    assert audit["source_registered"] is False
    assert audit["coordinate_provenance"]["source_registered"] is False
    assert audit["coordinate_provenance"]["kind"] == "user_manual_coordinates"
    assert audit["layout_refinement"]["reason"] == "user_manual_coordinates_preserved"
    assert audit["source"] is None
    assert audit["cofex_id"] is None
    assert audit["source_registration_history"]["scope"] == "historical_parent_only"
    assert audit["source_registration_history"]["explicit_graph"]["coordinate_provenance"]["source_registered"] is True


@pytest.mark.parametrize("name", ["f-cof-1-source-v1", "pyurea-source-layout-v1"])
def test_clear_source_registration_is_idempotent_and_keeps_chemistry_and_style(name):
    source = recipe(name)
    original = deepcopy(source)
    assert clear_source_registration(source) is source
    assert source["source_registered"] is False
    assert source["explicit_graph"]["source_coordinates_preserved"] is False
    assert source["explicit_graph"]["coordinate_provenance"]["source_registered"] is False
    assert source["source_style"]["source_registered"] is False
    assert source["source_style"]["exact_style_equivalence"] is False
    assert source["explicit_graph"]["atoms"] == original["explicit_graph"]["atoms"]
    assert source["explicit_graph"]["bonds"] == original["explicit_graph"]["bonds"]
    assert source["source_style"]["atom_colors"] == original["source_style"]["atom_colors"]
    if "layout_evidence" in source:
        assert source["layout_evidence"]["method"] == "user_manual_coordinates"
        assert source["source_registration_history"]["layout_evidence"] == original["layout_evidence"]
    once = deepcopy(source)
    clear_source_registration(source)
    assert source == once


def test_clear_source_registration_handles_assembled_source_anchors():
    source = json.loads((DATA / "tfb-tab-pore-v1.recipe.json").read_text(encoding="utf-8"))
    original_coordinates = deepcopy(source["atom_coordinates"])
    clear_source_registration(source)
    assert source["source_coordinate_provenance"]["source_registered"] is False
    assert source["source_registration_history"]["source_coordinate_provenance"]["source_registered"] is True
    assert source["atom_coordinates"] == original_coordinates


@pytest.mark.parametrize("filename,observed,hidden", [
    ("tfb-tab-pore-v1.recipe.json", 108, 0),
    ("wiley-r198-tppa-so3h-drawing-0002.json", 138, 18),
])
def test_original_source_anchor_layouts_are_used_without_refitting(filename, observed, hidden):
    source = json.loads((DATA / filename).read_text(encoding="utf-8"))
    mol, audit = build_patch(source)
    keyed = {(r["instance"], r["fragment_atom_index"]): r["xy"] for r in source["atom_coordinates"]}
    expected = [keyed[(a["instance"], a["fragment_atom_index"])] for a in audit["atom_origin"]]
    assert np.array_equal(mol.GetConformer().GetPositions()[:, :2], expected)
    assert source["source_coordinate_provenance"]["visible_atom_count"] == observed
    assert source["source_coordinate_provenance"]["unobserved_abbreviated_atom_count"] == hidden
    assert audit["layout_refinement"]["applied"] is False
