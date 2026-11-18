"""Validate admitted finite graphs, source bindings and usable design pairs."""
from collections import Counter
from copy import deepcopy
import hashlib
import importlib.util
import json
from pathlib import Path
import re

import pytest
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

from cof_vision.cofdraw_chemistry import decorated_monomers, prepare_template
from cof_vision.cofdraw_core import read_catalog
from cof_vision.framework_drawing import VERSION, build_patch

REPO = Path(__file__).resolve().parents[1]
DATA = REPO / "data/framework-drawing"
CATALOG = read_catalog(REPO / "data")
ROWS = {r["id"]: r for r in CATALOG["monomer_presets"]}
SELECTION = json.loads((DATA / "catalog/manual-monomer-selection.json").read_text("utf-8"))["selection"]
AUDIT = json.loads((DATA / "catalog/manual-monomer-admission.audit.json").read_text("utf-8"))
spec = importlib.util.spec_from_file_location("manual_monomers", REPO / "scripts/import_manual_monomers.py")
importer = importlib.util.module_from_spec(spec)
spec.loader.exec_module(importer)


@pytest.mark.parametrize("selected", SELECTION, ids=lambda x: x["id"])
def test_selected_graphs_are_complete_hash_bound_and_qualified(selected):
    row = ROWS[selected["id"]]
    mol = importer.finite_molecule(row["smiles"])
    assert Chem.MolToSmiles(mol, isomericSmiles=True) == selected["canonical_smiles"]
    assert rdMolDescriptors.CalcMolFormula(mol) == row["formula"]
    assert row["ground_truth"] is False and row["reviewer_type"] == "agent"
    assert row["source_evidence"][0]["annotation_file_sha256"] == selected["annotation_file_sha256"]
    assert row["source_evidence"][0]["source_crop_sha256"] == selected["source_crop_sha256"]
    assert row["source_evidence"][0]["annotation_id"] == selected["annotation_id"]
    assert row["source_dois"][0] == row["source_evidence"][0]["source_doi"]
    assert row["node_symmetry"] is None and row["qualification"]
    assert not re.search(r'(?:^|["\s])[A-Za-z]:[\\/]|/Users/|/home/|page_text', json.dumps(row))


def test_new_catalog_preserves_existing_graphs_and_all_source_bindings():
    old = json.loads((DATA / "catalog/history" / (AUDIT["previous_catalog_sha256"] + ".json")).read_text("utf-8"))
    current_sha = importer.digest(DATA / "catalog/current.json")
    admitted_path = DATA / "catalog/current.json" if current_sha == AUDIT["catalog_sha256"] else (
        DATA / "catalog/history" / (AUDIT["catalog_sha256"] + ".json"))
    assert importer.digest(admitted_path) == AUDIT["catalog_sha256"]
    admitted = json.loads(admitted_path.read_text("utf-8"))
    assert len(old["monomer_presets"]) == 49 and len(admitted["monomer_presets"]) == 63
    for row in old["monomer_presets"]:
        assert ROWS[row["id"]] == row
    # This audit is the historical first 14-entry admission, not the subsequent
    # master-key update. Preserve its exact rows across later additive updates.
    for row in admitted["monomer_presets"]:
        assert ROWS[row["id"]] == row
    for key in ("recipe_bindings", "source_metadata", "source_registry_sha256"):
        assert admitted[key] == old[key]
    graphs = [Chem.MolToSmiles(Chem.MolFromSmiles(r["smiles"]), isomericSmiles=True) for r in ROWS.values()]
    assert len(set(graphs)) == len(graphs)
    registry = DATA / "registry.json"
    if importer.digest(registry) != AUDIT["registry_sha256"]:
        registry = DATA / "history/registries" / (AUDIT["registry_sha256"] + ".json")
    assert importer.digest(registry) == AUDIT["registry_sha256"]


def test_portable_manifest_and_admission_audits_bind_every_file():
    manifest = json.loads((DATA / "current-bundle-manifest.json").read_text("utf-8"))
    assert manifest["catalog_sha256"] == importer.digest(DATA / "catalog/current.json")
    current = json.loads((DATA / "catalog/current.audit.json").read_text("utf-8"))
    assert current["monomer_count"] == len(ROWS) and current["source_count"] == 17
    assert current["manual_monomer_admission_sha256"] == importer.digest(DATA / "catalog/manual-monomer-admission.audit.json")
    for record in manifest["files"]:
        path = (DATA / record["path"]).resolve()
        assert path.is_relative_to(DATA.resolve())
        assert path.stat().st_size == record["bytes"]
        assert importer.digest(path) == record["sha256"]


@pytest.mark.parametrize("a,b,linkage,template", [
    ("tfppy_ome", "pa", "imine", "sql_4_2"),
    ("tfppy_f", "pa", "imine", "sql_4_2"),
    ("tfppy_13c", "pa", "imine", "sql_4_2"),
    ("tfb", "dabd", "imine", "hcb_3_2"),
    ("tfb", "thz", "acylhydrazone", "hcb_3_2"),
    ("ta", "tape", "imine", "sql_4_2"),
    ("tfb", "odh", "acylhydrazone", "hcb_3_2"),
    ("tfb_oh", "tapb", "imine", "hcb_3_3"),
    ("tfb_oh2", "tapb", "imine", "hcb_3_3"),
    ("tfb", "daaq", "imine", "hcb_3_2"),
])
def test_each_new_supported_monomer_builds_a_valid_explicit_pair(a, b, linkage, template):
    recipe = prepare_template(ROWS[a]["smiles"], ROWS[b]["smiles"], linkage, a, b, 6, template)
    recipe.update(schema_version=VERSION, recipe_id="monomer-admission-design", designed_variant=True)
    mol, audit = build_patch(recipe)
    # The builder validates valence while retaining each finite component's
    # aromatic assignment. Do not assign new aromatic circuits across the
    # assembled macrocycle (both chemically unjustified and combinatorial).
    mol.UpdatePropertyCache(strict=True)
    assert len(Chem.GetMolFrags(mol)) == 1
    assert audit["atoms"] == mol.GetNumAtoms()
    assert min(recipe["compatibility"]["reactive_site_counts"]) >= 2
    assert recipe.get("cofex_id") is None
    assembled = {(o["instance"], o["fragment_atom_index"]): i for i, o in enumerate(audit["atom_origin"])}
    for instance in recipe["instances"]:
        fragment = Chem.MolFromSmiles(recipe["components"][instance["component"]]["fragment_smiles"])
        for bond in fragment.GetBonds():
            left = assembled.get((instance["id"], bond.GetBeginAtomIdx()))
            right = assembled.get((instance["id"], bond.GetEndAtomIdx()))
            if left is not None and right is not None:
                actual = mol.GetBondBetweenAtoms(left, right)
                assert actual is not None and actual.GetBondType() == bond.GetBondType()
                assert actual.GetIsAromatic() == bond.GetIsAromatic()
    if a == "tfppy_13c":
        assert sum(atom.GetIsotope() == 13 for atom in mol.GetAtoms()) == 16


def test_reference_only_monomers_are_inspectable_without_automatic_roles():
    decorated = {row["id"]: row for row in decorated_monomers(CATALOG["monomer_presets"])}
    assert set(AUDIT["reference_only_ids"]) == {"ppcn", "tbpm_boronic", "allyl_tbpa", "tmc"}
    for key in AUDIT["reference_only_ids"]:
        row = decorated[key]
        assert row["catalog_use"] == "reference_only"
        assert row["automatic_prepare_supported"] is False
        assert not row["compatible_roles"]
    assert next(x["count"] for x in decorated["ppcn"]["reactive_groups"] if x["kind"] == "nitrile") == 4
    assert next(x["count"] for x in decorated["tbpm_boronic"]["reactive_groups"] if x["kind"] == "boronic_acid") == 4
    assert next(x["count"] for x in decorated["allyl_tbpa"]["reactive_groups"] if x["kind"] == "boronic_acid") == 3


def test_isotopes_salt_scope_and_substituent_positions_are_retained():
    mol = Chem.MolFromSmiles(ROWS["tfppy_13c"]["smiles"])
    isotopes = [a for a in mol.GetAtoms() if a.GetIsotope()]
    assert len(isotopes) == 4 and all(a.GetAtomicNum() == 6 and a.GetIsotope() == 13 for a in isotopes)
    carbonyls = {hit[1] for hit in mol.GetSubstructMatches(Chem.MolFromSmarts("[c][CH1](=O)"))}
    assert {a.GetIdx() for a in isotopes} == carbonyls
    assert Chem.MolToSmiles(mol) != ROWS["tfppy"]["canonical_smiles"]
    assert "dihydrochloride" in ROWS["dabd"]["qualification"]
    assert Counter(a.GetSymbol() for a in Chem.MolFromSmiles(ROWS["dabd"]["smiles"]).GetAtoms())["Cl"] == 0
    for key, sub in (("tfppy_ome", "O[CH3]"), ("tfppy_f", "F")):
        molecule = Chem.MolFromSmiles(ROWS[key]["smiles"])
        assert len(molecule.GetSubstructMatches(Chem.MolFromSmarts(sub))) == 4


@pytest.mark.parametrize("smiles", ["*c1ccccc1", "[CH3]", "CC.O", "not-a-smiles"])
def test_importer_rejects_boundary_radical_disconnected_and_invalid_graphs(smiles):
    with pytest.raises(ValueError):
        importer.finite_molecule(smiles)


def test_importer_rejects_changed_selected_graph_before_reading_source_files():
    selected = deepcopy(SELECTION[0])
    row = {"task": "molecule", "review_status": "agent_reviewed", "molecules": [{"smiles": "CC"}]}
    with pytest.raises(ValueError, match="changed after visual review"):
        importer.make_row(selected, row, annotation_sha256="0" * 64, annotation_line=1, curation_sha256="0" * 64)
