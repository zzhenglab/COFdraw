"""Source/installed-master binding, stereochemistry and actual drawing use."""
from copy import deepcopy
import importlib.util
import json
from pathlib import Path
import re
import sqlite3
import sys

import pytest
from rdkit import Chem

from cof_vision.cofdraw_chemistry import prepare_template
from cof_vision.cofdraw_core import read_catalog
from cof_vision.framework_drawing import VERSION, build_patch

REPO = Path(__file__).resolve().parents[1]
DATA = REPO / "data/framework-drawing"
CATALOG = read_catalog(REPO / "data")
ROWS = {row["id"]: row for row in CATALOG["monomer_presets"]}
AUDIT = json.loads((DATA / "catalog/master-key-monomer-admission.audit.json").read_text("utf-8"))
SELECTION = json.loads((DATA / "catalog/master-key-monomer-selection.json").read_text("utf-8"))["selection"]
sys.path.insert(0, str(REPO / "scripts"))
spec = importlib.util.spec_from_file_location("master_monomers", REPO / "scripts/import_master_key_monomers.py")
importer = importlib.util.module_from_spec(spec)
spec.loader.exec_module(importer)


@pytest.mark.parametrize("selection", SELECTION, ids=lambda r: r["id"])
def test_selected_complete_graph_and_source_master_bindings(selection):
    row = ROWS[selection["id"]]
    assert Chem.MolToSmiles(importer.finite_molecule(row["smiles"]), isomericSmiles=True) == selection["canonical_smiles"]
    evidence = row["source_evidence"][0]
    assert evidence["source_crop_sha256"] == selection["source_crop_sha256"]
    assert evidence["annotation_file_sha256"] == selection["annotation_file_sha256"]
    assert evidence["source_doi"] == selection["source_doi"]
    assert row["master_binding"]["record_ids"] == sorted(selection["master_record_ids"])
    assert row["master_binding"]["component_ids"] == sorted(selection["master_component_ids"])
    assert row["master_binding"]["reference_index_sha256"] == AUDIT["reference_index_sha256"]
    assert row["master_binding"]["installation_receipt_sha256"] == selection["installation_receipt_sha256"]
    assert row["ground_truth"] is False and row["master_binding"]["ground_truth"] is False
    assert row["reviewer_type"] == "agent" and "not independent human verification" in row["master_binding"]["qualification"]
    assert row["node_symmetry"] is None


def test_additive_history_and_unchanged_source_recipes():
    old = json.loads((DATA / "catalog/history" / (AUDIT["previous_catalog_sha256"] + ".json")).read_text("utf-8"))
    admitted_path=DATA / 'catalog/current.json'
    if importer.digest(admitted_path)!=AUDIT['catalog_sha256']:
        admitted_path=DATA / 'catalog/history' / (AUDIT['catalog_sha256']+'.json')
    admitted=json.loads(admitted_path.read_text('utf-8'))
    assert importer.digest(admitted_path)==AUDIT['catalog_sha256']
    assert len(old["monomer_presets"]) == 63 and len(admitted['monomer_presets']) == 81
    for row in admitted['monomer_presets']:
        assert ROWS[row['id']]==row
    assert len(AUDIT["added_ids"]) == 18
    assert len(AUDIT["automatic_drawing_candidate_ids"]) == 15
    assert AUDIT["reference_only_ids"] == ["aph", "flt_boronic", "ash"]
    for row in old["monomer_presets"]:
        assert ROWS[row["id"]] == row
    for key in ("recipe_bindings", "source_metadata", "source_registry_sha256"):
        assert CATALOG[key] == old[key]
    assert importer.digest(DATA / "registry.json") == AUDIT["registry_sha256"]
    current = json.loads((DATA / "catalog/current.audit.json").read_text("utf-8"))
    assert current["master_key_monomer_admission_sha256"] == importer.digest(DATA / "catalog/master-key-monomer-admission.audit.json")


def test_ground_truth_review_is_qualified_and_new_names_do_not_shadow_old_names():
    report = json.loads((DATA / "catalog/master-key-ground-truth-review.audit.json").read_text("utf-8"))
    assert report["catalog_sha256"] == AUDIT["catalog_sha256"]
    assert report["admission_audit_sha256"] == importer.digest(DATA / "catalog/master-key-monomer-admission.audit.json")
    assert report["annotation_scan"]["rows"] == report["annotation_scan"]["ground:False"] == 2489
    assert report["ground_truth"] is False and report["human_gold_found_in_inspected_records"] is False
    assert [r["unreviewed_count"] for r in report["human_review_queues"]] == [20, 100]
    assert len(report["excluded_master_matched_candidates"]) == 10
    assert len(report["existing_monomer_reconciliation"]) == 63
    assert not re.search(r'(?:^|["\s])[A-Za-z]:[\\/]|/Users/|/home/|source_page_text', json.dumps(report))
    normalize = lambda s: re.sub(r"[^a-z0-9]+", "", s.lower())
    aliases = lambda row: {normalize(s) for s in [row["abbreviation"], row["full_name"], *(row.get("literature_names") or [])]}
    old_names = set().union(*(aliases(row) for row in CATALOG["monomer_presets"][:63]))
    assert not set().union(*(aliases(ROWS[key]) for key in AUDIT["added_ids"])) & old_names


@pytest.mark.parametrize("a,b,linkage,template", [
    ("tfpc", "pa", "imine", "hcb_3_2"),
    ("aminal_a1", "pa", "imine", "sql_4_2"),
    ("ta", "tabpb", "imine", "hcb_3_2"),
    ("dba12_cho", "pa", "imine", "hcb_3_2"),
    ("bfbapy", "tapb", "imine", "hcb_3_2"),
    ("tfb", "nda15", "imine", "hcb_3_2"),
    ("ta", "tbt", "imine", "hcb_3_2"),
    ("aminal_a2", "pa", "imine", "sql_4_2"),
    ("bpp", "tapb", "imine", "hcb_3_2"),
    ("tfb", "dhzdpr", "acylhydrazone", "hcb_3_2"),
    ("bpf", "tapb", "imine", "hcb_3_2"),
    ("tfb", "btbba", "imine", "hcb_3_2"),
    ("tfb", "btbbav", "imine", "hcb_3_2"),
    ("tfb", "hth_s", "acylhydrazone", "hcb_3_2"),
    ("tfb", "hth_r", "acylhydrazone", "hcb_3_2"),
])
def test_each_supported_monomer_builds_a_scoped_design(a, b, linkage, template):
    recipe = prepare_template(ROWS[a]["smiles"], ROWS[b]["smiles"], linkage, a, b, 6, template)
    recipe.update(schema_version=VERSION, recipe_id="master-reviewed-design", designed_variant=True)
    mol, audit = build_patch(recipe)
    mol.UpdatePropertyCache(strict=True)
    assert len(Chem.GetMolFrags(mol)) == 1 and audit["atoms"] == mol.GetNumAtoms()
    assert min(recipe["compatibility"]["reactive_site_counts"]) >= 2
    assert not recipe.get("cofex_id")


def test_chiral_pair_is_not_collapsed_and_mixed_precursors_remain_reference_only():
    for key, expected in (("hth_s", "S"), ("hth_r", "R")):
        mol = Chem.MolFromSmiles(ROWS[key]["smiles"])
        assert [label for _, label in Chem.FindMolChiralCenters(mol)] == [expected, expected]
    assert ROWS["hth_s"]["smiles"] != ROWS["hth_r"]["smiles"]
    assert ROWS["hth_r"]["source_evidence"][0]["upstream_image_hold_reasons"] == ["unmasked_neighbor_reaction_arrow_context"]
    for key in ("aph", "ash", "flt_boronic"):
        assert ROWS[key]["automatic_prepare_supported"] is False
        assert max((r["reactive_sites"] for r in ROWS[key]["compatible_roles"]), default=0) < 2
    assert "tfpyripy_ac" not in ROWS


@pytest.mark.parametrize("field,value", [("identity_held", True), ("chemical_identity_held", True),
    ("source_name_structure_conflict", "unresolved"), ("ground_truth", True), ("human_verified", True)])
def test_identity_holds_and_human_claims_cannot_enter_provisional_lane(field, value):
    row = {"ground_truth": False, "molecules": [{"visibility": "complete"}], field: value}
    with pytest.raises(ValueError):
        importer.eligible_source(row)


def test_master_graph_match_is_stereo_and_doi_specific():
    connection = sqlite3.connect(":memory:")
    connection.executescript("CREATE TABLE records(record_id INTEGER,doi TEXT); CREATE TABLE record_monomers(record_id INTEGER,smiles TEXT); CREATE TABLE components(component_no TEXT,canonical_smiles TEXT);")
    selection = deepcopy(next(r for r in SELECTION if r["id"] == "hth_s"))
    connection.execute("INSERT INTO records VALUES (?,?)", (322, selection["source_doi"]))
    connection.execute("INSERT INTO record_monomers VALUES (?,?)", (322, selection["canonical_smiles"]))
    connection.execute("INSERT INTO components VALUES (?,?)", ("1143", selection["canonical_smiles"]))
    assert importer.master_binding(connection, selection, {})["record_ids"] == [322]
    selection["canonical_smiles"] = ROWS["hth_r"]["smiles"]
    with pytest.raises(ValueError, match="precursor graph"):
        importer.master_binding(connection, selection, {})
    selection["canonical_smiles"] = ROWS["hth_s"]["smiles"]
    selection["source_doi"] = "10.0000/wrong-paper"
    with pytest.raises(ValueError, match="precursor graph"):
        importer.master_binding(connection, selection, {})
