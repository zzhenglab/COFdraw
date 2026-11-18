"""Local retrieval preserves evidence and does not turn references into support."""
from copy import deepcopy
import hashlib
import json
from pathlib import Path
import shutil

import pytest

from cof_vision.cofdraw_knowledge import capability_inventory, query_catalog
from cof_vision.drawing_api import KnowledgeQuery, PowderDrawing


DATA = Path(__file__).resolve().parents[1] / "data"


def test_name_punctuation_unicode_and_doi_find_exact_source_payload():
    catalog = json.loads((DATA / "framework-drawing/catalog/current.json").read_text(encoding="utf-8"))
    answers = [query_catalog(DATA, q, ["source"]) for q in
               ["Pyurea–COF", "show me Pyurea COF", "https://doi.org/10.1039/d1qm00231g", "DOI: 10.1039/d1qm00231g"]]
    for answer in answers:
        assert answer["matched_count"] == 1
        assert not answer["errors"]
        row = answer["records"][0]
        assert row["id"] == "source:pyurea-source-layout-v1"
        assert row["data"]["metadata"] == catalog["source_metadata"]["pyurea-source-layout-v1"]
        assert row["support"]["source_registered"] is False
        assert row["data"]["metadata"]["cofex_id"] is None
        assert row["provenance"]["recipe_sha256"] == row["data"]["registry"]["sha256"]
        assert row["retrieval"]["path"] == "/api/drawing/sources/pyurea-source-layout-v1"


def test_linkage_support_uses_runtime_and_alias_is_not_automatic_equivalence():
    beta = query_catalog(DATA, "β-ketoenamine", ["linkage"])["records"][0]
    assert beta["record_id"] == "beta_ketoenamine"
    assert beta["support"]["automatic_prepare_supported"] is True
    boronate = query_catalog(DATA, "boronate ester", ["linkage"])["records"][0]
    assert boronate["support"]["status"] == "reference_only"
    assert boronate["support"]["automatic_prepare_supported"] is False
    hydrazone = next(r for r in query_catalog(DATA, "hydrazone", ["linkage"])["records"] if r["record_id"] == "hydrazone")
    assert hydrazone["support"]["automatic_prepare_supported"] is False


def test_topology_templates_and_crystals_keep_distinct_scope():
    hcb = query_catalog(DATA, "hcb", ["topology"])["records"][0]
    assert set(hcb["support"]["template_ids"]) == {"hcb_3_2", "hcb_3_3"}
    assert hcb["support"]["periodic_net_assignment"] is False
    crystals = query_catalog(DATA, "", ["crystal"])
    assert crystals["matched_count"] == 4
    assert all(r["support"]["status"] == "unrelaxed_periodic_design" for r in crystals["records"])
    assert all(r["support"]["experimentally_validated"] is False for r in crystals["records"])


def test_query_is_deterministic_paginates_and_unknown_query_has_no_invented_answer():
    first = query_catalog(DATA, "", ["monomer"], limit=7)
    second = query_catalog(DATA, "", ["monomer"], limit=7, offset=7)
    assert first == query_catalog(DATA, "", ["monomer"], limit=7)
    catalog = json.loads((DATA / "framework-drawing/catalog/current.json").read_text(encoding="utf-8"))
    assert first["matched_count"] == len(catalog["monomer_presets"])
    assert first["next_offset"] == 7
    assert not {r["id"] for r in first["records"]} & {r["id"] for r in second["records"]}
    missing = query_catalog(DATA, "unlistedmaterialxyz987")
    assert missing["records"] == [] and missing["matched_count"] == 0
    assert "not evidence" in missing["no_match_meaning"]


def test_returned_records_are_independent_and_capability_schemas_match_runtime():
    row = query_catalog(DATA, "tfb", ["monomer"])["records"][0]
    original = deepcopy(row)
    row["data"]["smiles"] = "invalid"
    assert query_catalog(DATA, "tfb", ["monomer"])["records"][0] == original
    capabilities = {r["record_id"]: r for r in capability_inventory()}
    assert capabilities["local-query"]["data"]["request_schema"] == KnowledgeQuery.model_json_schema()
    assert capabilities["powder-diffraction"]["data"]["request_schema"] == PowderDrawing.model_json_schema()
    assert capabilities["saved-report"]["data"]["method"] == "GET"


@pytest.mark.parametrize("kwargs", [{"query": "x" * 501}, {"limit": 51}, {"limit": True},
                                   {"offset": -1}, {"offset": 10001}, {"record_types": ["made_up"]}])
def test_invalid_work_bounds_are_rejected(kwargs):
    with pytest.raises(ValueError):
        query_catalog(DATA, **kwargs)


def test_tampered_source_is_excluded_with_integrity_error(tmp_path):
    shutil.copytree(DATA / "framework-drawing", tmp_path / "framework-drawing")
    shutil.copytree(DATA / "crystal-models", tmp_path / "crystal-models")
    path = tmp_path / "framework-drawing/recovered/pyurea-source-layout-v1.recipe.json"
    path.write_text("{}", encoding="utf-8")
    answer = query_catalog(tmp_path, "Pyurea", ["source"])
    assert answer["records"] == []
    assert any(r.get("record_id") == "pyurea-source-layout-v1" and "hash" in r["reason"] for r in answer["errors"])


def test_registry_escape_is_not_read_even_with_matching_external_hash(tmp_path):
    shutil.copytree(DATA / "framework-drawing", tmp_path / "framework-drawing")
    external = tmp_path / "outside.json"
    external.write_text('{"private":"not returned"}', encoding="utf-8")
    registry_path = tmp_path / "framework-drawing/registry.json"
    registry = json.loads(registry_path.read_text(encoding="utf-8"))
    row = next(r for r in registry["recipes"] if r["recipe_id"] == "pyurea-source-layout-v1")
    row.update(path="../outside.json", sha256=hashlib.sha256(external.read_bytes()).hexdigest())
    registry_path.write_text(json.dumps(registry), encoding="utf-8")
    catalog_path = tmp_path / "framework-drawing/catalog/current.json"
    catalog = json.loads(catalog_path.read_text(encoding="utf-8"))
    catalog["source_registry_sha256"] = hashlib.sha256(registry_path.read_bytes()).hexdigest()
    catalog_path.write_text(json.dumps(catalog), encoding="utf-8")
    answer = query_catalog(tmp_path, "Pyurea", ["source"])
    assert answer["records"] == []
    assert any("escapes" in r["reason"] for r in answer["errors"])
    assert "not returned" not in json.dumps(answer)
