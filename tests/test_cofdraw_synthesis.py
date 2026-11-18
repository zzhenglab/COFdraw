"""Published SA Score reference, conservative inputs and offline package resources."""
from copy import deepcopy
import importlib.util
import json
from pathlib import Path
import shutil
import socket

import pytest
from rdkit.Contrib.SA_Score import sascorer

from cof_vision import cofdraw_synthesis as synthesis
from cof_vision.cofdraw_core import read_catalog


REPO = Path(__file__).resolve().parents[1]
DATA = REPO / "data/framework-drawing"
CATALOG = read_catalog(REPO / "data")
# Exact example and score published by the RDKit maintainer on 2023-12-01:
# https://greglandrum.github.io/rdkit-blog/posts/2023-12-01-using_sascore_and_npscore.html
REFERENCE_SMILES = "FC1=CN=C(NC2=NC=C(C(N3CCN(C(C4CC4)=O)CC3)=O)C=C2)N=C1C5=CC=C6N=C(N(CC)CC)SC6=C5"


def test_published_rdkit_reference_value_and_qualified_metadata():
    result = synthesis.synthetic_accessibility(REFERENCE_SMILES)
    assert result["available"] is True
    assert result["score"] == pytest.approx(2.8716389090191434, abs=1e-12)
    assert result["method_id"] == "rdkit_sa_score"
    assert result["scale"] == {"min": 1, "max": 10, "direction": "lower_is_easier"}
    assert result["fragment_table_verified"] is True
    assert result["fragment_scores_sha256"] == synthesis.FRAGMENT_SCORES_SHA256
    assert result["rdkit_version"] and result["canonical_smiles"]
    assert "not a probability" in " ".join(result["limitations"])
    assert "Not ASKCOS" in " ".join(result["limitations"])


@pytest.mark.parametrize("smiles,code", [
    (None, "invalid_input"), ("", "invalid_input"), ("C" * 4001, "invalid_input"),
    ("not a SMILES", "invalid_graph"), ("C(C)(C)(C)(C)C", "invalid_graph"),
    ("CC ethanol", "invalid_graph"), ("C |$foo$|", "invalid_graph"),
    ("C" * 251, "atom_limit"), ("*c1ccccc1", "continuation_atoms"),
    ("CCO.CCO", "disconnected_graph"), ("[Na+].[Cl-]", "disconnected_graph"),
    ("[Fe]", "unsupported_elements"), ("[CH3]", "radical_graph"),
    ("N->[B-](F)(F)F", "unsupported_bonds"),
])
def test_unsupported_inputs_have_no_numeric_score_or_silent_graph_repair(smiles, code):
    result = synthesis.synthetic_accessibility(smiles)
    assert result["available"] is False and result["score"] is None
    assert result["reason_code"] == code and result["reason"]


def test_chemically_equivalent_smiles_and_caller_mutation_do_not_change_cached_scores():
    first = synthesis.synthetic_accessibility("CCO")
    assert first["score"] == pytest.approx(1.9802570386349831, abs=1e-12)
    first["score"] = 900
    first["limitations"].clear()
    first["scale"]["min"] = -999
    second = synthesis.synthetic_accessibility("OCC")
    assert second["score"] == pytest.approx(1.9802570386349831, abs=1e-12)
    assert second["limitations"] and second["scale"]["min"] == 1
    assert second["canonical_smiles"] == "CCO"


def test_isotope_and_formal_charge_are_retained_without_cost_prediction():
    labeled = synthesis.synthetic_accessibility("[13CH3]CO")
    charged = synthesis.synthetic_accessibility("C[NH3+]")
    assert labeled["available"] and "13" in labeled["canonical_smiles"]
    assert labeled["isotope_labels_present"] is True
    assert "Isotope-label" in " ".join(labeled["limitations"])
    assert charged["available"] and charged["formal_charge"] == 1
    assert "+" in charged["canonical_smiles"]


def test_all_catalog_monomers_score_offline_and_keep_qualified_reference_entries(monkeypatch):
    def no_network(*args, **kwargs):
        pytest.fail("An offline SA Score must not use the network")
    monkeypatch.setattr(socket, "create_connection", no_network)
    rows = deepcopy(CATALOG["monomer_presets"])
    result = synthesis.catalog_scores(rows)
    assert result["count"] == result["available_count"] == len(rows) == 118
    assert {r["id"] for r in result["monomers"]} == {r["id"] for r in rows}
    assert all(1 <= r["synthetic_accessibility"]["score"] <= 10 for r in result["monomers"])
    original = {r["id"]: r for r in rows}
    for row in result["monomers"]:
        assert row["smiles"] == original[row["id"]]["smiles"]
        assert row.get("compatible_roles") == original[row["id"]].get("compatible_roles")
    assert CATALOG["monomer_presets"] == rows
    assert "not reaction support" in result["scope"]


def test_assembled_precursor_summary_is_unweighted_and_ignores_product_fragments():
    recipe = {"components": {
        "a": {"monomer_smiles": "CCO", "fragment_smiles": "*CCO"},
        "a_alias": {"monomer_smiles": "OCC"},
        "b": {"monomer_smiles": "c1ccccc1"}},
        "instances": [{"component": "a"}] * 20 + [{"component": "b"}],
        "patch_smiles": "*CC*"}
    before = deepcopy(recipe)
    result = synthesis.precursor_scores(recipe)
    assert recipe == before
    assert result["component_count"] == 3
    assert result["available_count"] == result["total_count"] == 2 and result["complete"]
    assert result["precursors"][0]["components"] == ["a", "a_alias"]
    assert result["mean"] == pytest.approx((1.9802570386349831 + 1) / 2)
    assert result["max"] == pytest.approx(1.9802570386349831)
    assert "not a score of the drawing" in result["scope"]


def test_explicit_source_uses_declared_precursors_not_drawn_atoms():
    path = DATA / "recovered/brick-wall-pore-v1.recipe.json"
    recipe = json.loads(path.read_text("utf-8"))
    assert recipe["explicit_graph"]
    result = synthesis.precursor_scores(recipe)
    assert result["complete"] and result["total_count"] == 2
    assert {r["component"] for r in result["precursors"]} == set(recipe["components"])
    recipe["explicit_graph"] = {"atoms": [{"symbol": "*"}], "bonds": []}
    assert synthesis.precursor_scores(recipe) == result


@pytest.mark.parametrize("components", [
    {}, {"a": {"fragment_smiles": "CCO"}},
    {"a": {"monomer_smiles": "CCO"}, "b": {"monomer_smiles": "*N"}},
    {"a": {"monomer_smiles": "CCO"}, "b": {"monomer_smiles": "[Na+].[Cl-]"}},
])
def test_incomplete_precursors_never_produce_partial_numeric_aggregate(components):
    result = synthesis.precursor_scores({"components": components})
    assert result["complete"] is False and result["mean"] is None and result["max"] is None
    assert result["reason"]
    assert result["available_count"] < result["total_count"] or result["total_count"] == 0


def test_installed_table_verification_and_tampered_or_missing_table_rejected(tmp_path):
    result = synthesis.resource_information()
    assert result["sha256"] == synthesis.FRAGMENT_SCORES_SHA256
    assert result["bytes"] == synthesis.FRAGMENT_SCORES_BYTES and result["verified"]
    with pytest.raises(RuntimeError, match="missing or failed verification"):
        synthesis.resource_information(tmp_path)
    # The correct length alone cannot admit another pickle payload.
    table = tmp_path / "fpscores.pkl.gz"
    table.write_bytes(b"x" * synthesis.FRAGMENT_SCORES_BYTES)
    with pytest.raises(RuntimeError, match="missing or failed verification"):
        synthesis.resource_information(tmp_path)


def test_unverified_table_is_never_deserialized_and_score_becomes_unavailable(monkeypatch):
    synthesis._score_canonical.cache_clear()
    monkeypatch.setattr(synthesis, "_SCORER", None)
    def unavailable():
        raise RuntimeError("test resource mismatch")
    def forbidden_load(*args, **kwargs):
        pytest.fail("Unverified table must never reach pickle deserialization")
    monkeypatch.setattr(synthesis, "resource_information", unavailable)
    monkeypatch.setattr(sascorer, "readFragmentScores", forbidden_load)
    result = synthesis.synthetic_accessibility("CCO")
    assert result["available"] is False and result["score"] is None
    assert result["reason_code"] == "resource_unavailable"
    assert "test resource mismatch" not in json.dumps(result)


def test_build_gate_checks_staged_scorer_and_table_and_distributes_license(tmp_path):
    spec = importlib.util.spec_from_file_location("synthesis_build_gate", REPO / "scripts/build_cofdraw_desktop.py")
    builder = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(builder)
    with pytest.raises(SystemExit, match="verified RDKit"):
        builder.validate_synthesis_resources(tmp_path)
    package = Path(sascorer.__file__).parent
    for name in ("sascorer.py", "fpscores.pkl.gz"):
        shutil.copy2(package / name, tmp_path / name)
    assert builder.validate_synthesis_resources(tmp_path)["verified"] is True
    notice = (REPO / "docs/third-party-licenses/RDKit-SA-Score-BSD.txt").read_text("utf-8")
    assert "Novartis Institutes for BioMedical Research Inc." in notice
    assert "Redistributions in binary form must reproduce" in notice
    assert "THIS SOFTWARE IS PROVIDED" in notice


@pytest.mark.parametrize("function,arg", [
    (synthesis.catalog_scores, {}), (synthesis.catalog_scores, [None]),
    (synthesis.catalog_scores, [{}] * 1001), (synthesis.precursor_scores, []),
    (synthesis.precursor_scores, {"components": []}),
    (synthesis.precursor_scores, {"components": {str(i): {} for i in range(101)}}),
])
def test_public_collection_inputs_are_bounded(function, arg):
    with pytest.raises(ValueError):
        function(arg)
