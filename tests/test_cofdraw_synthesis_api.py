"""SA API results bind the actual finite graphs and persist with exported drawings."""
from pathlib import Path
import shutil

from fastapi.testclient import TestClient
import pytest

from cof_vision.cofdraw_desktop import create_desktop_app
from cof_vision.cofdraw_knowledge import capability_inventory
from cof_vision.cofdraw_synthesis import synthetic_accessibility


@pytest.fixture(scope="module")
def client(tmp_path_factory):
    root = tmp_path_factory.mktemp("sa-api")
    shutil.copytree(Path(__file__).resolve().parents[1] / "data/framework-drawing", root / "framework-drawing")
    with TestClient(create_desktop_app(root)) as active:
        yield active


def post(client, path, value):
    response = client.post("/api/drawing/" + path, json=value)
    assert response.status_code == 200, response.text
    return response.json()


def test_catalog_and_science_library_score_all_actual_monomers_with_same_support(client):
    library = client.get("/api/drawing/library").json()
    report = client.get("/api/drawing/science/synthesis-library").json()
    rows = {r["id"]: r for r in report["monomers"]}
    assert report["count"] == report["available_count"] == len(library["monomer_presets"]) == 81
    for monomer in library["monomer_presets"]:
        score = monomer["synthetic_accessibility"]
        assert score == rows[monomer["id"]]["synthetic_accessibility"]
        assert score["available"] and 1 <= score["score"] <= 10
        assert score["canonical_smiles"] == synthetic_accessibility(monomer["smiles"])["canonical_smiles"]
    # Having a score must not promote reference-only boronate chemistry.
    assert library["linkages"]["boronate_ester"]["supported"] is False
    assert rows["tbpm_boronic"]["synthetic_accessibility"]["available"]


def test_monomer_api_scores_current_graph_and_reports_unsupported_graphs_honestly(client):
    ethanol = post(client, "monomer", {"smiles": "CCO"})
    benzene = post(client, "monomer", {"smiles": "c1ccccc1"})
    assert ethanol["synthetic_accessibility"]["score"] == pytest.approx(1.9802570386349831)
    assert benzene["synthetic_accessibility"]["score"] == 1
    assert ethanol["synthetic_accessibility"]["canonical_smiles"] != benzene["synthetic_accessibility"]["canonical_smiles"]
    incomplete = post(client, "monomer", {"smiles": "*c1ccccc1"})
    assert incomplete["synthetic_accessibility"]["available"] is False
    assert incomplete["synthetic_accessibility"]["score"] is None


def test_fresh_design_render_scores_declared_precursors_and_exports_identical_audit(client):
    request = {"monomer_a": "O=Cc1cc(C=O)cc(C=O)c1", "monomer_b": "Nc1ccc(N)cc1",
               "label_a": "TFB", "label_b": "Pa", "linkage": "imine", "template_id": "hcb_3_2"}
    recipe = post(client, "prepare", request)["recipe"]
    result = post(client, "render", {"recipe": recipe, "size": 500})
    audit = result["audit"]["precursor_synthetic_accessibility"]
    assert audit["complete"] and audit["total_count"] == audit["available_count"] == 2
    scores = [synthetic_accessibility(request[key])["score"] for key in ("monomer_a", "monomer_b")]
    assert audit["mean"] == pytest.approx(sum(scores) / 2)
    assert audit["max"] == max(scores)
    assert {r["smiles"] for r in audit["precursors"]} == {request["monomer_a"], request["monomer_b"]}
    assert "infinite COF" in audit["scope"]
    exported = client.get(result["downloads"]["json"]).json()
    assert exported["audit"]["precursor_synthetic_accessibility"] == audit
    assert exported["recipe"]["components"] == recipe["components"]


def test_explicit_source_has_organic_precursor_scores_without_chemical_reassignment(client):
    result = post(client, "render", {"source_recipe_id": "cof-bp-ni1-network-v1", "size": 500})
    assert result["status"] == "source_bound_provisional"
    score = result["audit"]["precursor_synthetic_accessibility"]
    assert score["complete"] and score["total_count"] == 2
    assert all("Ni" not in row["smiles"] for row in score["precursors"])
    assert "Precursor SA Score summary only" in score["scope"]
    assert "polymerization prediction" in " ".join(score["limitations"])


def test_scientific_monomer_report_persists_score_with_inputs(client):
    result = post(client, "science/monomer", {"smiles": "CCO", "conformers": 1})
    assert result["synthetic_accessibility"]["score"] == pytest.approx(1.9802570386349831)
    assert result["inputs"]["smiles"] == "CCO"
    assert client.get(result["downloads"]["json"]).json()["synthetic_accessibility"] == result["synthetic_accessibility"]
    assert client.get("/api/drawing/science/reports/" + result["report_id"]).json() == result


def test_capability_discovery_distinguishes_geometry_and_finite_sa_scope():
    rows = {r["record_id"]: r for r in capability_inventory()}
    geometry = rows["periodic-geometry"]["data"]
    assert geometry["endpoint"] == "/api/drawing/science/geometry"
    assert "grid_spacing" in geometry["request_schema"]["properties"]
    assert "Coarse grids can leave transport unresolved" in geometry["scope"]
    sa = rows["precursor-sa-library"]["data"]
    assert sa["method"] == "GET" and sa["endpoint"].endswith("/synthesis-library")
    assert "not ASKCOS" in sa["scope"]
