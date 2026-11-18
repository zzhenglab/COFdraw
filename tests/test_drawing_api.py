import hashlib
import json
from copy import deepcopy

import pytest
from fastapi.testclient import TestClient
from rdkit import Chem

from cof_vision.cofdraw_desktop import create_desktop_app as create_app


@pytest.fixture
def client(tmp_path):
    with TestClient(create_app(tmp_path)) as active:
        yield active


def prepare(client, keto=False):
    response = client.post("/api/drawing/prepare", json={
        "linkage": "beta_ketoenamine" if keto else "imine",
        "monomer_a": "O=Cc1c(O)c(C=O)c(O)c(C=O)c1O" if keto else "O=Cc1cc(C=O)cc(C=O)c1",
        "monomer_b": "Nc1ccc(N)cc1", "instances": 12 if keto else 6,
        "label_a": "Tp" if keto else "TFB", "label_b": "Pa",
    })
    assert response.status_code == 200, response.text
    return response.json()["recipe"]


@pytest.mark.parametrize("keto", [False, True])
def test_actual_prepare_render_and_four_downloads(client, keto):
    recipe = prepare(client, keto)
    assert recipe["layout_refinement"] == "bond_angle_pose_v2"
    response = client.post("/api/drawing/render", json={"recipe": recipe, "size": 500})
    assert response.status_code == 200, response.text
    result = response.json()
    assert result["status"] == "designed_variant"
    assert result["audit"]["graph_scope"] == "finite_designed_pore_patch_with_continuation_ports"
    assert result["audit"]["cofex_id"] is None
    assert result["audit"]["ground_truth"] is False
    assert result["audit"]["layout_refinement"]["method"] == "bond_angle_pose_v2"
    assert len(result["audit"]["continuation_ports"]) == (6 if keto else 3)
    assert result["audit"]["building_block_instances"] == (12 if keto else 6)
    outputs = {ext: client.get(url) for ext, url in result["downloads"].items()}
    assert all(response.status_code == 200 for response in outputs.values())
    assert "<svg" in outputs["svg"].text
    assert outputs["png"].content.startswith(b"\x89PNG")
    mol = Chem.MolFromMolBlock(outputs["mol"].text)
    assert mol is not None and len(Chem.GetMolFrags(mol)) == 1
    assert outputs["json"].json()["recipe"]["components"]["A"]["label"] == recipe["components"]["A"]["label"]
    assert outputs["json"].json()["audit"]["whole_cof_smiles"] is None


def test_registry_source_identity_is_hash_bound_and_edits_lose_identity(client, tmp_path):
    recipe = prepare(client)
    recipe.pop("layout_refinement")  # A pre-existing source recipe is not upgraded.
    recipe.update(recipe_id="registered-pore", cofex_id="COFEXM-SOURCE", source={"doi": "10.1000/example"}, designed_variant=False)
    directory = tmp_path / "framework-drawing"
    directory.mkdir()
    path = directory / "registered-pore.json"
    path.write_text(json.dumps(recipe), encoding="utf-8")
    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    registry = {"recipes": [{"recipe_id": recipe["recipe_id"], "cofex_id": recipe["cofex_id"], "path": path.name, "sha256": digest}]}
    (directory / "registry.json").write_text(json.dumps(registry), encoding="utf-8")
    assert client.get("/api/drawing/library").json()["source_recipes"][0]["sha256"] == digest
    source = client.post("/api/drawing/render", json={"source_recipe_id": "registered-pore", "size": 500}).json()
    assert source["audit"]["cofex_id"] == "COFEXM-SOURCE"
    assert source["audit"]["registered_recipe_sha256"] == digest
    assert source["status"] == "source_bound_provisional"
    assert source["audit"]["layout_refinement"] == {"method": None, "applied": False}
    edited = client.post("/api/drawing/render", json={"recipe": recipe, "size": 500}).json()
    assert edited["audit"]["cofex_id"] is None
    assert edited["audit"]["source"] is None
    assert edited["status"] == "designed_variant"
    assert hashlib.sha256(path.read_bytes()).hexdigest() == digest
    path.write_text("{}", encoding="utf-8")
    bad = client.post("/api/drawing/render", json={"source_recipe_id": "registered-pore"})
    assert bad.status_code == 422 and "hash mismatch" in bad.text


def test_bad_geometry_ports_and_size_are_rejected(client, tmp_path):
    recipe = prepare(client)
    for mutation in ["coincident", "bad_port", "duplicate_port", "too_large", "bad_instance"]:
        value = deepcopy(recipe)
        if mutation == "coincident":
            value["instances"][1]["layout_position"] = value["instances"][0]["layout_position"]
        elif mutation == "bad_port":
            value["instances"][0]["incoming_port"] = 999
        elif mutation == "duplicate_port":
            value["connections"].append(value["connections"][0])
        elif mutation == "too_large":
            value["components"]["A"]["fragment_smiles"] = "C" * 400 + "C([*:1])C[*:2]"
        else:
            value["instances"][0] = "invalid"
        response = client.post("/api/drawing/render", json={"recipe": value, "size": 500})
        assert response.status_code == 422, (mutation, response.text)
    assert client.post("/api/drawing/render", json={"recipe": recipe, "source_recipe_id": "x"}).status_code == 422
    assert client.post("/api/drawing/render", json={"recipe": recipe, "size": 5000}).status_code == 422
    assert not list(tmp_path.rglob("*.svg"))


def test_user_filenames_cannot_escape_export_folder(client, tmp_path):
    recipe = prepare(client)
    recipe["recipe_id"] = "../../../escaped"
    response = client.post("/api/drawing/render", json={"recipe": recipe, "size": 500})
    assert response.status_code == 200, response.text
    assert not (tmp_path / "escaped.svg").exists()
    assert len(list((tmp_path / "framework-drawing" / "exports").rglob("*.svg"))) == 1
    assert client.get("/api/drawing/exports/not-an-id/svg").status_code == 404
    assert client.get("/api/drawing/exports/" + "a" * 32 + "/exe").status_code == 422


def test_drawing_assets_available_and_cross_site_writes_blocked(client):
    for path in ("/cofdraw.html", "/cofdraw.js", "/cofdraw.css"):
        assert client.get(path).status_code == 200
    assert client.post("/api/drawing/prepare", json={}, headers={"Origin": "https://elsewhere.example"}).status_code == 403


def test_registry_cannot_read_outside_data_root(client, tmp_path):
    directory = tmp_path / "framework-drawing"
    directory.mkdir()
    (directory / "registry.json").write_text(json.dumps({"recipes": [{
        "recipe_id": "outside", "cofex_id": None, "path": "../../outside.json", "sha256": "0" * 64,
    }]}), encoding="utf-8")
    response = client.get("/api/drawing/sources/outside")
    assert response.status_code == 422
    assert "outside this data workspace" in response.text
    library = client.get("/api/drawing/library").json()
    assert library["source_recipes"] == []
    assert len(library["source_registry_errors"]) == 1


@pytest.mark.parametrize("template,instances,b", [
    ("sql_4_2", 8, "Nc1ccc(N)cc1"),
    ("sql_4_4", 4, "Nc1cc(N)c(N)cc1N"),
])
def test_square_template_instance_counts_reach_chemical_engine(client, template, instances, b):
    # Regression: the HTTP model previously rejected valid 4/8-block templates
    # before the chemistry helper could validate their reactive-site counts.
    response = client.post("/api/drawing/prepare", json={
        "linkage": "imine", "monomer_a": "O=Cc1cc(C=O)c(C=O)cc1C=O",
        "monomer_b": b, "template_id": template, "instances": instances,
    })
    assert response.status_code == 200, response.text
    recipe = response.json()["recipe"]
    assert len(recipe["instances"]) == instances
    assert recipe["design_topology"] == "sql"
    assert recipe["cofex_id"] is None and recipe["topology"] is None
