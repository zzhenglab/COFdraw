import copy
import json

import pytest
from rdkit import Chem
from fastapi.testclient import TestClient

from cof_vision.cofdraw_core import normalize_monomer, STYLES, PRESETS
from cof_vision.cofdraw_desktop import create_desktop_app as create_app


@pytest.mark.parametrize("smiles", ["COc1cc(C=O)c(OC)cc1C=O", "N[C@@H](C)C(=O)O", "F/C=C/F", "Nc1ccc(N)cc1S(=O)(=O)O"])
def test_editor_graph_roundtrip_preserves_structure_and_stereo(smiles):
    imported = normalize_monomer(smiles=smiles)
    graph = copy.deepcopy(imported["graph"])
    for atom in graph["atoms"]:
        atom["x"] += 2.5
        atom["y"] -= 1
    edited = normalize_monomer(graph=graph)
    assert edited["smiles"] == imported["smiles"]
    assert normalize_monomer(molblock=edited["molblock"])["smiles"] == imported["smiles"]


def test_editor_rejects_invalid_valence_and_unconnected_graph():
    with pytest.raises(ValueError):
        normalize_monomer(smiles="CC.CC")
    value = normalize_monomer(smiles="CO")["graph"]
    value["bonds"][0]["order"] = 3
    with pytest.raises(ValueError):
        normalize_monomer(graph=value)


def test_editor_clean_and_mol_preserve_donor_direction_and_metal_charge():
    graph = {
        "atoms": [{"element": "N", "x": 0, "y": 0}, {"element": "Zn", "x": 5, "y": 1, "charge": 2}],
        "bonds": [{"a": 0, "b": 1, "order": 1, "kind": "dative"}],
    }
    clean = normalize_monomer(graph=graph, clean_layout=True)
    assert clean["smiles"] == "[NH3]->[Zn+2]"
    assert clean["graph"]["bonds"][0]["kind"] == "dative"
    assert clean["graph"]["atoms"][1]["charge"] == 2
    assert normalize_monomer(graph=clean["graph"])["smiles"] == clean["smiles"]
    assert normalize_monomer(molblock=clean["molblock"])["smiles"] == clean["smiles"]


def test_every_style_exports_same_chemical_graph_and_random_is_reproducible(tmp_path):
    with TestClient(create_app(tmp_path)) as client:
        # Sulfonated Tp/Pa exercises condensed labels without repeating the
        # expensive symmetric methoxy macrocycle for every presentation style.
        pre = {k: v for k, v in PRESETS[3].items() if k not in ("id", "label")}
        recipe = client.post("/api/drawing/prepare", json=pre).json()["recipe"]
        identities, images = set(), set()
        for style in STYLES:
            result = client.post("/api/drawing/render", json={"recipe": recipe, "size": 600, "style_id": style["id"]})
            assert result.status_code == 200, result.text
            result = result.json()
            assert result["audit"]["cofex_string"] is None
            # A finite port-cut macrocycle is not a molecule-wide aromaticity target.
            mol = Chem.MolFromMolBlock(client.get(result["downloads"]["mol"]).text, sanitize=False, removeHs=False)
            mol.UpdatePropertyCache(strict=True)
            atoms = tuple((a.GetAtomicNum(), a.GetFormalCharge(), a.GetIsotope(), int(a.GetChiralTag())) for a in mol.GetAtoms())
            bonds = tuple((b.GetBeginAtomIdx(), b.GetEndAtomIdx(), str(b.GetBondType()), int(b.GetStereo())) for b in mol.GetBonds())
            identities.add((atoms, bonds))
            images.add(client.get(result["downloads"]["svg"]).text)
            assert client.get(result["downloads"]["png"]).content.startswith(b"\x89PNG")
            if style["id"] == "compact_labels":
                groups = result["audit"]["drawing_style"]["rules"]["collapsed_groups"]
                assert {g["label"] for g in groups} == {"SO3H", "HO3S"}
                assert all(g["attachment_atom"] == "S" for g in groups)
                assert all(g["label"].endswith("S") if g["attachment_side"] == "right" else g["label"].startswith("S") for g in groups)
        assert len(identities) == 1
        assert len(images) == len(STYLES)
        choices = [client.post("/api/drawing/render", json={"recipe": recipe, "size": 500, "style_id": "random", "style_seed": 42}).json()["audit"]["drawing_style"] for _ in range(2)]
        assert choices[0] == choices[1]


def test_catalog_identity_requires_current_recipe_and_exact_named_precursor_set(tmp_path):
    import hashlib
    with TestClient(create_app(tmp_path)) as client:
        pre = {k: v for k, v in PRESETS[0].items() if k not in ("id", "label")}
        recipe = client.post("/api/drawing/prepare", json=pre).json()["recipe"]
        recipe.update(recipe_id="source", cofex_id="COFEXM-TEST")
        directory = tmp_path / "framework-drawing"; directory.mkdir()
        p = directory / "source.json"; p.write_text(json.dumps(recipe))
        digest = hashlib.sha256(p.read_bytes()).hexdigest()
        (directory / "registry.json").write_text(json.dumps({"recipes": [{"recipe_id": "source", "cofex_id": "COFEXM-TEST", "path": "source.json", "sha256": digest}]}))
        notation = pre["monomer_a"] + "." + pre["monomer_b"] + " COFex.unk.cat0;Example"
        metadata = {"cofex_id": "COFEXM-TEST", "source_recipe_sha256": digest, "cofex_notation": notation, "monomers": [{"smiles": pre["monomer_a"]}, {"smiles": pre["monomer_b"]}]}
        catalog = directory / "catalog-v2" / "catalog.json"; catalog.parent.mkdir()
        catalog.write_text(json.dumps({"source_metadata": {"source": metadata}}))
        result = client.post("/api/drawing/render", json={"source_recipe_id": "source", "size": 500}).json()
        assert result["audit"]["cofex_string"] == notation
        recipe.update(cofex_string=notation, cofex_notation=notation)
        edited = client.post("/api/drawing/render", json={"recipe": recipe, "size": 500}).json()
        assert edited["audit"]["cofex_string"] is None
        metadata["monomers"][0]["smiles"] = "CO"
        catalog.write_text(json.dumps({"source_metadata": {"source": metadata}}))
        assert client.get("/api/drawing/sources/source").status_code == 422
def test_editor_preserves_radical_electrons_across_graph_round_trip():
    from cof_vision.cofdraw_core import normalize_monomer
    radical = normalize_monomer(smiles="[O]c1ccccc1")
    assert radical["graph"]["atoms"][0]["radical_electrons"] == 1
    assert normalize_monomer(graph=radical["graph"])["smiles"] == radical["smiles"]
