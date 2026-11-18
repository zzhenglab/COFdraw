import csv
import io
from pathlib import Path
import shutil

from fastapi.testclient import TestClient
import pytest

from cof_vision.cofdraw_desktop import create_desktop_app


DATA = Path(__file__).resolve().parents[1] / "data"


def test_lookup_estimates_and_literature_load_use_distinct_amounts(tmp_path):
    shutil.copytree(DATA / "solubility", tmp_path / "solubility")
    catalog = tmp_path / "framework-drawing" / "catalog"
    catalog.mkdir(parents=True)
    shutil.copy2(DATA / "framework-drawing/catalog/current.json", catalog)
    with TestClient(create_desktop_app(tmp_path)) as client:
        capabilities = client.get("/api/drawing/science/capabilities").json()
        path = capabilities["operations"]["solubility"]["path"]
        assert "dioxane" in capabilities["solubility"]["solvents"]
        response = client.post(path, json={"monomer_id": "tfb", "solvent_ids": ["dioxane", "dmf"], "solvent_ml": 3})
        assert response.status_code == 200, response.text
        match = response.json()
        reported, estimate = match["solvents"]
        assert reported["status"] == "reported" and reported["source_doi"].startswith("10.")
        assert 293.15 <= reported["temperature_k"] <= 303.15
        assert reported["mass_mg"] == pytest.approx(3 * reported["mol_per_liter_solvent"] * match["molecular_weight_g_mol"])
        assert estimate["status"] == "estimated"
        assert estimate["screening_band"]["lower"]["mass_mg"] < estimate["mass_mg"] < estimate["screening_band"]["upper"]["mass_mg"]
        assert len(estimate["analogs"]) >= 3 and all(analog["smiles"] != match["canonical_smiles"] for analog in estimate["analogs"])
        saved = client.get(match["downloads"]["json"]).json()
        assert saved["solvents"] == match["solvents"]

        example = client.post(path, json={"monomer_id": "dhp_cho4", "solvent_ids": ["dioxane"], "solvent_ml": 5}).json()
        assert example["solvents"][0]["status"] != "reported"
        loading, = example["preparation_loadings"]
        assert loading["status"] == "reported_loading"
        assert loading["temperature_k"] is None
        assert loading["mass_mg"] == 50 and loading["solvent_ml"] == 2
        assert loading["mg_per_ml_solvent"] == 25
        csv_response = client.get(example["downloads"]["csv"])
        rows = list(csv.DictReader(io.StringIO(csv_response.text)))
        assert rows[-1]["status"] == "reported_loading"
        assert rows[-1]["solvent_ml"] == "2.0"

        new = client.post(path, json={"smiles": "Nc1ccc(-c2cccc(N)c2)cc1", "solvent_ids": ["dmf", "mesitylene"]})
        assert new.status_code == 200, new.text
        first, second = new.json()["solvents"]
        assert first["status"] == "estimated" and new.json()["preparation_loadings"] == []
        assert second["status"] == "unavailable" and "mass_mg" not in second
        metal = client.post(path, json={"smiles": "N->[Zn+2]", "solvent_ids": ["dmf"]}).json()
        assert metal["solvents"][0]["status"] == "unavailable"
        assert client.post(path, json={"smiles": "CC.O"}).status_code == 422
