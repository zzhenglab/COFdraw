import json
from pathlib import Path

import numpy as np
import pytest

from cof_vision.cofdraw_crystal_figures import crystal_details
from cof_vision.cofdraw_network import inspect_periodic_network
from cof_vision.cofdraw_science import import_cif


@pytest.mark.parametrize("filename,has_circuit", [
    ("cof-300.cif",True), ("curated/13040N3.cif",True), ("curated/07011N3.cif",False)])
def test_coordinate_illustrations(filename,has_circuit):
    root=Path(__file__).resolve().parents[1]/"data/crystal-models"
    structure=import_cif((root/filename).read_text())
    original=json.dumps(structure,sort_keys=True)
    details=crystal_details(structure,inspect_periodic_network(structure)["edges"])
    assert len(details["junction"]["atoms"])==25
    assert ("ring" in details)==has_circuit
    if has_circuit:
        assert len(details["overlay"]["guides"])==4
        assert len(details["net"]["bonds"])==16
        ring=details["ring"]
        xyz=np.array([[a[k] for k in ("x","y","z")] for a in ring["atoms"]])
        for a,b in ring["bonds"]:
            assert 1 < np.linalg.norm(xyz[a]-xyz[b]) < 1.9
        assert np.linalg.norm(ring["direction"])==pytest.approx(1)
    assert json.dumps(structure,sort_keys=True)==original
