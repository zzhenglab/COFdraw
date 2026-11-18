"""Bundled reticular designs must remain hash-bound and geometrically usable."""
from collections import Counter
import hashlib
import itertools
import json
from pathlib import Path
import re

import numpy as np
import pytest

from cof_vision.cofdraw_science import import_cif, structure_exports
from cof_vision.cofdraw_network import inspect_periodic_network

ROOT = Path(__file__).resolve().parents[1] / "data/crystal-models"
MANIFEST = json.loads((ROOT / "manifest.json").read_text("utf-8"))

# Independent expected product stoichiometry and periodic adjacency. Keep these
# separate from manifest metadata so a coordinated metadata/CIF error fails.
EXPECTED_PRODUCTS = {
    "hcb-imine-aa": ({"C": 36, "H": 24, "N": 6}, 2, 1, 66, 72),
    "hcb-imine-ab": ({"C": 72, "H": 48, "N": 12}, 2, 2, 132, 144),
    "hcb-expanded-aa": ({"C": 72, "H": 48, "N": 6}, 2, 1, 126, 138),
    "dia-imine": ({"C": 82, "H": 56, "N": 8}, 3, 1, 146, 160),
    "hcb-tfb-tapb-aa": ({"C": 33, "H": 21, "N": 3}, 2, 1, 57, 63),
}


@pytest.mark.parametrize("record", MANIFEST["templates"], ids=lambda row: row["id"])
def test_crystal_template_product_graph_formula_and_supercell_roundtrip(record):
    path = ROOT / record["path"]
    assert path.resolve().is_relative_to(ROOT.resolve())
    assert hashlib.sha256(path.read_bytes()).hexdigest() == record["sha256"]
    structure = import_cif(path.read_text("utf-8"))
    assert len(structure["atoms"]) == record["atoms"]
    counts = Counter(atom["element"] for atom in structure["atoms"])
    formula = {element: int(count or 1) for element, count in re.findall(r"([A-Z][a-z]?)(\d*)", record["formula"])}
    assert counts == formula
    expected_formula, dimension, components, atoms, edges = EXPECTED_PRODUCTS[record["id"]]
    assert counts == expected_formula
    assert counts["H"] > 0
    network = inspect_periodic_network(structure)
    assert network["dimensions"] == [dimension]
    assert network["quotient_component_count"] == components
    assert network["atom_count"] == atoms and network["edge_count"] == edges
    assert network["short_contact_count"] == 0 and network["isolated_atom_indices"] == []
    assert not network["threshold_sensitivity"]["classification_changed"]
    assert all(row["coordination"] == 1 for row in network["coordination"]["per_atom"] if row["element"] == "H")
    assert network["coordination"]["unusual_atoms"] == []
    cell = np.asarray(structure["cell"])
    fraction = np.array([atom["fractional"] for atom in structure["atoms"]])
    delta = fraction[:, None, :] - fraction[None, :, :]
    delta -= np.rint(delta)
    minimum = np.inf
    for image in itertools.product(range(-2, 3), repeat=3):
        distances = np.linalg.norm((delta + image) @ cell, axis=2)
        if image == (0, 0, 0):
            np.fill_diagonal(distances, np.inf)
        minimum = min(minimum, float(distances.min()))
    assert minimum > .7  # severe atom overlaps, not an energy/relaxation criterion
    assert minimum == pytest.approx(record["minimum_atom_distance_angstrom"], abs=2e-5)
    exported = structure_exports(structure, repeats=(2, 2, 2))
    assert len(exported["structure"]["atoms"]) == record["atoms"]
    assert len(exported["viewer"]["atoms"]) == record["atoms"] * 8
    roundtrip = import_cif(exported["exports"]["cif"])
    assert len(roundtrip["atoms"]) == record["atoms"] * 8
    assert Counter(atom["element"] for atom in roundtrip["atoms"]) == {element: count * 8 for element, count in counts.items()}
    assert roundtrip["cell_volume_angstrom3"] == pytest.approx(structure["cell_volume_angstrom3"] * 8, rel=2e-9)
