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

ROOT = Path(__file__).resolve().parents[1] / "data/crystal-models"
MANIFEST = json.loads((ROOT / "manifest.json").read_text("utf-8"))


@pytest.mark.parametrize("record", MANIFEST["templates"], ids=lambda row: row["id"])
def test_crystal_template_hash_formula_no_clashes_and_supercell_roundtrip(record):
    path = ROOT / record["path"]
    assert path.resolve().is_relative_to(ROOT.resolve())
    assert hashlib.sha256(path.read_bytes()).hexdigest() == record["sha256"]
    structure = import_cif(path.read_text("utf-8"))
    assert len(structure["atoms"]) == record["atoms"]
    counts = Counter(atom["element"] for atom in structure["atoms"])
    formula = {element: int(count or 1) for element, count in re.findall(r"([A-Z][a-z]?)(\d*)", record["formula"])}
    assert counts == formula
    assert counts["H"] > 0
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
