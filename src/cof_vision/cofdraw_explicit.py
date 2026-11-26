"""Validated source graph tables with an explicit, unchanged 2D conformer."""
from __future__ import annotations

from copy import deepcopy
import hashlib
import json
import math

import numpy as np
from rdkit import Chem
from rdkit.Geometry import Point3D

from .depiction_rings import symmetric_depiction_rings


ORDERS = {1: Chem.BondType.SINGLE, 2: Chem.BondType.DOUBLE,
          3: Chem.BondType.TRIPLE, 1.5: Chem.BondType.AROMATIC,
          "DATIVE": Chem.BondType.DATIVE}
STEREO = {name: getattr(Chem.BondStereo, name) for name in (
    "STEREONONE", "STEREOANY", "STEREOE", "STEREOZ", "STEREOCIS", "STEREOTRANS")}
CHIRAL = {name: getattr(Chem.ChiralType, name) for name in (
    "CHI_UNSPECIFIED", "CHI_TETRAHEDRAL_CW", "CHI_TETRAHEDRAL_CCW")}


def clear_source_registration(recipe):
    """Mutate registration metadata for an edited recipe; keep history once.

    Return the same dictionary. Atom/bond tables, coordinates and display colors
    are untouched. Historical evidence is explicitly parent-only and must never
    be interpreted as registration of the current user-edited coordinates.
    This also handles assembled recipes with source_coordinate_provenance.
    """
    graph = recipe.get("explicit_graph")
    graph = graph if isinstance(graph, dict) else None
    style = recipe.get("source_style")
    style = style if isinstance(style, dict) else None
    if "source_registration_history" not in recipe:
        history = {"scope": "historical_parent_only"}
        for key in ("source_registered", "source_coordinates_preserved", "source_coordinate_provenance",
                    "coordinate_provenance", "layout_evidence", "source_representation_kind", "source_scope"):
            if key in recipe:
                history[key] = deepcopy(recipe[key])
        if graph:
            history["explicit_graph"] = {key: deepcopy(graph[key]) for key in
                ("source_coordinates_preserved", "coordinate_provenance") if key in graph}
        if style:
            history["source_style"] = {key: deepcopy(style[key]) for key in
                ("source_registered", "exact_style_equivalence", "qualification",
                 "coordinate_evidence_sha256", "review_evidence_sha256") if key in style}
        recipe["source_registration_history"] = history
    qualification = "User-edited coordinates; no active source registration or exact source placement is asserted."
    provenance = {"kind": "user_manual_coordinates", "source_registered": False,
        "source_coordinates_preserved": False, "qualification": qualification,
        "historical_parent_provenance_field": "source_registration_history"}
    recipe["source_registered"] = False
    if "source_coordinates_preserved" in recipe:
        recipe["source_coordinates_preserved"] = False
    for key in ("source_coordinate_provenance", "coordinate_provenance"):
        if key in recipe:
            if key == "coordinate_provenance" and not graph:
                basis = recipe[key]
                if (isinstance(basis, dict) and basis.get('kind') == 'fitted_sheet_face'
                        and (recipe.get('selected_template') or {}).get('id') == 'sql_4_4'
                        and recipe.get('atom_coordinates')):
                    from .cofdraw_chemistry import sheet_face_signature
                    if basis.get('layout_sha256') == sheet_face_signature(recipe):
                        continue
            recipe[key] = deepcopy(provenance)
    if "layout_evidence" in recipe:
        recipe["layout_evidence"] = {"method": "user_manual_coordinates", "exact_atom_registration": False,
                                    "qualification": qualification}
    if "source_representation_kind" in recipe:
        recipe["source_representation_kind"] = "user_edited_finite_patch"
    if "source_scope" in recipe:
        recipe["source_scope"] = "User-edited finite chemical patch; parent source scope is historical only."
    if graph is not None:
        graph["source_coordinates_preserved"] = False
        graph["coordinate_provenance"] = deepcopy(provenance)
    if style is not None:
        style.update(source_registered=False, exact_style_equivalence=False,
            qualification="Inherited presentation styling on a user-edited drawing; no source registration is asserted.")
        style.pop("coordinate_evidence_sha256", None)
        style.pop("review_evidence_sha256", None)
    return recipe


def validate_explicit_recipe(recipe):
    from .framework_drawing import VERSION
    if recipe.get("schema_version") != VERSION:
        raise ValueError("Unsupported explicit drawing recipe version")
    if len(json.dumps(recipe, allow_nan=False)) > 1_500_000:
        raise ValueError("Explicit recipe exceeds 1.5 MB")
    graph = recipe.get("explicit_graph")
    if not isinstance(graph, dict) or graph.get("coordinate_system") != "cartesian_y_up":
        raise ValueError("Explicit graph needs Cartesian y-up coordinates")
    atoms, bonds = graph.get("atoms"), graph.get("bonds")
    if not isinstance(atoms, list) or not 2 <= len(atoms) <= 1000:
        raise ValueError("Explicit graph requires 2–1,000 atoms")
    if not isinstance(bonds, list) or not 1 <= len(bonds) <= 2000:
        raise ValueError("Explicit graph requires 1–2,000 bonds")
    for index, atom in enumerate(atoms):
        if not isinstance(atom, dict) or atom.get("index") != index:
            raise ValueError("Explicit atoms must have consecutive zero-based indices")
        symbol = atom.get("element")
        if not isinstance(symbol, str) or (symbol != "*" and Chem.GetPeriodicTable().GetAtomicNumber(symbol) <= 0):
            raise ValueError("Invalid explicit atom element")
        xy = atom.get("xy")
        if not isinstance(xy, list) or len(xy) != 2 or any(
                isinstance(v, bool) or not isinstance(v, (int, float)) or not math.isfinite(v) or abs(v) > 10000 for v in xy):
            raise ValueError("Explicit atom coordinates must be finite and bounded")
        for field, lower, upper in (("formal_charge", -8, 8), ("isotope", 0, 300),
                                    ("explicit_hs", 0, 8), ("atom_map", 0, 100000)):
            value = atom.get(field, 0)
            if isinstance(value, bool) or not isinstance(value, int) or not lower <= value <= upper:
                raise ValueError(f"Invalid explicit atom {field}")
        if atom.get("chiral_tag", "CHI_UNSPECIFIED") not in CHIRAL:
            raise ValueError("Unsupported explicit atom chirality")
    seen = set()
    for bond in bonds:
        if not isinstance(bond, dict):
            raise ValueError("Explicit bonds must be objects")
        a, b = bond.get("a"), bond.get("b")
        if (type(a) is not int or type(b) is not int or not 0 <= a < len(atoms)
                or not 0 <= b < len(atoms) or a == b):
            raise ValueError("Explicit bond endpoints must identify distinct atoms")
        edge = tuple(sorted((a, b)))
        if edge in seen:
            raise ValueError("Duplicate explicit bond")
        seen.add(edge)
        if bond.get("order") not in ORDERS:
            raise ValueError("Unsupported explicit bond order")
        if math.dist(atoms[a]["xy"], atoms[b]["xy"]) < 1e-6:
            raise ValueError("Bonded atoms cannot occupy the same coordinate")
        if bond.get("stereo", "STEREONONE") not in STEREO:
            raise ValueError("Unsupported explicit bond stereo")
        selectors = bond.get("stereo_atoms", [])
        if not isinstance(selectors, list) or len(selectors) not in (0, 2) or any(
                type(v) is not int or not 0 <= v < len(atoms) for v in selectors):
            raise ValueError("Invalid explicit stereo selectors")
        if bond.get("stereo", "STEREONONE") not in {"STEREONONE", "STEREOANY"} and len(selectors) != 2:
            raise ValueError("Defined double-bond stereo needs two selectors")


def build_explicit_patch(recipe):
    """Construct only the supplied graph; never infer chemistry from its layout."""
    from .framework_drawing import VERSION, _double_bond_geometry_audit
    validate_explicit_recipe(recipe)
    designed = bool(recipe.get("designed_variant"))
    if designed:
        recipe = clear_source_registration(deepcopy(recipe))
    graph = recipe["explicit_graph"]
    rw = Chem.RWMol()
    origins = []
    for row in graph["atoms"]:
        atom = Chem.Atom(0 if row["element"] == "*" else row["element"])
        atom.SetFormalCharge(row.get("formal_charge", 0))
        atom.SetIsotope(row.get("isotope", 0))
        atom.SetIsAromatic(row.get("aromatic", False))
        atom.SetNumExplicitHs(row.get("explicit_hs", 0))
        atom.SetNoImplicit(row.get("no_implicit", False))
        atom.SetChiralTag(CHIRAL[row.get("chiral_tag", "CHI_UNSPECIFIED")])
        atom.SetAtomMapNum(row.get("atom_map", 0))
        rw.AddAtom(atom)
        origin = row.get("origin") or {}
        origins.append({"instance": str(origin.get("instance", "graph")),
                        "component": str(origin.get("component", "graph")),
                        "fragment_atom_index": int(origin.get("fragment_atom_index", row["index"])),
                        "explicit_atom_index": row["index"]})
    for row in graph["bonds"]:
        rw.AddBond(row["a"], row["b"], ORDERS[row["order"]])
        rw.GetBondBetweenAtoms(row["a"], row["b"]).SetIsAromatic(row.get("aromatic", row["order"] == 1.5))
    mol = rw.GetMol()
    for row in graph["bonds"]:
        bond = mol.GetBondBetweenAtoms(row["a"], row["b"])
        if row.get("stereo_atoms"):
            bond.SetStereoAtoms(*row["stereo_atoms"])
        bond.SetStereo(STEREO[row.get("stereo", "STEREONONE")])
    mol.UpdatePropertyCache(strict=True)
    symmetric_depiction_rings(mol)
    if len(Chem.GetMolFrags(mol)) != 1:
        raise ValueError("Explicit source patch must be connected")
    conformer = Chem.Conformer(mol.GetNumAtoms())
    conformer.Set3D(False)
    ports = []
    for row, atom in zip(graph["atoms"], mol.GetAtoms()):
        conformer.SetAtomPosition(row["index"], Point3D(*row["xy"], 0))
        if atom.GetAtomicNum() == 0:
            ports.append({"atom_index": row["index"], "component_port": atom.GetAtomMapNum(), **origins[row["index"]]})
            atom.SetAtomMapNum(0)
    mol.AddConformer(conformer)
    if any(b.GetStereo() not in (Chem.BondStereo.STEREONONE, Chem.BondStereo.STEREOANY) for b in mol.GetBonds()):
        Chem.SetDoubleBondNeighborDirections(mol)
    lengths = np.asarray([math.dist(graph["atoms"][b.GetBeginAtomIdx()]["xy"], graph["atoms"][b.GetEndAtomIdx()]["xy"])
                          for b in mol.GetBonds()])
    median = float(np.median(lengths))
    serialized = json.dumps(graph, sort_keys=True, separators=(",", ":"))
    audit = {
        "schema_version": VERSION, "cofex_id": None if designed else recipe.get("cofex_id"),
        "recipe_id": recipe["recipe_id"],
        "recipe_sha256": hashlib.sha256(json.dumps(recipe, sort_keys=True, separators=(",", ":")).encode()).hexdigest(),
        "graph_scope": "finite_explicit_designed_patch_with_continuation_ports" if designed else "finite_explicit_source_diagram_with_continuation_ports",
        "whole_cof_smiles": None, "fragment_smiles": None, "fragment_smiles_canonical": False,
        "fragment_serialization_sha256": hashlib.sha256(serialized.encode()).hexdigest(),
        "fragment_serialization_format": "explicit_atom_bond_table",
        "atoms": mol.GetNumAtoms(), "bonds": mol.GetNumBonds(),
        "building_block_instances": len({o["instance"] for o in origins}),
        "connections": sum(origins[b.GetBeginAtomIdx()]["instance"] != origins[b.GetEndAtomIdx()]["instance"] for b in mol.GetBonds()),
        "continuation_ports": ports, "atom_origin": origins, "topology": None if designed else recipe.get("topology"),
        "topology_evidence": None if designed else recipe.get("topology_evidence"), "reference": None if designed else recipe.get("reference"),
        "source": None if designed else recipe.get("source"), "review_status": "user_design" if designed else recipe.get("review_status", "source_bound_provisional"),
        "ground_truth": False, "model_prediction": False, "synthetic_training_admitted": False,
        "layout_refinement": {"method": None, "applied": False, "reason": "user_manual_coordinates_preserved" if designed else "explicit_source_coordinates_preserved"},
        "depiction_geometry": {"median_bond_length": median, "minimum_bond_length": float(lengths.min()),
            "maximum_bond_length": float(lengths.max()), "max_to_median_ratio": float(lengths.max() / median),
            "bond_length_outliers_35percent": int(np.sum(abs(lengths / median - 1) > .35)),
            "interpretation": "Preserved 2D drawing coordinates, not physical geometry or strain energies."},
        "coordinate_provenance": graph.get("coordinate_provenance"),
        "source_registered": False if designed else bool(graph.get("coordinate_provenance", {}).get("source_registered", False)),
        "source_registration_history": recipe.get("source_registration_history"),
        "aromaticity_scope": "Explicit source atom/bond aromatic masks preserved; strict valence checked, no global aromaticity reassignment.",
        "stereochemistry": "Only explicit source table stereo retained; no new assignment from coordinates.",
        "double_bond_stereo_selector_remapping": {"defined_bonds": sum(bool(b.GetStereoAtoms()) for b in mol.GetBonds()), "consumed_selectors_rebound": 0},
        "defined_double_bond_geometry": _double_bond_geometry_audit(mol),
        "boundary": "Dummy atoms are outgoing continuation/crop interfaces, not terminal atoms or complete-material formulae."}
    from .cofdraw_layout import layout_diagnostics
    audit["layout_diagnostics"] = layout_diagnostics(mol, conformer.GetPositions(), origins)
    return mol, audit
