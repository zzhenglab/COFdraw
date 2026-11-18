"""Admit reviewed coordinate-only revisions, preserving chemical/source records.

Candidate recipes and their source-overlay evidence are supplied explicitly.
This does not run an optimizer, inspect source pixels or grant source registration.
"""
from __future__ import annotations

import argparse
from copy import deepcopy
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import re
import sys

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))
sys.path.insert(0, str(REPO / "scripts"))
from cof_vision.framework_drawing import build_patch
from import_manual_monomers import digest, read, write, refresh_manifest

QUALIFICATION = (
    "Same reviewed finite chemical graph; revised rigid-block planar coordinates were visually compared "
    "with the original source. This is a source-informed reconstruction, not exact atomwise registration, "
    "experimental E/Z assignment, crystal geometry or human ground truth."
)


def graph_signature(mol, audit):
    """Atom-order signature avoids canonicalizing a large aromatic macrocycle."""
    graph = {"origins": audit["atom_origin"],
             "atoms": [(a.GetAtomicNum(), a.GetIsotope(), a.GetFormalCharge(), a.GetIsAromatic(),
                        a.GetNumExplicitHs(), a.GetTotalNumHs(), int(a.GetChiralTag())) for a in mol.GetAtoms()],
             "bonds": [(b.GetBeginAtomIdx(), b.GetEndAtomIdx(), str(b.GetBondType()), b.GetIsAromatic(),
                        int(b.GetStereo()), list(b.GetStereoAtoms())) for b in mol.GetBonds()],
             "ports": audit["continuation_ports"]}
    return hashlib.sha256(json.dumps(graph, sort_keys=True, separators=(",", ":")).encode()).hexdigest()


def revised_recipe(original, candidate, evidence):
    if original["recipe_id"] != candidate.get("recipe_id") or candidate.get("designed_variant"):
        raise ValueError("Expected a candidate for the existing source recipe")
    for key in ("components", "connections", "linkage", "cofex_id", "topology"):
        if original.get(key) != candidate.get(key):
            raise ValueError("Coordinate-only admission cannot change chemistry or source identity")
    # Placement helpers may choose another port orientation to obtain the same
    # incidence. Only final keyed coordinates are transferred; original graph
    # construction and all original metadata remain authoritative.
    old_mol, old_audit = build_patch(original)
    candidate_mol, candidate_audit = build_patch(candidate)
    signature = graph_signature(old_mol, old_audit)
    if graph_signature(candidate_mol, candidate_audit) != signature:
        raise ValueError("Candidate assembled graph or atom-origin mapping changed")
    provenance = candidate.get("coordinate_provenance", {})
    if provenance.get("source_registered") is not False or evidence.get("source_registered") is not False:
        raise ValueError("This workflow cannot admit exact source registration")
    if not candidate.get("atom_coordinates") or evidence.get("graph_unchanged") is not True:
        raise ValueError("Explicit coordinates and graph-preservation evidence are required")
    result = deepcopy(original)
    result.update(atom_coordinates=deepcopy(candidate["atom_coordinates"]), layout_refinement="manual",
                  coordinate_provenance=deepcopy(provenance), source_registered=False,
                  representation_kind="source_informed_pore", preferred_style_id="source",
                  review_status="agent_reviewed_source_informed_layout")
    result["coordinate_provenance"]["qualification"] = QUALIFICATION
    result.setdefault("portable_provenance", {}).setdefault("subsequent_revisions", []).append({
        "kind": "source_informed_coordinate_revision", "previous_recipe_sha256": evidence["recipe_before_sha256"],
        "chemistry_unchanged": True, "source_registered": False, "qualification": QUALIFICATION})
    if "source_style" in result:
        result["source_style"]["source_registered"] = False
    new_mol, new_audit = build_patch(result)
    if graph_signature(new_mol, new_audit) != signature or not new_audit["defined_double_bond_geometry"]["valid"]:
        raise ValueError("Revised coordinates conflict with the retained graph or defined stereochemistry")
    return result, signature, new_audit["layout_diagnostics"]


def immutable(path, raw):
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists() and path.read_bytes() != raw:
        raise ValueError("Hash-addressed history contains conflicting bytes")
    if not path.exists():
        path.write_bytes(raw)


def admit(candidates, evidence_files, destination):
    destination = Path(destination).resolve()
    catalog_path, registry_path = destination / "catalog/current.json", destination / "registry.json"
    catalog, registry = read(catalog_path), read(registry_path)
    old_catalog_sha, old_registry_sha = digest(catalog_path), digest(registry_path)
    evidence, evidence_hashes = {}, {}
    for path in evidence_files:
        document = read(path)
        records = {r["recipe_id"]: r for r in document["accepted"]} if "accepted" in document else document
        for key, row in records.items():
            if key in evidence:
                raise ValueError("Duplicate recipe evidence")
            row = deepcopy(row)
            row.setdefault("recipe_before_sha256", row.get("original_recipe_sha256"))
            evidence[key], evidence_hashes[key] = row, digest(path)
    pending = []
    for path in candidates:
        candidate = read(path)
        rid = candidate["recipe_id"]
        rows = [row for row in registry["recipes"] if row["recipe_id"] == rid]
        if len(rows) != 1 or rid not in evidence:
            raise ValueError("Candidate needs one existing registry binding and reviewed evidence")
        entry = rows[0]
        recipe_path = (destination / entry["path"]).resolve()
        if not recipe_path.is_relative_to(destination) or digest(recipe_path) != entry["sha256"]:
            raise ValueError("Existing source recipe binding is invalid")
        record = evidence[rid]
        if record.get("recipe_before_sha256") != entry["sha256"]:
            raise ValueError("Source comparison did not use the current recipe")
        if record.get("selected_candidate_sha256") != digest(path):
            raise ValueError("Candidate differs from the reviewed selection")
        recipe, graph_sha, diagnostics = revised_recipe(read(recipe_path), candidate, record)
        raw = (json.dumps(recipe, ensure_ascii=False, indent=2) + "\n").encode("utf-8")
        if re.search(r'(?:^|["\s])[A-Za-z]:[\\/]|/Users/|/home/|private_alignment_input', raw.decode()):
            raise ValueError("Candidate contains private filesystem metadata")
        sha = hashlib.sha256(raw).hexdigest()
        pending.append((recipe_path, raw, {"recipe_id": rid, "previous_recipe_sha256": entry["sha256"],
            "recipe_sha256": sha, "candidate_sha256": digest(path), "graph_signature_sha256": graph_sha,
            "review_evidence_sha256": evidence_hashes[rid],
            "source_image_sha256": candidate["coordinate_provenance"]["review_crop_file_sha256"],
            "source_registered": False, "graph_unchanged": True, "qualification": QUALIFICATION,
            "source_overlay_evidence": record, "layout_diagnostics": diagnostics}))
        entry.update(sha256=sha, preferred_style_id="source", representation_kind="source_informed_pore")
        for bound in catalog["recipe_bindings"]:
            if bound["recipe_id"] == rid:
                bound.update(sha256=sha, portable_sha256=sha)
        metadata = catalog["source_metadata"][rid]
        metadata.update(source_recipe_sha256=sha, source_registered=False,
                        preferred_style_id="source", representation_kind="source_informed_pore",
                        coordinate_provenance=recipe["coordinate_provenance"])
        if QUALIFICATION not in metadata["qualifications"]:
            metadata["qualifications"].append(QUALIFICATION)
    if not pending:
        raise ValueError("No candidate revisions supplied")
    immutable(destination / f"catalog/history/{old_catalog_sha}.json", catalog_path.read_bytes())
    immutable(destination / f"catalog/history/{old_catalog_sha}.audit.json", (destination / "catalog/current.audit.json").read_bytes())
    immutable(destination / f"history/registries/{old_registry_sha}.json", registry_path.read_bytes())
    for path, raw, row in pending:
        immutable(destination / f"history/recipes/{row['previous_recipe_sha256']}.recipe.json", path.read_bytes())
        path.write_bytes(raw)
    write(registry_path, registry)
    catalog.update(catalog_version="cofdraw-source-layouts-7", source_registry_sha256=digest(registry_path),
                   source_layout_revision={"previous_catalog_sha256": old_catalog_sha,
                       "recipe_ids": [r[2]["recipe_id"] for r in pending], "source_registered": False})
    write(catalog_path, catalog)
    report = {"schema_version": "cofdraw-source-layout-revision-1", "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "previous_catalog_sha256": old_catalog_sha, "catalog_sha256": digest(catalog_path),
        "previous_registry_sha256": old_registry_sha, "registry_sha256": digest(registry_path),
        "monomer_count": len(catalog["monomer_presets"]), "source_count": len(registry["recipes"]),
        "monomers_unchanged": True, "reviewer_type": "agent", "ground_truth": False,
        "source_pixels_distributed": False, "revisions": [r[2] for r in pending]}
    report_path = destination / "catalog/source-layout-revisions.audit.json"
    write(report_path, report)
    audit = read(destination / "catalog/current.audit.json")
    audit.update(catalog_sha256=digest(catalog_path), registry_sha256=digest(registry_path),
                 source_layout_revisions_sha256=digest(report_path))
    write(destination / "catalog/current.audit.json", audit)
    refresh_manifest(destination)
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--candidate", action="append", type=Path, required=True)
    parser.add_argument("--evidence", action="append", type=Path, required=True)
    parser.add_argument("--destination", type=Path, default=REPO / "data/framework-drawing")
    args = parser.parse_args()
    report = admit(args.candidate, args.evidence, args.destination)
    print(json.dumps({k: v for k, v in report.items() if k != "revisions"}, indent=2))
