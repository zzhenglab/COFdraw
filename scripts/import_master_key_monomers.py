"""Admit an explicit source-reviewed selection against a read-only installed master.

The master is corroboration of a finite precursor slot, not ground truth. Original
annotations, PDFs, crops, master files and installation receipts are never edited
or copied to the portable bundle. Only exact hashes and qualified identifiers are
exported. Run from a source checkout with the private research inputs available.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import re
import sqlite3

from import_manual_monomers import checked_file, digest, finite_molecule, make_row, read, refresh_manifest, write
from rdkit import Chem


def eligible_source(row):
    if row.get("ground_truth") is not False or row.get("human_verified") is True:
        raise ValueError("This admission lane requires explicitly provisional agent evidence")
    if row.get("identity_held") or row.get("chemical_identity_held") or row.get("source_name_structure_conflict"):
        raise ValueError("Unresolved identity or source conflicts cannot be admitted")
    molecules = row.get("molecules", [])
    if len(molecules) != 1 or molecules[0].get("visibility", "complete") != "complete":
        raise ValueError("A complete single source molecule is required")


def master_binding(connection, selected, snapshot):
    """Verify exact DOI, isomeric finite graph, record IDs and component IDs."""
    smiles, doi = selected["canonical_smiles"], selected["source_doi"].lower()
    matched = connection.execute(
        "SELECT r.record_id FROM records r JOIN record_monomers m "
        "ON r.record_id=m.record_id WHERE lower(r.doi)=? AND m.smiles=?", (doi, smiles)
    ).fetchall()
    record_ids = sorted(row[0] for row in matched)
    components = sorted(row[0] for row in connection.execute(
        "SELECT component_no FROM components WHERE canonical_smiles=?", (smiles,)))
    if not record_ids or record_ids != sorted(selected["master_record_ids"]):
        raise ValueError("Installed master precursor graph/DOI/record IDs changed")
    if not components or components != sorted(selected["master_component_ids"]):
        raise ValueError("Installed master component graph/IDs changed")
    return {
        "status": "source_doi_and_finite_precursor_graph_match",
        "record_ids": record_ids, "component_ids": components,
        **snapshot, "installation_receipt_sha256": selected["installation_receipt_sha256"],
        "ground_truth": False,
        "qualification": "Source-reviewed finite precursor graph agrees with the installed master at the stated DOI. Master component IDs and record IDs are reference identifiers, not independent human verification, a verified material identity, topology, or crystal structure."
    }


def extra_review_hashes(row):
    result = []
    def visit(value, role):
        if isinstance(value, list):
            for item in value:
                visit(item, role)
        elif isinstance(value, dict):
            if value.get("path") and value.get("sha256"):
                checked_file(value["path"], value["sha256"])
                result.append({"role": role, "sha256": value["sha256"]})
    for role in ("independent_peers", "reference_resolution", "qualified_source_references",
                 "current_reference_evidence", "component_reference_qualification", "admission_authorization"):
        visit(row.get(role), role)
    return result


def import_selection(research_root, destination):
    research_root, destination = Path(research_root).resolve(), Path(destination).resolve()
    pointer_path = research_root / "curation/current.json"
    continuation_path = research_root / "curation/continuation-current.json"
    pointer, continuation = read(pointer_path), read(continuation_path)
    curation_ref = pointer["catalog"]
    curation_path = checked_file(curation_ref["path"], curation_ref["sha256"])
    curation = read(curation_path)
    files = {Path(r["path"]).name: r for r in curation["annotation_files"]}
    reference = continuation["current_reference_index"]
    index = checked_file(reference["path"], reference["sha256"])
    selection_path = destination / "catalog/master-key-monomer-selection.json"
    selections = read(selection_path)["selection"]
    catalog_path = destination / "catalog/current.json"
    original_bytes = catalog_path.read_bytes()
    original_audit = (destination / "catalog/current.audit.json").read_bytes()
    previous_sha = digest(catalog_path)
    catalog = read(catalog_path)
    existing = {r["id"]: r for r in catalog["monomer_presets"]}
    graphs = {Chem.MolToSmiles(finite_molecule(r["smiles"]), isomericSmiles=True) for r in existing.values()}
    receipts = {}
    required_receipts = {sha for r in selections for sha in r["installation_receipt_sha256"]}
    for path in research_root.rglob("installation.json"):
        sha = digest(path)
        if sha in required_receipts:
            if read(path).get("status") != "installed_and_hash_verified":
                raise ValueError("A selected installation receipt is not an installed transaction")
            receipts[sha] = True
    if set(receipts) != required_receipts:
        raise ValueError("A pinned installation receipt is missing")
    additions = []
    with sqlite3.connect(index.as_uri() + "?mode=ro", uri=True) as connection:
        metadata = json.loads(connection.execute("SELECT value FROM metadata WHERE key='summary'").fetchone()[0])
        for kind in ("master", "mapped"):
            checked_file(metadata[kind + "_file"], metadata[kind + "_sha256"])
        snapshot = {"reference_index_sha256": reference["sha256"],
                    "master_sha256": metadata["master_sha256"], "mapped_sha256": metadata["mapped_sha256"]}
        for selected in selections:
            eligible_file = files.get(selected["annotation_file"])
            if not eligible_file or eligible_file["sha256"] != selected["annotation_file_sha256"]:
                raise ValueError("Selection is not bound to the current sealed annotation catalog")
            annotation_path = checked_file(eligible_file["path"], eligible_file["sha256"])
            matches = [(n, json.loads(line)) for n, line in enumerate(annotation_path.read_text("utf-8-sig").splitlines(), 1)
                       if line.strip() and json.loads(line).get("annotation_id") == selected["annotation_id"]]
            if len(matches) != 1:
                raise ValueError("Selected annotation occurrence is absent or ambiguous")
            line, row = matches[0]
            eligible_source(row)
            new = make_row(selected, row, annotation_sha256=eligible_file["sha256"],
                           annotation_line=line, curation_sha256=curation_ref["sha256"])
            new["master_binding"] = master_binding(connection, selected, snapshot)
            new["source_evidence"][0]["review_artifacts"].extend(extra_review_hashes(row))
            new["source_evidence"][0]["upstream_image_hold_reasons"] = row.get("hold_reasons", [])
            if new["id"] in existing:
                if existing[new["id"]] != new:
                    raise ValueError("An existing selected row differs; review the changed evidence")
                continue
            if new["canonical_smiles"] in graphs:
                raise ValueError("Selection duplicates an existing finite graph")
            additions.append(new)
            existing[new["id"]] = new
            graphs.add(new["canonical_smiles"])
    if not additions:
        return {"added": 0, "monomer_count": len(existing), "unchanged": True}
    previous_count = len(catalog["monomer_presets"])
    catalog["monomer_presets"].extend(additions)
    catalog["catalog_version"] = "cofdraw-master-reviewed-monomers-7"
    catalog["master_key_monomer_update"] = {"selection_sha256": digest(selection_path),
        "previous_catalog_sha256": previous_sha, "added_ids": [r["id"] for r in additions],
        "reviewer_type": "agent", "ground_truth": False, "source_drawings_unchanged": True}
    if re.search(r'(?:^|["\s])[A-Za-z]:[\\/]|/Users/|/home/|source_page_text', json.dumps(catalog)):
        raise ValueError("Portable catalog contains private paths or source page text")
    for suffix, content in ((".json", original_bytes), (".audit.json", original_audit)):
        path = destination / "catalog/history" / (previous_sha + suffix)
        if path.exists() and path.read_bytes() != content:
            raise ValueError("Historical snapshot would be overwritten")
        path.parent.mkdir(exist_ok=True)
        if not path.exists():
            path.write_bytes(content)
    # Recheck the inputs after the read phase, before publishing the local bundle.
    checked_file(index, snapshot["reference_index_sha256"])
    checked_file(curation_path, curation_ref["sha256"])
    write(catalog_path, catalog)
    report = {"schema_version": "cofdraw-master-key-monomer-admission-1",
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "previous_catalog_sha256": previous_sha, "catalog_sha256": digest(catalog_path),
        "selection_sha256": digest(selection_path), "upstream_curation_catalog_sha256": curation_ref["sha256"],
        "upstream_current_pointer_sha256": digest(pointer_path),
        "upstream_continuation_pointer_sha256": digest(continuation_path), **snapshot,
        "original_monomer_count": previous_count, "monomer_count": len(existing),
        "added_ids": [r["id"] for r in additions],
        "automatic_drawing_candidate_ids": [r["id"] for r in additions if r["automatic_prepare_supported"]],
        "reference_only_ids": [r["id"] for r in additions if not r["automatic_prepare_supported"]],
        "installation_receipt_sha256": sorted(required_receipts),
        "source_count": len(read(destination / "registry.json")["recipes"]),
        "registry_sha256": digest(destination / "registry.json"),
        "reviewer_type": "agent", "ground_truth": False,
        "research_annotations_modified": False, "source_pixels_distributed": False,
        "scope": "Selected source-reviewed finite precursors matched to the current installed master by DOI and exact isomeric molecular graph. No bulk master import, training admission, human-gold promotion, or inferred periodic identity."}
    audit_path = destination / "catalog/master-key-monomer-admission.audit.json"
    write(audit_path, report)
    audit = json.loads(original_audit)
    audit.update(catalog_sha256=report["catalog_sha256"], monomer_count=len(existing),
                 master_key_monomer_admission_sha256=digest(audit_path))
    write(destination / "catalog/current.audit.json", audit)
    refresh_manifest(destination)
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--research-root", type=Path, required=True)
    parser.add_argument("--destination", type=Path, default=Path(__file__).resolve().parents[1] / "data/framework-drawing")
    args = parser.parse_args()
    print(json.dumps(import_selection(args.research_root, args.destination), indent=2))
