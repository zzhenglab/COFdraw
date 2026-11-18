"""Admit the explicitly reviewed finite monomer selection from read-only annotations.

Requires the original annotation directory and its sealed curation catalog. No
source images, PDF text, filesystem paths, or training permissions are exported.
Selection hashes pin the exact records reviewed for this drawing-library update.
"""
from __future__ import annotations

import argparse
from collections import Counter
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import re
import sys

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))
from cof_vision.cofdraw_chemistry import inspect_monomer


def read(path):
    return json.loads(Path(path).read_text(encoding="utf-8-sig"))


def write(path, value):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    Path(path).write_text(json.dumps(value, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def finite_molecule(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None or len(Chem.GetMolFrags(mol)) != 1:
        raise ValueError("Monomer must be one valid finite molecular graph")
    if any(a.GetAtomicNum() == 0 or a.GetNumRadicalElectrons() for a in mol.GetAtoms()):
        raise ValueError("Boundary placeholders and radical fragments are not admitted")
    return mol


def checked_file(path, expected):
    path = Path(path)
    if not re.fullmatch(r"[a-f0-9]{64}", str(expected)) or not path.is_file() or digest(path) != expected:
        raise ValueError("Upstream evidence file is missing or its pinned hash changed")
    return path


def review_artifacts(row):
    artifacts = []
    def inspect(value, role):
        if not isinstance(value, dict):
            return
        path = value.get("path") or value.get("manifest_path")
        sha = value.get("sha256") or value.get("manifest_sha256")
        if path and sha:
            checked_file(path, sha)
            artifacts.append({"role": role, "sha256": sha})
        else:
            for key, child in value.items():
                if isinstance(child, dict):
                    inspect(child, role + "/" + key)
    for key in ("source_review", "context_review", "paired_audit", "graph_ref", "reference_graph_binding",
                "independent_graph_review", "final_crop_review"):
        inspect(row.get(key), key)
    return artifacts


def make_row(selection, row, *, annotation_sha256, annotation_line, curation_sha256):
    if row.get("task") != "molecule" or row.get("review_status") != "agent_reviewed":
        raise ValueError("Only explicitly agent-reviewed finite molecule records are selected")
    molecules = row.get("molecules", [])
    if len(molecules) != 1:
        raise ValueError("Selected occurrence must contain exactly one finite monomer")
    source_molecule = molecules[0]
    mol = finite_molecule(source_molecule.get("smiles", ""))
    smiles = Chem.MolToSmiles(mol, isomericSmiles=True)
    if smiles != selection["canonical_smiles"]:
        raise ValueError("Selected molecular graph changed after visual review")
    source = row.get("source") or source_molecule.get("source", {})
    doi = source.get("doi") or row.get("doi")
    if not isinstance(doi, str) or not doi.startswith("10."):
        raise ValueError("Selected source lacks a DOI")
    checked_file(source["path"], source["file_sha256"])
    checked_file(source["crop_path"], selection["source_crop_sha256"])
    if source.get("crop_sha256") != selection["source_crop_sha256"]:
        raise ValueError("Annotation and reviewed crop do not agree")
    graph_ref = source_molecule.get("graph_ref") or row.get("graph_ref")
    if graph_ref:
        graph_path = checked_file(graph_ref["path"], graph_ref["sha256"])
        if graph_path.suffix.lower() == ".mol":
            graph = Chem.MolFromMolFile(str(graph_path), removeHs=True)
            if graph is None or Chem.MolToSmiles(graph, isomericSmiles=True) != smiles:
                raise ValueError("Reviewed MOL graph and selected SMILES disagree")
    functionality = inspect_monomer(smiles)
    roles = functionality["compatible_roles"]
    automatic = any(item["reactive_sites"] >= 2 for item in roles)
    if automatic != (selection["catalog_use"] == "scoped_drawing_candidate"):
        raise ValueError("Selected drawing/reference scope disagrees with runtime reactive sites")
    formula = rdMolDescriptors.CalcMolFormula(mol)
    expected_formula = source_molecule.get("expected_formula") or source_molecule.get("formula")
    if expected_formula and formula != expected_formula:
        raise ValueError("Source molecular formula and admitted graph disagree")
    isotopes = Counter((a.GetSymbol(), a.GetIsotope()) for a in mol.GetAtoms() if a.GetIsotope())
    evidence = {"annotation_id": row["annotation_id"], "annotation_file_sha256": annotation_sha256,
        "annotation_line": annotation_line, "upstream_curation_catalog_sha256": curation_sha256,
        "annotation_row_sha256": hashlib.sha256(json.dumps(row, sort_keys=True, separators=(",", ":")).encode()).hexdigest(),
        "source_doi": doi, "source_pdf_sha256": source["file_sha256"], "source_page": source["page"],
        "source_bbox": source["bbox"], "source_crop_sha256": selection["source_crop_sha256"],
        "source_label": source_molecule.get("label"), "review_status": row["review_status"],
        "upstream_admission_status": row.get("admission_status", "agent_reviewed_occurrence"),
        "review_artifacts": review_artifacts(row), "ground_truth": False,
        "visual_review": selection["visual_review"]}
    if graph_ref:
        evidence["source_graph_sha256"] = graph_ref["sha256"]
    names = list(dict.fromkeys([selection["abbreviation"], selection["full_name"],
                                source_molecule.get("name") or source_molecule.get("full_name") or selection["full_name"]]))
    result = {"id": selection["id"], "abbreviation": selection["abbreviation"], "full_name": selection["full_name"],
        "label": selection["abbreviation"] + " · " + selection["full_name"], "smiles": smiles,
        "canonical_smiles": smiles, "formula": formula,
        "isotopic_formula": rdMolDescriptors.CalcMolFormula(mol, separateIsotopes=True),
        "isotope_counts": [{"element": key[0], "isotope": key[1], "count": count} for key, count in sorted(isotopes.items())],
        "literature_names": names, "source_dois": [doi], "source_evidence": [evidence],
        "review_status": "agent_reviewed_source_monomer", "reviewer_type": "agent", "ground_truth": False,
        "reactive_groups": functionality["reactive_groups"], "compatible_roles": roles,
        "connection_count": max((r["reactive_sites"] for r in roles), default=None),
        "geometry_hint": "tetrahedral" if selection["id"] in {"tbpm_boronic", "allyl_tbpa"} else functionality["geometry_hint"],
        "node_symmetry": None, "node_symmetry_qualification": "Not assigned from functional-group counts or a 2D depiction.",
        "catalog_use": selection["catalog_use"], "automatic_prepare_supported": automatic,
        "qualification": selection["qualification"], "limitations": [selection["qualification"]],
        "compatibility_scope": "Recognized groups support explicit drawing designs only. They do not establish the source reaction, synthesis feasibility, topology or a periodic crystal.",
        "name_scope": "Names describe this exact finite graph. Source compound numbers and abbreviations remain paper-scoped."}
    if selection["id"] == "tmc":
        result["reactive_groups"].append({"kind": "acyl_chloride", "label": "Acyl chloride (reference-only)", "count": 3})
    return result


def refresh_manifest(destination):
    manifest = read(destination / "current-bundle-manifest.json")
    manifest.update(created_at_utc=datetime.now(timezone.utc).isoformat(),
                    catalog_sha256=digest(destination / manifest["catalog"]),
                    registry_sha256=digest(destination / manifest["registry"]))
    manifest["files"] = [{"path": p.relative_to(destination).as_posix(), "sha256": digest(p), "bytes": p.stat().st_size}
                         for p in sorted(destination.rglob("*")) if p.is_file()
                         and p.name != "current-bundle-manifest.json" and "exports" not in p.parts]
    write(destination / "current-bundle-manifest.json", manifest)


def import_selected(annotations, curation_path, destination):
    annotations, destination = Path(annotations).resolve(), Path(destination).resolve()
    selection_path = destination / "catalog/manual-monomer-selection.json"
    selections = read(selection_path)["selection"]
    curation = read(curation_path)
    admitted_files = {Path(row["path"]).name: row["sha256"] for row in curation["annotation_files"]}
    catalog_path = destination / "catalog/current.json"
    catalog = read(catalog_path)
    original_catalog_bytes = catalog_path.read_bytes()
    original_audit_bytes = (destination / "catalog/current.audit.json").read_bytes()
    previous_sha = digest(catalog_path)
    old_rows = catalog["monomer_presets"]
    existing_ids = {row["id"]: row for row in old_rows}
    existing_graphs = {Chem.MolToSmiles(finite_molecule(row["smiles"])): row["id"] for row in old_rows}
    additions = []
    for selected in selections:
        name = selected["annotation_file"]
        if Path(name).name != name or admitted_files.get(name) != selected["annotation_file_sha256"]:
            raise ValueError("Selected annotation is not hash-bound by the supplied curation catalog")
        path = checked_file(annotations / name, selected["annotation_file_sha256"])
        matches = [(i, json.loads(line)) for i, line in enumerate(path.read_text(encoding="utf-8-sig").splitlines(), 1)
                   if line.strip() and json.loads(line).get("annotation_id") == selected["annotation_id"]]
        if len(matches) != 1:
            raise ValueError("Selected annotation ID is absent or ambiguous")
        line, row = matches[0]
        new = make_row(selected, row, annotation_sha256=selected["annotation_file_sha256"], annotation_line=line,
                       curation_sha256=digest(curation_path))
        if new["id"] in existing_ids:
            if new != existing_ids[new["id"]]:
                raise ValueError("Existing selected ID differs; review the update explicitly")
            continue
        if new["canonical_smiles"] in existing_graphs:
            raise ValueError("Selected graph duplicates an existing monomer under another ID")
        additions.append(new)
        existing_ids[new["id"]] = new
        existing_graphs[new["canonical_smiles"]] = new["id"]
    if not additions:
        return {"added": 0, "monomer_count": len(old_rows), "unchanged": True}
    original_count = len(old_rows)
    catalog["monomer_presets"].extend(additions)
    catalog.update(catalog_version="cofdraw-manual-monomers-6", manual_monomer_update={
        "selection_sha256": digest(selection_path), "previous_catalog_sha256": previous_sha,
        "added_ids": [row["id"] for row in additions], "reviewer_type": "agent", "ground_truth": False,
        "source_drawings_unchanged": True, "source_images_distributed": False})
    text = json.dumps(catalog, ensure_ascii=False)
    if re.search(r'(?:^|["\s])[A-Za-z]:[\\/]|/Users/|/home/|source_page_text', text):
        raise ValueError("Portable catalog contains a private path or source page text")
    history = destination / "catalog/history"
    history.mkdir(exist_ok=True)
    for suffix, data in ((".json", original_catalog_bytes), (".audit.json", original_audit_bytes)):
        target = history / (previous_sha + suffix)
        if target.exists() and target.read_bytes() != data:
            raise ValueError("Existing historical catalog bytes conflict")
        if not target.exists():
            target.write_bytes(data)
    write(catalog_path, catalog)
    report = {"schema_version": "cofdraw-manual-monomer-admission-1", "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "previous_catalog_sha256": previous_sha, "catalog_sha256": digest(catalog_path),
        "upstream_curation_catalog_sha256": digest(curation_path), "selection_sha256": digest(selection_path),
        "original_monomer_count": original_count, "monomer_count": len(catalog["monomer_presets"]),
        "added_ids": [r["id"] for r in additions],
        "automatic_drawing_candidate_ids": [r["id"] for r in additions if r["automatic_prepare_supported"]],
        "reference_only_ids": [r["id"] for r in additions if not r["automatic_prepare_supported"]],
        "source_count": len(read(destination / "registry.json")["recipes"]),
        "registry_sha256": digest(destination / "registry.json"), "reviewer_type": "agent", "ground_truth": False,
        "source_pixels_distributed": False, "research_annotations_modified": False,
        "scope": "Selected complete molecular precursors, individually inspected against existing manual source/graph evidence. No automatic training admission or source reaction implementation.",
        "excluded_scope": "Model compounds, postsynthetic reagents, catalytic products, ambiguous precursor labels, boundary fragments and unreviewed records were not bulk imported."}
    write(destination / "catalog/manual-monomer-admission.audit.json", report)
    current_audit = json.loads(original_audit_bytes)
    current_audit.update(catalog_sha256=digest(catalog_path), monomer_count=report["monomer_count"],
                         manual_monomer_admission_sha256=digest(destination / "catalog/manual-monomer-admission.audit.json"))
    write(destination / "catalog/current.audit.json", current_audit)
    refresh_manifest(destination)
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--annotations", type=Path, required=True)
    parser.add_argument("--curation-catalog", type=Path, required=True)
    parser.add_argument("--destination", type=Path, default=REPO / "data/framework-drawing")
    args = parser.parse_args()
    print(json.dumps(import_selected(args.annotations, args.curation_catalog, args.destination), indent=2))


if __name__ == "__main__":
    main()
