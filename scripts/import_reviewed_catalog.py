"""Recover selected reviewed drawing tables; omit source pixels and private paths.

The external review archive is read-only. Run with its manual-wave3/full-paper-reviews
directory as --reviews and the repository's data/framework-drawing as --destination.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
import math
from pathlib import Path
import re
import sys

from rdkit import Chem
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from cof_vision.cofdraw_explicit import build_explicit_patch
from cof_vision.framework_drawing import build_patch


SELECTION = [
    ("acs-r327/f-cof1-source-registered-v1.audit.json", "f-cof-1-source-v1", "F-COF-1", "source_registered_pore"),
    ("root-r584/toc-network-regular/COF300_TOC.json", "cof-300-projection-v1", "COF-300 · projected network", "source_informed_projected_network"),
    ("acs-r169/BTA-COF1-local-v1.json", "bta-cof-1-local-v1", "BTA-COF1 · local boronate fragment", "source_informed_local_fragment"),
    ("acs-r169/BTA-COF2-local-v1.json", "bta-cof-2-local-v1", "BTA-COF2 · local boronate fragment", "source_informed_local_fragment"),
    ("acs-r169/BTA-COF3-local-v1.json", "bta-cof-3-local-v1", "BTA-COF3 · local boronate fragment", "source_informed_local_fragment"),
    ("acs-r169/BTA-COF4-pore-v1.json", "bta-cof-4-pore-v1", "BTA-COF4 · boroxine pore", "source_informed_pore"),
    ("acs-r187/brick-wall-pore-v1.json", "brick-wall-pore-v1", "Brick-wall COF", "source_informed_pore"),
    ("acs-r187/herringbone-hypothetical-pore-v1.json", "herringbone-proposal-v1", "Herringbone · rejected source proposal", "rejected_source_hypothesis"),
    ("acs-r381/COF-BP-parent-network-v1.json", "cof-bp-network-v1", "COF-BP · parent network", "source_informed_network"),
    ("acs-r381/COF-BP-Ni1-network-v1.json", "cof-bp-ni1-network-v1", "COF-BP-Ni1 · drawn coordination sites", "source_informed_coordination_patch"),
    ("acs-r381/COF-BP-Ni2-network-v1.json", "cof-bp-ni2-network-v1", "COF-BP-Ni2 · drawn coordination sites", "source_informed_coordination_patch"),
]


def read(path):
    return json.loads(path.read_text(encoding="utf-8"))


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")


def portable_source(record):
    source = record.get("source_pair", {}).get("source", record.get("source", {}))
    allowed = ("doi", "file_sha256", "sha256", "filename", "page", "page_count", "figure",
               "bbox", "bbox_points", "document_role", "crop_sha256", "dpi", "width", "height")
    result = {k: source[k] for k in allowed if k in source}
    result["doi"] = record.get("doi") or source.get("doi")
    result["source_pixels_distributed"] = False
    return result


def verify_review(path, record):
    """Verify peer-bound draft and review hashes without distributing those files."""
    evidence = [{"artifact": path.name, "sha256": digest(path), "role": "reviewed_graph_table"}]
    refs = []
    for key in ("original_peer_bound_draft", "independent_bounded_review", "full_paper_review"):
        if record.get(key): refs.append(record[key])
    refs.extend(record.get("independent_bounded_reviews", []))
    for ref in refs:
        if not isinstance(ref, dict) or not ref.get("path") or not ref.get("sha256"):
            continue
        original = Path(ref["path"])
        # The archived basename is sufficient when moved to another workstation.
        source = path.parent / original.name
        if not source.is_file():
            source = original
        if not source.is_file() or digest(source) != ref["sha256"]:
            raise ValueError(f"Upstream review hash mismatch for {path.name}: {original.name}")
        evidence.append({"artifact": source.name, "sha256": ref["sha256"], "role": "upstream_review_or_bound_draft"})
    if path.name == "COF300_TOC.json":
        peer = path.with_name("independent-web-toc-network-review-0001.json")
        reviewed = read(peer)
        if reviewed.get("status") != "accepted_provisionally_with_scope_limits" or digest(path) not in reviewed.get("files", {}).values():
            raise ValueError("COF-300 review does not bind this graph revision")
        evidence.append({"artifact": peer.name, "sha256": digest(peer), "role": "independent_bounded_review"})
    elif not ("reviewed" in record.get("review_status", "") or "reviewed" in record.get("visual_audit_status", "")):
        raise ValueError(f"Unreviewed candidate: {path.name}")
    return evidence


def convert_graph(record):
    raw_atoms, raw_bonds = record["atom_table"], record["bond_table"]
    anchor_map = {a["atom_index"]: a for a in record.get("source_coordinate_anchors", [])}
    xy = [a.get("xy", anchor_map.get(a["index"], {}).get("xy_pixels")) for a in raw_atoms]
    if not all(p is not None for p in xy):
        raise ValueError("An admitted drawing must provide coordinates for every atom")
    orders = {"SINGLE": 1, "DOUBLE": 2, "TRIPLE": 3, "AROMATIC": 1.5, "DATIVE": "DATIVE"}
    edges = [(b.get("a", b.get("begin")), b.get("b", b.get("end"))) for b in raw_bonds]
    lengths = sorted(math.dist(xy[a], xy[b]) for a, b in edges)
    scale = 1.5 / lengths[len(lengths) // 2]
    center = [sum(p[i] for p in xy) / len(xy) for i in range(2)]
    atoms = []
    for a, point in zip(raw_atoms, xy):
        origin = a.get("origin") or {"instance": a.get("component", "graph"),
                                     "component": re.sub(r"\d+$", "", a.get("component", "graph")),
                                     "fragment_atom_index": a.get("parent_atom_index", a["index"])}
        row = {"index": a["index"], "element": a["element"], "formal_charge": a.get("formal_charge", 0),
               "aromatic": a.get("aromatic", False), "atom_map": a.get("atom_map", 0),
               "isotope": a.get("isotope", 0), "chiral_tag": a.get("chiral_tag", "CHI_UNSPECIFIED"),
               "origin": origin, "xy": [(point[0] - center[0]) * scale, -(point[1] - center[1]) * scale]}
        if "hydrogens" in a:
            row.update(explicit_hs=a["hydrogens"], no_implicit=True)
        atoms.append(row)
    stereo_map = {tuple(sorted(b["bond"])): b for b in record.get("stereo_geometry", [])}
    bonds = []
    for b, (a, c) in zip(raw_bonds, edges):
        order = orders.get(b["order"], b["order"])
        selectors = b.get("stereo_atoms", [])
        if not selectors and tuple(sorted((a, c))) in stereo_map:
            s = stereo_map[tuple(sorted((a, c)))]; selectors = s["selectors"]
            if s["bond"] != [a, c]: selectors = list(reversed(selectors))
        bonds.append({"a": a, "b": c, "order": order, "aromatic": b.get("aromatic", order == 1.5),
                      "stereo": b.get("stereo", "STEREONONE"), "stereo_atoms": selectors})
    return {"atoms": atoms, "bonds": bonds, "coordinate_system": "cartesian_y_up",
            "source_coordinates_preserved": True,
            "coordinate_provenance": {"upstream_coordinate_system": "drawing_pixels_y_down",
                "source_registered": bool(anchor_map), "uniform_scale": scale, "translation_before_scale": center,
                "y_axis_reversed": True, "geometry_optimization_performed": False,
                "qualification": "All upstream atom coordinates retained up to a uniform affine viewport conversion; generic reconstructions remain generic."}}


def reviewed_monomers(reviews, folder, ids):
    doc = read(reviews / folder / "small-graphs-v1.json")
    if "reviewed" not in doc.get("review_status", ""):
        raise ValueError("Monomer list has not been reviewed")
    return [{"label": row["id"], "full_name": row["label"], "smiles": row["smiles"],
             "upstream_monomer_review_sha256": digest(reviews / folder / "small-graphs-v1.json")}
            for row in doc["data"] if row["id"] in ids]


def monomers_for(reviews, spec, record, existing):
    path, recipe_id, _, _ = spec
    if recipe_id == "f-cof-1-source-v1":
        return [{"label": key, "full_name": "1,3,5-Tris(4-formylphenyl)benzene" if key == "FTDA" else "2,4,6-Trimethyl-1,3,5-triazine",
                 "smiles": sm} for key, sm in record["parent_graphs"].items()]
    if recipe_id == "cof-300-projection-v1":
        return [{"label": m["id"], "full_name": m["full_name"], "smiles": m["smiles"]} for m in existing if m["id"] in {"tapm", "ta"}]
    if path.startswith("acs-r169"):
        partner = {"bta-cof-1-local-v1": "HHTP", "bta-cof-2-local-v1": "THDMA", "bta-cof-3-local-v1": "THDEA"}.get(recipe_id)
        return reviewed_monomers(reviews, "acs-r169", {"compound2", partner})
    if path.startswith("acs-r187"):
        return reviewed_monomers(reviews, "acs-r187", {"T-node", "PDA"})
    return reviewed_monomers(reviews, "acs-r381", {"TFBP", "BD"})


def source_style(reviews, spec, record, graph):
    path, recipe_id, _, kind = spec
    style = {"preferred_style_id": "source", "source_registered": kind == "source_registered_pore",
             "exact_style_equivalence": False, "atom_colors": {}, "component_colors": {}, "component_fills": {},
             "observed": record.get("style", {}), "qualification": "Reconstructed vector display style; exact source fonts are not known."}
    folder = path.split("/")[0]
    observed = reviews / folder / "source-style-v1.json"
    if observed.exists():
        style["observed"] = read(observed)
        style["evidence_sha256"] = digest(observed)
    if recipe_id == "f-cof-1-source-v1":
        style["atom_colors"] = {str(a["atom_index"]): "#171717" if a["color"] == "black" else a["color"] for a in record["source_coordinate_anchors"]}
        style["halo"] = {"color": "#fff17a", "opacity": .55, "radius_bonds": .22,
                         "scope": "Decorative fluorescence halo, not chemical density or orbital data."}
        style["edge_style"] = "wavy_cross"
    elif recipe_id == "cof-300-projection-v1":
        style["atom_colors"] = {str(a["index"]): "#c82828" for a in graph["atoms"] if a["element"] == "N"}
        style.update(faded_atoms=record.get("faded_atoms", []), fade_opacity=.28, edge_style="crossbar")
    elif path.startswith("acs-r187"):
        style["component_colors"] = {"T-node": "#171717", "PDA": "#2768c6"}
        style["atom_colors"] = {str(a["index"]): "#2768c6" for a in graph["atoms"] if a["element"] == "N" and a["origin"]["component"] == "PDA"}
        style["edge_style"] = "wavy_cross"
    elif path.startswith("acs-r381"):
        style["component_fills"] = {"TFBP": "#f2a3c1", "BD": "#9fdfdf", "salicylideneaniline": "#f6ce91"}
        style["atom_colors"] = {str(a["index"]): {"N": "#2222cc", "O": "#e34c32", "S": "#e78b21", "Ni": "#26843b"}[a["element"]]
                                for a in graph["atoms"] if a["element"] in {"N", "O", "S", "Ni"}}
        style["edge_style"] = "dashed"
    return style


def refresh_manifest(destination):
    catalog = destination / "catalog/current.json"
    registry = destination / "registry.json"
    manifest = {"schema_version": "cofdraw-current-portable-bundle-1",
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "catalog": "catalog/current.json", "catalog_sha256": digest(catalog),
        "registry": "registry.json", "registry_sha256": digest(registry),
        "files": [{"path": p.relative_to(destination).as_posix(), "sha256": digest(p), "bytes": p.stat().st_size}
                  for p in sorted(destination.rglob("*")) if p.is_file()
                  and p.name != "current-bundle-manifest.json" and "exports" not in p.parts],
        "source_pdf_or_image_included": False, "private_path_scan_passed": True, "historical_files_retained": True}
    write(destination / "current-bundle-manifest.json", manifest)


def restore_known_registrations(reviews, destination, catalog, registry):
    """Retain previously reviewed source anchors omitted from the desktop bundle."""
    base = reviews.parent
    selections = [
        ("tfb-tab-pore-v1", "framework-redraw-pilot/source-registration/tfb-tab-source-anchors-v1.json",
         "framework-redraw-pilot/source-registration/tfb-tab-source-registered-v1.audit.json"),
        ("wiley-r198-tppa-so3h-pore-0002", "periodic-fragments/wiley-r198-source-depiction-0002.json",
         "periodic-fragments/wiley-r198-source-depiction-render-v4-visual-review.json"),
    ]
    restored = []
    for rid, relative, review_relative in selections:
        anchor_path, review_path = base / relative, base / review_relative
        annotation, review = read(anchor_path), read(review_path)
        anchor_hash = digest(anchor_path)
        bound = (review.get("source_anchor_sha256") == anchor_hash
                 or any(r.get("sha256") == anchor_hash for r in review.get("inputs", [])))
        if not bound:
            raise ValueError(f"Registration review does not bind {rid}")
        entry = next(r for r in registry["recipes"] if r["recipe_id"] == rid)
        target = destination / entry["path"]
        recipe = read(target)
        recipe.pop("atom_coordinates", None)
        mol, audit = build_patch(recipe)
        origins = audit["atom_origin"]
        keys = [(o["instance"], o["fragment_atom_index"]) for o in origins]
        index = {f"{a}:{b}": i for i, (a, b) in enumerate(keys)}
        anchors = {(r["instance"], r["fragment_atom_index"]): np.asarray(r["xy"], dtype=float) for r in annotation["atom_coordinates"]}
        if not set(anchors).issubset(set(keys)):
            raise ValueError(f"Source anchor keys changed: {rid}")
        pixel_positions = np.asarray(mol.GetConformer().GetPositions())[:, :2].copy()
        pixel_positions[:, 1] *= -1
        placed = np.zeros_like(pixel_positions)
        hidden = []
        for instance in sorted({o["instance"] for o in origins}):
            members = [i for i, o in enumerate(origins) if o["instance"] == instance]
            observed = [i for i in members if keys[i] in anchors]
            if len(observed) < 3:
                raise ValueError(f"Insufficient source anchors for {rid}/{instance}")
            x = pixel_positions[observed]
            y = np.asarray([anchors[keys[i]] for i in observed])
            xc, yc = x - x.mean(0), y - y.mean(0)
            u, s, vt = np.linalg.svd(xc.T @ yc)
            rotation = u @ vt
            local_scale = float(s.sum() / np.sum(xc * xc))
            for i in members:
                if keys[i] in anchors:
                    placed[i] = anchors[keys[i]]
                else:
                    # Source-condensed SO3H oxygen positions are unobserved.
                    # Keep chemical export atoms near their fitted parent group.
                    placed[i] = (pixel_positions[i] - x.mean(0)) @ rotation * local_scale + y.mean(0)
                    hidden.append(i)
        observed_lengths = [math.dist(placed[b.GetBeginAtomIdx()], placed[b.GetEndAtomIdx()]) for b in mol.GetBonds()
                            if keys[b.GetBeginAtomIdx()] in anchors and keys[b.GetEndAtomIdx()] in anchors]
        viewport_scale = 1.5 / float(np.median(observed_lengths))
        center = placed.mean(0)
        normalized = (placed - center) * viewport_scale
        normalized[:, 1] *= -1
        recipe["atom_coordinates"] = [{"instance": o["instance"], "fragment_atom_index": o["fragment_atom_index"], "xy": list(normalized[i])}
                                      for i, o in enumerate(origins)]
        observed_style = annotation["style"]
        atom_colors = {str(index[k]): color for k, color in observed_style.get("atom_colors", {}).items() if k in index}
        for group in observed_style.get("collapsed_groups", []):
            if group["anchor"] in index and group.get("color"):
                atom_colors[str(index[group["anchor"]])] = group["color"]
        source_style = {"preferred_style_id": "source", "source_registered": True, "exact_style_equivalence": False,
            "atom_colors": atom_colors, "component_colors": observed_style.get("component_colors", {}),
            "edge_style": observed_style.get("port_style", "dashed"),
            "abbreviation_mode": "condensed" if hidden else "explicit", "font_family": observed_style.get("font_family", "Arial"),
            "qualification": "Previously reviewed manual source anchors. Fonts and stroke details remain reconstructed substitutes.",
            "coordinate_evidence_sha256": anchor_hash, "review_evidence_sha256": digest(review_path)}
        recipe["source_style"] = source_style
        recipe["source_coordinate_provenance"] = {"source_registered": True, "visible_atom_count": len(anchors),
            "unobserved_abbreviated_atom_count": len(hidden), "unobserved_atom_indices": hidden,
            "source_anchor_sha256": anchor_hash, "review_sha256": digest(review_path),
            "uniform_viewport_scale": viewport_scale, "pixel_center": list(center),
            "y_axis_reversed": True, "source_anchor_positions_changed": False,
            "unobserved_coordinate_method": "Rigid similarity fit of the parent component; illustrative chemical export positions only." if hidden else None}
        checked_mol, checked_audit = build_patch(recipe)
        if checked_mol.GetNumAtoms() != mol.GetNumAtoms() or not checked_audit["defined_double_bond_geometry"]["valid"]:
            raise ValueError(f"Registration changed graph or invalidated stereo: {rid}")
        write(target, recipe)
        bound_hash = digest(target)
        entry.update(sha256=bound_hash, preferred_style_id="source", representation_kind="source_registered_pore")
        metadata = catalog["source_metadata"][rid]
        metadata.update(source_recipe_sha256=bound_hash, preferred_style_id="source", source_registered=True,
                        coordinate_provenance=recipe["source_coordinate_provenance"])
        for binding in catalog.get("recipe_bindings", []):
            if binding["recipe_id"] == rid:
                binding.update(sha256=bound_hash, portable_sha256=bound_hash)
        restored.append({"recipe_id": rid, "visible_atoms": len(anchors), "unobserved_atoms": len(hidden),
                         "source_anchor_sha256": anchor_hash, "recipe_sha256": bound_hash})
    return restored


def graph_signature(atoms, bonds):
    """Coordinate-independent, explicit graph identity for reviewed revisions."""
    orders = {"SINGLE": 1, "DOUBLE": 2, "TRIPLE": 3, "AROMATIC": 1.5, "DATIVE": "DATIVE"}
    atom_rows = [(a["index"], a["element"], a.get("formal_charge", a.get("charge", 0)),
                  a.get("aromatic", False), a.get("isotope", 0), a.get("atom_map", 0),
                  a.get("chiral_tag", "CHI_UNSPECIFIED"),
                  a.get("explicit_hs", a.get("hydrogens", a.get("total_hydrogens", 0)))) for a in atoms]
    bond_rows = [(b.get("a", b.get("begin")), b.get("b", b.get("end")),
                  orders.get(b["order"], b["order"]), b.get("aromatic", False),
                  b.get("stereo", "STEREONONE"), tuple(b.get("stereo_atoms", []))) for b in bonds]
    return atom_rows, bond_rows


def admit_layout_revisions(reviews, revision_dir, dest, catalog, registry):
    """Admit source-informed coordinates only when reviewed chemistry is unchanged.

    The supplied revision must carry its manual placement evidence. The public
    recipe is itself the portable coordinate revision; the private source crop
    and the intermediate recipe are identified by hashes, never copied.
    """
    specs = [
        ("pyurea", "rsc-r049", "Pyurea-pore.json", "Pyurea-COF",
         "independent-rsc-source-incidence-review-0001.json", {"Py", "Urea"}),
        ("ucof", "rsc-r055", "UCOF-pore-graph.json", "UCOF",
         "independent-acs-urethane-and-finite-review-0001.json", {"M01", "M02"}),
    ]
    admitted, new_monomers = [], []
    for short, folder, graph_name, label, peer_name, monomer_ids in specs:
        upstream = reviews / folder / "chemistry-drafts-v1" / graph_name
        record = read(upstream)
        revision_path = revision_dir / f"{short}-source-layout-adjusted.recipe.json"
        recipe = read(revision_path)
        graph = recipe["explicit_graph"]
        if recipe["reference"]["upstream_graph_sha256"] != digest(upstream):
            raise ValueError(f"Layout revision does not bind upstream graph: {short}")
        if graph_signature(record["atom_table"], record["bond_table"]) != graph_signature(graph["atoms"], graph["bonds"]):
            raise ValueError(f"Layout revision changed explicit chemistry: {short}")
        peer_path = reviews / folder / peer_name
        peer = read(peer_path)
        if not peer.get("decision", peer.get("status", "")).startswith("accepted"):
            raise ValueError(f"Bounded chemical review not accepted: {short}")
        monomer_path = upstream.with_name("finite-graphs-draft.json")
        bound_refs = list(peer.get("graph_inputs", []))
        if short == "ucof":
            bound_refs += [peer["finite_manifest"], peer["pore_manifest"]]
        for ref in bound_refs:
            p = upstream.parent / Path(ref["path"]).name
            if digest(p) != ref["sha256"]:
                raise ValueError(f"Bounded review artifact changed: {p.name}")
        if short == "ucof":
            manifest = read(upstream.with_name("connected-graph-manifest-draft.json"))
            if digest(upstream) not in [r["graph"]["sha256"] for r in manifest["records"]]:
                raise ValueError("UCOF reviewed manifest does not bind graph")
        elif digest(upstream) not in [r["sha256"] for r in bound_refs]:
            raise ValueError("Pyurea review does not bind graph")
        layout = recipe["layout_evidence"]
        if layout.get("exact_atom_registration") is not False or not layout.get("block_anchor_records"):
            raise ValueError("Manual source-informed revision requires qualified placement evidence")
        if layout["source_image_sha256"] != recipe["source"]["crop_sha256"]:
            raise ValueError("Layout and graph source crop hashes differ")
        evidence = [{"artifact": p.name, "sha256": digest(p), "role": role} for p, role in [
            (upstream, "reviewed_graph_table"), (peer_path, "independent_bounded_chemical_review"),
            (monomer_path, "reviewed_monomer_table"), (revision_path, "source_informed_layout_revision")]]
        monomers = [{"label": m["label"], "full_name": m["name"], "smiles": m["smiles"],
                     "upstream_monomer_review_sha256": digest(monomer_path)}
                    for m in read(monomer_path) if m["id"] in monomer_ids]
        recipe.pop("private_alignment_input", None)
        rid = recipe["recipe_id"]
        kind = "source_informed_manual_pore"
        scope = ("Finite source-visible chemical pore with manually source-informed block placement; "
                 "not exact atom registration or a periodic crystal model. "
                 + ("Eight directly para-aryl-bound boundary cuts." if short == "pyurea" else
                    "Organic urethane-linked graph only; source photograph and measurement overlays excluded."))
        recipe.update(review_status="agent_reviewed_source_informed_layout", source_representation_kind=kind,
                      source_scope=scope,
                      components={m["label"]: {"label": m["label"], "monomer_smiles": m["smiles"]} for m in monomers})
        recipe["reference"]["review_evidence"] = evidence
        graph["source_coordinates_preserved"] = False
        graph["coordinate_provenance"] = {
            "source_registered": False, "upstream_generic_coordinates_preserved": False,
            "saved_coordinate_revision_preserved_on_render": True,
            "method": "Manual source-informed rigid block placement and junction-spacing correction",
            "physical_geometry_optimization_performed": False,
            "layout_evidence_field": "layout_evidence",
            "qualification": "Source ring centers guide this new coordinate revision; atomwise pixel correspondence is not claimed."}
        recipe["source_style"].update(exact_style_equivalence=False,
            qualification="Source-informed colors, block arrangement and boundary notation; fonts and exact atom placement are reconstructed.")
        mol, audit = build_explicit_patch(recipe)
        if not audit["defined_double_bond_geometry"]["valid"]:
            raise ValueError(f"Defined stereo contradicted by source-informed placement: {short}")
        target = dest / "recovered" / f"{rid}.recipe.json"
        write(target, recipe)
        sha = digest(target)
        src = recipe["source"]
        metadata = {"literature_name": label, "display_name": label, "source_doi": src["doi"],
            "source_dois": [src["doi"]], "monomers": monomers, "cofex_id": None, "cofex_notation": None,
            "identity_qualification": "Source-qualified finite drawing; no complete periodic identity asserted.",
            "representation_kind": kind, "preferred_style_id": "source", "ground_truth": False,
            "source_registered": False, "review_status": recipe["review_status"], "scope": scope,
            "review_evidence": evidence, "source_recipe_sha256": sha}
        catalog["source_metadata"][rid] = metadata
        registry["recipes"] = [r for r in registry["recipes"] if r["recipe_id"] != rid]
        registry["recipes"].append({"recipe_id": rid, "cofex_id": None, "label": label,
            "path": target.relative_to(dest).as_posix(), "sha256": sha, "upstream_sha256": digest(upstream),
            "preferred_style_id": "source", "representation_kind": kind})
        catalog["recipe_bindings"] = [r for r in catalog["recipe_bindings"] if r["recipe_id"] != rid]
        catalog["recipe_bindings"].append({"recipe_id": rid, "filename": target.relative_to(dest).as_posix(),
            "sha256": sha, "portable_sha256": sha, "upstream_sha256": digest(upstream)})
        style_id = "recovered-" + rid
        catalog["source_styles"] = [s for s in catalog["source_styles"] if s["id"] != style_id]
        catalog["source_styles"].append({"id": style_id, "source_dois": [src["doi"]], "source_recipe_ids": [rid],
            "observed": {}, "suggested_runtime_styles": ["source"], "exact_style_equivalence": False,
            "source_registered": False, "qualification": scope})
        admitted.append({"recipe_id": rid, "label": label, "doi": src["doi"], "atoms": mol.GetNumAtoms(),
            "bonds": mol.GetNumBonds(), "ports": len(audit["continuation_ports"]), "representation_kind": kind,
            "upstream_sha256": digest(upstream), "recipe_sha256": sha,
            "stereo": audit["defined_double_bond_geometry"]})
        for m in monomers:
            key = Chem.MolToSmiles(Chem.MolFromSmiles(m["smiles"]))
            if any(Chem.MolToSmiles(Chem.MolFromSmiles(old["smiles"])) == key for old in catalog["monomer_presets"]):
                continue
            mid = re.sub(r"[^a-z0-9]+", "_", "recovered_" + m["label"].lower()).strip("_")
            catalog["monomer_presets"].append({"id": mid, "label": m["label"] + " · " + m["full_name"],
                "abbreviation": m["label"], "full_name": m["full_name"], "smiles": m["smiles"],
                "source_doi": src["doi"], "source_dois": [src["doi"]], "review_status": "agent_reviewed_source_monomer",
                "ground_truth": False, "qualification": "Source-derived finite monomer; not an automatic linkage support claim."})
            new_monomers.append(mid)
    return admitted, new_monomers


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reviews", type=Path, required=True)
    parser.add_argument("--destination", type=Path, required=True)
    parser.add_argument("--layout-revisions", type=Path,
                        help="Optional reviewed Pyurea/UCOF manual placement revision directory")
    args = parser.parse_args()
    reviews, dest = args.reviews.resolve(), args.destination.resolve()
    catalog_path = dest / "catalog/current.json"
    catalog, registry = read(catalog_path), read(dest / "registry.json")
    restored = restore_known_registrations(reviews, dest, catalog, registry)
    admitted, new_monomers = [], []
    for spec in SELECTION:
        relative, recipe_id, label, kind = spec
        upstream = reviews / relative
        record = read(upstream)
        evidence = verify_review(upstream, record)
        graph = convert_graph(record)
        monomers = monomers_for(reviews, spec, record, catalog["monomer_presets"])
        src = portable_source(record)
        if not src.get("doi") and relative.startswith("root-r584"):
            src["doi"] = "10.1021/acsami.3c12826"
        metadata = {"literature_name": record.get("material_name", label), "display_name": label,
            "source_doi": src["doi"], "source_dois": [src["doi"]], "monomers": monomers,
            "cofex_id": record.get("cofex_reference", {}).get("key"), "cofex_notation": None,
            "identity_qualification": "Source-qualified drawing record. No complete periodic identity is asserted.",
            "representation_kind": kind, "preferred_style_id": "source", "ground_truth": False,
            "source_registered": kind == "source_registered_pore", "review_status": "agent_reviewed_source_drawing",
            "scope": record.get("scope", record.get("material_scope", record.get("site_scope", "Finite source-informed diagram with explicit boundary ports."))),
            "review_evidence": evidence}
        if kind == "rejected_source_hypothesis":
            metadata["scope"] = "Rejected herringbone proposal in the source paper; not a synthesized material claim."
        recipe = {"schema_version": "cof-chemical-patch-1", "recipe_id": recipe_id, "name": label,
            "cofex_id": metadata["cofex_id"], "ground_truth": False, "review_status": "agent_reviewed_source_drawing",
            "source": src, "reference": {"upstream_graph_sha256": digest(upstream), "review_evidence": evidence},
            "explicit_graph": graph, "components": {m["label"]: {"label": m["label"], "monomer_smiles": m["smiles"]} for m in monomers},
            "instances": [], "connections": [], "topology": None, "source_representation_kind": kind,
            "source_scope": metadata["scope"], "source_style": source_style(reviews, spec, record, graph)}
        mol, audit = build_explicit_patch(recipe)
        if not audit["defined_double_bond_geometry"]["valid"]:
            raise ValueError(f"Contradictory defined stereo in {recipe_id}: {audit['defined_double_bond_geometry']}")
        if mol.GetNumAtoms() != len(record["atom_table"]) or mol.GetNumBonds() != len(record["bond_table"]):
            raise ValueError(f"Graph size changed during import: {recipe_id}")
        target = dest / "recovered" / f"{recipe_id}.recipe.json"
        write(target, recipe)
        bound_hash = digest(target)
        metadata["source_recipe_sha256"] = bound_hash
        catalog["source_metadata"][recipe_id] = metadata
        registry["recipes"] = [r for r in registry["recipes"] if r["recipe_id"] != recipe_id]
        registry["recipes"].append({"recipe_id": recipe_id, "cofex_id": metadata["cofex_id"], "label": label,
            "path": target.relative_to(dest).as_posix(), "sha256": bound_hash, "upstream_sha256": digest(upstream),
            "preferred_style_id": "source", "representation_kind": kind})
        if isinstance(catalog.get("recipe_bindings"), dict):
            catalog["recipe_bindings"][recipe_id] = bound_hash
        else:
            catalog["recipe_bindings"] = [r for r in catalog.get("recipe_bindings", []) if r["recipe_id"] != recipe_id]
            catalog["recipe_bindings"].append({"recipe_id": recipe_id, "filename": target.relative_to(dest).as_posix(),
                "sha256": bound_hash, "portable_sha256": bound_hash, "upstream_sha256": digest(upstream)})
        style_id = "recovered-" + recipe_id
        catalog["source_styles"] = [s for s in catalog.get("source_styles", []) if s["id"] != style_id]
        catalog["source_styles"].append({"id": style_id, "source_dois": [src["doi"]],
            "source_recipe_ids": [recipe_id], "observed": recipe["source_style"].get("observed", {}),
            "suggested_runtime_styles": ["source"], "exact_style_equivalence": False,
            "source_registered": kind == "source_registered_pore", "qualification": metadata["scope"]})
        admitted.append({"recipe_id": recipe_id, "label": label, "doi": src["doi"], "atoms": mol.GetNumAtoms(),
            "bonds": mol.GetNumBonds(), "ports": len(audit["continuation_ports"]), "representation_kind": kind,
            "upstream_sha256": digest(upstream), "recipe_sha256": bound_hash,
            "stereo": audit["defined_double_bond_geometry"]})
        for m in monomers:
            key = Chem.MolToSmiles(Chem.MolFromSmiles(m["smiles"]))
            if any(Chem.MolToSmiles(Chem.MolFromSmiles(old["smiles"])) == key for old in catalog["monomer_presets"]):
                continue
            mid = re.sub(r"[^a-z0-9]+", "_", "recovered_" + m["label"].lower()).strip("_")
            if any(old["id"] == mid for old in catalog["monomer_presets"]):
                mid += "_" + hashlib.sha256(key.encode()).hexdigest()[:6]
            row = {"id": mid, "label": m["label"] + " · " + m["full_name"], "abbreviation": m["label"],
                "full_name": m["full_name"], "smiles": m["smiles"], "source_doi": src["doi"], "source_dois": [src["doi"]],
                "review_status": "agent_reviewed_source_monomer", "ground_truth": False,
                "qualification": "Source-derived finite monomer. Presence in the library does not imply an implemented automatic linkage transform."}
            catalog["monomer_presets"].append(row); new_monomers.append(mid)
    if args.layout_revisions:
        revisions, revised_monomers = admit_layout_revisions(reviews, args.layout_revisions.resolve(), dest, catalog, registry)
        admitted.extend(revisions)
        new_monomers.extend(revised_monomers)
    write(dest / "registry.json", registry)
    catalog.update(catalog_version="cofdraw-reviewed-drawings-5", source_registry_sha256=digest(dest / "registry.json"),
                   recovered_drawing_scope="Reviewed explicit graph drawings with retained coordinate-revision and source-style qualifications; no source pixels distributed.")
    write(catalog_path, catalog)
    report = {"schema_version": "cofdraw-reviewed-catalog-import-1", "imported_at_utc": datetime.now(timezone.utc).isoformat(),
              "entries": admitted, "new_monomer_ids": new_monomers, "original_source_count": 4,
              "source_count": len(registry["recipes"]), "monomer_count": len(catalog["monomer_presets"]),
              "restored_existing_source_registrations": restored,
              "human_verification_claimed": False, "source_images_distributed": False,
              "limits": ["Most recovered coordinates are organized source-informed reconstructions, not pixel registrations.",
                         "F-COF-1 and the restored existing TFB-TAB/SCOF entries retain source atomwise anchors; Pyurea/UCOF are manually source-informed revisions.",
                         "Herringbone is a rejected proposal. COF-300 is a projection. Ni site drawings do not establish bulk occupancy."]}
    write(dest / "catalog/recovered-drawings.audit.json", report)
    write(dest / "catalog/current.audit.json", {"schema_version": "cofdraw-current-catalog-audit-2",
        "catalog_sha256": digest(catalog_path), "registry_sha256": digest(dest / "registry.json"),
        "source_count": report["source_count"], "monomer_count": report["monomer_count"],
        "recovered_drawing_audit_sha256": digest(dest / "catalog/recovered-drawings.audit.json"),
        "reviewer_type": "agent", "ground_truth": False, "source_pixels_distributed": False})
    # The curated output must contain neither workstation paths nor source text.
    for p in (dest / "recovered").glob("*.json"):
        text = p.read_text(encoding="utf-8")
        if re.search(r"[A-Za-z]:[\\/]", text) or "page_text" in text:
            raise ValueError(f"Private path or publication text found in {p.name}")
    refresh_manifest(dest)
    print(json.dumps({"sources": report["source_count"], "monomers": report["monomer_count"], "added": [r["recipe_id"] for r in admitted]}, ensure_ascii=True))


if __name__ == "__main__":
    main()
