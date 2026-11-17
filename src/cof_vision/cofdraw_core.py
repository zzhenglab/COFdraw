"""Portable, offline COFdraw monomer editing and presentation styles.

Presentation changes never change the chemical graph or establish a COF identity.
"""
from __future__ import annotations

from copy import deepcopy
import hashlib
import json
import math
from pathlib import Path
import random

import numpy as np
from rdkit import Chem
from rdkit.Chem import rdDepictor, rdMolDescriptors
from rdkit.Chem.Draw import rdMolDraw2D

from .cofdraw_styles import STYLES, PALETTES, style_configuration, decorate_svg
MONOMERS = [
    {"id": "tfb", "label": "TFB · benzene-1,3,5-tricarbaldehyde", "smiles": "O=Cc1cc(C=O)cc(C=O)c1", "role": "aldehyde"},
    {"id": "tp", "label": "Tp · triformylphloroglucinol", "smiles": "O=Cc1c(O)c(C=O)c(O)c(C=O)c1O", "role": "aldehyde"},
    {"id": "pa", "label": "Pa · p-phenylenediamine", "smiles": "Nc1ccc(N)cc1", "role": "amine"},
    {"id": "tapb", "label": "TAPB · tris(4-aminophenyl)benzene", "smiles": "Nc1ccc(-c2cc(-c3ccc(N)cc3)cc(-c3ccc(N)cc3)c2)cc1", "role": "amine"},
    {"id": "dmtp", "label": "DMTP · dimethoxyterephthalaldehyde", "smiles": "COc1cc(C=O)c(OC)cc1C=O", "role": "aldehyde"},
    {"id": "pa-so3h", "label": "Diaminobenzenesulfonic acid", "smiles": "Nc1ccc(N)cc1S(=O)(=O)O", "role": "amine"},
]
PRESETS = [
    {"id": "tfb-tapb", "label": "TFB + TAPB · imine", "linkage": "imine", "monomer_a": MONOMERS[0]["smiles"], "monomer_b": MONOMERS[3]["smiles"], "label_a": "TFB", "label_b": "TAPB", "instances": 6},
    {"id": "tp-pa", "label": "Tp + Pa · keto-enamine", "linkage": "beta_ketoenamine", "monomer_a": MONOMERS[1]["smiles"], "monomer_b": MONOMERS[2]["smiles"], "label_a": "Tp", "label_b": "Pa", "instances": 12},
    {"id": "dmtp-tapb", "label": "DMTP + TAPB · imine", "linkage": "imine", "monomer_a": MONOMERS[4]["smiles"], "monomer_b": MONOMERS[3]["smiles"], "label_a": "DMTP", "label_b": "TAPB", "instances": 12},
    {"id": "tp-paso3h", "label": "Tp + Pa-SO3H · keto-enamine", "linkage": "beta_ketoenamine", "monomer_a": MONOMERS[1]["smiles"], "monomer_b": MONOMERS[5]["smiles"], "label_a": "Tp", "label_b": "Pa-SO3H", "instances": 12},
]


def read_catalog(data_root: Path) -> dict:
    path = data_root / "framework-drawing" / "catalog" / "current.json"
    if not path.is_file():
        path = data_root / "framework-drawing" / "catalog-v3" / "catalog.json"
    if not path.is_file():
        path = data_root / "framework-drawing" / "catalog-v2" / "catalog.json"
    if not path.is_file():
        return {}
    if path.stat().st_size > 4_000_000:
        raise ValueError("Drawing catalog is too large")
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise ValueError("Drawing catalog must be an object")
    return value


def _number(value, limit=10000):
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value) or abs(value) > limit:
        raise ValueError("Sketch coordinates must be finite numbers")
    return float(value)


def normalize_monomer(*, smiles=None, molblock=None, graph=None, clean_layout=False) -> dict:
    if sum(x is not None for x in (smiles, molblock, graph)) != 1:
        raise ValueError("Supply one SMILES, MOL block or sketch graph")
    if smiles is not None:
        if len(smiles) > 10000:
            raise ValueError("SMILES is too long")
        mol = Chem.MolFromSmiles(smiles)
    elif molblock is not None:
        if len(molblock) > 500000:
            raise ValueError("MOL file is too large")
        mol = Chem.MolFromMolBlock(molblock, sanitize=True, removeHs=False)
    else:
        if not isinstance(graph, dict):
            raise ValueError("Invalid sketch graph")
        atoms, bonds = graph.get("atoms", []), graph.get("bonds", [])
        if not 1 <= len(atoms) <= 400 or len(bonds) > 700:
            raise ValueError("Use 1–400 atoms and at most 700 bonds per monomer")
        rw = Chem.RWMol()
        conf = Chem.Conformer(len(atoms)); conf.Set3D(False)
        for i, row in enumerate(atoms):
            atom = Chem.Atom(row["element"])
            atom.SetFormalCharge(int(row.get("charge", 0)))
            atom.SetIsotope(int(row.get("isotope", 0)))
            atom.SetNumExplicitHs(int(row.get("explicit_hs", 0)))
            atom.SetNoImplicit(bool(row.get("no_implicit", False)))
            atom.SetIsAromatic(bool(row.get("aromatic", False)))
            atom.SetChiralTag(Chem.ChiralType.values[int(row.get("chiral_tag", 0))])
            rw.AddAtom(atom)
            conf.SetAtomPosition(i, (_number(row["x"]), _number(row["y"]), 0.0))
        kinds = {1: Chem.BondType.SINGLE, 2: Chem.BondType.DOUBLE, 3: Chem.BondType.TRIPLE, 1.5: Chem.BondType.AROMATIC}
        for row in bonds:
            a, b = row["a"], row["b"]
            if not isinstance(a, int) or not isinstance(b, int) or not 0 <= a < len(atoms) or not 0 <= b < len(atoms) or a == b:
                raise ValueError("Bond endpoint is missing")
            if row["order"] not in kinds or rw.GetBondBetweenAtoms(a, b):
                raise ValueError("Unsupported or duplicate bond")
            rw.AddBond(a, b, kinds[row["order"]])
        mol = rw.GetMol(); mol.AddConformer(conf)
        for row in bonds:
            bond = mol.GetBondBetweenAtoms(row["a"], row["b"])
            selectors = row.get("stereo_atoms", [])
            if selectors:
                if len(selectors) != 2 or any(not isinstance(x, int) or not 0 <= x < len(atoms) for x in selectors):
                    raise ValueError("Invalid double-bond stereochemistry selectors")
                bond.SetStereoAtoms(*selectors)
            bond.SetStereo(Chem.BondStereo.values[int(row.get("stereo", 0))])
            bond.SetBondDir(Chem.BondDir.values[int(row.get("direction", 0))])
        Chem.SanitizeMol(mol)
    if mol is None:
        raise ValueError("This structure could not be read; check atom valences and bond orders")
    if mol.GetNumAtoms() > 400 or len(Chem.GetMolFrags(mol)) != 1:
        raise ValueError("A monomer must be one connected structure with at most 400 atoms")
    if clean_layout or not mol.GetNumConformers():
        rdDepictor.Compute2DCoords(mol)
    Chem.AssignStereochemistry(mol, cleanIt=False, force=True)
    canonical = Chem.MolToSmiles(mol, isomericSmiles=True)
    # Explicit Kekulé bonds make ordinary benzene bond/atom edits unambiguous.
    editable = Chem.Mol(mol); Chem.Kekulize(editable, clearAromaticFlags=True)
    pos = editable.GetConformer().GetPositions()
    result_graph = {"atoms": [{"element": a.GetSymbol(), "x": float(pos[a.GetIdx()][0]), "y": float(pos[a.GetIdx()][1]), "charge": a.GetFormalCharge(), "isotope": a.GetIsotope(), "aromatic": a.GetIsAromatic(), "chiral_tag": int(a.GetChiralTag()), "explicit_hs": a.GetNumExplicitHs(), "no_implicit": a.GetNoImplicit()} for a in editable.GetAtoms()],
                    "bonds": [{"a": b.GetBeginAtomIdx(), "b": b.GetEndAtomIdx(), "order": b.GetBondTypeAsDouble(), "stereo": int(b.GetStereo()), "stereo_atoms": list(b.GetStereoAtoms()), "direction": int(b.GetBondDir())} for b in editable.GetBonds()]}
    drawer = rdMolDraw2D.MolDraw2DSVG(600, 400)
    drawer.drawOptions().useBWAtomPalette(); drawer.drawOptions().includeMetadata = False
    rdMolDraw2D.PrepareAndDrawMolecule(drawer, mol); drawer.FinishDrawing()
    from .cofdraw_chemistry import inspect_monomer
    return {"valid": True, "smiles": canonical, "molblock": Chem.MolToMolBlock(mol), "graph": result_graph, "svg": drawer.GetDrawingText(), "formula": rdMolDescriptors.CalcMolFormula(mol), **inspect_monomer(canonical)}


def chosen_style(style_id: str, seed=None):
    selected,seed,_=style_configuration(style_id,seed)
    return selected,seed


def render_styled(recipe, directory, size, style_id, seed=None, style_options=None, *, prepared_patch=None):
    from .framework_drawing import build_patch, render_recipe
    from .source_depiction import render_source_svg
    style_id, seed, options = style_configuration(style_id, seed, style_options)
    native_standard = style_id == "standard" and not style_options
    directory = Path(directory)
    # Custom SVG styles supply their own images; do not render and discard two
    # potentially expensive native macrocycle images first.
    if prepared_patch is None:
        prepared_patch = build_patch(recipe)
    audit = render_recipe(recipe, directory, size=size, draw_images=native_standard, prepared_patch=prepared_patch)
    if native_standard:
        audit["drawing_style"] = {"id": style_id, "seed": seed, "options": options, "category":"chemical", "chemistry_changed":False}
        return audit
    mol, _ = prepared_patch
    if any(a.GetChiralTag() != Chem.ChiralType.CHI_UNSPECIFIED for a in mol.GetAtoms()) or any(b.GetStereo() not in (Chem.BondStereo.STEREONONE, Chem.BondStereo.STEREOANY) for b in mol.GetBonds()):
        raise ValueError("Use Standard style for a stereospecified framework; custom literature styles do not yet render stereo marks")
    xy = np.array(mol.GetConformer().GetPositions())[:, :2]; xy[:, 1] *= -1
    span = np.ptp(xy, axis=0)
    scale = size * .82 / max(float(max(span)), 1)
    xy = (xy - (xy.min(0) + xy.max(0)) / 2) * scale + size / 2
    lengths = [np.linalg.norm(xy[b.GetBeginAtomIdx()] - xy[b.GetEndAtomIdx()]) for b in mol.GetBonds()]
    bond = float(np.median(lengths)); font = min(34, max(12, bond * .50))
    components = list(dict.fromkeys(o["component"] for o in audit["atom_origin"]))
    colors = PALETTES[options["palette"]]
    colored = style_id in ("red_blue","bold_color","stacked_eclipsed","stacked_slipped","perspective")
    style = {"font_family": options["font_family"], "atom_font_size": font*options["font_scale"],
             "bond_width": max(1, size / 800)*options["bond_scale"]*(1.6 if style_id=="bold_color" else 1),
             "double_offset": bond * .065, "ring_double_offset": bond * .16,
             "font_weight":"bold" if style_id=="bold_color" else "normal",
             "port_style": options["edge_style"],
             "component_colors": {c: colors[i % len(colors)] if colored else "#171717" for i, c in enumerate(components)}}
    def key(i):
        o = audit["atom_origin"][i]; return f"{o['instance']}:{o['fragment_atom_index']}"
    if style_id == "hetero_colors":
        atom_palette={7:"#245bd1",8:"#d82c2c",16:"#b18d0c",15:"#c05c22",5:"#bd598c",9:"#268a4b",17:"#268a4b",35:"#93552b",53:"#823c9c"}
        style["atom_colors"]={key(a.GetIdx()):atom_palette[a.GetAtomicNum()] for a in mol.GetAtoms() if a.GetAtomicNum() in atom_palette}
    if options["abbreviation_mode"] != "explicit":
        groups, used = [], set()
        for pattern, label in [("[OX2]-[CH2]-[CH3]", "OEt"), ("[OX2]-[CH3]", "OMe" if options["abbreviation_mode"]=="OMe" else "OCH3"), ("[SX4](=[OX1])(=[OX1])-[OX2H]", "SO3H")]:
            for hit in mol.GetSubstructMatches(Chem.MolFromSmarts(pattern)):
                if used.intersection(hit) or mol.GetAtomWithIdx(hit[0]).GetDegree() < 2:
                    continue
                parent = next((a.GetIdx() for a in mol.GetAtomWithIdx(hit[0]).GetNeighbors() if a.GetIdx() not in hit), None)
                if parent is None:
                    continue
                # Place the chemically attached S/O glyph next to the incoming
                # bond, including when the substituent points left.
                left = xy[hit[0], 0] < xy[parent, 0]
                displayed = {"SO3H":"HO3S", "OCH3":"H3CO", "OMe":"MeO", "OEt":"EtO"}.get(label, label) if left else label
                shift = (-1 if left else 1) * (len(displayed)-1) * style["atom_font_size"] * .56 / 2
                label_xy = [float(xy[hit[0], 0] + shift), float(xy[hit[0], 1])]
                groups.append({"anchor": key(hit[0]), "label": displayed, "label_xy": label_xy,
                               "attachment_atom": mol.GetAtomWithIdx(hit[0]).GetSymbol(),
                               "attachment_side": "right" if left else "left",
                               "hidden_atoms": [key(i) for i in hit[1:]]}); used.update(hit)
        style["collapsed_groups"] = groups
    annotation = {"canvas_size": [size, size], "atom_coordinates": [{"instance": o["instance"], "fragment_atom_index": o["fragment_atom_index"], "xy": list(p)} for o, p in zip(audit["atom_origin"], xy)], "style": style}
    stem = directory / recipe["recipe_id"]
    render_source_svg(mol, audit, annotation, stem.with_suffix(".svg"))
    svg=decorate_svg(stem.with_suffix(".svg").read_text(encoding="utf-8"),mol,audit,recipe,xy,style_id,options,size)
    stem.with_suffix(".svg").write_text(svg,encoding="utf-8")
    import resvg_py
    stem.with_suffix(".png").write_bytes(resvg_py.svg_to_bytes(svg_string=svg, background=None if options["transparent_background"] else "white"))
    descriptor=next(s for s in STYLES if s["id"]==style_id)
    audit["drawing_style"] = {"id": style_id, "seed": seed, "rules": style, "options":options,
                              "category":descriptor["category"],"representation_scope":descriptor["representation_scope"],
                              "source_registered": False, "chemistry_changed": False,
                              "graph_exports":"Underlying single finite chemical graph; visual layers and cartoon symbols are presentation only."}
    return audit


def export_chemical_formats(directory, stem, *, source_mol=None):
    """Add chemical files without re-perceiving aromaticity of a large pore.

    The MOL has already passed the source-scoped graph/geometry checks. A full
    sanitizer here can stall on symmetric macrocycles or invent a different
    aromatic circuit. Compare the explicit tables, hydrogen totals and stereo
    after raw parsing instead. Optional CDXML errors retain the existing outputs.
    """
    import xml.etree.ElementTree as ET
    directory = Path(directory)
    # A retry must not expose a stale CDXML from an earlier successful export.
    cdxml_target = directory / f"{stem}.cdxml"
    if cdxml_target.is_file():
        cdxml_target.unlink()
    mol_bytes = (directory / f"{stem}.mol").read_bytes()
    molblock = mol_bytes.decode("utf-8").replace("\r\n", "\n")
    (directory / f"{stem}.sdf").write_bytes(mol_bytes.rstrip() + b"\n$$$$\n")
    available = ["sdf"]
    report = {"mol_sha256": hashlib.sha256(mol_bytes).hexdigest(),
              "cdxml": "unavailable_in_this_rdkit_build",
              "validation_scope": "Direct CDXML chemical atom/bond/coordinate records, explicit external-connection points and local wedge stereo; no native CDXML reparse, global macrocycle aromaticity or canonical ranking.",
              "ChemDraw_desktop_import_tested": False,
              "style_scope": "Chemical graph/layout only; SVG presentation styles are not transferred."}

    def supplied_molecule_or_parse():
        # The renderer can pass the exact graph it just serialized. Bind that
        # private copy to the stored MOL before avoiding the costly parser.
        # Defined stereochemistry retains the established parser path; for
        # ordinary V2000 patches, transfer the literal unknown-double flags.
        if source_mol is not None:
            try:
                candidate = Chem.Mol(source_mol)
                serialized = Chem.MolToMolBlock(candidate).replace("\r\n", "\n")
                if serialized.rstrip() != molblock.rstrip():
                    raise ValueError("supplied molecule does not serialize to the stored MOL")
                if (candidate.GetStereoGroups()
                        or any(a.GetChiralTag() != Chem.ChiralType.CHI_UNSPECIFIED for a in candidate.GetAtoms())
                        or any(b.GetStereo() not in (Chem.BondStereo.STEREONONE, Chem.BondStereo.STEREOANY)
                               or b.GetBondDir() not in (Chem.BondDir.NONE, Chem.BondDir.EITHERDOUBLE)
                               for b in candidate.GetBonds())):
                    raise ValueError("defined atom/bond stereo retains the scoped MOL parser")
                lines = molblock.splitlines()
                if len(lines) < 4 or "V2000" not in lines[3]:
                    raise ValueError("non-V2000 input retains the scoped MOL parser")
                atom_count, bond_count = int(lines[3][:3]), int(lines[3][3:6])
                if (atom_count, bond_count) != (candidate.GetNumAtoms(), candidate.GetNumBonds()):
                    raise ValueError("stored MOL counts differ")
                Chem.Kekulize(candidate, clearAromaticFlags=True)
                for i, line in enumerate(lines[4 + atom_count:4 + atom_count + bond_count]):
                    a, b, order, stereo = (int(line[j:j+3]) for j in (0, 3, 6, 9))
                    bond = candidate.GetBondWithIdx(i)
                    if (a-1, b-1, order) != (bond.GetBeginAtomIdx(), bond.GetEndAtomIdx(), int(bond.GetBondTypeAsDouble())):
                        raise ValueError("stored MOL explicit bond table differs")
                    if stereo == 3 and order == 2:
                        bond.SetStereo(Chem.BondStereo.STEREOANY)
                        bond.SetBondDir(Chem.BondDir.EITHERDOUBLE)
                    elif stereo != 0:
                        raise ValueError("MOL stereo mark requires the scoped parser")
                # Use the file's rounded coordinates, making this path match
                # its parser fallback even for a high-precision source graph.
                conf = candidate.GetConformer()
                for i, line in enumerate(lines[4:4+atom_count]):
                    conf.SetAtomPosition(i, tuple(float(line[j:j+10]) for j in (0, 10, 20)))
                candidate.UpdatePropertyCache(strict=True)
                report["molecule_input"] = "verified_source_mol_fast_path"
                report["source_mol_serialization_matches_stored_MOL"] = True
                return candidate
            except Exception as exc:
                report["source_mol_fallback_reason"] = f"{type(exc).__name__}: {exc}"
        report["molecule_input"] = "scoped_MOL_parser"
        return Chem.MolFromMolBlock(molblock, sanitize=False, removeHs=False)

    def table(mol):
        mol.UpdatePropertyCache(strict=True)
        Chem.FastFindRings(mol)
        Chem.SetBondStereoFromDirections(mol)
        atoms = [(a.GetAtomicNum(), a.GetFormalCharge(), a.GetIsotope(),
                  a.GetTotalNumHs(includeNeighbors=True), int(a.GetChiralTag()))
                 for a in mol.GetAtoms()]
        # The native writer/parser preserves atom order. Fail safely if an
        # unsupported writer changes it instead of accepting a lossy checksum.
        bonds = [(b.GetBeginAtomIdx(), b.GetEndAtomIdx(), str(b.GetBondType()),
                  int(b.GetStereo()), tuple(b.GetStereoAtoms()), int(b.GetBondDir()))
                 for b in mol.GetBonds()]
        return atoms, bonds

    def normalize_continuation_points(data, mol):
        # RDKit's native writer serializes a dummy atom as Element="0". The
        # ChemDraw desktop rejects that as an invalid atomic number. In the
        # ChemDraw SDK these are a different node type, not chemical elements:
        # NodeType=ExternalConnectionPoint / ExternalConnectionType=Wavy.
        # Never repair this by inventing a carbon, hydrogen or R-group molecule.
        root = ET.fromstring(data)
        fragments = list(root.iter("fragment"))
        if len(fragments) != 1:
            raise ValueError("CDXML must contain exactly one chemical fragment")
        nodes = list(fragments[0].findall("n"))
        if len(nodes) != mol.GetNumAtoms():
            raise ValueError("CDXML node count changed before continuation conversion")
        converted = []
        for atom, row in zip(mol.GetAtoms(), nodes):
            if atom.GetAtomicNum() != 0:
                continue
            if (atom.GetDegree() != 1 or atom.GetFormalCharge() or atom.GetIsotope()
                    or atom.GetTotalNumHs(includeNeighbors=True)
                    or atom.GetChiralTag() != Chem.ChiralType.CHI_UNSPECIFIED):
                raise ValueError("Only neutral degree-one unlabelled boundary dummies can become CDXML connection points")
            if row.get("NodeType", "Element") != "Element" or row.get("Element") != "0":
                raise ValueError("Unexpected native writer boundary-node representation")
            # Remove fields that assert chemical-element identity. Their values
            # have already been checked on the original graph above.
            for key in ("Element", "NumHydrogens", "Isotope", "Charge", "Radical", "AS"):
                row.attrib.pop(key, None)
            for child in list(row):
                if child.tag != "t":
                    raise ValueError("Unsupported child of a CDXML boundary node")
                row.remove(child)  # Native '*' text must not override port display.
            row.set("NodeType", "ExternalConnectionPoint")
            row.set("ExternalConnectionType", "Wavy")
            converted.append({"atom_index": atom.GetIdx(), "cdxml_node_id": row.get("id"),
                              "neighbor_atom_index": atom.GetNeighbors()[0].GetIdx(),
                              "bond_order": atom.GetBonds()[0].GetBondTypeAsDouble()})
        if root.get("BondLength") == "":
            root.attrib.pop("BondLength")  # Native writer emits an invalid empty number.
        if any(row.get("Element") == "0" for row in root.iter("n")):
            raise ValueError("CDXML contains an invalid zero-element node")
        text = '<?xml version="1.0" encoding="UTF-8"?>\n'
        text += '<!DOCTYPE CDXML SYSTEM "https://static.chemistry.revvitycloud.com/cdxml/CDXML.dtd">\n'
        return text + ET.tostring(root, encoding="unicode"), converted

    def read_native_xml(data):
        # Native CDXML readers perform expensive stereo/ring perception even
        # with sanitize=False on some large symmetric pores. The writer emits
        # a bounded, simple graph subset; parse that subset directly instead.
        if len(data) > 5_000_000:
            raise ValueError("CDXML exceeds export validation size limit")
        root = ET.fromstring(data)
        fragments = list(root.iter("fragment"))
        if len(fragments) != 1:
            raise ValueError("CDXML must contain exactly one chemical fragment")
        nodes, edges = list(fragments[0].findall("n")), list(fragments[0].findall("b"))
        rw = Chem.RWMol(); coordinates = []; ids = {}
        for row in nodes:
            node_type = row.get("NodeType", "Element")
            if row.get("id") in ids or node_type not in ("Element", "ExternalConnectionPoint"):
                raise ValueError("Unsupported or duplicate CDXML node")
            if node_type == "ExternalConnectionPoint":
                if row.get("Element") is not None or row.get("ExternalConnectionType") != "Wavy":
                    raise ValueError("Invalid CDXML external-connection semantics")
                if any(key in row.attrib for key in ("Charge", "Isotope", "NumHydrogens", "Radical")) or list(row):
                    raise ValueError("Chemical atom properties attached to a boundary connection point")
                atom = Chem.Atom(0)
                atom.SetNoImplicit(True)
            else:
                element = int(row.get("Element", "6"))
                if not 1 <= element <= 118:
                    raise ValueError("Invalid CDXML chemical element number")
                atom = Chem.Atom(element)
            atom.SetFormalCharge(int(row.get("Charge", "0")))
            atom.SetIsotope(int(row.get("Isotope", "0")))
            if row.get("NumHydrogens") is not None:
                atom.SetNumExplicitHs(int(row.get("NumHydrogens")))
                atom.SetNoImplicit(True)
            if row.get("Radical"):
                raise ValueError("Radical CDXML export is not validated")
            ids[row.get("id")] = rw.AddAtom(atom)
            xy = [float(x) for x in row.get("p", "").split()]
            if len(xy) != 2 or not all(math.isfinite(x) for x in xy):
                raise ValueError("Invalid CDXML atom coordinates")
            coordinates.append((xy[0], -xy[1], 0.0))
        orders = {"1": Chem.BondType.SINGLE, "2": Chem.BondType.DOUBLE,
                  "3": Chem.BondType.TRIPLE, "1.5": Chem.BondType.AROMATIC,
                  "dative": Chem.BondType.DATIVE}
        for row in edges:
            a, b = ids[row.get("B")], ids[row.get("E")]
            order = orders[row.get("Order", "1")]
            display = row.get("Display", "Solid")
            if display in ("WedgeEnd", "WedgedHashEnd"):
                a, b = b, a
            rw.AddBond(a, b, order); bond = rw.GetBondBetweenAtoms(a, b)
            if display in ("WedgeBegin", "WedgeEnd"):
                bond.SetBondDir(Chem.BondDir.BEGINWEDGE)
            elif display in ("WedgedHashBegin", "WedgedHashEnd"):
                bond.SetBondDir(Chem.BondDir.BEGINDASH)
            elif display == "Wavy":
                bond.SetBondDir(Chem.BondDir.EITHERDOUBLE if order == Chem.BondType.DOUBLE else Chem.BondDir.UNKNOWN)
                if order == Chem.BondType.DOUBLE:
                    bond.SetStereo(Chem.BondStereo.STEREOANY)
            elif display != "Solid":
                raise ValueError("Unsupported CDXML bond display")
        result = rw.GetMol(); conf = Chem.Conformer(len(nodes)); conf.Set3D(False)
        for i, xyz in enumerate(coordinates):
            conf.SetAtomPosition(i, xyz)
        result.AddConformer(conf); result.UpdatePropertyCache(strict=True)
        Chem.FastFindRings(result)
        Chem.AssignChiralTypesFromBondDirs(result, replaceExistingTags=True)
        return result

    def fit_printable_page(data):
        # The native writer leaves coordinates centred around (0,0) and omits
        # page setup. ChemDraw then inherits the last printer/style (e.g. B4)
        # and can open with most of the structure outside the visible page.
        # Use one positive, uniform transform: no reflection or chemistry edit.
        # Letter paper is 612x792 pt; 36 pt margins leave a 540x720 drawing page.
        import struct
        root = ET.fromstring(data)
        pages = list(root.findall("page"))
        fragments = list(root.iter("fragment"))
        if len(pages) != 1 or len(fragments) != 1:
            raise ValueError("Page fitting requires one chemical drawing space")
        nodes = list(fragments[0].findall("n"))
        points = np.asarray([[float(x) for x in n.get("p", "").split()] for n in nodes])
        if points.shape != (len(nodes), 2) or not len(points) or not np.isfinite(points).all():
            raise ValueError("Page fitting requires finite 2D coordinates")
        ids = {n.get("id"): i for i, n in enumerate(nodes)}
        lengths = [float(np.linalg.norm(points[ids[b.get("B")]] - points[ids[b.get("E")]]))
                   for b in fragments[0].findall("b")]
        positive_lengths = [x for x in lengths if x > 1e-6]
        if lengths and len(positive_lengths) != len(lengths):
            raise ValueError("A collapsed chemical bond cannot be fitted to the page")
        native_bond = float(np.median(positive_lengths)) if positive_lengths else 14.4
        preferred_scale = 14.4 / native_bond
        page_size = np.asarray([540., 720.])
        padding = 24.  # label/charge/continuation-wave allowance inside print area
        extent = np.ptp(points, axis=0)
        fit_scale = min((page_size[i] - 2 * padding) / extent[i]
                        for i in range(2) if extent[i] > 1e-8) if np.any(extent > 1e-8) else preferred_scale
        scale = min(preferred_scale, fit_scale)
        offset = page_size / 2 - scale * (points.min(0) + points.max(0)) / 2
        for row in root.iter():
            if row.get("p"):
                xy = np.asarray([float(x) for x in row.get("p").split()])
                if xy.shape != (2,) or not np.isfinite(xy).all():
                    raise ValueError("Unsupported CDXML position while fitting page")
                row.set("p", " ".join(f"{v:.6f}" for v in xy * scale + offset))
            if row.get("BoundingBox"):
                box = np.asarray([float(x) for x in row.get("BoundingBox").split()])
                if box.shape != (4,) or not np.isfinite(box).all():
                    raise ValueError("Unsupported CDXML bounding box")
                row.set("BoundingBox", " ".join(f"{v:.6f}" for v in box * scale + np.tile(offset, 2)))
        label_scale = min(1., scale / preferred_scale)
        root.set("BondLength", f"{native_bond * scale:.6f}")
        root.set("LabelSize", f"{10 * label_scale:.6f}")
        root.set("CaptionSize", f"{10 * label_scale:.6f}")
        root.set("LabelFont", "3"); root.set("CaptionFont", "3")
        root.set("LineWidth", f"{.6 * label_scale:.6f}")
        root.set("BoldWidth", f"{2 * label_scale:.6f}")
        root.set("HashSpacing", f"{2.5 * label_scale:.6f}")
        root.set("MarginWidth", f"{1.6 * label_scale:.6f}")
        root.set("BondSpacing", "18")
        root.set("PrintMargins", "36 36 36 36")
        root.set("WindowIsZoomed", "yes")
        for field in ("WindowPosition", "WindowSize", "WinPrintInfo"):
            root.attrib.pop(field, None)
        # SDK-documented Macintosh TPrint fields also set Windows ChemDraw's
        # paper geometry; WinPrintInfo is documented as not read/written.
        printing = [0] * 60
        printing[2:4] = [72, 72]
        printing[4:8] = [0, 0, 792, 612]
        printing[8:12] = [0, 0, 792, 612]
        printing[12:15] = [-1, 792, 612]
        printing[24:26] = [100, 100]
        root.set("MacPrintInfo", struct.pack(">60h", *printing).hex().upper())
        fonts = root.find("fonttable")
        if fonts is None:
            fonts = ET.Element("fonttable"); root.insert(0, fonts)
        if not any(f.get("id") == "3" for f in fonts):
            ET.SubElement(fonts, "font", {"id": "3", "charset": "iso-8859-1", "name": "Arial"})
        for run in root.iter("s"):
            if run.get("size"):
                run.set("size", f"{float(run.get('size')) * scale:.6f}")
        page = pages[0]
        page.set("BoundingBox", "0 0 540 720")
        page.set("DrawingSpace", "Pages")
        page.set("WidthPages", "1"); page.set("HeightPages", "1")
        page.set("HeaderPosition", "36"); page.set("FooterPosition", "36")
        final_points = points * scale + offset
        lower = final_points.min(0) - padding / 2
        upper = final_points.max(0) + padding / 2
        if np.any(lower < 0) or np.any(upper > page_size):
            raise ValueError("Chemical drawing does not fit within printable page")
        bounds = " ".join(f"{v:.6f}" for v in np.r_[lower, upper])
        fragments[0].set("BoundingBox", bounds)
        root.set("BoundingBox", bounds)
        text = '<?xml version="1.0" encoding="UTF-8"?>\n'
        text += '<!DOCTYPE CDXML SYSTEM "https://static.chemistry.revvitycloud.com/cdxml/CDXML.dtd">\n'
        return text + ET.tostring(root, encoding="unicode"), {
            "paper": "US Letter portrait", "paper_size_points": [612, 792],
            "print_margins_points": [36, 36, 36, 36], "drawing_page_points": [540, 720],
            "width_pages": 1, "height_pages": 1, "uniform_scale": scale,
            "translation_points": offset.tolist(), "node_bounds_points": np.r_[final_points.min(0), final_points.max(0)].tolist(),
            "fragment_bounds_points": np.r_[lower, upper].tolist(),
            "label_size_points": 10 * label_scale,
            "scope": "Single-page placement with positive coordinates; original MOL/SDF geometry is unchanged."}

    if hasattr(Chem, "HasChemDrawCDXSupport") and Chem.HasChemDrawCDXSupport():
        try:
            mol = supplied_molecule_or_parse()
            if mol is None:
                raise ValueError("Exported chemical MOL could not be reloaded")
            table(mol)
            if mol.GetStereoGroups() or any(a.GetNumRadicalElectrons() for a in mol.GetAtoms()):
                raise ValueError("Enhanced stereo groups/radicals need a separately validated CDXML export")
            # The native writer rounds coordinates. Centre and normalize its
            # private copy first so tiny drawings or large absolute origins do
            # not collapse bonds before the final page transform can run.
            writer_mol = Chem.Mol(mol)
            writer_xy = np.asarray(writer_mol.GetConformer().GetPositions())[:, :2]
            writer_lengths = [float(np.linalg.norm(writer_xy[b.GetBeginAtomIdx()] - writer_xy[b.GetEndAtomIdx()]))
                              for b in writer_mol.GetBonds()]
            if any(x <= 1e-8 for x in writer_lengths):
                raise ValueError("Original chemical drawing has a collapsed bond")
            writer_scale = 1.5 / float(np.median(writer_lengths)) if writer_lengths else 1.
            writer_center = writer_xy.mean(0)
            for i, xy in enumerate((writer_xy - writer_center) * writer_scale):
                writer_mol.GetConformer().SetAtomPosition(i, (float(xy[0]), float(xy[1]), 0.))
            data, boundary_points = normalize_continuation_points(Chem.MolToCDXMLBlock(writer_mol), mol)
            data, page_geometry = fit_printable_page(data)
            page_geometry["writer_input_normalization"] = {"uniform_scale": writer_scale,
                                                            "subtracted_center": writer_center.tolist()}
            restored = read_native_xml(data)
            def atom_table(value):
                return [(a.GetAtomicNum(), a.GetFormalCharge(), a.GetIsotope(),
                         a.GetTotalNumHs(includeNeighbors=True), int(a.GetChiralTag())) for a in value.GetAtoms()]
            def bond_table(value):
                return sorted((b.GetBeginAtomIdx(), b.GetEndAtomIdx(), str(b.GetBondType())) if b.GetBondType() == Chem.BondType.DATIVE
                              else (*sorted((b.GetBeginAtomIdx(), b.GetEndAtomIdx())), str(b.GetBondType())) for b in value.GetBonds())
            if atom_table(restored) != atom_table(mol) or bond_table(restored) != bond_table(mol):
                raise ValueError("CDXML roundtrip changed atom, bond, hydrogen or stereo tables")
            # A uniform translation/positive scale is allowed by chemical
            # drawing units, but not reflection or coordinate collapse.
            original_xy = np.asarray(mol.GetConformer().GetPositions())[:, :2]
            returned_xy = np.asarray(restored.GetConformer().GetPositions())[:, :2]
            if len(original_xy) > 1:
                a = original_xy - original_xy.mean(0)
                b = returned_xy - returned_xy.mean(0)
                norm = float((a * a).sum())
                if norm > 1e-12:
                    scale = float((a * b).sum()) / norm
                    if scale <= 0 or not np.allclose(b, scale * a, atol=max(1e-3, np.ptp(b) * 1e-4)):
                        raise ValueError("CDXML roundtrip changed planar drawing geometry")
            # Defined alkene/imine geometry is preserved by the coordinate
            # check; Wavy instead means explicitly unknown and must not replace
            # a defined stereo bond (or disappear from an unknown bond).
            for bond in mol.GetBonds():
                other = restored.GetBondBetweenAtoms(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx())
                if bond.GetBondType() == Chem.BondType.DOUBLE:
                    if (bond.GetStereo() == Chem.BondStereo.STEREOANY) != (other.GetStereo() == Chem.BondStereo.STEREOANY):
                        raise ValueError("CDXML changed unknown/defined double-bond stereo")
                elif (bond.GetBondDir() == Chem.BondDir.UNKNOWN) != (other.GetBondDir() == Chem.BondDir.UNKNOWN):
                    raise ValueError("CDXML changed an explicitly unknown single-bond stereo mark")
            cdxml_target.write_text(data, encoding="utf-8")
            available.append("cdxml")
            report["cdxml"] = "validated_native_writer_explicit_XML_tables"
            report["boundary_representation"] = "ExternalConnectionPoint/Wavy; no Element=0 or invented C/H cap"
            report["boundary_connection_points"] = boundary_points
            report["page_geometry"] = page_geometry
            report["ChemDraw_desktop_regression"] = "Release cases passed hidden native ChemDraw graph/stereo and single-page bounds checks; see benchmarks/cofdraw-chemdraw-page-fit.json. This particular generated file has not been independently opened in ChemDraw."
        except Exception as exc:
            report["cdxml"] = "omitted_validation_or_writer_failure"
            report["reason"] = f"{type(exc).__name__}: {exc}"
    (directory / f"{stem}.chemical-export.json").write_text(
        json.dumps(report, indent=2) + "\n", encoding="utf-8")
    return available
