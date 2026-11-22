"""Portable, offline COFdraw monomer editing and presentation styles.

Presentation changes never change the chemical graph or establish a COF identity.
"""
from __future__ import annotations

from copy import deepcopy
import hashlib
import json
from functools import lru_cache
import math
from pathlib import Path
import random

import numpy as np
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem.Draw import rdMolDraw2D

from .cofdraw_styles import STYLES, PALETTES, ELEMENT_SCHEMES, style_configuration, decorate_svg, resolved_palette
from .cofdraw_depiction import compute_2d_coords
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


@lru_cache(maxsize=1)
def _reviewed_feeds():
    path=Path(__file__).with_name('web')/'cofdraw-reviewed-feeds.json'
    return json.loads(path.read_text(encoding='utf8'))


def catalog_monomers(catalog):
    from .cofdraw_chemistry import enable_reviewed_feed
    from .cofdraw_catalog_review import reviewed_row
    from .cofdraw_self import additional_precursors
    from .cofdraw_mixed import TPDA
    occurrences = catalog.get("recent_monomer_update", {}).get("additional_occurrences", {})
    rows = []
    for monomer in catalog.get("monomer_presets", MONOMERS):
        extra = occurrences.get(monomer["id"])
        if extra:
            monomer = {
                **monomer,
                "aliases": list(dict.fromkeys([*monomer.get("aliases", []), *(e["source_label"] for e in extra)])),
                "source_dois": sorted({*monomer.get("source_dois", []), *(e["source_doi"] for e in extra)}),
                "source_evidence": [*monomer.get("source_evidence", []), *extra],
            }
        rows.append(reviewed_row(monomer))
    seen = {row['id'] for row in rows}
    from .cofdraw_starters import tpb_precursors
    rows.extend(row for row in [*additional_precursors(),TPDA,*tpb_precursors()] if row['id'] not in seen)
    by_smiles={Chem.CanonSmiles(row['smiles']):i for i,row in enumerate(rows)}
    for feed in _reviewed_feeds():
        index=by_smiles.get(feed['smiles'])
        admission=dict(catalog_use='framework_feed',automatic_prepare_supported=True)
        if index is None:
            rows.append({**feed,**admission,'label':feed['abbreviation']})
        else:
            old=rows[index]
            rows[index]={**old,**admission,
                'source_dois':sorted(set(old.get('source_dois',[]))|set(feed['source_dois'])),
                'source_evidence':[*old.get('source_evidence',[]),*feed['source_evidence']]}
    return [enable_reviewed_feed(row) for row in rows]


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
            radicals = row.get("radical_electrons", 0)
            if type(radicals) is not int or not 0 <= radicals <= 8:
                raise ValueError("Radical electron count must be an integer between 0 and 8")
            atom.SetNumRadicalElectrons(radicals)
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
            kind = row.get("kind")
            if kind not in (None, "dative") or (kind == "dative" and row["order"] != 1):
                raise ValueError("Unsupported bond kind")
            rw.AddBond(a, b, Chem.BondType.DATIVE if kind == "dative" else kinds[row["order"]])
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
        compute_2d_coords(mol)
    Chem.AssignStereochemistry(mol, cleanIt=False, force=True)
    canonical = Chem.MolToSmiles(mol, isomericSmiles=True)
    # Explicit Kekulé bonds make ordinary benzene bond/atom edits unambiguous.
    editable = Chem.Mol(mol); Chem.Kekulize(editable, clearAromaticFlags=True)
    pos = editable.GetConformer().GetPositions()
    result_graph = {"atoms": [{"element": a.GetSymbol(), "x": float(pos[a.GetIdx()][0]), "y": float(pos[a.GetIdx()][1]), "charge": a.GetFormalCharge(), "isotope": a.GetIsotope(), "aromatic": a.GetIsAromatic(), "chiral_tag": int(a.GetChiralTag()), "explicit_hs": a.GetNumExplicitHs(), "no_implicit": a.GetNoImplicit(), **({"radical_electrons": a.GetNumRadicalElectrons()} if a.GetNumRadicalElectrons() else {})} for a in editable.GetAtoms()],
                    "bonds": [{"a": b.GetBeginAtomIdx(), "b": b.GetEndAtomIdx(), "order": b.GetBondTypeAsDouble(),
                               **({"kind": "dative"} if b.GetBondType() == Chem.BondType.DATIVE else {}),
                               "stereo": int(b.GetStereo()), "stereo_atoms": list(b.GetStereoAtoms()),
                               "direction": int(b.GetBondDir())} for b in editable.GetBonds()]}
    drawer = rdMolDraw2D.MolDraw2DSVG(600, 400)
    drawer.drawOptions().useBWAtomPalette(); drawer.drawOptions().includeMetadata = False
    rdMolDraw2D.PrepareAndDrawMolecule(drawer, mol); drawer.FinishDrawing()
    from .cofdraw_chemistry import inspect_monomer
    return {"valid": True, "smiles": canonical, "molblock": Chem.MolToMolBlock(mol), "graph": result_graph, "svg": drawer.GetDrawingText(), "formula": rdMolDescriptors.CalcMolFormula(mol), **inspect_monomer(canonical)}


def chosen_style(style_id: str, seed=None):
    selected,seed,_=style_configuration(style_id,seed)
    return selected,seed


def _native_styled_svg(mol, audit, recipe, style_id, options, rules, size, model_coordinates=None):
    """Use RDKit's chemical label and bond placement for every presentation.

    Native SVG paths retain wedge/stereo geometry and clip bonds at the actual
    glyphs. Decorations and colors never ask a second renderer to reinterpret
    the bond geometry. Collapsed aliases affect only a drawing copy.
    """
    import re
    import xml.etree.ElementTree as ET
    from rdkit.Geometry import Point2D
    from .depiction_rings import symmetric_depiction_rings
    tag = lambda name: '{http://www.w3.org/2000/svg}'+name
    ET.register_namespace('', 'http://www.w3.org/2000/svg')
    original_xy = np.asarray(mol.GetConformer().GetPositions())[:, :2]
    graph_keys = [f"{o['instance']}:{o['fragment_atom_index']}" for o in audit['atom_origin']]
    hidden = {key for group in rules.get('collapsed_groups', []) for key in group['hidden_atoms']}
    kept = [i for i, key in enumerate(graph_keys) if key not in hidden]
    old_to_new = {old: new for new, old in enumerate(kept)}
    draw_mol = Chem.Mol(mol)
    Chem.Kekulize(draw_mol, clearAromaticFlags=True)
    if hidden:
        editable = Chem.RWMol(draw_mol)
        for i in reversed(range(mol.GetNumAtoms())):
            if graph_keys[i] in hidden:
                editable.RemoveAtom(i)
        draw_mol = editable.GetMol()
        draw_mol.UpdatePropertyCache(strict=False)
    Chem.WedgeMolBonds(draw_mol, draw_mol.GetConformer())
    symmetric_depiction_rings(draw_mol)
    drawer = rdMolDraw2D.MolDraw2DSVG(size, size)
    opt = drawer.drawOptions()
    layered_front = style_id in ('gradient_stack','half_space_stack')
    opt.useBWAtomPalette(); opt.padding = .16 if layered_front else .10
    if layered_front:
        drawer.SetOffset(0, -round(size*.07))
    opt.bondLineWidth = (2.8 if style_id == 'bold_color' else 2.)*options['bond_scale']
    opt.minFontSize = round(12*options['font_scale']); opt.maxFontSize = round(30*options['font_scale'])
    opt.baseFontSize = .6*options['font_scale']
    import os
    font_names = {'Arial': 'arial.ttf', 'Times New Roman': 'times.ttf', 'Cambria': 'cambria.ttc',
                  'Calibri': 'calibri.ttf', 'Verdana': 'verdana.ttf', 'Courier New': 'cour.ttf', 'Segoe UI': 'segoeui.ttf'}
    native_font = Path(os.environ.get('WINDIR', 'C:/Windows'))/'Fonts'/font_names[options['font_family']]
    if native_font.is_file():
        opt.fontFile = str(native_font)
    opt.addStereoAnnotation = False; opt.includeMetadata = False
    opt.prepareMolsBeforeDrawing = False; opt.clearBackground = False
    # Built-in font files keep exports portable; system families remain in the
    # style metadata and schematic labels without a host-specific file path.
    if options['hydrogen_labels'] == 'compact':
        for old in kept:
            atom = mol.GetAtomWithIdx(old)
            if (atom.GetAtomicNum() in (7, 8, 15, 16) and atom.GetTotalNumHs() and
                    not (atom.GetFormalCharge() or atom.GetIsotope() or atom.GetNumRadicalElectrons()
                         or atom.GetIsAromatic() or atom.GetAtomMapNum() or atom.GetChiralTag())):
                opt.atomLabels[old_to_new[old]] = atom.GetSymbol()
    for group in rules.get('collapsed_groups', []):
        old = graph_keys.index(group['anchor'])
        opt.atomLabels[old_to_new[old]] = group['label']
    for i, atom in enumerate(draw_mol.GetAtoms()):
        if atom.GetAtomicNum() == 0 and style_id != 'standard':
            opt.atomLabels[i] = ''
    drawer.DrawMolecule(draw_mol)
    drawer.FinishDrawing()
    def project(points):
        pixels = [drawer.GetDrawCoords(Point2D(float(x), float(y))) for x, y in points]
        return np.array([(p.x, p.y) for p in pixels])
    xy = project(original_xy)
    root = ET.fromstring(drawer.GetDrawingText())
    defs = ET.Element(tag('defs')); root.insert(0, defs)
    atom_colors = {}
    source_style = recipe.get('source_style', {}) if style_id == 'source' else {}
    source_colors = source_style.get('atom_colors', {})
    for old in kept:
        atom_colors[old_to_new[old]] = source_colors.get(str(old), source_colors.get(old,
            rules.get('atom_colors', {}).get(graph_keys[old], source_style.get('component_colors', {}).get(audit['atom_origin'][old]['component'],
                rules['component_colors'].get(audit['atom_origin'][old]['component'], '#171717')))))
    def recolor(element, color):
        css = element.get('style', '')
        css = re.sub(r'(?<!-)stroke:(?!none)[^;]+', 'stroke:'+color, css)
        css = re.sub(r'(?<!-)fill:(?!none)[^;]+', 'fill:'+color, css)
        element.set('style', css)
        if element.get('fill') not in (None, 'none'): element.set('fill', color)
        if element.get('stroke') not in (None, 'none'): element.set('stroke', color)
    gradients = set()
    kept_xy = xy[kept]
    for element in list(root):
        classes = element.get('class', '')
        atoms = [int(v) for v in re.findall(r'atom-(\d+)', classes)]
        bond_match = re.search(r'bond-(\d+)', classes)
        if not atoms and element.tag.endswith('path'):
            # RDKit adds a short, unclassified polyline at carbon vertices to
            # close butt-ended bonds. Color that join with its exact nearest
            # atom too; leaving it black breaks otherwise colored rings.
            match = re.fullmatch(r'M\s+([-\d.]+)[ ,]+([-\d.]+)\s+L\s+([-\d.]+)[ ,]+([-\d.]+)\s+L\s+([-\d.]+)[ ,]+([-\d.]+)\s*', element.get('d', ''))
            if match:
                values = list(map(float, match.groups()))
                point = np.array(values[2:4])
                delta = kept_xy - point
                squared = np.einsum('ij,ij->i',delta,delta)
                position = int(np.argmin(squared))
                nearest = kept[position]
                if squared[position] < 1.:
                    atoms = [old_to_new[nearest]]
                    element.set('data-vertex-atom', str(nearest))
        if not atoms:
            continue
        if bond_match and len(atoms) >= 2:
            a, b = atoms[:2]; ca, cb = atom_colors[a], atom_colors[b]
            if ca != cb:
                name = 'bond-color-'+bond_match.group(1)
                if name not in gradients:
                    p, q = xy[kept[a]], xy[kept[b]]
                    gradient = ET.SubElement(defs, tag('linearGradient'), {'id': name, 'gradientUnits': 'userSpaceOnUse',
                        'x1': str(p[0]), 'y1': str(p[1]), 'x2': str(q[0]), 'y2': str(q[1])})
                    for stop, color in (('0%', ca), ('48%', ca), ('52%', cb), ('100%', cb)):
                        ET.SubElement(gradient, tag('stop'), {'offset': stop, 'stop-color': color})
                    gradients.add(name)
                recolor(element, 'url(#'+name+')')
            else:
                recolor(element, ca)
            if style_id != 'standard' and any(draw_mol.GetAtomWithIdx(i).GetAtomicNum() == 0 for i in atoms):
                if options['edge_style'] == 'dashed':
                    element.set('stroke-dasharray', f'{size*.003},{size*.003}')
        else:
            if style_id != 'standard' and draw_mol.GetAtomWithIdx(atoms[0]).GetAtomicNum() == 0:
                root.remove(element)
                continue
            recolor(element, atom_colors[atoms[0]])
        opacities = source_style.get('atom_opacities', {})
        values = [opacities.get(str(kept[a]), opacities.get(kept[a], source_style.get('fade_opacity', .28) if kept[a] in source_style.get('faded_atoms', []) else 1.)) for a in atoms]
        if values and min(values) < 1:
            element.set('opacity', str(max(0., min(1., float(min(values))))))
    if style_id != 'standard' and options['edge_style'] in ('wavy_cross', 'crossbar'):
        for port in audit['continuation_ports']:
            i = port['atom_index']; neighbors = mol.GetAtomWithIdx(i).GetNeighbors()
            if not neighbors: continue
            j = neighbors[0].GetIdx(); p, q = xy[j], xy[i]
            delta = q-p; length = np.linalg.norm(delta)
            if length < .01: continue
            tangent = delta/length; normal = np.array([-tangent[1], tangent[0]])
            center = p+.8*delta; span = length*.48
            samples = [center+(u-.5)*span*normal+(0 if options['edge_style'] == 'crossbar' else
                        .045*length*math.sin(6*math.pi*u))*tangent for u in np.linspace(0, 1, 35)]
            ET.SubElement(root, tag('path'), {'d': 'M '+' L '.join(f'{a:.3f},{b:.3f}' for a, b in samples),
                'fill': 'none', 'stroke': atom_colors.get(old_to_new.get(j), '#171717'),
                'stroke-width': str(opt.bondLineWidth), 'stroke-linecap': 'round', 'data-continuation': str(i)})
    return ET.tostring(root, encoding='unicode'), xy, xy if model_coordinates is None else project(model_coordinates)


def render_styled(recipe, directory, size, style_id, seed=None, style_options=None, *, prepared_patch=None, recent_styles=()):
    from .framework_drawing import build_patch, render_recipe
    from .cofdraw_boundary import aryl_cut_sites, aryl_figure
    if style_id == 'source':
        source_options = {key: recipe.get('source_style', {})[key] for key in ('abbreviation_mode','edge_style') if key in recipe.get('source_style', {})}
        source_options.update(style_options or {})
        style_options = source_options
    style_id, seed, options = style_configuration(style_id, seed, style_options, recent_styles)
    directory = Path(directory)
    # Custom SVG styles supply their own images; do not render and discard two
    # potentially expensive native macrocycle images first.
    if prepared_patch is None:
        prepared_patch = build_patch(recipe)
    audit = render_recipe(recipe, directory, size=size, draw_images=False, prepared_patch=prepared_patch)
    mol, _ = prepared_patch
    model_mol, model_audit = mol, audit
    boundary_mode = options['boundary_cut']
    shortened = 0
    if boundary_mode == 'aryl':
        cuts = aryl_cut_sites(mol)
        if cuts:
            mol, audit, recipe = aryl_figure(mol, audit, recipe, cuts)
            shortened = len(cuts)
    elif boundary_mode == 'hub':
        from .cofdraw_self_kagome import hub_figure
        figure = hub_figure(mol, audit, recipe)
        if figure is not None:
            mol, audit, recipe, shortened = figure
    model_coordinates = None
    editing_indices = None
    if shortened:
        model_coordinates = np.asarray(model_mol.GetConformer().GetPositions())[:, :2]
        kept = np.asarray(audit['model_indices'])
        positions = np.asarray(mol.GetConformer().GetPositions())[:, :2]
        editing_indices = kept[np.isclose(positions, model_coordinates[kept], atol=1e-8).all(axis=1)].tolist()
    xy = np.array(mol.GetConformer().GetPositions())[:, :2]; xy[:, 1] *= -1
    span = np.ptp(xy, axis=0)
    scale = size * .82 / max(float(max(span)), 1)
    xy = (xy - (xy.min(0) + xy.max(0)) / 2) * scale + size / 2
    lengths = [np.linalg.norm(xy[b.GetBeginAtomIdx()] - xy[b.GetEndAtomIdx()]) for b in mol.GetBonds()]
    bond = float(np.median(lengths)); font = min(34, max(12, bond * .50))
    components = list(dict.fromkeys(o["component"] for o in audit["atom_origin"]))
    colors = resolved_palette(options)
    colored = style_id in ("red_blue","bold_color","stacked_eclipsed","stacked_slipped","perspective","foreground_pore")
    style = {"font_family": options["font_family"], "atom_font_size": font*options["font_scale"],
             "bond_width": max(1, size / 800)*options["bond_scale"]*(1.6 if style_id=="bold_color" else 1),
             "double_offset": bond * .065, "ring_double_offset": bond * .16,
             "font_weight":"bold" if style_id=="bold_color" else "normal",
             "port_style": options["edge_style"],
             "component_colors": {c: colors[i % len(colors)] if colored else "#171717" for i, c in enumerate(components)}}
    def key(i):
        o = audit["atom_origin"][i]; return f"{o['instance']}:{o['fragment_atom_index']}"
    if style_id == "hetero_colors" or options["element_scheme"] in ELEMENT_SCHEMES:
        atom_palette=ELEMENT_SCHEMES.get(options["element_scheme"], ELEMENT_SCHEMES["conventional"])
        style["atom_colors"]={key(a.GetIdx()):atom_palette[a.GetAtomicNum()] for a in mol.GetAtoms() if a.GetAtomicNum() in atom_palette}
    if style_id=="foreground_pore" and recipe.get("repeat_faces"):
        foreground=set(recipe["repeat_faces"][0])
        ink=style.setdefault("atom_colors",{})
        for i,origin in enumerate(audit["atom_origin"]):
            if origin["instance"] not in foreground:
                color=ink.get(key(i),style["component_colors"][origin["component"]])
                ink[key(i)]='#'+''.join(f'{round(255*.72+int(color[j:j+2],16)*.28):02x}' for j in (1,3,5))
    if options["abbreviation_mode"] != "explicit":
        groups, used = [], set()
        patterns = [("[OX2]-[CH2]-[CH3]", "OEt"), ("[OX2]-[CH3]", "OMe" if options["abbreviation_mode"]=="OMe" else "OCH3"), ("[SX4](=[OX1])(=[OX1])-[OX2H]", "SO3H")]
        if options['abbreviation_mode'] == 'source_aldehyde':
            patterns = [("[CH1](=[OX1])", 'CHO')]
        for pattern, label in patterns:
            for hit in mol.GetSubstructMatches(Chem.MolFromSmarts(pattern)):
                if used.intersection(hit) or mol.GetAtomWithIdx(hit[0]).GetDegree() < 2:
                    continue
                parent = next((a.GetIdx() for a in mol.GetAtomWithIdx(hit[0]).GetNeighbors() if a.GetIdx() not in hit), None)
                if parent is None:
                    continue
                # Place the chemically attached S/O glyph next to the incoming
                # bond, including when the substituent points left.
                displacement = xy[hit[0]]-xy[parent]
                left = displacement[0] < -abs(displacement[1])*.25
                displayed = {"SO3H":"HO3S", "OCH3":"H3CO", "OMe":"MeO", "OEt":"EtO", 'CHO':'OHC'}.get(label, label) if left else label
                shift = (-1 if left else 1) * (len(displayed)-1) * style["atom_font_size"] * .56 / 2
                label_xy = [float(xy[hit[0], 0] + shift), float(xy[hit[0], 1])]
                groups.append({"anchor": key(hit[0]), "label": displayed, "label_xy": label_xy,
                               "attachment_atom": mol.GetAtomWithIdx(hit[0]).GetSymbol(),
                               "attachment_side": "right" if left else "left",
                               "hidden_atoms": [key(i) for i in hit[1:]]}); used.update(hit)
        style["collapsed_groups"] = groups
    stem = directory / recipe["recipe_id"]
    chemical_svg, xy, model_pixels = _native_styled_svg(mol, audit, recipe, style_id, options, style, size, model_coordinates)
    svg=decorate_svg(chemical_svg,mol,audit,recipe,xy,style_id,options,size)
    figure_atoms = mol.GetNumAtoms()
    # Geometry editing, pore measurements and all chemical files use the original
    # connected model. Project those atoms with the same drawer used for the figure.
    mol, audit, xy = model_mol, model_audit, model_pixels
    stem.with_suffix(".svg").write_text(svg,encoding="utf-8")
    import resvg_py
    stem.with_suffix(".png").write_bytes(resvg_py.svg_to_bytes(svg_string=svg, background=None if options["transparent_background"] else "white"))
    descriptor=next(s for s in STYLES if s["id"]==style_id)
    audit["drawing_style"] = {"id": style_id, "seed": seed, "rules": style, "options":options,
                              "category":descriptor["category"],"representation_scope":descriptor["representation_scope"],
                              "source_registered": bool(recipe.get('source_style',{}).get('source_registered')) if style_id=='source' else False, "chemistry_changed": False,
                              "graph_exports":"Underlying single finite chemical graph; visual layers and cartoon symbols are presentation only."}
    note = ''
    if shortened:
        location = 'aryl rings' if boundary_mode == 'aryl' else 'aromatic hubs'
        note = f"PNG and SVG cut at the {location}. Chemical files retain the full terminal groups."
    elif boundary_mode != 'linkage':
        note = "This patch has no eligible terminal groups for the selected cut."
    audit["drawing_style"]["boundary"] = {"mode": boundary_mode, "shortened_sites": shortened,
        "figure_atoms": figure_atoms, "model_atoms": mol.GetNumAtoms(), "note": note}
    original_xy = np.asarray(mol.GetConformer().GetPositions())[:, :2]
    atom_keys = [f"{o['instance']}:{o['fragment_atom_index']}" for o in audit['atom_origin']]
    names = list(dict.fromkeys(o['instance'] for o in audit['atom_origin']))
    centers = {name: xy[[i for i,o in enumerate(audit['atom_origin']) if o['instance']==name]].mean(0).tolist() for name in names}
    distances = [np.linalg.norm(original_xy[b.GetBeginAtomIdx()]-original_xy[b.GetEndAtomIdx()]) for b in mol.GetBonds()]
    pixels = [np.linalg.norm(xy[b.GetBeginAtomIdx()]-xy[b.GetEndAtomIdx()]) for b in mol.GetBonds()]
    audit['drawing_style']['renderer'] = 'shared-rdkit-chemical-svg'
    audit['drawing_style']['layout'] = {'canvas_size': [size,size], 'atom_keys': atom_keys,
        'atom_coordinates': original_xy.tolist(), 'atom_pixels': xy.tolist(),
        'atom_normalized': (xy/size).tolist(), 'instance_centers_pixels': centers,
        'instance_centers_normalized': {name: (np.asarray(p)/size).tolist() for name,p in centers.items()},
        'pixels_per_unit': float(np.median(np.asarray(pixels)/np.asarray(distances))),
        'coordinate_system': 'cartesian_y_up',
        'display_projection': descriptor['category'] == 'layered' and not descriptor.get('flat_foreground', False),
        'editing_scope': 'Atom pixels describe the unprojected chemical layer; use flat chemical styles for coordinate editing.'}
    if editing_indices is not None:
        audit['drawing_style']['layout']['editing_atom_indices'] = editing_indices
    return audit


def export_chemical_formats(directory, stem, *, source_mol=None, repeated_patch=False):
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
                if any(a.GetAtomicNum() == 7 and a.GetIsAromatic() and a.GetTotalNumHs() > 0
                       for a in candidate.GetAtoms()):
                    # Native CDXML writing of the direct Kekule copy can add
                    # unknown-double marks inside a porphyrin ring system.
                    # The established MOL parser retains the stored ring/stereo
                    # state; keep full atom/H/bond/stereo validation afterward.
                    raise ValueError("aromatic N-H graphs retain the scoped MOL parser")
                if (candidate.GetStereoGroups()
                        or any(a.GetChiralTag() != Chem.ChiralType.CHI_UNSPECIFIED for a in candidate.GetAtoms())
                        or any(b.GetStereo() not in (Chem.BondStereo.STEREONONE, Chem.BondStereo.STEREOANY)
                               or b.GetBondDir() not in (Chem.BondDir.NONE, Chem.BondDir.EITHERDOUBLE)
                               for b in candidate.GetBonds())):
                    raise ValueError("defined atom/bond stereo retains the scoped MOL parser")
                lines = molblock.splitlines()
                if len(lines)>3 and "V3000" in lines[3] and repeated_patch:
                    # Large MOL files switch to V3000. Re-perceiving every
                    # macrocycle in a multi-pore graph can take minutes even
                    # in the unsanitized native reader. Use the already bound
                    # graph and verify the writer's atom/bond rows directly.
                    atom_rows=lines[lines.index('M  V30 BEGIN ATOM')+1:lines.index('M  V30 END ATOM')]
                    bond_rows=lines[lines.index('M  V30 BEGIN BOND')+1:lines.index('M  V30 END BOND')]
                    if (len(atom_rows),len(bond_rows))!=(candidate.GetNumAtoms(),candidate.GetNumBonds()):
                        raise ValueError('V3000 atom or bond count differs')
                    Chem.Kekulize(candidate,clearAromaticFlags=True)
                    for atom,row in zip(candidate.GetAtoms(),atom_rows):
                        index,symbol,x,y,z,*_ = row.split()[2:]
                        if int(index)!=atom.GetIdx()+1 or symbol not in (atom.GetSymbol(), 'R' if atom.GetAtomicNum()==0 else atom.GetSymbol()):
                            raise ValueError('V3000 atom identity differs')
                        candidate.GetConformer().SetAtomPosition(atom.GetIdx(),(float(x),float(y),float(z)))
                    for bond,row in zip(candidate.GetBonds(),bond_rows):
                        idx,order,a,b,*flags=row.split()[2:]
                        if tuple(map(int,(idx,order,a,b))) != (bond.GetIdx()+1,int(bond.GetBondTypeAsDouble()),
                                                              bond.GetBeginAtomIdx()+1,bond.GetEndAtomIdx()+1):
                            raise ValueError('V3000 bond identity differs')
                        if flags==['CFG=2'] and order=='2':
                            bond.SetStereo(Chem.BondStereo.STEREOANY)
                            bond.SetBondDir(Chem.BondDir.EITHERDOUBLE)
                        elif flags:
                            raise ValueError('V3000 bond annotation requires independent verification')
                    candidate.UpdatePropertyCache(strict=True)
                    report.update(molecule_input='verified_source_mol_v3000',source_mol_serialization_matches_stored_MOL=True)
                    return candidate
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
                if repeated_patch:
                    raise ValueError('CDXML could not be verified from the repeated graph. MOL and SDF retain the displayed chemistry.') from exc
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
