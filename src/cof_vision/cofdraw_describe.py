"""Deterministic, catalog-bound prose planner; never an unrestricted chemist.

Only registered source IDs or existing monomer SMILES can reach a ready plan.
Unknown terms and explicit constraints survive optional UI/LLM ID selection.
"""
from __future__ import annotations

from copy import deepcopy
import re
import unicodedata

from .cofdraw_chemistry import prepare_template

SELECTION_KEYS = {"source_recipe_id", "monomer_a_id", "monomer_b_id", "linkage", "template_id", "style_id", "palette"}
SCOPE = "Finite catalog-bound COF drawing plan for review; not a synthesis prediction, automatic topology assignment, periodic crystal, or literature verification."
SOURCE_ALIASES = {
    "pyurea-source-layout-v1": ["Pyurea", "Pyurea COF"],
    "f-cof-1-source-v1": ["FCOF-1", "FCOF1", "F COF 1"],
    "cof-300-projection-v1": ["COF-300"],
}
# Names with verified precursor pairs. These create qualified finite designs.
NAMED_DESIGNS = {
    "COF-LZU1": ("TFB + Pa imine", "10.1021/ja206846p"),
    "Py-1P": ("TA + PyTTA imine", "10.1038/s41467-022-29086-x"),
    "COF-42": ("TFB + DETH acylhydrazone", "10.1021/ja204728y"),
    "COF-43": ("TFPB + DETH acylhydrazone", "10.1021/ja204728y"),
    "COF-76": ("TFPPy + TAPA imine", "10.1021/jacs.9b13971"),
    "4PE-1P": ("TA + TAPE imine kgm_4_2", "10.1021/ja5092936"),
    "Py-Azine": ("TFPPy + hydrazine azine", "10.1021/ja4103293"),
    "H2Por-DETH-COF": ("Por-CHO + DETH acylhydrazone", "10.1038/s41467-021-21527-3"),
}
NAMED_PAIRS = {
    "COF-LZU1": ("tfb","pa","imine","hcb_3_2"),
    "Py-1P": ("ta","pytta","imine","sql_4_2"),
    "COF-42": ("tfb","deth","acylhydrazone","hcb_3_2"),
    "COF-43": ("tfpb","deth","acylhydrazone","hcb_3_2"),
    "COF-76": ("tfppy","tapa","imine","alternating_hexagon"),
    "4PE-1P": ("ta","tape","imine","kgm_4_2"),
    "Py-Azine": ("tfppy","hydrazine","azine","sql_4_2"),
    "H2Por-DETH-COF": ("por_cho","deth","acylhydrazone","sql_4_2"),
}
CORE_NAMES = {
    "tapp": ["porphyrin tetraamine", "aminophenyl porphyrin", "porphyrin node"],
    "por_cho": ["porphyrin tetraaldehyde", "formylphenyl porphyrin", "porphyrin node"],
    "pytta": ["pyrene tetraamine", "aminophenyl pyrene", "pyrene node"],
    "tfppy": ["pyrene tetraaldehyde", "formylphenyl pyrene", "pyrene node"],
    "tfppy_13c": ["13C pyrene tetraaldehyde", "13C-labelled pyrene tetraaldehyde"],
}
def expand_material_names(text):
    notes = []
    for name, (building_blocks, doi) in NAMED_DESIGNS.items():
        pattern = r"(?<![\w-])" + re.escape(name).replace(r"\-", r"[- ]?") + r"(?![\w-])"
        text, count = re.subn(pattern, building_blocks, text, flags=re.I)
        if count:
            notes.append(f"{name} uses the published precursors ({building_blocks}). This is a proposed finite drawing. DOI {doi}.")
    return text, notes


LINKAGE_ALIASES = {
    "imine": ["imine", "C=N"],
    "beta_ketoenamine": ["beta ketoenamine", "beta-ketoenamine", "keto-enamine", "ketoenamine", "beta keto enamine"],
    "acylhydrazone": ["acylhydrazone", "acyl hydrazone"],
    "nitrone": ["nitrone", "nitrone linkage"],
    "keto_hydrazide": ["keto hydrazide", "keto-enamine-hydrazide"],
    "azine": ["azine"], "secondary_amine": ["reduced imine", "reduced-imine", "secondary amine"],
    "boronate_ester": ["boronate", "boronate ester"],
    "cyano_vinylene_sp2_carbon": ["cyano-vinylene", "cyano vinylene"],
    "vinylene_review_routes": ["vinylene", "sp2 carbon", "sp2-carbon"],
    "urea_containing_cof": ["urea linkage", "urea-linked", "urea linked"],
}
STYLE_ALIASES = {
    "gradient_stack": ["gradient layers", "fading layers", "stacked gradient"],
    "pore_perspective": ["through the pore", "pore perspective", "look down the channel"],
    "half_space_stack": ["half space filling", "half-space-filling layers", "half-space filling"],
    "element_compact": ["compact hydrogen labels", "compact heteroatom labels"],
    "sphere_cutaway": ["space filling cutaway", "space-filling cutaway", "spheres over skeletal"],
    "standard": ["standard", "standard chemical drawing", "black bonds", "plain chemical drawing"],
    "source": ["source style", "original colors", "original colours", "source colors", "source colours"],
    "mono_dashed": ["dashed edges", "dashed bonds"], "mono_wave": ["wavy edges", "wavy bonds"],
    "red_blue": ["colored building blocks", "coloured building blocks", "colored bonds", "coloured bonds", "red and blue", "red-blue"],
    "serif": ["serif", "serif labels"], "compact_labels": ["condensed labels", "condensed substituents"],
    "hetero_colors": ["colored heteroatoms", "coloured heteroatoms", "heteroatom colors", "color by element", "element colors", "element colours"],
    "node_highlight": ["highlights", "highlight", "highlighted building blocks", "building-block highlights", "building block highlights"],
    "ring_fill": ["filled aromatic rings", "ring fill"], "pore_fill": ["filled pore", "pore fill"],
    "pastel_blocks": ["pastel blocks", "pastel block regions", "soft block regions"],
    "bold_color": ["bold colors", "bold colours", "bold colored structure"],
    "stacked_eclipsed": ["stacked layers", "eclipsed layers", "transparent layers"],
    "stacked_slipped": ["offset layers", "slipped layers"], "perspective": ["tilted sheet", "perspective", "oblique projection"],
    "cartoon_nodes": ["cartoon", "cartoon nodes", "nodes and struts"],
    "cartoon_overlay": ["cartoon overlay", "chemical cartoon overlay"],
    "glow_outline": ["glow", "molecular glow"], "molecular_spheres": ["molecular spheres", "ball and stick illustration"],
}
PALETTE_ALIASES = {"pastel": ["pastel", "soft palette", "soft colors", "soft colours"],
                   "ink": ["monochrome", "black and white"], "colorblind": ["colorblind", "colorblind friendly"]}
GEOMETRIES = {"hexagonal": ["hexagonal", "hexagon", "hcb", "honeycomb lattice", "honeycomb"], "square": ["square", "sql", "square lattice", "square grid"],
              "triangular": ["triangular", "hxl"], "chain": ["open chain", "chain"],
              "custom": ["custom six-unit ring", "custom ring"]}
GEOMETRY_TEMPLATES = {"hexagonal": {"hcb_3_3", "hcb_3_2"}, "square": {"sql_4_2", "sql_4_4"},
                      "triangular": {"hxl_6_2"}, "chain": {"chain_2_2"}, "custom": {"alternating_hexagon"}}
FILLER = set("please draw drawing depict render make create show me a an the this that my our some cof cofs framework pore patch finite chemical structure design designed new from using use with and plus of for in on as into to it its is be would like want need can could you i we style styles palette colors colours color colour colored coloured soft pretty beautiful clean publication ready high quality nice two monomers monomer building blocks block linkage linked bonds bond geometry topology template layout source catalog registered original look illustration illustrative diagram model highlight region regions by select choose selected chemistry supported units alternating 2d review proposal proposed only".split())
FILLER.update("open load get retrieve display redraw sketch visualize visualise generate prepare export save download give produce add apply drawn drawing shown opening loading rendering figure image picture diagram schematic version please also then just kindly".split())
FILLER.update("png svg mol sdf cdxml format file files output outputs".split())
COUNT_WORDS = {"di": 2, "tri": 3, "tetra": 4, "hexa": 6, "two": 2, "three": 3, "four": 4, "six": 6,
               "eight": 8, "twelve": 12}
ROLE_GROUPS = {"amine": "primary_arylamine", "aldehyde": "aryl_aldehyde", "hydrazide": "acylhydrazide",
               "hydroxylamine": "aryl_hydroxylamine"}


def _normal(text):
    value = unicodedata.normalize("NFKC", text)
    return re.sub(r"\s+", " ", value.translate(str.maketrans({"–": "-", "—": "-", "−": "-", "β": "beta", "Β": "beta", "′": "'", "’": "'"}))).strip()


def _unique(values):
    result = []
    seen = set()
    for value in values:
        if isinstance(value, str) and value.strip() and _normal(value).casefold() not in seen:
            result.append(value.strip()); seen.add(_normal(value).casefold())
    return result


def _aliases(row, fields):
    values = [row.get(key) for key in fields] + list(row.get("aliases") or [])
    for field in fields:
        value = row.get(field)
        if not isinstance(value, str):
            continue
        # A label suffix such as "local boronate fragment" or "projected
        # network" qualifies a source; it is not another material name.
        name = re.split(r"\s+[·—]\s+", value, maxsplit=1)[0]
        values.extend(re.split(r"\s+/\s+", name))
        trailing = re.search(r"\s+\(([^()]*)\)$", value)
        if trailing:
            values.extend([value[:trailing.start()], trailing[1]])
    return _unique(values)


def describe_catalog(library):
    """Compact allowlists for UI/LLM selection; no executable recipes or SMILES."""
    if not isinstance(library, dict):
        raise ValueError("Supply the current drawing library object")
    from .framework_drawing import canonical
    monomers = []
    smiles_ids = {}
    for row in library.get("monomer_presets", []):
        identifier = row.get("id")
        if not identifier:
            continue
        aliases = _aliases(row, ("id", "abbreviation", "full_name", "label")) + CORE_NAMES.get(identifier, []) + list(row.get("literature_names") or [])
        if "_" in identifier:
            aliases.extend([identifier.replace("_", "-"), identifier.replace("_", " ")])
        monomers.append({"id": identifier, "label": row.get("label", identifier),
                         "abbreviation": row.get("abbreviation", identifier), "full_name": row.get("full_name", ""),
                         "aliases": _unique(aliases), "compatible_roles": deepcopy(row.get("binary_roles", row.get("compatible_roles", [])))
                         if row.get("automatic_prepare_supported") is not False else []})
        if row.get("smiles"):
            smiles_ids.setdefault(canonical(row["smiles"]), []).append(identifier)
    by_id = {row["id"]: row for row in monomers}
    sources = []
    for row in library.get("source_recipes", []):
        identifier = row["recipe_id"]
        aliases = _aliases(row, ("recipe_id", "label", "literature_name", "display_name")) + SOURCE_ALIASES.get(identifier, [])
        aliases += [re.sub(r"(?<=[A-Za-z])-(?=\d)", sep, a) for a in list(aliases) for sep in (" ", "")]
        monomer_ids = []
        for precursor in row.get("monomers", []):
            ids = ([precursor["monomer_preset_id"]] if precursor.get("monomer_preset_id") in by_id else
                   smiles_ids.get(canonical(precursor["smiles"]), []) if precursor.get("smiles") else [])
            for monomer_id in ids:
                monomer_ids.append(monomer_id)
                extras = _aliases(precursor, ("label", "abbreviation", "full_name")) + list(precursor.get("literature_names") or [])
                for match in precursor.get("master_component_matches", []):
                    extras += [match.get("name"), match.get("abbrev")]
                by_id[monomer_id]["aliases"] = _unique(by_id[monomer_id]["aliases"] + extras)
        qualifications = list(row.get("qualifications", []))
        for field in ("identity_qualification", "scope"):
            if row.get(field):
                qualifications.append(row[field])
        sources.append({"id": identifier, "recipe_id": identifier, "label": row.get("label", identifier),
                        "aliases": _unique(aliases), "qualifications": _unique(qualifications),
                        "scope": row.get("scope", "Source-bound provisional finite patch; not a complete crystal."),
                        "source_registered": bool(row.get("source_registered", False)),
                        "sha256": row.get("source_recipe_sha256") or row.get("sha256"),
                        "monomer_ids": list(dict.fromkeys(monomer_ids))})
    linkages = []
    for identifier, row in library.get("linkages", {}).items():
        aliases = LINKAGE_ALIASES.get(identifier, []) + _aliases(row, ("label",)) + [identifier, identifier.replace("_", " ")]
        # "Azine / paper-labelled hydrazine" is a reference record; don't make
        # its broad alias collide with the explicit supported azine transform.
        if identifier == "azine_paper_labelled_hydrazine_link":
            aliases = [identifier, row.get("label", identifier)]
        linkages.append({"id": identifier, "label": row.get("label", identifier), "aliases": _unique(aliases),
                         "supported": row.get("supported") is True, "roles": list(row.get("roles", []))})
    templates = [{**{k: deepcopy(row.get(k)) for k in ("id", "label", "allow_pairs", "instances", "topology_label", "description")},
                  "aliases": [row["id"], row["id"].replace("_", " ")]} for row in library.get("templates", [])]
    styles = [{"id": row["id"], "label": row.get("label", row["id"]), "category": row.get("category"),
               "description": row.get("description"), "aliases": _unique([row["id"], row.get("label")] + STYLE_ALIASES.get(row["id"], []))}
              for row in library.get("styles", [])]
    palettes = [{"id": value, "label": value.replace("_", " ").title(),
                 "aliases": _unique([value, value.replace("_", " "), value.replace("_", "-")] + PALETTE_ALIASES.get(value, []))}
                for value in library.get("style_controls", {}).get("palettes", [])]
    demos = [
        {"id": "named_py", "label": "Py-1P by name", "text": "Draw Py-1P with gradient layers", "expected_status": "ready"},
        {"id": "porphyrin", "label": "Porphyrin node", "text": "Draw a porphyrin tetraamine with terephthalaldehyde in a square lattice", "expected_status": "ready"},
        {"id": "imine", "label": "Imine design", "text": "Draw TFB and TAPB in pastel colors", "expected_status": "ready"},
        {"id": "tppa", "label": "Catalog TpPa", "text": "Draw TpPa in standard style", "expected_status": "ready"},
        {"id": "pyurea", "label": "Pyurea highlights", "text": "Draw Pyurea with building-block highlights and a soft palette", "expected_status": "ready"},
        {"id": "fcof", "label": "F-COF-1 highlights", "text": "Draw FCOF-1 with highlights and ocean coral colors", "expected_status": "ready"},
        {"id": "source_style", "label": "Recorded source colors", "text": "Draw F-COF-1 in source style", "expected_status": "ready"},
        {"id": "hydrazone", "label": "Acylhydrazone design", "text": "Draw a hexagonal acylhydrazone COF from TFB and DETH with teal orange highlights", "expected_status": "ready"},
        {"id": "square", "label": "Square design", "text": "Draw TFPPy and Pa with imine linkage in a square template and lavender mint colors", "expected_status": "ready"},
        {"id": "unsupported", "label": "Unsupported transform", "text": "Draw a boronate ester COF from BDBA and HHTP", "expected_status": "unsupported"},
        {"id": "choose", "label": "Choose a material", "text": "Draw a COF", "expected_status": "needs_input"},
    ]
    return {"sources": sources, "monomers": monomers, "linkages": linkages, "templates": templates,
            "named_designs": [{"name": name, "precursors": pair, "doi": doi} for name, (pair, doi) in NAMED_DESIGNS.items()],
            "styles": styles, "palettes": palettes, "demos": demos, "scope": SCOPE}


def _find(text, rows, kind, blocked=()):
    matches = []
    deliberate_short = bool(re.search(r"\bmonomers?\b|\+", text, re.I))
    for row in rows:
        for alias in row["aliases"]:
            value = _normal(alias)
            if not value or (kind == "monomer" and value.isdecimal()):
                continue
            for found in re.finditer(r"(?<![\w])" + re.escape(value) + r"(?![\w])", text, re.I):
                span = found.span()
                if any(span[0] < b and span[1] > a for a, b in blocked):
                    continue
                if kind == "monomer" and len(re.sub(r"\W", "", value)) <= 2:
                    # Pa/TA/Ha are also ordinary words, pressure symbols, etc.
                    abbreviation = row.get("abbreviation", "")
                    if found[0] != abbreviation and not deliberate_short:
                        continue
                matches.append({"kind": kind, "id": row["id"], "label": row["label"],
                                "matched_text": found[0], "span": list(span)})
    # Longest aliases own their span, but ties for different IDs stay ambiguous.
    retained = []
    for match in sorted(matches, key=lambda x: (-(x["span"][1]-x["span"][0]), x["span"][0], x["id"])):
        a, b = match["span"]
        if any(a >= x["span"][0] and b <= x["span"][1] and
               (match["span"] != x["span"] or match["id"] == x["id"]) for x in retained):
            continue
        retained.append(match)
    return sorted(retained, key=lambda x: (x["span"][0], x["id"]))


def _mask(text, matches):
    characters = list(text)
    for match in matches:
        a, b = match["span"]
        characters[a:b] = " " * (b-a)
    return "".join(characters)


def _choices(matches):
    return list(dict.fromkeys(match["id"] for match in matches))


def _explicit_requirements(text):
    """Recognize functional counts, never point-group labels such as C3/C6."""
    roles = r"(?P<role>amine|aldehyde|hydrazide|hydroxylamine)s?"
    modifier = r"(?:(?:aromatic|aryl|primary)\s+)*"
    patterns = [
        r"\b(?P<count>di|tri|tetra|hexa)[- ]?(?:topic|functional)\s+" + modifier + roles + r"\b",
        r"\b(?P<count>\d{1,2}|two|three|four|six)[- ](?:connected|functional|site)\s+" + modifier + roles + r"\b",
        r"\b(?P<count>di|tri|tetra|hexa)" + roles + r"\b",
        r"\b" + roles + r"\s+with\s+(?P<count>\d{1,2}|two|three|four|six)\s+(?:reactive\s+)?(?:sites|groups|connections)\b",
    ]
    functionality = []
    matches = []
    for pattern in patterns:
        for found in re.finditer(pattern, text, re.I):
            value = found["count"].casefold()
            count = COUNT_WORDS[value] if value in COUNT_WORDS else int(value)
            role = found["role"].casefold()
            if (role, count) not in functionality:
                functionality.append((role, count))
            matches.append({"kind": "requirement", "id": f"{role}_{count}_sites", "label": found[0],
                            "matched_text": found[0], "span": list(found.span())})
    instances = []
    pattern = r"\b(?P<count>\d{1,3}|two|three|four|six|eight|twelve)[- ](?:instances?|units?|building[- ]blocks?)\b"
    for found in re.finditer(pattern, text, re.I):
        value = found["count"].casefold()
        count = COUNT_WORDS[value] if value in COUNT_WORDS else int(value)
        if count not in instances:
            instances.append(count)
        matches.append({"kind": "requirement", "id": f"{count}_instances", "label": found[0],
                        "matched_text": found[0], "span": list(found.span())})
    return functionality, instances, matches


def _check_functionality(requirements, monomers, unresolved):
    for role, count in requirements:
        actual = sorted({int(group["count"]) for row in monomers for group in row.get("functionality", [])
                         if group.get("kind") == ROLE_GROUPS[role]})
        if count not in actual:
            unresolved.append(f"The explicit {count}-site {role} requirement does not match the selected catalog monomers (recognized counts: {actual or 'none'}).")


def describe_cof(text, library, selection=None, *, allow_unmatched_words=False):
    """Resolve supported prose and selected catalog IDs without changing intent.

    A selection can fill gaps or disambiguate a matched alias; it cannot replace
    an explicitly named source/monomer/linkage/template. Appearance overrides
    are disclosed. Explicitly reviewed or model-proposed selections may request
    residual wording as warnings; this is not a semantic-fidelity guarantee.
    """
    if not isinstance(text, str) or not text.strip() or len(text) > 4000:
        raise ValueError("Describe the drawing in 1-4000 characters")
    if selection is not None and not isinstance(selection, dict):
        raise ValueError("selection must be an object of catalog IDs")
    if not isinstance(allow_unmatched_words, bool):
        raise ValueError("allow_unmatched_words must be a boolean")
    catalog = describe_catalog(library)
    normalized, name_notes = expand_material_names(_normal(text))
    supplied = {key: value for key, value in (selection or {}).items() if value is not None}
    result = {"status": "needs_input", "summary": "Choose a source drawing or two compatible catalog monomers.",
              "selection": {key: None for key in sorted(SELECTION_KEYS)}, "prepare_request": None,
              "source_recipe_id": None, "style_id": "standard", "style_options": {}, "matches": [],
              "warnings": [], "unresolved": [], "alternatives": [], "assumptions": name_notes, "scope": SCOPE}
    unresolved, warnings, assumptions = result["unresolved"], result["warnings"], result["assumptions"]
    invalid = sorted(set(supplied) - SELECTION_KEYS)
    if invalid:
        unresolved.append("Unsupported selection fields: " + ", ".join(invalid))
    groups = {"source_recipe_id": "sources", "monomer_a_id": "monomers", "monomer_b_id": "monomers",
              "linkage": "linkages", "template_id": "templates", "style_id": "styles", "palette": "palettes"}
    for key, value in list(supplied.items()):
        if key not in groups:
            supplied.pop(key); continue
        allowed = {row["id"] for row in catalog[groups[key]]}
        if not isinstance(value, str) or value not in allowed:
            unresolved.append(f"Unknown catalog ID for {key}: {str(value)[:100]}")
            supplied.pop(key)
    if supplied.get("source_recipe_id") and any(supplied.get(k) for k in ("monomer_a_id", "monomer_b_id", "linkage", "template_id")):
        unresolved.append("Choose a registered source or a new monomer design; a source selection cannot also change its chemistry/template.")
    source_matches = _find(normalized, catalog["sources"], "source")
    blocked = [tuple(m["span"]) for m in source_matches]
    contextual_linkages = [match for match in _find(normalized, catalog["linkages"], "linkage", blocked)
                           if re.match(r"\s*-?\s*(?:linked|linkage|bonds?)\b", normalized[match["span"][1]:], re.I)]
    monomer_matches = _find(normalized, catalog["monomers"], "monomer", blocked + [tuple(m["span"]) for m in contextual_linkages])
    blocked += [tuple(m["span"]) for m in monomer_matches]
    linkage_matches = _find(normalized, catalog["linkages"], "linkage", blocked)
    template_matches = _find(normalized, catalog["templates"], "template", blocked)
    style_matches = _find(normalized, catalog["styles"], "style", blocked)
    # Palette words inside a complete style phrase (e.g. pastel block regions)
    # describe that style; its renderer default remains intact unless separate.
    palette_matches = _find(normalized, catalog["palettes"], "palette", blocked + [tuple(m["span"]) for m in style_matches])
    geometry_rows = [{"id": key, "label": key, "aliases": aliases} for key, aliases in GEOMETRIES.items()]
    geometry_matches = _find(normalized, geometry_rows, "geometry", blocked + [tuple(m["span"]) for m in template_matches])
    functionality, instance_counts, requirement_matches = _explicit_requirements(normalized)
    result["matches"] = all_matches = source_matches + monomer_matches + linkage_matches + template_matches + style_matches + palette_matches + geometry_matches + requirement_matches
    residual = _mask(normalized, all_matches)
    negated = re.search(r"\b(?:not|no|without|except|instead|rather than|avoid|don't|do not|never|exclude)\b", normalized, re.I)
    if negated:
        unresolved.append("Negation or replacement needs an explicit positive drawing choice; the planner will not ignore '" + negated[0] + "'.")
    crystal = re.search(r"\b(?:3[- ]?d|(?:three|3)[- ]dimensional|crystal|crystalline|cif|periodic|optimi[sz](?:ed|ation|e)|relaxed|DFT|BET|bandgap|band gap|synthesi[sz]e|synthesis|calculate|compute|PXRD|diffraction|surface area|pore volume|pore size|density)\b", normalized, re.I)
    coordinate_export = re.search(r"\b(?:XYZ|PDB)\b", normalized, re.I)
    physical_unit = r"(?:nm|nanomet(?:er|re)s?|angstroms?|Å|A|pm|picomet(?:er|re)s?)"
    physical_target = re.search(r"\b\d+(?:\.\d+)?[- ]*" + physical_unit + r"[- ]+(?:wide[- ]+)?(?:pores?|channels?|diameter|radius)\b", normalized, re.I)
    physical_target = physical_target or re.search(r"\b(?:pores?|channels?|diameter|radius)\b.{0,30}?\b\d+(?:\.\d+)?[- ]*" + physical_unit + r"\b", normalized, re.I)
    unsupported = bool(crystal or coordinate_export or physical_target)
    if crystal or coordinate_export:
        unresolved.append("Prose drawing cannot supply a validated crystal, 3D atomic structure, synthesis, or calculated material property. Use the structure laboratory with explicit inputs.")
    if physical_target:
        unresolved.append("A requested physical pore dimension cannot be guaranteed by a finite drawing template. Use an explicit structure and the qualified geometry tools to evaluate a dimensional target.")
    if len(instance_counts) > 1:
        unresolved.append("Conflicting instance counts were requested; choose one finite template size.")
    unknown = _unique(word for word in re.findall(r"[\w][\w'-]*", residual) if word.casefold() not in FILLER)
    if unknown:
        if allow_unmatched_words and supplied:
            warnings.append("The selected parameters interpret wording that local rules could not match; verify: " + ", ".join(unknown[:20]) + ". This does not establish that the interpretation matches every requested name or property.")
        else:
            unresolved.append("Unresolved words or names: " + ", ".join(unknown[:20]) + ". Use a listed name or clarify the constraint.")

    def choose(key, matches, default=None):
        ids = _choices(matches)
        selected = supplied.get(key)
        appearance_override = key in {"style_id", "palette"} and selected
        if appearance_override and ids and (selected not in ids or len(ids) > 1):
            assumptions.append(f"The explicitly selected {key} ({selected}) overrides the appearance wording ({', '.join(ids)}).")
        elif selected and ids and selected not in ids:
            unresolved.append(f"Selected {key} conflicts with the explicit text ({', '.join(ids)}).")
        same_alias = len({tuple(match["span"]) for match in matches}) == 1
        if len(ids) > 1 and not appearance_override and not (selected in ids and same_alias):
            unresolved.append(f"More than one {key} matches; choose explicitly: {', '.join(ids)}.")
        return selected or (ids[0] if len(ids) == 1 else default)

    style = choose("style_id", style_matches, "standard")
    palette = choose("palette", palette_matches)
    # A color request must produce visible color, not an inert palette on black
    # standard bonds. Explicit standard style stays black with a review warning.
    if palette and style == "standard" and not style_matches and not supplied.get("style_id"):
        wanted = "pastel_blocks" if palette == "pastel" else "red_blue"
        style = wanted if any(row["id"] == wanted for row in catalog["styles"]) else "standard"
        assumptions.append("Use pastel block regions behind dark chemical bonds for legibility." if style == "pastel_blocks" else
                           "Use the colored-building-block style so the requested palette is visible.")
    elif palette and style in {"standard", "mono_wave", "mono_dashed", "serif"}:
        warnings.append("This chemical line style uses black bonds; the palette may have no visible effect. Choose colored blocks or highlights for color.")
    result.update(style_id=style, style_options={"palette": palette} if palette else {})
    result["selection"].update(style_id=style, palette=palette)
    style_row = next((r for r in catalog["styles"] if r["id"] == style), {})
    if style_row.get("category") in {"layered", "schematic"}:
        warnings.append("Selected layers/projection/cartoon is an illustration of the finite drawing, not a periodic atomic structure.")
    source_id = choose("source_recipe_id", source_matches)
    linkage = choose("linkage", linkage_matches)
    if linkage and not next(row for row in catalog["linkages"] if row["id"] == linkage)["supported"]:
        unresolved.append("The requested linkage is catalogued but has no automatic monomer construction rule.")
        unsupported = True
    template_id = choose("template_id", template_matches)
    geometries = _choices(geometry_matches)
    if len(geometries) > 1:
        unresolved.append("Several incompatible drawing geometries were requested; choose one template.")
    if template_id and geometries and template_id not in GEOMETRY_TEMPLATES[geometries[0]]:
        unresolved.append("Selected template conflicts with the explicitly requested geometry.")
    runtime_monomers = {row["id"]: row for row in library["monomer_presets"]}
    # Distinct same-span aliases (notably BD) cannot be resolved by order.
    ambiguous_spans = {}
    for match in monomer_matches:
        ambiguous_spans.setdefault(tuple(match["span"]), set()).add(match["id"])
    selected_monomers = [supplied.get(k) for k in ("monomer_a_id", "monomer_b_id") if supplied.get(k)]
    resolved_slots = []
    for span, ids in sorted(ambiguous_spans.items()):
        if len(ids) > 1:
            chosen = ids.intersection(selected_monomers)
            if len(chosen) == 1:
                resolved_slots.append(next(iter(chosen)))
            else:
                unresolved.append("Ambiguous monomer alias '" + normalized[slice(*span)] + "': " + ", ".join(sorted(ids)) + ".")
                if len(ambiguous_spans) <= 2:
                    field = ("monomer_a_id", "monomer_b_id")[len(resolved_slots)]
                    for identifier in sorted(ids):
                        row = runtime_monomers[identifier]
                        result["alternatives"].append({"kind": "monomer_choice", "field": field, "id": identifier,
                                                       "label": row.get("full_name") or row["label"]})
                resolved_slots.append(None)
        else:
            resolved_slots.append(next(iter(ids)))
    monomer_ids = _unique(identifier for identifier in resolved_slots if identifier)
    if source_id:
        source = next(r for r in catalog["sources"] if r["id"] == source_id)
        if any(supplied.get(k) for k in ("monomer_a_id", "monomer_b_id")) or template_id or geometries or linkage:
            unresolved.append("A source drawing preserves its recorded chemistry and incidence. Remove design/template constraints or choose a new design without the source identity.")
        if monomer_ids and not set(monomer_ids).issubset(source["monomer_ids"]):
            unresolved.append("The named monomers change or cannot be verified against this source; a changed design cannot keep its source identity.")
        _check_functionality(functionality, [runtime_monomers[key] for key in source["monomer_ids"]], unresolved)
        if instance_counts:
            unresolved.append("Source metadata does not expose a validated instance-count constraint. Load the unchanged source without a new instance count, or choose an explicit design template.")
        warnings.extend(source["qualifications"])
        result["source_metadata"] = deepcopy(source)
        result["selection"]["source_recipe_id"] = source_id
        result["summary"] = f"Registered source drawing: {source['label']}; {style_row.get('label', style)}."
        if not unresolved:
            result.update(status="ready", source_recipe_id=source_id)
        elif unsupported:
            result["status"] = "unsupported"
        return _finish(result, catalog)

    if selected_monomers:
        if len(selected_monomers) == 2 and not set(monomer_ids).issubset(selected_monomers):
            unresolved.append("Selected monomers conflict with one or more explicitly named monomers.")
        for identifier in selected_monomers:
            if identifier not in monomer_ids:
                monomer_ids.append(identifier)
    if len(monomer_ids) != 2:
        unresolved.append("Choose exactly two unambiguous existing monomers for a new design, or one registered source drawing.")
        result["selection"].update(linkage=linkage, template_id=template_id)
        slots = resolved_slots if None in resolved_slots else monomer_ids
        for key, identifier in zip(("monomer_a_id", "monomer_b_id"), slots):
            result["selection"][key] = identifier
        result["status"] = "unsupported" if unsupported else "needs_input"
        return _finish(result, catalog)
    a, b = [runtime_monomers[key] for key in monomer_ids]
    _check_functionality(functionality, (a, b), unresolved)
    supported = {row["id"]: row for row in catalog["linkages"] if row["supported"]}
    if linkage and linkage not in supported:
        unresolved.append("The requested linkage is catalogued but has no automatic monomer construction rule.")
        unsupported = True
    # Resolve roles from the actual decorated monomers, then validate the final
    # proposal through prepare_template (metadata alone never creates readiness).
    def orientations(link):
        roles = supported[link]["roles"]
        if len(roles) != 2:
            return []
        def valid(row, role):
            return any(item.get("linkage") == link and item.get("role") == role and item.get("reactive_sites", 0) >= 2 for item in row.get("binary_roles",row.get("compatible_roles", [])))
        return [(left, right) for left, right in ((a, b), (b, a)) if valid(left, roles[0]) and valid(right, roles[1])]
    options = [key for key in supported if key != "secondary_amine" and orientations(key)]
    if not linkage:
        preferences = {row.get("linkage") for row in library.get("cof_presets", [])
                       if {row.get("monomer_a_id"), row.get("monomer_b_id")} == set(monomer_ids)
                       and row.get("linkage") in options}
        if len(options) == 1:
            linkage = options[0]
            assumptions.append("Use the unique supported conventional linkage for these recognized reactive groups: " + linkage + ".")
        elif len(preferences) == 1:
            linkage = next(iter(preferences))
            assumptions.append("Use the existing catalog starting-pair preference " + linkage + "; this creates a new design and does not load that source material.")
        else:
            unresolved.append("Choose the linkage explicitly; compatible supported choices: " + (", ".join(options) or "none") + ".")
    result["selection"].update(monomer_a_id=a["id"], monomer_b_id=b["id"], linkage=linkage, template_id=template_id)
    if not linkage or linkage not in supported:
        result["status"] = "unsupported" if unsupported else "needs_input"
        return _finish(result, catalog)
    ordered = orientations(linkage)
    if not ordered:
        unresolved.append("The named monomers do not provide the required reactive groups for this linkage.")
        return _finish(result, catalog)
    a, b = ordered[0]
    if [a["id"], b["id"]] != monomer_ids:
        assumptions.append("Reorder the two monomers to match the transform's aldehyde/node and amine/linker roles.")
    roles = supported[linkage]["roles"]
    counts = [next(item["reactive_sites"] for item in row["compatible_roles"] if item["linkage"] == linkage and item["role"] == role) for row, role in zip((a, b), roles)]
    candidates = [row for row in catalog["templates"] if row.get("allow_pairs") and counts in row["allow_pairs"]]
    if geometries:
        candidates = [row for row in catalog["templates"] if row["id"] in GEOMETRY_TEMPLATES[geometries[0]] and (not row.get("allow_pairs") or counts in row["allow_pairs"])]
    if linkage!='imine' or not (a.get('kgm_dialdehyde') and b.get('kgm_etta')):
        candidates = [row for row in candidates if row['id']!='kgm_4_2']
    if instance_counts:
        candidates = [row for row in candidates if row["instances"] == instance_counts[0]]
    if not template_id:
        if len(candidates) == 1:
            template_id = candidates[0]["id"]
            assumptions.append("Choose the compatible finite drawing layout " + template_id + "; this is a layout proposal for review, not a prediction of periodic topology.")
        else:
            unresolved.append("Choose a compatible finite drawing template; monomer degree alone does not determine a unique topology.")
    result["selection"].update(monomer_a_id=a["id"], monomer_b_id=b["id"], linkage=linkage, template_id=template_id)
    if template_id:
        template = next(row for row in catalog["templates"] if row["id"] == template_id)
        if instance_counts and any(count != template["instances"] for count in instance_counts):
            unresolved.append(f"The requested instance count does not match {template_id}, which uses {template['instances']} units.")
        request = {"monomer_a": a["smiles"], "monomer_b": b["smiles"], "linkage": linkage,
                   "label_a": a.get("abbreviation", a["id"]), "label_b": b.get("abbreviation", b["id"]),
                   "instances": template["instances"], "template_id": template_id}
        try:
            prepared = prepare_template(**request)
            result["compatibility"] = {"reactive_site_counts": prepared["compatibility"]["reactive_site_counts"],
                                       "roles": prepared["compatibility"]["roles"], "validated_by": "prepare_template"}
        except (ValueError, TypeError, KeyError, RuntimeError) as exc:
            unresolved.append("Chemical/template validation failed: " + str(exc))
        else:
            if not unresolved:
                result.update(status="ready", prepare_request=request)
                result["summary"] = f"New finite {linkage} design from {request['label_a']} + {request['label_b']}, using {template['label']}."
    if unsupported:
        result["status"] = "unsupported"
        result["prepare_request"] = None
    warnings.append("This is a new designed patch. Matching starting monomers does not transfer a literature name, source identity, crystal structure, or experimental topology.")
    return _finish(result, catalog)


def _finish(result, catalog):
    result["warnings"] = _unique(result["warnings"])
    result["unresolved"] = _unique(result["unresolved"])
    if result["status"] != "ready":
        result["prepare_request"] = None
        result["source_recipe_id"] = None
        result["summary"] = ("This request includes unsupported operations or chemistry." if result["status"] == "unsupported" else "Review the unresolved choices before creating a drawing.")
        if not result["alternatives"]:
            result["alternatives"] = [{"kind": "demo", **row} for row in catalog["demos"]]
    return result


def validate_selection(selection, library, text="Draw a COF", *, allow_unmatched_words=False):
    """Same deterministic constraint and chemistry checks for a UI/LLM proposal."""
    return describe_cof(text, library, selection=selection, allow_unmatched_words=allow_unmatched_words)
