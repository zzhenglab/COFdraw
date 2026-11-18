"""Stateless, bounded conversation context for the catalog drawing planner.

The client retains its conversation in memory. Context contains IDs, a local
drawing token and unresolved wording only: never recipes, coordinates or SMILES.
"""
from copy import deepcopy
import json
import re

from .cofdraw_describe import (SELECTION_KEYS, _choices, _find, _mask, _normal,
                              describe_catalog, describe_cof, expand_material_names)

FIELDS = {"source_recipe_id": "sources", "monomer_a_id": "monomers", "monomer_b_id": "monomers",
          "linkage": "linkages", "template_id": "templates", "style_id": "styles", "palette": "palettes"}
CHEMISTRY = ("source_recipe_id", "monomer_a_id", "monomer_b_id", "linkage", "template_id")
FOLLOWUP_WORDS = r"\b(?:now|again|same|keep|continue|current|previous|drawing|please|change|switch|set)\b"
SUGGESTIONS = ["Make it pastel", "Use wavy bonds", "Change B to Pa", "How do I check SA scores?"]


def _reject_structure_payload(text):
    cif = re.search(r'(?im)^\s*(?:data_\S+|_atom_site_\S+|_cell_length_\S+|loop_)\s*', text)
    structured = text.lstrip().startswith(('{', '[')) and re.search(r'"(?:recipe|components|explicit_graph|atoms|monomers|format|schema_version|cell)"\s*:', text)
    if cif or structured:
        raise ValueError('Use Open project or the structure import tools for molecular files. Chat does not forward project JSON or CIF contents to a model provider.')


def checked_context(context, library):
    if context is None:
        return {"selection": {}, "pending": [], "requirements": [], "origin": "plan", "token": None}
    if not isinstance(context, dict) or set(context) - {"selection", "pending", "requirements", "origin", "token"}:
        raise ValueError("Conversation context accepts catalog selections and local state only.")
    result = {"selection": {}, "pending": [], "requirements": [], "origin": context.get("origin", "plan"), "token": context.get("token")}
    if result["origin"] not in {"plan", "current_drawing"}:
        raise ValueError("Unknown conversation context origin.")
    if result["token"] is not None and (not isinstance(result["token"], str) or not re.fullmatch(r"[A-Za-z0-9_-]{1,80}", result["token"])):
        raise ValueError("Invalid local drawing token.")
    selected = context.get("selection") or {}
    if not isinstance(selected, dict) or set(selected) - SELECTION_KEYS:
        raise ValueError("Context must contain catalog IDs, not chemical or project data.")
    catalog = describe_catalog(library)
    for key, value in selected.items():
        if value is None:
            continue
        if not isinstance(value, str) or value not in {row["id"] for row in catalog[FIELDS[key]]}:
            raise ValueError("Unknown context catalog ID for " + key + ".")
        result["selection"][key] = value
    if result["selection"].get("source_recipe_id") and any(result["selection"].get(key) for key in CHEMISTRY[1:]):
        raise ValueError("Context cannot combine a source identity with a changed design.")
    pending = context.get("pending") or []
    if not isinstance(pending, list) or len(pending) > 12 or any(not isinstance(row, str) or not 1 <= len(row) <= 700 for row in pending):
        raise ValueError("Unresolved conversation constraints exceed their bound.")
    result["pending"] = list(dict.fromkeys(pending))
    requirements = context.get('requirements') or []
    if not isinstance(requirements, list) or len(requirements) > 8 or any(not isinstance(row, str) or not re.fullmatch(r'(?:(?:amine|aldehyde|hydrazide|hydroxylamine) with [0-9]{1,2} reactive sites|[0-9]{1,3} instances)', row) for row in requirements):
        raise ValueError('Context functional requirements must use bounded recognized count/role phrases.')
    result['requirements'] = list(dict.fromkeys(requirements))
    return result


def model_turn_text(text, context, library):
    """The optional provider sees the new utterance and catalog IDs only."""
    _reject_structure_payload(text)
    context = checked_context(context, library)
    payload = {"current_request": text, "current_catalog_selection": context["selection"],
               "instruction": "Interpret this follow-up against the selected catalog IDs. Return the complete intended selection. Keep unresolved constraints explicit. A current drawing may be custom; appearance-only changes need no replacement monomers."}
    value = json.dumps(payload, ensure_ascii=False)
    if len(value) > 4000:
        raise ValueError("Shorten the message so the catalog context fits the provider request.")
    return value


def _help(text):
    normalized = _normal(text).casefold()
    question = re.search(r"\b(?:help|how|where|what|explain)\b", normalized)
    if not question:
        return None
    topics = [
        (r"\b(?:sa|synthetic accessibility)\b", "SA Score is a finite-molecule heuristic from 1 easier to 10 harder. Open the monomer library to see actual scores; it is not a COF synthesis or yield prediction.", "monomer_library", "Open monomer library"),
        (r"\b(?:pores?|pld|lcd|porosity)\b", "Measure pore estimates the finite 2D drawing. Level 1 geometry instead needs an explicit periodic structure and reports qualified PLD/LCD and accessible geometry. No calculation has run from this message.", "science_porosity", "Open Level 1 geometry"),
        (r"\b(?:3[- ]?d|stack|stacking|cif|crystal)\b", "COF structure can illustrate layers of the current finite drawing or inspect an explicit periodic model/CIF. A layered illustration is not an optimized crystal. No structure has been generated from this message.", "science_model", "Open COF structure"),
        (r"\b(?:adsorption|gcmc|raspa|henry)\b", "Level 2 adsorption runs the bundled classical engine on an explicit periodic structure with reviewed force fields, charges and run settings. A short run is not convergence evidence. No simulation has run from this message.", "science_adsorption", "Open Level 2 adsorption"),
        (r"\bbet\b", "BET fitting needs supplied isotherm data and a chosen fit interval. Geometric area and simulated uptake are separate quantities. No BET value has been calculated from this message.", "science_bet", "Open BET from isotherm"),
    ]
    for pattern, answer, action, label in topics:
        if re.search(pattern, normalized):
            return answer, action, label
    if re.fullmatch(r"(?:help|what can you do|how does this work)[?.! ]*", normalized):
        return ("I can plan catalog drawings, keep a selected pair for follow-ups, change appearance, or guide you to the scientific tools. Try 'Draw TFB and TAPB', then 'Make it pastel'. Review each plan before applying it.", "examples", "Show drawing examples")
    return None


def _response(text, summary, context, *, operation="help", status="needs_input"):
    return {"status": status, "operation": operation, "summary": summary, "selection": deepcopy(context["selection"]),
            "prepare_request": None, "source_recipe_id": None, "style_id": context["selection"].get("style_id", "standard"),
            "style_options": {}, "matches": [], "warnings": [], "unresolved": [], "alternatives": [], "assumptions": [],
            "scope": "Reviewable local drawing assistance; no calculation or synthesis prediction is performed by this chat.",
            "conversation": {"context": deepcopy(context), "suggestions": SUGGESTIONS, "changes": [], "user_text": text}}


def describe_turn(text, library, context=None, selection=None, *, allow_unmatched_words=False):
    """Interpret one turn and return a normal drawing plan or a local intent.

    Replacements use exact catalog matching and existing preparation validation.
    Unknown wording stays unresolved by default, even when context supplies IDs.
    """
    if not isinstance(text, str) or not text.strip() or len(text) > 3000:
        raise ValueError("Write a message of 1–3000 characters.")
    _reject_structure_payload(text)
    text = text.strip(); context = checked_context(context, library)
    catalog = describe_catalog(library)
    if re.fullmatch(r"(?:start over|start again|new conversation|reset conversation)[.! ]*", text, re.I):
        empty = checked_context(None, library)
        return _response(text, "Started a new conversation. Describe two catalog monomers or choose a source drawing.", empty, operation="reset")
    help_result = _help(text)
    if help_result:
        answer, action, label = help_result
        result = _response(text, answer, context)
        result["actions"] = [{"id": action, "label": label}]
        return result
    normalized, labels = expand_material_names(_normal(text))
    if labels:
        result = describe_turn(normalized, library, context=context, selection=selection,
                               allow_unmatched_words=allow_unmatched_words)
        result.setdefault("assumptions", []).extend(labels)
        result["conversation"]["user_text"] = text
        return result
    # Existing matching rules retain alias ambiguity and exact name boundaries.
    initial = describe_cof(normalized, library)
    matches = initial["matches"]
    kinds = {kind: [row for row in matches if row["kind"] == kind] for kind in ("source", "monomer", "linkage", "template", "geometry", "style", "palette")}
    base = deepcopy(context["selection"]); before = deepcopy(base)
    extra = []; assumptions = []; parse_text = normalized
    monomers = _choices(kinds["monomer"]); sources = _choices(kinds["source"])
    changing = bool(monomers or kinds["linkage"] or kinds["template"] or kinds["geometry"] or any(row["kind"] == "requirement" for row in matches))
    # An appearance turn can address any local custom recipe without sending it.
    appearance_only = context["origin"] == "current_drawing" and not changing and not sources
    if appearance_only:
        if not context["token"]:
            raise ValueError("Current-drawing changes need a local revision token.")
        # Validate prose against a real catalog pair, used only to exercise the
        # lexical/constraint validator. It is never returned as the user's graph.
        probe_text = re.sub(FOLLOWUP_WORDS, " ", normalized, flags=re.I)
        appearance = {key: value for key, value in (selection or {}).items() if key in {"style_id", "palette"}}
        if selection and any(value and value != base.get(key) for key, value in selection.items() if key in CHEMISTRY):
            extra.append("A current-drawing appearance request cannot replace its monomers or source identity.")
        # The existing planner validates style and explicit unsupported wording;
        # its missing-pair diagnostic is the only expected omission here.
        checked = describe_cof(probe_text, library, selection=appearance or None, allow_unmatched_words=allow_unmatched_words)
        unresolved = [value for value in checked["unresolved"] if not value.startswith("Choose exactly two unambiguous")]
        unresolved += extra + context["pending"]
        if not kinds["style"] and not kinds["palette"] and not appearance:
            unresolved.append("Choose a drawing style or palette for the current drawing.")
        style = checked["style_id"] if kinds["style"] or kinds["palette"] or appearance else base.get("style_id", "standard")
        result = _response(text, "Apply the reviewed appearance to the current drawing, preserving its exact chemical graph and coordinates.", context, operation="style_only", status="ready" if not unresolved else checked["status"])
        result.update(style_id=style, style_options=checked["style_options"], unresolved=list(dict.fromkeys(unresolved)), warnings=checked["warnings"], assumptions=checked["assumptions"], matches=matches)
        result["selection"] = {"style_id": style, "palette": checked["selection"].get("palette")}
        result["context_token"] = context["token"]
        return _finish(result, context, before)
    if context['origin'] == 'current_drawing' and changing:
        assumptions.append('Changing chemistry or connectivity creates a new design; local coordinate edits and custom graph changes are not carried into that new construction.')
    if sources:
        base = {key: value for key, value in base.items() if key not in CHEMISTRY}
    elif changing and base.get("source_recipe_id"):
        source = next(row for row in catalog["sources"] if row["id"] == base["source_recipe_id"])
        ids = list(dict.fromkeys(source["monomer_ids"]))
        base = {key: value for key, value in base.items() if key not in CHEMISTRY}
        if len(ids) == 2:
            base.update(monomer_a_id=ids[0], monomer_b_id=ids[1])
        assumptions.append("Changing source chemistry starts a new design; the source identity and recorded layout will not be inherited.")
    # Explicit A/B commands are consumed only when they name a recognized ID.
    role = re.match(r"^(?:please\s+)?(?:change|set|switch|replace|use)\s+(?:monomer\s+|building block\s+)?([AB])\s+(?:to|with|as|=)\s+", normalized, re.I)
    if role and len(monomers) == 1:
        base["monomer_" + role[1].lower() + "_id"] = monomers[0]
        parse_text = normalized[role.end():]
        for key in ("linkage", "template_id"):
            base.pop(key, None)
        assumptions.append("Recheck linkage and finite template after changing a building block.")
    elif len(monomers) == 2:
        base = {key: value for key, value in base.items() if key not in CHEMISTRY}
        # A narrowly recognized 'replace old with new' retains the other block.
        replacement = re.match(r"^(?:please\s+)?replace\s+(.+?)\s+with\s+(.+?)[.!]*$", normalized, re.I)
        if replacement:
            left_matches = _find(replacement[1], catalog["monomers"], "monomer")
            left = _choices(left_matches)
            right = _choices(_find(replacement[2], catalog["monomers"], "monomer"))
            roles = [key for key in ("monomer_a_id", "monomer_b_id") if before.get(key) in left]
            old_name_only = not _mask(replacement[1], left_matches).strip(' .,!')
            if len(left) == len(right) == len(roles) == 1 and old_name_only:
                base.update({key: before[key] for key in ("monomer_a_id", "monomer_b_id") if before.get(key)})
                base[roles[0]] = right[0]; parse_text = replacement[2]
                assumptions.append("Replace the named current block and recheck chemistry; keep the other block.")
            else:
                extra.append("The replacement must name exactly one current block and one unambiguous catalog replacement.")
    elif len(monomers) == 1 and not role:
        identifier = monomers[0]
        current = [base.get(key) for key in ("monomer_a_id", "monomer_b_id")]
        if identifier not in current:
            missing = [key for key in ("monomer_a_id", "monomer_b_id") if not base.get(key)]
            if missing:
                base[missing[0]] = identifier
            else:
                candidates = []
                for key in ("monomer_a_id", "monomer_b_id"):
                    candidate = {k: v for k, v in base.items() if k not in {"linkage", "template_id"}}
                    candidate[key] = identifier
                    check = describe_cof("Draw", library, selection=candidate)
                    if check["status"] == "ready":
                        candidates.append(candidate)
                if len(candidates) == 1:
                    base = candidates[0]
                    assumptions.append("Replace the only chemically compatible building-block role and recheck linkage/template.")
                else:
                    extra.append("Specify which block to change: 'Change A to NAME' or 'Change B to NAME'.")
    for kind, key in (("linkage", "linkage"), ("template", "template_id"), ("style", "style_id"), ("palette", "palette")):
        if kinds[kind]:
            base.pop(key, None)
    if kinds["palette"] and not kinds["style"]:
        base.pop("style_id", None)  # Let an explicit color request choose a visible supported style.
    if kinds["geometry"]:
        base.pop("template_id", None)
    if context["selection"]:
        parse_text = re.sub(FOLLOWUP_WORDS, " ", parse_text, flags=re.I).strip() or "Draw"
    if selection:
        if set(selection) - SELECTION_KEYS:
            raise ValueError("Manual choices must contain catalog IDs only.")
        if selection.get("source_recipe_id"):
            base = {key: value for key, value in base.items() if key not in CHEMISTRY}
        base.update({key: value for key, value in selection.items() if value is not None})
    parse_text = parse_text + ' ' + ' '.join(context['requirements'])
    result = describe_cof(parse_text.strip(), library, selection=base or None, allow_unmatched_words=allow_unmatched_words)
    result["operation"] = "drawing"
    result["assumptions"] += assumptions
    result["unresolved"] = list(dict.fromkeys(result["unresolved"] + extra + context["pending"]))
    if result["unresolved"]:
        if result["status"] == "ready":
            result["status"] = "needs_input"
        result["source_recipe_id"] = None; result["prepare_request"] = None
    result["conversation"] = {"user_text": text}
    return _finish(result, context, before)


def _finish(result, context, before):
    selected = {key: value for key, value in result["selection"].items() if value is not None}
    # Missing slots can be completed in later turns. Unsupported/unknown claims
    # stay visible until the user explicitly starts a fresh conversation.
    pending = list(context["pending"])
    for message in result["unresolved"]:
        if result["status"] == "unsupported" or message.startswith(("Unresolved words or names:", "Negation or replacement", "A requested physical", "Prose drawing cannot")):
            if message not in pending:
                pending.append(message)
    requirements = list(context['requirements'])
    for match in result.get('matches', []):
        role = re.fullmatch(r'(amine|aldehyde|hydrazide|hydroxylamine)_([0-9]{1,2})_sites', match.get('id', '')) if match.get('kind') == 'requirement' else None
        instances = re.fullmatch(r'([0-9]{1,3})_instances', match.get('id', '')) if match.get('kind') == 'requirement' else None
        phrase = f'{role[1]} with {role[2]} reactive sites' if role else f'{instances[1]} instances' if instances else None
        if phrase and phrase not in requirements:
            requirements.append(phrase)
    if len(requirements) > 8:
        raise ValueError('Too many separate functional requirements. Start a focused new conversation.')
    context_selection = {**context["selection"], **selected} if result["operation"] == "style_only" else selected
    next_context = {"selection": context_selection, "pending": pending[:12], "requirements": requirements[:8], "origin": "current_drawing" if result["operation"] == "style_only" else "plan", "token": context["token"] if result["operation"] == "style_only" else None}
    result.setdefault("conversation", {}).update(context=next_context, suggestions=SUGGESTIONS,
        changes=[{"field": key, "before": before.get(key), "after": selected.get(key)} for key in FIELDS if before.get(key) != selected.get(key)])
    if context["pending"]:
        result["warnings"].append("An earlier unresolved constraint is still retained. Start a new conversation to explicitly discard it.")
    return result
