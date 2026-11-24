"""Catalog identity, explicit-constraint retention, and real preparation checks."""
from copy import deepcopy
import json
from pathlib import Path

from fastapi.testclient import TestClient
import pytest

from cof_vision.cofdraw_describe import describe_catalog, describe_cof, validate_selection
from cof_vision.cofdraw_desktop import create_desktop_app
from cof_vision.cofdraw_chemistry import prepare_template
from cof_vision.drawing_api import PrepareDrawing, validate_recipe
from cof_vision.framework_drawing import VERSION, build_patch


@pytest.fixture(scope="module")
def library():
    with TestClient(create_desktop_app(Path(__file__).resolve().parents[1] / "data")) as client:
        result = client.get("/api/drawing/library")
        assert result.status_code == 200
        return result.json()


def test_catalog_is_an_id_allowlist_without_smiles_and_keeps_qualifications(library):
    result = describe_catalog(library)
    for kind in ["sources", "monomers", "linkages", "templates", "styles", "palettes"]:
        assert result[kind]
        assert len({row["id"] for row in result[kind]}) == len(result[kind])
        assert all(row["label"] and row["aliases"] for row in result[kind])
    assert '"smiles"' not in json.dumps(result)
    source = next(row for row in result["sources"] if row["id"] == "pyurea-source-layout-v1")
    assert source["sha256"] and source["qualifications"]
    assert not source["source_registered"]
    assert len(result["demos"]) >= 6


def test_every_demo_has_honest_status_and_ready_designs_really_build(library):
    for demo in describe_catalog(library)["demos"]:
        result = describe_cof(demo["text"], library)
        assert result["status"] == demo["expected_status"], (demo, result["unresolved"])
        if result["prepare_request"]:
            request = PrepareDrawing(**result["prepare_request"])
            prepared = prepare_template(**request.model_dump())
            recipe = {"schema_version": VERSION, "recipe_id": "describe-test-design", **prepared}
            validate_recipe(recipe)
            graph, _ = build_patch(recipe)
            assert graph.GetNumAtoms() > 20
            assert result["source_recipe_id"] is None
            assert result["compatibility"]["validated_by"] == "prepare_template"


@pytest.mark.parametrize("text,identifier", [
    ("Draw TpPa", "tppa-pore-v1"),
    ("Draw TpPa-SO3H", "wiley-r198-tppa-so3h-pore-0002"),
    ("Draw SCOF", "wiley-r198-tppa-so3h-pore-0002"),
    ("Draw FCOF1 with highlights", "f-cof-1-source-v1"),
    ("Draw F-COF-1 in source style", "f-cof-1-source-v1"),
    ("Draw Pyurea", "pyurea-source-layout-v1"),
    ("Open the Pyurea source drawing", "pyurea-source-layout-v1"),
    ("Load the F-COF-1 drawing", "f-cof-1-source-v1"),
    ("Draw COF-BP-Ni2", "cof-bp-ni2-network-v1"),
    ("Draw BTA-COF4", "bta-cof-4-pore-v1"),
    ("Draw COF-300", "cof-300-projection-v1"),
])
def test_source_names_longest_match_and_aliases_preserve_identity(text, identifier, library):
    result = describe_cof(text, library)
    assert result["status"] == "ready", result["unresolved"]
    assert result["source_recipe_id"] == identifier
    assert result["prepare_request"] is None
    assert result["source_metadata"]["sha256"]
    assert result["warnings"]


def test_rejected_source_hypothesis_qualification_survives(library):
    result = describe_cof("Draw Herringbone", library)
    assert result["status"] == "ready"
    assert "rejected" in " ".join(result["warnings"]).lower()


@pytest.mark.parametrize("text", ["Draw a projected network", "Draw a local boronate fragment",
                                  "Draw the parent network", "Draw coordination sites"])
def test_generic_source_qualifiers_are_not_material_aliases(text, library):
    result = describe_cof(text, library)
    assert result["status"] != "ready"
    assert not any(row["kind"] == "source" for row in result["matches"])


def test_monomer_order_reversed_safely_and_names_are_exact(library):
    result = describe_cof("Draw TAPB and TFB with imine linkage in a hexagonal template", library)
    assert result["status"] == "ready", result["unresolved"]
    assert result["selection"]["monomer_a_id"] == "tfb"
    assert result["selection"]["monomer_b_id"] == "tapb"
    assert result["selection"]["template_id"] == "hcb_3_3"
    assert any("Reorder" in line for line in result["assumptions"])
    alias = describe_cof("Draw TFB and TAB", library)
    assert alias["status"] == "ready"
    assert alias["selection"]["monomer_b_id"] == "tapb"
    assert alias["source_recipe_id"] is None


def test_inferred_choices_are_disclosed_without_inheriting_source_identity(library):
    result = describe_cof("Draw TFB and TAPB in pastel colors", library)
    assert result["status"] == "ready"
    assert result["selection"]["linkage"] == "imine"
    assert result["selection"]["style_id"] == "pastel_blocks"
    assert result["style_options"] == {"palette": "pastel"}
    assert len(result["assumptions"]) >= 3
    assert result["source_recipe_id"] is None
    assert "source_metadata" not in result
    tp = describe_cof("Draw Tp + Pa", library)
    assert tp["status"] == "ready"
    assert tp["selection"]["linkage"] == "beta_ketoenamine"
    assert tp["selection"]["template_id"] == "hcb_3_2"
    assert any("catalog starting-pair preference" in line for line in tp["assumptions"])


def test_short_monomer_aliases_need_deliberate_case_or_context(library):
    assert describe_cof("Draw tp + pa", library)["status"] == "ready"
    assert describe_cof("Draw monomers ta and TAPB", library)["status"] == "ready"
    for text in ["Draw a COF at 10 Pa", "Draw ha and TFB", "Draw ta and TAPB", "Draw TFB and TAPBish"]:
        assert describe_cof(text, library)["status"] != "ready", text


def test_duplicate_bd_abbreviation_stays_ambiguous_until_selected(library):
    text = "Draw TFB and BD with imine linkage"
    result = describe_cof(text, library)
    assert result["status"] == "needs_input"
    assert any("Ambiguous monomer alias" in line for line in result["unresolved"])
    selected = describe_cof(text, library, {"monomer_a_id": "tfb", "monomer_b_id": "benzidine"})
    assert selected["status"] == "ready", selected["unresolved"]
    assert selected["selection"]["monomer_b_id"] == "benzidine"
    labelled = describe_cof("BD + 13C-Tfppy", library)
    assert labelled["selection"]["monomer_b_id"] == "tfppy_13c" and labelled["selection"]["monomer_a_id"] is None
    assert {row["id"] for row in labelled["alternatives"]} == {"benzidine", "recovered_bd"}
    for choice in labelled["alternatives"]:
        selected_bd = {choice["field"]: choice["id"]}
        confirmed = describe_cof("BD + 13C-Tfppy", library, selected_bd)
        assert confirmed["status"] == "ready"
        assert {confirmed["selection"][key] for key in ("monomer_a_id", "monomer_b_id")} == {choice["id"], "tfppy_13c"}
        assert confirmed["prepare_request"]["template_id"] == "sql_4_2"
        assert {row["id"] for row in confirmed["alternatives"]} == {"kgm_4_2", "nbo_4_2"}
        spatial = describe_cof("BD + 13C-Tfppy", library, {**selected_bd, "template_id": "nbo_4_2"})
        assert spatial["status"] == "ready"
        assert spatial["prepare_request"]["template_id"] == "nbo_4_2"


@pytest.mark.parametrize("text", [
    "Draw a COF", "Draw an imine COF", "Draw TFB and mysterydiamine", "Draw TFB and TAPB and XYZ",
    "Draw TFB and TAPB and DETH", "Draw TpPa with square geometry", "Draw TpPa using TFB",
    "Draw TpPa and FCOF-1", "Draw TFB and TAPB with imine and acylhydrazone linkages",
    "Draw TFB and TAPB without imine", "Do not draw TpPa", "Draw Pyurea but replace urea",
    "Draw TFB and TAPB with neon magenta colors", "Draw TFB and TAPB with square geometry",
    "Draw TAPM and TFB with imine linkage and a square template",
])
def test_ambiguous_unknown_conflicting_or_negative_prose_never_silently_draws(text, library):
    result = describe_cof(text, library)
    assert result["status"] != "ready", (text, result)
    assert result["unresolved"] and result["alternatives"]
    assert result["prepare_request"] is None and result["source_recipe_id"] is None


def test_boronate_request_without_precursors_needs_input(library):
    assert describe_cof("Draw a boronate COF", library)["status"] == "needs_input"


@pytest.mark.parametrize("text", [
    "Draw a dioxin-linked COF",
    "Draw azo-linked TFB and TAPB", "Generate a 3D crystal of TpPa", "Draw TpPa and export CIF",
    "Draw TFB and TAPB with DFT optimization", "Draw a COF with a BET area of 2000",
])
def test_unsupported_requests_are_explicit(text, library):
    result = describe_cof(text, library)
    assert result["status"] == "unsupported", result
    assert result["prepare_request"] is None and result["source_recipe_id"] is None


def test_id_selection_fills_missing_choices_and_real_validation_is_required(library):
    selection = {"monomer_a_id": "tfb", "monomer_b_id": "tapb", "linkage": "imine",
                 "template_id": "hcb_3_3", "style_id": "standard", "palette": None}
    result = validate_selection(selection, library)
    assert result["status"] == "ready"
    assert result["prepare_request"] and result["compatibility"]
    bad = {**selection, "template_id": "sql_4_2"}
    result = validate_selection(bad, library)
    assert result["status"] == "needs_input"
    assert any("validation failed" in text for text in result["unresolved"])


@pytest.mark.parametrize("text,selection", [
    ("Draw TpPa", {"source_recipe_id": "f-cof-1-source-v1"}),
    ("Draw TpPa and FCOF1", {"source_recipe_id": "tppa-pore-v1"}),
    ("Draw TFB and TAPB", {"monomer_a_id": "tp", "monomer_b_id": "pa"}),
    ("Draw TFB and TAPB with imine linkage", {"linkage": "secondary_amine"}),
    ("Draw TFB and TAPB with square geometry", {"template_id": "hcb_3_3"}),
    ("Draw TFB and UnknownAmine", {"monomer_a_id": "tfb", "monomer_b_id": "tapb"}),
    ("Do not draw TpPa", {"source_recipe_id": "tppa-pore-v1"}),
    ("Create a 3D crystal", {"source_recipe_id": "tppa-pore-v1"}),
    ("Draw a COF", {"source_recipe_id": "tppa-pore-v1", "monomer_a_id": "tfb"}),
    ("Draw a COF", {"source_recipe_id": "invented"}),
    ("Draw a COF", {"monomer_a_id": "tfb", "monomer_b_id": "tapb", "smiles": "NC=O"}),
])
def test_external_id_proposals_cannot_bypass_prose_or_allowlists(text, selection, library):
    result = describe_cof(text, library, selection)
    assert result["status"] != "ready", result
    assert result["unresolved"]
    assert result["prepare_request"] is None and result["source_recipe_id"] is None


def test_input_bounds_and_reproducibility_without_mutation(library):
    original = deepcopy(library)
    text = "Draw Pyurea with building-block highlights and a soft palette"
    assert describe_cof(text, library) == describe_cof(text, library)
    assert library == original
    for value in [None, "", " ", "x"*4001, 123]:
        with pytest.raises(ValueError):
            describe_cof(value, library)
    with pytest.raises(ValueError):
        describe_cof(text, library, selection=[])


def test_reviewed_appearance_changes_are_disclosed_and_preserve_source(library):
    result = describe_cof("Draw TpPa in standard style with pastel colors", library,
                          {"style_id": "node_highlight", "palette": "teal_orange"})
    assert result["status"] == "ready"
    assert result["source_recipe_id"] == "tppa-pore-v1"
    assert result["style_id"] == "node_highlight"
    assert result["style_options"]["palette"] == "teal_orange"
    assert len([line for line in result["assumptions"] if "overrides the appearance" in line]) == 2


def test_reviewed_selection_can_interpret_unmatched_language_with_visible_warning(library):
    text = "Draw TFB using the tritopic aromatic amine in an elegant scheme"
    selection = {"monomer_a_id": "tfb", "monomer_b_id": "tapb", "linkage": "imine", "template_id": "hcb_3_3"}
    assert describe_cof(text, library, selection)["status"] == "needs_input"
    result = describe_cof(text, library, selection, allow_unmatched_words=True)
    assert result["status"] == "ready"
    assert any("local rules could not match" in line for line in result["warnings"])
    assert any("elegant" in line for line in result["warnings"])
    assert any(row["kind"] == "requirement" and row["id"] == "amine_3_sites" for row in result["matches"])
    # A bare flag without an actual allowlisted proposal cannot make a guess.
    assert describe_cof(text, library, allow_unmatched_words=True)["status"] == "needs_input"


@pytest.mark.parametrize("text", [
    "Draw TFB and Pa", "Draw TpPa", "Draw TFB and TAPB without imine",
    "Draw TFB and TAPB using a boronate ester linkage", "Create a 3D crystal",
    "Calculate the pore size of TFB and TAPB", "Draw TFB and TAPB in a square template",
])
def test_relaxed_wording_never_relaxes_recognized_chemistry_or_operation_constraints(text, library):
    result = describe_cof(text, library, {"monomer_a_id": "tfb", "monomer_b_id": "tapb",
                          "linkage": "imine", "template_id": "hcb_3_3"}, allow_unmatched_words=True)
    assert result["status"] != "ready"
    assert result["prepare_request"] is None


@pytest.mark.parametrize("text", [
    "Create a 3-D structure of TpPa", "Generate a 3 dimensional model of TpPa",
    "Generate XYZ coordinates of TpPa", "Export TpPa as a PDB file",
    "Draw TpPa with a 2 nm pore", "Draw TpPa with a pore diameter of 20 angstrom",
    "Draw TpPa with a 2-nm-wide pore", "Draw TpPa with a 20 A pore",
])
def test_coordinate_and_quantitative_pore_requests_cannot_become_drawings(text, library):
    result = describe_cof(text, library, {"source_recipe_id": "tppa-pore-v1"}, allow_unmatched_words=True)
    assert result["status"] == "unsupported"
    assert result["source_recipe_id"] is None


@pytest.mark.parametrize("text,count", [
    ("Use a tetratopic amine with TFB", 4), ("Use a 4-connected amine with TFB", 4),
    ("Use a hexatopic aldehyde and TAPB", 6), ("Draw TFB and a diamine", 2),
    ("Draw TFB with an amine with four reactive sites", 4),
])
def test_explicit_functional_counts_are_hard_constraints(text, count, library):
    selection = {"monomer_a_id": "tfb", "monomer_b_id": "tapb", "linkage": "imine", "template_id": "hcb_3_3"}
    result = describe_cof(text, library, selection, allow_unmatched_words=True)
    assert result["status"] == "needs_input"
    assert any(f"{count}-site" in line for line in result["unresolved"])
    assert result["prepare_request"] is None


def test_matching_functional_count_and_instance_requirements_are_accepted(library):
    selection = {"monomer_a_id": "tfb", "monomer_b_id": "tapb", "linkage": "imine", "template_id": "hcb_3_3"}
    result = describe_cof("Use a tritopic amine with TFB in a six-unit pore", library, selection)
    assert result["status"] == "ready", result["unresolved"]
    assert result["prepare_request"]["instances"] == 6
    result = describe_cof("Draw TFB and TAPB with 12 instances", library, selection, allow_unmatched_words=True)
    assert result["status"] == "needs_input"
    assert any("instance count" in line for line in result["unresolved"])


def test_c3_is_not_interpreted_as_connectivity_and_illustrative_styles_stay_drawings(library):
    result = describe_cof("Draw TFB and TAPB with C3 symmetry", library)
    assert result["status"] == "needs_input"
    assert not any(row["kind"] == "requirement" for row in result["matches"])
    result = describe_cof("Draw TpPa with molecular spheres illustration", library)
    assert result["status"] == "ready"
    assert result["style_id"] == "molecular_spheres"
    assert any("not a periodic atomic structure" in line for line in result["warnings"])
    assert describe_cof("Export TpPa as SVG and PNG files", library)["status"] == "ready"


def test_explicit_colored_bonds_are_kept_with_pastel_palette(library):
    result = describe_cof("Draw TFB and TAPB with colored bonds and a pastel palette", library)
    assert result["status"] == "ready"
    assert result["style_id"] == "red_blue"
    assert result["style_options"]["palette"] == "pastel"


@pytest.mark.parametrize("pair", [("tapm", "ta"), ("ta", "tapm")])
def test_diamond_pair_selects_spatial_cage_and_keeps_canvas_projection(library, pair):
    import numpy as np
    from cof_vision.cofdraw_science import patch_3d

    plan = describe_cof("Draw an imine COF", library, {"monomer_a_id": pair[0], "monomer_b_id": pair[1]})
    assert plan["status"] == "ready", plan["unresolved"]
    assert plan["selection"]["template_id"] == "dia_4_2"
    recipe = dict(schema_version=VERSION, recipe_id="dia-test", **prepare_template(**plan["prepare_request"]))
    validate_recipe(recipe)
    mol, audit = build_patch(recipe)
    canvas = mol.GetConformer().GetPositions()
    assert np.all(canvas[:, 2] == 0)
    centers = {item["id"]: np.array(item["layout_position"]) for item in recipe["instances"]}
    link_ends = {}
    for edge in recipe["connections"]:
        link_ends.setdefault(edge["b"]["instance"], []).append(edge["a"]["instance"])
    sides = [centers[a] - centers[b] for a, b in link_ends.values()]
    assert sum(abs(dx) < .01 * abs(dy) for dx, dy in sides) == 3
    model = patch_3d(recipe)
    assert model["kind"] == "finite_dia_cage"
    graph = patch_3d(recipe, topology_only=True)
    assert graph["kind"] == "spatial_building_block_net" and graph["net"] == "dia"
    assert graph["viewer"]["units"] == "graph units"
    atoms = model["viewer"]["atoms"]
    assert len(atoms) == mol.GetNumAtoms() == 402
    points = np.array([[atom[k] for k in ("x", "y", "z")] for atom in atoms])
    assert np.linalg.matrix_rank(points - points.mean(axis=0)) == 3
    saved = {(r["instance"], r["fragment_atom_index"]): r["xyz"] for r in recipe["atom_coordinates"]}
    assert np.allclose(points, [saved[o["instance"], o["fragment_atom_index"]] for o in audit["atom_origin"]])
    assert model["viewer"]["bonds"] == [
        dict(a=b.GetBeginAtomIdx(), b=b.GetEndAtomIdx(), order=b.GetBondTypeAsDouble()) for b in mol.GetBonds()
    ]


@pytest.mark.parametrize("pair", [("tapm", "tfppy"), ("tfppy", "tapm"), ("por_cho", "tapm"), ("tfs_cho4", "tapp")])
def test_pts_from_selected_monomers(library, pair):
    from collections import Counter
    import numpy as np
    from cof_vision.cofdraw_science import patch_3d

    plan = describe_cof("Draw an imine COF", library, {"monomer_a_id": pair[0], "monomer_b_id": pair[1]})
    assert plan["status"] == "ready", plan["unresolved"]
    assert plan["selection"]["template_id"] == "pts_4_4"
    recipe = dict(schema_version=VERSION, recipe_id="pts-test", **prepare_template(**plan["prepare_request"]))
    validate_recipe(recipe)
    assert recipe["design_topology"] == "pts"
    used = [(edge[side]["instance"], edge[side]["port"]) for edge in recipe["connections"] for side in ("a", "b")]
    assert len(used) == len(set(used)) == 48
    degrees = Counter(instance for instance, port in used)
    assert len(degrees) == 22 and max(degrees.values()) == 4
    model = patch_3d(recipe)
    assert model["kind"] == "finite_pts_patch" and not model["periodic"]
    graph = patch_3d(recipe, topology_only=True)
    assert graph["kind"] == "spatial_building_block_net" and graph["net"] == "pts"
    assert graph["viewer"]["units"] == "graph units"
    atoms, bonds = model["viewer"]["atoms"], model["viewer"]["bonds"]
    mol, audit = build_patch(recipe)
    assert len(atoms) == mol.GetNumAtoms()
    assert sum(a["continuation_site"] for a in atoms) == 40
    assert np.all(mol.GetConformer().GetPositions()[:, 2] == 0)
    assert bonds == [dict(a=b.GetBeginAtomIdx(), b=b.GetEndAtomIdx(), order=b.GetBondTypeAsDouble())
                     for b in mol.GetBonds()]
    coords = np.array([[a[k] for k in ("x", "y", "z")] for a in atoms])
    assert np.linalg.matrix_rank(coords - coords.mean(axis=0)) == 3
    for b in bonds:
        a, c = b["a"], b["b"]
        if atoms[a]["instance"] != atoms[c]["instance"]:
            assert atoms[a]["component"] != atoms[c]["component"]
            assert b["order"] == 2
            assert 1.1 < np.linalg.norm(coords[a] - coords[c]) < 1.8


@pytest.mark.parametrize("a,b", [("tfppy", "pytta"), ("tfs_cho4", "tapm"), ("ta", "tapm")])
def test_pts_requires_both_geometries_and_four_sites(library, a, b):
    monomers = {r["id"]: r["smiles"] for r in library["monomer_presets"]}
    with pytest.raises(ValueError, match="requires|needs"):
        prepare_template(monomers[a], monomers[b], "imine", a, b, 22, "pts_4_4")
    if a == "tfppy":
        assert prepare_template(monomers[a], monomers[b], "imine", a, b, 4, "sql_4_4")["design_topology"] == "sql"


def test_pts_periodic_neighbor_geometry():
    import numpy as np
    from cof_vision.cofdraw_net_patch import _pts_neighbors

    metric = np.array([1., 1., 1 / np.sqrt(2)])
    for site in [(0, 0, 1), (0, 1, 0), (0, 0, 3), (1, 0, 2)]:
        neighbors = _pts_neighbors(site)
        assert len(set(neighbors)) == 4
        assert all(site in _pts_neighbors(n) for n in neighbors)
        arms = (np.array(neighbors) - site) * metric
        if site[2] % 2:
            arms /= np.linalg.norm(arms, axis=1)[:, None]
            assert np.allclose((arms @ arms.T)[np.triu_indices(4, 1)], -1 / 3)
        else:
            assert np.linalg.matrix_rank(arms) == 2
            assert np.allclose(arms[0], -arms[3])

@pytest.mark.parametrize("a,b", [("ta", "pytta"), ("tfppy", "pa"), ("bfbapy", "pytta")])
def test_kagome_has_six_triangles_around_a_hexagon(library, a, b):
    from collections import defaultdict

    monomers = {r["id"]: r["smiles"] for r in library["monomer_presets"]}
    template = next(t for t in library["templates"] if t["id"] == "kgm_4_2")
    request = PrepareDrawing(monomer_a=monomers[a], monomer_b=monomers[b], template_id=template["id"],
                             instances=template["instances"])
    recipe = dict(schema_version=VERSION, recipe_id="kgm-star", **prepare_template(**request.model_dump()))
    validate_recipe(recipe)
    assert len(recipe["instances"]) == template["instances"] == 30
    mol, _ = build_patch(recipe)
    assert sum(atom.GetAtomicNum() == 0 for atom in mol.GetAtoms()) == 12
    linker_ends = defaultdict(list)
    for edge in recipe["connections"]:
        linker_ends[edge["b"]["instance"]].append(edge["a"]["instance"])
    assert len(linker_ends) == 18
    neighbors = defaultdict(set)
    for left, right in linker_ends.values():
        neighbors[left].add(right)
        neighbors[right].add(left)
    ring = {n for n, adjacent in neighbors.items() if len(adjacent) == 4}
    tips = set(neighbors) - ring
    assert len(tips) == len(ring) == 6
    assert all(len(neighbors[n] & ring) == 2 for n in ring)
    ring_edges = {frozenset((n, m)) for n in ring for m in neighbors[n] & ring}
    assert {frozenset(neighbors[tip]) for tip in tips} == ring_edges


def test_dia_view_matches_cof300_cage_projection(library):
    import numpy as np

    monomers = {r["id"]: r["smiles"] for r in library["monomer_presets"]}
    recipe = prepare_template(monomers["ta"], monomers["tapm"], "imine", "TA", "TAPM", 22, "dia_4_2")
    centers = np.array([i["layout_position"] for i in recipe["instances"] if i["component"] == "B"])
    # Approximate junction positions read from Fig. 2A. Compare shape in the
    # displayed orientation, allowing translation and one overall scale factor.
    reference = np.array([(432,87), (166,425), (285,488), (613,394), (225,159),
                          (345,226), (675,126), (80,560), (406,466), (524,536)], dtype=float)
    reference[:, 1] *= -1
    for points in (centers, reference):
        points -= points.mean(axis=0)
        points /= np.linalg.norm(points)
    assert np.linalg.norm(centers - reference) < .24


@pytest.mark.parametrize("alias,node,spacer,extent", [
    ("dia", "tapm", "ta", "circuit"), ("pts", "tapm", "tfppy", "circuit"),
    ("kagome", "pytta", "ta", "hexagon"),
])
def test_request_one_face_of_full_molecular_drawing(library, alias, node, spacer, extent):
    from cof_vision.cofdraw_chemistry import design_signature
    from cof_vision.cofdraw_repeat import repeat_patch

    selections = {"monomer_a_id": node, "monomer_b_id": spacer}
    plan = describe_cof(f"Draw one pore with {alias} topology and imine linkage", library, selections)
    assert plan['status'] == 'ready', plan['unresolved']
    choice = plan['drawing_extent']
    assert choice['extent_mode'] == extent and choice['pores'] == 1
    recipe = prepare_template(**plan['prepare_request'])
    recipe['design_template_signature'] = design_signature(recipe)
    _, _, _, summary, _ = repeat_patch(recipe, choice['repeats'], extent_mode=extent)
    assert summary['pore_count'] == 1


@pytest.mark.parametrize("a,b,key,n,joins", [
    ("ta", "tapm", "dia_4_2", 22, 24),
    ("tfppy", "tapm", "pts_4_4", 22, 24),
    ("ta", "pytta", "kgm_4_2", 30, 36),
])
def test_molecular_net_imines_have_two_trigonal_centers(library, a, b, key, n, joins):
    import numpy as np
    from rdkit import Chem

    monomers = {r["id"]: r["smiles"] for r in library["monomer_presets"]}
    recipe = prepare_template(monomers[a], monomers[b], "imine", a, b, n, key)
    mol, audit = build_patch(dict(schema_version=VERSION, recipe_id=key, **recipe))
    positions = {(row["instance"], row["fragment_atom_index"]): row.get("xyz", row["xy"])
                 for row in recipe["atom_coordinates"]}
    xyz = np.array([positions[o["instance"], o["fragment_atom_index"]] for o in audit["atom_origin"]])
    imines = mol.GetSubstructMatches(Chem.MolFromSmarts("c-[CH]=N-c"))
    assert len(imines) == joins
    for match in imines:
        vectors = np.diff(xyz[list(match)], axis=0)
        lengths = np.linalg.norm(vectors, axis=1)
        assert lengths[1] < min(lengths[0], lengths[2])
        directions = vectors / lengths[:, None]
        cosines = -np.sum(directions[:-1] * directions[1:], axis=1)
        assert cosines == pytest.approx([-.5, -.5], abs=1e-4)


@pytest.mark.parametrize("node", ["tape", "pytta"])
def test_kagome_boundary_arms_continue_opposite_imines(library, node):
    import numpy as np
    from rdkit import Chem
    from cof_vision.framework_drawing import _ports

    smiles = {r["id"]: r["smiles"] for r in library["monomer_presets"]}
    recipe = prepare_template(smiles["ta"], smiles[node], "imine", "TA", node, 30, "kgm_4_2")
    mol, audit = build_patch(dict(schema_version=VERSION, recipe_id="kgm-boundary", **recipe))
    xy = mol.GetConformer().GetPositions()[:, :2]
    fragment = Chem.MolFromSmiles(recipe["components"]["B"]["fragment_smiles"])
    ports = _ports(fragment)
    checked = 0
    for block in recipe["instances"]:
        if block["component"] != "B":
            continue
        atoms = {o["fragment_atom_index"]: i for i, o in enumerate(audit["atom_origin"])
                 if o["instance"] == block["id"]}
        arms = []
        for dummy, terminal in ports.values():
            anchor = next(a.GetIdx() for a in fragment.GetAtomWithIdx(terminal).GetNeighbors() if a.GetIdx() != dummy)
            tip = next(a.GetIdx() for a in mol.GetAtomWithIdx(atoms[terminal]).GetNeighbors() if a.GetIdx() != atoms[anchor])
            arms.append(xy[[atoms[anchor], atoms[terminal], tip]])
            checked += dummy in atoms
        center = np.mean([arm[0] for arm in arms], axis=0)
        for arm in arms:
            assert any(np.allclose(2 * center - arm, opposite, atol=1e-4) for opposite in arms)
    assert checked == 12
