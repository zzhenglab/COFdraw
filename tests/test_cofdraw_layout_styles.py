import copy
import json
import math
from pathlib import Path
import xml.etree.ElementTree as ET

import numpy as np
import pytest
from rdkit import Chem

from cof_vision.cofdraw_core import render_styled
from cof_vision.cofdraw_layout import clean_layout
from cof_vision.cofdraw_styles import STYLE_FAMILIES, style_configuration
from cof_vision.framework_drawing import build_patch, linkage_fragment, VERSION


def test_random_style_changes_depiction_even_with_repeated_seed():
    history = ['stacked_slipped']
    for seed in [42] * 10 + list(range(30)):
        chosen, _, _ = style_configuration('random', seed, recent_styles=history)
        assert STYLE_FAMILIES.get(chosen, chosen) not in {
            STYLE_FAMILIES.get(s, s) for s in history[-4:]
        }
        history.append(chosen)


def pore():
    a = linkage_fragment('O=Cc1cc(C=O)cc(C=O)c1', 'imine', 'aryl_aldehyde')
    b = linkage_fragment('Nc1cc(N)cc(N)c1', 'imine', 'primary_arylamine')
    return {'schema_version': VERSION, 'recipe_id': 'layout-fixture', 'designed_variant': True,
            'components': {'a': a, 'b': b},
            'instances': [{'id': str(i), 'component': 'a' if i % 2 else 'b', 'incoming_port': 1, 'outgoing_port': 2,
                           'layout_position': [math.cos(i*math.pi/3)*15, math.sin(i*math.pi/3)*15]} for i in range(6)],
            'connections': [{'a': {'instance': str(i), 'port': 2}, 'b': {'instance': str((i+1)%6), 'port': 1}, 'bond_order': 2} for i in range(6)]}


def test_explicit_cleanup_preserves_rigid_blocks_and_graph_and_remains_manual():
    recipe = pore()
    old, before = build_patch(recipe)
    result = clean_layout(recipe)
    fitted, after = build_patch(result['recipe'])
    assert recipe.get('layout_refinement') is None  # Caller-owned recipe is intact.
    assert result['recipe']['layout_refinement'] == 'manual'
    assert before['fragment_serialization_sha256'] == after['fragment_serialization_sha256']
    a = old.GetConformer().GetPositions()[:, :2]
    b = fitted.GetConformer().GetPositions()[:, :2]
    for item in recipe['instances']:
        indices = [i for i, origin in enumerate(before['atom_origin']) if origin['instance'] == item['id']]
        assert np.allclose(np.linalg.norm(a[indices, None]-a[np.array(indices)[None]], axis=2),
                           np.linalg.norm(b[indices, None]-b[np.array(indices)[None]], axis=2), atol=1e-9)
    assert after['layout_diagnostics']['nonbonded_overlap_pairs'] == []
    assert result['audit']['layout_refinement']['quality'] == 'regular'
    assert result['audit']['layout_refinement']['final_max_junction_angle_deviation_degrees'] < 1
    assert after['depiction_geometry']['max_to_median_ratio'] < before['depiction_geometry']['max_to_median_ratio']


def test_manual_source_coordinates_survive_render_and_block_drag_is_explicit(tmp_path):
    recipe = pore()
    mol, audit = build_patch(recipe)
    xy = mol.GetConformer().GetPositions()[:, :2]
    # A source-aligned skew is supplied directly, not inferred from a template.
    xy[:, 0] += xy[:, 1]*.12
    recipe['atom_coordinates'] = [{**origin, 'xy': position.tolist()} for origin, position in zip(audit['atom_origin'], xy)]
    recipe['layout_refinement'] = 'rigid_block_pose_v1'
    fixed, fixed_audit = build_patch(recipe)
    assert np.allclose(fixed.GetConformer().GetPositions()[:, :2], xy)
    assert fixed_audit['layout_refinement']['method'] == 'explicit_atom_coordinates'
    rendered = render_styled(recipe, tmp_path, 700, 'red_blue')
    layout = rendered['drawing_style']['layout']
    assert np.allclose(layout['atom_coordinates'], xy)
    pixels = np.array(layout['atom_pixels'])
    assert np.allclose(pixels[:, 0]-pixels[0, 0], (xy[:, 0]-xy[0, 0])*layout['pixels_per_unit'])
    assert np.allclose(pixels[:, 1]-pixels[0, 1], -(xy[:, 1]-xy[0, 1])*layout['pixels_per_unit'])
    manual = copy.deepcopy(recipe)
    manual.pop('atom_coordinates')
    manual['layout_refinement'] = 'manual'
    manual['instances'][0]['layout_position'][0] += 2
    moved, moved_audit = build_patch(manual)
    assert not moved_audit['layout_refinement']['applied']
    assert not np.allclose(moved.GetConformer().GetPositions(), fixed.GetConformer().GetPositions())


def test_chemical_styles_share_geometry_and_color_native_carbon_joins(tmp_path):
    recipe = pore()
    recipe['layout_refinement'] = 'rigid_block_translation_v1'
    original, _ = build_patch(recipe)
    original_xyz = original.GetConformer().GetPositions()
    graph_exports = []
    for style in ('standard', 'red_blue', 'glow_outline', 'pastel_blocks', 'stacked_slipped', 'molecular_spheres','half_space_stack','pore_perspective','element_compact'):
        out = tmp_path/style
        audit = render_styled(recipe, out, 650, style)
        assert np.allclose(audit['drawing_style']['layout']['atom_coordinates'], original_xyz[:, :2])
        graph_exports.append((out/'layout-fixture.mol').read_bytes())
        assert (out/'layout-fixture.png').read_bytes().startswith(b'\x89PNG')
        if style == 'half_space_stack':
            svg = ET.parse(out/'layout-fixture.svg').getroot()
            assert len(svg.findall('.//{*}use')) == audit['drawing_style']['options']['layer_count']-1
            assert svg.find(".//*[@id='half-space-sheet']") is not None
            assert audit['drawing_style']['layout']['display_projection'] is False
        if style == 'pore_perspective':
            import re
            from PIL import Image
            svg = ET.parse(out/'layout-fixture.svg').getroot()
            layers = svg.findall(".//*[@data-layer='pore-perspective']/{*}use")
            scales = [float(re.search(r'scale\(([^)]+)', el.get('transform'))[1]) for el in layers]
            assert len(scales) == 5 and all(a < b for a,b in zip(scales, scales[1:])) and scales[-1] == 1
            assert Image.open(out/'layout-fixture.png').getpixel((0, 0))[-1] == 0
        if style == 'red_blue':
            svg = ET.parse(out/'layout-fixture.svg').getroot()
            joins = [p for p in svg.iter() if p.get('data-vertex-atom') is not None]
            assert joins
            assert all('#000000' not in p.get('style', '') for p in joins)
    assert all(content == graph_exports[0] for content in graph_exports)


@pytest.mark.parametrize('palette', [['red', 'blue'], ['#123456'], ['#123456', 'url(javascript:x)'], ['#123456']*9])
def test_custom_palette_rejects_invalid_or_non_color_values(palette):
    with pytest.raises(ValueError, match='custom_palette'):
        style_configuration('red_blue', options={'custom_palette': palette})


def test_custom_palette_is_used_in_actual_svg_and_transparent_highlights(tmp_path):
    recipe = pore()
    render_styled(recipe, tmp_path, 650, 'node_highlight', style_options={
        'custom_palette': ['#123456', '#abcdef'], 'highlight_width': 1.6, 'transparent_background': True})
    svg = (tmp_path/'layout-fixture.svg').read_text()
    assert '#123456' in svg and '#abcdef' in svg and 'stroke-linejoin="round"' in svg
    assert '<rect width="100%"' not in svg


def test_coordinate_override_requires_exact_graph_keys():
    recipe = pore(); mol, audit = build_patch(recipe)
    recipe['atom_coordinates'] = [{**o, 'xy': [0, 0]} for o in audit['atom_origin']]
    recipe['atom_coordinates'][0]['fragment_atom_index'] = 999
    with pytest.raises(ValueError, match='keys must match'):
        build_patch(recipe)
