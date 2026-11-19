"""Regression checks for composition evidence identity and bonded repeat boundaries."""
import json
from pathlib import Path

from cof_vision.cofdraw_chemistry import design_signature, prepare_template
from cof_vision.cofdraw_repeat import repeated_drawing
from cof_vision.framework_drawing import VERSION


def test_space_snapshot_resolves_every_cell_and_model():
    root = Path(__file__).resolve().parents[1]
    snapshot = json.loads((root/'src/cof_vision/web/chemical-space.json').read_text('utf8'))
    ids = set(map(int, snapshot['blocks']))
    stride = snapshot['stride']
    assert max(ids) < stride
    keys = set()
    for point in snapshot['points']:
        key = point['f'], tuple(point['ids'])
        assert point['f'] in snapshot['families'] and key not in keys
        assert key[1] == tuple(sorted(set(key[1]))) and ids.issuperset(key[1])
        assert point['refs'] and all(doi.startswith('10.') for _,doi in point['refs'])
        keys.add(key)
    for family in snapshot['families'].values():
        pairs = family['eligible']
        assert pairs == sorted(set(pairs))
        for code in pairs:
            a,b = divmod(code,stride)
            assert a in ids and b in ids and a < b
    manifest = json.loads((root/'data/crystal-models/manifest.json').read_text('utf8'))
    models = {m['id']:m['sha256'] for m in manifest['curated_models']+manifest['reported_models']}
    assert {r['id']:r['sha256'] for r in snapshot['properties']} == models


def test_honeycomb_repeat_joins_image_edges_and_keeps_only_external_ports():
    recipe = prepare_template('O=Cc1cc(C=O)cc(C=O)c1','Nc1ccc(N)cc1','imine','A','B',12,'hcb_3_2')
    recipe.update(schema_version=VERSION, recipe_id='repeat-test')
    recipe['design_template_signature'] = design_signature(recipe)
    for n in (2,3):
        svg, summary = repeated_drawing(recipe,n,600)
        assert '<svg' in svg
        assert summary['blocks'] == 5*n*n
        # Three inter-sublattice edges per cell, two cross lattice boundaries.
        assert summary['continuation_sites'] == 6*n-2
