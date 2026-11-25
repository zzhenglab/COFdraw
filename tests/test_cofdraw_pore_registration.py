"""A change in drawing extent preserves the molecule at each corner of the shared face."""
import json
from pathlib import Path

import numpy as np
import pytest
from rdkit import Chem
from rdkit.Chem import rdMolTransforms

from cof_vision.cofdraw_chemistry import design_signature, prepare_template
from cof_vision.cofdraw_explicit import clear_source_registration
from cof_vision.cofdraw_periodic_options import validated_design
from cof_vision.cofdraw_repeat import repeat_patch
from cof_vision.framework_drawing import VERSION, build_patch


ROOT = Path(__file__).resolve().parents[1]


def sql_design(aldehyde):
    catalog = json.loads((ROOT / 'data/framework-drawing/catalog/current.json').read_bytes())
    reviewed_feeds = json.loads((ROOT / 'src/cof_vision/web/cofdraw-reviewed-feeds.json').read_bytes())
    smiles = {r['id']: r['smiles'] for r in [*catalog['monomer_presets'], *reviewed_feeds]}
    recipe = prepare_template(smiles[aldehyde], smiles['pytta'], 'imine', aldehyde, 'PyTTA', 4, 'sql_4_4')
    recipe.update(schema_version=VERSION, recipe_id='sql-test', design_template_signature=design_signature(recipe))
    return recipe


def face_atoms_and_joins(mol, origins, face):
    """Identify atoms by their corner on this rim and their atom in the original fragment."""
    keys = {atom.GetIdx(): (face.index(owner['instance']), owner['fragment_atom_index'])
            for atom, owner in zip(mol.GetAtoms(), origins)
            if atom.GetAtomicNum() and owner['instance'] in face}
    xyz = mol.GetConformer().GetPositions()[:, :2]
    positions = {key: xyz[idx] for idx, key in keys.items()}
    bonds = set()
    for bond in mol.GetBonds():
        i, j = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        if i in keys and j in keys and keys[i][0] != keys[j][0]:
            bonds.add(frozenset((keys[i], keys[j])))
    return positions, bonds


@pytest.mark.parametrize('aldehyde', ['feed_97d130d567b8', 'tfppy'], ids=['BDPT', 'TFPy'])
def test_single_pore_reuses_the_first_square_of_the_fitted_sheet(aldehyde):
    recipe = sql_design(aldehyde)
    single, audit = build_patch(recipe)
    sheet, owners, _, summary, _ = repeat_patch(recipe, 2)
    assert summary['pore_count'] == 4
    rim = [row['id'] for row in recipe['instances']]
    one, one_joins = face_atoms_and_joins(single, audit['atom_origin'], rim)
    shared, shared_joins = face_atoms_and_joins(sheet, owners, summary['pore_faces'][0])
    assert len(one_joins) == 4 and one_joins == shared_joins
    assert one.keys() == shared.keys()
    # Removing the neighboring pores must keep the pose of each labelled
    # monomer. Centering permits recentering the canvas on the smaller crop.
    keys = sorted(one)
    a, b = (np.stack([points[k] for k in keys]) for points in (one, shared))
    np.testing.assert_allclose(a - a.mean(axis=0), b - b.mean(axis=0), atol=1e-5)
    if aldehyde == 'feed_97d130d567b8':
        imines = single.GetSubstructMatches(Chem.MolFromSmarts('c-[C]=[N]-c'))
        assert len(imines) == 4
        for ring_c, carbon, nitrogen, ring_n in imines:
            conf = single.GetConformer()
            assert rdMolTransforms.GetAngleDeg(conf, ring_c, carbon, nitrogen) == pytest.approx(120, abs=.1)
            assert rdMolTransforms.GetAngleDeg(conf, carbon, nitrogen, ring_n) == pytest.approx(120, abs=.1)


def test_saved_bdpt_pytta_incidence_and_hand_edited_coordinates_survive():
    recipe = sql_design('feed_97d130d567b8')
    # A project drawn before the shared-face layout used the same directed
    # port pair at both occurrences of each monomer. Preserve that incidence.
    ring = [('A1', 1, 3), ('B2', 2, 3), ('A3', 1, 3), ('B4', 2, 3)]
    old_positions = [[0, 25], [-25, 0], [0, -25], [25, 0]]
    recipe['instances'] = [dict(id=k, component=k[0], incoming_port=p, outgoing_port=q, layout_position=xy)
                           for (k, p, q), xy in zip(ring, old_positions)]
    recipe['connections'] = [dict(a=dict(instance=a, port=q), b=dict(instance=b, port=p), bond_order=2)
                             for (a, _, q), (b, p, _) in zip(ring, ring[1:] + ring[:1])]
    recipe.pop('coordinate_provenance', None)
    recipe.pop('atom_coordinates', None)
    recipe['design_template_signature'] = design_signature(recipe)
    legacy, audit = build_patch(recipe)
    edited = legacy.GetConformer().GetPositions()[:, :2].copy()
    edited[0] += [.3, -.2]
    recipe['atom_coordinates'] = [dict(instance=o['instance'], fragment_atom_index=o['fragment_atom_index'],
                                      xy=point.tolist()) for o, point in zip(audit['atom_origin'], edited)]
    restored, _ = validated_design(recipe)
    assert restored['connections'] == recipe['connections']
    reopened, _ = build_patch(recipe)
    np.testing.assert_allclose(reopened.GetConformer().GetPositions()[:, :2], edited)
    assert repeat_patch(recipe, 2)[3]['pore_count'] == 4


def test_sheet_geometry_is_recognized_on_reopen_and_after_a_manual_edit():
    recipe = json.loads(json.dumps(sql_design('feed_97d130d567b8')))
    clear_source_registration(recipe)
    assert recipe['coordinate_provenance']['kind'] == 'fitted_sheet_face'
    recipe['atom_coordinates'][0]['xy'][0] += .25
    clear_source_registration(recipe)
    assert recipe['coordinate_provenance']['kind'] == 'user_manual_coordinates'
    _, audit = build_patch(recipe)
    assert audit['layout_refinement']['source_coordinates_preserved'] is False
