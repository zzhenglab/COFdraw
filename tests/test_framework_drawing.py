import copy
import json
import math

import pytest
import numpy as np
from rdkit import Chem

from cof_vision.framework_drawing import VERSION, TP, linkage_fragment, build_patch, draw_cofex
from cof_vision.cofdraw_depiction import compute_2d_coords, normalize_linkage_ports


def test_reduced_linkage_symmetry_preserves_individual_drawing_angles():
    fragment = linkage_fragment('O=Cc1cc(C=O)cc(C=O)c1', 'secondary_amine', 'aryl_aldehyde')
    mol = Chem.MolFromSmiles(fragment['fragment_smiles'])
    compute_2d_coords(mol)
    xy = mol.GetConformer().GetPositions()[:, :2].copy()
    scaffold = [a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum()]
    original = xy[scaffold].copy()
    ports = [a for a in mol.GetAtoms() if not a.GetAtomicNum()]
    expected = []
    for port, angle, length in zip(ports, (40., -73., 54.), (1.4, 1.5, 1.6), strict=True):
        site = port.GetNeighbors()[0]
        anchor = next(a for a in site.GetNeighbors() if a.GetAtomicNum())
        axis = xy[site.GetIdx()] - xy[anchor.GetIdx()]
        axis /= np.linalg.norm(axis)
        normal = np.array([-axis[1], axis[0]])
        theta = math.radians(angle)
        xy[port.GetIdx()] = xy[site.GetIdx()] + length * (math.cos(theta) * axis + math.sin(theta) * normal)
        expected.append((port.GetIdx(), site.GetIdx(), axis, abs(angle), length))

    assert normalize_linkage_ports(mol, xy, 'secondary_amine')
    np.testing.assert_array_equal(xy[scaffold], original)
    for port, site, axis, angle, length in expected:
        bond = xy[port] - xy[site]
        assert np.linalg.norm(bond) == pytest.approx(length)
        assert math.degrees(math.atan2(np.linalg.det([axis, bond]), axis @ bond)) == pytest.approx(angle)


def pore():
    a = linkage_fragment('O=Cc1cc(C=O)cc(C=O)c1', 'imine', 'aryl_aldehyde')
    b = linkage_fragment('Nc1cc(N)cc(N)c1', 'imine', 'primary_arylamine')
    nodes = [{'id': str(i), 'component': 'a' if i%2 else 'b',
              'incoming_port': 1, 'outgoing_port': 2,
              'layout_position': [math.cos(i*math.pi/3), math.sin(i*math.pi/3)]} for i in range(6)]
    edges = [{'a': {'instance': str(i), 'port': 2}, 'b': {'instance': str((i+1)%6), 'port': 1}, 'bond_order': 2} for i in range(6)]
    return {'schema_version': VERSION, 'recipe_id': 'test', 'components': {'a': a, 'b': b},
            'instances': nodes, 'connections': edges, 'designed_variant': True}


def test_linkage_rules_preserve_expected_bond_chemistry():
    tp = linkage_fragment(TP, 'beta_ketoenamine', 'tp_trialdehyde')
    mol = Chem.MolFromSmiles(tp['fragment_smiles'])
    assert tp['connectivity'] == 3
    assert len(mol.GetSubstructMatches(Chem.MolFromSmarts('[C]=O'))) == 3
    assert len(mol.GetSubstructMatches(Chem.MolFromSmarts('[C]=[CH][#0]'))) == 3
    assert not any(a.GetIsAromatic() for a in mol.GetAtoms())
    amine = linkage_fragment('Nc1ccc(N)cc1', 'beta_ketoenamine', 'primary_arylamine')
    assert all(a.GetTotalNumHs() == 1 for a in Chem.MolFromSmiles(amine['fragment_smiles']).GetAtoms() if a.GetAtomicNum()==7)
    with pytest.raises(ValueError, match='exact Tp'):
        linkage_fragment('O=Cc1cc(C=O)cc(C=O)c1', 'beta_ketoenamine', 'tp_trialdehyde')


def test_pore_is_connected_with_open_ports_and_not_whole_cof_smiles():
    mol, audit = build_patch(pore())
    assert len(Chem.GetMolFrags(mol)) == 1
    assert len(audit['continuation_ports']) == 6
    assert audit['connections'] == 6 and audit['whole_cof_smiles'] is None
    assert len(mol.GetSubstructMatches(Chem.MolFromSmarts('[C]=[N]'))) == 6
    assert audit['graph_scope'].startswith('finite_designed')
    assert not audit['synthetic_training_admitted']


def test_double_use_and_bond_mismatch_are_rejected():
    r = pore(); r['connections'].append(copy.deepcopy(r['connections'][0]))
    with pytest.raises(ValueError, match='consumed more than once'):
        build_patch(r)
    r = pore(); r['connections'][0]['bond_order'] = 1
    with pytest.raises(ValueError, match='bond order disagrees'):
        build_patch(r)
    r = pore(); r['topology_basis'] = 'monomer_degree_only'
    with pytest.raises(ValueError, match='degree alone'):
        build_patch(r)


def test_registry_rejects_unknown_and_tampered_recipe(tmp_path):
    registry = tmp_path/'registry.json'
    registry.write_text(json.dumps({'recipes': []}))
    with pytest.raises(ValueError, match='found 0'):
        draw_cofex('missing', registry, tmp_path/'out')
    p = tmp_path/'recipe.json'; p.write_text(json.dumps(pore()))
    registry.write_text(json.dumps({'recipes': [{'cofex_id': 'known', 'recipe_id': 'test', 'path': p.name, 'sha256': 'bad'}]}))
    with pytest.raises(ValueError, match='hash mismatch'):
        draw_cofex('known', registry, tmp_path/'out')


def test_spacing_refinement_preserves_block_shapes_and_chemical_graph():
    recipe = pore()
    # A large canvas radius formerly stretched only the interblock bonds.
    for item in recipe['instances']:
        item['layout_position'] = [20*x for x in item['layout_position']]
    old, before = build_patch(recipe)
    recipe['layout_refinement'] = 'rigid_block_translation_v1'
    new, after = build_patch(recipe)
    assert after['layout_refinement']['applied']
    assert before['fragment_serialization_sha256'] == after['fragment_serialization_sha256']
    assert before['continuation_ports'] == after['continuation_ports']
    assert after['depiction_geometry']['max_to_median_ratio'] < 1.001
    a = np.asarray(old.GetConformer().GetPositions())
    b = np.asarray(new.GetConformer().GetPositions())
    for item in recipe['instances']:
        indices = [i for i, o in enumerate(after['atom_origin']) if o['instance'] == item['id']]
        # All within-block vectors, including substituents and ports, stay fixed.
        assert np.allclose(a[indices]-a[indices[0]], b[indices]-b[indices[0]], atol=1e-10)
    recipe['layout_refinement'] = 'unknown'
    with pytest.raises(ValueError, match='Unsupported layout refinement'):
        build_patch(recipe)


def test_reused_patch_cannot_inherit_a_different_recipe_identity(tmp_path):
    from cof_vision.framework_drawing import render_recipe
    recipe = pore()
    prepared = build_patch(recipe)
    edited = copy.deepcopy(recipe)
    edited['instances'][0]['layout_position'][0] += 1
    with pytest.raises(ValueError, match='does not match'):
        render_recipe(edited, tmp_path, prepared_patch=prepared)
