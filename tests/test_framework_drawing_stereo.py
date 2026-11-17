"""MOL output must retain chemical stereo independently of display geometry."""
import math
from pathlib import Path

from rdkit import Chem
import pytest

from cof_vision.framework_drawing import VERSION, build_patch, linkage_fragment, render_recipe


def stereo_pore():
    # Three aldehyde nodes, each with a defined tetrahedral pendant and an E
    # fluoroalkene, alternate with three triamine nodes. Inter-node imines have
    # intentionally unspecified E/Z: rendering may not invent defined geometry.
    aldehyde = linkage_fragment(
        'O=Cc1c([C@H](F)Cl)c(C=O)c(/C=C/F)c(C=O)c1',
        'imine', 'aryl_aldehyde',
    )
    amine = linkage_fragment('Nc1cc(N)cc(N)c1', 'imine', 'primary_arylamine')
    nodes = [
        {'id': str(i), 'component': 'aldehyde' if i % 2 else 'amine',
         'incoming_port': 1, 'outgoing_port': 2,
         'layout_position': [math.cos(i * math.pi / 3), math.sin(i * math.pi / 3)]}
        for i in range(6)
    ]
    return {
        'schema_version': VERSION, 'recipe_id': 'stereo-export-test',
        'components': {'aldehyde': aldehyde, 'amine': amine},
        'instances': nodes,
        'connections': [
            {'a': {'instance': str(i), 'port': 2},
             'b': {'instance': str((i + 1) % 6), 'port': 1}, 'bond_order': 2}
            for i in range(6)
        ],
        'designed_variant': True,
    }


def _defined_stereo(mol):
    tetrahedra = [a for a in mol.GetAtoms()
                  if a.GetChiralTag() != Chem.ChiralType.CHI_UNSPECIFIED]
    alkenes = [b for b in mol.GetBonds() if b.GetStereo() in
               (Chem.BondStereo.STEREOE, Chem.BondStereo.STEREOZ,
                Chem.BondStereo.STEREOCIS, Chem.BondStereo.STEREOTRANS)]
    return len(tetrahedra), len(alkenes)


def _same_aromaticity_key(mol):
    # The production builder intentionally retains component aromaticity only.
    # A standard MOL reader perceives global macrocycle aromaticity. Normalize
    # BOTH small test graphs alike before comparing their isomeric graph keys;
    # aromatic resonance flags are not a change of atom/bond stereochemistry.
    normalized = Chem.Mol(mol)
    Chem.SanitizeMol(normalized)
    return Chem.MolToSmiles(normalized, isomericSmiles=True)


def test_explicit_stereo_survives_mol_export(tmp_path):
    recipe = stereo_pore()
    original, _ = build_patch(recipe)
    assert _defined_stereo(original) == (3, 3)
    audit = render_recipe(recipe, tmp_path, size=700)
    restored = Chem.MolFromMolFile(audit['artifacts']['mol'], removeHs=False)
    assert restored is not None
    assert _defined_stereo(restored) == (3, 3)
    assert _same_aromaticity_key(restored) == _same_aromaticity_key(original)
    # Both output metadata and its actual coordinates identify a planar drawing.
    assert not restored.GetConformer().Is3D()
    assert '2D' in Path(audit['artifacts']['mol']).read_text(encoding='utf-8').splitlines()[1]


def test_display_orientation_retains_explicit_stereo(tmp_path):
    first = stereo_pore()
    reflected_positions = stereo_pore()
    for node in reflected_positions['instances']:
        node['layout_position'][0] *= -1
    first_mol, _ = build_patch(first)
    second_mol, _ = build_patch(reflected_positions)
    assert Chem.MolToSmiles(first_mol, isomericSmiles=True) == Chem.MolToSmiles(
        second_mol, isomericSmiles=True
    )
    # Reversing the global pore arrangement does not authorize reflecting each
    # stereodefined molecular fragment or silently changing a double bond.
    audit = render_recipe(reflected_positions, tmp_path, size=700)
    restored = Chem.MolFromMolFile(audit['artifacts']['mol'], removeHs=False)
    assert _defined_stereo(restored) == (3, 3)
    assert _same_aromaticity_key(restored) == _same_aromaticity_key(first_mol)


def consumed_selector_pore():
    """Regression derived from source-visible alternating cyanoalkene/imine ends.

    Both terminal dummy atoms of every bridge select explicit double-bond
    stereochemistry. The six three-connected nodes have alternating pairs of
    alkene and imine attachments, so a simple uniform linker direction is wrong.
    """
    nodes = []
    for i in range(12):
        is_node = i % 2 == 0
        j = i // 2
        incoming, outgoing = (1, 2) if is_node or j % 2 == 0 else (2, 1)
        angle = 2 * math.pi / 3 - i * math.pi / 6
        nodes.append({'id': ('T' if is_node else 'L') + str(j),
                      'component': 'T' if is_node else 'L',
                      'incoming_port': incoming, 'outgoing_port': outgoing,
                      'layout_position': [math.cos(angle), math.sin(angle)]})
    return {
        'schema_version': VERSION, 'recipe_id': 'consumed-selector-regression',
        'components': {
            'T': {'fragment_smiles': '[*:1]c1cc([*:2])cc([*:3])c1'},
            'L': {'fragment_smiles': r'[*:1]/C=C(\C#N)c1ccc(/N=C/[*:2])cc1'},
        },
        'instances': nodes,
        'connections': [
            {'a': {'instance': nodes[i]['id'], 'port': nodes[i]['outgoing_port']},
             'b': {'instance': nodes[(i + 1) % 12]['id'],
                   'port': nodes[(i + 1) % 12]['incoming_port']}, 'bond_order': 1}
            for i in range(12)
        ],
        'layout_refinement': 'rigid_block_translation_v1',
        'designed_variant': True,
    }


def test_consumed_dummy_selectors_retain_alkene_and_imine_stereo():
    mol, audit = build_patch(consumed_selector_pore())
    defined = [b for b in mol.GetBonds() if b.GetStereo() in
               (Chem.BondStereo.STEREOE, Chem.BondStereo.STEREOZ)]
    assert len(defined) == 12
    assert sum(b.GetStereo() == Chem.BondStereo.STEREOE for b in defined) == 6
    assert sum(b.GetStereo() == Chem.BondStereo.STEREOZ for b in defined) == 6
    for b in defined:
        selectors = tuple(b.GetStereoAtoms())
        assert len(selectors) == 2
        assert mol.GetBondBetweenAtoms(b.GetBeginAtomIdx(), selectors[0]) is not None
        assert mol.GetBondBetweenAtoms(b.GetEndAtomIdx(), selectors[1]) is not None
        assert all(mol.GetAtomWithIdx(i).GetAtomicNum() != 0 for i in selectors)
    # Actual stereochemical serialization must retain all twelve definitions;
    # checking only enums would miss RDKit's empty-selector failure.
    restored = Chem.MolFromSmiles(audit['fragment_smiles'], sanitize=False)
    restored.UpdatePropertyCache(strict=True)
    Chem.FastFindRings(restored)
    Chem.SetBondStereoFromDirections(restored)
    Chem.AssignStereochemistry(restored, cleanIt=False, force=True)
    assert _defined_stereo(restored) == (0, 12)
    assert audit['double_bond_stereo_selector_remapping'] == {
        'defined_bonds': 12, 'consumed_selectors_rebound': 12}


def test_export_rejects_layout_that_would_invert_defined_imines(tmp_path):
    recipe = consumed_selector_pore()
    mol, audit = build_patch(recipe)
    # The rigid block layout gives three joins the wrong selector-side geometry.
    # It is a valid chemical graph, but exporting those 2D coordinates as a MOL
    # would silently change E to Z. Reject before any misleading artifact exists.
    assert audit['defined_double_bond_geometry']['checked'] == 12
    assert not audit['defined_double_bond_geometry']['valid']
    assert len(audit['defined_double_bond_geometry']['violations']) == 3
    with pytest.raises(ValueError, match='contradicts explicit double-bond'):
        render_recipe(recipe, tmp_path, size=400)
    assert not list(tmp_path.glob('*.mol'))
    assert not list(tmp_path.glob('*.png'))
