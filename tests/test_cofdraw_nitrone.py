"""Bounded nitrone graph construction, not a crystal or synthesis prediction.

The tetratopic precursor below is the finite graph in Scheme 1 of
10.1002/anie.202307674. Its use in a chosen drawing template does not identify
that drawing as the paper's COF-115 or assign its experimentally reported net.
"""
from copy import deepcopy
from pathlib import Path
import xml.etree.ElementTree as ET

from fastapi.testclient import TestClient
import pytest
from rdkit import Chem

from cof_vision.cofdraw_chemistry import decorated_monomers, inspect_monomer, prepare_template
from cof_vision.cofdraw_core import export_chemical_formats, render_styled
from cof_vision.cofdraw_conversation import describe_turn
from cof_vision.cofdraw_describe import describe_cof
from cof_vision.cofdraw_desktop import create_desktop_app
from cof_vision.framework_drawing import VERSION, build_patch, linkage_fragment

TA = 'O=Cc1ccc(C=O)cc1'
DHA = 'ONc1ccc(NO)cc1'
TPE_HYDROXYLAMINE = 'ONc1ccc(C(=C(c2ccc(NO)cc2)c2ccc(NO)cc2)c2ccc(NO)cc2)cc1'


def design(b=DHA, template='chain_2_2'):
    return dict(prepare_template(TA, b, 'nitrone', 'TA', 'Hydroxylamine', 6, template),
                schema_version=VERSION, recipe_id='nitrone-regression', designed_variant=True,
                layout_refinement='bond_angle_pose_v2')


def assert_nitrone_graph(mol, n_sites, n_joins):
    mol.UpdatePropertyCache(strict=True)
    nitrogens = [a for a in mol.GetAtoms() if a.GetAtomicNum() == 7]
    oxygens = [a for a in mol.GetAtoms() if a.GetAtomicNum() == 8]
    assert len(nitrogens) == len(oxygens) == n_sites
    assert sum(a.GetFormalCharge() for a in mol.GetAtoms()) == 0
    for nitrogen in nitrogens:
        assert nitrogen.GetFormalCharge() == 1
        assert nitrogen.GetTotalNumHs(includeNeighbors=True) == 0
        attached_o = [a for a in nitrogen.GetNeighbors() if a.GetAtomicNum() == 8]
        assert len(attached_o) == 1
        oxygen = attached_o[0]
        assert oxygen.GetFormalCharge() == -1 and oxygen.GetDegree() == 1
        assert oxygen.GetTotalNumHs(includeNeighbors=True) == 0
        assert mol.GetBondBetweenAtoms(nitrogen.GetIdx(), oxygen.GetIdx()).GetBondType() == Chem.BondType.SINGLE
    joins = [b for b in mol.GetBonds() if b.GetBondType() == Chem.BondType.DOUBLE
             and sorted([b.GetBeginAtom().GetAtomicNum(), b.GetEndAtom().GetAtomicNum()]) == [6, 7]]
    assert len(joins) == n_joins
    for bond in joins:
        carbon = next(a for a in (bond.GetBeginAtom(), bond.GetEndAtom()) if a.GetAtomicNum() == 6)
        assert carbon.GetTotalNumHs() == 1
        assert bond.GetStereo() in (Chem.BondStereo.STEREONONE, Chem.BondStereo.STEREOANY)
    return joins


def test_only_hydroxylamine_n_reacts_and_n_o_bonds_are_retained():
    fragment = linkage_fragment(DHA, 'nitrone', 'aryl_hydroxylamine')
    assert fragment['connectivity'] == 2
    mol = Chem.MolFromSmiles(fragment['fragment_smiles'])
    assert_nitrone_graph(mol, 2, 0)
    for port in (a for a in mol.GetAtoms() if a.GetAtomicNum() == 0):
        nitrogen = port.GetNeighbors()[0]
        assert nitrogen.GetAtomicNum() == 7
        assert mol.GetBondBetweenAtoms(port.GetIdx(), nitrogen.GetIdx()).GetBondType() == Chem.BondType.DOUBLE
    roles = inspect_monomer(DHA)['compatible_roles']
    assert [(r['linkage'], r['role'], r['reactive_sites']) for r in roles] == [('nitrone', 'aryl_hydroxylamine', 2)]


@pytest.mark.parametrize('smiles', [
    'Nc1ccc(N)cc1',             # primary arylamine: no N-O bond
    'NOc1ccc(ON)cc1',           # aryl-O-NH2 has the wrong attachment atom
    'CONc1ccc(NOC)cc1',         # O-methyl protected hydroxylamine
    'CN(O)c1ccc(N(C)O)cc1',     # no N-H: N is already disubstituted
    '[O-][NH2+]c1ccccc1',       # zwitterionic starting-state is not the neutral rule
    'O[NH2+]c1ccccc1',          # protonated hydroxylamine
    'ONCCNO',                  # alkyl hydroxylamine: outside this aromatic rule
    '[2H]ONc1ccccc1',           # O-D cannot silently disappear
    '[2H]N(O)c1ccccc1',         # N-D cannot silently disappear
])
def test_related_or_protected_groups_cannot_masquerade_as_supported_precursors(smiles):
    assert not any(r['role'] == 'aryl_hydroxylamine' for r in inspect_monomer(smiles)['compatible_roles'])
    with pytest.raises(ValueError, match='No supported reactive sites'):
        linkage_fragment(smiles, 'nitrone', 'aryl_hydroxylamine')


def test_nitrogen_and_oxygen_isotopes_survive_reaction():
    fragment = linkage_fragment('[18OH][15NH]c1ccccc1', 'nitrone', 'aryl_hydroxylamine')
    mol = Chem.MolFromSmiles(fragment['fragment_smiles'])
    assert_nitrone_graph(mol, 1, 0)
    assert sorted((a.GetAtomicNum(), a.GetIsotope()) for a in mol.GetAtoms() if a.GetIsotope()) == [(7, 15), (8, 18)]


@pytest.mark.parametrize('b,template,sites,joins,ports', [
    (DHA, 'chain_2_2', 6, 5, 2),
    (TPE_HYDROXYLAMINE, 'sql_4_2', 16, 8, 8),
])
def test_connected_patch_removes_only_aldehyde_oxygen_and_preserves_precursor_bonds(b, template, sites, joins, ports):
    recipe = design(b, template)
    mol, audit = build_patch(recipe)
    assert len(Chem.GetMolFrags(mol)) == 1
    assert_nitrone_graph(mol, sites, joins)
    assert len(audit['continuation_ports']) == ports
    assert all(audit['atom_origin'][a.GetIdx()]['component'] == 'B' for a in mol.GetAtoms() if a.GetAtomicNum() == 8)
    assert audit['whole_cof_smiles'] is None
    assert audit['ground_truth'] is False
    assert recipe.get('cofex_id') is None and recipe.get('source') is None
    assert 'neither is experimental identity' in recipe['design_identity_qualification']
    origins = {(row['instance'], row['fragment_atom_index']): i for i, row in enumerate(audit['atom_origin'])}
    for instance in recipe['instances']:
        fragment = Chem.MolFromSmiles(recipe['components'][instance['component']]['fragment_smiles'])
        for bond in fragment.GetBonds():
            if not bond.GetBeginAtom().GetAtomicNum() or not bond.GetEndAtom().GetAtomicNum():
                continue
            a, b_index = (origins[(instance['id'], index)] for index in (bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()))
            assert mol.GetBondBetweenAtoms(a, b_index).GetBondType() == bond.GetBondType()


def test_mol_sdf_and_cdxml_retain_charges_and_do_not_invent_e_z(tmp_path):
    recipe = design()
    mol, _ = build_patch(recipe)
    render_styled(recipe, tmp_path, 700, 'standard')
    export_chemical_formats(tmp_path, recipe['recipe_id'], source_mol=mol)
    for suffix in ('mol', 'sdf'):
        block = (tmp_path / f'{recipe["recipe_id"]}.{suffix}').read_text(encoding='utf-8').split('$$$$')[0]
        restored = Chem.MolFromMolBlock(block, removeHs=False)
        assert restored is not None
        assert_nitrone_graph(restored, 6, 5)
        assert len(Chem.GetMolFrags(restored)) == 1
    if Chem.HasChemDrawCDXSupport():
        cdxml = ET.parse(tmp_path / f'{recipe["recipe_id"]}.cdxml')
        charges = [int(atom.get('Charge', '0')) for atom in cdxml.iter('n')]
        assert charges.count(1) == charges.count(-1) == 6


def test_planner_infers_nitrone_only_for_its_supported_roles_and_validates_graph():
    with TestClient(create_desktop_app(Path(__file__).resolve().parents[1] / 'data')) as client:
        library = deepcopy(client.get('/api/drawing/library').json())
    library['monomer_presets'].extend(decorated_monomers([
        {'id':'nitrone-test-aldehyde', 'label':'QA-dialdehyde', 'abbreviation':'QA-dialdehyde', 'smiles':TA},
        {'id':'nitrone-test-hydroxylamine', 'label':'QA-hydroxylamine', 'abbreviation':'QA-hydroxylamine', 'smiles':DHA},
    ]))
    for text in ('Draw QA-dialdehyde and QA-hydroxylamine', 'Draw QA-hydroxylamine and QA-dialdehyde with nitrone linkage'):
        plan = describe_cof(text, library)
        assert plan['status'] == 'ready', plan['unresolved']
        assert plan['prepare_request']['linkage'] == 'nitrone'
        assert plan['prepare_request']['monomer_a'] == TA
        assert plan['prepare_request']['monomer_b'] == DHA
        assert plan['source_recipe_id'] is None
        prepared = prepare_template(**plan['prepare_request'])
        prepared.update(schema_version=VERSION, recipe_id='planned-nitrone', designed_variant=True)
        graph, _ = build_patch(prepared)
        assert len(Chem.GetMolFrags(graph)) == 1
    wrong = describe_cof('Draw TFB and TAPB with nitrone linkage', library)
    assert wrong['status'] != 'ready' and wrong['prepare_request'] is None
    count = describe_turn('Draw QA-dialdehyde and QA-hydroxylamine with ditopic hydroxylamine', library)
    assert count['status'] == 'ready', count['unresolved']
    assert 'hydroxylamine with 2 reactive sites' in count['conversation']['context']['requirements']
    styled = describe_turn('Make it pastel', library, count['conversation']['context'])
    assert styled['status'] == 'ready' and styled['prepare_request']['linkage'] == 'nitrone'
    mismatch = describe_cof('Draw QA-dialdehyde and QA-hydroxylamine with tetratopic hydroxylamine', library)
    assert mismatch['status'] != 'ready'
    assert any('4-site hydroxylamine' in line for line in mismatch['unresolved'])


def test_actual_api_prepares_and_renders_qualified_nitrone_design(tmp_path):
    with TestClient(create_desktop_app(tmp_path)) as client:
        response = client.post('/api/drawing/prepare', json={
            'monomer_a':TA, 'monomer_b':TPE_HYDROXYLAMINE,
            'linkage':'nitrone', 'template_id':'sql_4_2', 'label_a':'TA', 'label_b':'TPE-hydroxylamine',
        })
        assert response.status_code == 200, response.text
        prepared = response.json()
        assert prepared['compatibility']['reactive_site_counts'] == [2, 4]
        assert 'E/Z' in prepared['limits'] and 'not assigned' in prepared['limits']
        rendered = client.post('/api/drawing/render', json={'recipe':prepared['recipe'], 'style_id':'standard', 'size':700})
        assert rendered.status_code == 200, rendered.text
        payload = rendered.json()
        assert payload['status'] == 'designed_variant'
        assert payload['audit']['registered_recipe_sha256'] is None
        assert payload['audit']['cofex_string'] is None
        assert payload['audit']['whole_cof_smiles'] is None
        assert {'svg', 'png', 'mol', 'sdf', 'json'} <= set(payload['downloads'])
        exported = client.get(payload['downloads']['mol'])
        assert exported.status_code == 200
        restored = Chem.MolFromMolBlock(exported.text, removeHs=False)
        assert_nitrone_graph(restored, 16, 8)
        wrong = client.post('/api/drawing/prepare', json={'monomer_a':TA, 'monomer_b':'Nc1ccc(N)cc1', 'linkage':'nitrone'})
        assert wrong.status_code == 422
