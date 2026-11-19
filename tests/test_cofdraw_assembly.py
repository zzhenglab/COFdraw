from pathlib import Path

import numpy as np
import pytest
from rdkit import Chem

from cof_vision.cofdraw_assembly import assemble
from cof_vision.cofdraw_chemistry import decorated_monomers
from cof_vision.cofdraw_core import catalog_monomers, read_catalog
from cof_vision.framework_drawing import VERSION, build_patch


@pytest.fixture(scope='module')
def catalog():
    rows = catalog_monomers(read_catalog(Path(__file__).resolve().parents[1]/'data'))
    return {r['id']:r for r in decorated_monomers(rows)}


@pytest.mark.parametrize('key,shape', [('bfbapy','closed_pore'), ('cv_4e35cbcca820334e','open_chain')])
def test_self_imine_has_complementary_ends(catalog,key,shape):
    recipe = assemble([catalog[key]], '')
    assert recipe['layout_kind'] == shape
    sites = {s['port']:s['site_role'] for s in recipe['components']['A']['reactive_sites']}
    for edge in recipe['connections']:
        assert (sites[edge['a']['port']],sites[edge['b']['port']]) == ('aldehyde','amine')
    mol, _ = build_patch(dict(recipe,schema_version=VERSION,recipe_id='self'))
    assert len(Chem.GetMolFrags(mol)) == 1


@pytest.mark.parametrize('key',['cv_05cbf00025a2d4a1','cv_e0efd03afb224b2a','cv_6bb0458f6411d793'])
def test_documented_masked_precursors_keep_feed_identity(catalog,key):
    row = catalog[key]
    recipe = assemble([row],'')
    comp = recipe['components']['A']
    assert comp['monomer_smiles'] == row['canonical_smiles'] and not row['binary_eligible']
    activated = Chem.MolFromSmiles(comp['activated_smiles'])
    for pattern in ('[c][CH1]=O','[c][NX3;H2;+0]'):
        assert len(activated.GetSubstructMatches(Chem.MolFromSmarts(pattern))) == 2
    build_patch(dict(recipe,schema_version=VERSION,recipe_id='masked'))


def test_mixed_linker_sequence_changes_geometry_and_bond_graph(catalog):
    rows = [catalog[k] for k in ('tapb','ta','bpda_cho')]
    def draw(pattern):
        recipe = assemble(rows,pattern)
        assert ''.join(x['component'] for x in recipe['instances'][1::2]) == pattern
        xy = np.array([x['layout_position'] for x in recipe['instances'][::2]])
        lengths = np.linalg.norm(xy-np.roll(xy,-1,axis=0),axis=1)
        assert lengths[pattern.index('C')] > lengths[pattern.index('B')]
        mol,_ = build_patch(dict(recipe,schema_version=VERSION,recipe_id='mixed'))
        return Chem.MolToSmiles(mol)
    assert draw('BCBCBC') != draw('BBBCCC')
    for order in ('BBBBBB','BCBC','BBCXCX'):
        with pytest.raises(ValueError):assemble(rows,order)
    with pytest.raises(ValueError):assemble([rows[0],rows[1],rows[1]],'BCBCBC')
    with pytest.raises(ValueError):assemble([catalog['tapm'],rows[1],rows[2]],'BCBC')
    with pytest.raises(ValueError):assemble([catalog['tapb']],'')
    mid=catalog['cv_58b48601a6c68864']
    assert mid['paper_label']=='S6' and not mid['binary_eligible']
    assert catalog['cv_3aac60550b8330ff']['paper_label']=='1'
