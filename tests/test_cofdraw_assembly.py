from pathlib import Path

import numpy as np
import pytest
from rdkit import Chem

from cof_vision.cofdraw_assembly import assemble
from cof_vision.cofdraw_chemistry import decorated_monomers
from cof_vision.cofdraw_core import catalog_monomers, read_catalog
from cof_vision.cofdraw_mixed import MIXED_EXAMPLES, canonical_sequence, mixed_pore_patterns, pore_patch
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


@pytest.mark.parametrize('example',MIXED_EXAMPLES,ids=lambda e:e['id'])
def test_mixed_examples_are_connected_edge_sharing_molecular_patches(catalog,example):
    rows=[catalog[k] for k in example['monomer_ids']]
    pattern=example['linker_sequence']
    layout=example.get('pore_layout')
    recipe=assemble(rows,pattern,pore_layout=layout)
    mol,summary=pore_patch(recipe['components'],pattern,2,2,layout)
    origins=summary['atom_origin']
    atom_blocks=[row['instance'] for row in origins]
    xy=mol.GetConformer().GetPositions()[:,:2]
    joins=set()
    for bond in mol.GetBonds():
        a,b=bond.GetBeginAtomIdx(),bond.GetEndAtomIdx()
        if atom_blocks[a]!=atom_blocks[b]:
            assert bond.GetBondType()==Chem.BondType.DOUBLE
            assert {mol.GetAtomWithIdx(x).GetSymbol() for x in (a,b)}=={'C','N'}
            assert 1.45 < np.linalg.norm(xy[a]-xy[b]) < 1.55
            joins.add(frozenset((atom_blocks[a],atom_blocks[b])))
    assert len(Chem.GetMolFrags(mol))==1 and len(origins)==mol.GetNumAtoms()
    assert summary['pores']==len(joins)-summary['blocks']+1
    for face,label in zip(summary['pore_faces'],summary['face_labels'],strict=True):
        assert all(frozenset((a,face[(i+1)%len(face)])) in joins for i,a in enumerate(face))
        assert canonical_sequence(''.join(s[0] for s in face[1::2]))==label
    if layout=='kgm':
        assert summary['pore_sequences']=={'BBB':3,'CCC':3,'BCBCBC':1}
        assert all(r['canonical_smiles'] for r in rows) and 'doi' in example
        triangles={}
        for blocks,label in zip(summary['pore_faces'],summary['face_labels']):
            if label=='BCBCBC':
                continue
            centers=[xy[[i for i,r in enumerate(origins) if r['instance']==node]].mean(axis=0) for node in blocks[::2]]
            a,b,c=centers
            triangles.setdefault(label,[]).append(abs(np.linalg.det(np.array([b-a,c-a])))/2)
        assert min(triangles['CCC'])>max(triangles['BBB'])
        for seq,letter in [('BBB','B'),('CCC','C')]:
            one=assemble(rows,seq,pore_layout=layout)
            assert [i['component'] for i in one['instances']]==['A',letter]*3
            ring,_=build_patch(dict(one,schema_version=VERSION,recipe_id='triangle'))
            assert len(Chem.GetMolFrags(ring))==1


def test_different_necklaces_are_selectable_at_fixed_feed_ratio(catalog):
    patterns=mixed_pore_patterns()
    assert {k:len(v) for k,v in patterns.items()}=={'4':4,'6':11}
    seqs={s['sequence'] for s in patterns['6'] if s['counts']=={'B':3,'C':3}}
    assert seqs=={'BBBCCC','BBCBCC','BCBCBC'}
    for seq in seqs:
        assert canonical_sequence(seq)==canonical_sequence(seq[1:]+seq[:1])==canonical_sequence(seq[::-1])
    rows=[catalog[k] for k in ('tapb','ta','bpda_cho')]
    with pytest.raises(ValueError,match='ETTA'):
        assemble(rows,'BBB',pore_layout='kgm')
