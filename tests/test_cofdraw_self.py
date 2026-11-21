from collections import Counter
from pathlib import Path

import numpy as np
import pytest
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

from cof_vision.cofdraw_assembly import assemble
from cof_vision.cofdraw_chemistry import decorated_monomers
from cof_vision.cofdraw_core import read_catalog,catalog_monomers
from cof_vision.cofdraw_self import self_examples
from cof_vision.cofdraw_self_kagome import kagome_patch
from cof_vision.framework_drawing import VERSION,build_patch


@pytest.fixture(scope='module')
def examples():
    rows = catalog_monomers(read_catalog(Path(__file__).resolve().parents[1]/'data'))
    keys = {e['monomer_ids'][0] for e in self_examples()}
    selected = {r['id']:r for r in decorated_monomers([r for r in rows if r['id'] in keys])}
    assert set(selected)==keys
    return selected


def build(recipe):
    return build_patch(dict(recipe,schema_version=VERSION,recipe_id='self-example',
                            designed_variant=True,layout_refinement='bond_angle_pose_v2'))


@pytest.mark.parametrize('key,counts,ring_elements',[
    ('bdba', {'C':18,'B':6,'O':6}, {'B','O'}),
    ('cv_d4982b6c8defb0a2', {'C':24,'N':6}, {'C','N'}),
])
def test_cyclotrimerization_from_one_feed(examples,key,counts,ring_elements):
    recipe=assemble([examples[key]],'')
    assert [p['smiles'] for p in recipe['assembly']['precursors']]==[examples[key]['canonical_smiles']]
    mol,audit=build(recipe)
    assert len(audit['continuation_ports'])==6 and len(Chem.GetMolFrags(mol))==1
    cell=Counter()
    for k,mult in [('A',2),('B',3)]:
        frag=Chem.MolFromSmiles(recipe['components'][k]['fragment_smiles'])
        for atom in frag.GetAtoms():
            if atom.GetAtomicNum():
                cell[atom.GetSymbol()]+=mult
        if k=='A':
            assert {a.GetSymbol() for a in frag.GetAtoms() if a.GetAtomicNum()}==ring_elements
            assert frag.GetRingInfo().NumRings()==1
            if 'B' in ring_elements:
                assert all(b.GetBondTypeAsDouble()==1 for b in frag.GetBonds())
    assert cell==counts


@pytest.mark.parametrize('key,formula',[('bfbaeb','C36H24N2O2'),('bfbaepy','C46H28N2O2')])
def test_kagome_faces_share_one_atomistic_tiling(examples,key,formula):
    row = examples[key]
    precursor=Chem.MolFromSmiles(row['smiles'])
    assert rdMolDescriptors.CalcMolFormula(precursor)==formula
    assert len(precursor.GetSubstructMatches(Chem.MolFromSmarts('[c]-C#C-[c]')))==2
    for pore,expected in [('triangle',3),('hexagon',6)]:
        recipe=assemble([row],'',pore_view=pore)
        mol,audit=build(recipe)
        assert len(recipe['connections'])==expected
        assert not audit['layout_refinement']['source_coordinates_preserved']
        assert len(audit['continuation_ports'])==expected*2
        assert recipe['assembly']['pore_view']==pore
        assert len(Chem.GetMolFrags(mol))==1
    mol,summary=kagome_patch(recipe,2)
    origins=summary['atom_origin']
    assert summary['pore_counts']=={'triangle':6,'hexagon':4}
    assert summary['blocks']==19 and len(origins)==mol.GetNumAtoms()
    xy=mol.GetConformer().GetPositions()
    for bond in mol.GetBonds():
        a,b=bond.GetBeginAtomIdx(),bond.GetEndAtomIdx()
        assert np.linalg.norm(xy[a]-xy[b])==pytest.approx(1.5,abs=.001)
        if origins[a]['instance']!=origins[b]['instance']:
            assert sorted([bond.GetBeginAtom().GetSymbol(),bond.GetEndAtom().GetSymbol()])==['C','N']
            assert bond.GetBondType()==Chem.BondType.DOUBLE


def test_pore_choice_and_feed_are_graph_specific(examples):
    pasted=decorated_monomers([dict(id='paste',label='Pasted Py linker',smiles=examples['bfbaepy']['smiles'])])[0]
    assert pasted['assembly_roles']['self_kind']=='self_kagome'
    with pytest.raises(ValueError,match='pore view'):
        assemble([examples['bfbapy']],'',pore_view='triangle')
    for smiles in ['CB(O)O','N#Cc1cc(C#N)ccc1','OB(O)c1ccc(B2OCCO2)cc1']:
        row=decorated_monomers([dict(id='edit',label='Edited linker',smiles=smiles)])[0]
        assert not row['assembly_roles']['self_kind']
