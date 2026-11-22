import pytest
from rdkit import Chem

from cof_vision.cofdraw_chemistry import design_signature, inspect_monomer, prepare_template
from cof_vision.cofdraw_cofex import parse_record, recipe_record
from cof_vision.cofdraw_repeat import repeat_patch
from cof_vision.framework_drawing import VERSION, build_patch, linkage_fragment


def design(a,b,link,template):
    recipe=prepare_template(a,b,link,'Node','Spacer',6,template)
    recipe.update(schema_version=VERSION,recipe_id='precursor-test',design_template_signature=design_signature(recipe))
    restored=parse_record(recipe_record(recipe)['text'],[])
    assert restored[0]['prepare_request']['template_id']==template
    return recipe


def test_extended_tetraaldehyde_and_aryl_diacetonitrile_generate_kagome_from_smiles():
    node='O=Cc1ccc(-c2ccc(C(=C(c3ccc(-c4ccc(C=O)cc4)cc3)c3ccc(-c4ccc(C=O)cc4)cc3)c3ccc(-c4ccc(C=O)cc4)cc3)cc2)cc1'
    spacer='N#CCc1ccc(CC#N)cc1'
    recipe=design(node,spacer,'cyano_vinylene_sp2_carbon','kgm_4_2')
    mol,_,ends,summary,_=repeat_patch(recipe,2)
    assert summary['pore_counts']=={'triangle':6,'hexagon':1}
    assert summary['component_counts']=={'A':12,'B':18}
    assert len(ends)==12
    assert len(mol.GetSubstructMatches(Chem.MolFromSmarts('[c][CH1]=[C]([C]#[N])[c]')))==36
    assert not inspect_monomer('N#Cc1ccc(C#N)cc1')['kagome_roles']
    assert not inspect_monomer('N#CCc1cccc(CC#N)c1')['kagome_roles']


@pytest.mark.parametrize('a,b,link,expected',[
    ('O=Cc1ccc(C=O)cc1','Cc1nc(C)nc(C)n1','vinylene_methyltriazine','c[CH]=[CH]c'),
    ('O=C(Cl)c1cc(C(=O)Cl)cc(C(=O)Cl)c1','Nc1ccc(N)cc1','amide','cC(=O)[NH]c'),
    ('O=C(O)c1cc(C(=O)O)cc(C(=O)O)c1','Nc1ccc(N)cc1','amide','cC(=O)[NH]c'),
    ('O=c1oc(=O)c2cc3c(=O)oc(=O)c3cc12','Nc1ccc(N(c2ccc(N)cc2)c2ccc(N)cc2)cc1','imide','[O]=[C,c]~[n,N](-c)~[c,C]=O'),
    ('OB(O)c1ccc(B(O)O)cc1','Oc1cc2c3cc(O)c(O)cc3c3cc(O)c(O)cc3c2cc1O','boronate_ester','c[B]1Oc2ccccc2O1'),
])
def test_local_product_graph_and_import(a,b,link,expected):
    recipe=design(a,b,link,'hcb_3_2')
    mol,_=build_patch(recipe)
    assert len(Chem.GetMolFrags(mol))==1
    assert len(mol.GetSubstructMatches(Chem.MolFromSmarts(expected)))==12


def test_regioselectivity_and_isotopic_transfer_are_not_silently_assigned():
    for s,link,role in [('Oc1cccc(O)c1O','boronate_ester','aryl_vicinal_diol'),
                        ('O[10B](O)c1ccccc1','boronate_ester','aryl_boronic_acid'),
                        ('[15NH2]c1ccc(N)cc1','imide','primary_arylamine'),
                        ('Cc1ccc(C)cc1','vinylene_methyltriazine','methyltriazine')]:
        with pytest.raises(ValueError):
            linkage_fragment(s,link,role)


def test_user_selected_local_ring_extent_does_not_assign_a_net():
    node='O=Cc1cc(C=O)cc(C=O)c1'
    amine='Nc1ccc(-c2cc(-c3ccc(N)cc3)cc(-c3ccc(N)cc3)c2)cc1'
    recipe=design(node,amine,'imine','local_hexagon')
    _,_,_,summary,_=repeat_patch(recipe,2)
    assert summary['pore_count']==2 and summary['blocks']==10
    assert summary['interblock_bonds']==11
    assert recipe_record(recipe)['compact'].endswith(' COFex.unk.cat0')

def test_source_comparison_preserves_substitution_and_distinguishes_a_boundary_cut():
    from scripts.audit_annotation_drawings import aryl_boundary, same_graph
    rim=Chem.MolFromSmiles('*=Cc1ccc(C=C(C#N)c2ccc(N*)cc2)cc1')
    clipped=aryl_boundary(rim)
    target=Chem.MolFromSmiles('*c1ccc(C=C(C#N)c2ccc(*)cc2)cc1')
    assert clipped.GetNumAtoms()==rim.GetNumAtoms()-2
    assert same_graph(clipped,Chem.RenumberAtoms(target,list(reversed(range(target.GetNumAtoms())))))
    assert not same_graph(clipped,Chem.MolFromSmiles('*c1ccc(C=C(C#N)c2cccc(*)c2)cc1'))
    assert not same_graph(target,Chem.MolFromSmiles('*c1ccc(C=C([13C]#N)c2ccc(*)cc2)cc1'))
