"""Generate cage and saturated-spacer designs from editable precursor SMILES."""
import pytest
from rdkit import Chem

from cof_vision.cofdraw_chemistry import design_signature, inspect_monomer, prepare_template
from cof_vision.cofdraw_cofex import parse_record, recipe_record
from cof_vision.cofdraw_repeat import repeat_patch
from cof_vision.framework_drawing import ALKYL_VARIANTS, VERSION, TP, build_patch, linkage_fragment


TAA = 'NC12CC3(N)CC(N)(C1)CC(N)(C2)C3'
TFB = 'O=Cc1cc(C=O)cc(C=O)c1'
BCO = 'O=CC12CCC(C=O)(CC1)CC2'
CHDA = 'N[C@H]1CC[C@H](N)CC1'


@pytest.mark.parametrize('aldehyde,donor,link,net,template,product', [
    (TFB, TAA, 'imine', 'ctn', 'ctn_4_3', '[c][CH]=[N]-[C;X4]'),
    (TP, TAA, 'beta_ketoenamine', 'ctn', 'ctn_4_3', '[O]=[C]~[C](=[CH]-[NH]-[C;X4])'),
    (BCO, TAA, 'imine', 'dia', 'dia_4_2', '[C;X4]-[CH]=[N]-[C;X4]'),
    (BCO, TAA, 'secondary_amine', 'dia', 'dia_4_2', '[C;X4]-[CH2]-[NH]-[C;X4]'),
    (TFB, CHDA, 'imine', 'hcb', 'hcb_3_2', '[c][CH]=[N]-[C;X4]'),
])
def test_product_and_editable_net_from_monomers(aldehyde, donor, link, net, template, product):
    recipe=prepare_template(aldehyde,donor,link,'A','B',6,template)
    recipe.update(recipe_id='alkyl-design',schema_version=VERSION,design_template_signature=design_signature(recipe))
    record=recipe_record(recipe)
    assert record['compact'].endswith(' COFex.'+net+'.cat0')
    assert record['linkage']==ALKYL_VARIANTS[link]
    restored=parse_record(record['text'],[])[0]['prepare_request']
    assert restored['template_id']==template and restored['linkage']==record['linkage']
    initial,_=build_patch(recipe)
    extended,_,_,summary,_=repeat_patch(recipe,2)
    pattern=Chem.MolFromSmarts(product)
    assert initial.HasSubstructMatch(pattern) and extended.GetNumAtoms()>initial.GetNumAtoms()
    assert len(Chem.GetMolFrags(extended))==1 and extended.HasSubstructMatch(pattern)
    if donor==CHDA:
        stereo=Chem.rdchem.ChiralType.CHI_UNSPECIFIED
        assert sum(atom.GetChiralTag()!=stereo for atom in initial.GetAtoms())==12
    if net=='dia':
        assert all('crystal_project' not in choice for choice in parse_record(record['compact'],[]))


def test_amine_geometry_and_chemical_role_constraints():
    for donor in ('NCCN','NC1CCC(N)CC1','N[C@H]1CC[C@@H](N)CC1'):
        assert inspect_monomer(donor)['site_shapes']['primary_alkylamine']=='unresolved_alkyl'
        with pytest.raises(ValueError,match='geometry'):
            prepare_template(TFB,donor,'imine','A','B',12,'hcb_3_2')
        assert prepare_template(TFB,donor,'imine','A','B',6,'local_hexagon')['linkage']=='alkyl_imine'
    for donor in ('NC(=O)CC(=O)N','NNC(=O)CC(=O)NN','CNCCNC','[NH3+]CC[NH3+]','Nc1ccc(N)cc1','[2H]N([2H])CCN([2H])[2H]'):
        with pytest.raises(ValueError):
            linkage_fragment(donor,'alkyl_imine','primary_alkylamine')
    for acceptor in ('CC(=O)C','NC=O','C=O'):
        with pytest.raises(ValueError):
            linkage_fragment(acceptor,'alkyl_imine','carbon_aldehyde')
