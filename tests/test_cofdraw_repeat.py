"""A larger drawing must export the same shared-edge graph shown on screen."""
from collections import Counter
from pathlib import Path
from xml.etree import ElementTree

import pytest
from fastapi.testclient import TestClient
from rdkit import Chem

from cof_vision.cofdraw_chemistry import design_signature, prepare_template
from cof_vision.cofdraw_desktop import create_desktop_app
from cof_vision.cofdraw_repeat import repeat_options, repeat_patch
from cof_vision.framework_drawing import VERSION


@pytest.fixture
def api(tmp_path):
    source = Path(__file__).resolve().parents[1]/'data'
    with TestClient(create_desktop_app(source,output_root=tmp_path)) as client:
        yield client


def test_binary_shared_faces_keep_isotopes_and_stereochemistry():
    recipe = prepare_template('O=Cc1cc([13CH]=O)cc(C=O)c1',
                              'Nc1cc([C@H](F)Cl)c(N)cc1/C=C/F','imine','node A','linker B',12,'hcb_3_2')
    recipe.update(schema_version=VERSION,recipe_id='test')
    recipe['design_template_signature']=design_signature(recipe)
    mol,origins,ports,summary,_ = repeat_patch(recipe,2)
    assert summary['pore_count']==4 and len(summary['pore_faces'])==4
    assert summary['blocks']==35 and summary['interblock_bonds']==38 and len(ports)==10
    assert len(Chem.GetMolFrags(mol))==1
    assert summary['component_counts']=={'A':16,'B':19} and len(origins)==mol.GetNumAtoms()
    assert sum(a.GetIsotope()==13 for a in mol.GetAtoms())==16
    assert sum(a.GetChiralTag()!=Chem.ChiralType.CHI_UNSPECIFIED for a in mol.GetAtoms())==19
    assert sum(b.GetStereo() in (Chem.BondStereo.STEREOE,Chem.BondStereo.STEREOTRANS)
               for b in mol.GetBonds())==19
    extents={o['value']:o for o in repeat_options(recipe)['extents']}
    assert extents[2]['atoms']==mol.GetNumAtoms()
    tri,origins,ports,summary,_=repeat_patch(recipe,3,extent_mode='triad')
    faces=[set(face) for face in summary['pore_faces']]
    # Three hexagons meet at a single junction. Every pair shares one
    # linker and the two adjacent junctions, each present once in the MOL.
    assert len(faces)==3 and len(set.intersection(*faces))==1
    assert all(len(faces[i]&faces[i-1])==3 for i in range(3))
    assert summary['component_counts']=={'A':13,'B':15}
    assert summary['interblock_bonds']==30 and len(ports)==9
    assert tri.GetNumAtoms()==extents['three-pores']['atoms']
    assert len({(o['instance'],o['fragment_atom_index']) for o in origins})==tri.GetNumAtoms()


def test_actual_kagome_pore_tints_and_v3000_exports(api):
    request = {'monomer_ids':['tape','ta','bpda_cho'],'linker_sequence':'BCBCBC','pore_layout':'kgm'}
    prepared=api.post('/api/drawing/prepare/assembly',json=request)
    assert prepared.status_code==200,prepared.text
    body={'recipe':prepared.json()['recipe'],'style_id':'standard','size':500}
    options=api.post('/api/drawing/repeat-options',json=body).json()
    assert [o['label'] for o in options['extents'][1:]]==['3 pores','7 pores','20 pores']
    result=api.post('/api/drawing/repeat',json={**body,'repeats':3,'color_pores':True})
    assert result.status_code==200,result.text
    report=result.json()
    summary,urls=report['summary'],report['downloads']
    assert summary['pore_sequences']=={'BBB':8,'CCC':8,'BCBCBC':4}
    colors=summary['pore_colors']
    svg=ElementTree.fromstring(api.get(urls['svg']).content)
    counts=Counter(p.get('fill') for p in svg.iter() if p.get('data-pore-fill'))
    assert counts=={colors[label]:n for label,n in summary['pore_sequences'].items()}
    mol=api.get(urls['mol']).content
    assert b'V3000' in mol and b'M  V30 COUNTS 1418 1625' in mol
    assert api.get(urls['png']).content.startswith(b'\x89PNG')
    assert api.get(urls['sdf']).content.startswith(mol.rstrip())
    audit=api.get(urls['json']).json()['audit']
    assert audit['atoms']==1418 and len(audit['atom_origin'])==1418
    if 'cdxml' in urls:
        assert audit['chemical_exports']['molecule_input']=='verified_source_mol_v3000'
        assert len(list(ElementTree.fromstring(api.get(urls['cdxml']).content).iter('n')))==1418
    plain=api.post('/api/drawing/repeat',json={**body,'repeats':2}).json()
    assert plain['summary']['pore_colors']=={} and 'data-pore-fill' not in api.get(plain['svg']).text


def test_expansion_uses_only_reviewed_source_boundaries(api):
    source='wiley-r198-tppa-so3h-pore-0002'
    checked=api.get('/api/drawing/sources/'+source).json()
    selected={'source_recipe_id':source,'size':500}
    response=api.post('/api/drawing/repeat',json={**selected,'repeats':2}).json()
    assert response['summary']['source_drawing_basis']['recipe_sha256']==checked['source_sha256']
    assert response['summary']['pore_count']==4
    other={'source_recipe_id':'pyurea-source-layout-v1','size':500}
    options=api.post('/api/drawing/repeat-options',json=other).json()
    assert options['extents']==[{'value':1,'label':'Source drawing'}] and not options['available']
    assert api.post('/api/drawing/repeat',json={**other,'repeats':2}).status_code==422
    triad=api.post('/api/drawing/repeat',json={**selected,'repeats':3,'extent_mode':'triad'}).json()
    assert triad['summary']['pore_count']==3
    assert triad['summary']['source_drawing_basis']['recipe_sha256']==checked['source_sha256']


@pytest.mark.parametrize(('inputs','expected'),[
    ({'monomer_ids':['bdba']},{'hexagon':3}),
    ({'monomer_ids':['bfbaeb'],'pore_view':'triangle'},{'triangle':2,'hexagon':1}),
    ({'monomer_ids':['tape','ta','bpda_cho'],'linker_sequence':'CCC','pore_layout':'kgm'},
     {'BBB':1,'CCC':1,'BCBCBC':1}),
    ({'monomer_ids':['tapp','ta','bpda_cho'],'linker_sequence':'BCBC'},{'BCBC':3}),
])
def test_three_pore_choice_returns_the_visible_molecule(api,inputs,expected):
    recipe=api.post('/api/drawing/prepare/assembly',json=inputs).json()['recipe']
    body={'recipe':recipe,'size':500}
    option=next(o for o in api.post('/api/drawing/repeat-options',json=body).json()['extents']
                if o['value']=='three-pores')
    response=api.post('/api/drawing/repeat',json={**body,'repeats':option['repeats'],
                      'extent_mode':option['extent_mode'],'color_pores':True})
    assert response.status_code==200,response.text
    result=response.json();summary=result['summary'];urls=result['downloads']
    assert summary.get('pore_counts',summary.get('pore_sequences'))==expected
    assert summary['pore_count']==summary['interblock_bonds']-summary['blocks']+1==3
    assert summary['atoms']==option['atoms']
    svg=ElementTree.fromstring(api.get(urls['svg']).content)
    assert sum(p.get('data-pore-fill') is not None for p in svg.iter())==3
    mol=Chem.MolFromMolBlock(api.get(urls['mol']).text,sanitize=False,removeHs=False)
    assert len(Chem.GetMolFrags(mol))==1 and mol.GetNumAtoms()==summary['atoms']
