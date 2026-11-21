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
    assert repeat_options(recipe)['extents'][1]['atoms']==mol.GetNumAtoms()


def test_actual_kagome_pore_tints_and_v3000_exports(api):
    request = {'monomer_ids':['tape','ta','bpda_cho'],'linker_sequence':'BCBCBC','pore_layout':'kgm'}
    prepared=api.post('/api/drawing/prepare/assembly',json=request)
    assert prepared.status_code==200,prepared.text
    body={'recipe':prepared.json()['recipe'],'style_id':'standard','size':500}
    options=api.post('/api/drawing/repeat-options',json=body).json()
    assert [o['label'] for o in options['extents'][1:]]==['7 pores','20 pores']
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
