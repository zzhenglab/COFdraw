"""A periodic design must match chemistry and unchanged incidence, not labels."""
from copy import deepcopy
from pathlib import Path
import json
import shutil

import pytest
from fastapi.testclient import TestClient
from cof_vision.cofdraw_chemistry import prepare_template, design_signature
from cof_vision.cofdraw_core import MONOMERS
from cof_vision.cofdraw_periodic_options import periodic_options


def recipe():
    r=prepare_template(MONOMERS[0]['smiles'],MONOMERS[2]['smiles'],'imine','untrusted A','untrusted B',12,'hcb_3_2')
    r['design_template_signature']=design_signature(r)
    return r


def manifest():
    return {'templates':[{'id':'exact','label':'Design','atoms':66,'formula':'C36H24N6','sha256':'a'*64,'description':'Unrelaxed',
      'precursor_match':{'canonical_smiles':[MONOMERS[0]['smiles'],MONOMERS[2]['smiles']],'linkage':'imine',
                         'drawing_template_ids':['hcb_3_2'],'preserves_drawing_coordinates':False}}]}


def test_match_canonical_graphs_and_not_labels_or_drawing_coordinates():
    r=recipe();r['components']['A']['label']='Anything';r['instances'][0]['layout_position']=[100,33]
    r['components']['A']['monomer_smiles']='c1(C=O)cc(C=O)cc(C=O)c1'
    r['design_template_signature']=design_signature(r)
    result=periodic_options(r,manifest())
    assert [x['id'] for x in result['options']]==['exact']
    assert result['preserves_drawing_coordinates'] is False


@pytest.mark.parametrize('change',['precursor','linkage','connectivity','fragment','signature','template','explicit_graph','forged_template'])
def test_modified_or_unknown_chemistry_never_acquires_matching_cell(change):
    r=recipe();m=manifest()
    if change=='precursor':m['templates'][0]['precursor_match']['canonical_smiles'][0]=MONOMERS[1]['smiles']
    if change=='linkage':m['templates'][0]['precursor_match']['linkage']='beta_ketoenamine'
    if change=='connectivity':r['connections'].pop()
    if change=='fragment':
        r['components']['A']['fragment_smiles']='CC';r['design_template_signature']=design_signature(r)
    if change=='signature':r.pop('design_template_signature')
    if change=='template':m['templates'][0]['precursor_match']['drawing_template_ids']=[]
    if change=='explicit_graph':r['explicit_graph']={'atoms':[{'element':'C','x':0,'y':0}],'bonds':[]}
    if change=='forged_template':
        r['selected_template']['id']='custom_pore';m['templates'][0]['precursor_match']['drawing_template_ids']=['custom_pore']
    result=periodic_options(r,m)
    assert result['options']==[] and result['reason']


def test_actual_api_match_is_revalidated_before_model_generation(tmp_path):
    from cof_vision.cofdraw_desktop import create_desktop_app
    shutil.copytree(Path(__file__).resolve().parents[1]/'data/crystal-models',tmp_path/'crystal-models')
    with TestClient(create_desktop_app(tmp_path)) as client:
        prepared=client.post('/api/drawing/prepare',json={'monomer_a':MONOMERS[0]['smiles'],'monomer_b':MONOMERS[2]['smiles'],'template_id':'hcb_3_2'}).json()
        request={'recipe':prepared['recipe']}
        response=client.post('/api/drawing/science/periodic-options',json=request)
        assert response.status_code==200,response.text
        assert 'hcb-imine-aa' in [r['id'] for r in response.json()['options']]
        wrong=client.post('/api/drawing/science/model',json={'template_id':'dia-imine','recipe_match':request})
        assert wrong.status_code==422
        valid=client.post('/api/drawing/science/model',json={'template_id':'hcb-imine-aa','recipe_match':request})
        assert valid.status_code==200,valid.text
        result=valid.json();match=result['structure']['provenance']['drawing_match']
        assert match['preserves_drawing_coordinates'] is False
        assert client.get(result['downloads']['cif']).status_code==200
