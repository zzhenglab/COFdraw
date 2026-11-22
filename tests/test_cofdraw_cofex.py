"""Round-trip the identities exposed under the active 2D and periodic models."""
from pathlib import Path
import json

import pytest
from fastapi.testclient import TestClient

from cof_vision.cofdraw_desktop import create_desktop_app


DATA = Path(__file__).resolve().parents[1] / 'data'


@pytest.fixture()
def client(tmp_path):
    with TestClient(create_desktop_app(DATA,output_root=tmp_path)) as active:
        yield active


def post(client,path,payload):
    result = client.post('/api/drawing/'+path,json=payload)
    assert result.status_code == 200, result.text
    return result.json()


@pytest.mark.parametrize(('ids','settings'),[
    (['bdba'],{}),
    (['bfbaeb'],{'pore_view':'triangle'}),
    (['tapb','ta','bpda_cho'],{'linker_sequence':'BCBBCC'}),
])
def test_assembly_text_restores_feeds_and_face(client,ids,settings):
    original = post(client,'prepare/assembly',dict(monomer_ids=ids,**settings))['recipe']
    record = post(client,'render',dict(recipe=original,size=500))['audit']['drawing_cofex']
    restored, = post(client,'cofex/resolve',dict(text=record['text']))['choices']
    assembly = post(client,'prepare/assembly',restored['assembly_request'])['recipe']['assembly']
    assert assembly['monomer_ids']==ids
    assert all(assembly[k]==v for k,v in settings.items())


def test_bare_notation_offers_compatible_products_and_preserves_isotopes(client):
    # COFex orders SMILES lexically. The preparer expects aldehyde first.
    text = 'Nc1ccc(N)cc1.O=[13CH]c1cc(C=O)cc(C=O)c1 COFex.hcb.cat0;example'
    choices = post(client,'cofex/resolve',dict(text=text))['choices']
    assert {p['prepare_request']['linkage'] for p in choices}=={'imine','secondary_amine'}
    reduced = next(p['prepare_request'] for p in choices if p['prepare_request']['linkage']=='secondary_amine')
    assert '[13CH]' in reduced['monomer_a']
    assert reduced['template_id']=='hcb_3_2'
    previous_record = json.dumps({'design_cofex_notation':text, 'linkage':'secondary_amine',
                                  'selected_template':{'id':'hcb_3_2'}})
    only, = post(client,'cofex/resolve',dict(text=previous_record))['choices']
    assert only['prepare_request']==reduced


def test_registered_id_and_notation_resolve_to_same_hash_bound_drawing(client):
    source = client.get('/api/drawing/sources/tppa-pore-v1').json()
    meta = source['metadata']
    for text in (meta['cofex_notation'],meta['cofex_id']):
        choices=post(client,'cofex/resolve',dict(text=text))['choices']
        assert any(p.get('source_recipe_id')=='tppa-pore-v1' for p in choices)


def test_diamond_phase_counts_and_reference_paste(client):
    for id_,count in [('cof-300-c5',5),('curated-13040N3',9)]:
        phase = post(client,'science/model',dict(template_id=id_,repeats=[1,1,1]))['structure']['provenance']['drawing_cofex']
        assert f' COFex.dia.cat{count};' in phase['text'] and phase['scope']=='source'
        only, = post(client,'cofex/resolve',dict(text=phase['text']))['choices']
        assert only['crystal_project']=={'template_id':id_,'version':2,'reference_only':True}
    pair = dict(monomer_a='Nc1ccc(C(c2ccc(N)cc2)(c2ccc(N)cc2)c2ccc(N)cc2)cc1', monomer_b='O=Cc1ccc(C=O)cc1',
                label_a='TAPM',label_b='TA',linkage='imine',network_count=2)
    product = post(client,'science/build-dia',pair)['structure']['provenance']['drawing_cofex']
    assert ' COFex.dia.cat2;' in product['text'] and product['scope']=='design'
    only, = post(client,'cofex/resolve',dict(text=product['text']))['choices']
    settings=only['crystal_project']
    assert settings['network_count']==2 and settings['topology']=='dia' and settings['linkage']=='imine'
    assert [m['label'] for m in settings['monomers']]==['TAPM','TA']


@pytest.mark.parametrize('model,code',[('curated-07011N3','ctn.cat0'),('curated-16250N3','pts.cat2'),
                                       ('curated-21270N3','pts.cat5'),('curated-21271N3','ljh.cat0')])
def test_reference_junction_and_copy_full_record(client,model,code):
    report=post(client,'science/model',dict(template_id=model,repeats=[1,1,1]))
    provenance=report['structure']['provenance']
    assert 'COFex.'+code+';' in provenance['drawing_cofex']['text']
    matched=post(client,'cofex/resolve',dict(text=json.dumps(provenance)))['choices']
    assert any(p.get('crystal_project',{}).get('template_id')==model for p in matched)
    endpoint=f"science/reports/{report['report_id']}/circuits/1"
    base=post(client,endpoint,{})
    edited=post(client,endpoint,{'arm_spread':16})
    assert edited['junction_svg']!=base['junction_svg']==report['chemical_junction_svg']
    assert edited['projection_svg']==base['projection_svg']


def test_reference_precursor_colors_follow_record_order_and_group_composition(client):
    from xml.etree import ElementTree as ET
    from cof_vision.cofdraw_crystal_figures import _reference_components
    report=post(client,'science/model',{'template_id':'curated-16250N3'})
    structure=report['structure']
    record=structure['provenance']['drawing_cofex']
    endpoint=f"science/reports/{report['report_id']}/circuits/2"
    def ownership(svg):
        return {el.get('data-block') for el in ET.fromstring(svg).iter() if el.get('data-block')}
    base=post(client,endpoint,{})
    assert ownership(base['scheme_svg'])==ownership(base['projection_svg'])=={'A','B'}
    tilted=post(client,endpoint,{'tilt':15})
    assert tilted['scheme_svg']!=base['scheme_svg'] and ownership(tilted['scheme_svg'])=={'A','B'}
    graph=[[] for atom in structure['atoms']]
    edges=report['display_edges']
    for e in edges:
        a,b,t=e['a'],e['b'],e['translation']
        graph[a].append((b,t));graph[b].append((a,tuple(-v for v in t)))
    elements=[a['element'] for a in structure['atoms']]
    original=_reference_components(record,graph,elements,edges)
    reversed_record={**record,'precursor_smiles':record['precursor_smiles'][::-1]}
    reverse=_reference_components(reversed_record,graph,elements,edges)
    assert set(original)=={'A','B'} and all(x!=y for x,y in zip(original,reverse))


@pytest.mark.parametrize('text',[
    '  ', 'C..N COFex.dia.cat2', 'C.C.C.C COFex.hcb.cat0', 'C.[*:1] COFex.sql.cat0',
    'C.N COFex.lon.cat0;COFdraw/1?mode=crystal&order=0%2C0&link=imine',
    'C.N COFex.dia.cat2;COFdraw/1?mode=crystal&order=0%2C1&link=imine&face=triangle',
    '{"source_recipe_id":"tppa-pore-v1","cofex_notation":42}', 'C COFex.unknown.cat9',
])
def test_bad_or_unsupported_cofex_never_switches_to_unrelated_structure(client,text):
    result=client.post('/api/drawing/cofex/resolve',json={'text':text})
    assert result.status_code==422
