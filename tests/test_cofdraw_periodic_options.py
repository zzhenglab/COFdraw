"""A periodic design must match chemistry and unchanged incidence, not labels."""
from collections import Counter
import hashlib
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
        automatic=client.post('/api/drawing/science/build-periodic',json=request).json()
        assert automatic['available'] is True
        assert automatic['cell_validation']['template_cif_sha256']==hashlib.sha256(
            (tmp_path/'crystal-models/hcb-imine-aa.cif').read_bytes()).hexdigest()
        assert automatic['structure']['atoms']==result['structure']['atoms']


def test_cof300_recorded_projection_opens_published_interpenetrated_crystal(tmp_path):
    from cof_vision.cofdraw_desktop import create_desktop_app
    data=Path(__file__).resolve().parents[1]/'data'
    for folder in ('framework-drawing','crystal-models'):
        shutil.copytree(data/folder,tmp_path/folder,ignore=shutil.ignore_patterns('exports'))
    with TestClient(create_desktop_app(tmp_path)) as client:
        source='cof-300-projection-v1'
        result=client.post('/api/drawing/science/build-periodic',json={
            'source_recipe_id':source,'repeats':[2,1,2]}).json()
        s=result['structure']
        assert result['available'] and len(s['atoms'])==292 and s['provenance']['network_dimension']==3
        atoms=result['viewer']['atoms']
        assert len({a['network'] for a in atoms})==5
        assert all(atoms[i]['network']==atoms[j]['network'] for i,j in result['viewer']['bonds'])
        recipe=client.get('/api/drawing/sources/'+source).json()['recipe']
        assert not client.post('/api/drawing/science/build-periodic',json={'recipe':recipe}).json()['available']


def test_registered_pores_build_fresh_audited_cells_and_keep_sulfonic_acid(tmp_path):
    from cof_vision.cofdraw_desktop import create_desktop_app
    from cof_vision.cofdraw_network import inspect_periodic_network
    from cof_vision.cofdraw_science import import_cif
    from cof_vision.cofdraw_periodic_build import edge_key
    data=Path(__file__).resolve().parents[1]/'data'
    for folder in ('framework-drawing','crystal-models'):
        shutil.copytree(data/folder,tmp_path/folder,ignore=shutil.ignore_patterns('exports'))
    expected={
        'tfb-tab-pore-v1':dict(C=33,H=21,N=3),
        'tppa-pore-v1':dict(C=36,H=24,N=6,O=6),
        'tpb-dmtp-cof-pore-v1':dict(C=78,H=60,N=6,O=6),
        'wiley-r198-tppa-so3h-pore-0002':dict(C=36,H=24,N=6,O=15,S=3),
    }
    with TestClient(create_desktop_app(tmp_path)) as client:
        for source,counts in expected.items():
            request={'source_recipe_id':source}
            entry=client.get('/api/drawing/sources/'+source).json()
            choices=client.post('/api/drawing/science/periodic-options',json=request).json()
            assert choices['source_drawing_basis']['recipe_sha256']==entry['source_sha256']
            result=client.post('/api/drawing/science/build-periodic',json=request).json()
            assert result.get('available'),result
            tags=result['structure']['provenance']
            assert tags['kind']=='cofdraw_periodic_starting_model'
            assert tags['source_drawing_basis']['recipe_id']==source
            assert tags['linkage']==entry['metadata']['periodic_design_inputs']['linkage']
            parsed=import_cif(client.get(result['downloads']['cif']).text)
            assert Counter(a['element'] for a in parsed['atoms'])==counts
            net=inspect_periodic_network(parsed)
            assert net['dimensions']==[2] and not net['short_contact_count']
            def identities(edges):
                return {edge_key(e['a'],e['b'],e['translation']) for e in edges}
            assert identities(net['edges'])==identities(result['cell_validation']['periodic_bonds'])
        assert tags['linkage']=='beta_ketoenamine'
        assert tags['spacing_angstrom']>=tags['interlayer_contact_height_angstrom']>3.4
        assert '3.40' in result['scope'] and tags['source_drawing_basis']['comparison_reference']['doi']=='10.1021/acs.chemmater.5b04947'
        rejected=client.post('/api/drawing/science/build-periodic',json={**request,'layer_spacing':3.4}).json()
        assert not rejected['available'] and 'contact screen' in rejected['reason']
        # A requested valid height should survive export. A replacement diamine
        # creates its own desulfonated design with no inherited paper identity.
        custom=client.post('/api/drawing/science/build-periodic',json={**request,'layer_spacing':6.1}).json()
        assert custom['structure']['provenance']['spacing_angstrom']==6.1
        drawing=client.post('/api/drawing/prepare',json=dict(
            monomer_a=tags['precursors'][0],monomer_b='Nc1ccc(N)cc1',linkage=tags['linkage'],template_id=tags['template'])).json()['recipe']
        edited=client.post('/api/drawing/science/build-periodic',json={'recipe':drawing}).json()['structure']
        assert 'source_drawing_basis' not in edited['provenance'] and not any(a['element']=='S' for a in edited['atoms'])


def test_exact_reported_ribbon_crystal_matches_the_pair_and_linkage(tmp_path):
    from cof_vision.cofdraw_desktop import create_desktop_app
    from cof_vision.cofdraw_network import inspect_periodic_network
    data=Path(__file__).resolve().parents[1]/'data/crystal-models'
    shutil.copytree(data,tmp_path/'crystal-models')
    model=json.loads((data/'manifest.json').read_text())['reported_models'][0]
    a,b=model['precursor_match']['canonical_smiles']
    with TestClient(create_desktop_app(tmp_path)) as client:
        for linkage in ('imine','secondary_amine'):
            payload=dict(monomer_a=a,monomer_b=b,linkage=linkage,template_id='alternating_hexagon')
            recipe=client.post('/api/drawing/prepare',json=payload).json()['recipe']
            result=client.post('/api/drawing/science/build-periodic',json={'recipe':recipe}).json()
            if linkage=='secondary_amine':
                assert result['available'] is False
                continue
            s=result['structure']
            assert s['provenance']['kind']=='reported_crystal'
            assert Counter(a['element'] for a in s['atoms'])=={'C':320,'H':216,'N':32}
            network=inspect_periodic_network(s)
            assert network['dimensions']==[1] and not network['short_contacts']
            assert model['citation_doi'] in client.get(result['downloads']['cif']).text


@pytest.mark.parametrize('a,b,linkage,template,formula,joins', [
    ('tp','pa','beta_ketoenamine','hcb_3_2',{'C':36,'N':6,'O':6,'H':24},6),
    ('tp','pa_so3h','beta_ketoenamine','hcb_3_2',{'C':36,'N':6,'O':15,'H':24,'S':3},6),
    ('ta','tapp','imine','sql_4_2',{'C':60,'N':8,'H':38},4),
])
def test_builds_an_unbundled_current_cell_with_real_bonds_and_cif(tmp_path,a,b,linkage,template,formula,joins):
    from cof_vision.cofdraw_desktop import create_desktop_app
    from cof_vision.cofdraw_network import inspect_periodic_network
    from cof_vision.cofdraw_science import import_cif
    catalog=Path(__file__).resolve().parents[1]/'data/framework-drawing/catalog/current.json'
    rows={row['id']:row['smiles'] for row in json.loads(catalog.read_text('utf8'))['monomer_presets']}
    with TestClient(create_desktop_app(tmp_path)) as client:
        prepared=client.post('/api/drawing/prepare',json=dict(monomer_a=rows[a],monomer_b=rows[b],linkage=linkage,template_id=template))
        assert prepared.status_code==200,prepared.text
        drawing=prepared.json()['recipe']
        response=client.post('/api/drawing/science/build-periodic',json={'recipe':drawing})
        assert response.status_code==200,response.text
        result=response.json()
        assert result['available'],result
        cif=client.get(result['downloads']['cif']).text
        assert 'unrelaxed AA starting model' in cif
        cell=import_cif(cif)
        assert Counter(atom['element'] for atom in cell['atoms'])==formula
        assert result['structure']['provenance']['kind']=='cofdraw_periodic_starting_model'
        assert result['structure']['provenance']['linkage']==linkage
        network=inspect_periodic_network(cell)
        assert network['dimensions']==[2] and network['quotient_component_count']==1
        assert not network['short_contacts'] and not network['isolated_atom_indices']
        assert network['edge_count']==result['cell_validation']['expected_atom_edges_verified']
        assert len(drawing['connections'])>joins
        if template=='sql_4_2':
            assert result['structure']['provenance']['tilted_phenyl_ring_count']==4
            shifted=client.post('/api/drawing/science/stacking',json={
                'structure':result['structure'],'stacking':'AB','spacing':5.,'corrugated_layer':True})
            assert shifted.status_code==200,shifted.text
            assert len(shifted.json()['structure']['atoms'])==2*sum(formula.values())
            pxrd=client.post('/api/drawing/science/pxrd',json={
                'structure':shifted.json()['structure'],'two_theta_max':20.})
            assert pxrd.status_code==200,pxrd.text
            assert pxrd.json()['peaks'] and client.get(pxrd.json()['downloads']['csv']).status_code==200
        drawing['connections'].pop()
        drawing['design_template_signature']=design_signature(drawing)
        rejected=client.post('/api/drawing/science/build-periodic',json={'recipe':drawing}).json()
        assert rejected['available'] is False and 'downloads' not in rejected


@pytest.mark.parametrize('a,b,linkage,template,method', [
    ('ta', 'pytta', 'imine', 'sql_4_2', 'phenyl_torsions'),
    ('tfppy_13c', 'benzidine', 'imine', 'sql_4_2', 'phenyl_torsions'),
    ('tfppy_13c', 'recovered_bd', 'imine', 'sql_4_2', 'restrained_mmff94'),
    ('tp', 'deth', 'keto_hydrazide', 'hcb_3_2', 'restrained_mmff94'),
])
def test_crowded_cells_retain_each_translated_bond_after_cif_export(tmp_path, a, b, linkage, template, method):
    from cof_vision.cofdraw_desktop import create_desktop_app
    from cof_vision.cofdraw_network import inspect_periodic_network
    from cof_vision.cofdraw_periodic_build import edge_key
    from cof_vision.cofdraw_science import import_cif
    path = Path(__file__).resolve().parents[1] / 'data/framework-drawing/catalog/current.json'
    rows = {r['id']: r['smiles'] for r in json.loads(path.read_text('utf8'))['monomer_presets']}
    with TestClient(create_desktop_app(tmp_path)) as client:
        prepared = client.post('/api/drawing/prepare', json={'monomer_a': rows[a], 'monomer_b': rows[b],
                               'linkage': linkage, 'template_id': template})
        assert prepared.status_code == 200, prepared.text
        response = client.post('/api/drawing/science/build-periodic', json={'recipe': prepared.json()['recipe']})
        assert response.status_code == 200, response.text
        result = response.json()
        assert result['available'] and result['structure']['provenance']['fragment_geometry'] == method
        cif = client.get(result['downloads']['cif']).text
        parsed = import_cif(cif)
        if a == 'tfppy_13c':
            from cof_vision.cofdraw_characterization import structure_summary
            from cof_vision.cofdraw_science import structure_exports, patch_3d
            labels = [site for site in parsed['atoms'] if site.get('isotope') == 13]
            assert len(labels) == 4 and all(site['element'] == 'C' for site in labels)
            natural = deepcopy(parsed)
            for site in natural['atoms']:
                site.pop('isotope', None)
            mass_gain = structure_summary(parsed)['cell_molar_mass_g_mol'] - structure_summary(natural)['cell_molar_mass_g_mol']
            assert mass_gain == pytest.approx(4 * (13.00335484 - 12.011), abs=.001)
            patch = patch_3d(prepared.json()['recipe'], layers=1)
            lines = patch['exports']['xyz'].splitlines()
            assert ':isotope:I:1' in lines[1]
            assert sum(row.split()[-1] == '13' for row in lines[2:]) == sum(s.get('isotope') == 13 for s in patch['viewer']['atoms']) > 0
            copies = structure_exports(parsed, (2, 1, 1))
            assert sum(s.get('isotope') == 13 for s in import_cif(copies['exports']['cif'])['atoms']) == 8
            assert ':isotope:I:1' in copies['exports']['xyz'] and 'pdb' not in copies['exports']
        network = inspect_periodic_network(parsed)
        assert network['dimensions'] == [2] and network['quotient_component_count'] == 1
        assert network['short_contact_count'] == 0 and not network['threshold_sensitivity']['classification_changed']
        def keys(edges):
            return {edge_key(e['a'], e['b'], e['translation']) for e in edges}
        assert keys(network['edges']) == keys(result['cell_validation']['periodic_bonds'])
