"""Level-1 geometry API contracts, saved provenance, and legacy compatibility."""
from copy import deepcopy
import hashlib
import json

from fastapi.testclient import TestClient
import pytest

from cof_vision.cofdraw_desktop import create_desktop_app


def structure():
    return {'periodic': True, 'cell': [[10.,0,0],[0,10.,0],[0,0,10.]],
            'atoms': [{'element': 'C', 'fractional': [.5,.5,.5]}],
            'provenance': {'kind': 'analytic numerical fixture'}}


def digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, allow_nan=False).encode()).hexdigest()


@pytest.fixture
def client(tmp_path):
    with TestClient(create_desktop_app(tmp_path), raise_server_exceptions=False) as instance:
        yield instance


def test_level1_reports_bind_all_effective_inputs_and_roundtrip_exactly(client):
    payload = {'structure': structure(), 'probe_radius': 1., 'grid_spacing': 1., 'surface_samples': 16, 'refinement': True}
    response = client.post('/api/drawing/science/geometry', json=payload)
    assert response.status_code == 200, response.text
    result = response.json()
    assert result['operation'] == 'geometry' and result['engine'] == 'cofdraw-grid-v1'
    assert result['inputs'] == payload and result['input_sha256'] == digest(payload)
    assert result['connectivity']['status'] == 'certified_channel'
    assert result['pld_bounds_angstrom'][0] == result['pld_angstrom']
    assert result['lcd_bounds_angstrom'][0] == result['lcd_angstrom']
    assert result['convergence']['performed'] is True
    assert result['transport_accessible_surface_m2_g'] == result['local_probe_exposed_surface_m2_g']
    assert result['bet_surface_area_m2_g'] is None
    assert client.get('/api/drawing/science/reports/'+result['report_id']).json() == result
    assert client.get(result['downloads']['json']).json() == result
    assert client.get('/api/drawing/science/reports?limit=1').json()['reports'][0]['operation'] == 'geometry'


def test_capabilities_expose_current_defaults_and_repeatable_numerical_results(client):
    capability = client.get('/api/drawing/science/capabilities').json()['operations']['geometry']
    assert capability['method'] == 'POST' and capability['path'] == '/api/drawing/science/geometry'
    properties = capability['request_schema']['properties']
    assert properties['grid_spacing']['default'] == .6
    assert properties['surface_samples']['default'] == 64
    assert properties['refinement']['default'] is False
    responses = [client.post(capability['path'], json={'structure': structure()}) for _ in range(2)]
    assert all(response.status_code == 200 for response in responses)
    a,b = [response.json() for response in responses]
    assert a['report_id'] != b['report_id']
    assert a['structure_sha256'] == b['structure_sha256'] and a['input_sha256'] == b['input_sha256']
    for key in ['pld_bounds_angstrom','lcd_bounds_angstrom','void_fraction','transport_accessible_surface_m2_g','method','connectivity']:
        assert a[key] == b[key]


@pytest.mark.parametrize('changes', [{'grid_spacing':0}, {'grid_spacing':float('inf')}, {'probe_radius':-1},
    {'surface_samples':True}, {'refinement':'yes'}, {'samples':4000}])
def test_invalid_geometry_controls_are_rejected(client, changes):
    if changes.get('grid_spacing') == float('inf'):
        response = client.post('/api/drawing/science/geometry', content=json.dumps({'structure':structure(),**changes}),
                               headers={'Content-Type':'application/json'})
    else:
        response = client.post('/api/drawing/science/geometry', json={'structure':structure(),**changes})
    assert response.status_code == 422


@pytest.mark.parametrize('mutation', ['finite','occupancy','coincident','bad_cell'])
def test_geometry_requires_resolved_periodic_full_cell_structure(client, mutation):
    value = deepcopy(structure())
    if mutation == 'finite': value['periodic'] = False
    elif mutation == 'occupancy': value['atoms'][0]['occupancy'] = .5
    elif mutation == 'coincident': value['atoms'].append(deepcopy(value['atoms'][0]))
    else: value['cell'][2][2] = -10
    response = client.post('/api/drawing/science/geometry', json={'structure':value})
    assert response.status_code == 422 and response.json()['detail']


def test_previous_sampling_porosity_endpoint_remains_separately_available(client):
    response = client.post('/api/drawing/science/porosity', json={'structure':structure(),'samples':1000})
    assert response.status_code == 200, response.text
    result = response.json()
    assert result['operation'] == 'porosity' and result['kind'] == 'periodic_hard_sphere_geometry'
    assert result['method']['requested_volume_samples'] == 1000
    assert 'largest_sampled_included_sphere_diameter_angstrom' in result
    assert 'pld_angstrom' not in result
