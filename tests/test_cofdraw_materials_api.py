"""Real materials API workflows, report provenance, and safe CSV interchange."""
import csv
from copy import deepcopy
import hashlib
import io
import json
import math
from pathlib import Path

from fastapi.testclient import TestClient
import pytest

from cof_vision.cofdraw_desktop import create_desktop_app
from cof_vision.cofdraw_science import import_cif, normalize_structure


@pytest.fixture
def materials(tmp_path):
    with TestClient(create_desktop_app(tmp_path), raise_server_exceptions=False) as client:
        yield client, tmp_path


def carbon_cell():
    return {'periodic': True, 'cell': [[5, 0, 0], [0, 5, 0], [0, 0, 5]],
            'atoms': [{'element': 'C', 'fractional': [.2, .3, .4]}],
            'provenance': {'kind': 'test_geometry', 'label': 'single-site numerical fixture'}}


def digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, allow_nan=False).encode()).hexdigest()


def request(client, path, payload):
    response = client.post('/api/drawing/science/'+path, json=payload)
    assert response.status_code == 200, response.text
    return response.json()


@pytest.mark.parametrize('name,dimension,components', [('hcb-imine-aa', 2, 1), ('hcb-imine-ab', 2, 2), ('dia-imine', 3, 1)])
def test_inspection_preserves_structure_hash_and_actual_periodic_dimensions(materials, name, dimension, components):
    client, _ = materials
    path = Path(__file__).resolve().parents[1]/'data/crystal-models'/f'{name}.cif'
    structure = import_cif(path.read_text('utf-8'))
    result = request(client, 'inspect', {'structure': structure})
    assert result['operation'] == 'inspect'
    assert result['summary']['structure_sha256'] == result['network']['structure_sha256'] == digest(normalize_structure(structure))
    assert result['network']['dimensions'] == [dimension]
    assert result['network']['quotient_component_count'] == components
    assert {row['periodic_rank'] for row in result['network']['components']} == {dimension}
    assert result['network']['short_contact_count'] == 0
    assert result['network']['explicit_hydrogen_count'] > 0
    assert any('inferred' in text.lower() for text in result['network']['warnings'])
    assert result['input_sha256'] == digest(result['inputs'])
    assert result['inputs']['structure'] == structure


def test_report_readback_export_and_listing_reproduce_exact_effective_inputs(materials):
    client, root = materials
    payload = {'structure': carbon_cell(), 'bond_scale': 1.1, 'tolerance_angstrom': .12}
    result = request(client, 'inspect', payload)
    assert result['inputs'] == payload
    assert result['input_sha256'] == digest(payload)
    assert result['network']['method']['bond_scale'] == 1.1
    assert result['network']['method']['tolerance_angstrom'] == .12
    assert result['report_schema_version'] == 'cofdraw-science-1'
    assert result['engine_version'] and result['created_at_utc'].endswith('+00:00')
    report_id = result['report_id']
    response = client.get('/api/drawing/science/reports/'+report_id)
    assert response.status_code == 200 and response.json() == result
    assert client.get(result['downloads']['json']).json() == result
    listed = client.get('/api/drawing/science/reports?limit=1').json()['reports']
    assert len(listed) == 1 and listed[0]['report_id'] == report_id
    assert listed[0]['input_sha256'] == digest(payload)
    changed = deepcopy(payload); changed['structure']['atoms'][0]['fractional'][0] += .1
    other = request(client, 'inspect', changed)
    assert other['input_sha256'] != result['input_sha256']
    assert other['report_id'] != report_id
    assert client.get('/api/drawing/science/reports/'+report_id).json() == result


def test_mixed_valid_invalid_batch_preserves_each_status_and_downloads_safe_table(materials):
    client, _ = materials
    invalid = carbon_cell(); invalid['cell'][2][2] = -5
    entries = [{'name': 'Model, "quoted"', 'structure': carbon_cell()},
               {'name': 'invalid geometry', 'structure': invalid}]
    result = request(client, 'screen', {'entries': entries, 'include_porosity': True, 'samples': 1000})
    assert result['operation'] == 'screen' and result['successful'] == 1
    assert [row['status'] for row in result['entries']] == ['ok', 'error']
    good, bad = result['entries']
    assert good['summary']['formula'] == 'C'
    assert good['inferred_network_dimensions'] == [0]
    assert good['porosity']['method']['requested_volume_samples'] == 1000
    assert 'right-handed' in bad['error']
    assert result['input_sha256'] == digest(result['inputs'])
    response = client.get(result['downloads']['csv'])
    assert response.status_code == 200 and response.headers['content-type'].startswith('text/csv')
    rows = list(csv.DictReader(io.StringIO(response.text)))
    assert len(rows) == 2
    assert rows[0]['name'] == entries[0]['name'] and rows[0]['formula'] == 'C'
    assert rows[0]['inferred_network_dimensions'] == '0'
    assert float(rows[0]['geometric_void_volume_cm3_g']) == pytest.approx(good['porosity']['geometric_void_volume_cm3_g'])
    assert rows[1]['status'] == 'error' and rows[1]['error'] == bad['error']
    assert not rows[1]['density_g_cm3']


def test_csv_names_are_inert_even_after_spreadsheet_whitespace_normalization(materials):
    client, _ = materials
    names = ['=1+1', '+2', '-3', '@SUM(1,2)', '\t=1+1', '\r=1+1', '\n=1+1', '  =1+1']
    result = request(client, 'screen', {'entries': [{'name': name, 'structure': carbon_cell()} for name in names]})
    rows = list(csv.DictReader(io.StringIO(client.get(result['downloads']['csv']).text)))
    assert result['successful'] == len(names)
    assert [row['name'] for row in result['entries']] == names  # JSON keeps the original evidence.
    assert len(rows) == len(names)
    for row, original in zip(rows, names):
        assert row['name'].startswith("'"), f'Unsafe spreadsheet name: {original!r}'
        assert row['name'][1:] == original
        assert row['formula'] == 'C'


def test_pxrd_route_has_a_real_profile_bragg_peak_method_and_numeric_csv(materials):
    client, _ = materials
    result = request(client, 'pxrd', {'structure': carbon_cell(), 'two_theta_min': 10, 'two_theta_max': 45,
                                     'wavelength': 1.5406, 'step': .02, 'fwhm': .15})
    assert result['kind'] == 'simulated_powder_xrd' and result['quality'] == 'model_simulation'
    assert result['operation'] == 'pxrd'
    assert result['structure_sha256'] == digest(normalize_structure(carbon_cell()))
    expected_first = math.degrees(2*math.asin(1.5406/(2*5)))
    assert result['peaks'][0]['two_theta_degrees'] == pytest.approx(expected_first, abs=1e-6)
    assert result['peaks'][0]['multiplicity'] == 6
    assert max(result['profile']['intensity']) == pytest.approx(100)
    assert len(result['method']['coefficient_sha256']) == 64
    assert 'neutral' in ' '.join(result['warnings']).lower()
    assert 'experimental' in ' '.join(result['warnings']).lower()
    rows = list(csv.DictReader(io.StringIO(client.get(result['downloads']['csv']).text)))
    assert len(rows) == len(result['profile']['two_theta_degrees'])
    for index in (0, len(rows)//2, len(rows)-1):
        assert float(rows[index]['two_theta_degrees']) == pytest.approx(result['profile']['two_theta_degrees'][index], abs=1e-7)
        assert float(rows[index]['relative_intensity']) == pytest.approx(result['profile']['intensity'][index], abs=1e-5)


def test_pxrd_bad_ranges_fail_clearly_and_narrow_width_works_with_resolved_sampling(materials):
    client, _ = materials
    base = {'structure': carbon_cell(), 'two_theta_min': 10, 'two_theta_max': 30}
    for options in ({'two_theta_min': 40}, {'fwhm': .02, 'step': .02}, {'wavelength': .1}):
        response = client.post('/api/drawing/science/pxrd', json={**base, **options})
        assert response.status_code == 422 and response.json()['detail']
    result = request(client, 'pxrd', {**base, 'fwhm': .02, 'step': .01})
    assert result['parameters']['actual_step_degrees'] <= result['parameters']['fwhm_degrees']/2+1e-12
    assert result['peak_count'] > 0


def test_capability_schemas_drive_real_requests_and_bound_batch_inputs(materials):
    client, _ = materials
    response = client.get('/api/drawing/science/capabilities')
    assert response.status_code == 200
    capabilities = response.json()
    assert capabilities['schema_version'] == 'cofdraw-capabilities-1'
    assert {'inspect', 'pxrd', 'screen', 'query', 'model', 'monomer', 'porosity', 'bet'} <= capabilities['operations'].keys()
    operation = capabilities['operations']['inspect']
    schema = operation['request_schema']
    assert operation['method'] == 'POST' and schema['additionalProperties'] is False
    assert schema['required'] == ['structure']
    payload = {name: field['default'] for name, field in schema['properties'].items() if 'default' in field}
    payload['structure'] = carbon_cell()
    result = client.request(operation['method'], operation['path'], json=payload)
    assert result.status_code == 200, result.text
    assert result.json()['inputs'] == payload
    assert client.post(operation['path'], json={**payload, 'unsupported_option': True}).status_code == 422
    field = schema['properties']['bond_scale']
    assert client.post(operation['path'], json={**payload, 'bond_scale': field['maximum']+.1}).status_code == 422
    screen = capabilities['operations']['screen']['request_schema']
    assert screen['properties']['entries']['maxItems'] == 8
    assert screen['properties']['entries']['items']['$ref'] == '#/$defs/ScreenEntry'
    assert screen['$defs']['ScreenEntry']['required'] == ['structure']
    assert client.post(capabilities['operations']['screen']['path'], json={'entries': []}).status_code == 422
    assert client.post(capabilities['operations']['screen']['path'], json={'entries': [{'structure': carbon_cell()}]*9}).status_code == 422
    assert client.get(capabilities['reports']['list']).status_code == 200
    assert client.get(capabilities['reports']['read'].replace('{report_id}', result.json()['report_id'])).status_code == 200


def test_recent_reports_skip_unfinished_nonobject_and_ordinary_drawing_files(materials):
    client, root = materials
    valid = request(client, 'inspect', {'structure': carbon_cell()})
    folder = root/'framework-drawing'/'exports'
    for identifier, content in [('a'*32, '{'), ('b'*32, '[]'), ('c'*32, '{"recipe":{}}')]:
        path = folder/identifier; path.mkdir()
        (path/(identifier+'.json')).write_text(content, encoding='utf-8')
    response = client.get('/api/drawing/science/reports')
    assert response.status_code == 200, response.text
    assert [row['report_id'] for row in response.json()['reports']] == [valid['report_id']]
    for identifier in ('a'*32, 'b'*32, 'c'*32, 'd'*32, 'not-a-report'):
        assert client.get('/api/drawing/science/reports/'+identifier).status_code in {404, 422}
    assert client.get('/api/drawing/science/reports?limit=0').status_code == 422
    assert client.get('/api/drawing/science/reports?limit=101').status_code == 422
