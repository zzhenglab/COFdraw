"""Description plans cross the same chemistry and identity boundaries as manual designs."""
import json
from pathlib import Path
import shutil

from fastapi.testclient import TestClient
import pytest

from cof_vision.cofdraw_desktop import create_desktop_app
from cof_vision.drawing_api import DescribeDrawing


@pytest.fixture
def describe_client(tmp_path):
    data = Path(__file__).resolve().parents[1] / 'data/framework-drawing'
    shutil.copytree(data, tmp_path / 'framework-drawing', ignore=shutil.ignore_patterns('exports'))
    with TestClient(create_desktop_app(tmp_path), raise_server_exceptions=False) as client:
        yield client, tmp_path


def plan(client, text, **kwargs):
    response = client.post('/api/drawing/describe', json={'text': text, **kwargs})
    assert response.status_code == 200, response.text
    return response.json()


def test_rules_plan_builds_real_graph_and_export_without_source_identity(describe_client):
    client, root = describe_client
    result = plan(client, 'Draw TFB and TAPB with an imine linkage in pastel colors')
    assert result['status'] == 'ready', result
    assert result['interpretation'] == {'mode': 'rules'}
    assert result['source_recipe_id'] is None
    prepared = client.post('/api/drawing/prepare', json=result['prepare_request'])
    assert prepared.status_code == 200, prepared.text
    recipe = prepared.json()['recipe']
    assert recipe['ground_truth'] is False and recipe['cofex_id'] is None
    rendered = client.post('/api/drawing/render', json={'recipe': recipe,
        'style_id': result['style_id'], 'style_options': result['style_options'], 'size': 600})
    assert rendered.status_code == 200, rendered.text
    drawing = rendered.json()
    assert drawing['audit']['atoms'] > 50
    for ext in ('png', 'svg', 'mol', 'sdf', 'json'):
        assert client.get(drawing['downloads'][ext]).status_code == 200


def test_source_plan_preserves_existing_hash_binding(describe_client):
    client, _ = describe_client
    result = plan(client, 'Open the Pyurea source drawing')
    assert result['status'] == 'ready', result
    assert result['prepare_request'] is None
    assert result['source_recipe_id'] == 'pyurea-source-layout-v1'
    original = client.get('/api/drawing/sources/'+result['source_recipe_id'])
    assert original.status_code == 200
    assert original.json()['source_sha256']
    assert original.json()['status'] == 'source_bound_provisional'


def test_catalog_and_capabilities_have_discoverable_actual_request_schema(describe_client):
    client, _ = describe_client
    catalog = client.get('/api/drawing/describe/catalog')
    assert catalog.status_code == 200, catalog.text
    data = catalog.json()
    assert len(data['demos']) >= 6 and len(data['monomers']) >= 49
    assert all(row['text'] and row['id'] for row in data['demos'])
    capabilities = client.get('/api/drawing/science/capabilities').json()
    assert capabilities['description']['request_schema'] == DescribeDrawing.model_json_schema()
    assert capabilities['description']['path'] == '/api/drawing/describe'
    assert any(row['record_id'] == 'describe-cof' for row in capabilities['workflows'])


def test_llm_mode_uses_provider_ids_then_actual_local_chemistry(describe_client, monkeypatch):
    from cof_vision import cofdraw_llm
    client, root = describe_client
    calls = []
    def selector(text, library, **kwargs):
        calls.append((text, kwargs))
        return {'selection': {'monomer_a_id': 'tfb', 'monomer_b_id': 'tapb', 'linkage': 'imine'},
                'explanation': 'Matched the two named library monomers.', 'unresolved': []}
    monkeypatch.setattr(cofdraw_llm, 'select_with_llm', selector)
    result = plan(client, 'Draw TFB and TAPB with an imine linkage', mode='llm',
                  llm={'model': 'test-model', 'api_key': 'private-test-secret'})
    assert len(calls) == 1 and calls[0][1]['api_key'] == 'private-test-secret'
    assert result['status'] == 'ready', result
    assert result['prepare_request']['linkage'] == 'imine'
    assert 'private-test-secret' not in json.dumps(result)
    assert 'api_key' not in json.dumps(result)
    assert not (root / 'framework-drawing/exports').exists()


def test_llm_unresolved_or_invented_ids_cannot_be_drawn(describe_client, monkeypatch):
    from cof_vision import cofdraw_llm
    client, _ = describe_client
    monkeypatch.setattr(cofdraw_llm, 'select_with_llm', lambda *a, **k: {
        'selection': {'source_recipe_id': 'pyurea-source-layout-v1'},
        'explanation': 'Ambiguous request.', 'unresolved': ['Confirm the source material.']})
    result = plan(client, 'Open Pyurea source drawing', mode='llm', llm={'model': 'test'})
    assert result['status'] != 'ready' and result['prepare_request'] is None
    assert result['source_recipe_id'] is None
    invalid = client.post('/api/drawing/describe', json={'text': 'Draw TFB and TAPB',
        'selection': {'monomer_a_id': 'invented-monomer'}})
    assert invalid.status_code == 422 or invalid.json()['status'] != 'ready'


def test_reviewed_interpretation_keeps_unmatched_words_visible_and_hard_constraints(describe_client, monkeypatch):
    from cof_vision import cofdraw_llm
    client, _ = describe_client
    selection = {'monomer_a_id': 'tfb', 'monomer_b_id': 'tapb', 'linkage': 'imine'}
    prompt = 'Draw TFB and TAPB with an imine linkage, beautifully luminous'
    initial = plan(client, prompt)
    assert initial['status'] == 'needs_input'
    reviewed = plan(client, prompt, selection=selection)
    assert reviewed['status'] == 'ready'
    assert any('luminous' in warning for warning in reviewed['warnings'])
    monkeypatch.setattr(cofdraw_llm, 'select_with_llm', lambda *a, **k: {
        'selection': selection, 'explanation': 'Matched the named monomers.', 'unresolved': []})
    interpreted = plan(client, prompt, mode='llm', llm={'model': 'test'})
    assert interpreted['status'] == 'ready'
    assert any('luminous' in warning for warning in interpreted['warnings'])
    contradictory = plan(client, 'Draw TFB and TAPB with a boronate ester linkage',
                         mode='llm', llm={'model': 'test'})
    assert contradictory['status'] != 'ready' and contradictory['prepare_request'] is None


@pytest.mark.parametrize('payload', [
    {'text': '   '}, {'text': 'x' * 4001},
    {'text': 'draw', 'mode': 'execute'},
    {'text': 'draw', 'mode': 'llm'},
    {'text': 'draw', 'llm': {'model': 'test'}},
    {'text': 'draw', 'mode': 'llm', 'llm': {'model': 'test'}, 'selection': {}},
    {'text': 'draw', 'selection': {'smiles': 'invented'}},
])
def test_invalid_requests_do_not_call_provider(describe_client, monkeypatch, payload):
    from cof_vision import cofdraw_llm
    client, _ = describe_client
    def forbidden(*a, **k):
        pytest.fail('An invalid description request contacted the provider')
    monkeypatch.setattr(cofdraw_llm, 'select_with_llm', forbidden)
    assert client.post('/api/drawing/describe', json=payload).status_code == 422


def test_untrusted_origin_cannot_trigger_paid_api(describe_client, monkeypatch):
    from cof_vision import cofdraw_llm
    client, _ = describe_client
    def forbidden(*a, **k):
        pytest.fail('Cross-origin write reached provider')
    monkeypatch.setattr(cofdraw_llm, 'select_with_llm', forbidden)
    response = client.post('/api/drawing/describe', json={'text': 'Draw a COF',
        'mode': 'llm', 'llm': {'model': 'test', 'api_key': 'private-test-secret'}},
        headers={'Origin': 'https://untrusted.example'})
    assert response.status_code == 403 and 'private-test-secret' not in response.text


def test_secret_not_serialized_or_represented_in_model():
    model = DescribeDrawing(text='Draw a COF', mode='llm', llm={
        'model': 'test', 'api_key': 'private-test-secret'})
    assert 'private-test-secret' not in repr(model)
    assert 'private-test-secret' not in model.model_dump_json()


@pytest.mark.parametrize('llm', [
    'private-test-secret', {'model': 'test', 'api_key': ['private-test-secret']},
    {'model': 'test', 'api_key': 'private-test-secret', 'unknown': 'private-test-secret'},
])
def test_malformed_credentials_are_not_echoed_by_validation(describe_client, llm):
    client, _ = describe_client
    response = client.post('/api/drawing/describe', json={
        'text': 'Draw a COF', 'mode': 'llm', 'llm': llm})
    assert response.status_code == 422
    assert 'private-test-secret' not in response.text
