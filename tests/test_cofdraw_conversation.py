"""Conversation continuity without weakening chemical or privacy constraints."""
import json
from pathlib import Path
from fastapi.testclient import TestClient
import pytest

from cof_vision.cofdraw_conversation import describe_turn, model_turn_text
from cof_vision.cofdraw_desktop import create_desktop_app
from cof_vision.cofdraw_chemistry import prepare_template


@pytest.fixture(scope='module')
def client():
    with TestClient(create_desktop_app(Path(__file__).resolve().parents[1] / 'data'), raise_server_exceptions=False) as value:
        yield value


@pytest.fixture(scope='module')
def library(client):
    return client.get('/api/drawing/library').json()


def follow(text, previous, library):
    return describe_turn(text, library, previous['conversation']['context'])


def test_followups_retain_pair_and_revalidate_actual_changed_chemistry(library):
    first = describe_turn('Draw TFB and TAPB', library)
    assert first['status'] == 'ready'
    pastel = follow('Make it pastel', first, library)
    assert pastel['status'] == 'ready', pastel['unresolved']
    assert pastel['selection']['monomer_a_id'] == 'tfb'
    assert pastel['style_id'] == 'pastel_blocks'
    changed = follow('Change B to DETH', pastel, library)
    assert changed['status'] == 'ready', changed['unresolved']
    assert changed['prepare_request']['linkage'] == 'acylhydrazone'
    assert changed['selection']['monomer_b_id'] == 'deth'
    assert prepare_template(**changed['prepare_request'])['instances']


def test_incremental_pair_and_explicit_replacement(library):
    partial = describe_turn('Use TFB', library)
    assert partial['status'] == 'needs_input'
    complete = follow('Use TAPB', partial, library)
    assert complete['status'] == 'ready', complete['unresolved']
    changed = follow('Replace TAPB with Pa', complete, library)
    assert changed['status'] == 'ready', changed['unresolved']
    assert changed['selection']['monomer_a_id'] == 'tfb'
    assert changed['selection']['monomer_b_id'] == 'pa'
    bad = follow('Replace DETH with Pa', complete, library)
    assert bad['status'] != 'ready' and bad['prepare_request'] is None
    constrained = follow('Replace TAPB without oxygen with Pa', complete, library)
    assert constrained['status'] != 'ready' and constrained['prepare_request'] is None


def test_source_appearance_preserves_identity_but_new_chemistry_cannot(library):
    source = describe_turn('Open TpPa', library)
    same = follow('Use wavy bonds', source, library)
    assert same['status'] == 'ready', same['unresolved']
    assert same['source_recipe_id'] == 'tppa-pore-v1'
    changed = follow('Change A to TFB', same, library)
    assert changed['status'] == 'ready', changed['unresolved']
    assert changed['source_recipe_id'] is None
    assert changed['prepare_request']
    assert any('source identity' in row for row in changed['assumptions'])


@pytest.mark.parametrize('text', ['Draw TFB and TAPB with a 2 nm pore', 'Draw TFB and TAPB with DFT', 'Use FooUnknown and TFB'])
def test_unknown_or_physical_constraints_survive_later_style_turn(text, library):
    first = describe_turn(text, library)
    assert first['status'] != 'ready'
    next_turn = follow('Make it pastel', first, library)
    assert next_turn['status'] != 'ready' and next_turn['prepare_request'] is None
    assert next_turn['conversation']['context']['pending']
    reset = follow('Start over', next_turn, library)
    assert reset['operation'] == 'reset' and not reset['conversation']['context']['pending']


def test_unknown_followup_not_downgraded_to_warning_by_inherited_context(library):
    base = describe_turn('Draw TFB and TAPB', library)
    unknown = follow('Use MythicalAmine', base, library)
    assert unknown['status'] == 'needs_input'
    assert any('MythicalAmine' in row for row in unknown['unresolved'])
    negation = follow('Use pastel without oxygen', base, library)
    assert negation['status'] != 'ready'


def test_explicit_functionality_can_be_completed_but_not_lost_between_turns(library):
    partial = describe_turn('Use a tritopic amine', library)
    assert partial['status'] != 'ready'
    assert partial['conversation']['context']['requirements'] == ['amine with 3 reactive sites']
    valid = follow('Draw TFB and TAPB', partial, library)
    assert valid['status'] == 'ready', valid['unresolved']
    changed = follow('Change B to Pa', valid, library)
    assert changed['status'] != 'ready'
    assert any('3-site amine' in row for row in changed['unresolved'])
    later = follow('Make it pastel', changed, library)
    assert later['status'] != 'ready'


def test_current_custom_drawing_appearance_is_a_bound_local_intent(library):
    context = {'origin': 'current_drawing', 'token': 'revision_12', 'selection': {'style_id': 'standard'}}
    result = describe_turn('Make the current drawing pastel', library, context)
    assert result['status'] == 'ready', result['unresolved']
    assert result['operation'] == 'style_only' and result['context_token'] == 'revision_12'
    assert result['style_id'] == 'pastel_blocks'
    assert result['prepare_request'] is None and result['source_recipe_id'] is None
    assert set(result['selection']) <= {'style_id', 'palette'}
    assert 'smiles' not in json.dumps(result).lower()
    invalid = describe_turn('Use wavy bonds with a tritopic amine', library, context)
    assert invalid['status'] != 'ready'
    crystal = describe_turn('Make a 3D crystal', library, context)
    assert crystal['status'] == 'unsupported'


@pytest.mark.parametrize('message,action', [('How do I check SA scores?', 'monomer_library'), ('How do I measure pores?', 'science_porosity'), ('Where can I make 3D layers?', 'science_model'), ('How can I run GCMC?', 'science_adsorption')])
def test_local_help_has_real_action_ids_and_no_invented_results(message, action, library):
    result = describe_turn(message, library)
    assert result['operation'] == 'help' and result['status'] == 'needs_input'
    assert result['actions'][0]['id'] == action
    assert result['prepare_request'] is None


def test_api_context_forbids_private_data_and_redacts_bad_credentials(client):
    endpoint = '/api/drawing/describe/chat'
    for context in [{'recipe': {'private': 'PRIVATE'}}, {'selection': {'smiles': 'PRIVATE'}}, {'origin': 'current_drawing', 'token': '../PRIVATE'}]:
        response = client.post(endpoint, json={'text': 'Make it pastel', 'context': context})
        assert response.status_code == 422 and 'PRIVATE' not in response.text
    response = client.post(endpoint, json={'text': 'Draw', 'mode': 'llm', 'llm': {'api_key': ['private-secret']}})
    assert response.status_code == 422 and 'private-secret' not in response.text
    response = client.post(endpoint, json={'text': 'private-secret', 'mode': 'llm', 'llm': {'model': 'fixture', 'api_key': 'private-secret'}})
    assert response.status_code == 422 and 'private-secret' not in response.text


def test_optional_model_gets_only_new_message_and_allowlisted_ids(client, library, monkeypatch):
    from cof_vision import cofdraw_llm
    calls = []
    def selector(text, catalog, **kwargs):
        calls.append(json.loads(text))
        assert '"smiles"' not in json.dumps(catalog).lower()
        return {'selection': {'monomer_a_id': 'tfb', 'monomer_b_id': 'tapb', 'style_id': 'pastel_blocks'}, 'unresolved': [], 'explanation': 'Fixture selection.'}
    monkeypatch.setattr(cofdraw_llm, 'select_with_llm', selector)
    context = {'selection': {'monomer_a_id': 'tfb', 'monomer_b_id': 'tapb'}, 'pending': ['PRIVATE local pending constraint']}
    response = client.post('/api/drawing/describe/chat', json={'text': 'Make it pastel', 'context': context, 'mode': 'llm', 'llm': {'model': 'fixture'}})
    assert response.status_code == 200, response.text
    assert len(calls) == 1
    assert 'PRIVATE' not in json.dumps(calls)
    assert calls[0]['current_catalog_selection']['monomer_a_id'] == 'tfb'
    assert response.json()['status'] != 'ready'  # Pending local constraint still blocks Apply.
    assert not any(key in model_turn_text('hello', context, library) for key in ['PRIVATE', 'token'])


def test_help_does_not_spend_an_optional_model_request(client, monkeypatch):
    from cof_vision import cofdraw_llm
    monkeypatch.setattr(cofdraw_llm, 'select_with_llm', lambda *a, **k: pytest.fail('Local help contacted provider'))
    response = client.post('/api/drawing/describe/chat', json={'text': 'How do I check SA scores?', 'mode': 'llm', 'llm': {'model': 'fixture'}})
    assert response.status_code == 200 and response.json()['interpretation']['mode'] == 'rules'


@pytest.mark.parametrize('text', ['data_private\n_cell_length_a 12\nloop_\n_atom_site_label C1', '{"recipe":{"private":"LOCAL_PROJECT"}}'])
def test_pasted_structure_payload_does_not_reach_provider(client, monkeypatch, text):
    from cof_vision import cofdraw_llm
    monkeypatch.setattr(cofdraw_llm, 'select_with_llm', lambda *a, **k: pytest.fail('Structured file contents reached provider'))
    response = client.post('/api/drawing/describe/chat', json={'text': text, 'mode': 'llm', 'llm': {'model': 'fixture'}})
    assert response.status_code == 422
    assert 'LOCAL_PROJECT' not in response.text
