"""Check chosen stacking transforms and periodic-copy hypotheses numerically."""
from pathlib import Path
import numpy as np
import pytest
from fastapi import FastAPI
from fastapi.testclient import TestClient

from cof_vision.cofdraw_arrangements import network_hypothesis
from cof_vision.cofdraw_chemistry import prepare_template
from cof_vision.cofdraw_science import import_cif, normalize_structure, patch_3d
from cof_vision.drawing_api import create_drawing_router
from cof_vision.framework_drawing import VERSION


def recipe():
    prepared = prepare_template('O=Cc1cc(C=O)cc(C=O)c1', 'Nc1ccc(N)cc1', 'imine', 'TFB', 'Pa', 12, 'hcb_3_2')
    return {'schema_version': VERSION, 'recipe_id': 'stack-test', 'layout_refinement': 'bond_angle_pose_v2', **prepared}


@pytest.mark.parametrize('mode,factors', [('AA',[0,0,0,0]), ('AB',[0,1,0,1]), ('ABC',[0,1,2,0]), ('slipped',[0,1,2,3])])
def test_layer_sequence_preserves_geometry(mode, factors):
    result = patch_3d(recipe(), layers=4, stacking=mode, slip=[1.2,.4])
    atoms = result['viewer']['atoms']
    count = len(atoms)//4
    base = np.array([[a[k] for k in ('x','y','z')] for a in atoms[:count]])
    for layer, factor in enumerate(factors):
        xyz = np.array([[a[k] for k in ('x','y','z')] for a in atoms[layer*count:(layer+1)*count]])
        assert xyz-base == pytest.approx(np.tile([1.2*factor,.4*factor,3.4*layer], (count,1)))
    assert all(a['instance'] and a['component'] for a in atoms)
    assert not result['periodic']


def test_twisted_stack_preserves_distances_and_uses_chosen_angle():
    result = patch_3d(recipe(), layers=3, stacking='twisted', twist_degrees=30)
    atoms = result['viewer']['atoms']; count = len(atoms)//3
    layers = [np.array([[a[k] for k in ('x','y','z')] for a in atoms[i*count:(i+1)*count]]) for i in range(3)]
    for i, xyz in enumerate(layers):
        angle = np.radians(i*30)
        matrix = np.array([[np.cos(angle), -np.sin(angle),0], [np.sin(angle),np.cos(angle),0],[0,0,1]])
        assert xyz == pytest.approx(layers[0]@matrix.T+[0,0,i*3.4])
        assert np.linalg.norm(xyz[:,None]-xyz[None,:], axis=-1) == pytest.approx(np.linalg.norm(layers[0][:,None]-layers[0][None,:], axis=-1))


def structure():
    return {'periodic':True, 'cell':[[10,0,0],[0,10,0],[0,0,10]],
            'atoms':[{'element':'C','fractional':[.1,.1,.1]}, {'element':'N','fractional':[.24,.1,.1]}]}


def test_network_copies_preserve_atoms_fractional_shifts_and_group_metadata():
    result = network_hypothesis(structure(), [[0,0,0],[.5,.5,.5]], repeats=[2,1,1])
    assert len(result['structure']['atoms']) == 4
    assert len(result['viewer']['atoms']) == 8
    assert [a['network'] for a in result['viewer']['atoms']] == [0,0,1,1]*2
    coordinates = np.array([a['fractional'] for a in result['structure']['atoms']])
    assert coordinates[2:]-coordinates[:2] == pytest.approx(np.full((2,3),.5))
    assert result['interpenetration_degree'] is None
    assert result['structure']['provenance']['interpenetration_degree'] is None
    assert result['minimum_intercopy_distance_angstrom'] == pytest.approx(np.sqrt(3.6**2+5**2+5**2))
    assert len(import_cif(result['exports']['cif'])['atoms']) == 8
    assert normalize_structure(result['structure'])['atoms'] == result['structure']['atoms']


@pytest.mark.parametrize('offsets', [[[0,0,0],[1,0,0]], [[0,0,0],[.01,0,0]], [[0,0]], [[0,0,float('nan')]]])
def test_network_bad_offsets_and_collisions_are_rejected(offsets):
    with pytest.raises(ValueError):
        network_hypothesis(structure(), offsets)


def test_network_api_exports_and_capability_contract(tmp_path):
    app=FastAPI(); app.include_router(create_drawing_router(tmp_path))
    client=TestClient(app)
    body={'structure':structure(),'offsets':[[0,0,0],[.5,.5,.5]],'repeats':[1,1,1]}
    response=client.post('/api/drawing/science/networks',json=body)
    assert response.status_code == 200, response.text
    result=response.json()
    assert result['operation']=='networks' and result['interpenetration_degree'] is None
    assert set(result['downloads'])=={'cif','xyz','pdb','json'}
    assert all(client.get(url).status_code==200 for url in result['downloads'].values())
    assert 'networks' in client.get('/api/drawing/science/capabilities').json()['operations']
    response=client.post('/api/drawing/science/networks',json={**body,'offsets':[[0,0,0],[.01,0,0]]})
    assert response.status_code==422 and 'collide' in response.text


def test_bundled_dia_hypothesis_and_colliding_offset():
    root = Path(__file__).resolve().parents[1]
    original = import_cif((root/'data/crystal-models/dia-imine.cif').read_text('utf-8'))
    result = network_hypothesis(original, [[0,0,0],[.5,.5,.5]])
    assert len(result['structure']['atoms']) == 292
    assert result['minimum_intercopy_distance_angstrom'] > 10
    assert result['interpenetration_degree'] is None
    with pytest.raises(ValueError, match='collide'):
        network_hypothesis(original, [[0,0,0],[.25,.25,.25]])
