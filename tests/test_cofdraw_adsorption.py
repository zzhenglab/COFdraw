import hashlib
import json
from pathlib import Path
import sys
import time
import zipfile

import pytest
from pydantic import ValidationError

from cof_vision.cofdraw_adsorption import AdsorptionJobs, AdsorptionRequest, mixture_selectivity, zero_loading_heat


class Backend:
    def __init__(self, delay=0, fail=False):
        self.delay, self.fail = delay, fail

    def engine_status(self):
        return {'available': True, 'version': 'test-worker'}

    def prepare_run(self, directory, structure, **settings):
        if structure.get('invalid'):
            raise ValueError('Unsupported structure')
        (directory/'engine-input.json').write_text(json.dumps(settings))
        return settings

    def command_for(self, directory):
        # Genuine isolated child process tests lifecycle; no invented RASPA results.
        return [sys.executable, '-c', f'import time;time.sleep({self.delay!r});raise SystemExit({int(self.fail)})']

    def parse_run(self, directory, manifest):
        return {'components': [{'gas': gas, 'loading_mol_kg': i+1., 'loading_confidence95_mol_kg': .1,
                                'henry_mol_kg_Pa': 1e-5*(i+1), 'henry_confidence95_mol_kg_Pa': 1e-6}
                               for i, gas in enumerate(manifest['gases'])]}


def await_job(manager, identifier):
    deadline = time.monotonic()+15
    while time.monotonic() < deadline:
        result = manager.status(identifier)
        if result['status'] in {'complete', 'failed', 'cancelled'}:
            return result
        time.sleep(.02)
    pytest.fail('Local test worker did not finish')


def test_jobs_run_distinct_widom_and_pressure_stages_and_preserve_provenance(tmp_path):
    manager = AdsorptionJobs(tmp_path, Backend())
    try:
        job = manager.start({'structure': {}, 'gases': ['CO2','N2'], 'mole_fractions':[.2,.8]})
        result = await_job(manager, job['job_id'])
        assert result['status'] == 'complete', result
        assert result['completed_stages'] == result['total_stages'] == 5
        directory = tmp_path/job['job_id']
        inputs = [json.loads(p.read_text()) for p in directory.glob('*/engine-input.json')]
        assert len(inputs) == 5
        assert len({p['seed'] for p in inputs}) == 5
        assert len([p for p in inputs if p['mode'] == 'widom']) == 3
        assert {p['temperature_K'] for p in inputs if p['mode'] == 'widom'} == {293.15,298.15,303.15}
        assert result['result']['henry_selectivity'][0]['selectivity'] == .5
        assert result['result']['isotherm'][0]['selectivity'][0]['selectivity'] == 2
        with zipfile.ZipFile(directory/f'{job["job_id"]}.zip') as archive:
            hashes=json.loads(archive.read('sha256.json'))
            assert 'widom/engine-input.json' in hashes
            assert 'pressure-01/worker.log' in hashes
            assert all(hashlib.sha256(archive.read(name)).hexdigest() == digest for name,digest in hashes.items())
    finally:
        manager.close()


def test_cancel_terminates_owned_worker_and_allows_next_run(tmp_path):
    manager=AdsorptionJobs(tmp_path, Backend(delay=5))
    try:
        job=manager.start({'structure':{},'calculate_heat':False,'pressures_Pa':[100000.]})
        with pytest.raises(RuntimeError, match='already running'):
            manager.start({'structure':{}})
        manager.cancel(job['job_id'])
        assert await_job(manager,job['job_id'])['status']=='cancelled'
        if manager._thread:
            manager._thread.join(2)
        assert manager._process is None
        manager.backend=Backend()
        new=manager.start({'structure':{},'calculate_heat':False,'pressures_Pa':[100000.]})
        assert await_job(manager,new['job_id'])['status']=='complete'
    finally:
        manager.close()


def test_failure_exports_logs_without_claiming_results(tmp_path):
    manager=AdsorptionJobs(tmp_path,Backend(fail=True))
    try:
        job=manager.start({'structure':{}})
        result=await_job(manager,job['job_id'])
        assert result['status']=='failed'
        assert result['result']['widom'] is None
        assert result['result']['isotherm']==[]
        assert result['downloads']['zip']
    finally:
        manager.close()


@pytest.mark.parametrize('change',[
    {'gases':['CO2','CO2'],'mole_fractions':[.5,.5]},
    {'gases':['CO2','N2'],'mole_fractions':[.1,.1]},
    {'fugacity_coefficients':[1,1]}, {'pressures_Pa':[100,100]},
    {'temperature_K':float('nan')}, {'cycles':0}, {'charge_mode':'guess'},
])
def test_invalid_run_settings_rejected(change):
    with pytest.raises(ValidationError):
        AdsorptionRequest.model_validate({'structure':{},**change})


def test_missing_sampling_does_not_report_infinite_selectivity():
    row=mixture_selectivity(['CO2','N2'],{'CO2':5,'N2':0},[.15,.85],'loading')[0]
    assert row['selectivity'] is None and row['reason']


def test_zero_loading_heat_recovers_analytic_vant_hoff_and_uncertainty():
    import math
    heat=25.
    def row(t):
        kh=1e-8*math.exp(heat/(.00831446261815324*t))
        return {'components':[{'gas':'CO2','henry_mol_kg_Pa':kh,'henry_confidence95_mol_kg_Pa':kh*.01}]}
    result=zero_loading_heat(row(295),row(305),295,305)['components'][0]
    assert result['qst0_kJ_mol']==pytest.approx(heat)
    assert result['input_interval_envelope_kJ_mol'][0]<heat<result['input_interval_envelope_kJ_mol'][1]


def test_stale_disk_job_is_interrupted_not_running(tmp_path):
    identifier='a'*32
    (tmp_path/identifier).mkdir()
    (tmp_path/identifier/'job-status.json').write_text(json.dumps({'job_id':identifier,'status':'running'}))
    manager=AdsorptionJobs(tmp_path,Backend())
    try:
        assert manager.status(identifier)['status']=='interrupted'
        with pytest.raises(KeyError):
            manager.status('../escape')
    finally:
        manager.close()


def test_cancel_during_terminal_thread_teardown_preserves_complete_result(tmp_path):
    # A completed result can be observed before the worker thread returns.
    from types import SimpleNamespace
    manager = AdsorptionJobs(tmp_path, Backend())
    identifier = 'b'*32
    manager._job = {'job_id': identifier, 'status': 'complete'}
    manager._thread = SimpleNamespace(is_alive=lambda: True, join=lambda timeout: None)
    try:
        assert manager.cancel(identifier)['status'] == 'complete'
        assert not manager._cancel.is_set()
    finally:
        manager.close()
