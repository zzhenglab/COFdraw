"""Bounded, cancellable local RASPA jobs; each stage retains its real inputs/output."""
from __future__ import annotations

import atexit
from copy import deepcopy
import csv
import hashlib
import json
import math
from pathlib import Path
import re
import subprocess
import threading
import time
import uuid
import zipfile
from datetime import datetime, timezone

from fastapi import APIRouter, HTTPException
from pydantic import BaseModel, ConfigDict, Field, model_validator
from typing import Annotated, Literal


class AdsorptionRequest(BaseModel):
    model_config = ConfigDict(extra='forbid', allow_inf_nan=False)
    structure: dict
    gases: list[Literal['CO2', 'N2', 'CH4', 'H2', 'Ar']] = Field(default=['CO2'], min_length=1, max_length=3)
    mole_fractions: list[Annotated[float, Field(ge=1e-8, le=1)]] = Field(default=[1.], min_length=1, max_length=3)
    temperature_K: float = Field(default=298.15, ge=20, le=600)
    pressures_Pa: list[Annotated[float, Field(ge=1, le=10_000_000)]] = Field(default=[10000., 100000.], min_length=1, max_length=8)
    cycles: int = Field(default=2000, ge=100, le=100000)
    init_cycles: int = Field(default=1000, ge=100, le=100000)
    widom_cycles: int = Field(default=2000, ge=100, le=100000)
    calculate_heat: bool = True
    seed: int = Field(default=42137, ge=1, le=2147483000)
    charge_mode: Literal['cif', 'zero'] = 'zero'
    fugacity_coefficients: list[Annotated[float, Field(ge=1e-6, le=10)]] | None = Field(default=None, min_length=1, max_length=3)

    @model_validator(mode='after')
    def consistent_mixture(self):
        if len(set(self.gases)) != len(self.gases):
            raise ValueError('Each gas must appear once')
        if len(self.gases) != len(self.mole_fractions) or abs(sum(self.mole_fractions)-1) > 1e-8:
            raise ValueError('Supply one positive mole fraction per gas, summing to one')
        if self.fugacity_coefficients is not None and len(self.fugacity_coefficients) != len(self.gases):
            raise ValueError('Supply one fugacity coefficient per gas')
        if len(set(self.pressures_Pa)) != len(self.pressures_Pa):
            raise ValueError('Pressure points must be distinct')
        if len(json.dumps(self.structure, allow_nan=False)) > 2_000_000:
            raise ValueError('Structure input limit is 2 MB')
        return self


def _write(path, value):
    temporary = path.with_suffix(path.suffix+'.tmp')
    temporary.write_text(json.dumps(value, indent=2, allow_nan=False), encoding='utf-8')
    temporary.replace(path)


def _now():
    return datetime.now(timezone.utc).isoformat()


def mixture_selectivity(components, values, fractions, metric):
    """Do not turn absent/zero component sampling into infinite selectivity."""
    result = []
    for i, a in enumerate(components):
        for j, b in enumerate(components[i+1:], i+1):
            x, y = values.get(a), values.get(b)
            valid = all(isinstance(v, (int, float)) and math.isfinite(v) and v > 0 for v in (x, y))
            value = ((x/y)/(fractions[i]/fractions[j]) if metric == 'loading' else x/y) if valid else None
            result.append({'numerator': a, 'denominator': b, 'selectivity': value,
                           'definition': '(q_i/q_j)/(y_i/y_j)' if metric == 'loading' else 'K_H,i / K_H,j',
                           'reason': None if valid else 'Both component estimates must be finite and positive; inspect sampling.'})
    return result


def zero_loading_heat(low, high, temperature_low, temperature_high):
    """Van't Hoff estimate at fixed framework: R d ln(K_H,p) / d(1/T)."""
    denominator = 1/temperature_low-1/temperature_high
    rows = []
    high_by_gas = {row['gas']: row for row in high.get('components', [])}
    for lower in low.get('components', []):
        upper = high_by_gas.get(lower['gas'], {})
        a, b = lower.get('henry_mol_kg_Pa'), upper.get('henry_mol_kg_Pa')
        valid = all(isinstance(v, (int, float)) and math.isfinite(v) and v > 0 for v in (a, b))
        heat = .00831446261815324*math.log(a/b)/denominator if valid else None
        bounds = None
        ea, eb = lower.get('henry_confidence95_mol_kg_Pa'), upper.get('henry_confidence95_mol_kg_Pa')
        if valid and all(isinstance(v, (int, float)) and math.isfinite(v) and v >= 0 for v in (ea, eb)) and a > ea and b > eb:
            bounds = [.00831446261815324*math.log((a-ea)/(b+eb))/denominator,
                      .00831446261815324*math.log((a+ea)/(b-eb))/denominator]
        rows.append({'gas': lower['gas'], 'qst0_kJ_mol': heat, 'input_interval_envelope_kJ_mol': bounds,
                     'reason': None if valid else 'Positive finite Henry estimates are required at both temperatures.'})
    return {'components': rows, 'temperature_low_K': temperature_low, 'temperature_high_K': temperature_high,
            'method': 'R d ln(K_H,p)/d(1/T), central finite difference of independent zero-loading Widom runs',
            'scope': 'Zero-loading isosteric heat estimate for this rigid classical model; not a finite-loading heat curve. Interval envelope uses the two input 95% intervals; it is not a joint 95% confidence interval and omits finite-difference bias.'}


class AdsorptionJobs:
    def __init__(self, output_root, backend=None, stage_timeout=3600):
        self.root = Path(output_root).resolve()
        self.backend = backend
        self.stage_timeout = stage_timeout
        self._lock = threading.RLock()
        self._job = None
        self._cancel = threading.Event()
        self._process = None
        self._thread = None
        self._closed = False
        atexit.register(self.close)

    def _engine(self):
        if self.backend is None:
            from . import cofdraw_raspa
            self.backend = cofdraw_raspa
        return self.backend

    def status(self, identifier):
        if not re.fullmatch('[0-9a-f]{32}', identifier):
            raise KeyError('Unknown adsorption job')
        with self._lock:
            if self._job and self._job['job_id'] == identifier:
                return deepcopy(self._job)
            path = self.root/identifier/'job-status.json'
            if not path.is_file() or path.stat().st_size > 30_000_000:
                raise KeyError('Unknown adsorption job')
            result = json.loads(path.read_text(encoding='utf-8'))
            if result.get('status') in {'queued', 'running', 'cancelling'}:
                result.update(status='interrupted', message='The application stopped before this job completed. Raw files remain in the job folder.')
            return result

    def _update(self, **values):
        with self._lock:
            self._job.update(values, updated_at_utc=_now())
            _write(self.root/self._job['job_id']/'job-status.json', self._job)

    def start(self, payload):
        request = AdsorptionRequest.model_validate(payload).model_dump()
        engine = self._engine()
        availability = engine.engine_status()
        if not availability.get('available'):
            raise ValueError(availability.get('reason', 'The bundled RASPA engine is unavailable'))
        with self._lock:
            if self._closed:
                raise RuntimeError('The application is shutting down')
            if self._thread and self._thread.is_alive():
                raise RuntimeError('An adsorption job is already running. Wait or cancel it first.')
            identifier = uuid.uuid4().hex
            directory = self.root/identifier
            directory.mkdir(parents=True, exist_ok=False)
            # Validate the complete model/force field before launching a background job.
            first = self._prepare(engine, directory/'widom', request, mode='widom', pressure=100000., index=0)
            _write(directory/'request.json', request)
            self._cancel.clear()
            self._job = {'job_id': identifier, 'report_id': identifier, 'operation': 'adsorption',
                         'status': 'queued', 'stage': 'Zero-loading Widom', 'completed_stages': 0,
                         'total_stages': 1+2*int(request['calculate_heat'])+len(request['pressures_Pa']), 'created_at_utc': _now(),
                         'input_sha256': hashlib.sha256(json.dumps(request, sort_keys=True).encode()).hexdigest(),
                         'engine': availability, 'message': 'Local simulation queued. Short runs are sampling checks, not converged predictions.'}
            self._update()
            self._thread = threading.Thread(target=self._run, args=(request, first), daemon=True, name='cofdraw-raspa')
            self._thread.start()
            return deepcopy(self._job)

    @staticmethod
    def _prepare(engine, directory, request, *, mode, pressure, index):
        directory.mkdir(parents=True, exist_ok=True)
        return engine.prepare_run(directory, request['structure'], gases=request['gases'],
            mole_fractions=request['mole_fractions'], temperature_K=request['temperature_K'], pressure_Pa=pressure,
            cycles=request['widom_cycles'] if mode == 'widom' else request['cycles'],
            init_cycles=request['init_cycles'], seed=request['seed']+index,
            charge_mode=request['charge_mode'], mode=mode, fugacity_coefficients=request['fugacity_coefficients'])

    def _execute(self, directory, manifest):
        if self._cancel.is_set():
            raise InterruptedError('Simulation cancelled')
        log_path = directory/'worker.log'
        with log_path.open('wb') as log:
            process = subprocess.Popen(self._engine().command_for(directory), cwd=directory, stdout=log,
                stderr=subprocess.STDOUT, stdin=subprocess.DEVNULL, creationflags=getattr(subprocess, 'CREATE_NO_WINDOW', 0))
            with self._lock:
                self._process = process
            deadline = time.monotonic()+self.stage_timeout
            try:
                while process.poll() is None:
                    if self._cancel.wait(.2):
                        raise InterruptedError('Simulation cancelled')
                    if time.monotonic() > deadline:
                        raise TimeoutError('Simulation exceeded the one-hour limit for one stage; retain inputs and use an external longer run.')
                    if log_path.stat().st_size > 50_000_000:
                        raise RuntimeError('Simulation log exceeded 50 MB; the run was stopped')
                if self._cancel.is_set():
                    raise InterruptedError('Simulation cancelled')
                if process.returncode:
                    raise RuntimeError(f'RASPA worker failed with exit code {process.returncode}. See worker.log in the saved job bundle.')
            finally:
                if process.poll() is None:
                    process.terminate()
                    try:
                        process.wait(timeout=5)
                    except subprocess.TimeoutExpired:
                        process.kill()
                        process.wait(timeout=5)
                with self._lock:
                    self._process = None
        result = self._engine().parse_run(directory, manifest)
        json.dumps(result, allow_nan=False)
        return result

    def _run(self, request, first_manifest):
        from . import __version__
        identifier = self._job['job_id']
        directory = self.root/identifier
        result = {'kind': 'raspa_adsorption', 'operation': 'adsorption', 'report_id': identifier,
                  'engine_version': __version__, 'report_schema_version': 'cofdraw-science-1',
                  'structure_sha256': first_manifest.get('structure_sha256'),
                  'inputs': request, 'input_sha256': self._job['input_sha256'],
                  'created_at_utc': self._job['created_at_utc'], 'widom': None, 'isotherm': [],
                  'scope': 'Rigid-framework classical simulations with the documented force field and charges. Numerical convergence and agreement with experiment require separate checks.',
                  'warnings': ['Short default runs support screening and workflow checks; they do not establish convergence.',
                               'Henry coefficients come from a separate zero-loading Widom calculation.',
                               'Selectivity ratios have no propagated confidence interval. Inspect component errors and rare insertion sampling.']}
        status, message = 'complete', 'Simulation finished. Review block errors, model assumptions and convergence before interpreting values.'
        completed = 0
        try:
            self._update(status='running', stage='Zero-loading Widom')
            result['widom'] = self._execute(directory/'widom', first_manifest)
            components = result['widom'].get('components', [])
            result['henry_selectivity'] = mixture_selectivity(request['gases'],
                {c['gas']: c.get('henry_mol_kg_Pa') for c in components}, request['mole_fractions'], 'henry')
            completed = 1
            result['zero_loading_heat'] = None
            if request['calculate_heat']:
                delta = min(5., .025*request['temperature_K'])
                neighbors = result['heat_widom_runs'] = []
                for offset, label in [(-delta, 'low'), (delta, 'high')]:
                    adjusted = {**request, 'temperature_K': request['temperature_K']+offset}
                    self._update(completed_stages=completed, stage=f'Zero-loading heat: {adjusted["temperature_K"]:g} K')
                    stage = directory/f'widom-{label}'
                    manifest = self._prepare(self._engine(), stage, adjusted, mode='widom', pressure=100000., index=completed)
                    neighbors.append(self._execute(stage, manifest))
                    completed += 1
                result['zero_loading_heat'] = zero_loading_heat(*neighbors, request['temperature_K']-delta, request['temperature_K']+delta)
            for index, pressure in enumerate(request['pressures_Pa'], 1):
                if self._cancel.is_set():
                    raise InterruptedError('Simulation cancelled')
                self._update(completed_stages=completed, stage=f'GCMC {pressure:g} Pa ({index}/{len(request["pressures_Pa"])})')
                stage = directory/f'pressure-{index:02d}'
                manifest = self._prepare(self._engine(), stage, request, mode='gcmc', pressure=pressure, index=completed)
                entry = self._execute(stage, manifest)
                entry['pressure_Pa'] = pressure
                entry['selectivity'] = mixture_selectivity(request['gases'],
                    {c['gas']: c.get('loading_mol_kg') for c in entry.get('components', [])}, request['mole_fractions'], 'loading')
                result['isotherm'].append(entry)
                completed += 1
                _write(directory/'partial-result.json', result)
        except InterruptedError as exc:
            status, message = 'cancelled', str(exc)
        except Exception as exc:
            status, message = 'failed', str(exc)
        try:
            result.update(status=status, message=message, finished_at_utc=_now())
            self._archive(directory, result)
            self._update(status=status, message=message, result=result,
                         completed_stages=completed,
                         downloads=result['downloads'])
        except Exception as exc:
            self._update(status='failed', message=f'{message} Report export failed: {exc}')

    @staticmethod
    def _archive(directory, result):
        identifier = directory.name
        result['downloads'] = {ext: f'/api/drawing/exports/{identifier}/{ext}' for ext in ('json', 'csv', 'zip')}
        _write(directory/f'{identifier}.json', result)
        with (directory/f'{identifier}.csv').open('w', newline='', encoding='utf-8') as stream:
            writer = csv.writer(stream)
            writer.writerow(['pressure_Pa', 'gas', 'loading_mol_kg', 'loading_confidence95_mol_kg'])
            for row in result['isotherm']:
                for gas in row.get('components', []):
                    writer.writerow([row['pressure_Pa'], gas['gas'], gas.get('loading_mol_kg'), gas.get('loading_confidence95_mol_kg')])
        files = [p for p in directory.rglob('*') if p.is_file() and not p.is_symlink() and p.suffix not in {'.zip', '.tmp'} and p.name != 'job-status.json']
        if sum(p.stat().st_size for p in files) > 100_000_000:
            raise ValueError('Run files exceed the 100 MB bundle limit; raw files remain in the local job folder')
        hashes = {p.relative_to(directory).as_posix(): hashlib.sha256(p.read_bytes()).hexdigest() for p in files}
        _write(directory/'sha256.json', hashes)
        with zipfile.ZipFile(directory/f'{identifier}.zip', 'w', zipfile.ZIP_DEFLATED) as archive:
            for path in [*files, directory/'sha256.json']:
                archive.write(path, path.relative_to(directory).as_posix())

    def cancel(self, identifier):
        with self._lock:
            if not self._job or self._job['job_id'] != identifier:
                raise KeyError('This job is not active')
            if self._job['status'] in {'queued', 'running', 'cancelling'} and self._thread and self._thread.is_alive():
                self._cancel.set()
                self._update(status='cancelling', message='Stopping the local worker and preserving completed stages…')
            return deepcopy(self._job)

    def close(self):
        with self._lock:
            self._closed = True
            self._cancel.set()
            process = self._process
            thread = self._thread
        if process and process.poll() is None:
            process.terminate()
        if thread and thread.is_alive() and thread is not threading.current_thread():
            thread.join(timeout=7)


def create_adsorption_router(output_root):
    router = APIRouter(prefix='/science/adsorption')
    jobs = AdsorptionJobs(output_root)
    router.add_event_handler('shutdown', jobs.close)

    @router.get('/engine')
    def engine():
        return jobs._engine().engine_status()

    @router.post('/jobs', status_code=202)
    def start(request: AdsorptionRequest):
        try:
            return jobs.start(request.model_dump())
        except RuntimeError as exc:
            raise HTTPException(409, str(exc)) from exc
        except (ValueError, KeyError, OSError) as exc:
            raise HTTPException(422, str(exc)) from exc

    @router.get('/jobs/{identifier}')
    def status(identifier: str):
        try:
            return jobs.status(identifier)
        except KeyError as exc:
            raise HTTPException(404, str(exc)) from exc

    @router.post('/jobs/{identifier}/cancel')
    def cancel(identifier: str):
        try:
            return jobs.cancel(identifier)
        except KeyError as exc:
            raise HTTPException(404, str(exc)) from exc

    return router
