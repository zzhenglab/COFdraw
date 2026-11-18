"""CPU-only native RASPA worker; execute as a child process, never in the API."""
from __future__ import annotations

import hashlib
import json
import math
import os
from pathlib import Path
import sys
import time

if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from cof_vision.cofdraw_raspa import ENGINE_VERSION, profile, resource_directory, _sha


def _write(path, value):
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(value, indent=2, allow_nan=False)+"\n", encoding="utf-8")
    temporary.replace(path)


def _cell_parameters(cell):
    import numpy as np
    matrix = np.asarray(cell, dtype=float)
    lengths = np.linalg.norm(matrix, axis=1)
    # Native 3.1.0 takes radians, despite the wheel's .pyi prose saying degrees.
    angles = [math.acos(float(np.clip(np.dot(matrix[i], matrix[j]) / (lengths[i] * lengths[j]), -1, 1)))
              for i,j in [(1,2),(0,2),(0,1)]]
    return [*lengths.tolist(), *angles]


def run(directory):
    """Construct the declared force field, run RASPA, and save actual statistics."""
    job = Path(directory).resolve()
    input_path = job / "engine-input.json"
    if not input_path.is_file() or input_path.stat().st_size > 5_000_000:
        raise ValueError("Missing or oversized prepared RASPA input")
    payload = json.loads(input_path.read_text("utf-8"))
    manifest, structure = payload["manifest"], payload["structure"]
    data = profile()
    if manifest["profile_sha256"] != _sha(resource_directory()/"profile.json"):
        raise ValueError("Prepared force-field profile changed")
    os.environ["OMP_NUM_THREADS"] = "1"
    os.environ["RASPALIB_QUIET"] = "1"  # CPU simulations do not require OpenCL.
    os.environ.pop("RASPALIB_FORCE", None)  # Use the library's checked CPU dispatch.
    import raspalib as raspa
    from rdkit import Chem
    if raspa.__version__ != ENGINE_VERSION:
        raise ValueError("Unexpected RASPA runtime version")
    pt = Chem.GetPeriodicTable()
    pseudo, parameters, index = [], [], {}
    for element in sorted({a["element"] for a in structure["atoms"]}):
        atom = data["framework"][element]
        index["framework:"+element] = len(pseudo)
        pseudo.append(raspa.PseudoAtom(name="framework_"+element, framework_type=True,
                                      mass=pt.GetAtomicWeight(element), charge=0., atomic_number=pt.GetAtomicNumber(element), source=atom["source"]))
        parameters.append(raspa.VDWParameters(atom["epsilon_K"], atom["sigma_angstrom"]))
    for gas in manifest["gases"]:
        for atom in data["gases"][gas]["sites"]:
            name = gas+":"+atom["name"]
            if name in index:
                continue
            index[name] = len(pseudo)
            pseudo.append(raspa.PseudoAtom(name=gas+"_"+atom["name"], framework_type=False, mass=atom["mass"], charge=atom["charge"],
                                          atomic_number=pt.GetAtomicNumber(atom["element"]), print_to_pdb=atom.get("print_to_output", True), source=atom["source"]))
            parameters.append(raspa.VDWParameters(atom["epsilon_K"], atom["sigma_angstrom"]))
    ff = raspa.ForceField(pseudo_atoms=pseudo, parameters=parameters,
                         cutoff_framework_vdw=data["cutoff_angstrom"], cutoff_molecule_vdw=data["cutoff_angstrom"],
                         cutoff_coulomb=data["cutoff_angstrom"], shifted=data["shifted"], tail_corrections=data["tail_corrections"],
                         use_charge=any(a["charge"] for a in structure["atoms"]) or any(a["charge"] for gas in manifest["gases"] for a in data["gases"][gas]["sites"]))
    components = []
    for k, gas in enumerate(manifest["gases"]):
        definition = data["gases"][gas]
        atoms = [raspa.Atom(position=tuple(atom["position_angstrom"]), charge=atom["charge"], type=index[gas+":"+atom["name"]]) for atom in definition["sites"]]
        if manifest["mode"] == "widom":
            moves = raspa.MCMoveProbabilities(widom_probability=1.)
        else:
            moves = raspa.MCMoveProbabilities(translation_probability=.5, rotation_probability=.5 if len(atoms)>1 else 0.,
                                             swap_cbmc_probability=1.)
        components.append(raspa.Component(force_field=ff, component_name=gas,
                                          critical_temperature=definition["critical_temperature_K"], critical_pressure=definition["critical_pressure_Pa"],
                                          acentric_factor=definition["acentric_factor"], defined_atoms=atoms, particle_probabilities=moves,
                                          fugacity_coefficient=manifest["internal_fugacity_coefficients"][k]))
    parameters = _cell_parameters(structure["cell"])
    framework = raspa.Framework(force_field=ff, component_name="COFdraw_explicit_framework",
                                simulation_box=raspa.SimulationBox(*parameters), space_group_hall_number=1,
                                defined_atoms=[raspa.Atom(position=tuple(a["fractional"]), charge=a["charge"], type=index["framework:"+a["element"]], is_fractional=True) for a in structure["atoms"]],
                                number_of_unit_cells=manifest["repeats"])
    system = raspa.System(force_field=ff, framework_components=framework, components=components,
                         external_temperature=manifest["temperature_K"], external_pressure=manifest["pressure_Pa"], initial_number_of_molecules=[0]*len(components))
    # Disable upstream formatted final output: its mixture enthalpy formatter can
    # terminate on a singular covariance. Read supported native statistics only.
    simulation = raspa.MonteCarlo(number_of_cycles=manifest["cycles"], number_of_initialization_cycles=manifest["init_cycles"],
                                  print_every=max(100,manifest["cycles"]), systems=[system], random_seed=manifest["seed"], number_of_blocks=5,
                                  output_to_files=False)
    os.chdir(job)
    _write(job/"engine-progress.json", {"phase":"running", "engine_version":raspa.__version__, "mode":manifest["mode"]})
    started = time.monotonic()
    simulation.run()
    actual = simulation.systems[0]
    mass = float(actual.framework_mass())
    if not math.isclose(mass, manifest["framework_mass_g_mol"], rel_tol=1e-9):
        raise RuntimeError("Native framework mass differs from prepared composition")
    volume = float(actual.simulation_box.volume)
    expected_volume = structure["cell_volume_angstrom3"] * math.prod(manifest["repeats"])
    if not math.isclose(volume, expected_volume, rel_tol=1e-9):
        raise RuntimeError("Native simulation cell volume differs from the supplied periodic structure")
    result = {"status":"completed", "engine":"RASPA3", "engine_version":raspa.__version__,
              "cpu_variant":getattr(raspa,"__cpu_level__","unknown"), "mode":manifest["mode"], "gases":manifest["gases"],
              "seed":manifest["seed"], "input_sha256":_sha(input_path), "components":[], "elapsed_seconds":time.monotonic()-started,
              "simulation_volume_angstrom3":volume, "framework_mass_g_mol":mass, "number_of_blocks":5,
              "uncertainty":"RASPA five-block 95% confidence interval; not a convergence guarantee", "warnings":[],
              "enthalpy_kJ_mol":None,"enthalpy_reason":manifest["enthalpy_status"]}
    for k,component in enumerate(actual.components):
        row = {"gas":manifest["gases"][k],"mole_fraction":manifest["mole_fractions"][k],
               "fugacity_coefficient":manifest["fugacity_coefficients"][k],"partial_fugacity_Pa":manifest["partial_fugacities_Pa"][k]}
        if manifest["mode"] == "widom":
            mean, confidence = component.average_rosenbluth_weights.result()
            conversion = volume*1e-27 / (1.380649e-23 * manifest["temperature_K"] * mass)
            row.update(henry_mol_kg_Pa=float(mean)*conversion,henry_confidence95_mol_kg_Pa=float(confidence)*conversion,
                       rosenbluth_weight=float(mean),rosenbluth_confidence95=float(confidence))
        else:
            mean, confidence = actual.average_loadings.average_loading_number_of_molecules(k)
            row.update(loading_mol_kg=float(mean)*1000/mass,loading_confidence95_mol_kg=float(confidence)*1000/mass,
                       molecules_per_simulation_cell=float(mean),molecules_confidence95=float(confidence))
        result["components"].append(row)
    if manifest["mode"] == "gcmc":
        # This accessor emits actual per-block loadings without evaluating the
        # separate unsafe mixture-enthalpy statistic.
        (job/"loading-statistics.txt").write_text(actual.average_loadings.write_averages_statistics(actual.components,mass,manifest["repeats"]),encoding="utf-8")
    _write(job/"engine-result.json",result)
    _write(job/"engine-progress.json",{"phase":"completed","elapsed_seconds":result["elapsed_seconds"]})
    return result


def main(argv=None):
    args = sys.argv[1:] if argv is None else argv
    if len(args)!=1:
        raise SystemExit("Usage: --raspa-worker JOB_DIRECTORY")
    try:
        run(args[0])
    except Exception as exc:
        job=Path(args[0]).resolve()
        if job.is_dir():
            _write(job/"engine-error.json",{"status":"failed","error":str(exc),"type":type(exc).__name__})
        print(str(exc),file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
