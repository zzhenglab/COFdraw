"""Validated inputs and results for an isolated, pinned RASPA3 CPU worker.

This module never imports the native simulation library into the API process.
Cell vectors are rows; framework atoms carry fractional coordinates and charges.
"""
from __future__ import annotations

import hashlib
import importlib.metadata
import importlib.util
import json
import math
from pathlib import Path
import platform
import sys

import numpy as np
from rdkit import Chem

from .cofdraw_science import normalize_structure, _finite, _integer
from .cofdraw_network import inspect_periodic_network

ENGINE_VERSION = "3.1.0"
MAX_SUPERCELL_ATOMS = 12000
MAX_WORK = 100_000_000_000


def resource_directory():
    if getattr(sys, "frozen", False):
        return Path(sys.executable).resolve().parent / "data" / "raspa"
    return Path(__file__).resolve().parents[2] / "data" / "raspa"


def _sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def _json_hash(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False).encode()).hexdigest()


def profile():
    root = resource_directory()
    manifest = json.loads((root / "engine-manifest.json").read_text("utf-8"))
    if _sha(root / "profile.json") != manifest["profile_sha256"]:
        raise RuntimeError("The bundled RASPA force-field profile failed its checksum")
    return json.loads((root / "profile.json").read_text("utf-8"))


def engine_status():
    """Report package/resource availability without loading a native extension."""
    result = {"engine": "RASPA3", "distribution": "raspalib", "required_version": ENGINE_VERSION,
              "available": False, "execution": "isolated CPU process", "gases": ["CO2", "N2", "CH4", "H2", "Ar"],
              "modes": ["gcmc", "widom"], "requires_gpu": False}
    try:
        result["version"] = importlib.metadata.version("raspalib")
        spec = importlib.util.find_spec("raspalib")
        if result["version"] != ENGINE_VERSION or not spec or not spec.origin:
            raise RuntimeError("Install the pinned raspalib==3.1.0 distribution")
        package = Path(spec.origin).parent
        variants = list(package.glob("raspalib_*.*"))
        if not any(p.suffix.lower() in (".pyd", ".so") for p in variants):
            raise RuntimeError("RASPA native CPU modules are missing")
        if sys.platform == "win32" and sys.version_info[:2] == (3, 11):
            manifest = json.loads((resource_directory()/"engine-manifest.json").read_text("utf-8"))
            for filename, expected in manifest["windows_cp311_native_files"].items():
                if not (package/filename).is_file() or _sha(package/filename) != expected:
                    raise RuntimeError("RASPA native module failed its official-wheel checksum: " + filename)
        data = profile()
        result.update(available=True, profile_id=data["id"], profile_name=data["name"],
                      profile_sha256=_sha(resource_directory() / "profile.json"), platform=platform.system())
    except (ImportError, importlib.metadata.PackageNotFoundError, OSError, ValueError, RuntimeError) as exc:
        result["reason"] = str(exc)
    return result


def prepare_run(directory, structure, gases, mole_fractions, temperature_K, pressure_Pa,
                cycles, init_cycles, seed, charge_mode, mode="gcmc", fugacity_coefficients=None):
    """Write a bounded immutable job input and return its reviewable manifest.

    phi defaults to one (an ideal-gas reservoir), never an implicit mixture EOS.
    The framework uses explicit CIF charges or an explicitly requested zero model.
    """
    status = engine_status()
    if not status["available"]:
        raise ValueError("RASPA3 is unavailable: " + status.get("reason", "missing runtime"))
    data = profile()
    if not isinstance(gases, list) or not 1 <= len(gases) <= 5 or any(not isinstance(g, str) or g not in data["gases"] for g in gases) or len(set(gases)) != len(gases):
        raise ValueError("Choose 1–5 distinct supported gases: CO2, N2, CH4, H2, Ar")
    if not isinstance(mole_fractions, (list, tuple)) or len(mole_fractions) != len(gases):
        raise ValueError("Provide one mole fraction for every gas")
    fractions = [_finite(y, "mole fraction", 1e-8, 1) for y in mole_fractions]
    if abs(sum(fractions) - 1) > 1e-6:
        raise ValueError("Physical gas mole fractions must sum to one")
    fractions = [y / sum(fractions) for y in fractions]
    temperature = _finite(temperature_K, "temperature_K", 19, 1000)
    pressure = _finite(pressure_Pa, "pressure_Pa", 1, 1e8)
    cycles = _integer(cycles, "cycles", 100, 2_000_000)
    initial = _integer(init_cycles, "init_cycles", 0, 1_000_000)
    seed = _integer(seed, "seed", 0, 2**31 - 1)
    if mode not in ("gcmc", "widom"):
        raise ValueError("mode must be gcmc or widom")
    if charge_mode not in ("zero", "cif"):
        raise ValueError("Select explicit CIF charges or the explicit zero-framework-charge approximation")
    if fugacity_coefficients is None:
        coefficients = [1.] * len(gases)
    else:
        if not isinstance(fugacity_coefficients, (list, tuple)) or len(fugacity_coefficients) != len(gases):
            raise ValueError("Supply one positive fugacity coefficient per gas")
        coefficients = [_finite(v, "fugacity coefficient", 1e-6, 100) for v in fugacity_coefficients]
    structure = normalize_structure(structure)
    missing = sorted({a["element"] for a in structure["atoms"]} - set(data["framework"]))
    if missing:
        raise ValueError("No reviewed framework force-field parameters for: " + ", ".join(missing))
    warnings = [data["scope"], "Five-block confidence intervals describe sampling uncertainty; they do not establish equilibration, convergence or force-field accuracy."]
    if charge_mode == "cif" and any("charge" not in a for a in structure["atoms"]):
        raise ValueError("CIF charge mode requires an explicit charge for every framework atom")
    if charge_mode == "zero":
        for a in structure["atoms"]:
            a["charge"] = 0.
        warnings.append("Framework charges are explicitly zero; CO2/N2 guest charges remain active. Electrostatic host–guest selectivity may be inaccurate.")
    net_charge = sum(a["charge"] for a in structure["atoms"])
    if abs(net_charge) > 1e-5:
        raise ValueError("Ewald adsorption requires a neutral framework cell; resolve charges/counterions (net charge %.8g e)" % net_charge)
    contacts = inspect_periodic_network(structure)
    if contacts["short_contact_count"]:
        raise ValueError("Framework has severe short contacts; resolve the periodic model before adsorption simulation")
    cell = np.asarray(structure["cell"])
    heights = 1 / np.linalg.norm(np.linalg.inv(cell), axis=0)
    repeats = np.maximum(1, np.ceil(2 * data["cutoff_angstrom"] / heights - 1e-10).astype(int))
    count = int(len(structure["atoms"]) * np.prod(repeats))
    if count > MAX_SUPERCELL_ATOMS or max(repeats) > 64:
        raise ValueError("The required >=24 Å cell heights exceed the bounded 12000-atom simulation limit")
    if count * (cycles + initial) * len(gases) > MAX_WORK:
        raise ValueError("Requested simulation exceeds the bounded atom-cycle work limit; reduce cycles or cell size")
    if "H2" in gases:
        warnings.append(data["gases"]["H2"]["qualification"])
    if fugacity_coefficients is None:
        warnings.append("Reservoir fugacity coefficients are explicitly one (ideal gas); supply validated coefficients for nonideal/high-pressure reservoirs.")
    if len(gases) > 1:
        warnings.append("Physical fractions and fugacities are preserved here. Internal equal RASPA fractions are compensated by N*y*phi; upstream excess loadings and mixture EOS are not reported.")
    pt = Chem.GetPeriodicTable()
    mass = sum(pt.GetAtomicWeight(a["element"]) for a in structure["atoms"]) * int(np.prod(repeats))
    manifest = {"schema": 1, "engine": "RASPA3", "engine_version": ENGINE_VERSION,
                "profile_id": data["id"], "profile_sha256": status["profile_sha256"],
                "mode": mode, "gases": gases, "mole_fractions": fractions,
                "fugacity_coefficients": coefficients, "partial_fugacities_Pa": [pressure*y*phi for y,phi in zip(fractions,coefficients)],
                "internal_fugacity_coefficients": [len(gases)*y*phi for y,phi in zip(fractions,coefficients)],
                "temperature_K": temperature, "pressure_Pa": pressure, "cycles": cycles,
                "init_cycles": initial if mode == "gcmc" else 0, "seed": seed,
                "charge_mode": charge_mode, "framework_net_charge_e": net_charge,
                "repeats": repeats.tolist(), "cell_heights_angstrom": heights.tolist(),
                "supercell_heights_angstrom": (heights*repeats).tolist(),
                "simulation_atom_count": count, "framework_mass_g_mol": mass,
                "structure_sha256": _json_hash(structure), "warnings": warnings,
                "enthalpy_status": "Finite-loading enthalpy is not reported by this adapter. Zero-loading Qst requires independent temperature-dependent Widom calculations.",
                "scope": "Model adsorption from an explicit periodic rigid framework; not an experimental isotherm or BET measurement."}
    resources = [("force-field-profile.json", (resource_directory()/"profile.json").read_bytes()),
                 *[(name, (resource_directory()/name).read_bytes()) for name in ("RASPA2-LICENSE.txt", "RASPA3-LICENSE.txt")]]
    if hashlib.sha256(resources[0][1]).hexdigest() != manifest["profile_sha256"]:
        raise RuntimeError("The RASPA force-field profile changed during job preparation")
    manifest["resource_files"] = [{"path": name, "sha256": hashlib.sha256(content).hexdigest()} for name, content in resources]
    payload = {"manifest": manifest, "structure": structure}
    job = Path(directory).resolve()
    if job.exists() and any(job.iterdir()):
        raise ValueError("Use a new empty job directory")
    job.mkdir(parents=True, exist_ok=True)
    for (name, content), record in zip(resources, manifest["resource_files"]):
        (job/name).write_bytes(content)
        if _sha(job/name) != record["sha256"]:
            raise RuntimeError("A saved RASPA resource failed its checksum: " + name)
    (job / "engine-input.json").write_text(json.dumps(payload, indent=2, allow_nan=False)+"\n", encoding="utf-8")
    manifest["input_sha256"] = _sha(job / "engine-input.json")
    (job / "manifest.json").write_text(json.dumps(manifest, indent=2, allow_nan=False)+"\n", encoding="utf-8")
    return manifest


def command_for(directory):
    job = str(Path(directory).resolve())
    if getattr(sys, "frozen", False):
        return [sys.executable, "--raspa-worker", job]
    return [sys.executable, str(Path(__file__).with_name("cofdraw_raspa_worker.py")), job]


def parse_run(directory, manifest):
    """Accept only a completed worker result tied to the actual input/version."""
    job = Path(directory).resolve()
    if _sha(job / "engine-input.json") != manifest["input_sha256"]:
        raise ValueError("RASPA job input changed after preparation")
    result_path = job / "engine-result.json"
    if not result_path.is_file():
        raise ValueError("RASPA did not produce a completed result")
    result = json.loads(result_path.read_text("utf-8"))
    if result.get("engine_version") != ENGINE_VERSION or result.get("input_sha256") != manifest["input_sha256"] or result.get("status") != "completed":
        raise ValueError("RASPA result is incomplete or does not match this input")
    if result.get("gases") != manifest["gases"] or result.get("mode") != manifest["mode"]:
        raise ValueError("RASPA result gas/mode mismatch")
    if not isinstance(result.get("components"), list) or [r.get("gas") for r in result["components"]] != manifest["gases"]:
        raise ValueError("RASPA result components do not match this input")
    for row in result["components"]:
        metric = "henry_mol_kg_Pa" if manifest["mode"] == "widom" else "loading_mol_kg"
        confidence = "henry_confidence95_mol_kg_Pa" if manifest["mode"] == "widom" else "loading_confidence95_mol_kg"
        for field in (metric, confidence):
            value = row.get(field)
            if not isinstance(value, (float, int)) or not math.isfinite(value) or value < 0:
                raise ValueError("RASPA produced a nonfinite or negative " + field)
    result["manifest"] = manifest
    result["warnings"] = list(manifest["warnings"]) + result.get("warnings", [])
    result["result_sha256"] = _sha(result_path)
    result["output_files"] = [{"path": str(p.relative_to(job)).replace("\\", "/"), "sha256": _sha(p)}
                              for p in sorted(job.rglob("*")) if p.is_file() and p.stat().st_size <= 20_000_000]
    return result
