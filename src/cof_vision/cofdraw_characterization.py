"""Composition and kinematic powder diffraction for explicit periodic models.

Cell vectors are rows in angstrom. These calculations do not identify a phase,
relax a model, estimate an experimental BET area, or assign accessible active sites.
"""
from __future__ import annotations

from collections import Counter
from functools import reduce
import hashlib
import json
import math

import numpy as np
from rdkit import rdBase

from .cofdraw_science import AVOGADRO, PT, site_mass, _cell_parameters, _finite, normalize_structure
from .cofdraw_xray_factors import DATA_SHA256, IT92, SOURCE_VERSION

MAX_HKL_CANDIDATES = 600_000
MAX_REFLECTIONS = 60_000
MAX_ATOM_REFLECTIONS = 24_000_000
MAX_PROFILE_POINTS = 20_001
MERGE_TOLERANCE_DEGREES = 1e-6


def _structure_hash(structure):
    return hashlib.sha256(json.dumps(structure, sort_keys=True, allow_nan=False).encode()).hexdigest()


def _hill_order(counts):
    if "C" in counts:
        return [e for e in ("C", "H") if e in counts] + sorted(set(counts) - {"C", "H"})
    return sorted(counts)


def _formula(counts, divisor=1):
    return "".join(e + (str(counts[e] // divisor) if counts[e] // divisor != 1 else "") for e in _hill_order(counts))


def structure_summary(structure):
    """Count the supplied full-cell sites, including all explicit guests/H atoms.

    Element-site densities are atom counts per mass of the supplied model,
    not accessible binding-site concentrations. No missing hydrogens are added.
    """
    normalized = normalize_structure(structure)
    counts = Counter(atom["element"] for atom in normalized["atoms"])
    masses = Counter()
    for atom in normalized["atoms"]:
        masses[atom["element"]] += site_mass(atom)
    mass = sum(masses.values())
    count = sum(counts.values())
    volume = normalized["cell_volume_angstrom3"]
    lengths_angles = _cell_parameters(np.asarray(normalized["cell"]))
    warnings = [
        "Formula, mass, density, and element-site concentrations count only the supplied full-cell atoms. Missing hydrogens, guests, counterions, and disorder are not reconstructed.",
        "Element-site density is a stoichiometric atom inventory, not an accessible adsorption or catalytic site count.",
        "Model geometry and composition do not establish experimental phase identity, stability, activation state, or synthesis feasibility.",
    ]
    if "H" not in counts:
        warnings.append("No explicit hydrogen sites are present. If hydrogens are missing, mass, density, and mass-normalized quantities are incomplete.")
    from .cofdraw_network import inspect_periodic_network
    try:
        network = inspect_periodic_network(normalized)
        contacts = {"available": True, **{key: network[key] for key in (
            "minimum_candidate_distance_angstrom", "short_contact_count", "short_contacts",
            "isolated_atom_indices", "borderline_contact_count")},
            "method": {key: network["method"][key] for key in ("short_contact_rule", "radii", "periodic_search", "bond_scale", "tolerance_angstrom")},
            "warnings": network["warnings"]}
        if network["short_contact_count"]:
            warnings.append("Implausibly short periodic contacts are present; inspect the geometry before interpreting model-derived properties.")
    except ValueError as exc:
        # Composition remains exact even when an expensive contact search exceeds
        # its independent limits. Never substitute an incomplete neighbor graph.
        contacts = {"available": False, "reason": str(exc)}
        warnings.append("Periodic contact checks are unavailable: " + str(exc))
    return {
        "kind": "periodic_structure_summary", "quality": "model_geometry",
        "formula": _formula(counts), "reduced_formula": _formula(counts, reduce(math.gcd, counts.values())),
        "isotope_site_counts": dict(Counter(f"{a['isotope']}{a['element']}" for a in normalized["atoms"] if a.get("isotope"))),
        "cell_molar_mass_g_mol": float(mass), "cell_volume_angstrom3": volume,
        "density_g_cm3": float(mass / (AVOGADRO * volume * 1e-24)), "atom_count": count,
        "elements": [{"element": e, "count": counts[e], "atomic_fraction": counts[e] / count,
                      "mass_fraction": masses[e] / mass,
                      "site_density_mmol_g": 1000 * counts[e] / mass} for e in _hill_order(counts)],
        "cell_parameters": dict(zip(("a", "b", "c", "alpha", "beta", "gamma"), map(float, lengths_angles))),
        "units": {"cell_lengths": "angstrom", "cell_angles": "degree", "site_density": "mmol_of_element_atoms_per_g_of_supplied_structure"},
        "contact_checks": contacts, "structure_sha256": _structure_hash(normalized),
        "method": {"atomic_weights": "RDKit isotope masses for labelled sites and standard atomic weights elsewhere",
                   "rdkit_version": rdBase.rdkitVersion, "occupancy": "fully occupied resolved sites only"},
        "warnings": warnings,
    }


def atomic_form_factor(element, stol_squared):
    """Neutral-atom IT92 f0, in electrons; argument is (sin(theta)/lambda)^2.

    Supports H through Cf and 0 <= sin(theta)/lambda <= 2 inverse angstrom.
    No anomalous dispersion, ion charge, or coefficient renormalization.
    """
    if element not in IT92:
        raise ValueError(f"No bundled neutral-atom IT92 form factor for {element}")
    values = np.asarray(stol_squared, dtype=float)
    if not np.isfinite(values).all() or np.any(values < 0) or np.any(values > 4 + 1e-12):
        raise ValueError("IT92 form factors require 0 <= sin(theta)/wavelength <= 2 inverse angstrom")
    coeff = np.asarray(IT92[element])
    return np.sum(coeff[:4] * np.exp(-values[..., None] * coeff[4:8]), axis=-1) + coeff[8]


def _reflections(cell, qmax, atom_count):
    # h_i = q dot a_i gives a rigorous enumeration bound for skewed cells too.
    bounds = np.floor(qmax * np.linalg.norm(cell, axis=1) + 1e-9).astype(int)
    candidates = math.prod(int(2 * n + 1) for n in bounds)
    if candidates > MAX_HKL_CANDIDATES:
        raise ValueError("Powder XRD reciprocal search exceeds the bounded workload; use a primitive/smaller cell, lower maximum angle, or longer wavelength")
    axes = [np.arange(-n, n + 1, dtype=np.int32) for n in bounds]
    hkls = np.stack(np.meshgrid(*axes, indexing="ij"), axis=-1).reshape(-1, 3)
    reciprocal = np.linalg.inv(cell).T
    lengths = np.linalg.norm(hkls @ reciprocal, axis=1)
    mask = (lengths > 1e-12) & (lengths <= qmax + 1e-12)
    hkls, lengths = hkls[mask], lengths[mask]
    if len(hkls) > MAX_REFLECTIONS or len(hkls) * atom_count > MAX_ATOM_REFLECTIONS:
        raise ValueError("Powder XRD structure-factor calculation exceeds the bounded workload; use a primitive/smaller cell or lower maximum angle")
    return hkls, lengths, candidates


def powder_xrd(structure, *, wavelength=1.5406, two_theta_min=2, two_theta_max=50,
               step=.02, fwhm=.15, b_iso=0.0):
    """Deterministic kinematic random-powder pattern, single wavelength in angstrom.

    Enumerate all reciprocal vectors, sum atomic complex amplitudes, apply the
    unpolarized Bragg-Brentano LP factor and a common isotropic Debye-Waller B.
    Gaussian FWHM is a chosen display width, not a measured crystallite size.
    """
    normalized = normalize_structure(structure)
    wave = _finite(wavelength, "wavelength", .1, 5)
    lower = _finite(two_theta_min, "two_theta_min", .1, 169.9)
    upper = _finite(two_theta_max, "two_theta_max", .2, 170)
    if upper <= lower:
        raise ValueError("two_theta_max must be greater than two_theta_min")
    spacing = _finite(step, "step", .005, 1)
    width = _finite(fwhm, "fwhm", .02, 3)
    bfactor = _finite(b_iso, "b_iso", 0, 20)
    if spacing > width / 2:
        raise ValueError("Profile step must be at most half the chosen FWHM to resolve the peak shape")
    grid_count = math.ceil((upper - lower) / spacing) + 1
    if grid_count > MAX_PROFILE_POINTS:
        raise ValueError("Powder XRD profile exceeds 20001 points; increase the step or narrow the angle range")
    elements = [atom["element"] for atom in normalized["atoms"]]
    missing = sorted(set(elements) - set(IT92))
    if missing:
        raise ValueError("No bundled neutral-atom IT92 form factors for " + ", ".join(missing))
    cell = np.asarray(normalized["cell"])
    fraction = np.asarray([atom["fractional"] for atom in normalized["atoms"]])
    sigma = width / math.sqrt(8 * math.log(2))
    # Include neighboring peaks for continuous Gaussian tails at display edges.
    search_low, search_high = max(.001, lower - 5 * sigma), min(179, upper + 5 * sigma)
    qmax = 2 * math.sin(math.radians(search_high / 2)) / wave
    if qmax / 2 > 2:
        raise ValueError("Requested diffraction range exceeds the IT92 sin(theta)/wavelength <= 2 inverse angstrom limit")
    hkls, q, candidates = _reflections(cell, qmax, len(elements))
    angles = np.degrees(2 * np.arcsin(np.clip(wave * q / 2, -1, 1)))
    keep = angles >= search_low - 1e-9
    hkls, q, angles = hkls[keep], q[keep], angles[keep]
    groups = {e: np.flatnonzero(np.asarray(elements) == e) for e in sorted(set(elements))}
    squared = np.zeros(len(q))
    for start in range(0, len(q), 256):
        stop = start + 256
        stol2 = (q[start:stop] / 2) ** 2
        phases = np.exp(2j * np.pi * (hkls[start:stop] @ fraction.T))
        amplitude = np.zeros(len(stol2), dtype=complex)
        for element, indices in groups.items():
            amplitude += atomic_form_factor(element, stol2) * phases[:, indices].sum(axis=1)
        amplitude *= np.exp(-bfactor * stol2)
        squared[start:stop] = np.abs(amplitude) ** 2
    # A relative numerical extinction threshold refers to the unattenuated
    # neutral electron inventory, not the largest LP-enhanced low-angle peak.
    electrons = sum(PT.GetAtomicNumber(e) for e in elements)
    keep = squared > electrons ** 2 * 1e-14
    hkls, q, angles, squared = hkls[keep], q[keep], angles[keep], squared[keep]
    theta = np.radians(angles / 2)
    lp = (1 + np.cos(2 * theta) ** 2) / (np.sin(theta) ** 2 * np.cos(theta))
    intensities = squared * lp
    peaks = []
    for index in np.argsort(angles, kind="stable"):
        angle = float(angles[index])
        hkl = [int(x) for x in hkls[index]]
        if peaks and angle - peaks[-1]["two_theta_degrees"] <= MERGE_TOLERANCE_DEGREES:
            peaks[-1]["integrated_intensity_arbitrary"] += float(intensities[index])
            peaks[-1]["hkls"].append(hkl)
            peaks[-1]["multiplicity"] += 1
        else:
            peaks.append({"two_theta_degrees": angle, "d_spacing_angstrom": float(1 / q[index]),
                          "integrated_intensity_arbitrary": float(intensities[index]),
                          "hkls": [hkl], "multiplicity": 1})
    visible = [p for p in peaks if lower - 1e-9 <= p["two_theta_degrees"] <= upper + 1e-9]
    maximum = max((p["integrated_intensity_arbitrary"] for p in visible), default=0)
    for peak in visible:
        peak["intensity"] = 100 * peak["integrated_intensity_arbitrary"] / maximum
    grid = np.linspace(lower, upper, grid_count)
    profile = np.zeros_like(grid)
    for peak in peaks:
        angle = peak["two_theta_degrees"]
        first, last = np.searchsorted(grid, [angle - 5 * sigma, angle + 5 * sigma])
        profile[first:last] += peak["integrated_intensity_arbitrary"] * np.exp(-.5 * ((grid[first:last] - angle) / sigma) ** 2) / (sigma * math.sqrt(2 * math.pi))
    if profile.size and profile.max() > 0:
        profile *= 100 / profile.max()
    warnings = [
        "Simulated kinematic diffraction from the supplied coordinates; this is not an experimental pattern or a phase identification.",
        "Neutral isolated-atom form factors; anomalous dispersion, ionic bonding corrections, absorption, preferred orientation, and multiple scattering are omitted.",
        "The default 1.5406 angstrom is a nominal single Cu K-alpha1 wavelength. No K-alpha doublet, instrument calibration, or spectral distribution is modeled.",
        "The chosen Gaussian width is constant in 2-theta. It does not estimate crystallite size, strain, stacking disorder, or instrument resolution.",
        "All supplied atoms scatter, including guests. Only fully occupied resolved sites and a common user-supplied isotropic B factor are supported; per-site CIF displacement parameters are not retained.",
        "Relative peak integrated intensities and the sampled profile are independently normalized to 100; neither is an absolute intensity or an experimental phase fraction.",
    ]
    if not visible:
        warnings.append("No non-extinct Bragg peak centers occur in the displayed range. Any nonzero profile consists only of neighboring peak tails.")
    csv = "two_theta_degrees,relative_intensity\n" + "".join(f"{angle:.8f},{value:.8g}\n" for angle, value in zip(grid, profile))
    return {"kind": "simulated_powder_xrd", "quality": "model_simulation", "structure_sha256": _structure_hash(normalized),
            "parameters": {"wavelength_angstrom": wave, "two_theta_min": lower, "two_theta_max": upper,
                           "step_degrees": spacing, "actual_step_degrees": float(grid[1] - grid[0]),
                           "fwhm_degrees": width, "b_iso_angstrom2": bfactor},
            "profile": {"two_theta_degrees": grid.tolist(), "intensity": profile.tolist()},
            "peaks": visible, "peak_count": len(visible),
            "method": {"scattering": SOURCE_VERSION, "coefficient_sha256": DATA_SHA256,
                       "form_factor": "f0(s) = sum_i a_i exp(-b_i*s^2) + c; s=sin(theta)/wavelength",
                       "structure_factor": "F(hkl) = sum_j f0_j(s) exp(2*pi*i*hkl.fractional_j) exp(-B*s^2)",
                       "lorentz_polarization": "(1 + cos(2*theta)^2) / (sin(theta)^2 * cos(theta)); unpolarized Bragg-Brentano",
                       "multiplicity": "All reciprocal vectors including Friedel pairs enumerated once; equal-angle intensities summed, no additional multiplicity factor",
                       "peak_merge_tolerance_degrees": MERGE_TOLERANCE_DEGREES,
                       "numerical_extinction_threshold_f_squared": float(electrons ** 2 * 1e-14),
                       "reciprocal_candidates": candidates, "non_extinct_reflections": len(hkls),
                       "profile": "Area-normalized Gaussians, +/-5 sigma support, neighboring peaks included at range edges"},
            "warnings": warnings, "exports": {"csv": csv}}
