"""Derive the compact CC BY reference table and benchmark the offline screen."""
from __future__ import annotations

import argparse
from collections import defaultdict
import csv
from dataclasses import asdict
from functools import lru_cache
import hashlib
import json
import math
from pathlib import Path
import sys

import numpy as np
from rdkit import Chem, rdBase
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem.Scaffolds import MurckoScaffold

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))
from cof_vision.cofdraw_solubility import (  # noqa: E402
    DATASET, DATASET_DOI, MIN_SOLUTES, SOLVENTS, TARGET_K, WINDOW_K, Observation, ReferenceSet,
    neighbor_estimate, solvent_trend,
)

SOURCE_SHA256 = {
    "BigSolDBv2.2.csv": "6b63f86601369f92530bedf6b93c3fb92a2731390b28a0aabf55ad005acba401",
    "BigSolDBv2.2_methods.csv": "1ed6e3aadfd6095c4ff3e5f9307087ba7ec4a688201b94e8b5efc4827a635eae",
}


@lru_cache(maxsize=None)
def canonical_molecule(smiles):
    with rdBase.BlockLogs():
        mol = Chem.MolFromSmiles(smiles)
    if not mol or len(Chem.GetMolFrags(mol)) != 1 or any(a.GetAtomicNum() == 0 for a in mol.GetAtoms()):
        return None
    return Chem.MolToSmiles(mol, isomericSmiles=True)


def observations(directory):
    for name, checksum in SOURCE_SHA256.items():
        if hashlib.sha256((directory / name).read_bytes()).hexdigest() != checksum:
            raise ValueError(f"Unexpected BigSolDB input hash for {name}")
    with (directory / "BigSolDBv2.2_methods.csv").open(encoding="utf-8-sig", newline="") as source:
        methods = {r["doi"].lower(): r["method"] for r in csv.DictReader(source)
                   if r["solubility_type"] == "equilibrium"}
    solvents = {canonical_molecule(smiles): key for key, _, smiles in SOLVENTS}
    selected = {}
    with (directory / "BigSolDBv2.2.csv").open(encoding="utf-8-sig", newline="") as source:
        for row in csv.DictReader(source):
            doi = row["Source"].strip().lower()
            temperature = float(row["Temperature_K"])
            distance = abs(temperature - TARGET_K)
            if distance > WINDOW_K + 1e-8 or doi not in methods:
                continue
            solvent = solvents.get(canonical_molecule(row["SMILES_Solvent"]))
            smiles = canonical_molecule(row["SMILES_Solute"])
            if not solvent or not smiles or smiles == canonical_molecule(row["SMILES_Solvent"]):
                continue
            fraction = float(row["Solubility(mole_fraction)"] or "nan")
            concentration = float(row["Solubility(mol/L)"] or "nan")
            if not 0 < fraction < 1 or not math.isfinite(concentration) or concentration <= 0:
                continue
            key = smiles, solvent, doi, row["solid_form_label"].strip()
            obs = Observation(smiles, solvent, temperature, fraction, concentration, doi, key[3], methods[doi])
            rank = distance, temperature, concentration
            if key not in selected or rank < selected[key][0]:
                selected[key] = rank, obs
    return [selected[key][1] for key in sorted(selected)]


def scaffold_folds(dataset):
    folds = []
    for smiles, mol in zip(dataset.smiles, dataset.mols):
        scaffold = MurckoScaffold.MurckoScaffoldSmiles(mol=mol) or smiles
        folds.append(int.from_bytes(hashlib.sha256(scaffold.encode()).digest()[:4], "big") % 5)
    return np.array(folds)


def statistics(errors, baselines):
    if not errors:
        return {"predicted": 0}
    return {
        "predicted": len(errors),
        "mae_log10": round(float(np.mean(errors)), 4),
        "median_baseline_mae_log10": round(float(np.mean(baselines)), 4),
        "p90_absolute_error_log10": round(float(np.quantile(errors, .90, method="higher")), 4),
    }


def benchmark(dataset):
    folds = scaffold_folds(dataset)
    # Two aromatic rings are a broader proxy than the small set of measured
    # exact COF precursors. The benchmark still spans ordinary organics.
    aromatic = [rdMolDescriptors.CalcNumAromaticRings(mol) >= 2 for mol in dataset.mols]
    errors, baselines, aromatic_errors, aromatic_baselines = (defaultdict(list) for _ in range(4))
    pools = {}
    for solvent, (indices, values) in dataset.pools.items():
        for fold in range(5):
            train = folds[indices] != fold
            if train.sum() >= MIN_SOLUTES:
                pool = indices[train], values[train]
                pools[solvent, fold] = *pool, solvent_trend(dataset.features, *pool)
    for index, mol in enumerate(dataset.mols):
        resemblance, distances = dataset.distances(mol)
        for solvent, (indices, values) in dataset.pools.items():
            where = np.searchsorted(indices, index)
            if where == len(indices) or indices[where] != index:
                continue
            pool = pools.get((solvent, folds[index]))
            if pool is None:
                continue
            training_indices, training_values, trend = pool
            prediction = neighbor_estimate(dataset.features[index], resemblance, distances, training_indices,
                                           training_values, dataset.features, trend)
            if prediction is None:
                continue
            actual = values[where]
            error = abs(actual - prediction[0])
            baseline = abs(actual - float(np.median(training_values)))
            errors[solvent].append(error)
            baselines[solvent].append(baseline)
            if aromatic[index]:
                aromatic_errors[solvent].append(error)
                aromatic_baselines[solvent].append(baseline)

    result = {}
    all_aromatic_errors = [error for errors_in_solvent in aromatic_errors.values() for error in errors_in_solvent]
    pooled_aromatic_band = float(np.quantile(all_aromatic_errors, .90, method="higher"))
    for solvent, (indices, _) in dataset.pools.items():
        stats = statistics(errors[solvent], baselines[solvent])
        aromatic_stats = statistics(aromatic_errors[solvent], aromatic_baselines[solvent])
        enabled = stats["predicted"] >= 25 and stats["mae_log10"] < stats["median_baseline_mae_log10"]
        if aromatic_stats["predicted"] >= 20:
            enabled &= aromatic_stats["mae_log10"] < aromatic_stats["median_baseline_mae_log10"]
        aromatic_band = (aromatic_stats["p90_absolute_error_log10"] if aromatic_stats["predicted"] >= 20
                         else pooled_aromatic_band)
        band = round(max(.75, stats["p90_absolute_error_log10"], aromatic_band), 4) if enabled else None
        result[solvent] = {
            "enabled": enabled, "reference_solutes": len(indices), **stats,
            "screening_half_width_log10": band,
            "two_or_more_aromatic_rings": {
                "reference_solutes": sum(aromatic[i] for i in indices), **aromatic_stats,
            },
        }
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source", type=Path, help="Directory with both pinned BigSolDB 2.2 files")
    parser.add_argument("--output", type=Path, default=REPO / "data/solubility")
    args = parser.parse_args()
    rows = observations(args.source)
    dataset = ReferenceSet(rows)
    validation = benchmark(dataset)
    args.output.mkdir(parents=True, exist_ok=True)
    table = args.output / "measurements.csv"
    with table.open("w", encoding="utf-8", newline="") as destination:
        writer = csv.DictWriter(destination, fieldnames=list(asdict(rows[0])), lineterminator="\n")
        writer.writeheader()
        writer.writerows(asdict(row) for row in rows)
    metadata = {
        "schema": "cofdraw-solubility-1",
        "dataset": DATASET, "dataset_doi": DATASET_DOI, "license": "CC-BY-4.0",
        "source_sha256": SOURCE_SHA256,
        "measurements_sha256": hashlib.sha256(table.read_bytes()).hexdigest(),
        "rows": len(rows), "unique_solutes": len({row.smiles for row in rows}),
        "neutral_organic_estimator_solutes": len(dataset.smiles),
        "preparation_rdkit_version": rdBase.rdkitVersion,
        "protocol": {"split": "Five folds by Bemis–Murcko scaffold; acyclic molecules grouped by canonical SMILES.",
                     "baseline": "Training-solute median in the same solvent, evaluated on the same held-out predictions.",
                     "band": "Larger of the full and multi-aromatic held-out 90th-percentile absolute errors, with a 0.75 log-unit floor. The multi-aromatic error is pooled across solvents if the solvent has fewer than 20 aromatic estimates.",
                     "eligibility": "At least 25 held-out predictions and lower MAE than the solvent median. Also require improvement for the multi-aromatic subset when it has 20 or more estimates."},
        "solvents": validation,
    }
    (args.output / "validation.json").write_text(json.dumps(metadata, indent=2, ensure_ascii=False) + "\n", encoding="utf-8", newline="\n")
    print(f"{len(rows)} measurements and {len(dataset.smiles)} estimator solutes")
    for key, stats in validation.items():
        print(key, stats)


if __name__ == "__main__":
    main()
