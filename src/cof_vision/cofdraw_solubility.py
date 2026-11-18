"""Room-temperature monomer solubility screening from equilibrium measurements."""
from __future__ import annotations

from collections import defaultdict
import csv
from dataclasses import dataclass
from functools import lru_cache
import hashlib
import io
import json
import math
from pathlib import Path

import numpy as np
from rdkit import Chem, DataStructs, rdBase
from rdkit.Chem import Crippen, Descriptors, rdFingerprintGenerator, rdMolDescriptors


# IDs are stable API inputs; canonical solvent graphs also disambiguate aliases.
SOLVENTS = (
    ("dmf", "DMF", "CN(C)C=O"),
    ("dmac", "DMAc", "CC(=O)N(C)C"),
    ("dmso", "DMSO", "CS(C)=O"),
    ("nmp", "NMP", "CN1CCCC1=O"),
    ("dioxane", "1,4-Dioxane", "C1COCCO1"),
    ("thf", "THF", "C1CCOC1"),
    ("mesitylene", "Mesitylene", "Cc1cc(C)cc(C)c1"),
    ("odcb", "o-Dichlorobenzene", "Clc1ccccc1Cl"),
    ("acoh", "Acetic acid", "CC(=O)O"),
    ("toluene", "Toluene", "Cc1ccccc1"),
    ("acetone", "Acetone", "CC(C)=O"),
    ("mecn", "Acetonitrile", "CC#N"),
    ("etoac", "Ethyl acetate", "CCOC(C)=O"),
    ("dcm", "Dichloromethane", "ClCCl"),
    ("chloroform", "Chloroform", "ClC(Cl)Cl"),
    ("ethanol", "Ethanol", "CCO"),
    ("methanol", "Methanol", "CO"),
    ("ipoh", "Isopropanol", "CC(C)O"),
    ("water", "Water", "O"),
    ("buth", "1-Butanol", "CCCCO"),
    ("proh", "1-Propanol", "CCCO"),
    ("iboh", "Isobutanol", "CC(C)CO"),
    ("sboh", "2-Butanol", "CCC(C)O"),
    ("tboh", "tert-Butanol", "CC(C)(C)O"),
    ("pentoh", "1-Pentanol", "CCCCCO"),
    ("hexoh", "1-Hexanol", "CCCCCCO"),
    ("octoh", "1-Octanol", "CCCCCCCCO"),
    ("eg", "Ethylene glycol", "OCCO"),
    ("mek", "2-Butanone", "CCC(C)=O"),
    ("cyone", "Cyclohexanone", "O=C1CCCCC1"),
    ("meoac", "Methyl acetate", "COC(C)=O"),
    ("buoac", "n-Butyl acetate", "CCCCOC(C)=O"),
    ("iproac", "Isopropyl acetate", "CC(C)OC(C)=O"),
    ("ether", "Diethyl ether", "CCOCC"),
    ("mtbe", "Methyl tert-butyl ether", "COC(C)(C)C"),
    ("cpme", "Cyclopentyl methyl ether", "COC1CCCC1"),
    ("hexane", "n-Hexane", "CCCCCC"),
    ("heptane", "n-Heptane", "CCCCCCC"),
    ("cyclohexane", "Cyclohexane", "C1CCCCC1"),
    ("dce", "1,2-Dichloroethane", "ClCCCl"),
    ("phcl", "Chlorobenzene", "Clc1ccccc1"),
    ("anisole", "Anisole", "COc1ccccc1"),
    ("oxylene", "o-Xylene", "Cc1ccccc1C"),
    ("mxylene", "m-Xylene", "Cc1cccc(C)c1"),
    ("pxylene", "p-Xylene", "Cc1ccc(C)cc1"),
    ("pc", "Propylene carbonate", "CC1COC(=O)O1"),
)
SOLVENT_NAMES = {key: name for key, name, _ in SOLVENTS}
DATASET = "BigSolDB 2.2"
DATASET_DOI = "10.5281/zenodo.22648301"
TARGET_K = 298.15
WINDOW_K = 5.0
NUM_NEIGHBORS = 8
MIN_SOLUTES = 25
MIN_FINGERPRINT = 0.18
MAX_DESCRIPTOR_DISTANCE = 1.15
_SCALE = np.array([110, 1.6, 45, 2, 3, 3, 1.5, 0.30])
_FP = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048, includeChirality=True)
_ORDINARY_ELEMENTS = {1, 5, 6, 7, 8, 9, 14, 15, 16, 17, 35, 53}


def prediction_domain(mol):
    """The estimator uses neutral, ordinary organic structures with defined valence."""
    return mol.GetNumHeavyAtoms() >= 3 and any(a.GetAtomicNum() == 6 for a in mol.GetAtoms()) and all(
        a.GetAtomicNum() in _ORDINARY_ELEMENTS
        and not (a.GetFormalCharge() or a.GetNumRadicalElectrons() or a.GetIsotope())
        for a in mol.GetAtoms()
    )


def molecular_features(mol):
    return np.array([
        Descriptors.MolWt(mol), Crippen.MolLogP(mol), rdMolDescriptors.CalcTPSA(mol),
        rdMolDescriptors.CalcNumHBD(mol), rdMolDescriptors.CalcNumHBA(mol),
        rdMolDescriptors.CalcNumRotatableBonds(mol), rdMolDescriptors.CalcNumAromaticRings(mol),
        rdMolDescriptors.CalcFractionCSP3(mol),
    ]) / _SCALE


@dataclass(frozen=True)
class Observation:
    smiles: str
    solvent: str
    temperature_k: float
    mole_fraction: float
    mol_per_liter_solvent: float
    doi: str
    solid_form: str
    method: str


class ReferenceSet:
    def __init__(self, observations):
        self.exact = defaultdict(list)
        for obs in observations:
            self.exact[obs.smiles, obs.solvent].append(obs)
        mols = {}
        for smiles, _ in self.exact:
            if smiles not in mols:
                mols[smiles] = Chem.MolFromSmiles(smiles)
        self.smiles = sorted(s for s, mol in mols.items() if mol is not None and prediction_domain(mol))
        self.mols = [mols[s] for s in self.smiles]
        self.weights = [Descriptors.MolWt(mol) for mol in self.mols]
        self.features = np.array([molecular_features(mol) for mol in self.mols])
        self.fingerprints = [_FP.GetFingerprint(mol) for mol in self.mols]
        self.pools = {}
        for solvent in SOLVENT_NAMES:
            indices, values = [], []
            for index, smiles in enumerate(self.smiles):
                # Combine the nearest room-temperature values from independent sources
                # and reported solid forms into one log-space value per solute.
                logs = [math.log10(o.mol_per_liter_solvent) for o in self.exact.get((smiles, solvent), ())]
                if logs:
                    indices.append(index)
                    values.append(float(np.median(logs)))
            self.pools[solvent] = np.asarray(indices, dtype=int), np.asarray(values)
        self.trends = {solvent: solvent_trend(self.features, indices, values)
                       for solvent, (indices, values) in self.pools.items() if len(indices) >= MIN_SOLUTES}
        self.shared_trends = shared_solvent_trends(self.features, self.pools)

    def distances(self, mol):
        resemblance = np.array(DataStructs.BulkTanimotoSimilarity(_FP.GetFingerprint(mol), self.fingerprints))
        distance = np.sqrt(np.mean((self.features - molecular_features(mol)) ** 2, axis=1))
        return resemblance, distance


def solvent_trend(features, indices, values):
    """Fit a small descriptor baseline from the measurements in memory."""
    training = features[indices]
    mean = training.mean(axis=0)
    centered = training - mean
    penalty = 8 * np.eye(centered.shape[1])
    coefficients = np.linalg.solve(centered.T @ centered + penalty, centered.T @ (values - values.mean()))
    return coefficients, float(values.mean() - mean @ coefficients)


def shared_solvent_trends(features, pools):
    """Share descriptor slopes across solvents, retaining each solvent's intercept."""
    gram = np.eye(features.shape[1]) * 8
    rhs = np.zeros(features.shape[1])
    centers = {}
    for key, (ids, values) in pools.items():
        if len(ids) < 8:
            continue
        x = features[ids]
        mean = x.mean(axis=0)
        x = x - mean
        y = values - values.mean()
        local_gram, local_rhs = x.T @ x, x.T @ y
        weight = 1 / math.sqrt(len(ids))
        gram += weight * local_gram
        rhs += weight * local_rhs
        centers[key] = mean, values.mean(), local_gram, local_rhs
    prior = np.linalg.solve(gram, rhs)
    trends = {}
    for key, (center, mean_y, local_gram, local_rhs) in centers.items():
        coef = np.linalg.solve(local_gram + 60*np.eye(len(prior)), local_rhs + 60*prior)
        trends[key] = coef, float(mean_y - center @ coef)
    return trends


def descriptor_estimate(features, distances, indices, trend):
    """Limit transfer to molecular sizes and descriptor neighborhoods seen in this solvent."""
    local = distances[indices]
    if trend is None or len(local) < 8 or np.partition(local, 2)[2] > .95:
        return None
    coef, intercept = trend
    return float(features @ coef + intercept)


def neighbor_estimate(features, resemblance, descriptor_distance, indices, values, training_features, trend):
    similarities = resemblance[indices]
    distances = descriptor_distance[indices]
    qualified = np.flatnonzero((similarities >= MIN_FINGERPRINT) & (distances <= MAX_DESCRIPTOR_DISTANCE))
    if len(qualified) < 3:
        return None
    ranked = (1 - similarities[qualified]) + 0.40 * distances[qualified]
    nearest = np.argsort(ranked, kind="stable")[:NUM_NEIGHBORS]
    positions = qualified[nearest]
    weights = np.exp(-(ranked[nearest] - ranked[nearest[0]]) / 0.17)
    coefficients, intercept = trend
    residuals = values[positions] - (training_features[indices[positions]] @ coefficients + intercept)
    estimate = float(features @ coefficients + intercept + .5 * np.average(residuals, weights=weights))
    return estimate, indices[positions], similarities[positions], weights / weights.sum()


@lru_cache(maxsize=4)
def load_reference(data_directory: Path):
    path = data_directory / "measurements.csv"
    metadata = json.loads((data_directory / "validation.json").read_text("utf-8"))
    raw = path.read_bytes()
    if hashlib.sha256(raw).hexdigest() != metadata["measurements_sha256"]:
        raise ValueError("Bundled solubility measurements failed their checksum")
    reader = csv.DictReader(io.StringIO(raw.decode("utf-8")))
    dataset = ReferenceSet([
        Observation(row["smiles"], row["solvent"], float(row["temperature_k"]), float(row["mole_fraction"]),
                    float(row["mol_per_liter_solvent"]), row["doi"], row["solid_form"], row["method"])
        for row in reader
    ])
    return dataset, metadata


def _amount(log_s, weight, volume):
    concentration = 10 ** log_s
    mg_ml = concentration * weight
    return {"mol_per_liter_solvent": concentration, "mg_per_ml_solvent": mg_ml, "mass_mg": mg_ml * volume}


def screen_solubility(smiles, data_directory, *, solvent_ids=None, solvent_ml=1.0, preparation_loadings=()):
    mol = Chem.MolFromSmiles(smiles)
    if not mol or not 1 <= mol.GetNumAtoms() <= 400 or len(Chem.GetMolFrags(mol)) != 1 or any(
        a.GetAtomicNum() == 0 for a in mol.GetAtoms()
    ):
        raise ValueError("Enter one complete connected monomer with at most 400 atoms")
    if not math.isfinite(solvent_ml) or not 0.01 <= solvent_ml <= 1000:
        raise ValueError("Solvent volume must be from 0.01 to 1000 mL")
    canonical = Chem.MolToSmiles(mol, isomericSmiles=True)
    weight = Descriptors.MolWt(mol)
    solvent_ids = list(SOLVENT_NAMES) if solvent_ids is None else list(solvent_ids)
    if not solvent_ids or len(set(solvent_ids)) != len(solvent_ids) or any(s not in SOLVENT_NAMES for s in solvent_ids):
        raise ValueError("Select distinct supported solvent IDs from the solubility catalog")
    references, metadata = load_reference(Path(data_directory).resolve())
    features = molecular_features(mol) if prediction_domain(mol) else None
    distances = references.distances(mol) if features is not None else None
    rows = []
    for solvent in solvent_ids:
        row = {"solvent_id": solvent, "solvent": SOLVENT_NAMES[solvent], "status": "unavailable"}
        indices, values = references.pools[solvent]
        row["reference_solute_count"] = len(indices)
        exact = references.exact.get((canonical, solvent))
        calibration = metadata["solvents"][solvent]
        if exact:
            observations = sorted(exact, key=lambda o: (abs(o.temperature_k - TARGET_K), o.doi, o.solid_form))
            row["status"] = "reported"
            row["measurements"] = [
                {"temperature_k": o.temperature_k, "mole_fraction": o.mole_fraction, "source_doi": o.doi,
                 "solid_form": o.solid_form or None, "method": o.method,
                 **_amount(math.log10(o.mol_per_liter_solvent), weight, solvent_ml)}
                for o in observations
            ]
            row.update(row["measurements"][0])
        elif distances is None:
            row["reason"] = "The estimator covers neutral organic molecules. No exact report matched this graph."
        elif not calibration["enabled"]:
            row["reason"] = "Insufficient held-out support for a numerical estimate in this solvent."
        else:
            # Existing records always take precedence, so the query is absent
            # from this solvent's prediction pool.
            predicted = neighbor_estimate(features, *distances, indices, values, references.features,
                                          references.trends[solvent])
            if predicted is None:
                row["reason"] = "Too few related measured solutes in this solvent."
            else:
                log_s, neighbors, similarity, fractions = predicted
                error = calibration["screening_half_width_log10"]
                row.update(status="estimated", evaluation=calibration, temperature_k=TARGET_K,
                           **_amount(log_s, weight, solvent_ml))
                row["screening_band"] = {
                    "lower": _amount(log_s - error, weight, solvent_ml),
                    "upper": _amount(log_s + error, weight, solvent_ml),
                }
                row["analogs"] = [
                    {"smiles": references.smiles[i], "morgan_similarity": float(s), "relative_weight": float(w)}
                    for i, s, w in zip(neighbors, similarity, fractions)
                ]
        if row["status"] == "estimated":
            row["estimate_method"] = "same_solvent_analogs"
        elif row["status"] == "unavailable" and distances is not None:
            transfer = calibration.get("descriptor_transfer", {})
            if transfer.get("enabled"):
                log_s = descriptor_estimate(features, distances[1], indices, references.shared_trends.get(solvent))
                if log_s is not None:
                    error = transfer["screening_half_width_log10"]
                    row.pop("reason")
                    row.update(status="estimated", estimate_method="shared_descriptor_regression",
                               evaluation=transfer, temperature_k=TARGET_K, **_amount(log_s, weight, solvent_ml))
                    row["screening_band"] = {"lower": _amount(log_s-error, weight, solvent_ml),
                                             "upper": _amount(log_s+error, weight, solvent_ml)}
            if row["status"] == "unavailable":
                # Only describe observations on the named reference molecule.
                # Its mass concentration must use its own molecular weight.
                sim, dist = distances
                neighbors = indices[sim[indices] >= .20]
                order = np.argsort(1-sim[neighbors] + .16*dist[neighbors], kind="stable")[:3]
                examples = []
                for i in neighbors[order]:
                    graph = references.smiles[i]
                    obs = min(references.exact[graph, solvent],
                              key=lambda o: (abs(o.temperature_k - TARGET_K), o.doi, o.solid_form))
                    examples.append({"smiles": graph, "morgan_similarity": float(sim[i]), "source_doi": obs.doi,
                                     "temperature_k": obs.temperature_k, "solid_form": obs.solid_form,
                                     "mg_per_ml_solvent": obs.mol_per_liter_solvent*references.weights[i]})
                if examples:
                    row["reference_compounds"] = examples
        rows.append(row)

    published_loadings = json.loads((Path(data_directory) / "preparation-loadings.json").read_text("utf-8"))
    loadings = [
        {**record, "solvent": SOLVENT_NAMES[record["solvent_id"]],
         "mg_per_ml_solvent": record["mass_mg"] / record["solvent_ml"],
         "status": "reported_loading"}
        for record in (*preparation_loadings, *published_loadings)
        if record["canonical_smiles"] == canonical and record["solvent_id"] in solvent_ids
    ]
    output = io.StringIO()
    columns = [
        "solvent_id", "solvent", "status", "temperature_k", "mol_per_liter_solvent", "mg_per_ml_solvent",
        "solvent_ml", "mass_mg", "low_mg_per_ml_solvent", "high_mg_per_ml_solvent", "mole_fraction", "source_doi",
        "source_locator", "solid_form", "method", "estimate_method", "reference_solute_count", "reason",
        "reference_smiles", "reference_mg_per_ml_solvent", "reference_similarity",
    ]
    writer = csv.DictWriter(output, fieldnames=columns, extrasaction="ignore")
    writer.writeheader()
    for row in rows:
        row = {**row, "solvent_ml": solvent_ml}
        if row["status"] == "reported":
            for measurement in row["measurements"]:
                writer.writerow({**row, **measurement})
        elif row["status"] == "estimated":
            writer.writerow({**row, "low_mg_per_ml_solvent": row["screening_band"]["lower"]["mg_per_ml_solvent"],
                             "high_mg_per_ml_solvent": row["screening_band"]["upper"]["mg_per_ml_solvent"]})
        else:
            writer.writerow(row)
            for ref in row.get("reference_compounds", ()):
                writer.writerow({"solvent": row["solvent"], "solvent_id": row["solvent_id"],
                                 "status": "reference_compound", "source_doi": ref["source_doi"],
                                 "temperature_k": ref["temperature_k"], "solid_form": ref["solid_form"],
                                 "reference_smiles": ref["smiles"], "reference_similarity": ref["morgan_similarity"],
                                 "reference_mg_per_ml_solvent": ref["mg_per_ml_solvent"]})
    writer.writerows(loadings)
    return {
        "scope": "Equilibrium solubility screening for a single monomer in neat solvents near room temperature.",
        "canonical_smiles": canonical, "formula": rdMolDescriptors.CalcMolFormula(mol),
        "molecular_weight_g_mol": weight, "solvent_ml": solvent_ml, "target_temperature_k": TARGET_K,
        "temperature_window_k": WINDOW_K, "solvents": rows, "preparation_loadings": loadings,
        "method": {"name": "Same-solvent analogs and scaffold-tested shared descriptor regression", "version": metadata["schema"],
                   "dataset": DATASET, "dataset_doi": DATASET_DOI, "license": "CC BY 4.0",
                   "measurements_sha256": metadata["measurements_sha256"], "rdkit_version": rdBase.rdkitVersion,
                   "calibration": metadata["protocol"]},
        "limitations": [
            "Concentrations use the starting volume of neat solvent. Reported mg/mL is converted from published mole fractions using pure-solvent densities.",
            "Analog bands summarize held-out errors. Crystal form, salt state, pH and experimental conditions can shift actual solubility beyond these bands.",
            "Prediction applies near room temperature in individual solvents. Mixed solvents, elevated temperatures and reaction mixtures require their own measurements.",
        ],
        "exports": {"csv": output.getvalue()},
    }
