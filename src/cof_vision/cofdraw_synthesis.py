"""Offline RDKit SA Score for complete finite precursors, never a whole COF.

This wraps the published RDKit Contrib implementation without retraining or
altering its score. Its bundled fragment table is hash checked before RDKit
loads the trusted pickle. No user-supplied tables, web services, or credentials
are accepted. Graph parsing and caches are bounded.
"""
from __future__ import annotations

from copy import deepcopy
from functools import lru_cache
import hashlib
import math
from pathlib import Path
import threading

from rdkit import Chem, rdBase

METHOD_ID = "rdkit_sa_score"
VERSION = "cofdraw-sa-1"
FRAGMENT_SCORES_SHA256 = "10dcef9340c873e7b987924461b0af5365eb8dd96be607203debe8ddf80c1e73"
FRAGMENT_SCORES_BYTES = 3_848_394
MAX_ATOMS = 250
MAX_SMILES_LENGTH = 4000
ALLOWED_ELEMENTS = frozenset((1, 5, 6, 7, 8, 9, 14, 15, 16, 17, 34, 35, 53))
LIMITATIONS = (
    "A fragment-frequency and molecular-complexity heuristic for finite molecules; "
    "1 is easier and 10 is harder on this scale, not a probability.",
    "Not ASKCOS, an ASKCOS-distilled model, a retrosynthetic route, commercial "
    "availability, cost, safety, reaction feasibility, or experimental yield.",
    "Not validated specifically for COF precursors; no whole-COF synthesis, "
    "crystallization, stability, topology, or polymerization prediction.",
    "Uses the supplied graph without salt stripping, neutralization, tautomer "
    "selection, or correction of uncertain source annotations.",
    "Isotope-label procurement/synthesis cost and reaction conditions are not modeled.",
)
_LOAD_LOCK = threading.Lock()
_SCORER = None


def _metadata():
    return {
        "method": "RDKit Contrib SA Score (Ertl-Schuffenhauer; RDKit implementation)",
        "method_id": METHOD_ID,
        "version": VERSION,
        "method_version": "RDKit SA_Score 1.1",
        "rdkit_version": rdBase.rdkitVersion,
        "scale": {"min": 1, "max": 10, "direction": "lower_is_easier"},
        "reference": "https://doi.org/10.1186/1758-2946-1-8",
        "limitations": list(LIMITATIONS),
    }


def resource_information(directory=None):
    """Verify the installed or build-staged, pinned fragment table without loading it.

    ``directory`` is for the local build validator only, never an API parameter.
    The returned record has no filesystem paths. A changed table requires an
    explicit reviewed update to the pinned digest and scientific tests.
    """
    if directory is None:
        from rdkit.Contrib.SA_Score import sascorer
        directory = Path(sascorer.__file__).parent
    table = Path(directory) / "fpscores.pkl.gz"
    try:
        if table.stat().st_size != FRAGMENT_SCORES_BYTES:
            raise ValueError("size")
        digest = hashlib.sha256(table.read_bytes()).hexdigest()
        if digest != FRAGMENT_SCORES_SHA256:
            raise ValueError("digest")
    except (OSError, ValueError):
        raise RuntimeError("The bundled SA Score fragment table is missing or failed verification") from None
    return {"filename": "fpscores.pkl.gz", "sha256": digest,
            "bytes": FRAGMENT_SCORES_BYTES, "verified": True}


def _load_scorer():
    global _SCORER
    with _LOAD_LOCK:
        if _SCORER is None:
            from rdkit.Contrib.SA_Score import sascorer
            resource_information()
            # Never delegate path lookup or deserialization to a user input.
            # Reload even if another module populated the upstream global first.
            sascorer.readFragmentScores(str(Path(sascorer.__file__).with_name("fpscores.pkl.gz")))
            _SCORER = sascorer
        return _SCORER


def _unavailable(code, reason, canonical=None):
    return {**_metadata(), "available": False, "score": None,
            "reason_code": code, "reason": reason, "canonical_smiles": canonical}


def _parse(smiles):
    if not isinstance(smiles, str) or not 1 <= len(smiles) <= MAX_SMILES_LENGTH:
        return None, _unavailable("invalid_input", "Supply a complete finite monomer SMILES of 1–4000 characters")
    params = Chem.SmilesParserParams()
    params.parseName = False
    params.allowCXSMILES = False
    with rdBase.BlockLogs():
        mol = Chem.MolFromSmiles(smiles, params)
    if mol is None or mol.GetNumAtoms() == 0:
        return None, _unavailable("invalid_graph", "The monomer SMILES does not parse as a valid molecular graph")
    if mol.GetNumAtoms() > MAX_ATOMS:
        return None, _unavailable("atom_limit", "SA Score is limited here to finite monomers with at most 250 atoms")
    if any(a.GetAtomicNum() == 0 for a in mol.GetAtoms()):
        return None, _unavailable("continuation_atoms", "Wildcard or continuation atoms are not complete finite precursors")
    if len(Chem.GetMolFrags(mol)) != 1:
        return None, _unavailable("disconnected_graph", "Disconnected salts or mixtures are unsupported; no component is silently selected")
    if any(a.GetAtomicNum() not in ALLOWED_ELEMENTS for a in mol.GetAtoms()):
        return None, _unavailable("unsupported_elements", "This molecular heuristic is limited here to H, B, C, N, O, F, Si, P, S, Se, Cl, Br and I")
    if any(a.GetNumRadicalElectrons() for a in mol.GetAtoms()):
        return None, _unavailable("radical_graph", "Open-shell radical precursors are outside this scoring scope")
    allowed_bonds = (Chem.BondType.SINGLE, Chem.BondType.DOUBLE, Chem.BondType.TRIPLE, Chem.BondType.AROMATIC)
    if any(b.GetBondType() not in allowed_bonds for b in mol.GetBonds()):
        return None, _unavailable("unsupported_bonds", "Only ordinary covalent finite molecular graphs are supported")
    return mol, None


@lru_cache(maxsize=4096)
def _score_canonical(canonical):
    scorer = _load_scorer()
    mol = Chem.MolFromSmiles(canonical)
    score = float(scorer.calculateScore(mol))
    if not math.isfinite(score) or not 1 <= score <= 10:
        raise ValueError("Invalid SA Score result")
    return {**_metadata(), "available": True, "score": score,
            "reason_code": None, "reason": None, "canonical_smiles": canonical,
            "fragment_scores_sha256": FRAGMENT_SCORES_SHA256,
            "fragment_table_verified": True, "heavy_atoms": mol.GetNumHeavyAtoms(),
            "isotope_labels_present": any(a.GetIsotope() for a in mol.GetAtoms()),
            "formal_charge": Chem.GetFormalCharge(mol)}


def synthetic_accessibility(smiles):
    """Return a qualified 1–10 finite-molecule score, or an explicit unavailable result."""
    mol, error = _parse(smiles)
    if error is not None:
        return error
    canonical = Chem.MolToSmiles(mol, isomericSmiles=True)
    try:
        # Cache only validated canonical graphs; caller edits cannot poison it.
        return deepcopy(_score_canonical(canonical))
    except (ImportError, OSError, RuntimeError):
        return _unavailable("resource_unavailable", "The local SA Score implementation or verified fragment table is unavailable", canonical)
    except (ValueError, TypeError, ZeroDivisionError):
        return _unavailable("calculation_failed", "SA Score could not evaluate this finite molecular graph", canonical)


def catalog_scores(rows):
    """Score actual supplied catalog monomers, retaining IDs and support qualifiers."""
    if not isinstance(rows, (list, tuple)) or len(rows) > 1000:
        raise ValueError("Supply at most 1000 catalog monomer records")
    result = []
    for row in rows:
        if not isinstance(row, dict):
            raise ValueError("Catalog monomers must be records")
        item = {"id": row.get("id"), "label": row.get("label", row.get("id")),
                "smiles": row.get("smiles"),
                "synthetic_accessibility": synthetic_accessibility(row.get("smiles"))}
        # A score does not turn reference-only chemistry into an implemented reaction.
        for key in ("qualification", "source_dois", "ground_truth", "reviewer_type", "compatible_roles"):
            if key in row:
                item[key] = deepcopy(row[key])
        result.append(item)
    return {**_metadata(), "monomers": result, "count": len(result),
            "available_count": sum(r["synthetic_accessibility"]["available"] for r in result),
            "scope": "Scores of the supplied finite catalog monomers, including reference-only entries; not reaction support or COF synthesis scores"}


def precursor_scores(recipe):
    """Summarize distinct finite precursor graphs explicitly declared in a recipe.

    Both fragment-built and explicit-source recipes use components[*].monomer_smiles.
    Product fragments, patch SMILES, and instances are never substituted. Repeated
    component aliases of the same canonical graph count once; all component names
    remain visible. An incomplete set has no numeric aggregate.
    """
    if not isinstance(recipe, dict):
        raise ValueError("Supply a drawing recipe record")
    components = recipe.get("components", {})
    if not isinstance(components, dict) or len(components) > 100:
        raise ValueError("Recipe components must be a mapping of at most 100 precursor records")
    precursors, seen = [], {}
    for name, row in components.items():
        row = row if isinstance(row, dict) else {}
        smiles = row.get("monomer_smiles")
        score = synthetic_accessibility(smiles)
        canonical = score.get("canonical_smiles")
        if score["available"] and canonical in seen:
            seen[canonical]["components"].append(str(name))
            continue
        item = {"component": str(name), "components": [str(name)],
                "label": row.get("label", str(name)), "smiles": smiles,
                "synthetic_accessibility": score}
        precursors.append(item)
        if score["available"]:
            seen[canonical] = item
    available = [p["synthetic_accessibility"]["score"] for p in precursors if p["synthetic_accessibility"]["available"]]
    complete = bool(precursors) and len(available) == len(precursors)
    return {**_metadata(), "precursors": precursors,
            "available_count": len(available), "total_count": len(precursors),
            "component_count": len(components), "complete": complete,
            "mean": sum(available) / len(available) if complete else None,
            "max": max(available) if complete else None,
            "aggregation": "unweighted distinct finite precursor graphs; repeated instances do not change the summary",
            "reason": None if complete else "Not all declared precursors have an available finite-molecule score" if precursors else "No finite precursor monomer SMILES are declared in this recipe",
            "scope": "Precursor SA Score summary only; not a score of the drawing, reaction, infinite COF or its synthesizability"}
