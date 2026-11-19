"""Build the offline explorer from the merged chemistry audit, mapped master and bundled CIFs.

Build environment additionally needs scikit-learn. The desktop app reads the resulting JSON.
"""
import argparse
import csv
import hashlib
import json
import re
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np
from rdkit import Chem
from rdkit.Chem import rdFingerprintGenerator
from rdkit.Chem.Draw import rdMolDraw2D
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
from threadpoolctl import threadpool_limits

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from cof_vision.cofdraw_science import import_cif, pore_analysis

# Representative product motifs. Broad master classes retain their original scope.
FAMILIES = {
    "imine": ("Imine", "*c1ccc(C=Nc2ccc(*)cc2)cc1"),
    "beta-ketoenamine": ("β-Ketoenamine", "*C1=CC(=O)C(=CNc2ccc(*)cc2)C=C1"),
    "c=c": ("C=C linked", "*c1ccc(C=Cc2ccc(*)cc2)cc1"),
    "hydrazone": ("Hydrazone", "*c1ccc(C=NNC(=O)c2ccc(*)cc2)cc1"),
    "azine": ("Azine", "*c1ccc(C=N/N=Cc2ccc(*)cc2)cc1"),
    "triazine": ("Triazine", "*c1nc(*)nc(*)n1"),
    "boronate-ester": ("Boronate ester", "*B1Oc2cc(*)ccc2O1"),
    "boroxine": ("Boroxine", "*B1OB(*)OB(*)O1"),
    "imide": ("Imide", "*N1C(=O)c2cc(*)ccc2C1=O"),
    "imidazole": ("Imidazole", "*c1nc(*)c[nH]1"),
    "benzoxazole": ("Benzoxazole", "*c1nc2ccc(*)cc2o1"),
    "benzobisoxazole": ("Benzobisoxazole", "*c1nc2cc3oc(*)nc3cc2o1"),
    "quinoline": ("Quinoline", "*c1ccc2ncc(*)cc2c1"),
    "amide": ("Amide", "*c1ccc(C(=O)Nc2ccc(*)cc2)cc1"),
    "squaraine": ("Squaraine", "*c1ccc([C+]2C(=O)C(c3ccc(*)cc3)=C2[O-])cc1"),
    "c-n": ("C–N linked", "*c1ccc(Nc2ccc(*)cc2)cc1"),
    "c-o": ("C–O linked", "*c1ccc(Oc2ccc(*)cc2)cc1"),
    "c-s": ("C–S linked", "*c1ccc(Sc2ccc(*)cc2)cc1"),
    "keto-hydrazide": ("Ketoenamine–hydrazide", "*C(=O)NNC=C(*)C=O"),
    "nitrone": ("Nitrone", "*c1ccc(C=[N+]([O-])c2ccc(*)cc2)cc1"),
    "reduced-imine": ("Reduced imine", "*c1ccc(CNc2ccc(*)cc2)cc1"),
    "phenazine": ("Phenazine", "*c1ccc2nc3ccc(*)cc3nc2c1"),
    "benzothiazole": ("Benzothiazole", "*c1nc2ccc(*)cc2s1"),
    "azo": ("Azo", "*c1ccc(N=Nc2ccc(*)cc2)cc1"),
    "dioxin": ("1,4-Dioxin", "*c1ccc2c(c1)Oc1ccc(*)cc1O2"),
    "pyridine": ("Pyridine", "*c1cc(*)ncc1*"),
}
ALIASES = {"beta_ketoenamine": "beta-ketoenamine", "acylhydrazone": "hydrazone", "keto_hydrazide": "keto-hydrazide"}
ROLE = {k: Chem.MolFromSmarts(v) for k, v in {
    "amine": "[NX3;H2]-[c]", "aldehyde": "[CX3H1](=O)", "boron": "[B]([OH])[OH]",
}.items()}


def rows(path):
    with path.open(encoding="utf-8-sig", newline="") as stream:
        yield from csv.DictReader(stream)


def canonical(smiles):
    mol = Chem.MolFromSmiles(smiles) if smiles else None
    if mol is None:
        return None
    for atom in mol.GetAtoms():
        atom.SetAtomMapNum(0)
    return Chem.MolToSmiles(mol)


def motif(smiles):
    mol = Chem.MolFromSmiles(smiles)
    draw = rdMolDraw2D.MolDraw2DSVG(180, 70)
    draw.drawOptions().clearBackground = False
    rdMolDraw2D.PrepareAndDrawMolecule(draw, mol)
    draw.FinishDrawing()
    return "<svg" + draw.GetDrawingText().split("<svg", 1)[1]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("inventory", type=Path)
    parser.add_argument("mapped_master", type=Path)
    parser.add_argument("--output", type=Path, default=ROOT / "src/cof_vision/web/chemical-space.json")
    args = parser.parse_args()
    blocks, by_smiles, fingerprints, roles = {}, {}, {}, {}
    generator = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=1024)
    for row in rows(args.inventory / "merged-building-blocks.csv"):
        idx = int(row["audit ID"])
        smiles = row["canonical isomeric SMILES"]
        mol = Chem.MolFromSmiles(smiles)
        blocks[idx] = row["label"]
        by_smiles[canonical(smiles)] = idx
        fingerprints[idx] = generator.GetFingerprintAsNumPy(mol)
        roles[idx] = {role for role, query in ROLE.items() if mol.HasSubstructMatch(query)}
    evidence = defaultdict(set)
    for row in rows(args.mapped_master):
        family = row["linkage"].strip().lower()
        doi = row["doi"].strip().lower()
        if family not in FAMILIES or row["article_trial_or_failure"].lower() != "no" or not re.fullmatch(r"10\.\d{4,9}/\S+", doi):
            continue
        count = int(row["bb_total"] or 0)
        if count < 1 or count > 4 or int(row["bb_mapped"] or 0) != count:
            continue
        ids = [by_smiles.get(canonical(row[f"bb_{i}_smiles"])) for i in range(1, count+1)]
        if None in ids or len(set(ids)) != count:
            continue
        evidence[family, tuple(sorted(ids))].add((row["cof_name"].strip() or "Unnamed COF", doi))

    family_ids = sorted(FAMILIES)
    points, vectors = [], []
    for (family, ids), refs in sorted(evidence.items()):
        points.append({"f": family, "ids": ids, "refs": sorted(refs)})
        # Mean pooling is invariant to the order of the precursor list.
        one_hot = np.array([2.0 * (family == k) for k in family_ids])
        vectors.append(np.r_[np.mean([fingerprints[i] for i in ids], axis=0), one_hot])
    with threadpool_limits(limits=1):
        matrix = PCA(n_components=min(48,len(vectors)-1), random_state=42).fit_transform(np.asarray(vectors))
        xy = TSNE(perplexity=min(30,len(points)-1), init="pca", learning_rate="auto", random_state=42).fit_transform(matrix)
    for point, coords in zip(points, xy):
        point["xy"] = [round(float(v), 3) for v in coords]

    stride = max(blocks) + 1
    eligible = defaultdict(set)
    for row in rows(args.inventory / "compatible-pairs.csv"):
        a, b = sorted(int(row[k]) for k in ("block 1 ID", "block 2 ID"))
        for key in row["all supported linkages"].split(";"):
            eligible[ALIASES.get(key, key)].add(a*stride+b)
    families = {k: {"name": label, "svg": motif(smiles), "eligible": sorted(eligible[k])}
                for k, (label, smiles) in FAMILIES.items()}
    models = ROOT / "data/crystal-models"
    manifest = json.loads((models / "manifest.json").read_text("utf8"))
    properties = []
    for row in manifest["curated_models"] + manifest["reported_models"]:
        result = pore_analysis(import_cif((models/row["path"]).read_text("utf8")), samples=2048, seed=42, grid_size=8, surface_samples=16)
        properties.append(dict(name=row["model_name"], doi=row["citation_doi"], id=row["id"],
                               density=round(result["density_g_cm3"],4), void=round(result["void_fraction"]*100,2),
                               diameter=round(result["largest_sampled_included_sphere_diameter_angstrom"],2),
                               samples=result["method"]["volume_samples"], sha256=row["sha256"]))
    metadata = json.loads((args.inventory / "counts.json").read_text("utf8"))
    # The literature index can advance independently of the compatibility audit.
    evidence_hash = hashlib.sha256(args.mapped_master.read_bytes()).hexdigest()
    result = dict(schema=1, stride=stride, blocks=blocks, roles={k:sorted(v) for k,v in roles.items()},
                  families=families, points=points, properties=properties, source_hash=evidence_hash,
                  compatibility_source_hash=metadata["mapped_file_sha256"])
    args.output.write_text(json.dumps(result, ensure_ascii=False, separators=(",",":")), encoding="utf8")
    print(f"{len(points)} compositions, {len(properties)} CIFs, {len(blocks)} building blocks")


if __name__ == "__main__":
    main()
