"""Reagent weigh-outs based on complementary functional groups of finite monomers."""
from math import gcd

from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors

from .cofdraw_chemistry import inspect_monomer
from .framework_drawing import LINKAGE_LIBRARY


def reagent_balance(monomer_a, monomer_b, linkage, mass_a_mg, purity_a_percent,
                    purity_b_percent, b_site_equivalents):
    if linkage not in LINKAGE_LIBRARY:
        raise ValueError("Choose an implemented linkage for the reagent balance.")
    roles = LINKAGE_LIBRARY[linkage]["roles"]
    rows = []
    functions = []
    for slot, smiles, purity in zip("AB", (monomer_a, monomer_b), (purity_a_percent, purity_b_percent)):
        mol = Chem.MolFromSmiles(smiles)
        if (not mol or not mol.GetNumAtoms() or len(Chem.GetMolFrags(mol)) != 1
                or any(atom.GetAtomicNum() == 0 for atom in mol.GetAtoms())):
            raise ValueError(f"Block {slot} needs one complete molecular structure.")
        functions.append({r["role"]: r["reactive_sites"] for r in inspect_monomer(smiles)["compatible_roles"]
                          if r["linkage"] == linkage and r["reactive_sites"] >= 2})
        rows.append({"block": slot, "formula": rdMolDescriptors.CalcMolFormula(mol, separateIsotopes=True),
                     "smiles": Chem.MolToSmiles(mol), "molar_mass_g_mol": Descriptors.MolWt(mol),
                     "purity_percent": purity})
    orientations = [pair for pair in dict.fromkeys((tuple(roles), tuple(reversed(roles))))
                    if all(role in counts for counts, role in zip(functions, pair))]
    if len(orientations) != 1:
        raise ValueError("The selected linkage needs one unambiguous pair of complementary reactive roles.")
    for row, role, counts in zip(rows, orientations[0], functions):
        row.update(role=role, sites_per_molecule=counts[role])
    a, b = rows
    a_mmol = mass_a_mg * purity_a_percent / 100 / a["molar_mass_g_mol"]
    b_mmol = a_mmol * a["sites_per_molecule"] * b_site_equivalents / b["sites_per_molecule"]
    for row, mmol in zip(rows, (a_mmol, b_mmol)):
        row.update(monomer_mmol=mmol, site_mmol=mmol * row["sites_per_molecule"],
                   weigh_mg=mmol * row["molar_mass_g_mol"] * 100 / row["purity_percent"])
    divisor = gcd(a["sites_per_molecule"], b["sites_per_molecule"])
    return {"reagents": rows, "balanced_molar_ratio_a_b": [b["sites_per_molecule"] // divisor, a["sites_per_molecule"] // divisor],
            "b_site_equivalents": b_site_equivalents, "linkage": linkage,
            "scope": "Nominal feed assuming every counted site reacts once. Purity is the mass fraction of the drawn species. Conversion and isolated yield are unknown."}
