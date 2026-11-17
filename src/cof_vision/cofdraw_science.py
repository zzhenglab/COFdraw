"""Bounded molecular and periodic-geometry calculations for COFdraw.

All lengths are angstrom, cell vectors are rows, fractional @ cell is Cartesian.
No drawing is silently promoted to a periodic crystal or an adsorption isotherm.
Public functions return JSON-compatible dictionaries and in-memory export text.
"""
from __future__ import annotations

from collections import deque
from fractions import Fraction
import hashlib
import itertools
import json
import math
import re
import shlex

import numpy as np
from rdkit import Chem, rdBase
from rdkit.Chem import AllChem, Crippen, Descriptors, Lipinski, rdMolDescriptors, rdPartialCharges

AVOGADRO = 6.02214076e23
MAX_ATOMS = 2000
PT = Chem.GetPeriodicTable()
ELEMENTS = {PT.GetElementSymbol(z): z for z in range(1, 119)}


def _finite(value, name, low=None, high=None):
    if isinstance(value, (bool, str)):
        raise ValueError(f"{name} must be a finite number")
    try:
        number = float(value)
    except (TypeError, ValueError):
        raise ValueError(f"{name} must be a finite number") from None
    if not math.isfinite(number) or (low is not None and number < low) or (high is not None and number > high):
        raise ValueError(f"{name} must be finite and in [{low}, {high}]")
    return number


def _integer(value, name, low, high):
    number = _finite(value, name, low, high)
    if int(number) != number:
        raise ValueError(f"{name} must be an integer")
    return int(number)


def _bounds(positions):
    xyz = np.asarray(positions, dtype=float)
    return {"min": xyz.min(axis=0).tolist(), "max": xyz.max(axis=0).tolist(),
            "center": ((xyz.min(axis=0) + xyz.max(axis=0)) / 2).tolist()}


def monomer_analysis(smiles, *, conformers=4, seed=42, max_iterations=300):
    """RDKit descriptors and seeded ETKDGv3/MMFF94s (or UFF) conformers.

    Energies are intramolecular force-field kcal/mol, meaningful for comparing
    conformers of this molecule only. No orbital gap, DFT, xTB, or stability claim.
    """
    if not isinstance(smiles, str) or not 1 <= len(smiles) <= 8000:
        raise ValueError("Supply a finite monomer SMILES of at most 8000 characters")
    count = _integer(conformers, "conformers", 1, 8)
    seed = _integer(seed, "seed", 0, 2**31 - 1)
    iterations = _integer(max_iterations, "max_iterations", 1, 1000)
    mol = Chem.MolFromSmiles(smiles)
    if mol is None or not 1 <= mol.GetNumAtoms() <= 250 or len(Chem.GetMolFrags(mol)) != 1:
        raise ValueError("Use one connected monomer with 1–250 atoms and valid valences")
    if any(atom.GetAtomicNum() == 0 for atom in mol.GetAtoms()):
        raise ValueError("Continuation/wildcard atoms are not finite monomers; supply a complete molecular graph")
    canonical = Chem.MolToSmiles(mol, isomericSmiles=True)
    descriptors = {"formula": rdMolDescriptors.CalcMolFormula(mol),
                   "molecular_weight_g_mol": Descriptors.MolWt(mol),
                   "exact_mass_da": rdMolDescriptors.CalcExactMolWt(mol),
                   "formal_charge": Chem.GetFormalCharge(mol),
                   "heavy_atoms": mol.GetNumHeavyAtoms(), "rings": rdMolDescriptors.CalcNumRings(mol),
                   "aromatic_rings": rdMolDescriptors.CalcNumAromaticRings(mol),
                   "rotatable_bonds": Lipinski.NumRotatableBonds(mol),
                   "hydrogen_bond_donors": Lipinski.NumHDonors(mol),
                   "hydrogen_bond_acceptors": Lipinski.NumHAcceptors(mol),
                   "tpsa_angstrom2": rdMolDescriptors.CalcTPSA(mol),
                   "crippen_logp": Crippen.MolLogP(mol),
                   "fraction_csp3": rdMolDescriptors.CalcFractionCSP3(mol)}
    embedded = Chem.AddHs(mol)
    params = AllChem.ETKDGv3()
    params.randomSeed = seed
    params.numThreads = 1
    params.maxIterations = 500
    params.enforceChirality = True
    params.pruneRmsThresh = 0.35
    if hasattr(params, "timeout"):
        params.timeout = 12
    ids = list(AllChem.EmbedMultipleConfs(embedded, numConfs=count, params=params))
    if not ids:
        raise ValueError("ETKDGv3 could not embed this monomer within the bounded search")
    forcefield = None
    if AllChem.MMFFHasAllMoleculeParams(embedded):
        forcefield = "MMFF94s"
        results = AllChem.MMFFOptimizeMoleculeConfs(embedded, numThreads=1, maxIters=iterations, mmffVariant=forcefield)
    elif AllChem.UFFHasAllMoleculeParams(embedded):
        forcefield = "UFF"
        results = AllChem.UFFOptimizeMoleculeConfs(embedded, numThreads=1, maxIters=iterations)
    else:
        results = [(None, None)] * len(ids)
    records = [{"id": int(cid), "converged": status == 0, "optimizer_status": status,
                "energy_kcal_mol": float(energy) if status in (0, 1) and energy is not None and math.isfinite(energy) else None}
               for cid, (status, energy) in zip(ids, results)]
    finite = [row for row in records if row["energy_kcal_mol"] is not None]
    converged = [row for row in finite if row["converged"]]
    selected = min(converged or finite, key=lambda row: row["energy_kcal_mol"]) if finite else records[0]
    cid = selected["id"]
    xyz = embedded.GetConformer(cid).GetPositions()
    charges = None
    try:
        rdPartialCharges.ComputeGasteigerCharges(embedded, nIter=12, throwOnParamFailure=True)
        candidate = [float(atom.GetProp("_GasteigerCharge")) for atom in embedded.GetAtoms()]
        if all(math.isfinite(charge) for charge in candidate):
            charges = candidate
    except (ValueError, RuntimeError):
        pass
    atoms = [{"element": atom.GetSymbol(), "x": float(pos[0]), "y": float(pos[1]), "z": float(pos[2]),
              "formal_charge": atom.GetFormalCharge(), "isotope": atom.GetIsotope(),
              "gasteiger_charge_e": charges[atom.GetIdx()] if charges is not None else None}
             for atom, pos in zip(embedded.GetAtoms(), xyz)]
    bonds = [{"a": b.GetBeginAtomIdx(), "b": b.GetEndAtomIdx(), "order": b.GetBondTypeAsDouble()}
             for b in embedded.GetBonds()]
    notes = ["ETKDG conformer search is finite and does not establish the global minimum.",
             "Force-field energies compare conformers of this same monomer only; these are not quantum energies or framework stability.",
             "No xTB calculation was run; RDKit descriptors and molecular mechanics are the bundled low-cost method."]
    if not selected["converged"]:
        notes.append("Selected geometry is not force-field converged; inspect it before further calculations.")
    notes.append("Gasteiger charges are empirical graph-based partial charges, not a quantum electron density or oxidation-state assignment.")
    if charges is None:
        notes.append("Gasteiger charges omitted because parameters or finite values were unavailable.")
    if any(flag == "?" for _, flag in Chem.FindMolChiralCenters(mol, includeUnassigned=True)):
        notes.append("Input contains unspecified stereocentres; the conformer is one realization, not a stereochemical assignment.")
    embedded.SetProp("_Name", "COFdraw monomer conformer")
    sdf = Chem.MolToMolBlock(embedded, confId=cid) + "\n>  <canonical_smiles>\n" + canonical + "\n\n$$$$\n"
    xyz_text = str(len(atoms)) + "\nCOFdraw ETKDGv3; angstrom; " + str(forcefield or "unoptimized") + "\n"
    xyz_text += "".join(f"{a['element']} {a['x']:.8f} {a['y']:.8f} {a['z']:.8f}\n" for a in atoms)
    return {"kind": "finite_monomer_conformer", "canonical_smiles": canonical, "descriptors": descriptors,
            "method": {"embedding": "ETKDGv3", "force_field": forcefield, "rdkit_version": rdBase.rdkitVersion,
                       "seed": seed, "requested_conformers": count, "embedded_conformers": len(ids), "max_iterations": iterations},
            "conformers": records, "selected_conformer": selected, "viewer": {"atoms": atoms, "bonds": bonds,
            "periodic": False, "bounds": _bounds(xyz), "units": "angstrom"}, "exports": {"xyz": xyz_text, "sdf": sdf},
            "partial_charges": {"method": "Gasteiger, 12 iterations, explicit hydrogens", "units": "elementary_charge",
                                "available": charges is not None, "values": charges,
                                "sum": float(sum(charges)) if charges is not None else None},
            "warnings": notes}


def cell_from_parameters(a, b, c, alpha=90, beta=90, gamma=90):
    lengths = [_finite(x, name, 0.5, 500) for x, name in zip((a, b, c), ("a", "b", "c"))]
    angles = [math.radians(_finite(x, name, 10, 170)) for x, name in zip((alpha, beta, gamma), ("alpha", "beta", "gamma"))]
    a, b, c = lengths
    ca, cb, cg = [math.cos(x) for x in angles]
    sg = math.sin(angles[2])
    z2 = 1 - cb * cb - ((ca - cb * cg) / sg) ** 2
    if z2 <= 1e-8:
        raise ValueError("Cell angles do not define a nondegenerate 3D cell")
    return [[a, 0.0, 0.0], [b * cg, b * sg, 0.0], [c * cb, c * (ca - cb * cg) / sg, c * math.sqrt(z2)]]


def normalize_structure(structure):
    """Validate explicitly supplied full-cell atoms. Occupancies must be exactly 1."""
    if not isinstance(structure, dict) or structure.get("periodic", True) is not True:
        raise ValueError("A periodic structure with an explicit 3D cell is required")
    cell = np.asarray(structure.get("cell"), dtype=float)
    if cell.shape != (3, 3) or not np.isfinite(cell).all():
        raise ValueError("Cell must contain three finite Cartesian lattice vectors in angstrom")
    volume = float(np.linalg.det(cell))
    if volume <= 0.1 or volume > 1.25e8 or np.linalg.cond(cell) > 100:
        raise ValueError("Cell must be right-handed, nondegenerate, and reasonably conditioned")
    lengths = np.linalg.norm(cell, axis=1)
    if min(lengths) < .5 or max(lengths) > 500:
        raise ValueError("Cell vector lengths must be between 0.5 and 500 angstrom")
    source_atoms = structure.get("atoms")
    if not isinstance(source_atoms, list) or not 1 <= len(source_atoms) <= MAX_ATOMS:
        raise ValueError(f"Provide 1–{MAX_ATOMS} atoms in the complete unit cell")
    inverse = np.linalg.inv(cell)
    atoms = []
    seen = set()
    for index, atom in enumerate(source_atoms):
        if not isinstance(atom, dict):
            raise ValueError(f"Atom {index + 1} must be an object with an element and coordinates")
        element = atom.get("element")
        if element not in ELEMENTS:
            raise ValueError(f"Unsupported or dummy element at atom {index + 1}")
        if abs(_finite(atom.get("occupancy", 1), "occupancy", 0, 1) - 1) > 1e-6:
            raise ValueError("Partial occupancies/disorder require an explicitly resolved structure")
        if "fractional" in atom:
            fractional = np.asarray(atom["fractional"], dtype=float)
        else:
            position = atom.get("position", [atom.get("x"), atom.get("y"), atom.get("z")])
            fractional = np.asarray(position, dtype=float) @ inverse
        if fractional.shape != (3,) or not np.isfinite(fractional).all() or np.max(np.abs(fractional)) > 1e6:
            raise ValueError("Every atom needs three finite fractional or Cartesian coordinates")
        fractional %= 1.0
        fractional[np.isclose(fractional, 1, atol=1e-10, rtol=0)] = 0
        key = tuple(np.round(fractional, 8))
        if key in seen:
            raise ValueError("Coincident atomic sites remain after symmetry expansion; resolve the structure")
        seen.add(key)
        atoms.append({"element": element, "fractional": fractional.tolist(), "occupancy": 1.0})
    provenance = structure.get("provenance", {})
    return {"cell": cell.tolist(), "atoms": atoms, "periodic": True, "units": "angstrom",
            "cell_volume_angstrom3": volume, "provenance": provenance if isinstance(provenance, dict) else {}}


def _cif_number(text):
    text = re.sub(r"\(\d+\)(?=([eE][+-]?\d+)?$)", "", str(text))
    try:
        result = float(text)
    except ValueError:
        raise ValueError("CIF contains a missing or unsupported numeric value") from None
    if not math.isfinite(result):
        raise ValueError("CIF values must be finite")
    return result


def _symmetry_operation(expression):
    """Parse crystallographic affine expressions; never use eval or exec."""
    pieces = str(expression).lower().replace(" ", "").split(",")
    if len(pieces) != 3 or len(expression) > 120:
        raise ValueError("Unsupported CIF symmetry operation")
    matrix = np.zeros((3, 3)); offset = np.zeros(3)
    for axis, piece in enumerate(pieces):
        terms = re.findall(r"[+-]?[^+-]+", piece)
        if "".join(terms) != piece:
            raise ValueError("Unsupported CIF symmetry expression")
        for term in terms:
            sign = -1 if term.startswith("-") else 1
            body = term.lstrip("+-")
            if body in ("x", "y", "z"):
                matrix[axis, "xyz".index(body)] += sign
            elif re.fullmatch(r"\d+(?:/\d+|\.\d+)?", body):
                try:
                    offset[axis] += sign * float(Fraction(body))
                except (ValueError, ZeroDivisionError):
                    raise ValueError("Invalid CIF symmetry translation") from None
            else:
                raise ValueError("Only affine x/y/z plus numeric translations are supported in CIF symmetry")
    if not math.isclose(abs(float(np.linalg.det(matrix))), 1, abs_tol=1e-8):
        raise ValueError("CIF symmetry operation is not invertible with determinant ±1")
    return matrix, offset


def import_cif(text):
    """Read one conventional core-CIF block, expand supplied affine operations.

    Space-group lookup is intentionally not guessed: non-P1 CIF requires an
    explicit symmetry-operation loop. Full occupancy and fractional sites only.
    """
    if not isinstance(text, str) or not 1 <= len(text) <= 2_000_000:
        raise ValueError("CIF must be text no larger than 2 MB")
    # CIF multiline text fields can contain tokens that look like data names.
    lines, multiline, block = [], False, []
    for line in text.splitlines():
        if line.startswith(";"):
            if multiline:
                lines.append("'COFdraw_multiline_text_omitted'"); block = []
            multiline = not multiline
        elif multiline:
            block.append(line)
        else:
            lines.append(line)
    if multiline:
        raise ValueError("Unterminated CIF multiline text field")
    lexer = shlex.shlex("\n".join(lines), posix=True)
    lexer.whitespace_split = True; lexer.commenters = "#"
    try:
        tokens = list(lexer)
    except ValueError as exc:
        raise ValueError("Malformed CIF quoting") from exc
    if len(tokens) > 180000:
        raise ValueError("CIF contains too many tokens")
    if sum(token.lower().startswith("data_") for token in tokens) != 1:
        raise ValueError("Import one CIF data block at a time")
    scalar, loops = {}, []
    i = 0
    while i < len(tokens):
        token = tokens[i].lower()
        if token == "loop_":
            i += 1; headers = []
            while i < len(tokens) and tokens[i].startswith("_"):
                headers.append(tokens[i].lower()); i += 1
            if not headers or len(set(headers)) != len(headers):
                raise ValueError("CIF loop needs distinct column names")
            values = []
            while i < len(tokens) and not (tokens[i].startswith("_") or tokens[i].lower() in ("loop_", "stop_") or tokens[i].lower().startswith(("data_", "save_"))):
                values.append(tokens[i]); i += 1
            if len(values) % len(headers):
                raise ValueError("CIF loop has an incomplete row")
            loops.append([dict(zip(headers, values[j:j + len(headers)])) for j in range(0, len(values), len(headers))])
        elif token.startswith("_"):
            if i + 1 >= len(tokens) or token in scalar:
                raise ValueError("Missing or duplicated CIF scalar")
            scalar[token] = tokens[i + 1]; i += 2
        else:
            i += 1
    try:
        values = [_cif_number(scalar["_cell_" + name]) for name in ("length_a", "length_b", "length_c", "angle_alpha", "angle_beta", "angle_gamma")]
    except KeyError:
        raise ValueError("CIF requires all six cell lengths/angles") from None
    cell = cell_from_parameters(*values)
    operations = []
    for rows in loops:
        for row in rows:
            for key in ("_space_group_symop_operation_xyz", "_symmetry_equiv_pos_as_xyz", "_space_group_symop.operation_xyz"):
                if key in row:
                    operations.append(_symmetry_operation(row[key]))
    for key in ("_space_group_symop_operation_xyz", "_symmetry_equiv_pos_as_xyz"):
        if key in scalar:
            operations.append(_symmetry_operation(scalar[key]))
    if not operations:
        number = scalar.get("_space_group_it_number", scalar.get("_symmetry_int_tables_number"))
        symbol = scalar.get("_space_group_name_h-m_alt", scalar.get("_symmetry_space_group_name_h-m", "")).replace(" ", "")
        if number not in (None, "1") or (symbol and symbol.lower() != "p1"):
            raise ValueError("Non-P1 CIF requires explicit symmetry operations; export a fully expanded P1 CIF")
        if number is None and not symbol:
            raise ValueError("CIF must declare P1 or provide explicit symmetry operations; asymmetric-unit expansion cannot be guessed")
        operations = [(np.eye(3), np.zeros(3))]
    if len(operations) > 192:
        raise ValueError("At most 192 explicit CIF symmetry operations are supported")
    metric = np.asarray(cell) @ np.asarray(cell).T
    for matrix, _ in operations:
        if not np.allclose(matrix.T @ metric @ matrix, metric, atol=1e-4, rtol=1e-4):
            raise ValueError("Supplied CIF symmetry is inconsistent with the cell metric")
    def operation_key(matrix, translation):
        shift = np.asarray(translation) % 1
        shift[np.isclose(shift, 1, atol=1e-8, rtol=0)] = 0
        return tuple(np.round(matrix, 8).ravel()) + tuple(np.round(shift, 8))
    keys = {operation_key(matrix, translation) for matrix, translation in operations}
    if operation_key(np.eye(3), np.zeros(3)) not in keys:
        raise ValueError("Explicit CIF symmetry operations must include the identity operation")
    for matrix_a, shift_a in operations:
        for matrix_b, shift_b in operations:
            if operation_key(matrix_a @ matrix_b, matrix_a @ shift_b + shift_a) not in keys:
                raise ValueError("CIF symmetry operations are incomplete or not closed under composition")
    atoms, seen = [], {}
    site_count = 0
    for rows in loops:
        if not rows or "_atom_site_fract_x" not in rows[0]:
            continue
        for row in rows:
            site_count += 1
            occupancy = _cif_number(row.get("_atom_site_occupancy", "1"))
            if abs(occupancy - 1) > 1e-6:
                raise ValueError("CIF contains partial occupancy; resolve disorder before geometric screening")
            raw = row.get("_atom_site_type_symbol", row.get("_atom_site_label", ""))
            match = re.match(r"([A-Z][a-z]?)", raw)
            element = match.group(1) if match else ""
            if element not in ELEMENTS or (match and len(raw) > match.end() and raw[match.end()].isalpha()):
                raise ValueError("Every CIF atom must have an unambiguous supported chemical element")
            try:
                fractional = np.array([_cif_number(row["_atom_site_fract_" + axis]) for axis in "xyz"])
            except KeyError:
                raise ValueError("CIF atom loop requires fractional x, y and z") from None
            for matrix, translation in operations:
                pos = (matrix @ fractional + translation) % 1
                pos[np.isclose(pos, 1, atol=1e-8, rtol=0)] = 0
                key = tuple(np.round(pos, 7))
                if key in seen:
                    if seen[key] != element:
                        raise ValueError("Coincident different elements require resolved disorder")
                    continue
                seen[key] = element
                atoms.append({"element": element, "fractional": pos.tolist()})
                if len(atoms) > MAX_ATOMS:
                    raise ValueError(f"Expanded CIF exceeds {MAX_ATOMS} atoms per cell")
    if not atoms:
        raise ValueError("No supported fractional atom-site loop found in CIF")
    return normalize_structure({"cell": cell, "atoms": atoms, "provenance": {
        "kind": "imported_cif", "input_sha256": hashlib.sha256(text.encode("utf-8")).hexdigest(),
        "asymmetric_sites": site_count, "supplied_symmetry_operations": len(operations),
        "expanded_sites": len(atoms), "identity": "user-supplied structure; not independently validated"}})


def _cell_parameters(cell):
    cell = np.asarray(cell, dtype=float)
    lengths = np.linalg.norm(cell, axis=1)
    angles = [math.degrees(math.acos(float(np.clip(np.dot(cell[a], cell[b]) / (lengths[a] * lengths[b]), -1, 1))))
              for a, b in ((1, 2), (0, 2), (0, 1))]
    return [*map(float, lengths), *angles]


def structure_exports(structure, repeats=(1, 1, 1)):
    normalized = normalize_structure(structure)
    if len(repeats) != 3:
        raise ValueError("Repeat counts need three integers")
    repeats = [_integer(n, "repeat", 1, 6) for n in repeats]
    count = len(normalized["atoms"]) * math.prod(repeats)
    if count > 12000:
        raise ValueError("Displayed/exported supercells are limited to 12000 atoms")
    # CIF/CRYST1 encode lengths and angles in conventional cell orientation.
    # Rotate coordinates consistently rather than pairing rotated Cartesian
    # positions with a differently oriented PDB unit cell.
    cell = np.asarray(cell_from_parameters(*_cell_parameters(normalized["cell"])))
    supercell = cell * np.array(repeats)[:, None]
    atoms, fractions = [], []
    for shift in itertools.product(*(range(n) for n in repeats)):
        for source in normalized["atoms"]:
            fractional = (np.asarray(source["fractional"]) + shift) / repeats
            xyz = fractional @ supercell
            atoms.append({"element": source["element"], "x": float(xyz[0]), "y": float(xyz[1]), "z": float(xyz[2])})
            fractions.append(fractional)
    params = _cell_parameters(supercell)
    cif = "data_COFdraw\n# Explicit full-cell geometry; no new experimental identity assigned.\n"
    cif += "".join(f"_cell_{name} {number:.9f}\n" for name, number in zip(("length_a", "length_b", "length_c", "angle_alpha", "angle_beta", "angle_gamma"), params))
    cif += "_space_group_name_H-M_alt 'P 1'\n_space_group_IT_number 1\nloop_\n_space_group_symop_operation_xyz\n'x,y,z'\n"
    cif += "loop_\n_atom_site_label\n_atom_site_type_symbol\n_atom_site_fract_x\n_atom_site_fract_y\n_atom_site_fract_z\n_atom_site_occupancy\n"
    for index, (atom, fractional) in enumerate(zip(atoms, fractions), 1):
        cif += f"{atom['element']}{index} {atom['element']} " + " ".join(f"{x:.10f}" for x in fractional) + " 1\n"
    lattice = " ".join(f"{x:.9f}" for x in supercell.ravel())
    xyz = f'{count}\nLattice="{lattice}" Properties=species:S:1:pos:R:3 pbc="T T T"\n'
    xyz += "".join(f"{atom['element']} {atom['x']:.8f} {atom['y']:.8f} {atom['z']:.8f}\n" for atom in atoms)
    a, b, c, alpha, beta, gamma = params
    pdb = "REMARK COFDRAW COORDINATES; BOND ORDERS AND FRAMEWORK IDENTITY NOT INFERRED\n"
    pdb += f"CRYST1{a:9.3f}{b:9.3f}{c:9.3f}{alpha:7.2f}{beta:7.2f}{gamma:7.2f} P 1           1\n"
    for index, atom in enumerate(atoms, 1):
        pdb += f"HETATM{index:5d} {atom['element']:>2s}   COF A   1    {atom['x']:8.3f}{atom['y']:8.3f}{atom['z']:8.3f}  1.00  0.00          {atom['element']:>2s}\n"
    pdb += "END\n"
    positions = [[atom[k] for k in ("x", "y", "z")] for atom in atoms]
    exports = {"cif": cif, "xyz": xyz, "pdb": pdb}
    warnings = ["Imported/proposed coordinates are not automatically relaxed or experimentally validated.",
                "CIF/extended XYZ retain the lattice; plain XYZ readers may ignore the extended lattice header.",
                "PDB is a coordinate interchange preview and does not encode framework bond orders."]
    if np.min(positions) < -999.999 or np.max(positions) > 9999.999:
        del exports["pdb"]
        warnings.append("PDB omitted because Cartesian coordinates exceed its fixed-column field width; use CIF or extended XYZ.")
    return {"kind": "explicit_periodic_structure", "structure": normalized, "repeats": repeats,
            "viewer": {"atoms": atoms, "bonds": [], "cell": supercell.tolist(), "periodic": True,
                       "bounds": _bounds(positions), "units": "angstrom"}, "exports": exports, "warnings": warnings}


class _PeriodicDistances:
    """Exact finite minimum-image search with a certified per-axis image bound.

    Orthogonal metric components separate; this avoids exploding c-axis image
    searches for long 2D-COF a/b vectors and short perpendicular layer spacing.
    Pathological unreduced cells fail explicitly instead of giving wrong pores.
    """
    def __init__(self, cell, fractional, radii):
        self.cell = np.asarray(cell); self.fractional = np.asarray(fractional); self.radii = np.asarray(radii)
        gram = self.cell @ self.cell.T
        remaining = set(range(3)); self.components = []
        while remaining:
            axes = {remaining.pop()}; todo = list(axes)
            while todo:
                i = todo.pop()
                linked = {j for j in remaining if abs(gram[i, j]) > 1e-10 * math.sqrt(gram[i, i] * gram[j, j])}
                remaining -= linked; axes |= linked; todo.extend(linked)
            axes = sorted(axes); metric = gram[np.ix_(axes, axes)]
            if len(axes) == 1:
                offsets = [np.zeros(1)]
            else:
                corners = np.array(list(itertools.product((-.5, .5), repeat=len(axes))))
                cover = math.sqrt(float(np.max(np.einsum("ni,ij,nj->n", corners, metric, corners))))
                bounds = np.floor(cover * np.sqrt(np.diag(np.linalg.inv(metric))) + .5 + 1e-10).astype(int)
                if math.prod(2 * bounds + 1) > 343:
                    raise ValueError("This highly skewed cell needs lattice reduction before bounded pore analysis")
                offsets = [np.asarray(x) for x in itertools.product(*(range(-int(b), int(b) + 1) for b in bounds))]
            self.components.append((axes, metric, offsets))

    def clearance(self, points):
        values = np.empty(len(points))
        for start in range(0, len(points), 128):
            delta = np.asarray(points[start:start + 128])[:, None, :] - self.fractional[None, :, :]
            delta -= np.rint(delta)
            square = np.zeros(delta.shape[:2])
            for axes, metric, offsets in self.components:
                d = delta[:, :, axes]
                smallest = np.full(delta.shape[:2], np.inf)
                for shift in offsets:
                    shifted = d - shift
                    ds = np.einsum("...i,ij,...j->...", shifted, metric, shifted, optimize=False)
                    np.minimum(smallest, ds, out=smallest)
                square += smallest
            values[start:start + len(delta)] = np.min(np.sqrt(np.maximum(square, 0)) - self.radii, axis=1)
        return values


def _wilson(successes, total):
    p = successes / total; z = 1.95996398454; denominator = 1 + z * z / total
    middle = (p + z * z / (2 * total)) / denominator
    half = z * math.sqrt(p * (1 - p) / total + z * z / (4 * total * total)) / denominator
    return [max(0.0, middle - half), min(1.0, middle + half)]


def _percolation(free):
    n = free.shape[0]; labels = np.full(free.shape, -1, dtype=int); positions = {}; components = []
    for start in zip(*np.where(free)):
        if labels[start] >= 0:
            continue
        identifier = len(components); queue = deque([start]); labels[start] = identifier
        positions[start] = np.array(start, dtype=int); size = 0; winding = np.zeros(3, dtype=bool)
        while queue:
            current = queue.popleft(); size += 1
            for axis in range(3):
                for direction in (-1, 1):
                    step = np.zeros(3, dtype=int); step[axis] = direction
                    neighbor = tuple((np.array(current) + step) % n)
                    if not free[neighbor]:
                        continue
                    expected = positions[current] + step
                    if labels[neighbor] < 0:
                        labels[neighbor] = identifier; positions[neighbor] = expected; queue.append(neighbor)
                    elif labels[neighbor] == identifier:
                        winding |= expected != positions[neighbor]
        components.append({"voxels": size, "percolating_axes": ["abc"[i] for i in range(3) if winding[i]]})
    connected = sum(c["voxels"] for c in components if c["percolating_axes"])
    return {"grid_size": n, "grid_points": n ** 3, "probe_free_fraction": float(free.mean()),
            "percolating_probe_center_fraction": connected / n ** 3, "components": components[:100],
            "component_count": len(components), "components_truncated": len(components) > 100,
            "method": "periodic six-neighbor voxel graph; winding cycles identify spanning channels",
            "limitation": "Coarse grid estimate; narrow or diagonal channels may be missed. Repeat at finer resolution or use Zeo++ for screening-grade accessibility."}


def pore_analysis(structure, *, probe_radius=1.82, samples=12000, surface_samples=96, seed=42, grid_size=16):
    """Seeded hard-sphere geometry of a supplied periodic full-cell structure.

    Probe-center free volume differs from zero-probe void volume. Local exposed
    sphere area includes enclosed cavities; percolation is separately estimated.
    No adsorption/BET, pore-limiting diameter, or helium Widom volume is inferred.
    """
    normalized = normalize_structure(structure)
    radius = _finite(probe_radius, "probe_radius", 0, 5)
    requested = _integer(samples, "samples", 512, 50000)
    per_atom = _integer(surface_samples, "surface_samples", 16, 512)
    seed = _integer(seed, "seed", 0, 2**31 - 1)
    grid_requested = _integer(grid_size, "grid_size", 8, 28)
    atoms = normalized["atoms"]; n = len(atoms)
    actual = min(requested, max(512, 10_000_000 // n))
    area_count = min(max(512, n * per_atom), max(512, 10_000_000 // n))
    grid = min(grid_requested, max(8, int((8_000_000 / n) ** (1 / 3))))
    elements = [a["element"] for a in atoms]
    radii = np.array([PT.GetRvdw(ELEMENTS[e]) for e in elements])
    if not np.isfinite(radii).all() or (radii <= 0).any():
        raise ValueError("RDKit has no usable van der Waals radius for an element")
    mass = float(sum(PT.GetAtomicWeight(ELEMENTS[e]) for e in elements))
    volume = normalized["cell_volume_angstrom3"]
    density = mass / AVOGADRO / (volume * 1e-24)
    cell = np.asarray(normalized["cell"]); fractional = np.array([a["fractional"] for a in atoms])
    distance = _PeriodicDistances(cell, fractional, radii)
    image_searches = sum(len(offsets) for _, _, offsets in distance.components)
    # Keep unusual, non-orthogonal imported cells bounded too. Effective counts
    # are exposed rather than silently reporting the requested precision.
    image_budget_samples = max(512, 100_000_000 // (n * image_searches))
    actual = min(actual, image_budget_samples)
    area_count = min(area_count, image_budget_samples)
    grid = min(grid, max(8, int((80_000_000 / (n * image_searches)) ** (1 / 3))))
    rng = np.random.default_rng(seed)
    points = rng.random((actual, 3)); clearance = distance.clearance(points)
    void = clearance >= 0; free = clearance >= radius
    inflated = radii + radius
    sphere_areas = 4 * math.pi * inflated * inflated
    chosen = rng.choice(n, size=area_count, p=sphere_areas / sphere_areas.sum())
    directions = rng.normal(size=(area_count, 3)); directions /= np.linalg.norm(directions, axis=1)[:, None]
    cartesian = fractional[chosen] @ cell + directions * inflated[chosen, None]
    surface = cartesian @ np.linalg.inv(cell)
    exposed = distance.clearance(surface) >= radius - 1e-7
    area = float(sphere_areas.sum() * exposed.mean())
    area_interval = [float(sphere_areas.sum() * x) for x in _wilson(int(exposed.sum()), area_count)]
    grid_points = np.array(list(itertools.product(range(grid), repeat=3)), dtype=float)
    grid_clearance = distance.clearance((grid_points + .5) / grid).reshape((grid, grid, grid))
    connectivity = _percolation(grid_clearance >= radius)
    diameters = 2 * clearance[void]
    hist, edges = np.histogram(diameters, bins=12) if len(diameters) else (np.zeros(12, dtype=int), np.linspace(0, 1, 13))
    notes = ["Hard-sphere geometry of the supplied coordinates, not adsorption, BET, or a relaxed-crystal prediction.",
             "Local probe-exposed surface includes isolated cavities; transport accessibility is only the separate coarse-grid estimate.",
             "Largest sampled included-sphere diameter is a sampling lower bound, not a pore-limiting diameter or a rigorous pore-size distribution.",
             "Density and mass-normalized metrics use only the listed full-occupancy atoms and standard atomic weights.",
             "Random-sampling intervals exclude radius-model, structural, and grid-resolution uncertainty."]
    if "H" not in elements:
        notes.append("No hydrogen atoms are listed. Missing hydrogens, solvent, defects, or unresolved disorder can change geometry and density.")
    return {"kind": "periodic_hard_sphere_geometry", "structure_sha256": hashlib.sha256(json.dumps(normalized, sort_keys=True).encode()).hexdigest(),
            "method": {"probe_radius_angstrom": radius, "radii": "RDKit periodic-table van der Waals radii",
                       "radii_angstrom": {e: PT.GetRvdw(ELEMENTS[e]) for e in sorted(set(elements))},
                       "seed": seed, "requested_volume_samples": requested, "volume_samples": actual,
                       "surface_samples": area_count, "requested_surface_samples_per_atom": per_atom,
                       "periodic_image_searches_per_pair": image_searches,
                       "rdkit_version": rdBase.rdkitVersion, "distance": "bounded exact periodic minimum-image enumeration"},
            "atom_count": n, "cell_volume_angstrom3": volume, "cell_molar_mass_g_mol": mass, "density_g_cm3": density,
            "void_fraction": float(void.mean()), "void_fraction_95_interval": _wilson(int(void.sum()), actual),
            "geometric_void_volume_cm3_g": float(void.mean() / density),
            "probe_center_free_fraction": float(free.mean()), "probe_center_free_fraction_95_interval": _wilson(int(free.sum()), actual),
            "probe_center_free_volume_cm3_g": float(free.mean() / density),
            "local_probe_exposed_surface_angstrom2": area, "local_probe_exposed_surface_95_interval_angstrom2": area_interval,
            "local_probe_exposed_surface_m2_g": area * 1e-20 * AVOGADRO / mass,
            "largest_sampled_included_sphere_diameter_angstrom": float(max(0, 2 * max(float(clearance.max()), float(grid_clearance.max())))),
            "local_clearance_diameter_histogram": {"edges_angstrom": edges.tolist(), "counts": hist.tolist(),
                "definition": "Twice distance from sampled void points to nearest atomic van der Waals surface; not a standard PSD."},
            "connectivity": connectivity, "bet_surface_area_m2_g": None, "warnings": notes}


def fit_bet(relative_pressures, uptake_mmol_g, *, fit_min=.05, fit_max=.30, cross_section_nm2=.162):
    """Fit supplied adsorption isotherm points; never synthesize uptake from CIF.

    User-selected range, ordinary least squares of x/[n(1-x)]. Return basic
    consistency checks without claiming an automatic valid BET range selection.
    """
    x = np.asarray(relative_pressures, dtype=float); uptake = np.asarray(uptake_mmol_g, dtype=float)
    if x.ndim != 1 or uptake.shape != x.shape or not 3 <= len(x) <= 5000:
        raise ValueError("Supply 3–5000 paired relative pressures and mmol/g uptake values")
    if not np.isfinite(x).all() or not np.isfinite(uptake).all() or np.any((x <= 0) | (x >= 1)) or np.any(uptake <= 0):
        raise ValueError("Relative pressures must be between 0 and 1 and uptake must be positive")
    if np.any(np.diff(x) <= 0):
        raise ValueError("Relative pressures must be strictly increasing; supply one adsorption branch")
    low = _finite(fit_min, "fit_min", 0, 1); high = _finite(fit_max, "fit_max", 0, 1)
    if high <= low:
        raise ValueError("BET fit_max must exceed fit_min")
    cross = _finite(cross_section_nm2, "cross_section_nm2", .01, 1)
    selected = (x >= low) & (x <= high)
    if selected.sum() < 3:
        raise ValueError("The selected BET range must include at least three measured points")
    px, n = x[selected], uptake[selected]
    y = px / (n * (1 - px))
    slope, intercept = map(float, np.polyfit(px, y, 1))
    if intercept <= 0 or slope + intercept <= 0:
        raise ValueError("Selected range gives nonphysical BET intercept or monolayer capacity")
    capacity = 1 / (slope + intercept); constant = 1 + slope / intercept
    if constant <= 0:
        raise ValueError("Selected range gives nonpositive BET C")
    predicted = slope * px + intercept
    variation = float(np.sum((y - y.mean()) ** 2))
    r2 = 1 - float(np.sum((y - predicted) ** 2)) / variation if variation > 0 else None
    p_mono = 1 / (math.sqrt(constant) + 1)
    capacity_in_range = bool(n.min() <= capacity <= n.max())
    observed_pressure = float(np.interp(capacity, n, px)) if capacity_in_range and np.all(np.diff(n) > 0) else None
    consistency = abs(observed_pressure - p_mono) / p_mono if observed_pressure is not None else None
    checks = {"positive_C": constant > 0, "n_times_one_minus_p_increases": bool(np.all(np.diff(n * (1 - px)) > 0)),
              "monolayer_uptake_within_selected_data": capacity_in_range,
              "monolayer_relative_pressure_in_fit_range": bool(px[0] <= p_mono <= px[-1]),
              "monolayer_pressure_agreement_within_20_percent": consistency is not None and consistency <= .20}
    area = capacity * 1e-3 * AVOGADRO * cross * 1e-18
    return {"kind": "BET_fit_from_supplied_isotherm", "surface_area_m2_g": area,
            "monolayer_capacity_mmol_g": capacity, "C": constant, "slope": slope, "intercept": intercept,
            "r_squared": r2, "fit_range": [float(px[0]), float(px[-1])], "points_used": int(selected.sum()),
            "molecular_cross_section_nm2": cross, "monolayer_relative_pressure": p_mono,
            "interpolated_monolayer_relative_pressure": observed_pressure, "consistency_checks": checks,
            "all_basic_checks_pass": all(checks.values()), "plot": {"relative_pressure": px.tolist(), "bet_y": y.tolist(), "fitted_y": predicted.tolist()},
            "warnings": ["Range is user-selected; passing these basic checks is not a complete automatic BET/Rouquerol validation.",
                         "A high regression R² alone does not establish a physically valid BET area, especially for microporous COFs.",
                         "Cross-section 0.162 nm² is a conventional nitrogen value; match adsorptive, temperature, and analysis assumptions."]}


def patch_3d(recipe, *, layer_spacing=3.4, layers=3, stacking="AA", slip=None, bond_length=1.4):
    """Chemically scaled finite drawing layers; explicitly not a periodic crystal."""
    from .framework_drawing import build_patch
    layer_count = _integer(layers, "layers", 1, 6)
    spacing = _finite(layer_spacing, "layer_spacing", 1, 20)
    target_bond = _finite(bond_length, "bond_length", 1, 2)
    if stacking not in ("AA", "AB", "slipped"):
        raise ValueError("Stacking must be AA, AB, or slipped")
    shift = np.asarray([0, 0] if slip is None else slip, dtype=float)
    if shift.shape != (2,) or not np.isfinite(shift).all() or np.max(np.abs(shift)) > 20:
        raise ValueError("Layer slip must contain two finite angstrom offsets between -20 and 20")
    if stacking == "AB" and slip is None:
        raise ValueError("AB illustration requires an explicit two-coordinate slip; it is not inferred from pore shape")
    mol, audit = build_patch(recipe)
    if mol.GetNumAtoms() * layer_count > 12000:
        raise ValueError("Finite layer preview exceeds 12000 atom/continuation sites")
    base = mol.GetConformer().GetPositions()
    real_bonds = [b for b in mol.GetBonds() if b.GetBeginAtom().GetAtomicNum() and b.GetEndAtom().GetAtomicNum()]
    distances = [np.linalg.norm(base[b.GetBeginAtomIdx()] - base[b.GetEndAtomIdx()]) for b in real_bonds]
    if not distances or min(distances) <= 0:
        raise ValueError("Patch has invalid drawing bond lengths")
    scale = target_bond / float(np.median(distances))
    base = (base - base.mean(axis=0)) * scale
    atoms, bonds = [], []
    for layer in range(layer_count):
        offset = np.zeros(3)
        offset[:2] = shift * (layer % 2 if stacking == "AB" else layer) if stacking != "AA" else 0
        offset[2] = layer * spacing
        for atom, position in zip(mol.GetAtoms(), base + offset):
            atoms.append({"element": atom.GetSymbol(), "x": float(position[0]), "y": float(position[1]), "z": float(position[2]),
                          "continuation_site": atom.GetAtomicNum() == 0, "layer": layer})
        for bond in mol.GetBonds():
            bonds.append({"a": layer * mol.GetNumAtoms() + bond.GetBeginAtomIdx(), "b": layer * mol.GetNumAtoms() + bond.GetEndAtomIdx(),
                          "order": bond.GetBondTypeAsDouble()})
    xyz = f"{len(atoms)}\nCOFdraw finite illustrative layers; angstrom; X denotes continuation placeholder; not periodic or relaxed\n"
    xyz += "".join(f"{'X' if a['continuation_site'] else a['element']} {a['x']:.8f} {a['y']:.8f} {a['z']:.8f}\n" for a in atoms)
    positions = [[atom[k] for k in ("x", "y", "z")] for atom in atoms]
    return {"kind": "finite_illustrative_layer_stack", "periodic": False,
            "viewer": {"atoms": atoms, "bonds": bonds, "periodic": False, "bounds": _bounds(positions), "units": "angstrom"},
            "exports": {"xyz": xyz}, "parameters": {"layers": layer_count, "spacing_angstrom": spacing, "stacking": stacking,
                        "slip_angstrom": shift.tolist(), "median_bond_length_angstrom": target_bond},
            "source_audit": audit, "warnings": ["Finite illustrative layers made from scaled 2D drawing coordinates; not an optimized periodic crystal.",
                "No interlayer covalent bonds or missing hydrogens were invented. X sites preserve continuation placeholders.",
                "Do not use this preview as input for periodic pore, adsorption, or BET calculations."]}
