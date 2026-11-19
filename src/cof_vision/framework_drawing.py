"""Auditable COF linkage rules and chemical pore-fragment reconstruction.

COFex keys identify reference records, not complete periodic atom graphs. A
versioned drawing recipe supplies the missing site incidence and depiction
scope. Output is a finite, connected, port-labelled patch, never a whole-COF
SMILES, a crystal structure, or an experimentally established topology.
"""
from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path

import numpy as np
from rdkit import Chem
from rdkit.Chem import rdDepictor
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Geometry import Point3D
from .depiction_rings import symmetric_depiction_rings
from .cofdraw_depiction import compute_2d_coords

VERSION = "cof-chemical-patch-1"
TP = "O=Cc1c(O)c(C=O)c(O)c(C=O)c1O"
# Ordinary explicit H is removed by the SMILES reader. Retained isotope-labelled
# hydrogen neighbours need a separate atom-balanced reaction rule, so exclude
# them rather than silently deleting D/T during the condensation.
ARYL_HYDROXYLAMINE_SMARTS = "[c]-[NX3;H1;+0;!$(*-[#1])]-[OX2;H1;+0;!$(*-[#1])]"
LINKAGE_LIBRARY = {
    "imine": {
        "version": 1, "status": "implemented_scoped",
        "roles": ["aryl_aldehyde", "primary_arylamine"],
        "site_smarts": {"aryl_aldehyde": "[c][CH1](=[O])", "primary_arylamine": "[c][NX3;H2;+0]"},
        "join_bond_order": 2,
        "bond_changes": "Remove aldehyde O; connect aldehyde C to amine N with C=N. Condensation water is not part of the polymer patch.",
        "limits": "Aryl aldehydes and primary arylamines only; source ring incidence required; E/Z is not assigned by geometry.",
    },
    "beta_ketoenamine": {
        "version": 1, "status": "implemented_scoped",
        "roles": ["tp_trialdehyde", "primary_arylamine"],
        "site_smarts": {"tp_trialdehyde": "exact canonical Tp graph", "primary_arylamine": "[c][NX3;H2;+0]"},
        "join_bond_order": 1,
        "bond_changes": "Tp OH sites become C=O; aldehyde carbon becomes exocyclic =CH-NH-Ar; six-membered Tp ring bonds become single.",
        "limits": "Exact unsubstituted Tp plus primary arylamines; keto tautomer must be source supported. Mixed keto/enol and partial conversion require explicit recipes.",
    },
}

# These are explicit product-graph construction rules, not claims that any
# arbitrary pair will form a crystalline framework under a universal recipe.
LINKAGE_LIBRARY.update({
    "nitrone": {"version":1,"status":"implemented_scoped","roles":["aryl_aldehyde","aryl_hydroxylamine"],
                "site_smarts":{"aryl_aldehyde":"[c][CH1](=[O])","aryl_hydroxylamine":ARYL_HYDROXYLAMINE_SMARTS},
                "join_bond_order":2,
                "bond_changes":"Remove aldehyde O; connect aldehyde C to hydroxylamine N as Ar-CH=[N+]([O-])-Ar. Retain N-O; remove N-H and O-H from the patch.",
                "limits":"Neutral directly aryl-bound N(H)-O(H) only; protected, charged and N-substituted hydroxylamines or reactive D/T need explicit recipes. New C=N E/Z, crystallization, topology and synthesis conditions are not assigned.",
                "source_doi":"10.1002/anie.202307674",
                "source_url":"https://yaghi.berkeley.edu/pdfPublications/23NitroneCOF.pdf",
                "source_scope":"Scheme 1 supports the local product motif; this rule does not assign COF-115 identity or its kgm net."},
    "acylhydrazone": {"version":1,"status":"implemented_scoped","roles":["aryl_aldehyde","acylhydrazide"],"join_bond_order":2,
                      "bond_changes":"Remove aldehyde O; connect aldehyde C to terminal hydrazide N as C=N-NH-C(=O).",
                      "limits":"Terminal acylhydrazide NH2 only. Acyl NH is retained. E/Z and synthesis conditions require source evidence."},
    "keto_hydrazide": {"version":1,"status":"implemented_scoped","roles":["tp_trialdehyde","acylhydrazide"],"join_bond_order":1,
                      "bond_changes":"Tp keto-enamine tautomer with CH-NH-NH-C(=O) connectivity; hydrazide carbonyl retained.",
                      "limits":"Exact Tp and terminal acylhydrazides; keto tautomer is a chosen product-state assumption."},
    "azine": {"version":1,"status":"implemented_scoped","roles":["aryl_aldehyde","hydrazine"],"join_bond_order":2,
              "bond_changes":"Hydrazine terminal nitrogens form C=N-N=C connections.",
              "limits":"Exact neutral hydrazine only; reagents/counterions, E/Z and synthetic feasibility are not inferred."},
    "secondary_amine": {"version":1,"status":"implemented_scoped","roles":["aryl_aldehyde","primary_arylamine"],"join_bond_order":1,
                        "bond_changes":"Reduced imine product connectivity Ar-CH2-NH-Ar.",
                        "limits":"Explicit reduced-imine product design; not a direct aldehyde/amine reaction prediction. Reduction evidence and conditions are required for a literature assignment."},
})


def canonical(value):
    mol = Chem.MolFromSmiles(value)
    if mol is None:
        raise ValueError("Invalid molecular SMILES")
    return Chem.MolToSmiles(mol, isomericSmiles=True)


def linkage_fragment(smiles: str, linkage: str, role: str) -> dict:
    """Convert a supported monomer to a reacted fragment with mapped ports.

    Port numbers are assigned by the input atom order and returned explicitly;
    they must not be interpreted as a canonical numbering across input SMILES.
    This is a conditional construction rule, not a reaction feasibility claim.
    """
    if linkage not in LINKAGE_LIBRARY or (role not in LINKAGE_LIBRARY[linkage]["roles"] and (linkage,role)!=("imine","self_imine")):
        raise ValueError("Unsupported linkage/role combination")
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError("Invalid monomer SMILES")
    if any(a.GetAtomMapNum() for a in mol.GetAtoms()):
        raise ValueError("Monomer atom maps must be absent; ports use their own map namespace")
    if role == "hydrazine" and canonical(smiles) != canonical("NN"):
        raise ValueError("The azine construction rule requires exact neutral hydrazine")
    if len(Chem.GetMolFrags(mol)) != 1 or any(a.GetAtomicNum()==0 for a in mol.GetAtoms()):
        raise ValueError("Choose one connected precursor with no attachment placeholders")
    rw = Chem.RWMol(mol)
    sites, removed = [], []
    if role == "primary_arylamine":
        sites = [hit[1] for hit in mol.GetSubstructMatches(Chem.MolFromSmarts("[c][NX3;H2;+0]"))]
    elif role == "aryl_hydroxylamine":
        for _, nitrogen, oxygen in mol.GetSubstructMatches(Chem.MolFromSmarts(ARYL_HYDROXYLAMINE_SMARTS)):
            sites.append(nitrogen)
            n_atom, o_atom = rw.GetAtomWithIdx(nitrogen), rw.GetAtomWithIdx(oxygen)
            n_atom.SetFormalCharge(1)
            n_atom.SetChiralTag(Chem.ChiralType.CHI_UNSPECIFIED)
            o_atom.SetFormalCharge(-1)
            o_atom.SetNumExplicitHs(0)
            o_atom.SetNoImplicit(True)
    elif role == "acylhydrazide":
        sites = [hit[0] for hit in mol.GetSubstructMatches(Chem.MolFromSmarts("[NX3;H2]-[NX3;H1]-[CX3](=[OX1])"))]
    elif role == "hydrazine":
        sites = [a.GetIdx() for a in mol.GetAtoms()]
    elif role in ("aryl_aldehyde","self_imine"):
        for _, carbon, oxygen in mol.GetSubstructMatches(Chem.MolFromSmarts("[c][CH1](=[O])")):
            sites.append(carbon); removed.append(oxygen)
        if role == "self_imine":
            aldehydes = set(sites)
            amines = [n for _, n in mol.GetSubstructMatches(Chem.MolFromSmarts("[c][NH2;+0]"))]
            if len(aldehydes) != len(amines) or len(amines) not in (1,2):
                raise ValueError("Choose a balanced aminoaldehyde with one or two of each free group")
            sites.extend(amines)
    elif role == "tp_trialdehyde":
        if canonical(smiles) != canonical(TP):
            raise ValueError("This ketoenamine rule supports exact Tp only")
        ring = mol.GetRingInfo().AtomRings()[0]
        for atom_id in ring:
            rw.GetAtomWithIdx(atom_id).SetIsAromatic(False)
        for bond in rw.GetBonds():
            if bond.GetBeginAtomIdx() in ring and bond.GetEndAtomIdx() in ring:
                bond.SetIsAromatic(False); bond.SetBondType(Chem.BondType.SINGLE)
        for aryl, carbon, oxygen in mol.GetSubstructMatches(Chem.MolFromSmarts("[c][CH1](=[O])")):
            sites.append(carbon); removed.append(oxygen)
            rw.GetBondBetweenAtoms(aryl, carbon).SetBondType(Chem.BondType.DOUBLE)
        for aryl, oxygen in mol.GetSubstructMatches(Chem.MolFromSmarts("[c][OH1]")):
            rw.GetBondBetweenAtoms(aryl, oxygen).SetBondType(Chem.BondType.DOUBLE)
            rw.GetAtomWithIdx(oxygen).SetNumExplicitHs(0)
    else:
        raise ValueError("No construction rule for this reactive role")
    if not sites:
        raise ValueError("No supported reactive sites")
    order = LINKAGE_LIBRARY[linkage]["join_bond_order"]
    site_records = []
    for port, atom_id in enumerate(sites, 1):
        # Explicit H counts from bracketed SMILES must respond to the new bond.
        # Otherwise [NH2] can retain two H after a double bond is added.
        rw.GetAtomWithIdx(atom_id).SetNumExplicitHs(0)
        rw.GetAtomWithIdx(atom_id).SetNoImplicit(False)
        dummy = Chem.Atom(0); dummy.SetAtomMapNum(port)
        index = rw.AddAtom(dummy)
        rw.AddBond(atom_id, index, Chem.BondType.DOUBLE if order == 2 else Chem.BondType.SINGLE)
        site_records.append({"port": port, "monomer_atom_index": atom_id, "join_bond_order": order})
        if role == "self_imine":
            site_records[-1]["site_role"] = "aldehyde" if atom_id in aldehydes else "amine"
    for atom_id in sorted(set(removed), reverse=True):
        rw.RemoveAtom(atom_id)
    result = rw.GetMol()
    Chem.SanitizeMol(result)
    return {"monomer_smiles": canonical(smiles), "fragment_smiles": Chem.MolToSmiles(result),
            "linkage": linkage, "role": role, "reactive_sites": site_records,
            "connectivity": len(sites), "library_version": VERSION}


def _ports(mol):
    out = {}
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 0:
            key = atom.GetAtomMapNum()
            if key < 1 or key in out or atom.GetDegree() != 1:
                raise ValueError("Each dummy port needs a unique positive map and exactly one neighbor")
            neighbor = atom.GetNeighbors()[0].GetIdx()
            out[key] = (atom.GetIdx(), neighbor)
    return out


def _fragment_layout(mol, incoming, outgoing, position, previous, following, flip_across_attachment_axis=False, endpoint_alignment=False, porphyrin_template=False):
    depict = compute_2d_coords if porphyrin_template else rdDepictor.Compute2DCoords
    depict(mol, canonOrient=False, forceRDKit=True)
    coords = np.array(mol.GetConformer().GetPositions())[:, :2]
    real = [a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum()]
    center = coords[real].mean(axis=0)
    coords -= center
    ports = _ports(mol)
    endpoint = 1 if endpoint_alignment else 0
    left, right = coords[ports[incoming][endpoint]], coords[ports[outgoing][endpoint]]
    target_left, target_right = np.asarray(previous)-position, np.asarray(following)-position
    # A planar reflection is safe only when no stereochemistry is specified.
    local = np.array([left/np.linalg.norm(left), right/np.linalg.norm(right)])
    target = np.array([target_left/np.linalg.norm(target_left), target_right/np.linalg.norm(target_right)])
    u, _, vt = np.linalg.svd(local.T @ target)
    stereo = any(a.GetChiralTag() != Chem.ChiralType.CHI_UNSPECIFIED for a in mol.GetAtoms()) or any(b.GetStereo() != Chem.BondStereo.STEREONONE for b in mol.GetBonds())
    correction = np.eye(2)
    if stereo:
        correction[-1, -1] = np.linalg.det(u @ vt)
    rotation = u @ correction @ vt
    if endpoint_alignment and np.linalg.norm(target[0]+target[1]) < 1e-8:
        # A linear linker has a rank-one fit. SVD's arbitrary null-space sign
        # used to reflect equivalent linkers differently around the same pore.
        # Its directed attachment chord fixes a unique proper rotation.
        chord = right-left
        desired = target_right-target_left
        theta = math.atan2(desired[1], desired[0])-math.atan2(chord[1], chord[0])
        rotation = np.array([[math.cos(theta), math.sin(theta)], [-math.sin(theta), math.cos(theta)]])
    fitted = coords @ rotation
    if flip_across_attachment_axis:
        if stereo:
            raise ValueError("Source layout reflection cannot be applied to a stereospecified fragment")
        anchor = fitted[ports[incoming][0]]
        axis = fitted[ports[outgoing][0]] - anchor
        axis /= np.linalg.norm(axis)
        relative = fitted - anchor
        fitted = anchor + 2*np.outer(relative @ axis, axis) - relative
    return fitted


def _refine_block_spacing(mol, coordinates, origins):
    """Translate rigid blocks to shorten overextended connection bonds.

    This is an optional 2D depiction operation. It never rotates, reflects,
    rescales or deforms a block, and does not change the molecular graph.
    A failed closure or a newly introduced atom overlap keeps the old layout.
    """
    original = np.asarray(coordinates, dtype=float)
    names = list(dict.fromkeys(o['instance'] for o in origins))
    block = np.array([names.index(o['instance']) for o in origins])
    bonds = [(b.GetBeginAtomIdx(), b.GetEndAtomIdx()) for b in mol.GetBonds()]
    joins = [(i, j) for i, j in bonds if block[i] != block[j]]
    internal = [np.linalg.norm(original[i]-original[j]) for i, j in bonds if block[i] == block[j]]
    target = float(np.median(internal))
    audit = {'method': 'rigid_block_translation_v1', 'applied': False,
             'target_bond_length': target, 'interpretation': '2D spacing only; not a physical energy minimization.'}
    if not joins or target <= 0 or not np.isfinite(original).all():
        raise ValueError('Cannot refine a degenerate block layout')
    incidence = np.zeros((len(joins)+1, len(names)))
    incidence[-1] = 1  # Preserve the centroid of the block translations.
    for n, (i, j) in enumerate(joins):
        incidence[n, block[i]] = -1; incidence[n, block[j]] = 1
    initial_lengths = np.array([np.linalg.norm(original[j]-original[i]) for i, j in joins])
    audit['initial_maximum_bond_error'] = float(max(abs(initial_lengths-target)))
    current = original.copy()
    for iteration in range(200):
        rhs = np.zeros((len(joins)+1, 2))
        for n, (i, j) in enumerate(joins):
            delta = current[j]-current[i]; length = np.linalg.norm(delta)
            if length < 1e-8:
                audit['reason'] = 'coincident_connection_endpoints'
                return original, audit
            rhs[n] = target*delta/length-delta
        shift = np.linalg.lstsq(incidence, rhs, rcond=None)[0]
        current += shift[block]
        lengths = np.array([np.linalg.norm(current[j]-current[i]) for i, j in joins])
        if max(abs(lengths-target)) < target/10000:
            break
    audit.update(iterations=iteration+1, final_maximum_bond_error=float(max(abs(lengths-target))))
    if audit['final_maximum_bond_error'] >= target/10000:
        audit['reason'] = 'connection_lengths_did_not_converge'
        return original, audit
    bonded = {tuple(sorted(pair)) for pair in bonds}
    def overlaps(xy):
        result = set()
        for i in range(len(xy)):
            for j in range(i+1, len(xy)):
                if block[i] != block[j] and (i, j) not in bonded and np.linalg.norm(xy[i]-xy[j]) < .35*target:
                    result.add((i, j))
        return result
    new_overlaps = overlaps(current)-overlaps(original)
    if new_overlaps:
        audit.update(reason='new_interblock_atom_overlap', new_overlap_pairs=sorted(new_overlaps))
        return original, audit
    translations = {name: (current-original)[np.flatnonzero(block==n)[0]].tolist() for n, name in enumerate(names)}
    audit.update(applied=True, block_translations=translations,
                 maximum_translation=float(max(np.linalg.norm(current-original, axis=1))))
    return current, audit


def build_patch(recipe: dict):
    """Assemble explicit port incidence; return RDKit graph + portable audit.

    The current layout helper accepts a single closed pore sequence. Every
    chemical connection is explicit; degrees do not select a topology.
    """
    if recipe.get('explicit_graph'):
        from .cofdraw_explicit import build_explicit_patch
        return build_explicit_patch(recipe)
    if recipe.get("schema_version") != VERSION:
        raise ValueError("Unsupported drawing recipe version")
    instances = recipe["instances"]
    if len(instances) < 3 or len({x["id"] for x in instances}) != len(instances):
        raise ValueError("A pore needs at least three uniquely identified instances")
    if recipe.get("topology_basis") == "monomer_degree_only":
        raise ValueError("Monomer degree alone cannot establish topology")
    components, parts, offsets, used, positions, coords = recipe["components"], {}, {}, set(), {}, {}
    joined = Chem.Mol()
    atom_origin = []
    for i, item in enumerate(instances):
        name = item["id"]; component = components[item["component"]]
        mol = Chem.MolFromSmiles(component["fragment_smiles"])
        if mol is None:
            raise ValueError("Invalid reacted component fragment")
        _ports(mol)
        position = np.asarray(item["layout_position"], dtype=float)
        prev = instances[(i-1) % len(instances)]["layout_position"]
        nxt = instances[(i+1) % len(instances)]["layout_position"]
        if recipe.get("layout_kind") == "open_chain":
            if i == 0:
                prev = (2*position-np.asarray(nxt)).tolist()
            if i == len(instances)-1:
                nxt = (2*position-np.asarray(prev)).tolist()
        coords[name] = _fragment_layout(mol, item["incoming_port"], item["outgoing_port"], position, prev, nxt,
                                       item.get("flip_across_attachment_axis", False),
                                       recipe.get('layout_refinement') == 'bond_angle_pose_v2',
                                       bool(recipe.get('designed_variant')) and not recipe.get('atom_coordinates'))
        angle = item.get('layout_rotation_degrees', 0)
        if isinstance(angle, bool) or not isinstance(angle, (int, float)) or not math.isfinite(angle):
            raise ValueError('Block rotation must be a finite angle')
        if angle:
            theta = math.radians(angle)
            coords[name] = coords[name] @ np.array([[math.cos(theta), math.sin(theta)], [-math.sin(theta), math.cos(theta)]])
        offsets[name] = joined.GetNumAtoms(); positions[name] = position; parts[name] = mol
        joined = Chem.CombineMols(joined, mol)
        atom_origin += [{"instance": name, "component": item["component"], "fragment_atom_index": a.GetIdx()} for a in mol.GetAtoms()]
    # Choose a common drawing scale from endpoint positions and bond lengths.
    numerator, denominator = 0.0, 0.0
    for edge in recipe["connections"]:
        a, b = edge["a"], edge["b"]
        pa = _ports(parts[a["instance"]])[a["port"]][1]
        pb = _ports(parts[b["instance"]])[b["port"]][1]
        delta = positions[b["instance"]] - positions[a["instance"]]
        offset = coords[b["instance"]][pb] - coords[a["instance"]][pa]
        numerator += 1.5*np.linalg.norm(delta) - np.dot(delta, offset)
        denominator += np.dot(delta, delta)
    if denominator <= 0:
        raise ValueError("Degenerate layout positions")
    scale = max(1.0, numerator/denominator)
    # A dummy port can itself be the neighboring atom that defines a double
    # bond's E/Z state. RDKit RemoveAtom keeps the stereo enum but clears those
    # selector indices. Capture them before joining, then replace a consumed
    # dummy selector by the real atom on the opposite side of that connection.
    specified_double_stereo = [
        (b.GetBeginAtomIdx(), b.GetEndAtomIdx(), b.GetStereo(), tuple(b.GetStereoAtoms()))
        for b in joined.GetBonds()
        if b.GetStereo() in (Chem.BondStereo.STEREOE, Chem.BondStereo.STEREOZ,
                             Chem.BondStereo.STEREOCIS, Chem.BondStereo.STEREOTRANS)
    ]
    rw = Chem.RWMol(joined); to_remove = []; consumed_selector_replacements = {}
    all_coords = []; port_joins = []
    for item in instances:
        name = item["id"]
        all_coords.extend(coords[name] + scale*positions[name])
    for edge in recipe["connections"]:
        endpoints = []; edge_dummies = []; directions = []
        for end in (edge["a"], edge["b"]):
            name, port = end["instance"], end["port"]
            if (name, port) in used:
                raise ValueError("Attachment port consumed more than once")
            used.add((name, port))
            dummy, neighbor = _ports(parts[name])[port]
            direction = coords[name][dummy]-coords[name][neighbor]
            directions.append(direction/np.linalg.norm(direction))
            dummy += offsets[name]; neighbor += offsets[name]
            source_order = rw.GetBondBetweenAtoms(dummy, neighbor).GetBondTypeAsDouble()
            if source_order != edge["bond_order"]:
                raise ValueError("Connection bond order disagrees with reacted fragment port")
            endpoints.append(neighbor); to_remove.append(dummy); edge_dummies.append(dummy)
        if edge["bond_order"] not in (1, 2, 3):
            raise ValueError("Unsupported connection bond order")
        rw.AddBond(*endpoints, {1: Chem.BondType.SINGLE, 2: Chem.BondType.DOUBLE, 3: Chem.BondType.TRIPLE}[edge["bond_order"]])
        consumed_selector_replacements[edge_dummies[0]] = endpoints[1]
        consumed_selector_replacements[edge_dummies[1]] = endpoints[0]
        port_joins.append({'a': endpoints[0], 'b': endpoints[1], 'direction_a': directions[0], 'direction_b': directions[1]})
    removed = set(to_remove)
    for idx in sorted(removed, reverse=True):
        rw.RemoveAtom(idx)
    mol = rw.GetMol()
    retained_atom_indices = {old: new for new, old in enumerate(
        i for i in range(joined.GetNumAtoms()) if i not in removed)}
    remapped_selector_count = 0
    for begin, end, state, selectors in specified_double_stereo:
        if len(selectors) != 2 or begin in removed or end in removed:
            raise ValueError("Cannot preserve an explicitly defined double bond at a removed endpoint")
        rewritten = []
        for selector in selectors:
            if selector in removed:
                selector = consumed_selector_replacements[selector]
                remapped_selector_count += 1
            rewritten.append(retained_atom_indices[selector])
        begin, end = retained_atom_indices[begin], retained_atom_indices[end]
        bond = mol.GetBondBetweenAtoms(begin, end)
        if bond.GetBeginAtomIdx() != begin:
            rewritten.reverse()
        bond.SetStereoAtoms(*rewritten)
        bond.SetStereo(state)
    # Component aromatic rings were sanitized independently. Global aromaticity
    # perception can enumerate an exponential number of fused macrocycles in
    # a large COF pore. Preserve those source component assignments, validate
    # valence/kekulization and every other sanitizer, and do not assert a new
    # macrocyclic aromatic circuit across the intercomponent linkages.
    aromatic_atoms = [a.GetIdx() for a in mol.GetAtoms() if a.GetIsAromatic()]
    aromatic_explicit_h = {i: mol.GetAtomWithIdx(i).GetNumExplicitHs() for i in aromatic_atoms}
    aromatic_bonds = [b.GetIdx() for b in mol.GetBonds() if b.GetIsAromatic()]
    for name in ("SANITIZE_CLEANUP", "SANITIZE_CLEANUP_ORGANOMETALLICS", "SANITIZE_PROPERTIES",
                 "SANITIZE_SYMMRINGS", "SANITIZE_KEKULIZE", "SANITIZE_FINDRADICALS",
                 "SANITIZE_SETCONJUGATION", "SANITIZE_SETHYBRIDIZATION", "SANITIZE_CLEANUPCHIRALITY",
                 "SANITIZE_CLEANUPATROPISOMERS", "SANITIZE_ADJUSTHS"):
        if name == "SANITIZE_SYMMRINGS":
            symmetric_depiction_rings(mol)
        else:
            Chem.SanitizeMol(mol, sanitizeOps=getattr(Chem.SanitizeFlags, name))
    for idx in aromatic_atoms:
        mol.GetAtomWithIdx(idx).SetIsAromatic(True)
        # ADJUSTHS sees the temporary Kekule representation and may turn [nH]
        # into a cached implicit H. Restoring aromatic bonds alone then loses
        # that H on the next property-cache refresh (e.g. a porphyrin core).
        # These atoms retain their source-fragment bonding and hydrogen state.
        mol.GetAtomWithIdx(idx).SetNumExplicitHs(aromatic_explicit_h[idx])
    for idx in aromatic_bonds:
        mol.GetBondWithIdx(idx).SetIsAromatic(True)
        mol.GetBondWithIdx(idx).SetBondType(Chem.BondType.AROMATIC)
    mol.UpdatePropertyCache(strict=True)
    # Each SanitizeMol call clears ring information, including operations that
    # do not calculate rings. Cache the selected basis after the final one so
    # stereo direction setup and drawing do not trigger RDL recomputation.
    symmetric_depiction_rings(mol)
    if len(Chem.GetMolFrags(mol)) != 1:
        raise ValueError("Chemical patch must be connected")
    kept_coords = [p for i, p in enumerate(all_coords) if i not in removed]
    kept_origin = [p for i, p in enumerate(atom_origin) if i not in removed]
    refinement = {'method': None, 'applied': False}
    if recipe.get('layout_refinement') not in (None, 'manual') and not recipe.get('atom_coordinates'):
        if recipe['layout_refinement'] in ('rigid_block_pose_v1', 'bond_angle_pose_v2'):
            from .cofdraw_layout import refine_rigid_poses, refine_bond_angles
            if specified_double_stereo or any(a.GetChiralTag() != Chem.ChiralType.CHI_UNSPECIFIED for a in mol.GetAtoms()):
                refinement = {'method': recipe['layout_refinement'], 'applied': False, 'reason': 'explicit_stereochemistry_preserved_without_automatic_reorientation'}
                if recipe['layout_refinement'] == 'bond_angle_pose_v2':
                    kept_coords, spacing = _refine_block_spacing(mol, kept_coords, kept_origin)
                    refinement.update(applied=spacing['applied'], angle_fit_applied=False, quality='needs_review',
                                      spacing=spacing, rigid_intrablock_geometry_preserved=True)
            else:
                joins = [{**j, 'a': retained_atom_indices[j['a']], 'b': retained_atom_indices[j['b']]} for j in port_joins]
                fit = refine_bond_angles if recipe['layout_refinement'] == 'bond_angle_pose_v2' else refine_rigid_poses
                kept_coords, refinement = fit(mol, kept_coords, kept_origin, joins)
        elif recipe['layout_refinement'] == 'rigid_block_translation_v1':
            kept_coords, refinement = _refine_block_spacing(mol, kept_coords, kept_origin)
        else:
            raise ValueError('Unsupported layout refinement')
    if recipe.get('atom_coordinates'):
        rows = recipe['atom_coordinates']
        if not isinstance(rows, list) or len(rows) != len(kept_origin):
            raise ValueError('Provide one coordinate for every retained graph atom')
        positions_by_atom = {(p['instance'], p['fragment_atom_index']): p['xy'] for p in rows}
        if len(positions_by_atom) != len(kept_origin):
            raise ValueError('Repeated atom coordinate key')
        try:
            kept_coords = np.asarray([positions_by_atom[(o['instance'], o['fragment_atom_index'])] for o in kept_origin], dtype=float)
        except (KeyError, TypeError, ValueError) as exc:
            raise ValueError('Atom coordinate keys must match the assembled graph') from exc
        if kept_coords.shape != (len(kept_origin), 2) or not np.isfinite(kept_coords).all() or np.max(abs(kept_coords)) > 100000:
            raise ValueError('Atom coordinates must be finite xy pairs')
        refinement = {'method': 'explicit_atom_coordinates', 'applied': False, 'source_coordinates_preserved': True}
    conf = Chem.Conformer(mol.GetNumAtoms())
    conf.Set3D(False)
    continuation_ports = []
    for i, (atom, xy) in enumerate(zip(mol.GetAtoms(), kept_coords)):
        conf.SetAtomPosition(i, Point3D(float(xy[0]), float(xy[1]), 0))
        if atom.GetAtomicNum() == 0:
            continuation_ports.append({"atom_index": i, "component_port": atom.GetAtomMapNum(), **kept_origin[i]})
            atom.SetAtomMapNum(0)
    mol.RemoveAllConformers(); mol.AddConformer(conf)
    # Recreate slash/backslash neighbor directions after port replacement.
    # Without this, RDKit can omit all stereo from SMILES despite nonempty
    # selector arrays and retained STEREOE/Z enums on the double bonds.
    if specified_double_stereo:
        Chem.SetDoubleBondNeighborDirections(mol)
    # Large symmetric pore graphs can make canonical ranking very expensive.
    # Atom-order serialization is deterministic for a recipe, not a graph key.
    has_stereo = any(a.GetChiralTag() != Chem.ChiralType.CHI_UNSPECIFIED or a.GetIsotope() for a in mol.GetAtoms()) or any(b.GetStereo() != Chem.BondStereo.STEREONONE for b in mol.GetBonds())
    fragment_smiles = Chem.MolToSmiles(mol, isomericSmiles=has_stereo, canonical=False)
    lengths=np.asarray([np.linalg.norm(np.asarray(kept_coords[b.GetBeginAtomIdx()])-kept_coords[b.GetEndAtomIdx()]) for b in mol.GetBonds()])
    median_length=float(np.median(lengths))
    audit = {"schema_version": VERSION, "cofex_id": recipe.get("cofex_id"),
        "recipe_sha256": hashlib.sha256(json.dumps(recipe, sort_keys=True, separators=(",", ":")).encode()).hexdigest(),
        "recipe_id": recipe["recipe_id"], "graph_scope": "finite_designed_pore_patch_with_continuation_ports" if recipe.get("designed_variant") else "finite_source_pore_patch_with_continuation_ports",
        "whole_cof_smiles": None, "fragment_smiles": fragment_smiles,
        "fragment_smiles_canonical": False,
        "fragment_serialization_sha256": hashlib.sha256(fragment_smiles.encode()).hexdigest(),
        "atoms": mol.GetNumAtoms(), "bonds": mol.GetNumBonds(), "building_block_instances": len(instances),
        "connections": len(recipe["connections"]), "continuation_ports": continuation_ports,
        "atom_origin": kept_origin, "topology": recipe.get("topology"),
        "depiction_geometry": {"median_bond_length":median_length,"minimum_bond_length":float(lengths.min()),
            "maximum_bond_length":float(lengths.max()),"max_to_median_ratio":float(lengths.max()/median_length),
            "bond_length_outliers_35percent":int(np.sum(abs(lengths/median_length-1)>.35)),
            "interpretation":"2D drawing spacing diagnostics, not physical strain energies."},
        "layout_refinement": refinement,
        "topology_evidence": recipe.get("topology_evidence"), "reference": recipe.get("reference"),
        "source": recipe.get("source"), "review_status": recipe.get("review_status", "agent_draft"),
        "ground_truth": False, "model_prediction": False, "synthetic_training_admitted": False,
        "aromaticity_scope": "Source component aromatic rings retained; no new global macrocycle aromaticity assignment. Symmetric rings use the per-molecule legacy SSSR basis; all other sanitizers passed.",
        "stereochemistry": "Only explicit source stereochemistry in fragments retained; consumed dummy selectors rebound to connected atoms. Drawing coordinates are not E/Z assignments.",
        "double_bond_stereo_selector_remapping": {"defined_bonds": len(specified_double_stereo),
                                                   "consumed_selectors_rebound": remapped_selector_count},
        "defined_double_bond_geometry": _double_bond_geometry_audit(mol),
        "boundary": "Stars denote outgoing continuation/crop ports, not terminal atoms or isolated molecules."}
    from .cofdraw_layout import layout_diagnostics
    audit['layout_diagnostics'] = layout_diagnostics(mol, kept_coords, kept_origin)
    return mol, audit


def _double_bond_geometry_audit(mol):
    """Compare explicit selector relationships with actual planar coordinates.

    A 2D MOL reader derives alkene/imine stereo from geometry. An E/Z graph
    tag alone cannot rescue a contradictory assembled drawing. This audit does
    not assign stereo to unspecified bonds or alter any graph/coordinates.
    """
    xy = np.asarray(mol.GetConformer().GetPositions())[:, :2]
    violations = []; checked = 0
    for b in mol.GetBonds():
        state = b.GetStereo()
        if state not in (Chem.BondStereo.STEREOE, Chem.BondStereo.STEREOZ,
                         Chem.BondStereo.STEREOCIS, Chem.BondStereo.STEREOTRANS):
            continue
        checked += 1
        selectors = tuple(b.GetStereoAtoms())
        record = {"bond": b.GetIdx(), "atoms": [b.GetBeginAtomIdx(), b.GetEndAtomIdx()],
                  "stereo": str(state), "selector_atoms": list(selectors)}
        if len(selectors) != 2:
            violations.append(record | {"reason": "missing_stereo_selectors"})
            continue
        a, c = record["atoms"]
        axis = xy[c]-xy[a]
        left, right = xy[selectors[0]]-xy[a], xy[selectors[1]]-xy[c]
        def sine(vector):
            scale = np.linalg.norm(axis)*np.linalg.norm(vector)
            return (axis[0]*vector[1]-axis[1]*vector[0])/scale if scale > 1e-12 else 0.0
        sl, sr = sine(left), sine(right)
        if min(abs(sl), abs(sr)) < 1e-4:
            violations.append(record | {"reason": "collinear_or_degenerate_stereo_geometry"})
            continue
        expected_same_side = state in (Chem.BondStereo.STEREOZ, Chem.BondStereo.STEREOCIS)
        if (sl*sr > 0) != expected_same_side:
            violations.append(record | {"reason": "coordinates_contradict_defined_stereo"})
    return {"checked": checked, "valid": not violations, "violations": violations,
            "scope": "Defined double bonds only; no new stereochemistry inferred from coordinates."}


def render_recipe(recipe: dict, output_dir: str | Path, *, size=1600, draw_images=True, prepared_patch=None) -> dict:
    """Write editable SVG, PNG, MOL and full graph/port/provenance JSON."""
    output = Path(output_dir); output.mkdir(parents=True, exist_ok=True)
    mol, audit = prepared_patch if prepared_patch is not None else build_patch(recipe)
    recipe_hash = hashlib.sha256(json.dumps(recipe, sort_keys=True, separators=(",", ":")).encode()).hexdigest()
    if audit.get("recipe_sha256") != recipe_hash:
        raise ValueError("Prepared chemical patch does not match the rendering recipe")
    audit = dict(audit)
    if not audit["defined_double_bond_geometry"]["valid"]:
        raise ValueError("Pore layout contradicts explicit double-bond stereochemistry; "
                         "supply a compatible source-supported layout before MOL/image export")
    draw_mol = Chem.Mol(mol)
    Chem.Kekulize(draw_mol, clearAromaticFlags=True)
    Chem.WedgeMolBonds(draw_mol, draw_mol.GetConformer())
    symmetric_depiction_rings(draw_mol)
    stem = recipe["recipe_id"]
    style=recipe.get('drawing_style',{})
    drawers = [("svg", rdMolDraw2D.MolDraw2DSVG(size, size)),
               ("png", rdMolDraw2D.MolDraw2DCairo(size, size))] if draw_images else []
    for suffix, drawer in drawers:
        opt = drawer.drawOptions(); opt.useBWAtomPalette(); opt.bondLineWidth = style.get('bond_width',3.0)
        opt.padding = .07; opt.minFontSize = style.get('min_font_size',12); opt.maxFontSize = style.get('max_font_size',30)
        opt.addStereoAnnotation = False; opt.includeMetadata = False
        opt.prepareMolsBeforeDrawing = False
        drawer.DrawMolecule(draw_mol); drawer.FinishDrawing()
        value = drawer.GetDrawingText(); target = output / f"{stem}.{suffix}"
        target.write_bytes(value.encode("utf-8") if isinstance(value, str) else value)
    Chem.MolToMolFile(mol, str(output / f"{stem}.mol"))
    audit["recipe_sha256"] = recipe_hash
    audit["artifacts"] = {ext: str(output / f"{stem}.{ext}") for ext in ["svg", "png", "mol"]}
    (output/f"{stem}.json").write_text(json.dumps(audit, indent=2)+"\n", encoding="utf-8")
    return audit


def draw_cofex(cofex_id: str, registry: str | Path, output_dir: str | Path, *, recipe_id=None):
    """Resolve only an existing audited recipe; never guess a net from an ID."""
    registry = Path(registry)
    data = json.loads(registry.read_text(encoding="utf-8"))
    candidates = [x for x in data["recipes"] if x["cofex_id"] == cofex_id and (recipe_id is None or x["recipe_id"] == recipe_id)]
    if len(candidates) != 1:
        raise ValueError(f"Expected one source-supported drawing recipe; found {len(candidates)}. Supply a recipe ID or annotate missing incidence.")
    record = candidates[0]
    path = (registry.parent / record["path"]).resolve()
    if hashlib.sha256(path.read_bytes()).hexdigest() != record["sha256"]:
        raise ValueError("Drawing recipe hash mismatch")
    recipe = json.loads(path.read_text(encoding="utf-8"))
    if recipe["cofex_id"] != cofex_id or recipe["recipe_id"] != record["recipe_id"]:
        raise ValueError("Drawing recipe identity mismatch")
    return render_recipe(recipe, output_dir)
