"""Scoped 2D layouts for porphyrin macrocycles and triptycene-based monomers."""
import math

import numpy as np
from rdkit import Chem
from rdkit.Chem import rdCoordGen, rdDepictor


# Porphine connectivity, independent of substituents, atom order, isotope and
# Kekule form. Atom/bond state checks below exclude reduced, charged, metalated
# and N-substituted cores; those keep the ordinary depiction path.
_PORPHINE = Chem.MolFromSmarts(
    "[#6]1~[#6]2~[#7]~[#6](~[#6]~[#6]3~[#6]~[#6]~[#6](~[#7]~3)"
    "~[#6]~[#6]3~[#7]~[#6](~[#6]~[#6]4~[#6]~[#6]~[#6]~1~[#7]~4)"
    "~[#6]~[#6]~3)~[#6]~[#6]~2"
)
_TRIPTYCENE = Chem.MolFromSmiles("c1ccc2c(c1)C1c3ccccc3C2c2ccccc21")


def port_symmetry_ranks(mol):
    """Rank chemically equivalent sites independently of their attachment-port numbers."""
    unlabeled = Chem.Mol(mol)
    for atom in unlabeled.GetAtoms():
        atom.SetAtomMapNum(0)
    return Chem.CanonicalRankAtoms(unlabeled, breakTies=False)


def straight_linker_axes(sites, axes):
    """Both outward scaffold bonds must follow opposite ends of the same chord."""
    if len(sites) != 2:
        return False
    chord = sites[1] - sites[0]
    lengths = np.linalg.norm(axes, axis=1)
    span = np.linalg.norm(chord)
    if min(lengths) < 1e-8 or span < 1e-8:
        return False
    alignment = (axes @ chord) / (lengths * span)
    return bool(min(-alignment[0], alignment[1]) > math.cos(math.radians(25)))


def normalize_imine_ports(mol, xy, *, preserve_uniform=False, straight_only=False):
    """Place imine bends around equivalent arms or on opposite sides of a straight linker.

    Only the replaceable dummy positions move. Mapped port numbers encode
    incidence; stereochemistry and isotopes still distinguish different arms.
    """
    ports = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 0]
    if len(ports) < 2 or any(atom.GetDegree() != 1 for atom in ports):
        return False
    stereo_selectors = {i for b in mol.GetBonds() for i in b.GetStereoAtoms()}
    ends = []
    for port in ports:
        site = port.GetNeighbors()[0]
        bond = port.GetBonds()[0]
        if (site.GetAtomicNum() not in (6, 7) or site.GetDegree() != 2 or
                bond.GetBondType() != Chem.BondType.DOUBLE or bond.GetStereo() != Chem.BondStereo.STEREONONE or
                port.GetIdx() in stereo_selectors):
            return False
        anchor = next(a for a in site.GetNeighbors() if a.GetIdx() != port.GetIdx())
        if anchor.GetAtomicNum() != 6:
            return False
        ends.append((port.GetIdx(), site.GetIdx(), anchor.GetIdx()))
    dummies, sites, anchors = np.asarray(ends).T
    axis = xy[sites, :2] - xy[anchors, :2]
    lengths = np.linalg.norm(axis, axis=1)
    if min(lengths) < 1e-8:
        return False
    # Test both aryl/alkyl exit bonds against the same end-to-end chord. A
    # substituent or a mixed aldehyde/amine linker may break atomic equivalence
    # without changing the direction in which its two ends point.
    straight = straight_linker_axes(xy[sites, :2], axis)
    if not straight:
        if straight_only or any(not mol.GetAtomWithIdx(int(i)).GetIsAromatic() for i in anchors):
            return False
        ranks = port_symmetry_ranks(mol)
        if len({ranks[p.GetIdx()] for p in ports}) != 1:
            return False
    outgoing = xy[dummies, :2] - xy[sites, :2]
    if preserve_uniform:
        turns = axis[:, 0]*outgoing[:, 1] - axis[:, 1]*outgoing[:, 0]
        if np.all(turns > 1e-8) or np.all(turns < -1e-8):
            return True
    port_lengths = np.linalg.norm(outgoing, axis=1)
    # Rotate the outward aryl bond by +60 degrees. The angle drawn at the
    # imine atom is then 120 degrees; a whole-block reflection reverses all arms.
    sine = math.sqrt(3.) / 2
    turn = np.array([[.5, sine], [-sine, .5]])
    xy[dummies, :2] = xy[sites, :2] + ((port_lengths / lengths)[:, None] * axis) @ turn
    return True


def free_base_porphyrin_core(mol):
    """Recognize one neutral, fully conjugated porphine core, not a catalog ID."""
    matches = mol.GetSubstructMatches(_PORPHINE, maxMatches=2)
    if len(matches) != 1:
        return ()
    atoms = [mol.GetAtomWithIdx(i) for i in matches[0]]
    if any(a.GetFormalCharge() or a.GetHybridization() != Chem.HybridizationType.SP2
           or a.GetNumRadicalElectrons() for a in atoms):
        return ()
    nitrogens = [a for a in atoms if a.GetAtomicNum() == 7]
    if (all(a.GetDegree() == 2 for a in nitrogens)
            and sorted(a.GetTotalNumHs() for a in nitrogens) == [0, 0, 1, 1]):
        return matches[0]
    return ()


def has_free_base_porphyrin(mol):
    return bool(free_base_porphyrin_core(mol))


def has_triptycene_core(mol):
    return mol.HasSubstructMatch(_TRIPTYCENE)


def coordgen_coords(mol):
    rdCoordGen.AddCoords(mol)
    conformer=mol.GetConformer()
    for i in range(mol.GetNumAtoms()):
        p=conformer.GetAtomPosition(i)
        conformer.SetAtomPosition(i,(1.5*p.x,1.5*p.y,0.))
    return conformer.GetId()


def compute_2d_coords(mol, **kwargs):
    """Generate new coordinates; callers retain authority over existing ones.

    Preserve the caller's original RDKit defaults for all other structures.
    Force RDKit only for the recognized core, because its ring templates are
    not used by CoordGen even if another host has set that global preference.
    """
    if has_free_base_porphyrin(mol):
        kwargs.update(useRingTemplates=True, forceRDKit=True)
    elif not kwargs and mol.HasSubstructMatch(_TRIPTYCENE):
        # CoordGen depicts the three blades separately. Rescale its ~1 Å
        # drawing bonds to the RDKit convention used by the sketch toolbar.
        return coordgen_coords(mol)
    return rdDepictor.Compute2DCoords(mol, **kwargs)
