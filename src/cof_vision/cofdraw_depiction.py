"""Scoped coordinate generation for conventional free-base porphyrin drawings.

RDKit's default macrocycle layout flattens the pyrroles. Its built-in ring
template gives the four rings regular pentagonal geometry without changing the
graph. Opt in per molecule, never through global depiction preferences. This
is a 2D drawing convention, not an optimized molecular geometry.
"""
from rdkit import Chem
from rdkit.Chem import rdDepictor


# Porphine connectivity, independent of substituents, atom order, isotope and
# Kekule form. Atom/bond state checks below exclude reduced, charged, metalated
# and N-substituted cores; those keep the ordinary depiction path.
_PORPHINE = Chem.MolFromSmarts(
    "[#6]1~[#6]2~[#7]~[#6](~[#6]~[#6]3~[#6]~[#6]~[#6](~[#7]~3)"
    "~[#6]~[#6]3~[#7]~[#6](~[#6]~[#6]4~[#6]~[#6]~[#6]~1~[#7]~4)"
    "~[#6]~[#6]~3)~[#6]~[#6]~2"
)


def has_free_base_porphyrin(mol):
    """Recognize one neutral, fully conjugated porphine core, not a catalog ID."""
    matches = mol.GetSubstructMatches(_PORPHINE, maxMatches=2)
    if len(matches) != 1:
        return False
    atoms = [mol.GetAtomWithIdx(i) for i in matches[0]]
    if any(a.GetFormalCharge() or a.GetHybridization() != Chem.HybridizationType.SP2
           or a.GetNumRadicalElectrons() for a in atoms):
        return False
    nitrogens = [a for a in atoms if a.GetAtomicNum() == 7]
    return (all(a.GetDegree() == 2 for a in nitrogens)
            and sorted(a.GetTotalNumHs() for a in nitrogens) == [0, 0, 1, 1])


def compute_2d_coords(mol, **kwargs):
    """Generate new coordinates; callers retain authority over existing ones.

    Preserve the caller's original RDKit defaults for all other structures.
    Force RDKit only for the recognized core, because its ring templates are
    not used by CoordGen even if another host has set that global preference.
    """
    if has_free_base_porphyrin(mol):
        kwargs.update(useRingTemplates=True, forceRDKit=True)
    return rdDepictor.Compute2DCoords(mol, **kwargs)
