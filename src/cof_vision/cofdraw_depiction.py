"""Scoped 2D layouts for porphyrin macrocycles and triptycene-based monomers."""
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
        rdCoordGen.AddCoords(mol)
        conformer = mol.GetConformer()
        for atom in mol.GetAtoms():
            index = atom.GetIdx()
            position = conformer.GetAtomPosition(index)
            conformer.SetAtomPosition(index, (1.5 * position.x, 1.5 * position.y, 0.))
        return conformer.GetId()
    return rdDepictor.Compute2DCoords(mol, **kwargs)
