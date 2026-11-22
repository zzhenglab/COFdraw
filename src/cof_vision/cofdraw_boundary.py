"""Presentation cuts at the aromatic ring adjoining an open linkage."""
from rdkit import Chem

from .depiction_rings import symmetric_depiction_rings


def aryl_cut_sites(mol):
    """Return (linkage atom, aryl atom, continuation port) in model indices."""
    cuts = []
    for port in mol.GetAtoms():
        if port.GetAtomicNum() != 0 or port.GetDegree() != 1:
            continue
        middle = port.GetNeighbors()[0]
        if (middle.GetAtomicNum() not in (6, 7) or middle.GetDegree() != 2
                or middle.GetIsAromatic() or middle.GetFormalCharge() or middle.GetIsotope()
                or middle.GetNumRadicalElectrons() or middle.GetChiralTag() != Chem.ChiralType.CHI_UNSPECIFIED):
            continue
        inside = next(a for a in middle.GetNeighbors() if a.GetIdx() != port.GetIdx())
        if inside.GetAtomicNum() != 6 or not inside.GetIsAromatic():
            continue
        attachment = mol.GetBondBetweenAtoms(inside.GetIdx(), middle.GetIdx())
        if attachment.GetBondType() != Chem.BondType.SINGLE:
            continue
        if middle.GetAtomicNum() == 6 and port.GetBonds()[0].GetBondType() != Chem.BondType.DOUBLE:
            continue
        cuts.append((middle.GetIdx(), inside.GetIdx(), port.GetIdx()))
    return cuts


def aryl_boundary(mol, *, cuts=None):
    """Return a cropped figure copy, or None when no eligible terminal stub exists."""
    cuts = aryl_cut_sites(mol) if cuts is None else cuts
    if not cuts:
        return None
    clipped = Chem.RWMol(mol)
    for middle, inside, port in cuts:
        clipped.AddBond(inside, port, Chem.BondType.SINGLE)
        for conf in clipped.GetConformers():
            conf.SetAtomPosition(port, conf.GetAtomPosition(middle))
    for i in sorted({row[0] for row in cuts}, reverse=True):
        clipped.RemoveAtom(i)
    clipped.UpdatePropertyCache(strict=True)
    symmetric_depiction_rings(clipped)
    return clipped.GetMol()


def aryl_figure(mol, audit, recipe, cuts):
    cropped = aryl_boundary(mol, cuts=cuts)
    removed = {middle for middle, _, _ in cuts}
    return remap_figure(mol, cropped, removed, audit, recipe)


def remap_figure(mol, cropped, removed, audit, recipe):
    """Remap display attributes onto a figure copy without editing the recipe."""
    kept = [i for i in range(mol.GetNumAtoms()) if i not in removed]
    new_index = {old: new for new, old in enumerate(kept)}
    figure_audit = {**audit, 'model_indices': kept, 'atom_origin': [audit['atom_origin'][i] for i in kept],
                    'continuation_ports': [{**p, 'atom_index': new_index[p['atom_index']]}
                                           for p in audit['continuation_ports']]}
    source = dict(recipe.get('source_style') or {})
    for key in ('atom_opacities', 'atom_colors'):
        if key in source:
            values = source[key]
            source[key] = {str(new): values.get(str(old), values.get(old))
                           for new, old in enumerate(kept) if str(old) in values or old in values}
    if 'faded_atoms' in source:
        source['faded_atoms'] = [new_index[i] for i in source['faded_atoms'] if i in new_index]
    return cropped, figure_audit, {**recipe, 'source_style': source}
