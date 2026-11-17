"""Ring optimization must preserve chemistry and stereochemical depiction."""
from rdkit import Chem
from rdkit.Chem import rdDepictor
from rdkit.Chem.Draw import rdMolDraw2D
import pytest

from cof_vision.depiction_rings import symmetric_depiction_rings
from cof_vision.cofdraw_chemistry import prepare_template
from cof_vision.framework_drawing import VERSION, build_patch


@pytest.mark.parametrize("smiles", [
    "[13CH3][C@H](F)CC/C=C/Cl", "C[C@H](F)[C@H](Cl)Br |&1:1,3|",
    "c1ccc2ccccc2c1", "C1C2CC3CC1CC(C2)C3", "C[C@H]1CC[C@@H](C)CC1",
])
def test_ring_basis_preserves_stereo_chemistry_and_native_pixels(smiles):
    mol = Chem.MolFromSmiles(smiles)
    rdDepictor.Compute2DCoords(mol)
    Chem.Kekulize(mol, clearAromaticFlags=True)
    Chem.WedgeMolBonds(mol, mol.GetConformer())
    original = Chem.MolToCXSmiles(mol)

    def draw(value):
        drawer = rdMolDraw2D.MolDraw2DCairo(500, 500)
        drawer.drawOptions().prepareMolsBeforeDrawing = False
        drawer.DrawMolecule(value)
        drawer.FinishDrawing()
        return drawer.GetDrawingText()

    before = draw(mol)
    global_before = Chem.GetUseLegacyRingFinding() if hasattr(Chem, "GetUseLegacyRingFinding") else None
    symmetric_depiction_rings(mol)
    assert Chem.MolToCXSmiles(mol) == original
    assert draw(mol) == before
    if global_before is not None:
        assert Chem.GetUseLegacyRingFinding() == global_before


def test_large_pore_has_bounded_symmetric_basis_and_all_local_rings():
    prepared = prepare_template(
        "COc1cc(C=O)c(OC)cc1C=O",
        "Nc1ccc(-c2cc(-c3ccc(N)cc3)cc(-c3ccc(N)cc3)c2)cc1",
        "imine", "DMTP", "TAPB", 12, "hcb_3_2",
    )
    recipe = {"schema_version": VERSION, "recipe_id": "large-ring-regression",
              "layout_refinement": "rigid_block_translation_v1", **prepared}
    mol, audit = build_patch(recipe)
    rings = list(mol.GetRingInfo().AtomRings())
    # All thirty phenyl rings remain available; there is no need to enumerate
    # the 262,174 alternative relevant cycles of this symmetric pore.
    assert sum(len(ring) == 6 for ring in rings) == 30
    assert len(rings) < 1000
    assert mol.GetNumAtoms() == 240
    assert len(audit["continuation_ports"]) == 6
    assert audit["defined_double_bond_geometry"]["valid"]
