"""Chemical export checks; these do not claim a ChemDraw desktop import test."""
import json
import hashlib
import struct
import xml.etree.ElementTree as ET
import pytest
from rdkit import Chem
from rdkit.Chem import rdDepictor

from cof_vision.cofdraw_core import export_chemical_formats


@pytest.mark.parametrize("smiles", [
    "[13CH3][C@H](F)/C=C/Cl",
    r"[13CH3][C@H](F)/C=C\Cl",
    "*c1ccc(C=N*)cc1",
    "[NH3+]CC(=O)[O-]",
    "[NH3]->[Cu+2]<-[NH3]",
])
def test_export_preserves_chemical_features(tmp_path, smiles):
    mol = Chem.MolFromSmiles(smiles)
    rdDepictor.Compute2DCoords(mol)
    Chem.MolToMolFile(mol, str(tmp_path / "fixture.mol"), forceV3000=True)
    expected = Chem.MolToSmiles(mol, isomericSmiles=True)
    formats = export_chemical_formats(tmp_path, "fixture")
    export_report = json.loads((tmp_path / "fixture.chemical-export.json").read_text())
    assert export_report["mol_sha256"] == hashlib.sha256((tmp_path / "fixture.mol").read_bytes()).hexdigest()
    assert (tmp_path / "fixture.sdf").read_bytes().startswith((tmp_path / "fixture.mol").read_bytes().rstrip())
    if Chem.HasChemDrawCDXSupport():
        assert "cdxml" in formats, (tmp_path / "fixture.chemical-export.json").read_text()
    sdf = list(Chem.SDMolSupplier(str(tmp_path / "fixture.sdf"), removeHs=False))
    assert len(sdf) == 1 and sdf[0] is not None
    assert Chem.MolToSmiles(sdf[0], isomericSmiles=True) == expected
    assert sum(a.GetAtomicNum() == 0 for a in sdf[0].GetAtoms()) == sum(
        a.GetAtomicNum() == 0 for a in mol.GetAtoms())
    if "cdxml" in formats:
        data = (tmp_path / "fixture.cdxml").read_text()
        root = ET.fromstring(data)
        assert all(n.get("Element") != "0" for n in root.iter("n"))
        ports = [n for n in root.iter("n") if n.get("NodeType") == "ExternalConnectionPoint"]
        assert len(ports) == sum(a.GetAtomicNum() == 0 for a in mol.GetAtoms())
        assert all(n.get("ExternalConnectionType") == "Wavy" and n.get("Element") is None for n in ports)
        if not ports:
            restored = Chem.MolsFromCDXML(data)
            assert len(restored) == 1
            assert Chem.MolToSmiles(restored[0], isomericSmiles=True) == expected
            assert not restored[0].GetConformer().Is3D()
        # RDKit's CDXML reader ignores ExternalConnectionPoint and substitutes
        # carbon. It is not an authority for ChemDraw boundary-node semantics.


@pytest.mark.parametrize("smiles", ["*NC=C*", "*N=CC=N*", "*=NC=N*", "*C#N"])
def test_continuation_points_use_documented_cdxml_node_type(tmp_path, smiles):
    if not Chem.HasChemDrawCDXSupport():
        pytest.skip("Native writer unavailable")
    mol = Chem.MolFromSmiles(smiles)
    rdDepictor.Compute2DCoords(mol)
    Chem.MolToMolFile(mol, str(tmp_path / "fixture.mol"))
    assert export_chemical_formats(tmp_path, "fixture") == ["sdf", "cdxml"]
    root = ET.fromstring((tmp_path / "fixture.cdxml").read_text())
    nodes = list(root.iter("n")); edges = list(root.iter("b"))
    assert len(nodes) == mol.GetNumAtoms() and len(edges) == mol.GetNumBonds()
    assert root.get("BondLength") != ""
    for atom, node in zip(mol.GetAtoms(), nodes):
        if atom.GetAtomicNum() == 0:
            assert node.get("NodeType") == "ExternalConnectionPoint"
            assert node.get("ExternalConnectionType") == "Wavy"
            assert node.get("Element") is None and node.get("NumHydrogens") is None
            attached = [e for e in edges if node.get("id") in [e.get("B"), e.get("E")]]
            assert len(attached) == 1
            assert float(attached[0].get("Order", "1")) == atom.GetBonds()[0].GetBondTypeAsDouble()


@pytest.mark.parametrize("smiles", ["C*C", "[1*]C", "[*+]C"])
def test_nonboundary_dummy_is_not_silently_converted_to_external_point(tmp_path, smiles):
    if not Chem.HasChemDrawCDXSupport():
        pytest.skip("Native writer unavailable")
    mol = Chem.MolFromSmiles("*C" if smiles == "[*+]C" else smiles)
    rdDepictor.Compute2DCoords(mol)
    if smiles == "[*+]C":
        # RDKit's writer itself cannot serialize a charged dummy. Supply the
        # documented MOL charge record to test the import/export safety guard.
        block = Chem.MolToMolBlock(mol).replace("M  END", "M  CHG  1   1   1\nM  END")
        (tmp_path / "fixture.mol").write_text(block)
    else:
        Chem.MolToMolFile(mol, str(tmp_path / "fixture.mol"), forceV3000=True)
    assert export_chemical_formats(tmp_path, "fixture") == ["sdf"]
    assert not (tmp_path / "fixture.cdxml").exists()


def test_unsupported_native_cdxml_keeps_sdf_available(tmp_path, monkeypatch):
    mol = Chem.MolFromSmiles("Oc1cc(O)cc(O)c1")
    rdDepictor.Compute2DCoords(mol)
    Chem.MolToMolFile(mol, str(tmp_path / "fixture.mol"))
    (tmp_path / "fixture.cdxml").write_text("stale earlier output")
    monkeypatch.setattr(Chem, "HasChemDrawCDXSupport", lambda: False)
    assert export_chemical_formats(tmp_path, "fixture") == ["sdf"]
    assert (tmp_path / "fixture.sdf").is_file()
    assert not (tmp_path / "fixture.cdxml").exists()


def test_cdxml_stereo_loss_is_rejected_without_losing_other_exports(tmp_path, monkeypatch):
    if not Chem.HasChemDrawCDXSupport():
        pytest.skip("Native CDXML writer unavailable")
    mol = Chem.MolFromSmiles("F[C@H](Cl)Br")
    rdDepictor.Compute2DCoords(mol)
    Chem.MolToMolFile(mol, str(tmp_path / "fixture.mol"))
    native = Chem.MolToCDXMLBlock

    def remove_stereo(value):
        changed = Chem.Mol(value)
        Chem.RemoveStereochemistry(changed)
        return native(changed)

    monkeypatch.setattr(Chem, "MolToCDXMLBlock", remove_stereo)
    assert export_chemical_formats(tmp_path, "fixture") == ["sdf"]
    assert not (tmp_path / "fixture.cdxml").exists()
    report = json.loads((tmp_path / "fixture.chemical-export.json").read_text())
    assert report["cdxml"] == "omitted_validation_or_writer_failure"
    assert "stereo" in report["reason"]


def test_export_avoids_full_sanitization_and_canonical_ranking(tmp_path, monkeypatch):
    mol = Chem.MolFromSmiles("*c1ccc(C=N*)cc1")
    rdDepictor.Compute2DCoords(mol)
    Chem.MolToMolFile(mol, str(tmp_path / "fixture.mol"))
    native_mol = Chem.MolFromMolBlock

    def scoped_mol(*args, **kwargs):
        assert kwargs.get("sanitize") is False
        return native_mol(*args, **kwargs)

    def scoped_cdxml(*args, **kwargs):
        raise AssertionError("Native CDXML parser must not run during export")

    def no_global_ranking(*args, **kwargs):
        raise AssertionError("Export must not canonicalize the whole macrocycle")

    monkeypatch.setattr(Chem, "MolFromMolBlock", scoped_mol)
    monkeypatch.setattr(Chem, "MolsFromCDXML", scoped_cdxml)
    monkeypatch.setattr(Chem, "MolToSmiles", no_global_ranking)
    formats = export_chemical_formats(tmp_path, "fixture")
    assert formats == (["sdf", "cdxml"] if Chem.HasChemDrawCDXSupport() else ["sdf"])


@pytest.mark.parametrize("smiles,make_unknown", [("F/C=C/Cl", True), ("FC=CCl", False)])
def test_export_rejects_changed_double_stereo_scope(tmp_path, monkeypatch, smiles, make_unknown):
    if not Chem.HasChemDrawCDXSupport():
        pytest.skip("Native CDXML writer unavailable")
    mol = Chem.MolFromSmiles(smiles)
    rdDepictor.Compute2DCoords(mol)
    Chem.MolToMolFile(mol, str(tmp_path / "fixture.mol"))
    native = Chem.MolToCDXMLBlock

    def change_unknown(value):
        root = ET.fromstring(native(value))
        for bond in root.iter("b"):
            if bond.get("Order") == "2":
                if make_unknown:
                    bond.set("Display", "Wavy")
                else:
                    assert bond.get("Display") == "Wavy"
                    bond.attrib.pop("Display")
        return ET.tostring(root, encoding="unicode")

    monkeypatch.setattr(Chem, "MolToCDXMLBlock", change_unknown)
    assert export_chemical_formats(tmp_path, "fixture") == ["sdf"]
    report = json.loads((tmp_path / "fixture.chemical-export.json").read_text())
    assert "double-bond stereo" in report["reason"]


@pytest.mark.parametrize("offset,factor", [((-10000., -6000.), 200.), ((3000., -900.), .02), ((0., 0.), 1.)])
def test_cdxml_page_fit_is_explicit_positive_and_preserves_graph(tmp_path, offset, factor):
    """Page fit must not depend on ChemDraw's last printer or molecule origin."""
    if not Chem.HasChemDrawCDXSupport():
        pytest.skip("Native CDXML writer unavailable")
    mol = Chem.MolFromSmiles("*c1cc(C=Nc2ccc(*)cc2)cc(*)c1")
    rdDepictor.Compute2DCoords(mol)
    conf = mol.GetConformer()
    for i in range(mol.GetNumAtoms()):
        p = conf.GetAtomPosition(i)
        conf.SetAtomPosition(i, (p.x * factor + offset[0], p.y * factor + offset[1], 0))
    Chem.MolToMolFile(mol, str(tmp_path / "fixture.mol"), forceV3000=True)
    original = (tmp_path / "fixture.mol").read_bytes()
    assert export_chemical_formats(tmp_path, "fixture") == ["sdf", "cdxml"]
    assert (tmp_path / "fixture.mol").read_bytes() == original
    root = ET.fromstring((tmp_path / "fixture.cdxml").read_bytes())
    page = root.find("page")
    assert page.get("BoundingBox") == "0 0 540 720"
    assert page.get("WidthPages") == page.get("HeightPages") == "1"
    assert page.get("DrawingSpace") == "Pages"
    assert root.get("PrintMargins") == "36 36 36 36"
    words = struct.unpack(">60h", bytes.fromhex(root.get("MacPrintInfo")))
    assert words[2:4] == (72, 72)
    assert words[4:12] == (0, 0, 792, 612, 0, 0, 792, 612)
    assert words[12:15] == (-1, 792, 612)
    assert words[24:26] == (100, 100)
    for n in root.iter("n"):
        x, y = map(float, n.get("p").split())
        assert 24 - 1e-5 <= x <= 516 + 1e-5
        assert 24 - 1e-5 <= y <= 696 + 1e-5
        assert n.get("Element") != "0"
    bbox = list(map(float, root.find("page/fragment").get("BoundingBox").split()))
    assert 0 < bbox[0] < bbox[2] < 540 and 0 < bbox[1] < bbox[3] < 720
    assert len(list(root.iter("n"))) == mol.GetNumAtoms()
    assert len(list(root.iter("b"))) == mol.GetNumBonds()
    report = json.loads((tmp_path / "fixture.chemical-export.json").read_text())
    assert report["page_geometry"]["uniform_scale"] > 0
    assert report["page_geometry"]["label_size_points"] > 0


def test_verified_source_mol_fast_path_matches_parser_and_unknown_stereo(tmp_path, monkeypatch):
    if not Chem.HasChemDrawCDXSupport():
        pytest.skip("Native writer unavailable")
    mol = Chem.MolFromSmiles("*c1ccc(C=N*)cc1")
    rdDepictor.Compute2DCoords(mol)
    before = Chem.MolToMolBlock(mol)
    (tmp_path / "fixture.mol").write_text(before)
    assert export_chemical_formats(tmp_path, "fixture") == ["sdf", "cdxml"]
    slow_xml = (tmp_path / "fixture.cdxml").read_text()

    def no_parse(*args, **kwargs):
        raise AssertionError("Verified ordinary molecule should avoid parser")

    monkeypatch.setattr(Chem, "MolFromMolBlock", no_parse)
    assert export_chemical_formats(tmp_path, "fixture", source_mol=mol) == ["sdf", "cdxml"]
    assert (tmp_path / "fixture.cdxml").read_text() == slow_xml
    assert Chem.MolToMolBlock(mol) == before
    report = json.loads((tmp_path / "fixture.chemical-export.json").read_text())
    assert report["molecule_input"] == "verified_source_mol_fast_path"
    assert any(b.get("Display") == "Wavy" and b.get("Order") == "2"
               for b in ET.fromstring(slow_xml).iter("b"))


@pytest.mark.parametrize("smiles", ["F[C@H](Cl)Br", "F/C=C/Cl", r"F/C=C\Cl"])
def test_source_mol_defined_stereo_uses_validated_parser_fallback(tmp_path, smiles):
    if not Chem.HasChemDrawCDXSupport():
        pytest.skip("Native writer unavailable")
    mol = Chem.MolFromSmiles(smiles); rdDepictor.Compute2DCoords(mol)
    Chem.MolToMolFile(mol, str(tmp_path / "fixture.mol"))
    assert export_chemical_formats(tmp_path, "fixture", source_mol=mol) == ["sdf", "cdxml"]
    report = json.loads((tmp_path / "fixture.chemical-export.json").read_text())
    assert report["molecule_input"] == "scoped_MOL_parser"
    assert "stereo" in report["source_mol_fallback_reason"]


def test_source_mol_mismatch_cannot_replace_stored_graph(tmp_path):
    if not Chem.HasChemDrawCDXSupport():
        pytest.skip("Native writer unavailable")
    stored = Chem.MolFromSmiles("NC(=O)O"); rdDepictor.Compute2DCoords(stored)
    wrong = Chem.MolFromSmiles("CCO"); rdDepictor.Compute2DCoords(wrong)
    Chem.MolToMolFile(stored, str(tmp_path / "fixture.mol"))
    assert export_chemical_formats(tmp_path, "fixture", source_mol=wrong) == ["sdf", "cdxml"]
    report = json.loads((tmp_path / "fixture.chemical-export.json").read_text())
    assert report["molecule_input"] == "scoped_MOL_parser"
    assert "does not serialize" in report["source_mol_fallback_reason"]
    restored = Chem.MolsFromCDXML((tmp_path / "fixture.cdxml").read_text())[0]
    assert Chem.MolToSmiles(restored) == Chem.MolToSmiles(stored)
