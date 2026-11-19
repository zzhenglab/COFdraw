import pytest
from rdkit import Chem
from cof_vision.cofdraw_template_sheet import read_templates


def test_cdxml_sheet_import_keeps_atoms_and_skips_nonmolecular_fragments():
    xml = '''<CDXML BondLength="30"><page><group>
      <fragment id="1" Name="Methanol"><n id="2" p="0 0"/><n id="3" p="30 0" Element="8"/>
      <b id="4" B="2" E="3"/></fragment>
      <fragment id="5"><n id="6" p="0 0" NodeType="GenericNickname"/></fragment>
      </group></page></CDXML>'''
    result = read_templates(xml)
    assert result['skipped'] == 1
    first, = result['templates']
    assert first['name'] == 'Methanol'
    assert Chem.MolToSmiles(Chem.MolFromMolBlock(first['molblock'])) == 'CO'


@pytest.mark.parametrize('text', ['<!DOCTYPE CDXML [<!ENTITY x "x">]><CDXML/>', '<svg/>', '<CDXML><page/></CDXML>'])
def test_template_sheet_rejects_nonmolecular_input(text):
    with pytest.raises(ValueError):
        read_templates(text)
