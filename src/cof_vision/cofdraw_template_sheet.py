"""Read finite molecular building blocks from a user-supplied ChemDraw sheet."""
import xml.etree.ElementTree as ET
import os
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors


COLLECTION_NAMES = ('Aromatics', 'Bicyclics', 'Cycloalkanes', 'Functional Groups',
                    'Ph Rings', 'Conformers', 'Hexoses', 'Amino Acids')


def installed_sheets():
    root = os.getenv('PROGRAMDATA')
    if not root or not Chem.HasChemDrawCDXSupport():
        return {}
    folder = Path(root) / 'RevvitySignalsSoftware/ChemDrawApplications/ChemDraw/ChemDraw Items'
    return {name: folder / (name + '.ctp') for name in COLLECTION_NAMES
            if (folder / (name + '.ctp')).is_file()}


def read_binary_templates(data):
    if not data.startswith(b'VjCD0100') or len(data) > 500000:
        raise ValueError('Choose a molecular ChemDraw CTP or CDX file up to 500 kB.')
    if not Chem.HasChemDrawCDXSupport():
        raise ValueError('This RDKit build reads XML sheets. Save the structure as CDXML in ChemDraw.')
    params = Chem.CDXMLParserParams()
    params.format = Chem.CDXMLFormat.CDX
    params.removeHs = False
    templates, skipped = [], 0
    molecules = Chem.MolsFromCDXML(data, params)
    if len(molecules) > 160:
        raise ValueError('Use a smaller molecular template sheet.')
    for index, mol in enumerate(molecules, 1):
        if (not 0 < mol.GetNumAtoms() <= 300 or len(Chem.GetMolFrags(mol)) != 1 or
                any(a.HasQuery() or not a.GetAtomicNum() for a in mol.GetAtoms()) or
                any(b.HasQuery() or b.GetBondTypeAsDouble() not in (1, 1.5, 2, 3) for b in mol.GetBonds())):
            skipped += 1
            continue
        templates.append({'name': mol.GetProp('_Name') if mol.HasProp('_Name') else f'Structure {index}',
                          'formula': rdMolDescriptors.CalcMolFormula(mol), 'molblock': Chem.MolToMolBlock(mol)})
    if not templates:
        raise ValueError('This palette has no supported finite molecular structures.')
    return {'templates': templates, 'skipped': skipped}


def read_templates(text):
    if '<!ENTITY' in text.upper():
        raise ValueError('External and custom entities are unsupported.')
    try:
        document = ET.fromstring(text)
    except ET.ParseError as exc:
        raise ValueError('The file does not contain readable CDXML.') from exc
    if document.tag != 'CDXML':
        raise ValueError('Choose a ChemDraw XML (.cdxml) file.')
    if sum(1 for _ in document.iter()) > 12000:
        raise ValueError('Use a smaller template sheet.')
    fragments, pending = [], list(reversed(document))
    while pending:
        obj = pending.pop()
        if obj.tag == 'fragment':
            fragments.append(obj)
        elif obj.tag in ('page', 'group'):
            pending.extend(reversed(obj))
    if len(fragments) > 64:
        raise ValueError('Import up to 64 structures per sheet.')
    templates, skipped = [], 0
    for index, fragment in enumerate(fragments, 1):
        nodes = fragment.findall('n')
        if (not 0 < len(nodes) <= 300 or fragment.find('.//fragment') is not None or
                any(n.get('NodeType', 'Element') not in ('Element', 'Unspecified') for n in nodes)):
            skipped += 1
            continue
        doc = ET.Element('CDXML', document.attrib)
        ET.SubElement(doc, 'page').append(fragment)
        try:
            molecules = Chem.MolsFromCDXML(ET.tostring(doc, encoding='unicode'), removeHs=False)
            if len(molecules) != 1:
                raise ValueError()
            mol = molecules[0]
            if (mol.GetNumAtoms() != len(nodes) or mol.GetNumBonds() != len(fragment.findall('b')) or
                    len(Chem.GetMolFrags(mol)) != 1 or any(a.HasQuery() or not a.GetAtomicNum() for a in mol.GetAtoms())
                    or any(b.HasQuery() for b in mol.GetBonds())):
                raise ValueError()
            templates.append({'name': fragment.get('Name') or f'Structure {index}',
                              'formula': rdMolDescriptors.CalcMolFormula(mol),
                              'molblock': Chem.MolToMolBlock(mol)})
        except (ValueError, RuntimeError):
            skipped += 1
    if not templates:
        raise ValueError('No supported connected molecular structures found. Expand abbreviations and export a molecular CDXML sheet.')
    return {'templates': templates, 'skipped': skipped}
