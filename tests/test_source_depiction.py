import hashlib
import numpy as np
import pytest
from PIL import Image
from rdkit import Chem
from cof_vision.source_depiction import render_source_svg, overlay_metrics


def test_collapsed_group_keeps_chemical_oxygen_atoms_but_not_fake_observed_coordinates(tmp_path):
    mol=Chem.MolFromSmiles('CS(=O)(=O)O')
    before=Chem.MolToSmiles(mol)
    audit={'atom_origin':[{'instance':'a','component':'A','fragment_atom_index':i} for i in range(5)],
           'fragment_serialization_sha256':hashlib.sha256(before.encode()).hexdigest()}
    annotation={'canvas_size':[100,100],'atom_coordinates':[
        {'instance':'a','fragment_atom_index':0,'xy':[20,50]},
        {'instance':'a','fragment_atom_index':1,'xy':[60,50]}],
        'style':{'collapsed_groups':[{'anchor':'a:1','hidden_atoms':['a:2','a:3','a:4'],'label':'SO3H'}]}}
    out=tmp_path/'graph.svg';result=render_source_svg(mol,audit,annotation,out)
    assert result['graph_atoms']==5 and result['source_anchored_atoms']==2
    assert result['graph_atoms_hidden_by_source_abbreviation']==3
    assert Chem.MolToSmiles(mol)==before
    assert 'baseline-shift="sub">3' in out.read_text()


def test_missing_source_coordinate_is_not_invented(tmp_path):
    mol=Chem.MolFromSmiles('CO');audit={'atom_origin':[{'instance':'a','component':'A','fragment_atom_index':i} for i in range(2)]}
    with pytest.raises(ValueError,match='Missing source coordinate'):
        render_source_svg(mol,audit,{'canvas_size':[100,100],'atom_coordinates':[]},tmp_path/'bad.svg')


def test_condensed_aldehyde_label_is_drawn_at_carbon_and_retains_oxygen(tmp_path):
    import xml.etree.ElementTree as ET
    mol=Chem.MolFromSmiles('CC=O');before=Chem.MolToSmiles(mol)
    audit={'atom_origin':[{'instance':'a','component':'A','fragment_atom_index':i} for i in range(3)],
           'fragment_serialization_sha256':hashlib.sha256(before.encode()).hexdigest()}
    ann={'canvas_size':[100,100],'atom_coordinates':[
        {'instance':'a','fragment_atom_index':0,'xy':[20,50]},
        {'instance':'a','fragment_atom_index':1,'xy':[70,50]}],
        'style':{'collapsed_groups':[{'anchor':'a:1','hidden_atoms':['a:2'],'label':'CHO'}]}}
    out=tmp_path/'aldehyde.svg';result=render_source_svg(mol,audit,ann,out)
    svg=ET.parse(out).getroot();ns='{http://www.w3.org/2000/svg}'
    assert [''.join(t.itertext()) for t in svg.findall(ns+'text')]==['CHO']
    assert len(svg.findall(ns+'path'))==1
    assert result['graph_atoms_hidden_by_source_abbreviation']==1
    assert Chem.MolToSmiles(mol)==before


def test_cross_wave_retains_straight_continuation_bond_and_graph(tmp_path):
    import re
    import xml.etree.ElementTree as ET
    mol=Chem.MolFromSmiles('C[*:1]');before=Chem.MolToSmiles(mol)
    audit={'atom_origin':[{'instance':'a','component':'A','fragment_atom_index':i} for i in range(2)],
           'fragment_serialization_sha256':hashlib.sha256(before.encode()).hexdigest()}
    annotation={'canvas_size':[100,100],'atom_coordinates':[
        {'instance':'a','fragment_atom_index':0,'xy':[20,50]},
        {'instance':'a','fragment_atom_index':1,'xy':[80,50]}],
        'style':{'port_style':'wavy_cross','wave_cross_span':20,'wave_amplitude':2}}
    out=tmp_path/'cross.svg';render_source_svg(mol,audit,annotation,out)
    paths=ET.parse(out).getroot().findall('{http://www.w3.org/2000/svg}path')
    coords=[np.array([float(v) for v in re.findall(r'-?\d+\.\d+',p.attrib['d'])]).reshape(-1,2) for p in paths]
    assert any(np.allclose(p,[[20,50],[80,50]]) for p in coords if len(p)==2)
    wave=next(p for p in coords if len(p)>2)
    assert np.ptp(wave[:,1])==pytest.approx(20)
    assert np.ptp(wave[:,0])<=4.01
    assert Chem.MolToSmiles(mol)==before


def test_overlay_tolerance_and_yellow_background_are_not_chemical_accuracy(tmp_path):
    source=np.full((20,20,3),255,np.uint8);target=source.copy()
    source[5:15,5]=[0,0,0];target[5:15,6]=[0,0,0]
    source[1:4,1:18]=[255,255,0]
    a=tmp_path/'a.png';b=tmp_path/'b.png';Image.fromarray(source).save(a);Image.fromarray(target).save(b)
    result=overlay_metrics(a,b,tolerance_pixels=1,pixel_threshold=200,ink_measure='luminance')
    assert result['ink_iou']==0 and result['tolerant_ink_f1']==1
    assert result['not_chemical_accuracy'] and result['fit_uses_source_pixels']


def test_portable_dashes_and_unicode_subscript_render_without_changing_graph(tmp_path):
    import xml.etree.ElementTree as ET
    mol=Chem.MolFromSmiles('N[*:1]');before=Chem.MolToSmiles(mol)
    audit={'atom_origin':[{'instance':'a','component':'A','fragment_atom_index':i} for i in range(2)],
           'fragment_serialization_sha256':hashlib.sha256(before.encode()).hexdigest()}
    ann={'canvas_size':[100,100],'atom_coordinates':[
        {'instance':'a','fragment_atom_index':0,'xy':[20,50]},
        {'instance':'a','fragment_atom_index':1,'xy':[80,50]}],
        'style':{'port_style':'dashed','port_dasharray':[10,10],'subscript_style':'unicode','atom_labels':{'a:0':'NH2'}}}
    out=tmp_path/'portable.svg';render_source_svg(mol,audit,ann,out)
    root=ET.parse(out).getroot();ns='{http://www.w3.org/2000/svg}'
    assert [p.attrib['d'] for p in root.findall(ns+'path')]==[
        'M 20.000,50.000 L 30.000,50.000','M 40.000,50.000 L 50.000,50.000','M 60.000,50.000 L 70.000,50.000']
    assert root.find(ns+'text').text=='NH₂'
    assert Chem.MolToSmiles(mol)==before


@pytest.mark.parametrize('style,expected_paths', [('bond_stub',3),('hidden',1)])
def test_open_imine_port_preserves_visible_multiplicity_and_hidden_continuation(tmp_path,style,expected_paths):
    import xml.etree.ElementTree as ET
    mol=Chem.MolFromSmiles('CN=[*:1]');before=Chem.MolToSmiles(mol)
    audit={'atom_origin':[{'instance':'a','component':'A','fragment_atom_index':i}for i in range(3)],
           'fragment_serialization_sha256':hashlib.sha256(before.encode()).hexdigest()}
    ann={'canvas_size':[120,100],'atom_coordinates':[
        {'instance':'a','fragment_atom_index':i,'xy':[20+40*i,50]}for i in range(3)],
        'style':{'port_style':style,'double_offset':2}}
    out=tmp_path/'imine.svg';render_source_svg(mol,audit,ann,out)
    paths=ET.parse(out).getroot().findall('{http://www.w3.org/2000/svg}path')
    assert len(paths)==expected_paths
    if style=='bond_stub':
        assert {p.attrib['d']for p in paths[1:]}=={'M 60.000,48.000 L 100.000,48.000','M 60.000,52.000 L 100.000,52.000'}
    assert Chem.MolToSmiles(mol)==before


def test_transverse_cut_mark_does_not_add_chemical_bond(tmp_path):
    import xml.etree.ElementTree as ET
    mol=Chem.MolFromSmiles('C=N[*:1]');before=Chem.MolToSmiles(mol)
    audit={'atom_origin':[{'instance':'a','component':'A','fragment_atom_index':i}for i in range(3)],
           'fragment_serialization_sha256':hashlib.sha256(before.encode()).hexdigest()}
    ann={'canvas_size':[120,100],'atom_coordinates':[
        {'instance':'a','fragment_atom_index':i,'xy':[20+40*i,50]}for i in range(3)],
        'style':{'port_style':'crossbar','crossbar_span':10}}
    out=tmp_path/'crossbar.svg';render_source_svg(mol,audit,ann,out)
    paths=ET.parse(out).getroot().findall('{http://www.w3.org/2000/svg}path')
    assert [p.attrib['d']for p in paths[-2:]]==['M 60.000,50.000 L 100.000,50.000','M 100.000,45.000 L 100.000,55.000']
    assert mol.GetNumBonds()==2 and Chem.MolToSmiles(mol)==before


def test_double_bond_crossbar_preserves_imine_and_separate_marker_color(tmp_path):
    import hashlib
    import xml.etree.ElementTree as ET
    mol=Chem.MolFromSmiles('CN=[*:1]');before=Chem.MolToSmiles(mol)
    audit={'atom_origin':[{'instance':'a','component':'A','fragment_atom_index':i}for i in range(3)],
           'fragment_serialization_sha256':hashlib.sha256(before.encode()).hexdigest()}
    ann={'canvas_size':[120,100],'atom_coordinates':[
        {'instance':'a','fragment_atom_index':i,'xy':[20+40*i,50]}for i in range(3)],
        'style':{'port_style':'crossbar','crossbar_span':10,'port_colors':{'a:2':'#ff0000'},'crossbar_color':'#000000'}}
    out=tmp_path/'imine-crossbar.svg';render_source_svg(mol,audit,ann,out)
    paths=ET.parse(out).getroot().findall('{http://www.w3.org/2000/svg}path')
    assert [p.attrib['d']for p in paths[-3:]]==['M 60.000,48.500 L 100.000,48.500','M 60.000,51.500 L 100.000,51.500','M 100.000,45.000 L 100.000,55.000']
    assert [p.attrib['stroke']for p in paths[-3:]]==['#ff0000','#ff0000','#000000']
    assert mol.GetBondWithIdx(1).GetBondTypeAsDouble()==2 and Chem.MolToSmiles(mol)==before


@pytest.mark.parametrize('smiles,strokes', [('CN=[*:1]',2),('CC#[*:1]',3)])
def test_wavy_continuation_retains_visible_boundary_bond_order(tmp_path,smiles,strokes):
    import xml.etree.ElementTree as ET
    mol=Chem.MolFromSmiles(smiles);before=Chem.MolToSmiles(mol)
    audit={'atom_origin':[{'instance':'a','component':'A','fragment_atom_index':i}for i in range(3)],
           'fragment_serialization_sha256':hashlib.sha256(before.encode()).hexdigest()}
    ann={'canvas_size':[120,100],'atom_coordinates':[
        {'instance':'a','fragment_atom_index':i,'xy':[20+40*i,50]}for i in range(3)],
        'style':{'port_style':'wavy_cross'}}
    out=tmp_path/'wave.svg';render_source_svg(mol,audit,ann,out)
    paths=ET.parse(out).getroot().findall('{http://www.w3.org/2000/svg}path')
    straight=[p.attrib['d'] for p in paths if p.attrib['d'].count(' L ')==1]
    # One internal bond plus the double/triple boundary, and a separate wave.
    assert len(straight)==1+strokes
    assert len(paths)==2+strokes
    assert len(set(straight))==1+strokes
    assert Chem.MolToSmiles(mol)==before


def test_protonated_secondary_amine_keeps_both_nitrogen_hydrogens(tmp_path):
    import hashlib
    mol=Chem.MolFromSmiles('C[NH2+]C.[Cl-]');before=Chem.MolToSmiles(mol)
    audit={'atom_origin':[{'instance':'a','component':'A','fragment_atom_index':i}for i in range(4)],
           'fragment_serialization_sha256':hashlib.sha256(before.encode()).hexdigest()}
    ann={'canvas_size':[160,100],'atom_coordinates':[
        {'instance':'a','fragment_atom_index':i,'xy':xy}for i,xy in enumerate([[20,50],[60,30],[100,50],[130,80]])]}
    out=tmp_path/'secondary-ammonium.svg';render_source_svg(mol,audit,ann,out)
    assert 'NH<tspan font-size="70%" baseline-shift="sub">2</tspan>+' in out.read_text()
    assert Chem.MolToSmiles(mol)==before and mol.GetAtomWithIdx(1).GetTotalNumHs()==2
