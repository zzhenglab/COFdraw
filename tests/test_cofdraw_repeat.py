"""A larger drawing must export the same shared-edge graph shown on screen."""
from collections import Counter
from pathlib import Path
from xml.etree import ElementTree
import math

import pytest
from fastapi.testclient import TestClient
from rdkit import Chem

from cof_vision.cofdraw_chemistry import design_signature, prepare_template
from cof_vision.cofdraw_desktop import create_desktop_app
from cof_vision.cofdraw_repeat import repeat_options, repeat_patch
from cof_vision.framework_drawing import VERSION, build_patch
from cof_vision.cofdraw_mixed import ETTA


@pytest.fixture
def api(tmp_path):
    source = Path(__file__).resolve().parents[1]/'data'
    with TestClient(create_desktop_app(source,output_root=tmp_path)) as client:
        yield client


@pytest.mark.parametrize(('endpoint', 'style', 'extra'), [
    ('render', 'ring_fill', {}), ('repeat', 'half_space_stack', {'repeats': 3, 'extent_mode': 'triad'})])
def test_aryl_figure_cut_keeps_the_substituted_connected_model(api, endpoint, style, extra):
    prepared = api.post('/api/drawing/prepare', json={
        'monomer_a': 'O=Cc1cc(C=O)cc(C=O)c1', 'monomer_b': 'Nc1ccc(N)c(F)c1',
        'linkage': 'imine', 'template_id': 'hcb_3_2', 'instances': 12,
    })
    assert prepared.status_code == 200, prepared.text
    body = {'recipe': prepared.json()['recipe'], 'size': 500, 'style_id': style, **extra}
    exports = {}
    for mode in ('linkage', 'aryl'):
        response = api.post('/api/drawing/' + endpoint, json={**body, 'style_options': {'boundary_cut': mode}})
        assert response.status_code == 200, response.text
        urls = response.json()['downloads']
        exports[mode] = {key: api.get(urls[key]) for key in ('json', 'svg', 'mol', 'sdf', 'png')}
    base, short = exports['linkage'], exports['aryl']
    doc = short['json'].json()
    audit = doc['audit']
    boundary = audit['drawing_style']['boundary']
    svg = ElementTree.fromstring(short['svg'].content)
    assert svg.get('data-boundary-cut') == 'aryl'
    assert int(svg.get('data-figure-atoms')) == audit['atoms'] - boundary['shortened_sites']
    assert boundary['shortened_sites'] > 0 and base['png'].content != short['png'].content
    assert short['mol'].content == base['mol'].content and short['sdf'].content == base['sdf'].content
    model = Chem.MolFromMolBlock(short['mol'].text)
    assert len(Chem.GetMolFrags(model)) == 1 and model.GetNumAtoms() == audit['atoms']
    assert any(atom.GetAtomicNum() == 9 for atom in model.GetAtoms())
    assert audit['fragment_smiles'] == base['json'].json()['audit']['fragment_smiles']
    assert len(audit['drawing_style']['layout']['atom_pixels']) == model.GetNumAtoms()
    assert len(audit['drawing_style']['layout']['editing_atom_indices']) < boundary['figure_atoms']
    if endpoint == 'render':
        assert audit['design_cofex_notation'] == base['json'].json()['audit']['design_cofex_notation']
    else:
        assert doc['summary']['figure_boundary'] == boundary
        assert doc['summary']['pore_count'] == 3


def test_kagome_star_cuts_only_its_open_arms_in_the_figure(api):
    prep = api.post('/api/drawing/prepare/assembly', json={'monomer_ids': ['feed_38117f8076cd']})
    assert prep.status_code == 200, prep.text
    recipe = prep.json()['recipe']
    choice = next(row for row in repeat_options(recipe)['extents'] if row['label'] == 'Hexagon + 6 triangles')
    response = api.post('/api/drawing/repeat', json={
        'recipe': recipe, 'repeats': choice['repeats'], 'extent_mode': choice['extent_mode'],
        'size': 500, 'style_id': 'foreground_pore', 'style_options': {'boundary_cut': 'hub'}})
    assert response.status_code == 200, response.text
    summary, urls = response.json()['summary'], response.json()['downloads']
    cut = summary['figure_boundary']
    assert cut['shortened_sites'] == 12 and cut['figure_atoms'] < cut['model_atoms'] - 12
    assert summary['pore_count'] == 7 and summary['continuation_sites'] == 12
    mol = Chem.MolFromMolBlock(api.get(urls['mol']).text)
    assert mol.GetNumAtoms() == cut['model_atoms'] and len(Chem.GetMolFrags(mol)) == 1
    svg = ElementTree.fromstring(api.get(urls['svg']).content)
    assert svg.get('data-boundary-cut') == 'hub' and int(svg.get('data-figure-atoms')) == cut['figure_atoms']


def test_binary_shared_faces_keep_isotopes_and_stereochemistry():
    recipe = prepare_template('O=Cc1cc([13CH]=O)cc(C=O)c1',
                              'Nc1cc([C@H](F)Cl)c(N)cc1/C=C/F','imine','node A','linker B',12,'hcb_3_2')
    recipe.update(schema_version=VERSION,recipe_id='test')
    recipe['design_template_signature']=design_signature(recipe)
    mol,origins,ports,summary,_ = repeat_patch(recipe,2)
    assert summary['pore_count']==4 and len(summary['pore_faces'])==4
    assert summary['blocks']==35 and summary['interblock_bonds']==38 and len(ports)==10
    assert len(Chem.GetMolFrags(mol))==1
    assert summary['component_counts']=={'A':16,'B':19} and len(origins)==mol.GetNumAtoms()
    assert sum(a.GetIsotope()==13 for a in mol.GetAtoms())==16
    assert sum(a.GetChiralTag()!=Chem.ChiralType.CHI_UNSPECIFIED for a in mol.GetAtoms())==19
    assert sum(b.GetStereo() in (Chem.BondStereo.STEREOE,Chem.BondStereo.STEREOTRANS)
               for b in mol.GetBonds())==19
    extents={o['value']:o for o in repeat_options(recipe)['extents']}
    assert extents[2]['atoms']==mol.GetNumAtoms()
    tri,origins,ports,summary,_=repeat_patch(recipe,3,extent_mode='triad')
    faces=[set(face) for face in summary['pore_faces']]
    # Three hexagons meet at a single junction. Every pair shares one
    # linker and the two adjacent junctions, each present once in the MOL.
    assert len(faces)==3 and len(set.intersection(*faces))==1
    assert all(len(faces[i]&faces[i-1])==3 for i in range(3))
    assert summary['component_counts']=={'A':13,'B':15}
    assert summary['interblock_bonds']==30 and len(ports)==9
    assert tri.GetNumAtoms()==extents['three-pores']['atoms']
    assert len({(o['instance'],o['fragment_atom_index']) for o in origins})==tri.GetNumAtoms()


def test_actual_kagome_pore_tints_and_v3000_exports(api):
    request = {'monomer_ids':['tape','ta','bpda_cho'],'linker_sequence':'BCBCBC','pore_layout':'kgm'}
    prepared=api.post('/api/drawing/prepare/assembly',json=request)
    assert prepared.status_code==200,prepared.text
    body={'recipe':prepared.json()['recipe'],'style_id':'standard','size':500}
    options=api.post('/api/drawing/repeat-options',json=body).json()
    assert [o['label'] for o in options['extents'][1:]]==['3 pores','Hexagon + 6 triangles','20 pores']
    result=api.post('/api/drawing/repeat',json={**body,'repeats':3,'color_pores':True})
    assert result.status_code==200,result.text
    report=result.json()
    summary,urls=report['summary'],report['downloads']
    assert summary['pore_sequences']=={'BBB':8,'CCC':8,'BCBCBC':4}
    colors=summary['pore_colors']
    svg=ElementTree.fromstring(api.get(urls['svg']).content)
    counts=Counter(p.get('fill') for p in svg.iter() if p.get('data-pore-fill'))
    assert counts=={colors[label]:n for label,n in summary['pore_sequences'].items()}
    mol=api.get(urls['mol']).content
    assert b'V3000' in mol and b'M  V30 COUNTS 1418 1625' in mol
    assert api.get(urls['png']).content.startswith(b'\x89PNG')
    assert api.get(urls['sdf']).content.startswith(mol.rstrip())
    audit=api.get(urls['json']).json()['audit']
    assert audit['atoms']==1418 and len(audit['atom_origin'])==1418
    if 'cdxml' in urls:
        assert audit['chemical_exports']['molecule_input']=='verified_source_mol_v3000'
        assert len(list(ElementTree.fromstring(api.get(urls['cdxml']).content).iter('n')))==1418
    plain=api.post('/api/drawing/repeat',json={**body,'repeats':2}).json()
    assert plain['summary']['pore_colors']=={} and 'data-pore-fill' not in api.get(plain['svg']).text


def test_expansion_uses_only_reviewed_source_boundaries(api):
    source='wiley-r198-tppa-so3h-pore-0002'
    checked=api.get('/api/drawing/sources/'+source).json()
    selected={'source_recipe_id':source,'size':500}
    response=api.post('/api/drawing/repeat',json={**selected,'repeats':2}).json()
    assert response['summary']['source_drawing_basis']['recipe_sha256']==checked['source_sha256']
    assert response['summary']['pore_count']==4
    other={'source_recipe_id':'pyurea-source-layout-v1','size':500}
    options=api.post('/api/drawing/repeat-options',json=other).json()
    assert options['extents']==[{'value':1,'label':'Source drawing'}] and not options['available']
    assert api.post('/api/drawing/repeat',json={**other,'repeats':2}).status_code==422
    triad=api.post('/api/drawing/repeat',json={**selected,'repeats':3,'extent_mode':'triad'}).json()
    assert triad['summary']['pore_count']==3
    assert triad['summary']['source_drawing_basis']['recipe_sha256']==checked['source_sha256']


@pytest.mark.parametrize(('inputs','expected'),[
    ({'monomer_ids':['bdba']},{'hexagon':3}),
    ({'monomer_ids':['bfbaeb'],'pore_view':'triangle'},{'triangle':2,'hexagon':1}),
    ({'monomer_ids':['tape','ta','bpda_cho'],'linker_sequence':'CCC','pore_layout':'kgm'},
     {'BBB':1,'CCC':1,'BCBCBC':1}),
    ({'monomer_ids':['tapp','ta','bpda_cho'],'linker_sequence':'BCBC'},{'BCBC':3}),
])
def test_three_pore_choice_returns_the_visible_molecule(api,inputs,expected):
    recipe=api.post('/api/drawing/prepare/assembly',json=inputs).json()['recipe']
    body={'recipe':recipe,'size':500}
    option=next(o for o in api.post('/api/drawing/repeat-options',json=body).json()['extents']
                if o['value']=='three-pores')
    response=api.post('/api/drawing/repeat',json={**body,'repeats':option['repeats'],
                      'extent_mode':option['extent_mode'],'color_pores':True})
    assert response.status_code==200,response.text
    result=response.json();summary=result['summary'];urls=result['downloads']
    assert summary.get('pore_counts',summary.get('pore_sequences'))==expected
    assert summary['pore_count']==summary['interblock_bonds']-summary['blocks']+1==3
    assert summary['atoms']==option['atoms']
    svg=ElementTree.fromstring(api.get(urls['svg']).content)
    assert sum(p.get('data-pore-fill') is not None for p in svg.iter())==3
    mol=Chem.MolFromMolBlock(api.get(urls['mol']).text,sanitize=False,removeHs=False)
    assert len(Chem.GetMolFrags(mol))==1 and mol.GetNumAtoms()==summary['atoms']
    measured=api.post('/api/drawing/repeat/'+result['drawing_id']+'/measure')
    assert measured.status_code==200,measured.text
    report=measured.json();aperture=report['measurement'];layout=report['audit']['drawing_style']['layout']
    assert aperture['focal_face']==summary['measurement_face']
    assert aperture['focal_face']['type'] in expected
    assert aperture['diameter_nm']>0 and report['drawing_id']==result['drawing_id']
    for anchor in aperture['annotation']['coordinate_anchors']:
        assert layout['atom_coordinates'][anchor['atom_index']]==pytest.approx(anchor['xy'],abs=1e-8)
    assert len(layout['atom_pixels'])==summary['atoms']


def test_etta_kagome_hexagon_and_each_shared_triangle():
    aldehyde='O=Cc1ccc(C=O)cc1'
    recipe=prepare_template(aldehyde,ETTA,'imine','TPA','ETTA',12,'kgm_4_2')
    recipe.update(schema_version=VERSION,recipe_id='dual-pore',design_template_signature=design_signature(recipe))
    assert recipe['design_topology']=='kgm'
    base,audit=build_patch(recipe)
    ring,_,_,census,_=repeat_patch(recipe,1,extent_mode='triangle')
    assert census['pore_counts']=={'triangle':1} and len(Chem.GetMolFrags(ring))==1
    mol,origins,ports,s,_=repeat_patch(recipe,2)
    assert s['pore_counts']=={'hexagon':1,'triangle':6}
    assert s['component_counts']=={'A':18,'B':12} and s['interblock_bonds']==36 and len(ports)==12
    faces=s['pore_faces']
    center=next(set(face) for face in faces if len(face)==12)
    triangles=[set(face) for face in faces if len(face)==6]
    shared=[center&rim for rim in triangles]
    assert len(triangles)==6 and len({frozenset(edge) for edge in shared})==6
    assert all(len(edge)==3 and sum(i.startswith('A') for i in edge)==1 for edge in shared)
    assert len({(o['instance'],o['fragment_atom_index']) for o in origins})==mol.GetNumAtoms()
    # Both the local face and full connected neighborhood must keep drawing
    # scale consistent between internal phenyl bonds and each new imine.
    for graph,rows in ((mol,origins),(base,audit['atom_origin'])):
        xy=graph.GetConformer().GetPositions()
        for bond in graph.GetBonds():
            a,b=bond.GetBeginAtomIdx(),bond.GetEndAtomIdx()
            if rows[a]['instance']!=rows[b]['instance']:
                assert bond.GetBondType()==Chem.BondType.DOUBLE
                assert {bond.GetBeginAtom().GetSymbol(),bond.GetEndAtom().GetSymbol()}=={'C','N'}
                assert 1.45<math.dist(xy[a],xy[b])<1.55
    for a,b in ((aldehyde,'Nc1cc(N)c(N)cc1N'),('O=Cc1cccc(C=O)c1',ETTA)):
        with pytest.raises(ValueError,match='ETTA and a linear'):
            prepare_template(a,b,'imine','A','B',12,'kgm_4_2')


def test_self_kagome_complete_shell_matches_mixed_kagome_incidence(api):
    for feed in (['bfbaeb'],['tape','ta','bpda_cho']):
        request={'monomer_ids':feed}
        if len(feed)==3:
            request.update(linker_sequence='CCC',pore_layout='kgm')
        recipe=api.post('/api/drawing/prepare/assembly',json=request).json()['recipe']
        choice=next(r for r in repeat_options(recipe)['extents'] if r['label']=='Hexagon + 6 triangles')
        mol,origins,_,s,_=repeat_patch(recipe,2,extent_mode=choice.get('extent_mode','pores'))
        faces=sorted(s['pore_faces'],key=len)
        triangles,center=faces[:-1],set(faces[-1])
        assert len(triangles)==6 and s['pore_count']==s['interblock_bonds']-s['blocks']+1==7
        assert all(len(t)*2==len(center) and len(center&set(t))==(2 if len(feed)==1 else 3) for t in triangles)
        assert len(Chem.GetMolFrags(mol))==1 and len(origins)==choice['atoms']==mol.GetNumAtoms()
