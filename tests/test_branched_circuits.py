import numpy as np
import pytest
from rdkit import Chem

from cof_vision.cofdraw_chemistry import design_signature,prepare_template
from cof_vision.cofdraw_net_graphs import graph_preview,net_keys,quotient
from cof_vision.cofdraw_repeat import repeat_options,repeat_patch


T_NODE='O=Cc1ccc(-c2nc3c(-c4ccc(C=O)cc4)ccc(-c4ccc(C=O)cc4)c3[nH]2)cc1'
TAPM='Nc1ccc(C(c2ccc(N)cc2)(c2ccc(N)cc2)c2ccc(N)cc2)cc1'


@pytest.mark.parametrize('aldehyde',[T_NODE,T_NODE.replace('O=Cc1ccc','O=Cc1c(F)cc')])
def test_substituted_circuits(aldehyde):
    r=prepare_template(aldehyde,TAPM,'imine','A','B',6,'branched_4_3')
    r['design_template_signature']=design_signature(r)
    assert r['design_topology'] is None
    choices=repeat_options(r)['extents']
    assert [c['value'] for c in choices]==[1,2,3]
    for count in (2,3):
        mol,_,boundaries,s,_=repeat_patch(r,count)
        assert len(Chem.GetMolFrags(mol))==1
        assert s['atoms']==choices[count-1]['atoms']
        assert [len(face) for face in s['pore_faces']]==[6,8,6][:count]
        assert s['interblock_bonds']-s['blocks']+1==count
        assert len(boundaries)==sum(s['component_counts'][k]*r['components'][k]['connectivity']
                                    for k in 'AB')-2*s['interblock_bonds']
        xyz=mol.GetConformer().GetPositions()
        lengths=[np.linalg.norm(xyz[b.GetBeginAtomIdx()]-xyz[b.GetEndAtomIdx()]) for b in mol.GetBonds()]
        assert min(lengths)>1.42 and max(lengths)<1.58


def test_gallery_draws_lifted_edges_and_junctions():
    for code in net_keys():
        g=graph_preview(code)
        assert len(g['edges'])>=len(g['nodes'])
        assert all(0<=i<j<len(g['nodes']) for i,j in g['edges'])
        assert np.isfinite(g['nodes']).all()
    assert sorted(set(map(len,quotient('bor')[0])))==[3,4]


def test_biphenyl_chairs_keep_distinct_junctions_clear(monkeypatch):
    from itertools import combinations
    from cof_vision.cofdraw_dia import build_dia
    from cof_vision.cofdraw_crystal_figures import adjacent_circuits,chemical_drawing,crystal_details
    from cof_vision import cofdraw_reticular_scheme as scheme

    model,audit=build_dia(TAPM,'O=Cc1ccc(-c2ccc(C=O)cc2)cc1')
    edges=audit['periodic_bonds']
    ring=crystal_details(model,edges)['ring']
    def closest(mol,positions,*args):
        junctions=[positions[a.GetIdx()] for a in mol.GetAtoms() if a.GetDegree()==4 and a.GetAtomicNum()==6]
        return min(np.linalg.norm(a-b) for a,b in combinations(junctions,2))
    monkeypatch.setattr(scheme,'vector_art',closest)
    for n in (2,3,4,7):
        chairs=adjacent_circuits(ring,model['cell'],n)
        # Two distinct four-aryl junctions require room for both sets of rings.
        assert chemical_drawing(model,edges,circuit=chairs,scheme=True)>12.


@pytest.mark.parametrize('key,a,b',[
    ('acs_6_2','hfptp_f','bpy55'),('acs_6_6','hfptp_h','hapt'),
    ('stp_6_4','hfptp_h','pytta'),('hea_6_4','hfptp_f','tapm'),
    ('ceq_6_3','hfptp_h','tapt'),('pts_4_4','tfppy','tapm'),
    ('nbo_4_2','ta','pytta'),('tbo_4_3','tfppy','tapb'),('dia_4_2','tfpa_adamantane','pa'),
    ('stp_6_4','recovered_tfbp','hapt'),('pts_4_4','tht_cho4','tapm'),('pts_4_4','tfs_cho4','tapp')])
def test_spatial_design_regenerates_with_selected_roles(key,a,b):
    import json
    from pathlib import Path
    rows={m['id']:m['smiles'] for m in json.loads((Path(__file__).resolve().parents[1]/
          'data/framework-drawing/catalog/current.json').read_bytes())['monomer_presets']}
    recipe=prepare_template(rows[a],rows[b],'imine',a,b,6,key)
    recipe['design_template_signature']=design_signature(recipe)
    assert recipe['design_topology'] is None
    assert repeat_options(recipe)['available']
    mol,orig,ports,summary,_=repeat_patch(recipe,3)
    assert len(Chem.GetMolFrags(mol))==1 and summary['pore_count']==3
    assert len(ports)==sum(recipe['components'][k]['connectivity']*n for k,n in
                           summary['component_counts'].items())-2*summary['interblock_bonds']
    assert len(orig)==mol.GetNumAtoms()
    assert all(k in orig[0] for k in ('component','instance','fragment_atom_index'))
    from cof_vision.cofdraw_science import patch_3d
    graph=patch_3d(recipe)
    assert graph['kind']=='spatial_building_block_net' and 'structure' not in graph
    vertices=graph['viewer']['atoms']
    edges=graph['viewer']['bonds']
    xy=np.array([[v[k] for k in ('x','y','z')] for v in vertices])
    assert np.linalg.matrix_rank(xy-xy.mean(0))==3
    assert {v['label'] for v in vertices}=={a,b}
    assert all(vertices[i]['label']!=vertices[j]['label'] for i,j in edges)
    for role in ('A','B'):
        assert max(sum(i in e for e in edges) for i,v in enumerate(vertices) if v['component']==role)==recipe['components'][role]['connectivity']


def test_opposed_imine_linkers_across_diamond_circuits():
    from collections import defaultdict
    from cof_vision.cofdraw_branched import branch_patch
    recipe = prepare_template('O=Cc1ccc(C=O)cc1', TAPM, 'imine', 'TA', 'TAPM', 12, 'dia_4_2')
    for size, count in ((1, 6), (3, 16)):
        mol, summary = branch_patch(recipe['components'], size, 'dia_4_2')
        owners = summary['atom_origin']
        ends = defaultdict(list)
        for bond in mol.GetBonds():
            a, b = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
            if owners[a]['instance'] == owners[b]['instance']:
                continue
            site, neighbor = (a, b) if owners[a]['component'] == 'A' else (b, a)
            ends[owners[site]['instance']].append((site, neighbor))
        assert len(ends) == count
        xy = mol.GetConformer().GetPositions()[:, :2]
        for (a, left), (b, right) in ends.values():
            axis = xy[b] - xy[a]
            normal = np.array([-axis[1], axis[0]]) / np.linalg.norm(axis)
            assert ((xy[left] - xy[a]) @ normal) * ((xy[right] - xy[b]) @ normal) < -.1
