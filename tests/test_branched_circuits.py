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
