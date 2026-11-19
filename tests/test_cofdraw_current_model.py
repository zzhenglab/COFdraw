"""Analytical finite-scattering checks and crystallographic stacking invariants."""
from copy import deepcopy
from pathlib import Path

import numpy as np
import pytest

from cof_vision.cofdraw_current_model import inspect_finite, finite_powder
from cof_vision.cofdraw_characterization import atomic_form_factor, powder_xrd
from cof_vision.cofdraw_science import import_cif
from cof_vision.cofdraw_stacking import periodic_stacking


@pytest.mark.parametrize("stereotag", ["@", "@@"])
def test_hth_periodic_cell_retains_all_side_chain_stereocentres(stereotag):
    from rdkit import Chem
    from cof_vision.cofdraw_chemistry import prepare_template
    from cof_vision.cofdraw_periodic_build import build_periodic_layer

    aldehyde = 'CCCCCCNC(=O)c1c(-c2ccc(C=O)cc2)c(C(=O)NCCCCCC)c(-c2ccc(C=O)cc2)c(C(=O)NCCCCCC)c1-c1ccc(C=O)cc1'
    partner = 'C[C@H](O)COc1cc(C(=O)NN)c(OC[C@H](C)O)cc1C(=O)NN'.replace('@', stereotag)
    recipe = prepare_template(aldehyde, partner, 'acylhydrazone', 'A', 'B', 6, 'hcb_3_2')
    cell, validation = build_periodic_layer(recipe)
    assert validation['short_contacts'] == 0 and validation['periodic_rank'] == 2
    sites = cell['atoms']
    assert len(sites) == 348
    fragment = Chem.AddHs(Chem.MolFromSmiles(recipe['components']['B']['fragment_smiles']))
    indices = [a.GetIdx() for a in fragment.GetAtoms() if a.GetAtomicNum()]
    # Recover local tetrahedra across periodic faces directly from the exported sites.
    basis = np.asarray(cell['cell'])
    fractional = np.array([a['fractional'] for a in sites])
    bond_images = {}
    for bond in validation['periodic_bonds']:
        i, j = bond['a'], bond['b']
        bond_images[i, j] = np.array(bond['translation'])
        bond_images[j, i] = -bond_images[i, j]
    starts = [a for a, site in enumerate(sites) if site['component'] == 'B'][::len(indices)]
    verified = 0
    for start in starts:
        for atom in fragment.GetAtoms():
            if atom.GetChiralTag() == Chem.ChiralType.CHI_UNSPECIFIED:
                continue
            i = start + indices.index(atom.GetIdx())
            probe = Chem.Mol(fragment)
            conf = Chem.Conformer(probe.GetNumAtoms())
            for neighbor in atom.GetNeighbors():
                j = start + indices.index(neighbor.GetIdx())
                vector = (fractional[j] + bond_images[i,j] - fractional[i]) @ basis
                conf.SetAtomPosition(neighbor.GetIdx(), vector.tolist())
            probe.AddConformer(conf)
            Chem.AssignAtomChiralTagsFromStructure(probe)
            assert probe.GetAtomWithIdx(atom.GetIdx()).GetChiralTag() == atom.GetChiralTag()
            verified += 1
    assert verified == 6


def finite(elements=('C', 'O'), positions=((0,0,0), (1.4,0,0))):
    return {'periodic': False, 'atoms': [dict(element=e, x=x, y=y, z=z) for e,(x,y,z) in zip(elements, positions)],
            'bonds': [{'a':0,'b':1,'order':1}] if len(elements)>1 else []}


def test_diatomic_absolute_debye_intensity_matches_analytic_equation():
    result = finite_powder(finite(), two_theta_max=12, step=.2)
    q = np.asarray(result['profile']['Q_inverse_angstrom'])
    f = atomic_form_factor('C', (q/(4*np.pi))**2)
    g = atomic_form_factor('O', (q/(4*np.pi))**2)
    expected = f*f + g*g + 2*f*g*np.sin(q*1.4)/(q*1.4)
    assert result['profile']['absolute_intensity_electron_squared'] == pytest.approx(expected, rel=1e-12)
    assert result['periodic'] is False and 'peaks' not in result
    assert 'intensity_electron_squared' in result['exports']['csv']


def test_single_atom_uses_only_self_scattering():
    result = finite_powder(finite(('Cl',), ((0,0,0),)))
    q = np.asarray(result['profile']['Q_inverse_angstrom'])
    assert result['profile']['absolute_intensity_electron_squared'] == pytest.approx(atomic_form_factor('Cl',(q/(4*np.pi))**2)**2)


def test_finite_scattering_rigid_motion_and_order_invariance():
    original = finite(('C','O','N'), ((0,0,0),(1.4,0,0),(0,2,1)))
    moved = finite(('N','O','C'), ((9,0,5),(11,1.4,4),(11,0,4)))
    # Rotation (x,y,z)->(-y,x,z), then translation (11,0,4).
    a=finite_powder(original);b=finite_powder(moved)
    assert a['profile']['absolute_intensity_electron_squared'] == pytest.approx(b['profile']['absolute_intensity_electron_squared'])


def test_continuation_not_carbon_and_no_fictitious_periodic_density():
    model=finite(('C','*'), ((0,0,0),(1.4,0,0)))
    report=inspect_finite(model)
    assert report['displayed_atom_inventory']=='C' and report['continuation_site_count']==1
    assert report['known_real_bond_count']==0 and report['known_graph_components']==1
    assert 'density' in report['unavailable'] and 'density_g_cm3' not in report
    result=finite_powder(model)
    assert result['parameters']['atom_count']==1
    assert result['profile']==finite_powder(finite(('C',), ((0,0,0),)))['profile']


@pytest.mark.parametrize('change', ['periodic','nan','element','bond','duplicate','null_bonds','real_continuation'])
def test_finite_rejects_malformed_models(change):
    model=finite()
    if change=='periodic':model['periodic']=True
    if change=='nan':model['atoms'][0]['x']=float('nan')
    if change=='element':model['atoms'][0]['element']='Invalid'
    if change=='bond':model['bonds'][0]['b']=3
    if change=='duplicate':model['bonds']*=2
    if change=='null_bonds':model['bonds']=None
    if change=='real_continuation':model['atoms'][0]['continuation_site']=True
    with pytest.raises(ValueError):inspect_finite(model)


def layer():
    root=Path(__file__).resolve().parents[1]
    return import_cif((root/'data/crystal-models/hcb-imine-aa.cif').read_text())


@pytest.mark.parametrize('mode,count', [('AA',1),('AB',2),('ABC',3),('slipped',1)])
def test_periodic_stack_composition_cell_and_charge_preserved(mode,count):
    source=layer()
    for i,atom in enumerate(source['atoms']):atom['charge']=(-1 if i%2 else 1)*.1
    before=deepcopy(source)
    result=periodic_stacking(source,stacking=mode,spacing=3.6,slip=[1.5,.2])
    out=result['structure']
    assert source==before
    assert len(out['atoms'])==len(source['atoms'])*count
    assert sum(a['charge'] for a in out['atoms'])==pytest.approx(0,abs=1e-12)
    c=np.asarray(out['cell']);normal=np.cross(c[0],c[1]);normal/=np.linalg.norm(normal)
    assert np.dot(c[2],normal)==pytest.approx(count*3.6)
    assert len(import_cif(result['exports']['cif'])['atoms'])==len(out['atoms'])
    if mode=='slipped':assert np.linalg.norm(c[2]-np.dot(c[2],normal)*normal)==pytest.approx(np.hypot(1.5,.2))


def test_zero_slip_ab_is_aa_supercell_with_identical_normalized_pxrd():
    source=layer()
    a=periodic_stacking(source,stacking='AA',spacing=3.4)['structure']
    b=periodic_stacking(source,stacking='AB',spacing=3.4,slip=[0,0])['structure']
    pa=powder_xrd(a,two_theta_max=25);pb=powder_xrd(b,two_theta_max=25)
    assert pa['profile']['intensity']==pytest.approx(pb['profile']['intensity'],abs=1e-6)


def test_restack_rejects_3d_multilayer_and_arbitrary_twist():
    root=Path(__file__).resolve().parents[1]
    for filename in ('dia-imine.cif','hcb-imine-ab.cif'):
        source=import_cif((root/'data/crystal-models'/filename).read_text())
        with pytest.raises(ValueError,match='one near-planar layer'):periodic_stacking(source)
    with pytest.raises(ValueError,match='arbitrary twist'):periodic_stacking(layer(),stacking='twisted')


def test_actual_viewer_model_api_checks_scattering_and_saved_reports(tmp_path):
    from fastapi.testclient import TestClient
    from cof_vision.cofdraw_desktop import create_desktop_app
    with TestClient(create_desktop_app(tmp_path)) as client:
        prefix='/api/drawing/science/'
        model=client.post(prefix+'monomer',json={'smiles':'CO','conformers':1}).json()['viewer']
        checked=client.post(prefix+'finite-inspect',json={'model':model})
        assert checked.status_code==200,checked.text
        report=checked.json()
        assert report['atom_count']==6 and report['displayed_atom_inventory']=='CH4O'
        assert client.get(report['downloads']['json']).json()['model_sha256']==report['model_sha256']
        response=client.post(prefix+'finite-powder',json={'model':model,'two_theta_max':10})
        assert response.status_code==200,response.text
        result=response.json()
        assert result['kind']=='finite_debye_powder_scattering'
        assert client.get(result['downloads']['csv']).text.startswith('two_theta_degrees,Q_inverse_angstrom')
        bad=client.post(prefix+'pxrd',json={'structure':model})
        assert bad.status_code==422
        stacked=client.post(prefix+'stacking',json={'structure':layer(),'stacking':'AB'})
        assert stacked.status_code==200,stacked.text
        assert len(stacked.json()['structure']['atoms'])==132
        assert client.get(stacked.json()['downloads']['cif']).status_code==200
        operations=client.get(prefix+'capabilities').json()['operations']
        assert {'finite-inspect','finite-powder','stacking'} <= set(operations)
