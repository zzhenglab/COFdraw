"""Geometry regressions from six-block and long aromatic-linker pore drawings.

These are finite user designs. Angles below are 2D drawing constraints, not
force-field energies, synthesized-material claims, or assigned imine E/Z.
"""
from copy import deepcopy
import json
from pathlib import Path

import numpy as np
import pytest
from rdkit import Chem

from cof_vision.cofdraw_chemistry import prepare_template
from cof_vision.cofdraw_core import render_styled
from cof_vision.framework_drawing import VERSION, build_patch


TP = 'O=Cc1c(O)c(C=O)c(O)c(C=O)c1O'
BPDA = 'O=Cc1ccc(-c2ccc(C=O)cc2)cc1'
TAPB = 'Nc1ccc(-c2cc(-c3ccc(N)cc3)cc(-c3ccc(N)cc3)c2)cc1'
ROOT = Path(__file__).resolve().parents[1]


def _design(smiles, template, count):
    return {'schema_version': VERSION, 'recipe_id': 'geometry-regression',
            'designed_variant': True, 'topology': None,
            **prepare_template(monomer_a=smiles, monomer_b=TAPB,
                               linkage='imine', label_a='aldehyde', label_b='TAPB',
                               instances=count, template_id=template)}


def _graph(mol):
    # Include properties that a coordinate cleanup must never silently change.
    atoms = [(a.GetAtomicNum(), a.GetIsotope(), a.GetFormalCharge(),
              a.GetNumExplicitHs(), a.GetNoImplicit(), a.GetIsAromatic(),
              int(a.GetChiralTag()), a.GetAtomMapNum()) for a in mol.GetAtoms()]
    bonds = [(b.GetBeginAtomIdx(), b.GetEndAtomIdx(), b.GetBondTypeAsDouble(),
              b.GetIsAromatic(), int(b.GetStereo()), tuple(b.GetStereoAtoms()))
             for b in mol.GetBonds()]
    return atoms, bonds


def _junction_geometry(mol, origins):
    xyz = mol.GetConformer().GetPositions()[:, :2]
    owner = [o['instance'] for o in origins]
    lengths, errors = [], []
    for bond in mol.GetBonds():
        a, b = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        lengths.append(np.linalg.norm(xyz[a]-xyz[b]))
        if owner[a] == owner[b]:
            continue
        for atom, external in ((a, b), (b, a)):
            for neighbor in mol.GetAtomWithIdx(atom).GetNeighbors():
                internal = neighbor.GetIdx()
                if internal == external or owner[internal] != owner[atom]:
                    continue
                u, v = xyz[internal]-xyz[atom], xyz[external]-xyz[atom]
                angle = np.degrees(np.arccos(np.clip(u @ v / np.linalg.norm(u) / np.linalg.norm(v), -1, 1)))
                errors.append(abs(angle-120.0))
    return np.asarray(lengths), np.asarray(errors)


@pytest.fixture(scope='module', params=[(TP, 'hcb_3_3', 6, 123), (BPDA, 'hcb_3_2', 12, 252)],
                ids=['tp-tapb-six-block', 'bpda-tapb-long-arms'])
def pore_geometry(request):
    smiles, template, count, atoms = request.param
    recipe = _design(smiles, template, count)
    original, original_audit = build_patch({**recipe, 'layout_refinement': 'manual'})
    fitted, fitted_audit = build_patch({**recipe, 'layout_refinement': 'bond_angle_pose_v2'})
    return recipe, original, original_audit, fitted, fitted_audit, atoms


def test_angle_cleanup_preserves_graph_and_all_real_intrablock_distances(pore_geometry):
    recipe, original, before, fitted, after, atom_count = pore_geometry
    assert fitted.GetNumAtoms() == original.GetNumAtoms() == atom_count
    assert _graph(fitted) == _graph(original)
    assert after['atom_origin'] == before['atom_origin']
    x = original.GetConformer().GetPositions()[:, :2]
    y = fitted.GetConformer().GetPositions()[:, :2]
    for instance in recipe['instances']:
        indices = [i for i, origin in enumerate(before['atom_origin'])
                   if origin['instance'] == instance['id'] and original.GetAtomWithIdx(i).GetAtomicNum()]
        a, b = x[indices], y[indices]
        assert np.allclose(np.linalg.norm(a[:, None]-a[None, :], axis=2),
                           np.linalg.norm(b[:, None]-b[None, :], axis=2), atol=1e-8)
    # Tp remains the user's enol/imine design; cleanup is not tautomer conversion.
    assert [a.GetTotalNumHs() for a in original.GetAtoms() if a.GetAtomicNum() == 8] == [
        a.GetTotalNumHs() for a in fitted.GetAtoms() if a.GetAtomicNum() == 8]


def test_screenshot_pores_have_regular_planar_junctions_without_stretched_bonds(pore_geometry):
    _, _, _, fitted, audit, _ = pore_geometry
    lengths, errors = _junction_geometry(fitted, audit['atom_origin'])
    assert len(errors) in (12, 24)
    assert errors.max() < 1.0
    assert np.max(abs(lengths/np.median(lengths)-1)) < 0.01
    assert audit['layout_diagnostics']['nonbonded_overlap_pairs'] == []
    assert audit['layout_diagnostics']['nonincident_bond_crossings'] == []
    assert np.isfinite(fitted.GetConformer().GetPositions()).all()


def _component_rings(recipe, origins):
    lookup = {(o['instance'], o['fragment_atom_index']): i for i, o in enumerate(origins)}
    rings = []
    # Per-fragment ring membership avoids macrocycle-wide ring perception.
    for instance in recipe['instances']:
        fragment = Chem.MolFromSmiles(recipe['components'][instance['component']]['fragment_smiles'])
        for ring in Chem.GetSymmSSSR(fragment):
            keys = [(instance['id'], int(i)) for i in ring]
            if len(ring) == 6 and all(k in lookup for k in keys):
                rings.append([lookup[k] for k in keys])
    return rings


def test_equivalent_long_edges_have_parallel_ring_axes_and_repeatable_offsets():
    recipe = _design(BPDA, 'hcb_3_2', 12)
    mol, audit = build_patch({**recipe, 'layout_refinement': 'bond_angle_pose_v2'})
    xy = mol.GetConformer().GetPositions()[:, :2]
    owner = [o['instance'] for o in audit['atom_origin']]
    rings = _component_rings(recipe, audit['atom_origin'])
    offsets = []
    for instance in recipe['instances']:
        if recipe['components'][instance['component']]['monomer_smiles'] != BPDA:
            continue
        block = instance['id']
        four = [r for r in rings if owner[r[0]] == block]
        for bond in mol.GetBonds():
            a, b = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
            if (owner[a] == block) == (owner[b] == block):
                continue
            outside = b if owner[a] == block else a
            neighbors = [n.GetIdx() for n in mol.GetAtomWithIdx(outside).GetNeighbors()
                         if owner[n.GetIdx()] == owner[outside]]
            four.extend(r for r in rings if any(n in r for n in neighbors))
        assert len(four) == 4
        axes = []
        for ring in four:
            ends = [i for i in ring if any(n.GetIdx() not in ring for n in mol.GetAtomWithIdx(i).GetNeighbors())]
            assert len(ends) == 2
            axis = xy[ends[1]]-xy[ends[0]]
            axes.append(axis/np.linalg.norm(axis))
        # Parallel phenyl axes, not falsely straightening the 120-degree imines.
        for a in axes:
            for b in axes:
                assert abs(a @ b) > np.cos(np.radians(1))
        centers = np.asarray([xy[r].mean(0) for r in four])
        normal = np.asarray([-axes[0][1], axes[0][0]])
        offsets.append(np.sqrt(np.mean(((centers-centers.mean(0)) @ normal)**2)))
    assert len(offsets) == 6
    assert np.ptp(offsets) < 0.01  # all equivalent edges, including left/right


def test_reviewed_tfb_tab_anchors_are_unchanged_by_automatic_layout_or_styles(tmp_path):
    recipe = json.loads((ROOT/'data/framework-drawing/tfb-tab-pore-v1.recipe.json').read_text(encoding='utf-8'))
    expected = {(row['instance'], row['fragment_atom_index']): row['xy'] for row in recipe['atom_coordinates']}
    assert len(expected) == 108
    requested = deepcopy(recipe)
    requested['layout_refinement'] = 'bond_angle_pose_v2'
    mol, audit = build_patch(requested)
    xy = np.asarray([expected[(o['instance'], o['fragment_atom_index'])] for o in audit['atom_origin']])
    assert np.array_equal(mol.GetConformer().GetPositions()[:, :2], xy)
    assert audit['layout_refinement']['method'] == 'explicit_atom_coordinates'
    signatures = []
    for style in ('standard', 'source', 'pastel_blocks'):
        folder = tmp_path/style
        rendered = render_styled(requested, folder, 600, style)
        assert np.array_equal(rendered['drawing_style']['layout']['atom_coordinates'], xy)
        signatures.append((folder/(recipe['recipe_id']+'.mol')).read_bytes())
    assert signatures[0] == signatures[1] == signatures[2]


@pytest.mark.parametrize('spacer,steps,defined_azo', [('azoda', 2, 6), ('oda44', 4, 0)])
def test_hcb_imines_follow_the_equivalent_linker_poses(spacer, steps, defined_azo):
    rows = json.loads((ROOT / 'data/framework-drawing/catalog/current.json').read_bytes())['monomer_presets']
    smiles = {row['id']: row['smiles'] for row in rows}
    recipe = dict(schema_version=VERSION, recipe_id='hcb-spacer',
                  **prepare_template(smiles['tfb'], smiles[spacer], 'imine', 'TFB', spacer, 12, 'hcb_3_2'))
    # A recipe without the optional layout setting and the normal Studio
    # request must produce the same usable molecular coordinates.
    mol, audit = build_patch(recipe)
    ui, _ = build_patch({**recipe, 'designed_variant': True, 'layout_refinement': 'bond_angle_pose_v2'})
    raw, _ = build_patch({**recipe, 'layout_refinement': 'manual'})
    assert _graph(mol) == _graph(raw)
    assert np.allclose(ui.GetConformer().GetPositions(), mol.GetConformer().GetPositions())
    assert sum(b.GetStereo() == Chem.BondStereo.STEREOE for b in mol.GetBonds()) == defined_azo
    assert audit['defined_double_bond_geometry']['valid']
    _, errors = _junction_geometry(mol, audit['atom_origin'])
    assert len(errors) == 24 and max(errors) < 1

    # Check corresponding atoms in every repeated species, including the
    # outward arms. The azo spacer repeats each edge, the ether every other.
    xy = mol.GetConformer().GetPositions()[:, :2]
    xy -= xy.mean(axis=0)
    owner = np.array([r['instance'] for r in audit['atom_origin']])
    blocks = [xy[owner == item['id']] for item in recipe['instances']]
    angle = np.radians(30 * steps)
    rotation = np.array([[np.cos(angle), np.sin(angle)], [-np.sin(angle), np.cos(angle)]])
    for current, following in zip(blocks, blocks[steps:] + blocks[:steps]):
        assert np.allclose(current @ rotation, following, atol=.01)


@pytest.mark.parametrize('a,b,template,coordination', [
    ('tfb','pa','hcb_3_2',3), ('ta','tapb','hcb_3_2',3),
    ('tfb','tapb','hcb_3_3',3), ('tfb','oda44','hcb_3_2',3),
    ('ta','pytta','sql_4_2',4), ('hexaformylbenzene','pa','hxl_6_2',6),
])
def test_all_imines_at_each_equivalent_node_turn_the_same_way(a, b, template, coordination):
    rows = json.loads((ROOT/'data/framework-drawing/catalog/current.json').read_bytes())['monomer_presets']
    smiles = {row['id']: row['smiles'] for row in rows}
    smiles['hexaformylbenzene'] = 'O=Cc1c(C=O)c(C=O)c(C=O)c(C=O)c1C=O'
    recipe = prepare_template(smiles[a], smiles[b], 'imine', a, b, 12, template)
    mol, audit = build_patch(dict(schema_version=VERSION, recipe_id='node-turn', **recipe))
    xy = mol.GetConformer().GetPositions()[:, :2]
    owners = [row['instance'] for row in audit['atom_origin']]
    nodes = 0
    for block in recipe['instances']:
        if recipe['components'][block['component']]['connectivity'] != coordination:
            continue
        nodes += 1
        arms = []
        for atom in mol.GetAtoms():
            i = atom.GetIdx()
            if owners[i] != block['id'] or atom.GetDegree() != 2:
                continue
            outside = [n for n in atom.GetNeighbors() if n.GetAtomicNum()==0 or owners[n.GetIdx()]!=owners[i]]
            if len(outside) != 1:
                continue
            j = outside[0].GetIdx()
            assert mol.GetBondBetweenAtoms(i,j).GetBondType() == Chem.BondType.DOUBLE
            k = next(n.GetIdx() for n in atom.GetNeighbors() if n.GetIdx()!=j)
            radial, bend = xy[i]-xy[k], xy[j]-xy[i]
            turn = np.degrees(np.arctan2(radial[0]*bend[1]-radial[1]*bend[0], radial@bend))
            arms.append((turn, outside[0].GetAtomicNum()==0))
        assert len(arms) == coordination
        assert sum(is_stub for _,is_stub in arms) == coordination - 2
        # A rectangular pyrene four-arm node can depart from the square angle.
        tolerance = 1 if coordination == 3 else 10
        assert all(abs(turn - 60) < tolerance for turn,_ in arms)
    assert nodes == {3: 6, 4: 4, 6: 3}[coordination]


@pytest.mark.parametrize('aldehyde,amine,net,blocks,rod_count', [
    ('ta', 'pa_so3h', 'chain_2_2', 6, 6),
    ('bpda_n3', 'pytta', 'sql_4_2', 8, 4),
    ('phenyl-ta', 'tapb', 'hcb_3_2', 12, 6),
    ('hexaaldehyde', 'cv_7dd0838a18ffc4a3', 'hxl_6_2', 6, 3),
])
def test_straight_imine_linker_ends_bend_to_opposite_sides(aldehyde, amine, net, blocks, rod_count):
    presets = json.loads((ROOT / 'data/framework-drawing/catalog/current.json').read_bytes())['monomer_presets']
    smiles = {row['id']: row['smiles'] for row in presets}
    smiles['phenyl-ta'] = 'O=Cc1c(-c2ccccc2)cc(C=O)cc1'
    smiles['hexaaldehyde'] = 'O=Cc1c(C=O)c(C=O)c(C=O)c(C=O)c1C=O'
    recipe = prepare_template(smiles[aldehyde], smiles[amine], 'imine', aldehyde, amine, blocks, net)
    mol, audit = build_patch(dict(schema_version=VERSION, recipe_id='straight-linker', **recipe))
    xy = mol.GetConformer().GetPositions()[:, :2]
    owners = [row['instance'] for row in audit['atom_origin']]
    rods = [row['id'] for row in recipe['instances'] if recipe['components'][row['component']]['connectivity'] == 2]
    assert len(rods) == rod_count
    for rod in rods:
        ends = []
        for atom in mol.GetAtoms():
            site = atom.GetIdx()
            if owners[site] != rod or atom.GetAtomicNum() == 0:
                continue
            for neighbor in atom.GetNeighbors():
                tip = neighbor.GetIdx()
                if owners[tip] != rod or neighbor.GetAtomicNum() == 0:
                    assert mol.GetBondBetweenAtoms(site, tip).GetBondType() == Chem.BondType.DOUBLE
                    ends.append((site, tip))
        # Use the same directed chord for both physical sides, even when one
        # terminal bond is an outward continuation at the end of an open chain.
        assert len(ends) == 2
        sites, tips = np.asarray(ends).T
        chord = xy[sites[1]] - xy[sites[0]]
        normal = np.array([-chord[1], chord[0]]) / np.linalg.norm(chord)
        lateral = (xy[tips] - xy[sites]) @ normal
        assert lateral.min() < -0.5 and lateral.max() > 0.5
    if net == 'hxl_6_2':
        # Bulky unreacted arms of these porphyrins cross in this 2D crop.
        # The fit recovers the imine bends and eliminates overlapping atoms.
        assert audit['layout_refinement']['quality'] == 'needs_review'
        assert audit['layout_diagnostics']['nonbonded_overlap_pairs'] == []
        assert len(audit['layout_diagnostics']['nonincident_bond_crossings']) == 12
    else:
        assert audit['layout_refinement']['quality'] == 'regular'
