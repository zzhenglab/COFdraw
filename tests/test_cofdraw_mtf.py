"""The COF-432 reference drawing has a large square surrounded by sixteen small squares."""
from collections import Counter, defaultdict
import json
from pathlib import Path

import numpy as np
import pytest
from fastapi.testclient import TestClient
from rdkit import Chem

from cof_vision.cofdraw_chemistry import design_signature, prepare_template
from cof_vision.cofdraw_desktop import create_desktop_app
from cof_vision.cofdraw_layout import layout_diagnostics
from cof_vision.cofdraw_periodic_options import validated_design
from cof_vision.cofdraw_repeat import repeat_options, repeat_patch


DATA = Path(__file__).resolve().parents[1] / 'data'
MONOMERS = {row['id']: row['smiles'] for row in
            json.loads((DATA / 'framework-drawing/catalog/current.json').read_text())['monomer_presets']}


def prepare(aldehyde, amine, instances):
    return prepare_template(MONOMERS[aldehyde], MONOMERS[amine], 'imine', aldehyde, amine, instances, 'mtf_4_3')


def attachment_centers(mol, origins):
    """The mean of each monomer's terminal C/N sites is independent of port numbers and cropping."""
    sites = defaultdict(list)
    xyz = mol.GetConformer().GetPositions()[:, :2]
    for atom, origin in zip(mol.GetAtoms(), origins):
        name = origin['instance']
        if atom.GetAtomicNum() and any(not n.GetAtomicNum() or origins[n.GetIdx()]['instance'] != name
                                       for n in atom.GetNeighbors()):
            sites[name].append(atom.GetIdx())
    return {name: xyz[ids].mean(axis=0) for name, ids in sites.items()}


@pytest.mark.parametrize('aldehyde,amine', [('tfb', 'tape'), ('tfb', 'tapp'), ('tfppy', 'tapb')])
def test_mtf_large_pore_has_four_corners_and_four_side_nodes(aldehyde, amine):
    # ETTA + TFB is the reported pair. TAPP and the reverse-degree pair
    # exercise applying the same net to other, explicitly designed molecules.
    recipe = prepare(aldehyde, amine, 32)
    mol, origins, _, summary, components = repeat_patch(recipe, 2)
    faces = summary['pore_faces']
    assert sorted(map(len, faces)) == [4, 8]
    assert len(set(faces[0]) & set(faces[1])) == 2
    large = max(faces, key=len)
    centers = attachment_centers(mol, origins)
    degrees = {o['instance']: components[o['component']]['connectivity'] for o in origins}
    points = np.array([centers[name] for name in large])
    edges = np.roll(points, -1, axis=0) - points
    edges /= np.linalg.norm(edges, axis=1)[:, None]
    angles = np.degrees(np.arccos(np.clip(np.sum(edges * np.roll(edges, 1, axis=0), axis=1), -1, 1)))
    corners = np.array([degrees[name] == 4 for name in large])
    assert np.all(corners != np.roll(corners, 1))
    assert np.max(abs(angles[corners] - 90)) < 15
    assert np.max(angles[~corners]) < 15
    if (aldehyde, amine) != ('tfb', 'tape'):
        # A designed analogue gets conventional imine angles. COF-432's
        # three measured projected angles are specific to its monomers.
        junctions = layout_diagnostics(mol, mol.GetConformer().GetPositions(), origins)['junction_angles']
        assert all(abs(a['angle_degrees'] - 120) < 5 for a in junctions)


def test_mtf_reference_patch_and_single_pore_selections():
    recipe = prepare('tfb', 'tape', 32)
    assert len(recipe['instances']) == 32 and len(recipe['connections']) == 48
    large_rim_distances = []
    # Euler rank distinguishes actual connected faces from overlapping rings.
    for n, mode, census in [(1, 'mtf_patch', {4: 16, 8: 1}), (1, 'small_pore', {4: 1}),
                            (1, 'large_pore', {8: 1}), (3, 'pores', {4: 2, 8: 1})]:
        mol, origins, boundaries, summary, _ = repeat_patch(recipe, n, extent_mode=mode)
        faces = summary['pore_faces']
        assert Counter(map(len, faces)) == census
        assert summary['pore_count'] == len(faces) == summary['interblock_bonds'] - summary['blocks'] + 1
        assert len(Chem.GetMolFrags(mol)) == 1 and len(origins) == mol.GetNumAtoms()
        centers = attachment_centers(mol, origins)
        for face in faces:
            if len(face) == 8:
                # Increasing drawing extent must preserve the already visible square.
                points = np.array([centers[name] for name in face])
                distances = np.linalg.norm(points[:, None] - points[None, :], axis=2)
                large_rim_distances.append(np.sort(distances[np.triu_indices(8, 1)]))
        if mode == 'mtf_patch':
            assert summary['component_counts'] == {'A': 16, 'B': 16}
            assert len(boundaries) == 16
            large = next(set(face) for face in faces if len(face) == 8)
            # Each of the eight perimeter edges borders one small square.
            assert sum(len(large.intersection(face)) == 2 for face in faces if len(face) == 4) == 8
            # Fig. 2b places four small-square clusters around a central diamond.
            # Compare areas in drawing coordinates, independent of canvas scale.
            areas = defaultdict(list)
            for face in faces:
                points = np.array([centers[name] for name in face])
                x, y = (points - points.mean(axis=0)).T
                areas[len(face)].append(abs(x @ np.roll(y, -1) - y @ np.roll(x, -1)) / 2)
                if len(face) == 8:
                    # The large square is shown on its point in the published view.
                    edges = np.roll(points, -1, axis=0) - points
                    directions = np.degrees(np.arctan2(edges[:, 1], edges[:, 0])) % 90
                    assert np.max(abs(directions - 45)) < 10
            assert 3 < areas[8][0] / np.mean(areas[4]) < 5
    for distances in large_rim_distances[1:]:
        np.testing.assert_allclose(distances, large_rim_distances[0], atol=1e-5)
    modes = {e.get('extent_mode') for e in repeat_options(recipe)['extents']}
    assert {'small_pore', 'large_pore'} <= modes
    # Older users explicitly selecting four monomers keep the single small square.
    legacy = prepare('tfb', 'tape', 4)
    assert len(legacy['instances']) == 4
    assert 'mtf_patch' in {e.get('extent_mode') for e in repeat_options(legacy)['extents']}


def test_saved_four_block_mtf_incidence_still_expands():
    # Instance order and port numbers from a project saved by the previous
    # mtf generator. A saved project keeps its local molecular graph.
    recipe = prepare('tfb', 'tape', 4)
    ring = [('B1', 3, 2), ('A2', 2, 1), ('B4', 2, 4), ('A3', 3, 2)]
    recipe['instances'] = [dict(id=name, component=name[0], incoming_port=p, outgoing_port=q)
                           for name, p, q in ring]
    recipe['connections'] = [
        dict(a=dict(instance=a, port=q), b=dict(instance=b, port=p), bond_order=2)
        for (a, _, q), (b, p, _) in zip(ring, ring[1:] + ring[:1])]
    recipe['design_template_signature'] = design_signature(recipe)
    restored, _ = validated_design(recipe)
    assert restored['connections'] == recipe['connections']
    _, _, _, summary, _ = repeat_patch(recipe, 1, extent_mode='large_pore')
    assert list(map(len, summary['pore_faces'])) == [8]


def test_desktop_mtf_exports_the_displayed_imine_graph(tmp_path):
    with TestClient(create_desktop_app(DATA, output_root=tmp_path)) as api:
        prepared = api.post('/api/drawing/prepare', json=dict(
            monomer_a=MONOMERS['tfb'], monomer_b=MONOMERS['tape'], linkage='imine',
            template_id='mtf_4_3', instances=32))
        prepared.raise_for_status()
        response = api.post('/api/drawing/render', json=dict(recipe=prepared.json()['recipe'], size=640))
        response.raise_for_status()
        urls = response.json()['downloads']
        mol = Chem.MolFromMolBlock(api.get(urls['mol']).text, sanitize=False, removeHs=False)
        audit = api.get(urls['json']).json()['audit']
        assert api.get(urls['png']).content.startswith(b'\x89PNG')
        origins = audit['atom_origin']
        assert mol.GetNumAtoms() == len(origins) == len(audit['drawing_style']['layout']['atom_pixels'])
        joins = [b for b in mol.GetBonds() if
                 origins[b.GetBeginAtomIdx()]['instance'] != origins[b.GetEndAtomIdx()]['instance']]
        assert len(joins) == 48
        assert all(b.GetBondType() == Chem.BondType.DOUBLE and
                   {b.GetBeginAtom().GetAtomicNum(), b.GetEndAtom().GetAtomicNum()} == {6, 7} for b in joins)
        xyz = mol.GetConformer().GetPositions()[:, :2]
        geometry = layout_diagnostics(mol, xyz, origins)
        assert not geometry['nonincident_bond_crossings']
        assert not geometry['nonbonded_overlap_pairs']
        assert all(.7 < ratio < 1.35 for ratio in geometry['junction_bond_ratios'])
        assert len(geometry['junction_angles']) == 2 * len(joins)
        assert all(95 < a['angle_degrees'] < 155 for a in geometry['junction_angles'])
        # The three imines on each TFB in the ab projection of COF-432
        # bend two ways (SI Table S1, C1/N1, C19/N2 and C6/N3).
        # Rotating all three alike, as in hcb, would erase its T-shaped exit.
        carbon_turns = defaultdict(list)
        projected_bends = defaultdict(list)
        angles = {a['atom']: a['angle_degrees'] for a in geometry['junction_angles']}
        for bond in joins:
            sites = (bond.GetBeginAtom(), bond.GetEndAtom())
            carbon, nitrogen = sorted(sites, key=lambda atom: atom.GetAtomicNum())
            i, j = carbon.GetIdx(), nitrogen.GetIdx()
            k = next(a.GetIdx() for a in carbon.GetNeighbors() if a.GetIdx() != j)
            u, v = xyz[i] - xyz[k], xyz[j] - xyz[i]
            turn = np.degrees(np.arctan2(np.linalg.det([u, v]), u @ v))
            carbon_turns[origins[i]['instance']].append(turn)
            projected_bends[origins[i]['instance']].append((angles[i], angles[j]))
        assert len(carbon_turns) == 16
        for name, turns in carbon_turns.items():
            assert len(turns) == 3 and min(turns) < -10 and max(turns) > 10
            np.testing.assert_allclose(sorted(projected_bends[name]),
                                       [(123, 121), (128, 142), (133, 143)], atol=2)
