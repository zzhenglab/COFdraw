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
