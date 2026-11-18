"""Conventional porphyrin geometry without graph edits or coordinate takeover."""
from copy import deepcopy
import json
from pathlib import Path

import numpy as np
import pytest
from rdkit import Chem
from rdkit.Chem import rdDepictor

from cof_vision.cofdraw_chemistry import prepare_template
from cof_vision.cofdraw_core import normalize_monomer
from cof_vision.cofdraw_depiction import compute_2d_coords, has_free_base_porphyrin
from cof_vision.cofdraw_layout import clean_layout
from cof_vision.depiction_rings import symmetric_depiction_rings
from cof_vision.framework_drawing import VERSION, build_patch


ROWS = {r['id']: r for r in json.loads((Path(__file__).resolve().parents[1]
        / 'data/framework-drawing/catalog/current.json').read_text())['monomer_presets']}


def graph_state(mol):
    return ([tuple((a.GetAtomicNum(), a.GetFormalCharge(), a.GetIsotope(),
                    a.GetTotalNumHs(), int(a.GetChiralTag()))) for a in mol.GetAtoms()],
            [(b.GetBeginAtomIdx(), b.GetEndAtomIdx(), b.GetBondTypeAsDouble(),
              int(b.GetStereo()), tuple(b.GetStereoAtoms())) for b in mol.GetBonds()])


def assert_pyrroles(mol, expected):
    """Geometric checks use the actual coordinates, not a template flag."""
    xy = mol.GetConformer().GetPositions()[:, :2]
    rings = [list(r) for r in symmetric_depiction_rings(mol) if len(r) == 5
             and sorted(mol.GetAtomWithIdx(i).GetAtomicNum() for i in r) == [6, 6, 6, 6, 7]]
    assert len(rings) == expected
    for ring in rings:
        for j, atom in enumerate(ring):
            u, v = xy[ring[j-1]]-xy[atom], xy[ring[(j+1) % 5]]-xy[atom]
            angle = np.degrees(np.arccos(np.clip(u @ v / (np.linalg.norm(u)*np.linalg.norm(v)), -1, 1)))
            assert 107 < angle < 109  # Default depiction included 44-degree pyrrole angles.
            assert 1.48 < np.linalg.norm(u) < 1.52


@pytest.mark.parametrize('identifier', ['tapp', 'por_cho', 'tbpp'])
def test_previews_have_regular_pyrroles_and_outward_phenyls(identifier):
    mol = Chem.MolFromSmiles(ROWS[identifier]['smiles'])
    # Isotopic substitution and a different input atom order must not require a
    # catalog identity. One representative variation suffices here.
    if identifier == 'tapp':
        mol.GetAtomWithIdx(1).SetIsotope(13)
        mol = Chem.RenumberAtoms(mol, list(reversed(range(mol.GetNumAtoms()))))
    before = graph_state(mol)
    compute_2d_coords(mol)
    assert graph_state(mol) == before
    assert_pyrroles(mol, 4)
    result = normalize_monomer(smiles=Chem.MolToSmiles(mol))
    restored = Chem.MolFromMolBlock(result['molblock'], removeHs=False)
    assert Chem.MolToSmiles(restored) == Chem.MolToSmiles(mol)
    assert_pyrroles(restored, 4)
    xy = restored.GetConformer().GetPositions()[:, :2]
    rings = [list(r) for r in symmetric_depiction_rings(restored)]
    core = set(i for r in rings if len(r) == 5 for i in r)
    center = xy[list(core)].mean(0)
    axes = []
    for ring in (r for r in rings if len(r) == 6):
        phenyl_center = xy[ring].mean(0)
        connectors = [(i, n.GetIdx()) for i in ring for n in restored.GetAtomWithIdx(i).GetNeighbors()
                      if n.GetIdx() not in ring and sum(x.GetIdx() in core for x in n.GetNeighbors()) == 2]
        assert len(connectors) == 1
        aryl, meso = connectors[0]
        radial, outward = phenyl_center-center, xy[aryl]-xy[meso]
        assert radial @ outward / np.linalg.norm(radial) / np.linalg.norm(outward) > .999
        axes.append(radial / np.linalg.norm(radial))
    assert len(axes) == 4
    # Four outward axes form an orthogonal cross (rotation is unconstrained).
    dots = sorted(abs(float(a @ b)) for i, a in enumerate(axes) for b in axes[i+1:])
    assert max(dots[:4]) < .02 and min(dots[4:]) > .999


def test_scope_excludes_other_macrocycles_and_n_substituted_cores():
    original = Chem.MolFromSmiles(ROWS['tapp']['smiles'])
    edited = Chem.RWMol(original)
    nitrogen = next(a for a in edited.GetAtoms() if a.GetAtomicNum() == 7 and a.GetIsAromatic() and a.GetTotalNumHs())
    nitrogen.SetNumExplicitHs(0)
    edited.AddBond(nitrogen.GetIdx(), edited.AddAtom(Chem.Atom('C')), Chem.BondType.SINGLE)
    substituted = edited.GetMol(); Chem.SanitizeMol(substituted)
    for mol in [substituted, Chem.MolFromSmiles('c1cc[nH]c1'),
                Chem.MolFromSmiles('NC1CCCC(N)CCCC(N)CCCC(N)CCC1')]:
        assert not has_free_base_porphyrin(mol)
        plain = Chem.Mol(mol)
        rdDepictor.Compute2DCoords(plain, canonOrient=False, forceRDKit=True)
        compute_2d_coords(mol, canonOrient=False, forceRDKit=True)
        np.testing.assert_array_equal(mol.GetConformer().GetPositions(), plain.GetConformer().GetPositions())


@pytest.mark.parametrize('a,b', [('ta', 'tapp'), ('por_cho', 'pa')])
def test_assembled_design_keeps_regular_cores_and_manual_coordinates(a, b):
    recipe = prepare_template(ROWS[a]['smiles'], ROWS[b]['smiles'], 'imine', a, b, 6, 'sql_4_2')
    recipe.update(schema_version=VERSION, recipe_id='porphyrin-depiction', designed_variant=True,
                  layout_refinement='bond_angle_pose_v2')
    mol, audit = build_patch(recipe)
    assert_pyrroles(mol, 16)
    assert audit['layout_refinement']['quality'] == 'regular'
    assert audit['layout_refinement']['final_max_junction_angle_deviation_degrees'] < 1
    assert not audit['layout_diagnostics']['nonbonded_overlap_pairs']
    assert not audit['layout_diagnostics']['nonincident_bond_crossings']
    assert sum(a.GetTotalNumHs() for a in mol.GetAtoms() if a.GetAtomicNum() == 7) == 8
    assert Chem.GetFormalCharge(mol) == 0
    # A user's deliberately different layout is authoritative for both source
    # and design recipes. The template must never silently clean these rows.
    xy = mol.GetConformer().GetPositions()[:, :2] @ np.array([[1.2, .2], [0, .8]]) + [4, -7]
    fixed = deepcopy(recipe)
    fixed['atom_coordinates'] = [dict(instance=o['instance'], fragment_atom_index=o['fragment_atom_index'], xy=p.tolist())
                                 for o, p in zip(audit['atom_origin'], xy)]
    for designed in (False, True):
        fixed['designed_variant'] = designed
        restored, report = build_patch(fixed)
        np.testing.assert_array_equal(restored.GetConformer().GetPositions()[:, :2], xy)
        assert graph_state(restored) == graph_state(mol)
        assert report['layout_refinement']['source_coordinates_preserved'] is True
    if a == 'ta':
        # The explicit Clean layout action is the supported upgrade path for
        # an old saved pore, while an ordinary render above retains its edits.
        fixed['designed_variant'] = False
        fixed['cofex_id'] = 'source-identity-must-not-follow-cleanup'
        upgraded = clean_layout(fixed)
        clean_mol, _ = build_patch(upgraded['recipe'])
        assert_pyrroles(clean_mol, 16)
        assert graph_state(clean_mol) == graph_state(mol)
        assert not upgraded['recipe'].get('cofex_id')
        assert upgraded['audit']['chemistry_changed'] is False


def test_imported_monomer_coordinates_change_only_on_explicit_cleanup():
    mol = Chem.MolFromSmiles(ROWS['por_cho']['smiles'])
    rdDepictor.Compute2DCoords(mol)  # An existing, deliberately non-template drawing.
    block = Chem.MolToMolBlock(mol)
    expected = Chem.MolFromMolBlock(block).GetConformer().GetPositions()
    untouched = normalize_monomer(molblock=block)
    np.testing.assert_array_equal(Chem.MolFromMolBlock(untouched['molblock']).GetConformer().GetPositions(), expected)
    graph = untouched['graph']
    from_sketch = normalize_monomer(graph=graph)
    np.testing.assert_array_equal(Chem.MolFromMolBlock(from_sketch['molblock']).GetConformer().GetPositions(), expected)
    cleaned = normalize_monomer(molblock=block, clean_layout=True)
    assert_pyrroles(Chem.MolFromMolBlock(cleaned['molblock']), 4)
