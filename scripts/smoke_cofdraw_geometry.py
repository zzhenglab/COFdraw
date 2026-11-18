"""Check drawing geometry through the running app, including frozen releases."""
from copy import deepcopy
import hashlib
import json
from pathlib import Path
import sys

import numpy as np
from rdkit import Chem

TP = 'O=Cc1c(O)c(C=O)c(O)c(C=O)c1O'
BPDA = 'O=Cc1ccc(-c2ccc(C=O)cc2)cc1'
TAPB = 'Nc1ccc(-c2cc(-c3ccc(N)cc3)cc(-c3ccc(N)cc3)c2)cc1'


def _mol_records(block):
    """Read V2000 records without rerunning RDKit stereo perception on a pore."""
    lines = block.splitlines()
    assert 'V2000' in lines[3]
    atoms, bonds = int(lines[3][:3]), int(lines[3][3:6])
    atom_rows, bond_rows = lines[4:4+atoms], lines[4+atoms:4+atoms+bonds]
    coordinates = np.asarray([[float(line[:10]), float(line[10:20])] for line in atom_rows])
    # Atom records retain element, charge/isotope/stereo/H/map flags. Bond
    # records and trailing M records retain connectivity, order and properties.
    graph = ([line[30:] for line in atom_rows], bond_rows, lines[4+atoms+bonds:])
    return atoms, coordinates, graph


def _long_edges(recipe, mol, audit, xy):
    origins = audit['atom_origin']
    lookup = {(o['instance'], o['fragment_atom_index']): i for i, o in enumerate(origins)}
    owner = [o['instance'] for o in origins]
    rings = []
    for instance in recipe['instances']:
        fragment = Chem.MolFromSmiles(recipe['components'][instance['component']]['fragment_smiles'])
        for ring in Chem.GetSymmSSSR(fragment):
            keys = [(instance['id'], int(i)) for i in ring]
            if len(ring) == 6 and all(k in lookup for k in keys):
                rings.append([lookup[k] for k in keys])
    offsets, errors = [], []
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
        errors.extend(float(np.degrees(np.arccos(np.clip(abs(a @ b), -1, 1))))
                      for a in axes for b in axes)
        centers = np.asarray([xy[r].mean(0) for r in four])
        normal = np.asarray([-axes[0][1], axes[0][0]])
        offsets.append(float(np.sqrt(np.mean(((centers-centers.mean(0)) @ normal)**2))))
    assert len(offsets) == 6 and max(errors) < .01 and np.ptp(offsets) < .01
    return {'long_edges_checked': 6, 'max_phenyl_axis_error_degrees': max(errors),
            'equivalent_edge_offset_spread': float(np.ptp(offsets))}


def check_geometry_examples(api, download, output):
    """Use API-generated coordinates; preserve graphs, rings and source anchors."""
    # The independent, unrefined graph oracle is built without rendering a
    # deliberately overlapping baseline. All fitted outputs come from the API.
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'src'))
    from cof_vision.framework_drawing import build_patch
    output = Path(output); output.mkdir(parents=True, exist_ok=True)
    reports = []
    for name, monomer, template, count, atoms in (
            ('tp-tapb', TP, 'hcb_3_3', 6, 123),
            ('bpda-tapb', BPDA, 'hcb_3_2', 12, 252)):
        directory = output/name; directory.mkdir(parents=True, exist_ok=True)
        recipe = api('/api/drawing/prepare', {'monomer_a': monomer, 'monomer_b': TAPB,
                     'linkage': 'imine', 'template_id': template, 'instances': count,
                     'label_a': 'aldehyde', 'label_b': 'TAPB'})['recipe']
        (directory/'prepared-recipe.json').write_text(json.dumps(recipe, indent=2)+'\n', encoding='utf8')
        assert recipe['layout_refinement'] == 'bond_angle_pose_v2'
        result = api('/api/drawing/render', {'recipe': recipe, 'style_id': 'standard', 'size': 1000})
        audit = result['audit']; xy = np.asarray(audit['drawing_style']['layout']['atom_coordinates'])
        assert audit['atoms'] == atoms and xy.shape == (atoms, 2) and np.isfinite(xy).all()
        assert audit['layout_refinement']['method'] == 'bond_angle_pose_v2'
        hashes = {}
        for ext in ('png', 'svg', 'mol', 'json'):
            payload = download(result['downloads'][ext])
            (directory/f'{name}.{ext}').write_bytes(payload)
            hashes[ext] = hashlib.sha256(payload).hexdigest()
        exported_atoms, exported_xy, exported_graph = _mol_records((directory/f'{name}.mol').read_text(encoding='utf8'))
        assert exported_atoms == atoms and np.allclose(exported_xy, xy, rtol=0, atol=.000051)
        manual = deepcopy(recipe); manual['layout_refinement'] = 'manual'
        unrefined, baseline_audit = build_patch(manual)
        assert exported_graph == _mol_records(Chem.MolToMolBlock(unrefined))[2]
        mol = unrefined
        assert audit['atom_origin'] == baseline_audit['atom_origin']
        old = unrefined.GetConformer().GetPositions()[:, :2]
        owner = [o['instance'] for o in audit['atom_origin']]
        for instance in recipe['instances']:
            indices = [i for i, block in enumerate(owner) if block == instance['id'] and mol.GetAtomWithIdx(i).GetAtomicNum()]
            a, b = old[indices], xy[indices]
            assert np.allclose(np.linalg.norm(a[:, None]-a[None, :], axis=2),
                               np.linalg.norm(b[:, None]-b[None, :], axis=2), atol=1e-8)
        lengths = np.asarray([np.linalg.norm(xy[b.GetBeginAtomIdx()]-xy[b.GetEndAtomIdx()]) for b in mol.GetBonds()])
        fraction_error = float(max(abs(lengths/np.median(lengths)-1)))
        errors = []
        for bond in mol.GetBonds():
            a, b = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
            if owner[a] == owner[b]: continue
            for center, external in ((a, b), (b, a)):
                for neighbor in mol.GetAtomWithIdx(center).GetNeighbors():
                    internal = neighbor.GetIdx()
                    if internal == external or owner[internal] != owner[center]: continue
                    u, v = xy[internal]-xy[center], xy[external]-xy[center]
                    angle = np.degrees(np.arccos(np.clip(u @ v/np.linalg.norm(u)/np.linalg.norm(v), -1, 1)))
                    errors.append(float(abs(angle-120)))
        assert len(errors) in (12, 24) and max(errors) < .01 and fraction_error < .001
        diagnostic = audit['layout_diagnostics']
        assert diagnostic['nonbonded_overlap_pairs'] == [] and diagnostic['nonincident_bond_crossings'] == []
        row = {'example': name, 'atoms': atoms, 'method': 'bond_angle_pose_v2',
               'max_junction_angle_error_degrees': max(errors), 'maximum_bond_fraction_error': fraction_error,
               'nonbonded_overlaps': 0, 'nonincident_crossings': 0, 'graph_preserved': True,
               'real_intrablock_distances_preserved': True, 'downloads_sha256': hashes}
        if monomer == BPDA: row.update(_long_edges(recipe, mol, audit, xy))
        (directory/'api-result.json').write_text(json.dumps(result, indent=2)+'\n', encoding='utf8')
        reports.append(row)
    source = api('/api/drawing/sources/tfb-tab-pore-v1')['recipe']
    expected = {(row['instance'], row['fragment_atom_index']): row['xy'] for row in source['atom_coordinates']}
    assert len(expected) == 108
    source['layout_refinement'] = 'bond_angle_pose_v2'
    for style in ('standard', 'source', 'pastel_blocks'):
        result = api('/api/drawing/render', {'recipe': source, 'style_id': style, 'size': 600})
        audit = result['audit']
        wanted = [expected[(o['instance'], o['fragment_atom_index'])] for o in audit['atom_origin']]
        assert audit['drawing_style']['layout']['atom_coordinates'] == wanted
        assert audit['layout_refinement']['method'] == 'explicit_atom_coordinates'
    reports.append({'example': 'tfb-tab-source-anchors', 'preserved_atom_coordinates': 108,
                    'styles_checked': ['standard', 'source', 'pastel_blocks']})
    (output/'geometry-report.json').write_text(json.dumps(reports, indent=2)+'\n', encoding='utf8')
    return reports
