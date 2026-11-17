"""Rigid 2D block layout and diagnostics; never a force field or crystal model."""
from __future__ import annotations

from copy import deepcopy
import math
import numpy as np


def layout_diagnostics(mol, coordinates, origins):
    xy = np.asarray(coordinates, dtype=float)[:, :2]
    block = [o['instance'] for o in origins]
    bonds = [(b.GetBeginAtomIdx(), b.GetEndAtomIdx()) for b in mol.GetBonds()]
    lengths = np.array([np.linalg.norm(xy[a]-xy[b]) for a, b in bonds])
    median = float(np.median(lengths)) if len(lengths) else 0.
    joined = [(a, b) for a, b in bonds if block[a] != block[b]]
    junction_angles = []
    for a, b in joined:
        for center, outside in ((a, b), (b, a)):
            atom = mol.GetAtomWithIdx(center)
            ideal = 180. if str(atom.GetHybridization()) == 'SP' else 109.5 if str(atom.GetHybridization()) == 'SP3' else 120.
            for neighbor in atom.GetNeighbors():
                inside = neighbor.GetIdx()
                if inside == outside or block[inside] != block[center]:
                    continue
                u, v = xy[inside]-xy[center], xy[outside]-xy[center]
                denominator = np.linalg.norm(u)*np.linalg.norm(v)
                if denominator < 1e-12:
                    continue
                angle = math.degrees(math.acos(float(np.clip(np.dot(u, v)/denominator, -1., 1.))))
                junction_angles.append({'atom': center, 'internal_neighbor': inside, 'external_neighbor': outside,
                                        'angle_degrees': angle, 'reference_degrees': ideal,
                                        'deviation_degrees': abs(angle-ideal)})
    bonded = {tuple(sorted(pair)) for pair in bonds}
    close = []
    if median > 0:
        for i in range(len(xy)):
            for j in range(i+1, len(xy)):
                if (i, j) not in bonded and np.linalg.norm(xy[i]-xy[j]) < median*.35:
                    close.append([i, j])
    return {'coordinate_scope': '2D depiction; not a strain or physical energy calculation',
            'finite_coordinates': bool(np.isfinite(xy).all()),
            'median_bond_length': median,
            'bond_length_outliers_35percent': int(np.sum(abs(lengths/median-1) > .35)) if median else 0,
            'junction_bond_ratios': [float(np.linalg.norm(xy[a]-xy[b])/median) for a, b in joined] if median else [],
            'junction_angles': junction_angles,
            'nonbonded_overlap_pairs': close,
            'overlap_threshold_bond_fraction': .35,
            'rigid_blocks': len(set(block))}


def refine_rigid_poses(mol, coordinates, origins, joins):
    """Fit port directions with only whole-block rotations and translations.

    Each join supplies actual atom indices and the original outward unit vectors
    of its removed ports. No intrablock bond, ring, or stereocentre is deformed.
    Reject a fit that introduces overlaps or does not improve the constraints.
    """
    original = np.asarray(coordinates, dtype=float)
    names = list(dict.fromkeys(o['instance'] for o in origins))
    owner = np.array([names.index(o['instance']) for o in origins])
    centers = np.array([original[owner == i].mean(0) for i in range(len(names))])
    local = original-centers[owner]
    internal = [np.linalg.norm(original[b.GetBeginAtomIdx()]-original[b.GetEndAtomIdx()])
                for b in mol.GetBonds() if owner[b.GetBeginAtomIdx()] == owner[b.GetEndAtomIdx()]]
    target = float(np.median(internal))
    report = {'method': 'rigid_block_pose_v1', 'applied': False,
              'interpretation': 'Rigid 2D drawing cleanup, not molecular mechanics.',
              'target_bond_length': target}
    if not joins or not math.isfinite(target) or target <= 0:
        report['reason'] = 'no_supported_block_junctions'
        return original, report
    count = len(names)
    parameters = np.zeros((count, 3))

    def evaluate(params, derivative=False):
        angles = params[:, 2]
        rotations = np.array([[[math.cos(t), -math.sin(t)], [math.sin(t), math.cos(t)]] for t in angles])
        rotated = np.einsum('nij,nj->ni', rotations[owner], local)
        xy = rotated+centers[owner]+params[owner, :2]
        residuals, jacobian = [], []
        turn = lambda vector: np.array([-vector[1], vector[0]])
        for join in joins:
            a, b = join['a'], join['b']; ia, ib = owner[a], owner[b]
            da = rotations[ia] @ join['direction_a']; db = rotations[ib] @ join['direction_b']
            residuals.extend(xy[b]-xy[a]-.5*target*(da-db))
            residuals.extend(.45*target*(da+db))
            if derivative:
                rows = np.zeros((4, count*3))
                rows[:2, ia*3:ia*3+2] = -np.eye(2)
                rows[:2, ib*3:ib*3+2] = np.eye(2)
                rows[:2, ia*3+2] = -turn(rotated[a])-.5*target*turn(da)
                rows[:2, ib*3+2] = turn(rotated[b])+.5*target*turn(db)
                rows[2:, ia*3+2] = .45*target*turn(da)
                rows[2:, ib*3+2] = .45*target*turn(db)
                jacobian.extend(rows)
        # Weak anchors retain the requested overall arrangement while allowing
        # local junction cleanup. Angles are radians, translations bond units.
        weights = np.tile([.008, .008, .025*target], count)
        residuals.extend(weights*params.ravel())
        if derivative:
            jacobian.extend(np.diag(weights))
        return xy, np.asarray(residuals), np.asarray(jacobian) if derivative else None

    _, initial, _ = evaluate(parameters)
    best_score = float(initial @ initial)
    report['initial_constraint_error'] = best_score
    for iteration in range(60):
        _, residual, jacobian = evaluate(parameters, True)
        delta = np.linalg.lstsq(jacobian, -residual, rcond=None)[0].reshape(count, 3)
        delta[:, 2] = np.clip(delta[:, 2], -.18, .18)
        accepted = False
        for rate in (1., .5, .2, .05):
            trial = parameters+rate*delta
            _, candidate, _ = evaluate(trial)
            score = float(candidate @ candidate)
            if score < best_score-1e-10:
                parameters, best_score = trial, score
                accepted = True
                break
        if not accepted or float(np.max(np.abs(delta))) < 1e-7:
            break
    fitted, _, _ = evaluate(parameters)
    report.update(iterations=iteration+1, final_constraint_error=best_score)
    before = layout_diagnostics(mol, original, origins)
    after = layout_diagnostics(mol, fitted, origins)
    new_overlaps = {tuple(p) for p in after['nonbonded_overlap_pairs']}-{tuple(p) for p in before['nonbonded_overlap_pairs']}
    if new_overlaps:
        report.update(reason='new_nonbonded_overlap', new_overlap_pairs=sorted(new_overlaps))
        return original, report
    if best_score >= report['initial_constraint_error']-1e-8:
        report['reason'] = 'layout_already_satisfies_constraints'
        return original, report
    report.update(applied=True, block_transforms={name: {'translation': parameters[i, :2].tolist(),
                  'rotation_degrees': float(np.degrees(parameters[i, 2]))} for i, name in enumerate(names)})
    return fitted, report


def clean_layout(recipe):
    """Return an explicit clean-layout design, with its exact fitted coordinates."""
    from .framework_drawing import build_patch
    if recipe.get('explicit_graph'):
        raise ValueError('This source-registered graph keeps its supplied atom coordinates; rigid cleanup requires an assembled monomer recipe')
    cleaned = deepcopy(recipe)
    cleaned.pop('atom_coordinates', None)
    cleaned['layout_refinement'] = 'rigid_block_pose_v1'
    cleaned['designed_variant'] = True
    for field in ('cofex_id', 'cofex_string', 'cofex_notation', 'cofex_source_records', 'literature_metadata'):
        cleaned.pop(field, None)
    mol, audit = build_patch(cleaned)
    cleaned['atom_coordinates'] = [{**origin, 'xy': [float(p[0]), float(p[1])]} for origin, p in
                                   zip(audit['atom_origin'], mol.GetConformer().GetPositions())]
    cleaned['layout_refinement'] = 'manual'
    return {'recipe': cleaned, 'audit': {'layout_refinement': audit['layout_refinement'],
                                        'layout_diagnostics': audit['layout_diagnostics'],
                                        'chemistry_changed': False,
                                        'scope': 'Explicit 2D coordinate cleanup; source identity cleared for the edited design.'}}
