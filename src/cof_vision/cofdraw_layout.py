"""Rigid 2D block layout and diagnostics; never a force field or crystal model."""
from __future__ import annotations

from copy import deepcopy
import math
import numpy as np


def _bond_crossings(xy, bonds):
    """Proper intersections between nonincident 2D bond segments."""
    if not bonds:
        return []
    endpoints = np.asarray(bonds, dtype=int)
    starts = xy[endpoints[:, 0]]
    vectors = xy[endpoints[:, 1]]-starts
    cross = lambda a, b: a[..., 0]*b[..., 1]-a[..., 1]*b[..., 0]
    found = []
    for start in range(0, len(bonds), 128):
        r = vectors[start:start+128, None]
        s = vectors[None]
        delta = starts[None]-starts[start:start+128, None]
        divisor = cross(r, s)
        valid = abs(divisor) > 1e-10
        safe = np.where(valid, divisor, 1.)
        t, u = cross(delta, s)/safe, cross(delta, r)/safe
        rows, columns = np.where(valid & (t > 1e-7) & (t < 1-1e-7) & (u > 1e-7) & (u < 1-1e-7))
        found.extend([int(start+i), int(j)] for i, j in zip(rows, columns)
                     if j > start+i and not (set(bonds[start+i]) & set(bonds[j])))
    return found


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
        # Row chunks bound temporary memory and avoid a Python distance call
        # for every atom pair during candidate layout evaluation.
        for start in range(0, len(xy), 128):
            delta = xy[start:start+128, None]-xy[None]
            rows, columns = np.where(np.einsum('ijk,ijk->ij', delta, delta) < (median*.35)**2)
            close.extend([int(start+i), int(j)] for i, j in zip(rows, columns)
                         if j > start+i and (int(start+i), int(j)) not in bonded)
    return {'coordinate_scope': '2D depiction; not a strain or physical energy calculation',
            'finite_coordinates': bool(np.isfinite(xy).all()),
            'median_bond_length': median,
            'bond_length_outliers_35percent': int(np.sum(abs(lengths/median-1) > .35)) if median else 0,
            'junction_bond_ratios': [float(np.linalg.norm(xy[a]-xy[b])/median) for a, b in joined] if median else [],
            'junction_angles': junction_angles,
            'nonbonded_overlap_pairs': close,
            'nonincident_bond_crossings': _bond_crossings(xy, bonds),
            'overlap_threshold_bond_fraction': .35,
            'rigid_blocks': len(set(block))}


def refine_rigid_poses(mol, coordinates, origins, joins, *, robust=False):
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
        if robust:
            # Damp nearly singular closed-loop modes and shorten the entire
            # step, rather than clipping only angles and losing descent.
            normal = jacobian.T @ jacobian
            delta = np.linalg.solve(normal + np.eye(count*3)*1e-5, -jacobian.T @ residual).reshape(count, 3)
            delta *= min(1., .18/max(float(np.max(abs(delta[:, 2]))), 1e-12))
        else:
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


def refine_bond_angles(mol, coordinates, origins, joins):
    """Choose unspecified port-side depictions, then fit rigid molecular blocks.

    A consumed dummy is a connection placeholder, not specified imine E/Z.
    Its two planar sides are allowed only at an endpoint with one internal
    neighbour. Every real intrablock distance remains fixed. This is a 2D
    depiction fit, never a force field or stereochemical assignment.
    """
    from .framework_drawing import _refine_block_spacing
    original = np.asarray(coordinates, dtype=float)
    spaced, spacing = _refine_block_spacing(mol, original, origins)
    spaced = np.asarray(spaced)
    original_diagnostics = layout_diagnostics(mol, original, origins)
    spacing_diagnostics = layout_diagnostics(mol, spaced, origins)
    new_spacing_crossings = ({tuple(p) for p in spacing_diagnostics['nonincident_bond_crossings']}
                            - {tuple(p) for p in original_diagnostics['nonincident_bond_crossings']})
    if new_spacing_crossings:
        spaced = original.copy()
        spacing = {**spacing, 'applied': False, 'reason': 'new_bond_crossing'}
    owners = [row['instance'] for row in origins]
    choices = []
    for join in joins:
        candidates = []
        for key, direction_key in (('a', 'direction_a'), ('b', 'direction_b')):
            endpoint = join[key]
            direction = np.asarray(join[direction_key], dtype=float)
            vectors = [direction]
            inside = [a.GetIdx() for a in mol.GetAtomWithIdx(endpoint).GetNeighbors()
                      if owners[a.GetIdx()] == owners[endpoint]]
            if len(inside) == 1:
                axis = spaced[inside[0]]-spaced[endpoint]
                axis /= np.linalg.norm(axis)
                reflected = 2*np.dot(direction, axis)*axis-direction
                if np.linalg.norm(reflected-direction) > 1e-8:
                    vectors.append(reflected)
            candidates.append(vectors)
        axis = spaced[join['b']]-spaced[join['a']]
        axis /= np.linalg.norm(axis)
        ranked = []
        for ia, da in enumerate(candidates[0]):
            for ib, db in enumerate(candidates[1]):
                score = 20*np.dot(da+db, da+db)+np.dot(da-axis, da-axis)+np.dot(db+axis, db+axis)
                ranked.append((round(float(score), 8), ia, ib, {**join, 'direction_a': da, 'direction_b': db}))
        choices.append([item[3] for item in sorted(ranked, key=lambda row: row[:3])])
    before = layout_diagnostics(mol, spaced, origins)
    target = before['median_bond_length']
    baseline_overlaps = {tuple(pair) for pair in before['nonbonded_overlap_pairs']}
    baseline_crossings = {tuple(pair) for pair in before['nonincident_bond_crossings']}
    def metrics(xy):
        diagnostic = layout_diagnostics(mol, xy, origins)
        angles = [row['deviation_degrees'] for row in diagnostic['junction_angles']]
        ratios = diagnostic['junction_bond_ratios']
        length_error = max((abs(value-1) for value in ratios), default=0.)
        return diagnostic, max(angles, default=0.), float(np.mean(angles)) if angles else 0., length_error
    best = spaced
    best_values = metrics(best)
    best_score = best_values[1]+best_values[2]+180*best_values[3]
    attempts = []
    selected = None
    for pattern in (0, 1, 2, 3):
        selected_joins = [options[min((pattern if pattern < 2 else (index+pattern) % 2), len(options)-1)]
                          for index, options in enumerate(choices)]
        fitted, fit = refine_rigid_poses(mol, spaced, origins, selected_joins, robust=True)
        values = metrics(fitted)
        new_overlaps = {tuple(pair) for pair in values[0]['nonbonded_overlap_pairs']}-baseline_overlaps
        new_crossings = {tuple(pair) for pair in values[0]['nonincident_bond_crossings']}-baseline_crossings
        score = values[1]+values[2]+180*values[3]
        accept = not new_overlaps and not new_crossings and values[3] <= .12 and score < best_score-1e-6
        attempts.append({'pattern': pattern, 'accepted': accept, 'max_angle_deviation_degrees': values[1],
                         'mean_angle_deviation_degrees': values[2], 'max_join_length_fraction_error': values[3],
                         'solver_applied': fit['applied'], 'solver_reason': fit.get('reason')})
        if accept:
            best, best_values, best_score, selected = fitted, values, score, fit
    report = {'method': 'bond_angle_pose_v2', 'applied': bool(spacing.get('applied') or selected),
              'interpretation': 'Rigid-block 2D depiction with unspecified connection-side choices; no energy or E/Z prediction.',
              'target_bond_length': target, 'angle_fit_applied': selected is not None,
              'initial_max_junction_angle_deviation_degrees': max((r['deviation_degrees'] for r in before['junction_angles']), default=0.),
              'final_max_junction_angle_deviation_degrees': best_values[1],
              'final_mean_junction_angle_deviation_degrees': best_values[2],
              'max_join_length_fraction_error': best_values[3], 'candidates': attempts,
              'quality': 'regular' if best_values[1] <= 15 and best_values[3] <= .12
                         and not best_values[0]['nonbonded_overlap_pairs']
                         and not best_values[0]['nonincident_bond_crossings'] else 'needs_review',
              'rigid_intrablock_geometry_preserved': True}
    if selected:
        report['block_transforms'] = selected['block_transforms']
    return best, report


def clean_layout(recipe):
    """Return an explicit clean-layout design, with its exact fitted coordinates."""
    from .framework_drawing import build_patch
    if recipe.get('explicit_graph'):
        raise ValueError('This source-registered graph keeps its supplied atom coordinates; rigid cleanup requires an assembled monomer recipe')
    cleaned = deepcopy(recipe)
    cleaned.pop('atom_coordinates', None)
    cleaned['layout_refinement'] = 'bond_angle_pose_v2'
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
