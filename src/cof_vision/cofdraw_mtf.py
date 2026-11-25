"""One chemical layer of the voided square grid, after COF-432 Fig. 2b–c."""
from collections import Counter
from functools import lru_cache
from itertools import permutations, product

import numpy as np
from rdkit import Chem

from .cofdraw_mixed import ETTA, _block, _turn, join_fragments
from .cofdraw_net_graphs import quotient

_TFB = Chem.MolToSmiles(Chem.MolFromSmiles('O=Cc1cc(C=O)cc(C=O)c1'))
_TETRAARYLETHENE = Chem.MolFromSmarts('[c]-[C;!R](=[C;!R]-[c])-[c]')


@lru_cache(maxsize=1)
def _motif():
    """Register the seven RCSR sites to the parent square grid in the paper's view."""
    adjacency, fractional, _ = quotient('mtf')
    # The first four-connected site is surrounded by four small faces. One
    # diagonal parent-grid step is (1,1); translations span four units in x/y.
    points = np.rint(4 * (fractional - fractional[0])).astype(int) + [2, 2]
    pairs, shifts, slots = [], [], []
    degree = Counter()
    for i, row in enumerate(adjacency):
        if len(row) != 3:
            continue
        for j, shift in row:
            pairs.append((i, j))
            shifts.append(shift)
            slots.append((degree[i], degree[j]))
            degree.update((i, j))
    return points, pairs, np.array(shifts), slots, [len(row) for row in adjacency]


def _faces():
    central = [(0, 2), (1, 1), (2, 0), (1, -1), (0, -2), (-1, -1), (-2, 0), (-1, 1)]
    faces = [central]
    for sx, sy in ((1, 1), (-1, 1), (-1, -1), (1, -1)):
        # Four little squares share each fully four-connected hub. Begin
        # with the two faces that also share edges with the central aperture.
        hub = (2*sx, 2*sy)
        for dx, dy in ((-sx, 0), (0, -sy), (sx, 0), (0, sy)):
            x, y = hub
            faces.append([hub, (x+dx-dy, y+dy+dx), (x+2*dx, y+2*dy), (x+dx+dy, y+dy-dx)])
    return faces


def graph(components, count):
    faces = _faces()
    if count == 0:
        faces = faces[1:2]
    elif count in (1, 2, 3, 17):
        faces = faces[:count]
    else:
        raise ValueError('Choose a supported voided square grid extent.')
    points = sorted({p for face in faces for p in face})
    index = {p: i for i, p in enumerate(points)}
    faces = [[index[p] for p in face] for face in faces]
    pairs = sorted({tuple(sorted((a, b))) for face in faces for a, b in zip(face, face[1:] + face[:1])})
    by_degree = {c['connectivity']: k for k, c in components.items()}
    kinds = [by_degree[3 if x % 2 else 4] for x, _ in points]
    return faces, pairs, np.array(points), kinds


def _poses(geometry, kinds, grid, edges, shifts):
    """Align the two tetratopic environments and direct one trigonal arm towards its hub."""
    incident = [[] for _ in grid]
    for (a, b), shift in zip(edges, shifts):
        ray = grid[b] + 4 * shift - grid[a]
        incident[a].append((b, ray))
        incident[b].append((a, -ray))
    orientations, port_numbers = [], []
    for i, row in enumerate(incident):
        mol, ports, xy, _ = geometry[kinds[i]]
        rays = xy[[s for _, s in ports]]
        phases = np.arctan2(rays[:, 1], rays[:, 0])
        directions = np.array([np.arctan2(v[1], v[0]) for _, v in row])
        if len(row) == 3:
            hub = next(j for j, (other, _) in enumerate(row) if np.all(grid[other] % 4 == (2, 2)))
            seeds = directions[hub] - phases
        else:
            # ETTA's central alkene is horizontal on the x-extending arms of
            # the vacancy and vertical on the y-extending arms (Fig. 2c).
            alkene = mol.GetSubstructMatch(_TETRAARYLETHENE)
            if alkene:
                axis = xy[alkene[2]] - xy[alkene[1]]
                baseline = np.arctan2(axis[1], axis[0])
            else:
                baseline = phases[0] - np.pi/4
            target = 0 if grid[i, 1] % 4 == 0 else np.pi/2
            seeds = [target - baseline, target - baseline + np.pi]
        candidates = []
        for seed in seeds:
            theta = np.arctan2(np.sin(seed), np.cos(seed))
            for slots in permutations(range(len(row))):
                error = np.sum(1 - np.cos(phases[list(slots)] + theta - directions))
                candidates.append(((round(error, 8), round(abs(theta), 8), slots), (theta, slots)))
        theta, slots = min(candidates)[1]
        orientations.append(theta)
        port_numbers.append(slots)
    return np.array(orientations), port_numbers


def _project_joins(geometry, kinds, reference, edges, shifts, incident, source_pair):
    initial, numbered = _poses(geometry, kinds, reference, edges, shifts)
    ports = np.array([[numbered[i][s] for i, s in zip(edge, pair)] for edge, pair in zip(edges, incident)])
    # The structural model has staggered layers; this is only its in-plane
    # drawing. Keep ring cores rigid and form each trans imine between its
    # two actual scaffold anchors. C and N are assigned by the selected feed.
    owners = np.array(edges)
    a, b = owners.T
    anchors = np.empty((len(edges), 2, 2))
    rays = np.empty_like(anchors)
    for e, (edge, pair) in enumerate(zip(edges, ports)):
        for side, (i, slot) in enumerate(zip(edge, pair)):
            mol, handles, xy, _ = geometry[kinds[i]]
            dummy, site = handles[slot]
            anchor = next(n.GetIdx() for n in mol.GetAtomWithIdx(site).GetNeighbors() if n.GetIdx() != dummy)
            anchors[e, side] = xy[anchor]
            rays[e, side] = (xy[site] - xy[anchor]) / np.linalg.norm(xy[site] - xy[anchor])
    arms = reference[b] + 4*shifts - reference[a]
    # Position the three- and four-connected blocks using complete incident
    # stars. Only the three-connected sites have an internal positional
    # parameter: the vacancy corners lie on the axes of a square cell.
    tri = sorted(set(a))
    nominal_span = np.hypot(3.75, 1.5*np.sqrt(3)/2)  # three standard 1.5-unit bonds in a zigzag
    offset = _turn(anchors, initial[owners])
    radial_span = np.sum(arms * (offset[:, 1] - offset[:, 0]), axis=1) / np.sqrt(2)
    scale = np.median(nominal_span - radial_span) / np.sqrt(2)
    start = np.r_[scale, np.zeros(2*len(tri)), initial]

    def geometry_at(parameters):
        centers = parameters[0] * reference.astype(float)
        centers[tri] += parameters[1:1+2*len(tri)].reshape(-1, 2)
        theta = parameters[1+2*len(tri):]
        positioned = _turn(anchors, theta[owners]) + centers[owners]
        left, right = positioned[:, 0], positioned[:, 1]
        right = right + shifts * (4*parameters[0])
        chord = right - left
        lengths = np.linalg.norm(chord, axis=1)
        axis = chord / lengths[:, None]
        original = _turn(rays, theta[owners])
        return theta, centers, left, right, lengths, axis, original
    # Fig. 2c is the mirror of the ab-coordinate convention of SI Table S1.
    # The table gives three distinct projected bends at the trigonal (C) and
    # tetratopic (N) end. Only transfer these angles to the reported feed.
    _, _, _, _, _, ax, ref_dirs = geometry_at(start)
    normal = np.column_stack((-ax[:, 1], ax[:, 0]))
    turns = np.full((len(edges), 2), 60.)
    for i in tri:
        ix = np.flatnonzero(a == i)
        if source_pair:
            hub = next(e for e in ix if np.all(reference[b[e]] % 4 == (2, 2)))
            for e in ix:
                left_turn = np.linalg.det([arms[hub], arms[e]]) > 0
                # Exterior turns at CH and N, rounded from C1/N1, C6/N3, C19/N2.
                if e == hub:
                    turns[e] = (57., 59.)
                elif left_turn:
                    turns[e] = (-47., -37.)
                else:
                    turns[e] = (52., 38.)
        else:
            # Pick the mtf 2:1 turn pattern closest to the selected monomers'
            # scaffold directions. No published coordinate is implied here.
            def cost(choice):
                direction = (3.75*ax[ix] - 1.5*np.sqrt(3)/2*np.array(choice)[:, None]*normal[ix])/nominal_span
                return np.sum((direction-ref_dirs[ix, 0])**2 + (direction+ref_dirs[ix, 1])**2)
            signs = min((v for v in product((-1, 1), repeat=3) if abs(sum(v)) == 1), key=cost)
            turns[ix] *= np.array(signs)[:, None]
    # Relative vectors for scaffold–CH=N–scaffold in an orientable local
    # frame. The two flanking bonds may have different projected angles.
    t1, t2 = np.radians(turns).T
    first_arm = np.array([1., 0.])
    bond = np.column_stack((np.cos(t1), np.sin(t1)))
    second_arm = np.column_stack((np.cos(t1-t2), np.sin(t1-t2)))
    # Use standard single/double skeletal lengths; only the small common
    # scale needed to close each three-bond span will vary by edge.
    mol, handles, _, _ = geometry[kinds[a[0]]]
    first_site = mol.GetAtomWithIdx(handles[ports[0, 0]][1])
    la, lb = (1.47, 1.42) if first_site.GetAtomicNum() == 6 else (1.42, 1.47)
    path = la*first_arm + 1.28*bond + lb*second_arm
    span = np.linalg.norm(path, axis=1)
    tilt = -np.arctan2(path[:, 1], path[:, 0])

    def arm_directions(axis):
        return _turn(axis, tilt), _turn(axis, tilt+t1-t2)

    def residual(parameters):
        theta, _, _, _, lengths, ax, original = geometry_at(parameters)
        first, second = arm_directions(ax)
        return np.r_[2.5*(lengths-span), (first-original[:, 0]).ravel(), (second+original[:, 1]).ravel(),
                     .65*parameters[1:1+2*len(tri)], theta-initial]

    # Small damped Gauss–Newton fit on seven poses; no crystal relaxation.
    parameters = start
    for _ in range(18):
        base = residual(parameters)
        jac = np.stack([(residual(parameters+d)-base)/1e-5 for d in 1e-5*np.eye(len(start))], axis=1)
        step = np.linalg.lstsq(jac, -base, rcond=None)[0]
        for factor in (1., .5, .25, .1):
            trial = parameters + factor*step
            if np.linalg.norm(residual(trial)) < np.linalg.norm(base):
                parameters = trial
                break
        else:
            break
        if np.max(abs(factor*step)) < 1e-5:
            break
    theta, centers, left, right, lengths, ax, _ = geometry_at(parameters)
    c_arm, n_arm = arm_directions(ax)
    scale = (lengths / span)[:, None]
    first = left + la*c_arm*scale
    second = right - lb*n_arm*scale
    tile = 4*parameters[0]*np.eye(2)
    local = [_turn(geometry[k][2], angle).copy() for k, angle in zip(kinds, theta)]
    for e, (i, j) in enumerate(edges):
        displacement = shifts[e]@tile
        for owner, s, site, tip, shift in ((i, 0, first[e], second[e], 0), (j, 1, second[e], first[e], displacement)):
            d, k = geometry[kinds[owner]][1][ports[e,s]]
            local[owner][k] = site - centers[owner] - shift
            local[owner][d] = tip - centers[owner] - shift
    return centers, tile, local, ports


@lru_cache(maxsize=16)
def _fitted_motif(fragments, source_pair):
    reference, edges, shifts, incident, degree = _motif()
    by_degree = {n: k for k, _, n in fragments}
    geometry = {k: _block(dict(fragment_smiles=s, linkage='imine')) for k, s, _ in fragments}
    poses = _project_joins(geometry, [by_degree[n] for n in degree], reference, edges, shifts, incident, source_pair)
    return geometry, poses


def patch(components, count):
    from .cofdraw_spatial import check_projection
    faces, pairs, points, kinds = graph(components, count)
    reference, edges, shifts, _, degree = _motif()
    monomers = {c['connectivity']: c['monomer_smiles'] for c in components.values()}
    source_pair = monomers == {3: _TFB, 4: ETTA}
    fragments = tuple((k, c['fragment_smiles'], c['connectivity']) for k, c in sorted(components.items()))
    geometry, (centers, tile, local, ports) = _fitted_motif(fragments, source_pair)
    # Place every crop from the same seven molecular poses and retain the
    # port assigned to every translated edge, including at the patch boundary.
    index = {tuple(p % 4): i for i, p in enumerate(reference)}
    site = np.array([index[tuple(p % 4)] for p in points])
    offset = (points - reference[site]) // 4
    translations = centers[site] + offset @ tile
    lookup = {(a, b, tuple(t)): port for (a, b), t, port in zip(edges, shifts, ports)}
    joined_ports = []
    for a, b in pairs:
        key = (site[a], site[b], tuple(offset[b] - offset[a]))
        if degree[site[a]] == 3:
            joined_ports.append(lookup[key])
        else:
            joined_ports.append(lookup[site[b], site[a], tuple(offset[a] - offset[b])][::-1])
    mol, origins = join_fragments(components, kinds, points, pairs, None, 2, geometry=geometry,
                                  placement=(translations, np.zeros(len(site)), joined_ports))
    conf = mol.GetConformer()
    blocks = {f'{kind}{i+1}': i for i, kind in enumerate(kinds)}
    for i, origin in enumerate(origins):
        block = blocks[origin['instance']]
        x, y = local[site[block]][origin['fragment_atom_index']] + translations[block]
        conf.SetAtomPosition(i, (float(x), float(y), 0.))
    check_projection(mol)
    rims = [[f'{kinds[i]}{i+1}' for i in face] for face in faces]
    labels = ['Large square (8 nodes)' if len(f) == 8 else 'Small square (4 nodes)' for f in faces]
    basis = 'One mtf layer in the orientation of COF-432 Figure 2 (DOI 10.1021/jacs.9b13094). '
    basis += ('Junction angles follow the projected model in Table S1.' if source_pair else
              'The selected monomers use schematic imine angles in a proposed structure.')
    return mol, dict(pores=len(faces), pore_faces=rims, face_labels=labels, pore_counts=dict(Counter(labels)),
                     atom_origin=origins, order_assumption=basis)
