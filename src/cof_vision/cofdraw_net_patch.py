"""Construct bounded molecular drawings on diamond, platinum sulfide and kagome net fragments."""
from itertools import combinations, permutations

import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem

from .cofdraw_depiction import compute_2d_coords, port_symmetry_ranks
from .framework_drawing import _ports


def _skeleton(diamond):
    if diamond:
        # The ten vertices of one adamantane cage in the diamond lattice.
        tetra = np.array([[1, 1, 1], [1, -1, -1], [-1, 1, -1], [-1, -1, 1]]) / np.sqrt(3)
        centers = list(tetra)
        edges = []
        for a, b in combinations(range(4), 2):
            middle = len(centers)
            centers.append(tetra[a] + tetra[b])
            edges.extend([(a, middle), (middle, b)])
        return np.array(centers), edges
    # One kagome star: every hexagon edge is shared with an exterior triangle.
    angles = np.arange(6) * np.pi / 3
    ring = np.column_stack((np.cos(angles), np.sin(angles), np.zeros(6)))
    centers = list(ring)
    ring_edges = [(i, (i + 1) % 6) for i in range(6)]
    edges = list(ring_edges)
    for i, j in ring_edges:
        tip = len(centers)
        centers.append(ring[i] + ring[j])
        edges.extend([(i, tip), (tip, j)])
    return np.array(centers), edges


def _coordinates(component, spatial):
    mol = Chem.MolFromSmiles(component["fragment_smiles"])
    ports = _ports(mol)
    if spatial:
        # Neutral carbon caps preserve port atom indices during embedding.
        capped = Chem.RWMol(mol)
        for dummy, _ in ports.values():
            capped.ReplaceAtom(dummy, Chem.Atom(6))
        Chem.SanitizeMol(capped)
        capped = Chem.AddHs(capped)
        params = AllChem.ETKDGv3()
        params.randomSeed = 300
        if AllChem.EmbedMolecule(capped, params) != 0:
            raise ValueError("Could not embed the selected tetrahedral building block")
        xyz = capped.GetConformer().GetPositions()[:mol.GetNumAtoms()].copy()
    else:
        compute_2d_coords(mol)
        xyz = mol.GetConformer().GetPositions()
    ends = [neighbor for _, neighbor in ports.values()]
    xyz -= xyz[ends].mean(axis=0)
    return xyz, ports


def _unit(vectors):
    return vectors / np.linalg.norm(vectors, axis=-1, keepdims=True)


def _orient_node(xyz, ports, directions):
    # Fit actual molecular arms to incident edges without reflecting stereocenters.
    best = None
    for keys in permutations(ports, len(directions)):
        arms = _unit(xyz[[ports[key][1] for key in keys]])
        u, _, vt = np.linalg.svd(arms.T @ directions)
        rotation = u @ np.diag([1, 1, np.linalg.det(u @ vt)]) @ vt
        error = np.linalg.norm(arms @ rotation - directions)
        if best is None or error < best[0] - 1e-9:
            best = error, keys, rotation
    return best[1], xyz @ best[2]


def _view(spatial):
    # Oblique, orthographic view. Store XYZ for other projections of the same graph.
    azimuth, elevation = np.radians([24, 27])
    ca, sa = np.cos(azimuth), np.sin(azimuth)
    ce, se = np.cos(elevation), np.sin(elevation)
    return np.array([[ca, -sa * se], [sa, ca * se], [0, ce]]) if spatial else np.eye(3)[:, :2]


def _add_block(items, atoms, view, spatial, kind, index, xyz, ports):
    name = f"{kind}{index + 1}"
    keys = sorted(ports)
    items.append(dict(id=name, component=kind, incoming_port=keys[0], outgoing_port=keys[1],
                      layout_position=(xyz.mean(axis=0) @ view).round(5).tolist()))
    for atom, pos in enumerate(xyz):
        row = dict(instance=name, fragment_atom_index=atom, xy=(pos @ view).round(5).tolist())
        if spatial:
            row["xyz"] = pos.round(5).tolist()
            row["depth"] = round(float(pos @ np.cross(view[:, 0], view[:, 1])), 5)
        atoms[name, atom] = row
    return name


def _cn_junctions(components, items, connections, atoms, view, spatial, *, reduced=False):
    """Place the C–N or C=N zigzag between the two fixed scaffold anchors."""
    ends = {}
    for kind, component in components.items():
        mol = Chem.MolFromSmiles(component["fragment_smiles"])
        for port, (dummy, joined) in _ports(mol).items():
            atom = mol.GetAtomWithIdx(joined)
            neighbors = [a.GetIdx() for a in atom.GetNeighbors() if a.GetIdx() != dummy]
            if len(neighbors) == 1:
                ends[kind, port] = (atom.GetSymbol(), joined, neighbors[0])
    kinds = {i["id"]: i["component"] for i in items}
    normal = np.cross(view[:, 0], view[:, 1])
    # Parallel scaffold bonds and a 60° turn give a conventional skeletal
    # zigzag. Scale nominal length ratios to the available anchor span.
    scaffold_c, junction, scaffold_n = (1.51, 1.47, 1.42) if reduced else (1.47, 1.28, 1.42)
    axial = scaffold_c + scaffold_n + junction * 0.5
    lateral = junction * np.sqrt(3.) / 2
    span_squared = axial**2 + lateral**2
    for edge in connections:
        groups = {}
        for end in (edge["a"], edge["b"]):
            name = end["instance"]
            info = ends.get((kinds[name], end["port"]))
            if info:
                symbol, joined, anchor = info
                pos = atoms[name, anchor]
                groups[symbol] = (atoms[name, joined], np.array(pos["xyz"] if spatial else [*pos["xy"], 0.]))
        if set(groups) != {"C", "N"}:
            continue
        carbon, left = groups["C"]
        nitrogen, right = groups["N"]
        delta = right - left
        distance = np.linalg.norm(delta)
        if distance < 1e-8:
            raise ValueError("The two C–N attachment sites coincide in this net drawing.")
        along = delta / distance
        # The azimuth around the attachment axis is free in this unrelaxed
        # schematic. Choose a transverse direction visible in the initial view.
        transverse = np.cross(normal, along)
        if np.linalg.norm(transverse) < 1e-8:
            transverse = view[:, 0]
        arm = (axial * along - lateral * _unit(transverse)) * distance / span_squared
        for row, xyz in ((carbon, left + scaffold_c * arm), (nitrogen, right - scaffold_n * arm)):
            row["xy"] = (xyz @ view).round(5).tolist()
            if spatial:
                row.update(xyz=xyz.round(5).tolist(), depth=round(float(xyz @ normal), 5))


def _complete_kagome_arms(component, kind, items, atoms, connections, linker_ports):
    """Continue a symmetry-equivalent opposite arm at the edge of the finite star."""
    mol = Chem.MolFromSmiles(component["fragment_smiles"])
    ports = _ports(mol)
    ranks = port_symmetry_ranks(mol)
    anchors = {}
    for port, (dummy, terminal) in ports.items():
        neighbors = [a.GetIdx() for a in mol.GetAtomWithIdx(terminal).GetNeighbors() if a.GetIdx() != dummy]
        if len(neighbors) != 1:
            return
        anchors[port] = neighbors[0]
    # net_patch stores each edge with its node first and its spacer second.
    partners = {}
    for edge in connections:
        node, spacer = edge["a"], edge["b"]
        other = linker_ports[spacer["port"]][1]
        partners[node["instance"], node["port"]] = atoms[spacer["instance"], other]["xy"]
    for item in items:
        if item["component"] != kind:
            continue
        name = item["id"]
        open_ports = [p for p, (dummy, _) in ports.items() if (name, dummy) in atoms]
        if not open_ports:
            continue
        xy = {p: np.array(atoms[name, anchor]["xy"]) for p, anchor in anchors.items()}
        center = np.mean(list(xy.values()), axis=0)
        for p in open_ports:
            dummy, terminal = ports[p]
            for q, (other_dummy, other_terminal) in ports.items():
                if q in open_ports or ranks[other_dummy] != ranks[dummy]:
                    continue
                if not np.allclose(xy[p] + xy[q], 2 * center, rtol=0, atol=0.02):
                    continue
                # Rotate the occupied arm by half a turn about the node, including
                # its first scaffold bond. The molecular backbone stays fixed.
                anchor_sum = xy[p] + xy[q]
                atoms[name, terminal]["xy"] = (anchor_sum - atoms[name, other_terminal]["xy"]).round(5).tolist()
                atoms[name, dummy]["xy"] = (anchor_sum - partners[name, q]).round(5).tolist()
                break


def net_patch(components, node, diamond, bond_order):
    linker = "B" if node == "A" else "A"
    centers, edges = _skeleton(diamond)
    node_xyz, node_ports = _coordinates(components[node], diamond)
    link_xyz, link_ports = _coordinates(components[linker], False)
    p, q = sorted(link_ports)
    left, right = (link_ports[k][1] for k in (p, q))
    axis = _unit(link_xyz[right] - link_xyz[left])
    across = np.array([-axis[1], axis[0], 0])
    radius = np.mean([np.linalg.norm(node_xyz[n]) for _, n in node_ports.values()])
    linkage = components[node]["linkage"]
    imine = linkage in ("imine", "alkyl_imine")
    reduced = not diamond and linkage in ("secondary_amine", "reduced_alkyl_imine")
    gap = 0.8 if imine or reduced else 1.4
    centers *= 2 * radius + np.linalg.norm(link_xyz[right] - link_xyz[left]) + 2 * gap
    view = _view(diamond)
    linker_plane_normal = np.cross(view[:, 0], view[:, 1])
    if diamond:
        # Orthographic sightline fitted to the ten COF-300 junctions in Ma et al.,
        # Science 361, Fig. 2A (doi:10.1126/science.aat7679). Roll the view so the
        # long [111] cage edges stand upright. Only the 2D projection changes.
        sightline = _unit(np.array([-.4312, .8936, -.1249]))
        horizontal = _unit(np.cross([1., 1., 1.], sightline))
        view = np.column_stack((horizontal, np.cross(sightline, horizontal)))
    items, atoms, connections, ends = [], {}, [], {}

    for i, center in enumerate(centers):
        neighbors = [b if a == i else a for a, b in edges if i in (a, b)]
        directions = _unit(centers[neighbors] - center)
        keys, oriented = _orient_node(node_xyz, node_ports, directions)
        coords = oriented + center
        name = _add_block(items, atoms, view, diamond, node, i, coords, node_ports)
        for j, key in zip(neighbors, keys):
            dummy, atom = node_ports[key]
            ends[i, j] = (dict(instance=name, port=key), coords[atom], dummy)

    for i, (a, b) in enumerate(edges):
        end_a, point_a, dummy_a = ends[a, b]
        end_b, point_b, dummy_b = ends[b, a]
        along = _unit(point_b - point_a)
        start, finish = point_a + gap * along, point_b - gap * along
        # Draw the linker in a plane containing its two real attachment atoms.
        # Its transverse basis stays in the picture plane for legible rings.
        sideways = _unit(np.cross(linker_plane_normal, along))
        offset = link_xyz - link_xyz[left]
        fraction = (offset @ axis) / np.linalg.norm(link_xyz[right] - link_xyz[left])
        xyz = start + fraction[:, None] * (finish - start) + np.outer(offset @ across, sideways)
        name = _add_block(items, atoms, view, diamond, linker, i, xyz, link_ports)
        for end, dummy, key in ((end_a, dummy_a, p), (end_b, dummy_b, q)):
            connections.append(dict(a=end, b=dict(instance=name, port=key), bond_order=bond_order))
            del atoms[end["instance"], dummy], atoms[name, link_ports[key][0]]
    if imine or reduced:
        _cn_junctions(components, items, connections, atoms, view, diamond, reduced=reduced)
        if not diamond:
            _complete_kagome_arms(components[node], node, items, atoms, connections, link_ports)
    return items, connections, list(atoms.values())


def _pts_neighbors(site):
    # PtS basis (AFLOW AB_tP4_131_c_e-001), in units of (a/2, a/2, c/4).
    # Pt sites have one horizontal axis; the S site's axis alternates above/below.
    x, y, z = site
    neighbors = []
    for dz in (-1, 1):
        planar_level = z if z % 2 == 0 else z + dz
        for side in (-1, 1):
            dx, dy = (0, side) if planar_level % 4 == 0 else (side, 0)
            neighbors.append((x + dx, y + dy, z + dz))
    return neighbors


def pts_patch(components, tetra, bond_order):
    """Two shells about one Pt–S edge, with exposed arms kept as continuation sites."""
    planar = "A" if tetra == "B" else "B"
    sites = {(0, 1, 0), (0, 0, 1)}
    for _ in range(2):
        sites |= {n for site in sites for n in _pts_neighbors(site)}
    sites = sorted(sites)
    index = {site: i for i, site in enumerate(sites)}
    fragments = {kind: _coordinates(component, spatial=kind == tetra)
                 for kind, component in components.items()}
    # Idealized c/a = sqrt(2) yields tetrahedral S and coplanar rectangular Pt.
    metric = np.array([1., 1., 1. / np.sqrt(2.)])
    imine = components[tetra]["linkage"] in ("imine", "alkyl_imine")
    gap = 0.8 if imine else 1.4
    blocks, edges, target_deltas = [], [], []
    for site in sites:
        kind = tetra if site[2] % 2 else planar
        xyz, ports = fragments[kind]
        neighbors = _pts_neighbors(site)
        directions = _unit((np.array(neighbors) - site) * metric)
        keys, oriented = _orient_node(xyz, ports, directions)
        blocks.append((kind, oriented, dict(zip(neighbors, keys))))

    for i, site in enumerate(sites):
        ka, xa, pa = blocks[i]
        for neighbor, key in pa.items():
            j = index.get(neighbor)
            if j is None or j <= i:
                continue
            kb, xb, pb = blocks[j]
            other_key = pb[site]
            end_a = fragments[ka][1][key][1]
            end_b = fragments[kb][1][other_key][1]
            direction = _unit((np.array(neighbor) - site) * metric)
            edges.append((i, key, j, other_key))
            target_deltas.append(xa[end_a] - xb[end_b] + gap * direction)

    # Fit translations jointly around all closed rings; keep each monomer rigid.
    incidence = np.zeros((len(edges), len(sites)))
    for k, (i, _, j, _) in enumerate(edges):
        incidence[k, i], incidence[k, j] = -1, 1
    centers = np.linalg.lstsq(incidence, target_deltas, rcond=None)[0]
    items, atoms, connections, names = [], {}, [], []
    view = _view(True)
    for i, (kind, xyz, _) in enumerate(blocks):
        names.append(_add_block(items, atoms, view, True, kind, i, xyz + centers[i], fragments[kind][1]))
    for i, a, j, b in edges:
        connections.append(dict(a=dict(instance=names[i], port=a), b=dict(instance=names[j], port=b),
                                bond_order=bond_order))
        del atoms[names[i], fragments[blocks[i][0]][1][a][0]]
        del atoms[names[j], fragments[blocks[j][0]][1][b][0]]
    if imine:
        _cn_junctions(components, items, connections, atoms, view, True)
    return items, connections, list(atoms.values())
