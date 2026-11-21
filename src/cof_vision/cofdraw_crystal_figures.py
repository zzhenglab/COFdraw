"""Finite illustrations unwrapped along the validated display graph's image shifts."""
from collections import deque
from functools import lru_cache

import numpy as np
from rdkit import Chem


def readable_direction(model, face=None):
    """Prefer views that keep six-membered rings open and distinct atoms apart."""
    atoms = model["atoms"]
    xyz = np.array([[a[k] for k in ("x","y","z")] for a in atoms])
    mol = Chem.RWMol()
    for _ in atoms:
        mol.AddAtom(Chem.Atom(6))
    for a,b in model["bonds"]:
        mol.AddBond(a,b,Chem.BondType.SINGLE)
    normals = []
    for ring in Chem.GetSymmSSSR(mol):
        if len(ring)==6:
            points = xyz[list(ring)]
            normals.append(np.linalg.svd(points-points.mean(axis=0),full_matrices=False)[2][-1])
    if not normals:
        return
    indices = np.arange(192)+.5
    z = 1-2*indices/192
    theta = indices*2.399963229728653
    r = np.sqrt(1-z*z)
    directions = np.column_stack([r*np.cos(theta),r*np.sin(theta),z])
    exposure = np.abs(directions @ np.asarray(normals).T)
    scores = .85*np.quantile(exposure,.15,axis=1) + .15*exposure.mean(axis=1)
    if face is not None:
        alignment = np.abs(directions @ face)
        scores += alignment
        scores[alignment < .80] = -np.inf
    # Discourage projected overlaps within a local neighborhood. Sample pairs
    # deterministically so large linkers cannot create a large work array.
    a,b = np.triu_indices(len(xyz),1)
    delta = xyz[a]-xyz[b]
    local = delta[np.linalg.norm(delta,axis=1)<8][::max(1,len(delta)//8000)]
    if len(local):
        perpendicular = (local*local).sum(axis=1)[:,None]-(local @ directions.T)**2
        scores -= 12*(perpendicular < .7**2).mean(axis=0)
    model["direction"] = directions[np.argmax(scores)].tolist()


def crystal_details(structure, edges):
    """Identify a four-aryl junction and trace linkers to adjacent junctions.

    A contraction is offered only when each arm of the selected motif terminates
    at exactly one other junction. Branched linkers and ambiguous graphs retain
    the atomistic representation.
    """
    sites = structure["atoms"]
    graph = [[] for _ in sites]
    for edge in edges:
        a, b = edge["a"], edge["b"]
        if "H" in (sites[a]["element"], sites[b]["element"]):
            continue
        t = tuple(edge["translation"])
        graph[a].append((b, t))
        graph[b].append((a, tuple(-v for v in t)))
    centers = {i for i, neighbors in enumerate(graph)
               if sites[i]["element"] in {"C", "Si"} and len(neighbors) == 4
               and all(sites[j]["element"] == "C" and len(graph[j]) == 3 for j, _ in neighbors)}
    if not centers:
        return {}
    cell = np.array(structure["cell"])
    def neighbors(site):
        i, x, y, z = site
        return [(j, x+a, y+b, z+c) for j, (a,b,c) in graph[i]]

    def model(nodes, connections):
        nodes = sorted(nodes)
        indices = {site: i for i, site in enumerate(nodes)}
        atoms = []
        for i, x, y, z in nodes:
            px, py, pz = (np.array(sites[i]["fractional"]) + [x,y,z]) @ cell
            atoms.append(dict(element=sites[i]["element"], x=px, y=py, z=pz))
        return dict(atoms=atoms, bonds=[[indices[a],indices[b]] for a,b in sorted(connections)])

    def molecular(nodes):
        bonds = {(a,b) for a in nodes for b in neighbors(a) if b in nodes and a < b}
        return model(nodes, bonds)

    root = (min(centers),0,0,0)
    nodes = frontier = {root}
    for _ in range(4):
        frontier = {other for a in frontier for other in neighbors(a)} - nodes
        nodes = nodes | frontier
    details = {"junction": molecular(nodes)}
    readable_direction(details["junction"])

    @lru_cache(maxsize=256)
    def arms(site):
        result = {}
        for first in neighbors(site):
            visited, ends, pending = {first}, set(), [first]
            while pending and len(visited) <= 160:
                for other in neighbors(pending.pop()):
                    if other[0] in centers:
                        ends.add(other)
                    elif other not in visited:
                        visited.add(other)
                        pending.append(other)
            ends.discard(site)
            if pending or len(ends) != 1:
                return {}
            result[ends.pop()] = visited
        return result if len(result) == 4 else {}

    adjacent = arms(root)
    if not adjacent or any(not arms(junction) for junction in adjacent):
        return details
    bonds = {tuple(sorted((a,b))) for a in [root,*adjacent] for b in arms(a)}
    details["net"] = model({site for edge in bonds for site in edge}, bonds)
    shell = {root,*adjacent}
    for group in adjacent.values():
        shell.update(group)
    overlay = molecular(shell)
    indices = {site:i for i,site in enumerate(sorted(shell))}
    overlay["guides"] = [[indices[root],indices[n]] for n in adjacent]
    details["overlay"] = overlay
    # A shortest closed junction circuit through the central junction. Its
    # atomistic drawing includes every atom of each intervening heavy-atom arm.
    first = min(adjacent)
    paths = deque([[first]])
    seen = {root,first}
    while paths and len(seen) <= 128:
        path = paths.popleft()
        if len(path) > 1 and path[-1] in adjacent:
            loop = [root,*path,root]
            atoms = set(loop)
            for a,b in zip(loop,loop[1:]):
                atoms.update(arms(a)[b])
            # Retain the other aryl arms at each junction as finite branches.
            for site in loop[:-1]:
                branch = {n for n in neighbors(site) if n not in atoms}
                for _ in range(4):
                    atoms.update(branch)
                    branch = {n for a in branch for n in neighbors(a)} - atoms
            details["ring"] = molecular(atoms)
            points = np.array([(np.asarray(sites[i]["fractional"]) + [x,y,z]) @ cell for i,x,y,z in path+[root]])
            _, _, vectors = np.linalg.svd(points-points.mean(axis=0), full_matrices=False)
            readable_direction(details["ring"],vectors[-1])
            break
        if len(path) < 8:
            for end in sorted(arms(path[-1])):
                if end not in seen:
                    seen.add(end)
                    paths.append([*path,end])
    return details
