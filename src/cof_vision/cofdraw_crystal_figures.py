"""Finite illustrations unwrapped along the validated display graph's image shifts."""
from collections import deque
from functools import lru_cache

import numpy as np
from xml.etree import ElementTree as ET


def _circuit_svg(mol, coordinates, elements, ports):
    """Project the traced circuit using depth for crossing order and explicit labels."""
    xy=coordinates[:,:2]
    low,high=xy.min(0),xy.max(0)
    scale=min(740/max(high[0]-low[0],1),530/max(high[1]-low[1],1))
    points=(xy-(low+high)/2)*[scale,-scale]+[400,295]
    root=ET.Element("svg",xmlns="http://www.w3.org/2000/svg",viewBox="0 0 800 590",width="800",height="590")
    ink={"N":"#185a9b","O":"#b8313f","F":"#267d64","S":"#8b6522","B":"#9a6144"}
    groups=[]
    for bond in mol.GetBonds():
        a,b=bond.GetBeginAtomIdx(),bond.GetEndAtomIdx()
        bond_group=ET.Element("g",stroke="#263842",attrib={"stroke-linecap":"round","stroke-width":"1.35"})
        p,q=points[a],points[b]
        def stroke(u,v,extra=None):
            ET.SubElement(bond_group,"line",x1=f"{u[0]:.2f}",y1=f"{u[1]:.2f}",x2=f"{v[0]:.2f}",y2=f"{v[1]:.2f}",attrib=extra or {})
        stroke(p,q)
        if bond.GetBondTypeAsDouble()>1:
            vector=q-p
            side=np.array([-vector[1],vector[0]])/max(np.linalg.norm(vector),.1)*2.1
            stroke(p+.15*vector+side,q-.15*vector+side,{"stroke-dasharray":"2 2"} if bond.GetIsAromatic() else {})
        groups.append(((coordinates[a,2]+coordinates[b,2])/2,bond_group))
    for _,g in sorted(groups,key=lambda item:item[0]):
        root.append(g)
    for i,(x,y) in enumerate(points):
        el=elements[i] if i<ports else "∗"
        if el=="C":
            continue
        label=ET.SubElement(root,"text",x=f"{x:.2f}",y=f"{y+3:.2f}",fill=ink.get(el,"#37454a"),
            stroke="white",attrib={"paint-order":"stroke","stroke-width":"2.4","font-size":"9",
            "font-family":"Arial, sans-serif","text-anchor":"middle"})
        label.text=el
    return ET.tostring(root,encoding="unicode")


def chemical_drawing(structure, edges, *, circuit=None, scheme=False):
    """Depict an unwrapped, connected local chemical graph with continuation ports."""
    from rdkit import Chem
    from rdkit.Chem import rdDepictor
    from rdkit.Chem.Draw import rdMolDraw2D

    atoms = structure["atoms"]
    elements = [a["element"] for a in atoms]
    graph = [[] for _ in atoms]
    for edge in edges:
        a, b, t = edge["a"], edge["b"], tuple(edge["translation"])
        graph[a].append((b, t))
        graph[b].append((a, tuple(-n for n in t)))
    def adjacent(site):
        i, image = site
        return [(j, tuple(a+b for a, b in zip(image, shift))) for j, shift in graph[i]
                if elements[j] != "H"]
    # Start at a chemically distinctive site on a framework component. Preserve
    # image shifts: neighbors across opposite cell faces may be different atoms.
    def priority(i):
        el = elements[i]
        neighbors = [elements[j] for j, _ in graph[i]]
        if el == "H" or not any(e != "H" for e in neighbors):
            return -1
        if el in {"F", "Cl", "Br"} and neighbors == ["C"]:
            return 6
        if el in {"O", "S"} and "C" in neighbors:
            return 5
        if el in {"Si","B","Cu","Zn","Co","Ni"}:
            return 4
        return 3 if el == "N" else 2 if len(neighbors) == 4 and "H" not in neighbors else 1
    seed = max(range(len(atoms)), key=priority)
    root = (seed, (0, 0, 0))
    nodes, frontier = {root}, {root}
    for _ in range(9):
        frontier = {b for a in frontier for b in adjacent(a)} - nodes
        if len(nodes | frontier) > 95:
            break
        nodes |= frontier
    # Complete rings crossing the crop where the missing path is short.
    for _ in range(2):
        outer = {b for a in nodes for b in adjacent(a)} - nodes
        bridges = {b for b in outer if len(set(adjacent(b)) & nodes) >= 2}
        for b in outer:
            for c in set(adjacent(b)) & outer:
                if len((set(adjacent(b)) | set(adjacent(c))) & nodes) >= 2:
                    bridges.update((b, c))
        nodes |= bridges
    if circuit:
        nodes = {(i, tuple(image)) for i,image in circuit["site_images"]}
    ordered = sorted(nodes)
    rw = Chem.RWMol()
    ids = {}
    for site in ordered:
        i = site[0]
        atom = Chem.Atom(elements[i])
        atom.SetNoImplicit(True)
        if elements[i] != "C":
            atom.SetNumExplicitHs(sum(elements[j] == "H" for j, _ in graph[i]))
        ids[site] = rw.AddAtom(atom)
    positions=list(ordered)
    for a in ordered:
        for b in adjacent(a):
            if b not in nodes:
                port = Chem.Atom(0)
                port.SetProp("atomLabel", "∗")
                rw.AddBond(ids[a], rw.AddAtom(port), Chem.BondType.SINGLE)
                positions.append(b)
            elif a < b:
                rw.AddBond(ids[a], ids[b], Chem.BondType.SINGLE)
    mol = rw.GetMol()
    mol.UpdatePropertyCache(strict=False)
    # CIF coordinates alone do not determine resonance or bond orders. Mark
    # six-carbon aromatic-like rings only when geometry supports that choice.
    xyz = np.array([(np.array(atoms[i]["fractional"])+image) @ structure["cell"] for i, image in ordered])
    carbon_graph = {ids[a]: [ids[b] for b in adjacent(a) if b in nodes and elements[b[0]]=="C"]
                    for a in ordered if elements[a[0]]=="C"}
    rings = []
    for start in carbon_graph:
        paths = [[start]]
        while paths:
            path=paths.pop()
            if len(path)==6:
                if start in carbon_graph[path[-1]] and path[1]<path[-1]:
                    rings.append(path)
            else:
                paths.extend(path+[v] for v in carbon_graph[path[-1]] if v>start and v not in path)
    for r in rings:
        points = xyz[r]
        lengths = np.linalg.norm(points-np.roll(points, 1, axis=0), axis=1)
        if min(lengths) < 1.29 or max(lengths) > 1.49 or np.linalg.svd(points-points.mean(0))[1][-1] > .24:
            continue
        for a, b in zip(r, r[1:]+r[:1]):
            mol.GetBondBetweenAtoms(a,b).SetBondType(Chem.BondType.AROMATIC)
            mol.GetBondBetweenAtoms(a,b).SetIsAromatic(True)
            mol.GetAtomWithIdx(a).SetIsAromatic(True)
    for bond in mol.GetBonds():
        a,b = bond.GetBeginAtomIdx(),bond.GetEndAtomIdx()
        if max(a,b) >= len(ordered) or bond.GetIsAromatic():
            continue
        i,j=ordered[a][0],ordered[b][0]
        if elements[j]=="C":
            i,j=j,i
        if elements[i]!="C" or len(graph[i])!=3:
            continue
        length=np.linalg.norm(xyz[a]-xyz[b])
        if (elements[j]=="O" and len(graph[j])==1 and 1.10<length<1.32
                or elements[j]=="N" and len(graph[j])==2 and 1.18<length<1.34
                and sorted(elements[k] for k,_ in graph[i])==["C","H","N"]):
            bond.SetBondType(Chem.BondType.DOUBLE)
    known=structure.get("chemical_bonds", [])
    if known:
        from .cofdraw_periodic_build import edge_key
        orders={edge_key(b["a"],b["b"],b["translation"]):b["order"] for b in known}
        types={1:Chem.BondType.SINGLE,1.5:Chem.BondType.AROMATIC,2:Chem.BondType.DOUBLE,3:Chem.BondType.TRIPLE}
        for bond in mol.GetBonds():
            a,b=bond.GetBeginAtomIdx(),bond.GetEndAtomIdx()
            if max(a,b)>=len(ordered):
                continue
            (i,u),(j,v)=ordered[a],ordered[b]
            order=orders.get(edge_key(i,j,np.subtract(v,u).tolist()))
            if order in types:
                bond.SetBondType(types[order])
                bond.SetIsAromatic(order==1.5)
                if order==1.5:
                    bond.GetBeginAtom().SetIsAromatic(True)
                    bond.GetEndAtom().SetIsAromatic(True)
    if scheme and circuit and circuit.get("backbone"):
        from .cofdraw_reticular_scheme import draw_scheme
        def index(site):
            return ids[site[0],tuple(site[1:])]
        backbone=[index(site) for site in circuit["backbone"]]
        segments=[[index(site) for site in group] for group in circuit["linkers"]]
        return draw_scheme(mol,backbone,segments)
    if circuit:
        # Preserve the source circuit geometry. A generic 2D macrocycle layout
        # can fold the pore shut and enumerate an enormous set of ring layouts.
        normal=np.asarray(circuit.get("direction",[0.,0.,1.]))
        horizontal=np.cross(np.eye(3)[np.argmin(np.abs(normal))],normal)
        horizontal/=np.linalg.norm(horizontal)
        vertical=np.cross(normal,horizontal)
        coords=np.array([(np.array(atoms[i]["fractional"])+image)@structure["cell"] for i,image in positions])
        return _circuit_svg(mol, coords@np.array([horizontal,vertical,normal]).T,
                            [elements[i] for i,_ in ordered],len(ordered))
    else:
        rdDepictor.Compute2DCoords(mol, canonOrient=True)
    draw = rdMolDraw2D.MolDraw2DSVG(800, 590)
    draw.drawOptions().prepareMolsBeforeDrawing = False
    draw.drawOptions().padding = .055
    draw.DrawMolecule(mol)
    draw.FinishDrawing()
    return draw.GetDrawingText()


def readable_direction(model, face=None):
    """Prefer views that keep six-membered rings open and distinct atoms apart."""
    atoms = model["atoms"]
    xyz = np.array([[a[k] for k in ("x","y","z")] for a in atoms])
    graph = [[] for _ in atoms]
    for a,b in model["bonds"]:
        graph[a].append(b)
        graph[b].append(a)
    normals = []
    # Only six-membered rings inform this camera. Enumerate short paths directly
    # to avoid finding a global cycle basis on large, branched molecular motifs.
    for start in range(len(graph)):
        paths = [(start,)]
        while paths:
            path = paths.pop()
            if len(path)==6:
                if start in graph[path[-1]] and path[1]<path[-1]:
                    points = xyz[list(path)]
                    normals.append(np.linalg.svd(points-points.mean(axis=0),full_matrices=False)[2][-1])
            else:
                paths.extend(path+(n,) for n in graph[path[-1]] if n>start and n not in path)
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


def precursor_scenes(structure, edges):
    """Unwrap reacted molecular residues using their assigned precursor and image shifts."""
    atoms=structure["atoms"]
    if not edges or not all(a.get("component") in ("A","B") for a in atoms):
        return {}
    adjacency=[[] for _ in atoms]
    for e in edges:
        i,j=e["a"],e["b"]
        t=np.asarray(e["translation"])
        adjacency[i].append((j,t))
        adjacency[j].append((i,-t))
    cell=np.asarray(structure["cell"])
    f=np.array([a["fractional"] for a in atoms])
    group_of={}
    offsets={}
    residues=[]
    details={}
    for i,atom in enumerate(atoms):
        if i in group_of:continue
        group=len(residues)
        family=atom["component"]
        group_of[i]=group
        offsets[i]=np.zeros(3,dtype=int)
        queue=[i]
        members=[]
        while queue:
            a=queue.pop()
            members.append(a)
            for b,t in adjacency[a]:
                if atoms[b]["component"]!=family:continue
                shift=offsets[a]+t
                if b not in group_of:
                    group_of[b]=group
                    offsets[b]=shift
                    queue.append(b)
                elif not np.array_equal(offsets[b],shift):
                    return {}
        residues.append((family,members))
    centers=[np.mean([f[a]+offsets[a] for a in members if atoms[a]["element"]!="H"],axis=0) for family,members in residues]
    connections=[set() for _ in residues]
    for e in edges:
        a,b=e["a"],e["b"]
        u,v=group_of[a],group_of[b]
        if u==v:continue
        t=offsets[a]+e["translation"]-offsets[b]
        connections[u].add((v,*t.tolist()))
        connections[v].add((u,*(-t).tolist()))
    def patched(groups,links,topology=False):
        out_atoms=[]
        bonds=[]
        origin={}
        centroids={}
        for group,x,y,z in sorted(groups):
            tag=(group,x,y,z)
            image=np.array([x,y,z])
            family,members=residues[group]
            centroids[tag]=((centers[group]+image)@cell).tolist()
            if topology:
                origin[tag]=len(out_atoms)
                px,py,pz=centroids[tag]
                out_atoms.append(dict(x=px,y=py,z=pz,component=family,element="*"))
                continue
            for i in members:
                origin[tag,i]=len(out_atoms)
                px,py,pz=(f[i]+offsets[i]+image)@cell
                out_atoms.append(dict(element=atoms[i]["element"],x=px,y=py,z=pz,component=family))
        if topology:
            bonds=[[origin[a],origin[b]] for a,b in links]
        else:
            for group,x,y,z in sorted(groups):
                tag=(group,x,y,z)
                for a in residues[group][1]:
                    for b,t in adjacency[a]:
                        other=(group_of[b],*(offsets[a]+[x,y,z]+t-offsets[b]).tolist())
                        left=origin[tag,a]
                        right=origin.get((other,b))
                        if right is not None and left<right:bonds.append([left,right])
        return dict(atoms=out_atoms,bonds=bonds,spatial_guides=[] if topology else [[centroids[a],centroids[b]] for a,b in links])
    for i,(family,members) in enumerate(residues):
        if "block"+family not in details:
            details["block"+family]=patched({(i,0,0,0)},[])
    root=(max(range(len(residues)),key=lambda i:len(connections[i])),0,0,0)
    groups={root}
    frontier={root}
    links=set()
    for _ in range(2):
        upcoming=set()
        for tag in sorted(frontier):
            i,x,y,z=tag
            for j,u,v,w in connections[i]:
                other=(j,x+u,y+v,z+w)
                links.add(tuple(sorted((tag,other))))
                upcoming.add(other)
        frontier=upcoming-groups
        groups.update(upcoming)
    if links:
        details["blocknet"]=patched(groups,sorted(links),True)
        details["block_overlay"]=patched(groups,sorted(links))
    return details


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
    centers.update(structure.get("junction_sites",[]))
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
    details = {"junction": molecular(nodes),**precursor_scenes(structure,edges)}
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
            details["ring"]["site_images"] = [[i,[x,y,z]] for i,x,y,z in sorted(atoms)]
            details["ring"]["backbone"]=loop[:-1]
            details["ring"]["linkers"]=[sorted(arms(a)[b]) for a,b in zip(loop,loop[1:])]
            ring_index={site:i for i,site in enumerate(sorted(atoms))}
            details["cage_overlay"]={k:details["ring"][k] for k in ("atoms","bonds")}
            details["cage_overlay"]["guides"]=[[ring_index[a],ring_index[b]] for a,b in zip(loop,loop[1:])]
            points = np.array([(np.asarray(sites[i]["fractional"]) + [x,y,z]) @ cell for i,x,y,z in path+[root]])
            _, _, vectors = np.linalg.svd(points-points.mean(axis=0), full_matrices=False)
            readable_direction(details["ring"],vectors[-1])
            details["cage_overlay"]["direction"]=details["ring"].get("direction",vectors[-1].tolist())
            break
        if len(path) < 8:
            for end in sorted(arms(path[-1])):
                if end not in seen:
                    seen.add(end)
                    paths.append([*path,end])
    return details
