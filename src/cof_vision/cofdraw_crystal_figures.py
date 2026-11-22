"""Finite illustrations unwrapped along the validated display graph's image shifts."""
from collections import deque
from functools import lru_cache

import numpy as np
from xml.etree import ElementTree as ET


def _circuit_svg(mol, coordinates, elements, ports, blocks, rings):
    """Project the traced circuit using depth for crossing order and explicit labels."""
    xy=coordinates[:,:2]
    low,high=xy.min(0),xy.max(0)
    scale=15.
    w,h=np.maximum((high-low)*scale,1)+32
    points=(xy-(low+high)/2)*[scale,-scale]+[w/2,h/2]
    root=ET.Element("svg",xmlns="http://www.w3.org/2000/svg",viewBox=f"0 0 {w:.2f} {h:.2f}",width=f"{w:.2f}",height=f"{h:.2f}")
    if any(blocks):
        root.set("data-unassigned-sites",str(sum(not blocks[a.GetIdx()] and a.GetAtomicNum()>1 for a in mol.GetAtoms())))
    ink={"N":"#185a9b","O":"#b8313f","F":"#267d64","S":"#8b6522","B":"#9a6144"}
    depths=coordinates[:,2];low_z,high_z=np.percentile(depths,[5,95])
    depths=np.clip((depths-low_z)/(high_z-low_z),0,1) if high_z-low_z>.001 else np.ones(len(depths))
    washes=ET.SubElement(root,"g",attrib={"data-ring-washes":"", "display":"none"})
    for ring in rings:
        if all(mol.GetAtomWithIdx(a).GetIsAromatic() for a in ring):
            coords=" ".join(f"{points[a,0]:.2f},{points[a,1]:.2f}" for a in ring)
            ET.SubElement(washes,"polygon",points=coords,attrib={"data-depth":f"{np.mean(depths[ring]):.3f}",
                "data-block":blocks[ring[0]] if len({blocks[a] for a in ring})==1 else ""})
    groups=[]
    for bond in mol.GetBonds():
        a,b=bond.GetBeginAtomIdx(),bond.GetEndAtomIdx()
        bond_group=ET.Element("g",stroke="#263842",attrib={"stroke-linecap":"round","stroke-width":"1.35"})
        bond_group.set("data-depth",f"{(depths[a]+depths[b])/2:.3f}")
        bond_group.set("data-block",blocks[a] if blocks[a]==blocks[b] else "")
        p,q=points[a],points[b]
        def stroke(u,v,extra=None):
            ET.SubElement(bond_group,"line",x1=f"{u[0]:.2f}",y1=f"{u[1]:.2f}",x2=f"{v[0]:.2f}",y2=f"{v[1]:.2f}",attrib=extra or {})
        stroke(p,q)
        if bond.GetBondTypeAsDouble()>1:
            vector=q-p
            side=np.array([-vector[1],vector[0]])/max(np.linalg.norm(vector),.1)*2.1
            stroke(p+.15*vector+side,q-.15*vector+side,{"stroke-dasharray":"2 2"} if bond.GetIsAromatic() else {})
            if bond.GetBondTypeAsDouble()==3:stroke(p+.15*vector-side,q-.15*vector-side)
        groups.append(((coordinates[a,2]+coordinates[b,2])/2,bond_group))
    for _,g in sorted(groups,key=lambda item:item[0]):
        root.append(g)
    for i,(x,y) in enumerate(points):
        el=elements[i] if i<ports else "∗"
        atom=mol.GetAtomWithIdx(i)
        if el=="C" and not atom.GetIsotope():
            continue
        label=ET.SubElement(root,"text",x=f"{x:.2f}",y=f"{y+3:.2f}",fill=ink.get(el,"#37454a"),
            stroke="white",attrib={"paint-order":"stroke","stroke-width":"2.4","font-size":"9",
            "font-family":"Arial, sans-serif","text-anchor":"middle"})
        hydrogens=atom.GetNumExplicitHs()
        label.set("data-depth",f"{depths[i]:.3f}");label.set("data-block",blocks[i])
        label.text=(str(atom.GetIsotope()) if atom.GetIsotope() else "")+el+("H"+(str(hydrogens) if hydrogens>1 else "") if hydrogens else "")
    return ET.tostring(root,encoding="unicode")


def junction_patch(structure, edges):
    """One tetraaryl junction plus its directly reacted neighbors, in lattice images."""
    atoms = structure["atoms"]
    graph = [[] for _ in atoms]
    for edge in edges:
        i, j, t = edge["a"], edge["b"], edge["translation"]
        if "H" not in (atoms[i]["element"], atoms[j]["element"]):
            graph[i].append((j, *t))
            graph[j].append((i, *(-v for v in t)))
    root = next((i for i, near in enumerate(graph)
                 if atoms[i]["element"] in ("C", "Si") and len(near) == 4
                 and all(atoms[j]["element"] == "C" and len(graph[j]) == 3 for j, *_ in near)), None)
    if root is None:
        return None
    assigned = all(a.get("component") in ("A", "B") for a in atoms)
    @lru_cache(None)
    def boundary(i, j):
        if assigned:
            return atoms[i]["component"] != atoms[j]["component"]
        elements = {atoms[i]["element"], atoms[j]["element"]}
        if elements == {"C", "B"}:
            return True
        if elements == {"C", "N"}:
            carbon, nitrogen = (i,j) if atoms[i]["element"] == "C" else (j,i)
            if len(graph[carbon]) != 2 or len(graph[nitrogen]) != 2:
                return False
            # Pyridine and pyrrole C-N contacts belong to a small ring.
            # Only an exocyclic CH=N bond can separate imine residues.
            front={k for k,*_ in graph[carbon] if k!=nitrogen}
            seen={carbon}
            for _ in range(5):
                if nitrogen in front:return False
                seen.update(front)
                front={neighbor for atom in front for neighbor,*_ in graph[atom]}-seen
            return True
        return False
    center = (root, 0, 0, 0)
    distance = {center: 0}
    queue = deque([center])
    while queue:
        site = queue.popleft()
        i,x,y,z = site
        for j,u,v,w in graph[i]:
            other = (j,x+u,y+v,z+w)
            cost = distance[site] + int(boundary(i,j))
            if cost <= 1 and cost < distance.get(other, 2):
                distance[other] = cost
                queue.append(other)
        if len(distance) > 450:
            return None
    if 1 not in distance.values():
        return None
    # End at the aryl single bond before the next imine so a wavy single
    # continuation never replaces a known C=N double bond.
    trim=[]
    for site,cost in distance.items():
        if cost!=1:continue
        i,x,y,z=site
        near=[(j,x+u,y+v,z+w) for j,u,v,w in graph[i]]
        if len(near)!=2:continue
        outside=[v for v in near if v not in distance]
        if len(outside)==1 and {atoms[i]['element'],atoms[outside[0][0]]['element']}=={'C','N'}:
            trim.append(site)
    for site in trim:del distance[site]
    return {"center": center, "site_images": [[i,[x,y,z]] for i,x,y,z in sorted(distance)],
            "outer_sites": [s for s,cost in distance.items() if cost], "direction": [0.,0.,1.]}


def _reference_components(record, graph, elements, edges):
    """Map atoms in a source phase to recorded precursors across the full cell."""
    from collections import Counter
    from rdkit import Chem
    linkage=record.get('linkage')
    if linkage=='boroxine':return ['A']*len(elements)
    if linkage not in ('boronate_ester','imine'):return ['']*len(elements)
    first,second=('B','O') if linkage=='boronate_ester' else ('N','C')
    letters={first:'A',second:'B'}
    formulas={}
    if linkage=='imine':
        smiles=record.get('precursor_smiles',[])
        if len(smiles)!=2:return ['']*len(elements)
        letters={}
        amine=Chem.MolFromSmarts('[c]-[NH2]')
        aldehyde=Chem.MolFromSmarts('[c]-[CH1]=O')
        for letter,s in zip(('A','B'),smiles):
            precursor=Chem.MolFromSmiles(s)
            if precursor is None:continue
            n=len(precursor.GetSubstructMatches(amine))
            c=len(precursor.GetSubstructMatches(aldehyde))
            if bool(n)==bool(c):continue
            letters['N' if n else 'C']=letter
            formula=Counter(a.GetSymbol() for a in precursor.GetAtoms() if a.GetSymbol()!='H')
            formula['O']-=c
            formulas[letter]=+formula
        if len(letters)!=2:return ['']*len(elements)
    boundaries=set()
    seeds={}
    for row in edges:
        i,j=row['a'],row['b']
        if {elements[i],elements[j]}!={first,second}:continue
        if linkage=='imine':
            carbon,nitrogen=(i,j) if elements[i]=='C' else (j,i)
            heavy=lambda k:[p for p,_ in graph[k] if elements[p]!='H']
            if (sorted(elements[p] for p in heavy(carbon))!=['C','N'] or
                    sorted(elements[p] for p in heavy(nitrogen))!=['C','C'] or
                    not 1.18<row.get('distance_angstrom',0.)<1.36):continue
            # Exclude endocyclic aromatic C=N, including pyridine and porphyrin.
            seen={carbon}
            front=set(heavy(carbon))-{nitrogen}
            for _ in range(6):
                if nitrogen in front:break
                seen.update(front)
                front={p for q in front for p in heavy(q)}-seen
            if nitrogen in front:continue
        boundaries.add(tuple(sorted((i,j))))
        seeds[i]=letters[elements[i]];seeds[j]=letters[elements[j]]
    owners={}
    for seed in range(len(elements)):
        if seed in owners:continue
        group={seed};queue=[seed]
        while queue:
            i=queue.pop()
            for j,_ in graph[i]:
                if j not in group and tuple(sorted((i,j))) not in boundaries:
                    group.add(j);queue.append(j)
        evidence={seeds[i] for i in group if i in seeds}
        letter=evidence.pop() if len(evidence)==1 else ''
        if letter and formulas and Counter(elements[i] for i in group if elements[i]!='H')!=formulas[letter]:
            letter=''
        owners.update(dict.fromkeys(group,letter))
    return [owners[i] for i in range(len(elements))]


def chemical_drawing(structure, edges, *, circuit=None, scheme=False, arm_spread=0):
    """Depict an unwrapped, connected local chemical graph with continuation ports."""
    from rdkit import Chem

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
        if atoms[i].get("isotope"):
            atom.SetIsotope(int(atoms[i]["isotope"]))
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
    blocks=[atoms[i].get("component","") for i,_ in positions]
    record=structure.get("provenance",{}).get("drawing_cofex",{})
    if not any(blocks) and record:
        owners=_reference_components(record,graph,elements,edges)
        blocks=[owners[i] for i,_ in positions]
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
    if circuit and "center" in circuit:
        from .cofdraw_reticular_scheme import junction_scheme
        i,*image = circuit["center"]
        return junction_scheme(mol, ids[i,tuple(image)], blocks, arm_spread)
    if scheme and circuit and circuit.get("backbone"):
        from .cofdraw_reticular_scheme import draw_scheme
        def index(site):
            return ids[site[0],tuple(site[1:])]
        backbone=[index(site) for site in circuit["backbone"]]
        segments=[[index(site) for site in group] for group in circuit["linkers"]]
        # Retain real crystal depths even when rigid linker fragments are laid
        # out with conventional bond angles in the plane of the drawing.
        normal=np.asarray(circuit["direction"])
        depths=[((np.asarray(atoms[i]["fractional"])+image)@structure["cell"])@normal for i,image in positions]
        rings=[([index(s) for s in r["backbone"]],[[index(s) for s in group] for group in r["linkers"]])
               for r in circuit.get("rings",[])]
        anchors=None
        if rings:
            u=np.cross(np.eye(3)[np.argmin(np.abs(normal))],normal);u/=np.linalg.norm(u)
            basis=np.column_stack([u,np.cross(normal,u)])
            anchors={j:((np.array(atoms[i]["fractional"])+t)@structure["cell"])@basis for (i,t),j in ids.items()}
        return draw_scheme(mol,backbone,segments,depths,blocks,rings,anchors)
    if circuit:
        # Preserve the source circuit geometry. A generic 2D macrocycle layout
        # can fold the pore shut and enumerate an enormous set of ring layouts.
        normal=np.asarray(circuit.get("direction",[0.,0.,1.]))
        horizontal=np.cross(np.eye(3)[np.argmin(np.abs(normal))],normal)
        horizontal/=np.linalg.norm(horizontal)
        vertical=np.cross(normal,horizontal)
        coords=np.array([(np.array(atoms[i]["fractional"])+image)@structure["cell"] for i,image in positions])
        projected=coords@np.array([horizontal,vertical,normal]).T
        if scheme:
            from .cofdraw_reticular_scheme import reference_scheme
            return reference_scheme(mol,projected,blocks,rings,circuit.get('scheme_tilt',0))
        return _circuit_svg(mol,projected,[elements[i] for i,_ in ordered],len(ordered),blocks,rings)
    from .cofdraw_reticular_scheme import fragment_layout,vector_art
    coords=fragment_layout(mol,range(mol.GetNumAtoms()))
    return vector_art(mol,coords,blocks=blocks)


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


def adjacent_circuits(ring, cell, count=3):
    """Grow a compact patch of circuits sewn together at real periodic edges."""
    if count not in (1,2,3,4,7):raise ValueError("Choose 1, 2, 3, 4 or 7 circuits")
    if count==1:return ring
    loop=[tuple(site) for site in ring["backbone"]]
    if ring.get("junction_arms"):
        return _chair_patch(ring, loop, count, np.asarray(cell))
    edges={tuple(sorted((a,b))) for a,b in zip(loop,loop[1:]+loop[:1])}
    def translate(site,shift):
        return (site[0],*(a+b for a,b in zip(site[1:],shift)))
    candidates={}
    for a in loop:
        for b in loop:
            if a==b or a[0]!=b[0]:continue
            shift=tuple(x-y for x,y in zip(a[1:],b[1:]))
            moved={tuple(sorted((translate(u,shift),translate(v,shift)))) for u,v in edges}
            if len(moved&edges)==1:candidates[shift]=moved
    if not candidates:return None
    shifts=[(0,0,0)]
    while len(shifts)<count:
        frontier={tuple(a+b for a,b in zip(p,q)) for p in shifts for q in candidates}-set(shifts)
        if not frontier:return None
        def score(q):
            shared=sum(tuple(a-b for a,b in zip(q,p)) in candidates for p in shifts)
            spread=sum(np.linalg.norm(np.subtract(q,p)@cell) for p in shifts)
            return -shared,spread,q
        shifts.append(min(frontier,key=score))
    sites={(i,*shift) for i,shift in ring["site_images"]}
    union={translate(s,d) for s in sites for d in shifts}
    rings=[{"backbone":[translate(s,d) for s in loop],
            "linkers":[[translate(s,d) for s in arm] for arm in ring["linkers"]]} for d in shifts]
    return {**rings[0],"site_images":[[site[0],list(site[1:])] for site in sorted(union)],
            "direction":ring["direction"],"circuits":count,"rings":rings}


def _chair_patch(ring, first, count, cell):
    """Follow different chair faces through the same real tetrahedral cage."""
    basis=ring["junction_arms"]
    def move(site,origin):
        return (site[0],*(site[k]+origin[k] for k in (1,2,3)))
    @lru_cache(None)
    def neighbors(site):
        return {move(end,site):arm for end,arm in basis[str(site[0])]}
    def edges(loop):
        return {tuple(sorted(pair)) for pair in zip(loop,loop[1:]+loop[:1])}
    def xyz(v):
        return (np.array(ring["junction_fractional"][str(v[0])])+v[1:])@cell
    normal=np.array(ring["direction"])
    focus=np.mean([xyz(v) for v in first],axis=0)
    used=edges(first)
    loops=[first]
    seen={frozenset(first)}
    while len(loops)<count:
        candidates={}
        for a,b in sorted(used):
            pending=[(a,b)]
            while pending:
                path=pending.pop()
                if len(path)==len(first):
                    key=frozenset(path)
                    if a in neighbors(path[-1]) and key not in seen and not edges(list(path))<=used:
                        candidates.setdefault(key,list(path))
                    continue
                pending.extend((*path,v) for v in neighbors(path[-1]) if v not in path)
        if not candidates:return None
        def rank(path):
            offset=np.mean([xyz(v) for v in path],axis=0)-focus
            shared=len(edges(path)&used)
            # Prefer the folded rear face of the same cage over an outward tile.
            return (-shared, round(np.linalg.norm(offset),5), round(offset@normal,5),path)
        chosen=min(candidates.values(),key=rank)
        loops.append(chosen);seen.add(frozenset(chosen));used.update(edges(chosen))
    sites={(i,*shift) for i,shift in ring["site_images"]}
    rings=[]
    for loop in loops:
        arms=[[move(v,a) for v in neighbors(a)[b]] for a,b in zip(loop,loop[1:]+loop[:1])]
        sites.update(loop)
        for node in loop:
            sites.update(move(v,node) for v in ring.get("junction_caps",{}).get(str(node[0]),[]))
        for arm in arms:sites.update(arm)
        rings.append(dict(backbone=loop,linkers=arms))
    # View between the normals of the front and rear faces. A face-on
    # projection of the first chair can hide a whole rear linker behind it.
    faces=[np.array([xyz(v) for v in loop]) for loop in loops]
    normals=np.array([np.linalg.svd(f-f.mean(0))[2][-1] for f in faces])
    vertex_ids=sorted({v for loop in loops for v in loop})
    vertices=np.array([xyz(v) for v in vertex_ids])
    i,j=np.triu_indices(len(vertices),1)
    delta=vertices[i]-vertices[j]
    length=max(np.linalg.norm(faces[0][1]-faces[0][0]),1.)
    k=np.arange(256)+.5
    z=1-2*k/256;r=np.sqrt(1-z*z)
    cameras=np.column_stack((r*np.cos(k*2.399963),r*np.sin(k*2.399963),z))
    gaps=np.sqrt(np.maximum((delta*delta).sum(1)[:,None]-(delta@cameras.T)**2,0))/length
    exposure=np.abs(normals@cameras.T)
    # Check angular separation at branch points as well as vertex distances.
    # Distinct endpoints may still project their phenyl groups on top of each other.
    cosine=np.full(len(cameras),-1.)
    for center in {v for edge in used for v in edge}:
        rays=[xyz(b if a==center else a)-xyz(center) for a,b in used if a==center or b==center]
        for i,left in enumerate(rays):
            for right in rays[:i]:
                lz,rz=cameras@left,cameras@right
                lengths=np.sqrt(np.maximum((left@left-lz**2)*(right@right-rz**2),1e-10))
                cosine=np.maximum(cosine,(left@right-lz*rz)/lengths)
    clearance=np.arccos(np.clip(cosine,-1,1))
    arm_span=np.ones(len(cameras))
    edge_clearance=np.ones(len(cameras))
    for a,b in used:
        d=xyz(b)-xyz(a)
        dz=cameras@d
        span=np.maximum(d@d-dz**2,1e-9)
        arm_span=np.minimum(arm_span,np.sqrt(span/(d@d)))
        points=np.array([xyz(v)-xyz(a) for v in vertex_ids if v not in (a,b)])
        pz=points@cameras.T
        along=points@d
        t=np.clip((along[:,None]-pz*dz)/span,0,1)
        squared=(points*points).sum(1)[:,None]-pz**2-2*t*(along[:,None]-pz*dz)+t**2*span
        edge_clearance=np.minimum(edge_clearance,np.sqrt(np.maximum(0,squared)).min(0)/length)
    scores=(exposure.min(0)+.6*np.minimum(gaps.min(0),.5)+.15*exposure.mean(0)
            +.9*np.minimum(clearance,1.25)+2*np.minimum(arm_span,.74)+2*np.minimum(edge_clearance,.43))
    from .cofdraw_crystal_circuits import bridge_overlap_cost
    links=[(vertex_ids.index(a),vertex_ids.index(b)) for a,b in sorted(used)]
    scores-=2*bridge_overlap_cost(vertices,links,cameras)
    view=cameras[scores.argmax()]
    if view@normal<0:view=-view
    return {**rings[0],"rings":rings,"circuits":count,"direction":view.tolist(),
            "site_images":[[s[0],list(s[1:])] for s in sorted(sites)]}


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
        # Keep fractional offsets so both skewed cells and supercell displays
        # can reconstruct the whole net without guessing across face boundaries.
        positions=[c-np.floor(c) for c in centers]
        spans=[]
        for i,neighbors in enumerate(connections):
            for j,x,y,z in sorted(neighbors):
                if (i,0,0,0)>=(j,x,y,z):continue
                shift=(np.floor(centers[j])-np.floor(centers[i])+[x,y,z]).astype(int)
                spans.append([i,j,*shift.tolist()])
        details["periodic_topology"]={"points":[p.tolist() for p in positions],"edges":spans,
                                     "kind":"Reacted building-block centers"}
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
    site_ids=sorted(centers)
    site_index={site:i for i,site in enumerate(site_ids)}
    spans=[]
    for i in site_ids:
        site=(i,0,0,0)
        ends=arms(site)
        if not ends:break
        spans.extend([site_index[i],site_index[j],x,y,z] for j,x,y,z in sorted(ends) if site<(j,x,y,z))
    else:
        details["periodic_topology"]={"points":[sites[i]["fractional"] for i in site_ids],
                                     "site_ids":site_ids,"edges":spans,"kind":"Traced junction-to-junction paths"}
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
            details["ring"]["junction_arms"]={str(i):[[end,sorted(arm)] for end,arm in arms((i,0,0,0)).items()] for i in site_ids}
            caps={}
            for i in site_ids:
                crop=front={(i,0,0,0)}
                for _ in range(4):
                    front={q for p in front for q in neighbors(p)}-crop
                    crop=crop|front
                caps[str(i)]=sorted(crop)
            details["ring"]["junction_caps"]=caps
            details["ring"]["junction_fractional"]={str(i):sites[i]["fractional"] for i in site_ids}
            details["ring"]["linkers"]=[sorted(arms(a)[b]) for a,b in zip(loop,loop[1:])]
            ring_index={site:i for i,site in enumerate(sorted(atoms))}
            details["cage_overlay"]={k:details["ring"][k] for k in ("atoms","bonds")}
            details["cage_overlay"]["guides"]=[[ring_index[a],ring_index[b]] for a,b in zip(loop,loop[1:])]
            points = np.array([(np.asarray(sites[i]["fractional"]) + [x,y,z]) @ cell for i,x,y,z in path+[root]])
            _, _, vectors = np.linalg.svd(points-points.mean(axis=0), full_matrices=False)
            readable_direction(details["ring"],vectors[-1])
            details["cage_overlay"]["direction"]=details["ring"].setdefault("direction",vectors[-1].tolist())
            break
        if len(path) < 8:
            for end in sorted(arms(path[-1])):
                if end not in seen:
                    seen.add(end)
                    paths.append([*path,end])
    return details
