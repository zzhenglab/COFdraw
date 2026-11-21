"""Construct diamond networks from reacted molecular graphs and a 3D edge conformer."""
from collections import Counter, defaultdict
from itertools import product
from math import acos, degrees, gcd

import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem

from .framework_drawing import LINKAGE_LIBRARY, _ports, canonical, linkage_fragment
from .cofdraw_network import inspect_periodic_network
from .cofdraw_periodic_build import edge_key
from .cofdraw_science import MAX_ATOMS, _PeriodicDistances, normalize_structure


SIGNS = np.array([[1,1,1], [1,-1,-1], [-1,1,-1], [-1,-1,1]])


def _components(smiles, linkage):
    if linkage not in ("imine", "secondary_amine", "acylhydrazone"):
        raise ValueError("The diamond builder supports imine, reduced imine and acylhydrazone products.")
    roles = LINKAGE_LIBRARY[linkage]["roles"]
    for ordered_roles in (roles, roles[::-1]):
        try:
            fragments = [linkage_fragment(s, linkage, role) for s,role in zip(smiles, ordered_roles)]
        except ValueError:
            continue
        counts = [f["connectivity"] for f in fragments]
        if sorted(counts) == [2,4]:
            return fragments, counts.index(4)
        if counts == [4,4]:
            for fragment in fragments:
                _symmetric_arm(Chem.MolFromSmiles(fragment['fragment_smiles']))
            return fragments, 0
    raise ValueError("Choose a four-arm tetrahedral precursor and a complementary two-ended or four-arm tetrahedral precursor.")


def _symmetric_arm(node):
    """Split only a true four-way articulation with four graph-identical reactive branches."""
    for a in node.GetAtoms():
        if a.GetAtomicNum() not in (6,14) or a.GetDegree() != 4 or a.IsInRing():
            continue
        center = a.GetIdx()
        arms = []
        for nb in a.GetNeighbors():
            sites, todo = {center}, [nb.GetIdx()]
            while todo:
                site = todo.pop()
                if site in sites:
                    continue
                sites.add(site)
                todo.extend(b.GetIdx() for b in node.GetAtomWithIdx(site).GetNeighbors() if b.GetIdx() not in sites)
            mol = Chem.RWMol(node)
            for i in reversed(range(mol.GetNumAtoms())):
                if i not in sites:
                    mol.RemoveAtom(i)
            arm = mol.GetMol()
            ports = _ports(arm)
            if len(ports) != 1:
                break
            dummy, _ = next(iter(ports.values()))
            arm.GetAtomWithIdx(dummy).SetAtomMapNum(1)
            c = arm.GetAtomWithIdx(sum(i < center for i in sites))
            c.SetBoolProp("_diaCenter", True)
            c.SetNumExplicitHs(3)
            c.SetNoImplicit(True)
            Chem.SanitizeMol(arm)
            arms.append(arm)
        if len(arms) == 4 and len({Chem.MolToSmiles(m) for m in arms}) == 1:
            return arms[0]
    raise ValueError("A diamond net requires four equivalent reactive arms on a tetrahedral carbon or silicon center. A planar four-site node needs another topology.")


def _edge_molecule(node, linker, order, node_slot):
    arm = _symmetric_arm(node)
    paired_hubs = len(_ports(linker)) == 4
    if paired_hubs:
        linker = _symmetric_arm(linker)
    for mol, component in ((arm,node_slot), (linker,1-node_slot)):
        for atom in mol.GetAtoms():
            atom.SetProp("_component", "AB"[component])
    arm_size, linker_size = arm.GetNumAtoms(), linker.GetNumAtoms()
    pa, ea = next(iter(_ports(arm).values()))
    paired_ports = list(_ports(linker).values())
    pair = Chem.CombineMols(arm,linker)
    edge = Chem.RWMol(pair if paired_hubs else Chem.CombineMols(pair,arm))
    excluded = []
    for off,(pl,el) in zip((0,arm_size+linker_size),paired_ports):
        a,b = off+ea, arm_size+el
        edge.AddBond(a,b,Chem.BondType.DOUBLE if order == 2 else Chem.BondType.SINGLE)
        edge.GetBondBetweenAtoms(a,b).SetBoolProp("_productJoin",True)
        excluded.extend((off+pa, arm_size+pl))
    for i in sorted(excluded,reverse=True):
        edge.RemoveAtom(i)
    edge = edge.GetMol()
    Chem.SanitizeMol(edge)
    # Choose extended E geometry at the new C=N junctions for this starting conformer.
    for bond in edge.GetBonds():
        if order == 2 and bond.HasProp("_productJoin"):
            a,b = bond.GetBeginAtom(), bond.GetEndAtom()
            left = next(n.GetIdx() for n in a.GetNeighbors() if n.GetIdx()!=b.GetIdx())
            right = next(n.GetIdx() for n in b.GetNeighbors() if n.GetIdx()!=a.GetIdx())
            bond.SetStereoAtoms(left,right)
            bond.SetStereo(Chem.BondStereo.STEREOE)
    Chem.SetDoubleBondNeighborDirections(edge)
    edge = Chem.AddHs(edge)
    for atom in edge.GetAtoms():
        if not atom.HasProp("_component"):
            atom.SetProp("_component",atom.GetNeighbors()[0].GetProp("_component"))
    return edge


def _unit(v):
    return v / np.linalg.norm(v)


def _basis(v):
    v = _unit(v)
    w = _unit(np.cross(v, np.eye(3)[np.argmin(abs(v))]))
    return np.array([v,w,np.cross(v,w)])


def _edge_orientations(xyz, ends, neighbours):
    start, end = ends
    v = _unit(xyz[neighbours[0]]-xyz[start])
    w = _unit(xyz[neighbours[1]]-xyz[end])
    displacement = xyz[end]-xyz[start]
    axis = _unit(v-w+_unit(displacement))
    basis = _basis(axis)
    target = _basis(np.ones(3))
    candidates = []
    for phi in np.arange(0, np.pi*2/3, np.pi/45):
        cs,sn = np.cos(phi),np.sin(phi)
        transform = basis.T @ np.array([[1,0,0],[0,cs,-sn],[0,sn,cs]]) @ target
        # Each junction sees four symmetry-related arm directions in the periodic cell.
        angle_error = max(abs(degrees(acos(np.clip((d@transform) @ ((d@transform)*s),-1,1)))-109.4712)
                          for s in SIGNS[1:] for d in (v,w))
        d = _unit(displacement@transform)
        if min(d) < .24:
            continue
        candidates.append((angle_error+2*np.ptp(d),angle_error,transform))
    for _,angle,transform in sorted(candidates,key=lambda c:c[0]):
        if angle <= 18:
            yield angle,(xyz-xyz[start])@transform


def _assemble_edges(edge, ends, cell, junctions, placements):
    """Share the tetrahedral atom at each lattice vertex; retain every other edge atom."""
    inv = np.linalg.inv(cell)
    cap_h = {h.GetIdx() for i in ends for h in edge.GetAtomWithIdx(i).GetNeighbors() if h.GetAtomicNum()==1}
    atom_rows, bonds = [], []
    def site(i, f):
        atom = edge.GetAtomWithIdx(i)
        atom_rows.append(dict(element=atom.GetSymbol(), fractional=f.tolist(), component=atom.GetProp("_component"),
                              **({"isotope":atom.GetIsotope()} if atom.GetIsotope() else {})))
        return len(atom_rows)-1
    for i,origin in enumerate(junctions):
        f=origin@inv
        site(ends[i%2],f-np.floor(f+1e-10))
    for start,end,coords in placements:
        frac = coords @ inv
        wraps = np.floor(frac+1e-10).astype(int)
        ids = {ends[0]:start,ends[1]:end}
        for i,f in enumerate(frac-wraps):
            if i not in cap_h and i not in ids:
                ids[i] = site(i,f)
        for b in edge.GetBonds():
            i,j=b.GetBeginAtomIdx(),b.GetEndAtomIdx()
            if i in ids and j in ids:
                bonds.append(dict(a=ids[i],b=ids[j],translation=(wraps[j]-wraps[i]).tolist(),order=b.GetBondTypeAsDouble()))
    return normalize_structure(dict(periodic=True, cell=cell.tolist(), atoms=atom_rows)), bonds


def _diamond_candidates(edge, ends, neighbours):
    xyz = edge.GetConformer().GetPositions()
    for angle,pos in _edge_orientations(xyz,ends,neighbours):
        vectors = SIGNS * pos[ends[1]]
        cell = vectors[0] - vectors[1:]
        structure,bonds = _assemble_edges(edge, ends, cell, [pos[ends[0]],pos[ends[1]]], [(0,1,pos*s) for s in SIGNS])
        yield angle, structure, bonds


def _junction_sequence(bonds, hubs):
    """Contract actual bonded bridges between hubs, carrying periodic image indices."""
    links = defaultdict(list)
    for b in bonds:
        t = tuple(b["translation"])
        links[b["a"]].append((b["b"],t))
        links[b["b"]].append((b["a"],tuple(-x for x in t)))
    graph = defaultdict(set)
    hubs = set(hubs)
    for h in hubs:
        start = (h,(0,0,0))
        visited, stack = {start},[start]
        while stack:
            i,offset = stack.pop()
            for j,step in links[i]:
                dest = (j,tuple(a+b for a,b in zip(offset,step)))
                if dest not in visited:
                    visited.add(dest)
                    if j in hubs:
                        graph[h].add(dest)
                    else:
                        stack.append(dest)
                if len(visited)>4*len(links):
                    raise ValueError("A bridge forms an unexpected periodic path without a tetrahedral junction.")
    if any(len(graph[h])!=4 for h in hubs):
        raise ValueError("An inferred junction does not connect four separate neighboring junctions.")
    shell = {(min(hubs),(0,0,0))}
    visited = shell.copy()
    sequence = []
    for _ in range(5):
        shell = {(j,tuple(a+b for a,b in zip(t,s))) for h,t in shell for j,s in graph[h]}-visited
        sequence.append(len(shell))
        visited.update(shell)
    return sequence


def _pack(structure, bonds, count, hubs):
    """Search distinct cyclic translation groups of the requested order in the primitive cell."""
    atoms, cell = structure["atoms"], np.array(structure["cell"])
    if count*len(atoms)>MAX_ATOMS:
        raise ValueError(f"The requested {count} nets exceed the {MAX_ATOMS}-atom primitive-cell limit.")
    f = np.array([a["fractional"] for a in atoms])
    radii = np.array([Chem.GetPeriodicTable().GetRvdw(a["element"]) for a in atoms])
    # Only close pairs can violate the vdW gate. Fractional pruning is exact
    # for that distance when the search box is narrower than half each cell height.
    cutoff = 2*radii.max()+1.
    widths = cutoff*np.linalg.norm(np.linalg.inv(cell),axis=0)
    delta = f[:,None,:]-f[None,:,:]
    pairs = .75*(radii[:,None]+radii[None,:])
    cache = {}
    def clearance(step):
        neg = tuple((-x)%count for x in step)
        k = min(step,neg)
        if k not in cache:
            shift = np.array(k)/count
            if max(widths) < .5:
                dx = delta+shift
                dx -= np.rint(dx)
                near = (np.abs(dx)<=widths).all(axis=2)
                cache[k] = float(np.min(np.linalg.norm(dx[near]@cell,axis=1)-pairs[near], initial=1.))
            else:
                finder=_PeriodicDistances(cell,f,.75*radii)
                cache[k] = float(np.min(finder.clearance(f+shift)-.75*radii))
        return cache[k]
    groups = set()
    for vector in product(range(count), repeat=3):
        if gcd(count,*vector)==1:
            groups.add(tuple(sorted(tuple(i*x%count for x in vector) for i in range(count))))
    ranked = []
    for group in sorted(groups):
        margins = []
        for step in group[1:]:
            margin = clearance(step)
            if margin < 0:
                break
            margins.append(margin)
        else:
            ranked.append((min(margins),group))
    for margin, group in sorted(ranked,reverse=True):
        sites, connected, junctions = [], [], []
        for copy, step in enumerate(group):
            frac = f+np.array(step)/count
            wrap = np.floor(frac).astype(int)
            off = len(sites)
            junctions.extend(off+hub for hub in hubs)
            sites.extend(dict(a,fractional=xyz.tolist(),network=copy) for a,xyz in zip(atoms,frac-wrap))
            for b in bonds:
                i,j=b["a"],b["b"]
                connected.append(dict(b,a=i+off,b=j+off,translation=(b["translation"]+wrap[j]-wrap[i]).tolist()))
        multi = normalize_structure(dict(periodic=True,cell=cell.tolist(),atoms=sites))
        check = inspect_periodic_network(multi)
        keys=lambda items:{edge_key(b["a"],b["b"],b["translation"]) for b in items}
        if (check["quotient_components_by_dimension"] == {"3":count} and not check["short_contact_count"]
            and not check["threshold_sensitivity"]["classification_changed"] and keys(connected)==keys(check["edges"])):
            return multi,connected,junctions,dict(fractional_offsets=(np.array(group)/count).tolist(),
                       vdw_fraction_gate=.75,intercopy_clearance_lower_bound_angstrom=float(np.floor(margin*10000)/10000))
    return None


def build_dia(monomer_a, monomer_b, linkage="imine", label_a="Monomer A", label_b="Monomer B", network_count=1, topology="dia"):
    if type(network_count) is not int or not 1<=network_count<=7:
        raise ValueError("Choose one to seven networks.")
    if topology not in ("dia","lon"):
        raise ValueError("Choose the cubic diamond (dia) or hexagonal diamond (lon) net.")
    hub_count = 2 if topology=="dia" else 4
    smiles = [canonical(s) for s in (monomer_a,monomer_b)]
    fragments, node_slot = _components(smiles,linkage)
    mols = [Chem.MolFromSmiles(f["fragment_smiles"]) for f in fragments]
    if any(len(Chem.GetMolFrags(m))!=1 or m.GetNumAtoms()>120 or
           any(a.GetFormalCharge() or a.GetNumRadicalElectrons() for a in m.GetAtoms()) for m in mols):
        raise ValueError("Each precursor must be one neutral closed-shell molecule with at most 120 fragment atoms.")
    two_nodes = all(f['connectivity']==4 for f in fragments)
    stoichiometry = [hub_count//2]*2 if two_nodes else [hub_count*(1 if i==node_slot else 2) for i in range(2)]
    counts = Counter()
    for i, mol in enumerate(mols):
        for a in Chem.AddHs(mol).GetAtoms():
            if a.GetAtomicNum():
                counts[a.GetSymbol(), a.GetIsotope()] += stoichiometry[i]
    if sum(counts.values()) > 650:
        raise ValueError("The diamond primitive cell exceeds the 650-atom build budget.")
    edge = _edge_molecule(mols[node_slot],mols[1-node_slot],LINKAGE_LIBRARY[linkage]["join_bond_order"],node_slot)
    ends = [a.GetIdx() for a in edge.GetAtoms() if a.HasProp("_diaCenter")]
    neighbours = [next(a.GetIdx() for a in edge.GetAtomWithIdx(i).GetNeighbors() if a.GetAtomicNum()!=1) for i in ends]
    best = "no extended tetrahedral conformer"
    tried_packing = False
    from .cofdraw_lon import lon_candidates
    cell_candidates = _diamond_candidates if topology=="dia" else lon_candidates
    for seed in (49, 511, 3289, 7717, 27011):
        params=AllChem.ETKDGv3()
        params.randomSeed=seed
        params.maxIterations=65
        if AllChem.EmbedMolecule(edge,params)<0:
            continue
        if AllChem.MMFFHasAllMoleculeParams(edge):
            converged=AllChem.MMFFOptimizeMolecule(edge,maxIters=450)
        elif AllChem.UFFHasAllMoleculeParams(edge):
            converged=AllChem.UFFOptimizeMolecule(edge,maxIters=450)
        else:
            raise ValueError("The molecular edge lacks supported conformer force-field parameters.")
        if converged!=0:
            continue
        for angle,structure,bonds in cell_candidates(edge,ends,neighbours):
            net = inspect_periodic_network(structure)
            expected = {edge_key(b["a"],b["b"],b["translation"]) for b in bonds}
            observed = {edge_key(b["a"],b["b"],b["translation"]) for b in net["edges"]}
            if (expected == observed and len(expected)==len(bonds) and net["dimensions"]==[3]
                    and not net["short_contact_count"] and not net["threshold_sensitivity"]["classification_changed"]):
                break
            best=f"{len(expected-observed)} missing bonds, {len(observed-expected)} extra bonds and {net['short_contact_count']} close contacts"
        else:
            continue
        if Counter((a["element"],a.get("isotope",0)) for a in structure["atoms"]) != counts:
            raise ValueError("The constructed cell failed the independent reacted-monomer composition check.")
        junctions=list(range(hub_count))
        sequence = _junction_sequence(net["edges"],junctions)
        if sequence != {"dia":[4,12,24,42,64],"lon":[4,12,25,44,67]}[topology]:
            raise ValueError("The periodic junction graph disagrees with the specified diamond net.")
        pack={}
        if network_count>1:
            tried_packing = True
            packed = _pack(structure,bonds,network_count,junctions)
            if packed is None:
                continue
            structure,bonds,junctions,pack=packed
        scope = (f"Unrelaxed {topology} hypothesis with {network_count} chosen network(s) generated from the selected precursors. "
                 "Tetrahedral junction angles, periodic bonds and rank three for each network were checked. Experimental packing and interpenetration are unspecified.")
        structure.update(chemical_bonds=bonds,junction_sites=junctions)
        structure["provenance"] = dict(kind="monomer_dia_design",model_name=f"{label_a} + {label_b}",net=topology,
                network_dimension=3,network_count=network_count,precursors=smiles,linkage=linkage,node_index=node_slot,
                precursor_connectivities=[f["connectivity"] for f in fragments],
                packing=pack,
                conformer_seed=seed,product_geometry="Extended E starting geometry at newly formed imines" if linkage!="secondary_amine" else "Reduced amine",
                qualification=scope)
        from .cofdraw_cofex import design_record
        feeds = [dict(smiles=s, label=l) for s,l in zip(smiles,(label_a,label_b))]
        structure["provenance"]["drawing_cofex"] = design_record(feeds, 'crystal', linkage, topology, network_count)
        return structure,dict(periodic_bonds=bonds,network_dimensions=net["dimensions"],short_contacts=0,
                              missing_bonds=0,extra_bonds=0,junction_sites=junctions,tetrahedral_max_angle_error_deg=round(angle,3),
                              network_count=network_count,topology=topology,junction_coordination_sequence=sequence,
                              packing=pack,cell_composition_verified=True,
                              precursor_stoichiometry=[network_count*n for n in stoichiometry])
    if tried_packing:
        raise ValueError(f"No {network_count}-network arrangement passed inter-network contact and periodic bond checks in the sampled conformers and translations. Try a different count or supply a CIF.")
    raise ValueError(f"A connected diamond net could not be fitted to this pair ({best}). The selected molecular geometry may require a different net or conformer.")
