"""Local pore environments and explicitly connected planar molecular patches."""
from collections import Counter
from functools import lru_cache
from itertools import permutations, product
import math

import numpy as np
from rdkit import Chem

from .cofdraw_depiction import compute_2d_coords, normalize_imine_ports, straight_linker_axes
from .depiction_rings import symmetric_depiction_rings
from .framework_drawing import _ports


PANG_DOI = '10.1021/jacs.6b01244'
ETTA = Chem.MolToSmiles(Chem.MolFromSmiles('Nc1ccc(C(=C(c2ccc(N)cc2)c2ccc(N)cc2)c2ccc(N)cc2)cc1'))
_ARYL_FORMYL = Chem.MolFromSmarts('[c]-[CH1;X3](=[O])')
TPDA = dict(id='tpda_terphenyl', abbreviation='TPDA', full_name='p-Terphenyl-4,4″-dicarbaldehyde',
            label='TPDA · p-terphenyl-4,4″-dicarbaldehyde', smiles='O=Cc1ccc(-c2ccc(-c3ccc(C=O)cc3)cc2)cc1',
            aliases=['4,4″-Diformyl-p-terphenyl'], source_dois=[PANG_DOI],
            source_evidence=[dict(source_doi=PANG_DOI,source_label='TPDA',location='Supporting Information pp. S2–S3, SIOC-COF-2')],
            geometry_hint='linear', qualification='Reported SIOC-COF-2 dialdehyde precursor.')

MIXED_EXAMPLES = [
    dict(id='sioc-cof-1',label='SIOC-COF-1 · three pore types',monomer_ids=['tape','ta','bpda_cho'],
         linker_sequence='BCBCBC',pore_layout='kgm',doi=PANG_DOI,
         description='ETTA, terephthalaldehyde and biphenyldialdehyde. The reported arrangement has two triangular pores and an alternating hexagon.'),
    dict(id='sioc-cof-2',label='SIOC-COF-2 · extended triple pore',monomer_ids=['tape','bpda_cho','tpda_terphenyl'],
         linker_sequence='BCBCBC',pore_layout='kgm',doi=PANG_DOI,
         description='Biphenyl and terphenyl dialdehydes line the two triangles; both occur in the hexagon.'),
    dict(id='mixed-tapb',label='Hexagonal design · alternating linkers',monomer_ids=['tapb','ta','bpda_cho'],
         linker_sequence='BCBCBC',description='A three-arm amine with alternating short and long dialdehydes.'),
    dict(id='mixed-tfb',label='Hexagonal design · adjacent diamines',monomer_ids=['tfb','pa','benzidine'],
         linker_sequence='BBBCCC',description='Three adjacent short diamines and three adjacent biphenyl diamines in the selected pore.'),
    dict(id='mixed-tapp',label='Porphyrin design · opposite linkers',monomer_ids=['tapp','ta','bpda_cho'],
         linker_sequence='BCBC',description='A free-base porphyrin tetraamine with a different dialdehyde along each sheet direction.'),
    dict(id='mixed-tfppy',label='Pyrene design · adjacent linkers',monomer_ids=['tfppy','pa','benzidine'],
         linker_sequence='BBCC',description='Two short and two extended diamines placed around a four-sided pore.'),
]


def canonical_sequence(seq):
    """Chemical B/C environments equivalent by rotation or reversal of a face."""
    return min(s[i:]+s[:i] for s in (seq,seq[::-1]) for i in range(len(s)))


def mixed_pore_patterns():
    patterns = {}
    for n in (4,6):
        unique = sorted({canonical_sequence(''.join(seq)) for seq in product('BC',repeat=n) if len(set(seq))==2})
        patterns[str(n)] = [dict(sequence=s,label=f'{s} · {s.count("B")} B + {s.count("C")} C',
                                 counts=dict(B=s.count('B'),C=s.count('C'))) for s in unique]
    return patterns


@lru_cache(maxsize=512)
def linear_dialdehyde(smiles):
    mol=Chem.MolFromSmiles(smiles)
    return linear_arms(mol, mol.GetSubstructMatches(_ARYL_FORMYL))


def linear_arms(mol, matches):
    """Require both outgoing attachment directions to follow the end-to-end rod."""
    if len(matches)!=2:
        return False
    mol=Chem.Mol(mol)
    compute_2d_coords(mol,canonOrient=False,forceRDKit=True)
    xy=mol.GetConformer().GetPositions()[:,:2]
    anchors, sites = np.array([match[:2] for match in matches]).T
    return straight_linker_axes(xy[sites], xy[sites] - xy[anchors])


_HEX = ((2,0),(1,1),(-1,1),(-2,0),(-1,-1),(1,-1))
_SQUARE = ((0,0),(1,0),(1,1),(0,1))


def _add(a,b):
    return (a[0]+b[0],a[1]+b[1])


def _faces(pattern,repeats,kgm,triad=False,single_linker=None):
    """Integer vertex coordinates make neighboring faces share exact graph edges."""
    n = len(pattern)
    if kgm:
        if repeats==1 and n==3:
            i = 0 if pattern=='BBB' else 1
            return [((_HEX[i],_HEX[(i+1)%6],_add(_HEX[i],_HEX[(i+1)%6])),(single_linker*3 if single_linker else pattern))]
        faces = {}
        radius=1 if triad else max(1,repeats-1)
        for u,v in product(range(radius),repeat=2):
            center = (4*u+2*v,2*v)
            verts = tuple(_add(center,p) for p in _HEX)
            faces[frozenset(verts)] = (verts,single_linker*6 if single_linker else 'BCBCBC')
            if repeats>1:
                for i,p in enumerate(_HEX):
                    if triad and i>=2:
                        break
                    q = _HEX[(i+1)%6]
                    tri = (_add(center,p),_add(center,q),_add(center,_add(p,q)))
                    faces[frozenset(tri)] = (tri,(single_linker or 'BC'[i%2])*3)
        return list(faces.values())
    corners = _HEX if n==6 else _SQUARE
    # Neighboring faces assign the same species to their common edge. Square
    # faces alternate by half a turn so each row spans both linker lengths.
    faces=[]
    cells=((0,0),(0,1),(1,0)) if triad else product(range(repeats),repeat=2)
    for u,v in cells:
        center = (3*u,u+2*v) if n==6 else (u,v)
        parity = (u,v,u) if n==6 else (u+v,u+v)
        seq = ''.join(pattern[(i+n//2*(parity[i%(n//2)]%2))%n] for i in range(n))
        faces.append((tuple(_add(center,p) for p in corners),seq))
    return faces


def _block(component, alternate=False, *, planar_pore=False):
    mol = Chem.MolFromSmiles(component['fragment_smiles'])
    from .cofdraw_depiction import has_triptycene_core,has_free_base_porphyrin,coordgen_coords
    if has_triptycene_core(mol):
        compute_2d_coords(mol)
    elif alternate and not has_free_base_porphyrin(mol):
        coordgen_coords(mol)
    else:
        compute_2d_coords(mol,canonOrient=False,forceRDKit=True)
    ports = list(_ports(mol).values())
    endpoints = [site for _,site in ports]
    xy = mol.GetConformer().GetPositions()[:,:2].copy()
    if component.get("linkage") in ("imine", "alkyl_imine") and normalize_imine_ports(mol, xy, straight_only=not planar_pore):
        mol.SetBoolProp("_fixedImineTurns", True)
    xy -= xy[endpoints].mean(axis=0)
    tangents = xy[[dummy for dummy,_ in ports]]-xy[endpoints]
    tangents /= np.linalg.norm(tangents,axis=1)[:,None]
    return mol,ports,xy,tangents


def _turn(v,theta):
    c,s = np.cos(theta),np.sin(theta)
    return np.stack((c*v[...,0]-s*v[...,1],s*v[...,0]+c*v[...,1]),axis=-1)


def _orient(axes,ports,xy):
    target = np.arctan2(axes[:,1],axes[:,0])
    # Several junctions (e.g. boroxine) can have short attachment radii.
    # The external dummy points unambiguously identify their outgoing arms.
    radial = xy[[p[0] for p in ports]]
    phase = np.arctan2(radial[:,1],radial[:,0])
    best=None
    for indices in permutations(range(len(ports)),len(axes)):
        diff=target-phase[list(indices)]
        theta=math.atan2(np.sin(diff).sum(),np.cos(diff).sum())
        error=float(np.sum(1-np.cos(diff-theta)))
        score=(round(error,9),round(abs(theta),9),indices)
        if best is None or score<best[0]:
            best=score,indices,theta
    return best[1:]


def _fit(geometry,kinds,centers,joins,slot_pairs,orientations=None, *, images=None, cell=None):
    """Fit rigid fragments to shared junctions with analytic rotation derivatives."""
    owners=np.array(joins)
    a,b=owners.T
    n=len(kinds)
    incident=[[] for _ in kinds]
    for edge,(x,y) in enumerate(joins):
        d=centers[y]-centers[x]
        if images is not None:
            d=d+images[edge]@cell
        d=d/np.linalg.norm(d)
        incident[x].append(d)
        incident[y].append(-d)
    placed=orientations or [_orient(np.array(axes),geometry[k][1],geometry[k][2]) for k,axes in zip(kinds,incident)]
    slots=np.array([[placed[k][0][slot] for k,slot in zip(pair,sp)] for pair,sp in zip(joins,slot_pairs)])
    angles=np.array([p[1] for p in placed])
    anchors=np.array([[geometry[kinds[k]][2][geometry[kinds[k]][1][s][1]] for k,s in zip(pair,ss)]
                      for pair,ss in zip(joins,slots)])
    tangents=np.array([[geometry[kinds[k]][3][s] for k,s in zip(pair,ss)] for pair,ss in zip(joins,slots)])
    # At a terminal trigonal reaction atom the two depiction sides are
    # interchangeable before an explicit E/Z assignment. Choose the side
    # facing the opposite molecule without reflecting its other atoms.
    for i,pair in enumerate(joins):
        choices=[]
        for k,s,d in zip(pair,slots[i],tangents[i]):
            mol,ports,xy,_=geometry[kinds[k]]
            dummy,end=ports[s]
            atom=mol.GetAtomWithIdx(end)
            others=[j.GetIdx() for j in atom.GetNeighbors() if j.GetIdx()!=dummy]
            variants=[d]
            if not mol.HasProp("_fixedImineTurns") and len(others)==1 and atom.GetHybridization()==Chem.HybridizationType.SP2 and all(
                    bond.GetStereo()==Chem.BondStereo.STEREONONE for bond in atom.GetBonds()):
                axis=xy[end]-xy[others[0]]
                axis/=np.linalg.norm(axis)
                variants.append(2*np.dot(axis,d)*axis-d)
            choices.append(variants)
        axis=incident[pair[0]][slot_pairs[i][0]]
        def penalty(pair_dirs):
            da,db=_turn(np.array(pair_dirs),angles[list(pair)])
            return np.linalg.norm(da+db)**2+.6*(2-(da-db)@axis)
        tangents[i]=min(product(*choices),key=penalty)
    half_junctions=anchors+.75*tangents
    incidence=np.zeros((len(joins),n))
    ix=np.arange(len(joins))
    incidence[ix,a]=-1
    incidence[ix,b]=1
    if images is not None:
        incidence=np.column_stack((incidence,images))
    solve=np.linalg.pinv(incidence)
    project=incidence@solve-np.eye(len(joins))
    initial=angles.copy()
    def evaluate(theta):
        ends=_turn(half_junctions,theta[owners])
        dirs=_turn(tangents,theta[owners])
        target=ends[:,0]-ends[:,1]
        residual=np.r_[(project@target).ravel(),(.7*(dirs[:,0]+dirs[:,1])).ravel(),.1*(theta-initial)]
        return residual,ends,dirs,solve@target
    for _ in range(18):
        residual,ends,dirs,translations=evaluate(angles)
        dt=np.zeros((len(joins),2,n))
        dt[ix,:,a]=np.stack((-ends[:,0,1],ends[:,0,0]),axis=1)
        dt[ix,:,b]=np.stack((ends[:,1,1],-ends[:,1,0]),axis=1)
        dd=np.zeros_like(dt)
        dd[ix,:,a]=.7*np.stack((-dirs[:,0,1],dirs[:,0,0]),axis=1)
        dd[ix,:,b]=.7*np.stack((-dirs[:,1,1],dirs[:,1,0]),axis=1)
        jac=np.vstack(((project@dt.reshape(len(joins),-1)).reshape(-1,n),dd.reshape(-1,n),.1*np.eye(n)))
        step=np.linalg.solve(jac.T@jac+np.eye(n)*1e-6,-jac.T@residual)
        step*=min(1.,.2/max(1e-10,float(np.abs(step).max())))
        if np.abs(step).max()<1e-6:
            break
        for scale in (1.,.5,.25,.1):
            trial=angles+step*scale
            if np.linalg.norm(evaluate(trial)[0])<np.linalg.norm(residual):
                angles=trial
                break
        else:
            break
    translations=evaluate(angles)[3]
    ends=_turn(anchors,angles[owners])
    # Different length linkers can leave short junctions at closed loops.
    # Correct inter-block distances by translations of complete molecules.
    for _ in range(60):
        bonds=incidence@translations+ends[:,1]-ends[:,0]
        lengths=np.linalg.norm(bonds,axis=1)
        if np.min(lengths)<1e-8 or np.max(abs(lengths-1.5))<.015:
            break
        translations=solve@(ends[:,0]-ends[:,1]+1.5*bonds/lengths[:,None])
    return translations,angles,slots


def join_fragments(components,kinds,centers,joins,slot_pairs,bond_order,orientations=None,geometry=None,placement=None):
    degrees=Counter(k for pair in joins for k in pair)
    geo=geometry or {k:_block(components[k]) for k in set(kinds)}
    if sum(geo[k][0].GetNumAtoms() for k in kinds)-2*len(joins)>2400:
        raise ValueError('This drawing exceeds 2,400 atoms. Use a smaller neighborhood.')
    if any(degrees[i]>len(geo[k][1]) for i,k in enumerate(kinds)):
        raise ValueError('The selected node has too few sites for this pore arrangement.')
    translations,angles,slots=placement or _fit(geo,kinds,centers,joins,slot_pairs,orientations)
    fragments=[]
    for i,k in enumerate(kinds):
        original,ports,local,_=geo[k]
        mol=Chem.Mol(original)
        xy=_turn(local,angles[i])+translations[i]
        for atom,(x,y) in zip(mol.GetAtoms(),xy):
            atom.SetProp('component',k)
            atom.SetProp('instance',f'{k}{i+1}')
            atom.SetIntProp('fragment_atom_index',atom.GetIdx())
            mol.GetConformer().SetAtomPosition(atom.GetIdx(),(float(x),float(y),0.))
            if atom.GetAtomicNum()==0:
                atom.SetAtomMapNum(0)
        fragments.append(mol)
    for tag,(pair,pair_slots) in enumerate(zip(joins,slots),1):
        for k,s in zip(pair,pair_slots):
            dummy,endpoint=geo[kinds[k]][1][s]
            mol=fragments[k]
            if mol.GetBondBetweenAtoms(dummy,endpoint).GetBondTypeAsDouble()!=bond_order:
                raise ValueError('The reacted fragments have incompatible bond orders.')
            mol.GetAtomWithIdx(dummy).SetAtomMapNum(tag)
    joined=fragments[0]
    for mol in fragments[1:]:
        joined=Chem.CombineMols(joined,mol)
    result=Chem.molzip(joined)
    if len(Chem.GetMolFrags(result))!=1 or result.GetNumAtoms()!=joined.GetNumAtoms()-2*len(joins):
        raise ValueError('The connected molecular patch failed its boundary and valence checks.')
    result.UpdatePropertyCache(strict=True)
    symmetric_depiction_rings(result)
    origins=[dict(component=a.GetProp('component'),instance=a.GetProp('instance'),
                  fragment_atom_index=a.GetIntProp('fragment_atom_index')) for a in result.GetAtoms()]
    return result,origins


def binary_kagome_roles(components):
    """Find the node and spacer from the active, reacted precursors."""
    return tuple(next(k for k,c in components.items() if c['connectivity']==n) for n in (4,2))


def _periodic_pore_placement(geometry, pattern, kgm, node, spacer, kinds, grid, joins, metric):
    """Fit the repeating molecular motif once, then crop it to the selected faces."""
    if kgm:
        cell = np.array([[4, 0], [2, 2]])
    elif len(pattern) == 6:
        cell = np.array([[3, 1], [0, 2]])
        # Face ordering in _faces toggles opposite edges along u or v.
        cell[0] *= 2 if pattern[0] != pattern[3] or pattern[2] != pattern[5] else 1
        cell[1] *= 2 if pattern[1] != pattern[4] else 1
    else:
        cell = np.eye(2, dtype=int)
        if pattern[:2] != pattern[2:]:
            # Alternating half-turns give each row and column the same inventory
            # of short and long linkers, with one shared molecule per edge.
            cell = np.array([[1, 1], [-1, 1]])
    inverse = np.linalg.inv(cell)
    def address(kind, point):
        image = np.floor(point @ inverse + 1e-8).astype(int)
        # The integer vertex grid also contains half-integer linker midpoints.
        x, y = np.rint(2 * (point - image @ cell)).astype(int)
        return (kind, x, y), image

    # Include complete incident stars even when the requested crop is one face.
    motif_edges = set()
    for vertices, sequence in _faces(pattern, 3, kgm, single_linker=spacer):
        for a, b, kind in zip(vertices, vertices[1:] + vertices[:1], sequence):
            mid, image = address(kind, (np.array(a) + b) / 2)
            for vertex in (a, b):
                corner, offset = address(node, np.array(vertex))
                motif_edges.add((corner, mid, tuple(image - offset)))
    ordered = sorted(motif_edges)
    sites = sorted({site for a, b, _ in ordered for site in (a, b)})
    index = {site: i for i, site in enumerate(sites)}
    pairs = [(index[a], index[b]) for a, b, _ in ordered]
    sides = []
    degree = Counter()
    for a, b in pairs:
        sides.append((degree[a], degree[b]))
        degree.update((a, b))
    points = np.array([site[1:] for site in sites]) * metric / 2
    shifts = np.array([shift for _, _, shift in ordered])
    positions, angles, ports = _fit(geometry, [s[0] for s in sites], points, pairs, sides,
                                    images=shifts, cell=cell * metric)
    addresses = [address(kind, point) for kind, point in zip(kinds, grid)]
    tile = positions[-2:]
    translations = np.array([positions[index[site]] + offset @ tile for site, offset in addresses])
    lookup = dict(zip(ordered, ports))
    slots = []
    for a, b in joins:
        source, start = addresses[a]
        target, end = addresses[b]
        slots.append(lookup[source, target, tuple(end - start)])
    return translations, np.array([angles[index[site]] for site, _ in addresses]), slots


def pore_patch(components,pattern,repeats,bond_order,layout='',triad=False):
    """Build a finite molecular graph with shared blocks, labeled outer ports, and pore census."""
    kgm=layout=='kgm'
    if repeats not in (1,2,3,4) or (kgm and pattern not in ('BBB','CCC','BCBCBC','CBCBCB')) or (
            not kgm and (len(pattern) not in (4,6) or not set(pattern)<=set('BC'))):
        raise ValueError('Choose a supported pore and a drawing extent from one to four.')
    binary_kgm=kgm and set(components)=={'A','B'}
    node,spacer=binary_kagome_roles(components) if binary_kgm else ('A',None)
    faces=_faces(pattern,repeats,kgm,triad=triad,single_linker=spacer)
    edge_types={}
    for vertices,sequence in faces:
        for i,letter in enumerate(sequence):
            edge=tuple(sorted((vertices[i],vertices[(i+1)%len(vertices)])))
            if edge_types.setdefault(edge,letter)!=letter:
                raise ValueError('Conflicting linker assignments on a shared edge.')
    nodes=sorted({vertex for edge in edge_types for vertex in edge})
    ix={point:i for i,point in enumerate(nodes)}
    kinds=[node]*len(nodes)
    centers=list(map(np.array,nodes))
    joins=[]; slot_pairs=[]; degrees=Counter()
    linker_ids={}
    for edge,letter in edge_types.items():
        mid=len(kinds)
        linker_ids[edge]=f'{letter}{mid+1}'
        kinds.append(letter)
        centers.append(.5*(np.array(edge[0])+edge[1]))
        for p in edge:
            junction=ix[p]
            joins.append((junction,mid))
            slot_pairs.append((degrees[junction],degrees[mid]))
            degrees[junction]+=1; degrees[mid]+=1
    grid = np.array(centers)
    metric = np.array([.5, math.sqrt(3)/2] if kgm or len(pattern)==6 else [1., 1.])
    imines = all(c.get("linkage") in ("imine", "alkyl_imine") for c in components.values())
    geometry = {k: _block(c, planar_pore=imines) for k, c in components.items()}
    placement = (_periodic_pore_placement(geometry, pattern, kgm, node, spacer, kinds, grid, joins, metric)
                 if imines else None)
    result,origins=join_fragments(components,kinds,grid*metric,joins,slot_pairs,bond_order,
                                  geometry=geometry,placement=placement)
    census=dict(sorted(Counter(canonical_sequence(s) for _,s in faces).items()))
    pore_faces=[]
    for vertices,_ in faces:
        rim=[]
        for i,p in enumerate(vertices):
            edge=tuple(sorted((p,vertices[(i+1)%len(vertices)])))
            rim.extend((f'{node}{ix[p]+1}',linker_ids[edge]))
        pore_faces.append(rim)
    return result,dict(extent=f'{len(faces)} complete pores',template='kgm' if kgm else 'hcb' if len(pattern)==6 else 'sql',
                       atoms=result.GetNumAtoms(),blocks=len(kinds),pores=len(faces),pore_sequences=census,
                       pore_faces=pore_faces,face_labels=[canonical_sequence(s) for _,s in faces],atom_origin=origins,
                       continuation_sites=sum(a.GetAtomicNum()==0 for a in result.GetAtoms()),
                       order_assumption='Separated B and C triangles share edges with alternating hexagons.' if kgm else
                       'A two-cell B/C ordering extends the selected pore. Neighboring pores share each molecular edge.')


def pore_face_geometry(components,pattern,layout):
    """Keep one editable pore in exactly the same rigid-fragment geometry as larger patches."""
    mol,summary=pore_patch(components,pattern,1,2,layout=layout)
    return circuit_geometry(components,mol,summary,'Single face of the repeating molecular drawing')


def circuit_geometry(components,mol,summary,provenance):
    rim=summary['pore_faces'][0]
    origins=summary['atom_origin']
    xy=mol.GetConformer().GetPositions()[:,:2]
    sites={k:{endpoint:port for port,(_,endpoint) in _ports(Chem.MolFromSmiles(c['fragment_smiles'])).items()}
           for k,c in components.items()}
    incidence={}
    orders={}
    for bond in mol.GetBonds():
        a=origins[bond.GetBeginAtomIdx()]; b=origins[bond.GetEndAtomIdx()]
        if a['instance']!=b['instance']:
            for here,there in ((a,b),(b,a)):
                incidence[here['instance'],there['instance']]=sites[here['component']][here['fragment_atom_index']]
                orders[here['instance'],there['instance']]=int(bond.GetBondTypeAsDouble())
    instances=[]; connections=[]
    for i,name in enumerate(rim):
        previous,following=rim[i-1],rim[(i+1)%len(rim)]
        indices=[j for j,o in enumerate(origins) if o['instance']==name]
        incoming,outgoing=incidence[name,previous],incidence[name,following]
        instances.append(dict(id=name,component=origins[indices[0]]['component'],incoming_port=incoming,
                              outgoing_port=outgoing,layout_position=xy[indices].mean(axis=0).tolist()))
        connections.append(dict(a=dict(instance=name,port=outgoing),
                                b=dict(instance=following,port=incidence[following,name]),bond_order=orders[name,following]))
    coords=[dict(instance=o['instance'],fragment_atom_index=o['fragment_atom_index'],xy=point.tolist())
            for o,point in zip(origins,xy)]
    return dict(instances=instances,connections=connections,atom_coordinates=coords,
                coordinate_provenance=provenance)
