"""Geometry-qualified molecular schematics traced through crystallographic net edges."""
from collections import deque
from functools import lru_cache
import math

import numpy as np
from rdkit import Chem

from .cofdraw_depiction import _TRIPTYCENE,compute_2d_coords
from .cofdraw_net_graphs import neighbors,quotient


# Finite incidence proposals. Dimensionality belongs to the selected parent graph.
DEFINITIONS = [
    ('dia_4_2','dia',('tetrahedral','linear'),(4,2),12,'Tetrahedral junction + ditopic spacer'),
    ('acs_6_2','acs',('prismatic','linear'),(6,2),8,'Prismatic node + ditopic spacer'),
    ('acs_6_6','acs',('prismatic','prismatic'),(6,6),4,'Two prismatic nodes'),
    ('ceq_6_3','ceq',('prismatic','triangular'),(6,3),4,'Prismatic + triangular'),
    ('stp_6_4','stp',('prismatic','quadrangular'),(6,4),4,'Prismatic + quadrangular'),
    ('hea_6_4','hea',('prismatic','tetrahedral'),(6,4),4,'Prismatic + tetrahedral'),
    ('pts_4_4','pts',('tetrahedral','quadrangular'),(4,4),4,'Tetrahedral + quadrangular'),
    ('tbo_4_3','tbo',('quadrangular','triangular'),(4,3),6,'Quadrangular + triangular'),
    ('pto_4_3','pto',('quadrangular','triangular'),(4,3),8,'Quadrangular + triangular'),
    ('fjh_4_3','fjh',('quadrangular','triangular'),(4,3),8,'Quadrangular + triangular'),
    ('ctn_4_3','ctn',('tetrahedral','triangular'),(4,3),8,'Tetrahedral + triangular'),
    ('nbo_4_2','nbo',('quadrangular','linear'),(4,2),12,'Crossed quadrangular nodes + ditopic spacer'),
]
SPATIAL={key:dict(code=code,shapes=shapes,degrees=degrees,instances=n,label=label)
         for key,code,shapes,degrees,n,label in DEFINITIONS}


def templates():
    return [dict(id=key,label=f'3D {r["code"]} · {r["label"]}',instances=r['instances'],
                 allow_pairs=[list(x) for x in dict.fromkeys((r['degrees'],r['degrees'][::-1]))],
                 shape_pair=list(r['shapes']),dimension=3,parent_net=r['code'],topology_label=None,
                 linkages=['imine','boronate_ester'] if key=='ctn_4_3' else ['imine'],
                 description=f'Unfolded molecular circuits traced from the {r["code"]} graph. The monomer pair and drawing do not establish a crystal phase.')
            for key,r in SPATIAL.items()]


_ADAMANTANE=Chem.MolFromSmiles('C1C2CC3CC1CC(C2)C3')
_ARYL_CENTER=Chem.MolFromSmarts('[C,Si;X4](-c)(-c)(-c)-c')


def site_shapes(mol):
    if len(Chem.GetMolFrags(mol)) != 1:
        return {}
    centers=mol.GetSubstructMatches(_ARYL_CENTER)
    cage=mol.GetSubstructMatch(_ADAMANTANE)
    bridgeheads={cage[a.GetIdx()] for a in _ADAMANTANE.GetAtoms() if a.GetDegree()==3} if cage else set()
    shapes={}
    from .framework_drawing import ARYL_BORONIC, catechol_pairs
    groups={role:mol.GetSubstructMatches(Chem.MolFromSmarts(pattern)) for role,pattern in
            [('aryl_aldehyde','[c][CH1](=O)'),('primary_arylamine','[c][NX3;H2;+0]'),('aryl_boronic_acid',ARYL_BORONIC)]}
    groups['aryl_vicinal_diol']=catechol_pairs(mol)
    for role,hits in groups.items():
        if not 2<=len(hits)<=6:continue
        if len(hits)==4:
            distinct_center=False
            for center,*arms in centers:
                paths=[Chem.GetShortestPath(mol,center,h[1]) for h in hits]
                if len({p[1] for p in paths})==4 and not mol.GetAtomWithIdx(center).IsInRing():
                    shapes[role]='tetrahedral';distinct_center=True;break
            if distinct_center:continue
            if bridgeheads:
                anchors=[min(bridgeheads,key=lambda b:len(Chem.GetShortestPath(mol,h[1],b))) for h in hits]
                if len(set(anchors))==4:
                    shapes[role]='tetrahedral'
                    continue
            if any(mol.GetAtomWithIdx(match[0]).IsInRing() for match in centers):
                shapes[role]='spiro'
                continue
        if len(hits)==6 and prismatic_sites(mol,hits):
            shapes[role]='prismatic'
            continue
        if len(hits)==6 and (centers or cage):
            shapes[role]='spatial_unspecified'
            continue
        # Inspect *outgoing aryl bonds*, since a pair of endpoint positions
        # alone is always collinear even for an ortho- or meta-substituted ring.
        work=Chem.Mol(mol)
        compute_2d_coords(work)
        xy=work.GetConformer().GetPositions()[:,:2]
        rays=(np.array([xy[h[0]]+xy[h[3]]-xy[h[1]]-xy[h[2]] for h in hits]) if role=='aryl_vicinal_diol'
              else np.array([xy[h[1]]-xy[h[0]] for h in hits]))
        angles=np.sort(np.arctan2(rays[:,1],rays[:,0]))
        gaps=np.degrees(np.diff(np.r_[angles,angles[0]+2*np.pi]))
        if len(hits)==2 and min(gaps)>144:shapes[role]='linear'
        if len(hits)==3 and min(gaps)>88 and max(gaps)<147:shapes[role]='triangular'
        if len(hits)==6 and min(gaps)>36 and max(gaps)<86:shapes[role]='hexagonal'
        if len(hits)==4 and min(gaps)>43 and max(gaps)<140:shapes[role]='quadrangular'
    return shapes


def prismatic_sites(mol,hits):
    core=mol.GetSubstructMatch(_TRIPTYCENE)
    if not core:return False
    blades=[set(ring) for ring in mol.GetRingInfo().AtomRings() if len(ring)==6 and set(ring)<=set(core)
            and all(mol.GetAtomWithIdx(i).GetIsAromatic() for i in ring)]
    if len(blades)!=3:return False
    attachment=[]
    for hit in hits:
        distances=[len(Chem.GetShortestPath(mol,hit[1],int(i))) for i in core]
        shortest=min(distances)
        sites={core[i] for i,n in enumerate(distances) if n==shortest}
        if len(sites)!=1:return False
        attachment.extend(sites)
    return len(set(attachment))==6 and all(len(set(attachment)&b)==2 for b in blades)


def matches(template, shapes):
    expected=template.get('shape_pair')
    return not expected or shapes==expected or shapes[::-1]==expected


def _edge(a,b):return tuple(sorted((a,b)))


def net_preview(recipe):
    """Build a finite 3D neighborhood of the chosen parent net in graph units."""
    from .cofdraw_chemistry import design_signature
    key=(recipe.get('selected_template') or {}).get('id')
    if key not in SPATIAL and key!='branched_4_3':
        return None
    if recipe.get('explicit_graph') or recipe.get('design_template_signature')!=design_signature(recipe):
        raise ValueError('Apply the modified connectivity before opening its topology in 3D.')
    rule=SPATIAL.get(key,dict(code='bor',shapes=('tetrahedral','tritopic'),degrees=(4,3)))
    components=recipe['components']
    first=components['A']
    degree_a=first['connectivity']
    if rule['degrees'][0]!=rule['degrees'][1]:
        primary='A' if degree_a==rule['degrees'][0] else 'B'
    else:
        shape=site_shapes(Chem.MolFromSmiles(first['monomer_smiles'])).get(first['role'])
        primary='A' if shape==rule['shapes'][0] else 'B'
    secondary='B' if primary=='A' else 'A'
    code=rule['code']
    adj,frac,cell=quotient(code)
    root=(0,0,0,0)
    def position(site):
        return (frac[site[0]]+np.array(site[1:]))@cell
    start_color=0
    if code=='pts':
        vectors=[position(v)-position(root) for v in neighbors(adj,root)]
        start_color=int(np.linalg.matrix_rank(vectors,tol=1e-6)==2)
    colors={root:start_color}
    frontier=[root]
    for _ in range(3):
        upcoming=[]
        for v in frontier:
            for w in neighbors(adj,v):
                if w not in colors:
                    colors[w]=1-colors[v]
                    upcoming.append(w)
        frontier=upcoming
    lookup={v:i for i,v in enumerate(colors)}
    points=[position(v) for v in colors]
    ditopic=rule['degrees'][1]==2
    groups=[]
    for v in colors:
        if ditopic:groups.append(primary)
        elif rule['degrees'][0]==rule['degrees'][1]:groups.append((primary,secondary)[colors[v]])
        else:groups.append(primary if len(adj[v[0]])==rule['degrees'][0] else secondary)
    edges=[]
    for v,i in lookup.items():
        for w in neighbors(adj,v):
            j=lookup.get(w)
            if j is None or i>=j:continue
            if ditopic:
                midpoint=len(points)
                points.append((points[i]+points[j])/2)
                groups.append(secondary)
                edges.extend([[i,midpoint],[midpoint,j]])
            else:edges.append([i,j])
    points=np.array(points)
    ij=np.array(edges)
    points*=4/np.median(np.linalg.norm(points[ij[:,0]]-points[ij[:,1]],axis=1))
    labels={k:c.get('label') or k for k,c in components.items()}
    nodes=[dict(component=g,label=labels[g],x=float(p[0]),y=float(p[1]),z=float(p[2])) for g,p in zip(groups,points)]
    return dict(kind='spatial_building_block_net',net=code,periodic=False,monomers={k:dict(label=labels[k],smiles=c['monomer_smiles']) for k,c in components.items()},
                viewer=dict(representation='building_block_net',units='graph units',atoms=nodes,bonds=edges,periodic=False),
                scope=f'Proposed {code} topology. Each vertex represents one building block. Edges show connectivity. Atom positions and a physical unit cell are undetermined.')


@lru_cache(maxsize=16)
def traced_faces(code):
    """Lift periodic edges. Trace adjacent shortest circuits and unfold across each shared bond."""
    adj,_,_=quotient(code)
    def alternate(a,b,exclude):
        if all(w in exclude or w==b for w in neighbors(adj,a)):return None
        queue=deque([b]);parent={b:None};distance={b:0}
        while queue:
            v=queue.popleft()
            if distance[v]>=10:continue
            for w in neighbors(adj,v):
                if w==a and v!=b:
                    path=[]
                    while v is not None:
                        path.append(v);v=parent[v]
                    return [a,*reversed(path)]
                if w not in parent and w not in exclude and w!=a:
                    parent[w]=v;distance[w]=distance[v]+1;queue.append(w)
    def polygon(face,start,edge,sign):
        angle=sign*2*np.pi/len(face)
        rot=np.array([[math.cos(angle),-math.sin(angle)],[math.sin(angle),math.cos(angle)]])
        pts={}
        for node in face:
            pts[node]=start
            start=start+edge
            edge=rot@edge
        return pts
    a=(0,0,0,0)
    start=alternate(a,next(neighbors(adj,a)),set())
    positions=polygon(start,np.zeros(2),np.array([1.,0]),1)
    faces=[start]
    while len(faces)<3:
        best=None
        for face in faces:
            for a,b in zip(face,face[1:]+face[:1]):
                new=alternate(a,b,set(positions)-{a,b})
                if not new:continue
                for sign in (1,-1):
                    pts=polygon(new,positions[a],positions[b]-positions[a],sign)
                    free=np.array([pts[i] for i in new[2:]])
                    old=np.array([pt for i,pt in positions.items() if i not in (a,b)])
                    gap=np.linalg.norm(free[:,None]-old,axis=2).min()
                    score=(round(min(gap,1.),5),-len(new),round(gap,5))
                    if best is None or score>best[0]:best=score,new,pts
        if best is None:break
        _,new,pts=best
        positions.update({k:v for k,v in pts.items() if k not in positions})
        faces.append(new)
    return faces,positions


def patch_graph(key,shape_a,count):
    rule=SPATIAL[key]
    faces,points=traced_faces(rule['code'])
    between=rule['degrees'][1]==2
    kinds={}
    xy={}
    adjacency=quotient(rule['code'])[0]
    root=faces[0][0]
    root_color=0
    if key=='pts_4_4':
        _,frac,cell=quotient('pts')
        rays=np.array([(frac[v]+image-frac[root[0]])@cell for v,*image in neighbors(adjacency,root)])
        root_color=int(np.linalg.matrix_rank(rays,tol=1e-6)==2)
    color={root:root_color}
    rims=[]
    for face in faces:
        anchor=next(v for v in face if v in color) if not between else face[0]
        rim=[]
        for i,v in enumerate(face):
            if between:index=0
            elif rule['degrees'][0]!=rule['degrees'][1]:index=rule['degrees'].index(len(adjacency[v[0]]))
            else:
                index=(color[anchor]+i-face.index(anchor))%2
                if v in color and color[v]!=index:raise ValueError('The proposed patch is not bipartite.')
            color[v]=index
            kinds[v]='A' if index==0 else 'B'
            xy[v]=points[v]
            rim.append(v)
            if between:
                w=face[(i+1)%len(face)]
                mid=('edge',*_edge(v,w));kinds[mid]='B';xy[mid]=(points[v]+points[w])/2
                rim.append(mid)
        rims.append(rim)
    if shape_a!=rule['shapes'][0]:
        kinds={v:('B' if k=='A' else 'A') for v,k in kinds.items()}
    sites=list(kinds)
    full_edges=sorted({tuple(sorted((sites.index(a),sites.index(b)))) for rim in rims for a,b in zip(rim,rim[1:]+rim[:1])})
    return [[sites.index(v) for v in r] for r in rims[:count]],full_edges,np.array([xy[v] for v in sites]),[kinds[v] for v in sites]



def check_projection(mol):
    """Reject severe interfragment overlaps and stretched joins in a flattened drawing."""
    xy=mol.GetConformer().GetPositions()[:,:2]
    owners=[a.GetProp('instance') for a in mol.GetAtoms()]
    # The finite neighborhoods have up to 2,400 atoms. Neighbor bins avoid
    # constructing an all-pairs matrix for a large user-defined monomer.
    bins={}
    tolerance=.5
    for atom in mol.GetAtoms():
        if not atom.GetAtomicNum():continue
        i=atom.GetIdx()
        x,y=np.floor(xy[i]/tolerance).astype(int)
        for col in range(x-1,x+2):
            for row in range(y-1,y+2):
                for j in bins.get((col,row),[]):
                    if owners[i]!=owners[j] and np.linalg.norm(xy[i]-xy[j])<tolerance:
                        raise ValueError('Separate monomers overlap in this projected circuit. Choose a different net or a smaller drawing extent.')
        bins.setdefault((x,y),[]).append(i)
    for b in mol.GetBonds():
        i,j=b.GetBeginAtomIdx(),b.GetEndAtomIdx()
        if owners[i]!=owners[j] and not 1.05<np.linalg.norm(xy[i]-xy[j])<1.95:
            raise ValueError('A linkage is excessively stretched in this circuit projection. Choose a different net or a smaller drawing extent.')
    # Internal perspective conventions belong to the depicted monomer. Edges
    # of independently placed blocks should never intersect on the page.
    bins={}
    for b in mol.GetBonds():
        i,j=b.GetBeginAtomIdx(),b.GetEndAtomIdx()
        if not b.GetBeginAtom().GetAtomicNum() or not b.GetEndAtom().GetAtomicNum():continue
        x,y=np.floor((xy[i]+xy[j])/4.).astype(int)
        u=xy[j]-xy[i]
        for cx in range(x-1,x+2):
            for cy in range(y-1,y+2):
                for p,q in bins.get((cx,cy),()):
                    if i in (p,q) or j in (p,q) or owners[i]==owners[j]==owners[p]==owners[q]:continue
                    v=xy[q]-xy[p]
                    det=u[0]*v[1]-u[1]*v[0]
                    if abs(det)<1e-7:continue
                    w=xy[p]-xy[i]
                    t=(w[0]*v[1]-w[1]*v[0])/det
                    s=(w[0]*u[1]-w[1]*u[0])/det
                    if .02<t<.98 and .02<s<.98:
                        raise ValueError('Bonds on separate building blocks cross in this projection. Choose a different net or a smaller extent.')
        bins.setdefault((x,y),[]).append((i,j))


def clear_junctions(mol):
    """Translate/rotate intact monomers to clear collisions while retaining the selected ports."""
    conformer=mol.GetConformer()
    start=conformer.GetPositions()[:,:2].copy()
    group_names=[a.GetProp('instance') for a in mol.GetAtoms()]
    _,groups=np.unique(group_names,return_inverse=True)
    centers=np.array([start[groups==k].mean(0) for k in range(groups.max()+1)])
    joins=np.array([(b.GetBeginAtomIdx(),b.GetEndAtomIdx()) for b in mol.GetBonds()
                    if groups[b.GetBeginAtomIdx()]!=groups[b.GetEndAtomIdx()]])
    heavy=np.array([a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum()])
    # Sample bond interiors as well, so moving atoms apart does not simply
    # push a neighboring aromatic ring through the middle of a bond.
    inside=np.array([(b.GetBeginAtomIdx(),b.GetEndAtomIdx()) for b in mol.GetBonds()
                     if groups[b.GetBeginAtomIdx()]==groups[b.GetEndAtomIdx()]
                     and b.GetBeginAtom().GetAtomicNum() and b.GetEndAtom().GetAtomicNum()])
    u,v=inside.T
    natoms=len(start)
    start=np.vstack((start,(2*start[u]+start[v])/3,(start[u]+2*start[v])/3))
    groups=np.r_[groups,groups[u],groups[u]]
    heavy=np.r_[heavy,np.arange(natoms,len(start))]
    local=start-centers[groups]
    joined={_edge(i,j) for i,j in joins}
    def contacts(xy):
        bins={}
        pairs=[]
        for i in heavy:
            x,y=np.floor(xy[i]/1.8).astype(int)
            for col in range(x-1,x+2):
                for row in range(y-1,y+2):
                    for j in bins.get((col,row),()):
                        if groups[i]!=groups[j] and _edge(i,j) not in joined and np.linalg.norm(xy[i]-xy[j])<1.8:
                            pairs.append((i,j))
            bins.setdefault((x,y),[]).append(i)
        return pairs
    # Penalize departure from the original pose and junction direction. The
    # clearance residual acts only on close heavy atoms of *different* blocks.
    pose=np.zeros((len(centers),3))
    xy=start.copy()
    for _ in range(30):
        collisions=contacts(xy)
        pairs=np.vstack((joins,np.array(collisions,dtype=int).reshape(-1,2)))
        i,j=pairs.T
        d=xy[i]-xy[j]
        lengths=np.linalg.norm(d,axis=1)
        units=d/np.maximum(lengths,1e-8)[:,None]
        units[lengths<1e-8]=[1.,0.]
        clash=np.arange(len(i))>=len(joins)
        active=~clash|(lengths<.9)
        weights=np.where(clash,2.,5.)*active
        n=len(pairs)
        rows=np.arange(n)
        jac=np.zeros((n+len(joins),pose.size))
        offsets=xy-centers[groups]-pose[groups,:2]
        for atoms,sign in ((i,1.),(j,-1.)):
            k=3*groups[atoms]
            jac[rows,k]=sign*units[:,0]*weights
            jac[rows,k+1]=sign*units[:,1]*weights
            jac[rows,k+2]=sign*(units[:,1]*offsets[atoms,0]-units[:,0]*offsets[atoms,1])*weights
        clearance=np.where((i>=natoms)|(j>=natoms),.59,.90)
        weights[clash&(lengths>=clearance)]=0
        jac[:n][weights==0]=0
        initial=start[joins[:,0]]-start[joins[:,1]]
        normal=np.column_stack((-initial[:,1],initial[:,0]))/np.linalg.norm(initial,axis=1)[:,None]
        tangent_weight=2.5
        for atoms,sign in ((joins[:,0],1.),(joins[:,1],-1.)):
            k=3*groups[atoms];ix=n+np.arange(len(joins))
            jac[ix,k]=sign*tangent_weight*normal[:,0]
            jac[ix,k+1]=sign*tangent_weight*normal[:,1]
            jac[ix,k+2]=sign*tangent_weight*(normal[:,1]*offsets[atoms,0]-normal[:,0]*offsets[atoms,1])
        residual=np.r_[weights*(lengths-np.where(clash,clearance,1.5)),tangent_weight*(normal*d[:len(joins)]).sum(1)]
        # Rigid-body gauge and a small restraint on rotations/translations.
        stiffness=np.tile([.12,.12,1.2],len(centers))
        step=np.linalg.solve(jac.T@jac+np.diag(stiffness**2+.02),
                             -jac.T@residual-stiffness**2*pose.ravel()).reshape(-1,3)
        pose+=step*min(1.,.45/max(1e-8,np.linalg.norm(step[:,:2],axis=1).max()),.07/max(1e-8,abs(step[:,2]).max()))
        c=np.cos(pose[groups,2]);s=np.sin(pose[groups,2])
        xy=centers[groups]+pose[groups,:2]+np.column_stack((local[:,0]*c-local[:,1]*s,local[:,0]*s+local[:,1]*c))
        if np.linalg.norm(step)<1e-5:break
    for i,(x,y) in enumerate(xy[:natoms]):
        conformer.SetAtomPosition(i,(float(x),float(y),0.))
    check_projection(mol)
