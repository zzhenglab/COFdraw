"""Rigid ring and tetraaryl placement for chemical projections of periodic circuits."""
import math

import numpy as np


def _ring_systems(adj):
    """Connected components of bonds in local rings, excluding framework circuits."""
    rings=[set() for _ in adj]
    for i,near in enumerate(adj):
        for j in near:
            if j<=i:continue
            visited={i}
            front=set(near)-{j}
            for _ in range(6):
                if any(j in adj[k] for k in front):
                    rings[i].add(j)
                    rings[j].add(i)
                    break
                visited.update(front)
                front={v for k in front for v in adj[k]}-visited
                if not front:break
    pending={i for i,near in enumerate(rings) if near}
    systems=[]
    while pending:
        seed=min(pending)
        system={seed}
        queue=[seed]
        while queue:
            new=rings[queue.pop()]-system
            system.update(new)
            queue.extend(new)
        pending-=system
        systems.append(sorted(system))
    return systems


def _normalize(v):
    return v/max(np.linalg.norm(v),1.e-8)


def _turn(v):
    return np.stack((-v[...,1],v[...,0]),axis=-1)


def nearby_pairs(points, owners, cutoff):
    """Enumerate spatially local contact pairs without an N² allocation."""
    buckets={}
    left=[];right=[]
    for i,(x,y) in enumerate(np.floor(points/cutoff).astype(int)):
        near=[j for dx in (-1,0,1) for dy in (-1,0,1)
              for j in buckets.get((x+dx,y+dy),()) if owners[i]!=owners[j]]
        if near:
            candidates=np.asarray(near)
            d=points[candidates]-points[i]
            candidates=candidates[(d*d).sum(1)<cutoff**2]
            left.extend(candidates);right.extend([i]*len(candidates))
        buckets.setdefault((x,y),[]).append(i)
    return np.array(left,dtype=int),np.array(right,dtype=int)


def crossing_gaps(mol,points,depth):
    """Intervals to leave open on the rear stroke at genuine projected crossings."""
    edges=np.array([(b.GetBeginAtomIdx(),b.GetEndAtomIdx()) for b in mol.GetBonds()],dtype=int)
    if len(edges)<2:return {}
    starts=points[edges[:,0]]
    delta=points[edges[:,1]]-starts
    lengths=np.linalg.norm(delta,axis=1)
    indices=np.arange(len(edges))
    a,b=nearby_pairs(starts+.5*delta,indices,max(1.,lengths.max()))
    separate=(edges[a,0]!=edges[b,0])&(edges[a,0]!=edges[b,1])&(edges[a,1]!=edges[b,0])&(edges[a,1]!=edges[b,1])
    ordinary=np.array([bond.GetBeginAtom().GetAtomicNum()>0 and bond.GetEndAtom().GetAtomicNum()>0
                       and not bond.HasProp("junctionDepth") for bond in mol.GetBonds()])
    separate&=ordinary[a]&ordinary[b]
    a,b=a[separate],b[separate]
    def cross(u,v):
        return u[:,0]*v[:,1]-u[:,1]*v[:,0]
    offset=starts[b]-starts[a]
    denom=cross(delta[a],delta[b])
    safe=np.where(np.abs(denom)>1.e-8,denom,np.inf)
    t=cross(offset,delta[b])/safe
    u=cross(offset,delta[a])/safe
    inside=(t>.04)&(t<.96)&(u>.04)&(u<.96)
    result={}
    for left,right,tl,tr,area in zip(a[inside],b[inside],t[inside],u[inside],np.abs(denom[inside])):
        za,zb=depth[edges[left]],depth[edges[right]]
        if za[0]+tl*(za[1]-za[0])<=zb[0]+tr*(zb[1]-zb[0]):
            rear,front,where=left,right,tl
        else:
            rear,front,where=right,left,tr
        half=min(.25,4.*lengths[front]/area)
        result.setdefault(rear,[]).append((max(0.,where-half),min(1.,where+half)))
    return result


def _remote_atoms(mol):
    """Sample pairs on separate parts of the graph for projection scoring."""
    n=mol.GetNumAtoms()
    adj=[{a.GetIdx() for a in atom.GetNeighbors()} for atom in mol.GetAtoms()]
    sampled=np.unique(np.linspace(0,n-1,min(n,240),dtype=int))
    pairs=[]
    for k,i in enumerate(sampled):
        seen={i};front={i}
        for _ in range(4):
            front={v for j in front for v in adj[j]}-seen
            seen.update(front)
        pairs.extend((i,j) for j in sampled[:k] if j not in seen)
    return np.asarray(pairs,dtype=int).reshape(-1,2)


def reference_layout(mol,xyz,fragment_layout,tilt_degrees=0):
    """Choose a readable view using the full molecular graph around the net."""
    xyz=xyz-xyz.mean(0)
    pairs=_remote_atoms(mol)
    if not len(pairs):
        return fit_ring_bodies(mol,dict(enumerate(xyz[:,:2])),fragment_layout),xyz[:,2]
    k=np.arange(256)+.5
    z=1-2*k/256
    r=np.sqrt(1-z*z)
    cameras=np.vstack(([0,0,1],np.column_stack((r*np.cos(k*2.399963),r*np.sin(k*2.399963),z))))
    scores=np.zeros(len(cameras))
    for chunk in np.array_split(pairs,max(1,(len(pairs)+1023)//1024)):
        delta=xyz[chunk[:,0]]-xyz[chunk[:,1]]
        depth=sum(delta[:,k,None]*cameras[:,k] for k in range(3))
        distance=np.sqrt(np.maximum(0.,(delta*delta).sum(1)[:,None]-depth**2))
        scores+=(np.maximum(2.9-distance,0.)**2).sum(0)
    choices=[]
    for idx in np.argsort(scores):
        if all(abs(cameras[idx]@cameras[other])<.985 for other in choices):choices.append(idx)
        if len(choices)==3:break
    best=None
    ends=np.array([(b.GetBeginAtomIdx(),b.GetEndAtomIdx()) for b in mol.GetBonds()])
    for idx in choices:
        normal=cameras[idx] if cameras[idx,2]>=0 else -cameras[idx]
        ex=_normalize(np.cross(normal,np.eye(3)[abs(normal).argmin()]))
        ey=np.cross(normal,ex)
        depth=xyz@normal
        flat=xyz@np.array([ex,ey]).T
        pos=fit_ring_bodies(mol,dict(enumerate(1.15*flat)),fragment_layout)
        points=np.array([pos[i] for i in range(len(xyz))])
        gap=np.linalg.norm(points[pairs[:,0]]-points[pairs[:,1]],axis=1)
        lengths=np.linalg.norm(points[ends[:,0]]-points[ends[:,1]],axis=1)
        crossings=sum(map(len,crossing_gaps(mol,points*16,depth).values()))
        score=(np.maximum(2.5-gap,0)**2).sum()+12*(np.maximum(abs(lengths-1.5)-.22,0)**2).sum()
        score+=18*crossings+90*max(0,1.85-gap.min())**2
        if best is None or score<best[0]:best=(score,pos,depth,normal,ex)
        if not crossings and gap.min()>1.8 and max(lengths)<2.15:break
    if not tilt_degrees:return best[1:3]
    normal,ex=best[3:]
    angle=math.radians(tilt_degrees)
    normal=math.cos(angle)*normal+math.sin(angle)*np.cross(normal,ex)
    flat=xyz@np.array([ex,np.cross(normal,ex)]).T
    return fit_ring_bodies(mol,dict(enumerate(1.15*flat)),fragment_layout),xyz@normal


def fit_ring_bodies(mol,positions,fragment_layout):
    """Keep cyclic residues rigid and reserve one drawing quadrant per aryl arm.

    Distance along the projection sets the initial ordering of different nodes.
    Only drawing coordinates are fitted. No atom or periodic bond is changed.
    """
    n=mol.GetNumAtoms()
    old=np.array([positions[i] for i in range(n)])
    adj=[[a.GetIdx() for a in atom.GetNeighbors()] for atom in mol.GetAtoms()]
    systems=_ring_systems(adj)
    owner={a:k for k,group in enumerate(systems) for a in group}
    ring_ids=set(owner)
    xy=old.copy()
    for group in systems:
        coords=fragment_layout(mol,group)
        local=np.array([coords[i] for i in group])
        local-=local.mean(0)
        points=old[group]
        u,_,vt=np.linalg.svd(local.T@(points-points.mean(0)))
        xy[group]=local@(u@vt)+points.mean(0)

    tetra=[]
    for i,near in enumerate(adj):
        if i in owner or len(near)!=4 or mol.GetAtomWithIdx(i).GetAtomicNum() not in (6,14):
            continue
        if not all(j in owner and mol.GetAtomWithIdx(j).GetAtomicNum()==6 for j in near):
            continue
        groups=[owner[j] for j in near]
        if len(set(groups))!=4 or any(len(systems[k])!=6 for k in groups):
            continue
        # Preserve the cyclic ordering from the viewing direction. A fixed
        # depiction atom order can swap two ports and cross a complete linker.
        angles=[math.atan2(*(old[systems[k]].mean(0)-old[i])[::-1]) for k in groups]
        order=np.argsort(angles)
        phases=np.array(angles)[order]-np.arange(4)*np.pi/2
        phase=np.angle(np.exp(1j*phases).sum())
        tetra.append(i)
        merged=[i]
        for rank,index in enumerate(order):
            group=systems[groups[index]]
            first=near[index]
            radial=_normalize(xy[group].mean(0)-xy[first])
            angle=phase+rank*np.pi/2
            target=np.array([math.cos(angle),math.sin(angle)])
            rotation=np.column_stack((target,_turn(target)))@np.array([radial,_turn(radial)])
            xy[group]=(xy[group]-xy[first])@rotation.T+xy[i]+1.5*target
            merged+=group
            systems[groups[index]]=[]
        systems.append(merged)
        owner.update((j,len(systems)-1) for j in merged)
    bodies=[group for group in systems if group]
    # Keep explicitly drawn acyclic bonds (e.g. an imine) together as well.
    remaining=set(range(n))-set(owner)
    while remaining:
        seed=min(remaining)
        queue=[seed];group={seed}
        while queue:
            new=(set(adj[queue.pop()])&remaining)-group
            queue.extend(new);group.update(new)
        remaining-=group
        bodies.append(sorted(group))
    ids=np.empty(n,dtype=int)
    centers=np.array([xy[g].mean(0) for g in bodies])
    for k,g in enumerate(bodies):ids[g]=k
    local=xy-centers[ids]
    bonds=[(i,j) for i,row in enumerate(adj) for j in row if i<j and ids[i]!=ids[j]]
    if not bonds:return dict(enumerate(xy))
    ends=np.array(bonds)
    sites=ids[ends]
    def direction(i,j):
        inside=[k for k in adj[i] if k in ring_ids and ids[k]==ids[i] and k!=j]
        if i in ring_ids and len(inside)>=2:
            return _normalize(sum((_normalize(xy[i]-xy[k]) for k in inside),start=np.zeros(2)))
        if len(adj[i])==2 and len(siblings:=[k for k in adj[i] if k!=j and ids[k]==ids[i]])==1:
            axis=_normalize(xy[siblings[0]]-xy[i])
            sign=1 if _turn(axis)@(old[j]-old[i])>=0 else -1
            return -.5*axis+sign*np.sqrt(.75)*_turn(axis)
        return _normalize(old[j]-old[i])
    rays=np.array([[direction(i,j),direction(j,i)] for i,j in bonds])
    angles=np.zeros(len(bodies))
    translations=centers.copy()
    count=len(bodies)
    # Relax contacts between separate rigid residues in the actual atom
    # drawing. Pure center-to-center separation misses a methoxy arm or a
    # boronate ring that cuts into an otherwise well separated neighbor.
    ring_sets={tuple(sorted(map(int,b.GetProp("drawingRing").split(",")))) for b in mol.GetBonds()
               if b.HasProp("drawingRing")}
    midpoints=[[i,j] for i,row in enumerate(adj) for j in row if i<j and ids[i]==ids[j]]
    markers=[[i] for i in range(n)]+[list(group) for group in sorted(ring_sets)]+midpoints
    marker_ids=np.array([ids[group[0]] for group in markers])
    offsets=np.array([local[g].mean(0) for g in markers])
    radii=np.r_[np.full(n,.86),[1.95 if len(g)==6 else 1.65 for g in sorted(ring_sets)],np.full(len(midpoints),.72)]
    initial=translations[marker_ids]+offsets*np.cos(angles[marker_ids,None])+_turn(offsets)*np.sin(angles[marker_ids,None])
    aa,bb=nearby_pairs(initial,marker_ids,11.)
    if tetra:
        cc,dd=np.triu_indices(len(tetra),1)
        aa=np.r_[aa,np.array(tetra)[cc]];bb=np.r_[bb,np.array(tetra)[dd]]
    clearance=radii[aa]+radii[bb]
    # Give fragments on unrelated branches a visible margin. Distinct aryl
    # edges or heteroatom labels spaced by one bond length can look bonded.
    attached={tuple(sorted(row)) for row in sites}
    unrelated=np.array([tuple(sorted((a,b))) not in attached
                        for a,b in zip(marker_ids[aa],marker_ids[bb])])
    clearance+=1.35*unrelated
    if len(tetra)>1:clearance[-len(cc):]=12.9
    target_translations=translations.copy()
    target_angles=angles.copy()
    rigid_weight=.06
    for iteration in range(44):
        angles_at=angles[marker_ids,None]
        lever=offsets*np.cos(angles_at)+_turn(offsets)*np.sin(angles_at)
        pos=translations[marker_ids]+lever
        diff=pos[aa]-pos[bb]
        dist=np.linalg.norm(diff,axis=1)
        active=dist<clearance
        contact_a,contact_b=aa[active],bb[active]
        # Distances protect individual ring perimeters and also their interiors.
        a=np.r_[ends[:,0],contact_a];b=np.r_[ends[:,1],contact_b]
        delta=pos[a]-pos[b];d=np.maximum(np.linalg.norm(delta,axis=1),1.e-7)
        v=delta/d[:,None]
        weights=np.r_[np.full(len(ends),3.8),np.full(len(contact_a),1.6)]
        v*=weights[:,None]
        res=(d-np.r_[np.full(len(ends),1.5),clearance[active]])*weights
        ja=np.column_stack([v,(_turn(lever[a])*v).sum(1)])
        jb=-np.column_stack([v,(_turn(lever[b])*v).sum(1)])
        # Exocyclic aryl bonds should point along the external bisector.
        angles_at=angles[sites][...,None]
        rot_rays=rays*np.cos(angles_at)+_turn(rays)*np.sin(angles_at)
        tips=lever[ends]+.75*rot_rays
        link=(pos[ends[:,0]]-pos[ends[:,1]]+.75*(rot_rays[:,0]-rot_rays[:,1]))
        def bond_jac(k):
            base=np.tile(np.eye(2),(len(ends),1))
            return np.column_stack((base,_turn(tips[:,k]).ravel()))
        ja=np.vstack((ja,.65*bond_jac(0)))
        jb=np.vstack((jb,-.65*bond_jac(1)))
        first=np.r_[marker_ids[a],np.repeat(sites[:,0],2)]
        second=np.r_[marker_ids[b],np.repeat(sites[:,1],2)]
        residual=np.r_[res,.65*link.ravel()]
        def scatter(at,values):
            return np.column_stack([np.bincount(at,col,minlength=count) for col in values.T])
        def jt(values):
            return scatter(first,ja*values[:,None])+scatter(second,jb*values[:,None])
        prior=np.column_stack((translations-target_translations,angles-target_angles))
        gradient=jt(residual)+rigid_weight**2*prior
        # Each constraint touches just two rigid bodies. A matrix-free
        # normal-equation solve avoids allocating a dense molecular Hessian.
        damp=.12+rigid_weight**2
        diag=scatter(first,ja**2)+scatter(second,jb**2)+damp
        remain=-gradient
        direction=remain/diag
        norm=np.sum(remain*direction)
        step=np.zeros_like(prior)
        for _ in range(24):
            if norm<1.e-10:break
            row=(ja*direction[first]+jb*direction[second]).sum(1)
            hess=jt(row)+damp*direction
            alpha=norm/max(np.sum(direction*hess),1.e-12)
            step+=alpha*direction
            remain=remain-alpha*hess
            next_norm=np.sum(remain**2/diag)
            direction=remain/diag+(next_norm/norm)*direction
            norm=next_norm
        step*=min(1.,.16/max(np.abs(step[:,2]).max(),1.e-9))
        translations+=step[:,:2];angles+=step[:,2]
        if np.abs(step).max()<.003:break
    cos=np.cos(angles[ids])[:,None];sin=np.sin(angles[ids])[:,None]
    return dict(enumerate(translations[ids]+local*cos+_turn(local)*sin))
