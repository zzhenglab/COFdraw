"""Lay out a traced periodic circuit as a synthetic chemist's finite bond-line scheme."""
import math
from itertools import product
from xml.etree import ElementTree as ET

import numpy as np
from rdkit import Chem
from rdkit.Chem import rdDepictor


def fragment_layout(mol, indices):
    indices=sorted(indices)
    local=Chem.RWMol()
    remap={old:local.AddAtom(Chem.Atom(mol.GetAtomWithIdx(old))) for old in indices}
    for i in indices:
        for b in mol.GetAtomWithIdx(i).GetBonds():
            j=b.GetOtherAtomIdx(i)
            if i<j and j in remap:
                local.AddBond(remap[i],remap[j],b.GetBondType())
    frag=local.GetMol()
    frag.UpdatePropertyCache(False)
    Chem.GetSymmSSSR(frag)
    rdDepictor.Compute2DCoords(frag,canonOrient=False,useRingTemplates=True)
    # Ring perception on the entire macrocycle can enumerate an enormous set
    # of equivalent paths through adjacent aromatic rings. The fragment contains
    # the complete aryl systems and suffices for bond placement.
    Chem.KekulizeIfPossible(frag,clearAromaticFlags=True)
    for ring in frag.GetRingInfo().AtomRings():
        if len(ring)>8:
            continue
        original=[indices[j] for j in ring]
        for i,j in zip(original,original[1:]+original[:1]):
            mol.GetBondBetweenAtoms(i,j).SetProp("drawingRing",",".join(map(str,original)))
    for b in frag.GetBonds():
        i,j=indices[b.GetBeginAtomIdx()],indices[b.GetEndAtomIdx()]
        mol.GetBondBetweenAtoms(i,j).SetDoubleProp("drawingOrder",b.GetBondTypeAsDouble())
    pos=frag.GetConformer().GetPositions()[:,:2]
    return {old:pos[k] for old,k in remap.items()}


def rotate_to(points, start, source, target, origin):
    a=math.atan2(source[1],source[0])
    b=math.atan2(target[1],target[0])-a
    rot=np.array([[math.cos(b),-math.sin(b)],[math.sin(b),math.cos(b)]])
    return {i:origin+rot@(p-start) for i,p in points.items()}


def segment_vertices(lengths):
    # Place complete chemically drawn arms around a convex cyclic polygon.
    # Its circumradius solves the measured chord lengths without stretching bonds.
    lo=max(lengths)/2
    hi=sum(lengths)
    def angles(radius):
        return 2*np.arcsin(np.clip(np.asarray(lengths)/(2*radius),0,1))
    if sum(angles(lo))<2*math.pi-1e-7:
        raise ValueError("One illustrated linker is too long to close this schematic circuit")
    for _ in range(55):
        radius=(lo+hi)/2
        if sum(angles(radius))>2*math.pi:
            lo=radius
        else:
            hi=radius
    t=np.r_[0,np.cumsum(angles(hi)[:-1])]+math.pi/2-math.pi/len(lengths)
    polygon=hi*np.column_stack((np.cos(t),np.sin(t)))
    # Draw the projected circuit slightly obliquely, like a framework scheme.
    # Solve for edge directions while keeping every chemical fragment rigid.
    projected=polygon@np.array([[1.,.09],[.13,.82]])
    steps=np.roll(projected,-1,axis=0)-projected
    theta=np.arctan2(steps[:,1],steps[:,0])
    lengths=np.array(lengths)
    for _ in range(12):
        displacements=np.column_stack((np.cos(theta),np.sin(theta)))*lengths[:,None]
        error=displacements.sum(0)
        if np.linalg.norm(error)<1e-8:
            vertices=np.vstack(([0.,0.],displacements.cumsum(0)[:-1]))
            return vertices-vertices.mean(0)
        jac=np.column_stack((-displacements[:,1],displacements[:,0])).T
        theta-=np.linalg.lstsq(jac,error,rcond=None)[0]
    return polygon


def draw_scheme(mol, backbone, segments, depths=None, blocks=None, rings=None, projected=None):
    mol=Chem.Mol(mol)
    layouts={}
    centers={}
    ring_centers=[]
    for loop,fragments in rings or [(backbone,segments)]:
        lengths=[]
        for a,b,fragment in zip(loop,loop[1:]+loop[:1],fragments):
            key=tuple(sorted((a,b)))
            if key not in layouts:
                layouts[key]=fragment_layout(mol,set(fragment)|{a,b})
            coords=layouts[key]
            lengths.append(np.linalg.norm(coords[b]-coords[a]))
        if projected is not None:
            centers.update((i,projected[i]) for i in loop)
            continue
        polygon=dict(zip(loop,segment_vertices(lengths)))
        if centers:
            shared=[a for a in loop if a in centers]
            a,b=max(((i,j) for i in shared for j in shared),key=lambda pair:np.linalg.norm(polygon[pair[0]]-polygon[pair[1]]))
            proposals=[]
            for reflection in (1,-1):
                source={i:xy*[1,reflection] for i,xy in polygon.items()}
                trial=rotate_to(source,source[a],source[b]-source[a],centers[b]-centers[a],centers[a])
                midpoint=np.mean(list(trial.values()),axis=0)
                # Shared edges have two sides. Put the new pore on the unoccupied side.
                mismatch=sum(np.linalg.norm(trial[i]-centers[i])**2 for i in shared)
                clearance=min(np.linalg.norm(midpoint-p) for p in ring_centers)
                proposals.append((mismatch-clearance*.001,trial))
            polygon=min(proposals,key=lambda p:p[0])[1]
        ring_centers.append(np.mean(list(polygon.values()),axis=0))
        for i,p in polygon.items():centers.setdefault(i,p)
    # Reconcile all shared junctions before rigidly posing their molecular linkers.
    # Length constraints use lengths of the chemically drawn fragments themselves.
    if projected is not None:
        # The camera locates the front and rear chairs, but a nearly end-on
        # arm can shrink a phenyl group to a dot. Fit schematic junctions to
        # natural linker spans with a soft restraint on their projected positions.
        keys=sorted(centers)
        indices={key:i for i,key in enumerate(keys)}
        ends=np.array([[indices[a],indices[b]] for a,b in layouts])
        spans=np.array([np.linalg.norm(p[b]-p[a]) for (a,b),p in layouts.items()])
        original=np.array([centers[key] for key in keys])
        original-=original.mean(axis=0)
        observed=np.linalg.norm(original[ends[:,0]]-original[ends[:,1]],axis=1)
        original*=np.median(spans/np.maximum(observed,.01))
        xy=original.copy()
        tether=.33
        identity=np.eye(xy.size)
        rows=np.arange(len(ends))
        turns=[]
        for center in range(len(keys)):
            adjacent=[b if a==center else a for a,b in ends if center in (a,b)]
            turns.extend((a,center,b) for i,a in enumerate(adjacent) for b in adjacent[:i])
        left,mid,right=np.array(turns).T
        angle_weight=.8*np.median(spans)
        # Two junctions from different chair faces can project to almost the
        # same point even when all incident arm lengths and angles look good.
        pairs=[(i,j) for i in range(len(keys)) for j in range(i)
               if tuple(sorted((keys[i],keys[j]))) not in layouts]
        clearance=min(15.,.78*np.median(spans))
        for _ in range(36):
            delta=xy[ends[:,0]]-xy[ends[:,1]]
            lengths=np.linalg.norm(delta,axis=1)
            directions=delta/np.maximum(lengths[:,None],1e-8)
            jac=np.zeros((len(ends),len(keys),2))
            jac[rows,ends[:,0]]=directions
            jac[rows,ends[:,1]]=-directions
            jac=jac.reshape(len(ends),-1)
            # Leave room for the aromatic rings on the arms around a junction.
            p,q=xy[left]-xy[mid],xy[right]-xy[mid]
            lp=np.maximum(np.linalg.norm(p,axis=1)[:,None],1e-8)
            lq=np.maximum(np.linalg.norm(q,axis=1)[:,None],1e-8)
            u,v=p/lp,q/lq
            dot=(u*v).sum(1)
            active=dot>math.cos(math.radians(78))
            angle=np.maximum(0,dot-math.cos(math.radians(78)))*angle_weight
            dl=angle_weight*(v-dot[:,None]*u)/lp*active[:,None]
            dr=angle_weight*(u-dot[:,None]*v)/lq*active[:,None]
            aj=np.zeros((len(turns),len(keys),2))
            turn_ids=np.arange(len(turns))
            aj[turn_ids,left]=dl
            aj[turn_ids,right]=dr
            aj[turn_ids,mid]=-dl-dr
            jac=np.vstack([jac,aj.reshape(len(turns),-1)])
            residual=np.r_[lengths-spans,angle]
            for i,j in pairs:
                ray=xy[i]-xy[j]
                distance=np.linalg.norm(ray)
                if distance>=clearance:continue
                direction=ray/distance if distance>1e-8 else np.array([1.,0.])
                row=np.zeros_like(xy)
                row[i],row[j]=2*direction,-2*direction
                jac=np.vstack((jac,row.ravel()))
                residual=np.append(residual,2*(distance-clearance))
            gradient=jac.T@residual+tether**2*(xy-original).ravel()
            step=np.linalg.solve(jac.T@jac+(tether**2+.1)*identity,gradient).reshape(-1,2)
            xy-=step
            if np.max(np.abs(step))<1e-5:break
        centers=dict(zip(keys,xy))
    error=0.
    for _ in range(0 if projected is not None else 150):
        error=0.
        for (a,b),coords in layouts.items():
            d=centers[b]-centers[a]
            length=np.linalg.norm(d)
            target=np.linalg.norm(coords[b]-coords[a])
            error=max(error,abs(length-target))
            correction=d*(length-target)/max(2*length,1e-8)
            centers[a]=centers[a]+correction
            centers[b]=centers[b]-correction
        if error<1e-6:break
    if error>.02:raise ValueError("Adjacent molecular circuits could not be flattened without stretching bonds")
    positions={}
    poses={}
    for (a,b),coords in layouts.items():
        if projected is not None:
            # Foreshorten an entire molecular arm uniformly in the oblique
            # view. Every aromatic ring retains its internal bond angles.
            ratio=np.linalg.norm(centers[b]-centers[a])/np.linalg.norm(coords[b]-coords[a])
            coords={i:p*ratio for i,p in coords.items()}
        layout=rotate_to(coords,coords[a],coords[b]-coords[a],centers[b]-centers[a],centers[a])
        if projected is not None:
            axis=centers[b]-centers[a]
            axis/=max(np.linalg.norm(axis),1e-8)
            mirror=2*np.outer(axis,axis)-np.eye(2)
            reverse={i:centers[a]+mirror@(p-centers[a]) for i,p in layout.items()}
            poses[a,b]=[layout,reverse]
        positions.update(layout)
    if poses:
        # Either planar orientation of a conventionally drawn arm can meet the
        # same junctions. Choose orientations together so adjacent aryls clear.
        clouds={edge:[np.array([p for i,p in pose.items() if i not in edge]) for pose in pair]
                for edge,pair in poses.items()}
        costs={}
        for n,e in enumerate(poses):
            for f in list(poses)[:n]:
                score=np.zeros((2,2))
                for a,p in enumerate(clouds[e]):
                    for b,q in enumerate(clouds[f]):
                        dist=np.linalg.norm(p[:,None,:]-q[None,:,:],axis=2)
                        score[a,b]=np.maximum(2.15-dist,0).sum()
                costs[e,f]=score
                costs[f,e]=score.T
        flips={e:0 for e in poses}
        for _ in range(6):
            changed=False
            for e in poses:
                energies=sum((costs[e,f][:,k] for f,k in flips.items() if f!=e),np.zeros(2))
                if energies[1-flips[e]]<energies[flips[e]]-1e-6:
                    flips[e]=1-flips[e]
                    changed=True
            if not changed:break
        for edge,k in flips.items():positions.update(poses[edge][k])
    # Keep both outgoing branches on the exterior. Each fragment has its own
    # conventional bond angles and ring coordinates and receives only a rotation.
    outside=set(range(mol.GetNumAtoms()))-set(positions)
    for center,vertex in centers.items():
        attachments=sorted(a.GetIdx() for a in mol.GetAtomWithIdx(center).GetNeighbors() if a.GetIdx() in outside)
        if not attachments:continue
        choices=[]
        fixed=[positions[a.GetIdx()]-vertex for a in mol.GetAtomWithIdx(center).GetNeighbors() if a.GetIdx() in positions]
        occupied=[d/np.linalg.norm(d) for d in fixed]
        angles=sorted(math.atan2(d[1],d[0]) for d in occupied)
        gap,angle_start=max((end-start,start) for start,end in zip(angles,angles[1:]+[angles[0]+2*math.pi]))
        radial=angle_start+gap/2
        for rank,first in enumerate(attachments):
            branch={first}
            queue=[first]
            while queue:
                i=queue.pop()
                for a in mol.GetAtomWithIdx(i).GetNeighbors():
                    j=a.GetIdx()
                    if j in outside and j not in branch:
                        branch.add(j)
                        queue.append(j)
            outside-=branch
            coords=fragment_layout(mol,branch|{center})
            angle=radial+math.radians((rank-(len(attachments)-1)/2)*84)
            nearby=np.array([p for i,p in positions.items() if i!=center and np.linalg.norm(p-vertex)<12])
            # A kink at the end of an imine means its aryl ipso bond may deviate
            # from the direction between node centers. Fit the two outward arms
            # against the actual drawing so no separate phenyl rings look fused.
            candidates=[]
            for offset in range(-40,41,4):
                theta=angle+math.radians(offset)
                target=np.array([math.cos(theta),math.sin(theta)])
                trial=rotate_to(coords,coords[center],coords[first]-coords[center],target,vertex)
                points=np.array([p for i,p in trial.items() if i!=center])
                distances=np.linalg.norm(points[:,None,:]-nearby[None,:,:],axis=2)
                penalty=(np.maximum(2.05-distances,0)**2).sum()+.00009*offset**2
                for direction in occupied:
                    separation=math.degrees(math.acos(np.clip(direction@target,-1,1)))
                    penalty+=max(0,73-separation)**2
                candidates.append((penalty,trial,target,points))
            choices.append(candidates)
        def clearance(combination):
            score=sum(c[0] for c in combination)
            for i,(_,_,v,p) in enumerate(combination):
                for _,_,w,q in combination[:i]:
                    separation=math.degrees(math.acos(np.clip(v@w,-1,1)))
                    score+=max(0,73-separation)**2
                    score+=(np.maximum(2.05-np.linalg.norm(p[:,None,:]-q[None,:,:],axis=2),0)**2).sum()
            return score
        for _,trial,_,_ in min(product(*choices),key=clearance):
            positions.update(trial)
    if len(positions)!=mol.GetNumAtoms():
        # Continuation marks clipped on a substituent along one of the linkers.
        for i in sorted(outside):
            atom=mol.GetAtomWithIdx(i)
            if atom.GetAtomicNum()!=0:
                raise ValueError("Unplaced heavy-atom fragment")
            anchor=atom.GetNeighbors()[0].GetIdx()
            other=[positions[a.GetIdx()]-positions[anchor] for a in mol.GetAtomWithIdx(anchor).GetNeighbors() if a.GetIdx() in positions]
            direction=-np.sum(other,axis=0) if other else np.array([0.,1.])
            positions[i]=positions[anchor]+1.45*direction/max(np.linalg.norm(direction),1e-6)
    return vector_art(mol,positions,depths,blocks)


def junction_scheme(mol, center, blocks, spread=0):
    """Draw complete arms at the same bond scale with a front/back convention."""
    mol=Chem.Mol(mol)
    arms=sorted(a.GetIdx() for a in mol.GetAtomWithIdx(center).GetNeighbors())
    positions={center:np.zeros(2)}
    depths=np.full(mol.GetNumAtoms(),.7)
    for k,first in enumerate(arms):
        members={center,first}
        pending=[first]
        while pending:
            for other in mol.GetAtomWithIdx(pending.pop()).GetNeighbors():
                j=other.GetIdx()
                if j not in members:
                    members.add(j)
                    pending.append(j)
        if len(members.intersection(arms))!=1:
            raise ValueError("Junction arms reconnect within this fragment")
        coords=fragment_layout(mol,members)
        angle=math.radians(12+90*k+(-1)**k*spread)
        target=np.array([math.cos(angle),math.sin(angle)])
        positions.update(rotate_to(coords,coords[center],coords[first]-coords[center],target,np.zeros(2)))
        for j in members-{center}:
            depths[j]=(.94,.72,.38,.72)[k]
        if k in (0,2):
            bond=mol.GetBondBetweenAtoms(center,first)
            bond.SetProp("junctionDepth","front" if k==0 else "rear")
            bond.SetIntProp("junctionCenter",center)
    return vector_art(mol,positions,depths,blocks)


def reference_scheme(mol, xyz, blocks, rings):
    """Keep fused aromatic rings readable in a projection of a reference CIF."""
    mol=Chem.Mol(mol)
    xy=xyz[:,:2].copy()
    groups=[]
    for ring in rings:
        if not all(mol.GetAtomWithIdx(i).GetIsAromatic() for i in ring):continue
        merged=set(ring)
        for group in groups.copy():
            if group & merged:
                merged.update(group)
                groups.remove(group)
        groups.append(merged)
    constraints={}
    for group in groups:
        ids=sorted(group)
        layout=fragment_layout(mol,ids)
        local=np.array([layout[i] for i in ids])
        local-=local.mean(0)
        points=xy[ids]
        u,_,vt=np.linalg.svd(local.T@(points-points.mean(0)))
        xy[ids]=local@(u@vt)+points.mean(0)
        for i in ids:
            near=[a.GetIdx() for a in mol.GetAtomWithIdx(i).GetNeighbors() if a.GetIdx() in group]
            for j in near:
                for k in [i,*near]:
                    if k!=j:constraints[tuple(sorted((k,j)))]=np.linalg.norm(layout[k]-layout[j])
    for bond in mol.GetBonds():
        a,b=sorted((bond.GetBeginAtomIdx(),bond.GetEndAtomIdx()))
        constraints.setdefault((a,b),np.linalg.norm(xyz[a]-xyz[b]))
    ends=np.array(list(constraints));lengths=np.array(list(constraints.values()))
    # Simultaneous distance relaxation prevents the locally opened rings
    # from leaving a long bond at their attachment to the rest of the net.
    original=xy.copy()
    degree=np.bincount(ends.ravel(),minlength=len(xy)).clip(1)
    for _ in range(80):
        d=xy[ends[:,0]]-xy[ends[:,1]]
        dist=np.linalg.norm(d,axis=1)
        correction=d*((lengths-dist)/np.maximum(dist,1e-7))[:,None]
        update=np.zeros_like(xy)
        np.add.at(update,ends[:,0],correction)
        np.add.at(update,ends[:,1],-correction)
        xy+=.72*update/degree[:,None]+.035*(original-xy)
    return vector_art(mol,dict(enumerate(xy)),xyz[:,2],blocks)


def vector_art(mol,positions,depths=None,blocks=None):
    pts=np.array([positions[i] for i in range(mol.GetNumAtoms())])
    # Use a fixed bond scale so exports retain consistent strokes and labels.
    pts=(pts-np.min(pts,axis=0))*[16,-16]
    pts-=np.min(pts,axis=0)-14
    w,h=np.max(pts,axis=0)+14
    root=ET.Element("svg",xmlns="http://www.w3.org/2000/svg",viewBox=f"0 0 {w:.2f} {h:.2f}",
                    width="840",height="660",role="img")
    ET.SubElement(root,"title").text="Framework circuit with outward continuation bonds"
    if depths is None:
        depths=np.ones(len(pts))
    else:
        low,high=np.percentile(depths,[5,95])
        depths=np.clip((np.asarray(depths)-low)/(high-low),0,1) if high-low>.001 else np.ones(len(pts))
    blocks=blocks or [""]*len(pts)
    def attributes(indices):
        block=blocks[indices[0]]
        return {"data-depth":f"{np.mean(depths[list(indices)]):.3f}",
                "data-block":block if block and all(blocks[i]==block for i in indices) else ""}
    def point(p):
        return f"{p[0]:.2f},{p[1]:.2f}"
    shades=ET.SubElement(root,"g",attrib={"data-ring-washes":"", "display":"none"})
    ring_seen=set()
    for bond in mol.GetBonds():
        if not bond.HasProp("drawingRing"):continue
        ring=tuple(int(k) for k in bond.GetProp("drawingRing").split(","))
        if frozenset(ring) in ring_seen:continue
        ring_seen.add(frozenset(ring))
        ET.SubElement(shades,"polygon",points=" ".join(point(pts[k]) for k in ring),attrib=attributes(ring))
    bonds=ET.SubElement(root,"g",fill="none",stroke="#262a30",attrib={"stroke-linecap":"round","stroke-linejoin":"round","stroke-width":"1.65"})
    def line(a,b,**attrs):
        ET.SubElement(group,"path",d="M"+point(a)+" L"+point(b),attrib=attrs)
    for b in mol.GetBonds():
        i,j=b.GetBeginAtomIdx(),b.GetEndAtomIdx()
        p,q=pts[i],pts[j]
        d=q-p
        size=np.linalg.norm(d)
        if size<.01:continue
        group=ET.SubElement(bonds,"g",attrib=attributes([i,j]))
        normal=np.array([-d[1],d[0]])/size
        if b.HasProp("junctionDepth"):
            tip,end=(p,q) if i==b.GetIntProp("junctionCenter") else (q,p)
            if b.GetProp("junctionDepth")=="front":
                ET.SubElement(group,"path",d="M"+point(tip)+" L"+point(end+normal*2.5)+" L"+point(end-normal*2.5)+" Z",fill="currentColor",stroke="none")
            else:
                for t in np.linspace(.14,.96,7):
                    mid=tip+(end-tip)*t
                    line(mid-normal*2.5*t,mid+normal*2.5*t)
        elif not mol.GetAtomWithIdx(i).GetAtomicNum() or not mol.GetAtomWithIdx(j).GetAtomicNum():
            end,start=(q,p) if mol.GetAtomWithIdx(i).GetAtomicNum() else (p,q)
            stem=start+(end-start)*.65
            line(start,stem)
            axis=(end-start)/size
            wave=[stem+normal*t+axis*1.35*math.sin(t*1.5) for t in np.linspace(-4,4,25)]
            ET.SubElement(group,"polyline",points=" ".join(point(x) for x in wave))
        else:
            line(p,q)
            order=b.GetDoubleProp("drawingOrder") if b.HasProp("drawingOrder") else b.GetBondTypeAsDouble()
            if order>=1.5:
                if b.HasProp("drawingRing"):
                    center=pts[[int(k) for k in b.GetProp("drawingRing").split(",")]].mean(0)
                    if normal@(center-(p+q)/2)<0:normal=-normal
                line(p+.15*d+normal*3.6,q-.15*d+normal*3.6)
            if order==3:
                line(p+.15*d-normal*3.6,q-.15*d-normal*3.6)
    labels=ET.SubElement(root,"g",attrib={"text-anchor":"middle","font-family":"Arial, sans-serif","font-size":"14"})
    for atom in mol.GetAtoms():
        symbol=atom.GetSymbol()
        if symbol=="*" or symbol=="C" and not atom.GetIsotope():continue
        i=atom.GetIdx()
        h=atom.GetNumExplicitHs()
        label=(str(atom.GetIsotope()) if atom.GetIsotope() else "")+symbol+("H"+(str(h) if h>1 else "") if h else "")
        charge=atom.GetFormalCharge()
        if charge:label+=(str(abs(charge)) if abs(charge)>1 else "")+("⁺" if charge>0 else "⁻")
        x,y=pts[i]
        halfwidth=4.8*len(label)+1.2
        ET.SubElement(labels,"rect",x=f"{x-halfwidth:.2f}",y=f"{y-7.4:.2f}",width=f"{halfwidth*2:.2f}",height="14.8",fill="white")
        ET.SubElement(labels,"text",x=f"{x:.2f}",y=f"{y+4.8:.2f}",fill={"N":"#17679b","O":"#b32e35","S":"#986619","F":"#247656"}.get(symbol,"#292a32"),attrib=attributes([i])).text=label
    return ET.tostring(root,encoding="unicode")
