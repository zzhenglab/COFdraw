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


def draw_scheme(mol, backbone, segments):
    mol=Chem.Mol(mol)
    n=len(backbone)
    layouts=[]
    lengths=[]
    for k,fragment in enumerate(segments):
        a,b=backbone[k],backbone[(k+1)%n]
        coords=fragment_layout(mol,set(fragment)|{a,b})
        layouts.append(coords)
        lengths.append(np.linalg.norm(coords[b]-coords[a]))
    vertices=segment_vertices(lengths)
    positions={}
    for k,coords in enumerate(layouts):
        a,b=backbone[k],backbone[(k+1)%n]
        moved=rotate_to(coords,coords[a],coords[b]-coords[a],vertices[(k+1)%n]-vertices[k],vertices[k])
        positions.update(moved)
    # Keep both outgoing branches on the exterior. Each fragment has its own
    # conventional bond angles and ring coordinates and receives only a rotation.
    outside=set(range(mol.GetNumAtoms()))-set(positions)
    for k,center in enumerate(backbone):
        attachments=sorted(a.GetIdx() for a in mol.GetAtomWithIdx(center).GetNeighbors() if a.GetIdx() in outside)
        side1=vertices[k]-vertices[(k-1)%n]
        side2=vertices[k]-vertices[(k+1)%n]
        bisector=side1/np.linalg.norm(side1)+side2/np.linalg.norm(side2)
        radial=math.atan2(bisector[1],bisector[0])
        choices=[]
        fixed=[positions[a.GetIdx()]-vertices[k] for a in mol.GetAtomWithIdx(center).GetNeighbors() if a.GetIdx() in positions]
        occupied=[d/np.linalg.norm(d) for d in fixed]
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
            nearby=np.array([p for i,p in positions.items() if i!=center and np.linalg.norm(p-vertices[k])<12])
            # A kink at the end of an imine means its aryl ipso bond may deviate
            # from the direction between node centers. Fit the two outward arms
            # against the actual drawing so no separate phenyl rings look fused.
            candidates=[]
            for offset in range(-40,41,4):
                theta=angle+math.radians(offset)
                target=np.array([math.cos(theta),math.sin(theta)])
                trial=rotate_to(coords,coords[center],coords[first]-coords[center],target,vertices[k])
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
    return vector_art(mol,positions)


def vector_art(mol,positions):
    pts=np.array([positions[i] for i in range(mol.GetNumAtoms())])
    # Use a fixed bond scale so exports retain consistent strokes and labels.
    pts=(pts-np.min(pts,axis=0))*[16,-16]
    pts-=np.min(pts,axis=0)-14
    w,h=np.max(pts,axis=0)+14
    root=ET.Element("svg",xmlns="http://www.w3.org/2000/svg",viewBox=f"0 0 {w:.2f} {h:.2f}",
                    width="840",height="660",role="img")
    ET.SubElement(root,"title").text="Framework circuit with outward continuation bonds"
    bonds=ET.SubElement(root,"g",fill="none",stroke="#262a30",attrib={"stroke-linecap":"round","stroke-linejoin":"round","stroke-width":"1.65"})
    def point(p):
        return f"{p[0]:.2f},{p[1]:.2f}"
    def line(a,b,**attrs):
        ET.SubElement(bonds,"path",d="M"+point(a)+" L"+point(b),attrib=attrs)
    for b in mol.GetBonds():
        i,j=b.GetBeginAtomIdx(),b.GetEndAtomIdx()
        p,q=pts[i],pts[j]
        d=q-p
        size=np.linalg.norm(d)
        if size<.01:continue
        normal=np.array([-d[1],d[0]])/size
        if not mol.GetAtomWithIdx(i).GetAtomicNum() or not mol.GetAtomWithIdx(j).GetAtomicNum():
            end,start=(q,p) if mol.GetAtomWithIdx(i).GetAtomicNum() else (p,q)
            stem=start+(end-start)*.65
            line(start,stem)
            axis=(end-start)/size
            wave=[stem+normal*t+axis*1.35*math.sin(t*1.5) for t in np.linspace(-4,4,25)]
            ET.SubElement(bonds,"polyline",points=" ".join(point(x) for x in wave))
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
        if symbol in ("C","*"):continue
        i=atom.GetIdx()
        h=atom.GetNumExplicitHs()
        label=symbol+("H"+(str(h) if h>1 else "") if h else "")
        charge=atom.GetFormalCharge()
        if charge:label+=(str(abs(charge)) if abs(charge)>1 else "")+("⁺" if charge>0 else "⁻")
        x,y=pts[i]
        halfwidth=4.8*len(label)+1.2
        ET.SubElement(labels,"rect",x=f"{x-halfwidth:.2f}",y=f"{y-7.4:.2f}",width=f"{halfwidth*2:.2f}",height="14.8",fill="white")
        ET.SubElement(labels,"text",x=f"{x:.2f}",y=f"{y+4.8:.2f}",fill={"N":"#17679b","O":"#b32e35","S":"#986619","F":"#247656"}.get(symbol,"#292a32")).text=label
    return ET.tostring(root,encoding="unicode")
