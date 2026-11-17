"""Drawing-plane opening estimate. Never an experimental or periodic pore size."""
from __future__ import annotations

import numpy as np
from rdkit import Chem

from .depiction_rings import symmetric_depiction_rings


def _polygon_area(polygon):
    p = polygon - polygon.mean(axis=0)
    return abs(float(np.sum(p[:, 0] * np.roll(p[:, 1], -1) - p[:, 1] * np.roll(p[:, 0], -1)))) * .5


def _simple_polygon(polygon):
    """Reject zero edges and nonadjacent intersections in a ring projection."""
    ends = np.roll(polygon, -1, axis=0)
    if np.any(np.linalg.norm(ends - polygon, axis=1) < 1e-8):
        return False
    cross = lambda a, b: a[..., 0] * b[..., 1] - a[..., 1] * b[..., 0]
    for i, (a, b) in enumerate(zip(polygon, ends)):
        js = np.arange(i + 2, len(polygon))
        if i == 0:
            js = js[js != len(polygon) - 1]
        if not len(js):
            continue
        c, d = polygon[js], ends[js]
        overlaps = np.all(np.maximum(np.minimum(a, b), np.minimum(c, d)) <=
                          np.minimum(np.maximum(a, b), np.maximum(c, d)) + 1e-10, axis=1)
        ab_c, ab_d = cross(b - a, c - a), cross(b - a, d - a)
        cd_a, cd_b = cross(d - c, a - c), cross(d - c, b - c)
        if np.any(overlaps & (ab_c * ab_d <= 1e-12) & (cd_a * cd_b <= 1e-12)):
            return False
    return True


def _edge_clearance(points, polygon):
    """Euclidean distance to the polygon boundary; centres are already inside."""
    starts = polygon
    vectors = np.roll(polygon, -1, axis=0) - starts
    relative = points[:, None, :] - starts[None, :, :]
    projection = np.sum(relative * vectors[None, :, :], axis=2) / np.sum(vectors * vectors, axis=1)
    projection = np.clip(projection, 0, 1)
    nearest = starts[None, :, :] + projection[:, :, None] * vectors[None, :, :]
    return np.linalg.norm(points[:, None, :] - nearest, axis=2).min(axis=1)


def measure_drawing(recipe):
    from .framework_drawing import build_patch
    mol, _ = build_patch(recipe)
    if mol.GetNumAtoms() > 1000:
        raise ValueError("Drawing measurement supports at most 1,000 atoms")
    xy = np.asarray(mol.GetConformer().GetPositions())[:, :2]
    if not np.isfinite(xy).all():
        raise ValueError("Drawing coordinates must be finite")
    bonds = [b for b in mol.GetBonds() if b.GetBeginAtom().GetAtomicNum() and b.GetEndAtom().GetAtomicNum()]
    aromatic = [b for b in bonds if b.GetIsAromatic()]
    reference = aromatic or bonds
    lengths = [np.linalg.norm(xy[b.GetBeginAtomIdx()]-xy[b.GetEndAtomIdx()]) for b in reference]
    if not lengths or min(lengths) < 1e-8 or np.median(lengths) < .01:
        raise ValueError("No usable bond-length calibration in this drawing")
    # The depiction coordinates have no physical unit. This is an explicit
    # conventional calibration, not a crystallographic determination.
    reference_length = 1.40 if aromatic else 1.50
    scale = reference_length / float(np.median(lengths))
    xy = xy * scale
    rings = [list(r) for r in symmetric_depiction_rings(mol) if len(r) >= 12]
    if not rings:
        raise ValueError("No closed macrocyclic opening was found. Use a periodic CIF for pore analysis or the 3D atom-distance ruler.")
    def area(indices):
        return _polygon_area(xy[indices])
    ring = max(rings, key=area)
    polygon=xy[ring];lo=polygon.min(0);hi=polygon.max(0)
    if not _simple_polygon(polygon) or area(ring) <= 1e-8:
        raise ValueError("The largest macrocycle has a crossed or degenerate projection; clean the layout before measuring")
    grid=np.stack(np.meshgrid(np.linspace(lo[0],hi[0],90),np.linspace(lo[1],hi[1],90)),axis=-1).reshape(-1,2)
    inside=np.zeros(len(grid),dtype=bool)
    for i in range(len(polygon)):
        a,b=polygon[i],polygon[(i+1)%len(polygon)]
        crossing=(a[1]>grid[:,1])!=(b[1]>grid[:,1])
        x=a[0]+(grid[:,1]-a[1])*(b[0]-a[0])/(b[1]-a[1]+1e-15)
        inside ^= crossing & (grid[:,0]<x)
    points=grid[inside]
    if not len(points):raise ValueError("The projected ring has no usable interior")
    indices=[a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum()>0]
    table=Chem.GetPeriodicTable()
    radii=np.array([table.GetRvdw(mol.GetAtomWithIdx(i).GetAtomicNum()) for i in indices])
    best=(-np.inf,None)
    for chunk in np.array_split(points,max(1,len(points)//256)):
        distances=np.linalg.norm(chunk[:,None,:]-xy[indices][None,:,:],axis=2)-radii
        clearance=np.minimum(distances.min(axis=1),_edge_clearance(chunk,polygon));i=int(clearance.argmax())
        if clearance[i]>best[0]:best=(float(clearance[i]),chunk[i])
    if best[0]<=0:raise ValueError("No positive opening remains after atomic radii are included")
    return {"diameter_angstrom":2*best[0],"center_angstrom":best[1].tolist(),
            "projected_ring_area_angstrom2":area(ring),"ring_atoms":len(ring),
            "calibration_bond_angstrom":reference_length,"grid_points":len(points),
            "grid_diagonal_angstrom":float(np.linalg.norm((hi-lo)/89)),
            "sampling_lower_bound":True,"explicit_hydrogens_included":sum(a.GetAtomicNum()==1 for a in mol.GetAtoms()),
            "method":"Largest sampled circle constrained inside the largest projected macrocycle in the bounded symmetric ring basis, excluding all explicit-atom van der Waals disks",
            "scope":"2D drawing estimate with conventional bond-length scaling; implicit H, stacking, and out-of-plane distortion are omitted. Sampling gives a lower bound for this projection, not a physical pore diameter. Use a periodic 3D model for physical pore screening."}
