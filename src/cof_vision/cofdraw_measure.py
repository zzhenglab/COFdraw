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
    aperture=prepare_aperture(mol)
    return measure_aperture(mol.GetConformer().GetPositions()[:,:2],aperture)


def prepare_aperture(mol,owners=None,face=None):
    """Record calibration and a specific molecular rim independently of display style."""
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
    rings = [list(r) for r in symmetric_depiction_rings(mol) if len(r) >= 12]
    if face is not None:
        target=set(face)
        rings=[r for r in rings if {owners[i] for i in r}==target]
    if not rings:
        raise ValueError("No closed macrocyclic opening was found. Use a periodic CIF for pore analysis or the 3D atom-distance ruler.")
    ring=max(rings,key=lambda r:_polygon_area(xy[r]))
    return dict(scale=scale,reference_length=reference_length,ring=ring,
                elements=[a.GetAtomicNum() for a in mol.GetAtoms()])


def measure_aperture(coordinates,aperture):
    scale=aperture['scale']; ring=aperture['ring']; elements=aperture['elements']
    xy=np.asarray(coordinates,dtype=float)*scale
    polygon=xy[ring];lo=polygon.min(0);hi=polygon.max(0)
    area=_polygon_area(polygon)
    if not _simple_polygon(polygon) or area <= 1e-8:
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
    indices=[i for i,z in enumerate(elements) if z>0]
    table=Chem.GetPeriodicTable()
    radii=np.array([table.GetRvdw(elements[i]) for i in indices])
    # More distant sites cannot reach the pore even at its maximum possible radius.
    margin=.5*min(hi-lo)+max(radii)
    nearby=np.all((xy[indices]>=lo-margin)&(xy[indices]<=hi+margin),axis=1)
    indices=np.asarray(indices)[nearby];radii=radii[nearby]
    best=(-np.inf,None)
    for chunk in np.array_split(points,max(1,len(points)//256)):
        distances=np.linalg.norm(chunk[:,None,:]-xy[indices][None,:,:],axis=2)-radii
        clearance=np.minimum(distances.min(axis=1),_edge_clearance(chunk,polygon));i=int(clearance.argmax())
        if clearance[i]>best[0]:best=(float(clearance[i]),chunk[i])
    if best[0]<=0:raise ValueError("No positive opening remains after atomic radii are included")
    center = best[1] / scale
    radius = best[0] / scale
    # This line spans the sampled inscribed circle, not a wall-to-wall chord.
    # Retain the uncalibrated coordinate frame for an exact rendered overlay.
    endpoints = [center + [-radius, 0.], center + [radius, 0.]]
    anchor_indices = sorted({0, len(elements) // 2, len(elements) - 1})
    annotation = {"coordinate_system": "drawing_xy_cartesian_y_up",
                  "center": center.tolist(), "radius": float(radius),
                  "endpoints": [point.tolist() for point in endpoints],
                  "angstrom_per_drawing_unit": scale,
                  "coordinate_anchors": [{"atom_index": i, "xy": (xy[i] / scale).tolist()}
                                         for i in anchor_indices],
                  "meaning": "Horizontal diameter of the sampled inscribed circle; not a wall-to-wall width",
                  "label_qualification": "Pore opening"}
    return {"diameter_angstrom":2*best[0],"diameter_nm":best[0]/5,
            "center_angstrom":best[1].tolist(),"annotation":annotation,
            "projected_ring_area_angstrom2":area,"ring_atoms":len(ring),
            "calibration_bond_angstrom":aperture['reference_length'],"grid_points":len(points),
            "grid_diagonal_angstrom":float(np.linalg.norm((hi-lo)/89)),
            "sampling_lower_bound":True,"explicit_hydrogens_included":elements.count(1),
            "method":"Largest sampled circle constrained inside the largest projected macrocycle in the bounded symmetric ring basis, excluding all explicit-atom van der Waals disks",
            "scope":"Largest sampled opening in the drawing plane, calibrated by conventional bond lengths. The value is a lower bound for the projected opening. Periodic pore calculations use a cell and its explicit atomic coordinates."}
