"""Place molecular bridges on the hexagonal diamond (ideal wurtzite) quotient graph."""
from collections import defaultdict
from itertools import combinations
import math

import numpy as np
from rdkit import Chem

from .cofdraw_dia import _assemble_edges, _basis


# 60-degree in-plane basis, c/a=sqrt(8/3). Bond length is one in this basis.
# A and B sublattices are occupied by chemically identical four-connected hubs.
HUBS = np.array([[0,0,0], [0,0,3/8], [1/3,1/3,1/2], [1/3,1/3,7/8]])
EDGES = [(0,1,(0,0,0)), (2,3,(0,0,0)),
         (0,3,(0,0,-1)), (0,3,(-1,0,-1)), (0,3,(0,-1,-1)),
         (2,1,(0,0,0)), (2,1,(1,0,0)), (2,1,(0,1,0))]


STACKS = {"4h":"ABAC","6h":"ABCACB"}

def polytype_graph(stack):
    layers = len(stack)
    hubs = np.array([[x,x,(i+d)/layers] for i,letter in enumerate(stack)
                     for x in ["ABC".index(letter)/3] for d in (0,.75)])
    scale = np.array([[math.sqrt(8/3),0,0],[math.sqrt(2/3),math.sqrt(2),0],[0,0,4*layers/3]])
    edges=[]
    for a in range(0,len(hubs),2):
        for b in range(1,len(hubs),2):
            for shift in np.ndindex(3,3,3):
                image=np.asarray(shift)-1
                if abs(np.linalg.norm((hubs[b]+image-hubs[a])@scale)-1)<1e-8:
                    edges.append((a,b,tuple(map(int,image))))
    return hubs,edges


def _roll_states(edge, ends, neighbours, base, cell, hubs, edges):
    """Score pairs of bridge rolls on each actual junction in the periodic graph."""
    axis = _basis(base[ends[1]])
    positions, at_hub = [], defaultdict(list)
    for i,(a,b,t) in enumerate(edges):
        target = (hubs[b]+t-hubs[a]) @ cell
        dest = _basis(target)
        states = []
        for phi in np.arange(18)*math.pi/9:
            c,s = math.cos(phi),math.sin(phi)
            roll = np.array([[1,0,0],[0,c,-s],[0,s,c]])
            states.append(base @ (axis.T @ roll @ dest))
        coords = np.asarray(states)
        positions.append(coords)
        for hub,end,j in zip((a,b), ends, neighbours):
            centered = coords - coords[:,end:end+1]
            dist = np.linalg.norm(centered[0],axis=1)
            cap = {n.GetIdx() for n in edge.GetAtomWithIdx(end).GetNeighbors() if n.GetAtomicNum()==1}
            arm = [k for k in np.flatnonzero((dist > .01) & (dist < 5.)) if k not in cap]
            radii = np.array([Chem.GetPeriodicTable().GetRvdw(edge.GetAtomWithIdx(int(k)).GetAtomicNum()) for k in arm])
            at_hub[hub].append((i, centered[:,j]/dist[j], centered[:,arm], radii))
    costs, angles = {}, {}
    for arms in at_hub.values():
        for left,right in combinations(sorted(arms,key=lambda a:a[0]),2):
            i,u,x,r = left
            j,v,y,s = right
            err = abs(np.degrees(np.arccos(np.clip(u@v.T,-1,1))) - 109.4712)
            dif = x[:,None,:,None,:] - y[None,:,None,:,:]
            overlap = np.maximum(.78*(r[:,None]+s[None,:])-np.linalg.norm(dif,axis=4),0)
            penalty = (err/10)**2 + 2*np.maximum(err-16,0)**2 + 30*np.sum(overlap**2,axis=(2,3))
            costs[i,j] = costs.get((i,j),0) + penalty
            angles[i,j] = np.maximum(angles.get((i,j),0),err)
    return positions,costs,angles


def lon_candidates(edge, ends, neighbours, topology="lon"):
    hubs,edges=(HUBS,EDGES) if topology=="lon" else polytype_graph(STACKS[topology])
    xyz = edge.GetConformer().GetPositions()
    base = xyz-xyz[ends[0]]
    length = np.linalg.norm(base[ends[1]])
    a = math.sqrt(8/3)*length
    cell = np.array([[a,0,0],[a/2,a*math.sqrt(3)/2,0],[0,0,2*len(hubs)*length/3]])
    positions,costs,angles = _roll_states(edge,ends,neighbours,base,cell,hubs,edges)
    adjacency = defaultdict(list)
    for (i,j), cost in costs.items():
        adjacency[i].append((j,cost))
        adjacency[j].append((i,cost.T))
    rng = np.random.default_rng(7013)
    trials = {}
    for _ in range(40):
        rolls = rng.integers(0,len(positions[0]),size=len(edges))
        for _ in range(12):
            prev = rolls.copy()
            for i in rng.permutation(len(rolls)):
                rolls[i] = np.argmin(sum(cost[:,rolls[j]] for j,cost in adjacency[i]))
            if np.array_equal(rolls,prev):
                break
        angle = max(tab[rolls[i],rolls[j]] for (i,j),tab in angles.items())
        if angle <= 18:
            score = sum(c[rolls[i], rolls[j]] for (i,j),c in costs.items())
            trials[tuple(rolls)] = (score,angle)
    for rolls, (_,angle) in sorted(trials.items(),key=lambda t:t[1])[:6]:
        centers = hubs@cell
        copies = [(a,b,centers[a]+positions[i][rolls[i]]) for i,(a,b,_) in enumerate(edges)]
        structure,bonds = _assemble_edges(edge, ends, cell, centers, copies)
        yield angle,structure,bonds
