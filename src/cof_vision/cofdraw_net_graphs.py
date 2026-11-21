"""Local neighborhoods of RCSR nets from the Gavrog Systre archive keys."""
from collections import deque
from functools import lru_cache
import json
from pathlib import Path

import numpy as np


@lru_cache(maxsize=1)
def net_keys():
    path = Path(__file__).with_name('web') / 'cofdraw-net-keys.json'
    return json.loads(path.read_text(encoding='utf8'))


@lru_cache(maxsize=None)
def quotient(code):
    dim, *data = map(int, net_keys()[code].split())
    rows = np.array(data).reshape(-1, dim + 2)
    count = int(rows[:, :2].max())
    incidence = np.zeros((len(rows), count))
    adjacency = [[] for _ in range(count)]
    for i, (a, b, *t) in enumerate(rows):
        a, b = a-1, b-1
        incidence[i, a] -= 1
        incidence[i, b] += 1
        adjacency[a].append((b, tuple(t)))
        adjacency[b].append((a, tuple(-v for v in t)))
    fractional = np.linalg.lstsq(incidence, -rows[:, 2:], rcond=None)[0]
    vectors = incidence @ fractional + rows[:, 2:]
    eigen, axes = np.linalg.eigh(vectors.T @ vectors)
    cell = (axes / np.sqrt(eigen)) @ axes.T
    return adjacency, fractional, cell


def neighbors(adjacency, site):
    k, *image = site
    for j, offset in adjacency[k]:
        yield (j, *(a+b for a,b in zip(image,offset)))


@lru_cache(maxsize=None)
def graph_preview(code):
    adjacency, fractional, cell = quotient(code)
    dim = cell.shape[0]
    origin = (0, *(0 for _ in range(dim)))
    pending = deque([origin])
    sites = {origin:0}
    links=0
    while pending and (len(sites)<(41 if dim==2 else 36) or links-len(sites)+1<3):
        for other in neighbors(adjacency, pending.popleft()):
            if other not in sites:
                links+=sum(q in sites for q in neighbors(adjacency,other))
                sites[other] = len(sites)
                pending.append(other)
    edges = [[i, sites[q]] for p,i in sites.items() for q in neighbors(adjacency,p)
             if q in sites and i < sites[q]]
    points = np.array([(fractional[k]+t) @ cell for k,*t in sites])
    if dim == 3:
        # Favor a view where separate nodes remain visible and edges retain length.
        i,j = np.triu_indices(len(points),1)
        delta = points[i]-points[j]
        bond = points[np.array(edges)[:,0]]-points[np.array(edges)[:,1]]
        index = np.arange(180)+.5
        z = 1-2*index/len(index)
        radius = np.sqrt(1-z*z)
        views = np.column_stack((radius*np.cos(index*2.399963),radius*np.sin(index*2.399963),z))
        def clearances(vectors):
            return np.sqrt(np.maximum(0,(vectors*vectors).sum(1)[:,None]-(vectors@views.T)**2))
        view = views[np.argmax(clearances(delta).min(0)+.5*clearances(bond).min(0))]
        horizontal = np.cross(view,[0,0,1])
        horizontal /= np.linalg.norm(horizontal)
        points = points @ np.array([horizontal,np.cross(view,horizontal),view]).T
    xy = points[:,:2] - (np.max(points[:,:2],axis=0)+np.min(points[:,:2],axis=0))/2
    xy *= min(175/max(np.ptp(xy[:,0]),1e-8),100/max(np.ptp(xy[:,1]),1e-8))
    depth = points[:,2] if dim == 3 else np.zeros(len(xy))
    depth = (depth-depth.min())/max(np.ptp(depth),1e-8)
    return dict(nodes=[[round(float(x+100),1),round(float(y+58),1),round(float(d),2),len(adjacency[k])]
                       for (x,y),d,(k,*_) in zip(xy,depth,sites)],edges=edges,dimension=dim)
