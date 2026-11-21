"""Shape-preserving aryl linker substitutions within an explicit parent cell."""
from collections import defaultdict
from copy import deepcopy
import numpy as np

from .cofdraw_network import inspect_periodic_network


def fluorinated_model(parent, parent_id):
    network = inspect_periodic_network(parent)
    atoms, cell = parent["atoms"], np.asarray(parent["cell"])
    graph = [[] for _ in atoms]
    for e in network["edges"]:
        a, b, t = e["a"], e["b"], np.array(e["translation"])
        graph[a].append((b, t))
        graph[b].append((a, -t))
    def edge_set(result):
        return {(e["a"], e["b"], tuple(e["translation"])) for e in result["edges"]}
    if network["short_contact_count"] or network["dimensions"] != [3]:
        raise ValueError("This crystal needs contact or network review before aryl substitution.")
    elements = [a["element"] for a in atoms]
    colors = elements
    for _ in range(4):
        signatures = [(colors[i], tuple(sorted(colors[j] for j, _ in neighbors))) for i, neighbors in enumerate(graph)]
        keys = {sig: i for i, sig in enumerate(sorted(set(signatures)))}
        colors = [keys[s] for s in signatures]
    candidates = defaultdict(list)
    xyz = np.array([a["fractional"] for a in atoms])
    for c, neighbors in enumerate(graph):
        if elements[c] != "C" or sorted(elements[j] for j, _ in neighbors) != ["C", "C", "H"]:
            continue
        h, shift = next((j,t) for j,t in neighbors if elements[j] == "H")
        if len(graph[h]) != 1:
            continue
        carbons = [(j,t) for j,t in neighbors if elements[j] == "C"]
        # Require a closed six-carbon ring in the unwrapped graph.
        first, last = carbons
        paths = [(first[0], first[1], (c,first[0]))]
        ring = False
        while paths:
            site, image, visited = paths.pop()
            if len(visited) == 6:
                if site == last[0] and np.array_equal(image,last[1]):
                    ring = True
                    break
                continue
            for j, t in graph[site]:
                if elements[j] == "C" and j not in visited and 1.28 <= np.linalg.norm((xyz[j]+t-xyz[site]) @ cell) <= 1.49:
                    paths.append((j,image+t,(*visited,j)))
        if ring:
            direction = (xyz[h]+shift-xyz[c]) @ cell
            candidates[colors[c]].append((h, c, direction / np.linalg.norm(direction)))
    expected = edge_set(network)
    inverse = np.linalg.inv(cell)
    for sites in sorted(candidates.values(), key=lambda group:(len(group),group[0][1]))[:6]:
        trial = deepcopy(parent)
        for h,c,direction in sites:
            trial["atoms"][h] = {"element":"F", "fractional":((xyz[c]+1.35*direction@inverse)%1).tolist(), "occupancy":1.0}
        # Parent partial charges cannot be inherited by a changed composition.
        for a in trial["atoms"]:
            a.pop("charge",None)
        check = inspect_periodic_network(trial)
        if check["short_contact_count"] or edge_set(check) != expected or check["threshold_sensitivity"]["classification_changed"]:
            continue
        name = parent["provenance"]["model_name"]
        trial["provenance"] = {"kind":"linker_design", "model_name":name+" · aryl-F design",
            "parent_model":parent_id, "parent_doi":parent["provenance"].get("doi"), "network_dimension":3,
            "substituted_atom_indices":[h for h,_,_ in sites], "carbon_atom_indices":[c for _,c,_ in sites],
            "qualification":"Unrelaxed substitution at matching local aryl C–H environments in the parent cell. "
                            "Original adjacency and periodic rank checked. Site classes use local connectivity; packing is inherited."}
        return trial
    raise ValueError("No aryl substitution passed the explicit-hydrogen, periodic-bond and contact checks for this cell.")
