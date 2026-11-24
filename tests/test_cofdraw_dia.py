"""The 3D design route must derive each single-net crystal from both precursor graphs."""
from collections import Counter, defaultdict

import numpy as np
import pytest
from fastapi.testclient import TestClient
from rdkit import Chem

from cof_vision.cofdraw_desktop import create_desktop_app
from cof_vision.cofdraw_network import inspect_periodic_network
from cof_vision.cofdraw_periodic_build import edge_key
from cof_vision.cofdraw_science import import_cif


TAPM = "Nc1ccc(C(c2ccc(N)cc2)(c2ccc(N)cc2)c2ccc(N)cc2)cc1"
TETRA_CHO = "O=Cc1ccc(C(c2ccc(C=O)cc2)(c2ccc(C=O)cc2)c2ccc(C=O)cc2)cc1"


def assert_opposed_terephthalaldehyde_ends(structure, bonds, expected_linkers):
    """Measure every para-linker's two outgoing C=N vectors in its lattice images."""
    atoms = structure['atoms']
    linker = 'B' if structure['provenance']['node_index'] == 0 else 'A'
    remaining = {i for i, atom in enumerate(atoms) if atom['component'] == linker}
    neighbors, imines = defaultdict(set), {}
    fractional = np.array([atom['fractional'] for atom in atoms])
    for bond in bonds:
        a, b = bond['a'], bond['b']
        if a in remaining and b in remaining:
            neighbors[a].add(b)
            neighbors[b].add(a)
        elif (a in remaining) != (b in remaining):
            assert bond['order'] == 2
            vector = (fractional[b] + bond['translation'] - fractional[a]) @ structure['cell']
            imines[a if a in remaining else b] = vector if a in remaining else -vector
    count = 0
    while remaining:
        component = {remaining.pop()}
        front = component
        while front:
            front = {j for i in front for j in neighbors[i]} - component
            component |= front
        remaining -= component
        left, right = [imines[i] for i in component if i in imines]
        assert left @ right / (np.linalg.norm(left) * np.linalg.norm(right)) < -.98
        count += 1
    assert count == expected_linkers


@pytest.mark.parametrize("net,reverse,cells", [("dia", False, 1), ("lon", True, 2)])
def test_two_distinct_tetrahedral_nodes_keep_si_formula_and_input_order(tmp_path, net, reverse, cells):
    silicon = TETRA_CHO.replace("ccc(C(", "ccc([Si](")
    feed = [silicon, TAPM] if reverse else [TAPM, silicon]
    with TestClient(create_desktop_app(tmp_path)) as client:
        result = client.post("/api/drawing/science/build-reticular",
                             json=dict(monomer_a=feed[0], monomer_b=feed[1], topology=net))
        assert result.status_code == 200, result.text
        model = result.json()
        structure = model["structure"]
        assert structure["provenance"]["precursor_connectivities"] == [4,4]
        assert structure["provenance"]["precursors"] == [Chem.MolToSmiles(Chem.MolFromSmiles(s)) for s in feed]
        assert model["cell_validation"]["precursor_stoichiometry"] == [cells,cells]
        assert model["cell_validation"]["network_dimensions"] == [3]
        assert Counter(a["element"] for a in structure["atoms"]) == {
            element: n*cells for element,n in {"C":53,"H":36,"N":4,"Si":1}.items()}
        assert model["cell_validation"]["cell_composition_verified"]


@pytest.mark.parametrize("a,b,linkage,counts", [
    (TAPM,"O=Cc1ccc(C=O)cc1","imine", dict(C=82,H=56,N=8)),
    ("O=Cc1c(F)cc(C=O)cc1",TAPM,"imine", dict(C=82,H=52,N=8,F=4)),
    (TETRA_CHO,"NNC(=O)c1ccc(C(=O)NN)cc1","acylhydrazone", dict(C=90,H=64,N=16,O=8)),
    (TAPM,"O=Cc1ccc(C=O)cc1","secondary_amine", dict(C=82,H=72,N=8)),
])
def test_draw_custom_periodic_product_without_template(tmp_path,a,b,linkage,counts):
    with TestClient(create_desktop_app(tmp_path)) as client:
        route = "/api/drawing/science"
        resp = client.post(route+"/build-dia",json=dict(monomer_a=a,monomer_b=b,linkage=linkage))
        assert resp.status_code==200, resp.text
        data = resp.json()
        m = data["structure"]
        assert data["available"] and m["provenance"]["kind"]=="monomer_dia_design"
        assert Counter(x["element"] for x in m["atoms"])==counts
        cif = client.get(data["downloads"]["cif"])
        net = inspect_periodic_network(import_cif(cif.text))
        assert net["dimensions"]==[3] and net["quotient_component_count"]==1 and not net["short_contact_count"]
        assert not net["threshold_sensitivity"]["classification_changed"]
        explicit = data["cell_validation"]["periodic_bonds"]
        keys=lambda edges:{edge_key(e["a"],e["b"],e["translation"]) for e in edges}
        assert keys(explicit)==keys(net["edges"])
        assert sum(e["order"]==2 for e in explicit) == (0 if linkage=="secondary_amine" else 16 if linkage=="acylhydrazone" else 8)
        if linkage == 'imine':
            assert_opposed_terephthalaldehyde_ends(m, explicit, 4)
        repeated = client.post(route+"/reports/"+data["report_id"]+"/view",json={"repeats":[2,1,1]})
        assert repeated.status_code==200 and len(repeated.json()["viewer"]["atoms"])==2*sum(counts.values())


def test_planar_node_rejected(tmp_path):
    with TestClient(create_desktop_app(tmp_path)) as client:
        resp=client.post("/api/drawing/science/build-dia", json=dict(monomer_a="Nc1c(N)cc(N)c(N)c1",monomer_b="O=Cc1ccc(C=O)cc1"))
        assert resp.status_code==422 and "tetrahedral" in resp.text


def test_interwoven_networks_match_requested_count_and_have_no_crosslinks(tmp_path):
    with TestClient(create_desktop_app(tmp_path)) as client:
        resp=client.post("/api/drawing/science/build-dia",json={"monomer_a":TAPM,
                         "monomer_b":"O=Cc1ccc(C=O)cc1","network_count":3})
        assert resp.status_code==200,resp.text
        data=resp.json()
        a=data["cell_validation"]
        assert a["network_count"]==3 and len(a["packing"]["fractional_offsets"])==3
        assert a["packing"]["intercopy_clearance_lower_bound_angstrom"]>=0
        atoms=data["structure"]["atoms"]
        assert all(atoms[b["a"]]["network"]==atoms[b["b"]]["network"] for b in a["periodic_bonds"])
        cif=client.get(data["downloads"]["cif"])
        assert inspect_periodic_network(import_cif(cif.text))["quotient_components_by_dimension"]=={"3":3}


@pytest.mark.parametrize('topology,linkers', [('4h',16), ('6h',24)])
def test_opposed_linkers_in_diamond_polytypes(topology, linkers):
    from cof_vision.cofdraw_dia import build_dia
    model, validation = build_dia(TAPM, 'O=Cc1ccc(C=O)cc1', topology=topology)
    assert_opposed_terephthalaldehyde_ends(model, validation['periodic_bonds'], linkers)


def test_dia_and_lon_have_distinct_periodic_junction_graphs_for_same_custom_linker(tmp_path):
    from cof_vision.cofdraw_dia import _junction_sequence
    with TestClient(create_desktop_app(tmp_path)) as client:
        structures = {}
        for net,expected,sites in [("dia",[4,12,24,42,64],146),("lon",[4,12,25,44,67],292)]:
            resp = client.post("/api/drawing/science/build-reticular",json={"monomer_a":TAPM,
                               "monomer_b":"O=Cc1c(F)cc(C=O)cc1","topology":net})
            assert resp.status_code==200,resp.text
            data=resp.json()
            structures[net]=data["structure"]
            assert len(data["structure"]["atoms"])==sites
            assert_opposed_terephthalaldehyde_ends(data["structure"], data["cell_validation"]["periodic_bonds"],
                                               4 if net == "dia" else 8)
            assert data["cell_validation"]["junction_coordination_sequence"]==expected
            imported=inspect_periodic_network(import_cif(client.get(data["downloads"]["cif"]).text))
            assert imported["short_contact_count"]==0 and imported["dimensions"]==[3]
            assert _junction_sequence(imported["edges"],data["cell_validation"]["junction_sites"])==expected
            assert data["chemical_framework_svg"] and data["chemical_circuit_svg"]
        assert structures["dia"]["cell"]!=structures["lon"]["cell"]
