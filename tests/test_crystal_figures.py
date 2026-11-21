import json
from pathlib import Path

import numpy as np
import pytest

from cof_vision.cofdraw_crystal_figures import crystal_details
from cof_vision.cofdraw_network import inspect_periodic_network
from cof_vision.cofdraw_science import import_cif


@pytest.mark.parametrize("filename,has_circuit", [
    ("cof-300.cif",True), ("curated/13040N3.cif",True), ("curated/07011N3.cif",False)])
def test_coordinate_illustrations(filename,has_circuit):
    root=Path(__file__).resolve().parents[1]/"data/crystal-models"
    structure=import_cif((root/filename).read_text())
    original=json.dumps(structure,sort_keys=True)
    details=crystal_details(structure,inspect_periodic_network(structure)["edges"])
    assert len(details["junction"]["atoms"])==25
    assert ("ring" in details)==has_circuit
    if has_circuit:
        assert len(details["overlay"]["guides"])==4
        assert len(details["net"]["bonds"])==16
        ring=details["ring"]
        xyz=np.array([[a[k] for k in ("x","y","z")] for a in ring["atoms"]])
        for a,b in ring["bonds"]:
            assert 1 < np.linalg.norm(xyz[a]-xyz[b]) < 1.9
        assert np.linalg.norm(ring["direction"])==pytest.approx(1)
    else:
        from cof_vision.cofdraw_crystal_circuits import joined_reference_circuits
        edges=inspect_periodic_network(structure)["edges"]
        patches=joined_reference_circuits(structure,edges)
        assert list(patches)==[1,2,3,4]
        sites=set()
        for n,patch in patches.items():
            expanded={(i,tuple(image)) for i,image in patch["site_images"]}
            assert sites < expanded and len(expanded)==len(patch["site_images"])
            assert patch["circuits"]==n and np.linalg.norm(patch["direction"])==pytest.approx(1)
            sites=expanded
    assert json.dumps(structure,sort_keys=True)==original


def test_linker_derivative_keeps_periodic_graph_and_has_distinct_identity():
    from collections import Counter
    from cof_vision.cofdraw_crystal_linkers import fluorinated_model
    root=Path(__file__).resolve().parents[1]/"data/crystal-models"
    parent=import_cif((root/"cof-300.cif").read_text())
    parent["provenance"].update(model_name="COF-300",doi="10.1021/ja8096256")
    original=json.dumps(parent,sort_keys=True)
    design=fluorinated_model(parent,"cof-300-c5")
    provenance=design["provenance"]
    assert provenance["kind"]=="linker_design" and "doi" not in provenance
    assert inspect_periodic_network(design)["dimensions"]==[3]
    assert design["cell"]==parent["cell"] and json.dumps(parent,sort_keys=True)==original
    before=Counter(a["element"] for a in parent["atoms"])
    after=Counter(a["element"] for a in design["atoms"])
    count=len(provenance["substituted_atom_indices"])
    assert count==after["F"]-before["F"]==before["H"]-after["H"] and count>0


def test_generated_circuit_preserves_linker_substitution_and_block_views():
    from xml.etree import ElementTree as ET
    from cof_vision.cofdraw_dia import build_dia
    from cof_vision.cofdraw_crystal_figures import adjacent_circuits, chemical_drawing
    tetraamine="Nc1ccc(C(c2ccc(N)cc2)(c2ccc(N)cc2)c2ccc(N)cc2)cc1"
    model,audit=build_dia(tetraamine,"O=Cc1cc(F)c(C=O)cc1")
    edges=audit["periodic_bonds"]
    details=crystal_details(model,edges)
    ring=details["ring"]
    assert len(ring["backbone"])==len(ring["linkers"])==len(details["cage_overlay"]["guides"])==6
    svg=ET.fromstring(chemical_drawing(model,edges,circuit=ring,scheme=True))
    ns={"s":"http://www.w3.org/2000/svg"}
    labels=[t.text for t in svg.findall(".//s:text",ns)]
    assert labels.count("F")==6 and labels.count("N")==12
    assert svg.findall(".//s:polyline",ns)
    assert all(a["component"]=="A" for a in details["blockA"]["atoms"])
    assert any(a["element"]=="F" for a in details["blockB"]["atoms"])
    assert len(details["blocknet"]["bonds"])==len(details["block_overlay"]["spatial_guides"])>4
    topology=details["periodic_topology"]
    assert set(topology["site_ids"])==set(audit["junction_sites"])
    assert len(topology["edges"])==2*len(topology["site_ids"])
    assert any(any(edge[2:]) for edge in topology["edges"])
    joined=adjacent_circuits(ring,model["cell"])
    site_images={(i,tuple(image)) for i,image in joined["site_images"]}
    assert len(joined["site_images"])==len(site_images)
    assert {(i,tuple(image)) for i,image in ring["site_images"]}<site_images
    assert len(site_images)<3*len(ring["site_images"])
    depths=[float(el.attrib["data-depth"]) for el in svg.iter() if "data-depth" in el.attrib]
    assert min(depths)<.1 and max(depths)>.9
    prior=set()
    for n in (1,2,3,4,7):
        patch=adjacent_circuits(ring,model["cell"],n)
        sites={(i,tuple(image)) for i,image in patch["site_images"]}
        assert prior<sites
        prior=sites
        if n==1:continue
        assert len(patch["rings"])==n
        loops=[{tuple(sorted((a,b))) for a,b in zip(r["backbone"],r["backbone"][1:]+r["backbone"][:1])} for r in patch["rings"]]
        assert all(any(bonds & earlier for earlier in loops[:i]) for i,bonds in enumerate(loops) if i)
        figure=ET.fromstring(chemical_drawing(model,edges,circuit=patch,scheme=True))
        assert all(a.get('data-block') in ('','A','B') for a in figure.iter() if 'data-depth' in a.attrib)
        expected_f=sum(model["atoms"][i]["element"]=="F" for i,_ in sites)
        assert sum(t.text=="F" for t in figure.findall(".//s:text",ns))==expected_f
