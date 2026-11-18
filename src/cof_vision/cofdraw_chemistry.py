"""Reactive-site inspection and explicit, user-chosen drawing templates."""
from __future__ import annotations
import math
import re
import hashlib
import json
import numpy as np
from rdkit import Chem
from .cofdraw_depiction import compute_2d_coords
from .framework_drawing import ARYL_HYDROXYLAMINE_SMARTS, LINKAGE_LIBRARY, TP, canonical, linkage_fragment, _ports

ROLE_LABELS = {"aryl_aldehyde":"Aromatic aldehyde", "primary_arylamine":"Primary aromatic amine",
               "tp_trialdehyde":"Tp trialdehyde", "acylhydrazide":"Terminal acylhydrazide", "hydrazine":"Hydrazine",
               "aryl_hydroxylamine":"Aromatic hydroxylamine · Ar–NH–OH"}
LINKAGE_LABELS = {"imine":"Imine · C=N", "beta_ketoenamine":"β-Ketoenamine", "acylhydrazone":"Acylhydrazone · C=N–NH–CO",
                  "keto_hydrazide":"Tp keto-enamine–hydrazide", "azine":"Azine · C=N–N=C", "secondary_amine":"Reduced imine · CH₂–NH",
                  "nitrone":"Nitrone · C=N⁺(O⁻)"}
TEMPLATES = [
    {"id":"hcb_3_3","label":"Hexagonal · 3c + 3c","instances":6,"allow_pairs":[[3,3]],"topology_label":"hcb","description":"Alternating trigonal vertices; one representative pore."},
    {"id":"hcb_3_2","label":"Hexagonal · 3c + 2c","instances":12,"allow_pairs":[[3,2],[2,3]],"topology_label":"hcb","description":"Trigonal nodes at six corners; linear linkers on the six edges."},
    {"id":"sql_4_2","label":"Square · 4c + 2c","instances":8,"allow_pairs":[[4,2],[2,4]],"topology_label":"sql","description":"Four tetratopic corners joined by four linear linkers."},
    {"id":"sql_4_4","label":"Square · 4c + 4c","instances":4,"allow_pairs":[[4,4]],"topology_label":"sql","description":"Alternating tetratopic vertices; unused sites remain continuation ports."},
    {"id":"hxl_6_2","label":"Triangular · 6c + 2c","instances":6,"allow_pairs":[[6,2],[2,6]],"topology_label":"hxl","description":"Three hexatopic corners and three edge linkers."},
    {"id":"chain_2_2","label":"Open chain · 2c + 2c","instances":6,"allow_pairs":[[2,2]],"topology_label":None,"description":"Six alternating units with two open ends; no 2D net is assigned."},
    {"id":"alternating_hexagon","label":"Custom six-unit ring","instances":6,"allow_pairs":[],"topology_label":None,"description":"Explicit ring design with any supported functionality ≥2; topology remains unknown."},
]
for template in TEMPLATES:
    template.update(linkages=list(LINKAGE_LIBRARY),periodic_topology_inferred=False,
                    required_functionalities={"allow_pairs":template["allow_pairs"]})

def inspect_monomer(smiles):
    mol=Chem.MolFromSmiles(smiles)
    if mol is None: raise ValueError("Invalid monomer SMILES")
    patterns={"aryl_aldehyde":"[c][CH1](=[O])", "primary_arylamine":"[c][NH2]",
              "aryl_hydroxylamine":ARYL_HYDROXYLAMINE_SMARTS,
              "acylhydrazide":"[NX3;H2]-[NX3;H1]-[CX3](=[OX1])",
              "boronic_acid":"[B]([OH1])[OH1]", "nitrile":"[C]#[N]"}
    groups=[]
    for kind,pattern in patterns.items():
        # Symmetric OH permutations must not count a boronic acid twice.
        hits=mol.GetSubstructMatches(Chem.MolFromSmarts(pattern))
        count=len({hit[0] for hit in hits})
        if count: groups.append({"kind":kind,"label":ROLE_LABELS.get(kind,kind.replace('_',' ').title()),"count":count})
    if canonical(smiles)==canonical(TP): groups.append({"kind":"tp_trialdehyde","label":"Tp trialdehyde","count":3})
    if canonical(smiles)==canonical("NN"): groups.append({"kind":"hydrazine","label":"Hydrazine terminal nitrogens","count":2})
    roles=[]
    for linkage,rule in LINKAGE_LIBRARY.items():
        for role in rule["roles"]:
            count=next((g["count"] for g in groups if g["kind"]==role),0)
            if count:
                roles.append({"linkage":linkage,"role":role,"reactive_sites":count,"compatible":True,
                              "reason":"Recognized functional group; crystallization and reaction conditions are not predicted."})
    tetra=mol.HasSubstructMatch(Chem.MolFromSmarts("[CX4](-[c])(-[c])(-[c])-[c]"))
    return {"functionality":groups,"reactive_groups":groups,"compatible_roles":roles,
            "geometry_hint":"tetrahedral" if tetra else "unspecified",
            "geometry_note":"Tetrahedral tetraaryl-carbon node: a planar template does not establish a 3D net." if tetra else None}

def library_linkages(catalog):
    result={}
    for row in catalog.get("linkage_catalog",[]):
        result[row["id"]]={**row,"supported":False,"roles":[],"role_labels":[],"category":"Other reported linkages",
                           "limits":row.get("qualification") or "Catalogued chemistry; no automatic product transform implemented."}
    for key,rule in LINKAGE_LIBRARY.items():
        result[key]={**result.get(key,{}),**rule,"label":LINKAGE_LABELS[key],"supported":True,
                     "automatic_prepare_supported":True,"runtime_status":"bounded_monomer_transform",
                     "category":"Nitrogen linkages","role_labels":[ROLE_LABELS[r] for r in rule["roles"]]}
    return result

def decorated_monomers(rows):
    output=[]
    for row in rows:
        try:
            metadata=inspect_monomer(row["smiles"])
            output.append({**row,**metadata,"catalog_reactive_groups":row.get("reactive_groups", []),
                           "abbreviation":row.get("abbreviation") or row.get("id", ""),
                           "full_name":row.get("full_name") or row.get("label", ""),
                           "geometry_hint":row.get("geometry_hint") or metadata["geometry_hint"]})
        except (ValueError,KeyError):
            continue
    return output

def best_port_pair(component, angle):
    mol=Chem.MolFromSmiles(component["fragment_smiles"]); compute_2d_coords(mol)
    xy=np.asarray(mol.GetConformer().GetPositions())[:,:2]
    center=xy[[a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum()]].mean(0)
    ports=_ports(mol); candidates=[]
    for incoming,(i,_) in ports.items():
        for outgoing,(j,_) in ports.items():
            if incoming==outgoing: continue
            left,right=xy[i]-center,xy[j]-center
            scale=np.linalg.norm(left)*np.linalg.norm(right)
            if scale<1e-10: continue
            actual=math.degrees(math.acos(float(np.clip(np.dot(left,right)/scale,-1,1))))
            candidates.append((abs(actual-angle),incoming,outgoing,actual))
    if not candidates: raise ValueError("Cannot orient the selected reactive sites")
    error,incoming,outgoing,actual=min(candidates)
    return incoming,outgoing,{"target_angle":angle,"depicted_site_angle":actual,"absolute_angle_difference":error,
                              "basis":"Closest pair in a monomer 2D drawing; inspect incidence before assigning source identity."}

def design_signature(recipe):
    value={"components":{k:{f:c.get(f) for f in ("fragment_smiles","monomer_smiles","role")} for k,c in recipe["components"].items()},
           "instances":[{k:x.get(k) for k in ("id","component","incoming_port","outgoing_port")} for x in recipe["instances"]],
           "connections":recipe["connections"],"linkage":recipe.get("linkage"),"layout_kind":recipe.get("layout_kind")}
    return hashlib.sha256(json.dumps(value,sort_keys=True,separators=(',',':')).encode()).hexdigest()

def refresh_design_identity(recipe):
    selected=recipe.get("selected_template") or {}
    matched=recipe.get("design_template_signature")==design_signature(recipe)
    known=next((t for t in TEMPLATES if t["id"]==selected.get("id")),None) if matched else None
    net=known["topology_label"] if known else None
    components=recipe["components"]
    smiles=sorted({canonical(c["monomer_smiles"]) for c in components.values() if c.get("monomer_smiles")})
    name=re.sub(r"[^A-Za-z0-9_+-]","","-".join(c.get("label",k) for k,c in components.items()))[:80] or "user-design"
    recipe.update(design_topology=net,selected_template=known,
                  design_cofex_notation=".".join(smiles)+f" COFex.{net or 'unk'}.cat0;{name}" if smiles else None,
                  design_identity_qualification="Proposed design notation. Net is the chosen unchanged template, cat0 is assumed; neither is a literature identity or a CIF.")

def prepare_template(monomer_a,monomer_b,linkage,label_a,label_b,instances,template_id=None):
    if linkage not in LINKAGE_LIBRARY: raise ValueError("This linkage is searchable but has no automatic construction rule")
    rule=LINKAGE_LIBRARY[linkage]
    components={key:{**linkage_fragment(smiles,linkage,role),"label":label}
                for key,smiles,role,label in zip(("A","B"),(monomer_a,monomer_b),rule["roles"],(label_a,label_b))}
    counts=[components[k]["connectivity"] for k in ("A","B")]
    if min(counts)<2: raise ValueError("Each component needs at least two compatible reactive sites")
    template=next((dict(t) for t in TEMPLATES if t["id"]==template_id),None)
    if template_id and not template: raise ValueError("Unknown drawing template")
    if template and template["allow_pairs"] and counts not in template["allow_pairs"]:
        raise ValueError(f"{template['label']} requires {template['allow_pairs']}; these monomers have {counts} reactive sites")
    inspection=[inspect_monomer(s) for s in (monomer_a,monomer_b)]
    if template and template["topology_label"] and any(x["geometry_hint"]=="tetrahedral" for x in inspection):
        raise ValueError("This tetrahedral monomer needs a source-supported 3D incidence recipe; choose a custom projection without assigning a planar net")
    positions=[]; kinds=[]; angles=[]
    n=template["instances"] if template else instances
    open_chain=template_id=="chain_2_2"
    if open_chain:
        positions=[[12*(i-(n-1)/2),0] for i in range(n)]; kinds=["A" if i%2==0 else "B" for i in range(n)]; angles=[180]*n
    elif template_id in ("hcb_3_2","sql_4_2","hxl_6_2"):
        corners=n//2; node="A" if counts[0]>counts[1] else "B"; edge="B" if node=="A" else "A"
        vertices=[np.array([25*math.cos(math.pi/2+2*math.pi*i/corners),25*math.sin(math.pi/2+2*math.pi*i/corners)]) for i in range(corners)]
        for i,point in enumerate(vertices):
            positions.extend([point.tolist(),((point+vertices[(i+1)%corners])/2).tolist()]); kinds.extend([node,edge]); angles.extend([180-360/corners,180])
    else:
        positions=[[25*math.cos(math.pi/2+2*math.pi*i/n),25*math.sin(math.pi/2+2*math.pi*i/n)] for i in range(n)]
        kinds=["A" if i%2==0 else "B" for i in range(n)]; angles=[180-360/n]*n
    items=[]; selected=[]
    for i,(position,kind,angle) in enumerate(zip(positions,kinds,angles)):
        incoming,outgoing,fit=best_port_pair(components[kind],angle)
        items.append({"id":f"{kind}{i+1}","component":kind,"incoming_port":incoming,"outgoing_port":outgoing,"layout_position":position})
        selected.append({"instance":items[-1]["id"],**fit})
    connections=[{"a":{"instance":items[i]["id"],"port":items[i]["outgoing_port"]},
                  "b":{"instance":items[(i+1)%n]["id"],"port":items[(i+1)%n]["incoming_port"]},"bond_order":rule["join_bond_order"]}
                 for i in range(n-1 if open_chain else n)]
    net=template["topology_label"] if template else None
    design_name=re.sub(r"[^A-Za-z0-9_+-]","",f"{label_a}-{label_b}")[:80] or "user-design"
    notation=".".join(sorted({canonical(monomer_a),canonical(monomer_b)}))+f" COFex.{net or 'unk'}.cat0;{design_name}"
    return {"components":components,"instances":items,"connections":connections,"linkage":linkage,
            "layout_kind":"open_chain" if open_chain else "closed_pore", "selected_template":template,
            "design_topology":net,"design_cofex_notation":notation,"design_catenation":"cat0",
            "design_identity_qualification":"Proposed design notation. Net is the chosen template, cat0 is assumed; neither is experimental identity or a CIF.",
            "compatibility":{"reactive_site_counts":counts,"roles":rule["roles"],"port_selection":selected,"monomers":inspection}}
