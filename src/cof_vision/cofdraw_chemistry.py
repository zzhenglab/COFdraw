"""Reactive-site inspection and explicit, user-chosen drawing templates."""
from __future__ import annotations
from copy import deepcopy
import math
import re
import hashlib
import json
from functools import lru_cache
import numpy as np
from rdkit import Chem
from .cofdraw_depiction import compute_2d_coords
from .cofdraw_linkage_groups import linkage_functional_groups
from .framework_drawing import ALKYL_ALDEHYDE, ALKYL_AMINE, ALKYL_VARIANTS, ARYL_HYDROXYLAMINE_SMARTS, ARYL_BORONIC, CARBON_SITES, CYCLIC_ANHYDRIDE, LINKAGE_LIBRARY, TP, canonical, catechol_pairs, linkage_fragment, _ports

ROLE_LABELS = {"aryl_aldehyde":"Aromatic aldehyde", "primary_arylamine":"Primary aromatic amine",
               "primary_alkylamine":"Primary aliphatic amine","alkyl_aldehyde":"Aliphatic aldehyde","carbon_aldehyde":"Aromatic or aliphatic aldehyde",
               "tp_trialdehyde":"Tp trialdehyde", "acylhydrazide":"Terminal acylhydrazide", "hydrazine":"Hydrazine",
               "aryl_hydroxylamine":"Aromatic hydroxylamine · Ar–NH–OH","aryl_acetonitrile":"Arylacetonitrile · Ar–CH₂–CN",
               "methyltriazine":"Methyl on 1,3,5-triazine","aryl_carboxyl":"Aromatic carboxylic acid or acyl chloride",
               "aryl_cyclic_anhydride":"Aromatic cyclic anhydride","aryl_vicinal_diol":"Aromatic ortho-diol pair","aryl_boronic_acid":"Arylboronic acid"}
LINKAGE_LABELS = {"imine":"Imine · C=N", "beta_ketoenamine":"β-Ketoenamine", "acylhydrazone":"Acylhydrazone · C=N–NH–CO",
                  "alkyl_imine":"Alkylamine-derived imine","alkyl_ketoenamine":"Alkylamine-derived β-ketoenamine",
                  "reduced_alkyl_imine":"Reduced alkylamine-derived imine",
                  "keto_hydrazide":"Tp keto-enamine–hydrazide", "azine":"Azine · C=N–N=C", "secondary_amine":"Reduced imine · CH₂–NH",
                  "nitrone":"Nitrone · C=N⁺(O⁻)","cyano_vinylene_sp2_carbon":"Cyanovinylene · CH=C(CN)",
                  "vinylene_methyltriazine":"Vinylene · methyltriazine condensation","amide":"Amide · C(O)–NH",
                  "imide":"Imide · cyclized dianhydride–amine","boronate_ester":"Boronate ester · cyclic B(O)₂"}
from . import cofdraw_spatial as spatial
from . import cofdraw_planar as planar

TEMPLATES = [
    {"id":"branched_4_3","label":"3D bor · tetrahedral + tritopic","dimension":3,"parent_net":"bor","instances":6,"allow_pairs":[[3,4],[4,3]],"topology_label":None,"linkages":["imine","boronate_ester"],"description":"Finite projection of one to three circuits sharing junctions along a proposed bor incidence graph. The pair alone does not identify its crystal net."},
    {"id":"hcb_3_3","label":"Hexagonal · 3c + 3c","instances":6,"allow_pairs":[[3,3]],"topology_label":"hcb","description":"Alternating trigonal vertices; one representative pore."},
    {"id":"hcb_3_2","label":"Hexagonal · 3c + 2c","instances":12,"allow_pairs":[[3,2],[2,3]],"topology_label":"hcb","description":"Trigonal nodes at six corners; linear linkers on the six edges."},
    {"id":"sql_4_2","label":"sql · Square · 4c + 2c","instances":8,"allow_pairs":[[4,2],[2,4]],"topology_label":"sql","description":"Four tetratopic corners joined by four linear linkers."},
    {"id":"kgm_4_2","label":"kgm · 2D kagome","instances":30,"allow_pairs":[[2,4],[4,2]],"topology_label":"kgm","description":"Hexagonal pore surrounded by six triangular pores. Choose a single face or connected neighborhood under Drawing extent."},
    {"id":"sql_4_4","label":"Square · 4c + 4c","instances":4,"allow_pairs":[[4,4]],"topology_label":"sql","description":"Alternating tetratopic vertices; unused sites remain continuation ports."},
    {"id":"hxl_6_2","label":"Triangular · 6c + 2c","instances":6,"allow_pairs":[[6,2],[2,6]],"topology_label":"hxl","description":"Three hexatopic corners and three edge linkers."},
    {"id":"chain_2_2","label":"Open chain · 2c + 2c","instances":6,"allow_pairs":[[2,2]],"topology_label":None,"description":"Six alternating units with two open ends; no 2D net is assigned."},
    {"id":"alternating_hexagon","label":"Custom six-unit ring","instances":6,"allow_pairs":[],"topology_label":None,"description":"Explicit ring design with any supported functionality ≥2; topology remains unknown."},
]
TEMPLATES[-1:-1]=planar.templates()+spatial.templates()
# The studio opens the full molecular cage/star by default. Shorter circuit
# extents are also valid for these parent nets.
for template in TEMPLATES:
    if template["id"] in ("dia_4_2", "pts_4_4"):
        template["instances"] = 22
        template["linkages"] = list(LINKAGE_LIBRARY)
TEMPLATES.extend(dict(id=key,label=f'Local drawing · {n}-block ring',instances=n,allow_pairs=[],topology_label=None,
                      description='An explicitly chosen finite ring. Adjacent rings can share an edge when both junctions have a spare site.')
                 for key,n in planar.LOCAL_FACES.items())
for template in TEMPLATES:
    template.update(linkages=template.get("linkages",list(LINKAGE_LIBRARY)),periodic_topology_inferred=False,
                    required_functionalities={"allow_pairs":template["allow_pairs"]})

@lru_cache(maxsize=512)
def _four_arm_geometry(smiles, pattern):
    """Screen the directions of the four reactive sites of this specific role."""
    from rdkit.Chem import AllChem

    mol = Chem.MolFromSmiles(smiles)
    sites = sorted({hit[0] for hit in mol.GetSubstructMatches(Chem.MolFromSmarts(pattern))})
    if len(sites) != 4 or mol.GetNumHeavyAtoms() > 140:
        return "unassigned"
    mol = Chem.AddHs(mol)
    params = AllChem.ETKDGv3()
    for seed in (300, 113, 42):
        params.randomSeed = seed
        params.maxIterations = 30
        if AllChem.EmbedMolecule(mol, params) < 0:
            continue
        arms = mol.GetConformer().GetPositions()[list(sites)]
        arms -= arms.mean(axis=0)
        lengths = np.linalg.norm(arms, axis=1)
        if min(lengths) < 1e-6:
            continue
        arms /= lengths[:, None]
        cosines = np.sort((arms @ arms.T)[np.triu_indices(4, 1)])
        singular = np.linalg.svd(arms, compute_uv=False)
        if singular[2] < .20 * singular[0] and cosines[1] < -.85 and max(abs(cosines[2:])) < .75:
            return "square_planar"
        if cosines[0] > -.78 and cosines[-1] < .25 and singular[2] > .48 * singular[0]:
            return "tetrahedral"
    return "unassigned"


_TETRAARYLETHENE = Chem.MolFromSmarts('[c]-[C;!R](-[c])=[C;!R](-[c])-[c]')


def kagome_roles(mol, patterns):
    """Each of the four arms must terminate in a different reactive site."""
    from .cofdraw_mixed import linear_arms
    cores=mol.GetSubstructMatches(_TETRAARYLETHENE)
    roles={}
    for role in ('aryl_aldehyde','primary_arylamine','aryl_acetonitrile','aryl_hydroxylamine'):
        hits=mol.GetSubstructMatches(Chem.MolFromSmarts(patterns[role]))
        if len(hits)==2 and linear_arms(mol,hits):
            roles[role]='spacer'
        if len(hits)!=4:
            continue
        for arm1,left,arm2,right,arm3,arm4 in cores:
            arms={arm1,arm2,arm3,arm4}
            starts=[]
            for hit in hits:
                path=min((Chem.GetShortestPath(mol,center,hit[1]) for center in (left,right)),key=len)
                starts.append(path[1])
            if set(starts)==arms:
                roles[role]='node'
                break
    return roles


def kagome_match(rows,roles):
    return {row.get('kagome_roles',{}).get(role) for row,role in zip(rows,roles)}=={'node','spacer'}


def inspect_monomer(smiles):
    return deepcopy(_inspect_monomer(smiles))


@lru_cache(maxsize=2048)
def _inspect_monomer(smiles):
    mol=Chem.MolFromSmiles(smiles)
    if mol is None: raise ValueError("Invalid monomer SMILES")
    patterns={"aryl_aldehyde":"[c][CH1](=[O])", "primary_arylamine":"[c][NX3;H2;+0]",
              "primary_alkylamine":ALKYL_AMINE,"alkyl_aldehyde":ALKYL_ALDEHYDE,
              "aryl_hydroxylamine":ARYL_HYDROXYLAMINE_SMARTS,
              "acylhydrazide":"[NX3;H2]-[NX3;H1]-[CX3](=[OX1])",
              "boronic_acid":"[B]([OH1])[OH1]", "aryl_boronic_acid":ARYL_BORONIC, "nitrile":"[C]#[N]", "aryl_cyclic_anhydride":CYCLIC_ANHYDRIDE, **CARBON_SITES}
    groups=[]
    directions = {}
    for kind,pattern in patterns.items():
        # Symmetric OH permutations must not count a boronic acid twice.
        hits=mol.GetSubstructMatches(Chem.MolFromSmarts(pattern))
        # Distinct formyl carbons or nitrogens can share a saturated carbon.
        count=len({hit[1 if kind in ('primary_alkylamine','alkyl_aldehyde') else 0] for hit in hits})
        if count == 4:
            directions[kind] = _four_arm_geometry(canonical(smiles), pattern)
        if count == 2 and kind in ("aryl_aldehyde", "primary_arylamine", "aryl_hydroxylamine"):
            if not mol.GetNumConformers():
                compute_2d_coords(mol)
            xy = mol.GetConformer().GetPositions()
            arms = [xy[h[1]] - xy[h[0]] for h in hits]
            cosine = np.dot(*arms) / (np.linalg.norm(arms[0]) * np.linalg.norm(arms[1]))
            directions[kind] = "linear" if cosine <= math.cos(math.radians(145)) else "bent"
        if count: groups.append({"kind":kind,"label":ROLE_LABELS.get(kind,kind.replace('_',' ').title()),"count":count})
    diols=catechol_pairs(mol)
    if diols:groups.append(dict(kind='aryl_vicinal_diol',label=ROLE_LABELS['aryl_vicinal_diol'],count=len(diols)))
    if canonical(smiles)==canonical(TP): groups.append({"kind":"tp_trialdehyde","label":"Tp trialdehyde","count":3})
    if canonical(smiles)==canonical("NN"): groups.append({"kind":"hydrazine","label":"Hydrazine terminal nitrogens","count":2})
    roles=[]
    for linkage,rule in LINKAGE_LIBRARY.items():
        for role in rule["roles"]:
            count=next((g["count"] for g in groups if g["kind"]==role),0)
            if role=="carbon_aldehyde":
                count=sum(g["count"] for g in groups if g["kind"] in ("aryl_aldehyde","alkyl_aldehyde"))
            if count:
                roles.append({"linkage":linkage,"role":role,"reactive_sites":count,"compatible":True,
                              "reason":"Recognized functional group; crystallization and reaction conditions are not predicted."})
    shapes=spatial.site_shapes(mol)
    if any(g['kind']=='tp_trialdehyde' for g in groups):
        shapes['tp_trialdehyde']='triangular'
    kagome=kagome_roles(mol,patterns)
    if 'aryl_aldehyde' in kagome and not any(g['kind']=='alkyl_aldehyde' for g in groups):
        kagome['carbon_aldehyde']=kagome['aryl_aldehyde']
    for role in ('carbon_aldehyde','primary_alkylamine'):
        if shapes.get(role)=='linear':
            kagome[role]='spacer'
    tetra='tetrahedral' in shapes.values() and mol.HasSubstructMatch(Chem.MolFromSmarts("[C,Si;X4;!R](-[c])(-[c])(-[c])-[c]"))
    from .cofdraw_mixed import ETTA, linear_dialdehyde
    return {"attachment_directions":directions,"site_shapes":shapes,"functionality":groups,"reactive_groups":groups,"compatible_roles":roles,"kagome_roles":kagome,
            "kgm_etta":any(g["kind"]=="primary_arylamine" and g["count"]==4 for g in groups) and Chem.MolToSmiles(mol)==ETTA,
            "kgm_dialdehyde":any(g["kind"]=="aryl_aldehyde" and g["count"]==2 for g in groups) and linear_dialdehyde(smiles),
            "geometry_hint":"tetrahedral" if tetra else "tetrahedral_cage" if "tetrahedral" in shapes.values() else "unspecified",
            "geometry_note":("Four independently attached aryl arms on a tetrahedral center." if tetra else
                             "Four substituted adamantane bridgeheads form a tetrahedral node." if "tetrahedral" in shapes.values() else
                             "Six reactive arms occupy the three blades of a triptycene prism." if "prismatic" in shapes.values() else
                             "The nonplanar core needs an independently chosen connectivity." if {"spiro","spatial_unspecified"}&set(shapes.values()) else None)}

def library_topologies(catalog):
    from .cofdraw_starters import topology_catalog
    rows = topology_catalog(catalog.get("topology_catalog", []))
    for row in rows:
        if row["id"] in ("dia", "pts"):
            key = "dia_4_2" if row["id"] == "dia" else "pts_4_4"
            row.update(compatible_template_ids=[key], construction_supported=True, chemical_template_supported=True,
                       runtime_scope="Finite molecular cage or fragment from the selected monomers.")
    return rows


def library_linkages(catalog):
    result={}
    for row in catalog.get("linkage_catalog",[]):
        result[row["id"]]={**row,"supported":False,"roles":[],"role_labels":[],"category":"Other reported linkages",
                           "limits":row.get("qualification") or "Catalogued chemistry; no automatic product transform implemented."}
    for key,rule in LINKAGE_LIBRARY.items():
        result[key]={**result.get(key,{}),**rule,"label":LINKAGE_LABELS[key],"supported":True,
                     "automatic_prepare_supported":True,"runtime_status":"bounded_monomer_transform",
                     "category":"Boron linkages" if key=="boronate_ester" else "Carbon–carbon linkages" if key in ("cyano_vinylene_sp2_carbon","vinylene_methyltriazine") else "Nitrogen linkages","role_labels":[ROLE_LABELS[r] for r in rule["roles"]]}
    return result

# Library facets describe the molecule, independently of its construction roles.
_GROUP_PATTERNS = {label: Chem.MolFromSmarts(pattern) for label, pattern in {
    "Aldehyde": "[CX3H1](=O)[#6]",
    "Primary amine": "[NX3;H2;+0]-[#6;!$([C]=[O,S,N])]",
    "Secondary amine": "[NX3;H1;+0](-[#6;!$([C]=[O,S,N])])-[#6;!$([C]=[O,S,N])]",
    "Azide": "[#6]-[$(N=[N+]=[N-]),$([N-]-[N+]#N)]",
    "Ether": "[#6;!$([C]=[O,S])]-[OX2]-[#6;!$([C]=[O,S])]",
    "Thiol": "[#6]-[SH1]",
    "Hydroxylamine": "[NX3H1]([OH1])-[#6]",
    "Acyl hydrazide": "[NH2]-[NH]-C(=O)",
    "Hydrazine": "[NH2]-[NH2]",
    "Boronic acid": "[#6]-B([OH1])[OH1]",
    "Nitrile": "[#6]#N",
    "Phenol": "c-[OH1]",
    "Catechol": "[OH1]-c:c-[OH1]",
    "Alcohol": "[CX4]-[OH1]",
    "Methoxy": "[CH3]-[O;X2]-[#6;!$([C]=O)]",
    "Carboxylic acid / carboxylate": "[CX3](=O)[O;H1,-1]",
    "Ketone": "[#6][CX3](=O)[#6]",
    "Ester": "[#6][CX3](=O)O[#6;!$([C]=O)]",
    "Amide": "[#6][CX3](=O)[NX3]",
    "Anhydride": "C(=O)OC(=O)",
    "Acyl chloride": "C(=O)Cl",
    "Nitro": "[N+](=O)[O-]",
    "Sulfonic acid / sulfonate": "[#6]S(=O)(=O)[O;H1,-1]",
    "Sulfone": "[#6]S(=O)(=O)[#6]",
    "Alkyne": "[#6]#[#6]",
    "Organic halide": "[#6;!$([C]=O)]-[F,Cl,Br,I]",
}.items()}




def enable_reviewed_feed(row):
    """Re-evaluate an intact archived feed when its reaction family becomes drawable."""
    if not str(row.get('catalog_use','')).startswith('reference_only') or not row.get('source_dois') or row.get('annotation_only'):
        return row
    donor_roles={'aryl_acetonitrile','methyltriazine','aryl_carboxyl','aryl_cyclic_anhydride','aryl_boronic_acid','aryl_vicinal_diol','primary_alkylamine'}
    inspected=_inspect_monomer(row['smiles'])
    if any(g['kind']=='alkyl_aldehyde' for g in inspected['reactive_groups']):
        donor_roles.add('carbon_aldehyde')
    for r in inspected['compatible_roles']:
        if r['role'] not in donor_roles or r['reactive_sites']<2:
            continue
        try:
            donor=linkage_fragment(row['smiles'],r['linkage'],r['role'])
        except (ValueError,RuntimeError):
            continue
        if donor['connectivity']==r['reactive_sites']:
            return {**row,'catalog_use':'framework_feed','automatic_prepare_supported':True}
    return row

def decorated_monomers(rows):
    from .cofdraw_self import self_kind, can_tile_rhombus
    from .cofdraw_mixed import linear_dialdehyde
    output=[]
    for row in rows:
        try:
            metadata=inspect_monomer(row["smiles"])
            mol = Chem.MolFromSmiles(row['smiles'])
            complete = len(Chem.GetMolFrags(mol)) == 1 and all(
                a.GetAtomicNum() and not a.GetAtomMapNum() and not a.GetNumRadicalElectrons() for a in mol.GetAtoms())
            eligible = complete and not str(row.get('catalog_use','')).startswith('reference_only') and row.get('automatic_prepare_supported') is not False and not row.get('annotation_only')
            candidates = []
            if eligible and not row.get('self_activation'):
                for entry in metadata['compatible_roles']:
                    if entry['reactive_sites'] < 2:
                        continue
                    try:
                        fragment = linkage_fragment(row['smiles'], entry['linkage'], entry['role'])
                    except (ValueError, RuntimeError):
                        continue
                    if fragment['connectivity'] == entry['reactive_sites']:
                        fragment_atoms = Chem.MolFromSmiles(fragment['fragment_smiles']).GetNumAtoms()
                        candidates.append({**entry,'fragment_atoms':fragment_atoms})
            groups = {g['kind']:g['count'] for g in metadata['functionality']}
            aldehydes, amines = groups.get('aryl_aldehyde',0), groups.get('primary_arylamine',0)
            self_imine = complete and row.get('catalog_use') != 'synthetic_intermediate' and (bool(row.get('self_activation')) or aldehydes == amines and amines in (1,2))
            roles = {k:next((e['reactive_sites'] for e in candidates if e['linkage']=='imine' and e['role']==k),0)
                     for k in ('aryl_aldehyde','primary_arylamine')}
            kind = self_kind(row,mol,self_imine) if complete else None
            pores = ([{'id':'hexagon','label':'Kagome · hexagonal pore'}, {'id':'triangle','label':'Kagome · triangular pore'}]
                     if kind=='self_kagome' else [])
            if pores and can_tile_rhombus(canonical(row['smiles'])):
                pores.append({'id':'rhombus','label':'Square lattice · rhombic pore'})
            functional_groups = sorted(label for label, pattern in _GROUP_PATTERNS.items() if mol.HasSubstructMatch(pattern))
            output.append({**row,**metadata,"functional_groups":functional_groups,
                           "linkage_functional_groups":linkage_functional_groups(row, mol, candidates, kind),
                           'binary_eligible':bool(candidates),'binary_roles':candidates,
                           'canonical_smiles':canonical(row['smiles']), 'assembly_roles':{**roles,'self_imine':self_imine,
                           'self_kind':kind,'self_pores':pores,
                           'kgm_linear_dialdehyde':roles['aryl_aldehyde']==2 and linear_dialdehyde(row['smiles'])},
                           "catalog_reactive_groups":row.get("reactive_groups", []),
                           "abbreviation":row.get("abbreviation") or row.get("id", ""),
                           "full_name":row.get("full_name") or row.get("label", ""),
                           "geometry_hint":metadata["geometry_hint"] if metadata["geometry_hint"].startswith("tetrahedral")
                           else row.get("geometry_hint") or metadata["geometry_hint"]})
        except (ValueError,KeyError):
            continue
    return output

def best_port_pair(component, angle):
    roles = tuple((s["port"],s["site_role"]) for s in component.get("reactive_sites",[]) if "site_role" in s)
    error,incoming,outgoing,actual = _port_fit(component["fragment_smiles"], angle, roles)
    return incoming,outgoing,{"target_angle":angle,"depicted_site_angle":actual,"absolute_angle_difference":error,
                              "basis":"Closest pair in a monomer 2D drawing; inspect incidence before assigning source identity."}

@lru_cache(maxsize=1024)
def _port_fit(fragment_smiles, angle, site_roles):
    mol=Chem.MolFromSmiles(fragment_smiles); compute_2d_coords(mol)
    xy=np.asarray(mol.GetConformer().GetPositions())[:,:2]
    center=xy[[a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum()]].mean(0)
    ports=_ports(mol); candidates=[]
    roles = dict(site_roles)
    for incoming,(i,_) in ports.items():
        for outgoing,(j,_) in ports.items():
            if incoming==outgoing or roles and (roles[incoming],roles[outgoing])!=("amine","aldehyde"): continue
            left,right=xy[i]-center,xy[j]-center
            scale=np.linalg.norm(left)*np.linalg.norm(right)
            if scale<1e-10: continue
            actual=math.degrees(math.acos(float(np.clip(np.dot(left,right)/scale,-1,1))))
            candidates.append((abs(actual-angle),incoming,outgoing,actual))
    if not candidates: raise ValueError("Cannot orient the selected reactive sites")
    return min(candidates)

def design_signature(recipe):
    value={"components":{k:{f:c.get(f) for f in ("fragment_smiles","monomer_smiles","role")} for k,c in recipe["components"].items()},
           "instances":[{k:x.get(k) for k in ("id","component","incoming_port","outgoing_port")} for x in recipe["instances"]],
           "connections":recipe["connections"],"linkage":recipe.get("linkage"),"layout_kind":recipe.get("layout_kind")}
    return hashlib.sha256(json.dumps(value,sort_keys=True,separators=(',',':')).encode()).hexdigest()

def refresh_design_identity(recipe):
    if recipe.get('assembly'):
        from .cofdraw_cofex import recipe_record
        record=recipe_record(recipe) if recipe.get('design_template_signature')==design_signature(recipe) else None
        recipe.update(selected_template=None,design_topology=None,design_cofex_notation=None,
                      design_identity_qualification='Local assembly. Periodic connectivity and catenation are unspecified.')
        if record:
            recipe['design_cofex_notation']=record['text']
        return
    selected=recipe.get("selected_template") or {}
    matched=recipe.get("design_template_signature")==design_signature(recipe)
    known=next((t for t in TEMPLATES if t["id"]==selected.get("id")),None) if matched else None
    if known and known["id"] in ("kgm_4_2", "dia_4_2", "pts_4_4", "mtf_4_3"):
        known = {**known, "instances": len(recipe["instances"])}
        if recipe.get("layout_kind") == "net_patch" and known["id"] in spatial.SPATIAL:
            known["topology_label"] = known["parent_net"]
    net=known["topology_label"] if known else None
    parent=net or (known or {}).get('parent_net') or 'unk'
    components=recipe["components"]
    smiles=sorted({canonical(c["monomer_smiles"]) for c in components.values() if c.get("monomer_smiles")})
    name=re.sub(r"[^A-Za-z0-9_+-]","","-".join(c.get("label",k) for k,c in components.items()))[:80] or "user-design"
    recipe.update(design_topology=net,selected_template=known,
                  design_cofex_notation=".".join(smiles)+f" COFex.{parent}.cat0;{name}" if smiles else None,
                  design_identity_qualification="Proposed design notation. Net is the chosen unchanged template, cat0 is assumed; neither is a literature identity or a CIF.")

def net_geometry_error(template, counts, inspection, roles):
    """Check candidate geometry after the linkage roles and site counts match."""
    key = template["id"]
    shapes = [row.get("site_shapes", {}).get(role) for row, role in zip(inspection, roles)]
    directions = [row.get("attachment_directions", {}).get(role) for row, role in zip(inspection, roles)]
    if key == "pts_4_4":
        if set(directions) != {"tetrahedral", "square_planar"} and not spatial.matches(template, shapes):
            return "The pts candidate needs one tetrahedral and one square-planar set of four reactive sites."
    elif key == "dia_4_2":
        node = counts.index(4)
        if not spatial.matches(template, shapes) and not (
                directions[1-node] == "linear" and
                (directions[node] == "tetrahedral" or shapes[node] == "tetrahedral")):
            return "The diamond candidate requires a tetrahedral node and a linear ditopic spacer."
    elif key == "kgm_4_2":
        if not kagome_match(inspection, roles):
            node = counts.index(4)
            if (directions[1-node] != "linear" or shapes[node] in ("tetrahedral","spiro","spatial_unspecified")
                    or directions[node] == "tetrahedral" or "unresolved_alkyl" in shapes):
                return "The Kagome drawing requires a planar four-arm node such as ETTA and a linear ditopic spacer."
    elif not spatial.matches(template, shapes):
        return "The attachment geometry does not fit this template."
    elif key == "branched_4_3" and inspection[counts.index(4)].get("geometry_hint") != "tetrahedral":
        return "The branched design requires a tetrahedral tetratopic monomer and a complementary tritopic monomer."
    elif template.get("topology_label") and ("tetrahedral" in shapes or "tetrahedral" in directions):
        return "A tetrahedral node requires 3D incidence; select a compatible spatial template."
    return None


def prepare_template(monomer_a,monomer_b,linkage,label_a,label_b,instances,template_id=None):
    if linkage not in LINKAGE_LIBRARY: raise ValueError("This linkage is searchable but has no automatic construction rule")
    if linkage in ALKYL_VARIANTS:
        donors={row['kind']:row['count'] for row in _inspect_monomer(monomer_b)['reactive_groups']}
        if donors.get('primary_alkylamine',0)>=2 and donors.get('primary_arylamine',0)<2:
            linkage=ALKYL_VARIANTS[linkage]
    rule=LINKAGE_LIBRARY[linkage]
    components={key:{**linkage_fragment(smiles,linkage,role),"label":label}
                for key,smiles,role,label in zip(("A","B"),(monomer_a,monomer_b),rule["roles"],(label_a,label_b))}
    return _assemble_template(components, linkage, instances, template_id)


MTF_PATCH_ATOM_LIMIT = 1500


def _assemble_template(components, linkage, instances, template_id, *, legacy_mtf=False):
    rule=LINKAGE_LIBRARY[linkage]
    monomer_a, monomer_b=(components[k]['monomer_smiles'] for k in ('A','B'))
    label_a, label_b=(components[k]['label'] for k in ('A','B'))
    counts=[components[k]["connectivity"] for k in ("A","B")]
    if min(counts)<2: raise ValueError("Each component needs at least two compatible reactive sites")
    template=next((dict(t) for t in TEMPLATES if t["id"]==template_id),None)
    if template_id and not template: raise ValueError("Unknown drawing template")
    if template and linkage not in template["linkages"]:
        raise ValueError("This connectivity template has no construction rule for the selected linkage.")
    if template and template["allow_pairs"] and counts not in template["allow_pairs"]:
        raise ValueError(f"{template['label']} requires {template['allow_pairs']}; these monomers have {counts} reactive sites")
    inspection=[inspect_monomer(s) for s in (monomer_a,monomer_b)]
    if template:
        reason = net_geometry_error(template, counts, inspection, rule["roles"])
        if reason:
            raise ValueError(reason)

    if template_id == 'mtf_4_3':
        from .cofdraw_branched import branch_patch
        from .cofdraw_mixed import circuit_geometry
        # Every size has equal numbers of the two alternating monomers.
        # The complete finite patch has four groups around one large aperture.
        size = instances if instances in (4, 8) else 32
        pair_atoms = sum(Chem.MolFromSmiles(c['fragment_smiles']).GetNumAtoms() for c in components.values())
        size_note = None
        if size == 32 and 16 * pair_atoms > MTF_PATCH_ATOM_LIMIT:
            size = 8 if 4 * pair_atoms <= 1000 else 4
            size_note = 'Showing one square because the complete mtf patch exceeds the drawing atom limit.'
        face_count = {4: 0, 8: 1, 32: 17}[size]
        if legacy_mtf:
            mol, summary = branch_patch(components, 1, 'legacy_mtf')
        else:
            mol, summary = branch_patch(components, face_count, template_id)
        geometry = circuit_geometry(components, mol, summary, summary['order_assumption'])
        template['instances'] = size
        recipe = dict(components=components, linkage=linkage, selected_template=template,
                      layout_kind='net_patch' if size == 32 else 'closed_pore', designed_variant=True,
                      compatibility=dict(reactive_site_counts=counts, roles=rule['roles'], monomers=inspection),
                      **geometry)
        if size == 32:
            recipe['repeat_faces'] = summary['pore_faces']
        if size_note:
            recipe['drawing_extent_note'] = size_note
        recipe['design_template_signature'] = design_signature(recipe)
        refresh_design_identity(recipe)
        return recipe

    # The molecular cage and star have their own fixed sizes. Smaller saved
    # requests use the corresponding single-circuit / single-face depiction.
    full_patch = template_id in ("dia_4_2", "pts_4_4", "kgm_4_2") and instances == template["instances"]
    directions = [row["attachment_directions"].get(role) for row, role in zip(inspection,rule["roles"])]
    if template_id == "dia_4_2" and full_patch:
        full_patch = "linear" in directions
    if template_id == "pts_4_4" and full_patch:
        full_patch = set(directions) == {"tetrahedral", "square_planar"}
    if full_patch:
        from .cofdraw_net_patch import net_patch
        if template_id == "pts_4_4":
            from .cofdraw_net_patch import pts_patch
            tetra = "A" if inspection[0]["attachment_directions"][rule["roles"][0]] == "tetrahedral" else "B"
            items, connections, atoms = pts_patch(components, tetra, rule["join_bond_order"])
        else:
            node = "A" if counts[0] == 4 else "B"
            items, connections, atoms = net_patch(components, node, template_id == "dia_4_2", rule["join_bond_order"])
        recipe = dict(components=components, instances=items, connections=connections, atom_coordinates=atoms,
                      linkage=linkage, layout_kind="net_patch", selected_template=template, designed_variant=True,
                      compatibility=dict(reactive_site_counts=counts, roles=rule["roles"], monomers=inspection))
        recipe["design_template_signature"] = design_signature(recipe)
        refresh_design_identity(recipe)
        return recipe
    if template_id == "kgm_4_2":
        if not kagome_match(inspection, rule["roles"]):
            raise ValueError("This single-pore drawing requires a four-arm tetraarylethene such as ETTA and a linear ditopic spacer.")
        template["instances"] = 12
    if template_id in ("dia_4_2", "pts_4_4"):
        template["instances"] = spatial.SPATIAL[template_id]["instances"]
    positions=[]; kinds=[]; angles=[]
    n=template["instances"] if template else instances
    open_chain=template_id=="chain_2_2"
    if open_chain:
        positions=[[12*(i-(n-1)/2),0] for i in range(n)]; kinds=["A" if i%2==0 else "B" for i in range(n)]; angles=[180]*n
    elif template_id in ("hcb_3_2","sql_4_2","hxl_6_2","kgm_4_2"):
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
    geometry={}
    if template_id=="branched_4_3" or template_id in spatial.SPATIAL or template_id in planar.PLANAR:
        from .cofdraw_branched import branch_patch
        from .cofdraw_mixed import circuit_geometry
        mol,summary=branch_patch(components,1,template_id)
        geometry=circuit_geometry(components,mol,summary,summary['order_assumption'])
        selected=[dict(instance=item['id'],incoming_port=item['incoming_port'],outgoing_port=item['outgoing_port'],
                       basis=summary['order_assumption']) for item in geometry['instances']]
    if template_id=="kgm_4_2":
        from .cofdraw_mixed import pore_face_geometry
        geometry=pore_face_geometry(components,'BCBCBC','kgm')
        selected=[{'instance':item['id'],'incoming_port':item['incoming_port'],'outgoing_port':item['outgoing_port'],
                   'basis':'Kagome edge incidence'} for item in geometry['instances']]
    net=template["topology_label"] if template else None
    design_name=re.sub(r"[^A-Za-z0-9_+-]","",f"{label_a}-{label_b}")[:80] or "user-design"
    parent=net or (template or {}).get('parent_net') or 'unk'
    notation=".".join(sorted({canonical(monomer_a),canonical(monomer_b)}))+f" COFex.{parent}.cat0;{design_name}"
    return {"components":components,"instances":items,"connections":connections,"linkage":linkage,
            "layout_kind":"open_chain" if open_chain else "closed_pore", "selected_template":template,
            "design_topology":net,"design_cofex_notation":notation,"design_catenation":"cat0",
            "design_identity_qualification":"Proposed design notation. Net is the chosen template, cat0 is assumed; neither is experimental identity or a CIF.",
            "compatibility":{"reactive_site_counts":counts,"roles":rule["roles"],"port_selection":selected,"monomers":inspection},**geometry}
