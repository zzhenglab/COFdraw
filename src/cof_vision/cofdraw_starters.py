"""Cited precursor/net assignments and design-only alternatives for the picker."""
from copy import deepcopy

# These are inputs to the molecular builder. A generated single net is a
# separate design from a published interpenetrated crystal with the same feeds.
SINGLE_CRYSTAL_DOI = "10.1126/science.aat7679"
GENERATED = [
    ("starter-cof-303", "COF-303", "tfpm", "pa", "dia", SINGLE_CRYSTAL_DOI),
    ("starter-lzu-111", "LZU-111", "tapm", "tfs_cho4", "lon", SINGLE_CRYSTAL_DOI),
    ("starter-carbon-dia", None, "tfpm", "tapm", "dia", None),
    ("starter-carbon-lon", None, "tfpm", "tapm", "lon", None),
    ("starter-tapm-ta-lon", None, "tapm", "ta", "lon", None),
    ("starter-tapm-bpda-lon", None, "tapm", "bpda_cho", "lon", None),
    ("starter-tfpm-pa-lon", None, "tfpm", "pa", "lon", None),
]

# Original reported or proposed nets for specific bundled crystal models.
# Crystallographic assignment stays with this phase when a user edits a feed.
PHASES = {
    "cof-300-c5": ("dia", ["tapm","ta"], "generated"),
    "curated-13040N3": ("dia", ["tapm","bpda_cho"], "generated"),
    # El-Kaderi et al., Science 316 (2007), 268–272, Fig. 1 and Fig. 4.
    "curated-07010N3": ("ctn", ["tbpm_boronic"], "reference"),
    "curated-07011N3": ("ctn", ["tbps_boronic"], "reference"),
    "curated-07012N3": ("ctn", ["tbps_boronic","hhtp"], "reference"),
    "curated-07013N3": ("bor", ["tbpm_boronic","hhtp"], "reference"),
    # Lin et al., JACS (2016), DOI 10.1021/jacs.6b00652; Lin et al.,
    # JACS (2017), DOI 10.1021/jacs.7b04141. Both nets are structural models.
    "curated-16250N3": ("pts", ["tapm","tfppy"], "reference"),
    "curated-17163N3": ("pts", ["tapm","por_cho"], "reference"),
    "curated-19240N3": ("pts", ["tapm","tpb_h"], "reference"),
    "curated-19241N3": ("pts", ["tapm","tpb_me"], "reference"),
    "curated-19242N3": ("pts", ["tapm","tpb_f"], "reference"),
    "curated-21270N3": ("pts", ["tapm","tpb_ome"], "reference"),
    "curated-21271N3": ("ljh", ["tapm","tpb_ph"], "reference"),
}

def tpb_precursors():
    # The groups occupy the two remaining positions on the central benzene.
    formylphenyl="(-c2ccc(C=O)cc2)"
    substitutions=[("H","", "hydrogen"),("Me","(C)","methyl"),("F","(F)","fluoro"),
                   ("OMe","(OC)","methoxy"),("Ph","(-c2ccccc2)","phenyl")]
    rows=[]
    for code,branch,name in substitutions:
        label="TPB-"+code
        doi="10.1021/jacs.1c03042" if code in ("OMe","Ph") else "10.1002/anie.201905591"
        smiles="c1"+branch+"c"+formylphenyl+"c"+formylphenyl+"c"+branch+"c"+formylphenyl+"c1-c1ccc(C=O)cc1"
        rows.append(dict(id="tpb_"+code.lower(),label=label,abbreviation=label,smiles=smiles,
            full_name="1,2,4,5-Tetrakis(4-formylphenyl)benzene" if code=="H" else
            "3,6-Di"+name+"-1,2,4,5-tetrakis(4-formylphenyl)benzene",
            geometry_hint="planar_four_arm",source_dois=[doi],literature_names=["3D-TPB-COF-"+code]))
    return rows


REFERENCE_NETWORKS = {
    "curated-07010N3":1, "curated-07011N3":1, "curated-07012N3":1, "curated-07013N3":1,
    "curated-16250N3":2, "curated-17163N3":2,
    "curated-19240N3":5, "curated-19241N3":5, "curated-19242N3":5,
    "curated-21270N3":5, "curated-21271N3":1,
}

def crystal_starters(catalog):
    rows = [*catalog.get('reported_models',[]), *catalog.get('curated_models',[]),
            *catalog.get('templates',[]), *catalog.get('linker_models',[])]
    for row in rows:
        phase = PHASES.get(row['id'])
        if phase:
            net, ids, kind = phase
            row['starter'] = dict(net=net,monomer_ids=ids,kind=kind,common_name=row['model_name'],
                                  doi=row['citation_doi'],linkage='boroxine' if len(ids)==1 else
                                  'boronate_ester' if net in ('ctn','bor') else 'imine')
            if row['id'] in REFERENCE_NETWORKS:
                note = {"curated-16250N3":"The paper proposes two interpenetrating nets. The bundled coordinate graph resolves one network.",
                        "curated-21271N3":"One self-penetrated ljh net; no second independent net."}.get(row['id'])
                row['cofex_phase']=dict(monomer_ids=ids,net=net,network_count=REFERENCE_NETWORKS[row['id']],
                                        linkage=row['starter']['linkage'],qualification=note)
        elif row.get('linker_label') or row['id']=='dia-imine':
            row['starter'] = dict(net='dia',monomer_ids=row.get('monomer_ids') or ['tapm','ta'],kind='generated')
    return [dict(id=key,model_name=name or '',network_dimension=3,
                 starter=dict(net=net,monomer_ids=[a,b],kind='generated',common_name=name,doi=doi))
            for key,name,a,b,net,doi in GENERATED]


def topology_catalog(rows):
    rows = deepcopy(rows)
    if rows:
        for code,stack in (('4h','ABAC'),('6h','ABCACB')):
            rows.append(dict(id=code,label=f'{code.upper()} tetrahedral polytype', dimension=3,
                kind='polytype',connectivity=[[4,2],[4,4]],construction_supported=True,
                description=f'Build a tetrahedral network with the {stack} stacking sequence.',
                qualification='Ideal tetrahedral design. The sequence specifies the junction net across the cell.',
                source_refs=[dict(title='Tetrahedral stacking sequences',url='https://journals.iucr.org/j/issues/2013/01/00/kk5117/')]))
    capability = {
        'kgm':'Finite binary, self-condensed and mixed-linker Kagome drawings for supported precursor geometries.',
        'dia':'Editable periodic model for tetrahedral C/Si nodes with complementary ditopic spacers or a second tetrahedral node.',
        'lon':'Editable hexagonal-diamond model for tetrahedral C/Si nodes with complementary ditopic spacers or a second tetrahedral node.',
    }
    descriptions = {
        'kgm':('Kagome','ETTA-based sheets can contain both triangular and hexagonal pores.',[4,2],'10.1039/C6SC05673C'),
        'dia':('Diamond','COF-300 forms a diamond network with tetrahedral TAPM junctions and terephthalaldehyde spacers.',[4,2],'10.1021/ja8096256'),
        'lon':('Hexagonal diamond','Single-crystal diffraction identified a threefold interlocked lon-b network for LZU-111.',[4,4],SINGLE_CRYSTAL_DOI),
        'pts':('Platinum sulfide','The structural model for 3D-Py-COF pairs tetrahedral TAPM with tetratopic pyrene.',[4,4],'10.1021/jacs.6b00652'),
        'ctn':('Cubic carbon nitride','COF-102, COF-103 and COF-105 were modeled with tetrahedral and triangular junctions.',[4,3],'10.1126/science.1139915'),
        'bor':('Boracite','COF-108 has tetrahedral boronic-acid-derived sites and planar HHTP junctions.',[4,3],'10.1126/science.1139915'),
        'hea':('hea','3D-hea-COF couples a six-connected triptycene aldehyde with a four-connected tetrahedral amine.',[6,4],'10.1021/acs.chemmater.1c03156'),
    }
    if rows and not any(row['id']=='hea' for row in rows):
        rows.append(dict(id='hea',net_code='hea',dimension=3,kind='net',construction_supported=False))
    for row in rows:
        code = row['id']
        if code not in descriptions:
            continue
        name,description,counts,doi = descriptions[code]
        row.update(label=f'{code} · {name}' if code!='hea' else 'hea · 3D network', description=description,
                   connectivity=[counts], source_refs=[{'title':'Source paper','url':'https://doi.org/'+doi}],
                   qualification=description,limitations=[capability.get(code,'A literature assignment applies to the cited phase.')])
        if code in capability:
            row.update(construction_supported=True, construction_note=capability[code])
    return rows
