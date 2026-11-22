"""Single-feed chemistry for finite COF drawings."""
from functools import lru_cache
import math

import numpy as np
from rdkit import Chem

from .cofdraw_depiction import compute_2d_coords
from .framework_drawing import VERSION, _ports, linkage_fragment


LI_DOI = '10.1021/jacs.9b03463'
CHEN_DOI = '10.1021/acs.chemmater.4c01292'
# Li et al., Figure 6a and Supporting Information, section 21. Each amine
# arm has one ethynylene spacer; each aldehyde arm is attached directly.
KAGOME_PRECURSORS = {
    'bfbaeb': ('BFBAEB', '1,4-Bis(4-formylphenyl)-2,5-bis[(4-aminophenyl)ethynyl]benzene',
               'Nc1ccc(C#Cc2cc(-c3ccc(C=O)cc3)c(C#Cc3ccc(N)cc3)cc2-c2ccc(C=O)cc2)cc1'),
    'bfbaepy': ('BFBAEPy', '1,6-Bis(4-formylphenyl)-3,8-bis[(4-aminophenyl)ethynyl]pyrene',
                'Nc1ccc(C#Cc2cc(-c3ccc(C=O)cc3)c3ccc4c(C#Cc5ccc(N)cc5)cc(-c5ccc(C=O)cc5)c5ccc2c3c45)cc1'),
}


def additional_precursors():
    return [dict(id=k, abbreviation=abbr, full_name=name, label=name, smiles=smi,
                 source_dois=[LI_DOI], literature_names=[abbr+'-COF'],
                 geometry_hint='planar_four_arm', catalog_use='self_condensing_precursor', automatic_prepare_supported=False,
                 source_evidence=[{'source_doi':LI_DOI, 'source_label':abbr, 'location':'Figure 6a, SI section 21'}])
            for k,(abbr,name,smi) in KAGOME_PRECURSORS.items()]


SELF_EXAMPLES = (
    ('cof-1', 'COF-1 · boroxine', 'bdba', 'Three boronic acid ends condense into each B₃O₃ ring.', '10.1126/science.1120411'),
    ('ctf-1', 'CTF-1 · triazine', 'cv_d4982b6c8defb0a2', 'Three nitrile groups form each triazine in the idealized sheet.', '10.1002/anie.200705710'),
    ('py-cof', 'Py-COF · rhombic imine', 'bfbapy', 'Four two-in-one pyrene precursors surround a single local pore.', LI_DOI),
    ('bfbaeb-cof', 'BFBAEB-COF · triangle and hexagon', 'bfbaeb', 'Self-condensed benzene nodes in the reported Kagome arrangement.', LI_DOI),
    ('bfbaepy-cof', 'BFBAEPy-COF · triangle and hexagon', 'bfbaepy', 'Self-condensed pyrene nodes in the reported Kagome arrangement.', LI_DOI),
    ('cof-953', 'COF-953 · Kagome', 'feed_a633d27d4963', 'Choose Hexagon + 6 triangles to see the surrounding pores.', '10.1002/adma.202308251'),
    ('cof-954', 'COF-954 · Kagome', 'feed_f955f8b23f7a', 'Self-condensing aminoaldehyde with vinylene spacers.', '10.1002/adma.202308251'),
    ('latent-imine', 'Latent aminoaldehyde · pyridyl imine', 'cv_e0efd03afb224b2a', 'Two dimethyl acetals unmask the aldehydes of this reported precursor. The pore drawn is a design.', CHEN_DOI),
    ('aba-chain', '4-Aminobenzaldehyde · linear chain', 'cv_4e35cbcca820334e', 'One aldehyde and one amine make an open-chain model with two continuation sites.', None),
)


def self_examples():
    return [dict(id=id_,label=label,monomer_ids=[monomer],linker_sequence='',description=description,
                 **({'doi':doi} if doi else {})) for id_,label,monomer,description,doi in SELF_EXAMPLES]


TRIMER = {
    'self_boroxine': ('[c]-[B;X3;+0]([O;H1;X2;+0])[O;H1;X2;+0]', '[*:1]B1OB([*:2])OB([*:3])O1', 'Boroxine'),
    'self_triazine': ('[c]-[C;X2;+0]#[N;X1;+0]', '[*:1]c1nc([*:2])nc([*:3])n1', 'Triazine'),
}
TRIMER_QUERY = {k:Chem.MolFromSmarts(rule[0]) for k,rule in TRIMER.items()}


def _spacer(mol, matches):
    rw = Chem.RWMol(mol)
    remove = []
    for idx, (anchor, *end) in enumerate(matches, 1):
        dummy = Chem.Atom(0)
        dummy.SetAtomMapNum(idx)
        rw.AddBond(anchor, rw.AddAtom(dummy), Chem.BondType.SINGLE)
        remove.extend(end)
    for idx in sorted(remove, reverse=True):
        rw.RemoveAtom(idx)
    spacer = rw.GetMol()
    Chem.SanitizeMol(spacer)
    return spacer


@lru_cache(maxsize=256)
def can_tile_kagome(smiles):
    from .cofdraw_self_kagome import _cell
    try:
        part=linkage_fragment(smiles,'imine','self_imine')
        roles=tuple(sorted((r['port'],r['site_role']) for r in part['reactive_sites']))
        _cell(part['fragment_smiles'],roles)
    except (ValueError,RuntimeError):
        return False
    return True


@lru_cache(maxsize=128)
def can_tile_rhombus(smiles):
    from .cofdraw_repeat import _fit_sheet
    part=linkage_fragment(smiles,'imine','self_imine')
    fragments=(('A',part['fragment_smiles'],part['connectivity']),)
    return part['connectivity']==4 and _fit_sheet('self_square',fragments,2) is not None


def rhombus_recipe(row,component):
    """Take one face from the same molecular SQL tiling used for larger extents."""
    from .cofdraw_repeat import _sheet
    recipe={'components':{'A':component},'connections':[{'bond_order':2}]}
    mol,summary=_sheet(recipe,'self_square',1,'pores')
    face,=summary['pore_faces']
    port_for_atom={a:n for n,(_,a) in _ports(Chem.MolFromSmiles(component['fragment_smiles'])).items()}
    positions=mol.GetConformer().GetPositions()[:,:2]
    atoms={p:[] for p in face}
    names={p:f'A{i+1}' for i,p in enumerate(face)}
    coordinates=[]
    for atom,xy in zip(mol.GetAtoms(),positions):
        owner=atom.GetProp('instance')
        atoms[owner].append(xy)
        coordinates.append(dict(instance=names[owner],fragment_atom_index=atom.GetIntProp('fragment_atom_index'),
                                xy=xy.tolist()))
    slots={}
    for bond in mol.GetBonds():
        a,b=bond.GetBeginAtom(),bond.GetEndAtom()
        p,q=a.GetProp('instance'),b.GetProp('instance')
        if p!=q:
            slots[p,q]=port_for_atom[a.GetIntProp('fragment_atom_index')]
            slots[q,p]=port_for_atom[b.GetIntProp('fragment_atom_index')]
    instances=[]
    connections=[]
    for i,p in enumerate(face):
        previous,q=face[i-1],face[(i+1)%len(face)]
        instances.append(dict(id=names[p],component='A',incoming_port=slots[p,previous],outgoing_port=slots[p,q],
                              layout_position=np.mean(atoms[p],axis=0).tolist()))
        connections.append(dict(a=dict(instance=names[p],port=slots[p,q]),b=dict(instance=names[q],port=slots[q,p]),
                                bond_order=2))
    return dict(components=recipe['components'],instances=instances,connections=connections,linkage='imine',
                layout_kind='closed_pore',atom_coordinates=coordinates,selected_template=None,design_topology=None,
                coordinate_provenance='generated from the fitted self-condensed square lattice',
                assembly=dict(monomer_ids=[row['id']],precursors=[feed_metadata(row)],kind='self_imine',linker_sequence='',
                              pore_view='rhombus',scope='One rhombic aperture of a proposed self-condensed sql sheet.'))


def self_kind(row, mol, self_imine):
    if self_imine:
        return 'self_kagome' if not row.get('self_activation') and can_tile_kagome(Chem.MolToSmiles(mol)) else 'self_imine'
    if row.get('catalog_use') == 'synthetic_intermediate' or row.get('annotation_only'):
        return None
    for kind,query in TRIMER_QUERY.items():
        ends = mol.GetSubstructMatches(query)
        if len(ends) != 2 or any(mol.GetAtomWithIdx(i).GetIsotope() for end in ends for i in end[1:]):
            continue
        spacer = _spacer(mol,ends)
        if len(Chem.GetMolFrags(spacer)) != 1:
            continue
        compute_2d_coords(spacer)
        xy = spacer.GetConformer().GetPositions()[:, :2]
        (d1,a1),(d2,a2) = _ports(spacer).values()
        tangents = np.array([xy[d1]-xy[a1], xy[d2]-xy[a2]])
        tangents /= np.linalg.norm(tangents,axis=1)[:,None]
        axis = xy[a1]-xy[a2]
        if np.linalg.norm(axis) < 1:
            continue
        axis /= np.linalg.norm(axis)
        if (np.dot(*tangents) <= -math.cos(math.radians(22.5))
                and min(tangents[0]@axis,-tangents[1]@axis) >= math.cos(math.radians(30))):
            return kind
    return None


def trimer_assembly(row, kind):
    from .cofdraw_chemistry import best_port_pair

    mol = Chem.MolFromSmiles(row['smiles'])
    ends = mol.GetSubstructMatches(TRIMER_QUERY[kind])
    spacer = _spacer(mol,ends)
    _,ring_smiles,ring_label = TRIMER[kind]
    linkage = kind.removeprefix('self_')
    feed = row['canonical_smiles']
    def fragment(fragment_smiles, label, role):
        ports = _ports(Chem.MolFromSmiles(fragment_smiles))
        return dict(monomer_smiles=feed,fragment_smiles=fragment_smiles,label=label,role=role,
                    feed_id=row['id'],linkage=linkage,connectivity=len(ports), library_version=VERSION,
                    reactive_sites=[dict(port=p,join_bond_order=1) for p in sorted(ports)])
    components = {'A':fragment(ring_smiles, ring_label+' junction', 'product_junction'),
                  'B':fragment(Chem.MolToSmiles(spacer),'Arylene spacer','product_spacer')}
    coord = np.arange(6)*math.pi/3 + math.pi/6
    corners = np.column_stack((np.cos(coord),np.sin(coord)))
    instances = []
    for i,p in enumerate(corners):
        for key,point,angle in [('A',p,120),('B',(p+corners[(i+1)%6])/2,180)]:
            incoming,outgoing,_=best_port_pair(components[key],angle)
            instances.append(dict(id=f'{key}{i+1}',component=key, incoming_port=incoming,
                                  outgoing_port=outgoing,layout_position=point.tolist()))
    edges = [dict(a=dict(instance=a['id'],port=a['outgoing_port']),b=dict(instance=b['id'],port=b['incoming_port']),
                  bond_order=1) for a,b in zip(instances,instances[1:]+instances[:1])]
    scope = ('Dehydration of boronic acid ends forms B₃O₃ junctions.' if kind=='self_boroxine'
             else 'Cyclotrimerization of nitrile ends forms C₃N₃ junctions.')
    return dict(components=components,instances=instances,connections=edges,linkage=linkage,layout_kind='closed_pore',
                selected_template=None,design_topology=None,
                assembly=dict(monomer_ids=[row['id']],linker_sequence='',kind=kind,pore_view='hexagon',
                              scope=scope+' One finite hexagonal pore from a single precursor.'))


def feed_metadata(row):
    return dict(id=row['id'],label=row.get('abbreviation') or row['label'],smiles=row['canonical_smiles'])
