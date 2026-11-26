"""Compact COFex exchange for the precursor-based drawing and crystal builders."""
import hashlib
import itertools
import re
from functools import lru_cache
from urllib.parse import parse_qs, urlencode

from rdkit import Chem

from .cofdraw_chemistry import LINKAGE_LABELS, TEMPLATES, decorated_monomers, inspect_monomer, kagome_match, net_geometry_error
from .framework_drawing import LINKAGE_LIBRARY, canonical, matches_linkage_family


# Earlier exchange records for these templates describe a single circuit.
_CIRCUIT_SIZE = {'dia_4_2': 12, 'pts_4_4': 4, 'kgm_4_2': 12, 'mtf_4_3': 4}


def design_record(feeds, mode, linkage, net=None, count=1, **settings):
    """Keep canonical COFex intact and carry drawing choices in its optional comment."""
    smiles = [canonical(row['smiles']) for row in feeds]
    sorted_smiles = sorted(smiles)
    options = [('mode', mode), ('link', linkage), ('order', ','.join(str(sorted_smiles.index(s)) for s in smiles))]
    if mode=='crystal' and net in ('4h','6h'):
        settings['polytype']=net
        net='unk'
    options.extend((k, v) for k, v in settings.items() if v)
    options.extend(('label', str(row.get('label') or chr(65+i))[:100]) for i, row in enumerate(feeds))
    compact = '.'.join(sorted_smiles) + f" COFex.{net or 'unk'}.cat{0 if count == 1 else count}"
    return {'text': compact + ';COFdraw/1?' + urlencode(options), 'compact': compact,
            'scope': 'design', 'net': net or 'unk', 'network_count': count, 'linkage': linkage}


def phase_record(item, catalog_rows):
    phase = item.get('cofex_phase')
    if not phase:
        return None
    by_id = {row['id']:row['smiles'] for row in catalog_rows}
    if any(key not in by_id for key in phase['monomer_ids']):
        return None
    precursors=[canonical(by_id[key]) for key in phase['monomer_ids']]
    joined = '.'.join(sorted(precursors))
    count=phase['network_count']
    text = joined + f" COFex.{phase['net']}.cat{0 if count==1 else count};{item['model_name']}"
    return {'text':text, 'compact':text, 'scope':'source', 'doi':item['citation_doi'],
            'net':phase['net'], 'network_count':count, 'linkage':phase['linkage'], 'precursor_smiles':precursors,
            'qualification':phase.get('qualification')}


def recipe_record(recipe, metadata=None, source=False):
    if source:
        notation = (metadata or {}).get('cofex_notation')
        if notation or recipe.get('cofex_id'):
            return {'text': notation or recipe['cofex_id'], 'compact': notation or recipe['cofex_id'],
                    'scope': 'source', 'id': recipe.get('cofex_id')}
        return None
    assembly = recipe.get('assembly')
    if assembly:
        feeds = assembly.get('precursors')
        if not feeds:
            return None
        kind = assembly.get('kind')
        net = 'kgm' if kind == 'self_kagome' or assembly.get('pore_layout') == 'kgm' else None
        if kind in ('self_boroxine','self_triazine'):
            net='hcb'
        elif kind=='self_imine' and assembly.get('pore_view')=='rhombus':
            net='sql'
        return design_record(feeds, 'assembly', recipe['linkage'], net,
                             sequence=assembly.get('linker_sequence'), layout=assembly.get('pore_layout'),
                             face=assembly.get('pore_view') if kind in ('self_imine','self_kagome') else None)
    template = recipe.get('selected_template') or {}
    if not template or not all(k in recipe.get('components', {}) for k in ('A','B')):
        # Unregistered explicit graphs and edited connections have no matching builder recipe.
        return None
    feeds = [dict(smiles=c['monomer_smiles'], label=c.get('label')) for c in (recipe['components'][k] for k in ('A','B'))]
    view = 'molecular_patch' if recipe.get('layout_kind') == 'net_patch' else None
    if template['id'] == 'mtf_4_3' and len(recipe.get('instances', [])) == 8:
        view = 'large_square'
    return design_record(feeds, 'binary', recipe['linkage'], template.get('topology_label') or template.get('parent_net'),
                         template=template['id'], view=view)


@lru_cache(maxsize=3000)
def _catalog_smiles(smiles):
    try:
        return canonical(smiles)
    except (RuntimeError, ValueError):
        return ''


def parse_record(text, catalog_rows):
    match = re.fullmatch(r'(\S+)\s+COFex\.([A-Za-z0-9_-]{1,32})\.cat(0|[1-9][0-9]?)(?:;([^\r\n]*))?', text.strip())
    if not match:
        raise ValueError('Paste precursor SMILES followed by COFex.net.catN, or a registered source COFex ID.')
    fragments, net, cat, comment = match.groups()
    net = net.lower()
    strings = fragments.split('.')
    if len(strings) > 3:
        raise ValueError('This builder supports one, two or three connected precursors.')
    smiles = []
    for s in strings:
        mol = Chem.MolFromSmiles(s) if len(s) <= 4000 else None
        if mol is None or not 1 <= mol.GetNumAtoms() <= 250 or any(not a.GetAtomicNum() or a.GetAtomMapNum() for a in mol.GetAtoms()):
            raise ValueError('Each precursor needs valid, fully specified SMILES with up to 250 atoms.')
        smiles.append(Chem.MolToSmiles(mol))
    if len(smiles) != len(set(smiles)):
        raise ValueError('Enter each distinct precursor once in COFex.')
    settings = {}
    if (comment or '').startswith('COFdraw/'):
        if not comment.startswith('COFdraw/1?'):
            raise ValueError('Unsupported COFdraw exchange version.')
        fields = parse_qs(comment.split('?',1)[1], strict_parsing=True, max_num_fields=12)
        allowed = {'mode','order','label','link','template','face','sequence','layout','polytype','view'}
        if (fields.keys() - allowed or any(len(v) != 1 for k,v in fields.items() if k != 'label')
                or not {'mode','order','link'}.issubset(fields)):
            raise ValueError('Invalid or repeated COFdraw exchange fields.')
        settings = {k:v[0] for k,v in fields.items() if k != 'label'}
        if not re.fullmatch(r'[0-2](?:,[0-2]){0,2}', settings['order']):
            raise ValueError('Invalid precursor order.')
        order = [int(i) for i in settings['order'].split(',')]
        if sorted(order) != list(range(len(smiles))):
            raise ValueError('Precursor order disagrees with the SMILES count.')
        smiles = [smiles[i] for i in order]
        if 'label' in fields and (len(fields['label']) != len(smiles) or any(len(x)>100 or not x.isprintable() for x in fields['label'])):
            raise ValueError('Invalid precursor labels.')
        settings['labels'] = fields.get('label', [])
    relevant = {}
    wanted = set(smiles)
    for row in catalog_rows:
        s = row.get('canonical_smiles') or _catalog_smiles(row['smiles'])
        if s in wanted and s not in relevant:
            relevant[s] = row
    labels = settings.get('labels') or [relevant.get(s, {}).get('abbreviation') or chr(65+i) for i,s in enumerate(smiles)]
    if 'polytype' in settings:
        if net!='unk' or settings['polytype'] not in ('4h','6h'):
            raise ValueError('Invalid tetrahedral stacking sequence.')
        net=settings['polytype']
    mode = 'assembly' if len(smiles) != 2 else 'crystal' if net in ('dia','lon','4h','6h') else 'binary'
    if mode == 'crystal' and settings.get('mode') == 'binary' and settings.get('template'):
        mode = 'binary'
    if settings.get('mode', mode) != mode:
        raise ValueError('Builder mode disagrees with the COFex net or precursor count.')
    specific = set(settings) - {'mode','order','link','labels'}
    allowed = {'binary':{'template','view'},'assembly':{'face'} if len(smiles)==1 else {'sequence','layout'},'crystal':{'polytype'}}
    if specific - allowed[mode]:
        raise ValueError('The COFdraw exchange fields disagree with the precursor count or net.')
    if 'view' in settings:
        valid_patch = settings['view'] == 'molecular_patch' and settings.get('template') in _CIRCUIT_SIZE
        large_square = settings['view'] == 'large_square' and settings.get('template') == 'mtf_4_3'
        if not (valid_patch or large_square):
            raise ValueError('This drawing view does not match the selected template.')
    if mode == 'assembly':
        return _assembly_choices(smiles, labels, relevant, net, int(cat), settings)
    if int(cat) > (7 if mode == 'crystal' else 0):
        raise ValueError('This net and catenation count cannot be generated by the selected builder.')
    choices = []
    roles = [inspect_monomer(s) for s in smiles]
    for link, rule in LINKAGE_LIBRARY.items():
        if not matches_linkage_family(settings.get('link',link),link):
            continue
        orientation = next((order for order in ((0,1),(1,0)) if all(
            any(r['linkage']==link and r['role']==role and r['reactive_sites']>=2 for r in roles[i]['compatible_roles'])
            for i,role in zip(order,rule['roles']))), None)
        if orientation is None:
            continue
        if mode == 'crystal':
            if link in ('imine','secondary_amine','acylhydrazone') and any(x['geometry_hint']=='tetrahedral' for x in roles):
                from .cofdraw_dia import _components
                try:
                    _components(smiles,link)
                except ValueError:
                    pass
                else:
                    choices.append({'label': f'{LINKAGE_LABELS[link]} · {net} · {max(int(cat),1)} net(s)', 'crystal_project': {
                        'template_id':'custom-dia', 'monomers':[dict(smiles=s,label=l) for s,l in zip(smiles,labels)],
                        'linkage':link, 'topology':net, 'network_count':max(int(cat),1)}})
            if settings or cat != '0' or net != 'dia':
                continue
        counts = [next(r['reactive_sites'] for r in roles[i]['compatible_roles'] if r['linkage']==link and r['role']==role)
                  for i,role in zip(orientation,rule['roles'])]
        a,b = orientation
        for template in TEMPLATES:
            key = template['id']
            if link not in template['linkages'] or template['allow_pairs'] and counts not in template['allow_pairs']:
                continue
            ordered_roles = [roles[a], roles[b]]
            if net_geometry_error(template, counts, ordered_roles, rule['roles']):
                continue
            if settings.get('template') != key and template.get('topology_label') and any(roles[i]['site_shapes'].get(r) in ('spiro','spatial_unspecified','prismatic') for i,r in zip(orientation,rule['roles'])):
                continue
            instances = template['instances']
            if key in _CIRCUIT_SIZE and settings.get('template') == key and 'view' not in settings:
                instances = _CIRCUIT_SIZE[key]
                if key == 'kgm_4_2' and not kagome_match(ordered_roles, rule['roles']):
                    continue
            elif key == 'mtf_4_3' and settings.get('view') == 'large_square':
                instances = 8
            requested = settings.get('template')
            legacy_ring = requested == 'alternating_hexagon' and net == 'unk'
            if ((requested and requested != key and not legacy_ring)
                    or net != 'unk' and (template['topology_label'] or template.get('parent_net')) != net
                    or template['topology_label'] and any(x['geometry_hint']=='tetrahedral' for x in roles)):
                continue
            title = 'Restore original single ring' if legacy_ring and key == requested else template['label']
            choices.append({'label': f"{LINKAGE_LABELS[link]} · {title}", 'prepare_request': dict(
                monomer_a=smiles[a], monomer_b=smiles[b], label_a=labels[a], label_b=labels[b],
                linkage=link, template_id=key, instances=instances)})
    if not choices:
        raise ValueError('No supported drawing rule matches these reactive groups, net and linkage.')
    return choices


def _assembly_choices(smiles, labels, catalog, net, cat, settings):
    if cat or net not in ('unk','kgm','hcb','sql') or (net=='sql' and len(smiles)!=1):
        raise ValueError('This assembly cannot establish the requested net or catenation.')
    rows, custom = [], []
    for s,label in zip(smiles,labels):
        if s in catalog:
            row = catalog[s]
        else:
            key = 'custom-' + hashlib.sha256(s.encode()).hexdigest()[:24]
            row = {'id':key, 'label':label, 'abbreviation':label, 'smiles':s}
            custom.append(dict(smiles=s,label=label))
        rows.append(row)
    rows = decorated_monomers(rows)
    choices = []
    if len(rows) == 1:
        kind = rows[0]['assembly_roles'].get('self_kind')
        if not kind and not rows[0]['assembly_roles'].get('self_imine'):
            raise ValueError('This precursor has no supported self-condensation rule.')
        link = {'self_boroxine':'boroxine','self_triazine':'triazine'}.get(kind,'imine')
        if settings.get('link',link) != link or net=='kgm' and kind!='self_kagome' or net=='hcb' and kind not in ('self_boroxine','self_triazine'):
            raise ValueError('The self-condensation rule disagrees with the linkage or net.')
        faces=[p['id'] for p in rows[0]['assembly_roles']['self_pores']] if kind=='self_kagome' else [None]
        if net=='kgm':
            faces=[face for face in faces if face in ('hexagon','triangle')]
        if net=='sql':
            faces=[face for face in faces if face=='rhombus']
        for face in faces:
            if settings.get('face',face) == face:
                choices.append({'label':link.capitalize() + (' · '+face if face else ''),'assembly_request':dict(
                    monomer_ids=[rows[0]['id']], custom_precursors=custom, pore_view=face)})
    else:
        for a,b,c in itertools.permutations(range(3)):
            if settings and (a,b,c)!=(0,1,2) or not settings and b>c:
                continue
            role = next((role for role,other in (('primary_arylamine','aryl_aldehyde'),('aryl_aldehyde','primary_arylamine'))
                         if rows[a]['assembly_roles'][role] in (3,4) and all(rows[i]['assembly_roles'][other]==2 for i in (b,c))),None)
            if not role or settings.get('link','imine')!='imine' or net=='hcb':
                continue
            n = 6 if rows[a]['assembly_roles'][role]==3 else 4
            layout = settings.get('layout', 'kgm' if net=='kgm' else None)
            if layout not in (None,'kgm') or bool(layout)!=(net=='kgm'):
                continue
            if any(r['assembly_roles']['self_imine'] or r.get('geometry_hint')=='tetrahedral' for r in rows):
                continue
            edge_role='aryl_aldehyde' if role=='primary_arylamine' else 'primary_arylamine'
            if layout and (rows[a]['kagome_roles'].get(role)!='node' or not all(rows[i]['kagome_roles'].get(edge_role)=='spacer' for i in (b,c))):
                continue
            patterns = ('BBB','CCC','BCBCBC','CBCBCB') if layout else (('BC'*(n//2)), ('B'*(n//2)+'C'*(n//2)))
            for pattern in ([settings['sequence']] if 'sequence' in settings else patterns):
                if (layout and pattern not in patterns) or (not layout and (len(pattern)!=n or set(pattern)!={'B','C'})):
                    continue
                choices.append({'label':f"Mixed imine · {pattern}",'assembly_request':dict(
                    monomer_ids=[rows[i]['id'] for i in (a,b,c)], custom_precursors=custom, linker_sequence=pattern, pore_layout=layout)})
    if not choices:
        raise ValueError('No supported local pore arrangement matches this precursor set.')
    return choices
