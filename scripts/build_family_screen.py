"""Precompute reacted repeat compositions and a bounded sample of explicit AA cells."""
import argparse
import csv
import hashlib
import json
import math
import os
import sys
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor
from functools import lru_cache
from pathlib import Path

import numpy as np
from rdkit import Chem
from rdkit import rdBase

ROOT = Path(__file__).resolve().parents[1]
WEB = ROOT / 'src/cof_vision/web'
sys.path.insert(0,str(ROOT/'src'))
from cof_vision.framework_drawing import linkage_fragment, LINKAGE_LIBRARY
from cof_vision.cofdraw_chemistry import TEMPLATES, prepare_template
from cof_vision.cofdraw_periodic_build import build_periodic_layer, SUPPORTED_TEMPLATES
from cof_vision.cofdraw_porosity import geometric_porosity
from cof_vision.cofdraw_science import structure_exports, import_cif

RULES = {'imine':'imine','beta-ketoenamine':'beta_ketoenamine','hydrazone':'acylhydrazone',
         'keto-hydrazide':'keto_hydrazide','azine':'azine','nitrone':'nitrone'}
NETS = {tuple(pair):t['id'] for t in TEMPLATES for pair in t['allow_pairs']}
COLUMNS = ['pair','first_is_left_role','template','repeat_mass','nitrogen_wt','oxygen_wt','junction_density','aromatic_c_pct']


@lru_cache(None)
def fragment(smiles, rule, role):
    try:
        part=linkage_fragment(smiles,rule,role)
    except (ValueError, RuntimeError):
        return None
    if part['connectivity'] < 2:
        return None
    mol=Chem.AddHs(Chem.MolFromSmiles(part['fragment_smiles']))
    atoms=[a for a in mol.GetAtoms() if a.GetAtomicNum()]
    carbon=[a for a in atoms if a.GetAtomicNum()==6]
    return part['connectivity'], np.array([sum(a.GetMass() for a in atoms),
        sum(a.GetMass() for a in atoms if a.GetAtomicNum()==7),sum(a.GetMass() for a in atoms if a.GetAtomicNum()==8),
        len(carbon),sum(a.GetIsAromatic() for a in carbon)])


def pair_row(code,stride,blocks,rule):
    a,b=divmod(code,stride)
    left,right=LINKAGE_LIBRARY[rule]['roles']
    for first,second in ((a,b),(b,a)):
        fa=fragment(blocks[first][0],rule,left)
        fb=fragment(blocks[second][0],rule,right)
        if not fa or not fb or (fa[0],fb[0]) not in NETS:
            continue
        na,nb=fa[0],fb[0]
        divisor=math.gcd(na,nb)
        mass,n,o,c,aromatic=fa[1]*(nb//divisor)+fb[1]*(na//divisor)
        return [code,int(first==a),NETS[na,nb],round(mass,3),round(100*n/mass,3),
                round(100*o/mass,3),round(1000*na*nb/divisor/mass,4),round(100*aromatic/c,3) if c else None]
    raise ValueError(f'Cannot reproduce compatible pair {code} under {rule}')


def geometry(structure):
    g=geometric_porosity(structure,probe_radius=1.82,grid_spacing=.9,surface_samples=16)
    def rounded(key):
        return round(g[key],5) if g[key] is not None else None
    return dict(density=rounded('density_g_cm3'),void=round(100*g['void_fraction'],3),
        lcd=rounded('lcd_angstrom'),pld=rounded('pld_angstrom'),cell_volume=rounded('cell_volume_angstrom3'),
        pore_volume=rounded('geometric_void_volume_cm3_g'),free_volume=rounded('probe_center_free_volume_cm3_g'),
        connected_volume=rounded('transport_accessible_probe_center_volume_cm3_g'),
        area=rounded('transport_accessible_surface_m2_g'),local_area=rounded('local_probe_exposed_surface_m2_g'),
        lcd_bounds=g['lcd_bounds_angstrom'],pld_bounds=g['pld_bounds_angstrom'],
        channel=g['connectivity']['status'],grid=g['method']['grid_shape'])


def periodic(job):
    key,code,template,a,b=job
    try:
        design=prepare_template(a,b,RULES[key],'A','B',6,template)
        model,audit=build_periodic_layer(design)
        cif=structure_exports(model,periodic_bonds=audit['periodic_bonds'])['exports']['cif']
        return key,code,dict(**geometry(import_cif(cif)),cif=cif,sha256=hashlib.sha256(cif.encode()).hexdigest())
    except (ValueError, RuntimeError) as exc:
        return key,code,{'reason':str(exc)}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('inventory',type=Path)
    parser.add_argument('--workers',type=int,default=4)
    parser.add_argument('--attempts-per-family',type=int,default=36)
    parser.add_argument('--output',type=Path,default=WEB/'chemical-screen.json')
    args=parser.parse_args()
    atlas_path=WEB/'chemical-space.json'
    atlas=json.loads(atlas_path.read_text('utf8'))
    blocks={}
    with (args.inventory/'merged-building-blocks.csv').open(encoding='utf-8-sig',newline='') as stream:
        for rec in csv.DictReader(stream):
            s=rec['canonical isomeric SMILES'];i=int(rec['audit ID'])
            assert atlas['blocks'][str(i)]==rec['label']
            blocks[i]=[s,hashlib.sha256(s.encode()).hexdigest()[:24]]
    assert len(set(v[1] for v in blocks.values()))==len(atlas['blocks'])==len(blocks)
    stride=atlas['stride'];families={};jobs=[]
    for key,info in atlas['families'].items():
        if not info['eligible']:continue
        rows=[pair_row(code,stride,blocks,RULES[key]) for code in info['eligible']]
        families[key]=rows
        reported={tuple(p['ids']) for p in atlas['points'] if p['f']==key}
        pools=defaultdict(list)
        for row in rows:
            code,forward,t,mass=row[:4];a,b=divmod(code,stride)
            if (a,b) in reported or t not in SUPPORTED_TEMPLATES or mass>1350:continue
            left,right=(a,b) if forward else (b,a)
            # Bound the initial offline run; round-robin over nets and repeat-mass bands.
            if '@' in blocks[left][0]+blocks[right][0]:continue
            pools[t,int(mass//350)].append((key,code,t,blocks[left][0],blocks[right][0]))
        for pool in pools.values():
            pool.sort(key=lambda j:hashlib.sha256(repr(j).encode()).digest())
        chosen=0
        while chosen<args.attempts_per_family and any(pools.values()):
            for group in sorted(pools):
                if pools[group] and chosen<args.attempts_per_family:
                    jobs.append(pools[group].pop());chosen+=1
        print(key,len(rows),'candidates;',chosen,'periodic attempts',flush=True)
    models={};failed=defaultdict(dict)
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        for k,code,result in pool.map(periodic,jobs,chunksize=1):
            if 'reason' in result:failed[k][code]=result['reason']
            else:models.setdefault(k,{})[code]=result
    print('Periodic sample',sum(map(len,models.values())),'accepted of',len(jobs),flush=True)
    model_dir=ROOT/'data/crystal-models'
    catalog=json.loads((model_dir/'manifest.json').read_text('utf8'))
    references=[]
    for group in ['curated_models','reported_models','templates']:
        for r in catalog[group]:
            path=model_dir/r['path'];assert hashlib.sha256(path.read_bytes()).hexdigest()==r['sha256']
            references.append(dict(**geometry(import_cif(path.read_text('utf8'))),id=r['id'],name=r['label'],
                doi=r.get('citation_doi'),source='template' if group=='templates' else 'literature',sha256=r['sha256']))
    result=dict(schema=1,atlas_sha256=hashlib.sha256(atlas_path.read_bytes()).hexdigest(),
        columns=COLUMNS,blocks=blocks,families=families,rules=RULES,models=models,failed=failed,references=references,
        rdkit=rdBase.rdkitVersion,grid_spacing_angstrom=.9,probe_radius_angstrom=1.82,surface_samples_per_atom=16)
    args.output.write_text(json.dumps(result,ensure_ascii=False,separators=(',',':'),allow_nan=False),'utf8')
    print('Saved',args.output,flush=True)


if __name__=='__main__':
    os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
    main()
