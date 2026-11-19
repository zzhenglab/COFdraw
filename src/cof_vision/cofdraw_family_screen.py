"""Offline family screening assets, molecular depictions and batch exports."""
import csv
import io
import json
import zipfile
from functools import lru_cache
from pathlib import Path

from rdkit import Chem
from rdkit.Chem.Draw import rdMolDraw2D

from .cofdraw_depiction import compute_2d_coords

WEB = Path(__file__).with_name('web')
GEOMETRY = ['density','void','lcd','pld','cell_volume','pore_volume','free_volume','connected_volume','area','local_area']
PENDING = ['band_gap_eV','co2_uptake_mmol_g_298K_1bar','ch4_uptake_mmol_g_298K_1bar']


@lru_cache(maxsize=1)
def screen():
    return json.loads((WEB/'chemical-screen.json').read_text('utf8'))


@lru_cache(maxsize=256)
def block_svg(block):
    smiles=screen()['blocks'][str(block)][0]
    mol=Chem.MolFromSmiles(smiles)
    compute_2d_coords(mol)
    draw=rdMolDraw2D.MolDraw2DSVG(320,170)
    draw.drawOptions().padding=.08
    rdMolDraw2D.PrepareAndDrawMolecule(draw,mol)
    draw.FinishDrawing()
    return draw.GetDrawingText()


def export_family(family, path):
    payload=screen()
    rows=payload['families'][family]
    atlas=json.loads((WEB/'chemical-space.json').read_text('utf8'))
    reported={tuple(p['ids']):p['refs'] for p in atlas['points'] if p['f']==family and len(p['ids'])==2}
    fields=['candidate_id','family','drawing_rule','template','block_a_smiles','block_b_smiles','indexed_report','dois',
        'repeat_mass_g_mol','nitrogen_wt_pct','oxygen_wt_pct','junctions_mmol_g','aromatic_fraction_of_carbon_pct',
        'model_status','model_cif_sha256','channel_status']+GEOMETRY+[col for k in PENDING for col in (k,k+'_status')]
    with zipfile.ZipFile(path,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=3) as archive:
        with archive.open('candidates.csv','w') as raw, io.TextIOWrapper(raw,encoding='utf8',newline='') as stream:
            writer=csv.writer(stream);writer.writerow(fields)
            for code,forward,template,mass,n,o,links,aromatic in rows:
                a,b=divmod(code,atlas['stride'])
                ids=[a,b] if forward else [b,a]
                (sa,ha),(sb,hb)=[payload['blocks'][str(i)] for i in ids]
                model=payload['models'].get(family,{}).get(str(code),{})
                failed=str(code) in payload['failed'].get(family,{})
                status='computed_AA' if model else 'fit_failed' if failed else 'chain' if template=='chain_2_2' else 'not_sampled'
                refs=reported.get((a,b),[])
                writer.writerow([family+'-'+ha+'-'+hb,family,payload['rules'][family],template,sa,sb,bool(refs),'|'.join(sorted({doi for _,doi in refs})),
                    mass,n,o,links,aromatic,status,model.get('sha256',''),model.get('channel',''),
                    *[model.get(k) if model.get(k) is not None else ('' if model else 0) for k in GEOMETRY],
                    *[v for _ in PENDING for v in (0,'pending')]])
        for code,m in payload['models'].get(family,{}).items():
            archive.writestr('models/'+m['sha256']+'.cif',m['cif'])
        archive.writestr('methods.json',json.dumps(dict(family=family, rows=len(rows), atlas_sha256=payload['atlas_sha256'],
            rdkit=payload['rdkit'],candidate_id='Family and the two SHA-256 SMILES prefixes ordered by reaction role.',
            composition='Dry, fully reacted ideal repeat. Smallest integer block ratio balances reactive sites. Chain templates have no assigned 2D pores.',
            geometry='Unrelaxed generated AA cell, intended periodic edges and steric contacts checked. RDKit vdW sphere grid; 1.82 angstrom probe, 0.9 angstrom requested spacing, 16 surface samples per atom. LCD and PLD are grid lower estimates.',
            units={'density':'g/cm3','void':'percent','lcd':'angstrom','pld':'angstrom','cell_volume':'angstrom3',
                   'pore_volume':'geometric void cm3/g','free_volume':'probe-center cm3/g','connected_volume':'transport-connected probe-center cm3/g','area':'transport-accessible m2/g','local_area':'local probe-exposed m2/g'},
            missing='Geometry zeros require model_status=computed_AA to represent results. Empty surface fields within a computed model mean unavailable. Check channel_status for resolution limits. For electronic and adsorption targets use only rows with a computed target status. Gas placeholders are pure-component, 298 K, 1 bar.',
            literature='https://doi.org/10.1021/acsmaterialslett.4c02594'),indent=2))
