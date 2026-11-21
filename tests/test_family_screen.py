import csv
import hashlib
import io
import json
import zipfile
import runpy
from pathlib import Path
import pytest

from rdkit import Chem
from rdkit.Chem import Descriptors

from cof_vision.cofdraw_family_screen import WEB, block_svg, export_family, screen


def test_family_rows_and_periodic_models_match_their_source():
    s=screen()
    path=WEB/'chemical-space.json'
    assert hashlib.sha256(path.read_bytes()).hexdigest()==s['atlas_sha256']
    atlas=json.loads(path.read_text('utf8'))
    for family,rows in s['families'].items():
        assert [r[0] for r in rows]==atlas['families'][family]['eligible']
        positive={p['ids'][0]*atlas['stride']+p['ids'][1] for p in atlas['points'] if p['f']==family and len(p['ids'])==2}
        for key,model in s['models'].get(family,{}).items():
            assert int(key) not in positive and key not in s['failed'].get(family,{})
            assert hashlib.sha256(model['cif'].encode()).hexdigest()==model['sha256']
            assert model['pld']<=model['lcd']+.0001 and 0<=model['void']<=100


def test_imine_repeat_mass_accounts_for_condensation():
    s=screen()
    by_smiles={v[0]:int(k) for k,v in s['blocks'].items()}
    inputs=['O=Cc1cc(C=O)cc(C=O)c1','Nc1ccc(N)cc1']
    a,b=sorted(by_smiles[x] for x in inputs)
    stride=max(by_smiles.values())+1
    row=next(r for r in s['families']['imine'] if r[0]==a*stride+b)
    mw=lambda sm:Descriptors.MolWt(Chem.MolFromSmiles(sm))
    assert abs(row[3]-(2*mw(inputs[0])+3*mw(inputs[1])-6*mw('O')))<.01
    assert row[5]==0 and '<svg' in block_svg(a)


def test_export_masks_pending_targets_and_contains_exact_cifs(tmp_path):
    path=tmp_path/'test.zip'
    export_family('azine',path)
    with zipfile.ZipFile(path) as archive:
        rows=list(csv.DictReader(io.StringIO(archive.read('candidates.csv').decode())))
        assert len(rows)==len(screen()['families']['azine'])
        assert all(r['band_gap_eV']=='0' and r['band_gap_eV_status']=='pending' for r in rows)
        assert len({r['candidate_id'] for r in rows})==len(rows)
        for r in rows:
            if r['model_status']=='computed_AA':
                cif=archive.read('models/'+r['model_cif_sha256']+'.cif')
                assert hashlib.sha256(cif).hexdigest()==r['model_cif_sha256'] and float(r['density'])>0
            else:
                assert r['density']=='0' and not r['model_cif_sha256']


def test_campaign_confirmations_require_the_same_cell_and_complete_triplicates(tmp_path):
    merge = runpy.run_path(str(Path(__file__).resolve().parents[1]/'scripts/build_campaign_snapshot.py'))['merge_confirmations']
    table = tmp_path/'confirmed.csv'
    head = 'id,paired_seeds,delivery_mmol_g,delivery_vv,delivery_sd_mmol_g,mixture_seeds,selectivity,sha256\n'
    checksum = 'a'*64
    row = {'id':'imine-10','sha256':checksum,'values':{}}
    table.write_text(head+f'imine-10,3,0,0,0,2,12,{checksum}\n',encoding='utf8')
    merge([row],table,'confirmation')
    assert row['values']['confirmation_delivery_mmol_g']==0
    assert 'confirmation_selectivity' not in row['values']
    assert row['confirmation_seeds']=={'delivery':3,'mixture':2}
    table.write_text(head+f'imine-10,3,1,1,1,3,4,{"b"*64}\n',encoding='utf8')
    with pytest.raises(ValueError,match='CIF mismatch'):
        merge([row],table,'deep')
