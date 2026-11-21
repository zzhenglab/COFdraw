"""Pack completed COFdraw campaign receipts for the local chemical-space explorer."""
import argparse
import base64
import csv
from datetime import datetime, timezone
import gzip
import hashlib
import json
import math
import os
from pathlib import Path
import statistics
import tempfile

WEB = Path(__file__).resolve().parents[1] / 'src/cof_vision/web'
GEOMETRY = {
    'density_g_cm3': 'Density (g/cm³)', 'pld_angstrom': 'Limiting diameter (Å)',
    'lcd_angstrom': 'Largest cavity (Å)', 'void_fraction': 'Void fraction (0–1)',
    'transport_accessible_probe_center_volume_cm3_g': 'Connected pore volume (cm³/g)',
    'transport_accessible_surface_m2_g': 'Connected surface (m²/g)',
}
GAS = {
    'delivery_mmol_g': 'CH₄ delivery (mmol/g)', 'delivery_cm3_STP_cm3': 'CH₄ delivery (cm³ STP/cm³)',
    'dilute_CH4_H2': 'CH₄/H₂ Henry ratio',
}
CONFIRM = {
    'delivery_mmol_g': 'Long-run CH₄ delivery (mmol/g)', 'delivery_vv': 'Long-run CH₄ delivery (cm³ STP/cm³)',
    'delivery_sd_mmol_g': 'Long-run delivery seed SD (mmol/g)', 'selectivity': 'Mixture CH₄/H₂ selectivity',
}
SEEDS = [7127, 17807]


def read(path):
    return json.loads(path.read_bytes())


def digest(data):
    return hashlib.sha256(data).hexdigest()


def number(value):
    if value in (None, ''):
        return None
    result = float(value)
    if not math.isfinite(result):
        raise ValueError('Non-finite screening value')
    return result


def build_record(path, receipt, recipe, cutoff):
    g = read(path/'geometry.json')
    row = dict(id=recipe['id'], family=recipe['family'], pair=recipe['pair'], status=receipt['status'],
               template=recipe['template'], indexed=bool(recipe['atlas_same_family_refs']), values={})
    model_sha = g.get('files', {}).get('structure.cif')
    if not model_sha:
        if receipt['status'] != 'build_rejected':
            raise ValueError(f'Missing CIF identity for {row["id"]}')
        row['reason'] = g.get('reason', 'Cell fit failed')
        return row
    cif = (path/'structure.cif').read_bytes()
    if digest(cif) != model_sha:
        raise ValueError(f'Changed CIF for {row["id"]}')
    geom = g['geometry']
    row.update(sha256=model_sha, cif_gzip=base64.b64encode(gzip.compress(cif, compresslevel=5, mtime=0)).decode(),
               channel=geom['connectivity']['status'], pld_bounds=geom['pld_bounds_angstrom'],
               lcd_bounds=geom['lcd_bounds_angstrom'])
    row['values'].update({key: number(geom[key]) for key in GEOMETRY})
    if row['status'] == 'screened':
        if g.get('status') != 'gas_eligible':
            raise ValueError(f'Gas results failed geometry qualification for {row["id"]}')
        gas = read(path/'gas.json')
        if gas['seeds'] != SEEDS:
            raise ValueError(f'Unexpected seeds for {row["id"]}')
        replicates = {key: [] for key in GAS}
        for seed in SEEDS:
            values = []
            for mode, species, pressure in [('gcmc','CH4',65), ('gcmc','CH4',5.8), ('widom','CH4',1), ('widom','H2',1)]:
                job = read(path/f'{mode}-{species}-{pressure}-{seed}.json')
                if (job['cif_sha256'], job['seed'], job['bar'], job['gas'], job['mode']) != (model_sha, seed, pressure, species, mode):
                    raise ValueError(f'Native receipt mismatch for {row["id"]}')
                values.append(number(job['value']))
            delivery = values[0] - values[1]
            replicates['delivery_mmol_g'].append(delivery)
            replicates['delivery_cm3_STP_cm3'].append(delivery*geom['density_g_cm3']*22.41397)
            if values[3] > 0:
                replicates['dilute_CH4_H2'].append(values[2]/values[3])
        for key in GAS:
            complete = len(replicates[key]) == len(SEEDS)
            for field, fn, suffix in [('mean',statistics.mean,''),('seed_sd',statistics.stdev,'_sd')]:
                value = number(gas[key][field])
                if (complete and (value is None or not math.isclose(value, fn(replicates[key]), abs_tol=1e-10))) or (not complete and value is not None):
                    raise ValueError(f'Summary and native seeds disagree for {row["id"]} {key}')
                row['values'][key+suffix] = value
    scan_done = path/'sensitivity-done.json'
    if scan_done.exists() and read(scan_done)['time'] <= cutoff:
        scan = read(path/'sensitivity.json')
        if scan['cif_sha256'] != model_sha:
            raise ValueError(f'Spacing test CIF mismatch for {row["id"]}')
        volumes = [s['transport_accessible_probe_center_volume_cm3_g'] for s in scan['cases'] if s['graph_accepted']]
        row['valid_spacings'] = len(volumes)
        if len(volumes) == 3:
            row['values']['spacing_span'] = max(volumes) - min(volumes)
    return row


def merge_confirmations(rows, csv_file, stage):
    by_id = {r['id']: r for r in rows}
    seen = set()
    with csv_file.open(encoding='utf-8-sig', newline='') as stream:
        for item in csv.DictReader(stream):
            key = item['id']
            if key in seen:
                raise ValueError(f'Duplicate confirmation for {key}')
            seen.add(key)
            row = by_id.get(key)
            if row is None:
                continue
            if not row.get('sha256') or item['sha256'] != row['sha256']:
                raise ValueError(f'Confirmation CIF mismatch for {key}')
            paired, mixture = int(item['paired_seeds']), int(item['mixture_seeds'])
            if not (0 <= paired <= 3 and 0 <= mixture <= 3):
                raise ValueError(f'Invalid {stage} seed count for {key}')
            row[stage+'_seeds'] = {'delivery': paired, 'mixture': mixture}
            for field in CONFIRM:
                if (mixture if field == 'selectivity' else paired) >= 3:
                    row['values'][stage+'_'+field] = number(item[field])


def breadth_inputs(directories, campaign):
    output = []
    for directory in directories:
        if read(directory/'campaign.json') != campaign:
            raise ValueError('Breadth run uses a different recipe inventory or screen protocol')
        manifest = read(directory/'input-receipt.json')
        if manifest['protocol'] != campaign['protocol']:
            raise ValueError('Breadth run protocol mismatch')
        with (directory/'catalog.csv').open(encoding='utf-8-sig', newline='') as stream:
            catalog = {}
            for entry in csv.DictReader(stream):
                if entry['id'] in catalog:
                    raise ValueError(f'Duplicate breadth catalog entry {entry["id"]}')
                catalog[entry['id']] = entry
        output.append((directory, manifest, catalog))
    return output


def breadth_record(directory, recipe, entry, protocol):
    encoded = json.dumps(recipe, sort_keys=True, ensure_ascii=False, separators=(',', ':')).encode()
    if (entry['family'], int(entry['pair']), entry['recipe_sha256']) != (recipe['family'], recipe['pair'], digest(encoded)):
        raise ValueError(f'Frozen breadth recipe mismatch for {recipe["id"]}')
    if entry['status'] == 'operational_error':
        return dict(id=recipe['id'], family=recipe['family'], pair=recipe['pair'], template=recipe['template'],
                    indexed=bool(recipe['atlas_same_family_refs']), status='operational_error', values={})
    if entry['protocol'] != protocol:
        raise ValueError(f'Breadth job protocol mismatch for {recipe["id"]}')
    path = directory/'records'/recipe['family']/recipe['id']
    result = build_record(path, entry, recipe, float('inf'))
    if (entry.get('sha256') or None) != result.get('sha256'):
        raise ValueError(f'Breadth catalog CIF mismatch for {recipe["id"]}')
    return result


def build(args):
    if args.as_of and (args.extra_results or args.confirmation_csv or args.deep_csv):
        raise ValueError('Use the live snapshot for extra campaigns; their catalog rows have no cutoff timestamps')
    campaign = read(args.results/'campaign.json')
    snapshot = campaign['snapshot']
    screen = read(WEB/'chemical-screen.json')
    if screen['atlas_sha256'] != snapshot['atlas_sha256'] or digest((WEB/'chemical-screen.json').read_bytes()) != snapshot['screen_sha256']:
        raise ValueError('Campaign and installed explorer use different building-block catalogs')
    if digest(args.recipes.read_bytes()) != snapshot['recipes_sha256']:
        raise ValueError('Frozen recipe inventory digest mismatch')
    at = datetime.fromisoformat(args.as_of.replace('Z', '+00:00')) if args.as_of else datetime.now(timezone.utc)
    if at.tzinfo is None:
        raise ValueError('Specify a timezone for the cutoff')
    cutoff = at.timestamp()
    rows = []
    breadth = breadth_inputs(args.extra_results, campaign)
    with gzip.open(args.recipes, 'rt', encoding='utf8') as stream:
        for recipe in map(json.loads, stream):
            if recipe['flags']:
                continue
            path = args.results/'records'/recipe['family']/recipe['id']
            marker = path/'screen-done.json'
            record = None
            if marker.is_file():
                receipt = read(marker)
                if receipt['time'] <= cutoff:
                    record = build_record(path, receipt, recipe, cutoff)
                    record['cohorts'] = ['primary']
            for i, (directory, _, catalog) in enumerate(breadth, start=1):
                entry = catalog.pop(recipe['id'], None)
                if entry is None:
                    continue
                extra = breadth_record(directory, recipe, entry, campaign['protocol'])
                if record is None:
                    record = extra
                    record['cohorts'] = []
                elif any(record.get(k) != extra.get(k) for k in ('family','pair','sha256','status','values')):
                    raise ValueError(f'Conflicting outcomes for overlapping recipe {recipe["id"]}')
                record['cohorts'].append('breadth_'+str(i))
            if record:
                rows.append(record)
    if any(catalog for _,_,catalog in breadth):
        raise ValueError('A breadth catalog contains a held or unrecognized recipe')
    methods = {'screen': 'Unrelaxed AA cell, checked against intended bonds and close contacts. Geometry uses a 1.82 Å probe, 0.65 Å requested grid and 64 surface samples per atom. At 298 K the rigid zero-charge RASPA3 screen uses two independent seeds. CH₄ delivery is the 65 to 5.8 bar swing (500 initialization and 1,000 production cycles per pressure). Dilute CH₄/H₂ is the ratio of Widom Henry coefficients (3,000 cycles per gas). Volumetric delivery uses ideal-crystal density. SD is between seeds. Spacing sensitivity spans three graph-preserving heights at ±0.25 Å.'}
    metrics = {**GEOMETRY, **GAS, **{k+'_sd': label+' · seed SD' for k,label in GAS.items()},
               'spacing_span': 'Spacing sensitivity, volume span (cm³/g)'}
    for i,(_,manifest,_) in enumerate(breadth, start=1):
        methods['breadth_'+str(i)] = manifest
    for stage,table,method in [('confirmation',args.confirmation_csv,args.confirmation_method),
                               ('deep',args.deep_csv,args.deep_method)]:
        if bool(table) != bool(method):
            raise ValueError(f'Supply both the {stage} table and its method metadata')
        if table:
            merge_confirmations(rows, table, stage)
            methods[stage] = read(method)
            prefix = 'Extended run · ' if stage == 'deep' else ''
            metrics.update({stage+'_'+k: prefix+v for k,v in CONFIRM.items()})
    result = dict(schema_version=1, format='cofdraw-screen-results', snapshot_utc=at.astimezone(timezone.utc).isoformat(),
                  protocol=campaign['protocol'], atlas_sha256=snapshot['atlas_sha256'], recipes_sha256=snapshot['recipes_sha256'],
                  methods=methods, metrics=metrics, records=rows)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(dir=args.output.parent, delete=False) as raw:
        temp = Path(raw.name)
    try:
        with gzip.open(temp, 'wt', compresslevel=5, encoding='utf8', newline='\n') as archive:
            json.dump(result, archive, ensure_ascii=False, allow_nan=False, separators=(',', ':'))
        os.replace(temp, args.output)
    finally:
        temp.unlink(missing_ok=True)
    print(f'{len(rows)} candidate outcomes, {sum(bool(r.get("sha256")) for r in rows)} checked CIFs in {args.output}')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    for arg in ('results','recipes','output'):
        parser.add_argument('--'+arg, type=Path, required=True)
    parser.add_argument('--as-of', help='ISO 8601 cutoff with timezone for primary and spacing jobs')
    parser.add_argument('--confirmation-csv', type=Path)
    parser.add_argument('--confirmation-method', type=Path)
    parser.add_argument('--deep-csv', type=Path)
    parser.add_argument('--deep-method', type=Path)
    parser.add_argument('--extra-results', type=Path, action='append', default=[], help='Additional frozen breadth run with catalog and completed receipts')
    build(parser.parse_args())
