"""Refresh the drawing library from a sealed COFvision catalog and its finite-feed census."""
import argparse
from collections import Counter
from pathlib import Path

from rdkit import Chem
from import_manual_monomers import (REPO, checked_file, digest, finite_molecule, inspect_monomer,
                                   make_row, read, refresh_manifest, write)
import hashlib
import json


def refresh(curation, destination):
    pointer = read(curation / 'current.json')
    catalog_ref = pointer['catalog']
    annotations = read(checked_file(catalog_ref['path'], catalog_ref['sha256']))
    report_ref = pointer['report']
    inventory = checked_file(report_ref['path'], report_ref['sha256']).parent
    census_file = inventory / 'source-linked-monomer-census.json'
    census = read(census_file)
    roles = set(census['counted_source_roles'])
    directory = destination / 'catalog'
    current = directory / 'current.json'
    catalog = read(current)
    previous = digest(current)
    if catalog.get('cofvision_census', {}).get('census_sha256') == digest(census_file):
        return catalog['cofvision_census']
    graphs = {Chem.MolToSmiles(finite_molecule(r['smiles'])): r for r in catalog['monomer_presets']}
    by_id = {}
    for ref in annotations['annotation_files']:
        path = checked_file(ref['path'], ref['sha256'])
        for line_number, line in enumerate(path.read_text('utf-8-sig').splitlines(), 1):
            if line.strip():
                row = json.loads(line)
                if row.get('task') == 'molecule' and row.get('review_status') == 'agent_reviewed':
                    by_id[row['annotation_id']] = (row, ref['sha256'], line_number)
    added, excluded = [], Counter()
    for entry in census['records']:
        graph = entry['canonical_smiles']
        if graph in graphs:
            continue
        for identifier in entry['reviewed_ids']:
            if identifier not in by_id:
                continue
            row, file_hash, line = by_id[identifier]
            molecules = row.get('molecules', [])
            if len(molecules) != 1:
                continue
            m = molecules[0]
            if not roles.intersection((m.get('role'), row.get('molecule_role'))):
                continue
            name = m.get('full_name') or m.get('name') or m.get('label') or m.get('abbreviation')
            source = row.get('source') or m.get('source', {})
            if not name or not source.get('crop_sha256') or row.get('doi') not in entry['dois']:
                continue
            try:
                mol = finite_molecule(m['smiles'])
                if Chem.MolToSmiles(mol) != graph:
                    continue
                info = inspect_monomer(graph)
                supported = any(r['reactive_sites'] >= 2 for r in info['compatible_roles'])
                selection = dict(canonical_smiles=graph, source_crop_sha256=source['crop_sha256'],
                    id='cv_'+hashlib.sha256(graph.encode()).hexdigest()[:16], abbreviation=m.get('abbreviation') or m.get('label') or name,
                    full_name=name, catalog_use='scoped_drawing_candidate' if supported else 'reference_only',
                    qualification='Finite framework feed in the cited source. Automatic drawing depends on the selected partner and supported linkage.',
                    visual_review='COFvision agent-reviewed finite molecule crop, matched to its recorded isomeric SMILES.')
                new = make_row(selection, row, annotation_sha256=file_hash, annotation_line=line,
                               curation_sha256=catalog_ref['sha256'])
            except (ValueError, KeyError, OSError) as error:
                excluded[type(error).__name__ + ': ' + str(error).split(':')[0][:72]] += 1
                continue
            catalog['monomer_presets'].append(new)
            graphs[graph] = new
            added.append(new['id'])
            break
    if not added:
        return {'added': 0}
    info = dict(inventory=inventory.name, census_sha256=digest(census_file), upstream_catalog_sha256=catalog_ref['sha256'],
                reviewed_graphs=len(census['records']), added_ids=added, ground_truth=False)
    catalog.update(catalog_version='cofdraw-'+inventory.name, cofvision_census=info)
    body = json.dumps(catalog, ensure_ascii=False)
    if any(v in body for v in (':\\\\', '/Users/', '/home/')):
        raise ValueError('A nonportable path was found')
    audit_path = directory / 'current.audit.json'
    for path, suffix in ((current,'.json'),(audit_path,'.audit.json')):
        historical = directory / 'history' / (previous+suffix)
        if historical.exists() and historical.read_bytes() != path.read_bytes():
            raise ValueError('Historical catalog conflict')
        historical.write_bytes(path.read_bytes())
    write(current, catalog)
    audit = read(audit_path)
    audit.update(catalog_sha256=digest(current), previous_catalog_sha256=previous,
                 monomer_count=len(graphs), cofvision_census_sha256=info['census_sha256'])
    write(audit_path, audit)
    refresh_manifest(destination)
    return {'added':len(added),'total':len(catalog['monomer_presets']),'skipped_occurrence_checks':dict(excluded)}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--curation',type=Path,required=True)
    parser.add_argument('--destination',type=Path,default=REPO/'data/framework-drawing')
    args = parser.parse_args()
    print(json.dumps(refresh(args.curation,args.destination), indent=2))
