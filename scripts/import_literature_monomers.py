"""Admit reviewed finite precursor graphs; keep private publications outside Git.

The portable selection is reproducible catalog data. Admission also requires
the exact source files and two-agent graph-review receipts in an evidence root.
This is an agent-reviewed lane, never a promotion to human ground truth.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import json
from pathlib import Path
import re

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

from import_manual_monomers import digest, finite_molecule, read, refresh_manifest, write
from cof_vision.cofdraw_chemistry import inspect_monomer


def validate_row(row):
    if not re.fullmatch(r'[a-z][a-z0-9_]{1,79}', row.get('id', '')):
        raise ValueError('Invalid monomer ID')
    if row.get('ground_truth') is not False or row.get('reviewer_type') != 'agent':
        raise ValueError('Literature admission must retain provisional agent review')
    mol=finite_molecule(row['smiles'])
    canonical=Chem.MolToSmiles(mol,isomericSmiles=True)
    if canonical!=row['canonical_smiles'] or row['smiles']!=canonical:
        raise ValueError('Canonical finite graph mismatch')
    if rdMolDescriptors.CalcMolFormula(mol)!=row['formula']:
        raise ValueError('Formula disagrees with the finite graph')
    if rdMolDescriptors.CalcMolFormula(mol,separateIsotopes=True)!=row['isotopic_formula']:
        raise ValueError('Isotopic formula disagrees with graph')
    inspection=inspect_monomer(canonical)
    # This importer verifies the signed 37-monomer selection against the
    # reaction vocabulary used at its review date. Live capabilities evolve.
    links={'imine','beta_ketoenamine','keto_hydrazide','acylhydrazone','azine','secondary_amine','nitrone'}
    kinds={'aryl_aldehyde','primary_arylamine','tp_trialdehyde','acylhydrazide','hydrazine','aryl_hydroxylamine','boronic_acid','nitrile'}
    roles=[r for r in inspection['compatible_roles'] if r['linkage'] in links]
    groups=[g for g in inspection['reactive_groups'] if g['kind'] in kinds]
    if row['compatible_roles']!=roles or row['reactive_groups']!=groups:
        raise ValueError('Runtime functional groups/roles differ from reviewed selection')
    supported=any(role['reactive_sites']>=2 for role in roles)
    if row['automatic_prepare_supported'] is not supported:
        raise ValueError('Automatic construction scope disagrees with runtime roles')
    if row['catalog_use']!=('scoped_drawing_candidate' if supported else 'reference_only'):
        raise ValueError('Drawing/reference scope mismatch')
    if not row.get('qualification') or not row.get('source_evidence'):
        raise ValueError('Source evidence and qualification are required')
    for evidence in row['source_evidence']:
        if evidence.get('ground_truth') is not False or evidence.get('review_status')!='two_agent_source_graph_review':
            raise ValueError('Each source occurrence requires two-agent graph review')
        if evidence.get('source_doi') not in row['source_dois']:
            raise ValueError('Source DOI is not attached to the monomer')
        asset_key='source_asset_sha256' if evidence.get('source_document_kind')=='image' else 'source_pdf_sha256'
        if asset_key=='source_asset_sha256' and not str(evidence.get('source_asset_url','')).startswith('https://'):
            raise ValueError('Image-only evidence requires its explicit public source URL')
        for key in (asset_key,'candidate_selection_sha256','independent_review_receipt_sha256'):
            if not re.fullmatch(r'[0-9a-f]{64}',evidence.get(key,'')):
                raise ValueError('Missing exact source/selection/review hash')
    serialized=json.dumps(row)
    if re.search(r'(?:^|["\s])[A-Za-z]:[\\/]|/Users/|/home/|source_page_text|source_evidence_private_path',serialized):
        raise ValueError('Portable row contains a private path or source text field')
    return canonical


def verify_evidence(rows,evidence_root):
    required=set()
    for row in rows:
        for evidence in row['source_evidence']:
            for key in ('source_pdf_sha256','source_asset_sha256','source_page_image_sha256','si_pdf_sha256',
                        'candidate_selection_sha256','independent_review_receipt_sha256'):
                if evidence.get(key):required.add(evidence[key])
            required.update(page['sha256'] for page in evidence.get('reviewed_page_images',[]))
    found={}
    for path in Path(evidence_root).rglob('*'):
        if path.is_file() and path.suffix.lower() in ('.json','.pdf','.png'):
            sha=digest(path)
            if sha in required:found[sha]=path
    if required-set(found):raise ValueError(f'{len(required-set(found))} pinned evidence files are absent or changed')
    for row in rows:
        for evidence in row['source_evidence']:
            selection=read(found[evidence['candidate_selection_sha256']])
            candidates=selection if isinstance(selection,list) else selection['candidates']
            matches=[r for r in candidates if r['candidate_id']==row['id']]
            if len(matches)!=1 or (matches[0].get('canonical_smiles') or matches[0].get('canonical_isomeric_smiles'))!=row['canonical_smiles']:
                raise ValueError('Pinned candidate does not match admitted ID and graph')
            receipt=read(found[evidence['independent_review_receipt_sha256']])
            if receipt.get('ground_truth') is not False or receipt.get('selection_sha256')!=evidence['candidate_selection_sha256']:
                raise ValueError('Independent review is not bound to this candidate selection')
            matches=[r for r in receipt['reviews'] if r['candidate_id']==row['id']]
            if len(matches)!=1 or matches[0].get('canonical_smiles')!=row['canonical_smiles'] or matches[0].get('result')!='finite_graph_agrees_with_viewed_source':
                raise ValueError('Independent review does not admit this exact graph')
            reviewed_hashes={r['sha256'] for r in matches[0]['viewed_pages']}
            if not reviewed_hashes or not reviewed_hashes<={p['sha256'] for p in evidence.get('reviewed_page_images',[])}:
                raise ValueError('Portable evidence omits independently reviewed source pages')
    return len(required)


def import_selection(selection_path,evidence_root,destination,linkages_path=None):
    selection_path=Path(selection_path);destination=Path(destination)
    selection=read(selection_path);rows=selection['monomers']
    if selection.get('ground_truth') is not False:raise ValueError('Selection must be provisional')
    ids=set();graphs=set()
    for row in rows:
        graph=validate_row(row)
        if row['id'] in ids or graph in graphs:raise ValueError('Selection has a duplicate ID or graph')
        ids.add(row['id']);graphs.add(graph)
    evidence_count=verify_evidence(rows,evidence_root)
    path=destination/'catalog/current.json';audit_path=destination/'catalog/current.audit.json'
    original=path.read_bytes();original_audit=audit_path.read_bytes();previous_sha=digest(path)
    catalog=read(path);existing={r['id']:r for r in catalog['monomer_presets']}
    linkage_additions=[]
    if linkages_path:
        linkage_package=read(linkages_path)
        if linkage_package.get('ground_truth') is not False:raise ValueError('Linkage references must remain provisional')
        old_linkages={r['id']:r for r in catalog['linkage_catalog']}
        for row in linkage_package['linkage_catalog']:
            if row.get('ground_truth') is not False or row.get('automatic_prepare_supported') is not False or row.get('candidate_transform') is not False:
                raise ValueError('Reference admission cannot add an executable linkage')
            if not row.get('qualification') or not row.get('primary_sources'):raise ValueError('Linkage source and scope required')
            if row['id'] in old_linkages:
                if row!=old_linkages[row['id']]:raise ValueError('Existing linkage changed')
            else:linkage_additions.append(row)
    existing_graphs={Chem.MolToSmiles(finite_molecule(r['smiles']),isomericSmiles=True):r['id'] for r in existing.values()}
    additions=[]
    for row in rows:
        if row['id'] in existing:
            if existing[row['id']]!=row:raise ValueError('Existing reviewed ID changed; explicit re-review required')
            continue
        if row['canonical_smiles'] in existing_graphs:raise ValueError('Existing graph has a different ID; merge its source evidence explicitly')
        additions.append(row);existing_graphs[row['canonical_smiles']]=row['id']
    if not additions and not linkage_additions:return {'added':0,'monomer_count':len(existing),'unchanged':True,'verified_evidence_files':evidence_count}
    for suffix,content in (('.json',original),('.audit.json',original_audit)):
        saved=destination/'catalog/history'/(previous_sha+suffix)
        if saved.exists() and saved.read_bytes()!=content:raise ValueError('Conflicting historical snapshot')
        if not saved.exists():saved.write_bytes(content)
    catalog['monomer_presets'].extend(additions)
    catalog['linkage_catalog'].extend(linkage_additions)
    catalog['catalog_version']='cofdraw-literature-reviewed-monomers-8'
    catalog['literature_monomer_update']={'selection_sha256':digest(selection_path),'previous_catalog_sha256':previous_sha,
         'added_ids':[r['id'] for r in additions],'linkage_reference_ids':[r['id'] for r in linkage_additions],
         'linkage_reference_selection_sha256':digest(linkages_path) if linkages_path else None,
         'reviewer_type':'agent','ground_truth':False,'source_drawings_unchanged':True}
    coverage=destination/'catalog/yaghi-literature-coverage.json'
    if coverage.is_file():
        catalog['literature_coverage']={'path':'catalog/yaghi-literature-coverage.json','sha256':digest(coverage),'ground_truth':False}
    write(path,catalog)
    report={'schema_version':'cofdraw-literature-monomer-admission-1','created_at_utc':datetime.now(timezone.utc).isoformat(),
      'previous_catalog_sha256':previous_sha,'catalog_sha256':digest(path),'selection_sha256':digest(selection_path),
      'original_monomer_count':len(existing),'monomer_count':len(catalog['monomer_presets']),
      'added_ids':[r['id'] for r in additions],'automatic_drawing_candidate_ids':[r['id'] for r in additions if r['automatic_prepare_supported']],
      'reference_only_ids':[r['id'] for r in additions if not r['automatic_prepare_supported']],
      'verified_evidence_files':evidence_count,'registry_sha256':digest(destination/'registry.json'),
      'linkage_reference_ids':[r['id'] for r in linkage_additions],
      'linkage_reference_selection_sha256':digest(linkages_path) if linkages_path else None,
      'reviewer_type':'agent','ground_truth':False,'source_pixels_distributed':False,'research_annotations_modified':False,
      'scope':'Selected finite graphs checked against source figures and separately reviewed by a second agent. Source chemistry, periodic identity and human ground truth are not inferred.'}
    admission=destination/'catalog/literature-monomer-admission.audit.json';write(admission,report)
    audit=json.loads(original_audit);audit.update(catalog_sha256=digest(path),monomer_count=report['monomer_count'],literature_monomer_admission_sha256=digest(admission))
    write(audit_path,audit);refresh_manifest(destination)
    return report


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--selection',type=Path,required=True)
    parser.add_argument('--evidence-root',type=Path,required=True)
    parser.add_argument('--destination',type=Path,required=True)
    parser.add_argument('--linkages',type=Path)
    args=parser.parse_args()
    print(json.dumps(import_selection(args.selection,args.evidence_root,args.destination,args.linkages),indent=2))
