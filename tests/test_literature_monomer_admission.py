"""Independent invariants for the additive 37-precursor literature admission.

Public tests validate chemistry, portable bindings and importer failure paths.
They do not replace the private source-figure reviews or claim human gold.
"""
from collections import Counter
from copy import deepcopy
import hashlib
import importlib.util
import json
import math
from pathlib import Path
import re
import sys
import xml.etree.ElementTree as ET

from fastapi.testclient import TestClient
import pytest
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

from cof_vision.cofdraw_chemistry import inspect_monomer, prepare_template
from cof_vision.cofdraw_core import export_chemical_formats, read_catalog, render_styled
from cof_vision.cofdraw_desktop import create_desktop_app
from cof_vision.cofdraw_knowledge import query_catalog
from cof_vision.cofdraw_synthesis import synthetic_accessibility
from cof_vision.framework_drawing import VERSION, build_patch, linkage_fragment

REPO = Path(__file__).resolve().parents[1]
DATA = REPO / 'data/framework-drawing'
SELECTION_PATH = DATA / 'catalog/literature-monomer-selection.json'
SELECTION = json.loads(SELECTION_PATH.read_text(encoding='utf-8'))
SELECTED = {row['id']: row for row in SELECTION['monomers']}
BASELINE_81_SHA256 = '05f190c93835e687797f6ce3011a2879a6884275aa2e28ef2c26ff0e745993f0'

sys.path.insert(0, str(REPO / 'scripts'))
spec = importlib.util.spec_from_file_location('literature_monomer_importer', REPO / 'scripts/import_literature_monomers.py')
importer = importlib.util.module_from_spec(spec)
spec.loader.exec_module(importer)

# Deliberately chosen drawing incidences, not assignments to the papers' COFs.
# Nonplanar central methane/alcohol and tetraaryl-silicon examples use the
# unassigned finite projection rather than claiming an experimental planar net.
DRAWING_PAIRS = [
    ('tapp', 'ta', 'tapp', 'imine', 'sql_4_2'),
    ('ttf_cho4', 'ttf_cho4', 'pa', 'imine', 'sql_4_2'),
    ('tpe_hydroxylamine4', 'ta', 'tpe_hydroxylamine4', 'nitrone', 'sql_4_2'),
    ('tabta', 'tp', 'tabta', 'beta_ketoenamine', 'hcb_3_3'),
    ('dhbp_mom', 'dhbp_mom', 'pa', 'imine', 'chain_2_2'),
    ('dapy25', 'tp', 'dapy25', 'beta_ketoenamine', 'hcb_3_2'),
    ('oda44', 'tfb', 'oda44', 'imine', 'hcb_3_2'),
    ('tapm3', 'tfb', 'tapm3', 'imine', 'alternating_hexagon'),
    ('tapmol3', 'tfb', 'tapmol3', 'imine', 'alternating_hexagon'),
    ('daba44', 'tfb', 'daba44', 'imine', 'hcb_3_2'),
    ('pcba26', 'pcba26', 'pa', 'imine', 'chain_2_2'),
    ('tfs_cho4', 'tfs_cho4', 'pa', 'imine', 'alternating_hexagon'),
    ('ttfb_mes', 'ttfb_mes', 'pa', 'imine', 'hcb_3_2'),
    ('hapb', 'ta', 'hapb', 'imine', 'hxl_6_2'),
    ('hfpb', 'hfpb', 'pa', 'imine', 'hxl_6_2'),
    ('tftda', 'tftda', 'pa', 'imine', 'chain_2_2'),
    ('terphenyl_cooh2_cho2', 'terphenyl_cooh2_cho2', 'pa', 'imine', 'chain_2_2'),
    ('taps_e', 'ta', 'taps_e', 'imine', 'sql_4_2'),
    ('bpda_n3', 'bpda_n3', 'pa', 'imine', 'chain_2_2'),
    ('eba_n3', 'eba_n3', 'pa', 'imine', 'chain_2_2'),
]


def read(path):
    return json.loads(Path(path).read_text(encoding='utf-8'))


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


@pytest.fixture(scope='module')
def catalog():
    return read_catalog(REPO / 'data')


@pytest.fixture(scope='module')
def rows(catalog):
    return {row['id']: row for row in catalog['monomer_presets']}


@pytest.fixture(scope='module')
def library():
    with TestClient(create_desktop_app(REPO / 'data')) as client:
        response = client.get('/api/drawing/library')
        assert response.status_code == 200
        return response.json()


@pytest.mark.parametrize('identifier', list(SELECTED))
def test_every_admitted_graph_formula_isotope_stereo_and_sa_score(identifier, rows):
    row = rows[identifier]
    assert row == SELECTED[identifier]
    mol = Chem.MolFromSmiles(row['smiles'])
    assert mol is not None and len(Chem.GetMolFrags(mol)) == 1
    assert mol.GetNumAtoms() > 0
    assert all(a.GetAtomicNum() and not a.GetNumRadicalElectrons() and not a.GetAtomMapNum() for a in mol.GetAtoms())
    assert Chem.MolToSmiles(mol, isomericSmiles=True) == row['canonical_smiles'] == row['smiles']
    assert rdMolDescriptors.CalcMolFormula(mol) == row['formula']
    assert rdMolDescriptors.CalcMolFormula(mol, separateIsotopes=True) == row['isotopic_formula']
    restored = Chem.MolFromMolBlock(Chem.MolToMolBlock(mol), removeHs=False)
    assert restored is not None
    assert Chem.MolToSmiles(restored, isomericSmiles=True) == row['canonical_smiles']
    score = synthetic_accessibility(row['smiles'])
    assert score['available'] and score['fragment_table_verified'], score
    assert math.isfinite(score['score']) and 1 <= score['score'] <= 10
    assert score['canonical_smiles'] == row['canonical_smiles']
    assert score['formal_charge'] == Chem.GetFormalCharge(mol)
    assert score['limitations'] and row['ground_truth'] is False and row['reviewer_type'] == 'agent'
    assert row['node_symmetry'] is None
    assert row.get('cofex_id') is None and row.get('topology') is None
    assert not re.search(r'(?:^|["\s])[A-Za-z]:[\\/]|/Users/|/home/|source_page_text|source_evidence_private_path', json.dumps(row))


def test_exact_additive_history_and_manifest_bindings(catalog, rows):
    audit_path = DATA / 'catalog/literature-monomer-admission.audit.json'
    audit = read(audit_path)
    baseline_path = DATA / 'catalog/history' / (BASELINE_81_SHA256 + '.json')
    assert digest(baseline_path) == BASELINE_81_SHA256
    baseline = read(baseline_path)
    assert audit['previous_catalog_sha256'] == BASELINE_81_SHA256
    assert len(baseline['monomer_presets']) == audit['original_monomer_count'] == 81
    assert len(SELECTED) == len(audit['added_ids']) == 37
    assert len(rows) == audit['monomer_count'] == 118
    assert len({row['canonical_smiles'] for row in SELECTED.values()}) == 37
    assert set(audit['added_ids']) == set(SELECTED)
    assert set(rows) == {r['id'] for r in baseline['monomer_presets']} | set(SELECTED)
    for old in baseline['monomer_presets']:
        assert rows[old['id']] == old
    for field in ('recipe_bindings', 'source_metadata', 'source_registry_sha256'):
        assert catalog[field] == baseline[field]
    assert audit['registry_sha256'] == digest(DATA / 'registry.json')
    assert audit['catalog_sha256'] == digest(DATA / 'catalog/current.json')
    assert audit['selection_sha256'] == digest(SELECTION_PATH)
    assert audit['ground_truth'] is False and audit['source_pixels_distributed'] is False
    assert audit['research_annotations_modified'] is False
    assert audit['verified_evidence_files'] > 0
    current_audit = read(DATA / 'catalog/current.audit.json')
    assert current_audit['literature_monomer_admission_sha256'] == digest(audit_path)
    manifest = read(DATA / 'current-bundle-manifest.json')
    assert manifest['catalog_sha256'] == digest(DATA / 'catalog/current.json')
    indexed = {record['path']: record for record in manifest['files']}
    for relative in ('catalog/current.json', 'catalog/current.audit.json', 'catalog/literature-monomer-selection.json',
                     'catalog/literature-monomer-admission.audit.json', f'catalog/history/{BASELINE_81_SHA256}.json'):
        assert indexed[relative]['sha256'] == digest(DATA / relative)
        assert indexed[relative]['bytes'] == (DATA / relative).stat().st_size


def test_actual_runtime_roles_and_support_split(library):
    actual = {r['id']: r for r in library['monomer_presets']}
    automatic, reference = set(), set()
    for key, selected in SELECTED.items():
        inspection = inspect_monomer(selected['smiles'])
        assert actual[key]['compatible_roles'] == selected['compatible_roles'] == inspection['compatible_roles']
        eligible = [r for r in actual[key]['compatible_roles'] if r['reactive_sites'] >= 2]
        assert bool(eligible) == actual[key]['automatic_prepare_supported']
        (automatic if eligible else reference).add(key)
        for role in eligible:
            assert library['linkages'][role['linkage']]['supported'] is True
            assert role['role'] in library['linkages'][role['linkage']]['roles']
        if not eligible:
            assert actual[key]['catalog_use'] == 'reference_only'
    assert len(automatic) == 20 and len(reference) == 17
    assert automatic == {case[0] for case in DRAWING_PAIRS}
    # Constructor availability is deliberately separate from a paper's route.
    for key in ('dhbp_mom', 'bpda_n3', 'eba_n3'):
        assert actual[key]['automatic_prepare_supported'] is True
        assert actual[key]['reported_source_reaction_supported'] is False
    assert library['linkages']['nitrone']['supported'] is True
    assert actual['tpe_hydroxylamine4']['compatible_roles'][0]['role'] == 'aryl_hydroxylamine'


def test_literature_query_returns_exact_coverage_and_live_new_monomer_roles(catalog):
    coverage_path = DATA / 'catalog/yaghi-literature-coverage.json'
    coverage = read(coverage_path)
    assert catalog['literature_coverage']['sha256'] == digest(coverage_path)
    assert coverage['counts']['selected_index_entries'] == len(coverage['papers']) == 100
    assert coverage['counts']['unique_selected_index_entries'] == 99
    assert coverage['counts']['duplicate_book_entries'] == 1
    assert coverage['counts']['full_main_and_si_audits_claimed'] == 0
    assert coverage['counts']['independently_checked_new_finite_monomers'] == 37
    actual = []
    for offset in (0, 50):
        answer = query_catalog(REPO / 'data', '', ['literature'], limit=50, offset=offset)
        assert answer['matched_count'] == 100 and not answer['errors']
        actual.extend(answer['records'])
    assert len({r['record_id'] for r in actual}) == 100
    expected = {paper['review_id']: paper for paper in coverage['papers']}
    for record in actual:
        assert record['data'] == expected[record['record_id']]
        assert record['provenance']['coverage_sha256'] == digest(coverage_path)
        assert record['support']['ground_truth'] is False
        assert record['support']['full_main_and_si_audit'] is False
        assert record['qualifications']
    for key in ('tpe_hydroxylamine4', 'tcpb_cyanomethyl', 'nbpda', 'bpda_n3'):
        answer = query_catalog(REPO / 'data', key, ['monomer'])
        assert not answer['errors']
        record = next(r for r in answer['records'] if r['record_id'] == key)
        assert record['data'] == SELECTED[key]
        assert record['support']['runtime_compatible_roles'] == inspect_monomer(SELECTED[key]['smiles'])['compatible_roles']
        assert record['support']['automatic_prepare_supported'] == SELECTED[key]['automatic_prepare_supported']
    nitrone = next(r for r in query_catalog(REPO / 'data', 'nitrone', ['linkage'])['records'] if r['record_id'] == 'nitrone')
    assert nitrone['support']['automatic_prepare_supported'] is True


def test_porphyrin_inner_n_h_survives_cache_refresh_and_real_chemical_exports(rows, tmp_path):
    recipe = prepare_template(rows['ta']['smiles'], rows['tapp']['smiles'], 'imine', 'TA', 'TAPP', 6, 'sql_4_2')
    recipe.update(schema_version=VERSION, recipe_id='porphyrin-hydrogen-regression', designed_variant=True)
    mol, audit = build_patch(recipe)
    # Repeated refreshes expose a cached implicit-H artifact that can otherwise
    # look correct immediately after restoring component aromatic flags.
    for _ in range(2):
        mol.UpdatePropertyCache(strict=True)
        inner = [a for a in mol.GetAtoms() if a.GetAtomicNum() == 7 and a.GetIsAromatic()]
        assert len(inner) == 16
        assert sum(a.GetTotalNumHs() for a in inner) == 8
        assert sum(a.GetNumExplicitHs() for a in inner) == 8
        assert not any(a.GetNumRadicalElectrons() for a in mol.GetAtoms())
    render_styled(recipe, tmp_path, 700, 'standard', prepared_patch=(mol, audit))
    formats = export_chemical_formats(tmp_path, recipe['recipe_id'], source_mol=mol)
    assert 'cdxml' in formats
    export_report = read(tmp_path / f'{recipe["recipe_id"]}.chemical-export.json')
    assert export_report['molecule_input'] == 'scoped_MOL_parser'
    assert 'aromatic N-H' in export_report['source_mol_fallback_reason']
    assert export_report['cdxml'] == 'validated_native_writer_explicit_XML_tables'
    xml = ET.parse(tmp_path / f'{recipe["recipe_id"]}.cdxml')
    nodes = list(xml.iter('n'))
    assert sum(n.get('NodeType') == 'ExternalConnectionPoint' for n in nodes) == 8
    node_index = {n.get('id'): i for i, n in enumerate(nodes)}
    # CDXML may use ordinary implicit valence rather than NumHydrogens for
    # pyrrolic N. Parse its bounded literal atom/bond table; treating a missing
    # NumHydrogens attribute as zero would incorrectly reject that valid form.
    xml_graph = Chem.RWMol()
    for node in nodes:
        boundary = node.get('NodeType') == 'ExternalConnectionPoint'
        atom = Chem.Atom(0 if boundary else int(node.get('Element', '6')))
        atom.SetFormalCharge(int(node.get('Charge', '0')))
        if boundary:
            atom.SetNoImplicit(True)
        if node.get('NumHydrogens') is not None:
            atom.SetNumExplicitHs(int(node.get('NumHydrogens')))
            atom.SetNoImplicit(True)
        xml_graph.AddAtom(atom)
    for bond in xml.iter('b'):
        a, b = node_index[bond.get('B')], node_index[bond.get('E')]
        order = {'1': Chem.BondType.SINGLE, '2': Chem.BondType.DOUBLE, '3': Chem.BondType.TRIPLE}[bond.get('Order', '1')]
        xml_graph.AddBond(a, b, order)
        if audit['atom_origin'][a]['instance'] != audit['atom_origin'][b]['instance']:
            assert bond.get('Order') == '2' and bond.get('Display') == 'Wavy'
    xml_mol = xml_graph.GetMol()
    xml_mol.UpdatePropertyCache(strict=True)
    assert sum(a.GetTotalNumHs() for a in xml_mol.GetAtoms() if a.GetAtomicNum() == 7) == 8
    assert Chem.GetFormalCharge(xml_mol) == 0
    for suffix in ('mol', 'sdf'):
        block = (tmp_path / f'{recipe["recipe_id"]}.{suffix}').read_text(encoding='utf-8').split('$$$$')[0]
        restored = Chem.MolFromMolBlock(block, removeHs=False)
        assert restored is not None
        assert rdMolDescriptors.CalcMolFormula(restored) == rdMolDescriptors.CalcMolFormula(mol)
        assert sum(a.GetTotalNumHs() for a in restored.GetAtoms() if a.GetAtomicNum() == 7) == 8
        assert not any(a.GetNumRadicalElectrons() for a in restored.GetAtoms())


@pytest.mark.parametrize('identifier,a,b,linkage,template', DRAWING_PAIRS, ids=[r[0] for r in DRAWING_PAIRS])
def test_each_supported_precursor_constructs_connected_scoped_patch(identifier, a, b, linkage, template, rows):
    recipe = prepare_template(rows[a]['smiles'], rows[b]['smiles'], linkage, a, b, 6, template)
    recipe.update(schema_version=VERSION, recipe_id='literature-admission-design', designed_variant=True)
    mol, audit = build_patch(recipe)
    mol.UpdatePropertyCache(strict=True)
    assert len(Chem.GetMolFrags(mol)) == 1
    assert audit['atoms'] == mol.GetNumAtoms() and audit['connections'] == len(recipe['connections'])
    assert min(recipe['compatibility']['reactive_site_counts']) >= 2
    assert not recipe.get('source') and not recipe.get('cofex_id')
    assert audit['source'] is None and audit['reference'] is None and audit['cofex_id'] is None
    assert audit['whole_cof_smiles'] is None and audit['ground_truth'] is False
    assert audit['synthetic_training_admitted'] is False
    assert 'neither is experimental identity' in recipe['design_identity_qualification']
    if template == 'alternating_hexagon':
        assert recipe['design_topology'] is None

    # Atom-origin mapping checks every surviving internal atom and bond, not
    # just aggregate formula. This catches pendant azide charges or source E
    # alkene stereo disappearing while a connected patch still looks valid.
    origins = {(r['instance'], r['fragment_atom_index']): i for i, r in enumerate(audit['atom_origin'])}
    expected_charge = 0
    for instance in recipe['instances']:
        fragment = Chem.MolFromSmiles(recipe['components'][instance['component']]['fragment_smiles'])
        for atom in fragment.GetAtoms():
            if not atom.GetAtomicNum():
                continue
            actual = mol.GetAtomWithIdx(origins[(instance['id'], atom.GetIdx())])
            expected_charge += atom.GetFormalCharge()
            assert (actual.GetAtomicNum(), actual.GetFormalCharge(), actual.GetIsotope(), actual.GetTotalNumHs()) == (
                atom.GetAtomicNum(), atom.GetFormalCharge(), atom.GetIsotope(), atom.GetTotalNumHs())
        for bond in fragment.GetBonds():
            if not bond.GetBeginAtom().GetAtomicNum() or not bond.GetEndAtom().GetAtomicNum():
                continue
            left, right = (origins[(instance['id'], i)] for i in (bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()))
            actual = mol.GetBondBetweenAtoms(left, right)
            assert actual is not None and actual.GetBondType() == bond.GetBondType()
            assert actual.GetIsAromatic() == bond.GetIsAromatic()
            assert actual.GetStereo() == bond.GetStereo()
    assert Chem.GetFormalCharge(mol) == expected_charge
    if identifier == 'taps_e':
        assert sum(b.GetStereo() == Chem.BondStereo.STEREOE for b in mol.GetBonds()) == 4
    if identifier in ('bpda_n3', 'eba_n3'):
        azides = mol.GetSubstructMatches(Chem.MolFromSmarts('[N-]=[N+]=N'))
        assert len(azides) == 6
        assert Counter(a.GetFormalCharge() for a in mol.GetAtoms())[1] == 6
        assert Counter(a.GetFormalCharge() for a in mol.GetAtoms())[-1] == 6


@pytest.mark.parametrize('identifier,expected_free_amines', [
    ('tcpb_cyanomethyl', 0), ('tfpn_dioxin', 0), ('tfpc_dioxin', 0),
    ('bpda_n3', 0), ('eba_n3', 0), ('dmbdu', 0),
    ('nbpda', 1), ('paba', 1), ('tabta', 3), ('daba44', 2), ('tpe_hydroxylamine4', 0),
])
def test_nitrile_azide_amide_boc_and_hydroxylamine_are_not_false_primary_amines(identifier, expected_free_amines, rows):
    inspected = inspect_monomer(rows[identifier]['smiles'])
    count = next((g['count'] for g in inspected['reactive_groups'] if g['kind'] == 'primary_arylamine'), 0)
    assert count == expected_free_amines
    if expected_free_amines == 0:
        with pytest.raises(ValueError, match='No supported reactive sites'):
            linkage_fragment(rows[identifier]['smiles'], 'imine', 'primary_arylamine')
    elif expected_free_amines == 1:
        with pytest.raises(ValueError, match='at least two compatible reactive sites'):
            prepare_template(rows['tfb']['smiles'], rows[identifier]['smiles'], 'imine', 'TFB', identifier, 6, 'alternating_hexagon')
    if identifier in ('tabta', 'daba44'):
        fragment = Chem.MolFromSmiles(linkage_fragment(rows[identifier]['smiles'], 'imine', 'primary_arylamine')['fragment_smiles'])
        assert len(fragment.GetSubstructMatches(Chem.MolFromSmarts('[C](=O)-[NH1]-[c]'))) == (3 if identifier == 'tabta' else 1)


@pytest.mark.parametrize('identifier', list(SELECTED))
def test_portable_selection_passes_importer_and_keeps_exact_review_qualifications(identifier):
    row = SELECTED[identifier]
    assert importer.validate_row(row) == row['canonical_smiles']
    for evidence in row['source_evidence']:
        assert evidence['ground_truth'] is False
        assert evidence['review_status'] == 'two_agent_source_graph_review'
        assert evidence['source_doi'] in row['source_dois']
        assert evidence['reviewed_page_images']
        assert all(re.fullmatch('[0-9a-f]{64}', p['sha256']) for p in evidence['reviewed_page_images'])
        assert 'human ground truth' in evidence['scope']


@pytest.mark.parametrize('mutation,reason', [
    (lambda r: r.update(ground_truth=True), 'provisional agent review'),
    (lambda r: r.update(reviewer_type='human'), 'provisional agent review'),
    (lambda r: r.update(formula='C1H1'), 'Formula disagrees'),
    (lambda r: r.update(canonical_smiles='CC'), 'Canonical finite graph mismatch'),
    (lambda r: r.update(automatic_prepare_supported=True), 'Automatic construction scope'),
    (lambda r: r['source_evidence'][0].update(ground_truth=True), 'two-agent graph review'),
    (lambda r: r['source_evidence'][0].update(source_doi='10.0000/unrelated'), 'Source DOI'),
    (lambda r: r['source_evidence'][0].update(candidate_selection_sha256='bad'), 'exact source/selection/review hash'),
    (lambda r: r.update(source_evidence_private_path='C:/private/source.pdf'), 'private path'),
])
def test_invalid_graph_identity_claims_and_source_bindings_are_rejected(mutation, reason):
    row = deepcopy(SELECTED['btba'])
    mutation(row)
    with pytest.raises(ValueError, match=reason):
        importer.validate_row(row)


@pytest.fixture
def proof_fixture(tmp_path):
    """Synthetic files exercise hash/receipt validation, not source admission."""
    row = deepcopy(SELECTED['btba'])
    pdf = tmp_path / 'synthetic-source.pdf'
    page = tmp_path / 'synthetic-page.png'
    selection = tmp_path / 'synthetic-candidates.json'
    receipt = tmp_path / 'synthetic-review.json'
    pdf.write_bytes(b'SYNTHETIC test source bytes; not a publication')
    page.write_bytes(b'SYNTHETIC test page bytes; not visual evidence')
    candidates = [{'candidate_id': row['id'], 'canonical_smiles': row['canonical_smiles']}]
    importer.write(selection, candidates)
    review = {'ground_truth': False, 'selection_sha256': digest(selection), 'reviews': [{
        'candidate_id': row['id'], 'canonical_smiles': row['canonical_smiles'],
        'result': 'finite_graph_agrees_with_viewed_source', 'viewed_pages': [{'sha256': digest(page)}],
    }]}
    importer.write(receipt, review)
    row['source_evidence'] = [{
        'source_pdf_sha256': digest(pdf), 'source_page_image_sha256': digest(page),
        'candidate_selection_sha256': digest(selection), 'independent_review_receipt_sha256': digest(receipt),
        'reviewed_page_images': [{'sha256': digest(page)}],
    }]
    return row, tmp_path, pdf, selection, receipt


def test_evidence_files_must_still_match_the_reviewed_bytes(proof_fixture):
    row, root, pdf, _, _ = proof_fixture
    assert importer.verify_evidence([row], root) == 4
    pdf.write_bytes(pdf.read_bytes() + b' changed after review')
    with pytest.raises(ValueError, match='pinned evidence files are absent or changed'):
        importer.verify_evidence([row], root)


@pytest.mark.parametrize('wrong_field,wrong_value', [('candidate_id', 'unrelated'), ('canonical_smiles', 'CC')])
def test_hash_valid_candidate_must_bind_exact_selected_identity_and_graph(proof_fixture, wrong_field, wrong_value):
    row, root, _, selection, receipt = proof_fixture
    candidates = read(selection)
    candidates[0][wrong_field] = wrong_value
    importer.write(selection, candidates)
    evidence = row['source_evidence'][0]
    evidence['candidate_selection_sha256'] = digest(selection)
    review = read(receipt)
    review['selection_sha256'] = digest(selection)
    importer.write(receipt, review)
    evidence['independent_review_receipt_sha256'] = digest(receipt)
    with pytest.raises(ValueError, match='Pinned candidate does not match admitted ID and graph'):
        importer.verify_evidence([row], root)


@pytest.mark.parametrize('change,reason', [
    (lambda r: r.update(selection_sha256='0' * 64), 'not bound to this candidate selection'),
    (lambda r: r.update(ground_truth=True), 'not bound to this candidate selection'),
    (lambda r: r['reviews'][0].update(canonical_smiles='CC'), 'does not admit this exact graph'),
    (lambda r: r['reviews'][0].update(result='needs_review'), 'does not admit this exact graph'),
    (lambda r: r['reviews'][0].update(viewed_pages=[{'sha256':'0' * 64}]), 'omits independently reviewed source pages'),
])
def test_hash_valid_receipt_still_requires_exact_graph_review_and_page_binding(proof_fixture, change, reason):
    row, root, _, _, receipt = proof_fixture
    review = read(receipt)
    change(review)
    importer.write(receipt, review)
    row['source_evidence'][0]['independent_review_receipt_sha256'] = digest(receipt)
    with pytest.raises(ValueError, match=reason):
        importer.verify_evidence([row], root)
