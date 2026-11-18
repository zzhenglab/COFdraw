"""Generate a TFB/TAPB cell and audit exact precursor mapping of bundled CIFs.

Development only: pycofbuilder==0.0.9.4, pymatgen, networkx, numpy and RDKit.
Run after generate_crystal_templates.py if regenerating the original four CIFs.
The runtime consumes the resulting CIFs/manifest, not these dependencies.
"""
from collections import Counter
from pathlib import Path
import hashlib
import importlib.metadata
import json
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src'))

from rdkit import Chem
from cof_vision.cofdraw_science import import_cif
from cof_vision.cofdraw_network import inspect_periodic_network

DEST = ROOT / 'data/crystal-models'
NEW_ID = 'hcb-tfb-tapb-aa'
NEW_MODEL = 'T3_BENZ_CHO_H-T3_TPBZ_NH2_H-HCB-AA'
PAIRS = {
    'hcb-imine-aa': ('tfb', 'pa', 'hcb_3_2', 'hcb'),
    'hcb-imine-ab': ('tfb', 'pa', 'hcb_3_2', 'hcb'),
    'hcb-expanded-aa': ('ta', 'tapb', 'hcb_3_2', 'hcb'),
    'dia-imine': ('ta', 'tapm', None, 'dia'),
    NEW_ID: ('tfb', 'tapb', 'hcb_3_3', 'hcb'),
}


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _fragment(smiles, role):
    """Expected explicit-H precursor fragment after imine water elimination."""
    import networkx as nx
    mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    if role == 'aldehyde':
        remove = {hit[1] for hit in mol.GetSubstructMatches(Chem.MolFromSmarts('[CH1]=[O]'))}
    else:
        nitrogens = [hit[1] for hit in mol.GetSubstructMatches(Chem.MolFromSmarts('[c][NH2]'))]
        remove = {n.GetIdx() for i in nitrogens for n in mol.GetAtomWithIdx(i).GetNeighbors() if n.GetSymbol() == 'H'}
    graph = nx.Graph()
    for atom in mol.GetAtoms():
        if atom.GetIdx() not in remove:
            graph.add_node(atom.GetIdx(), element=atom.GetSymbol())
    for bond in mol.GetBonds():
        a, b = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        if a in graph and b in graph:
            graph.add_edge(a, b)
    return graph


def validate_imine_cell(text, aldehyde_smiles, amine_smiles):
    """Verify all atoms/adjacencies and explicit H against precursor fragments.

    CIF has coordinates, not bond orders. The independent geometric graph is
    cut only at N–CH(aryl) imine connections. Every resulting finite component
    must match one full precursor after the stated water elimination. The
    matched signatures include explicit H, so enol-imine cannot pass as keto.
    """
    import networkx as nx
    model = import_cif(text)
    network = inspect_periodic_network(model)
    if network['short_contact_count'] or network['isolated_atom_indices']:
        raise ValueError('Invalid short contacts or isolated atoms in candidate cell')
    if network['threshold_sensitivity']['classification_changed']:
        raise ValueError('Periodic connectivity is cutoff-sensitive')
    graph = nx.Graph()
    for index, atom in enumerate(model['atoms']):
        graph.add_node(index, element=atom['element'])
    for edge in network['edges']:
        a, b = edge['a'], edge['b']
        if a == b or graph.has_edge(a, b):
            raise ValueError('Ambiguous quotient adjacency; larger explicit cell required for fragment audit')
        graph.add_edge(a, b, translation=edge['translation'], distance=edge['distance_angstrom'])
    cuts = []
    for index, row in graph.nodes(data=True):
        if row['element'] != 'N':
            continue
        if sorted(graph.nodes[n]['element'] for n in graph[index]) != ['C', 'C']:
            raise ValueError('Imine cell nitrogen must have exactly two carbon neighbors and no hydrogen')
        candidates = [c for c in graph[index] if sorted(graph.nodes[n]['element'] for n in graph[c]) == ['C', 'H', 'N']]
        if len(candidates) != 1:
            raise ValueError('Imine N–CH junction is not uniquely identified')
        carbon = candidates[0]
        if not 1.15 < graph.edges[index, carbon]['distance'] < 1.38:
            raise ValueError('Expected short imine C=N coordinate distance is absent')
        cuts.append((index, carbon))
    if not cuts:
        raise ValueError('No imine junctions present')
    for a, b in cuts:
        graph.remove_edge(a, b)
    targets = {'aldehyde': _fragment(aldehyde_smiles, 'aldehyde'), 'amine': _fragment(amine_smiles, 'amine')}
    counts = Counter()
    match = nx.algorithms.isomorphism.categorical_node_match('element', None)
    for component in nx.connected_components(graph):
        subgraph = graph.subgraph(component)
        # Quotient adjacency alone could confuse a winding periodic chain with
        # a finite ring. Require every fragment cycle to have zero lattice
        # translation before comparing it with a molecular precursor graph.
        seed = next(iter(component)); potential = {seed: (0, 0, 0)}; queue = [seed]
        while queue:
            a = queue.pop()
            for b in subgraph[a]:
                t = subgraph.edges[a, b]['translation']
                shift = t if a < b else [-v for v in t]
                predicted = tuple(x+y for x, y in zip(potential[a], shift))
                if b in potential:
                    if potential[b] != predicted:
                        raise ValueError('A proposed molecular fragment winds through the periodic lattice')
                else:
                    potential[b] = predicted; queue.append(b)
        matches = [role for role, target in targets.items() if nx.is_isomorphic(subgraph, target, node_match=match)]
        if len(matches) != 1:
            raise ValueError('Periodic product fragment does not match the specified full precursor graph')
        counts[matches[0]] += 1
    # Independent atom balance of full precursors minus one H2O per imine.
    expected = Counter()
    for role, smiles in [('aldehyde', aldehyde_smiles), ('amine', amine_smiles)]:
        expected.update({element: n * counts[role] for element, n in Counter(a.GetSymbol() for a in Chem.AddHs(Chem.MolFromSmiles(smiles)).GetAtoms()).items()})
    expected.subtract({'H': 2*len(cuts), 'O': len(cuts)})
    actual = Counter(a['element'] for a in model['atoms'])
    if +expected != actual or any(v < 0 for v in expected.values()):
        raise ValueError('Condensation stoichiometry mismatch')
    return {'atom_count': len(model['atoms']), 'element_counts': dict(sorted(actual.items())),
            'precursor_multiplicities': dict(counts), 'water_eliminations': len(cuts),
            'all_product_fragments_match': True, 'explicit_hydrogens_validated': True,
            'all_fragment_cycle_translations_zero': True,
            'inferred_periodic_dimensions': network['dimensions'], 'short_contact_count': 0,
            'minimum_candidate_distance_angstrom': network['minimum_candidate_distance_angstrom'],
            'method': 'All-element explicit-H adjacency isomorphism after cutting uniquely identified short N–CH imine junctions, plus complete precursor-minus-water atom balance.',
            'bond_order_limit': 'CIF does not contain chemical bond orders; coordinate-inferred adjacency and explicit H are checked against the chosen imine product. No experimental structure assignment.'}


def main():
    import numpy as np
    import pycofbuilder as pcb
    from pymatgen.core import Structure
    from pymatgen.io.cif import CifWriter
    if importlib.metadata.version('pycofbuilder') != '0.0.9.4':
        raise ValueError('Use the reviewed pyCOFBuilder version 0.0.9.4')
    catalog = json.loads((ROOT/'data/framework-drawing/catalog/current.json').read_text(encoding='utf8'))
    monomers = {r['id']: r for r in catalog['monomer_presets']}
    framework = pcb.Framework(NEW_MODEL, save_bb=False, log_level='error', out_path=str(ROOT/'build/model-generator'))
    crystal = Structure(framework.cellMatrix, framework.atom_types, framework.atom_pos, coords_are_cartesian=True)
    text = '# COFdraw unrelaxed TFB/TAPB reticular hypothesis; not an experimental structure assignment.\n' + str(CifWriter(crystal, symprec=None))
    # Validate before changing the bundled file or metadata.
    validate_imine_cell(text, monomers['tfb']['smiles'], monomers['tapb']['smiles'])
    (DEST/(NEW_ID+'.cif')).write_text(text, encoding='utf8', newline='\n')
    manifest = json.loads((DEST/'manifest.json').read_text(encoding='utf8'))
    rows = [r for r in manifest['templates'] if r['id'] != NEW_ID]
    distances = crystal.distance_matrix.copy(); np.fill_diagonal(distances, np.inf)
    rows.append({'id': NEW_ID, 'label': 'TFB + TAPB imine · AA layers', 'model_name': NEW_MODEL,
                 'path': NEW_ID+'.cif', 'sha256': sha(DEST/(NEW_ID+'.cif')), 'atoms': len(crystal),
                 'formula': crystal.composition.formula, 'minimum_atom_distance_angstrom': float(distances.min())})
    audited = []
    for row in rows:
        aldehyde, amine, drawing_template, topology = PAIRS[row['id']]
        path = DEST/row['path']; assert sha(path) == row['sha256']
        smiles = [Chem.MolToSmiles(Chem.MolFromSmiles(monomers[k]['smiles']), isomericSmiles=True) for k in (aldehyde, amine)]
        validation = validate_imine_cell(path.read_text(encoding='utf8'), *smiles)
        if validation['inferred_periodic_dimensions'] != ([3] if topology == 'dia' else [2]):
            raise ValueError('Unexpected periodic graph dimension')
        row['precursor_match'] = {'monomer_ids': [aldehyde, amine], 'canonical_smiles': smiles,
                                  'roles': ['aryl_aldehyde', 'primary_arylamine'], 'linkage': 'imine',
                                  'drawing_template_ids': [drawing_template] if drawing_template else [],
                                  'generated_topology': topology, 'requires_explicit_topology_choice': drawing_template is None,
                                  'preserves_drawing_coordinates': False,
                                  'qualification': 'Exact precursor pair and chosen imine chemistry; an independent unrelaxed periodic hypothesis, not conversion of finite drawing coordinates or proof of source/crystal identity.'}
        row['description'] = f"{row['model_name']}. Full-cell explicit-H reticular hypothesis from the stated precursor pair. Topology and stacking are chosen; no experimental assignment or geometry optimization."
        audited.append({'template_id': row['id'], 'cif_sha256': row['sha256'], 'precursor_match': row['precursor_match'], **validation})
    library = Path(pcb.__file__).parent
    parts = ['framework.py','tools.py','data/topology.py','data/core/T3/BENZ.cjson','data/core/T3/TPBZ.cjson','data/conector/CHO.cjson','data/conector/NH2.cjson','data/func_groups/H.cjson']
    proof = {'schema_version': 1, 'generator': 'pyCOFBuilder 0.0.9.4', 'generator_library_files': [{'path': p, 'sha256': sha(library/p)} for p in parts],
             'templates': audited, 'ground_truth': False,
             'precursor_smiles_limit': 'pyCOFBuilder 0.0.9.4 emits unresolved [R] placeholders and malformed leading connector SMILES; these strings are not treated as chemical identity. Identity is checked independently from all generated coordinate fragments against the catalog full precursor graphs.',
             'excluded_beta_ketoenamine_candidate': {'model_name':'T3_BENZ_CHO_OH-L2_BENZ_NH2_H-HCB_A-AA',
                    'reason':'The stock builder retains six O–H bonds and no N–H bonds in its 72-atom Tp/Pa cell. This is enol-imine, not the beta-ketoenamine product. No beta-ketoenamine match is provided.'}}
    audit_path = DEST/'precursor-mapping.audit.json'; audit_path.write_text(json.dumps(proof,indent=2)+'\n',encoding='utf8',newline='\n')
    manifest.update(templates=rows, precursor_mapping_audit={'path': audit_path.name, 'sha256': sha(audit_path)},
                    matched_model_generator='scripts/generate_matched_crystal_templates.py')
    (DEST/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n',encoding='utf8',newline='\n')
    print(json.dumps({'template_count':len(rows),'new_template':NEW_ID,'new_cif_sha256':sha(DEST/(NEW_ID+'.cif')),'audit_sha256':sha(audit_path)}))


if __name__ == '__main__':
    main()
