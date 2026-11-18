"""Independent runtime checks for chemically matched periodic hypotheses."""
from pathlib import Path
import hashlib
import json

from rdkit import Chem

from cof_vision.cofdraw_network import inspect_periodic_network
from cof_vision.cofdraw_science import import_cif

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT/'data/crystal-models'


def test_match_manifest_hashes_and_exact_precursor_scope():
    manifest = json.loads((DATA/'manifest.json').read_text())
    ref = manifest['precursor_mapping_audit']
    assert hashlib.sha256((DATA/ref['path']).read_bytes()).hexdigest() == ref['sha256']
    audit = json.loads((DATA/ref['path']).read_text())
    assert len(audit['templates']) == len(manifest['templates']) == 5
    proofs = {r['template_id']:r for r in audit['templates']}
    catalog = json.loads((ROOT/'data/framework-drawing/catalog/current.json').read_text(encoding='utf8'))
    monomers = {m['id']:m for m in catalog['monomer_presets']}
    for row in manifest['templates']:
        assert hashlib.sha256((DATA/row['path']).read_bytes()).hexdigest() == row['sha256'] == proofs[row['id']]['cif_sha256']
        binding = row['precursor_match']
        assert binding['linkage'] == 'imine'
        assert binding['canonical_smiles'] == [Chem.MolToSmiles(Chem.MolFromSmiles(monomers[k]['smiles']),isomericSmiles=True) for k in binding['monomer_ids']]
        assert binding['preserves_drawing_coordinates'] is False
        assert proofs[row['id']]['all_product_fragments_match']
        assert proofs[row['id']]['all_fragment_cycle_translations_zero']
        if row['id'] == 'dia-imine':
            assert binding['drawing_template_ids'] == [] and binding['requires_explicit_topology_choice']
    assert 'beta-ketoenamine' in audit['excluded_beta_ketoenamine_candidate']['reason']


def _adjacency_signature(elements, edges):
    # RDKit canonical graph signature is an independent checker from the
    # generator's NetworkX isomorphism. Bond orders are deliberately excluded:
    # the artifact is a coordinate CIF and this check is of explicit-H adjacency.
    rw = Chem.RWMol()
    for element in elements:
        atom = Chem.Atom(element); atom.SetNoImplicit(True); rw.AddAtom(atom)
    for a,b in edges:
        rw.AddBond(a,b,Chem.BondType.SINGLE)
    mol = rw.GetMol(); mol.UpdatePropertyCache(strict=False)
    return Chem.MolToSmiles(mol, canonical=True)


def _precursor_fragment(smiles, role):
    mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    remove = set()
    if role == 'aldehyde':
        remove = {h[1] for h in mol.GetSubstructMatches(Chem.MolFromSmarts('[CH1]=O'))}
    else:
        for atom in mol.GetAtoms():
            if atom.GetSymbol() == 'N':
                remove.update(a.GetIdx() for a in atom.GetNeighbors() if a.GetSymbol()=='H')
    selected = [a.GetIdx() for a in mol.GetAtoms() if a.GetIdx() not in remove]
    mapping = {old:new for new,old in enumerate(selected)}
    edges = [(mapping[b.GetBeginAtomIdx()],mapping[b.GetEndAtomIdx()]) for b in mol.GetBonds() if b.GetBeginAtomIdx() in mapping and b.GetEndAtomIdx() in mapping]
    return _adjacency_signature([mol.GetAtomWithIdx(i).GetSymbol() for i in selected], edges)


def test_tfb_tapb_cell_contains_exact_full_precursor_fragments():
    model = import_cif((DATA/'hcb-tfb-tapb-aa.cif').read_text())
    net = inspect_periodic_network(model)
    elements = [a['element'] for a in model['atoms']]
    adjacency = [set() for _ in elements]
    for edge in net['edges']:
        a,b = edge['a'],edge['b']; adjacency[a].add(b); adjacency[b].add(a)
    cuts=[]
    for i,element in enumerate(elements):
        if element=='N':
            assert sorted(elements[j] for j in adjacency[i]) == ['C','C']
            carbons = [j for j in adjacency[i] if sorted(elements[k] for k in adjacency[j]) == ['C','H','N']]
            assert len(carbons)==1; cuts.append((i,carbons[0]))
    assert len(cuts)==3
    for a,b in cuts:adjacency[a].remove(b);adjacency[b].remove(a)
    unseen=set(range(len(elements)));signatures=[]
    while unseen:
        component={min(unseen)};queue=list(component)
        while queue:
            a=queue.pop()
            for b in adjacency[a]-component:component.add(b);queue.append(b)
        unseen-=component;indices=sorted(component);mapping={old:new for new,old in enumerate(indices)}
        edges=[(mapping[a],mapping[b]) for a in indices for b in adjacency[a] if a<b]
        signatures.append(_adjacency_signature([elements[a] for a in indices],edges))
    tfb='O=Cc1cc(C=O)cc(C=O)c1'
    tapb='Nc1ccc(-c2cc(-c3ccc(N)cc3)cc(-c3ccc(N)cc3)c2)cc1'
    assert sorted(signatures)==sorted([_precursor_fragment(tfb,'aldehyde'),_precursor_fragment(tapb,'amine')])
