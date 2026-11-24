"""Chemical invariants for explicit drawing designs, not material predictions.

Finite inputs below are copied from the reviewed source library. Tests inspect
local bonds and hydrogen/functional-group counts; they deliberately do not
canonicalize or globally re-aromatize the assembled macrocycles.
"""
from collections import Counter
from copy import deepcopy

import pytest
from rdkit import Chem

from cof_vision.cofdraw_chemistry import (
    design_signature,
    inspect_monomer,
    prepare_template,
    refresh_design_identity,
)
from cof_vision.framework_drawing import VERSION, TP, build_patch, linkage_fragment
from cof_vision.cofdraw_repeat import repeat_options


# Catalog v3 source occurrences (provisional agent reviews, not human gold):
# TA: acs-wave3-001-006, DOI 10.1021/jacs.3c10312;
# Pa: acs-context-r187-002, DOI 10.1021/acsmacrolett.6b00805;
# Tp: acs-context-r089-004, DOI 10.1021/acsami.5c01437;
# DETH: wiley-wave3-r241-01, DOI 10.1002/anie.202319876.
TA = "O=Cc1ccc(C=O)cc1"
PA = "Nc1ccc(N)cc1"
DETH = "CCOc1cc(C(=O)NN)c(OCC)cc1C(=O)NN"

# Planar tetratopic cores: wiley-wave3-r251-05 / r122-01, respectively;
# DOI 10.1002/anie.201902543 / 10.1002/chem.202303556.
PY_CHO = (
    "O=Cc1ccc(-c2cc(-c3ccc(C=O)cc3)c3ccc4c(-c5ccc(C=O)cc5)"
    "cc(-c5ccc(C=O)cc5)c5ccc2c3c54)cc1"
)
PY_NH2 = (
    "Nc1ccc(-c2cc(-c3ccc(N)cc3)c3ccc4c(-c5ccc(N)cc5)"
    "cc(-c5ccc(N)cc5)c5ccc2c3c54)cc1"
)

# Tetrahedral controls: rsc-molecule-r048-02 / acs-wave3-001-004;
# DOI 10.1039/d0gc02670k / 10.1021/jacs.0c12505.
TFPM = "O=Cc1ccc(C(c2ccc(C=O)cc2)(c2ccc(C=O)cc2)c2ccc(C=O)cc2)cc1"
TAPM = "Nc1ccc(C(c2ccc(N)cc2)(c2ccc(N)cc2)c2ccc(N)cc2)cc1"

# Six-site triptycene, wiley-wave3-r255-04, DOI 10.1002/anie.202117101.
# This is NOT evidence of an hxl material or a planar conformation. The hxl
# test below checks an explicitly chosen projection's incidence/intent only.
HFPTP_H = (
    "O=Cc1ccc(-c2cc3c(cc2-c2ccc(C=O)cc2)C2c4cc(-c5ccc(C=O)cc5)"
    "c(-c5ccc(C=O)cc5)cc4C3c3cc(-c4ccc(C=O)cc4)"
    "c(-c4ccc(C=O)cc4)cc32)cc1"
)


def _design(a, b, linkage, template):
    recipe = prepare_template(a, b, linkage, "A", "B", 6, template)
    recipe.update(schema_version=VERSION, recipe_id="chemistry-regression", designed_variant=True)
    return recipe


def _count(mol, smarts):
    return len(mol.GetSubstructMatches(Chem.MolFromSmarts(smarts)))


def _assert_internal_bonds_retained(recipe, mol, audit):
    """Every surviving source-fragment bond keeps its order/aromatic flag."""
    origins = audit["atom_origin"]
    indices = {(o["instance"], o["fragment_atom_index"]): i for i, o in enumerate(origins)}
    expected_internal_bonds = 0
    for instance in recipe["instances"]:
        fragment = Chem.MolFromSmiles(recipe["components"][instance["component"]]["fragment_smiles"])
        for atom in fragment.GetAtoms():
            key = instance["id"], atom.GetIdx()
            if key in indices and atom.GetAtomicNum():
                actual = mol.GetAtomWithIdx(indices[key])
                assert actual.GetAtomicNum() == atom.GetAtomicNum()
                assert actual.GetFormalCharge() == atom.GetFormalCharge()
                assert actual.GetIsAromatic() == atom.GetIsAromatic()
                assert actual.GetTotalNumHs() == atom.GetTotalNumHs()
        for bond in fragment.GetBonds():
            left = indices.get((instance["id"], bond.GetBeginAtomIdx()))
            right = indices.get((instance["id"], bond.GetEndAtomIdx()))
            if left is not None and right is not None:
                expected_internal_bonds += 1
                actual = mol.GetBondBetweenAtoms(left, right)
                assert actual is not None
                assert actual.GetBondType() == bond.GetBondType()
                assert actual.GetIsAromatic() == bond.GetIsAromatic()
    internal = [b for b in mol.GetBonds() if origins[b.GetBeginAtomIdx()]["instance"] == origins[b.GetEndAtomIdx()]["instance"]]
    assert len(internal) == expected_internal_bonds


@pytest.mark.parametrize(
    "linkage,a,b,template,carbonyls,nn_bonds,n_atoms,n_hydrogens,join_order,carbon_hydrogens",
    [
        ("imine", TA, PA, "chain_2_2", 0, 0, 6, 0, 2, 1),
        ("beta_ketoenamine", TP, PA, "hcb_3_2", 18, 0, 12, 12, 1, 1),
        ("acylhydrazone", TA, DETH, "chain_2_2", 6, 6, 12, 6, 2, 1),
        ("keto_hydrazide", TP, DETH, "hcb_3_2", 30, 12, 24, 24, 1, 1),
        # Neutral hydrazine: acs-context-r364-003, DOI 10.1021/jacs.4c08231.
        ("azine", TA, "NN", "chain_2_2", 0, 3, 6, 0, 2, 1),
        # A selected reduced product state; not a one-step reaction claim.
        ("secondary_amine", TA, PA, "chain_2_2", 0, 0, 6, 6, 1, 2),
    ],
    ids=["imine", "ketoenamine", "acylhydrazone", "keto-hydrazide", "azine", "reduced-imine"],
)
def test_six_product_states_preserve_local_chemistry(
    linkage, a, b, template, carbonyls, nn_bonds, n_atoms, n_hydrogens, join_order, carbon_hydrogens
):
    recipe = _design(a, b, linkage, template)
    mol, audit = build_patch(recipe)
    mol.UpdatePropertyCache(strict=True)
    assert len(Chem.GetMolFrags(mol)) == 1
    assert _count(mol, "[C]=[O]") == carbonyls
    assert _count(mol, "[#7]-[#7]") == nn_bonds
    nitrogens = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7]
    assert len(nitrogens) == n_atoms
    assert sum(atom.GetTotalNumHs() for atom in nitrogens) == n_hydrogens
    assert sum(atom.GetFormalCharge() for atom in mol.GetAtoms()) == 0

    origins = audit["atom_origin"]
    joins = [bond for bond in mol.GetBonds() if origins[bond.GetBeginAtomIdx()]["instance"] != origins[bond.GetEndAtomIdx()]["instance"]]
    assert len(joins) == len(recipe["connections"])
    for bond in joins:
        atoms = [bond.GetBeginAtom(), bond.GetEndAtom()]
        assert sorted(atom.GetAtomicNum() for atom in atoms) == [6, 7]
        assert bond.GetBondTypeAsDouble() == join_order
        carbon = next(atom for atom in atoms if atom.GetAtomicNum() == 6)
        assert carbon.GetTotalNumHs() == carbon_hydrogens
        assert bond.GetStereo() == Chem.BondStereo.STEREONONE
    _assert_internal_bonds_retained(recipe, mol, audit)
    assert audit["whole_cof_smiles"] is None
    assert audit["fragment_smiles_canonical"] is False
    assert not audit["synthetic_training_admitted"]


@pytest.mark.parametrize("linkage,terminal_h", [("acylhydrazone", 0), ("keto_hydrazide", 1)])
def test_only_terminal_hydrazide_n_reacts_and_acyl_nh_survives(linkage, terminal_h):
    # DETH has two terminal NH2 groups and two distinct acyl NH groups.
    fragment = linkage_fragment(DETH, linkage, "acylhydrazide")
    mol = Chem.MolFromSmiles(fragment["fragment_smiles"])
    assert fragment["connectivity"] == 2
    assert _count(mol, "[C](=[O])-[NH1]-[#7]") == 2
    assert _count(mol, "[C]=[O]") == 2
    for port in (a for a in mol.GetAtoms() if not a.GetAtomicNum()):
        terminal = port.GetNeighbors()[0]
        assert terminal.GetAtomicNum() == 7
        assert terminal.GetTotalNumHs() == terminal_h
        assert sorted(n.GetAtomicNum() for n in terminal.GetNeighbors()) == [0, 7]


@pytest.mark.parametrize("smiles", ["NC(=O)c1ccccc1", "CNC(=O)c1ccccc1", "CNNC(=O)c1ccccc1"])
def test_amide_or_n_substituted_hydrazide_is_not_terminal_hydrazide(smiles):
    roles = inspect_monomer(smiles)["compatible_roles"]
    assert not any(role["role"] == "acylhydrazide" for role in roles)
    for linkage in ("acylhydrazone", "keto_hydrazide"):
        with pytest.raises(ValueError, match="No supported reactive sites"):
            linkage_fragment(smiles, linkage, "acylhydrazide")


@pytest.mark.parametrize("smiles", ["NNC(=O)c1ccccc1", "NN.O", "[NH3+]N.[Cl-]"])
def test_azine_does_not_silently_strip_hydrazide_hydrate_or_salt(smiles):
    with pytest.raises(ValueError, match="exact neutral hydrazine"):
        linkage_fragment(smiles, "azine", "hydrazine")


@pytest.mark.parametrize("a,b", [(TFPM, PA), (TA, TAPM)])
def test_tetrahedral_four_site_nodes_cannot_imply_planar_sql(a, b):
    assert "tetrahedral" in [inspect_monomer(s)["geometry_hint"] for s in (a, b)]
    with pytest.raises(ValueError, match="tetrahedral.*3D incidence"):
        _design(a, b, "imine", "sql_4_2")
    projection = _design(a, b, "imine", "alternating_hexagon")
    assert projection["design_topology"] is None
    assert " COFex.unk.cat0;" in projection["design_cofex_notation"]


@pytest.mark.parametrize(
    "template,a,b,counts,instances,joins,ports,net",
    [
        ("sql_4_2", PY_CHO, PA, [4, 2], 8, 8, 8, "sql"),
        ("sql_4_2", TA, PY_NH2, [2, 4], 8, 8, 8, "sql"),
        ("sql_4_4", PY_CHO, PY_NH2, [4, 4], 4, 4, 8, "sql"),
        ("hxl_6_2", HFPTP_H, PA, [6, 2], 6, 6, 12, "hxl"),
        ("chain_2_2", TA, PA, [2, 2], 6, 5, 2, None),
    ],
    ids=["sql-4-2", "sql-2-4", "sql-4-4", "chosen-hxl-projection", "open-chain"],
)
def test_template_counts_and_ports_are_design_intent_not_source_identity(
    template, a, b, counts, instances, joins, ports, net
):
    recipe = _design(a, b, "imine", template)
    assert recipe["compatibility"]["reactive_site_counts"] == counts
    assert len(recipe["instances"]) == instances
    assert len(recipe["connections"]) == joins
    assert recipe["selected_template"]["periodic_topology_inferred"] is False
    assert recipe["design_topology"] == net
    assert f" COFex.{net or 'unk'}.cat0;" in recipe["design_cofex_notation"]
    assert "assumed" in recipe["design_identity_qualification"]
    assert not recipe.get("cofex_id") and not recipe.get("cofex_notation")

    mol, audit = build_patch(recipe)
    mol.UpdatePropertyCache(strict=True)
    assert len(Chem.GetMolFrags(mol)) == 1
    assert len(audit["continuation_ports"]) == ports
    assert sum(not atom.GetAtomicNum() for atom in mol.GetAtoms()) == ports
    assert audit["building_block_instances"] == instances
    assert audit["connections"] == joins
    assert audit["topology"] is None
    assert audit["whole_cof_smiles"] is None
    _assert_internal_bonds_retained(recipe, mol, audit)
    if template == "hxl_6_2":
        # The nonplanar triptycene in this projection cannot tile a flat
        # triangular net while preserving its attachment directions.
        choices = repeat_options(dict(recipe,design_template_signature=design_signature(recipe)))
        assert not choices['available'] and len(choices['extents']) == 1
        assert 'align' in choices['reason']
    if template == "chain_2_2":
        assert recipe["layout_kind"] == "open_chain"
        endpoints = Counter(p["instance"] for p in audit["continuation_ports"])
        assert endpoints == {recipe["instances"][0]["id"]: 1, recipe["instances"][-1]["id"]: 1}


def test_changed_attachment_incidence_clears_chosen_net_identity():
    recipe = _design(PY_CHO, PA, "imine", "sql_4_2")
    recipe["design_template_signature"] = design_signature(recipe)
    refresh_design_identity(recipe)
    assert recipe["design_topology"] == "sql"
    changed = deepcopy(recipe)
    changed["instances"][0]["incoming_port"], changed["instances"][0]["outgoing_port"] = (
        changed["instances"][0]["outgoing_port"], changed["instances"][0]["incoming_port"]
    )
    refresh_design_identity(changed)
    assert changed["selected_template"] is None
    assert changed["design_topology"] is None
    assert " COFex.unk.cat0;" in changed["design_cofex_notation"]


def test_library_functional_groups_are_independent_of_reaction_sites():
    from cof_vision.cofdraw_chemistry import decorated_monomers

    rows = decorated_monomers([{"id": "tp", "smiles": TP},
                               {"id": "hydrazide", "smiles": "NNC(=O)c1ccc(C(=O)NN)cc1"},
                               {"id": "ester", "smiles": "COC(=O)c1ccccc1"}])
    tp, hydrazide, ester = rows
    assert {"Phenol", "Aldehyde"} <= set(tp["functional_groups"])
    assert "Phenol" not in {g["label"] for g in tp["functionality"]}
    assert "Acyl hydrazide" in hydrazide["functional_groups"]
    assert "Primary amine" not in hydrazide["functional_groups"]
    assert set(ester["functional_groups"]) == {"Ester"}


def test_library_linkage_groups_across_reviewed_monomer_types():
    from pathlib import Path

    from cof_vision.cofdraw_chemistry import decorated_monomers
    from cof_vision.cofdraw_core import catalog_monomers, read_catalog

    expected = {
        'tp': {'aldehyde'}, 'dmtp': {'aldehyde'}, 'bpda_n3': {'aldehyde'},
        'pdan': {'activated_cyanomethyl'}, 'hhtp': {'catechol'}, 'tmc': {'acyl_chloride'},
        'bfbaepy': {'aldehyde', 'primary_amine'},
        'thpe': {'phenol'}, 'thpb': {'phenol'}, 'dpt': {'activated_ester'}, 'dpbp': {'activated_ester'},
        'tfpc_dioxin': {'vicinal_aryl_fluorides'}, 'ppcn': {'nitrile'},
        'recovered_monoformylphloroglucinol': {'phenol'},
        'recovered_1_4_phenylene_diisocyanate': {'isocyanate'}, 'recovered_urea': {'terminal_urea'},
        'cv_b1cb056e6150b738': {'benzyl_bromide'}, 'cv_2a97fff010220244': {'pyridyl_nitrogen'},
        'cv_58b48601a6c68864': set(),  # S6, an intermediate with free amino groups.
        'cv_05cbf00025a2d4a1': {'primary_amine', 'masked_aldehyde'},
        'cv_6bb0458f6411d793': {'aldehyde', 'masked_amine'},
        'terphenyl_cooh2_cho2': {'aldehyde', 'carboxylic_acid'},
    }
    catalog = read_catalog(Path(__file__).resolve().parents[1] / 'data')
    feeds = [r for r in catalog_monomers(catalog) if r['id'] in expected]
    by_id = {r['id']: r for r in decorated_monomers(feeds)}
    assert {k: {g['id'] for g in r['linkage_functional_groups']} for k, r in by_id.items()} == expected

    for row in by_id.values():
        candidates = {(r['linkage'], r['role']) for r in row['binary_roles']}
        facet_roles = {(r['linkage'], r['role']) for g in row['linkage_functional_groups'] for r in g['binary_roles']}
        assert candidates == facet_roles
    groups = {g['id']: g for g in by_id['terphenyl_cooh2_cho2']['linkage_functional_groups']}
    assert ('imine', 'aryl_aldehyde') in {(r['linkage'], r['role']) for r in groups['aldehyde']['binary_roles']}
    assert groups['carboxylic_acid']['binary_roles'] == [{'linkage': 'amide', 'role': 'aryl_carboxyl'}]
    assert by_id['dpt']['linkage_functional_groups'][0]['binary_roles'] == []


def test_linkage_facet_requires_a_reaction_context_for_reference_only_motifs():
    from cof_vision.cofdraw_chemistry import decorated_monomers

    unknown = {'smiles': 'C#Cc1ccc(C#C)cc1', 'catalog_use': 'reference_only',
               'source_dois': ['10.example/unassigned']}
    same_paper = dict(unknown, source_dois=['10.1021/acsapm.4c00560'])
    click = dict(same_paper, reported_source_chemistry=['azide_alkyne_cycloaddition'])
    unknown_row, same_paper_row, click_row = decorated_monomers([unknown, same_paper, click])
    assert unknown_row['functional_groups'] == ['Alkyne']
    assert not unknown_row['linkage_functional_groups']
    assert not same_paper_row['linkage_functional_groups']
    assert click_row['linkage_functional_groups'] == [
        {'id': 'terminal_alkyne', 'label': 'Terminal alkyne', 'binary_roles': []}]
