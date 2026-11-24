"""Monomer library facets for groups that supply framework-forming sites."""
from rdkit import Chem

from .framework_drawing import CARBON_SITES


_ROLE_GROUP = {
    'aryl_aldehyde': 'aldehyde', 'carbon_aldehyde': 'aldehyde', 'tp_trialdehyde': 'aldehyde',
    'primary_arylamine': 'primary_amine', 'primary_alkylamine': 'primary_amine',
    'acylhydrazide': 'acyl_hydrazide', 'hydrazine': 'hydrazine',
    'aryl_hydroxylamine': 'hydroxylamine', 'aryl_acetonitrile': 'activated_cyanomethyl',
    'methyltriazine': 'methyltriazine', 'aryl_cyclic_anhydride': 'cyclic_anhydride',
    'aryl_vicinal_diol': 'catechol', 'aryl_boronic_acid': 'boronic_acid',
}
_CARBOXYL = Chem.MolFromSmarts(CARBON_SITES['aryl_carboxyl'])

_LABELS = {
    'activated_cyanomethyl': 'Activated cyanomethyl (–CH₂–CN)',
    'methyltriazine': 'Triazine-bound methyl',
    'catechol': 'Catechol (aromatic ortho-diol)',
    'masked_aldehyde': 'Masked aldehyde (acetal)',
    'masked_amine': 'Masked primary amine (Boc)',
    'benzophenone_imine': 'Masked primary amine (benzophenone imine)',
    'pinacol_boronate': 'Pinacol boronate ester',
    'lactone_methylene': 'Activated methylene (lactone)',
    'vicinal_aryl_fluorides': 'Adjacent aryl fluorides',
    'activated_ester': 'Activated pyridyl ester',
    'benzyl_alcohol': 'Benzylic alcohol',
    'vicinal_diol': 'Vicinal diol (aliphatic)',
    'chlorophosphazene': 'Chlorophosphazene (P–Cl)',
}

# These source-only motifs need both a reported reaction family and the atoms
# it consumes. For example, F reacts in the dioxin feed while CN is retained.
_REPORTED_GROUPS = {
    'borosilicate': {'silanol': '[Si]-[OH1]'},
    'polyphenylene_aldol_cyclotrimerization': {'aryl_methyl_ketone': 'c-C(=O)-[CH3]'},
    'dioxin': {'vicinal_aryl_fluorides': 'F-c:c-F', 'catechol': '[OH1]-c:c-[OH1]'},
    'ester_transesterification': {'phenol': 'c-[OH1]', 'activated_ester': 'C(=O)Oc1ccccn1'},
    'urea_condensation': {'terminal_urea': '[NH2]C(=O)[NH]'},
    'protected_amine_deprotection_then_imine': {'primary_amine': 'c-[NH2]', 'masked_amine': 'c-NC(=O)OC(C)(C)C'},
    'hybrid_titanium_cluster_imine': {'primary_amine': 'c-[NH2]'},
    'azide_alkyne_cycloaddition': {'azide': '[#6]-[$(N=[N+]=[N-]),$([N-]-[N+]#N)]', 'terminal_alkyne': 'C#[CH]'},
    'porphyrazine': {'nitrile': 'N#C-c:c-C#N'},
    'imine_boroxine': {'aldehyde': 'c-[CH1]=O', 'boronic_acid': 'c-B([OH1])[OH1]'},
    'benzophenone_imine_activation': {'benzophenone_imine': 'c-N=C(-c)-c'},
    'pinacol_boronate_activation': {'pinacol_boronate': 'c-B1OC(C)(C)C(C)(C)O1'},
    'lactone_aldol': {'lactone_methylene': 'c1[CH2]C(=O)Oc1'},
    'mannich': {'aryl_methyl_ketone': 'c-C(=O)-[CH3]', 'primary_amine': 'c-[NH2]', 'aldehyde': 'c-[CH1]=O'},
    'urethane': {'isocyanate': '[#6]-N=C=O', 'phenol': 'c-[OH1]'},
    'pyridinium_alkylation': {'benzyl_bromide': 'c-[CH2]-Br', 'pyridyl_nitrogen': 'n1ccccc1'},
    'aryl_bromide_coupling': {'aryl_bromide': 'c-Br'},
    'terminal_alkyne_coupling': {'terminal_alkyne': 'C#[CH]'},
    'disulfide': {'thiol': '[#6]-[SH1]'},
    'nitrile_cyclotrimerization': {'nitrile': '[#6]-C#N'},
    'phosphazene_substitution': {'chlorophosphazene': 'Cl-P=N', 'primary_amine': 'c-[NH2]'},
    'oxidative_alcohol_amidine': {'benzyl_alcohol': 'c-[CH2]-[OH1]', 'amidine': 'c-C(=[NH])-[NH2]'},
    'aliphatic_vicinal_diol_boronate': {'vicinal_diol': '[OH1]-[CX4]-[CX4]-[OH1]',
                                         'boronic_acid': 'c-B([OH1])[OH1]'},
    'mixed_imine_acylhydrazone': {'primary_amine': 'c-[NH2]', 'acyl_hydrazide': '[NH2][NH]C=O'},
}
_REPORTED_QUERIES = {family: {name: Chem.MolFromSmarts(pattern) for name, pattern in groups.items()}
                     for family, groups in _REPORTED_GROUPS.items()}
# Some earlier catalog admissions record the primary source without a family
# field. Match its documented feed and family before checking the molecular graph.
_SOURCE_FAMILIES = {
    '10.1021/acsapm.4c00560': ('azide_alkyne_cycloaddition', {
        'cv_7117e3cd0c229ae6', 'cv_cbe0d434259731af',
        'cv_e95042e72d14d974'}),
    '10.1002/anie.202203250': ('porphyrazine', {'ppcn'}),
    '10.1021/ja510926w': ('imine_boroxine', {'cv_ecb3aaeac976ee9e'}),
    '10.1021/jacs.7b00925': ('mixed_imine_acylhydrazone', {'aph', 'ash'}),
    '10.1021/acs.chemmater.2c01401': ('benzophenone_imine_activation', {'cv_46a108241ae9d94d'}),
    '10.1021/jacs.9b08147': ('benzophenone_imine_activation', {'cv_eba16000faee0e1a'}),
    '10.1039/c2ce26343b': ('pinacol_boronate_activation', {'cv_7758080fa6dc23dc'}),
    '10.1021/acs.nanolett.4c00343': ('lactone_aldol', {'cv_822e07eb3181998d'}),
    '10.1039/d4sc04358h': ('mannich', {'cv_7c59097771241e56'}),
    '10.1039/c6cc02578a': ('urethane', {'cv_1660a9e85fdbb415'}),
    '10.1039/d1qm00231g': ('urea_condensation', {'recovered_urea'}),
    '10.1039/d2dt04057c': ('urethane', {
        'recovered_monoformylphloroglucinol', 'recovered_1_4_phenylene_diisocyanate'}),
    '10.1039/c4ra10594j': ('pyridinium_alkylation', {
        'cv_b1cb056e6150b738', 'cv_2a97fff010220244',
        'cv_b5da5378ba04f5a3', 'cv_a64ea7f70834dd08'}),
    '10.1021/acs.chemmater.0c03901': ('aryl_bromide_coupling', {'cv_67660b54a57b8322'}),
    '10.1021/acs.chemmater.1c02564': ('terminal_alkyne_coupling', {'cv_ffbe3ec526ed840f'}),
    '10.1039/c7ra04775d': ('disulfide', {'cv_f38a5fec3421bf57', 'cv_151f80e20fa1712b'}),
    '10.1021/acsanm.0c00155': ('nitrile_cyclotrimerization', {
        'cv_8fa21acdb6eb6828', 'cv_3efc00770386bcb6',
        'cv_f0b622349d0f7705'}),
    '10.1021/acsami.7b01783': ('nitrile_cyclotrimerization', {
        'cv_d2f14c65352cc54b', 'cv_f6bd8bed9a6dc1e5',
        'cv_6ebf2b0bfae4bb95', 'cv_0cd5203f029ba8cd',
        'cv_6867529ef0532007'}),
    '10.1039/c8ta11722e': ('nitrile_cyclotrimerization', {'cv_b1a1c700ee34de1a'}),
    '10.1021/acs.chemmater.8b05365': ('nitrile_cyclotrimerization', {'cv_cd7462f82d3942a1'}),
    '10.1039/d1nj01456k': ('phosphazene_substitution', {'cv_beeda257ae9686b9'}),
    '10.1002/anie.201806664': ('oxidative_alcohol_amidine', {
        'cv_5d89750f2acfdf79', 'cv_b65f882454199693',
        'cv_4b4d89701a4f8685', 'cv_1bc85ff84e5c8f9e'}),
    '10.1021/acssuschemeng.9b02772': ('aliphatic_vicinal_diol_boronate', {'cv_6af65fcccdeba8bd'}),
}
_REFERENCE_ROLE = {
    'activated_cyanomethyl': 'activated_cyanomethyl', 'aryl_methyl_ketone': 'aryl_methyl_ketone',
    'secondary_amine': 'secondary_amine', 'vicinal_diol_pair': 'catechol',
    'boronic_acid': 'boronic_acid', 'cyclic_anhydride': 'cyclic_anhydride',
}


def linkage_functional_groups(row, mol, binary_roles, self_kind):
    """Give each facet the exact binary roles it supports on this precursor."""
    if row.get('annotation_only') or row.get('catalog_use') == 'synthetic_intermediate':
        return []
    groups = {}
    for candidate in binary_roles:
        role = candidate['role']
        if role == 'aryl_carboxyl':
            ids = {'acyl_chloride' if mol.GetAtomWithIdx(leaving).GetAtomicNum() == 17 else 'carboxylic_acid'
                   for _, _, _, leaving in mol.GetSubstructMatches(_CARBOXYL)}
        elif role in _ROLE_GROUP:
            ids = (_ROLE_GROUP[role],)
        else:
            continue
        for group in ids:
            groups.setdefault(group, []).append({'linkage': candidate['linkage'], 'role': role})

    supplemental = set()
    if self_kind in ('self_imine', 'self_kagome'):
        activation = row.get('self_activation')
        supplemental.update(('masked_aldehyde' if activation == 'dimethyl_acetal' else 'aldehyde',
                             'masked_amine' if activation == 'Boc' else 'primary_amine'))
    elif self_kind == 'self_boroxine':
        supplemental.add('boronic_acid')
    elif self_kind == 'self_triazine':
        supplemental.add('nitrile')
    if row.get('source_dois'):
        families = set(row.get('reported_source_chemistry', []))
        for doi in row['source_dois']:
            source = _SOURCE_FAMILIES.get(doi)
            if source and row.get('id') in source[1]:
                families.add(source[0])
        for family in families:
            supplemental.update(name for name, query in _REPORTED_QUERIES.get(family, {}).items()
                                if mol.HasSubstructMatch(query))
        if row.get('catalog_use') == 'reference_only_family':
            supplemental.update(_REFERENCE_ROLE[g['kind']] for g in row.get('reactive_groups', [])
                                if g.get('count', 0) and g['kind'] in _REFERENCE_ROLE)
        # Reviewed high-connectivity entries can have an explicit single role
        # even when a matching drawing template is unavailable.
        if row.get('review_status') and row.get('role') in _ROLE_GROUP:
            supplemental.add(_ROLE_GROUP[row['role']])
    for group in supplemental:
        groups.setdefault(group, [])
    facets = [dict(id=group, label=_LABELS.get(group, group.replace('_', ' ').capitalize()),
                   binary_roles=roles) for group, roles in groups.items()]
    return sorted(facets, key=lambda f: f['label'])
