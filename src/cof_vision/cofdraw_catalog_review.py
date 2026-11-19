"""Display names and source roles, kept separate from signed upstream records."""
import re

# Chen et al., 10.1021/acs.chemmater.4c01292, SI pp. S3–S9. Descriptions
# summarize distinguishing groups in the drawn structures; IDs remain searchable.
INTERMEDIATES = {
    'S1':'4-Bromopyridine-2,6-dicarbaldehyde',
    'S2':'4-Bromo-2,6-bis(dimethoxymethyl)pyridine',
    'S3':'Pyridyl bisacetal pinacol boronate',
    'S4':'Bromo-bis(hexyloxy)phenyl pyridyl bisacetal',
    'S5':'5′-Bromo-m-terphenyl-4,4″-diamine',
    'S6':'Diaminoterphenyl pinacol boronate',
    'S7':'Bromo-dimethoxyphenyl pyridyl bisacetal',
    'S8':'Bromo-dimethoxyphenyl terphenyl diamine',
    'S9':'Bis(Boc)-protected bromo-dimethoxyphenyl terphenyl diamine',
}
NAMES = {
    'cv_3aac60550b8330ff':'Tributyltriazatruxene tricarbaldehyde',
    'cv_5141c71c7cba5e76':'p-Terphenyl-4,4″-diamine',
    'cv_05cbf00025a2d4a1':'Bis(hexyloxy) diamine–pyridyl bisacetal',
    'cv_e0efd03afb224b2a':'Dimethoxy diamine–pyridyl bisacetal',
    'cv_6bb0458f6411d793':'Bis(Boc)-protected diamine–dialdehyde',
}
MASKED_SELF = {'cv_05cbf00025a2d4a1':'dimethyl_acetal','cv_e0efd03afb224b2a':'dimethyl_acetal',
               'cv_6bb0458f6411d793':'Boc'}


def reviewed_row(row):
    alias = row.get('abbreviation', '')
    intermediate = '10.1021/acs.chemmater.4c01292' in row.get('source_dois',[]) and alias in INTERMEDIATES
    name = INTERMEDIATES[alias] if intermediate else NAMES.get(row['id'])
    if not name and re.fullmatch(r'(?:compound\s*\d+|r\d+-\d+)',alias,re.I):
        name = row.get('full_name')
    if name:
        previous = row.get('full_name') or alias
        full = name if re.fullmatch(r'(?:monomer|compound)\s+[SM]?\d+',previous,re.I) else previous
        row = {**row,'abbreviation':name,'full_name':full,'label':name,'paper_label':alias,
               'aliases':list(dict.fromkeys([alias,previous,*row.get('aliases',[])]))}
    if intermediate:
        row = {**row,'catalog_use':'synthetic_intermediate','automatic_prepare_supported':False,
               'qualification':'Synthetic intermediate in the preparation of M1–M3 (supporting information, S3–S9).'}
    if row['id'] in MASKED_SELF:
        row = {**row,'self_activation':MASKED_SELF[row['id']], 'catalog_use':'masked_self_precursor',
               'qualification':'Reported latent aminoaldehyde. The self-condensation design includes the source-specific deprotection.'}
    return row
