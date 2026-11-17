"""Regenerate the bundled reticular model examples with pyCOFBuilder 0.0.9.4.

Development-only: pip install pycofbuilder==0.0.9.4 pymatgen numpy scipy
simplejson ase gemmi. These packages are not needed to use the bundled models.
"""
from pathlib import Path
import hashlib
import importlib.metadata
import json
import numpy as np
import pycofbuilder as pcb
from pymatgen.core import Structure
from pymatgen.io.cif import CifWriter
from pycofbuilder.building_block import BuildingBlock

# 0.0.9.4 appends tetrahedral connector provenance labels as nested lists.
# Normalize this metadata before substitution; atom types and coordinates are
# untouched. Without this, NumPy rejects the ragged labels array for D4 cores.
_add_r_group=BuildingBlock.add_R_group
def _add_r_group_flat_labels(self,*args,**kwargs):
    self.atom_labels=[label for value in self.atom_labels for label in (value if isinstance(value,list) else [value])]
    if len(self.atom_labels)!=len(self.atom_types):
        raise ValueError('Connector provenance labels do not match atom count')
    return _add_r_group(self,*args,**kwargs)
BuildingBlock.add_R_group=_add_r_group_flat_labels

ROOT=Path(__file__).resolve().parents[1]
DEST=ROOT/'data/crystal-models'
MODELS=[
    ('hcb-imine-aa','Hexagonal imine · AA layers','T3_BENZ_CHO_H-L2_BENZ_NH2_H-HCB_A-AA'),
    ('hcb-imine-ab','Hexagonal imine · offset layers','T3_BENZ_CHO_H-L2_BENZ_NH2_H-HCB_A-AB1'),
    ('hcb-expanded-aa','Extended trigonal imine · AA layers','T3_TPBZ_NH2_H-L2_BENZ_CHO_H-HCB_A-AA'),
    ('dia-imine','Diamond-like imine · 3D network','D4_TTPM_NH2_H-L2_BENZ_CHO_H-DIA_A-1'),
]

def main():
    DEST.mkdir(parents=True,exist_ok=True)
    rows=[]
    for identifier,label,name in MODELS:
        framework=pcb.Framework(name,save_bb=False,log_level='error',out_path=str(ROOT/'build/model-generator'))
        structure=Structure(framework.cellMatrix,framework.atom_types,framework.atom_pos,coords_are_cartesian=True)
        distances=structure.distance_matrix.copy();np.fill_diagonal(distances,np.inf)
        shortest=float(distances.min())
        if shortest<.75 or any(str(site.specie) in {'He','X','Q','R'} for site in structure):
            raise ValueError(f'Invalid atom/collision in {identifier}: {shortest}')
        text='# COFdraw unrelaxed reticular template; not an experimental structure assignment.\n'+str(CifWriter(structure,symprec=None))
        path=DEST/(identifier+'.cif');path.write_text(text,encoding='utf-8',newline='\n')
        rows.append({'id':identifier,'label':label,'model_name':name,'path':path.name,
                     'sha256':hashlib.sha256(path.read_bytes()).hexdigest(),'atoms':len(structure),
                     'formula':structure.composition.formula,'minimum_atom_distance_angstrom':shortest,
                     'description':f'{name}. Generated, unrelaxed reticular model with explicit cell and hydrogen atoms. Topology/stacking are chosen assumptions; not a match to the active drawing or an experimental crystal.'})
        print(json.dumps({'model':identifier,'atoms':len(structure),'minimum_distance':shortest}))
    manifest={'generator':'pyCOFBuilder '+importlib.metadata.version('pycofbuilder'),
              'generator_metadata_workaround':'Flatten nested D4 connector provenance labels before functional-group substitution; no coordinate or element changes.',
              'source':'https://github.com/lipelopesoliveira/pyCOFBuilder','citation_doi':'10.1021/acs.jcim.3c01918',
              'qualification':'Explicit unrelaxed reticular designs generated for COFdraw; no experimental identity or stability asserted.',
              'templates':rows}
    (DEST/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n',encoding='utf-8',newline='\n')

if __name__=='__main__':main()
