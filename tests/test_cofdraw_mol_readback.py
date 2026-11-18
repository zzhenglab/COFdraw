"""Bounded import checks for para-phenyl macrocycles, isolated from app flags."""
import hashlib
import json
import subprocess
import sys

import numpy as np
import pytest
from rdkit import Chem

from cof_vision.cofdraw_chemistry import prepare_template
from cof_vision.framework_drawing import VERSION, build_patch


READER = r'''
import json,sys
from pathlib import Path
from rdkit import Chem
if hasattr(Chem, 'SetUseLegacyRingFinding'):
    Chem.SetUseLegacyRingFinding(True)
mol=Chem.MolFromMolBlock(Path(sys.argv[1]).read_text(encoding='utf8'),
                         sanitize=False,removeHs=False)
assert mol is not None
mol.UpdatePropertyCache(strict=True)
print(json.dumps({
    'atoms':[[a.GetAtomicNum(),a.GetIsotope(),a.GetFormalCharge(),a.GetTotalNumHs(),int(a.GetChiralTag())] for a in mol.GetAtoms()],
    'bonds':[[b.GetBeginAtomIdx(),b.GetEndAtomIdx(),b.GetBondTypeAsDouble()] for b in mol.GetBonds()],
    'xy':mol.GetConformer().GetPositions()[:,:2].tolist(),
    'defined_double_stereo':sum(b.GetStereo() not in (Chem.BondStereo.STEREONONE,Chem.BondStereo.STEREOANY) for b in mol.GetBonds()),
    'unknown_double_stereo':sum(b.GetStereo()==Chem.BondStereo.STEREOANY for b in mol.GetBonds())
}))
'''


@pytest.mark.parametrize('layout', ['manual', 'bond_angle_pose_v2'])
@pytest.mark.parametrize('v3000', [False, True])
def test_symmetric_pore_mol_can_be_imported_without_changing_graph_or_coordinates(tmp_path, layout, v3000):
    # The fully symmetric 252-atom pore triggers RDKit2026.09's documented
    # expensive symmetric-ring enumeration even with sanitize=False. Use its
    # documented legacy import switch in a separate process, never globally.
    recipe = {'schema_version': VERSION, 'recipe_id': 'readback-pore',
              **prepare_template(monomer_a='O=Cc1ccc(-c2ccc(C=O)cc2)cc1',
                 monomer_b='Nc1ccc(-c2cc(-c3ccc(N)cc3)cc(-c3ccc(N)cc3)c2)cc1',
                 linkage='imine', label_a='BPDA', label_b='TAPB', instances=12,
                 template_id='hcb_3_2'), 'layout_refinement': layout}
    mol, _ = build_patch(recipe)
    assert mol.GetNumAtoms() == 252
    path = tmp_path/'pore.mol'
    Chem.MolToMolFile(mol, str(path), forceV3000=v3000)
    before = hashlib.sha256(path.read_bytes()).hexdigest()
    completed = subprocess.run([sys.executable, '-c', READER, str(path)],
                               capture_output=True, text=True, timeout=10,
                               creationflags=getattr(subprocess, 'CREATE_NO_WINDOW', 0))
    assert completed.returncode == 0, completed.stderr
    restored = json.loads(completed.stdout)
    expected = Chem.Mol(mol); Chem.Kekulize(expected, clearAromaticFlags=True)
    assert restored['atoms'] == [[a.GetAtomicNum(), a.GetIsotope(), a.GetFormalCharge(),
                                 a.GetTotalNumHs(), int(a.GetChiralTag())] for a in expected.GetAtoms()]
    assert restored['bonds'] == [[b.GetBeginAtomIdx(), b.GetEndAtomIdx(), b.GetBondTypeAsDouble()]
                                for b in expected.GetBonds()]
    assert np.allclose(restored['xy'], mol.GetConformer().GetPositions()[:, :2], rtol=0, atol=.000051)
    assert restored['defined_double_stereo'] == 0
    assert restored['unknown_double_stereo'] == 12  # unspecified imines remain unspecified
    assert hashlib.sha256(path.read_bytes()).hexdigest() == before
