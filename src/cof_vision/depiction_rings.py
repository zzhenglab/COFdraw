"""Per-molecule symmetric ring basis for finite COF drawing patches.

Large pores can have hundreds of thousands of relevant macrocycles under the
new RDKit RDL default. The legacy symmetric SSSR basis retains the local rings
needed for depiction without enumerating every alternative macrocycle. This
does not alter atom/bond chemistry, aromaticity assignments or global settings.
"""
from rdkit import Chem


def symmetric_depiction_rings(mol):
    algorithms = getattr(Chem, "SymmetrizeSSSRAlgorithm", None)
    if algorithms is not None:
        return Chem.GetSymmSSSR(mol, algorithm=algorithms.LEGACY)
    # Older supported RDKit builds use this algorithm by default and do not
    # expose the per-call selector. Never change a process-wide ring setting.
    return Chem.GetSymmSSSR(mol)
