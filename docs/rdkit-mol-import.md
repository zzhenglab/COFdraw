# Importing large COF pore files with RDKit

COFdraw's MOL and SDF exports retain the finite chemical graph, continuation atoms, coordinates and specified or unspecified stereochemistry. A slow RDKit import does not necessarily mean the file is malformed.

RDKit 2026.09 changed its symmetric-ring enumeration. Its release notes identify large macrocycles containing several para-substituted rings as a case where the new algorithm can be expensive. The documented legacy switch is faster for these graphs, but its symmetric ring set can be incomplete. This is a reader limitation; changing the output to V3000 does not resolve it. [RDKit 2026.09 release notes](https://github.com/rdkit/rdkit/releases/tag/Release_2026_09_1).

The 252-atom BPDA/TAPB test pore reproduced this behavior with RDKit 2026.09.2: default imports exceeded a 10-second child-process limit for both V2000 and V3000, with sanitization enabled or disabled. Raw and fitted coordinates behaved the same. Changing the stereo-perception setting did not help. The documented ring switch reduced the measured parsing time to about 15–16 ms on the development machine. These are diagnostic timings, not a performance guarantee.

For a separate RDKit workflow that needs the literal exported chemical table, run this in a short-lived importer process:

```python
from pathlib import Path
from rdkit import Chem

# This affects the current Python process only. Keep it out of a shared
# scientific service: the legacy symmetric ring set may be incomplete.
if hasattr(Chem, "SetUseLegacyRingFinding"):
    Chem.SetUseLegacyRingFinding(True)

mol = Chem.MolFromMolBlock(
    Path("pore.mol").read_text(encoding="utf-8"),
    sanitize=False,
    removeHs=False,
)
if mol is None:
    raise ValueError("MOL import failed")
mol.UpdatePropertyCache(strict=True)
```

Retain `sanitize=False` when the purpose is to preserve the exported table. Whole-pore sanitization can reassign aromaticity across a macrocycle; it is a separate chemical interpretation and is not equivalent to preserving the source-local aromatic mask. The unsanitized MOL representation uses its explicit Kekulé bond orders. Keep the project JSON for the original graph provenance and aromaticity scope.

This workaround is for reading and editing the literal graph. It does not validate a complete symmetric ring set, COF topology, crystal structure or whole-framework aromaticity. COFdraw does not change the application's global ring-finding policy. Its renderer passes the already checked graph directly to the chemical exporter, bound to the written MOL, so ordinary export avoids this parser cost.

The application’s monomer importer expects a molecular building block. Do not import an entire large pore there: its direct RDKit normalization may encounter the same ring-enumeration slowdown. Open the saved project JSON to continue editing a COFdraw pore, or use the separate reader above for MOL inspection.

Regression checks cover both raw and fitted BPDA/TAPB layouts in V2000 and V3000, with each reader running in a child process under a 10-second limit. They verify atom/isotope/charge/hydrogen records, literal bond orders, rounded coordinates, unchanged file hashes and all 12 unspecified imine bonds. Actual frozen-executable Tp/TAPB and BPDA/TAPB MOL files also passed these checks. This RDKit check is separate from the [native ChemDraw validation](chemdraw-export-validation.md).
