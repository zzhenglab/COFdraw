# RDKit large-pore ring calculation: independent probe

The installed RDKit **2026.09.2** expands the 240-atom, 270-bond DMTP test patch to **262,174 relevant cycles** under its default ring algorithm. Explicit per-molecule `LEGACY` symmetrization returns 84 rings and keeps all 30 six-member ring bond sets. This is a ring-enumeration bottleneck; the probe does not justify disabling stereochemistry or setting `_StereochemDone` manually.

Upstream `FindRings.cpp` selects RDL for the current default and materializes its relevant cycles; the explicit legacy branch retains an SSSR plus selected extra rings. These are different ring representations, not different atom/bond graphs. The public Python API exposes `SymmetrizeSSSRAlgorithm.LEGACY`. [Ring implementation](https://github.com/rdkit/rdkit/blob/master/Code/GraphMol/FindRings.cpp), [Python API](https://www.rdkit.org/docs/source/rdkit.Chem.rdmolops.html#rdkit.Chem.rdmolops.GetSymmSSSR).

`MolDraw2D::drawMolecule()` requires ring information marked as symmetrized and recalculates it when that condition fails, independently of `prepareMolsBeforeDrawing`. An ordinary `GetSSSR` or `FastFindRings` cache does not satisfy this guard. `DrawMolecules([mol])` bypassed that particular guard in the inspected source and bounded probe, but explicit legacy symmetrization is the less intrusive choice. [Drawing implementation](https://github.com/rdkit/rdkit/blob/master/Code/GraphMol/MolDraw2D/MolDraw2D.cpp), [RingInfo types](https://github.com/rdkit/rdkit/blob/master/Code/GraphMol/RingInfo.h).

`SanitizeMol` always clears computed properties first. Therefore, computing the legacy rings before the remaining individual sanitizer calls loses that cache. Some operations may create a replacement ring basis; this does not establish that the final cache is symmetrized. On the explicit-Kekule DMTP test, every individual sanitizer tested, including `SANITIZE_NONE`, left RingInfo uninitialized. Direct molecule copying, `Kekulize`, `WedgeMolBonds`, `SetDoubleBondNeighborDirections`, and `UpdatePropertyCache` retained the final legacy cache. [Sanitization implementation](https://github.com/rdkit/rdkit/blob/master/Code/GraphMol/MolOps.cpp), [sanitization documentation](https://www.rdkit.org/docs/RDKit_Book.html#molecular-sanitization).

Recommended placement, without changing shared process settings:

```python
# Complete required sanitizer operations and source-local aromatic restoration.
if hasattr(Chem, "SymmetrizeSSSRAlgorithm"):
    Chem.GetSymmSSSR(mol, algorithm=Chem.SymmetrizeSSSRAlgorithm.LEGACY)
else:
    Chem.GetSymmSSSR(mol)  # Compatibility with older supported RDKit versions.
# Then preserve/assign the already specified stereo directions, copy and draw.
# A later SanitizeMol requires repeating the final ring-cache preparation.
```

One local run gave the following diagnostic measurements; these are not a comprehensive performance benchmark.

| Operation | Time | Rings |
|---|---:|---:|
| Default `MolFromMolFile(sanitize=False)` | 3.93 s | 262,174 |
| Draw already parsed default molecule | 0.950 s | 262,174 |
| Replace large default cache with explicit legacy cache | 0.038 s | 84 |
| Draw with legacy cache | 0.0031 s | 84 |
| Final legacy calculation after nonaromatic sanitizers | 0.0016 s | 84 |

The default and legacy 700×700 PNGs were pixel-identical. Atom elements, isotopes, charges, explicit H fields, chiral tags, bond orders/directions/stereo selectors, and stereo-group membership were unchanged. Six bounded controls also produced identical pixels: isotope/tetrahedral/alkene, enhanced AND stereogroup, fused aromatic, bridged cage, ring tetrahedral, and unknown double-bond cases. This does not establish universal equivalence for atropisomers, all enhanced-group arrangements, unusual query graphs, or changed global aromaticity. Preserve their existing validation paths.

Parsing remains separate: `sanitize=False` still calls bond-stereo detection, which calls `setDoubleBondNeighborDirections`, which demands symmetrized rings. No per-call parser ring-algorithm option was found in the installed Python API. Avoid reparsing an internally available graph; do not remove stereo processing from arbitrary imported MOL files. [MOL parser](https://github.com/rdkit/rdkit/blob/master/Code/GraphMol/FileParsers/MolFileParser.cpp), [bond-stereo implementation](https://github.com/rdkit/rdkit/blob/master/Code/GraphMol/Chirality.cpp).

Probe artifacts are under `cof-vision-data/qa/cofdraw-performance/`: `independent_ring_draw_probe.py`, `independent_ring_cache_probe.py`, and `independent-ring-probe/{report.json,cache-lifetime-report.json}`. Input MOL SHA256: `bb36d9b37cb5f2d2ea3be9f41e1cb6653924c0ff67f3234555bb8779ff24d89e`. The inspected upstream links reference mutable master; installed-version measurements are the direct evidence. No app code, global ring/stereo setting, source annotation, or master data was modified by this investigation.
