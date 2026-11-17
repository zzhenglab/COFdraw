# Chemical export validation

COFdraw provides MOL and SDF chemical graphs, CDXML chemical drawings when the installed RDKit build supports its writer, and separate SVG/PNG presentation images. The JSON recipe/audit carries the source, attachment ports and qualified COFex reference. A pore patch is not a finite representation of the entire periodic material.

The user exposed a real import defect: the native writer emitted boundary dummies as `Element="0"`, which ChemDraw rejected as an invalid atomic number. The exporter now converts only neutral, unlabelled, degree-one boundary dummies into `NodeType="ExternalConnectionPoint"` with `ExternalConnectionType="Wavy"`. It retains their positions, neighbor and bond order, removes element-only attributes, and rejects any remaining zero-element node. It does not replace a boundary with carbon or hydrogen.

This representation follows the original ChemDraw SDK: external connection points are a distinct node type, and `Wavy` denotes the perpendicular wave at a bond end. An ordinary element node without `Element` defaults to carbon; simply deleting `Element="0"` would therefore be chemically wrong. [Node types](https://iupac.github.io/IUPAC-FAIRSpec/cdx_sdk/properties/Node_Type.htm), [external-connection display](https://iupac.github.io/IUPAC-FAIRSpec/cdx_sdk/properties/Atom_ExternalConnectionType.htm), [element semantics](https://iupac.github.io/IUPAC-FAIRSpec/cdx_sdk/properties/Node_Element.htm).

Revvity documents desktop read/save support for CDX/CDXML, MOL and SDF, including V2000 and V3000. ChemDraw 23.0 made V3000 its default MOL export format. These establish format support, not a test of our particular exported files in ChemDraw. [Official format matrix](https://revvitysignals.com/sites/default/files/2025-06/Full%20Dot%20Chart%20ChemDraw%20Suite%20of%20Products.pdf), [23.0 release notes](https://revvitysignals.com/sites/default/files/2024-02/rs-release-notes-DT-23.0-chemdraw.pdf).

The current RDKit 2026.09.2 build includes its optional native ChemDraw writer. RDKit explicitly describes its CDXML implementation as covering molecular/reaction basics, rather than the complete drawing format. [RDKit file-format API](https://www.rdkit.org/docs/source/rdkit.Chem.rdmolfiles.html#rdkit.Chem.rdmolfiles.MolToCDXMLBlock).

| Export | What it carries | Limit |
| --- | --- | --- |
| MOL / SDF | Atoms, bonds, charge, isotope, stereo and 2D positions | Continuation dummies can display as R or stars; they are boundary ports. No literature font/color guarantee. |
| CDXML | Chemical objects, layout and explicit wavy external connection points, checked against the MOL tables | Current exporter does not transfer the SVG font, component colors or collapsed labels. |
| SVG / PNG | Selected presentation style | SVG paths are visual drawing objects, not a ChemDraw chemical graph. |
| JSON | Recipe, source, reference qualification and port metadata | Keep with the chemical file; a precursor-derived key does not establish a unique periodic structure. |

The export check loads MOL without full sanitization, then validates the emitted CDXML's XML records directly. It reconstructs only a local graph for strict valence/hydrogen and wedge-stereo checks with FastFindRings. Atom/bond tables, isotopes, charges, tetrahedral tags and planar coordinates are compared; positive scale and translation are allowed, reflection is rejected. Explicitly unknown double bonds must retain their wavy flag, and a defined E/Z bond must not acquire one. Full macrocycle aromaticity assignment, canonical ranking and native CDXML reparsing are excluded from this production path.

This separation fixes two observed stalls: the default MOL reader on the 234-atom TPB–DMTP pore, and the native CDXML reader even with sanitization disabled on a larger automatically generated pore. Small fixtures without boundary points still use the native reader as an independent round-trip check. The tested RDKit CDXML reader does not preserve `ExternalConnectionPoint` semantics and can reinterpret these nodes as carbon, so it is not an authority for port round trips. The production export diagnostic records direct XML validation, not a native-reader round trip.

An unavailable or failed optional CDXML export preserves the existing MOL/SDF/SVG/PNG. The adjacent `*.chemical-export.json` records the reason; a stale prior CDXML is removed before retry. No file is called CDXML merely because it contains an SVG image.

The regression suite checks an isotope-bearing stereocenter with E and Z alkenes, single/double/triple boundary bonds, rejection of charged/isotopic/nonterminal dummies, a zwitterion, a dative copper complex, unavailable native support, deliberate stereo loss, both unknown/defined double-stereo changes, and exclusion of the native CDXML reader, full sanitization and canonical ranking. Existing framework tests also verify stereo-selector remapping and reject geometries that would invert explicitly specified double-bond stereochemistry. Run:

```powershell
..\.venv\Scripts\python.exe -m pytest tests/test_cofdraw_exports.py tests/test_framework_drawing_stereo.py -q
```

The earlier local registry round-trip audit is stored outside the repository at `cof-vision-data/qa/cofdraw-export-audit-0001/`. It covers TFB–TAB (108 atoms), TpPa (126), TPB–DMTP (234), and TpPa–SO3H (156), each with six continuation ports. Audit `0002` adds the actual automatic DMTP–TAPB API preset (240 atoms, six ports), and `0003` rechecks all five MOL files with the final direct XML validator and raw-byte provenance hashes. The `0002` diagnostic's `mol_sha256` used normalized text rather than raw file bytes; `0003` corrects that reporting convention without changing the input graphs. The four source standard PNGs were visually inspected; this does not validate CDXML appearance inside ChemDraw.

## Actual ChemDraw import and save checks

The fixed exports were opened by the installed **ChemDraw 23.1.2.7** through its `ChemDraw_x64.Application` automation interface, then saved by ChemDraw to CDXML, MOL and PNG. Nine isolated, initially empty, hidden instances were used, closing only their temporary test documents; an existing user process was preserved. All nine native PNGs were actually inspected.

| Cases | Atoms including ports | Ports | Result |
| --- | --- | --- | --- |
| TFB–TAB / TPB–DMTP / TpPa / TpPa–SO3H | 108 / 234 / 126 / 156 | 6 each | Open/save passed; CDXML atom, bond, hydrogen and boundary tables retained |
| Automatic DMTP–TAPB preset | 240 | 6 | Same checks passed |
| Isotope + tetrahedral center + E/Z pair | 6 each | 0 | Saved MOL isomeric graph matched each input exactly |
| Single-boundary fixture / double-boundary fixture | 10 / 5 | 2 each | Correct external points and boundary bond orders retained in CDXML |

All saved CDXMLs retained the chemical-node IDs, element/charge/isotope records, bond orders/displays, hydrogen totals and connection types. Coordinates differed only by translation, within the serialized precision. These checks use explicit tables and local valence, avoiding whole-pore aromaticity or canonicalization. The hash-bound private report is `cof-vision-data/qa/chemdraw-native-port-fix-0001/native-roundtrip-audit.json` (SHA-256 `f013bed99ce750c537ab98626b29eb74ad3b4e6effbbefcac4bca6ab398d9027`). Its portable copy is [the native compatibility report](../benchmarks/cofdraw-chemdraw-native.json).

**ChemDraw's own SaveAs MOL removed every external connection point in these tests.** Keep edited boundary structures as CDXML; retain COFdraw's original MOL/SDF and JSON as the explicit boundary-graph record. Do not use that ChemDraw MOL save-back as a port-preserving training target or silently interpret its new terminal hydrogens as part of the framework.

These results cover the tested files and version, not arbitrary ChemDraw content, automatic cleanup or all presentation styles. The rebuilt frozen executable still needs its packaging smoke. Per-export `ChemDraw_desktop_import_tested:false` remains correct until that particular generated file is tested; the separate QA report identifies the files actually opened. Future releases should repeat the version-specific import, chemical comparison and stereo checks. Fonts, colors and line styles need separate validation. [Official settings documentation](https://support.revvitysignals.com/hc/en-us/articles/4408239581204-ChemDraw-Document-and-Object-Settings-of-a-ChemDraw-Document).
