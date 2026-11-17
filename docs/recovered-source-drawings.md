# Recovered source drawings

The catalog includes 17 source drawings and 49 monomers, including 13 reviewed chemical diagrams that were absent from the original four-drawing desktop bundle. Their chemical graphs, continuation interfaces, coordinates, source qualifications, and presentation metadata are now portable. Original publication images, PDFs, page text, and workstation paths are excluded.

The existing TFB–TAB and sulfonated TpPa/SCOF entries also regain their previously reviewed manual source layouts. TFB–TAB retains all 108 anchored atom positions. SCOF retains 138 visible positions; the 18 oxygen atoms hidden inside source SO3H abbreviations receive explicitly qualified illustrative coordinates for chemical export. They are not reported as observed source positions.

| Drawing | Atoms, including boundaries | Boundary ports | Scope |
| --- | ---: | ---: | --- |
| F-COF-1 | 90 | 6 | Atomwise source-registered vinylene pore; six explicit E bonds and decorative fluorescence halo |
| COF-300 | 402 | 16 | Source-informed projected network; crossing lines are not extra bonds |
| BTA-COF1 | 59 | 8 | Local boronate fragment, not a full pore |
| BTA-COF2 | 53 | 6 | Local boronate fragment, not a full pore |
| BTA-COF3 | 55 | 6 | Local boronate fragment, not a full pore |
| BTA-COF4 | 120 | 9 | Boroxine pore with explicit single B–O bonds |
| Brick-wall COF | 240 | 6 | Source-informed T-shaped-node pore |
| Herringbone proposal | 240 | 6 | Rejected source hypothesis, clearly labeled as such |
| COF-BP | 630 | 12 | Source-informed parent network with component highlights |
| COF-BP-Ni1 | 666 | 12 | Four drawn coordination sites, not a bulk occupancy claim |
| COF-BP-Ni2 | 694 | 12 | Four drawn coordination sites, not a bulk occupancy claim |
| Pyurea-COF | 232 | 8 | Compact source-informed pyrene pore; manually placed blocks, not exact atom registration |
| UCOF | 144 | 6 | Source-informed urethane pore; source angular arrangement and red urethane links |

F-COF-1 is the name attached to the manually registered red-and-black, yellow-halo drawing in its source review, DOI [10.1021/acsami.2c11052](https://doi.org/10.1021/acsami.2c11052). The halo is a decorative illustration, not calculated electron density. Most other recovered diagrams use organized, source-informed geometry; they are not claimed to be pixel registrations. Their saved coordinates are retained without an automatic layout replacement.

The COF-300 diagram is a 2D projection of a source network illustration, not a calculated 3D structure or a crystallographic unit cell. The independent review verifies its 10 nodes, 12 linkers, 24 internal imines, and 16 boundary cuts. Its faded branch remains distinct from the foreground and overlapping projected phenyls remain unconnected.

The nickel drawings retain explicit charge-balanced coordination models: four Ni atoms and 16 dative bonds in each example. This normalization records the source's drawn contacts; it does not establish a unique crystal, uniform loading, bulk coordination number, or transfer the parent COFex identity to a metal-containing material.

Pyurea-COF and UCOF are new, manually source-informed coordinate revisions. Their original archived graphs used generic layouts. The restored versions use source block centers and orientations, preserve every aromatic block without shear or stretching, and correct junction spacing by rigid translations. UCOF also changes terminal carbonyl display orientation where needed to remove overlapping oxygen positions. All atom identities, hydrogen counts, charges, bonds and stereo fields match the reviewed graphs. Their `layout_evidence` records the source-image hash, anchors, adjustments and limitations; exact pixel registration is not claimed. The original source photograph inside the UCOF pore is excluded.

Eight additional source-derived monomers accompany these examples: 2,4,6-trimethyl-1,3,5-triazine, the triazine triboronic acid used by BTA-COFs, the T-shaped trisaldehyde used by brick-wall COF, tetraformylbiphenyldiol, benzothiadiazole dianiline, 1,4-phenylenediurea, 2,4,6-trihydroxybenzaldehyde, and 1,4-phenylene diisocyanate. Pyurea's pyrene precursor was already present. An entry in the monomer library does not imply that an automatic linkage transform is implemented.

## Provenance and reproducibility

Each recovered recipe contains `explicit_graph` atom/bond tables, a 2D conformer, source identifiers and hashes, and `source_style` display metadata. Graph aromaticity, formal charges, hydrogen counts, isotope/stereo fields, and boundary bonds are retained explicitly. Complete-material SMILES and unverified topology assignments are omitted.

The original graph revisions and bounded review artifacts are identified by SHA-256 hashes in each recipe. The registry, current catalog metadata, recipe bindings, and portable bundle manifest bind the exact distributed recipe bytes. See [`recovered-drawings.audit.json`](../data/framework-drawing/catalog/recovered-drawings.audit.json) for the import inventory and bounds of the review. Agent review remains distinct from independent human verification.

Maintainers with the original read-only evidence archive can reproduce the curated import:

```powershell
python scripts/import_reviewed_catalog.py --reviews path/to/manual-wave3/full-paper-reviews --destination data/framework-drawing
```

The importer selects reviewed revisions, verifies review/draft hashes, validates the explicit graphs and defined double-bond geometry, omits private paths and publication text, and refreshes all current registry/catalog/manifest bindings. It does not modify the evidence archive.

To admit the Pyurea/UCOF manual revisions, also provide `--layout-revisions path/to/revisions`, containing `pyurea-source-layout-adjusted.recipe.json` and `ucof-source-layout-adjusted.recipe.json`. These are the intermediate placement records with `layout_evidence`; the importer verifies unchanged chemistry against the peer-bound upstream graphs before copying their coordinates into the public recipes. The two distributed recipes in `data/framework-drawing/recovered/` are the final portable coordinate revisions.

For the other recovered explicit diagrams, source coordinate conversion applies only a shared translation, a uniform scale, and a y-axis convention change. Rendering any saved source recipe retains its admitted coordinates. Future layout corrections must be admitted as a new reviewed coordinate revision with their scope recorded.
