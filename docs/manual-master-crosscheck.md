# Manual master-key and ground-truth cross-check

This update adds **18 individually reviewed finite monomers**, bringing the library from 63 to **81**. Fifteen additions can be used with the existing scoped drawing templates; three are reference-only. All 81 finite entries have an available offline RDKit SA Score. Neither a score nor a master-key match establishes synthetic feasibility.

The review used the sealed annotation catalog **0042** and the installed reference index identified in the current continuation pointer. The index hash, current master/mapped-file hashes, source DOI, exact isomeric graph, source-crop hash, annotation hash and applicable installation-receipt hashes are retained in the portable evidence. Private paths, source publication images and PDF text are excluded.

## What was admitted

| Entries | Source DOI | Drawing scope |
| --- | --- | --- |
| TFPC; 1,5-diaminonaphthalene | 10.1039/c8ta05070h | Imine design candidates |
| A1 and A2 tetraaldehydes | 10.1021/jacs.9b08017 | Imine design candidates; source aminal reaction is not implemented |
| TABPB; BPP | 10.1021/jacs.4c12973 | Imine design candidates |
| DBA[12]-CHO | 10.1021/jacs.0c01026 | Imine design candidate with corrected annulene graph |
| BFBAPy | 10.1021/jacs.9b03463 | One explicitly selected reaction role; both aldehyde and amine groups remain in the precursor |
| TBT | 10.1021/acsmaterialslett.2c00669 | Imine design candidate; source imide reaction is not implemented |
| DHzDPr | 10.1039/c7cc05779b | Acylhydrazone design candidate |
| BPF | 10.1038/s44160-022-00071-y | Imine design candidate; no complete source material key is asserted |
| BTBBA; BTBBAv | 10.1039/d4tc02995j | Diamine design candidates; acetals stay protected |
| Hth (S,S); Hth (R,R) | 10.1016/j.cclet.2020.11.063 | Acylhydrazone candidates retaining both source-assigned stereocenters |
| APH; ASH | 10.1021/jacs.7b00925 | Reference-only: mixed amine/hydrazide chemistry is outside the current templates |
| FLT diboronic acid | 10.1021/jacs.7b05555 | Reference-only: boronic-acid framework chemistry is not implemented |

Names and abbreviations are tied to exact graphs and source papers. The new entries have no normalized name/alias collisions or duplicate graphs with the previous 63 entries. In particular, new **TFPC** has phenoxy bridges; existing **TFPT** has direct phenyl–triazine bonds. A1/A2 names are explicitly scoped to the Aminal-COF paper. The 63 existing rows and all 17 source recipe bindings remain unchanged.

## What “ground truth” means here

The inspected sealed files contain 2,489 annotation rows, all explicitly marked `ground_truth: false`. Two available human-review templates contain 20 and 100 unreviewed decisions, respectively. A separate 25-case coarse audit identifies its reviewer as an agent. These artifacts do **not** constitute completed independent human adjudication.

The newer source reviews and installed corrections strengthen the evidence while remaining **provisional agent reviews**. Every new row retains `reviewer_type: "agent"` and `ground_truth: false`. Its `master_binding` confirms only a source DOI and exact finite precursor graph in the installed reference. It does not certify a whole COF, topology, catenation, crystal structure or experimental result.

The historical PSM v9 and installed-master v8b reports were checked, but they predate the current sealed catalog. Their older inventory totals and master hashes are not substituted for the current snapshot. APH and ASH are the finite parent precursors; grafted PSM states are not imported as monomers.

## Selection and exclusions

The complete inventory scan found 689 distinct finite graphs absent from the previous drawing library. This includes reagents, catalysts, model compounds, intermediates and products. Twenty-eight finite molecule graphs had an exact match at the same DOI in the installed master. Nineteen candidates received individual source-crop review; 18 were admitted.

TFPyriPy-Ac was excluded because its current source occurrence carries an identity hold. Ten of the 28 master-matched candidates were not admitted in this update: a master match alone is insufficient to establish a reviewed monomer role or drawing readiness. This bounded selection does not claim that every remaining molecule is unusable or incorrect.

The R,R-Hth occurrence has a neighboring reaction-arrow hold for unmasked image training. Its complete molecular graph and stereochemistry were independently bound to the source MOL and installed correction. That image-training hold is preserved in its evidence; no training permission or eligibility is granted by this library update.

## Reproducible records and checks

- [Explicit admission selection](../data/framework-drawing/catalog/master-key-monomer-selection.json)
- [Admission audit and installed snapshot hashes](../data/framework-drawing/catalog/master-key-monomer-admission.audit.json)
- [Ground-truth review, exclusions and existing-entry reconciliation](../data/framework-drawing/catalog/master-key-ground-truth-review.audit.json)
- [Read-only research-input importer](../scripts/import_master_key_monomers.py)
- [Graph, binding, stereo and actual-pair regression tests](../tests/test_master_key_monomer_admission.py)

The importer checks the current pointer, sealed annotation files, source PDFs/crops, graph files where available, installed master/mapped hashes and pinned correction receipts. It rejects held identities, source conflicts, boundary atoms, disconnected graphs, radicals, graph duplicates and changed source bindings. It snapshots the previous catalog and preserves the earlier [14-entry admission](manual-monomer-update.md) as a separate historical event.

Validation builds an actual connected patch for each of the 15 supported additions, preserves both Hth enantiomers, tests DOI/stereochemical binding failures, verifies unchanged earlier rows and source bindings, and checks every bundled file against its manifest. These checks establish the stated software behavior and evidence consistency; they do not convert agent labels into human gold.
