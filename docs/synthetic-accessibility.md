# Finite-precursor synthetic accessibility

COFdraw calculates the offline **RDKit Contrib SA Score** for complete finite monomers: **1 means easier and 10 means harder on this heuristic scale**. It combines fragment-frequency contributions with penalties for molecular complexity and the adjustments in RDKit's published implementation. This is the existing SA Score, not a model trained on or distilled from ASKCOS. See the [original Ertl–Schuffenhauer paper](https://doi.org/10.1186/1758-2946-1-8), [RDKit implementation](https://github.com/rdkit/rdkit/blob/master/Contrib/SA_Score/sascorer.py) and [RDKit usage example](https://greglandrum.github.io/rdkit-blog/posts/2023-12-01-using_sascore_and_npscore.html).

The score is useful as a rough molecular-complexity screening aid. It is not a synthesis probability, route, reaction validation, purchasing recommendation, price, safety assessment, yield, or prediction of COF polymerization or crystallization. The original method was developed for small-molecule synthesis assessment; COFdraw has not validated its accuracy for COF precursors. A low score does not establish that a particular COF can be made.

## Inputs and results

`synthetic_accessibility(smiles)` in [cofdraw_synthesis.py](../src/cof_vision/cofdraw_synthesis.py) accepts one connected, valence-valid finite molecular graph with at most 250 atoms and 4000 SMILES characters. This wrapper supports H, B, C, N, O, F, Si, P, S, Se, Cl, Br and I with ordinary covalent bonds. Wildcards/continuation atoms, disconnected salts or mixtures, radicals, unsupported elements and invalid graphs return `available: false`, `score: null` and a specific reason. It does not silently remove counterions, infer missing atoms, neutralize a molecule or replace the graph.

An available result includes the full-precision score, canonical isomeric SMILES, RDKit version, method ID/version, scale direction, fragment-table hash and limitations. Isotopes and formal charges remain in the reported graph. Isotope-label synthesis or procurement cost is not modeled. Source-annotation uncertainty remains applicable: computing a score does not verify the annotated identity or imply an implemented COF linkage.

`catalog_scores(rows)` evaluates the actual supplied monomer records, including reference-only precursors. Each row retains its ID, label, SMILES and available catalog support/provenance qualifiers. The catalog currently contains 81 eligible finite monomers. The result states how many were evaluated successfully; callers must use that count rather than assume every future entry is supported.

`precursor_scores(recipe)` reads only `components[*].monomer_smiles`, for both assembled designs and explicit source graphs. It never scores the product fragment, pore patch or infinite framework. It returns each distinct declared finite precursor graph and an **unweighted mean and maximum** of those scores. Repeated component aliases are grouped; the number of repeated pore instances does not change this summary. If any declared precursor is unavailable, both aggregate values are `null`; individual available results remain visible. No declared precursors also produces no aggregate.

## Offline resources and packaging

The scorer and [fragment table](https://github.com/rdkit/rdkit/blob/master/Contrib/SA_Score/fpscores.pkl.gz) are supplied by the installed RDKit package. The desktop build already collects RDKit's contribution modules and data, and now checks that these resources are present before packaging and in the staged application. No network request, model download, API key or ASKCOS service is used at runtime.

The resource is a trusted bundled compressed pickle, not a user-import format. COFdraw verifies its size and SHA-256 before allowing RDKit to load it. Unknown or missing tables make the score unavailable and fail the desktop build. A table update therefore requires an explicit reviewed hash change and regression checks.

| Resource | Expected value |
|---|---|
| File | `rdkit/Contrib/SA_Score/fpscores.pkl.gz` |
| Bytes | 3,848,394 |
| SHA-256 | `10dcef9340c873e7b987924461b0af5365eb8dd96be607203debe8ddf80c1e73` |

The [Novartis BSD notice](third-party-licenses/RDKit-SA-Score-BSD.txt) is included in source documentation and the Windows distribution's `third-party-licenses` directory. RDKit's original scoring implementation is used unchanged. Regression checks cover the published RDKit example, all catalog monomers, unsupported inputs, resource verification, precursor aggregation and operation without network access. Packaged release checks must also exercise an actual frozen score before distribution.
