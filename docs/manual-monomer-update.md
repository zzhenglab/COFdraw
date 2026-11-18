# Reviewed manual monomer update

This first admission expanded the drawing catalog from 49 to **63 distinct finite monomer graphs**, adding 14 entries from newer local manual annotations. Ten additions have recognized groups that work with existing drawing transforms; four are reference-only. The subsequent [manual master-record cross-check](manual-master-crosscheck.md) adds another 18 entries, bringing the current catalog to **81**. All entries remain agent-reviewed and provisional. These updates do not claim human verification or authorize training.

| Added IDs | Source | Drawing scope |
| --- | --- | --- |
| `tfppy_ome`, `tfppy_f`, `tfppy_13c` | [10.1002/anie.202208833](https://doi.org/10.1002/anie.202208833) | Four aryl aldehydes; methoxy, fluorine and four formyl-carbon-13 labels retained |
| `dabd` | [10.1021/acsmacrolett.6b00486](https://doi.org/10.1021/acsmacrolett.6b00486) | Two arylamines; neutral drawing, with experimental dihydrochloride qualification |
| `thz` | [10.1039/d1qm00231g](https://doi.org/10.1039/d1qm00231g) | Terephthalic dihydrazide; two terminal acylhydrazide groups |
| `tape` | [10.1021/acs.chemmater.1c04051](https://doi.org/10.1021/acs.chemmater.1c04051) | Four arylamines; central C=C retained |
| `odh` | [10.1021/acsami.8b18502](https://doi.org/10.1021/acsami.8b18502) | Oxalyldihydrazide; two terminal acylhydrazide groups |
| `tfb_oh`, `tfb_oh2` | [10.1021/acsami.5c01437](https://doi.org/10.1021/acsami.5c01437) | Three aldehydes with one or two phenolic hydroxyl groups |
| `daaq` | [10.1039/d4sc00241e](https://doi.org/10.1039/d4sc00241e) | 2,6-Diaminoanthraquinone; both quinone carbonyls retained |
| `ppcn` | [10.1002/anie.202203250](https://doi.org/10.1002/anie.202203250) | Reference-only nitrile precursor |
| `tbpm_boronic`, `allyl_tbpa` | [10.1039/c3cc40358k](https://doi.org/10.1039/c3cc40358k) | Reference-only boronic-acid precursors with tetrahedral carbon centers |
| `tmc` | [10.1002/aoc.5172](https://doi.org/10.1002/aoc.5172) | Reference-only trimesoyl chloride |

Functional-group compatibility supports an explicitly chosen new design. It does not reproduce the source synthesis automatically. In particular, source thiazole, benzobisoxazole, polyimide, boronic-acid and acyl-chloride reactions remain outside the implemented transforms. Source compound numbers are preserved in the evidence rather than promoted to globally unique chemical names. Site counts do not establish symmetry, a net, or a periodic structure.

Each selected source crop or paired source/graph image was inspected for connectivity, substituent placement, completeness and relevant isotope labels. Annotation files were checked against their sealed upstream curation catalog. Source PDF/crop and review-artifact hashes were verified, and available MOL graphs matched the admitted isomeric SMILES. Model compounds, catalytic products, postsynthetic reagents, ambiguous labels and boundary fragments were excluded from this selected update.

The selection is pinned in [manual-monomer-selection.json](../data/framework-drawing/catalog/manual-monomer-selection.json). The [admission audit](../data/framework-drawing/catalog/manual-monomer-admission.audit.json) records scope, counts and hashes. This monomer admission preserved the existing 49 entries and all 17 source-drawing bindings. The subsequent [source-layout revision](recovered-source-drawings.md) updates two drawing coordinate sets under a separate audit. Both catalog steps remain under hash-addressed `catalog/history`; the current portable manifest binds the updated data. No publication images, source text or private filesystem paths are distributed.

The importer can reproduce this selected update when the original local evidence is available:

```shell
python scripts/import_manual_monomers.py --annotations /path/to/annotations --curation-catalog /path/to/sealed/catalog.json
```

It rejects altered selected graphs, changed source hashes, missing reviewed records and duplicate graph identities. An unchanged repeat is a no-op. The application reads the bundled catalog directly; it does not require the original research annotations.
