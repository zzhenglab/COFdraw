# Yaghi literature review and catalog evidence

The November 2026 review screened the [Yaghi laboratory's official publication index](https://yaghi.berkeley.edu/publications.html) for COF research, molecular weaving and related reticular reviews. It produced **37 additional finite monomer graphs**, each checked by two agents against specified source drawings. These are agent reviews, not human ground truth or experimental verification.

The portable [coverage ledger](../data/framework-drawing/catalog/yaghi-literature-coverage.json) records titles, DOIs, official links, source hashes, actually viewed pages, existing monomers, checked candidates and unresolved structures. Publication PDFs and images are not distributed with COFdraw.

| Review scope | Count |
| --- | ---: |
| Official-index entries screened by title, 2005–2026 | 294 |
| Selected research/review/background entries | 100 |
| Unique selected entries | 99 |
| Title-only exclusions, without full-text review | 194 |
| New finite graphs independently checked by two agents | 37 |
| Exhaustive main-paper-plus-SI audits claimed | 0 |

The repeated entry is the same introductory textbook excerpt. Its contents/glossary pages were inspected; the full book's COF chapters were not. Ten early background entries have abstract/text scope rather than a full figure review. For COF-999, the public publisher Figure 1 was available and inspected, while the main PDF and SI remained unavailable. Each ledger entry lists its own boundaries; a download or text extraction is never counted as a figure inspection.

Unresolved structures remain outside admission, including the THAn positional-name discrepancy, metal coordination complexes and counterions, isotopic mixtures, several extended or protected precursors, and figures requiring further exact graph work. The 194 title-only exclusions may contain incidental COF discussion. This review therefore supports the listed additions and references; it does not establish exhaustive chemical coverage of the entire bibliography.

The [linkage reference supplement](../data/framework-drawing/catalog/literature-linkage-reference.json) adds ester transesterification, borosilicate, borophosphonate and methyltriazine–vinylene references with primary-figure evidence. A borazine reference has reviewed secondary figures and primary-abstract corroboration, with the original primary figures/SI explicitly unreviewed. **All five supplement records are references, without automatic construction rules.** Borophosphonate cubane and rod products remain distinct source-specific constitutions.

Runtime compatibility is determined separately by COFdraw's installed reaction registry and templates. A valid finite SMILES, reactive-site match or SA score does not establish a reported COF's identity, synthesis feasibility, topology, interpenetration or crystal coordinates. Source-specific post-synthetic modification and polymers are not collapsed into a generic two-monomer reaction.
