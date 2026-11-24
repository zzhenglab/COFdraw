# COF screening checkpoint · November 24, 2026

This snapshot freezes the committed calculation results when the campaign was paused at 11:46 a.m. Central. It includes **39,096 design outcomes**, **21,429 unique generated periodic CIFs**, and **19,090 first-pass gas screens**. The updated release adds temperature-dependent Henry coefficients and paired force-field results that passed the existing reporting rules. The earlier snapshot remains available for reproducibility. [Paused dataset and checksums](https://github.com/zzhenglab/COFdraw/releases/tag/screening-2026-11-24-paused) · [Earlier snapshot](screening-2026-11-24.md)

| Study | Results at pause | Remaining in the frozen plan |
| --- | ---: | ---: |
| Confirmation | 11,349 / 11,349 completed | 0 |
| Long GCMC, controls, and packing | 33,585 completed and 354 geometry holds / 33,939 | 0 |
| Additional first-pass construction and gas queue | 13,439 / 13,439 resolved | 0 |
| Temperature-dependent Widom sampling | 101,961 / 108,000 completed | 6,039 |
| Paired host force-field comparison | 2,247 / 3,246 completed | 999 |
| Remaining eligible coordinate builds | 0 / 6,343 attempted in this queue | 6,343 |

The data include **1,133 seven-pressure methane isotherms and equimolar CH₄/H₂ screens** plus **737** additional endpoint delivery results. Both gases have resolved apparent dilute adsorption heats for **5,664** candidates. **14** candidate CIFs have complete, matched force-field triplicates at both methane delivery pressures. The seven experimental controls form a separate 35 bar comparison. Incomplete groups are excluded from these comparisons.

## Why 19,090 of 39,096 designs have gas results

| Outcome | Designs |
| --- | ---: |
| Passed the geometry and accessibility checks and completed the first-pass gas screen | 19,090 |
| CIF available, pocket accessibility needs review | 2,193 |
| CIF available, geometry only | 127 |
| CIF available, force-field coverage needs review | 19 |
| Construction rejected | 17,667 |
| **Total recorded outcomes** | **39,096** |

The frozen inventory contains **79,064 hypotheses**. The current two-dimensional workflow admits **45,439** to construction and holds 30,222 chain hypotheses, 2,378 unsupported topologies, and 1,025 unresolved precursor identities. The public CIFs are generated, unrelaxed AA-stacked models. Even for reported precursor pairs, each modeled topology and packing remains a hypothesis.

The remaining breadth work starts with the **6,343 unbuilt eligible designs**, then reviews the **2,339 existing CIFs** that lack first-pass gas results. Rejected models need a changed structural hypothesis or a validated replacement CIF. After model checks, broad inexpensive screening and a stratified set of longer calculations can answer different questions. Blank fields in the release represent unavailable properties.

## Reuse the data

Use `results.csv` in the CIF ZIP to match precursors, net template, geometry, and adsorption results by candidate ID and CIF SHA‑256. The companion `research-screen.json.gz` opens in COFdraw under **Explore COF chemical space → Simulation results → Open results file**. Both archives contain the same outcomes and method metadata. For citation, record the release tag and exact CIF hash. [Adsorption methods](../docs/raspa-methods.md) · [Pore calculations](../docs/geometric-properties.md)
