# COF screening data · November 24, 2026

The 08:44 Central checkpoint contains 39,096 design outcomes from the 45,439 eligible two-dimensional construction attempts. Use the ZIP for ordinary CIF files and a CSV property index. The compressed JSON opens in COFdraw 0.10.7 on Windows and in the corresponding Mac beta. Both files come from the same snapshot. [Dataset downloads and SHA-256 checksums](https://github.com/zzhenglab/COFdraw/releases/tag/screening-2026-11-24)

| Outcome | Designs |
| --- | ---: |
| Periodic CIF with a completed first-pass gas screen | 19,090 |
| Periodic CIF held for pocket accessibility review | 2,193 |
| Periodic CIF with geometry only | 127 |
| Periodic CIF held for force-field coverage | 19 |
| Construction rejected | 17,667 |

The 21,429 unique CIFs are generated, unrelaxed AA-stacked sheets. Each CSV record identifies the two precursor SMILES, the proposed net template and the exact CIF hash where coordinates exist. Empty property fields indicate unavailable results. Reported precursor pairs can lead to multiple products and packings.

## Adsorption results

At 298 K, 1,133 models have three-seed, seven-pressure absolute methane isotherms and equimolar CH₄/H₂ mixture simulations at 10 bar. Median methane delivery over 65 to 5.8 bar is **111.6 cm³(STP)/cm³ of ideal crystal**, and median mixture selectivity is **4.43**. These medians describe the selected models. An additional 737 designs have three-seed delivery results under a shorter endpoint protocol. The temperature study has resolved apparent dilute adsorption heats for both gases in 5,215 models. The explorer and CSV carry the property units and simulation settings.

The long study completed 33,585 native jobs, with 354 planned jobs held on geometry grounds and zero execution errors. It yielded 325 matched packing comparisons. The temperature and two-force-field campaigns are still running at this checkpoint. The explorer admits a candidate force-field comparison once both host profiles have all three matched seeds at both endpoint pressures.

At 35 bar and 298 K, the seven experimental controls already have three matched seeds per force field. For COF-5, the baseline and UFF-host models give excess methane uptakes of 95.8 and 123.0 mg/g, respectively, against 89 mg/g measured. Both used the same CIF, methane guest and fugacity convention. These are rigid-host predictions. Packing, inaccessible pockets and sample activation affect their interpretation. [Experimental measurements](https://doi.org/10.1021/ja9015765)
