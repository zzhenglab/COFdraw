# Five-hour COF screening checkpoint

Local screening from 2026-11-20 18:31 UTC through 2026-11-20 23:31 UTC. Counts include jobs with a recorded outcome by the cutoff. The schedule prioritizes published precursor precedents and samples all six families. These proportions describe the selected sample.

| Linkage family | Recorded outcomes | Accepted cells | Completed gas screens | Build rejections | Execution errors or retries |
| --- | ---: | ---: | ---: | ---: | ---: |
| azine | 122 | 97 | 82 | 25 | 0 |
| beta-ketoenamine | 217 | 147 | 134 | 70 | 0 |
| hydrazone | 668 | 470 | 408 | 197 | 1 |
| imine | 672 | 437 | 392 | 235 | 0 |
| keto-hydrazide | 39 | 34 | 27 | 5 | 0 |
| nitrone | 122 | 93 | 90 | 29 | 0 |
| All | 1840 | 1278 | 1133 | 561 | 1 |

## Preliminary model results

| Observable | Median across completed gas screens |
| --- | ---: |
| CH4 delivery, mmol/g | 16.4 (n = 1133) |
| CH4 delivery, cm³ STP/cm³ of ideal crystal | 112 (n = 1133) |
| Infinite-dilution CH4/H2 Henry ratio | 4.77 (n = 1133) |

282 completed gas screens have a same-family precursor-pair reference in the frozen atlas. This evidence does not establish the phase, registry or performance of the proposed model.

114 selected designs completed the layer-spacing sensitivity test. The median span in connected probe-center volume across the three graph-preserving heights is 0.223 (n = 114) cm³/g. The selection includes controls and the top and middle methane-delivery candidates in successive groups of 20 completions.

## Calculation conditions

The initial structures are unrelaxed AA sheet models. CIFs are reimported and checked for periodic connectivity, close contacts and repeat mass. Initial pore geometry uses a 1.82 Å probe, 0.65 Å requested grid spacing and 64 surface samples per atom. Qualified cells receive two independent RASPA3 seeds per observable at 298 K. Pure methane delivery uses the 65 to 5.8 bar pressure swing with 500 initialization and 1,000 production cycles per seed. Dilute selectivity uses separate 3,000-cycle Widom runs for methane and hydrogen. Volumetric delivery uses ideal crystal density. The spacing test repeats geometry at the original AA height and at ±0.25 Å with a 0.45 Å requested grid and 128 surface samples per atom.

These are preliminary rigid-model predictions. Experimental calibration, polymorphism, activation and uncertainty in the force field affect interpretation. The public update contains aggregate results. Individual hypotheses and simulation archives remain in the local research collection.
