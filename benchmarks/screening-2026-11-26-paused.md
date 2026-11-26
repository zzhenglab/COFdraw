# COF screening checkpoint, November 26, 2026

The release covers 45,439 hypothetical designs and 24,760 unique, unrelaxed periodic CIFs. First-pass methane and hydrogen properties are available for **22,050** models. The final coordinate cohort adds **2,780** newly screened models to the November 25 release. It also resolves 31 pending records as force-field coverage holds. [Dataset and checksums](https://github.com/zzhenglab/COFdraw/releases/tag/screening-2026-11-26-paused)

| Outcome | Models |
| --- | ---: |
| First-pass gas screen complete | 22,050 |
| CIF available, accessibility review needed | 2,510 |
| CIF available, geometry only | 146 |
| CIF available, force-field coverage hold | 54 |
| Construction rejected | 20,679 |
| **Total** | **45,439** |

## Methane delivery and hydrogen competition

All 2,995 models admitted from the final coordinate cohort have a disposition. Of these, 2,960 passed the force-field checks. For each passing model the rigid, zero-charge RASPA3 calculation at 298 K used two independent seeds. Pure methane delivery is the absolute loading at 65 bar minus that at 5.8 bar. The dilute CH₄/H₂ ratio is the mean of two seedwise ratios of Widom Henry coefficients. The CIF archive reports loading per mass, methane delivery per ideal-crystal volume, and the dimensionless Henry ratio. Files also report seed standard deviations and reasons for holds. [Earlier protocol and coverage](screening-2026-11-25.md)

Longer methane calculations repeated both pressures with three new seeds for 192 imine COFs. These included 64 controls selected independently of the initial delivery ranking and 32 models from each of four ranking quartiles. Across the 192, the mean absolute difference in delivery was **0.71 cm³(STP)/cm³** and the Spearman rank correlation was **0.999**. These figures describe numerical repeatability for the selected models under the same host and gas model.

| Equimolar CH₄/H₂ at 298 K | Completed native runs | Planned |
| --- | ---: | ---: |
| 10 bar, short sampling | 8,880 | 8,880 |
| 10 bar, long sampling | 8,880 | 8,880 |
| 1 bar, long sampling | 1,079 | 8,880 |
| 30 bar, long sampling | 1,059 | 8,880 |

At 10 bar, the short mixture calculation and the dilute Henry ratio placed 192 of the same COFs in their respective top 296 among 2,960 imine models. Their Spearman correlation was 0.953. The independent long-reference labels remain sealed for prospective evaluation.

The 1- and 30-bar series paused on November 26 with 2,138 verified native archives and zero recorded failures. Each model has three independent seeds per pressure. At the pause, 352 models had all six runs. The public snapshot reports aggregate progress for this series. Blank material-level gas fields denote unavailable observations.

## Download and use

The CIF ZIP contains `results.csv`, `methods.json`, and `coordinate-first-pass.csv.gz` with all 2,995 dispositions in the final coordinate cohort. Join external calculations using the candidate ID and exact CIF SHA-256. Import `research-screen.json.gz` under **Explore COF chemical space → Simulation results → Open results file**. The earlier release supplies the original force-field and temperature comparisons. [Snapshot downloads](https://github.com/zzhenglab/COFdraw/releases/tag/screening-2026-11-26-paused) · [November 25 checkpoint](screening-2026-11-25.md)
