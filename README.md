<p align="center"><img src="docs/assets/cofdraw-logo.png" alt="" width="72"></p>

# COFdraw

Draw covalent organic frameworks from their building blocks, explore periodic crystal models, and compare simulated properties. The desktop app opens in a 2D chemical drawing workspace. You can start with monomers, a COF name, a COFex record, or a description. The drawing engine and catalog run locally.

**[Download for Windows](https://github.com/zzhenglab/COFdraw/releases/download/v0.10.7/COFdraw-Windows.zip)** · **[Download for Mac (beta)](https://github.com/zzhenglab/COFdraw/releases/tag/v0.10.7-mac-beta.21)** · **[Download the screening dataset](https://github.com/zzhenglab/COFdraw/releases/tag/screening-2026-11-24)**

![COFdraw 0.10.7 with the COF-953 precursor in the left editor and a kagome patch on the canvas. The central hexagonal pore is purple and its six triangular neighbors are green.](docs/assets/cofdraw-workspace-v0107.png)

*Try this drawing:* In **Build manually**, choose **One monomer · self-condensation**, select **COF-953 · Kagome**, then choose **Hexagon + 6 triangles** under **Drawing extent**. Each triangular pore shares an edge with the central hexagon.

## Run the desktop app

| Platform | First launch |
| --- | --- |
| Windows x64 | Extract the **entire ZIP** to a writable folder and open `COFdraw.exe`. Keep `COFdraw.exe.config`, `_internal`, and `data` beside it. The package includes Python. [Windows startup help](docs/windows-startup.md) |
| macOS beta | Choose **arm64** for Apple silicon or **x86_64** for Intel. Extract the ZIP, move **COFdraw.app** to Applications, and follow the first-open instructions for your build. [Mac beta downloads and tested macOS versions](https://github.com/zzhenglab/COFdraw/releases/tag/v0.10.7-mac-beta.21) |

## From building blocks to a figure

| If you want to… | In COFdraw |
| --- | --- |
| Draw or edit a COF | Choose a starting combination or draw monomer A and pick a compatible B. Automatic matching suggests a linkage and a layout. **Drawing extent** shows the neighboring pores supported by that topology. |
| Work from a name or a description | Open **Describe a COF** and try “Draw TFB and TAPB with an imine linkage in pastel colors.” Review the matched structures before applying. Local rules work offline. GPT and Claude are optional. [API setup and data handling](docs/llm-api.md) |
| Prepare a figure | Add precursor structures and a reaction arrow, label the COF, measure the drawn opening, or adjust the palette and geometry. Export SVG or PNG for a figure, CDXML for ChemDraw, MOL/SDF for the finite chemical graph, or project JSON to resume editing. [Chemical export guide](docs/chemdraw-export-validation.md) |
| Examine a crystal | Open **3D & properties**. Build an unrelaxed periodic model for a supported design or import a CIF. Inspect the cell, compare stacking and PXRD, and screen pore dimensions, accessible volume, and surface area. [Materials workflows](docs/materials-workflows.md) |
| Compare chemistry or adsorption | Inspect monomer descriptors and synthetic accessibility, or run local RASPA adsorption calculations with documented force-field and charge settings. BET fitting uses a supplied adsorption isotherm. [Scientific methods](docs/scientific-methods.md) |

<details>
<summary>See the geometry and layer-stacking examples</summary>

![Before-and-after examples of drawing geometry with rigid building blocks and adjusted linkage angles](docs/assets/cofdraw-geometry-v062.png)

![Five illustrations of layer stacking, colored by layer](docs/assets/cofdraw-stacking-v062.png)

</details>

## Explore the screening dataset

The **November 24, 2026** release records **39,096 design outcomes**. It includes **21,429 unique periodic CIFs** for generated, unrelaxed AA-stacked sheets. **19,090** have a first-pass gas screen. Each coordinate file is linked to its precursor SMILES and proposed net through the CSV index. The index also records rejected builds and held calculations. Empty property fields mean the result is unavailable. [Coverage, methods, and validation](benchmarks/screening-2026-11-24.md)

| Download | How to use it |
| --- | --- |
| [CIFs and property index (ZIP)](https://github.com/zzhenglab/COFdraw/releases/download/screening-2026-11-24/COFdraw-screening-CIFs.zip) | Open `results.csv` and use the `cif` column to find a model. The `cif_sha256` column identifies its exact coordinates. |
| [COFdraw results file (JSON.GZ)](https://github.com/zzhenglab/COFdraw/releases/download/screening-2026-11-24/research-screen.json.gz) | In **Explore COF chemical space → Simulation results → Open results file**, select the compressed file directly. |

![Methane simulations for the 1,133-model subset. A scatter plot relates crystal density to volumetric delivery. A second panel shows median uptake and the 10th to 90th percentile at seven pressures.](docs/assets/methane-screen-2026-11-24.png)

*Reading the figure.* Each point on the left is one of the **1,133 selected models** with a three-seed, seven-pressure methane isotherm at 298 K. Median predicted delivery from 65 to 5.8 bar is **111.6 cm³(STP) per cm³ of ideal crystal**. On the right, the shaded band shows variation across frameworks at each pressure. Both panels use the downloadable CSV.

**Interpreting a design.** A drawing records a finite patch and its continuation points. A COFex string identifies precursors and carries the assigned net and catenation when known. Periodic properties depend on the particular CIF, packing, and simulation settings. A literature precedent for a precursor pair does not establish that a proposed crystal has been synthesized. [Pore definitions](docs/geometric-properties.md) · [Adsorption assumptions](docs/raspa-methods.md)

## Use with an AI agent

Paste this request into an agent with access to GitHub and local tools:

> Use https://github.com/zzhenglab/COFdraw and read its SKILL.md. Draw the catalog TpPa pore in Standard style. Save an SVG, a ChemDraw file, and an editable project with its provenance.

The [COFdraw skill](skills/cofdraw/SKILL.md) covers catalog queries, headless drawing, exports, periodic structures, and retrieving reports.

## Run from source

For Windows development, use Python 3.11 from the repository root.

```powershell
py -3.11 -m venv .venv
.\.venv\Scripts\python.exe -m pip install -r requirements-cofdraw-desktop.txt
.\.venv\Scripts\python.exe scripts\launch_cofdraw.py --data-dir data
```

[User guide](docs/cofdraw-desktop-release.md) · [Contributing and tests](CONTRIBUTING.md) · [All releases](https://github.com/zzhenglab/COFdraw/releases) · [Report an issue](https://github.com/zzhenglab/COFdraw/issues)

## License and acknowledgments

COFdraw application code is distributed under the [MIT License](LICENSE). The desktop package includes third-party notices. The Gemmi IT92 coefficients retain their [MPL-2.0 license](docs/third-party-licenses/Gemmi-MPL-2.0.txt), the RDKit SA Score retains its [Novartis BSD notice](docs/third-party-licenses/RDKit-SA-Score-BSD.txt), and periodic templates credit [pyCOFBuilder](https://github.com/lipelopesoliveira/pyCOFBuilder). Catalog records retain their publication identifiers. [RASPA and force-field attribution](docs/raspa-methods.md)
