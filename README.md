<p align="center">
  <img src="docs/assets/cofdraw-logo.png" alt="COFdraw" width="100">
</p>

<h1 align="center">COFdraw</h1>

<p align="center"><strong>A desktop drawing studio for covalent organic frameworks.</strong><br>
Draw source-aligned COF pores, explore 3D structures, and screen geometric pore properties offline.</p>

<p align="center">
  <a href="https://github.com/zzhenglab/COFdraw/releases/latest"><img src="https://img.shields.io/github/v/release/zzhenglab/COFdraw?color=207b68" alt="Latest release"></a>
  <img src="https://img.shields.io/badge/desktop-Windows%20x64-0078d4" alt="Windows x64">
  <a href="LICENSE"><img src="https://img.shields.io/badge/license-MIT-207b68" alt="MIT license"></a>
</p>

<p align="center">
  <a href="https://github.com/zzhenglab/COFdraw/releases/latest/download/COFdraw-Windows.zip"><strong>Download for Windows</strong></a>
  &nbsp; · &nbsp; <a href="docs/cofdraw-desktop-release.md">User guide</a>
  &nbsp; · &nbsp; <a href="skills/cofdraw/SKILL.md">Agent skill</a>
  &nbsp; · &nbsp; <a href="https://github.com/zzhenglab/COFdraw/issues">Report an issue</a>
</p>

![COFdraw drawing interface with a finite imine pore and export controls](docs/assets/cofdraw-desktop.png)

## Download and run

1. Download **[COFdraw-Windows.zip](https://github.com/zzhenglab/COFdraw/releases/latest/download/COFdraw-Windows.zip)** from the latest release.
2. Extract the **entire ZIP** into a writable folder.
3. Open **`COFdraw.exe`** inside the extracted `COFdraw` folder.

**No Python installation is needed.** Keep the executable, `_internal`, and `data` folders together. The ZIP is the ready-to-run desktop application; GitHub's automatically generated “Source code” archives are for developers.

The Windows application requires [Microsoft WebView2 Runtime](https://developer.microsoft.com/en-us/microsoft-edge/webview2/). Install it if the desktop window cannot open. Once the runtime is present, the bundled drawing tools and catalog work offline. The application runs its chemistry service locally on `127.0.0.1`.

Release assets include a SHA-256 checksum file. Troubleshooting and export details are in the [user guide](docs/cofdraw-desktop-release.md).

## What you can do

| Feature | Included in v0.4.0 |
| --- | --- |
| Monomer library | 49 presets, structure cards, and automatic compatible linkage/template selection |
| Reviewed drawings | 17 source-associated patches, including recovered manual alignments, Pyurea and UCOF |
| Geometry editing | Direct atom drag, block-position editing, undo, rigid-block cleanup, and diagnostics |
| Presentation | 20 styles, 13 palettes, custom colors, source highlights, halos, and illustrative layers |
| Figure / chemical export | PNG, SVG, MOL, SDF, CDXML, and editable JSON |
| 3D structures | Rotatable viewer, four periodic model templates, CIF import, supercells, CIF/XYZ/PDB export |
| Pore screening | Drawing-plane opening; periodic sampled pore size, void volume, and probe-exposed surface |
| BET | Fit supplied adsorption isotherms with consistency diagnostics |
| Monomer properties | RDKit descriptors, Gasteiger charges, ETKDG conformers, MMFF/UFF optimization, XYZ/SDF |
| Automation | Local HTTP API and reusable GitHub-discoverable agent skill |

Source drawings retain their recorded coordinates and presentation. The restored Pyurea and UCOF layouts use source-informed block positions rather than forcing a rectangular or regular-hexagonal arrangement. These are qualified reconstructions, not exact image registrations. See the [recovered drawing catalog](docs/recovered-source-drawings.md).

![Periodic COF model in the 3D structure laboratory](docs/assets/cofdraw-3d.png)

The six implemented transforms cover imine, beta-ketoenamine, acylhydrazone, Tp keto-enamine–hydrazide, azine, and a selected reduced-imine secondary-amine product state. Broader linkage and topology catalogs also contain **reference-only entries**; those entries do not imply an available constructor.

<details>
<summary>Preview the monomer library</summary>

![Searchable monomer library with chemical structure cards](docs/assets/cofdraw-library.png)

</details>

<details>
<summary>Recovered source alignment and highlights</summary>

![Source-aligned F-COF-1 with its recorded yellow halo and red/black bonds](docs/assets/cofdraw-source-style.png)

</details>

## Use with an AI agent

Paste this into an agent that can read GitHub and run local tools:

```text
Use https://github.com/zzhenglab/COFdraw to draw a COF pore.
Read the repository's SKILL.md and follow the COFdraw skill.
Start with the bundled TFB–TAB source drawing and export PNG, SVG,
MOL, SDF, CDXML, and an editable JSON project to a new output folder.
```

The repository provides a root [`SKILL.md`](SKILL.md), [`AGENTS.md`](AGENTS.md), and the installable [`skills/cofdraw`](skills/cofdraw/SKILL.md) skill. The skill explains setup, library inspection, headless drawing, export, and the limits of chemical interpretation. A pasted link makes these instructions discoverable; automatic installation and execution depend on the agent's capabilities and permissions.

## Run from source

Use Python 3.11 on Windows x64. From a clone of this repository:

```powershell
py -3.11 -m venv .venv
.\.venv\Scripts\python.exe -m pip install -r requirements-cofdraw-desktop.txt
.\.venv\Scripts\python.exe scripts\launch_cofdraw.py --data-dir data
```

The reviewed portable catalog is included in `data/framework-drawing`. Source publications, private datasets, model weights, and training environments are not needed for this drawing application.

See [CONTRIBUTING.md](CONTRIBUTING.md) for tests and packaging. The Python import namespace remains `cof_vision` because COFdraw originated in that project; this repository contains the desktop drawing application.

## Scientific scope

COFdraw creates **finite molecular drawings with explicit continuation points**. A chosen template or a regular pore does not establish a unique periodic net, crystal structure, stacking arrangement, synthetic feasibility, or whole-material SMILES. Layered drawing styles and finite stacked patches are illustrations. The separate 3D laboratory imports explicit periodic structures and includes four unrelaxed model templates generated with pyCOFBuilder; these do not automatically identify the active 2D drawing as a crystal. Source identity is qualified by the supplied reference and is cleared when chemistry or connectivity changes.

**Pore screening is geometric, not an adsorption simulation.** Sampled pore size and probe-exposed surface depend on the input model, atom radii, probe and sampling. BET requires supplied isotherm data; geometric area is never labeled BET. Monomer calculations use inexpensive RDKit methods; xTB and quantum orbital energies are not bundled. Read the [scientific methods and limits](docs/scientific-methods.md).

CDXML preserves the checked chemical graph and continuation points, but does not reproduce every presentation style. Compatibility was checked with **ChemDraw 23.1.2.7**. Keep ChemDraw-edited framework boundaries in CDXML: ChemDraw's MOL save omitted external connection points in the recorded tests. Keep original MOL/SDF and project JSON files as well.

See [v0.4 validation](docs/v0.4-validation.json), the [catalog scope](docs/cofdraw-catalog-v3.md), [chemical export validation](docs/chemdraw-export-validation.md), [page-fit validation](docs/cofdraw-page-fit-validation.md), [original release validation](docs/cofdraw-release-validation.json), and [public-package verification](docs/publication-validation.json) for evidence and limitations.

## License and acknowledgments

Original project code is available under the [MIT License](LICENSE). COFdraw uses RDKit, pywebview, FastAPI, and other open-source components; the Windows download includes their license notices. Generated periodic templates cite [pyCOFBuilder](https://github.com/lipelopesoliveira/pyCOFBuilder) and retain its attribution in `data/crystal-models`. Source publications and external materials retain their own rights. Catalog records preserve their source identifiers and provenance.
