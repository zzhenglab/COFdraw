# COFdraw user guide

COFdraw is a Windows drawing application for finite COF pore patches and their monomers. Its drawing engine, catalog and rule-based description workflow work offline. An optional, explicitly configured LLM API can interpret descriptions; no language model is bundled. It includes an atom/bond editor, source-associated drawings, presentation styles and chemical exports. It performs no PDF OCR or training.

The application and executable icons use the supplied faceted teal logo with a transparent background.

## Run the release

1. Extract the **entire** `COFdraw-Windows.zip` to a writable folder. Keep `COFdraw.exe`, `COFdraw.exe.config`, its companion libraries and `data` together; copying the executable alone is insufficient.
2. Open `COFdraw.exe`. No Python installation is needed for the packaged release. The desktop window uses a local service bound to `127.0.0.1` and the included assets.
3. If the window cannot start, inspect `data/cofdraw.log`. This Windows build uses Microsoft WebView2 Runtime. For a computer without it, obtain the official standalone installer before moving the package offline. Edge browser installation alone is not the same prerequisite. [Microsoft distribution guidance](https://learn.microsoft.com/en-us/microsoft-edge/webview2/concepts/distribution), [pywebview requirements](https://pywebview.flowrl.com/guide/installation.html).

## Describe, review, draw

If an older downloaded release reports `Failed to resolve Python.Runtime.Loader.Initialize`, replace it with v0.6.1 or later. Windows can retain an Internet-zone marking on extracted DLLs; the included application configuration supports loading the bundled assemblies without removing those markings or changing machine settings. See [Windows startup troubleshooting](windows-startup.md).

Start in **Describe a COF**. Type a sentence such as `Draw TFB and TAPB with an imine linkage in pastel colors`, or choose a built-in example. Analyze the description, review the selected monomers, linkage, finite layout and style, then create the drawing. Names and words are matched to existing catalog records. The review shows assumptions and unresolved instructions; an unsupported or ambiguous request cannot become a ready drawing silently. Source-name matches open the qualified source recipe. A proposed design does not inherit the identity of a named material.

Rules mode needs no key or network. Expand the optional API settings to use an OpenAI-compatible language model with your endpoint, model and API key. Only an explicit API-mode analysis sends the description and compact catalog. The key remains transient and is excluded from saved projects and calculation reports. Review or edit the returned selections locally without another provider request. See [LLM API setup](llm-api.md) for endpoint requirements and supported behavior.

Use **Build manually** to choose or edit monomers directly. Advanced drawing and source controls remain available without crowding the initial description workflow.

## Draw and save

The guided library keeps **Building block A** and **Compatible partner B** visible at the top. Choosing A advances to compatible B choices. Use **Change** to revise either block, or **Single monomer** to inspect/edit one molecule. Search and a compact site-count filter stay visible; detailed explanations and source/master evidence expand when needed. See the [manual master-record review](manual-master-crosscheck.md) for the latest admission evidence.

The monomer library displays a scrollable grid of chemical structure cards, with three to six columns at desktop widths. Previews load as they come into view and are cached during the session. Use the panel to search abbreviations or full names and filter by 2c, 3c, 4c or 6c reactive sites. Load one monomer into the editor, or select A and then a compatible B. Incompatible choices remain visible with a reason and a disabled selection button. Reactive-site counts are distinct from molecular point-group symmetry. Choose a starting combination or a source drawing. Enter SMILES, import a single-molecule MOL, or use **Draw / edit** to change a monomer. Update its preview, select the supported linkage and connectivity template, then create the drawing. Inspect attachment ports and the connected structure before saving.

The v0.8.0 catalog contains **17 source-associated drawings and 118 monomers**, including recovered manual alignments and source-informed Pyurea/UCOF layouts. The latest 37 additions were checked against specific published figures by two agents; 20 have supported drawing roles and 17 remain references. The original 81 entries and their master/manual evidence are preserved. See the [literature coverage and limits](yaghi-literature-review.md). Choose A, then a compatible B; the visible A/B order stays fixed. **Show all** reveals unavailable choices with their reasons. Seven scoped transforms are available, including nitrone; a reference entry never implies an implemented constructor.

Monomer cards and checked structures display an offline **SA Score**, from 1 (easier) to 10 (harder). Each fresh drawing also shows its declared finite precursor scores and their unweighted mean/maximum. Editing the inputs clears the stale summary until a new drawing is rendered. This is RDKit's molecular heuristic, not ASKCOS, a synthesis route, or COF reaction/yield prediction. See [SA Score scope](synthetic-accessibility.md).

**Linkage and connectivity** starts in **Auto — match my monomers**. Validated monomer changes update the supported linkage and connectivity template. Choose **Manual — choose and modify** to retain your own compatible product-state and layout selections. An unsupported pair remains unsupported; selecting a reference example does not invent a constructor for it. The full COFex record, including precursor SMILES and qualified net/cat fields, appears above the canvas and can be copied.

Twenty styles use the same native chemical bond/label renderer. **Source drawing style** retains recorded colors, highlights and faded regions. Choose among 13 palettes or your own colors; adjust line weight, fonts, halo strength, highlight width, sphere size, and illustrative layer spacing/tilt. Styles preserve the underlying graph and coordinates.

Use **Move / adjust** to drag individual atoms on the displayed drawing, then **Render changes**. **Undo layout** restores the previous recipe. The block-position panel moves whole building blocks; **Clean bond geometry** optimizes their poses without deforming rings. Editing creates a design and clears the active source identity. Cleanup does not overwrite a reviewed source; explicit source graphs without block metadata retain their manual coordinates. Geometry diagnostics report bond-length and junction-angle concerns, not proof of a relaxed crystal.

New designs now fit both bond lengths and connection angles while preserving rigid monomer geometry. A visible review notice appears when the geometry checks are not satisfied. TpPa and TPB-DMTP use explicit coordinates chosen by comparison with their original figures. See [drawing geometry and comparisons](drawing-geometry.md).

**Measure pore** reports a calibrated 2D drawing-plane opening estimate where a suitable noncrossing macrocycle exists. A double-ended arrow marks the sampled inscribed-circle diameter on a flat drawing, with a nanometer label and **2D estimate** qualification. Use **Save annotated SVG** or **Save annotated PNG** in the measurement card to include it in a figure; ordinary exports remain unchanged. Editing or changing style clears the annotation. The SVG also retains the measurement method and scope. This is not a crystallographic pore diameter; see [drawing geometry and measurement](drawing-geometry.md).

## 3D and properties

Open **3D & properties** to generate layers from the current drawing. Changing the drawing invalidates the previous drawing-derived model; generate again to view the new structure. Representation and color controls change appearance independently of the coordinates. Drag to rotate, scroll to zoom, and select two atoms for a coordinate distance.

Compare AA, AB, ABC, progressive-slip and twisted layer illustrations using the same finite drawing and spacing. Slip is entered in Å and twist is a rotation per layer. Continuation boundaries remain explicit. These are illustrations and cannot supply periodic pore-volume analysis.

Quick checks and finite Debye scattering work directly on the displayed finite atoms, without a CIF. In the check, diffraction or geometry tab, **Find a periodic design for these monomers** checks exact precursors, linkage, product fragments and unchanged connectivity. Supported pairs include default TFB/TAPB, TFB/Pa and TAPB/terephthalaldehyde. Choosing a match generates a separate unrelaxed periodic design with a unit cell; it does not convert source drawing coordinates into a crystal. Unsupported designs keep finite checks and can import a corresponding CIF.

For independent periodic work, select one of five imine templates or import a CIF. The templates include a 3D diamond-like network. Set repeats along a/b/c and export CIF, extended XYZ, PDB, JSON or a view PNG. Single-layer 2D periodic cells can compare AA, AB, ABC and progressive-slip PXRD hypotheses; arbitrary twist requires a commensurate cell. Finite illustrations also offer a twist comparison. Interpenetration hypotheses use explicitly translated periodic copies with group colors and collision diagnostics. A copy count alone does not establish entanglement or a COFex `cat` value. Keep the JSON for copy assignments and inspect contacts.

In **Level 1 · geometry**, set the probe radius and grid resolution for the loaded periodic structure. The deterministic method reports PLD/LCD bounds, transport-connected probe-center volume and a sampled connected surface. Read the effective resolution and diagnostics. **Unresolved** means the coarse lower graph has not found a channel but the upper bound still permits one; it does not mean zero physical accessibility. These are hard-sphere geometric properties, separate from gas adsorption or experimental BET. See [geometric definitions and limitations](geometric-properties.md).

In **Level 2 · adsorption**, a separate pinned RASPA3 CPU worker simulates an explicit periodic model. Choose gases, composition, temperature, pressures, cycle counts and charge assumptions, then run and monitor the job. Cancellation preserves completed stages. Save the JSON/CSV and ZIP of actual simulation inputs/output. Default short runs need convergence checks; Henry coefficients use separate Widom sampling, and optional zero-loading heat uses finite differences at nearby temperatures. A rigid classical model does not predict chemisorption, framework flexibility or experimental performance. See [RASPA methods and assumptions](raspa-methods.md).

In **BET from isotherm**, paste two columns: relative pressure `p/p0` and uptake in `mmol/g`. Choose a fit range and adsorbate cross-section. Review the fitted area, C, capacity, R² and consistency checks before reporting. The app does not turn geometric area into BET.

In **Monomer descriptors**, choose A or B and a conformer count. RDKit supplies molecular descriptors and empirical Gasteiger charges, generates ETKDG conformers, and minimizes with MMFF94s or UFF when supported. Export XYZ/SDF and the method report. Energies compare conformers of the same monomer; they are not xTB, DFT or orbital energies. See [scientific methods](scientific-methods.md) for detailed assumptions and limitations.

## Materials inspection and evidence

Use **Quick check**, **PXRD**, and **Compare CIF files** to inspect composition, density, cell geometry and inferred periodic networks, compare candidate CIFs, and calculate model diffraction. Lattice-axis views, hydrogen visibility and atom tables support inspection. Catalog information and saved reports remain available to agents through the local API; the former Find information panel is removed. Follow the [materials workflow guide](materials-workflows.md) for instructions and interpretation.

## Exports and identity

Save the project JSON to reopen monomer inputs, recipe, layout and style. Use PNG/SVG for presentation, and MOL/SDF/CDXML for chemical editing. CDXML depends on the packaged native writer; optional failure leaves the other outputs available and is recorded in the chemical-export diagnostic. Current CDXML carries the chemical graph, coordinates and explicit wavy external connection points, **not all SVG fonts, colors or aliases**.

The fixed CDXML exports were actually opened and saved with **ChemDraw 23.1.2.7** in nine hidden, isolated document tests: four source pores, an automatic pore, isotope/E–Z fixtures and boundary fixtures. The saved CDXML retained the checked chemical tables and connection points; all nine native PNGs were visually inspected. Production export validates XML directly to avoid a native-reader stall on larger pores. The fix replaces invalid `Element="0"` nodes with documented external connection points, preserving boundary bond order without adding carbon or hydrogen. See [Chemical export validation](chemdraw-export-validation.md) for evidence and limits.

The current update also fixes the reported corner-only/B4 display: exports are centered inside one explicitly sized Letter page with room for labels and boundary waves. Seven native ChemDraw cases retained their full graph inside the page, including the default, large DMTP and four-connected examples. See [page-fit validation](cofdraw-page-fit-validation.md).

**Keep ChemDraw-edited framework boundaries in CDXML.** ChemDraw's own MOL save removed external connection points in these tests. Retain COFdraw's original MOL/SDF plus JSON for the explicit boundary graph. A continuation dummy in those files may display as `R` or `*`; it is an outgoing boundary, not a terminal atom or a complete material formula. Direct compatibility checks cover the tested ChemDraw version and files; frozen-executable and browser checks are recorded in [the release validation record](cofdraw-release-validation.json).

For external RDKit workflows, large para-phenyl pores can trigger the RDKit 2026.09 symmetric-ring slowdown. The [scoped MOL reader workaround](rdkit-mol-import.md) preserves the literal graph without changing COFdraw's scientific defaults. Reopen a pore with its project JSON; the application's MOL input is intended for individual monomers.

The full expanded COFex notation, compact key, literature name and precursor structures are shown together when a matching source record is available. The identity is provisional and source-qualified. `unk` remains unknown, and the reference convention `cat0` is not independent evidence of noninterpenetration. Changing chemistry/connectivity creates a design and clears the source identity; style-only changes preserve it. Seven explicit construction templates cover 3c+2c, 3c+3c, 4c+2c, 4c+4c, 6c+2c, an open 2c+2c chain and a custom six-unit ring. The 25-row topology/reference browser also records rhombic/kagome, multicomponent, staged, cage and 3D examples; reference entries do not enable an unimplemented chemical constructor. A regular pore does not establish a unique net, stacking, crystal structure, synthetic feasibility or whole-COF SMILES. See [catalog identity rules](cofdraw-catalog.md).

## Development and automation

For source setup, tests, and Windows packaging, see [CONTRIBUTING.md](../CONTRIBUTING.md). For scripted drawing and export, read the [COFdraw agent skill](../skills/cofdraw/SKILL.md).
