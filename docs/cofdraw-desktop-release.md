# COFdraw user guide

COFdraw is an offline Windows drawing application for finite COF pore patches and their monomers. It includes a basic atom/bond editor, source-associated starting drawings, presentation styles and chemical exports. It performs no PDF OCR, model inference or training.

The application and executable icons use the supplied faceted teal logo with a transparent background.

## Run the release

1. Extract the **entire** `COFdraw-Windows.zip` to a writable folder. Keep `COFdraw.exe`, its companion libraries and `data` together; copying the executable alone is insufficient.
2. Open `COFdraw.exe`. No Python installation is needed for the packaged release. The desktop window uses a local service bound to `127.0.0.1` and the included assets.
3. If the window cannot start, inspect `data/cofdraw.log`. This Windows build uses Microsoft WebView2 Runtime. For a computer without it, obtain the official standalone installer before moving the package offline. Edge browser installation alone is not the same prerequisite. [Microsoft distribution guidance](https://learn.microsoft.com/en-us/microsoft-edge/webview2/concepts/distribution), [pywebview requirements](https://pywebview.flowrl.com/guide/installation.html).

## Draw and save

The monomer library displays a scrollable grid of chemical structure cards, with three to six columns at desktop widths. Previews load as they come into view and are cached during the session. Use the panel to search abbreviations or full names and filter by 2c, 3c, 4c or 6c reactive sites. Load one monomer into the editor, or select A and then a compatible B. Incompatible choices remain visible with a reason and a disabled selection button. Reactive-site counts are distinct from molecular point-group symmetry. Choose a starting combination or a source drawing. Enter SMILES, import a single-molecule MOL, or use **Draw / edit** to change a monomer. Update its preview, select the supported linkage and connectivity template, then create the drawing. Inspect attachment ports and the connected structure before saving.

The current catalog supplies four source-associated drawings: TFB–TAB, TpPa, TPB–DMTP-COF and TpPa–SO3H/SCOF, plus eight design starting pairs. Its 41 monomers retain full names, paper-specific abbreviations and source-bound graphs. Six scoped automatic transforms are available: imine, beta-ketoenamine, acylhydrazone, Tp keto-enamine–hydrazide, azine and reduced-imine secondary amine. The last is a selected product state, not a universal one-step condensation. The 30 broader catalog entries include historical aliases and reference-only chemistry; an entry does not imply an implemented transform. See [v3 catalog scope](cofdraw-catalog-v3.md).

Seventeen presentation styles include fine or bold chemical lines, wavy/dashed edges, condensed substituents, colored atoms, building-block highlights, aromatic-ring and pore fills, pastel blocks, illustrative stacked layers, tilted sheets, and cartoons. Adjust fonts, palette, line weight, labels, edge marks, layer count and transparency. Random style retains its seed and resolved options. Layered views are illustrations, not calculated 3D crystals. Styles change the depiction, not the chemical graph. Stereospecified structures use Standard where custom styles cannot represent their stereo marks.

Save the project JSON to reopen monomer inputs, recipe, layout and style. Use PNG/SVG for presentation, and MOL/SDF/CDXML for chemical editing. CDXML depends on the packaged native writer; optional failure leaves the other outputs available and is recorded in the chemical-export diagnostic. Current CDXML carries the chemical graph, coordinates and explicit wavy external connection points, **not all SVG fonts, colors or aliases**.

The fixed CDXML exports were actually opened and saved with **ChemDraw 23.1.2.7** in nine hidden, isolated document tests: four source pores, an automatic pore, isotope/E–Z fixtures and boundary fixtures. The saved CDXML retained the checked chemical tables and connection points; all nine native PNGs were visually inspected. Production export validates XML directly to avoid a native-reader stall on larger pores. The fix replaces invalid `Element="0"` nodes with documented external connection points, preserving boundary bond order without adding carbon or hydrogen. See [Chemical export validation](chemdraw-export-validation.md) for evidence and limits.

The current update also fixes the reported corner-only/B4 display: exports are centered inside one explicitly sized Letter page with room for labels and boundary waves. Seven native ChemDraw cases retained their full graph inside the page, including the default, large DMTP and four-connected examples. See [page-fit validation](cofdraw-page-fit-validation.md).

**Keep ChemDraw-edited framework boundaries in CDXML.** ChemDraw's own MOL save removed external connection points in these tests. Retain COFdraw's original MOL/SDF plus JSON for the explicit boundary graph. A continuation dummy in those files may display as `R` or `*`; it is an outgoing boundary, not a terminal atom or a complete material formula. Direct compatibility checks cover the tested ChemDraw version and files; frozen-executable and browser checks are recorded in [the release validation record](cofdraw-release-validation.json).

The full expanded COFex notation, compact key, literature name and precursor structures are shown together when a matching source record is available. The identity is provisional and source-qualified. `unk` remains unknown, and the reference convention `cat0` is not independent evidence of noninterpenetration. Changing chemistry/connectivity creates a design and clears the source identity; style-only changes preserve it. Seven explicit construction templates cover 3c+2c, 3c+3c, 4c+2c, 4c+4c, 6c+2c, an open 2c+2c chain and a custom six-unit ring. The 25-row topology/reference browser also records rhombic/kagome, multicomponent, staged, cage and 3D examples; reference entries do not enable an unimplemented chemical constructor. A regular pore does not establish a unique net, stacking, crystal structure, synthetic feasibility or whole-COF SMILES. See [catalog identity rules](cofdraw-catalog.md).

## Development and automation

For source setup, tests, and Windows packaging, see [CONTRIBUTING.md](../CONTRIBUTING.md). For scripted drawing and export, read the [COFdraw agent skill](../skills/cofdraw/SKILL.md).
