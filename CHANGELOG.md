# Changelog

## 0.10.0

- Use the Framework canvas, rendering inspector and single Drawing extent menu for both planar and three-dimensional COFs.
- Display a calibrated pore arrow together with the precursors and reaction arrow. Add a common name or building-block formula to the exported figure, and double-click the label to edit it.
- Explore topology-specific neighborhoods and depth-faded circuits. Browse starting structures by chemistry and dimension.

## 0.9.8

- Add the selected precursors and a horizontal or vertical reaction arrow to the current drawing. SVG and PNG exports include the scheme.
- Choose exactly three connected pores on supported frameworks. The same graph supplies the drawing and its molecular exports.
- Build a periodic framework from eligible tetrahedral pairs directly from the standard editor. Draw depth-faded chemical circuits and inspect networks and their topology in the periodic viewer.
- Copy the precursor-based COFex below a drawing or paste a supported record to reconstruct a design. The rendering inspector uses a compact list.

## 0.9.7

- Draw one pore or its connected neighbors directly on the 2D canvas. SVG, PNG and molecular exports follow the visible extent. Project files remember the extent and pore colors.
- Explore seven one-precursor examples, including boroxine COF-1, triazine CTF-1 and two imine frameworks with both triangular and hexagonal pores.
- Explore six three-precursor examples. SIOC-COF-1 and SIOC-COF-2 expose both triangular environments and their alternating-linker hexagon. Change a library precursor or draw your own to create a design.

## 0.9.6

- Generate periodic layers from four reviewed 2D source drawings. The CIF and model report identify the selected precursors, reaction, periodic template and assumed stacking.
- Generate diamond and hexagonal-diamond frameworks from editable tetrahedral and ditopic precursors with imine, reduced-imine or acylhydrazone bonds. Named starting combinations use the same builder as custom SMILES. Choose a network count and inspect the resulting contact checks.
- Draw a closed framework circuit with conventional bond-line geometry. Rotate either reacted building block, highlight a network, superimpose the building-block graph or display a probe-center pore surface in the 3D viewer.
- Explore locally imported screening results with model identifiers, CIF hashes, and separate first-pass and confirmation results.

## 0.9.5

- Crystal examples open on the 2D drawing canvas. Open the same periodic coordinates in 3D, rotate the packing, highlight one network, or export the unit cell.
- Switch between reported linkers in four literature series or choose a checked TAPM diamond-net design. Qualified aryl-F derivatives retain their parent cell and have a separate design identity.
- Apple silicon and Intel Mac betas use the shared desktop application and ship separately from the Windows ZIP.

## 0.9.4

- Open COF-300, COF-320 and COF-103 crystal figures from Starting combination. Switch between a vector projection and the rotatable coordinates. Inspect traced junctions, molecular circuits and periodic packing where supported.
- Added four source-reviewed framework precursors and excluded radical fragments from automatic pairing.

## 0.9.3

- Explore linkage families, DOI-linked precursor pairs and their chemical structures. Open compatible pairs in the drawing studio or export a family screening batch.
- Compare ideal repeat chemistry across 79,064 family–pair combinations and geometric properties for 190 periodic models, including 143 screened hypothetical AA cells. Electronic and adsorption targets carry explicit pending flags.
- Draw connected repeat units and self-condensed or sequence-defined mixed-linker pores. Explore the active structure in the integrated 3D workspace.

## 0.9.2

- Added 128 finite building blocks from COFvision inventory 0110, with source figures identified by DOI, page and crop hash.
- Added translucent half-space-filling back layers, a half-space-filling 3D view, and element-colored drawings with selectable functional-group hydrogen labels.
- Free rotation and reset update the axis menu. Atom measurements require two distinct atoms and clear when measurement mode is switched off.

## 0.9.1

- Drag the horizontal separation slider to move from overlaid stacking patterns to individual 2θ axes. The exported SVG keeps the chosen layout.
- Search starting combinations by COF, building block or linkage. Auto includes compatible designs across linkage families and gives recorded COF drawings their own entries.
- Comparing again after viewing an imported model’s AB or ABC arrangement reuses the original layer.

## 0.9.0 - 2026-11-18

- Fuse five- or six-membered rings onto an existing outer bond with a placement preview. Ring double bonds use a shorter inner stroke. The sketch shows terminal amine and hydroxyl hydrogens and preserves radical electrons. Hover to change atom charge or delete an item.
- Open drawing and analysis tools from the workspace header. Use the arrow keys in the structure laboratory to switch tabs.
- Calculate a two-monomer weigh-out using linkage-specific functional groups, purity and chosen site equivalents.
- Compare measured PXRD positions with the current periodic model and export candidate reflections, unmatched centers and model assumptions. Inspect the pattern at full canvas width and save an editable SVG.
- Malformed BET input now identifies the affected line.

## 0.8.6 - 2026-11-18

- Build periodic isotope-labelled designs and organic sulfur, phosphorus, silicon and halogen fragments with checked geometry. Preserve isotope masses in crystal analysis, CIF and extended XYZ.
- Open in the manual builder. Describe materials by name or molecular core and optionally review a proposal using an OpenAI, Anthropic or custom API key.
- Vary element colors in randomized styles and view pore faces, layer profiles and mixed representations with verified periodic bonds.

## 0.8.5 - 2026-11-18

- Recover crowded periodic starting cells using phenyl rotations and restrained fragment geometry. Report the remaining clashes and missing or unintended bonds if validation fails.
- Screen 46 solvents. Add a scaffold-tested descriptor fallback, measured reference compounds, and published preparation amounts for PyTTA, TAPP and terephthalaldehyde.

## 0.8.4 - 2026-11-18

- Include the offline solubility measurements when the agent helper creates a scratch runtime.

## 0.8.3 - 2026-11-18

- Screen room-temperature monomer solubility across 19 solvents with reported measurements, estimates with error bands, and mass per chosen solvent volume. Show literature preparation loads when available. Export CSV and JSON.
- Add DHP, HAPT, and the figure-based HBC-LA12 reference. Find existing monomers by TAM and Trip-CHO, with source evidence in their details.
- Add a periodic picker for common elements, directed coordination bonds, bond-order shortcuts, snapped chain drawing, and undoable Ctrl+Shift+K cleanup. Use a three-blade triptycene depiction in monomer previews.

## 0.8.2 - 2026-11-18

- Move Appearance above Geometry and put Randomize style beside Measure pore. The pore arrow shows its diameter in nanometers and the method is expandable.
- Add Gradient layers and Spheres and bonds, four palettes, and coordinated random style combinations.
- Open a periodic CIF for supported current drawings. Assemble unrelaxed hexagonal and square models from the chosen precursors and connectivity when there is no bundled match. Keep a finite preview available for drawings that need a separate crystal model.

## 0.8.1 - 2026-11-18

- Use symmetric porphyrin ring depictions in monomer previews and generated COF drawings, preserving inner N-H atoms and saved source coordinates.
- Remove unused screenshots and superseded validation files; consolidate duplicate crystal-model tests while retaining their chemical and network checks.

## 0.8.0 - 2026-11-17

- Check the displayed finite atoms and calculate finite Debye scattering without a CIF. Keep periodic density, pore geometry and Bragg PXRD distinct.
- Match exact precursors, product fragments and unchanged connectivity to validated unrelaxed periodic designs. Add the default TFB/TAPB cell, and reject stale results or explicit edited graphs that no longer match.
- Compare AA/AB/ABC/slipped periodic PXRD hypotheses with contact and dimensionality checks; preserve skew cells and atom charges. Finite stacks additionally support illustrative twist.
- Enlarge the square drawing canvas, keep full COFex notation visible, collapse secondary facts, and preserve the user's visible A/B order during automatic chemistry selection.
- Add local conversational follow-ups, suggestions, review/apply actions and explicit current-drawing context; optional API mode remains catalog-bound and credential-safe.
- Expand to 118 monomers with 37 finite graphs checked against published figures by two agents. Twenty additions have supported drawing roles; 17 remain references. Preserve the previous 81 entries and all source recipes.
- Add a scoped nitrone constructor with retained N-O connectivity, explicit N+/O- charges and unassigned new E/Z stereochemistry.
- Preserve internal porphyrin N-H atoms through aromaticity restoration, property-cache refresh and chemical export.
- Add five reference linkage families and a searchable literature coverage ledger with viewed-page evidence, ambiguous cases and access limits. The review does not claim exhaustive main-paper/SI coverage.

## 0.7.0 - 2026-11-17

- Add deterministic periodic pore geometry with PLD/LCD bounds, transport-accessible surface and probe-center volume, density, void fraction and resolution diagnostics.
- Integrate a local RASPA3 CPU worker for five gases, gas mixtures, pressure-series GCMC and separate zero-loading Widom calculations. Save force fields, charges, reservoir settings, seeds, block errors and actual run files.
- Estimate zero-loading isosteric heat from temperature-dependent Henry coefficients; distinguish this estimate from finite-loading heats and experimental adsorption.
- Add offline RDKit SA Score for every catalog monomer and declared finite precursors in new drawings; scores are molecular heuristics, not synthesis routes or whole-COF feasibility predictions.

- Correct distorted linkage angles and inconsistent aromatic-edge alignment in new designs with a rigid-block depiction fit; preserve each monomer's internal geometry and the molecular graph.
- Compare TpPa and TPB-DMTP reconstructions against their original figures and retain the selected source-informed coordinates. Existing exact source coordinates remain unchanged.
- Filter monomer B by the selected A, refresh linkage/layout choices after edits, and prevent incompatible pairs from creating a drawing.
- Expand the reviewed monomer catalog from 49 to 81: 25 additional buildable monomers and seven reference-only entries, with source and installed-master evidence. Keep provisional agent review distinct from completed human ground truth.
- Show geometry review notices when cleanup cannot satisfy the drawing checks. Retain manual adjustment, source qualifiers, and the downloaded-Windows startup fix.
- Refresh 3D views from the active drawing, expose full COFex notation, and make automatic linkage/connectivity selection explicit with manual overrides.
- Draw pore-size arrows on the figure; add more 3D representations, layer coloring and AA/AB/ABC/slipped/twisted comparison options.
- Add explicit translated-network hypotheses with collision diagnostics and group coloring. Copy counts do not assign interpenetration or COFex cat values.
- Remove the Find information panel from the desktop; catalog querying remains available through the agent API.

## 0.6.1 - 2026-11-17

- Fix downloaded Windows releases failing at `Python.Runtime.Loader.Initialize`: include the process-scoped .NET loading configuration beside the executable.
- Preserve Windows download markings and machine security settings; no blanket file unblocking or administrator launch is required.
- Add native startup validation from a freshly extracted, download-marked package under a normal Windows user token.
- Return a controlled startup error and exit status instead of a second unhandled PyInstaller traceback dialog.

## 0.6.0 - 2026-11-17

- Add a description-to-drawing workflow with built-in sentence examples, local catalog matching and an editable parameter review.
- Resolve existing monomer and source names, supported linkage/layout choices and appearance words; show assumptions, ambiguity and unsupported instructions explicitly.
- Add optional OpenAI-compatible API parameter selection with local chemical validation and transient credentials.
- Simplify the desktop starting flow into Describe and Build manually modes with progressive disclosure of advanced controls.
- Expose description catalogs and request schemas to agents and extend the headless helper and skill.

## 0.5.0 - 2026-11-17

- Add quick periodic-model composition, density, cell and element-site inventories.
- Infer periodic network dimensionality with explicit translated-cell bonds, coordination, contact warnings and cutoff sensitivity.
- Calculate model powder XRD with IT92 form factors, physical reciprocal-cell geometry, peak tables and CSV exports.
- Compare up to eight CIFs using common screening settings.
- Add lattice-axis views, hydrogen visibility and fractional atom tables.
- Search local COF evidence by names, monomers, linkage, DOI and workflow; preserve support status and source qualifiers.
- Expose live API request schemas and saved reports with input hashes, parameters and engine versions; expand the agent skill and headless client.
- Preserve CSV payloads and protect spreadsheet names; skip unreadable historical report files.
- Include Gemmi coefficient provenance, MPL-2.0 license and editable coefficient source in the portable package.

## 0.4.0 — 2026-11-17

- Preserve source geometry and use a shared native chemical renderer across styles. Restore compact, skewed Pyurea and asymmetric UCOF layouts without distorting aromatic rings.
- Recover reviewed manual annotations: 17 source drawings and 49 monomer presets.
- Add source colors/halos, molecular-sphere illustration, 13 palettes and custom colors; 20 styles total.
- Add direct atom adjustment, block editing, undo, rigid-block cleanup and geometry diagnostics.
- Automatically update supported linkage and connectivity choices after monomer edits.
- Add local 3D viewer, four explicit periodic templates, CIF import, supercells and CIF/XYZ/PDB exports.
- Add geometric pore-size/volume/surface screening and separate BET fitting from supplied isotherms.
- Add RDKit monomer descriptors, Gasteiger charges and optimized 3D conformer export.
- Update agent skill, methods documentation and Windows package; retain explicit scientific scope.

## 0.3.0 — 2026-11-17

First public GitHub release of COFdraw's Windows desktop application.

- Offline drawing with 41 monomer presets, four source-associated pores, seven construction templates, and 17 presentation styles.
- PNG, SVG, MOL, SDF, CDXML, and editable JSON exports.
- ChemDraw continuation-point and page-fit fixes, with recorded native compatibility checks.
- Faster large-pore depiction for the measured DMTP–TAPB example.
- Public desktop source, usage documentation, screenshots, and a portable agent skill.
- Portable Windows package with required libraries, reviewed catalog data, third-party notices, and SHA-256 checksum.

The public package preserves the validated 0.3.0 executable. Packaging omits prior generated drawings and unrelated research handoff files. Detailed validation and limits are recorded in `docs/cofdraw-release-validation.json`.
