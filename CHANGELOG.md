# Changelog

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
