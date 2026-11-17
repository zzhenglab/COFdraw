# Changelog

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
