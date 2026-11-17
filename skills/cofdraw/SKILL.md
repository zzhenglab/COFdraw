---
name: cofdraw
description: Use COFdraw to draw source-informed or designed finite COF patches, preserve manually placed coordinates, export editable chemistry, generate molecular conformers, inspect periodic CIF models, and run qualified geometry or supplied-isotherm analyses. Use when asked to use COFdraw or its GitHub link for these drawings, structures, or calculations.
---

# COFdraw

Create actual PNG/SVG and editable MOL/SDF/CDXML drawings with the local COFdraw application. Preserve the accompanying JSON and its provenance. Its separate scientific tools handle molecular conformers, explicit periodic models, geometric pore estimates and supplied adsorption data. A finite drawing does not establish a crystal structure.

## Start from a GitHub link

The project is [zzhenglab/COFdraw](https://github.com/zzhenglab/COFdraw). A pasted link makes these instructions readable; automatic skill installation depends on the agent host. If the repository is available, use that checkout. Otherwise obtain the source or, on Windows, the complete `COFdraw-Windows.zip` from its Releases page. The EXE requires its `_internal` and `data` folders. Do not download a lone EXE or assume a browser-only agent can run it.

For an installable skill, select the **`skills/cofdraw` directory**, including its scripts and references. Installing this directory does not install the application. The helper accepts an explicit checkout or executable location and works when the skill is installed separately.

Read [local-api.md](references/local-api.md) when setting up a runtime, preparing custom monomers, editing a recipe, or changing presentation options. For 3D or calculations, also read [science-api.md](references/science-api.md). The helper [cofdraw_client.py](scripts/cofdraw_client.py) starts a temporary headless service, copies catalog/model data into a fresh output directory, saves actual returned exports, and stops its own service.

## Drawing workflow

1. Query `library` to discover real source recipe IDs, supported linkages, monomers, templates, and styles. Check `source_registry_errors`; do not invent catalog IDs or treat a reference-only linkage as implemented.
2. For a named catalog drawing, use `source` with its exact recipe ID and preferred style (`source` is the helper default). Preserve admitted coordinates, including a compact, tilted or irregular source arrangement. Do not automatically replace it with a regular polygon or run cleanup. For a new design, choose explicit finite monomer SMILES, supported linkage, and an appropriate template, then use `design` with the JSON prepare request. For an edited recipe or advanced style controls, use `render` with the full render request.
3. Inspect the returned `audit`, including continuation ports, source/design status, selected template, export diagnostics, and any style limits. Visually inspect the generated PNG/SVG when an image tool is available. Preserve isotope, charge, stereo, and attachment information when editing.
4. Return links to the files actually produced. Keep the JSON alongside chemical exports. Say when CDXML was omitted and use the reason in `chemical_exports`; never claim a successful format without a returned file.

From a checkout, with its Python environment active:

```shell
python skills/cofdraw/scripts/cofdraw_client.py --repo . --output outputs/library library
python skills/cofdraw/scripts/cofdraw_client.py --repo . --output outputs/tppa source tppa-pore-v1 --style standard --size 1200
python skills/cofdraw/scripts/cofdraw_client.py --repo . --output outputs/pyurea source pyurea-source-layout-v1 --style source
```

On Windows, the helper itself needs only Python's standard library when using the packaged application:

```powershell
python path\to\cofdraw_client.py --exe path\to\COFdraw.exe --output outputs\tppa source tppa-pore-v1
```

Use a new output directory for each call. The helper refuses overwrites. It does not open a desktop window or require WebView2 in headless mode.

## Chemistry and identity

- A source catalog entry is **source-bound provisional**, not independently verified chemistry. Preserve DOI, source recipe hashes, qualifiers, and unknown values.
- A custom or edited chemical recipe is a **designed variant**. Do not transfer the parent COFex identity, literature name, or experimental topology to it. Style-only rendering of an unchanged registered source can retain its qualified source identity.
- Template labels describe chosen connectivity. Monomer degree or a regular pore does not establish a unique topology, stacking, feasibility, crystal structure, or complete-material formula/SMILES. `cat0` is a reference convention, not independent evidence of noninterpenetration.
- Continuation atoms (`*`/`R`) and wavy external connection points are boundaries. Do not replace them with terminal carbons or hydrogens. Keep ChemDraw-edited boundary graphs in CDXML; ChemDraw MOL saves can remove external connection points.
- Layered, tilted, and cartoon styles are illustrations. Style changes must preserve the graph. Use Standard when another style cannot represent stereo marks.
- Exact source registration and source-informed reconstruction are different. Read `source_registered`, `coordinate_provenance` and `layout_evidence`. Pyurea/UCOF preserve manually source-informed revisions, not exact atomwise pixel registrations. COF-300 is a projection, herringbone is a rejected source proposal, and nickel-site drawings do not establish bulk metal occupancy.
- Unsupported chemistry should remain explicit. Report the limitation instead of fabricating a reaction transform or silently changing the user's monomers.

## Structures and calculations

Query `science library` for the actual periodic model IDs. Use `science monomer` for a complete finite molecule, `science model` for an explicit CIF or a separate built-in reticular design, and `science patch` for an illustrative stack of the active finite drawing. Their scope is different; retain the returned qualification and JSON with XYZ/SDF/CIF/PDB exports. Never attach the active drawing's source identity to an unrelated periodic template.

`science porosity` requires an explicit periodic structure. It estimates hard-sphere geometry, not adsorption or BET area. `science bet` requires the user's adsorption data and a chosen fit interval. Do not invent an isotherm to supply a requested material property. No bundled calculation runs DFT or xTB, establishes a band gap, or proves a stable crystal. Report units, parameters, convergence/sampling limits and failed checks from the actual result.

If this agent cannot execute local programs or access files, explain that limitation and provide the desktop download or exact runnable command. Do not claim that reading the repository generated a drawing.
