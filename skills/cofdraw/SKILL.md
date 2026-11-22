---
name: cofdraw
description: Use COFdraw or its GitHub link to query local chemistry, draw and export finite COF patches, preserve source layouts, inspect periodic models, screen monomer accessibility, calculate geometric properties or run qualified local adsorption jobs, and return actual reports.
---

# COFdraw

Create actual PNG/SVG and editable MOL/SDF/CDXML drawings with the local COFdraw application. Preserve the accompanying JSON and its provenance. Its separate scientific tools query local evidence and handle molecular conformers, periodic models, contact networks, ideal diffraction, geometric screening, classical adsorption jobs and supplied adsorption data. A finite drawing does not establish a crystal structure.

## Start from a GitHub link

The project is [zzhenglab/COFdraw](https://github.com/zzhenglab/COFdraw). A pasted link makes these instructions readable; automatic skill installation depends on the agent host. If the repository is available, use that checkout. Otherwise obtain the source or, on Windows, the complete `COFdraw-Windows.zip` from its Releases page. The EXE requires adjacent `COFdraw.exe.config`, `_internal`, and `data`. Do not download a lone EXE or assume a browser-only agent can run it. For downloaded Windows startup failures, read [startup troubleshooting](../../docs/windows-startup.md); do not change machine-wide policy or unblock unrelated folders.

For an installable skill, select the **`skills/cofdraw` directory**, including its scripts and references. Installing this directory does not install the application. The helper accepts an explicit checkout or executable location and works when the skill is installed separately.

Read [local-api.md](references/local-api.md) when setting up a runtime, preparing custom monomers, editing a recipe, or changing presentation options. For 3D or calculations, also read [science-api.md](references/science-api.md). For information retrieval, screening and saved reports, read [agent-workflows.md](references/agent-workflows.md). The helper [cofdraw_client.py](scripts/cofdraw_client.py) starts a temporary headless service, copies catalog/model data into a fresh output directory, saves actual returned exports, and stops its own service.

## Drawing workflow

For pasted precursor notation, POST `{"text":"… COFex.unk.cat0"}` to `/api/drawing/cofex/resolve`. Read the returned choices before preparing. An unknown net leaves the choice open. For Por-CHO + TAPA, compare the proposed planar `mtf_4_3` sheet with the spatial `tbo_4_3` and `fjh_4_3` projections. A legacy six-unit ring can also be restored. On a monomer change resolve the new pair again. Read the active recipe’s offered extents before requesting a repeat. A description plan can return `drawing_extent`. Use its `repeats` (or numeric `value`) and `extent_mode` (default `pores`) at `/api/drawing/repeat`. Return `drawing_cofex.compact` for exchange and preserve `drawing_cofex.text` when the exact linkage or drawing recipe matters.

For a prose request, use the local description workflow in [agent-workflows.md](references/agent-workflows.md) to obtain a reviewable, catalog-bound plan. Rules mode is the default. Use an optional provider API only when the user has explicitly chosen that mode and supplied its settings; never discover or reuse ambient credentials. An LLM proposal selects IDs and still passes the local chemistry checks.

The session conversation also supports constrained follow-ups through `/describe/chat`. Carry only its returned allowlisted context. Review each current plan before Create/Apply; a `style_only` intent must match the caller's current drawing revision and be applied to that local recipe without reloading or replacing its custom monomers. Do not send saved recipes, SMILES, CIFs, project files or history to a provider. See [drawing conversations](../../docs/conversation-methods.md).

1. Query `library` to discover real source recipe IDs, supported linkages, monomers, templates, and styles. Check `source_registry_errors`; do not invent catalog IDs or treat a reference-only linkage as implemented.
2. For a named catalog drawing, use `source` with its exact recipe ID and preferred style (`source` is the helper default). Preserve admitted coordinates, including a compact, tilted or irregular source arrangement. Do not automatically replace it with a regular polygon or run cleanup. For a new design, choose explicit finite monomer SMILES, supported linkage, and an appropriate template, then use `design` with the JSON prepare request. For an edited recipe or advanced style controls, use `render` with the full render request.
3. Inspect the returned `audit`, including continuation ports, source/design status, selected template, export diagnostics, and any style limits. Visually inspect the generated PNG/SVG when an image tool is available. Preserve isotope, charge, stereo, and attachment information when editing.
4. Return links to the files actually produced. Keep the JSON alongside chemical exports. Say when CDXML was omitted and use the reason in `chemical_exports`; never claim a successful format without a returned file.

For new designs, check the A/B pair together: an individually supported monomer may have no supported reaction with its partner. Inspect `audit.layout_refinement.quality` and `layout_diagnostics`; `needs_review` means the fit needs inspection or adjustment. The angle-aware fit preserves rigid molecular blocks, but is not an energy minimization or E/Z prediction. Keep source-coordinate revisions fixed rather than rerunning a generic fit. See [geometry review](references/local-api.md#geometry-review) for the relevant fields and source-overlay procedure.

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
- The scoped `nitrone` rule accepts aromatic aldehydes and neutral directly aryl-bound `Ar–NH–OH`. It preserves charged N–O in `Ar–CH=[N+]([O−])–Ar`; ordinary amines are not interchangeable. New E/Z, COF-115 identity, kgm topology and the paper's photo-rearrangement are not inferred. Read [nitrone limits](../../docs/nitrone-methods.md).

## Structures and calculations

Start with `science capabilities` for the running API's request schemas. Use `science query` to find exact local names, IDs, DOIs, monomers, reference-only linkages, topologies, models and workflows. Keep retrieved records and qualifiers; no match means absent from this bundle, not absent from chemistry. Distinguish catalog information from a calculation result.

Query `science library` for the actual periodic model IDs. Use `science monomer` for a complete finite molecule, `science model` for an explicit CIF or a separate built-in reticular design, and `science patch` for an illustrative stack of the active finite drawing. Their scope is different; retain the returned qualification and JSON with XYZ/SDF/CIF/PDB exports. Never attach the active drawing's source identity to an unrelated periodic template.

For the actual finite drawing, pass `patch`'s returned `viewer` to `science finite-inspect` or `science finite-powder` as `model`. These work without a CIF and report displayed-atom checks or finite-cluster Debye scattering. They do not supply cell density, periodic pore geometry or hkl peaks. Never substitute an unrelated periodic template merely to obtain a calculation.

`science build-periodic` uses a current drawing's declared precursor graphs, implemented linkage and unchanged connectivity to return a qualified unrelaxed CIF. Automatic spacing uses a bundled AA match if available. For supported hexagonal and square designs it can also build a new cell and verify the intended periodic bonds. Retain the model's assumptions and validation report. If it returns `available: false`, continue with finite checks or use a corresponding supplied CIF. `science periodic-options` and `model.recipe_match` let agents select a specific bundled match.

Use `science inspect` for cell composition/density and distance-inferred periodic contacts, `science pxrd` for an ideal periodic powder pattern, and `science screen` for a bounded comparison of supplied structures. Contact cutoff changes can change inferred dimensionality; it is not automatic RCSR identification. PXRD is a calculation from the supplied structure, not phase identification or refinement of an experimental pattern. Element counts per gram are not accessible catalytic-site densities. `science stacking` builds AA/AB/ABC/slipped periodic hypotheses only from a screened single planar 2D layer; arbitrary twist is available only as a finite illustration unless a commensurate periodic cell is supplied separately.

Use `science synthesis-library` for offline RDKit SA Scores of eligible finite monomers (1 easier–10 harder). Drawing audits summarize declared precursors only. These are heuristic scores, not ASKCOS, retrosynthesis, COF reaction feasibility or yield.

`science geometry` requires an explicit periodic structure and reports hard-sphere PLD/LCD and volume bounds plus a sampled transport-connected surface. Read `connectivity.status`: a coarse-grid zero can mean unresolved rather than closed. The older `science porosity` is a separate sampling method. Neither calculates adsorption or BET area.

For physical adsorption, first check `science adsorption-engine`, then read the asynchronous workflow in [science-api.md](references/science-api.md#local-adsorption-jobs). Use `science adsorption --request ...` only with an explicit periodic model and reviewed gas, charge, force-field and reservoir assumptions. Keep the service alive until its job finishes; the helper polls and cancels its own job on interruption. Retain JSON/CSV and the actual input/output ZIP. Short runs are not convergence evidence; zero-loading heat is a finite-difference Widom estimate, not a finite-loading heat curve.

`science bet` requires the user's supplied adsorption data and a chosen fit interval. Do not invent an isotherm or relabel simulated/geometric results as an experiment. No bundled calculation runs DFT or xTB, establishes a band gap, or proves a stable crystal. Report units, parameters, convergence/sampling limits and failed checks from the actual result.

Return the saved JSON report and any actual CSV/structure/image files, plus concise findings grounded in their values. Preserve `report_id`, `operation`, engine version, input/structure hashes, units, method and limitations where returned. A persistent local service can list and reread its previous reports; the helper's isolated runs also keep their complete result files on disk.

If this agent cannot execute local programs or access files, explain that limitation and provide the desktop download or exact runnable command. Do not claim that reading the repository generated a drawing.
