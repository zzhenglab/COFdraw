# COF materials workflows

COFdraw combines finite chemical drawings with explicit periodic-model inspection. A useful screening workflow starts by checking the model, then compares quantities calculated under the same assumptions. The structure, calculation and evidence remain queryable through the local API and agent skill.

## A quick check before using a CIF

Open **3D & properties**, import a CIF or generate a periodic template, and select **Quick check**. **Check this structure** reports the full-cell formula, density, cell parameters, elemental mass fractions and elemental atom inventories in mmol/g. Missing hydrogen, guests or counterions affect these quantities. An elemental N inventory is useful for comparing compositions; it is not a count of accessible catalytic or adsorption sites.

The same report inspects periodic contacts, isolated atoms, coordination and distance-inferred network dimensionality. A layer is 2D; a three-dimensional spanning network is 3D; a finite molecule crossing a cell boundary remains 0D. Warnings and cutoff sensitivity identify cases requiring chemical review. The graph does not assign a named net or an interpenetration multiplicity. See [periodic connectivity](periodic-network-inspection.md).

Turn on **Include quick pore screening** to add sampled pore opening, void volume and local probe-exposed surface estimates. These use the full periodic model and a common default probe; they retain the limitations in [scientific methods](scientific-methods.md). The local exposed area may include inaccessible cavities. It is never labeled BET.

## Compare candidate models

**Compare CIF files** accepts up to eight files. The comparison uses the same settings for every imported structure and exports CSV plus JSON. Review geometry errors and composition completeness before ranking models. The viewer keeps the currently loaded model; the table identifies each compared file separately.

Comparing AA/offset stacking models can reveal changes in geometric accessibility and calculated diffraction while leaving chemistry unchanged. A favorable pore metric does not establish an experimentally realizable stacking arrangement or energetic stability.

## Check expected diffraction

The **PXRD** tab calculates a kinematic powder pattern from the current cell and atoms, with an adjustable wavelength, 2-theta range and Gaussian width. The preview and strongest-peak table show positions, d spacings, relative integrated intensities and contributing Miller indices. Export CSV for experimental comparison or publication plotting.

This is useful for comparing candidate cells and stacking hypotheses. The model does not perform phase identification or refinement, and its width does not estimate crystallite size or strain. It omits preferred orientation, stacking disorder and several instrument effects. The exact formula, scattering coefficients, verification and bounds are in [characterization methods](characterization-methods.md).

## Inspect the structure and retrieve evidence

View along **a**, **b** or **c**, toggle hydrogen visibility, inspect the fractional atom table and measure distances between displayed atoms. Display-only hydrogen hiding does not remove hydrogen from calculations. The atom table shows the calculation's unit cell; the displayed model may be a repeated supercell.

The agent query API searches bundled monomer, linkage, topology, source, crystal and workflow records by normalized names, keywords or DOI. It reports available versus reference-only chemistry and retains source qualifications. Unknown queries return an empty result; the app does not invent an answer or silently query the internet. The desktop Find information panel has been removed.

The reports API retrieves local saved scientific reports. Each new report has an operation, engine version, timestamp, original request, input hash, method and units. Agents use these interfaces:

```text
GET  /api/drawing/science/capabilities
POST /api/drawing/science/query
POST /api/drawing/science/inspect
POST /api/drawing/science/pxrd
POST /api/drawing/science/screen
POST /api/drawing/science/networks
GET  /api/drawing/science/reports
GET  /api/drawing/science/reports/{report_id}
```

Read the returned request schemas instead of guessing arguments. Retain stable record IDs, source qualifiers and warnings when answering a user. The [agent workflow guide](../skills/cofdraw/references/agent-workflows.md) provides complete examples.

## Next scientific integrations

The most valuable further integrations are periodic geometry relaxation with validated parameters, a dedicated pore-accessibility backend, and adsorption simulation with explicit force fields and conditions. They would add physical validation beyond the current geometric screening. [Zeo++](https://www.zeoplusplus.org/) provides established pore-network and accessibility analysis; it is not bundled in this release. The existing monomer force fields do not relax a periodic COF, and no gas uptake, selectivity or electronic band gap is inferred from a drawing.

The design direction draws on structure inspection, cell measurements and diffraction workflows documented by [CrystalMaker](https://crystalmaker.com/crystalmaker/specs/index.html), with COF-specific evidence and calculation scopes retained. COFdraw is independently developed and has no affiliation with CrystalMaker.
