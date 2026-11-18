# Scientific API and scope

Use the same loopback service as drawing requests. The helper's `science ENDPOINT --request input.json` sends the supplied body, saves `result.json`, and downloads returned files as `structure.EXT`. `library`, `capabilities`, `synthesis-library` and `adsorption-engine` are read-only and need no request. `adsorption-engine` maps to `/api/drawing/science/adsorption/engine`; `adsorption` manages the asynchronous jobs below. Other endpoints map directly to `/api/drawing/science/ENDPOINT`. `artifacts.json` hashes actual files. A successful calculation is not experimental verification.

From a checkout with its environment active:

```shell
python skills/cofdraw/scripts/cofdraw_client.py --repo . --output outputs/models science library
python skills/cofdraw/scripts/cofdraw_client.py --repo . --output outputs/monomer science monomer --request monomer.json
```

## Choose the right object

| Endpoint | JSON body | Meaning and outputs |
| --- | --- | --- |
| `GET library` | None | Available independent periodic model templates and methods |
| `GET capabilities` | None | Running operation paths and actual JSON request schemas |
| `GET synthesis-library` | None | Qualified offline RDKit SA Score for each eligible finite catalog monomer, including reference-only entries |
| `GET adsorption/engine` | None | Actual local RASPA3 availability, version and profile; helper name `adsorption-engine` |
| `POST query` | `{"query":"Pyurea","record_types":["source"],"limit":10}` | Exact local catalog records with stable IDs, evidence and support status |
| `POST monomer` | `{"smiles":"Nc1ccc(N)cc1","conformers":2}` | Complete finite molecule; RDKit descriptors, ETKDG conformers, MMFF94s/UFF optimization, empirical charges, XYZ/SDF/JSON |
| `POST model` | `{"template_id":"hcb-imine-aa","repeats":[2,2,1]}` or `{"cif_text":"...","repeats":[1,1,1]}` | Explicit periodic structure and actual supercell; CIF/XYZ/PDB/JSON when representable |
| `POST patch` | `{"source_recipe_id":"ucof-source-layout-v1","layers":2,"layer_spacing":3.4,"stacking":"AB","slip":[1.7,0]}` | Finite layer illustration from admitted drawing coordinates; XYZ/JSON, continuation sites remain `X` |
| `POST finite-inspect` | `{"model":{...}}`, where `model` is the actual `patch` result's `viewer` | Displayed atom inventory, supplied bonds, finite components, extents and contact checks; no CIF needed |
| `POST finite-powder` | `{"model":{...},"wavelength":1.5406,"two_theta_min":2,"two_theta_max":50,"step":0.1}` | Finite-atom Debye scattering from that viewer model; CSV/JSON, no hkl assignment |
| `POST periodic-options` | `{"recipe":{...}}` or `{"source_recipe_id":"..."}` | Exact precursor/linkage/unchanged-template matches to separate bundled periodic designs, or an explicit no-match reason |
| `POST build-periodic` | `{"recipe":{...},"layer_spacing":null,"repeats":[1,1,1]}` | Matched or assembled unrelaxed AA cell, CIF and validation. `available:false` explains an unsupported drawing or fit |
| `POST stacking` | `{"structure":{...},"stacking":"AB","spacing":3.4,"slip":[1.5,0]}` | Periodic stacking hypothesis from one screened planar layer; CIF/XYZ/PDB/JSON and provenance |
| `POST networks` | `{"structure":{...},"offsets":[[0,0,0],[0.5,0.5,0.5]],"repeats":[1,1,1]}` | Explicit translated periodic-network hypothesis; CIF/XYZ/PDB/JSON, copy groups and collision diagnostics; does not assign interpenetration or COFex `cat` |
| `POST measure` | `{"source_recipe_id":"pyurea-source-layout-v1"}` | Conventionally scaled 2D opening estimate and sampling qualification; no structure export |
| `POST porosity` | `{"structure":{...},"probe_radius":1.82,"samples":12000}` | Hard-sphere geometric estimates for an explicit periodic cell; JSON |
| `POST geometry` | `{"structure":{...},"probe_radius":1.82,"grid_spacing":0.6,"surface_samples":64,"refinement":false}` | Deterministic periodic grid, PLD/LCD and transport-volume bounds, connected surface estimate and resolution status; JSON |
| `POST adsorption/jobs` | See below | Starts a local asynchronous RASPA3 job; polling returns actual JSON/CSV/input-output ZIP on termination |
| `POST bet` | `{"relative_pressures":[...],"uptake_mmol_g":[...],"fit_min":0.05,"fit_max":0.30,"cross_section_nm2":0.162}` | Fit of supplied adsorption data in a chosen interval; JSON and consistency checks |
| `POST inspect` | `{"structure":{...},"bond_scale":1.15,"tolerance_angstrom":0.15}` | Cell/composition/density summary and distance-inferred periodic contact network; JSON |
| `POST pxrd` | `{"structure":{...},"wavelength":1.5406,"two_theta_min":2,"two_theta_max":50,"step":0.02,"fwhm":0.15,"b_iso":0}` | Ideal calculated powder pattern and peak table; CSV/JSON |
| `POST screen` | `{"entries":[{"name":"Input A","structure":{...}}],"include_porosity":false,"probe_radius":1.82,"samples":4000}` | Bounded comparison of up to eight supplied periodic models; CSV/JSON |

Use IDs from the running `library`; the four bundled periodic designs are not crystal assignments to the active 2D source. `model` requires exactly one of `template_id` and `cif_text`. Repeats are three integers from 1 to 4 and are further limited by total atom count. `patch` and `measure` accept either a registered `source_recipe_id` or an editable `recipe`, never both. Layers are 1–6, spacing is 2.5–12 Å; stacking supports `AA`, `AB`, `ABC`, `slipped` and `twisted`. Offset modes require an explicit two-component slip in Å; twisted mode uses `twist_degrees` per layer (−180 to 180). These finite transforms do not assign crystallographic stacking. `monomer` accepts 1–8 requested conformers and rejects wildcard-containing incomplete molecules.

`finite-inspect` and `finite-powder` take the complete returned `viewer`, not the enclosing patch result or a periodic `structure`. It must declare `periodic: false`, with Cartesian atom fields `element`, `x`, `y`, `z` in Å and supplied `bonds`. Continuation `*`/`X`/`R` sites are omitted from atom physics and counted separately; missing hydrogens are not added. Inspection accepts up to 12,000 displayed sites, with pair-contact checks available only up to 2,000 real atoms. Scattering requires at most 2,000 real atoms and 300 million atom-pair/angle operations, so use fewer layers, a narrower range or a larger step if rejected. It returns `model_sha256`, the actual pair-sum profile and qualifications, not periodic Bragg peaks.

`periodic-options` accepts the drawing request used for rendering and returns `options`, `reason`, `scope` and `preserves_drawing_coordinates: false`. A match requires two exact isomeric precursor graphs, the supported product transform, and an unchanged compatible design template. Source IDs without that declared design signature do not silently acquire a periodic model. To select an offered option, call `model` with `template_id`, `repeats` and `recipe_match: {"recipe":...}` (or the same source-ID request). The backend rechecks the match before returning the hash-bound cell. `recipe_match` is not allowed with CIF imports. The new coordinates are a separate unrelaxed design hypothesis; source drawing identity is not an experimental validation.

`stacking` requires an explicit periodic input that can be screened as one near-planar layer parallel to cell vectors a,b: thickness at most 0.35 Å, a robust distance-inferred 2D covalent graph, no short contacts or cross-layer covalent/borderline contacts. Molecular cells, 3D networks and multi-layer cells are rejected. Use `AA`, `AB`, `ABC` or `slipped`, `spacing` 2.5–12 Å and a two-component `slip` in Å along the orthonormal in-plane frame defined by a and the layer normal. AA has one layer per cell, AB two and ABC three; progressive slip uses an oblique c vector. There is no periodic `twisted` mode. Compare each arrangement from the same original single layer, then call `pxrd` on each returned `structure` with identical parameters. These are unrelaxed hypotheses; inspect resulting contacts before interpreting them.

`networks` requires an explicit periodic structure and one to four distinct fractional offsets. It rejects coincident copies and severe atom collisions; inspect returned short-contact warnings and retain the JSON for copy grouping. The copy count is a user hypothesis, not proof of entanglement or a COFex `cat` assignment. Reuse the original structure when comparing offsets, rather than recursively copying an already combined candidate.

For `porosity`, take the returned model's explicit `structure` object or supply a full periodic cell. Cell vectors are rows in Å; atoms use `fractional: [x,y,z]` or Cartesian `position: [x,y,z]`, with full occupancies. Do not pass a finite patch illustration. Probe radius is in Å; the API accepts 0–3 and 1,000–30,000 requested volume samples, subject to computational limits. Check effective sample counts. CIF import requires explicit symmetry operations when needed and rejects unsupported disorder/partial occupancy rather than silently guessing.

`geometry` uses the same explicit structure. Probe radius is 0–5 Å, requested grid spacing 0.1–3 Å, surface samples 16–256 per atom, and `refinement` is a boolean. Inspect the effective grid/work limits and both diameter bounds. `certified_channel` means a continuous lower-grid path exists; `unresolved_at_this_resolution` means an upper enclosing graph can wind but no lower path was found; `no_channel_at_this_probe` means even that upper graph cannot wind. Zero lower accessible volume/surface in an unresolved result does not establish a closed pore. Surface estimates use finite quadrature and have no reported statistical confidence interval.

## Local adsorption jobs

The helper performs the complete job lifecycle:

```shell
python skills/cofdraw/scripts/cofdraw_client.py --repo . --output outputs/engine science adsorption-engine
python skills/cofdraw/scripts/cofdraw_client.py --repo . --output outputs/adsorption science adsorption --request adsorption.json
```

The request includes the actual `structure` returned by `model` or supplied from a validated full-cell CIF; do not substitute a finite drawing. A minimal settings example, with `structure` added from that actual result, is:

```json
{
  "structure": {},
  "gases": ["CO2", "N2"],
  "mole_fractions": [0.15, 0.85],
  "temperature_K": 298.15,
  "pressures_Pa": [10000, 100000],
  "cycles": 2000,
  "init_cycles": 1000,
  "widom_cycles": 2000,
  "calculate_heat": true,
  "seed": 42137,
  "charge_mode": "cif",
  "fugacity_coefficients": [1, 1]
}
```

The empty structure above is a placeholder and is rejected. Gas choices are `CO2`, `N2`, `CH4`, `H2`, `Ar`; select one to three distinct gases and positive fractions summing to one. Up to eight distinct positive pressures are accepted, in Pa, at 20–600 K. Each cycle count is 100–100,000, further bounded by model/work limits. These short defaults are workflow settings, not established convergence criteria. Use a fresh output directory for each comparison or seed.

`charge_mode: "cif"` requires a charge on every framework atom and a neutral cell. If charges are missing, either supply defensible external charges or explicitly choose `"zero"` as an approximation; never silently convert to zero or use finite-monomer Gasteiger values as validated framework charges. Fugacity coefficients default to one, an ideal-gas reservoir assumption. Supplied coefficients must correspond to the chosen conditions; no mixture equation of state is silently inferred. Preserve force-field/profile hashes and warnings. See `docs/raspa-methods.md` in the checkout for the pinned rigid-framework DREIDING/UFF and gas models, shifted 12 Å cutoff, Ewald treatment and limitations, including classical H2.

Direct clients use `POST /api/drawing/science/adsorption/jobs` (HTTP 202), then `GET /api/drawing/science/adsorption/jobs/{job_id}`. Active states are `queued`, `running`, `cancelling`; terminal states are `complete`, `failed`, `cancelled`, `interrupted`. Cancel only your own returned ID with `POST .../{job_id}/cancel`, then continue polling until terminal. One job runs per service. Each simulation stage has a one-hour limit; the native worker is isolated from the API process.

The helper saves `job-status.json`, retains the terminal result in `result.json`, downloads actual `structure.json`, `structure.csv`, `structure.zip` when returned, and hashes them in `artifacts.json`. The ZIP retains stage inputs, output and hashes; cancelled/failed jobs may have only partial data. On interruption the helper requests cancellation and awaits confirmation before stopping its service. If creation/cancellation cannot be confirmed, it leaves its local service running and reports `server.json` for inspection instead of killing an uncertain active worker. A failed or cancelled run is not a successful isotherm.

Henry coefficients come from a separate zero-loading Widom calculation, in mol kg⁻¹ Pa⁻¹. GCMC loadings are in mol kg⁻¹ (numerically mmol g⁻¹); mixture selectivity is `(q_i/q_j)/(y_i/y_j)`, while Henry selectivity is `K_H,i/K_H,j`. Missing/nonpositive component estimates yield no ratio, not infinity. Ratios have no propagated confidence interval.

When `calculate_heat` is true, two additional independent Widom runs at `T ± min(5 K, 0.025T)` estimate `Qst0 = R Δln(K_H,p)/Δ(1/T)` in kJ mol⁻¹. This is a finite-difference **zero-loading** estimate, with sampling noise and temperature-step bias; it is not a finite-loading heat curve. The envelope formed from the two input intervals is not a joint 95% confidence interval. Neither a low block error nor successful execution establishes convergence or agreement with experiment.

## Interpret the result

- Molecular descriptors and force-field conformers describe a finite molecular graph. Keep optimizer status, chosen force field, seed, requested/generated conformer counts and energies. Energies compare conformers of the same molecule under the same force field; they are not framework stability, electronic energies or heats of formation. Charges are empirical Gasteiger values, not quantum density.
- SA Score is the published offline RDKit finite-molecule heuristic, 1 easier–10 harder. It is not ASKCOS, a distilled model, a synthetic route or COF reaction/yield prediction. Drawing audits summarize distinct declared finite precursors only, with no aggregate if any is unavailable; no infinite framework is scored.
- Periodic model templates are unrelaxed reticular designs. CIF imports retain supplied geometry; export does not optimize or validate a crystal. CIF/PDB contain coordinates rather than the drawing's chemical bond-order annotations. Preserve JSON provenance and the source CIF.
- Finite stacks retain the drawing's arrangement after a conventional bond-length scale. They omit physical relaxation and do not infer a cell, interlayer bonds, missing hydrogens or stacking from a pore's appearance.
- `finite-inspect` reports only displayed real atoms, supplied real-atom bonds and finite graph components. Its atom inventory is not a complete COF formula. Cell density and periodic dimensionality are unavailable. `finite-powder` performs the independent-atom Debye pair sum without periodic images, thermal factors, polarization or instrumental broadening; it is not a crystal PXRD prediction. UI comparison curves are individually normalized to their own maxima, so they compare shapes and positions, not phase fractions or absolute scattering strength.
- `measure` samples a circle inside a projected macrocycle using conventional bond scaling and explicit-atom van der Waals disks. The result describes the opening in the drawing plane. Projected crossings or missing rings may be rejected. Use `geometry` on a periodic model for pore-limiting diameter.
- `porosity` reports volume/surface sampling intervals, atomic radii, listed-atom density and a separate voxel connectivity estimate. Local probe-exposed surface includes closed cavities; it is not BET area or certified transport-accessible area. Largest sampled included-sphere diameter is a lower bound, not pore-limiting diameter. The local clearance histogram is not a conventional pore-size distribution. Structural/radius/discretization errors are outside sampling intervals.
- `bet` needs actual supplied adsorption-branch data, positive uptake in mmol/g and increasing relative pressure between 0 and 1. At least three points must fall inside the chosen fit range. Preserve all consistency diagnostics; high R² alone does not establish a valid microporous-material BET area. Never substitute a geometric surface estimate or fabricated data.
- `inspect` reports composition and listed-atom crystal density plus distance-inferred contacts. Element inventories in mmol/g are not counts of accessible or catalytic sites. Network dimensionality comes from cell-translated contact cycles and can change with the contact cutoff; it does not assign chemical bond orders, formal valence, an RCSR net or interpenetration multiplicity.
- `pxrd` computes ideal intensities from the supplied cell/positions with neutral-atom IT92 form factors, Lorentz/polarization correction, a chosen common isotropic B value and Gaussian FWHM. Wavelength is in Å; angular values are degrees 2θ. No experimental phase identification/refinement, preferred orientation, anomalous dispersion, Kα doublet or instrumental profile is inferred. A chosen line width does not measure crystallite size or strain.
- `screen` compares the supplied models under recorded parameters and retains per-entry warnings/errors. A screening table is not a prediction of adsorption performance, selectivity, catalytic activity, stability or synthesizability.

No bundled endpoint runs DFT or xTB, predicts a band gap, or certifies an experimental crystal. For detailed method equations and bounded validation, read `docs/scientific-methods.md` in the repository. Return the actual files with their method, units, assumptions, warnings and unresolved uncertainties.

See [agent-workflows.md](agent-workflows.md) for precise local catalog queries, moving an actual model's `structure` into inspection/PXRD/screening requests, and rereading saved reports by `report_id`. Current schemas from `/science/capabilities` take precedence over examples if the application version changes.
