# Scientific API and scope

Use the same loopback service as drawing requests. The helper's `science ENDPOINT --request input.json` sends the file verbatim to `/api/drawing/science/ENDPOINT`, saves `result.json`, and downloads returned files as `structure.EXT`. `science library` needs no request. `artifacts.json` hashes actual files. A successful calculation is not experimental verification.

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
| `POST query` | `{"query":"Pyurea","record_types":["source"],"limit":10}` | Exact local catalog records with stable IDs, evidence and support status |
| `POST monomer` | `{"smiles":"Nc1ccc(N)cc1","conformers":2}` | Complete finite molecule; RDKit descriptors, ETKDG conformers, MMFF94s/UFF optimization, empirical charges, XYZ/SDF/JSON |
| `POST model` | `{"template_id":"hcb-imine-aa","repeats":[2,2,1]}` or `{"cif_text":"...","repeats":[1,1,1]}` | Explicit periodic structure and actual supercell; CIF/XYZ/PDB/JSON when representable |
| `POST patch` | `{"source_recipe_id":"ucof-source-layout-v1","layers":2,"layer_spacing":3.4,"stacking":"AB","slip":[1.7,0]}` | Finite layer illustration from admitted drawing coordinates; XYZ/JSON, continuation sites remain `X` |
| `POST measure` | `{"source_recipe_id":"pyurea-source-layout-v1"}` | Conventionally scaled 2D opening estimate and sampling qualification; no structure export |
| `POST porosity` | `{"structure":{...},"probe_radius":1.82,"samples":12000}` | Hard-sphere geometric estimates for an explicit periodic cell; JSON |
| `POST bet` | `{"relative_pressures":[...],"uptake_mmol_g":[...],"fit_min":0.05,"fit_max":0.30,"cross_section_nm2":0.162}` | Fit of supplied adsorption data in a chosen interval; JSON and consistency checks |
| `POST inspect` | `{"structure":{...},"bond_scale":1.15,"tolerance_angstrom":0.15}` | Cell/composition/density summary and distance-inferred periodic contact network; JSON |
| `POST pxrd` | `{"structure":{...},"wavelength":1.5406,"two_theta_min":2,"two_theta_max":50,"step":0.02,"fwhm":0.15,"b_iso":0}` | Ideal calculated powder pattern and peak table; CSV/JSON |
| `POST screen` | `{"entries":[{"name":"Input A","structure":{...}}],"include_porosity":false,"probe_radius":1.82,"samples":4000}` | Bounded comparison of up to eight supplied periodic models; CSV/JSON |

Use IDs from the running `library`; the four bundled periodic designs are not crystal assignments to the active 2D source. `model` requires exactly one of `template_id` and `cif_text`. Repeats are three integers from 1 to 4 and are further limited by total atom count. `patch` and `measure` accept either a registered `source_recipe_id` or an editable `recipe`, never both. Layers are 1–6, spacing is 2.5–12 Å, stacking is `AA` or `AB`; AB requires an explicit two-component slip in Å. `monomer` accepts 1–8 requested conformers and rejects wildcard-containing incomplete molecules.

For `porosity`, take the returned model's explicit `structure` object or supply a full periodic cell. Cell vectors are rows in Å; atoms use `fractional: [x,y,z]` or Cartesian `position: [x,y,z]`, with full occupancies. Do not pass a finite patch illustration. Probe radius is in Å; the API accepts 0–3 and 1,000–30,000 requested volume samples, subject to computational limits. Check effective sample counts. CIF import requires explicit symmetry operations when needed and rejects unsupported disorder/partial occupancy rather than silently guessing.

## Interpret the result

- Molecular descriptors and force-field conformers describe a finite molecular graph. Keep optimizer status, chosen force field, seed, requested/generated conformer counts and energies. Energies compare conformers of the same molecule under the same force field; they are not framework stability, electronic energies or heats of formation. Charges are empirical Gasteiger values, not quantum density.
- Periodic model templates are unrelaxed reticular designs. CIF imports retain supplied geometry; export does not optimize or validate a crystal. CIF/PDB contain coordinates rather than the drawing's chemical bond-order annotations. Preserve JSON provenance and the source CIF.
- Finite stacks retain the drawing's arrangement after a conventional bond-length scale. They omit physical relaxation and do not infer a cell, interlayer bonds, missing hydrogens or stacking from a pore's appearance.
- `measure` samples a circle inside a projected macrocycle using conventional bond scaling and explicit-atom van der Waals disks. Implicit H, stacking and out-of-plane shape are omitted. It is a 2D drawing estimate, not experimental pore size or a periodic pore-limiting diameter. Projected crossings or missing rings may be rejected.
- `porosity` reports volume/surface sampling intervals, atomic radii, listed-atom density and a separate voxel connectivity estimate. Local probe-exposed surface includes closed cavities; it is not BET area or certified transport-accessible area. Largest sampled included-sphere diameter is a lower bound, not pore-limiting diameter. The local clearance histogram is not a conventional pore-size distribution. Structural/radius/discretization errors are outside sampling intervals.
- `bet` needs actual supplied adsorption-branch data, positive uptake in mmol/g and increasing relative pressure between 0 and 1. At least three points must fall inside the chosen fit range. Preserve all consistency diagnostics; high R² alone does not establish a valid microporous-material BET area. Never substitute a geometric surface estimate or fabricated data.
- `inspect` reports composition and listed-atom crystal density plus distance-inferred contacts. Element inventories in mmol/g are not counts of accessible or catalytic sites. Network dimensionality comes from cell-translated contact cycles and can change with the contact cutoff; it does not assign chemical bond orders, formal valence, an RCSR net or interpenetration multiplicity.
- `pxrd` computes ideal intensities from the supplied cell/positions with neutral-atom IT92 form factors, Lorentz/polarization correction, a chosen common isotropic B value and Gaussian FWHM. Wavelength is in Å; angular values are degrees 2θ. No experimental phase identification/refinement, preferred orientation, anomalous dispersion, Kα doublet or instrumental profile is inferred. A chosen line width does not measure crystallite size or strain.
- `screen` compares the supplied models under recorded parameters and retains per-entry warnings/errors. A screening table is not a prediction of adsorption performance, selectivity, catalytic activity, stability or synthesizability.

No bundled endpoint runs DFT or xTB, predicts a band gap, or certifies an experimental crystal. For detailed method equations and bounded validation, read `docs/scientific-methods.md` in the repository. Return the actual files with their method, units, assumptions, warnings and unresolved uncertainties.

See [agent-workflows.md](agent-workflows.md) for precise local catalog queries, moving an actual model's `structure` into inspection/PXRD/screening requests, and rereading saved reports by `report_id`. Current schemas from `/science/capabilities` take precedence over examples if the application version changes.
