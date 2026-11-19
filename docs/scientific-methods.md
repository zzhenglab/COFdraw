# Molecular and periodic structure calculations

The 3D rendering panel displays an interpolated probe-center surface from the periodic clearance field used for geometric porosity. The channel view keeps samples connected to a certified periodic winding path. Rendering samples at most 22 centers per axis and repeats the surface over the displayed supercell. The section slider clips the full supercell along fractional c and fills the opening in the section plane. The exported geometry report retains the full calculation grid of the unit cell.

For tetrahedral stereochemistry, the planar-net builder embeds capped fragments with ETKDGv3, fixes junction positions during MMFF relaxation and checks the specified atom tags against the resulting coordinates. Placement uses proper 3D rotations. The assembled structure must pass the same periodic contact and intended-bond checks as other generated cells. Flexible substituents and the predicted layer spacing describe one starting conformation.

The external-analysis controls export the active periodic coordinates as a CIF for [Crystal Builder](https://github.com/JulesOpp/Crystal-Builder) or as an executable notebook for [Mofography](https://mofography.github.io/mofography/quick_start.html). The notebook runs the upstream segmentation in an environment with its Rust backend installed. Its radii, pocket definition and resulting descriptors follow Mofography.

COFdraw separates molecular conformers, explicitly periodic structures, finite drawing illustrations, hard-sphere geometry, classical adsorption simulations, and adsorption-data fitting. The result records the method, parameters, units, and limitations. A finite drawing never becomes a validated crystal merely by adding a cell around it.

## Monomer properties and conformers

`monomer_analysis(smiles, conformers=4, seed=42, max_iterations=300)` accepts one complete connected molecular graph, with at most 250 input atoms. Wildcards and continuation points are rejected. Molecular formula, molecular weight, exact mass, formal charge, ring counts, rotatable bonds, hydrogen-bond counts, topological polar surface area, Crippen logP, and fraction of sp3 carbons are calculated with RDKit. These are molecular descriptors, not predictions of a COF's electronic or adsorption properties. See the [RDKit descriptor reference](https://www.rdkit.org/docs/source/rdkit.Chem.rdMolDescriptors.html).

Explicit hydrogens are added before seeded, single-threaded ETKDGv3 conformer generation. Available MMFF94s parameters are preferred; UFF is the fallback. Each conformer's optimizer status and energy is retained. Selection prefers the lowest-energy converged conformer; if none converged, the result explicitly reports that limitation. Requested conformer count, pruning, seed, iteration limit, force field, and RDKit version are recorded. The conformer search is bounded and does not establish a global minimum. [RDKit conformer guide](https://www.rdkit.org/docs/GettingStartedInPython.html), [force-field status definitions](https://www.rdkit.org/docs/source/rdkit.Chem.rdForceFieldHelpers.html).

XYZ and SDF contain actual three-dimensional Cartesian coordinates in ångströms. SDF retains the molecular graph, stereochemistry, charge, and isotope information; XYZ is a coordinate format. Force-field energies are in kcal/mol and should only compare conformers of the same molecular graph under the same force field. They are not heats of formation, binding energies, orbital energies, or framework stability. The bundled calculation does **not** run xTB or DFT. XYZ is suitable as an initial geometry for separately installed quantum-chemistry programs after inspection.

The monomer response also provides empirical Gasteiger partial charges for all explicit atoms, including hydrogens, using 12 iterations. Missing parameters cause explicit omission instead of zero-valued stand-ins. These graph-based charges are not quantum electron density, oxidation states, HOMO/LUMO energies, or a band gap. [RDKit partial-charge implementation](https://www.rdkit.org/docs/source/rdkit.Chem.rdPartialCharges.html).

## Monomer solubility screening

**Screen solubility** takes a library monomer ID or a complete molecular SMILES. It runs locally without conformer generation. The screen compares 46 neat solvents near 25 °C. Exact canonical graph matches return the closest equilibrium observation within 5 K, with other publications or solid forms available in the result. Each observation carries the measured mole fraction, temperature, method and source DOI. BigSolDB converts mole fraction `x` to moles per liter of pure solvent using `x/(1 − x)` times the molarity of the pure solvent. COFdraw multiplies that concentration by the molecular weight to get mg/mL of solvent, then by the selected starting solvent volume to get mg. [BigSolDB 2.2](https://doi.org/10.5281/zenodo.22648301)

For queries without a matching observation, the estimator fits a small ridge regression for each solvent from the bundled equilibrium data. It uses molecular weight, Crippen logP, TPSA, hydrogen-bond counts, rotatable bonds, aromatic rings and fraction sp³. It adjusts the descriptor trend by half the weighted mean residual of up to eight measured analogs, using Morgan fingerprint similarity and scaled descriptor distance. At least three sufficiently close analogs are required. When this route has insufficient support, a second ridge model shares descriptor slopes across solvents and fits solvent-specific deviations and intercepts. It requires three training solutes within 0.95 scaled descriptor RMS distance in the requested solvent. Each test scaffold is withheld across all solvents. Only residual cases missed by the first method are used to benchmark the fallback. At least one method qualifies in 34 solvents. An unavailable query can still show measured concentrations and sources for up to three related compounds, each identified by its own SMILES. It withholds estimates for salts, metals and other graphs outside the neutral-organic prediction domain.

The bundled `data/solubility/validation.json` records five-fold Bemis–Murcko scaffold holdout with comparisons against a same-solvent median and a separate multi-aromatic subgroup. The displayed bands are based on the 90th-percentile held-out log errors. These can span orders of magnitude because the screen lacks a query-specific crystal form or melting point. The benchmark has limited coverage of large COF precursors. Mixed solvents, pH adjustments and elevated temperatures are outside the screen.

When a paper explicitly reports dissolving a catalog monomer in a neat solvent, COFdraw can list the actual mass and solvent volume as a separate preparation load. The [PyTTA and TA dioxane preparations](https://doi.org/10.1038/s41467-022-29086-x) and the [TAPP preparation in DMF](https://doi.org/10.3390/gels8040202) record the monomer dissolution step before mixing reagents. The DHP and TAM dioxane loads preserve the publication's quantities and unspecified preparation temperature. They are excluded from the equilibrium estimator and are saved separately in JSON and CSV. [Tavani et al., SI §S1.1](https://doi.org/10.1021/jacs.6c17514)

## Finite-precursor SA Score

The offline RDKit Contrib SA Score reports a molecular heuristic from 1 (easier) to 10 (harder) for eligible complete finite monomers. It is the published scorer, not ASKCOS or an ASKCOS-distilled model. Drawing audits give an unweighted mean/maximum of distinct declared precursor graphs only; an unsupported or missing precursor prevents a numeric aggregate. They do not score a product patch, polymerization, crystallization or infinite COF. Reference-only catalog entries retain their reaction-support restrictions. See [synthetic-accessibility inputs, resources and limitations](synthetic-accessibility.md).

## Explicit periodic structure import and export

`import_cif(text)` reads one conventional core-CIF data block, cell lengths/angles, fully occupied fractional atom sites, and explicit affine symmetry operations. Symmetry is expanded into a complete unit cell; identical symmetry images of special-position sites are deduplicated. The supplied operations must preserve the cell metric, include identity, and be closed under composition. A non-P1 space-group name alone is insufficient: the file must contain its operations, or be exported as fully expanded P1. Missing symmetry, partial occupancies, unresolved mixed sites, unsupported elements, invalid cells, and more than 2,000 expanded atoms are rejected with a reason. This intentionally limited parser is not a replacement for general crystallographic software. [IUCr core-CIF definitions](https://www.iucr.org/__data/iucr/cif/dictionaries/cif_core_2.4.2.dic.pdf).

`normalize_structure(structure)` also accepts explicit full-cell data:

```json
{
  "periodic": true,
  "cell": [[10, 0, 0], [0, 10, 0], [0, 0, 10]],
  "atoms": [{"element": "C", "fractional": [0.5, 0.5, 0.5]}],
  "provenance": {"kind": "user_supplied_geometry"}
}
```

This one-atom example is a geometry test fixture, not a COF. Cell vectors are rows, in ångströms, and `fractional @ cell` gives Cartesian coordinates. Cartesian atom input can instead use a `position: [x, y, z]` field. Cell orientation must be right-handed. Occupancy is exactly one; no chemical identity, bond order, or hydrogen placement is inferred from coordinates.

`structure_exports(structure, repeats=(1,1,1))` returns viewer atoms, cell, bounds, and CIF, extended XYZ, and PDB text. Integer repetition produces a real supercell containing the expected number of atoms and scaled lattice vectors; exports are limited to 12,000 sites. CIF is written as full-cell P1. Coordinates are rotated consistently into conventional cell orientation before PDB/CIF export. PDB is omitted if its fixed-width coordinate fields cannot represent the result. PDB does not carry framework bond orders, and ordinary XYZ readers may discard the extended lattice header. Neither an import nor an export performs crystal relaxation or establishes experimental provenance.

The bundled periodic examples were generated with **pyCOFBuilder 0.0.9.4**: HCB imine AA (66 atoms), its AB1 model (132), an extended HCB imine AA design (126), and a DIA imine network (146), including explicit hydrogens. Their manifest preserves generator identifiers, hashes, formula counts, and the metadata workaround used during generation. They are separate, unrelaxed reticular designs with chosen connectivity/stacking; selecting one does not transfer the active drawing's identity or turn that drawing into the same crystal. The generator is a development tool, not a runtime dependency. [pyCOFBuilder source](https://github.com/lipelopesoliveira/pyCOFBuilder), [method publication](https://doi.org/10.1021/acs.jcim.3c01918).

### Matching a drawing to a separate periodic design

`/science/periodic-options` compares an unchanged drawing-template signature, implemented linkage, and the two exact canonical isomeric precursor graphs with explicit match records in the bundled model manifest. It also reconstructs the expected product fragments and rejects edited chemistry that no longer follows those declarations. A match returns a bounded choice of separately generated models; it does not fit the finite drawing's coordinates to a cell. An empty result includes a reason and leaves finite-model checks available.

Selecting a returned option through `/science/model` with `recipe_match` revalidates that binding and the model file's manifest hash. The CIF and its provenance record the unrelaxed coordinates and the precursor match.

The viewer offers an ab-plane pore face, an inclined layer profile and an oblique mixed display. Color keys identify elements and labelled sites. Newly constructed cells draw bonds from the checked periodic image translations. Face and profile views are useful for comparing in-plane pore outlines with stacking, as illustrated in the Py-1P crystal structure. [Py-1P, Fig. 4h](https://pmc.ncbi.nlm.nih.gov/articles/PMC8927472/figure/Fig4/)

### Building a periodic cell from the current drawing

`/science/build-periodic` accepts the declared recipe. It verifies the same precursor chemistry and template incidence. With automatic spacing it first tries an exact bundled AA match. Otherwise it assembles the reacted building blocks on the chosen HCB 3+3, HCB 3+2, SQL 4+2 or SQL 4+4 connection graph. It uses each reactive port once across explicit cell translations. Building blocks are kept rigid while their orientations, positions and the cell vectors are fitted to the port directions and connection lengths.

The initial coordinates use a 1.40 Å drawing bond scale. Inter-block bonds target 1.29 Å for double bonds and 1.42 Å for single bonds. Hydrogens are added to the reacted blocks before their connection placeholders are removed. Supported neutral fragments can contain H, B, C, N, O, F, Si, P, S, Cl, Br and I. The builder requires converged restrained fragment geometry for the heavier elements. Tabulated isotope labels survive cell generation, replication, and mass-normalized characterization. CIF uses conventional element symbols plus the local site field `_atom_site_cofdraw_mass_number`. Extended XYZ uses an integer isotope column with zero for natural abundance. Other CIF readers may ignore the local field. Specified stereochemistry, radicals and unsupported topologies require supplied coordinates. Para-substituted meso phenyl rings on supported porphyrins are tilted 60° to reduce steric overlap. This is a chosen starting geometry inspired by observed aryl tilt in porphyrins. [Porphyrin diffraction study](https://journals.iucr.org/e/issues/2014/10/00/lx2293/)

The builder checks the periodic atomic edges against the intended graph, tests short contacts and verifies a connected rank-two layer. Automatic spacing uses 3.4 Å for the planar construction and starts at 5.0 Å for tilted phenyl arms. After a rejected fit, the builder tries additional unfused para-phenylene torsions and then restrained MMFF94 minimization of individually capped fragments. Cap and attachment positions stay fixed and cap atoms are excluded from the cell. The recovered layer uses its explicit thickness plus 3.4 angstrom for default spacing. Each candidate must pass the same periodic edge, contact and dimensionality checks. An explicit spacing is checked as supplied. The CIF is P1, with the unrelaxed model and its spacing identified in the header. The JSON stores the cell, precursor record and validation, including intended bonds and any formal charges for newly assembled models. Partial charges and experimental crystal identity require their own evidence.

Drawings without a supported declared template or a valid fit return `available: false` and a reason. The finite preview remains available. Phase assignment, stability and geometry relaxation require separate work. [Current-model workflow for agents](../skills/cofdraw/references/agent-workflows.md#a-periodic-structure-inspection-workflow)

## Finite layer illustrations

`patch_3d(recipe, layers=3, layer_spacing=3.4, stacking="AA", slip=None, bond_length=1.4, twist_degrees=15)` makes a separate **nonperiodic illustration** from an existing finite drawing. Its median real-atom bond length is scaled to the requested ångström value. AA aligns layers; AB uses alternating offsets `0, s, 0, s`; ABC uses `0, s, 2s, 0`; progressive slip uses `0, s, 2s, 3s`; twisted layers rotate successively around the drawing centroid by the chosen angle. Offset modes require an explicit slip vector in ångströms; pore shape does not infer it. These are chosen illustration transforms, not crystallographic stacking assignments such as AB1 or ABC1. Layer translations/rotations are recorded individually. Bond connectivity and internal distances are retained within each finite layer. No interlayer bonds, missing hydrogens, or optimized atomic geometry are invented. XYZ uses `X` for continuation placeholders. This output is not accepted as a periodic pore-analysis structure.

### Checking the actual displayed finite atoms

`/science/finite-inspect` accepts the patch result's `viewer` as `model`, with `periodic: false`, Cartesian coordinates in Å and supplied bonds. It reports real-atom inventory, extents, supplied real-atom bonds, finite graph components and bond-length range. `*`, `X` and `R` continuation placeholders are excluded from the physical inventory and counted separately. Missing hydrogens remain missing, so this inventory is not a complete molecular or COF formula. No cell density or periodic dimensionality is reported.

Up to 12,000 displayed sites can be inspected; a pair-contact check is available for at most 2,000 real atoms. It reports minimum pair distance and contacts shorter than 0.65 times the sum of RDKit covalent radii, including supplied bonds; only the first 50 are listed, while the total count remains available. This catches obvious geometric collisions, not chemical valence, stability or synthetic feasibility. Each result is bound to the submitted model's SHA-256. The UI invalidates old checks when the current drawing or generated model changes.

### Finite-cluster powder scattering

`/science/finite-powder` takes the same finite viewer model and evaluates the independent-atom Debye equation directly over every real-atom pair:

`I(Q) = Σᵢ fᵢ(Q)² + 2 Σᵢ<ⱼ fᵢ(Q) fⱼ(Q) sin(Qrᵢⱼ)/(Qrᵢⱼ)`, with `Q = 4π sin(θ)/λ`.

The zero-distance limit is one. This is an orientational average of the supplied finite object and does not require periodicity; it contains no reciprocal-lattice indexing or hkl assignments. [IUCr discussion of the Debye equation](https://journals.iucr.org/a/issues/2016/06/00/me0628/).

COFdraw uses its pinned neutral-atom IT92 form-factor table, coordinates in Å, wavelength in Å, and degrees 2θ. The response includes Q in Å⁻¹, intensity in electron-squared units, a profile normalized to a maximum of 100, CSV, model hash and coefficient hash. Default wavelength is 1.5406 Å, range 2–50°, step 0.1°. The work limit is 2,000 real atoms, 15,001 angular points and 300 million pair/angle products. A rejected calculation needs fewer layers, a narrower range or a larger step rather than a hidden approximation.

No periodic images, missing hydrogens, thermal motion, polarization, absorption, anomalous/incoherent scattering or instrumental profile are added. Continuation placeholders scatter no intensity. An illustrative drawing's finite scattering is therefore a comparison of its supplied geometry, not a periodic crystal PXRD prediction or experimental fit. The five finite stacking comparisons use the same drawing and layer count. Each plotted curve is independently normalized; peak heights do not measure phase fractions or relative total scattering power.

## Explicit periodic stacking comparisons

`/science/stacking` constructs AA, AB, ABC or progressively slipped periodic hypotheses from a supplied single near-planar layer parallel to a,b. The unwrapped layer thickness must be at most 0.35 Å. A distance-inferred covalent graph must have robust rank-two components, no flagged short contacts and no covalent or borderline contacts between layers along c. Boundary-wrapped atoms are lifted into the same physical layer before these checks. Molecular, chain, 3D-network, multilayer or cutoff-sensitive inputs are rejected rather than flattened into a different structure.

The two slip components are Å along the in-plane orthonormal frame formed by the a vector and layer normal. Internal Cartesian geometry is preserved, including the original oblique cell's in-plane c contribution. AA has one layer per cell; AB uses `0,s` with a two-layer repeat; ABC uses `0,s,2s` with a three-layer repeat; progressive slip uses one layer and an oblique c vector. Spacing is explicit, 2.5–12 Å. Arbitrary twist is not accepted because a valid commensurate periodic cell has not been supplied.

The returned coordinates, unit cell and provenance support actual CIF/XYZ/PDB export and the existing periodic `pxrd` calculation. Compare every arrangement from the same original layer and use identical diffraction settings. These cells are unrelaxed assumptions; the input connectivity screen does not certify favorable interlayer packing or experimental stacking. Inspect resulting contacts before calculating or interpreting pores and adsorption. The [periodic characterization method](characterization-methods.md) describes the separate reciprocal-cell powder calculation.

## Interpenetration hypotheses

`network_hypothesis(structure, offsets, repeats=(1,1,1))` creates one to four explicitly translated copies of a supplied periodic geometry. Fractional offsets are chosen by the user and reduced modulo the cell; coincident offsets and intercopy atom separations below 0.8 Å are rejected. Intercopy minimum distance uses the same metric-aware periodic distance engine as pore analysis. Contacts below 2 Å produce a review warning. These thresholds catch obvious collisions; they do not certify chemically sensible packing.

The JSON and viewer retain copy-group labels, offsets, input hash and parent provenance. CIF/XYZ/PDB contain the actual generated coordinates; retain JSON for group assignments. The viewer does not invent bonds between explicitly labeled copies. Parent geometry may itself contain more than one infinite net: a chosen copy count is therefore not a measured interpenetration degree. Neither copying a geometry nor counting quotient-graph components determines entanglement, COFex `cat`, or a stable crystal. Inspect contacts and supplied connectivity with an appropriate topology tool before assigning them. The distinction between independent nets and their interpenetration is discussed in the [RCSR methods paper](https://yaghi.berkeley.edu/pdfPublications/08resource.pdf); the broader reticular construction approach is documented by [pyCOFBuilder](https://lipelopesoliveira.github.io/pyCOFBuilder/).

## Periodic hard-sphere geometry

`pore_analysis(structure, probe_radius=1.82, samples=12000, surface_samples=96, seed=42, grid_size=16)` measures the explicitly supplied unit cell using RDKit's element-specific van der Waals radii. The complete radius table used is included in each result. Probe radius is an input in ångströms; choose it for the intended comparison rather than assuming all published screening datasets use the same convention.

For a point **x**, clearance is the smallest distance to any periodically repeated atom centre minus that atom's van der Waals radius. The implementation enumerates a finite set of lattice images with a bound obtained from the cell metric and its inverse. Orthogonal metric components are separated, which handles short layer spacing efficiently. Extremely skewed, unreduced cells that exceed the image-search bound are rejected rather than approximated with an incorrect nearest image.

Seeded uniform cell samples estimate two different volumes:

| Result | Definition |
| --- | --- |
| Void fraction | Fraction with clearance ≥ 0 |
| Probe-center free fraction | Fraction with clearance ≥ specified probe radius |
| Geometric void volume, cm³/g | Void fraction divided by listed-atom crystal density |
| Probe-center free volume, cm³/g | Probe-center free fraction divided by listed-atom crystal density |
| Largest sampled included-sphere diameter, Å | Twice the largest sampled clearance, including the comparison grid |

The last quantity is a **sampling lower bound** on the largest included sphere. It is not a pore-limiting diameter. The local-clearance histogram is not a conventional pore-size distribution. Free volume is not helium Widom volume, adsorption pore volume, or probe-occupiable volume.

Local probe-exposed surface is the exposed area of the union of atomic spheres expanded by the probe radius. Surface points are drawn from those spheres with probability proportional to their area; a point is exposed if it lies outside every other periodically repeated expanded sphere. The result is reported in Å² per cell and m²/g. This geometric quantity includes closed cavities. It must not be relabeled as a transport-accessible area or BET area. Established packages distinguish local probe overlap from channel accessibility; use [Zeo++](https://www.maciejharanczyk.info/Zeopp/examples.html) for dedicated Voronoi-based pore and accessibility analysis.

A separate periodic voxel graph estimates connectivity. Free voxel centres connect to their six axial neighbors, including across cell boundaries. Winding cycles identify components that span one or more lattice directions. The result reports percolating probe-center fraction, component sizes, and spanning axes. Finite grid resolution may miss narrow or diagonal channels; the connectivity result is not used to silently relabel the separately sampled surface. Increase resolution and compare against a specialist calculation before using these values for screening claims.

The density calculation is `M_cell / (N_A × V_cell × 10⁻²⁴)` in g/cm³, using standard atomic weights. Area conversion is `A_Å² × 10⁻²⁰ × N_A / M_cell` in m²/g. Missing hydrogens, solvent, defects, disorder, and incorrect cell content change mass-normalized results. Input structures are not independently chemically validated.

Volume and surface fractions include 95% Wilson sampling intervals. These intervals omit radius-model, structural, and voxel-discretization uncertainty. Random seeds make a run repeatable with the same numerical/software versions. Requested work is capped by atom count and periodic-image search complexity; the effective volume/surface sample counts and grid size are reported. No network service, trained model, SciPy, or external executable is needed for these calculations.

## Deterministic transport geometry

`geometric_porosity(structure, probe_radius=1.82, grid_spacing=0.6, surface_samples=64, refinement=False)` is the separate deterministic `/science/geometry` endpoint. It uses the same periodic hard-sphere clearance, with continuous-path lower bounds and enclosing-voxel upper bounds for pore-limiting diameter, largest-cavity diameter and transport-connected probe-center volume. The bounds cover grid discretization within that supplied model, not structural or atomic-radius uncertainty. The sampled inflated-sphere surface counts only points with collision-free connections to certified winding channels; unresolved attachments are excluded, and quadrature has no statistical confidence interval.

Check `connectivity.status`: `unresolved_at_this_resolution` means the lower graph found no channel while the enclosing upper graph still can wind. A zero lower volume or surface is then not proof of an inaccessible pore. Refinement compares grids and does not establish convergence. Effective grid/sample counts can be capped. The [geometric-properties method](geometric-properties.md) gives definitions, bound reasoning, API examples and analytical tests. The older `pore_analysis` above remains a distinct sampling calculation and does not acquire these transport certificates.

## Local classical adsorption

The asynchronous `/science/adsorption/jobs` service runs a pinned RASPA3 CPU worker in an isolated process. It requires an explicit full periodic structure and a documented rigid-framework force field, gas composition, temperature, pressure and charge model. Force-field support is checked; missing charges are rejected in CIF-charge mode. Explicit zero-framework-charge mode is an approximation, not inferred neutrality or validated host–guest electrostatics. Reservoir fugacity coefficients default to one (ideal gas), rather than an inferred mixture equation of state. Read the [RASPA model, resources and validation limits](raspa-methods.md).

A separate zero-loading Widom run provides Henry coefficients; GCMC stages provide loadings at each requested pressure. Mixture selectivity uses `(q_i/q_j)/(y_i/y_j)`; Henry selectivity uses the coefficient ratio. Missing or nonpositive estimates produce no ratio. Short default runs and five-block error estimates are not proof of equilibration, convergence or agreement with experimental adsorption. No flexibility, chemisorption or framework synthesis is modeled.

Optional zero-loading heat uses two extra independent Widom runs at `T ± min(5 K, 0.025T)`, evaluating `Qst0 = R Δln(K_H,p)/Δ(1/T)` with `R = 0.00831446261815324 kJ mol⁻¹ K⁻¹`. This finite difference has temperature-step bias and sampling noise. Its input-interval envelope is not a joint 95% confidence interval, and it does not supply a finite-loading heat curve.

Retain the actual JSON/CSV and ZIP of inputs, outputs and hashes. Jobs are queued/running until a terminal result is recorded; cancellation waits for worker shutdown and preserves completed stages. Failed or cancelled jobs may contain partial results. The [agent API workflow](../skills/cofdraw/references/science-api.md#local-adsorption-jobs) documents polling, cancellation and the helper. Geometric surfaces and simulated isotherms must remain distinct from experimental BET analysis.

## BET fitting requires adsorption data

`fit_bet(relative_pressures, uptake_mmol_g, fit_min=0.05, fit_max=0.30, cross_section_nm2=0.162)` requires actual supplied adsorption-branch data. Relative pressure must increase strictly between zero and one, and uptake is positive in mmol/g. The selected interval must contain at least three points. The default interval is an editable starting range, not a claim that every microporous COF should be fitted over it.

The linear fit uses `y = x / [n(1−x)] = intercept + slope × x`, giving `n_m = 1/(slope+intercept)` and `C = 1+slope/intercept`. Area is `n_m × 10⁻³ × N_A × molecular_cross_section_nm² × 10⁻¹⁸` in m²/g. The fit returns its points, fitted values, R², monolayer capacity, C, and basic consistency diagnostics: increasing `n(1−x)`, positive C, monolayer uptake/pressure inside the selected interval, and pressure agreement within 20%. Nonphysical intercept/capacity is rejected. A failed diagnostic remains visible.

This is a transparent fit of a **user-selected interval**, not a complete automatic Rouquerol selection procedure or certification of an experimental BET result. High R² is insufficient, particularly for microporous materials. The conventional nitrogen cross-section of 0.162 nm² is configurable and must match the adsorptive and conditions. Consult the [IUPAC physisorption report](https://doi.org/10.1515/pac-2014-1117). No BET value is calculated from a 2D picture or geometric CIF surface area.

## Validation coverage

Tests check a real chiral monomer's descriptor formula, conformer reproducibility, stereochemical SDF roundtrip, explicit nonconvergence, CIF symmetry and special positions, supercell atom/volume counts, fixed-column PDB coordinates, rejection of unresolved occupancy/symmetry, and triclinic distances against wide explicit image enumeration. A single sphere in a large cubic cell supplies analytic void-volume, surface-area, density, and unit-conversion targets; a small cell checks self-image occlusion. Separate voxel fixtures distinguish an isolated pocket from a wrapping channel. Synthetic BET data with known monolayer capacity and C verify the fit and m²/g conversion. These are numerical/software checks, not validation against experimental COF adsorption data.
