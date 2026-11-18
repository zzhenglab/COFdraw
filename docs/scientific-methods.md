# Molecular and periodic structure calculations

COFdraw separates molecular conformers, explicitly periodic structures, finite drawing illustrations, hard-sphere geometry, classical adsorption simulations, and adsorption-data fitting. The result records the method, parameters, units, and limitations. A finite drawing never becomes a validated crystal merely by adding a cell around it.

## Monomer properties and conformers

`monomer_analysis(smiles, conformers=4, seed=42, max_iterations=300)` accepts one complete connected molecular graph, with at most 250 input atoms. Wildcards and continuation points are rejected. Molecular formula, molecular weight, exact mass, formal charge, ring counts, rotatable bonds, hydrogen-bond counts, topological polar surface area, Crippen logP, and fraction of sp3 carbons are calculated with RDKit. These are molecular descriptors, not predictions of a COF's electronic or adsorption properties. See the [RDKit descriptor reference](https://www.rdkit.org/docs/source/rdkit.Chem.rdMolDescriptors.html).

Explicit hydrogens are added before seeded, single-threaded ETKDGv3 conformer generation. Available MMFF94s parameters are preferred; UFF is the fallback. Each conformer's optimizer status and energy is retained. Selection prefers the lowest-energy converged conformer; if none converged, the result explicitly reports that limitation. Requested conformer count, pruning, seed, iteration limit, force field, and RDKit version are recorded. The conformer search is bounded and does not establish a global minimum. [RDKit conformer guide](https://www.rdkit.org/docs/GettingStartedInPython.html), [force-field status definitions](https://www.rdkit.org/docs/source/rdkit.Chem.rdForceFieldHelpers.html).

XYZ and SDF contain actual three-dimensional Cartesian coordinates in ångströms. SDF retains the molecular graph, stereochemistry, charge, and isotope information; XYZ is a coordinate format. Force-field energies are in kcal/mol and should only compare conformers of the same molecular graph under the same force field. They are not heats of formation, binding energies, orbital energies, or framework stability. The bundled calculation does **not** run xTB or DFT. XYZ is suitable as an initial geometry for separately installed quantum-chemistry programs after inspection.

The monomer response also provides empirical Gasteiger partial charges for all explicit atoms, including hydrogens, using 12 iterations. Missing parameters cause explicit omission instead of zero-valued stand-ins. These graph-based charges are not quantum electron density, oxidation states, HOMO/LUMO energies, or a band gap. [RDKit partial-charge implementation](https://www.rdkit.org/docs/source/rdkit.Chem.rdPartialCharges.html).

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

## Finite layer illustrations

`patch_3d(recipe, layers=3, layer_spacing=3.4, stacking="AA", slip=None, bond_length=1.4, twist_degrees=15)` makes a separate **nonperiodic illustration** from an existing finite drawing. Its median real-atom bond length is scaled to the requested ångström value. AA aligns layers; AB uses alternating offsets `0, s, 0, s`; ABC uses `0, s, 2s, 0`; progressive slip uses `0, s, 2s, 3s`; twisted layers rotate successively around the drawing centroid by the chosen angle. Offset modes require an explicit slip vector in ångströms; pore shape does not infer it. These are chosen illustration transforms, not crystallographic stacking assignments such as AB1 or ABC1. Layer translations/rotations are recorded individually. Bond connectivity and internal distances are retained within each finite layer. No interlayer bonds, missing hydrogens, or optimized atomic geometry are invented. XYZ uses `X` for continuation placeholders. This output is not accepted as a periodic pore-analysis structure.

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
