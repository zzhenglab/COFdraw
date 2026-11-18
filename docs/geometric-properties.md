# Level 1: deterministic geometric properties

`geometric_porosity(structure, probe_radius=1.82, grid_spacing=0.6, surface_samples=64, refinement=False)` analyzes an explicit periodic full-cell structure accepted by the CIF/structure validator. It uses elemental hard spheres and the supplied coordinates; validation does not establish that the structure is experimentally correct, relaxed or free of missing atoms. All listed hydrogens are included. The existing `pore_analysis` remains a separate, compatible sampling calculation.

The engine is **COFdraw's bounded periodic-grid approximation**, not Zeo++ or a Voronoi calculation. Zeo++ defines largest included/free spheres and determines accessible channels from a periodic Voronoi network. PoreBlazer provides a primary precedent for lattice percolation and explicitly distinguishes total from network-accessible quantities. No external porosity executable is bundled or invoked by this calculation. [Zeo++ methods overview](https://www.zeoplusplus.org/about.html), [PoreBlazer source and method documentation](https://github.com/SarkisovGitHub/PoreBlazer).

In the desktop, open **3D & properties**, generate a periodic model or import a CIF, then select **Level 1 · geometry**. The default calculation uses a requested 0.6 Å grid step and 64 area samples per atom before work caps. **Compare a second, finer grid** adds a resolution check. Save the JSON report to retain all inputs, bounds, effective counts and the structure hash. The HTTP equivalent is `POST /api/drawing/science/geometry` with `{"structure": <full-cell structure>, "probe_radius": 1.82, "grid_spacing": 0.6, "surface_samples": 64, "refinement": false}`. The older `/science/porosity` endpoint retains its existing sampling behavior.

| Output | Meaning |
| --- | --- |
| `density_g_cm3` | Listed full-cell atomic mass divided by cell volume |
| `void_fraction` | Grid-center fraction outside all atomic van der Waals spheres, for a point probe |
| `lcd_angstrom` | Largest cavity diameter: twice the largest sampled clearance; a lower estimate |
| `pld_angstrom` | Pore-limiting diameter: twice a certified periodic-path clearance; a lower estimate |
| `transport_accessible_probe_center_volume_cm3_g` | Grid estimate of probe-center space connected to a periodic winding channel |
| `transport_accessible_surface_m2_g` | Deterministically sampled surface of probe-inflated atomic spheres whose probe centers have certified connections to winding channels |
| `local_probe_exposed_surface_m2_g` | All locally exposed probe-inflated surface samples, including closed cavities |

Accessible probe-center volume is not probe-occupiable volume, helium Widom volume or adsorption pore volume. Neither surface field is a measured BET area. The probe radius is a geometric input in Å; it is not a temperature-dependent gas force field. Atom radii, RDKit version, requested/effective resolution, and the normalized input hash accompany every result.

## Periodic paths and numerical bounds

The grid follows the three cell vectors, including skew cells. Exact bounded minimum-image enumeration provides each grid point's clearance from periodically repeated atomic surfaces. For endpoints with clearances `ca` and `cb` separated by Cartesian distance `L`, the whole line segment has clearance at least `(ca + cb - L) / 2`, capped by its endpoint values. This follows from the 1-Lipschitz property of distance-to-surface. Six-neighbor segments with these capacities form the conservative path graph.

A union-find tracks **integer cell-translation cycles**. Touching opposite faces alone does not establish transport. The rank of the independent winding vectors gives the graph's periodic dimensionality; a diagonal winding vector can involve two cell axes while still having rank one. The highest edge threshold at which any winding cycle exists gives the reported PLD lower estimate.

For an enclosing bound, every fractional voxel lies within a known Cartesian covering radius `h` of its center. Its clearance cannot exceed `c + h` or fall below `c - h`. A continuous periodic path must pass through a winding sequence of touching voxels, so a 26-neighbor graph of the enclosing voxel capacities supplies a PLD upper bound. LCD lies between `2 max(c, 0)` and `2 max(c + h, 0)`. The reported `*_bounds_angstrom` enclose these diameters for the specified hard-sphere model, within floating-point tolerance. They are not physical uncertainty intervals.

Void-volume bounds use fully free and possibly free voxels. Transport-volume bounds use winding components of fully free six-neighbor voxels and possibly free 26-neighbor voxels. The volume estimate itself counts centers connected by certified paths. Bounds can be broad at coarse resolution; a zero lower PLD or accessible volume does not exclude a narrow unresolved channel when the upper bound remains positive.

`connectivity.status` distinguishes `certified_channel`, `no_channel_at_this_probe` and `unresolved_at_this_resolution`. The last status means that the lower graph found no transport path but the enclosing graph permits one. The desktop labels the resulting zero accessible estimates **Unresolved**. For example, a coarse axis grid can miss an open diagonal FCC channel; changing the cell basis or refining the grid can resolve it. Grid estimates need not be invariant to basis or origin at finite resolution, while their conservative bounds continue to apply.

## Accessible surface and reproducibility

Surface quadrature uses a fixed three-dimensional Halton sequence: area-weighted atom selection and directions on each probe-inflated sphere. It has no random seed or Monte Carlo confidence interval. Samples occluded by another atom or periodic image are removed. A remaining sample contributes to accessible area only if an exact atom-to-segment check connects it to one of eight neighboring centers belonging to a winding component. The periodic image box includes every atom image that could intersect that segment. Tangencies use a numerical tolerance of 10⁻⁷ Å.

Samples without a certified connection are excluded and reported as unresolved exposed samples. This count can include both genuinely isolated surfaces and surfaces missed by the grid. The accessible-area estimate therefore retains quadrature and connectivity-resolution error; it is not a rigorous integral lower bound. Local exposed area remains separately named and is never relabeled accessible. The area conversion is `A_Å² × 10⁻²⁰ × N_A / M_cell`; volume per gram is the relevant volume fraction divided by density.

Work is bounded to at most 65,536 grid points and 96 points per axis, with tighter atom/image-dependent limits. Effective counts and a warning identify coarsening. Surface samples also have an explicit work cap; unsupported image searches produce an unavailable surface result rather than an invented value. `refinement=True` computes a second grid at half the requested step when the cap permits and returns both descriptor sets. This is a convergence diagnostic, not an automatic claim of convergence; it does not refine the surface quadrature.

The [regression tests](../tests/test_cofdraw_porosity.py) cover analytic LCD/PLD bounds, sphere surface and unit conversion, an isolated cavity with positive local area but zero accessible area, orthogonal and skew one-dimensional channels, boundary-touching versus winding, diagonal winding rank, deterministic repetition, resource caps, and segment collisions against independent periodic image enumeration. These tests validate algorithms on geometric fixtures, not agreement with experimental COF adsorption measurements. For high-accuracy screening, refine the grid and compare with a dedicated Voronoi or independently validated porosity engine using matched radii and definitions.
