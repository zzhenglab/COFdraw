# Periodic model characterization

These tools work on the explicitly supplied full unit cell. A generated COF model is a geometric hypothesis; its calculated formula, density, or diffraction pattern does not make it an experimentally identified phase. The broader [scientific methods](scientific-methods.md) describe CIF import, conformers, pore geometry, and the separate BET fit from measured isotherm data.

## Composition, cell, and contact summary

`cof_vision.cofdraw_characterization.structure_summary(structure)` accepts the same normalized structure as the structure laboratory: three Cartesian cell-vector rows in angstrom and a list of fully occupied atomic sites with element symbols and fractional or Cartesian coordinates. The cell must be right-handed and nondegenerate. Finite layered patches are not periodic inputs.

The summary contains the full-cell Hill formula, its greatest-common-divisor reduced formula, atom counts, and standard atomic-weight sum `cell_molar_mass_g_mol`. This mass is per mole of the supplied cells, not necessarily per mole of an independently defined COF repeat unit. Element rows report atomic fraction, mass fraction, and `site_density_mmol_g = 1000 * count / cell_molar_mass_g_mol`. This last quantity is an **elemental atom inventory**, useful for comparing N/O/B content; it is not a concentration of accessible binding sites or catalytically active sites.

Cell volume is the determinant of the row-vector cell matrix. Lengths are vector norms; alpha, beta, and gamma are the angles between b/c, a/c, and a/b, respectively. Density is

```text
rho [g/cm3] = cell_molar_mass [g/mol] / (N_A [1/mol] * volume [angstrom3] * 1e-24)
N_A = 6.02214076e23 mol^-1
```

All supplied sites, including guests and counterions, contribute. No missing hydrogen atoms are added. Resolve partial occupancy/disorder before import; no silent fractional-site mass approximation is made. Hydrogen-free models receive an explicit completeness warning. Density is model crystal density, not a powder's tapped or bulk density.

Contact checks use the periodic covalent-radius search in `cofdraw_network.py`, with its default bond scale and tolerance recorded in `contact_checks.method`. Implausibly short distances are warnings, not an energy calculation. `minimum_candidate_distance_angstrom` is the smallest distance inside that algorithm's finite candidate envelope; `null` means no candidates, not infinite atom separation. It is **not a general global nearest-distance result**. If resource bounds prevent a contact calculation, composition and cell properties remain available and the failed contact check is marked unavailable. The separate full network inspection may use user-selected cutoffs.

## Simulated powder XRD

```python
from cof_vision.cofdraw_characterization import powder_xrd

result = powder_xrd(
    structure, wavelength=1.5406,
    two_theta_min=2, two_theta_max=50,
    step=0.02, fwhm=0.15, b_iso=0.0,
)
```

Angles are degrees; wavelength is angstrom; the common isotropic `b_iso` is in angstrom squared. Default wavelength 1.5406 is a nominal single Cu K-alpha1 approximation. It does not model the Cu emission spectrum or K-alpha doublet. Real source spectra have finite line structure; see the [NIST Cu K-alpha measurement](https://www.nist.gov/publications/high-precision-measurement-x-ray-cu-k-alpha-spectrum).

For row-vector cell matrix A, reciprocal vectors are rows of `inverse(A).T`, without a 2-pi factor. Every nonzero integer hkl within the reciprocal sphere is considered. For `q = hkl @ inverse(A).T`, `d = 1 / |q|` and `2 theta = 2 asin(wavelength / (2d))`. The bound `|h_i| <= q_max * |a_i|` ensures skewed cells are covered. No space group, preferred net, or symmetry is guessed from a drawing.

Atomic scattering uses the neutral-atom IT92 four-Gaussian-plus-constant coefficients distributed by Gemmi 0.7.5. The bundled H-through-Cf coefficient table is exported without renormalization; its SHA-256 accompanies each result. Gemmi is not a runtime dependency. With `s = sin(theta)/wavelength`, the calculation is

```text
f_element(s) = sum_i a_i * exp(-b_i * s^2) + c
F_hkl = sum_j f_element_j(s) * exp(2*pi*i*hkl.fractional_j) * exp(-B*s^2)
I_hkl = |F_hkl|^2 * (1 + cos(2*theta)^2) / (sin(theta)^2 * cos(theta))
```

The scattering parameterization and squared-s argument follow the [Gemmi scattering documentation](https://gemmi.readthedocs.io/en/stable/scattering.html). The powder Lorentz-polarization factor assumes unpolarized Bragg-Brentano geometry. Enumerating all reciprocal vectors already includes Friedel pairs and multiplicity; the algorithm does not multiply by multiplicity again. Equal-angle contributions are summed. This direct-summation powder formalism is also documented by [pymatgen's diffraction implementation](https://pymatgen.org/pymatgen.analysis.diffraction.html#pymatgen.analysis.diffraction.xrd.XRDCalculator).

Peak centers merge within 0.000001 degree, independently of the selected display width. Numerical extinctions use `|F|^2 <= (sum atomic_number)^2 * 1e-14`; the exact threshold is recorded. Each peak reports every contributing hkl, multiplicity, d-spacing, its arbitrary integrated intensity, and its integrated intensity relative to the strongest displayed peak. Those arbitrary values are not absolute calibrated scattering intensities, and must not be interpreted as phase weight fractions.

The profile adds area-normalized Gaussians with the chosen constant FWHM, truncated at plus/minus five sigma. Reflections just outside the display window contribute their tails. Grid spacing is at most the requested step, with the actual value recorded; both endpoints are included. Profile height and peak integrated intensity are independently normalized to a maximum of 100. A window with only neighboring tails is explicitly flagged. The CSV export contains `two_theta_degrees,relative_intensity`; the JSON response retains the peak table, parameters, coefficient hash, and input-structure hash.

This model omits anomalous dispersion, ionic/bonding charge corrections, absorption, texture/preferred orientation, multiple scattering, and instrumental spectrum. It supports one common B factor; per-site CIF displacement tensors are not preserved by the current importer. The chosen width is not a measurement of crystallite size, strain, or stacking disorder. A pattern match by eye is not a phase identification or refinement. For disorder, texture, calibrated intensities, or quantitative fitting, export the CIF for a dedicated diffraction/refinement workflow.

## Workload limits and verification

The existing periodic validator allows at most 2000 fully occupied sites. XRD additionally limits the rectangular reciprocal enumeration to 600000 candidates, accepted vectors to 60000, atom-reflection evaluations to 24 million, and the profile to 20001 points. Over-budget requests fail with a clear error instead of returning a truncated lattice. Use a primitive cell, a narrower angular range, or a longer wavelength to reduce work. This implementation bounds `sin(theta)/wavelength` to at most 2 inverse angstrom. Step must be at most half the chosen FWHM. These are numerical workload/resolution limits, not a certification of model accuracy.

Tests cover cubic Bragg positions and multiplicities, body- and face-centered systematic absences, two-element amplitude interference, triclinic reciprocal metrics, origin and rigid-rotation invariance, supercell extinction/invariance, the squared Debye-Waller attenuation, Gaussian width, CSV fidelity, reproducibility, missing form factors, and bounded failures. The four bundled 2D/3D reticular models are exercised using their actual coordinates. Reference scattering values were independently obtained from Gemmi 0.7.5; CI does not require Gemmi.

## Measured peak positions

The PXRD panel compares user-picked 2θ centers against the current periodic model at the selected wavelength. A symmetric angular window and a cutoff on relative integrated reflection intensity control inclusion. The CSV lists every candidate family with its signed observed-minus-calculated deviation, measured peaks outside the model range, unmatched observations and calculated-only centers. More than one observation may overlap the same candidate. The procedure leaves the cell and zero point fixed. Preferred orientation, stacking disorder and a limited number of reflections can leave several structures consistent with the data. [Statistical representation of stacking disorder in layered COFs](https://doi.org/10.1021/acs.chemmater.1c04365)


## Coefficient provenance and license

Application algorithms and tests retain the repository's MIT license. The generated [coefficient module](https://github.com/zzhenglab/COFdraw/blob/main/src/cof_vision/cofdraw_xray_factors.py) retains its upstream MPL-2.0 notice. It comes from [Gemmi 0.7.5 `include/gemmi/it92.hpp`](https://github.com/project-gemmi/gemmi/blob/v0.7.5/include/gemmi/it92.hpp), copyright Global Phasing Ltd. The complete upstream license is included in [Gemmi-MPL-2.0.txt](third-party-licenses/Gemmi-MPL-2.0.txt). The source module remains available in this repository alongside the executable distribution's documentation.

The export operation is `gemmi.Element(z).it92.get_coefs()` for z=1..98 in an unmodified fresh Gemmi 0.7.5 process; values are a1..a4, b1..b4, c. The canonical sorted compact JSON coefficient mapping has SHA-256 `cd6ba74e08c44657c951f3498965d92399ada0bf17ca9b1eec4a1fd10848783a`. No coefficient normalization or interpolation was applied.
