# RASPA3 adsorption calculations

COFdraw uses the real **RASPA3 CPU Monte Carlo engine**, distributed as the pinned
`raspalib==3.1.0` Python package. The Windows application includes its native CPU
modules and MIT license; an installation of RASPA, WSL, a compiler, GPU or network
connection is not required for a calculation. Source installations obtain the
same dependency through the normal project requirements. Each calculation runs
in a separate process so cancellation or a native failure does not terminate the
drawing application.

The engine requires an **explicit periodic structure**. An illustrated stack or
finite drawing is not automatically a periodic adsorption model. Full atom
occupancies, valid cell geometry, supported elements, absence of severe atomic
collisions and explicit charge handling are checked before starting. The
framework stays rigid; generic parameters do not establish that a particular
COF is stable, experimentally accessible or accurately described by this model.

## Models, units and charges

The immutable, checksum-bound profile is
[`data/raspa/profile.json`](../data/raspa/profile.json). It records every gas site,
position, mass, charge, Lennard–Jones parameter and attribution. Original RASPA
force-field files and licenses are included beside it.

| Species | Interaction model |
| --- | --- |
| CO₂ | Three rigid sites, García-Sánchez et al. (2009); charged C/O sites |
| N₂ | Two N sites and a central charge site, Martín-Calvo et al. (2011) |
| CH₄ | One united-atom site, Martin et al. (2001), as supplied in RASPA's examples |
| H₂ | One neutral spherical Buch site, ε/kB = 34.2 K, σ = 2.96 Å |
| Ar | One neutral site, ε/kB = 124.070 K, σ = 3.38 Å, RASPA ExampleMOFsForceField |
| Framework | RASPA ExampleMOFsForceField DREIDING organic-element parameters and UFF parameters for its listed metals |

These are the specified published/example models, **not a blanket TraPPE
parameterization**. H₂ is a classical spherical approximation without a
quadrupole, Feynman–Hibbs correction, ortho/para distinction or nuclear quantum
effects. Cryogenic H₂ predictions need particular caution. The explicit H₂
parameters are documented in the [published correction](https://doi.org/10.1021/jp110139d)
to a molecular adsorption study, citing [Buch's potential](https://doi.org/10.1063/1.466854).

The profile uses Lorentz–Berthelot mixing, a **12 Å shifted** Lennard–Jones cutoff
and no tail correction. Unit-cell repeats are increased until every simulation
cell face-to-face height is at least 24 Å, including for skew cells. The adapter
checks native simulation volume and framework mass against the supplied structure.
Calculations are bounded to 12,000 framework atoms after replication and a finite
atom-cycle work limit; the job manager also controls runtime and cancellation.

Temperature is in K, pressure and fugacity in Pa, distances in Å, absolute loading
in **mol/kg-framework** (numerically equal to mmol/g), and Henry coefficients in
**mol/kg-framework/Pa**. Excess uptake is not reported: its conversion additionally
requires an appropriately defined accessible pore volume and fluid density.

Choose either complete, explicit CIF partial charges or the explicitly labeled
**zero framework charge** approximation. CIF mode refuses missing charges and
non-neutral framework cells. It does not silently neutralize a cell or guess a
charge assignment. Guest CO₂/N₂ charges remain active even with zero framework
charges; the engine uses electrostatics when any site is charged. These are
classical partial charges, not a quantum calculation.

## GCMC, mixtures, Henry coefficients and heat

GCMC uses rigid-molecule translations, rotations where applicable, and actual
configurational-bias insertion/deletion. Species loading and uncertainty are read
from native engine statistics. Physical gas mole fractions must sum to one.
The reservoir input is the explicit partial fugacity `f_i = y_i * phi_i * P`;
`phi_i` defaults to **1**, an ideal-gas reservoir approximation. Supply validated
coefficients for a nonideal reservoir. The adapter does not silently claim a
mixture equation of state.

RASPA's current Python constructor starts swappable components at equal internal
fractions. The adapter supplies `N * y_i * phi_i` as the internal coefficient so
that the insertion/deletion activity is exactly `y_i * phi_i * P`. Both the
physical and internal values are recorded. Upstream mixture composition labels,
EOS-derived values and excess loadings are not presented as the physical input.
An independent ideal-gas single-component and 20/80 binary calculation checks
`<N_i> = y_i * P * V / (kB * T)` against native sampling uncertainty.

**Henry mode is a separate zero-loading calculation:** it creates no adsorbate
molecules and performs Widom insertions only. It is not a finite-loading Widom
statistic relabeled as a Henry coefficient. The native Rosenbluth mean and 95%
confidence half-width are converted using the actual simulation-cell volume and
framework mass:

`KH = <W> * V[Å³] * 1e-27 / (kB[J/K] * T[K] * M_framework[g/mol])`.

The adsorption workflow can estimate **zero-loading Qst** from independent Henry
calculations at neighboring temperatures using `R * d ln(KH) / d(1/T)`. Its result
retains the temperature interval, seeds, Henry uncertainties and qualification
as a finite-difference estimate. A poorly resolved or zero Henry coefficient
cannot support a reliable logarithm or heat estimate. The adapter does not
report finite-loading mixture enthalpy: upstream covariance/ideal-gas reference
handling is not sufficiently robust for this integration.

Mixture separation factors use actual loadings divided by the specified reservoir
fractions, with zero/poorly sampled loadings left unresolved. Henry ratios are
**infinite-dilution selectivities**; they do not establish finite-pressure mixture
selectivity. No reported result is an experimental isotherm, BET measurement,
chemical selectivity, membrane permeability or kinetic separation prediction.

## Reproducibility and numerical checks

Each job saves its structure, force-field profile hash, temperature, pressure,
fractions, fugacity coefficients, cell repeats, cycle counts, random seed, engine
version and actual statistics. GCMC also saves native per-block loading text.
Five-block 95% confidence intervals describe finite sampling only. Short smoke
runs verify execution and units; they do not establish equilibration or
production-quality adsorption predictions. Compare longer runs and independent
seeds when interpreting material trends.

The selected library auto-dispatches among base, AVX, AVX2 and AVX512 CPU variants.
The build includes all four Windows modules and package metadata. A public native
RASPA3 3.0.24 Windows CLI was evaluated and rejected after its eigensolver produced
incorrect CO₂/N₂ inertia and invalid values; that executable is **not bundled**.
The chosen 3.1.0 library passes the molecular geometry, skew-cell volume,
ideal-gas activity, five-gas, reproducibility and charged-structure checks.
The adapter uses the native statistics API with formatted final output disabled;
this avoids a separate upstream mixture-enthalpy formatting failure. Reinsertion
is omitted after an upstream mixture failure; insertion/deletion and ordinary
translation/rotation still sample the declared grand-canonical ensemble.

The wheel's Python stub describes six-argument cell angles as degrees, whereas
the native constructor takes radians. The adapter explicitly supplies radians
and independently checks the resulting volume.

## Sources and redistribution

- [RASPA3 primary paper](https://doi.org/10.1063/5.0226249),
  [official documentation](https://iraspa.github.io/RASPA3/), and
  [official Python distribution](https://pypi.org/project/raspalib/3.1.0/).
- [RASPA3 example gas definitions](https://github.com/iRASPA/RASPA3/tree/v3.0.24/examples/basic/12_mc_henry_coefficient_of_co2_n2_methane_in_mfi).
- [RASPA2 ExampleMOFsForceField](https://github.com/iRASPA/RASPA2/tree/master/forcefield/ExampleMOFsForceField).
- [CO₂ force field](https://doi.org/10.1021/jp810871f) and
  [N₂/Ar example attribution](https://doi.org/10.1039/C1CP20168A).
- [DREIDING](https://doi.org/10.1021/j100389a010) and
  [UFF](https://doi.org/10.1021/ja00051a040).

RASPA2 and RASPA3 use MIT licenses; their notices are preserved. The exact Windows
wheel URL/SHA256 and profile hash are in
[`engine-manifest.json`](../data/raspa/engine-manifest.json). No global runtime
installation, machine policy change or download-time unblocking is part of this
calculation integration.
