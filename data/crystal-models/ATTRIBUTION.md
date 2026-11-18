# Crystal model attribution

These five COFdraw model examples were generated with **pyCOFBuilder 0.0.9.4**, by Felipe Lopes de Oliveira, using its building-block library and reticular construction methods.

- Upstream project: [lipelopesoliveira/pyCOFBuilder](https://github.com/lipelopesoliveira/pyCOFBuilder)
- Method publication: [DOI 10.1021/acs.jcim.3c01918](https://doi.org/10.1021/acs.jcim.3c01918)
- Upstream license: [MIT License, copyright © 2023 Felipe Lopes de Oliveira](third-party-licenses/pyCOFBuilder-LICENSE.txt), copied unchanged from the installed 0.0.9.4 distribution.
- Regeneration scripts in the source checkout: `scripts/generate_crystal_templates.py` for the original four models, then `scripts/generate_matched_crystal_templates.py` for the TFB–TAPB model and precursor audit.
- Exact model identifiers, atom/formula counts, hashes, and generation notes: [manifest.json](manifest.json).

The examples are generated, **unrelaxed reticular designs**. They are not experimental crystal assignments, optimized structures, or conversions of the user's active 2D drawing coordinates. Exact precursor-pair and imine-chemistry matches are recorded in the manifest and [precursor audit](precursor-mapping.audit.json). Topology and stacking are explicit model choices. Missing experimental agreement or a large geometric surface area must not be interpreted as demonstrated stability, synthesis feasibility, or measured BET area.

The audit checks all generated atoms, explicit hydrogens, finite precursor-fragment connectivity, zero lattice winding of each molecular fragment, and precursor-minus-water stoichiometry. CIF coordinates do not supply bond orders; distance-inferred connectivity is compared against the selected imine product. The generator's unresolved precursor SMILES metadata is not used as chemical identity. The stock hydroxylated Tp–Pa example retains O–H and lacks N–H, so it is excluded from β-ketoenamine matching.

For the tetrahedral D4 example, the generation script flattens nested connector-provenance labels before functional-group substitution to address a metadata-array incompatibility in version 0.0.9.4. The workaround does not alter atom types or coordinates. Pymatgen writes the resulting full-cell coordinates to P1 CIF. Runtime use of these bundled CIFs does not require pyCOFBuilder or Pymatgen.
