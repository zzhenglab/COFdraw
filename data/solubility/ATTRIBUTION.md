# BigSolDB 2.2 reference measurements

The derived table in this directory uses the equilibrium solubility and methods data from Lev Krasnov, Dmitry Malikov, Marina Kiseleva, Sergei Tatarin, Sergey Sosnin, and Stanislav Bezzubov. It is distributed under [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/). [Archived BigSolDB 2.2 release](https://doi.org/10.5281/zenodo.22648301) · [Dataset paper](https://doi.org/10.1038/s41597-025-05559-8)

COFdraw canonicalizes solute and solvent SMILES, selects 19 neat solvents, keeps equilibrium records within 5 K of 298.15 K, and retains the closest observation for each solute, solvent, publication and solid-form label. Original study DOIs, solid-form labels, recorded mole fractions and the dataset's converted concentrations are retained. The concentration basis is moles of solute per liter of pure solvent.

Run `scripts/prepare_cofdraw_solubility.py` to regenerate the table and benchmark. See `validation.json` for source hashes, the derived CSV checksum and scaffold-held-out errors.
