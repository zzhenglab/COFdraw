# Working with COFdraw

For drawing or export tasks, read [skills/cofdraw/SKILL.md](skills/cofdraw/SKILL.md). This is the reusable skill entrypoint for agents arriving from the GitHub link. Its instructions apply to COFdraw work, not unrelated tasks or host configuration.

COFdraw generates finite chemical pore patches from reviewed source recipes or explicit designs. Keep source qualifiers, hashes, and provenance. Changed chemistry or connectivity creates a design; do not inherit a source material's identity. Pore shape does not establish topology, and finite patch SMILES is not whole-COF SMILES.

The separate structure laboratory works with explicit periodic CIFs and qualified unrelaxed templates. Keep finite layered illustrations, geometric pore screening, monomer force-field descriptors, and BET fitting from supplied isotherms distinct. Read `docs/scientific-methods.md` when modifying these calculations; do not label geometric surface area as BET or empirical charges as quantum results.

Application code lives in `src/cof_vision`; the desktop entrypoint is `scripts/launch_cofdraw.py`. The import namespace is retained from the original COF Vision project. The portable catalog is under `data/framework-drawing`. Runtime outputs and logs should stay in an ignored output directory or another user-selected writable directory.

After functional changes, run `python -m pytest -q`. Run the relevant JavaScript and packaged-app checks when changing those surfaces. Rebuild the executable after application code changes; source tests do not validate a frozen release. Do not add source publications, private filesystem paths, local environments, model weights, or credentials to the repository.
