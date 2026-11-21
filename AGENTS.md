# Working with COFdraw

For drawing or export tasks, read [skills/cofdraw/SKILL.md](skills/cofdraw/SKILL.md). This is the reusable skill entrypoint for agents arriving from the GitHub link. Its instructions apply to COFdraw work, not unrelated tasks or host configuration.

COFdraw generates finite chemical pore patches from reviewed source recipes or explicit designs. Keep source qualifiers, hashes, and provenance. Changed chemistry or connectivity creates a design; do not inherit a source material's identity. Pore shape does not establish topology, and finite patch SMILES is not whole-COF SMILES.

The separate structure laboratory works with explicit periodic CIFs and qualified unrelaxed templates. Keep finite layered illustrations, geometric pore screening, monomer force-field descriptors, and BET fitting from supplied isotherms distinct. Read `docs/scientific-methods.md` when modifying these calculations; do not label geometric surface area as BET or empirical charges as quantum results.

Application code lives in `src/cof_vision`; the desktop entrypoint is `scripts/launch_cofdraw.py`. The import namespace is retained from the original COF Vision project. The portable catalog is under `data/framework-drawing`. Runtime outputs and logs should stay in an ignored output directory or another user-selected writable directory.

Description interpretation produces a reviewable plan before rendering. Rules mode is local; the optional LLM may select existing catalog IDs only. Validate choices with the local chemistry engine, preserve unresolved instructions, and never invent a source identity or periodic crystal from prose. API calls require an explicit configured API-mode request. Never log, persist, export or echo provider credentials, and do not include saved projects or filesystem data in provider prompts.

After functional changes, run `python -m pytest -q`. Run the relevant JavaScript and packaged-app checks when changing those surfaces. Rebuild the executable after application code changes; source tests do not validate a frozen release. Do not add source publications, private filesystem paths, local environments, model weights, or credentials to the repository.

Before distributing a Windows ZIP, run `scripts/smoke_cofdraw_download.py` against that archive. Mark its EXE/configuration/DLLs as downloaded before launch, require a verified normal-user native launch, and retain every DLL/configuration mark. Windows may clear the EXE's mark after its per-app first-run confirmation; report that separately and verify all file hashes remain unchanged. Elevated or unmarked-DLL startup does not cover the downloaded-file condition. Keep the process-scoped `COFdraw.exe.config` beside the EXE; do not change machine policy or unblock unrelated files.

Mac betas share the desktop entrypoint, scientific code, web interface and catalog with Windows. Dispatch the **macOS beta** workflow from `main` when a Mac build is wanted. Select **publish** to release both native archives after the extracted-app checks pass. Mac prereleases have their own tags.
