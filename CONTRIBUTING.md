# Contributing

COFdraw welcomes reproducible bug reports, drawing improvements, and reviewed chemistry additions. Describe the user-visible problem and include a minimal input or saved project that you can share publicly.

## Development setup

Use Windows x64 and Python 3.11 for the desktop application and packaging.

```powershell
py -3.11 -m venv .venv
.\.venv\Scripts\python.exe -m pip install -r requirements-cofdraw-desktop.txt
.\.venv\Scripts\python.exe -m pip install pytest httpx
.\.venv\Scripts\python.exe scripts\launch_cofdraw.py
```

The launcher uses this repository's portable `data` directory by default. Use `--data-dir` to select a separate writable data directory. For automation, follow the [COFdraw skill](skills/cofdraw/SKILL.md).

## Validate changes

```powershell
.\.venv\Scripts\python.exe -m pytest -q
node scripts\test_cofdraw_library.cjs
```

Node.js is needed only for the JavaScript library checks. Exercise the desktop window and affected exports when changing the UI, save bridge, or packaging. Browser/API checks alone do not establish that a native Save dialog works.

For chemistry changes, preserve stereochemistry, source-qualified identity, provenance, and explicit continuation points. A successful parse is not a chemistry validation. A template is a construction choice, not evidence of a periodic topology. Keep reference-only catalog entries distinct from implemented transforms.

## Build a Windows release

```powershell
.\.venv\Scripts\python.exe scripts\build_cofdraw_desktop.py
```

The builder produces `dist/COFdraw/`, `dist/COFdraw-Windows.zip`, and a SHA-256 checksum file. It refuses an existing output directory; use `--output-dir dist/another-build` for a new build. Keep the EXE, `_internal`, and `data` together. The build environment and executable output are excluded from Git.

The first public v0.3.0 release repackages the existing validated executable with a clean portable catalog and documentation. The standalone builder lets contributors produce a fresh executable from this desktop source. A rebuilt binary can have a different hash; validate that binary before distributing it.

## Catalog and documentation

The current catalog, registry, and recipes are bound together by SHA-256 hashes in `data/framework-drawing/current-bundle-manifest.json`. Keep those records consistent. Preserve source identifiers, and do not include private publications, local source paths, datasets, model weights, or generated exports in commits.

For detailed scientific scope and recorded validation, see the [catalog guide](docs/cofdraw-catalog-v3.md) and [release record](docs/cofdraw-release-validation.json). That record describes the original 0.3.0 application's broader validation run; the public desktop repository contains its focused drawing tests.
