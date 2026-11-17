# Contributing

COFdraw welcomes reproducible bug reports, drawing improvements, and reviewed chemistry additions. Describe the user-visible problem and include a minimal input or saved project that you can share publicly.

## Development setup

Use Windows x64 and Python 3.11 for the desktop application and packaging.

```powershell
py -3.11 -m venv .venv
.\.venv\Scripts\python.exe -m pip install -r requirements-cofdraw-desktop.txt
.\.venv\Scripts\python.exe -m pip install pytest httpx "scipy>=1.11"
.\.venv\Scripts\python.exe scripts\launch_cofdraw.py
```

The launcher uses this repository's portable `data` directory by default. Use `--data-dir` to select a separate writable data directory. For automation, follow the [COFdraw skill](skills/cofdraw/SKILL.md).

## Validate changes

```powershell
.\.venv\Scripts\python.exe -m pytest -q
node scripts\test_cofdraw_library.cjs
```

Node.js is needed only for the JavaScript library checks. SciPy is a development dependency for the image-overlay comparison test; it is not required by the drawing application. Exercise the desktop window and affected exports when changing the UI, save bridge, or packaging. Browser/API checks alone do not establish that a native Save dialog works.

For chemistry changes, preserve stereochemistry, source-qualified identity, provenance, and explicit continuation points. A successful parse is not a chemistry validation. A template is a construction choice, not evidence of a periodic topology. Keep reference-only catalog entries distinct from implemented transforms.

## Build a Windows release

```powershell
.\.venv\Scripts\python.exe scripts\build_cofdraw_desktop.py
```

The builder produces `dist/COFdraw/`, `dist/COFdraw-Windows.zip`, and a SHA-256 checksum file. It refuses an existing output directory; use `--output-dir dist/another-build` for a new build. Keep the EXE, `COFdraw.exe.config`, `_internal`, and `data` together. The build environment and executable output are excluded from Git.

The builder validates and includes the reviewed drawing catalog and periodic model templates. It retains its build cache and copies into a fresh destination. A rebuilt binary has a new hash; run `scripts/smoke_cofdraw_release.py` against that executable and its bundled data before distributing it. Browser smoke checks require development Playwright and installed Chrome or Edge.

Also run `scripts/smoke_cofdraw_download.py --zip <release.zip> --output-dir <fresh-directory>` on Windows. It marks the archive's extracted EXE, configuration, and DLLs as downloaded, then verifies a normal-user native process and retained DLL/configuration markings. Windows may clear only the EXE's mark when its first-run confirmation is accepted; the report distinguishes that from file-content changes. An elevated or unmarked-DLL local startup is insufficient evidence for downloaded-release readiness. See [Windows startup validation](docs/windows-startup.md).

## Catalog and documentation

The current catalog, registry, and recipes are bound together by SHA-256 hashes in `data/framework-drawing/current-bundle-manifest.json`. Keep those records consistent. Preserve source identifiers, and do not include private publications, local source paths, datasets, model weights, or generated exports in commits.

For detailed scientific scope and recorded validation, see the [catalog guide](docs/cofdraw-catalog-v3.md) and [release record](docs/cofdraw-release-validation.json). That record describes the original 0.3.0 application's broader validation run; the public desktop repository contains its focused drawing tests.

Current structure methods are documented in [scientific methods](docs/scientific-methods.md). Description interpretation and API data handling are documented in [description methods](docs/describe-methods.md) and [optional API setup](docs/llm-api.md). The [v0.6 validation record](docs/v0.6-validation.json) identifies the tested executable and browser workflows. The optional provider path is tested with a localhost fixture; paid external model calls are not required for tests.
