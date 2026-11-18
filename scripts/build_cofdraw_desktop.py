"""Build a self-contained Windows ZIP in dist/ using the pinned desktop environment.

Existing output is never overwritten. Use --output-dir dist/another-build to keep
multiple builds. The application requires its entire directory, not only the EXE.
"""
from __future__ import annotations

import argparse
import hashlib
import importlib.metadata as metadata
import json
from pathlib import Path
import shutil
import subprocess
import sys
import uuid
import zipfile

REPO = Path(__file__).resolve().parents[1]
VERSION = "0.8.1"


def validate_synthesis_resources(directory=None):
    """Fail the build if the offline SA Score table was omitted or changed."""
    sys.path.insert(0, str(REPO / "src"))
    from cof_vision.cofdraw_synthesis import resource_information
    try:
        if directory is not None and not (directory / "sascorer.py").is_file():
            raise RuntimeError("SA Score source was not collected")
        return resource_information(directory)
    except (ImportError, RuntimeError):
        raise SystemExit("The build requires the verified RDKit SA Score implementation and fragment table") from None


def validate_raspa_resources():
    sys.path.insert(0, str(REPO / "src"))
    from cof_vision.cofdraw_raspa import engine_status
    status = engine_status()
    if not status["available"]:
        raise SystemExit("The build requires the pinned RASPA3 CPU runtime: " + status.get("reason", "missing package"))
    return status


def sha(path):
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def inside_repository(path):
    resolved = path.resolve()
    if resolved == REPO or not resolved.is_relative_to(REPO):
        raise SystemExit("Build outputs must be inside this repository")
    return resolved


def validate_catalog(data):
    framework = data / "framework-drawing"
    manifest = json.loads((framework / "current-bundle-manifest.json").read_text("utf-8"))
    for row in manifest["files"]:
        item = (framework / row["path"]).resolve()
        if not item.is_relative_to(framework.resolve()) or not item.is_file() or sha(item) != row["sha256"]:
            raise SystemExit(f"Portable catalog hash mismatch: {row['path']}")
    for key in ("catalog", "registry"):
        if sha(framework / manifest[key]) != manifest[key + "_sha256"]:
            raise SystemExit(f"Portable {key} hash mismatch")
    return manifest


def copy_licenses(destination):
    destination.mkdir()
    packages = []
    for dist in metadata.distributions():
        name = dist.metadata.get("Name", "unknown")
        packages.append({"name": name, "version": dist.version,
                         "license": dist.metadata.get("License-Expression") or dist.metadata.get("License", "See included files")})
        for file in dist.files or []:
            if not any(word in file.name.lower() for word in ("license", "copying", "notice")):
                continue
            source = Path(dist.locate_file(file))
            target = (destination / name / str(file)).resolve()
            if (target.is_relative_to(destination.resolve()) and source.is_file()
                    and source.stat().st_size < 2_000_000):
                target.parent.mkdir(parents=True, exist_ok=True)
                shutil.copy2(source, target)
    python_license = Path(sys.base_prefix) / "LICENSE.txt"
    if python_license.is_file():
        shutil.copy2(python_license, destination / "Python-LICENSE.txt")
    (destination / "packages.json").write_text(json.dumps(sorted(packages, key=lambda p: p["name"]), indent=2), "utf-8")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--catalog-data", type=Path, default=REPO / "data")
    parser.add_argument("--output-dir", type=Path, default=REPO / "dist")
    args = parser.parse_args()
    if sys.platform != "win32":
        raise SystemExit("Build the Windows executable on Windows x64 with Python 3.11")
    destination = inside_repository(args.output_dir)
    target = inside_repository(destination / "COFdraw")
    archive = inside_repository(destination / "COFdraw-Windows.zip")
    if target.exists() or archive.exists():
        raise SystemExit("Output already exists; choose a fresh --output-dir")
    data = args.catalog_data.resolve()
    catalog = validate_catalog(data)
    validate_synthesis_resources()
    run = inside_repository(REPO / "build" / uuid.uuid4().hex)
    run.mkdir(parents=True)
    destination.mkdir(parents=True, exist_ok=True)
    from PIL import Image
    logo = REPO / "src/cof_vision/web/cofdraw-logo.png"
    icon = run / "COFdraw.ico"
    with Image.open(logo) as original:
        rgba = original.convert("RGBA")
        edge = max(rgba.size)
        square = Image.new("RGBA", (edge, edge), (255, 255, 255, 0))
        square.alpha_composite(rgba, ((edge - rgba.width) // 2, (edge - rgba.height) // 2))
        square.save(icon, format="ICO", sizes=[(x, x) for x in (16, 24, 32, 48, 64, 128, 256)])
    validate_raspa_resources()
    command = [sys.executable, "-m", "PyInstaller", "--noconfirm", "--onedir", "--windowed",
               "--name", "COFdraw", "--distpath", str(run / "dist"),
               "--workpath", str(run / "work"), "--specpath", str(run),
               "--paths", str(REPO / "src"), "--icon", str(icon),
               "--add-data", f"{REPO / 'src/cof_vision/web'};cof_vision/web",
               "--collect-all", "rdkit", "--collect-all", "resvg_py",
               "--collect-all", "raspalib", "--copy-metadata", "raspalib",
               "--hidden-import", "uvicorn.logging", "--hidden-import", "uvicorn.loops.asyncio",
               "--hidden-import", "uvicorn.protocols.http.h11_impl", "--hidden-import", "uvicorn.lifespan.on"]
    for module in ("torch", "onnx", "onnxruntime", "pymupdf", "fitz", "matplotlib", "IPython",
                   "tkinter", "pytest", "pandas", "scipy", "sklearn", "notebook", "PyQt5", "PyQt6", "PySide2", "PySide6"):
        command += ["--exclude-module", module]
    command.append(str(REPO / "scripts/launch_cofdraw.py"))
    subprocess.run(command, check=True)
    staged = inside_repository(run / "dist/COFdraw")
    if not (staged / "COFdraw.exe").is_file():
        raise SystemExit("PyInstaller did not produce COFdraw.exe")
    synthesis_resource = validate_synthesis_resources(staged / "_internal/rdkit/Contrib/SA_Score")
    # Windows preserves the Internet-zone marking when Explorer extracts a
    # downloaded ZIP. The process-scoped CLR setting lets Python.NET/WebView2
    # load the supplied managed DLLs without unblocking files or changing the OS.
    shutil.copy2(REPO / "packaging/COFdraw.exe.config", staged / "COFdraw.exe.config")
    portable = staged / "data/framework-drawing"
    portable.mkdir(parents=True)
    # Copy only the hash-bound catalog. Local exports and logs stay private.
    for relative in ["current-bundle-manifest.json", *[row["path"] for row in catalog["files"]]]:
        output = portable / relative
        output.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(data / "framework-drawing" / relative, output)
    validate_catalog(staged / "data")
    crystal_data = data / "crystal-models"
    crystal_manifest = json.loads((crystal_data / "manifest.json").read_text("utf-8"))
    for row in crystal_manifest["templates"]:
        crystal = (crystal_data / row["path"]).resolve()
        if not crystal.is_relative_to(crystal_data.resolve()) or sha(crystal) != row["sha256"]:
            raise SystemExit("Periodic template hash mismatch")
    shutil.copytree(crystal_data, staged / "data/crystal-models")
    # Immutable engine/profile resources belong to the application, independent
    # of a user's writable catalog-data location.
    shutil.copytree(REPO / "data/raspa", staged / "data/raspa")
    for cpu in ("base", "avx", "avx2", "avx512"):
        if not list((staged / "_internal/raspalib").glob(f"raspalib_{cpu}*.pyd")):
            raise SystemExit(f"Bundled RASPA CPU variant missing: {cpu}")
    shutil.copy2(icon, staged / "COFdraw.ico")
    shutil.copy2(REPO / "LICENSE", staged / "LICENSE.txt")
    for name in ("README.md", "LICENSE", "CONTRIBUTING.md", "CHANGELOG.md", "SKILL.md", "AGENTS.md"):
        if (REPO / name).is_file():
            shutil.copy2(REPO / name, staged / name)
    for name in ("docs", "skills", "benchmarks"):
        if (REPO / name).is_dir():
            shutil.copytree(REPO / name, staged / name,
                            ignore=shutil.ignore_patterns("__pycache__", "*.pyc", "*.log"))
    copy_licenses(staged / "third-party-licenses")
    shutil.copy2(REPO / "docs/third-party-licenses/RDKit-SA-Score-BSD.txt",
                 staged / "third-party-licenses/RDKit-SA-Score-BSD.txt")
    # The generated Gemmi coefficient table retains MPL-2.0. Supply its exact
    # editable source beside the executable, as well as the upstream license.
    coefficient_source = staged / "third-party-sources/Gemmi"
    coefficient_source.mkdir(parents=True)
    shutil.copy2(REPO / "src/cof_vision/cofdraw_xray_factors.py", coefficient_source)
    shutil.copy2(REPO / "docs/third-party-licenses/Gemmi-MPL-2.0.txt", coefficient_source)
    manifest = {"application_id": "cofdraw", "version": VERSION,
                "executable_sha256": sha(staged / "COFdraw.exe"),
                "catalog_sha256": catalog["catalog_sha256"],
                "synthetic_accessibility": {"method_id": "rdkit_sa_score", "fragment_table": synthesis_resource},
                "files": [{"path": p.relative_to(staged).as_posix(), "bytes": p.stat().st_size, "sha256": sha(p)}
                          for p in sorted(staged.rglob("*")) if p.is_file()]}
    (staged / "release-manifest.json").write_text(json.dumps(manifest, indent=2), "utf-8")
    # Both checked locations are inside this repository. Keep the build cache
    # and copy into a fresh target; Windows scanners may hold the staging EXE.
    shutil.copytree(staged, target)
    with zipfile.ZipFile(archive, "x", compression=zipfile.ZIP_DEFLATED, compresslevel=6) as output:
        for item in sorted(target.rglob("*")):
            if item.is_file():
                output.write(item, "COFdraw/" + item.relative_to(target).as_posix())
    checksum = f"{sha(archive)}  {archive.name}\n"
    (destination / "SHA256SUMS.txt").write_text(checksum, "utf-8")
    print(json.dumps({"archive": str(archive), "sha256": sha(archive), "bytes": archive.stat().st_size}, indent=2))


if __name__ == "__main__":
    main()
