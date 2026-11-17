"""Run the local COFdraw API and retain real outputs without altering its catalog."""
from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time
from urllib.error import HTTPError, URLError
from urllib.parse import urlparse
from urllib.request import Request, build_opener, ProxyHandler


def write_json(path: Path, value):
    path.write_text(json.dumps(value, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    runtime = parser.add_mutually_exclusive_group(required=True)
    runtime.add_argument("--repo", type=Path, help="COFdraw source checkout")
    runtime.add_argument("--exe", type=Path, help="Extracted Windows COFdraw.exe")
    parser.add_argument("--data-dir", type=Path, help="Reviewed catalog template; never modified")
    parser.add_argument("--output", type=Path, required=True, help="New directory for this invocation")
    commands = parser.add_subparsers(dest="operation", required=True)
    commands.add_parser("library", help="Save the actual available catalog")
    source = commands.add_parser("source", help="Render a registered source recipe")
    source.add_argument("recipe_id")
    design = commands.add_parser("design", help="Prepare monomers and render the resulting design")
    design.add_argument("--request", type=Path, required=True, help="JSON /prepare body")
    render = commands.add_parser("render", help="Render a complete JSON /render request")
    render.add_argument("--request", type=Path, required=True)
    for sub in (source, design):
        sub.add_argument("--style", default="standard")
        sub.add_argument("--size", type=int, default=1200, choices=range(500, 2401), metavar="500..2400")
        sub.add_argument("--seed", type=int)
    args = parser.parse_args()

    if args.repo:
        root = args.repo.resolve()
        launcher = root / "scripts/launch_cofdraw.py"
        command = [sys.executable, str(launcher)]
    else:
        launcher = args.exe.resolve()
        root = launcher.parent
        command = [str(launcher)]
    template = (args.data_dir or root / "data").resolve()
    if not launcher.is_file():
        parser.error(f"Launcher not found: {launcher}")
    if not (template / "framework-drawing/registry.json").is_file():
        parser.error(f"Portable catalog registry not found under {template}")
    output = args.output.resolve()
    if output.exists():
        parser.error("Output already exists; choose a new directory to preserve previous work")
    if output == template or template in output.parents:
        parser.error("Output must be outside the catalog template")
    request_body = None
    if args.operation in {"design", "render"}:
        request_body = json.loads(args.request.read_text(encoding="utf-8-sig"))
        if not isinstance(request_body, dict):
            parser.error("Request JSON must be an object")

    output.mkdir(parents=True)
    data = output / "runtime-data"
    shutil.copytree(template / "framework-drawing", data / "framework-drawing",
                    ignore=shutil.ignore_patterns("exports", "*.log", "__pycache__"))
    status_file = output / "server.json"
    command += ["--headless", "--data-dir", str(data), "--status-file", str(status_file)]
    opener = build_opener(ProxyHandler({}))  # Loopback calls must not use an external proxy.
    process = None

    def fetch(base, path, value=None, timeout=120):
        if not path.startswith("/api/") or path.startswith("//"):
            raise ValueError("COFdraw returned a nonlocal API path")
        payload = None if value is None else json.dumps(value).encode("utf-8")
        request = Request(base + path, data=payload, headers={"Content-Type": "application/json"})
        try:
            with opener.open(request, timeout=timeout) as response:
                if not response.url.startswith(base + "/"):
                    raise ValueError("Unexpected API redirect")
                return response.read()
        except HTTPError as exc:
            detail = exc.read().decode("utf-8", errors="replace")
            raise RuntimeError(f"COFdraw HTTP {exc.code}: {detail}") from exc

    try:
        with (output / "server-output.log").open("wb") as log:
            process = subprocess.Popen(command, stdout=log, stderr=subprocess.STDOUT,
                                       creationflags=getattr(subprocess, "CREATE_NO_WINDOW", 0))
            deadline = time.monotonic() + 45
            base = None
            while time.monotonic() < deadline:
                if process.poll() is not None:
                    raise RuntimeError(f"COFdraw exited with code {process.returncode}; see server-output.log")
                try:
                    status = json.loads(status_file.read_text(encoding="utf-8"))
                    candidate = status["url"].rstrip("/")
                    parsed = urlparse(candidate)
                    # Windows venv launchers can spawn a child Python process,
                    # so its server PID can differ from Popen.pid. This fresh
                    # private status file and copied data root identify our run.
                    if (not isinstance(status.get("pid"), int)
                            or Path(status["data_directory"]).resolve() != data.resolve()
                            or parsed.scheme != "http" or parsed.hostname != "127.0.0.1"):
                        raise ValueError("Unexpected local service identity")
                    about = json.loads(fetch(candidate, "/api/desktop", timeout=2))
                    if about["offline"]:
                        base = candidate
                        break
                except (OSError, ValueError, KeyError, URLError):
                    time.sleep(0.2)
            if base is None:
                raise RuntimeError("Local COFdraw service did not become ready within 45 seconds")
            library = json.loads(fetch(base, "/api/drawing/library"))
            write_json(output / "library.json", library)
            if library.get("source_registry_errors"):
                raise RuntimeError("Catalog errors: " + "; ".join(library["source_registry_errors"]))
            if args.operation == "library":
                print(json.dumps({"library": str(output / "library.json"), "source_recipes": len(library["source_recipes"])}))
                return
            if args.operation == "source":
                body = {"source_recipe_id": args.recipe_id}
            elif args.operation == "design":
                prepared = json.loads(fetch(base, "/api/drawing/prepare", request_body))
                write_json(output / "prepared.json", prepared)
                body = {"recipe": prepared["recipe"]}
            else:
                body = request_body
            if args.operation in {"source", "design"}:
                body.update(style_id=args.style, size=args.size)
                if args.seed is not None:
                    body["style_seed"] = args.seed
            result = json.loads(fetch(base, "/api/drawing/render", body))
            write_json(output / "result.json", result)
            artifacts = []
            for extension, path in result["downloads"].items():
                if extension not in {"svg", "png", "mol", "sdf", "cdxml", "json"}:
                    raise ValueError(f"Unexpected export format: {extension}")
                content = fetch(base, path)
                filename = f"drawing.{extension}"
                (output / filename).write_bytes(content)
                artifacts.append({"file": filename, "bytes": len(content), "sha256": hashlib.sha256(content).hexdigest()})
            write_json(output / "artifacts.json", artifacts)
            print(json.dumps({"output": str(output), "status": result["status"], "formats": list(result["downloads"])}))
    finally:
        if process is not None and process.poll() is None:
            if os.name == "nt":
                # Include the child interpreter of a Windows venv launcher.
                subprocess.run(["taskkill", "/PID", str(process.pid), "/T", "/F"],
                               capture_output=True, timeout=10,
                               creationflags=getattr(subprocess, "CREATE_NO_WINDOW", 0))
            else:
                process.terminate()
            try:
                process.wait(timeout=10)
            except subprocess.TimeoutExpired:
                process.kill()
                process.wait(timeout=10)


if __name__ == "__main__":
    try:
        main()
    except (OSError, ValueError, RuntimeError) as error:
        print(f"COFdraw: {error}", file=sys.stderr)
        raise SystemExit(1)
