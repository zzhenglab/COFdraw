"""Run the local COFdraw API and retain real outputs without altering its catalog."""
from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import time
from urllib.error import HTTPError, URLError
from urllib.parse import urlparse
from urllib.request import Request, build_opener, ProxyHandler


def write_json(path: Path, value):
    path.write_text(json.dumps(value, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")


class LocalAPIError(RuntimeError):
    def __init__(self, status, message):
        super().__init__(message)
        self.status = status


class AdsorptionJob:
    """Own one asynchronous job; only a confirmed terminal state permits shutdown."""
    terminal = frozenset({"complete", "failed", "cancelled", "interrupted"})

    def __init__(self, fetch, base, output):
        self.fetch, self.base, self.output = fetch, base, output
        self.identifier = None
        self.status = None
        self.start_pending = False

    @property
    def safe_to_stop(self):
        return not self.start_pending and (self.identifier is None or (self.status or {}).get("status") in self.terminal)

    def _retain(self, value):
        if not isinstance(value, dict) or value.get("job_id") != self.identifier:
            raise RuntimeError("Unexpected adsorption job identity")
        if value.get("status") not in self.terminal | {"queued", "running", "cancelling"}:
            raise RuntimeError("Unexpected adsorption job state")
        self.status = value
        write_json(self.output / "job-status.json", value)
        return value

    def start(self, request):
        self.start_pending = True
        try:
            value = json.loads(self.fetch(self.base, "/api/drawing/science/adsorption/jobs", request))
        except LocalAPIError as exc:
            if exc.status in {400, 409, 422}:
                self.start_pending = False  # Request was rejected before starting a job.
            raise
        identifier = value.get("job_id") if isinstance(value, dict) else None
        if not isinstance(identifier, str) or not re.fullmatch("[0-9a-f]{32}", identifier):
            raise RuntimeError("Missing adsorption job identity; local server retained for inspection")
        self.identifier = identifier
        self.start_pending = False
        self._retain(value)

    def poll(self):
        return self._retain(json.loads(self.fetch(self.base, "/api/drawing/science/adsorption/jobs/" + self.identifier, timeout=10)))

    def wait(self, interval=1):
        last_stage = None
        while self.status["status"] not in self.terminal:
            stage = (self.status.get("status"), self.status.get("stage"))
            if stage != last_stage:
                print(json.dumps({"job_id": self.identifier, "status": stage[0], "stage": stage[1]}), flush=True)
                last_stage = stage
            time.sleep(interval)
            self.poll()
        return self.result()

    def result(self):
        if not self.safe_to_stop or self.status is None:
            raise RuntimeError("Adsorption job is not confirmed finished")
        return self.status.get("result") or self.status

    def cancel_if_active(self, timeout=45, interval=.25):
        if self.start_pending:
            raise RuntimeError("Job creation outcome is unknown; inspect the retained server before stopping it")
        if self.safe_to_stop:
            return self.result() if self.status else None
        # This is our returned ID only; never cancel unrelated jobs or processes.
        self._retain(json.loads(self.fetch(self.base, "/api/drawing/science/adsorption/jobs/" + self.identifier + "/cancel", {}, timeout=10)))
        deadline = time.monotonic() + timeout
        while not self.safe_to_stop:
            if time.monotonic() >= deadline:
                raise RuntimeError("Cancellation is not confirmed; local server retained until its job finishes")
            time.sleep(interval)
            self.poll()
        return self.result()


def save_result(output, result, fetch, base, prefix):
    write_json(output / "result.json", result)
    artifacts = []
    for extension, path in result.get("downloads", {}).items():
        if extension not in {"svg", "png", "mol", "sdf", "cdxml", "json", "cif", "xyz", "pdb", "csv", "zip"}:
            raise ValueError(f"Unexpected export format: {extension}")
        content = fetch(base, path)
        filename = f"{prefix}.{extension}"
        (output / filename).write_bytes(content)
        artifacts.append({"file": filename, "bytes": len(content), "sha256": hashlib.sha256(content).hexdigest()})
    write_json(output / "artifacts.json", artifacts)
    return artifacts


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    runtime = parser.add_mutually_exclusive_group(required=True)
    runtime.add_argument("--repo", type=Path, help="COFdraw source checkout")
    runtime.add_argument("--exe", type=Path, help="Extracted Windows COFdraw.exe")
    parser.add_argument("--data-dir", type=Path, help="Reviewed catalog template; never modified")
    parser.add_argument("--output", type=Path, required=True, help="New directory for this invocation")
    commands = parser.add_subparsers(dest="operation", required=True)
    commands.add_parser("library", help="Save the actual available catalog")
    commands.add_parser("describe-catalog", help="Read description examples and selectable catalog IDs")
    describe = commands.add_parser("describe", help="Interpret a description into a reviewable plan; does not draw")
    describe.add_argument("--request", type=Path, required=True, help="JSON /describe body; rules mode is the default")
    source = commands.add_parser("source", help="Render a registered source recipe")
    source.add_argument("recipe_id")
    design = commands.add_parser("design", help="Prepare monomers and render the resulting design")
    design.add_argument("--request", type=Path, required=True, help="JSON /prepare body")
    render = commands.add_parser("render", help="Render a complete JSON /render request")
    render.add_argument("--request", type=Path, required=True)
    layout = commands.add_parser("layout", help="Request an explicit layout revision; saves the returned recipe")
    layout.add_argument("--request", type=Path, required=True)
    science = commands.add_parser("science", help="Run a documented local structure/analysis endpoint")
    science.add_argument("endpoint", choices=["library", "capabilities", "query", "model", "patch", "networks", "monomer", "solubility", "porosity", "geometry", "synthesis-library", "adsorption-engine", "adsorption", "bet", "measure", "inspect", "pxrd", "screen", "finite-inspect", "finite-powder", "stacking", "periodic-options", "build-periodic"])
    science.add_argument("--request", type=Path)
    for sub in (source, design):
        sub.add_argument("--style", default="standard")
        sub.add_argument("--size", type=int, default=1200, choices=range(500, 2401), metavar="500..2400")
        sub.add_argument("--seed", type=int)
    source.set_defaults(style="source")
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
    science_get = {"library": "library", "capabilities": "capabilities", "synthesis-library": "synthesis-library", "adsorption-engine": "adsorption/engine"}
    if args.operation == "science" and args.endpoint not in science_get and args.request is None:
        parser.error("A science calculation requires --request with its JSON request body")
    if args.operation == "science" and args.endpoint in science_get and args.request is not None:
        parser.error("This read-only science endpoint takes no request body")
    if args.operation in {"design", "render", "layout", "describe"} or (args.operation == "science" and args.request):
        request_body = json.loads(args.request.read_text(encoding="utf-8-sig"))
        if not isinstance(request_body, dict):
            parser.error("Request JSON must be an object")

    output.mkdir(parents=True)
    data = output / "runtime-data"
    shutil.copytree(template / "framework-drawing", data / "framework-drawing",
                    ignore=shutil.ignore_patterns("exports", "*.log", "__pycache__"))
    if (template / "crystal-models").is_dir():
        shutil.copytree(template / "crystal-models", data / "crystal-models")
    status_file = output / "server.json"
    command += ["--headless", "--data-dir", str(data), "--status-file", str(status_file)]
    opener = build_opener(ProxyHandler({}))  # Loopback calls must not use an external proxy.
    process = None
    adsorption_job = None

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
            if path == "/api/drawing/describe":
                # Even a framework validation response could echo a malformed
                # credential input. Do not read or print that response body.
                code = exc.code
                exc.close()
                raise LocalAPIError(code, f"COFdraw description request failed (HTTP {code}); check the settings and supported inputs.") from None
            detail = exc.read().decode("utf-8", errors="replace")
            raise LocalAPIError(exc.code, f"COFdraw HTTP {exc.code}: {detail}") from exc

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
            if args.operation == "science":
                if args.endpoint == "adsorption":
                    adsorption_job = AdsorptionJob(fetch, base, output)
                    adsorption_job.start(request_body)
                    result = adsorption_job.wait()
                else:
                    result = json.loads(fetch(base, "/api/drawing/science/" + science_get.get(args.endpoint, args.endpoint), request_body))
            elif args.operation == "layout":
                result = json.loads(fetch(base, "/api/drawing/layout", request_body))
            elif args.operation == "describe-catalog":
                result = json.loads(fetch(base, "/api/drawing/describe/catalog"))
            elif args.operation == "describe":
                result = json.loads(fetch(base, "/api/drawing/describe", request_body))
                # Never persist the request or credentials. The server returns
                # a transient, catalog-bound plan, not a raw provider response.
            if args.operation == "source":
                body = {"source_recipe_id": args.recipe_id}
            elif args.operation == "design":
                prepared = json.loads(fetch(base, "/api/drawing/prepare", request_body))
                write_json(output / "prepared.json", prepared)
                body = {"recipe": prepared["recipe"]}
            elif args.operation == "render":
                body = request_body
            if args.operation in {"source", "design"}:
                body.update(style_id=args.style, size=args.size)
                if args.seed is not None:
                    body["style_seed"] = args.seed
            if args.operation in {"source", "design", "render"}:
                result = json.loads(fetch(base, "/api/drawing/render", body))
            save_result(output, result, fetch, base, "structure" if args.operation == "science" else "drawing")
            print(json.dumps({"output": str(output), "status": result.get("status", result.get("scope", args.operation)),
                              "formats": list(result.get("downloads", {}))}))
            if adsorption_job and result.get("status") != "complete":
                raise RuntimeError("Adsorption job did not complete successfully; retained its status and available partial exports")
    finally:
        safe_to_stop = True
        if adsorption_job is not None:
            try:
                was_active = not adsorption_job.safe_to_stop
                terminal_result = adsorption_job.cancel_if_active()
                if was_active and terminal_result:
                    save_result(output, terminal_result, fetch, base, "structure")
            except BaseException as exc:
                # A second interrupt or transport failure cannot justify
                # killing a server whose native job may still be running.
                safe_to_stop = adsorption_job.safe_to_stop
                print(f"Adsorption cleanup: {exc}. See job-status.json and server.json in {output}.", file=sys.stderr)
                if not safe_to_stop:
                    print("The owned local server was left running for job completion/inspection; no other process was stopped.", file=sys.stderr)
        if safe_to_stop and process is not None and process.poll() is None:
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
