"""Test a fresh Windows ZIP with download marks and a normal-user native launch.

This requires DLL/configuration download marks and reports the EXE mark separately.
It never unblocks files, changes machine policy, or modifies the supplied ZIP. Requires Windows,
NTFS, and an interactive Explorer desktop running at medium integrity.
"""
from __future__ import annotations

import argparse
import ctypes
from ctypes import wintypes
import hashlib
import json
import os
from pathlib import Path, PurePosixPath
import subprocess
import sys
import time
from urllib.error import URLError
from urllib.request import urlopen
import zipfile
import xml.etree.ElementTree as ET


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def extract_fresh(archive, destination):
    destination = Path(destination).resolve()
    if destination.exists():
        raise ValueError("Extraction directory must be new")
    with zipfile.ZipFile(archive) as bundle:
        seen = set()
        for entry in bundle.infolist():
            name = PurePosixPath(entry.filename.replace("\\", "/"))
            if (name.is_absolute() or ".." in name.parts or not name.parts
                    or name.parts[0] != "COFdraw" or any(":" in x for x in name.parts)
                    or (entry.external_attr >> 16) & 0o170000 == 0o120000):
                raise ValueError("Archive contains an unsafe path")
            normalized = str(name).casefold()
            if normalized in seen:
                raise ValueError("Archive contains duplicate paths")
            seen.add(normalized)
        destination.mkdir(parents=True)
        bundle.extractall(destination)
    return destination / "COFdraw"


def mark_downloads(application):
    application = Path(application).resolve()
    executable = application / "COFdraw.exe"
    config = application / "COFdraw.exe.config"
    managed = application / "_internal/pythonnet/runtime/Python.Runtime.dll"
    core = application / "_internal/webview/lib/Microsoft.Web.WebView2.Core.dll"
    forms = application / "_internal/webview/lib/Microsoft.Web.WebView2.WinForms.dll"
    loader = application / "_internal/webview/lib/runtimes/win-x64/native/WebView2Loader.dll"
    for required in (executable, config, managed, core, forms, loader):
        if not required.is_file():
            raise ValueError("Release is missing a required Windows desktop file: " + required.name)
    paths = sorted({executable, config, *application.rglob("*.dll")})
    records = []
    zone = b"[ZoneTransfer]\r\nZoneId=3\r\nHostUrl=https://github.com/zzhenglab/COFdraw/releases\r\n"
    for path in paths:
        if not path.resolve().is_relative_to(application):
            raise ValueError("Download-mark target leaves the extracted application")
        before = sha256(path)
        with open(str(path) + ":Zone.Identifier", "wb") as stream:
            stream.write(zone)
        records.append({"path": path.relative_to(application).as_posix(), "sha256": before})
    verify_marks(application, records)
    return records


def validate_packaged_config(application):
    application = Path(application)
    config = application / "COFdraw.exe.config"
    tree = ET.fromstring(config.read_bytes())
    if (tree.tag != "configuration" or tree.attrib or len(tree) != 1
            or tree[0].tag != "runtime" or tree[0].attrib or len(tree[0]) != 1
            or tree[0][0].tag != "loadFromRemoteSources"
            or tree[0][0].attrib != {"enabled": "true"} or len(tree[0][0])):
        raise ValueError("Unexpected executable-scoped CLR configuration")
    manifest = json.loads((application / "release-manifest.json").read_text("utf-8"))
    matches = [row for row in manifest["files"] if row["path"] == "COFdraw.exe.config"]
    if len(matches) != 1 or matches[0]["sha256"] != sha256(config):
        raise ValueError("Release manifest does not bind the executable configuration")
    return matches[0]["sha256"]


def verify_marks(application, records):
    application = Path(application).resolve()
    for record in records:
        path = (application / record["path"]).resolve()
        if not path.is_relative_to(application) or sha256(path) != record["sha256"]:
            raise AssertionError("A marked application file changed")
        with open(str(path) + ":Zone.Identifier", "rb") as stream:
            if b"ZoneId=3" not in stream.read(8192).splitlines():
                raise AssertionError("Internet download mark was removed")


def verify_after_launch_marks(application, records):
    """Require unchanged contents and all DLL/config marks; report EXE separately.

    Windows may remove the EXE mark after its ordinary first-run confirmation.
    This function does not remove or restore any mark and never excuses a
    missing DLL or configuration mark.
    """
    application = Path(application).resolve()
    executable_after = None
    for record in records:
        path = (application / record["path"]).resolve()
        if not path.is_relative_to(application) or sha256(path) != record["sha256"]:
            raise AssertionError("A marked application file changed")
        try:
            with open(str(path) + ":Zone.Identifier", "rb") as stream:
                if b"ZoneId=3" not in stream.read(8192).splitlines():
                    raise AssertionError("A file no longer has the expected Internet zone")
        except FileNotFoundError:
            if record["path"] != "COFdraw.exe":
                raise AssertionError("A DLL or configuration download mark was removed") from None
            executable_after = False
        else:
            if record["path"] == "COFdraw.exe":
                executable_after = True
    if executable_after is None:
        raise AssertionError("Executable mark evidence is missing")
    return {"marked_file_hashes_unchanged": True,
            "download_marks_retained": executable_after,
            "dll_config_download_marks_retained": True,
            "executable_download_mark_after": executable_after}


def _windows_api():
    if sys.platform != "win32":
        raise RuntimeError("The marked-download smoke test requires Windows")
    kernel = ctypes.WinDLL("kernel32", use_last_error=True)
    security = ctypes.WinDLL("advapi32", use_last_error=True)
    kernel.OpenProcess.argtypes = [wintypes.DWORD, wintypes.BOOL, wintypes.DWORD]
    kernel.OpenProcess.restype = wintypes.HANDLE
    kernel.CloseHandle.argtypes = [wintypes.HANDLE]
    kernel.WaitForSingleObject.argtypes = [wintypes.HANDLE, wintypes.DWORD]
    kernel.WaitForSingleObject.restype = wintypes.DWORD
    kernel.TerminateProcess.argtypes = [wintypes.HANDLE, wintypes.UINT]
    kernel.QueryFullProcessImageNameW.argtypes = [wintypes.HANDLE, wintypes.DWORD, wintypes.LPWSTR, ctypes.POINTER(wintypes.DWORD)]
    security.OpenProcessToken.argtypes = [wintypes.HANDLE, wintypes.DWORD, ctypes.POINTER(wintypes.HANDLE)]
    security.GetTokenInformation.argtypes = [wintypes.HANDLE, ctypes.c_int, wintypes.LPVOID, wintypes.DWORD, ctypes.POINTER(wintypes.DWORD)]
    security.GetSidSubAuthorityCount.argtypes = [wintypes.LPVOID]
    security.GetSidSubAuthorityCount.restype = ctypes.POINTER(ctypes.c_ubyte)
    security.GetSidSubAuthority.argtypes = [wintypes.LPVOID, wintypes.DWORD]
    security.GetSidSubAuthority.restype = ctypes.POINTER(wintypes.DWORD)
    return kernel, security


def _process_evidence(pid, executable):
    kernel, security = _windows_api()
    process = kernel.OpenProcess(0x1000 | 0x100000 | 0x0001, False, pid)
    if not process:
        raise OSError("Cannot inspect the launched process; normal-user launch is unverified")
    token = wintypes.HANDLE()
    try:
        image = ctypes.create_unicode_buffer(32768)
        size = wintypes.DWORD(len(image))
        if not kernel.QueryFullProcessImageNameW(process, 0, image, ctypes.byref(size)):
            raise OSError("Cannot verify the launched executable")
        if Path(image.value).resolve() != Path(executable).resolve():
            raise AssertionError("Status PID does not belong to the tested executable")
        if not security.OpenProcessToken(process, 8, ctypes.byref(token)):
            raise OSError("Cannot inspect the launched process token")
        needed = wintypes.DWORD()
        security.GetTokenInformation(token, 25, None, 0, ctypes.byref(needed))
        data = ctypes.create_string_buffer(needed.value)
        if not security.GetTokenInformation(token, 25, data, needed.value, ctypes.byref(needed)):
            raise OSError("Cannot read the launched process integrity level")
        class SidAndAttributes(ctypes.Structure):
            _fields_ = [("sid", wintypes.LPVOID), ("attributes", wintypes.DWORD)]
        sid = ctypes.cast(data, ctypes.POINTER(SidAndAttributes)).contents.sid
        count = security.GetSidSubAuthorityCount(sid)[0]
        integrity = security.GetSidSubAuthority(sid, count - 1)[0]
        return process, {"pid": pid, "integrity_rid": integrity, "medium_integrity": integrity == 8192}
    except BaseException:
        kernel.CloseHandle(process)
        raise
    finally:
        if token:
            kernel.CloseHandle(token)


def _powershell_literal(value):
    return "'" + str(value).replace("'", "''") + "'"


def launch_from_explorer(executable, arguments, output):
    """Create a transparent Shell Link and ask ordinary Explorer to open it.

    Returns the Explorer request's Popen object, not the native application's
    process. Callers must independently verify the actual application's token.
    `arguments` is an explicit list of strings, not shell command text.
    """
    _windows_api()
    executable = Path(executable).resolve()
    output = Path(output).resolve()
    if not executable.is_file() or not isinstance(arguments, (list, tuple)):
        raise ValueError("Expected an existing executable and an argument list")
    if any(not isinstance(arg, str) or "\0" in arg for arg in arguments):
        raise ValueError("Shortcut arguments must be strings without null characters")
    output.mkdir(parents=True, exist_ok=True)
    shortcut = output / "COFdraw startup check.lnk"
    helper = output / "create-startup-shortcut.ps1"
    if shortcut.exists() or helper.exists():
        raise ValueError("Choose a fresh Explorer launch output directory")
    script = ("$ErrorActionPreference='Stop'\n"
              "$shellLinks=New-Object -ComObject WScript.Shell\n"
              "$link=$shellLinks.CreateShortcut(" + _powershell_literal(shortcut) + ")\n"
              "$link.TargetPath=" + _powershell_literal(executable) + "\n"
              "$link.Arguments=" + _powershell_literal(subprocess.list2cmdline(arguments)) + "\n"
              "$link.WorkingDirectory=" + _powershell_literal(executable.parent) + "\n"
              "$link.Description='Launch the reviewed COFdraw application'\n"
              "$link.WindowStyle=1\n$link.Save()\n")
    helper.write_text(script, "utf-8-sig")
    system = Path(os.environ.get("SystemRoot", r"C:\Windows"))
    powershell = system / "System32/WindowsPowerShell/v1.0/powershell.exe"
    # Default Windows execution policy can disable .ps1 files while permitting
    # ordinary interactive commands. Keep their literal text for review and
    # pass those same readable commands directly; do not alter policy, encode
    # content, or load/evaluate a script file.
    completed = subprocess.run([str(powershell), "-NoProfile", "-NonInteractive", "-Command", script],
                               capture_output=True, timeout=30, creationflags=subprocess.CREATE_NO_WINDOW)
    if completed.returncode != 0 or not shortcut.is_file():
        (output / "shortcut-helper-error.txt").write_bytes(completed.stdout + completed.stderr)
        raise RuntimeError("Windows shortcut creation failed; see the private helper error file")
    return subprocess.Popen([str(system / "explorer.exe"), str(shortcut)],
                            creationflags=subprocess.CREATE_NO_WINDOW)


def launch_medium_smoke(executable, output, *, expected_status="desktop_smoke_passed", timeout=180):
    """Launch an existing isolated copy; do not modify its DLLs or configuration."""
    kernel, _ = _windows_api()
    executable = Path(executable).resolve()
    output = Path(output).resolve()
    output.mkdir(parents=True, exist_ok=False)
    status_file = output / "native-status.json"
    arguments = ["--smoke-window", "--status-file", str(status_file)]
    report = {"method": "ordinary_explorer_reviewed_shell_link", "expected_status": expected_status,
              "native": None, "process": None, "server_closed": False}
    process = None
    helper = None
    status = None
    try:
        print("Opening the verified COFdraw shortcut in Explorer. If Windows shows a first-run "
              "confirmation, review the exact COFdraw.exe identity before choosing Run. "
              "Do not change security settings or unblock DLLs.", flush=True)
        helper = launch_from_explorer(executable, arguments, output)
        report["explorer_request_pid"] = helper.pid
        deadline = time.monotonic() + timeout
        while time.monotonic() < deadline:
            try:
                status = json.loads(status_file.read_text("utf-8"))
            except (OSError, ValueError):
                status = None
            if status and process is None:
                process, evidence = _process_evidence(int(status["pid"]), executable)
                report["process"] = evidence
                if not evidence["medium_integrity"]:
                    raise AssertionError("Explorer did not launch COFdraw at medium integrity (8192)")
            if process and status and status.get("status") == expected_status == "failed":
                # Old releases rethrow into PyInstaller's windowed exception
                # dialog. Capture the negative evidence, then close only this
                # owned test process instead of leaving a modal dialog open.
                if kernel.WaitForSingleObject(process, 0) != 0:
                    kernel.TerminateProcess(process, 1)
                    kernel.WaitForSingleObject(process, 10000)
                    report["negative_process_stopped_after_failure"] = True
                break
            if process and kernel.WaitForSingleObject(process, 0) == 0:
                break
            time.sleep(.05)
        else:
            raise TimeoutError("Native marked-download launch did not complete")
        status = json.loads(status_file.read_text("utf-8"))
        report["native"] = status
        if status["status"] != expected_status:
            raise AssertionError("Unexpected native startup status: " + str(status["status"]))
        if expected_status == "desktop_smoke_passed" and status.get("native_save_bridge") is not True:
            raise AssertionError("Native save bridge did not become ready")
        try:
            with urlopen(status["url"] + "api/desktop", timeout=2):
                raise AssertionError("Local server remained available after the native window closed")
        except URLError:
            report["server_closed"] = True
        report["passed"] = True
        return report
    except BaseException as exc:
        report["passed"] = False
        report["failure"] = type(exc).__name__ + ": " + str(exc)
        if status:
            report["native"] = status
        raise
    finally:
        if process:
            if kernel.WaitForSingleObject(process, 0) != 0:
                kernel.TerminateProcess(process, 1)
                kernel.WaitForSingleObject(process, 10000)
            kernel.CloseHandle(process)
        # Explorer can hand off to an existing desktop or remain alive. Never
        # terminate Explorer as test cleanup; only the verified test PID above.
        (output / "native-launch-report.json").write_text(json.dumps(report, indent=2) + "\n", "utf-8")


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--zip", dest="archive", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args(argv)
    _windows_api()
    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=False)
    report = {"archive_sha256": sha256(args.archive), "passed": False}
    try:
        application = extract_fresh(args.archive, output / "Downloaded COFdraw with spaces")
        report["runtime_config_sha256"] = validate_packaged_config(application)
        records = mark_downloads(application)
        report["marked_files"] = records
        report["executable_download_mark_before"] = True
        report["launch"] = launch_medium_smoke(application / "COFdraw.exe", output / "normal-user-native")
        report.update(verify_after_launch_marks(application, records))
        report["passed"] = True
    except BaseException as exc:
        report["failure"] = type(exc).__name__ + ": " + str(exc)
        raise
    finally:
        (output / "download-smoke-report.json").write_text(json.dumps(report, indent=2) + "\n", "utf-8")
        print(json.dumps({"passed": report["passed"], "report": str(output / "download-smoke-report.json")}))


if __name__ == "__main__":
    main()
