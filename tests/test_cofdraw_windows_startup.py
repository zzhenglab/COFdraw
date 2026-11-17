"""Startup failure handling and the packaged CLR sidecar contract."""
import importlib.util
import json
from pathlib import Path
import socket
import sys
from types import SimpleNamespace
from urllib.parse import urlsplit
import zipfile

import pytest

from cof_vision import cofdraw_desktop


REPO = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("cofdraw_download_smoke", REPO / "scripts/smoke_cofdraw_download.py")
download_smoke = importlib.util.module_from_spec(spec)
spec.loader.exec_module(download_smoke)


def _assert_local_service_stopped(status):
    address = urlsplit(status["url"])
    with socket.socket() as probe:
        probe.settimeout(1)
        assert probe.connect_ex((address.hostname, address.port)) != 0


def test_native_loader_failure_returns_error_and_stops_real_server(monkeypatch, tmp_path):
    def fail(*args, **kwargs):
        raise RuntimeError("Cannot initialize the test native runtime")
    fake = SimpleNamespace(create_window=lambda *a, **kw: SimpleNamespace(), start=fail)
    monkeypatch.setitem(sys.modules, "webview", fake)
    status_file = tmp_path / "status.json"
    result = cofdraw_desktop.main(["--smoke-window", "--data-dir", str(tmp_path / "data"),
                                   "--status-file", str(status_file)])
    status = json.loads(status_file.read_text("utf-8"))
    assert result == 1
    assert status["status"] == "failed"
    assert "test native runtime" in status["error"]
    _assert_local_service_stopped(status)


def test_failed_native_page_smoke_returns_error_and_destroys_window(monkeypatch, tmp_path):
    destroyed = []
    def fail_script(script):
        raise RuntimeError("Native page did not load")
    window = SimpleNamespace(evaluate_js=fail_script, destroy=lambda: destroyed.append(True))
    fake = SimpleNamespace(create_window=lambda *a, **kw: window, start=lambda callback, **kw: callback())
    monkeypatch.setitem(sys.modules, "webview", fake)
    status_file = tmp_path / "status.json"
    result = cofdraw_desktop.main(["--smoke-window", "--data-dir", str(tmp_path / "data"),
                                   "--status-file", str(status_file)])
    status = json.loads(status_file.read_text("utf-8"))
    assert result == 1
    assert status["status"] == "desktop_smoke_failed"
    assert destroyed == [True]
    _assert_local_service_stopped(status)


def test_sidecar_has_only_app_runtime_setting_and_is_manifest_bound(tmp_path):
    config = tmp_path / "COFdraw.exe.config"
    config.write_bytes((REPO / "packaging/COFdraw.exe.config").read_bytes())
    digest = download_smoke.sha256(config)
    manifest = tmp_path / "release-manifest.json"
    manifest.write_text(json.dumps({"files": [{"path": config.name, "sha256": digest}]}), "utf-8")
    assert download_smoke.validate_packaged_config(tmp_path) == digest
    manifest.write_text(json.dumps({"files": []}), "utf-8")
    with pytest.raises(ValueError, match="manifest"):
        download_smoke.validate_packaged_config(tmp_path)


def test_mismatched_config_hash_is_not_validated(tmp_path):
    config = tmp_path / "COFdraw.exe.config"
    config.write_bytes((REPO / "packaging/COFdraw.exe.config").read_bytes())
    (tmp_path / "release-manifest.json").write_text(json.dumps({"files": [
        {"path": config.name, "sha256": "0" * 64}]}), "utf-8")
    with pytest.raises(ValueError, match="manifest"):
        download_smoke.validate_packaged_config(tmp_path)


@pytest.mark.parametrize("entry", ["../outside", "COFdraw/../../outside", "/COFdraw/file", "COFdraw/C:/file", "another/app.exe"])
def test_download_extraction_rejects_unsafe_paths(tmp_path, entry):
    archive = tmp_path / "bad.zip"
    with zipfile.ZipFile(archive, "w") as bundle:
        bundle.writestr(entry, "data")
    with pytest.raises(ValueError, match="unsafe"):
        download_smoke.extract_fresh(archive, tmp_path / "fresh extraction")
    assert not (tmp_path / "fresh extraction").exists()


def test_download_extraction_preserves_existing_output(tmp_path):
    output = tmp_path / "existing"
    output.mkdir()
    sentinel = output / "keep.txt"
    sentinel.write_text("unchanged", "utf-8")
    with pytest.raises(ValueError, match="new"):
        download_smoke.extract_fresh(tmp_path / "unused.zip", output)
    assert sentinel.read_text("utf-8") == "unchanged"


def _marked_fixture(tmp_path):
    records = []
    for name in ("COFdraw.exe", "COFdraw.exe.config", "Python.Runtime.dll"):
        path = tmp_path / name
        path.write_bytes(b"unchanged application contents")
        Path(str(path) + ":Zone.Identifier").write_bytes(b"[ZoneTransfer]\r\nZoneId=3\r\n")
        records.append({"path": name, "sha256": download_smoke.sha256(path)})
    return records


def test_after_launch_marks_reports_executable_absence_without_claiming_all_retained(tmp_path):
    records = _marked_fixture(tmp_path)
    Path(str(tmp_path / "COFdraw.exe") + ":Zone.Identifier").unlink()
    result = download_smoke.verify_after_launch_marks(tmp_path, records)
    assert result["marked_file_hashes_unchanged"] is True
    assert result["dll_config_download_marks_retained"] is True
    assert result["executable_download_mark_after"] is False
    assert result["download_marks_retained"] is False
    with pytest.raises(FileNotFoundError):
        download_smoke.verify_marks(tmp_path, records)


def test_after_launch_marks_reports_retained_executable(tmp_path):
    result = download_smoke.verify_after_launch_marks(tmp_path, _marked_fixture(tmp_path))
    assert result["executable_download_mark_after"] is True
    assert result["download_marks_retained"] is True


@pytest.mark.parametrize("name", ["COFdraw.exe.config", "Python.Runtime.dll"])
def test_after_launch_missing_library_or_config_mark_always_fails(tmp_path, name):
    records = _marked_fixture(tmp_path)
    Path(str(tmp_path / name) + ":Zone.Identifier").unlink()
    with pytest.raises(AssertionError, match="DLL or configuration"):
        download_smoke.verify_after_launch_marks(tmp_path, records)


def test_after_launch_changed_file_always_fails(tmp_path):
    records = _marked_fixture(tmp_path)
    (tmp_path / "COFdraw.exe").write_bytes(b"changed")
    with pytest.raises(AssertionError, match="file changed"):
        download_smoke.verify_after_launch_marks(tmp_path, records)


def test_after_launch_changed_zone_always_fails(tmp_path):
    records = _marked_fixture(tmp_path)
    Path(str(tmp_path / "COFdraw.exe") + ":Zone.Identifier").write_bytes(b"ZoneId=2\r\n")
    with pytest.raises(AssertionError, match="Internet zone"):
        download_smoke.verify_after_launch_marks(tmp_path, records)


def test_explorer_launch_uses_transparent_shortcut_and_plain_helper(monkeypatch, tmp_path):
    executable = tmp_path / "COFdraw.exe"
    executable.write_bytes(b"test fixture, not executable")
    output = tmp_path / "launch"
    calls = []
    def create_shortcut(command, **kwargs):
        assert "-Command" in command and "-EncodedCommand" not in command
        assert "-File" not in command and "-ExecutionPolicy" not in command
        helper = output / "create-startup-shortcut.ps1"
        script = helper.read_text("utf-8-sig")
        assert command[-1] == script
        assert "Invoke-Expression" not in script and "Set-ExecutionPolicy" not in script
        assert "WScript.Shell" in script and ".CreateShortcut(" in script
        assert "Shell.Application" not in script and "FindWindowSW" not in script
        assert ".TargetPath=" in script and ".WorkingDirectory=" in script
        assert ".Arguments=" in script and "--smoke-window" in script
        (output / "COFdraw startup check.lnk").write_bytes(b"mock shell link")
        calls.append(command)
        return SimpleNamespace(returncode=0, stdout=b"", stderr=b"")
    def open_explorer(command, **kwargs):
        calls.append(command)
        assert Path(command[0]).name == "explorer.exe"
        assert command[1] == str(output / "COFdraw startup check.lnk")
        return SimpleNamespace(pid=123)
    monkeypatch.setattr(download_smoke, "_windows_api", lambda: None)
    monkeypatch.setattr(download_smoke.subprocess, "CREATE_NO_WINDOW", 0x08000000, raising=False)
    monkeypatch.setattr(download_smoke.subprocess, "run", create_shortcut)
    monkeypatch.setattr(download_smoke.subprocess, "Popen", open_explorer)
    result = download_smoke.launch_from_explorer(executable, ["--smoke-window"], output)
    assert result.pid == 123 and len(calls) == 2
