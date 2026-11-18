"""The standalone agent helper must await/cancel only its own native job."""
from copy import deepcopy
import hashlib
import importlib.util
import json
from pathlib import Path

import pytest

PATH = Path(__file__).resolve().parents[1] / "skills/cofdraw/scripts/cofdraw_client.py"
spec = importlib.util.spec_from_file_location("cofdraw_skill_client", PATH)
helper = importlib.util.module_from_spec(spec)
spec.loader.exec_module(helper)
JOB = "a" * 32


class Service:
    def __init__(self, polls=()):
        self.polls = list(polls)
        self.calls = []

    def __call__(self, base, path, value=None, timeout=120):
        self.calls.append((path, deepcopy(value)))
        if path.endswith("/cancel"):
            return json.dumps({"job_id": JOB, "status": "cancelling"}).encode()
        if path.endswith("/jobs"):
            return json.dumps({"job_id": JOB, "status": "queued"}).encode()
        result = self.polls.pop(0)
        if isinstance(result, BaseException):
            raise result
        return json.dumps({"job_id": JOB, **result}).encode()


def test_helper_waits_for_actual_terminal_result_not_accepted_job(tmp_path):
    result = {"status": "complete", "report_id": JOB, "isotherm": [{"pressure_Pa": 100000}]}
    service = Service([{"status": "running", "stage": "Widom"}, {"status": "complete", "result": result}])
    job = helper.AdsorptionJob(service, "http://127.0.0.1:1234", tmp_path)
    job.start({"structure": {}})
    assert not job.safe_to_stop
    assert job.wait(interval=0) == result and job.safe_to_stop
    assert job.cancel_if_active() == result
    assert not any(path.endswith("/cancel") for path, _ in service.calls)
    assert json.loads((tmp_path / "job-status.json").read_text())["status"] == "complete"


def test_interrupt_cancels_own_job_then_waits_for_worker_terminal_status(tmp_path):
    service = Service([KeyboardInterrupt(), {"status": "cancelled", "result": {"status": "cancelled", "downloads": {}}}])
    job = helper.AdsorptionJob(service, "http://127.0.0.1:1234", tmp_path)
    job.start({"structure": {}})
    with pytest.raises(KeyboardInterrupt):
        job.wait(interval=0)
    assert not job.safe_to_stop
    assert job.cancel_if_active(interval=0)["status"] == "cancelled"
    assert job.safe_to_stop
    assert [(p, v) for p, v in service.calls if p.endswith("/cancel")] == [(f"/api/drawing/science/adsorption/jobs/{JOB}/cancel", {})]


def test_uncertain_cancellation_never_allows_server_shutdown(tmp_path):
    service = Service()
    job = helper.AdsorptionJob(service, "http://127.0.0.1:1234", tmp_path)
    job.start({"structure": {}})
    with pytest.raises(RuntimeError, match="not confirmed"):
        job.cancel_if_active(timeout=0)
    assert job.status["status"] == "cancelling" and not job.safe_to_stop


@pytest.mark.parametrize("error,safe", [(helper.LocalAPIError(422, "invalid"), True), (OSError("connection lost"), False)])
def test_rejected_request_is_safe_but_unknown_creation_outcome_retains_server(tmp_path, error, safe):
    def unavailable(*args, **kwargs):
        raise error
    job = helper.AdsorptionJob(unavailable, "http://127.0.0.1:1234", tmp_path)
    with pytest.raises(type(error)):
        job.start({})
    assert job.safe_to_stop is safe
    if safe:
        assert job.cancel_if_active() is None
    else:
        with pytest.raises(RuntimeError, match="outcome is unknown"):
            job.cancel_if_active()


def test_mismatched_response_cannot_adopt_another_job(tmp_path):
    service = Service()
    job = helper.AdsorptionJob(service, "http://127.0.0.1:1234", tmp_path)
    job.start({})
    with pytest.raises(RuntimeError, match="job identity"):
        job._retain({"job_id": "b" * 32, "status": "complete"})
    assert job.identifier == JOB and not job.safe_to_stop


def test_failed_job_keeps_actual_partial_exports_without_claiming_success(tmp_path):
    result = {"status": "failed", "downloads": {ext: "/api/drawing/exports/" + JOB + "/" + ext for ext in ("json", "csv", "zip")}}
    content = {"json": json.dumps(result).encode(), "csv": b"pressure_Pa,gas\n", "zip": b"PK\x03\x04test-bundle"}
    def download(base, path):
        return content[path.rsplit("/", 1)[1]]
    artifacts = helper.save_result(tmp_path, result, download, "http://127.0.0.1:1234", "structure")
    assert json.loads((tmp_path / "result.json").read_text())["status"] == "failed"
    for row in artifacts:
        payload = (tmp_path / row["file"]).read_bytes()
        assert len(payload) == row["bytes"] and hashlib.sha256(payload).hexdigest() == row["sha256"]
    assert {r["file"] for r in artifacts} == {"structure.json", "structure.csv", "structure.zip"}
