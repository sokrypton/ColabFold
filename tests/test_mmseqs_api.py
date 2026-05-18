"""Regression tests for the run_mmseqs2 retry loops.

Each of the four retry loops in run_mmseqs2 (submit, status, download, templates)
caps at a fixed number of attempts. status() and the templates loop used to
reset their counter inside the loop, which silently turned them into
unbounded retries — these tests pin the bounded behavior.
"""

import pytest
import requests

from colabfold.colabfold import run_mmseqs2

# Safety bound so a regression to unbounded retries fails fast instead of
# spinning forever (time.sleep is monkeypatched to a no-op below).
MAX_RETRY_CALLS = 50


class JsonResponse:
    def __init__(self, payload):
        self.payload = payload
        self.text = str(payload)
        self.content = b""

    def json(self):
        return self.payload


def _failing_get(calls, message):
    def get(*args, **kwargs):
        calls["get"] += 1
        if calls["get"] > MAX_RETRY_CALLS:
            raise AssertionError(
                f"retry loop exceeded {MAX_RETRY_CALLS} calls — "
                "the bound on retries has regressed"
            )
        raise requests.exceptions.ConnectionError(message)

    return get


def test_status_retries_are_bounded(monkeypatch, tmp_path):
    calls = {"get": 0}

    monkeypatch.setattr(
        "colabfold.colabfold.requests.post",
        lambda *a, **kw: JsonResponse({"status": "RUNNING", "id": "job-id"}),
    )
    monkeypatch.setattr(
        "colabfold.colabfold.requests.get",
        _failing_get(calls, "status failed"),
    )
    monkeypatch.setattr("colabfold.colabfold.time.sleep", lambda _s: None)

    with pytest.raises(requests.exceptions.ConnectionError, match="status failed"):
        run_mmseqs2(
            "ACDE",
            str(tmp_path / "msa"),
            use_env=False,
            user_agent="colabfold/test",
        )

    assert calls["get"] == 6


def test_download_retries_are_bounded(monkeypatch, tmp_path):
    calls = {"get": 0}

    monkeypatch.setattr(
        "colabfold.colabfold.requests.post",
        lambda *a, **kw: JsonResponse({"status": "COMPLETE", "id": "job-id"}),
    )
    monkeypatch.setattr(
        "colabfold.colabfold.requests.get",
        _failing_get(calls, "download failed"),
    )
    monkeypatch.setattr("colabfold.colabfold.time.sleep", lambda _s: None)

    with pytest.raises(requests.exceptions.ConnectionError, match="download failed"):
        run_mmseqs2(
            "ACDE",
            str(tmp_path / "msa"),
            use_env=False,
            user_agent="colabfold/test",
        )

    assert calls["get"] == 6


def test_templates_retries_are_bounded(monkeypatch, tmp_path):
    # Pre-create the MSA-stage artifacts so run_mmseqs2 skips submit/status/
    # download (tar_gz_file exists) and a3m extraction (uniref.a3m exists),
    # going straight to the templates fetch we want to exercise.
    msa_dir = tmp_path / "msa_all"
    msa_dir.mkdir()
    (msa_dir / "out.tar.gz").write_bytes(b"")
    (msa_dir / "uniref.a3m").write_text(">seq1\nACDE\n")
    (msa_dir / "pdb70.m8").write_text(
        "101\t1abc_A\t50.0\t100\t10\t1\t1\t100\t1\t100\t1e-10\t100.0\n"
    )

    calls = {"get": 0}
    monkeypatch.setattr(
        "colabfold.colabfold.requests.get",
        _failing_get(calls, "template failed"),
    )
    monkeypatch.setattr("colabfold.colabfold.time.sleep", lambda _s: None)

    with pytest.raises(requests.exceptions.ConnectionError, match="template failed"):
        run_mmseqs2(
            "ACDE",
            str(tmp_path / "msa"),
            use_env=False,
            use_templates=True,
            user_agent="colabfold/test",
        )

    assert calls["get"] == 6
