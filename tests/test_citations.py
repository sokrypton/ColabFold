import importlib.metadata
import json
import shutil
import sys
from pathlib import Path
from types import SimpleNamespace
from zipfile import ZipFile

import pytest

from colabfold.citations import write_bibtex, write_rgi_citation


@pytest.mark.parametrize("existing", [False, True])
def test_rgi_citation_preserves_existing_bibliography_on_rerun(tmp_path, existing):
    original = b""
    if existing:
        path = write_bibtex("alphafold2_ptm", True, True, True, True, tmp_path)
        original = path.read_bytes()
        assert b"hori2026rgi" not in original

    path = write_rgi_citation(tmp_path)
    first = path.read_bytes()
    assert first.startswith(original)
    assert first.count(b"@article{hori2026rgi,") == 1
    assert b"doi     = {10.64898/2026.10.05.756905}" in first
    assert write_rgi_citation(tmp_path) == path
    assert path.read_bytes() == first


def test_rgi_citation_recognizes_existing_bibtex_key(tmp_path):
    path = tmp_path / "cite.bibtex"
    original = "@ARTICLE { hori2026rgi ,\n  title = {Existing citation}\n}\n"
    path.write_text(original, encoding="utf-8")
    write_rgi_citation(tmp_path)
    assert path.read_text(encoding="utf-8") == original


def test_runtime_output_citation_is_scoped_to_guided_jobs(tmp_path, monkeypatch):
    pytest.importorskip("jax")
    pytest.importorskip("rgi_toolkit")
    from colabfold.rgi import ALPHAFOLD3_COLABFOLD_VERSION, runtime

    class Input:
        @classmethod
        def from_json(cls, json_str, json_path=None):
            result = cls()
            result.name = json.loads(json_str)["name"]
            return result

        def sanitised_name(self):
            return self.name

    def write_outputs(results, output_dir, job_name):
        output_dir.mkdir()
        write_bibtex("", False, False, False, False, output_dir)

    runner = SimpleNamespace(
        folding_input=SimpleNamespace(Input=Input),
        predict_structure=lambda *args, **kwargs: [],
        write_fold_input_json=lambda *args: None,
        write_outputs=write_outputs,
    )
    model_runner = SimpleNamespace(
        model_name="boltz2", _rgi_apply=None, run_inference=lambda *args: None
    )
    original_version = importlib.metadata.version
    monkeypatch.setattr(
        importlib.metadata,
        "version",
        lambda name: (
            ALPHAFOLD3_COLABFOLD_VERSION
            if name == "alphafold3-colabfold"
            else original_version(name)
        ),
    )
    runtime.install(runner)
    for name, guided in (("guided", True), ("vanilla", False)):
        raw = {"name": name}
        if guided:
            raw["restraints_config"] = {}
        fold_input = Input.from_json(json.dumps(raw))
        results = runner.predict_structure(fold_input, model_runner)
        output_dir = tmp_path / name
        runner.write_outputs(results, output_dir=output_dir, job_name=name)
        bibliography = (output_dir / "cite.bibtex").read_text()
        assert "@article{Mirdita2022," in bibliography
        assert ("@article{hori2026rgi," in bibliography) is guided
        assert (output_dir / "rgi_report.json").exists() is guided


@pytest.mark.parametrize(
    ("notebook", "download_id", "output_dir"),
    [
        ("ColabFold2_preview.ipynb", "download", "outputs/test"),
        ("AlphaFold3_of3.ipynb", "download", "outputs/test"),
        ("Boltz1.ipynb", "jdSBSTOpaULF", "test/boltz_results_test"),
    ],
)
def test_downloaded_results_include_rgi_citation(
    notebook, download_id, output_dir, tmp_path, monkeypatch
):
    if notebook != "Boltz1.ipynb" and shutil.which("zip") is None:
        pytest.skip("The AF3 notebook download cell requires zip.")
    root = Path(__file__).resolve().parents[1]
    cells = json.loads((root / notebook).read_text())["cells"]
    source = next(
        "".join(cell["source"])
        for cell in cells
        if cell["metadata"]["id"] == download_id
    )
    monkeypatch.chdir(tmp_path)
    output = Path(output_dir)
    output.mkdir(parents=True)
    citation = write_rgi_citation(output).read_bytes()
    input_yaml = tmp_path / "test.yaml"
    input_yaml.write_text("version: 1\n")
    downloaded = []
    monkeypatch.setitem(
        sys.modules,
        "google.colab",
        SimpleNamespace(files=SimpleNamespace(download=downloaded.append)),
    )
    exec(
        source,
        {"jobname": "test", "OUTPUT_DIR": "outputs", "input_yaml": input_yaml},
    )
    assert len(downloaded) == 1
    with ZipFile(downloaded[0]) as archive:
        names = [name for name in archive.namelist() if name.endswith("cite.bibtex")]
        assert len(names) == 1
        assert archive.read(names[0]) == citation
