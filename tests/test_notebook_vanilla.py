"""Keep the notebook's vanilla command independent of RGI state and dependencies."""

import builtins
import copy
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path
from types import SimpleNamespace

import pytest


@pytest.mark.parametrize(
    ("notebook", "model", "extra_args"),
    [
        ("ColabFold2_preview.ipynb", "openbind0", []),
        ("ColabFold2_preview.ipynb", "af2_ptm", ["--model_dir=af2_weights"]),
        ("ColabFold2_preview.ipynb", "chai1", ["--use_esm_embeddings"]),
        ("ColabFold2_preview.ipynb", "esmfold2", ["--use_esm_embeddings"]),
        ("ColabFold2_preview.ipynb", "esmfold2_lm300m", ["--use_esm_embeddings"]),
        ("AlphaFold3_of3.ipynb", "alphafold3", ["--model_dir=af3_weights"]),
    ],
)
def test_vanilla_prediction_ignores_previous_rgi_state(
    notebook, model, extra_args, tmp_path, monkeypatch
):
    root = Path(__file__).resolve().parents[1]
    cells = json.loads((root / notebook).read_text())["cells"]
    source = next(
        "".join(cell["source"]) for cell in cells if cell["metadata"]["id"] == "run"
    )
    monkeypatch.chdir(tmp_path)
    monkeypatch.setenv("AF3_NB_OVERRIDES", "{}")
    Path("inputs").mkdir()
    original = {
        "name": "previous_run",
        "modelSeeds": [42],
        "sequences": [
            {
                "protein": {
                    "id": "A",
                    "sequence": "ACDEFGHIK",
                    "unpairedMsa": ">query\nACDEFGHIK\n",
                    "pairedMsa": "",
                    "templates": [],
                }
            }
        ],
    }
    guided = copy.deepcopy(original)
    guided["restraints_config"] = {"invalid_previous_config": True}
    guided["sequences"][0]["protein"]["conformer_restraints"] = True
    namespace = {
        "model": model,
        "sys": sys,
        "json": json,
        "fold_input": guided,
        "basejob": "test",
        "INPUT_DIR": "inputs",
        "OUTPUT_DIR": "outputs",
        "CACHE_DIR": "cache",
        "AF2_DIR": "af2_weights",
        "NATIVE_DIR": "af3_weights",
        "IS_AF2": model.startswith("af2_"),
        "IS_AF3": model == "alphafold3",
        "msa_mode": "single_sequence",
        "on_existing": "overwrite",
        "use_rgi": False,
        "restraints_config": "invalid",
        "ref_pdb": "missing.pdb",
    }
    commands = []

    def run(command, **kwargs):
        assert command[0] == "nvidia-smi", command
        return SimpleNamespace(stdout="7.5\n", returncode=0)

    def popen(command, **kwargs):
        commands.append(command)
        output = Path(namespace["job_dir"])
        output.mkdir(parents=True)
        (output / "result.cif").touch()
        return SimpleNamespace(stdout=iter([]), wait=lambda: 0)

    real_import = builtins.__import__

    def without_rgi_toolkit(name, *args, **kwargs):
        if name == "rgi_toolkit" or name.startswith("rgi_toolkit."):
            raise AssertionError("Vanilla prediction must not import rgi-toolkit")
        return real_import(name, *args, **kwargs)

    monkeypatch.setattr(subprocess, "run", run)
    monkeypatch.setattr(subprocess, "Popen", popen)
    monkeypatch.setattr(builtins, "__import__", without_rgi_toolkit)
    exec(source, namespace)

    assert len(commands) == 1
    command = commands[0]
    assert command[:2] == [sys.executable, "run_alphafold.py"]
    assert f"--model={model}" in command
    assert "--num_recycles=10" in command
    assert "--num_diffusion_samples=5" in command
    assert [
        arg
        for arg in command
        if arg == "--use_esm_embeddings" or arg.startswith("--model_dir=")
    ] == extra_args
    assert not any("rgi" in argument for argument in command[1:])
    written = json.loads(Path(namespace["json_path"]).read_text())
    assert written["sequences"] == original["sequences"]
    assert written["modelSeeds"] == original["modelSeeds"]
    assert "restraints_config" not in written
    assert namespace["rgi_config"] is None
    assert not (Path(namespace["job_dir"]) / "cite.bibtex").exists()


@pytest.mark.parametrize(
    ("rgi_cached", "produce_structure"), [(False, True), (True, True), (True, False)]
)
def test_native_boltz_vanilla_run_all_switches_back_to_upstream(
    rgi_cached, produce_structure, tmp_path, monkeypatch
):
    yaml = pytest.importorskip(
        "yaml", reason="The native Boltz notebook requires PyYAML."
    )

    root = Path(__file__).resolve().parents[1]
    cells = json.loads((root / "Boltz1.ipynb").read_text())["cells"]
    sources = {cell["metadata"]["id"]: "".join(cell["source"]) for cell in cells}
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(shutil, "which", lambda name: name)
    Path("test").mkdir()
    Path("test/query.a3m").write_text(">A\nACDEFGHIK\n")
    stale_output = Path("test/boltz_results_test")
    (stale_output / "processed").mkdir(parents=True)
    (stale_output / "processed/manifest.json").write_text("{}")
    stale_predictions = stale_output / "predictions/test"
    stale_predictions.mkdir(parents=True)
    (stale_predictions / "stale.cif").touch()
    if rgi_cached:
        (stale_output / ".rgi-enabled").touch()
        (stale_output / "cite.bibtex").write_text("@article{hori2026rgi}\n")
    namespace = {
        "os": os,
        "use_rgi": False,
        "boltz_python": ".cache/boltz-rgi/bin/python",
        "restraints_config": "invalid",
        "ref_pdb": "missing.pdb",
        "fasta_entries": [(">A|protein", "ACDEFGHIK"), (">D|smiles", "CCO")],
        "msa_mode": "single_sequence",
        "jobname": "test",
    }
    commands = []

    def run(command, **kwargs):
        commands.append(command)
        if "predict" in command:
            if rgi_cached:
                assert not stale_output.exists()
            else:
                assert (stale_output / "processed/manifest.json").is_file()
                assert (stale_predictions / "stale.cif").is_file()
            assert Path("test/query.a3m").read_text() == ">A\nACDEFGHIK\n"
            if produce_structure:
                output = stale_output / "predictions/test"
                output.mkdir(parents=True, exist_ok=True)
                (output / "result.cif").touch()
        return SimpleNamespace(returncode=0)

    real_import = builtins.__import__

    def without_rgi_toolkit(name, *args, **kwargs):
        if name == "rgi_toolkit" or name.startswith("rgi_toolkit."):
            raise AssertionError("Vanilla prediction must not import rgi-toolkit")
        return real_import(name, *args, **kwargs)

    monkeypatch.setattr(subprocess, "run", run)
    monkeypatch.setattr(builtins, "__import__", without_rgi_toolkit)
    exec(sources["4eXNO1JJHYrB"], namespace)
    if produce_structure:
        exec(sources["bgaBXxXtIAu9"], namespace)
    else:
        with pytest.raises(
            RuntimeError, match="Prediction produced no structure files"
        ):
            exec(sources["bgaBXxXtIAu9"], namespace)
    assert not (stale_output / ".rgi-enabled").exists()
    assert not (stale_output / "cite.bibtex").exists()

    assert any(command[-1] == "boltz==2.2.1" for command in commands)
    assert not any(
        "boltz_restr" in arg or "rgi_toolkit" in arg
        for command in commands
        for arg in command
    )
    command = commands[-1]
    assert command[:4] == [
        ".cache/boltz-vanilla/bin/python",
        "-m",
        "boltz.main",
        "predict",
    ]
    assert command[command.index("--model") + 1] == "boltz1"
    assert command[command.index("--recycling_steps") + 1] == "3"
    assert command[command.index("--diffusion_samples") + 1] == "1"
    written = yaml.safe_load(Path("test/test.yaml").read_text())
    assert "restraints_config" not in written
    assert all(
        "conformer_restraints" not in next(iter(entity.values()))
        for entity in written["sequences"]
    )


def test_native_boltz_rgi_reprocesses_when_target_changes(tmp_path, monkeypatch):
    pytest.importorskip("yaml", reason="The native Boltz notebook requires PyYAML.")
    root = Path(__file__).resolve().parents[1]
    cells = json.loads((root / "Boltz1.ipynb").read_text())["cells"]
    source = next(
        "".join(cell["source"])
        for cell in cells
        if cell["metadata"]["id"] == "bgaBXxXtIAu9"
    )
    monkeypatch.chdir(tmp_path)
    output = Path("test/boltz_results_test")
    (output / "processed").mkdir(parents=True)
    (output / "processed/manifest.json").write_text("{}")
    namespace = {
        "Path": Path,
        "subprocess": subprocess,
        "use_rgi": True,
        "boltz_python": "boltz-rgi-python",
        "fasta_entries": [(">A|protein", "ACDEFGHIK")],
        "msa_mode": "single_sequence",
        "jobname": "test",
    }

    def config_from_fields(fields):
        return {"target_distance": fields["target_distance"]}, ""

    monkeypatch.setitem(
        sys.modules,
        "rgi_toolkit.notebook_colab",
        SimpleNamespace(config_from_fields=config_from_fields),
    )
    monkeypatch.setattr(
        "colabfold.rgi.boltz.notebook_input",
        lambda entries, *, config, **kwargs: {"restraints_config": config},
    )
    targets = []

    def run(command, **kwargs):
        assert (output / ".rgi-enabled").is_file()
        assert not (output / "processed").exists()
        targets.append(namespace["rgi_config"]["target_distance"])
        (output / "processed").mkdir()
        (output / "processed/manifest.json").write_text("{}")
        prediction = output / "predictions/test"
        prediction.mkdir(parents=True)
        (prediction / "result.cif").touch()
        return SimpleNamespace(returncode=0)

    monkeypatch.setattr(subprocess, "run", run)
    for target in (25, 35):
        namespace["target_distance"] = target
        exec(source, namespace)
        bibliography = (output / "cite.bibtex").read_text()
        assert bibliography.count("@article{hori2026rgi,") == 1
    assert targets == [25, 35]
