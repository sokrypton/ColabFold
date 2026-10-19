"""CLI forwarding tests that do not require the AlphaFold runtime or weights."""

import argparse
import runpy
import sys
from pathlib import Path
from unittest import mock

import pytest

from colabfold.caat_options import add_prediction_options, prediction_kwargs


@pytest.mark.parametrize("flag", ["--num-recycles", "--num-recycle", "--num_recycles"])
def test_recycle_aliases_and_defaults(flag):
    parser = argparse.ArgumentParser()
    names = add_prediction_options(parser)
    assert prediction_kwargs(parser.parse_args([]), names) == {}
    assert prediction_kwargs(parser.parse_args([flag, "0"]), names) == {"num_recycles": 0}


@pytest.mark.parametrize("script", ["run_e2e_pipeline.py", "run_attention_heads.py"])
def test_prediction_options_reach_every_run(script, tmp_path, monkeypatch):
    batch = mock.MagicMock()
    batch.get_queries.return_value = ([], False)
    for name in [
        "colabfold.download", "colabfold.utils", "alphafold.analysis.utils",
        "alphafold.analysis.attention_pipeline", "alphafold.model.modules",
    ]:
        monkeypatch.setitem(sys.modules, name, mock.MagicMock())
    monkeypatch.setitem(sys.modules, "colabfold.batch", batch)
    query = tmp_path / "query.a3m"
    query.touch()
    argv = [script, "--query-seq-path", str(query), "--query-name", "query",
            "--result-dir", str(tmp_path / "results"),
            "--num-recycles", "7", "--num-ensemble", "2",
            "--recycle-early-stop-tolerance", "0.0", "--model-order", "3,1",
            "--max-seq", "64", "--max-extra-seq", "128", "--use-dropout",
            "--msa-mode", "single_sequence", "--pair-strategy", "complete",
            "--templates", "--custom-template-path", "templates",
            "--disable-cluster-profile", "--rank", "ptm", "--save-recycles"]
    if script == "run_e2e_pipeline.py":
        argv += ["--target-seq-path", str(query), "--target-name", "target"]
    monkeypatch.setattr(sys, "argv", argv)
    runpy.run_path(str(Path(__file__).resolve().parents[1] / "scripts" / script), run_name="__main__")
    assert batch.run.call_count == (2 if script == "run_e2e_pipeline.py" else 1)
    for call in batch.run.call_args_list:
        for key, value in {
            "num_recycles": 7, "num_ensemble": 2, "recycle_early_stop_tolerance": 0.0,
            "model_order": [3, 1], "max_seq": 64, "max_extra_seq": 128,
            "use_dropout": True, "msa_mode": "single_sequence",
            "pairing_strategy": "complete", "use_templates": True,
            "custom_template_path": "templates", "use_cluster_profile": False,
            "rank_by": "ptm", "save_recycles": True,
        }.items():
            assert call.kwargs[key] == value
