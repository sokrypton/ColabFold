import logging
import pickle
from pathlib import Path
from typing import Any, Mapping, List, Tuple
from unittest import mock

import numpy
from absl import logging as absl_logging
from alphafold.model.features import FeatureDict

from colabfold.batch import run, get_queries
from colabfold.download import download_alphafold_params


def test_get_queries(caplog, tmp_path):
    tmp_path.joinpath("5AWL_1.fasta").write_text(">5AWL_1\nYYDPETGTWY")
    tmp_path.joinpath("6A5J.fasta").write_text(">6A5J\nIKKILSKIKKLLK")

    queries = get_queries(tmp_path)
    assert queries == [("5AWL_1", "YYDPETGTWY", None), ("6A5J", "IKKILSKIKKLLK", None)]
    assert caplog.messages == []


def predict(
    known_inputs: List[Tuple[FeatureDict, Mapping[str, Any]]], feat: FeatureDict
) -> Mapping[str, Any]:
    for input_fix, prediction in known_inputs:
        try:
            numpy.testing.assert_equal(feat, input_fix)
            return prediction
        except AssertionError:
            continue
    else:
        raise AssertionError("input not stored")


def load_known_input(root_path: Path) -> List[Tuple[FeatureDict, Mapping[str, Any]]]:
    known_inputs = []
    for input_fix_file in root_path.joinpath("test-data/run_batch").glob(
        "*/*_input_fix.pkl"
    ):
        prediction_file = input_fix_file.with_name(
            input_fix_file.name.replace("input_fix", "prediction_result")
        )
        with input_fix_file.open("rb") as input_fix_fp, prediction_file.open(
            "rb"
        ) as prediction_fp:
            known_inputs.append((pickle.load(input_fix_fp), pickle.load(prediction_fp)))
    return known_inputs


def mock_run_mmseqs2(
    x,
    prefix,
    use_env=True,
    use_filter=True,
    use_templates=False,
    filter=None,
    use_pairing=False,
    host_url="https://a3m.mmseqs.com",
):
    assert (
        prefix
        and use_filter
        and not use_templates
        and use_env
        and not filter
        and host_url == "https://a3m.mmseqs.com"
    )
    if x == "YYDPETGTWY":
        return ">101\nYYDPETGTWY\n"
    elif x == "IKKILSKIKKLLK":
        return ">101\nIKKILSKIKKLLK"
    else:
        assert False, x


def test_batch(pytestconfig, caplog, tmp_path):
    caplog.set_level(logging.INFO)
    # otherwise jax will tell us about its search for devices
    absl_logging.set_verbosity("error")

    data_dir = pytestconfig.rootpath

    # We'll also want to mock that out later
    download_alphafold_params(data_dir.joinpath("params"))

    known_inputs = load_known_input(pytestconfig.rootpath)
    queries = [("5AWL_1", "YYDPETGTWY", None), ("6A5J", "IKKILSKIKKLLK", None)]

    with mock.patch(
        "alphafold.model.model.RunModel.predict",
        lambda self, feat: predict(known_inputs, feat),
    ), mock.patch("colabfold.batch.run_mmseqs2", mock_run_mmseqs2):
        run(
            queries,
            tmp_path,
            use_templates=False,
            use_amber=False,
            msa_mode="MMseqs2 (UniRef+Environmental)",
            num_models=1,
            homooligomer=1,
            data_dir=data_dir,
            do_not_overwrite_results=False,
            rank=0,
            pair_mode="unpaired+paired",
        )

    # Very simple test, it would be better to check coordinates
    assert (
        len(tmp_path.joinpath("5AWL_1_unrelaxed_model_1.pdb").read_text().splitlines())
        == 92
    )
    assert (
        len(tmp_path.joinpath("6A5J_unrelaxed_model_1.pdb").read_text().splitlines())
        == 108
    )

    assert caplog.messages == [
        "Found 5 citations for tools or databases",
        "Predicting 2 structures",
        "Running: 5AWL_1",
        "running model_1",
        "reranking models based on avg. predicted lDDT",
        "model_1 94.2",
        "Running: 6A5J",
        "running model_1",
        "reranking models based on avg. predicted lDDT",
        "model_1 89.5",
    ]
