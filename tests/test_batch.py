import logging
import pickle
from pathlib import Path
from tempfile import TemporaryDirectory
from typing import Any, Mapping, List, Tuple
from unittest import mock

import numpy
from alphafold.model.features import FeatureDict

from colabfold.batch import run
from colabfold.download import download_alphafold_params


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
    host_url="https://a3m.mmseqs.com",
):
    assert (
        prefix
        and use_filter
        and not use_templates
        and not use_env
        and not filter
        and host_url == "https://a3m.mmseqs.com"
    )
    if x == "YYDPETGTWY":
        return ">101\nYYDPETGTWY\n"
    elif x == "IKKILSKIKKLLK":
        return ">101\nIKKILSKIKKLLK"
    else:
        assert False, x


def test_batch(pytestconfig, caplog):
    caplog.set_level(logging.INFO)
    data_dir = pytestconfig.rootpath

    # We'll also want to mock that out later
    download_alphafold_params(data_dir.joinpath("params"))

    known_inputs = load_known_input(pytestconfig.rootpath)

    with TemporaryDirectory() as input_dir, TemporaryDirectory() as result_dir:
        Path(input_dir).joinpath("7_to_be_renamed.fasta").write_text(
            ">5AWL_1\nYYDPETGTWY"
        )
        Path(input_dir).joinpath("6BGN.fasta").write_text(">6A5J\nIKKILSKIKKLLK")

        with mock.patch(
            "alphafold.model.model.RunModel.predict",
            lambda self, feat: predict(known_inputs, feat),
        ), mock.patch(
            "colabfold_.batch.run_mmseqs2",
            mock_run_mmseqs2,
        ):
            run(
                input_dir,
                result_dir,
                use_templates=False,
                use_amber=False,
                use_env=False,
                num_models=1,
                homooligomer=1,
                data_dir=data_dir,
                do_not_overwrite_results=False,
            )

    assert caplog.messages == [
        "Found 4 citations for tools or databases",
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
