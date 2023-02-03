import json
import lzma
import os
import pickle
from pathlib import Path
from typing import List, Tuple, Mapping, Any, Dict

import numpy

from alphafold.model.features import FeatureDict
from alphafold.model.model import RunModel
from colabfold.colabfold import run_mmseqs2

# Copy the original method before mocking
original_run_model = RunModel.predict


class MockRunModel:
    """Mocks FeatureDict -> prediction

    The class is stateful, i.e. predictions need to be done in the given order

    msa_feat is a) large and b) has some variance between machines, so we ignore it
    """

    fixture_dir: Path
    predictions: List[str]
    pos: int

    def __init__(self, fixture_dir: Path, predictions: List[str]):
        self.fixture_dir = fixture_dir
        self.predictions = predictions
        self.pos = 0

    def predict(
        self, model_runner: RunModel, feat: FeatureDict, random_seed: int, prediction_callback: Any = None
    ) -> Mapping[str, Any]:
        """feat["msa"] or feat["msa_feat"] for normal/complexes is non-deterministic, so we remove it before storing,
        but we keep it for predicting or returning, where we need it for plotting"""
        feat_no_msa = dict(feat)
        if "msa_feat" in feat_no_msa.keys():
            del feat_no_msa["msa_feat"]
        elif "msa" in feat_no_msa.keys():
            del feat_no_msa["msa"]
        else:
            raise AssertionError("neither msa nor msa_feat in feat")

        prediction_file = self.fixture_dir.joinpath(
            self.predictions[self.pos]
        ).joinpath("model_prediction_result.pkl.xz")
        input_fix_file = self.fixture_dir.joinpath(self.predictions[self.pos]).joinpath(
            "model_input_fix.pkl.xz"
        )
        self.pos += 1

        if (
            not prediction_file.is_file() or not input_fix_file.is_file()
        ) and os.environ.get("UPDATE_SNAPSHOTS"):
            print("Running new prediction")
            with lzma.open(input_fix_file) as fp:
                pickle.dump(feat_no_msa, fp)
            prediction, _ = original_run_model(model_runner, feat)
            del prediction["distogram"]
            del prediction["experimentally_resolved"]
            del prediction["masked_msa"]
            del prediction["aligned_confidence_probs"]
            with lzma.open(prediction_file) as fp:
                pickle.dump(prediction, fp)

        with lzma.open(input_fix_file) as input_fix_fp:
            input_fix = pickle.load(input_fix_fp)
        with lzma.open(prediction_file) as prediction_fp:
            prediction = pickle.load(prediction_fp)

        is_same = True
        for key in input_fix:
            if (
                key not in feat_no_msa
                or feat_no_msa[key].shape != input_fix[key].shape
                or not numpy.allclose(feat_no_msa[key], input_fix[key])
            ):
                is_same = False
                break

        if is_same:
            return prediction, 3

        if os.environ.get("UPDATE_SNAPSHOTS"):
            print("Running new prediction")
            with lzma.open(input_fix_file, "wb") as fp:
                pickle.dump(feat_no_msa, fp)
            prediction, _ = original_run_model(model_runner, feat)
            with lzma.open(prediction_file, "wb") as fp:
                pickle.dump(prediction, fp)
            return prediction, _
        else:
            for key in input_fix:
                # Generate a more helpful error message
                assert feat_no_msa[key].shape != input_fix[
                    key
                ].shape and numpy.allclose(feat_no_msa[key], input_fix[key]), key


class MMseqs2Mock:
    """Mocks out the call to the mmseqs2 api

    Each test has its own json file which contains the run_mmseqs2 input data in the
    config field and the saved response. To update responses or to add new tests,
    set the UPDATE_SNAPSHOTS env var (e.g. `UPDATE_SNAPSHOTS=1 pytest`
    """

    data_file: Path
    saved_responses: List[Dict[str, Any]]

    def __init__(self, rootpath: Path, name: str):
        self.data_file = (
            rootpath.joinpath("test-data/mmseqs-api-reponses")
            .joinpath(name)
            .with_suffix(".json")
        )
        if os.environ.get("UPDATE_SNAPSHOTS") and not self.data_file.is_file():
            self.data_file.write_text("[]")
        with self.data_file.open() as fp:
            self.saved_responses = []
            for saved_response in json.load(fp):
                # Join lines we've split before
                response = join_lines(saved_response["response"])
                self.saved_responses.append(
                    {"config": saved_response["config"], "response": response}
                )

    def mock_run_mmseqs2(
        self,
        query,
        prefix,
        use_env=True,
        use_filter=True,
        use_templates=False,
        filter=None,
        use_pairing=False,
        host_url="https://a3m.mmseqs.com",
    ):
        assert prefix
        config = {
            "query": query,
            "use_env": use_env,
            "use_filter": use_filter,
            "use_templates": use_templates,
            "filter": filter,
            "use_pairing": use_pairing,
        }

        for saved_response in self.saved_responses:
            if saved_response["config"] == config:
                return saved_response["response"]

        if os.environ.get("UPDATE_SNAPSHOTS"):
            print(f"\nrun_mmseqs2 with {config}")
            response = run_mmseqs2(
                x=config["query"],
                prefix=prefix,
                use_env=config["use_env"],
                use_filter=config["use_filter"],
                use_templates=config["use_templates"],
                filter=config["filter"],
                use_pairing=config["use_pairing"],
                host_url=host_url,
            )
            # Split lines so we get a readable json file
            response = split_lines(response)
            self.saved_responses.append({"config": config, "response": response})
            self.data_file.write_text(json.dumps(self.saved_responses, indent=2))
        else:
            assert False, config


def split_lines(x):
    """Split each files into a list of lines"""
    if isinstance(x, list):
        return [split_lines(i) for i in x]
    elif isinstance(x, str):
        return x.splitlines()
    else:
        raise TypeError(f"{type(x)} {str(x)[:20]}")


def join_lines(x):
    """Inverse of split_lines"""
    if all(isinstance(i, str) for i in x):
        return "\n".join(x)
    elif all(isinstance(i, list) for i in x):
        return [join_lines(i) for i in x]
    else:
        raise TypeError(f"{[type(i) for i in x]} {str(x)[:20]}")
