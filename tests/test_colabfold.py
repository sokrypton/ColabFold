import json
import logging
import lzma
import os
import pickle
import re
from pathlib import Path
from typing import Any, Mapping, List, Tuple, Dict
from unittest import mock

import numpy
from absl import logging as absl_logging
from alphafold.model.features import FeatureDict
from alphafold.model.model import RunModel
from alphafold.model.tf import utils

from colabfold.batch import run
from colabfold.colabfold import run_mmseqs2
from colabfold.download import download_alphafold_params

# Copy the original method before mocking
original_run_model = RunModel.predict


class MockRunModel:
    """Mocks FeatureDict -> prediction

    msa_feat is a) large and b) has some variance between machines, so we ignore it
    """

    known_inputs: List[Tuple[FeatureDict, Mapping[str, Any]]]
    fixtures: Path

    def __init__(self, fixtures: Path):
        self.fixtures = fixtures
        self.known_inputs = []
        for prediction_file in fixtures.glob("*/*_prediction_result.pkl.xz"):
            input_fix_file = prediction_file.with_name(
                prediction_file.name.replace("prediction_result", "input_fix")
            )
            with lzma.open(input_fix_file) as input_fix_fp, lzma.open(
                prediction_file, "rb"
            ) as prediction_fp:
                input_fix = pickle.load(input_fix_fp)
                self.known_inputs.append((input_fix, pickle.load(prediction_fp)))

    def predict(
        self, model_runner: RunModel, feat: FeatureDict
    ) -> Tuple[Mapping[str, Any], Tuple[Any, Any]]:
        """feat["msa"] or feat["msa_feat"] for normal/complexes is non-deterministic, so we remove it before storing,
        but add it back before prediction or returning, as we need it for plotting"""
        is_complex = False
        if "msa_feat" in feat.keys():
            msa_feat = feat["msa_feat"]
            # noinspection PyUnresolvedReferences
            del feat["msa_feat"]
        elif "msa" in feat.keys():
            msa_feat = feat["msa"]
            # noinspection PyUnresolvedReferences
            del feat["msa"]
            is_complex = True
        else:
            raise AssertionError("neither msa nor msa_feat in feat")

        for input_fix, prediction in self.known_inputs:
            try:
                numpy.testing.assert_equal(feat, input_fix)
                # TODO: Also mock (recycles,tol) from the patches
                if is_complex:
                    # noinspection PyUnresolvedReferences
                    feat["msa"] = msa_feat
                else:
                    # noinspection PyUnresolvedReferences
                    feat["msa_feat"] = msa_feat
                return prediction, (None, None)
            except AssertionError:
                continue

        if os.environ.get("UPDATE_SNAPSHOTS"):
            print("Running new prediction")
            # This is real hacky and you'll have to rename the stuff yourself, sorry
            counter = 0
            while self.fixtures.joinpath(str(counter)).is_dir():
                counter += 1
            folder = self.fixtures.joinpath(str(counter))
            folder.mkdir(exist_ok=True, parents=True)
            with lzma.open(folder.joinpath(f"model_input_fix.pkl.xz"), "wb") as fp:
                pickle.dump(feat, fp)
            # Put msa_feat back, we need it for the prediction
            if is_complex:
                # noinspection PyUnresolvedReferences
                feat["msa"] = msa_feat
            else:
                # noinspection PyUnresolvedReferences
                feat["msa_feat"] = msa_feat
            prediction, (_, _) = original_run_model(model_runner, feat)
            self.known_inputs.append((feat, prediction))
            with lzma.open(
                folder.joinpath(f"model_prediction_result.pkl.xz"), "wb"
            ) as fp:
                pickle.dump(prediction, fp)
            return prediction, (None, None)
        else:
            for input_fix, prediction in self.known_inputs:
                try:
                    numpy.testing.assert_equal(feat, input_fix)
                    return prediction, (None, None)
                except AssertionError as e:
                    print(e)
            raise AssertionError("input not stored")


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


class MMseqs2Mock:
    """Mocks out the call to the mmseqs2 api

    Each test has its own json file which contains the run_mmseqs2 input data in the
    config field and the saved response. To update responses or to add new tests,
    set the UPDATE_SNAPSHOTS env var (e.g. `UPDATE_SNAPSHOTS=1 pytest`
    """

    data_file: Path
    saved_responses: List[Dict[str, Any]]

    def __init__(self, rootpath: Path, name: str):
        self.data_file = rootpath.joinpath(f"test-data/mmseqs-api-reponses/{name}.json")
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


def prepare_prediction_test(caplog):
    caplog.set_level(logging.INFO)
    # otherwise jax will tell us about its search for devices
    absl_logging.set_verbosity("error")
    # We'll also want to mock that out later
    download_alphafold_params(True)
    download_alphafold_params(False)
    # alphafold uses a method called `make_random_seed`, which deterministically starts with a seed
    # of zero and increases it by one for each protein. This means the input features would become
    # dependent on the number and order of tests. Here we just reset the seed to 0
    utils.seed_maker = utils.SeedMaker()


def test_batch(pytestconfig, caplog, tmp_path):
    prepare_prediction_test(caplog)

    queries = [("5AWL_1", "YYDPETGTWY", None), ("6A5J", "IKKILSKIKKLLK", None)]

    mock_run_model = MockRunModel(pytestconfig.rootpath.joinpath("test-data/batch"))
    mock_run_mmseqs = MMseqs2Mock(pytestconfig.rootpath, "batch").mock_run_mmseqs2
    with mock.patch(
        "alphafold.model.model.RunModel.predict",
        lambda model_runner, feat: mock_run_model.predict(model_runner, feat),
    ), mock.patch("colabfold.batch.run_mmseqs2", mock_run_mmseqs):
        run(
            queries,
            tmp_path,
            use_templates=False,
            use_amber=False,
            msa_mode="MMseqs2 (UniRef+Environmental)",
            model_type="auto",
            num_models=1,
            num_recycles=3,
            model_order=[1, 2, 3, 4, 5],
            is_complex=False,
            keep_existing_results=False,
            rank_mode="auto",
            pair_mode="unpaired+paired",
            stop_at_score=100,
        )

    assert caplog.messages == [
        "Found 5 citations for tools or databases",
        "Query 1/2: 5AWL_1 (length 10)",
        "Running model_1",
        "model_1 took 0.0s with pLDDT 94.2",
        "reranking models based on avg. predicted lDDT",
        "Query 2/2: 6A5J (length 13)",
        "Running model_1",
        "model_1 took 0.0s with pLDDT 90.8",
        "reranking models based on avg. predicted lDDT",
        "Done",
    ]

    # Very simple test, it would be better to check coordinates
    assert (
        len(
            tmp_path.joinpath("5AWL_1_unrelaxed_model_1_rank_1.pdb")
            .read_text()
            .splitlines()
        )
        == 96
    )
    assert (
        len(
            tmp_path.joinpath("6A5J_unrelaxed_model_1_rank_1.pdb")
            .read_text()
            .splitlines()
        )
        == 112
    )
    assert tmp_path.joinpath("config.json").is_file()


def test_single_sequence(pytestconfig, caplog, tmp_path):
    prepare_prediction_test(caplog)

    queries = [("5AWL_1", "YYDPETGTWY", None)]

    mock_run_model = MockRunModel(pytestconfig.rootpath.joinpath("test-data/batch"))
    mock_run_mmseqs = MMseqs2Mock(pytestconfig.rootpath, "batch").mock_run_mmseqs2
    with mock.patch(
        "alphafold.model.model.RunModel.predict",
        lambda model_runner, feat: mock_run_model.predict(model_runner, feat),
    ), mock.patch("colabfold.batch.run_mmseqs2", mock_run_mmseqs):
        run(
            queries,
            tmp_path,
            use_templates=False,
            use_amber=False,
            msa_mode="single_sequence",
            model_type="auto",
            num_models=1,
            num_recycles=3,
            model_order=[1, 2, 3, 4, 5],
            is_complex=False,
            keep_existing_results=False,
            rank_mode="auto",
            pair_mode="unpaired+paired",
            stop_at_score=100,
        )

    assert caplog.messages == [
        "Found 2 citations for tools or databases",
        "Query 1/1: 5AWL_1 (length 10)",
        "Running model_1",
        "model_1 took 0.0s with pLDDT 94.2",
        "reranking models based on avg. predicted lDDT",
        "Done",
    ]

    # Very simple test, it would be better to check coordinates
    assert (
        len(
            tmp_path.joinpath("5AWL_1_unrelaxed_model_1_rank_1.pdb")
            .read_text()
            .splitlines()
        )
        == 96
    )
    assert tmp_path.joinpath("config.json").is_file()


def test_complex(pytestconfig, caplog, tmp_path):
    prepare_prediction_test(caplog)

    pdb_3g50_A = "MRILPISTIKGKLNEFVDAVSSTQDQITITKNGAPAAVLVGADEWESLQETLYWLAQPGIRESIAEADADIASGRTYGEDEIRAEFGVPRRPH"
    pdb_3g50_B = "MPYTVRFTTTARRDLHKLPPRILAAVVEFAFGDLSREPLRVGKPLRRELAGTFSARRGTYRLLYRIDDEHTTVVILRVDHRADIYRR"
    queries = [("3G5O_A_3G5O_B", [pdb_3g50_A, pdb_3g50_B], None)]

    mock_run_model = MockRunModel(pytestconfig.rootpath.joinpath("test-data/complex"))
    mock_run_mmseqs2 = MMseqs2Mock(pytestconfig.rootpath, "complex").mock_run_mmseqs2
    with mock.patch(
        "alphafold.model.model.RunModel.predict",
        lambda model_runner, feat: mock_run_model.predict(model_runner, feat),
    ), mock.patch("colabfold.batch.run_mmseqs2", mock_run_mmseqs2):
        run(
            queries,
            tmp_path,
            use_templates=False,
            use_amber=False,
            msa_mode="MMseqs2 (UniRef+Environmental)",
            model_type="auto",
            num_models=1,
            num_recycles=3,
            model_order=[1, 2, 3, 4, 5],
            is_complex=True,
            keep_existing_results=False,
            rank_mode="auto",
            pair_mode="unpaired+paired",
            stop_at_score=100,
        )

    messages = list(caplog.messages)
    # noinspection PyUnresolvedReferences
    messages[3] = re.sub(r"\d+\.\d+s", "0.0s", messages[3])
    assert messages == [
        "Found 5 citations for tools or databases",
        "Query 1/1: 3G5O_A_3G5O_B (length 180)",
        "Running model_1",
        "model_1 took 0.0s with pLDDT 94.4",
        "reranking models based on avg. predicted lDDT",
        "Done",
    ]


def test_complex_ptm(pytestconfig, caplog, tmp_path):
    prepare_prediction_test(caplog)

    pdb_3g50_A = "MRILPISTIKGKLNEFVDAVSSTQDQITITKNGAPAAVLVGADEWESLQETLYWLAQPGIRESIAEADADIASGRTYGEDEIRAEFGVPRRPH"
    pdb_3g50_B = "MPYTVRFTTTARRDLHKLPPRILAAVVEFAFGDLSREPLRVGKPLRRELAGTFSARRGTYRLLYRIDDEHTTVVILRVDHRADIYRR"
    queries = [("3G5O_A_3G5O_B", [pdb_3g50_A, pdb_3g50_B], None)]

    mock_run_model = MockRunModel(
        pytestconfig.rootpath.joinpath("test-data/complex_ptm")
    )
    mock_run_mmseqs2 = MMseqs2Mock(pytestconfig.rootpath, "complex").mock_run_mmseqs2
    with mock.patch(
        "alphafold.model.model.RunModel.predict",
        lambda model_runner, feat: mock_run_model.predict(model_runner, feat),
    ), mock.patch("colabfold.batch.run_mmseqs2", mock_run_mmseqs2):
        run(
            queries,
            tmp_path,
            use_templates=False,
            use_amber=False,
            msa_mode="MMseqs2 (UniRef+Environmental)",
            model_type="AlphaFold2",
            num_models=1,
            num_recycles=3,
            model_order=[1, 2, 3, 4, 5],
            is_complex=True,
            keep_existing_results=False,
            rank_mode="auto",
            pair_mode="unpaired+paired",
            stop_at_score=100,
        )

    messages = list(caplog.messages)
    # noinspection PyUnresolvedReferences
    messages[3] = re.sub(r"\d+\.\d+s", "0.0s", messages[3])
    assert messages == [
        "Found 5 citations for tools or databases",
        "Query 1/1: 3G5O_A_3G5O_B (length 180)",
        "Running model_1",
        "model_1 took 0.0s with pLDDT 91.6",
        "reranking models based on avg. predicted lDDT",
        "Done",
    ]


def test_complex_monomer_ptm(pytestconfig, caplog, tmp_path):
    prepare_prediction_test(caplog)

    A = "PIAQIHILEGRSDEQKETLIREVSEAISRSLDAPLTSVRVIITEMAKGHFGIGGELASK"
    queries = [("A_A", [A, A], None)]

    mock_run_model = MockRunModel(
        pytestconfig.rootpath.joinpath("test-data/complex_monomer_ptm")
    )
    mock_run_mmseqs2 = MMseqs2Mock(
        pytestconfig.rootpath, "complex_monomer"
    ).mock_run_mmseqs2
    with mock.patch(
        "alphafold.model.model.RunModel.predict",
        lambda model_runner, feat: mock_run_model.predict(model_runner, feat),
    ), mock.patch("colabfold.batch.run_mmseqs2", mock_run_mmseqs2):
        run(
            queries,
            tmp_path,
            use_templates=False,
            use_amber=False,
            msa_mode="MMseqs2 (UniRef+Environmental)",
            model_type="AlphaFold2",
            num_models=1,
            num_recycles=3,
            model_order=[1, 2, 3, 4, 5],
            is_complex=True,
            keep_existing_results=False,
            rank_mode="auto",
            pair_mode="unpaired+paired",
            stop_at_score=100,
        )

    messages = list(caplog.messages)
    # noinspection PyUnresolvedReferences
    messages[3] = re.sub(r"\d+\.\d+s", "0.0s", messages[3])
    assert messages == [
        "Found 5 citations for tools or databases",
        "Query 1/1: A_A (length 118)",
        "Running model_1",
        "model_1 took 0.0s with pLDDT 95.5",
        "reranking models based on avg. predicted lDDT",
        "Done",
    ]


def test_complex_monomer(pytestconfig, caplog, tmp_path):
    prepare_prediction_test(caplog)

    A = "PIAQIHILEGRSDEQKETLIREVSEAISRSLDAPLTSVRVIITEMAKGHFGIGGELASK"
    queries = [("A_A", [A, A], None)]

    mock_run_model = MockRunModel(
        pytestconfig.rootpath.joinpath("test-data/complex_monomer")
    )
    mock_run_mmseqs2 = MMseqs2Mock(
        pytestconfig.rootpath, "complex_monomer"
    ).mock_run_mmseqs2
    with mock.patch(
        "alphafold.model.model.RunModel.predict",
        lambda model_runner, feat: mock_run_model.predict(model_runner, feat),
    ), mock.patch("colabfold.batch.run_mmseqs2", mock_run_mmseqs2):
        run(
            queries,
            tmp_path,
            use_templates=False,
            use_amber=False,
            msa_mode="MMseqs2 (UniRef+Environmental)",
            model_type="auto",
            num_models=1,
            num_recycles=3,
            model_order=[1, 2, 3, 4, 5],
            is_complex=True,
            keep_existing_results=False,
            rank_mode="auto",
            pair_mode="unpaired+paired",
            stop_at_score=100,
        )

    messages = list(caplog.messages)
    # noinspection PyUnresolvedReferences
    messages[3] = re.sub(r"\d+\.\d+s", "0.0s", messages[3])
    assert messages == [
        "Found 5 citations for tools or databases",
        "Query 1/1: A_A (length 118)",
        "Running model_1",
        "model_1 took 0.0s with pLDDT 95.3",
        "reranking models based on avg. predicted lDDT",
        "Done",
    ]
