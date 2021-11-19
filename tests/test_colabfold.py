import logging
import re
from functools import lru_cache
from unittest import mock

import haiku
import pytest
import numpy
from absl import logging as absl_logging

from alphafold.model.data import get_model_haiku_params
from alphafold.model.tf import utils
from colabfold.batch import run
from colabfold.download import download_alphafold_params
from tests.mock import MockRunModel, MMseqs2Mock


# Without this, we're reading the params each time again which is slow
@lru_cache(maxsize=None)
def get_model_haiku_params_cached(model_name: str, data_dir: str) -> haiku.Params:
    return get_model_haiku_params(model_name, data_dir)

@pytest.fixture
def prediction_test(caplog):
    caplog.set_level(logging.INFO)
    # otherwise jax will tell us about its search for devices
    absl_logging.set_verbosity("error")
    # We'll also want to mock that out later
    download_alphafold_params("AlphaFold2-multimer")
    download_alphafold_params("AlphaFold2-ptm")
    # alphafold uses a method called `make_random_seed`, which deterministically starts with a seed
    # of zero and increases it by one for each protein. This means the input features would become
    # dependent on the number and order of tests. Here we just reset the seed to 0
    utils.seed_maker = utils.SeedMaker()

    # This works because it's used as `data.get_model_haiku_params`
    with mock.patch(
        "alphafold.model.data.get_model_haiku_params", get_model_haiku_params_cached
    ):
        yield



def test_batch(pytestconfig, caplog, tmp_path, prediction_test):
    queries = [("5AWL_1", "YYDPETGTWY", None), ("6A5J", "IKKILSKIKKLLK", None)]

    mock_run_model = MockRunModel(
        pytestconfig.rootpath.joinpath("test-data/batch"), ["5AWL_1", "6A5J"]
    )
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
        "model_1 took 0.0s with pLDDT 94.3",
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


def test_single_sequence(pytestconfig, caplog, tmp_path, prediction_test):
    queries = [("5AWL_1", "YYDPETGTWY", None)]

    mock_run_model = MockRunModel(
        pytestconfig.rootpath.joinpath("test-data/batch"), ["5AWL_1"]
    )
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
        "model_1 took 0.0s with pLDDT 94.3",
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


def test_complex(pytestconfig, caplog, tmp_path, prediction_test):
    pdb_3g50_A = "MRILPISTIKGKLNEFVDAVSSTQDQITITKNGAPAAVLVGADEWESLQETLYWLAQPGIRESIAEADADIASGRTYGEDEIRAEFGVPRRPH"
    pdb_3g50_B = "MPYTVRFTTTARRDLHKLPPRILAAVVEFAFGDLSREPLRVGKPLRRELAGTFSARRGTYRLLYRIDDEHTTVVILRVDHRADIYRR"
    queries = [("3G5O_A_3G5O_B", [pdb_3g50_A, pdb_3g50_B], None)]

    mock_run_model = MockRunModel(
        pytestconfig.rootpath.joinpath("test-data/complex"), ["3G5O_A_3G5O_B"]
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


def test_complex_ptm(pytestconfig, caplog, tmp_path, prediction_test):
    pdb_3g50_A = "MRILPISTIKGKLNEFVDAVSSTQDQITITKNGAPAAVLVGADEWESLQETLYWLAQPGIRESIAEADADIASGRTYGEDEIRAEFGVPRRPH"
    pdb_3g50_B = "MPYTVRFTTTARRDLHKLPPRILAAVVEFAFGDLSREPLRVGKPLRRELAGTFSARRGTYRLLYRIDDEHTTVVILRVDHRADIYRR"
    queries = [("3G5O_A_3G5O_B", [pdb_3g50_A, pdb_3g50_B], None)]

    mock_run_model = MockRunModel(
        pytestconfig.rootpath.joinpath("test-data/complex_ptm"), ["3G5O_A_3G5O_B"]
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
            model_type="AlphaFold2-ptm",
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
        "model_1 took 0.0s with pLDDT 91.9",
        "reranking models based on avg. predicted lDDT",
        "Done",
    ]


def test_complex_monomer_ptm(pytestconfig, caplog, tmp_path, prediction_test):
    A = "PIAQIHILEGRSDEQKETLIREVSEAISRSLDAPLTSVRVIITEMAKGHFGIGGELASK"
    queries = [("A_A", [A, A], None)]

    mock_run_model = MockRunModel(
        pytestconfig.rootpath.joinpath("test-data/complex_monomer_ptm"), ["A_A"]
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
            model_type="AlphaFold2-ptm",
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


def test_complex_monomer(pytestconfig, caplog, tmp_path, prediction_test):
    A = "PIAQIHILEGRSDEQKETLIREVSEAISRSLDAPLTSVRVIITEMAKGHFGIGGELASK"
    queries = [("A_A", [A, A], None)]

    mock_run_model = MockRunModel(
        pytestconfig.rootpath.joinpath("test-data/complex_monomer"), ["A_A"]
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


def test_msa_serialization(pytestconfig, caplog, tmp_path):
    from colabfold.batch import msa_to_str, unserialize_msa
    # heteromer
    unpaired_alignment = [">101\nAAAAAAAA\n>UP1\nAACCcccVVAA\n", ">102\nCCCC\n>UP1\nCCCC\n>UP2\nCaCaCC\n"]
    paired_alignment = [">101\nAAAAAAAA\n>UP1\nVVaVVAAAA\n", ">102\nCCCC\n>UP2\nGGGG\n"]
    query_sequence = ["AAAAAAAA", "AAAAAAAA", "CCCC"]
    query_sequence_unique = ["AAAAAAAA", "CCCC"]
    query_sequence_cardinality = [2,1]
    msa = msa_to_str(unpaired_alignment, paired_alignment, query_sequence_unique, query_sequence_cardinality)
    (unpaired_alignment_ret,
     paired_alignment_ret,
     query_sequence_unique_ret,
     query_sequence_cardinality_ret,
     template) = unserialize_msa(msa, query_sequence)
    numpy.testing.assert_equal(numpy.array(unpaired_alignment_ret), numpy.array(unpaired_alignment))
    numpy.testing.assert_equal(numpy.array(paired_alignment_ret), numpy.array(paired_alignment))
    numpy.testing.assert_equal(numpy.array(query_sequence_unique_ret), numpy.array(query_sequence_unique))
    numpy.testing.assert_equal(numpy.array(query_sequence_cardinality), numpy.array(query_sequence_cardinality_ret))

    # heteromer three complex
    unpaired_alignment = [">101\nAAAAAAAA\n>UP1\nAACCcccVVAA\n",
                          ">102\nCCCC\n>UP1\nCCCC\n>UP2\nCaCaCC\n",
                          ">103\nGGGG\n>UP1\nR--R\n",
                          ">104\nW\n"]
    paired_alignment = [">101\nAAAAAAAA\n>UP1\nVVaVVAAAA\n",
                        ">102\nCCCC\n>UP2\nGGGG\n",
                        ">103\nGGGG\n>UP3\nGGgGG\n",
                        ">104\nW\n>UP4\nW\n"]
    query_sequence = ["AAAAAAAA", "CCCC", "GGGG", "W", "W"]
    query_sequence_unique = ["AAAAAAAA", "CCCC", "GGGG", "W"]
    query_sequence_cardinality = [1, 1, 1, 2]
    msa = msa_to_str(unpaired_alignment, paired_alignment, query_sequence_unique, query_sequence_cardinality)
    (unpaired_alignment_ret,
     paired_alignment_ret,
     query_sequence_unique_ret,
     query_sequence_cardinality_ret,
     template) = unserialize_msa(msa, query_sequence)
    numpy.testing.assert_equal(numpy.array(unpaired_alignment_ret), numpy.array(unpaired_alignment))
    numpy.testing.assert_equal(numpy.array(paired_alignment_ret), numpy.array(paired_alignment))
    numpy.testing.assert_equal(numpy.array(query_sequence_unique_ret), numpy.array(query_sequence_unique))
    numpy.testing.assert_equal(numpy.array(query_sequence_cardinality), numpy.array(query_sequence_cardinality_ret))

    # heteromer with unpaired
    unpaired_alignment = [">101\nAAAAAAAA\n>UP1\nAACCcccVVAA\n", ">102\nCCCC\n>UP1\nCCCC\n>UP2\nCaCaCC\n"]
    paired_alignment = [">101\nAAAAAAAA\n", ">102\nCCCC\n"]
    query_sequence = ["AAAAAAAA", "CCCC", "CCCC"]
    query_sequence_unique = ["AAAAAAAA", "CCCC"]
    query_sequence_cardinality = [1, 2]
    msa = msa_to_str(unpaired_alignment, paired_alignment, query_sequence_unique, query_sequence_cardinality)
    (unpaired_alignment_ret,
     paired_alignment_ret,
     query_sequence_unique_ret,
     query_sequence_cardinality_ret,
     template) = unserialize_msa(msa, query_sequence)
    numpy.testing.assert_equal(numpy.array(unpaired_alignment_ret), numpy.array(unpaired_alignment))
    numpy.testing.assert_equal(numpy.array(paired_alignment_ret), numpy.array(paired_alignment))
    numpy.testing.assert_equal(numpy.array(query_sequence_unique_ret), numpy.array(query_sequence_unique))
    numpy.testing.assert_equal(numpy.array(query_sequence_cardinality), numpy.array(query_sequence_cardinality_ret))

    # homooligomer
    unpaired_alignment = [">101\nAAAAAAAA\n>UP2\nAAAVVAAA\n>UP1\nA-CCcccVV-A\n"]
    paired_alignment = [">101\nAAAAAAAA\n", ">102\nAAAAAAAA\n"]
    query_sequence = ["AAAAAAAA", "AAAAAAAA"]
    query_sequence_unique = ["AAAAAAAA"]
    query_sequence_cardinality = [2]
    msa = msa_to_str(unpaired_alignment, paired_alignment, query_sequence_unique, query_sequence_cardinality)
    (unpaired_alignment_ret,
     paired_alignment_ret,
     query_sequence_unique_ret,
     query_sequence_cardinality_ret,
     template) = unserialize_msa(msa, query_sequence)

    numpy.testing.assert_equal(numpy.array(unpaired_alignment_ret), numpy.array(unpaired_alignment))
    numpy.testing.assert_equal(numpy.array(paired_alignment_ret), numpy.array(paired_alignment))
    numpy.testing.assert_equal(numpy.array(query_sequence_unique_ret), numpy.array(query_sequence_unique))
    numpy.testing.assert_equal(numpy.array(query_sequence_cardinality), numpy.array(query_sequence_cardinality_ret))

    # a3m without header
    unpaired_alignment = ">101\nAAAAAAAA\n>UP2\nAAAVVAAA\n>UP1\nA-CCcccVV-A\n"
    paired_alignment = None
    query_sequence = ["AAAAAAAA"]
    query_sequence_unique = ["AAAAAAAA"]
    query_sequence_cardinality = [1]
    (unpaired_alignment_ret,
     paired_alignment_ret,
     query_sequence_unique_ret,
     query_sequence_cardinality_ret,
     template) = unserialize_msa(unpaired_alignment, query_sequence)

    numpy.testing.assert_equal(numpy.array(unpaired_alignment_ret), numpy.array(unpaired_alignment))
    numpy.testing.assert_equal(numpy.array(paired_alignment_ret), numpy.array(paired_alignment))
    numpy.testing.assert_equal(numpy.array(query_sequence_unique_ret), numpy.array(query_sequence_unique))
    numpy.testing.assert_equal(numpy.array(query_sequence_cardinality), numpy.array(query_sequence_cardinality_ret))

    msa = "#10\t1\n>101\nYYDPETGTWY"
    (unpaired_alignment_ret,
     paired_alignment_ret,
     query_sequence_unique_ret,
     query_sequence_cardinality_ret,
     template) = unserialize_msa(msa, "YYDPETGTWY")
    numpy.testing.assert_equal(numpy.array(unpaired_alignment_ret), numpy.array(['>101\nYYDPETGTWY\n']))
    assert paired_alignment_ret == None
    numpy.testing.assert_equal(numpy.array(query_sequence_unique_ret), numpy.array(["YYDPETGTWY"]))
    numpy.testing.assert_equal(numpy.array(query_sequence_cardinality), numpy.array([1]))

