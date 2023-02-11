from unittest import mock

import haiku
import logging
import pytest
import re
from absl import logging as absl_logging
from functools import lru_cache
from zipfile import ZipFile

from alphafold.model.data import get_model_haiku_params
from alphafold.model.tf import utils
from colabfold.batch import msa_to_str, unserialize_msa, get_queries
from colabfold.batch import run
from colabfold.download import download_alphafold_params
from tests.mock import MockRunModel, MMseqs2Mock


# Without this, we're reading the params each time again which is slow
@lru_cache(maxsize=None)
def get_model_haiku_params_cached(model_name: str, data_dir: str, fuse: bool = True) -> haiku.Params:
    return get_model_haiku_params(model_name, data_dir, fuse)


@pytest.fixture
def prediction_test(caplog):
    caplog.set_level(logging.INFO)
    # otherwise jax will tell us about its search for devices
    absl_logging.set_verbosity("error")
    # We'll also want to mock that out later
    download_alphafold_params("alphafold2_multimer_v1")
    download_alphafold_params("alphafold2_ptm")
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
        lambda model_runner, feat, random_seed, return_representations, callback: \
        mock_run_model.predict(model_runner, feat, random_seed, return_representations, callback),
    ), mock.patch("colabfold.colabfold.run_mmseqs2", mock_run_mmseqs):
        run(
            queries,
            tmp_path,
            num_models=1,
            num_recycles=3,
            model_order=[1, 2, 3, 4, 5],
            is_complex=False,
        )

    messages = [re.sub(r"\d+\.\d+s", "0.0s", i) for i in caplog.messages]
    expected = [
      #'Running on GPU',
      'Found 4 citations for tools or databases',
      'Query 1/2: 5AWL_1 (length 10)',
      'Padding length to 13',
      'alphafold2_ptm_model_1_seed_000 took 0.0s (3 recycles)',
      "reranking models by 'plddt' metric",
      'rank_001_alphafold2_ptm_model_1_seed_000 pLDDT=94.2 pTM=0.0567',
      'Query 2/2: 6A5J (length 13)',
      'alphafold2_ptm_model_1_seed_000 took 0.0s (3 recycles)',
      "reranking models by 'plddt' metric",
      'rank_001_alphafold2_ptm_model_1_seed_000 pLDDT=90.8 pTM=0.0455',
      'Done'
    ]
    for x in expected:
      assert x in messages

    # Very simple test, it would be better to check coordinates
    assert (
        len(
            tmp_path.joinpath("5AWL_1_unrelaxed_rank_001_alphafold2_ptm_model_1_seed_000.pdb")
            .read_text()
            .splitlines()
        )
        == 96
    )
    assert (
        len(
            tmp_path.joinpath("6A5J_unrelaxed_rank_001_alphafold2_ptm_model_1_seed_000.pdb")
            .read_text()
            .splitlines()
        )
        == 112
    )
    assert tmp_path.joinpath("config.json").is_file()


def test_zip(pytestconfig, caplog, tmp_path, prediction_test):
    queries = [("5AWL_1", "YYDPETGTWY", None), ("6A5J", "IKKILSKIKKLLK", None)]

    mock_run_model = MockRunModel(
        pytestconfig.rootpath.joinpath("test-data/batch"), ["5AWL_1", "6A5J"]
    )
    mock_run_mmseqs = MMseqs2Mock(pytestconfig.rootpath, "batch").mock_run_mmseqs2
    with mock.patch(
        "alphafold.model.model.RunModel.predict",
        lambda model_runner, feat, random_seed, return_representations, callback: \
        mock_run_model.predict(model_runner, feat, random_seed, return_representations, callback),
    ), mock.patch("colabfold.colabfold.run_mmseqs2", mock_run_mmseqs):
        run(
            queries,
            tmp_path,
            num_models=1,
            num_recycles=3,
            model_order=[1, 2, 3, 4, 5],
            is_complex=False,
            zip_results=True,
        )

    # Ensure that the correct files are packaged and that they do not contain the dir prefix
    expect_zip = [
      '5AWL_1_unrelaxed_rank_001_alphafold2_ptm_model_1_seed_000.pdb', 
      '5AWL_1_scores_rank_001_alphafold2_ptm_model_1_seed_000.json', 
      '5AWL_1_coverage.png',
      '5AWL_1_predicted_aligned_error_v1.json', 
      '5AWL_1_pae.png', 
      '5AWL_1_plddt.png', 
      '5AWL_1.a3m', 
      'cite.bibtex', 
      'config.json'
    ]
    with ZipFile(tmp_path.joinpath("5AWL_1.result.zip")) as result_zip:
        actual_zip = [i.filename for i in result_zip.infolist()]
    
    for x in expect_zip:
      assert x in actual_zip

def test_single_sequence(pytestconfig, caplog, tmp_path, prediction_test):
    queries = [("5AWL_1", "YYDPETGTWY", None)]

    mock_run_model = MockRunModel(
        pytestconfig.rootpath.joinpath("test-data/single"), ["5AWL_1"]
    )
    mock_run_mmseqs = MMseqs2Mock(pytestconfig.rootpath, "single").mock_run_mmseqs2
    with mock.patch(
        "alphafold.model.model.RunModel.predict",
        lambda model_runner, feat, random_seed, return_representations, callback: \
        mock_run_model.predict(model_runner, feat, random_seed, return_representations, callback),
    ), mock.patch("colabfold.colabfold.run_mmseqs2", mock_run_mmseqs):
        run(
            queries,
            tmp_path,
            msa_mode="single_sequence",
            num_models=1,
            num_recycles=3,
            model_order=[1, 2, 3, 4, 5],
            is_complex=False,
            stop_at_score=100,
        )

    messages = [re.sub(r"\d+\.\d+s", "0.0s", i) for i in caplog.messages]
    expected = [
      #'Running on GPU',
      'Found 1 citations for tools or databases', 
      'Query 1/1: 5AWL_1 (length 10)', 
      'alphafold2_ptm_model_1_seed_000 took 0.0s (3 recycles)',
      "reranking models by 'plddt' metric",
      'rank_001_alphafold2_ptm_model_1_seed_000 pLDDT=94.2 pTM=0.0567',
      'Done']

    for x in expected:
      assert x in messages

    # Very simple test, it would be better to check coordinates
    assert (
        len(
            tmp_path.joinpath("5AWL_1_unrelaxed_rank_001_alphafold2_ptm_model_1_seed_000.pdb")
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
        lambda model_runner, feat, random_seed, return_representations, callback: \
        mock_run_model.predict(model_runner, feat, random_seed, return_representations, callback),
    ), mock.patch("colabfold.colabfold.run_mmseqs2", mock_run_mmseqs2):
        run(
            queries,
            tmp_path,
            num_models=1,
            model_type="alphafold2_multimer_v1",
            num_recycles=3,
            model_order=[1, 2, 3, 4, 5],
            is_complex=True,
            stop_at_score=100,
        )

    messages = [re.sub(r"\d+\.\d+s", "0.0s", i) for i in caplog.messages]
    expected = [
      #'Running on GPU',
      'Found 4 citations for tools or databases',
      'Query 1/1: 3G5O_A_3G5O_B (length 180)',
      'Setting max_seq=252, max_extra_seq=1152',
      'alphafold2_multimer_v1_model_1_seed_000 took 0.0s (3 recycles)',
      "reranking models by 'multimer' metric",
      'rank_001_alphafold2_multimer_v1_model_1_seed_000 pLDDT=94.4 pTM=0.884 ipTM=0.878',
      'Done'
    ]
    for x in expected:
      assert x in messages

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
        lambda model_runner, feat, random_seed, return_representations, callback: \
        mock_run_model.predict(model_runner, feat, random_seed, return_representations, callback),
    ), mock.patch("colabfold.colabfold.run_mmseqs2", mock_run_mmseqs2):
        run(
            queries,
            tmp_path,
            model_type="alphafold2_ptm",
            num_models=1,
            num_recycles=3,
            model_order=[1, 2, 3, 4, 5],
            is_complex=True,
            stop_at_score=100,
        )

    messages = [re.sub(r"\d+\.\d+s", "0.0s", i) for i in caplog.messages]
    expected = [
      #'Running on GPU',
      'Found 4 citations for tools or databases',
      'Query 1/1: 3G5O_A_3G5O_B (length 180)',
      'Setting max_seq=512, max_extra_seq=5120',
      'alphafold2_ptm_model_1_seed_000 took 0.0s (3 recycles)',
      "reranking models by 'multimer' metric",
      'rank_001_alphafold2_ptm_model_1_seed_000 pLDDT=92 pTM=0.846 ipTM=0.849',
      'Done'
    ]
    for x in expected:
      assert x in messages

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
        lambda model_runner, feat, random_seed, return_representations, callback: \
        mock_run_model.predict(model_runner, feat, random_seed, return_representations, callback),
    ), mock.patch("colabfold.colabfold.run_mmseqs2", mock_run_mmseqs2):
        run(
            queries,
            tmp_path,
            model_type="alphafold2_ptm",
            num_models=1,
            num_recycles=3,
            model_order=[1, 2, 3, 4, 5],
            is_complex=True,
            stop_at_score=100,
        )

    messages = [re.sub(r"\d+\.\d+s", "0.0s", i) for i in caplog.messages]
    expected = [
      #'Running on GPU', 
      'Found 4 citations for tools or databases', 
      'Query 1/1: A_A (length 118)', 
      'Setting max_seq=512, max_extra_seq=5120', 
      'alphafold2_ptm_model_1_seed_000 took 0.0s (3 recycles)', 
      "reranking models by 'multimer' metric", 
      'rank_001_alphafold2_ptm_model_1_seed_000 pLDDT=95.6 pTM=0.867 ipTM=0.864', 
      'Done'
    ]
    for x in expected:
      assert x in messages

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
        lambda model_runner, feat, random_seed, return_representations, callback: \
        mock_run_model.predict(model_runner, feat, random_seed, return_representations, callback),
    ), mock.patch("colabfold.colabfold.run_mmseqs2", mock_run_mmseqs2):
        run(
            queries,
            tmp_path,
            num_models=1,
            model_type="alphafold2_multimer_v1",
            num_recycles=3,
            model_order=[1, 2, 3, 4, 5],
            is_complex=True,
            stop_at_score=100,
        )

    messages = [re.sub(r"\d+\.\d+s", "0.0s", i) for i in caplog.messages]
    expected = [
      #'Running on GPU', 
      'Found 4 citations for tools or databases', 
      'Query 1/1: A_A (length 118)', 
      'Setting max_seq=252, max_extra_seq=1152', 
      'alphafold2_multimer_v1_model_1_seed_000 took 0.0s (3 recycles)', 
      "reranking models by 'multimer' metric", 
      'rank_001_alphafold2_multimer_v1_model_1_seed_000 pLDDT=95.3 pTM=0.866 ipTM=0.861', 
      'Done'
    ]
    for x in expected:
      assert x in messages

def test_msa_serialization(pytestconfig):
    # heteromer
    unpaired_alignment = [
        ">101\nAAAAAAAA\n>UP1\nAACCcccVVAA\n",
        ">102\nCCCC\n>UP1\nCCCC\n>UP2\nCaCaCC\n",
    ]
    paired_alignment = [">101\nAAAAAAAA\n>UP1\nVVaVVAAAA\n", ">102\nCCCC\n>UP2\nGGGG\n"]
    query_sequence = ["AAAAAAAA", "AAAAAAAA", "CCCC"]
    query_sequence_unique = ["AAAAAAAA", "CCCC"]
    query_sequence_cardinality = [2, 1]
    msa = msa_to_str(
        unpaired_alignment,
        paired_alignment,
        query_sequence_unique,
        query_sequence_cardinality,
    )
    (
        unpaired_alignment_ret,
        paired_alignment_ret,
        query_sequence_unique_ret,
        query_sequence_cardinality_ret,
        template,
    ) = unserialize_msa([msa], query_sequence)
    assert unpaired_alignment_ret == unpaired_alignment
    assert paired_alignment_ret == paired_alignment
    assert query_sequence_unique_ret == query_sequence_unique
    assert query_sequence_cardinality == query_sequence_cardinality_ret

    # heteromer three complex
    unpaired_alignment = [
        ">101\nAAAAAAAA\n>UP1\nAACCcccVVAA\n",
        ">102\nCCCC\n>UP1\nCCCC\n>UP2\nCaCaCC\n",
        ">103\nGGGG\n>UP1\nR--R\n",
        ">104\nW\n",
    ]
    paired_alignment = [
        ">101\nAAAAAAAA\n>UP1\nVVaVVAAAA\n",
        ">102\nCCCC\n>UP2\nGGGG\n",
        ">103\nGGGG\n>UP3\nGGgGG\n",
        ">104\nW\n>UP4\nW\n",
    ]
    query_sequence = ["AAAAAAAA", "CCCC", "GGGG", "W", "W"]
    query_sequence_unique = ["AAAAAAAA", "CCCC", "GGGG", "W"]
    query_sequence_cardinality = [1, 1, 1, 2]
    msa = msa_to_str(
        unpaired_alignment,
        paired_alignment,
        query_sequence_unique,
        query_sequence_cardinality,
    )
    (
        unpaired_alignment_ret,
        paired_alignment_ret,
        query_sequence_unique_ret,
        query_sequence_cardinality_ret,
        template,
    ) = unserialize_msa([msa], query_sequence)
    assert unpaired_alignment_ret == unpaired_alignment
    assert paired_alignment_ret == paired_alignment
    assert query_sequence_unique_ret == query_sequence_unique
    assert query_sequence_cardinality == query_sequence_cardinality_ret

    # heteromer with unpaired
    unpaired_alignment = [
        ">101\nAAAAAAAA\n>UP1\nAACCcccVVAA\n",
        ">102\nCCCC\n>UP1\nCCCC\n>UP2\nCaCaCC\n",
    ]
    paired_alignment = [">101\nAAAAAAAA\n", ">102\nCCCC\n"]
    query_sequence = ["AAAAAAAA", "CCCC", "CCCC"]
    query_sequence_unique = ["AAAAAAAA", "CCCC"]
    query_sequence_cardinality = [1, 2]
    msa = msa_to_str(
        unpaired_alignment,
        paired_alignment,
        query_sequence_unique,
        query_sequence_cardinality,
    )
    (
        unpaired_alignment_ret,
        paired_alignment_ret,
        query_sequence_unique_ret,
        query_sequence_cardinality_ret,
        template,
    ) = unserialize_msa([msa], query_sequence)
    assert unpaired_alignment_ret == unpaired_alignment
    assert paired_alignment_ret == paired_alignment
    assert query_sequence_unique_ret == query_sequence_unique
    assert query_sequence_cardinality == query_sequence_cardinality_ret

    # homooligomer
    unpaired_alignment = [">101\nAAAAAAAA\n>UP2\nAAAVVAAA\n>UP1\nA-CCcccVV-A\n"]
    paired_alignment = [">101\nAAAAAAAA\n", ">102\nAAAAAAAA\n"]
    query_sequence = ["AAAAAAAA", "AAAAAAAA"]
    query_sequence_unique = ["AAAAAAAA"]
    query_sequence_cardinality = [2]
    msa = msa_to_str(
        unpaired_alignment,
        paired_alignment,
        query_sequence_unique,
        query_sequence_cardinality,
    )
    (
        unpaired_alignment_ret,
        paired_alignment_ret,
        query_sequence_unique_ret,
        query_sequence_cardinality_ret,
        template,
    ) = unserialize_msa([msa], query_sequence)

    assert unpaired_alignment_ret == unpaired_alignment
    assert paired_alignment_ret == paired_alignment
    assert query_sequence_unique_ret == query_sequence_unique
    assert query_sequence_cardinality == query_sequence_cardinality_ret

    # a3m without header
    unpaired_alignment = ">101\nAAAAAAAA\n>UP2\nAAAVVAAA\n>UP1\nA-CCcccVV-A"
    paired_alignment = None
    query_sequence = "AAAAAAAA"
    query_sequence_unique = ["AAAAAAAA"]
    query_sequence_cardinality = [1]
    (
        unpaired_alignment_ret,
        paired_alignment_ret,
        query_sequence_unique_ret,
        query_sequence_cardinality_ret,
        template,
    ) = unserialize_msa([unpaired_alignment], query_sequence)

    assert unpaired_alignment_ret == [unpaired_alignment]
    assert paired_alignment_ret is None
    assert query_sequence_unique_ret == query_sequence_unique
    assert query_sequence_cardinality == query_sequence_cardinality_ret

    msa = "#10\t1\n>101\nYYDPETGTWY"
    (
        unpaired_alignment_ret,
        paired_alignment_ret,
        query_sequence_unique_ret,
        query_sequence_cardinality_ret,
        template,
    ) = unserialize_msa([msa], "YYDPETGTWY")
    assert unpaired_alignment_ret == [">101\nYYDPETGTWY\n"]
    assert paired_alignment_ret is None
    assert query_sequence_unique_ret == ["YYDPETGTWY"]
    assert query_sequence_cardinality == [1]

    # non-complex a3m files
    a3m_file = pytestconfig.rootpath.joinpath("test-data/a3m/5AWL1.a3m")
    [(_, query_sequence, _)], is_complex = get_queries(a3m_file)
    assert not is_complex
    msa = a3m_file.read_text()
    (
        unpaired_alignment_ret,
        paired_alignment_ret,
        query_sequence_unique_ret,
        query_sequence_cardinality_ret,
        template,
    ) = unserialize_msa([msa], query_sequence)
    assert unpaired_alignment_ret
    assert not paired_alignment
    assert query_sequence_unique_ret == [query_sequence]
    assert query_sequence_cardinality_ret == [1]
