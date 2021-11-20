from pathlib import Path

import pytest

from colabfold.batch import get_queries, run


def test_get_queries_fasta_dir(pytestconfig, caplog):
    dir_path = pytestconfig.rootpath.joinpath("test-data/batch/input")
    queries, is_complex = get_queries(dir_path)
    assert queries == [("5AWL_1", "YYDPETGTWY", None), ("6A5J", "IKKILSKIKKLLK", None)]
    assert not is_complex
    assert caplog.messages == [f"{dir_path}/empty.fasta is empty"]


def test_get_queries_empty_a3m(pytestconfig, caplog):
    with pytest.raises(ValueError, match="a3m/empty.a3m is empty"):
        get_queries(pytestconfig.rootpath.joinpath("test-data/a3m/empty.a3m"))
    assert caplog.messages == []


def test_get_queries_csv(pytestconfig, caplog, tmp_path):
    queries, is_complex = get_queries(
        pytestconfig.rootpath.joinpath("test-data/complex/input.csv")
    )

    assert queries == [
        (
            "3G5O_A_3G5O_B",
            [
                "MRILPISTIKGKLNEFVDAVSSTQDQITITKNGAPAAVLVGADEWESLQETLYWLAQPGIRESIAEADADIASGRTYGEDEIRAEFGVPRRPH",
                "MPYTVRFTTTARRDLHKLPPRILAAVVEFAFGDLSREPLRVGKPLRRELAGTFSARRGTYRLLYRIDDEHTTVVILRVDHRADIYRR",
            ],
            None,
        ),
        ("5AWL_1", "YYDPETGTWY", None),
    ]
    assert is_complex
    assert caplog.messages == []


def test_a3m_input(pytestconfig, caplog, tmp_path):
    queries, is_complex = get_queries(pytestconfig.rootpath.joinpath("test-data/a3m"))

    assert queries == [
        ("5AWL1", "YYDPETGTWY", [">101\nYYDPETGTWY"]),
        ("6A5J", "IKKILSKIKKLLK", [">101\nIKKILSKIKKLLK\n>101\nIKKILSKIKKLLK"]),
    ]
    assert not is_complex

    queries, is_complex = get_queries(
        pytestconfig.rootpath.joinpath("test-data/a3m/6A5J.a3m")
    )

    assert queries == [
        ("6A5J", "IKKILSKIKKLLK", [">101\nIKKILSKIKKLLK\n>101\nIKKILSKIKKLLK"])
    ]
    assert not is_complex

    assert caplog.messages == [
        f"{pytestconfig.rootpath}/test-data/a3m/empty.a3m is empty"
    ]
