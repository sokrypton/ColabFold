import pytest

from colabfold.batch import get_queries, convert_pdb_to_mmcif, validate_and_fix_mmcif


def test_get_queries_fasta_dir(pytestconfig, caplog):
    dir_path = pytestconfig.rootpath.joinpath("test-data/batch/input")
    queries, is_complex = get_queries(dir_path)
    assert queries == [("5AWL_1", "YYDPETGTWY", None, None), ("6A5J", "IKKILSKIKKLLK", None, None)]
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
        ("5AWL_1", "YYDPETGTWY", None, None),
        (
            "3G5O_A_3G5O_B",
            [
                "MRILPISTIKGKLNEFVDAVSSTQDQITITKNGAPAAVLVGADEWESLQETLYWLAQPGIRESIAEADADIASGRTYGEDEIRAEFGVPRRPH",
                "MPYTVRFTTTARRDLHKLPPRILAAVVEFAFGDLSREPLRVGKPLRRELAGTFSARRGTYRLLYRIDDEHTTVVILRVDHRADIYRR",
            ],
            None,
            None,
        ),
    ]
    assert is_complex
    assert caplog.messages == []


def test_a3m_input(pytestconfig, caplog, tmp_path):
    queries, is_complex = get_queries(pytestconfig.rootpath.joinpath("test-data/a3m"))

    assert queries == [
        ("5AWL1", "YYDPETGTWY", [">101\nYYDPETGTWY"], None),
        ("6A5J", "IKKILSKIKKLLK", [">101\nIKKILSKIKKLLK\n>101\nIKKILSKIKKLLK"], None),
    ]
    assert not is_complex

    queries, is_complex = get_queries(
        pytestconfig.rootpath.joinpath("test-data/a3m/6A5J.a3m")
    )

    assert queries == [
        ("6A5J", "IKKILSKIKKLLK", [">101\nIKKILSKIKKLLK\n>101\nIKKILSKIKKLLK"], None)
    ]
    assert not is_complex

    assert caplog.messages == [
        f"{pytestconfig.rootpath}/test-data/a3m/empty.a3m is empty"
    ]


def test_convert_pdb_to_mmcif(pytestconfig, tmp_path):
    base_name = "ERR550519_2213899_unrelaxed_model_1"
    tmp_path.joinpath(f"{base_name}.pdb").write_text(
        pytestconfig.rootpath.joinpath(f"test-data/{base_name}.pdb").read_text()
    )

    convert_pdb_to_mmcif(tmp_path.joinpath(f"{base_name}.pdb"))

    actual = tmp_path.joinpath(f"{base_name}.cif").read_text()
    expected = pytestconfig.rootpath.joinpath(f"test-data/{base_name}.cif").read_text()
    assert actual == expected


def test_validate_and_fix_mmcif(pytestconfig, tmp_path):
    base_name = "ERR550519_2213899_unrelaxed_model_1"
    tmp_path.joinpath(f"{base_name}.cif").write_text(
        pytestconfig.rootpath.joinpath(f"test-data/{base_name}.cif").read_text()
    )

    validate_and_fix_mmcif(tmp_path.joinpath(f"{base_name}.cif"))

    actual = tmp_path.joinpath(f"{base_name}.cif").read_text()
    expected = pytestconfig.rootpath.joinpath(f"test-data/{base_name}.cif").read_text()
    assert actual == expected


def test_validate_and_fix_mmcif_add_revision_date(pytestconfig, tmp_path):
    base_name = "ERR550519_2213899_unrelaxed_model_1"
    cif = pytestconfig.rootpath.joinpath(f"test-data/{base_name}.cif").read_text()
    cif_without_revision_date = "\n".join(cif.split("\n")[:-9]) + "\n"

    tmp_path.joinpath(f"{base_name}.cif").write_text(cif_without_revision_date)

    validate_and_fix_mmcif(tmp_path.joinpath(f"{base_name}.cif"))

    actual = tmp_path.joinpath(f"{base_name}.cif").read_text()
    expected = pytestconfig.rootpath.joinpath(f"test-data/{base_name}.cif").read_text()
    assert actual == expected


def test_validate_and_fix_mmcif_missing_poly_seq(pytestconfig, tmp_path):
    base_name = "ERR550519_2213899_unrelaxed_model_1"
    cif_lines = (
        pytestconfig.rootpath.joinpath(f"test-data/{base_name}.cif")
        .read_text()
        .split("\n")
    )
    del cif_lines[2:138]
    cif_without_poly_seq = "\n".join(cif_lines) + "\n"

    cif_file = tmp_path.joinpath(f"{base_name}.cif")
    cif_file.write_text(cif_without_poly_seq)

    with pytest.raises(ValueError) as e:
        validate_and_fix_mmcif(cif_file)
        assert (
            e.message
            == f"mmCIF file {cif_file} is missing required field _entity_poly_seq.mon_id."
        )


def test_convert_pdb_to_mmcif_parse_with_alphafold(pytestconfig, tmp_path):
    base_name = "ERR550519_2213899_unrelaxed_model_1"
    tmp_path.joinpath(f"{base_name}.pdb").write_text(
        pytestconfig.rootpath.joinpath(f"test-data/{base_name}.pdb").read_text()
    )

    convert_pdb_to_mmcif(tmp_path.joinpath(f"{base_name}.pdb"))

    from alphafold.data.templates import mmcif_parsing

    parsing_result = mmcif_parsing.parse(
        file_id=base_name,
        mmcif_string=tmp_path.joinpath(f"{base_name}.cif").read_text(),
        catch_all_errors=False,
    )

    assert len(parsing_result.errors) == 0


def test_af3_merge_paired_msa_into_unpaired():
    from colabfold.utils import AF3Utils

    query_seqs_unique = ["AAAA", "GGGG"]
    query_seqs_cardinality = [1, 1]
    # Per-chain unpaired A3Ms (first entry is the query itself)
    unpaired_msa = [
        ">query\nAAAA\n>u1\nAAta\n",
        ">query\nGGGG\n>u2\nGGgg\n",
    ]
    # Per-chain paired A3Ms: row-aligned across chains, gap-only dummy row
    # for chains without a hit for that species
    paired_msa = [
        ">query\nAAAA\n>p1\nAA--\n>p2\n--AA\n",
        ">query\nGGGG\n>p1\n----\n>p2\nGGGG\n",
    ]

    # Without a paired MSA, nothing changes and pairedMsa stays empty
    af3 = AF3Utils("test", query_seqs_unique, query_seqs_cardinality, unpaired_msa, None)
    for i, seq in enumerate(af3.content["sequences"]):
        assert seq["protein"]["unpairedMsa"] == unpaired_msa[i]
        assert seq["protein"]["pairedMsa"] == ""

    # With a paired MSA, the paired rows are merged into unpairedMsa (paired
    # block first, query sequence as the first entry) and pairedMsa is left
    # empty so AF3 does not re-pair the MSA by header name matching
    af3 = AF3Utils("test", query_seqs_unique, query_seqs_cardinality, unpaired_msa, paired_msa)
    expected_unpaired = [
        paired_msa[0] + unpaired_msa[0],
        paired_msa[1] + unpaired_msa[1],
    ]
    for i, seq in enumerate(af3.content["sequences"]):
        assert seq["protein"]["unpairedMsa"] == expected_unpaired[i]
        # First MSA entry must be exactly the query sequence (AF3 requirement)
        a3m_lines = seq["protein"]["unpairedMsa"].splitlines()
        assert a3m_lines[1] == query_seqs_unique[i]
        assert seq["protein"]["pairedMsa"] == ""


def test_af3_merge_paired_msa_into_unpaired_paired_only():
    from colabfold.utils import AF3Utils

    query_seqs_unique = ["AAAA", "GGGG"]
    query_seqs_cardinality = [1, 1]
    paired_msa = [
        ">query\nAAAA\n>p1\nAA--\n",
        ">query\nGGGG\n>p1\n----\n",
    ]

    # paired-only mode: no unpaired MSA available
    af3 = AF3Utils("test", query_seqs_unique, query_seqs_cardinality, None, paired_msa)
    for i, seq in enumerate(af3.content["sequences"]):
        assert seq["protein"]["unpairedMsa"] == paired_msa[i]
        assert seq["protein"]["pairedMsa"] == ""
