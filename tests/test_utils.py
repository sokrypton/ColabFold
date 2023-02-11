import pytest

from colabfold.batch import get_queries, convert_pdb_to_mmcif, validate_and_fix_mmcif


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
        ("5AWL_1", "YYDPETGTWY", None),
        (
            "3G5O_A_3G5O_B",
            [
                "MRILPISTIKGKLNEFVDAVSSTQDQITITKNGAPAAVLVGADEWESLQETLYWLAQPGIRESIAEADADIASGRTYGEDEIRAEFGVPRRPH",
                "MPYTVRFTTTARRDLHKLPPRILAAVVEFAFGDLSREPLRVGKPLRRELAGTFSARRGTYRLLYRIDDEHTTVVILRVDHRADIYRR",
            ],
            None,
        ),
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
