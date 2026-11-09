import shutil

import pytest

from colabfold.input import get_queries


STRUCTURE_SEQUENCE = (
    "CPSGWVGFSDNCYQWNSQKLTWSDAREFCNNHEADLVSIHSDLENNFISYQLSPKKQRFWIGTQRKAT"
    "GSTEWLWTDGSSFDYENFCQGKLDNWLGNQDCLDIGYTWGSQCSRDTWDDQPCDHQLPFVCK"
)


@pytest.mark.parametrize("suffix", [".pdb", ".cif", ".PDB", ".CIF"])
def test_get_queries_structure_dir(pytestconfig, tmp_path, suffix):
    source = pytestconfig.rootpath.joinpath(
        f"test-data/ERR550519_2213899_unrelaxed_model_1{suffix.lower()}"
    )
    shutil.copyfile(source, tmp_path.joinpath(f"structure{suffix}"))

    queries, is_complex = get_queries(tmp_path)

    assert queries == [("structure", [STRUCTURE_SEQUENCE], None, None)]
    assert is_complex


@pytest.mark.parametrize("suffix", [".pdb", ".cif"])
@pytest.mark.parametrize("fasta_first", [True, False])
def test_get_queries_mixed_structure_dir(pytestconfig, tmp_path, suffix, fasta_first):
    source = pytestconfig.rootpath.joinpath(
        f"test-data/ERR550519_2213899_unrelaxed_model_1{suffix}"
    )
    shutil.copyfile(source, tmp_path.joinpath(f"1_structure{suffix}"))
    fasta_name = "0_fasta" if fasta_first else "2_fasta"
    tmp_path.joinpath(f"{fasta_name}.fasta").write_text(">query\nACDEFG\n")

    queries, is_complex = get_queries(tmp_path)

    assert queries == [
        (fasta_name, "ACDEFG", None, None),
        ("1_structure", [STRUCTURE_SEQUENCE], None, None),
    ]
    assert is_complex
