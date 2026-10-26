from unittest import mock

from colabfold.batch import (
    generate_input_feature,
    get_msa_and_templates,
    unserialize_msa,
)
from tests.mock import MMseqs2Mock


def test_get_msa_and_templates(pytestconfig, caplog, tmp_path):
    Q60262 = "MEIIALLIEEGIIIIKDKKVAERFLKDLESSQGMDWKEIRERAERAKKQLEEGIEWAKKTKL"

    for msa_mode, tag, lines in [
        ("mmseqs2_uniref_env", "uniref_env", 12),
        ("mmseqs2_uniref", "uniref", 8),
        ("single_sequence", "single_sequence", 2),
    ]:
        mmseqs2mock = MMseqs2Mock(pytestconfig.rootpath, f"get_msa_{tag}")
        with mock.patch("colabfold.colabfold.run_mmseqs2", mmseqs2mock.mock_run_mmseqs2):
            (
                unpaired_msa,
                paired_msa,
                query_seqs_unique,
                query_seqs_cardinality,
                template_features,
            ) = get_msa_and_templates(
                "test",
                Q60262,
                None,
                tmp_path,
                msa_mode,
                False,
                None,
                "unpaired_paired",
                "colabfold/test"
            )

        assert len(unpaired_msa[0].splitlines()) == lines
        assert paired_msa is None
        assert query_seqs_unique == [Q60262]
        assert query_seqs_cardinality == [1]

    assert caplog.messages == []

def test_paired_only_custom_msa():
    """A custom a3m without unpaired sequences must still build features (#770)"""
    chain_a = "MEIIALLIEEGIIIIKDKK"
    chain_b = "VAERFLKDLESSQGMDWKEI"
    a3m_lines = "\n".join(
        [
            f"#{len(chain_a)},{len(chain_b)}\t1,1",
            ">101\t102",
            chain_a + chain_b,
            ">hit_1\thit_1",
            chain_a.replace("M", "L", 1) + chain_b.replace("V", "I", 1),
        ]
    )

    (
        unpaired_msa,
        paired_msa,
        query_seqs_unique,
        query_seqs_cardinality,
        template_features,
    ) = unserialize_msa([a3m_lines], [chain_a, chain_b])

    # there are no unpaired sequences in this a3m, so the unpaired msa is empty
    assert unpaired_msa == ["", ""]
    assert len(paired_msa) == 2

    input_feature, _ = generate_input_feature(
        query_seqs_unique,
        query_seqs_cardinality,
        unpaired_msa,
        paired_msa,
        template_features,
        True,
        "alphafold2_multimer_v3",
        512,
    )

    assert input_feature["msa"].shape[0] > 0
    assert input_feature["aatype"].shape[0] == len(chain_a) + len(chain_b)
