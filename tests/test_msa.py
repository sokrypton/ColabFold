from unittest import mock

from colabfold.batch import get_msa_and_templates
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
                tmp_path,
                msa_mode,
                False,
                None,
                "unpaired_paired",
            )

        assert len(unpaired_msa[0].splitlines()) == lines
        assert paired_msa is None
        assert query_seqs_unique == [Q60262]
        assert query_seqs_cardinality == [1]

    assert caplog.messages == []