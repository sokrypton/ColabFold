import sys

import pytest

from colabfold.mmseqs import search


def test_use_templates_without_db2(tmp_path, monkeypatch):
    """--use-templates must use the template database of setup_databases.sh (#748)"""
    dbbase = tmp_path.joinpath("colabfold_db")
    dbbase.mkdir()
    for suffix in [".dbtype", ".idx"]:
        dbbase.joinpath(f"uniref30_2302_db{suffix}").touch()
    query = tmp_path.joinpath("query.fasta")
    query.write_text(">test\nMEIIALLIEEGIIIIKDKKVAERFLKDLESSQGMDWKEIRERAERAKKQLEEGIEWAKKTKL\n")

    monkeypatch.setattr(search, "run_mmseqs", lambda *args, **kwargs: None)
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "colabfold_search",
            str(query),
            str(dbbase),
            str(tmp_path.joinpath("out")),
            "--use-env",
            "0",
            "--use-templates",
            "1",
        ],
    )

    # the template database is not downloaded here, but the error has to name it
    with pytest.raises(FileNotFoundError, match="pdb100_230517"):
        search.main()
