from pathlib import Path

from colabfold.mmseqs import search


def _touch_db(dbbase: Path, name: str) -> None:
    dbbase.joinpath(f"{name}.dbtype").touch()
    dbbase.joinpath(f"{name}.idx").touch()


def _capture_mmseqs_calls(monkeypatch):
    calls = []

    def fake_run_mmseqs(mmseqs, params):
        calls.append([str(param) for param in params])

    monkeypatch.setattr(search, "run_mmseqs", fake_run_mmseqs)
    return calls


def _assert_default_k_score_is_unquoted(calls):
    search_calls = [call for call in calls if call[0] == "search"]
    assert len(search_calls) == 1

    search_call = search_calls[0]
    assert "'seq:96,prof:80'" not in search_call
    assert search_call[search_call.index("--k-score") + 1] == "seq:96,prof:80"


def test_monomer_search_default_k_score_is_passed_without_shell_quotes(
    monkeypatch, tmp_path
):
    dbbase = tmp_path.joinpath("db")
    base = tmp_path.joinpath("work")
    dbbase.mkdir()
    base.mkdir()
    base.joinpath("tmp").mkdir()
    _touch_db(dbbase, "uniref30_2302_db")
    calls = _capture_mmseqs_calls(monkeypatch)

    search.mmseqs_search_monomer(
        dbbase=dbbase,
        base=base,
        use_env=False,
        use_templates=False,
        s=None,
        threads=1,
        unpack=False,
    )

    _assert_default_k_score_is_unquoted(calls)


def test_pair_search_default_k_score_is_passed_without_shell_quotes(
    monkeypatch, tmp_path
):
    dbbase = tmp_path.joinpath("db")
    base = tmp_path.joinpath("work")
    dbbase.mkdir()
    base.mkdir()
    base.joinpath("tmp").mkdir()
    _touch_db(dbbase, "uniref30_2302_db")
    calls = _capture_mmseqs_calls(monkeypatch)

    search.mmseqs_search_pair(
        dbbase=dbbase,
        base=base,
        pair_env=False,
        s=None,
        threads=1,
        unpack=False,
    )

    _assert_default_k_score_is_unquoted(calls)
