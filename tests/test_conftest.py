"""Tests for ``tests/conftest.py``'s import-probe classifier.

The probe must distinguish two cases:

1. "Optional extra not installed" — expected, and the affected test module
   should be skipped from collection with a helpful header message.
2. "Optional extra is installed but importing it broke" — a real bug that
   must surface as a collection error, not be silently hidden behind the
   same "missing extra" message.

These tests pin that discrimination so a regression to a broader
``except Exception:`` / ``except ImportError:`` would fail loudly here.
"""

import pytest

from tests.conftest import _import_failure


@pytest.fixture
def fake_import(monkeypatch):
    """Replace ``importlib.import_module`` with a stub that raises ``exc``."""

    def install(exc):
        def _stub(name):
            raise exc

        monkeypatch.setattr("tests.conftest.importlib.import_module", _stub)

    return install


# -- happy path --------------------------------------------------------------


def test_returns_none_when_import_succeeds(monkeypatch):
    monkeypatch.setattr(
        "tests.conftest.importlib.import_module", lambda name: object()
    )
    assert _import_failure("colabfold.batch") is None


# -- expected missing-extras shapes -----------------------------------------


def test_module_not_found_for_known_optional_root_skips(fake_import):
    # Bare ``import haiku`` when dm-haiku isn't installed.
    fake_import(ModuleNotFoundError("No module named 'haiku'", name="haiku"))
    result = _import_failure("haiku")
    assert result is not None
    assert "ModuleNotFoundError" in result
    assert "haiku" in result


def test_module_not_found_for_biopython_skips(fake_import):
    # ``from Bio import ...`` failing in a bare minimal venv.
    fake_import(ModuleNotFoundError("No module named 'Bio'", name="Bio"))
    result = _import_failure("colabfold.batch")
    assert result is not None
    assert "ModuleNotFoundError" in result


def test_batch_alphafold_runtime_error_skips(fake_import):
    # Exact shape colabfold/batch.py raises when alphafold extra is missing.
    fake_import(
        RuntimeError(
            "\n\nalphafold is not installed. "
            "Please run `pip install colabfold[alphafold]`\n"
        )
    )
    result = _import_failure("colabfold.batch")
    assert result is not None
    assert "RuntimeError" in result
    assert "alphafold is not installed" in result


# -- adversarial cases: must NOT be hidden as missing-extras ---------------


def test_module_not_found_for_unknown_root_re_raises(fake_import):
    # A test-only helper module not in the optional set — must surface,
    # never be silently skipped under the missing-extras banner.
    err = ModuleNotFoundError(
        "No module named 'colabfold.unreleased_helper'",
        name="colabfold.unreleased_helper",
    )
    fake_import(err)
    with pytest.raises(ModuleNotFoundError):
        _import_failure("colabfold.unreleased_helper")


def test_module_not_found_for_submodule_of_known_root_re_raises(fake_import):
    # alphafold is installed but ``alphafold.common`` is somehow missing —
    # that's API drift / a layout problem, not a missing extra.
    err = ModuleNotFoundError(
        "No module named 'alphafold.common'", name="alphafold.common"
    )
    fake_import(err)
    with pytest.raises(ModuleNotFoundError):
        _import_failure("alphafold.model.data")


def test_import_error_cannot_import_name_re_raises(fake_import):
    # The reviewer's exact scenario: package installed, but a symbol is
    # missing because of version skew. Must surface — installing
    # ``-E alphafold`` should not produce silent test-collection skips.
    fake_import(
        ImportError(
            "cannot import name 'data' from 'alphafold.model' "
            "(/opt/venv/lib/.../alphafold/model/__init__.py)"
        )
    )
    with pytest.raises(ImportError):
        _import_failure("alphafold.model.data")


def test_runtime_error_unrelated_to_extras_re_raises(fake_import):
    # A genuine bug in batch.py top-level code — must NOT be hidden by the
    # batch-specific RuntimeError allowance.
    fake_import(RuntimeError("some unrelated runtime failure"))
    with pytest.raises(RuntimeError, match="some unrelated runtime failure"):
        _import_failure("colabfold.batch")


def test_runtime_error_alphafold_message_from_non_batch_re_raises(fake_import):
    # The batch-specific allowance must not extend to other modules, even
    # if someone elsewhere happens to raise a similarly-worded RuntimeError.
    fake_import(RuntimeError("alphafold is not installed"))
    with pytest.raises(RuntimeError):
        _import_failure("alphafold.model.data")


def test_attribute_error_re_raises(fake_import):
    # Classic "code bug in a module's top-level execution" — must surface.
    fake_import(AttributeError("module 'jax' has no attribute 'tree_util'"))
    with pytest.raises(AttributeError):
        _import_failure("colabfold.batch")


def test_type_error_re_raises(fake_import):
    fake_import(TypeError("metaclass conflict"))
    with pytest.raises(TypeError):
        _import_failure("colabfold.batch")
