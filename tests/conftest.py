"""Skip optional-dependency-heavy test modules when their imports aren't
available, so the suite can still run a useful subset on a minimal install.

Without the optional ``alphafold`` extra (and its transitive deps biopython,
jax, haiku, absl, etc.) installed, pytest used to abort at collection time
because ``tests/test_msa.py``, ``tests/test_utils.py``, and
``tests/test_colabfold.py`` all import ``colabfold.batch`` at module load,
which in turn imports ``Bio`` and the rest.

This conftest probes the required imports up front and adds any unmet
modules to ``collect_ignore``, reporting what was skipped and why in the
pytest header. Install with ``poetry install -E alphafold`` (or the
equivalent pip extras) to enable the full suite.
"""

import importlib


# Module-level imports each test file performs at collection time. If any of
# these can't be imported, that test file is skipped from collection. We try
# the actual import rather than just locating the module because
# ``colabfold.batch`` raises ``RuntimeError`` (not ``ImportError``) when the
# alphafold extra is missing, which a spec lookup would not catch.
_OPTIONAL_REQUIREMENTS = {
    "test_msa.py": ["colabfold.batch"],
    "test_utils.py": ["colabfold.batch"],
    "test_colabfold.py": [
        "colabfold.batch",
        "alphafold.model.data",
        "haiku",
        "absl",
    ],
}


# Top-level package names whose absence we treat as "this optional extra was
# not installed". Anything else — a missing submodule of an installed
# package, or a non-ModuleNotFoundError ImportError — points at API drift or
# a real import-time bug and must surface, not be hidden behind a skip.
#
# Sources:
#   - alphafold extra in pyproject.toml: alphafold-colabfold (→ ``alphafold``),
#     jax, absl-py (→ ``absl``), dm-tree (→ ``tree``), dm-haiku (→ ``haiku``),
#     tensorflow / tensorflow-cpu (→ ``tensorflow``), py3Dmol.
#   - Base deps that are nevertheless commonly absent in minimal local venvs:
#     biopython (→ ``Bio``), importlib_metadata. We accept these so that
#     ``pytest`` is usable on a bare checkout without ``poetry install``.
_OPTIONAL_TOP_LEVEL_PACKAGES = frozenset({
    "Bio",
    "absl",
    "alphafold",
    "haiku",
    "importlib_metadata",
    "jax",
    "jaxlib",
    "py3Dmol",
    "tensorflow",
    "tree",
})

# colabfold.batch deliberately raises RuntimeError (not ModuleNotFoundError)
# when the alphafold extra is missing — see colabfold/batch.py. We treat that
# exact signal as a missing-extras failure; any other RuntimeError indicates
# a real bug that must surface, not be swallowed by collection.
_BATCH_MISSING_EXTRA_MARKER = "alphafold is not installed"


def _format_failure(exc):
    """Render an exception as a single, capped line for the pytest header."""
    message = " ".join(str(exc).split())
    if not message:
        return type(exc).__name__
    if len(message) > 160:
        message = message[:157] + "..."
    return f"{type(exc).__name__}: {message}"


def _import_failure(module_name):
    """Return a short failure description if importing fails with a known
    missing-extras signal, or ``None`` on success.

    Only two failure shapes are treated as "missing optional extra":

    1. ``ModuleNotFoundError`` whose ``name`` attribute is exactly one of the
       packages in :data:`_OPTIONAL_TOP_LEVEL_PACKAGES`. A submodule miss
       (``name == "alphafold.common"``) implies the root is installed but
       its layout doesn't match what the code expects — that's API drift.
    2. ``RuntimeError`` raised by ``colabfold.batch`` with the canonical
       "alphafold is not installed" marker.

    Every other exception — including plain ``ImportError``
    ("cannot import name 'X' from 'Y'"), ``AttributeError``, ``TypeError``,
    or unexpected ``RuntimeError`` — propagates so genuine bugs surface
    during collection instead of being misreported as a missing extra.
    """
    try:
        importlib.import_module(module_name)
        return None
    except ModuleNotFoundError as e:
        if e.name in _OPTIONAL_TOP_LEVEL_PACKAGES:
            return _format_failure(e)
        raise
    except RuntimeError as e:
        if (
            module_name == "colabfold.batch"
            and _BATCH_MISSING_EXTRA_MARKER in str(e)
        ):
            return _format_failure(e)
        raise


_missing_by_file = {}
collect_ignore = []

for _filename, _requirements in _OPTIONAL_REQUIREMENTS.items():
    _failures = {}
    for _module in _requirements:
        _reason = _import_failure(_module)
        if _reason is not None:
            _failures[_module] = _reason
    if _failures:
        collect_ignore.append(_filename)
        _missing_by_file[_filename] = _failures


def pytest_report_header(config):
    if not _missing_by_file:
        return None
    lines = [
        "Skipping test modules that depend on the optional `alphafold` extra "
        "(install with `poetry install -E alphafold`):"
    ]
    for filename in sorted(_missing_by_file):
        for module, reason in _missing_by_file[filename].items():
            lines.append(f"  tests/{filename}: cannot import {module} ({reason})")
    return "\n".join(lines)
