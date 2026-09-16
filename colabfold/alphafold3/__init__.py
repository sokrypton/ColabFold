"""
AlphaFold3 backend, on alphafold3-colabfold (needs python >= 3.12).
"""
from pathlib import Path
from typing import Optional

INSTALL = "pip install alphafold3-colabfold"


def require(data_dir: Optional[Path] = None) -> None:
    """Raise unless alphafold3 is importable, fetching what it needs to import."""
    import importlib.util

    # find_spec does not run the package, so it answers before alphafold3.cpp
    # can raise over a components.cif that is not there yet
    if importlib.util.find_spec("alphafold3") is None:
        raise RuntimeError(
            f"\n\nalphafold3-colabfold is not installed. Install it with:\n\n"
            f"    {INSTALL}\n"
        )

    from colabfold.alphafold3.weights import ensure_ccd

    ensure_ccd(data_dir)
    import alphafold3  # noqa: F401
