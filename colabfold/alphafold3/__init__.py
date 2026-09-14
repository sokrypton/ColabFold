"""
AlphaFold3 backend, on sokrypton's alphafold3-open (needs python >= 3.12).
"""
INSTALL = (
    "pip install 'alphafold3-open @ "
    "git+https://github.com/sokrypton/alphafold3.git@af3-any-model'"
)


def require() -> None:
    """Raise unless alphafold3-open is importable."""
    try:
        import alphafold3  # noqa: F401
    except ModuleNotFoundError:
        raise RuntimeError(
            f"\n\nalphafold3-open is not installed. It is not on PyPI, so install it with:\n\n"
            f"    {INSTALL}\n"
        )
