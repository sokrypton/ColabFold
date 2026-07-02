"""
AlphaFold2 backend
Importing this subpackage requires the ``alphafold`` extra.
"""
try:
    import alphafold  # noqa: F401
except ModuleNotFoundError:
    raise RuntimeError(
        "\n\nalphafold is not installed. Please run `pip install colabfold[alphafold]`\n"
    )
