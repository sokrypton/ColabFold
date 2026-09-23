"""Stand in for tokamax before any test imports alphafold3, as the CLI does."""
from colabfold.alphafold3.tokamax_shim import install

install()
