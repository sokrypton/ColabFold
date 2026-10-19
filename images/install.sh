#!/bin/sh
# Shared by the Apptainer and Docker builds. Run as root inside the image.
set -eu
export DEBIAN_FRONTEND=noninteractive
export PIP_NO_CACHE_DIR=1
export POETRY_CACHE_DIR=/tmp/caat-poetry-cache
export POETRY_VIRTUALENVS_IN_PROJECT=true

apt-get update
apt-get install -y --no-install-recommends ca-certificates libgomp1
rm -rf /var/lib/apt/lists/*

python -m venv /opt/poetry
/opt/poetry/bin/pip install 'poetry==2.2.1'
cd /opt/CAAT
/opt/poetry/bin/poetry install --no-interaction -E alphafold

# Install this checkout's AlphaFold implementation, not a published substitute.
.venv/bin/python -m pip install --no-deps ./alphafold
# Match the CUDA plugins to the JAX version selected by poetry.lock.
CAAT_JAX_VERSION=$(.venv/bin/python -c 'from importlib.metadata import version; print(version("jax"))')
.venv/bin/python -m pip install "jax[cuda12]==${CAAT_JAX_VERSION}"
.venv/bin/python -m pip check
# Source checkouts on shared clusters can have group-private directories.
# The installed application must also be readable by arbitrary runtime UIDs.
chmod -R a+rX /opt/CAAT
chmod +x /opt/CAAT/images/entrypoint.sh
mkdir -p /work
rm -rf /opt/poetry "$POETRY_CACHE_DIR" /root/.cache
