# CAAT

**CAAT** is an end-to-end pipeline designed to bridge the gap between protein structure prediction and functional biological insight. It automates the generation of AlphaFold2 structures and performs deep-dive analysis on the raw attention heads to identify residues of high structural importance to the model.

## Using CAAT via Colab Notebook

The quickest way to get up and running with CAAT is to use the publicly available notebook with the GPU runtime. You can find that [here](https://colab.research.google.com/drive/11UVndoYaP5cQD7762o8rT0DRlyzSnTnn?usp=sharing).

## CAAT Docs
For an in-depth look into CAAT, you can access the docs [here](https://prameshsharma25.github.io/CAAT/).

## Coverage Report

View coverage report [here](https://app.codecov.io/github/prameshsharma25/caat/tree/main).

## Quick Start with Poetry

Follow these steps to set up and run your local version of ColabFold using **Poetry**.

### Prerequisites

You need **Python 3.11 or newer**. The setup script installs Poetry into the
project environment, or you can install Poetry yourself for the manual setup.

* **Python:** Install [Python 3.11 or newer](https://www.python.org/downloads/).
* **Poetry:** Install it by following the official [Poetry installation guide](https://python-poetry.org/docs/#installation).
---

## Running Custom Attention Head Analysis Pipeline

This section guides you through installing the environment and executing the custom script for collecting attention heads.

## 1. Installation

To set up the project, you must first clone the repository and then install all dependencies using **Poetry**.

1.  **Clone the Repository:** Get the project source code from the remote repository.
    ```bash
    git clone https://github.com/prameshsharma25/CAAT.git
    cd CAAT
    ```

2.  **Automated setup:** Install CPU dependencies with:
    ```bash
    bash setup_env.sh
    ```

    For an NVIDIA GPU with CUDA 12, use:
    ```bash
    bash setup_env.sh --gpu
    ```

    The bundled AlphaFold package is installed normally into `.venv`; the
    source checkout is not linked in editable mode.

3.  **Manual setup (optional):** If Poetry is already installed, run:
    ```bash
    poetry install -E alphafold
    poetry run python -m pip install ./alphafold
    ```

    For CUDA 12 support, also run:
    ```bash
    poetry run python -m pip install --upgrade 'jax[cuda12]>=0.5.2,<0.11'
    ```

#### JAX version policy

CAAT no longer pins JAX to one patch release. It accepts
`jax>=0.5.2,<0.11`; `poetry.lock` records the concrete version tested by this
checkout. Install the CUDA extra as one package specification so `jax`,
`jaxlib`, `jax-cuda12-plugin`, and `jax-cuda12-pjrt` remain compatible. Do not
install those components at unrelated versions. See the official JAX
[installation guide](https://docs.jax.dev/en/latest/installation.html) for CUDA
driver and platform requirements.

Verify the resolved version and accelerator before a run:

```bash
poetry run python -c "import jax; print(jax.__version__); print(jax.devices())"
```
---

## Running the container

Both container recipes package the local checkout, including its bundled
AlphaFold code and CAAT CLI options. Build from the repository root:

```bash
sudo apptainer build caat.sif images/caat.def
```

Run from a directory containing your input. Mount that directory for inputs,
outputs, and a persistent model-weight cache:

```bash
apptainer run --cleanenv --nv --bind "$PWD:/work" --pwd /work caat.sif \
  scripts/run_e2e_pipeline.py \
  --query-seq-path query.a3m --query-name query \
  --num-recycles 3 --data-dir /work/weights
```

The entrypoint preserves the working directory. Relative output paths are
written there. `--data-dir` (alias `--data`) selects the directory containing
`params/`; downloads and prediction use the same directory. If omitted,
ColabFold's normal user cache is used. When binding existing weights, mount
the actual directory rather than a parent containing a symlink to an unmounted
location. Model weights are downloaded on demand and are not bundled in the image.

For analysis on a CPU, omit `--nv` and add `--env JAX_PLATFORMS=cpu`:

```bash
apptainer run --cleanenv --env JAX_PLATFORMS=cpu \
  --bind "$PWD:/work" --pwd /work caat.sif \
  scripts/run_analysis_pipeline.py --query-seq-path query.a3m \
  --query-name query --query-attn-dir attention_outputs/query
```

The Docker image uses the same installer and entrypoint. It defaults to the
non-root `caat` user (UID/GID `10001:10001`), with writable `/work` and
`/home/caat` directories. For bind-mounted inputs and outputs, override the user
to match your host UID/GID and put caches in the mounted working directory:

```bash
docker build -t caat .
docker run --rm --gpus all --user "$(id -u):$(id -g)" \
  -v "$PWD:/work" -e HOME=/work -e XDG_CACHE_HOME=/work/.cache \
  -e MPLCONFIGDIR=/work/.mplconfig caat \
  scripts/run_e2e_pipeline.py --query-seq-path query.a3m --query-name query \
  --num-recycles 3 --data-dir /work/weights
```

Without `--data-dir`, the container entrypoint uses the home cache when writable,
then `.cache/colabfold` in the working directory. If neither is writable, it
creates a temporary cache and prints its location. Set `--data-dir` for persistent
weights, or `XDG_CACHE_HOME` to explicitly choose the cache root; explicit cache
settings are preserved.

Container builds run an offline smoke test covering imports, CLI options,
relative paths, and real attention analysis. Run it again with
`apptainer test --cleanenv caat.sif`. CI also runs the packaged CAAT tests
that do not require model weights before publishing; GPU inference still
requires a CUDA-compatible host and `--nv` (Apptainer) or `--gpus all` (Docker).
The CUDA plugin versions match the JAX version selected by `poetry.lock`.

Published images use `sha-<full-commit-SHA>-sif-<SIF-SHA256>` tags. The publishing
job summary includes the exact URI and file checksum. Different SIF contents
receive different tags, including rebuilds of the same commit; the workflow no
longer writes short commit-only tags. `latest` remains a moving tag for `main`.
Base images and system packages can change between builds, so rebuilding a
commit is not guaranteed to reproduce identical bytes. Use the recorded image
URI to retrieve a particular build and verify the downloaded SIF with
`sha256sum` against the recorded checksum.


## 2. Running the Pipeline

CAAT offers three entry points depending on whether you need a full structural run or just specific analysis components.

**Option A: Full End-to-End Run**

```bash
poetry run python scripts/run_e2e_pipeline.py [OPTIONS]
```

**Option B: Generate Attention Heads Only**

If you only require attention heads for your own analysis, run the following script:

```bash
poetry run python scripts/run_attention_heads.py [OPTIONS]
```

Additionally, attention heads can be retrieved from CAAT in the same way as ColabFold with an additional custom flag for outputting attention heads locally:

```bash
poetry run colabfold_batch [OPTIONS] --output-attention-dir 'PATH/TO/HEADS'
```

Both prediction scripts accept additional ColabFold settings, including
`--num-recycles` (aliases: `--num-recycle`, `--num_recycles`),
`--recycle-early-stop-tolerance`, `--num-ensemble`, `--model-order`,
`--max-seq`, `--max-extra-seq`, `--msa-mode`, `--pair-mode`, `--pair-strategy`,
`--templates`, `--custom-template-path`, `--use-dropout`,
`--disable-cluster-profile`, `--rank`, `--stop-at-score`, `--save-recycles`,
and `--save-all`. Omitted settings retain ColabFold's defaults. In the
end-to-end pipeline, these settings apply to both query and target predictions.

```bash
poetry run python scripts/run_e2e_pipeline.py \
  --query-seq-path query.a3m --query-name query \
  --num-recycles 10 --recycle-early-stop-tolerance 0.0
```

**Option C: Run Analysis Only**

If you already have attention head .npy files and simply need to generate the plots and difference maps:

```bash
poetry run python scripts/run_analysis_pipeline.py [OPTIONS]
```

---

## 3. Output Overview

### 1. Prediction Results (`results/`)
This directory contains the standard structural outputs from the AlphaFold engine.

### 2. Raw Attention Data (`attention_outputs/`)
This folder contains the raw numerical matrices extracted during the model's forward pass.

### 3. Visualizations (`attention_visualizations/`)
This is the primary directory for human-readable insights.

#### **A. Average Attention Plots** (`*_average_attention.png`)
* The mean attention score for each amino acid residue across all layers and heads.

#### **B. Difference Maps** (`*_attention_difference.png`)
* The residue-by-residue delta between the Query and the Target.

## Visualization Dashboard

To view outputs of CAAT in depth, run the streamlit dashboard locally:

```bash
poetry run streamlit run dashboard/app.py
```


## References

- Mirdita M, Schütze K, Moriwaki Y, Heo L, Ovchinnikov S and Steinegger M. ColabFold: Making protein folding accessible to all. <br />
  Nature Methods (2022) doi: [10.1038/s41592-022-01488-1](https://www.nature.com/articles/s41592-022-01488-1)
- If you’re using **AlphaFold**, please also cite: <br />
  Jumper et al. "Highly accurate protein structure prediction with AlphaFold." <br />
  Nature (2021) doi: [10.1038/s41586-021-03819-2](https://doi.org/10.1038/s41586-021-03819-2)
- If you’re using **AlphaFold-multimer**, please also cite: <br />
  Evans et al. "Protein complex prediction with AlphaFold-Multimer." <br />
  biorxiv (2021) doi: [10.1101/2021.10.04.463034v1](https://www.biorxiv.org/content/10.1101/2021.10.04.463034v1)
- If you are using **RoseTTAFold**, please also cite: <br />
  Minkyung et al. "Accurate prediction of protein structures and interactions using a three-track neural network." <br />
  Science (2021) doi: [10.1126/science.abj8754](https://doi.org/10.1126/science.abj8754)
