# CAAT

**CAAT** is an end-to-end pipeline designed to bridge the gap between protein structure prediction and functional biological insight. It automates the generation of AlphaFold2 structures and performs deep-dive analysis on the raw attention heads to identify residues of high structural importance to the model.

## Using CAAT via Colab Notebook

The quickest way to get up and running with CAAT is to use the publicly available notebook with the GPU runtime. You can find that [here](https://colab.research.google.com/drive/11UVndoYaP5cQD7762o8rT0DRlyzSnTnn?usp=sharing).

View coverage report [here](https://app.codecov.io/github/prameshsharma25/caat/tree/main).

## Quick Start with Poetry

Follow these steps to set up and run your local version of ColabFold using **Poetry**.

### Prerequisites

You need **Python 3.11** and **Poetry** installed on your system.

* **Python:** Install [Python 3.11](https://www.python.org/downloads/).
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

2.  **Install Dependencies:** Run the following command using Poetry.
    ```bash
    poetry install
    ```

    **Install AlphaFold Dependencies**:
    ```bash
    poetry install -E alphafold
    ```

    This command installs all dependencies defined in `pyproject.toml`, including your custom `alphafold-colabfold` package from PyPI.

3.  **Activate Environment (Optional):**
    ```bash
    poetry shell
    ```
    This prepares your terminal session for direct script execution.

---

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
