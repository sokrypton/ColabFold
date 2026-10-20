"""Offline container check: installed dependencies, CLIs, and real analysis."""

import os
from pathlib import Path
import subprocess
import sys
import tempfile

import jax
import numpy as np
import pandas as pd
from PIL import Image
from packaging.version import Version

ROOT = Path(__file__).resolve().parents[1]


def main():
    subprocess.run([sys.executable, "-m", "pip", "check"], check=True)
    print(f"Python: {sys.version}; JAX: {jax.__version__}; devices: {jax.devices()}")
    assert sys.version_info >= (3, 11) and sys.version_info.releaselevel == "final"
    assert Version("0.5.2") <= Version(jax.__version__) < Version("0.11")
    np.testing.assert_array_equal(np.asarray(jax.numpy.ones(2) + 1), [2, 2])

    with tempfile.TemporaryDirectory(prefix="caat-smoke-") as directory:
        work = Path(directory)
        env = dict(os.environ, MPLCONFIGDIR=str(work / "mpl"))
        for script in ("run_e2e_pipeline.py", "run_attention_heads.py", "run_analysis_pipeline.py"):
            help_text = subprocess.check_output(
                [str(ROOT / "images/entrypoint.sh"), f"scripts/{script}", "--help"],
                cwd=work, env=env, text=True,
            )
            if script != "run_analysis_pipeline.py":
                assert "--num-recycles" in help_text and "--data-dir" in help_text

        query = work / "query.a3m"
        query.write_text(">smoke\nACDEFGHIKLMNPQRSTVWY\n")
        attention = work / "attention"
        attention.mkdir()
        rng = np.random.default_rng(7)
        for head in range(2):
            np.save(
                attention / f"seed_000_model_1_recycle_0_main_evoformer_loop_1_head_{head}.npy",
                rng.normal(size=(20, 4, 20, 20)).astype(np.float16),
            )
        # Relative inputs and outputs deliberately test preservation of cwd.
        subprocess.run(
            [str(ROOT / "images/entrypoint.sh"), "scripts/run_analysis_pipeline.py",
             "--query-seq-path", "query.a3m", "--query-attn-dir", "attention",
             "--query-name", "smoke", "--output-dir", "analysis", "--save-attention-npy"],
            cwd=work, env=env, check=True,
        )
        output = work / "analysis/smoke"
        ranking = pd.read_csv(output / "smoke_residue_ranking.csv")
        scores = ranking["Attention score"].to_numpy()
        assert len(scores) == 20 and np.isfinite(scores).all()
        # Seeded, nonuniform inputs must retain meaningful residue variation.
        assert np.ptp(scores) > 1e-6, "Attention scores lost their variation"
        with Image.open(output / "smoke_average_attention.png") as plot:
            plot.verify()
    print("CAAT container smoke test passed.")


if __name__ == "__main__":
    main()
