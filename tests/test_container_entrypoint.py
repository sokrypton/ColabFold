"""Check command dispatch and relative paths without a container runtime."""

import json
import os
from pathlib import Path
import subprocess
import sys

import pytest


@pytest.mark.parametrize(
    "arguments, expected",
    [
        ([], ["/opt/CAAT/scripts/run_e2e_pipeline.py", "--help"]),
        (["--num-recycles", "3"], ["/opt/CAAT/scripts/run_e2e_pipeline.py", "--num-recycles", "3"]),
        (["scripts/run_e2e_pipeline.py", "--data-dir", "weights with spaces"],
         ["/opt/CAAT/scripts/run_e2e_pipeline.py", "--data-dir", "weights with spaces"]),
        (["run_attention_heads.py", "--help"], ["/opt/CAAT/scripts/run_attention_heads.py", "--help"]),
        (["scripts/run_analysis_pipeline.py", "--help"], ["/opt/CAAT/scripts/run_analysis_pipeline.py", "--help"]),
        (["python", "-m", "pip", "check"], ["-m", "pip", "check"]),
        (["custom script.py", "input file.a3m"], ["custom script.py", "input file.a3m"]),
        (["colabfold_batch", "--help"], ["--help"]),
    ],
)
def test_entrypoint_preserves_working_directory_and_arguments(tmp_path, arguments, expected):
    bin_dir = tmp_path / "bin"
    bin_dir.mkdir()
    stub = (
        f"#!{sys.executable}\n"
        "import json, os, sys\n"
        "print(json.dumps({'cwd': os.getcwd(), 'args': sys.argv[1:]}))\n"
    )
    for command in ("python", "colabfold_batch"):
        executable = bin_dir / command
        executable.write_text(stub)
        executable.chmod(0o755)
    entrypoint = Path(__file__).resolve().parents[1] / "images/entrypoint.sh"
    result = subprocess.check_output(
        [str(entrypoint), *arguments], cwd=tmp_path,
        env=dict(os.environ, PATH=f"{bin_dir}:{os.environ['PATH']}"), text=True,
    )
    assert json.loads(result) == {"cwd": str(tmp_path.resolve()), "args": expected}
