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


@pytest.mark.parametrize("cache_mode", ["home", "blocked_home", "no_home", "explicit", "temporary"])
def test_entrypoint_selects_writable_cache(tmp_path, cache_mode):
    home = tmp_path / "home"
    home.mkdir()
    work = tmp_path / "work"
    work.mkdir()
    env = dict(os.environ, HOME=str(home), TMPDIR=str(tmp_path))
    env.pop("XDG_CACHE_HOME", None)
    if cache_mode in ("blocked_home", "temporary"):
        (home / ".cache").write_text("not a directory")
    if cache_mode == "no_home":
        env.pop("HOME")
    if cache_mode == "temporary":
        (work / ".cache").write_text("not a directory")
    if cache_mode == "explicit":
        env["XDG_CACHE_HOME"] = str(tmp_path / "chosen")
    entrypoint = Path(__file__).resolve().parents[1] / "images/entrypoint.sh"
    # Use the real interpreter to check appdirs' actual cache resolution.
    script = tmp_path / "check.py"
    script.write_text(
        "import appdirs, json, pathlib\n"
        "p = pathlib.Path(appdirs.user_cache_dir('colabfold'))\n"
        "(p / 'params').mkdir(parents=True)\n"
        "(p / 'params' / 'probe').write_text('ok')\n"
        "print(json.dumps(str(p)))\n"
    )
    env["PATH"] = str(Path(sys.executable).parent) + os.pathsep + env["PATH"]
    result = subprocess.check_output([str(entrypoint), str(script)], cwd=work, env=env, text=True)
    cache = Path(json.loads(result))
    if cache_mode == "temporary":
        assert cache.parent.parent == tmp_path
        assert cache.parent.name.startswith("caat-cache.")
    else:
        base = {"home": home / ".cache", "explicit": tmp_path / "chosen"}.get(cache_mode, work / ".cache")
        assert cache == base / "colabfold"
    assert (cache / "params" / "probe").read_text() == "ok"
