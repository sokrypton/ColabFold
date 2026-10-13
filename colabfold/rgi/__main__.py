"""Run the pinned upstream CLI with optional input-driven RGI hooks."""

import importlib.util
import sys
from pathlib import Path


def main():
    if len(sys.argv) < 2 or sys.argv[1] in {"-h", "--help"}:
        print("Usage: python -m colabfold.rgi run_alphafold.py [upstream flags]")
        return
    script = Path(sys.argv[1]).resolve()
    if not script.is_file():
        raise FileNotFoundError(script)
    sys.argv = [str(script), *sys.argv[2:]]
    spec = importlib.util.spec_from_file_location("_rgi_colabfold_runner", script)
    runner = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = runner
    spec.loader.exec_module(runner)
    from colabfold.rgi.runtime import install

    install(runner)
    runner.app.run(runner.main)


if __name__ == "__main__":
    main()
