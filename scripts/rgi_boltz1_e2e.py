"""Exercise native Boltz-1 using the exact notebook input builder."""

import argparse
import json
import os
import subprocess
import sys
from pathlib import Path

import yaml
from rgi_e2e import measure
from rgi_toolkit.notebook import make_config

from colabfold.rgi.boltz import notebook_input


def run(args):
    work = Path(args.work_dir).resolve()
    work.mkdir(parents=True, exist_ok=True)
    summary = []
    entries = [(f">{chain}|protein|empty", "ACDEFGHIK") for chain in ("A", "B")]
    for arm in ("vanilla", "rgi"):
        config = (
            None
            if arm == "vanilla"
            else make_config(
                {
                    "verbose": True,
                    "distance_restraints_config": [
                        {
                            "atom_selection1": "chain A",
                            "atom_selection2": "chain B",
                            "harmonic": {"target_distance": 25},
                        }
                    ],
                }
            )
        )
        data = notebook_input(entries, config=config, single_sequence=True)
        path = work / f"boltz1_{arm}.yaml"
        path.write_text(yaml.safe_dump(data, sort_keys=False))
        cmd = [
            args.python,
            "-m",
            "boltz.main",
            "predict",
            str(path),
            "--model",
            "boltz1",
            "--out_dir",
            str(work),
            "--recycling_steps",
            "1",
            "--diffusion_samples",
            "1",
            "--seed",
            "42",
            "--override",
            "--num_workers",
            "0",
        ]
        log_path = work / f"boltz1_{arm}.log"
        with log_path.open("w") as log:
            proc = subprocess.run(
                cmd, stdout=log, stderr=subprocess.STDOUT, check=False
            )
        if proc.returncode:
            raise RuntimeError(log_path.read_text()[-12000:])
        cifs = list(
            (work / f"boltz_results_{path.stem}" / "predictions" / path.stem).glob(
                "*.cif"
            )
        )
        assert len(cifs) == 1, cifs
        distance = measure(cifs[0])
        if arm == "rgi":
            assert "distances=1" in log_path.read_text()
            assert abs(distance - 25) < 0.5, distance
        summary.append({"model": "boltz1", "arm": arm, "distance_angstrom": distance})
        print(json.dumps(summary[-1]), flush=True)
    assert abs(summary[0]["distance_angstrom"] - summary[1]["distance_angstrom"]) > 1
    (work / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--python", default=os.environ.get("RGI_BOLTZ_PYTHON", sys.executable)
    )
    parser.add_argument("--work-dir", default=".cache/rgi-e2e/boltz1")
    run(parser.parse_args())
