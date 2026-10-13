"""Verify vanilla/guided transitions and config isolation in one real JAX process."""

import argparse
import json
import os
import subprocess
import sys
from pathlib import Path

from rgi_e2e import measure, prepare_runtime
from rgi_toolkit.notebook import make_config


def run(work):
    work = Path(work).resolve()
    work.mkdir(parents=True, exist_ok=True)
    runner = prepare_runtime(work)
    inputs = work / "inputs"
    inputs.mkdir(exist_ok=True)
    targets = {
        "a_vanilla": None,
        "b_distance25": 25,
        "c_distance35": 35,
        "d_custom30": 30,
        "e_vanilla": None,
    }
    for name, target in targets.items():
        data = {
            "name": name,
            "dialect": "alphafold3",
            "version": 1,
            "modelSeeds": [42, 43],
            "sequences": [
                {
                    "protein": {
                        "id": chain,
                        "sequence": "ACDEFGHIK",
                        "templates": [],
                        "unpairedMsa": ">query\nACDEFGHIK\n",
                        "pairedMsa": "",
                    }
                }
                for chain in ("A", "B")
            ],
        }
        if name == "d_custom30":
            data["restraints_config"] = {
                "verbose": True,
                "custom_restraints_config": [
                    {
                        "name": "centroid_target",
                        "selections": {"A": "chain A", "B": "chain B"},
                        "energy": "(distance(A, B) - 30)**2",
                    }
                ],
            }
        elif target is not None:
            data["restraints_config"] = make_config(
                {
                    "distance_restraints_config": [
                        {
                            "atom_selection1": "chain A",
                            "atom_selection2": "chain B",
                            "harmonic": {"target_distance": target},
                        }
                    ]
                }
            )
        (inputs / f"{name}.json").write_text(json.dumps(data, indent=2))
    os.environ.setdefault("XLA_FLAGS", "--xla_gpu_enable_triton_gemm=false")
    os.environ.setdefault("XLA_PYTHON_CLIENT_PREALLOCATE", "false")
    command = [
        sys.executable,
        "-m",
        "colabfold.rgi",
        str(runner),
        f"--input_dir={inputs}",
        "--model=boltz2",
        "--norun_data_pipeline",
        f"--output_dir={work / 'outputs'}",
        "--force_output_dir",
        "--flash_attention_implementation=xla",
        "--num_recycles=1",
        "--num_diffusion_samples=2",
        "--buckets=32",
    ]
    with (work / "batch.log").open("w") as log:
        subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, check=True)
    summary = {}
    for name, target in targets.items():
        job = work / "outputs" / name
        cifs = sorted(job.glob("seed-*_sample-*/*.cif"))
        assert len(cifs) == 4, (name, cifs)
        values = [measure(path) for path in cifs]
        if target is not None:
            assert all(abs(value - target) < 0.5 for value in values), (name, values)
            report = json.loads((job / "rgi_report.json").read_text())
            assert [record["seed"] for record in report["seeds"]] == [42, 43]
            term = "custom" if name == "d_custom30" else "distance"
            assert all(record["inventory"][term] == 1 for record in report["seeds"])
        else:
            assert not (job / "rgi_report.json").exists()
        summary[name] = values
    assert summary["a_vanilla"] == summary["e_vanilla"], summary
    (work / "summary.json").write_text(json.dumps(summary, indent=2))
    print("BATCH E2E PASSED", json.dumps(summary), flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--work-dir", default=".cache/rgi-e2e/batch")
    run(parser.parse_args().work_dir)
