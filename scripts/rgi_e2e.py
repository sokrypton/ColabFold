"""Run paired real predictions and validate written structures and RGI evidence.

Run with an environment containing the notebook dependencies:
  uv run --no-project python scripts/rgi_e2e.py --models boltz2
Use --work-dir to keep downloaded parameters and outputs between invocations.
"""

from __future__ import annotations

import argparse
import importlib.metadata
import json
import os
import subprocess
import sys
import urllib.request
from pathlib import Path


def prepare_runtime(work):
    from alphafold3.constants import ccd_fetch

    runner = work / "run_alphafold.py"
    if not runner.exists():
        urllib.request.urlretrieve(
            "https://raw.githubusercontent.com/sokrypton/alphafold3/v3.1.11/run_alphafold.py",
            runner,
        )
    root = Path(importlib.metadata.distribution("alphafold3-colabfold").locate_file(""))
    conv = root / "alphafold3/constants/converters"
    conv.mkdir(parents=True, exist_ok=True)
    if not (conv / "ccd.pickle").exists():
        ccd_fetch.write_pickles(
            ccd_fetch.codes_for_input(extra=["ATP"]),
            str(conv / "ccd.pickle"),
            str(conv / "chemical_component_sets.pickle"),
            libcifpp_dir=str(root / "share/libcifpp"),
        )
    return runner


def measure(path):
    import gemmi
    import numpy as np

    structure = gemmi.read_structure(str(path))
    groups = {}
    for chain in structure[0]:
        xyz = np.asarray(
            [
                [atom.pos.x, atom.pos.y, atom.pos.z]
                for residue in chain
                for atom in residue
            ]
        )
        if not len(xyz) or not np.isfinite(xyz).all():
            raise AssertionError(f"Empty or non-finite chain in {path.name}")
        groups[chain.name] = xyz.mean(axis=0)
    if "B" not in groups:
        residues = list(structure[0]["A"])
        groups = {
            key: np.asarray(
                [[a.pos.x, a.pos.y, a.pos.z] for r in subset for a in r]
            ).mean(axis=0)
            for key, subset in (("A", residues[:9]), ("B", residues[9:]))
        }
    return float(np.linalg.norm(groups["A"] - groups["B"]))


def run(args):
    import numpy as np
    from rgi_toolkit.notebook import make_config

    work = Path(args.work_dir).resolve()
    work.mkdir(parents=True, exist_ok=True)
    runner = prepare_runtime(work)
    os.environ.setdefault("XLA_FLAGS", "--xla_gpu_enable_triton_gemm=false")
    os.environ.setdefault("XLA_PYTHON_CLIENT_PREALLOCATE", "false")
    os.environ.setdefault("AF3_WEIGHTS_DIR", str(work / "weights"))
    summary = []
    for model in args.models:
        distances = {}
        for arm in args.arms:
            name = f"{model}_{arm}"
            sequence = "ACDEFGHIK"
            raw = {
                "name": name,
                "dialect": "alphafold3",
                "version": 1,
                "modelSeeds": [42],
                "sequences": [
                    {
                        "protein": {
                            "id": chain,
                            "sequence": sequence,
                            "templates": [],
                            "unpairedMsa": f">query\n{sequence}\n",
                            "pairedMsa": "",
                        }
                    }
                    for chain in ("A", "B")
                ],
            }
            if model.startswith("esmfold2"):
                protein = raw["sequences"][0]["protein"]
                protein["sequence"] = sequence * 2
                protein["unpairedMsa"] = f">query\n{sequence * 2}\n"
                raw["sequences"] = raw["sequences"][:1]
            if args.fixture == "combined":
                raw["sequences"].append(
                    {"ligand": {"id": "C", "smiles": "O=C(O)/C=C/C(=O)O"}}
                )
            if arm != "vanilla":
                target = 25.0 if arm == "rgi" else 35.0
                config = {
                    "verbose": True,
                    "distance_restraints_config": [
                        {
                            "atom_selection1": "chain A and resid 1 to 9"
                            if model.startswith("esmfold2")
                            else "chain A",
                            "atom_selection2": "chain A and resid 10 to 18"
                            if model.startswith("esmfold2")
                            else "chain B",
                            "harmonic": {"target_distance": target},
                        }
                    ],
                }
                if args.fixture == "combined":
                    config["conformer_restraints_config"] = {"plane": {"weight": 1.0}}
                    raw["sequences"][-1]["ligand"]["conformer_restraints"] = True
                raw["restraints_config"] = make_config(config)
            path = work / f"{name}.json"
            path.write_text(json.dumps(raw, indent=2) + "\n")
            out = work / "outputs"
            command = [
                sys.executable,
                *([] if arm == "vanilla" else ["-m", "colabfold.rgi"]),
                str(runner),
                f"--json_path={path}",
                f"--model={model}",
                "--norun_data_pipeline",
                f"--output_dir={out}",
                "--force_output_dir",
                "--flash_attention_implementation=xla",
                "--num_recycles=1",
                "--num_diffusion_samples=1",
                "--buckets=32",
                "--weights_precision=fp32"
                if model == "alphafold3"
                else "--weights_precision=int8",
                f"--cache_dir={work / 'compile_cache'}",
            ]
            if model == "chai1" or model.startswith("esmfold2"):
                command.append("--use_esm_embeddings")
            if model == "alphafold3":
                weights = work / "af3_weights"
                weights.mkdir(exist_ok=True)
                if not (weights / "af3.bin.zst").exists():
                    urllib.request.urlretrieve(
                        "https://storage.googleapis.com/alphafold3/af3.bin.zst",
                        weights / "af3.bin.zst",
                    )
                command.append(f"--model_dir={weights}")
            log_path = work / f"{name}.log"
            print(f"Running {model} / {arm}; log: {log_path.name}", flush=True)
            with log_path.open("w") as log:
                proc = subprocess.run(
                    command, stdout=log, stderr=subprocess.STDOUT, check=False
                )
            if proc.returncode:
                print(log_path.read_text()[-12000:], flush=True)
                raise RuntimeError(f"{model}/{arm} exited {proc.returncode}")
            job = out / name
            cifs = sorted(job.glob("seed-*_sample-*/*.cif"))
            if not cifs:
                raise AssertionError(f"No structures for {model}/{arm}")
            values = [measure(cif) for cif in cifs]
            distances[arm] = values
            report = job / "rgi_report.json"
            if arm != "vanilla":
                evidence = json.loads(report.read_text())
                assert evidence["seeds"][0]["inventory"]["distance"] == 1
                if args.fixture == "combined":
                    inventory = evidence["seeds"][0]["inventory"]
                    assert all(
                        inventory[key] > 0
                        for key in ("bond", "angle", "plane", "cistrans")
                    ), inventory
                assert all(abs(value - target) < 0.5 for value in values), values
            elif report.exists():
                raise AssertionError("Vanilla must not produce an RGI report")
            summary.append(
                {
                    "model": model,
                    "arm": arm,
                    "samples": len(cifs),
                    "distance_angstrom": values,
                }
            )
            (work / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
            print(json.dumps(summary[-1]), flush=True)
        if "rgi" in distances and "vanilla" in distances:
            assert not np.allclose(distances["rgi"], distances["vanilla"], atol=1.0)
    print("E2E PASSED", json.dumps(summary), flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--models", nargs="+", default=["boltz2"])
    parser.add_argument(
        "--arms",
        nargs="+",
        choices=["vanilla", "rgi", "rgi_other"],
        default=["vanilla", "rgi"],
    )
    parser.add_argument("--work-dir", default=".cache/rgi-e2e/run")
    parser.add_argument(
        "--fixture", choices=["distance", "combined"], default="distance"
    )
    run(parser.parse_args())
