"""Scoped glue between the ColabFold JAX runner and the shared RGI engine."""

from __future__ import annotations

import copy
import functools
import inspect
import json
from contextlib import contextmanager
from pathlib import Path

from colabfold.citations import write_rgi_citation
from colabfold.rgi import ALPHAFOLD3_COLABFOLD_VERSION, require_supported_model


@contextmanager
def _replace(obj, name, value):
    original = getattr(obj, name)
    setattr(obj, name, value)
    try:
        yield
    finally:
        setattr(obj, name, original)


def split_input(raw, base_dir=None):
    """Strip only RGI extensions before calling the unchanged upstream parser."""
    from rgi_toolkit.config import RestraintsConfig, resolve_restraints_config

    clean = copy.deepcopy(raw)
    config = clean.pop("restraints_config", None)
    opted = {}
    for entity in clean.get("sequences", []):
        for molecule in entity.values():
            flag = molecule.pop("conformer_restraints", False)
            if not isinstance(flag, bool):
                raise TypeError("conformer_restraints must be true or false.")
            ids = molecule.get("id", [])
            for chain in [ids] if isinstance(ids, str) else ids:
                opted[chain] = flag
    if config is not None:
        config = resolve_restraints_config(config, base_dir=base_dir)
        RestraintsConfig.from_dict(config)
    elif any(opted.values()):
        raise ValueError("Conformer opt-in requires restraints_config.")
    return clean, config, opted


def build_adapter(fold_input, example, opted):
    """Resolve predictor chemistry; reuse the existing framework-free AF3 adapter."""
    from alphafold3.common import folding_input
    from alphafold3.constants import decoded_ccd, residue_names
    from alphafold3.data.tools import rdkit_utils
    from rdkit import Chem
    from rgi_toolkit.alphafold3.adapter import AF3RestraintAdapter

    mols = []
    ccd = None
    for chain in fold_input.chains:
        if not opted.get(chain.id) or not isinstance(chain, folding_input.Ligand):
            continue
        if chain.smiles is not None:
            mol = Chem.MolFromSmiles(chain.smiles)
            if mol is None:
                raise ValueError(f"Invalid ligand SMILES for chain {chain.id}.")
            mols.append((chain.id, mol, True))
        else:
            if len(chain.ccd_ids or []) != 1:
                raise ValueError(
                    "Ligand geometry currently needs one CCD per ligand chain."
                )
            if ccd is None:
                ccd = decoded_ccd.get_ccd(user_ccd=fold_input.user_ccd)
            mol = rdkit_utils.mol_from_ccd_cif(
                ccd[chain.ccd_ids[0]], sort_alphabetically=False, remove_hydrogens=True
            )
            mols.append((chain.id, mol, False))
    mapping = {chain.id: i + 1 for i, chain in enumerate(fold_input.chains)}
    return AF3RestraintAdapter(
        example,
        mapping,
        residue_names.POLYMER_TYPES,
        mols,
        {mapping[chain]: flag for chain, flag in opted.items() if chain in mapping},
    )


def _guided_apply(runner, model_runner):
    import haiku as hk
    import jax
    from alphafold3.model.network import diffusion_head

    from colabfold.rgi.sampler import coordinate_minimizer, make_guided_sampler

    original = diffusion_head.sample

    @hk.transform
    def forward(batch, minimize, gather):
        sampler = make_guided_sampler(original, coordinate_minimizer(minimize, gather))
        with _replace(diffusion_head, "sample", sampler):
            return runner.model.Model(model_runner._model_config)(
                batch, use_dropout=model_runner._use_dropout
            )

    model_runner._preinit_tokamax_context()
    return forward.apply if runner._NOJIT.value else jax.jit(forward.apply)


def install(runner):
    """Install input and per-job hooks in one CLI process, preserving vanilla."""
    import importlib.metadata

    if (
        importlib.metadata.version("alphafold3-colabfold")
        != ALPHAFOLD3_COLABFOLD_VERSION
    ):
        raise RuntimeError(
            "RGI requires alphafold3-colabfold==3.1.11; re-run installation."
        )

    input_class = runner.folding_input.Input
    original_parse = input_class.from_json
    original_predict = runner.predict_structure
    original_write_input = runner.write_fold_input_json
    original_write_outputs = runner.write_outputs
    jobs, reports = {}, {}

    @classmethod
    def parse(cls, json_str, json_path=None):
        raw = json.loads(json_str)
        clean, config, opted = split_input(
            raw, Path(json_path).parent if json_path is not None else None
        )
        result = original_parse(json.dumps(clean), json_path)
        settings = (config, opted)
        if result.name in jobs and jobs[result.name] != settings:
            raise ValueError("Give jobs with different RGI settings different names.")
        jobs[result.name] = settings
        return result

    def write_input(fold_input, output_dir):
        original_write_input(fold_input, output_dir)
        # A failed rerun must never leave the previous run's success report.
        (Path(output_dir) / "rgi_report.json").unlink(missing_ok=True)
        config, opted = jobs.get(fold_input.name, (None, {}))
        if config is not None:
            path = Path(output_dir) / f"{fold_input.sanitised_name()}_data.json"
            raw = json.loads(path.read_text())
            raw["restraints_config"] = config
            entities = []
            for entity in raw["sequences"]:
                kind = next(iter(entity))
                molecule = next(iter(entity.values()))
                ids = molecule["id"]
                ids = [ids] if isinstance(ids, str) else ids
                for chain in ids:
                    item = dict(molecule, id=chain)
                    if opted.get(chain):
                        item["conformer_restraints"] = True
                    entities.append({kind: item})
            raw["sequences"] = entities
            path.write_text(json.dumps(raw, indent=2) + "\n")

    def predict(fold_input, model_runner, **kwargs):
        import jax.numpy as jnp
        import numpy as np
        from rgi_toolkit.combined import CombinedRestraints
        from rgi_toolkit.notebook import distance_report, restraint_inventory

        config, opted = jobs.get(fold_input.name, (None, {}))
        if config is None:
            return original_predict(fold_input, model_runner, **kwargs)
        require_supported_model(model_runner.model_name)
        records = []
        original_inference = model_runner.run_inference
        if not hasattr(model_runner, "_rgi_apply"):
            model_runner._rgi_apply = _guided_apply(runner, model_runner)

        def inference(example, rng_key):
            restraints = CombinedRestraints()
            adapter = build_adapter(fold_input, example, opted)
            restraints.setup(adapter, config=config)
            if not restraints.is_active():
                raise ValueError(
                    "RGI built no active restraints. "
                    "Check selections and ligand opt-ins."
                )
            if config.get("conformer_restraints_config") is not None and not any(
                opted.values()
            ):
                raise ValueError(
                    "Ligand/polymer geometry requires "
                    "conformer_restraints: true on a chain."
                )
            inventory = restraint_inventory(restraints)
            print(
                "RGI built restraints:",
                json.dumps(inventory, sort_keys=True),
                flush=True,
            )
            minimize = restraints.get_minimizer()
            gather = example.get("structbook/residue_atom_gather")
            if gather is not None:
                gather = np.asarray(gather)
                mask = np.asarray(example["ref_mask"], bool)
                if gather.shape != mask.shape or np.any(mask & (gather < 0)):
                    raise ValueError(
                        "OpenDDE structural atom map does not cover the RGI atoms."
                    )
                valid = gather[gather >= 0]
                if len(np.unique(valid)) != len(valid):
                    raise ValueError(
                        "OpenDDE structural atom map contains duplicate atoms."
                    )
                gather = jnp.asarray(gather)
            apply = model_runner._rgi_apply
            had_model = "_model" in vars(model_runner)
            cached_model = vars(model_runner).get("_model")
            model_runner._model = lambda key, batch: apply(
                model_runner.model_params, key, batch, minimize, gather
            )
            try:
                result = original_inference(example, rng_key)
            finally:
                if had_model:
                    model_runner._model = cached_model
                else:
                    del model_runner._model
            coords = np.asarray(result["diffusion_samples"]["atom_positions"])
            flat = coords.reshape(coords.shape[0], -1, 3)
            if not np.isfinite(flat).all():
                raise RuntimeError(
                    "The predictor produced non-finite atom coordinates."
                )
            restraints.finalize(jnp.asarray(flat), 0)
            records.append(
                {
                    "inventory": inventory,
                    "finite": True,
                    "distances": distance_report(restraints, flat),
                }
            )
            print(
                "RGI measured distances:",
                json.dumps(records[-1]["distances"]),
                flush=True,
            )
            return result

        with _replace(model_runner, "run_inference", inference):
            results = original_predict(fold_input, model_runner, **kwargs)
        for record, result in zip(records, results, strict=True):
            record["seed"] = result.seed
        reports[fold_input.sanitised_name()] = {
            "model": model_runner.model_name,
            "restraints_config": config,
            "seeds": records,
        }
        return results

    @functools.wraps(original_write_outputs)
    def write_outputs(*args, **kwargs):
        original_write_outputs(*args, **kwargs)
        bound = inspect.signature(original_write_outputs).bind(*args, **kwargs)
        report = reports.pop(bound.arguments["job_name"], None)
        if report is not None:
            path = Path(bound.arguments["output_dir"]) / "rgi_report.json"
            path.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
            write_rgi_citation(path.parent)

    input_class.from_json = parse
    runner.predict_structure = predict
    runner.write_fold_input_json = write_input
    runner.write_outputs = write_outputs
