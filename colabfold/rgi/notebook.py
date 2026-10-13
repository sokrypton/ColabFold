"""Input preparation shared by the notebook forms and headless E2E runs."""

from __future__ import annotations

import copy
import hashlib
import json
import re

from colabfold.rgi import require_supported_model


def chain_table(fold_input):
    rows = []
    for entity in fold_input["sequences"]:
        kind, molecule = next(iter(entity.items()))
        ids = molecule["id"]
        for chain in [ids] if isinstance(ids, str) else ids:
            rows.append(
                {
                    "chain": chain,
                    "type": kind,
                    "residues": len(molecule["sequence"])
                    if "sequence" in molecule
                    else None,
                }
            )
    return rows


def prepare_input(
    fold_input, model, *, use_rgi=False, config=None, conformer_chains="ligands"
):
    """Return a fresh input; disabling RGI removes every RGI opt-in."""
    raw = copy.deepcopy(fold_input)
    raw.pop("restraints_config", None)
    entities = []
    for entity in raw["sequences"]:
        kind, molecule = next(iter(entity.items()))
        molecule.pop("conformer_restraints", None)
        ids = molecule["id"]
        for chain in [ids] if isinstance(ids, str) else ids:
            entities.append({kind: dict(molecule, id=chain)})
    raw["sequences"] = entities
    if model.startswith("esmfold2"):
        proteins = [row for row in chain_table(raw) if row["type"] == "protein"]
        if len(proteins) != 1:
            raise ValueError(
                "ESMFold2 language-model inputs need exactly one protein chain. "
                "Use one sequence or choose another predictor for a protein complex."
            )
    if not use_rgi:
        return raw
    require_supported_model(model)
    from rgi_toolkit.config import RestraintsConfig

    if not config:
        raise ValueError("RGI is enabled but no restraints were configured.")
    RestraintsConfig.from_dict(config)
    if config.get("conformer_restraints_config") is not None:
        table = chain_table(raw)
        if conformer_chains.strip() == "ligands":
            chosen = {row["chain"] for row in table if row["type"] == "ligand"}
        else:
            chosen = {
                item.strip() for item in conformer_chains.split(",") if item.strip()
            }
        known = {row["chain"] for row in table}
        if not chosen or not chosen <= known:
            raise ValueError(
                "Choose existing conformer chains from the table, or add a ligand."
            )
        for entity in entities:
            molecule = next(iter(entity.values()))
            if molecule["id"] in chosen:
                molecule["conformer_restraints"] = True
    raw["restraints_config"] = copy.deepcopy(config)
    return raw


def name_job(fold_input, model, *, base, num_recycles, num_diffusion_samples):
    """Keep changed model, sampling, and restraint settings out of stale results."""
    material = {key: value for key, value in fold_input.items() if key != "name"}
    material.update(
        model=model,
        num_recycles=num_recycles,
        num_diffusion_samples=num_diffusion_samples,
    )
    digest = hashlib.sha256(json.dumps(material, sort_keys=True).encode()).hexdigest()[
        :10
    ]
    base = re.sub(r"\W+", "", base).lower() or "job"
    if "restraints_config" in fold_input:
        return f"{base}_{model}_rgi_{digest}"
    return f"{base}_{model}_{digest}"
