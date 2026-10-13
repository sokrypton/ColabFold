"""Preserve legacy notebook molecules and MSAs in Boltz's RGI-capable YAML."""

import copy


def notebook_input(
    entries, *, config=None, conformer_chains="ligands", single_sequence=False
):
    sequences = []
    ligands = set()
    for header, sequence in entries:
        fields = header.lstrip(">").split("|")
        chain, kind = fields[:2]
        molecule = {"id": chain}
        if kind == "protein":
            molecule.update(
                sequence=sequence, msa="empty" if single_sequence else fields[2]
            )
        elif kind.lower() in {"dna", "rna"}:
            kind = kind.lower()
            molecule["sequence"] = sequence
        elif kind in {"smiles", "ccd"}:
            molecule[kind] = sequence
            kind = "ligand"
            ligands.add(chain)
        else:
            raise ValueError(f"Unsupported molecule type {kind!r}.")
        sequences.append({kind: molecule})
    raw = {"version": 1, "sequences": sequences}
    if config is not None:
        from rgi_toolkit.config import RestraintsConfig

        RestraintsConfig.from_dict(config)
        raw["restraints_config"] = copy.deepcopy(config)
        if config.get("conformer_restraints_config") is not None:
            chosen = (
                ligands
                if conformer_chains == "ligands"
                else {
                    chain.strip()
                    for chain in conformer_chains.split(",")
                    if chain.strip()
                }
            )
            known = {next(iter(entity.values()))["id"] for entity in sequences}
            if not chosen or not chosen <= known:
                raise ValueError("Choose existing conformer chains or add a ligand.")
            for entity in sequences:
                molecule = next(iter(entity.values()))
                if molecule["id"] in chosen:
                    molecule["conformer_restraints"] = True
    return raw
