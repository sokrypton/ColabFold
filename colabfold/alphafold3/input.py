"""
Build an alphafold3 ``folding_input.Input`` from ColabFold queries and MSAs.
"""
import logging
from pathlib import Path
from typing import Any, List, Optional, Sequence, Tuple

from colabfold.utils import MolType

logger = logging.getLogger(__name__)


def chain_id(index: int) -> str:
    """A, B, ... Z, AA, AB, ... in the order alphafold3 expects."""
    out = []
    while True:
        out.append(chr(index % 26 + ord("A")))
        index = index // 26 - 1
        if index < 0:
            break
    return "".join(reversed(out))


def paired_for_af3(paired_msa, num_chains: int, pairing: str = "colabfold"):
    """``colabfold`` keeps our row pairing; ``uniprot`` lets alphafold3 pair by species."""
    from colabfold.msa_pairing import PairingError, rewrite_paired_descriptions

    if not paired_msa or pairing == "uniprot":
        return paired_msa or None
    try:
        return rewrite_paired_descriptions(paired_msa, num_chains)
    except PairingError as e:
        logger.warning(f"keeping the original descriptions, alphafold3 will pair by species: {e}")
        return paired_msa


def _protein_chains(query_seqs_unique, query_seqs_cardinality, unpaired_msa, paired_msa,
                    start=0, templates=None):
    from alphafold3.common import folding_input

    chains, index = [], start
    for i, sequence in enumerate(query_seqs_unique):
        # "" is no MSA; None would make alphafold3 build one itself
        unpaired = unpaired_msa[i] if unpaired_msa else ""
        paired = paired_msa[i] if paired_msa else ""
        for _ in range(query_seqs_cardinality[i]):
            chains.append(folding_input.ProteinChain(
                id=chain_id(index),
                sequence=sequence,
                ptms=[],
                unpaired_msa=unpaired,
                paired_msa=paired,
                templates=list(templates[i]) if templates else [],
            ))
            index += 1
    return chains, index


def _molecule_chains(molecules, start=0):
    from alphafold3.common import folding_input

    chains, index = [], start
    for moltype, payload, copies in molecules or ():
        for _ in range(copies):
            cid = chain_id(index)
            if moltype is MolType.SMILES:
                chains.append(folding_input.Ligand(id=cid, smiles=payload))
            elif moltype is MolType.CCD:
                chains.append(folding_input.Ligand(id=cid, ccd_ids=[payload]))
            elif moltype is MolType.RNA:
                chains.append(folding_input.RnaChain(id=cid, sequence=payload,
                                                     modifications=[], unpaired_msa=""))
            elif moltype is MolType.DNA:
                chains.append(folding_input.DnaChain(id=cid, sequence=payload, modifications=[]))
            else:
                raise ValueError(f"cannot pass {moltype} to alphafold3")
            index += 1
    return chains, index


def build_fold_input(
    name: str,
    query_seqs_unique: List[str],
    query_seqs_cardinality: List[int],
    unpaired_msa: Optional[List[str]],
    paired_msa: Optional[List[str]],
    molecules: Optional[Sequence[Tuple[MolType, str, int]]] = None,
    seeds: Sequence[int] = (1,),
    pairing: str = "colabfold",
    templates=None,
):
    from alphafold3.common import folding_input

    paired_msa = paired_for_af3(paired_msa, len(query_seqs_unique), pairing)
    chains, index = _protein_chains(
        query_seqs_unique, query_seqs_cardinality, unpaired_msa, paired_msa, templates=templates
    )
    extra, _ = _molecule_chains(molecules, start=index)
    return no_af3_search(
        folding_input.Input(name=name, chains=chains + extra, rng_seeds=list(seeds))
    )


def load_fold_inputs(path: Path) -> List[Any]:
    """Read one or more fold inputs from an alphafold3 JSON file."""
    from alphafold3.common import folding_input

    return list(folding_input.load_fold_inputs_from_path(Path(path)))


def is_json_query(query) -> bool:
    extras = query[3] if len(query) > 3 else None
    return getattr(extras, "fold_input", None) is not None


def msa_state(fold_input) -> str:
    """``provided``, ``missing`` or ``mixed`` across the protein chains."""
    states = set()
    for chain in fold_input.protein_chains:
        has = chain.unpaired_msa is not None or chain.paired_msa is not None
        states.add("provided" if has else "missing")
    if not states:
        return "provided"
    return states.pop() if len(states) == 1 else "mixed"


def msas_of(fold_input) -> Tuple[List[str], Optional[List[str]]]:
    """The file's own MSAs, one per unique protein sequence."""
    unique, _ = sequences_of(fold_input)
    unpaired, paired = [""] * len(unique), [""] * len(unique)
    for chain in fold_input.protein_chains:
        i = unique.index(chain.sequence)
        unpaired[i] = chain.unpaired_msa or unpaired[i]
        paired[i] = chain.paired_msa or paired[i]
    return unpaired, (paired if any(paired) else None)


def with_msas(fold_input, unpaired_msa, paired_msa, pairing: str = "colabfold", templates=None):
    """Return a copy of ``fold_input`` with ColabFold's MSAs and templates on its protein chains."""
    import dataclasses

    from alphafold3.common import folding_input

    unique, _ = sequences_of(fold_input)
    paired_msa = paired_for_af3(paired_msa, len(unique), pairing)

    chains = []
    for chain in fold_input.chains:
        if isinstance(chain, folding_input.ProteinChain):
            # one MSA and one template list per unique sequence, so copies of a chain share them
            i = unique.index(chain.sequence)
            fetched = list(templates[i]) if templates else []
            # an MSA already in the file is the user's, and "" is a deliberate empty one:
            # only a field left unset asks ColabFold to fill it in
            unpaired, paired = chain.unpaired_msa, chain.paired_msa
            if unpaired is None:
                unpaired = unpaired_msa[i] if unpaired_msa else ""
            if paired is None:
                paired = paired_msa[i] if paired_msa else ""
            # ProteinChain is not a dataclass, so it is rebuilt rather than replaced
            chain = folding_input.ProteinChain(
                id=chain.id,
                sequence=chain.sequence,
                ptms=chain.ptms,
                description=chain.description,
                unpaired_msa=unpaired,
                paired_msa=paired,
                templates=chain.templates or fetched,
            )
        chains.append(chain)
    return no_af3_search(dataclasses.replace(fold_input, chains=chains))


def no_af3_search(fold_input):
    """Fill every unset MSA and template field: alphafold3 reads ``None`` as "search"."""
    return fold_input.fill_missing_fields()


def polymer_lengths(fold_input) -> List[int]:
    """Residues per protein, RNA and DNA chain, in chain order, for the plots."""
    from alphafold3.common import folding_input

    polymers = (folding_input.ProteinChain, folding_input.RnaChain, folding_input.DnaChain)
    return [len(c.sequence) for c in fold_input.chains if isinstance(c, polymers)]


def sequences_of(fold_input) -> Tuple[List[str], List[int]]:
    """Unique protein sequences and their copy counts, for the MSA search."""
    unique, counts = [], []
    for chain in fold_input.protein_chains:
        if chain.sequence in unique:
            counts[unique.index(chain.sequence)] += 1
        else:
            unique.append(chain.sequence)
            counts.append(1)
    return unique, counts
