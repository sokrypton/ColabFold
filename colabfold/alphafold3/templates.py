"""
Turn ColabFold's template hits into alphafold3 ``folding_input.Template`` objects.

alphafold3 consumes templates but never searches for them: it wants one mmCIF per hit,
filtered to a single polymer chain, plus a query-residue -> template-residue map. The
map indexes the chain's full SEQRES *including residues missing from the structure*,
which is not what hhsearch reports, so the hit sequence is aligned to the SEQRES here.
"""
import logging
from pathlib import Path
from typing import Dict, List, Tuple

logger = logging.getLogger(__name__)

# Not wired into the backend yet: the index composition is tested, a real hhsearch hit is not.


class TemplateError(ValueError):
    pass


def single_letter_sequence(mmcif_text: str, chain_id: str) -> str:
    """The SEQRES of one label chain, from ``_entity_poly_seq``/``_pdbx_poly_seq_scheme``."""
    from alphafold3 import structure

    struc = structure.from_mmcif(mmcif_text, include_water=False, include_bonds=False)
    sequences = struc.chain_single_letter_sequence()
    if chain_id in sequences:
        return sequences[chain_id]
    raise TemplateError(f"chain {chain_id!r} is not in the template ({sorted(sequences)})")


def map_hit_to_seqres(hit_sequence: str, seqres: str) -> Dict[int, int]:
    """Ungapped-hit index -> SEQRES index. The hit is a contiguous run of the SEQRES."""
    ungapped = hit_sequence.replace("-", "")
    if not ungapped:
        raise TemplateError("the hit aligns nothing")
    start = seqres.find(ungapped)
    if start < 0:
        raise TemplateError("the hit sequence is not a substring of the template SEQRES")
    return {i: start + i for i in range(len(ungapped))}


def query_to_hit_map(hit, query_sequence: str) -> Dict[int, int]:
    """Query index -> ungapped-hit index, the way AlphaFold2 reads an hhr hit."""
    hit_indices = [i for i in hit.indices_hit if i > -1]
    query_indices = [i for i in hit.indices_query if i > -1]
    if not hit_indices or not query_indices:
        raise TemplateError("the hit has no aligned residue")
    hit_min, query_min = min(hit_indices), min(query_indices)
    offset = query_sequence.find(hit.query.replace("-", ""))
    if offset < 0:
        raise TemplateError("the hit's query sequence is not part of the query")

    ungapped_hit = len(hit.hit_sequence.replace("-", ""))
    mapping = {}
    for q, h in zip(hit.indices_query, hit.indices_hit):
        if q == -1 or h == -1:
            continue
        q, h = q - query_min + offset, h - hit_min
        if h < ungapped_hit and q < len(query_sequence):
            mapping[q] = h
    return mapping


def build_templates(hits, query_sequence: str, cif_dir: Path,
                    max_templates: int = 20) -> List:
    """``folding_input.Template`` per hhsearch hit that maps cleanly onto its mmCIF."""
    from alphafold3.common import folding_input

    from colabfold.alphafold3 import require

    require()
    templates = []
    for hit in hits:
        if len(templates) >= max_templates:
            break
        pdb_id, _, chain_id = hit.name.partition("_")
        chain_id = (chain_id or "A").split()[0]
        cif = Path(cif_dir).joinpath(f"{pdb_id.lower()}.cif")
        if not cif.is_file():
            logger.warning(f"{hit.name}: no {cif.name} beside the hits, skipping")
            continue
        try:
            text, label_id = filter_to_chain(cif.read_text(), chain_id)
            seqres = single_letter_sequence(text, label_id)
            hit_to_seqres = map_hit_to_seqres(hit.hit_sequence, seqres)
            mapping = {q: hit_to_seqres[h]
                       for q, h in query_to_hit_map(hit, query_sequence).items()
                       if h in hit_to_seqres}
            if not mapping:
                raise TemplateError("nothing mapped onto the SEQRES")
            templates.append(folding_input.Template(mmcif=text, query_to_template_map=mapping))
        except (TemplateError, KeyError, ValueError) as e:
            logger.warning(f"{hit.name}: dropping this template, {e}")
    return templates


def filter_to_chain(mmcif_text: str, chain_id: str) -> Tuple[str, str]:
    """Keep one polymer chain, and say which label id it is keyed by.

    pdb70 names a hit by its auth chain, the sequences below are keyed by label.
    """
    from alphafold3 import structure

    struc = structure.from_mmcif(mmcif_text, include_water=False, include_bonds=False)
    chains = struc.polymer_auth_asym_id_to_label_asym_id()
    label = chains.get(chain_id, chain_id)
    filtered = struc.filter(chain_id=[c for c in (chain_id, label) if c])
    return filtered.to_mmcif(), label
