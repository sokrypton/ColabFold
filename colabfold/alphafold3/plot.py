"""
MSA coverage plot for a fold input, from the a3m blocks on its protein chains.
"""
from typing import List

import numpy as np

from colabfold.plot import plot_msa

# the alphabet plot_msa expects, where 21 is a gap
_RESTYPES = "ARNDCQEGHILKMFPSTWYVX-"
_INDEX = {c: min(i, 21) for i, c in enumerate(_RESTYPES)}


def _rows(a3m: str) -> List[str]:
    rows, current = [], None
    for line in (a3m or "").splitlines():
        if line.startswith(">"):
            if current is not None:
                rows.append(current)
            current = ""
        elif current is not None:
            current += line.strip()
    if current:
        rows.append(current)
    return rows


_LOOKUP = np.full(256, 20, dtype=np.int8)
for _c, _i in _INDEX.items():
    _LOOKUP[ord(_c)] = _i


def _encode(rows: List[str], length: int) -> np.ndarray:
    out = np.full((len(rows), length), 21, dtype=np.int8)
    try:
        codes = np.frombuffer("".join(rows).encode("ascii"), dtype=np.uint8)
    except UnicodeEncodeError:
        return _encode_loop(rows, length, out)
    # a3m lower case is an insertion relative to the query
    keep = ~((codes >= ord("a")) & (codes <= ord("z")))
    row = np.repeat(np.arange(len(rows)), [len(r) for r in rows])[keep]
    kept = np.flatnonzero(keep)
    first = np.searchsorted(row, np.arange(len(rows)))
    col = np.arange(len(kept)) - first[row]
    inside = col < length
    out[row[inside], col[inside]] = _LOOKUP[codes[kept[inside]]]
    return out


def _encode_loop(rows: List[str], length: int, out: np.ndarray) -> np.ndarray:
    for i, row in enumerate(rows):
        ungapped = [c for c in row if not c.islower()]
        for j, c in enumerate(ungapped[:length]):
            out[i, j] = _INDEX.get(c.upper(), 20)
    return out


def plot_msa_coverage(fold_input, dpi: int = 200):
    chains = list(fold_input.protein_chains)
    if not chains:
        return None
    lengths = [len(c.sequence) for c in chains]
    total = sum(lengths)
    query = np.concatenate([_encode([c.sequence], len(c.sequence))[0] for c in chains])

    blocks = []
    offset = 0
    for chain, length in zip(chains, lengths):
        rows = _rows(chain.unpaired_msa) + _rows(chain.paired_msa)
        if not rows:
            offset += length
            continue
        block = np.full((len(rows), total), 21, dtype=np.int8)
        block[:, offset:offset + length] = _encode(rows, length)
        blocks.append(block)
        offset += length
    if not blocks:
        return None
    return plot_msa(np.concatenate(blocks), query, lengths, total, dpi=dpi)
