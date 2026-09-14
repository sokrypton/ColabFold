"""
Carry ColabFold's pairing through alphafold3's species-id pairing.

alphafold3 pairs rows by a UniProt mnemonic species id parsed from each description,
which only really works for bacteria, and drops rows that do not match the pattern.
ColabFold has already paired by taxid, so row *i* of every chain belongs together:
writing the row index in as the species id makes alphafold3 keep that pairing.
"""
import logging
from typing import List, Sequence

logger = logging.getLogger(__name__)

# alphafold3's species id is [A-Z0-9]{1,5}, so five digits is the ceiling
MAX_PAIRED_ROWS = 100_000

_SYNTHETIC = "tr|AAAAAA|AAAAAA_{row:05d}"


class PairingError(ValueError):
    pass


def synthetic_description(row: int) -> str:
    if not 0 < row < MAX_PAIRED_ROWS:
        raise PairingError(f"row {row} is outside 1..{MAX_PAIRED_ROWS - 1}")
    return _SYNTHETIC.format(row=row)


def _rows(a3m: str) -> List[List[str]]:
    rows: List[List[str]] = []
    for line in (a3m or "").splitlines():
        if line.startswith(">"):
            rows.append([line[1:], ""])
        elif line.strip() and rows:
            rows[-1][1] += line.strip()
    return rows


def check_pairing_preconditions(paired_msas: Sequence[str], num_chains: int) -> int:
    """Return the shared depth, or refuse if row-order pairing is undefined."""
    if len(paired_msas) != num_chains:
        raise PairingError(
            f"{len(paired_msas)} paired block(s) for {num_chains} chain(s); "
            f"row-order pairing needs one per chain"
        )
    blank = [i for i, m in enumerate(paired_msas) if not (m and m.strip())]
    if blank:
        raise PairingError(
            f"chain(s) {', '.join(map(str, blank))} have no paired block, so every paired "
            f"row would be dropped; re-run the MSA search or use --pair-mode unpaired"
        )
    depths = {len(_rows(m)) for m in paired_msas}
    if len(depths) != 1:
        raise PairingError(
            f"paired blocks have unequal depths {sorted(depths)}; row *i* of one block has "
            f"no partner in another and pairing them would be wrong"
        )
    return depths.pop()


def rewrite_paired_descriptions(paired_msas: Sequence[str], num_chains: int) -> List[str]:
    """Give row *i* of every block the same synthetic species id."""
    depth = check_pairing_preconditions(paired_msas, num_chains)
    if depth > MAX_PAIRED_ROWS:
        logger.warning(f"keeping the first {MAX_PAIRED_ROWS} paired rows of {depth}")
        depth = MAX_PAIRED_ROWS

    out = []
    for a3m in paired_msas:
        lines = []
        for row, (_description, sequence) in enumerate(_rows(a3m)[:depth]):
            # row 0 is the query, which must keep an empty species id so it is not paired
            lines.append(f">{101 + row}" if row == 0 else f">{synthetic_description(row)}")
            lines.append(sequence)
        out.append("\n".join(lines) + "\n")
    return out
