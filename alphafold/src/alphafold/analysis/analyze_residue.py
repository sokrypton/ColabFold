import logging
import numpy as np
import pandas as pd

from Bio import Align
from Bio.Align import substitution_matrices
from typing import Tuple
from pathlib import Path


logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
logger = logging.getLogger(__name__)


def consecutive(data: np.ndarray, stepsize: int = 1) -> list:
    """Split a sorted 1D array into runs of consecutive values.

    Consecutive is a small convenience that groups indices (or values) that
    differ by exactly `stepsize` into contiguous runs.

    Args:
        data: 1D numpy array of integer-like values (assumed sorted).
        stepsize: step size that determines consecutiveness (default 1).

    Returns:
        List of numpy.ndarray objects, each containing a consecutive run.
    """
    return np.split(data, np.where(np.diff(data) != stepsize)[0] + 1)


def find_highest_attention(
    attention: np.ndarray, sequence: str, output_dir, protein: str
) -> None:
    """Write a CSV ranking residues by attention score.

    Produces a CSV with columns: Rank, Residue number (1-based), Amino acid, and
    Attention score. Rows are ordered from highest to lowest attention.

    Args:
        attention: 1D numpy array of per-residue attention scores.
        sequence: amino-acid sequence string whose length matches attention.
        output_dir: directory path (str or Path) where CSV will be written.
        protein: base name used for the output filename.

    Side effects:
        Creates output_dir (if needed) and writes "{protein}_residue_ranking.csv".
    """
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)

    descending_indices = (-attention).argsort()
    sorted_attention = attention[descending_indices]
    sorted_residue_numbers = descending_indices + 1
    sorted_sequence = np.array(list(sequence))[descending_indices]

    df = pd.DataFrame(
        {
            "Rank": np.arange(1, len(sequence) + 1),
            "Residue number": sorted_residue_numbers,
            "Amino acid": sorted_sequence,
            "Attention score": sorted_attention,
        }
    )

    out_file = out / f"{protein}_residue_ranking.csv"
    df.to_csv(out_file, index=False)
    logger.info("Wrote residue ranking to %s", out_file)


def find_important(attention: np.ndarray, zscores: np.ndarray = None) -> np.ndarray:
    """Identify a small set of top residues by attention.

    Current heuristic: select the top 5 residues by attention score and return
    an array with non-zero values only at those positions.

    Args:
        attention: 1D numpy array of per-residue attention scores.
        zscores: optional array of z-scores (currently unused, reserved for future).

    Returns:
        1D numpy array same shape as `attention` with original attention values
        at the selected top-5 indices and zeros elsewhere.
    """
    important = np.zeros_like(attention)
    top_indices = np.argsort(attention)[-5:]
    important[top_indices] = attention[top_indices]
    return important


def blosum_scores(sequence1: str, sequence2: str) -> Tuple[np.ndarray, np.ndarray]:
    """Compute per-position BLOSUM62 substitution scores for two aligned sequences.

    The function expects two sequences of equal length (alignment including gaps).
    Gap characters are handled by mapping '-' to '*' for lookup into the matrix.

    Args:
        sequence1: aligned amino-acid sequence (may contain '-' for gaps).
        sequence2: aligned amino-acid sequence (same length as sequence1).

    Returns:
        Tuple of two numpy arrays (scores1, scores2) containing the per-position
        substitution score for the corresponding residue in each sequence.
    """
    aligner = Align.PairwiseAligner()
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    scores1 = np.zeros([len(sequence1)])
    scores2 = np.zeros([len(sequence2)])

    for i, (s1, s2) in enumerate(zip(sequence1, sequence2)):
        if s1 == "-" and s2 == "-":
            s1_temp = "*"
            s2_temp = "*"
            scores1[i] = aligner.substitution_matrix[s1_temp, s2_temp]
            scores2[i] = aligner.substitution_matrix[s1_temp, s2_temp]
        elif s1 == "-":
            s1_temp = "*"
            scores1[i] = aligner.substitution_matrix[s1_temp, s2]
        elif s2 == "-":
            s2_temp = "*"
            scores2[i] = aligner.substitution_matrix[s1, s2_temp]
        else:
            scores1[i] = aligner.substitution_matrix[s1, s2]
            scores2[i] = aligner.substitution_matrix[s1, s2]
    return scores1, scores2


def calculate_differences(
    scores1: np.ndarray,
    scores2: np.ndarray,
    attention1: np.ndarray,
    attention2: np.ndarray,
    gaps1,
    gaps2,
    bool_alignment: bool = False,
) -> Tuple[list, list]:
    """Compute adjusted attention-difference arrays with edge handling and BLOSUM weighting.

    This function computes two importance-weighted difference arrays (for sequence1
    and sequence2). When `bool_alignment` is True, positions corresponding to gaps
    in the opposite sequence are zeroed out before further processing.

    The output arrays are scaled by the provided BLOSUM scores, and edge regions
    (long runs of zeros at the beginning or end) are accounted for via offsets to
    avoid spurious signals near alignment ends.

    Args:
        scores1: per-position BLOSUM (or similarity) scores for sequence1.
        scores2: per-position BLOSUM (or similarity) scores for sequence2.
        attention1: per-residue attention values for sequence1.
        attention2: per-residue attention values for sequence2.
        gaps1: iterable of indices (0-based) marking gaps in sequence1 (used when bool_alignment True).
        gaps2: iterable of indices (0-based) marking gaps in sequence2 (used when bool_alignment True).
        bool_alignment: if True treat gaps by zeroing positions in the opposite-difference arrays.

    Returns:
        A tuple (important_diff_blosum1, important_diff_blosum2) where each element
        is a list (same length as the input differences) containing the BLOSUM-weighted
        importance scores with edge offsets applied.
    """
    if bool_alignment:
        difference1 = attention1 - attention2
        important_diff1 = np.where(difference1 > 0, difference1, 0)
        for g in gaps2:
            important_diff1[g] = 0
        difference2 = attention2 - attention1
        important_diff2 = np.where(difference2 > 0, difference2, 0)
        for g in gaps1:
            important_diff2[g] = 0
    else:
        difference1 = attention1 - attention2
        important_diff1 = np.where(difference1 > 0, difference1, 0)

        difference2 = attention2 - attention1
        important_diff2 = np.where(difference2 > 0, difference2, 0)

    # Account for Edge Effects LLP
    beginning_offset = 0
    ending_offset = len(difference1)

    consecutive1 = consecutive(difference1)
    consecutive2 = consecutive(difference2)

    logger.info(
        "Beginning offset: %d, Ending offset: %d", beginning_offset, ending_offset
    )

    if [len(x) for x in consecutive1 if 0 in x]:
        beginning_offset = [len(x) for x in consecutive1 if 0 in x][0] + 5

    if [len(x) for x in consecutive2 if 0 in x]:
        beginning_offset = [len(x) for x in consecutive2 if 0 in x][0] + 5

    if [len(x) for x in consecutive1 if len(difference1) - 1 in x]:
        ending_offset = [len(x) for x in consecutive1 if len(difference1) - 1 in x][
            0
        ] - 4

    if [len(x) for x in consecutive2 if len(difference2) - 1 in x]:
        ending_offset = [max(x) for x in consecutive2 if len(difference2) - 1 in x][
            0
        ] - 4

    logger.info(
        "Beginning offset: %d, Ending offset: %d", beginning_offset, ending_offset
    )

    # Update blosum with edge effects LLP
    important_diff_blosum1 = (
        [0] * beginning_offset
        + [diff * s for diff, s in zip(important_diff1, scores1)][
            beginning_offset:ending_offset
        ]
        + [0] * (len(difference1) - ending_offset)
    )
    important_diff_blosum2 = (
        [0] * beginning_offset
        + [diff * s for diff, s in zip(important_diff2, scores2)][
            beginning_offset:ending_offset
        ]
        + [0] * (len(difference1) - ending_offset)
    )
    return important_diff_blosum1, important_diff_blosum2
