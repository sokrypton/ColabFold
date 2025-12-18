import logging
import numpy as np
import pandas as pd

from Bio import Align
from Bio.Align import substitution_matrices
from typing import Tuple
from pathlib import Path


logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
logger = logging.getLogger(__name__)


def find_highest_attention(
    attention: np.ndarray, sequence: str, output_dir, protein: str
) -> None:
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
    important = np.zeros_like(attention)
    top_indices = np.argsort(attention)[-5:]
    important[top_indices] = attention[top_indices]

    return important


def blosum_scores(sequence1: str, sequence2: str) -> Tuple[np.ndarray, np.ndarray]:
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

    important_diff_blosum1 = [diff * s for diff, s in zip(important_diff1, scores1)]
    important_diff_blosum2 = [diff * s for diff, s in zip(important_diff2, scores2)]

    return important_diff_blosum1, important_diff_blosum2
