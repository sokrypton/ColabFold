import numpy as np
import matplotlib.pyplot as plt

from pathlib import Path
from typing import Optional, List, Union


def create_custom_xticks(
    residue_indices: np.ndarray, sequence: str, label_interval: int = 5
) -> List[str]:
    """Generates formatted labels containing the amino acid and its index."""
    labels = [
        f"{sequence[i-1]}\n{i}" if i % label_interval == 0 else sequence[i - 1]
        for i in residue_indices
    ]
    return labels


def plot_attention(
    attention_scores: np.ndarray,
    highlighted_scores: np.ndarray,
    protein_name: str,
    output_dir: Union[str, Path],
    sequence: Optional[str] = None,
) -> None:
    """Plots average attention per residue with important residues highlighted."""
    residue_indices = np.arange(1, attention_scores.size + 1)
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    plt.figure(figsize=(8, 6))
    plt.bar(residue_indices, attention_scores, color="gray", zorder=1)
    plt.bar(residue_indices, highlighted_scores, color="#AE0639", zorder=2)

    if sequence:
        plt.xticks(residue_indices, create_custom_xticks(residue_indices, sequence))

    plt.xlabel("Amino Acid Residue")
    plt.ylabel("Average Attention Score")
    plt.title(f"Attention Analysis: {protein_name}")

    plt.savefig(output_path / f"{protein_name}_average_attention.png", dpi=600)
    plt.close()


def plot_difference(
    attn_diff_scores: np.ndarray,
    protein_name: str,
    output_dir: Union[str, Path],
    sequence: Optional[str] = None,
) -> None:
    """Plots the attention difference between two proteins or states."""
    residue_indices = np.arange(1, len(attn_diff_scores) + 1)
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    plt.figure(figsize=(8, 6))
    plt.bar(residue_indices, attn_diff_scores, color="gray")

    if sequence:
        plt.xticks(residue_indices, create_custom_xticks(residue_indices, sequence))

    plt.xlabel("Amino Acid Residue")
    plt.ylabel("Attention Difference")
    plt.title(f"Attention Difference: {protein_name}")

    plt.savefig(output_path / f"{protein_name}_attention_difference.png", dpi=600)
    plt.close()
