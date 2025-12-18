import numpy as np
import matplotlib.pyplot as plt

from pathlib import Path
from typing import Optional, List, Union

#LLP convert positive values to 0
def negative_only(x):
    if x < 0:
        return x
    else:
        return 0

#LLP custom ticks for difference plots
def create_custom_xticks_top(x, sequence, interval=5):
    # Use sequence[i-1]because x is 1-indexed
    labels = [f'{i}\n{sequence[i-1]}' if i % interval == 0 else sequence[i-1] for i in x]
    return labels

def create_custom_xticks(
    residue_indices: np.ndarray, sequence: str, label_interval: int = 5
) -> List[str]:
    """Generates formatted labels containing the amino acid and its index."""
    labels = [
        f"{sequence[i-1]}\n{i}" if i % label_interval == 0 else sequence[i - 1]
        for i in residue_indices
    ]
    return labels


def color_xtick_labels(
    residue_indices: np.ndarray,
    sequence: str,
    ax: plt.Axes,
    query_highlight_positions: Optional[List[int]] = None,
    target_highlight_positions: Optional[List[int]] = None,
    query_highlight_color: str = "#AE0639",
    target_highlight_color: str = "#1f77b4",
) -> None:
    """Set x-tick label text and color amino-acid letters at provided 1-based indices.

    If a position is present in both highlight lists, target_highlight_color takes precedence.
    """
    ax.set_xticks(residue_indices)
    ax.set_xticklabels(create_custom_xticks(residue_indices, sequence))
    tick_labels = ax.get_xticklabels()

    # default color
    for lbl in tick_labels:
        lbl.set_color("black")

    if query_highlight_positions:
        for pos in query_highlight_positions:
            if 1 <= pos <= len(tick_labels):
                tick_labels[pos - 1].set_color(query_highlight_color)

    if target_highlight_positions:
        for pos in target_highlight_positions:
            if 1 <= pos <= len(tick_labels):
                tick_labels[pos - 1].set_color(target_highlight_color)


def plot_attention(
    attention_scores: np.ndarray,
    highlighted_scores: np.ndarray,
    protein_name: str,
    output_dir: Union[str, Path],
    sequence: Optional[str] = None,
    query_highlight_positions: Optional[List[int]] = None,
    target_highlight_positions: Optional[List[int]] = None,
    query_highlight_color: str = "#AE0639",
    target_highlight_color: str = "#1f77b4",
) -> None:
    """Plots average attention per residue with important residues highlighted.

    Note: amino-acid letters on the x-axis are color-coded (not the bars).
    """
    residue_indices = np.arange(1, attention_scores.size + 1)
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    fig, ax = plt.subplots(figsize=(8, 6))
    ax.bar(residue_indices, attention_scores, color="gray", zorder=1)

    if highlighted_scores is not None:
        ax.bar(
            residue_indices, highlighted_scores, color="#AE0639", zorder=2, alpha=0.6
        )

    if sequence:
        color_xtick_labels(
            residue_indices,
            sequence,
            ax,
            query_highlight_positions=query_highlight_positions,
            target_highlight_positions=target_highlight_positions,
            query_highlight_color=query_highlight_color,
            target_highlight_color=target_highlight_color,
        )
    else:
        ax.set_xticks(residue_indices)

    ax.set_xlabel("Amino Acid Residue")
    ax.set_ylabel("Average Attention Score")
    ax.set_title(f"Attention Analysis: {protein_name}")

    fig.savefig(output_path / f"{protein_name}_average_attention.png", dpi=600)
    plt.close(fig)


def plot_difference(
    attn_diff_scores: np.ndarray,
    protein_name: str,
    output_dir: Union[str, Path],
    sequence: Optional[str] = None,
    query_highlight_positions: Optional[List[int]] = None,
    target_highlight_positions: Optional[List[int]] = None,
    query_highlight_color: str = "#AE0639",
    target_highlight_color: str = "#1f77b4",
) -> None:
    """Plots the attention difference between two proteins or states.

    Note: amino-acid letters on the x-axis are color-coded (not the bars).
    """

    #Save only negative values; set positive values to 0 LLP
    negative_attention_diff_scores = [negative_only(x) for x in attn_diff_scores]
    
    residue_indices = np.arange(1, len(attn_diff_scores) + 1)
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    fig, ax = plt.subplots(figsize=(8, 6))

    #Use negative values (LLP)
    ax.bar(residue_indices, negative_attention_diff_scores, color="gray")

    if sequence:
        #Use create_custom_xticks_top to make proper ticks for this function
        color_xtick_labels(
            residue_indices,
            sequence,
            ax,
            query_highlight_positions=query_highlight_positions,
            target_highlight_positions=target_highlight_positions,
            query_highlight_color=query_highlight_color,
            target_highlight_color=target_highlight_color,
        )
    else:
        ax.set_xticks(residue_indices)

    #Get rid of bottom spine LLP
    #ax.set_xlabel("Amino Acid Residue")
    ax.set_ylabel("Attention Difference")
    ax.set_title(f"Attention Difference: {protein_name}")

    # Hide the right and bottom spines LLP
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(False)

    fig.savefig(output_path / f"{protein_name}_attention_difference.png", dpi=600)
    plt.close(fig)
