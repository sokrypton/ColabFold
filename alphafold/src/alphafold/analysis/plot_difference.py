import logging
import numpy as np
import matplotlib.pyplot as plt

from pathlib import Path
from typing import Optional, List, Union


logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
logger = logging.getLogger(__name__)


def negative_only(x: float) -> float:
    """Return the value if negative, otherwise return 0.

    Purpose:
    - Zero out positive attention-difference values so plots show only decreases.

    Args:
    - x: A single float attention-difference value.

    Returns:
    - float: x if x < 0 else 0.
    """
    if x < 0:
        return x
    else:
        return 0


def create_custom_xticks_top(x, sequence, interval=5):
    """Build top x-tick labels showing indices above residue letters.

    Labels are formatted so the numeric residue index appears on the first
    (top) line every `interval` residues and the single-letter residue code
    appears on the second line for every tick.

    Args:
    - x: iterable of 1-based residue indices.
    - sequence: string of amino-acid one-letter codes (1-based indexing).
    - interval: int, show numeric index only for residues where index % interval == 0.

    Returns:
    - List[str]: labels for use with ax.set_xticklabels(...).
    """
    labels = [
        f"{i}\n{sequence[i-1]}" if i % interval == 0 else sequence[i - 1] for i in x
    ]
    return labels


def create_custom_xticks_bottom(residue_indices, sequence, label_interval=5):
    """Build bottom x-tick labels with residue letter and periodic numeric index.

    Each label contains the residue letter on the first line. The numeric index is
    placed on the second line only for residues where index % label_interval == 0.
    Non-index ticks include a blank second line to keep vertical alignment.

    Args:
    - residue_indices: iterable of 1-based residue indices.
    - sequence: amino-acid sequence string.
    - label_interval: int, how often to display the numeric index on the second line.

    Returns:
    - List[str]: labels formatted as "AA\\nindex" or "AA\\n".
    """
    labels = []
    for i in residue_indices:
        res = sequence[i - 1]
        labels.append(f"{res}\n{i}" if i % label_interval == 0 else f"{res}\n")
    return labels


def compute_bar_colors(
    residue_indices: np.ndarray,
    query_highlight_positions: Optional[List[int]] = None,
    target_highlight_positions: Optional[List[int]] = None,
    query_highlight_color: str = "#AE0639",
    target_highlight_color: str = "#1f77b4",
    default_color: str = "#D3D3D3",
) -> List[str]:
    """Return a per-residue list of colors for bar plotting.

    Args:
        residue_indices: 1-based residue indices (array-like).
        query_highlight_positions: optional list of 1-based indices to mark with
            query_highlight_color.
        target_highlight_positions: optional list of 1-based indices to mark with
            target_highlight_color.
        query_highlight_color: color string used for query-highlighted positions.
        target_highlight_color: color string used for target-highlighted positions.
        default_color: color string used for all other positions.

    Returns:
        List[str]: color for each index in residue_indices. If a position appears
        in both target and query highlight lists, target_highlight_color takes
        precedence.
    """
    colors = []
    for pos in residue_indices:
        if target_highlight_positions and pos in target_highlight_positions:
            colors.append(target_highlight_color)
        elif query_highlight_positions and pos in query_highlight_positions:
            colors.append(query_highlight_color)
        else:
            colors.append(default_color)
    return colors


def plot_attention(
    attention_scores: np.ndarray,
    highlighted_scores: np.ndarray,
    protein_name: str,
    output_dir: Union[str, Path],
    sequence: Optional[str] = None,
    query_highlight_positions: Optional[List[int]] = None,
    target_highlight_positions: Optional[List[int]] = None,
    query_highlight_color: str = None,
    target_highlight_color: str = None,
) -> None:
    """Create and save an average-attention bar plot with optional highlights.

    The function builds a bar plot of attention_scores (one bar per residue),
    overlays highlighted_scores if provided, and optionally labels the x-axis
    with residue letters and periodic indices (bottom ticks).

    Args:
        attention_scores: 1D numpy array of per-residue attention values.
        highlighted_scores: 1D numpy array (same length) plotted on top of the
            base bars to indicate special regions (may be None).
        protein_name: used for the plot title and output filename.
        output_dir: directory (str or Path) where the PNG will be written.
        sequence: optional amino-acid sequence (string). When provided, x-tick
            labels show residue letters (and periodic numeric indices).
        query_highlight_positions: optional list of 1-based indices to color as query.
        target_highlight_positions: optional list of 1-based indices to color as target.
        query_highlight_color: color string for query highlights (fallback used if None).
        target_highlight_color: color string for target highlights (fallback used if None).

    Side effects:
        Saves a PNG named "{protein_name}_average_attention.png" at 600 dpi in output_dir.
    """
    residue_indices = np.arange(1, attention_scores.size + 1)
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    bar_colors = compute_bar_colors(
        residue_indices,
        query_highlight_positions=query_highlight_positions,
        target_highlight_positions=target_highlight_positions,
        query_highlight_color=query_highlight_color or "#AE0639",
        target_highlight_color=target_highlight_color or "#1f77b4",
        default_color="#D3D3D3",
    )

    fig, ax = plt.subplots(figsize=(8, 6))
    ax.bar(residue_indices, attention_scores, color=bar_colors, zorder=1, alpha=0.8)

    if highlighted_scores is not None:
        ax.bar(
            residue_indices, highlighted_scores, color="#AE0639", zorder=2, alpha=0.6
        )

    if sequence:
        ax.set_xticks(residue_indices)
        ax.set_xticklabels(create_custom_xticks_bottom(residue_indices, sequence))

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
    """Create and save an attention-difference bar plot showing decreases only.

    Positive differences are zeroed (kept only negative values) so the plot
    highlights decreases in attention. Bars can be colored using the same
    highlight parameters as plot_attention. X-axis tick labels are placed on
    the top of the plot when a sequence is provided.

    Args:
        attn_diff_scores: 1D numpy array of per-residue attention differences.
            Positive values are ignored (set to 0) before plotting.
        protein_name: used for the plot title and output filename.
        output_dir: directory (str or Path) where the PNG will be written.
        sequence: optional amino-acid sequence (string). When provided, x-tick
            labels are drawn on the top with periodic numeric indices.
        query_highlight_positions: optional list of 1-based indices to color as query.
        target_highlight_positions: optional list of 1-based indices to color as target.
        query_highlight_color: color string for query highlights.
        target_highlight_color: color string for target highlights.

    Side effects:
        Saves a PNG named "{protein_name}_attention_difference.png" at 600 dpi in output_dir.
    """
    # Save only negative values; set positive values to 0 LLP
    negative_attention_diff_scores = [negative_only(x) for x in attn_diff_scores]
    logger.info(
        "Negative attention difference scores: %s", negative_attention_diff_scores
    )

    residue_indices = np.arange(1, len(attn_diff_scores) + 1)
    logger.info("Residue indices: %s", residue_indices)
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    fig, ax = plt.subplots(figsize=(8, 6))

    ax.xaxis.set_ticks_position("top")
    ax.xaxis.set_label_position("top")
    ax.tick_params(
        axis="x", which="both", top=True, labeltop=True, bottom=False, labelbottom=False
    )

    bar_colors = compute_bar_colors(
        residue_indices,
        query_highlight_positions=query_highlight_positions,
        target_highlight_positions=target_highlight_positions,
        query_highlight_color=query_highlight_color,
        target_highlight_color=target_highlight_color,
        default_color="gray",
    )

    # Use negative values (LLP)
    ax.bar(residue_indices, negative_attention_diff_scores, color=bar_colors)

    if sequence:
        ax.set_xticks(residue_indices)
        ax.set_xticklabels(create_custom_xticks_top(residue_indices, sequence))
    else:
        ax.set_xticks(residue_indices)

    ax.set_ylabel("Attention Difference")
    ax.set_title(f"Attention Difference: {protein_name}")

    ax.spines["right"].set_visible(False)
    ax.spines["bottom"].set_visible(False)

    fig.savefig(output_path / f"{protein_name}_attention_difference.png", dpi=600)
    plt.close(fig)
