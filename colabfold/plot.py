from pathlib import Path

import numpy as np
from matplotlib import pyplot as plt


def plot_predicted_alignment_error(
    jobname: str, num_models: int, outs: dict, result_dir: Path, show: bool = False
):
    plt.figure(figsize=(3 * num_models, 2), dpi=100)
    for n, (model_name, value) in enumerate(outs.items()):
        plt.subplot(1, num_models, n + 1)
        plt.title(model_name)
        plt.imshow(value["pae"], label=model_name, cmap="bwr", vmin=0, vmax=30)
        plt.colorbar()
    plt.savefig(result_dir.joinpath(jobname + "_PAE.png"))
    if show:
        plt.show()
    plt.close()


def plot_msa(msa, query_sequence, seq_len_list, total_seq_len, dpi=100):
    # gather MSA info
    prev_pos = 0
    lines = []
    Ln = np.cumsum(np.append(0, [len for len in seq_len_list]))

    for l in seq_len_list:
        chain_seq = np.array(query_sequence[prev_pos : prev_pos + l])
        chain_msa = np.array(msa[:, prev_pos : prev_pos + l])
        seqid = np.array(
            [
                np.count_nonzero(chain_seq == msa_line[prev_pos : prev_pos + l])
                / len(chain_seq)
                for msa_line in msa
            ]
        )
        non_gaps = (chain_msa != 21).astype(float)
        non_gaps[non_gaps == 0] = np.nan
        lines.append(non_gaps[:] * seqid[:, None])
        prev_pos += l

    # Nn = np.cumsum(np.append(0, Nn))
    lines = np.concatenate(lines, 1)

    plt.figure(figsize=(8, 5), dpi=dpi)
    plt.title("Sequence coverage")
    plt.imshow(
        lines[::-1],
        interpolation="nearest",
        aspect="auto",
        cmap="rainbow_r",
        vmin=0,
        vmax=1,
        origin="lower",
        extent=(0, lines.shape[1], 0, lines.shape[0]),
    )
    for i in Ln[1:-1]:
        plt.plot([i, i], [0, lines.shape[0]], color="black")
    # for i in Ln_dash[1:-1]:
    #    plt.plot([i, i], [0, lines.shape[0]], "--", color="black")
    # for j in Nn[1:-1]:
    #    plt.plot([0, lines.shape[1]], [j, j], color="black")

    plt.plot((np.isnan(lines) == False).sum(0), color="black")
    plt.xlim(0, lines.shape[1])
    plt.ylim(0, lines.shape[0])
    plt.colorbar(label="Sequence identity to query")
    plt.xlabel("Positions")
    plt.ylabel("Sequences")

    return plt
