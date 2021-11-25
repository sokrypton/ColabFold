from pathlib import Path

import numpy as np
from matplotlib import pyplot as plt
from matplotlib import collections as mcoll
import matplotlib
from string import ascii_uppercase, ascii_lowercase

from colabfold.colabfold import (
    pymol_color_list,
    pymol_cmap,
    alphabet_list,
    kabsch,
    plot_pseudo_3D,
    add_text,
    plot_ticks,
)


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


def plot_protein_confidence(
    plot_path,
    protein=None,
    pos=None,
    plddt=None,
    pae=None,
    Ls=None,
    dpi=200,
    best_view=True,
    line_w=2.0,
    show=False,
):

    use_ptm = pae is not None

    if protein is not None:
        pos = np.asarray(protein.atom_positions[:, 1, :])
        plddt = np.asarray(protein.b_factors[:, 0])

    # get best view
    if best_view:
        if plddt is not None:
            weights = plddt / 100
            pos = pos - (pos * weights[:, None]).sum(0, keepdims=True) / weights.sum()
            pos = pos @ kabsch(pos, pos, weights, return_v=True)
        else:
            pos = pos - pos.mean(0, keepdims=True)
            pos = pos @ kabsch(pos, pos, return_v=True)

    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(
        nrows=2, ncols=2, gridspec_kw={"height_ratios": [1.5, 1]}
    )
    fig.set_figwidth(7)
    fig.set_figheight(7)
    ax = {"prot_chain": ax1, "prot_plddt": ax2, "IDDT": ax3, "Ali_error": ax4}

    fig.set_dpi(dpi)
    fig.subplots_adjust(top=0.9, bottom=0.1, right=1, left=0.1, hspace=0, wspace=0.05)

    # 3D PLOT:
    xy_min = pos[..., :2].min() - line_w
    xy_max = pos[..., :2].max() + line_w
    for a in [ax["prot_chain"], ax["prot_plddt"]]:
        a.set_xlim(xy_min, xy_max)
        a.set_ylim(xy_min, xy_max)
        a.axis(False)

    if Ls is None or len(Ls) == 1:
        # color N->C
        c = np.arange(len(pos))[::-1]
        plot_pseudo_3D(pos, line_w=line_w, ax=ax1)
        add_text("colored by N→C", ax1)
    else:
        # color by chain
        c = np.concatenate([[n] * L for n, L in enumerate(Ls)])
        if len(Ls) > 40:
            plot_pseudo_3D(pos, c=c, line_w=line_w, ax=ax1)
        else:
            plot_pseudo_3D(
                pos, c=c, cmap=pymol_cmap, cmin=0, cmax=39, line_w=line_w, ax=ax1
            )
        add_text("colored by chain", ax1)

    if plddt is not None:
        # color by pLDDT
        plot_pseudo_3D(pos, c=plddt, cmin=50, cmax=90, line_w=line_w, ax=ax2)
        add_text("colored by pLDDT", ax2)

    # Conf plot:
    ax["IDDT"].set_title("Predicted lDDT")
    ax["IDDT"].plot(plddt)
    if Ls is not None:
        L_prev = 0
        for L_i in Ls[:-1]:
            L = L_prev + L_i
            L_prev += L_i
            ax["IDDT"].plot([L, L], [0, 100], color="black")
    ax["IDDT"].set_ylim(0, 100)
    ax["IDDT"].set_ylabel("plDDT")
    ax["IDDT"].set_xlabel("position")

    if use_ptm:
        ax["Ali_error"].set_title("Predicted Aligned Error")
        Ln = pae.shape[0]
        ax["Ali_error"].imshow(pae, cmap="bwr", vmin=0, vmax=30, extent=(0, Ln, Ln, 0))
        if Ls is not None and len(Ls) > 1:
            plot_ticks(Ls)
        # ax['Ali_error'].colorbar()
        ax["Ali_error"].set_xlabel("Scored residue")
        ax["Ali_error"].set_ylabel("Aligned residue")

    plt.savefig(plot_path)
    if show:
        plt.show()
    plt.close()
