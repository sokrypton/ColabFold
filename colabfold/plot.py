from pathlib import Path

import numpy as np
from matplotlib import pyplot as plt
from matplotlib import collections as mcoll
import matplotlib
from string import ascii_uppercase, ascii_lowercase

pymol_color_list = [
    "#33ff33",
    "#00ffff",
    "#ff33cc",
    "#ffff00",
    "#ff9999",
    "#e5e5e5",
    "#7f7fff",
    "#ff7f00",
    "#7fff7f",
    "#199999",
    "#ff007f",
    "#ffdd5e",
    "#8c3f99",
    "#b2b2b2",
    "#007fff",
    "#c4b200",
    "#8cb266",
    "#00bfbf",
    "#b27f7f",
    "#fcd1a5",
    "#ff7f7f",
    "#ffbfdd",
    "#7fffff",
    "#ffff7f",
    "#00ff7f",
    "#337fcc",
    "#d8337f",
    "#bfff3f",
    "#ff7fff",
    "#d8d8ff",
    "#3fffbf",
    "#b78c4c",
    "#339933",
    "#66b2b2",
    "#ba8c84",
    "#84bf00",
    "#b24c66",
    "#7f7f7f",
    "#3f3fa5",
    "#a5512b",
]

pymol_cmap = matplotlib.colors.ListedColormap(pymol_color_list)
alphabet_list = list(ascii_uppercase + ascii_lowercase)


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


def plot_lddt(
    jobname: str, msa, outs: dict, query_sequence, result_dir: Path, show: bool = False
):
    # gather MSA info
    seqid = (query_sequence == msa).mean(-1)
    seqid_sort = seqid.argsort()  # [::-1]
    non_gaps = (msa != 21).astype(float)
    non_gaps[non_gaps == 0] = np.nan

    plt.figure(figsize=(14, 4), dpi=100)

    plt.subplot(1, 2, 1)
    plt.title("Sequence coverage")
    plt.imshow(
        non_gaps[seqid_sort] * seqid[seqid_sort, None],
        interpolation="nearest",
        aspect="auto",
        cmap="rainbow_r",
        vmin=0,
        vmax=1,
        origin="lower",
    )
    plt.plot((msa != 21).sum(0), color="black")
    plt.xlim(-0.5, msa.shape[1] - 0.5)
    plt.ylim(-0.5, msa.shape[0] - 0.5)
    plt.colorbar(label="Sequence identity to query")
    plt.xlabel("Positions")
    plt.ylabel("Sequences")

    plt.subplot(1, 2, 2)
    plt.title("Predicted lDDT per position")
    for model_name, value in outs.items():
        plt.plot(value["plddt"], label=model_name)

    plt.legend()
    plt.ylim(0, 100)
    plt.ylabel("Predicted lDDT")
    plt.xlabel("Positions")
    plt.savefig(str(result_dir.joinpath(jobname + "_coverage_lDDT.png")))
    if show:
        plt.show()
    plt.close()


def kabsch(a, b, weights=None, return_v=False):
    a = np.asarray(a)
    b = np.asarray(b)
    if weights is None:
        weights = np.ones(len(b))
    else:
        weights = np.asarray(weights)
    B = np.einsum("ji,jk->ik", weights[:, None] * a, b)
    u, s, vh = np.linalg.svd(B)
    if np.linalg.det(u @ vh) < 0:
        u[:, -1] = -u[:, -1]
    if return_v:
        return u
    else:
        return u @ vh


def plot_pseudo_3D(
    xyz,
    c=None,
    ax=None,
    chainbreak=5,
    cmap="gist_rainbow",
    line_w=2.0,
    cmin=None,
    cmax=None,
    zmin=None,
    zmax=None,
):
    def rescale(a, amin=None, amax=None):
        a = np.copy(a)
        if amin is None:
            amin = a.min()
        if amax is None:
            amax = a.max()
        a[a < amin] = amin
        a[a > amax] = amax
        return (a - amin) / (amax - amin)

    # make segments
    xyz = np.asarray(xyz)
    seg = np.concatenate([xyz[:-1, None, :], xyz[1:, None, :]], axis=-2)
    seg_xy = seg[..., :2]
    seg_z = seg[..., 2].mean(-1)
    ord = seg_z.argsort()

    # set colors
    if c is None:
        c = np.arange(len(seg))[::-1]
    else:
        c = (c[1:] + c[:-1]) / 2
    c = rescale(c, cmin, cmax)

    if isinstance(cmap, str):
        if cmap == "gist_rainbow":
            c *= 0.75
        colors = matplotlib.cm.get_cmap(cmap)(c)
    else:
        colors = cmap(c)

    if chainbreak is not None:
        dist = np.linalg.norm(xyz[:-1] - xyz[1:], axis=-1)
        colors[..., 3] = (dist < chainbreak).astype(np.float)

    # add shade/tint based on z-dimension
    z = rescale(seg_z, zmin, zmax)[:, None]
    tint, shade = z / 3, (z + 2) / 3
    colors[:, :3] = colors[:, :3] + (1 - colors[:, :3]) * tint
    colors[:, :3] = colors[:, :3] * shade

    set_lim = False
    if ax is None:
        fig, ax = plt.subplots()
        fig.set_figwidth(5)
        fig.set_figheight(5)
        set_lim = True
    else:
        fig = ax.get_figure()
        if ax.get_xlim() == (0, 1):
            set_lim = True

    if set_lim:
        xy_min = xyz[:, :2].min() - line_w
        xy_max = xyz[:, :2].max() + line_w
        ax.set_xlim(xy_min, xy_max)
        ax.set_ylim(xy_min, xy_max)

    ax.set_aspect("equal")

    # determine linewidths
    width = fig.bbox_inches.width * ax.get_position().width
    linewidths = line_w * 72 * width / np.diff(ax.get_xlim())

    lines = mcoll.LineCollection(
        seg_xy[ord],
        colors=colors[ord],
        linewidths=linewidths,
        path_effects=[matplotlib.patheffects.Stroke(capstyle="round")],
    )

    return ax.add_collection(lines)


def add_text(text, ax):
    return plt.text(
        0.5,
        1.01,
        text,
        horizontalalignment="center",
        verticalalignment="bottom",
        transform=ax.transAxes,
    )


def plot_ticks(Ls):
    Ln = sum(Ls)
    L_prev = 0
    for L_i in Ls[:-1]:
        L = L_prev + L_i
        L_prev += L_i
        plt.plot([0, Ln], [L, L], color="black")
        plt.plot([L, L], [0, Ln], color="black")
    ticks = np.cumsum([0] + Ls)
    ticks = (ticks[1:] + ticks[:-1]) / 2
    plt.yticks(ticks, alphabet_list[: len(ticks)])


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


# fig = plot_protein_confidence(unrelaxed_protein,
#                              pae=prediction_result['predicted_aligned_error'],
#                              Ls=query_sequence_len_array)
