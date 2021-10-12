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


def plot_lddt(
    jobname: str,
    msa,
    outs: dict,
    query_sequence,
    result_dir: Path,
    show: bool = False,
):
    # gather MSA info
    deduped_full_msa = list(dict.fromkeys(msa))
    msa_arr = np.array([list(seq) for seq in deduped_full_msa])

    if isinstance(query_sequence, str):
        query_str = query_sequence
        query_len = len(query_sequence)
    else:
        query_str = "".join(query_sequence)
        query_len = sum(len(s) for s in query_sequence)

    seqid = (np.array(list(query_str)) == msa_arr).mean(-1)
    seqid_sort = seqid.argsort()  # [::-1]
    non_gaps = (msa_arr != "-").astype(float)
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
    plt.plot((msa_arr != "-").sum(0), color="black")
    plt.xlim(-0.5, msa_arr.shape[1] - 0.5)
    plt.ylim(-0.5, msa_arr.shape[0] - 0.5)
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
