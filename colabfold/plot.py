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


def plot_msa_v2(feature_dict, sort_lines=True, dpi=100):
    seq = feature_dict["msa"][0]
    if "asym_id" in feature_dict:
      Ls = [0]
      k = feature_dict["asym_id"][0]
      for i in feature_dict["asym_id"]:
        if i == k: Ls[-1] += 1
        else: Ls.append(1)
        k = i
    else:
      Ls = [len(seq)]    
    Ln = np.cumsum([0] + Ls)

    try:
        N = feature_dict["num_alignments"][0]
    except:
        N = feature_dict["num_alignments"] 
    
    msa = feature_dict["msa"][:N]
    gap = msa != 21
    qid = msa == seq
    gapid = np.stack([gap[:,Ln[i]:Ln[i+1]].max(-1) for i in range(len(Ls))],-1)
    lines = []
    Nn = []
    for g in np.unique(gapid, axis=0):
        i = np.where((gapid == g).all(axis=-1))
        qid_ = qid[i]
        gap_ = gap[i]
        seqid = np.stack([qid_[:,Ln[i]:Ln[i+1]].mean(-1) for i in range(len(Ls))],-1).sum(-1) / (g.sum(-1) + 1e-8)
        non_gaps = gap_.astype(float)
        non_gaps[non_gaps == 0] = np.nan
        if sort_lines:
            lines_ = non_gaps[seqid.argsort()] * seqid[seqid.argsort(),None]
        else:
            lines_ = non_gaps[::-1] * seqid[::-1,None]
        Nn.append(len(lines_))
        lines.append(lines_)
    
    Nn = np.cumsum(np.append(0,Nn))
    lines = np.concatenate(lines,0)
    plt.figure(figsize=(8,5), dpi=dpi)
    plt.title("Sequence coverage")
    plt.imshow(lines,
              interpolation='nearest', aspect='auto',
              cmap="rainbow_r", vmin=0, vmax=1, origin='lower',
              extent=(0, lines.shape[1], 0, lines.shape[0]))
    for i in Ln[1:-1]:
        plt.plot([i,i],[0,lines.shape[0]],color="black")
    for j in Nn[1:-1]:
        plt.plot([0,lines.shape[1]],[j,j],color="black")
    
    plt.plot((np.isnan(lines) == False).sum(0), color='black')
    plt.xlim(0,lines.shape[1])
    plt.ylim(0,lines.shape[0])
    plt.colorbar(label="Sequence identity to query")
    plt.xlabel("Positions")
    plt.ylabel("Sequences")
    return plt

def plot_msa(msa, query_sequence, seq_len_list, total_seq_len, dpi=100):
    # gather MSA info
    prev_pos = 0
    msa_parts = []
    Ln = np.cumsum(np.append(0, [len for len in seq_len_list]))
    for id, l in enumerate(seq_len_list):
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
        msa_parts.append((non_gaps[:] * seqid[:, None]).tolist())
        prev_pos += l
    lines = []
    lines_to_sort = []
    prev_has_seq = [True] * len(seq_len_list)
    for line_num in range(len(msa_parts[0])):
        has_seq = [True] * len(seq_len_list)
        for id in range(len(seq_len_list)):
            if np.sum(~np.isnan(msa_parts[id][line_num])) == 0:
                has_seq[id] = False
        if has_seq == prev_has_seq:
            line = []
            for id in range(len(seq_len_list)):
                line += msa_parts[id][line_num]
            lines_to_sort.append(np.array(line))
        else:
            lines_to_sort = np.array(lines_to_sort)
            lines_to_sort = lines_to_sort[np.argsort(-np.nanmax(lines_to_sort, axis=1))]
            lines += lines_to_sort.tolist()
            lines_to_sort = []
            line = []
            for id in range(len(seq_len_list)):
                line += msa_parts[id][line_num]
            lines_to_sort.append(line)
        prev_has_seq = has_seq
    lines_to_sort = np.array(lines_to_sort)
    lines_to_sort = lines_to_sort[np.argsort(-np.nanmax(lines_to_sort, axis=1))]
    lines += lines_to_sort.tolist()

    # Nn = np.cumsum(np.append(0, Nn))
    # lines = np.concatenate(lines, 1)
    xaxis_size = len(lines[0])
    yaxis_size = len(lines)

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
        extent=(0, xaxis_size, 0, yaxis_size),
    )
    for i in Ln[1:-1]:
        plt.plot([i, i], [0, yaxis_size], color="black")
    # for i in Ln_dash[1:-1]:
    #    plt.plot([i, i], [0, lines.shape[0]], "--", color="black")
    # for j in Nn[1:-1]:
    #    plt.plot([0, lines.shape[1]], [j, j], color="black")

    plt.plot((np.isnan(lines) == False).sum(0), color="black")
    plt.xlim(0, xaxis_size)
    plt.ylim(0, yaxis_size)
    plt.colorbar(label="Sequence identity to query")
    plt.xlabel("Positions")
    plt.ylabel("Sequences")

    return plt
