from pathlib import Path
import numpy as np

def plot_predicted_alignment_error(
    jobname: str, num_models: int, outs: dict, result_dir: Path, show: bool = False
):
    from matplotlib import pyplot as plt
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
    from matplotlib import pyplot as plt
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

def _msa_lines(msa, query_sequence, seq_len_list):
    """Rows of per-chain identity to the query (nan at gaps), grouped and sorted for display.

    Consecutive rows that cover the same set of chains form a group, and each group is
    sorted by its best identity.
    """
    msa = np.asarray(msa)
    query_sequence = np.asarray(query_sequence)
    starts = np.cumsum(np.append(0, seq_len_list))[:-1]
    parts, has_seq = [], []
    for start, length in zip(starts, seq_len_list):
        chain_msa = msa[:, start:start + length]
        seqid = np.count_nonzero(chain_msa == query_sequence[start:start + length], axis=1) / length
        non_gaps = (chain_msa != 21).astype(float)
        non_gaps[non_gaps == 0] = np.nan
        part = non_gaps * seqid[:, None]
        parts.append(part)
        has_seq.append(~np.isnan(part).all(axis=1))
    parts = np.concatenate(parts, axis=1)
    has_seq = np.stack(has_seq, axis=1)
    # a group ends where the covered set changes; the first row is compared against all chains
    previous = np.vstack([np.ones((1, has_seq.shape[1]), dtype=bool), has_seq[:-1]])
    bounds = [i for i in np.flatnonzero((has_seq != previous).any(axis=1)) if i > 0]
    groups = np.split(parts, bounds)
    return np.concatenate([g[np.argsort(-np.nanmax(g, axis=1))] for g in groups if len(g)])


def plot_msa(msa, query_sequence, seq_len_list, total_seq_len, dpi=100):
    from matplotlib import pyplot as plt
    Ln = np.cumsum(np.append(0, [len for len in seq_len_list]))
    lines = _msa_lines(msa, query_sequence, seq_len_list)

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
