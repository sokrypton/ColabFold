import os
import sys
import jax
import logging
import linecache
import numpy as np

from typing import List, Tuple

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
logger = logging.getLogger(__name__)


def get_n(folder_path: str) -> int:
    shape_counts = {}

    for fname in os.listdir(folder_path):
        if fname.endswith(".npy"):
            try:
                arr = np.load(os.path.join(folder_path, fname))
                shape = arr.shape
                shape_counts[shape] = shape_counts.get(shape, 0) + 1
            except Exception as e:
                logger.warning("Could not load %s: %s", fname, e)

    if not shape_counts:
        logger.error("No .npy files found in %s", folder_path)
        sys.exit(1)

    sorted_shapes = sorted(shape_counts.items(), key=lambda item: item[1], reverse=True)

    logger.info("Unique shapes found (sorted by frequency):")
    for shape, count in sorted_shapes:
        logger.info("%s: %d files", shape, count)

    most_frequent_shape = sorted_shapes[0][0]

    return most_frequent_shape[0]


def get_attention(folder_path: str, n: int) -> np.ndarray:
    file_list = [fname for fname in os.listdir(folder_path) if fname.endswith(".npy")]
    file_roots = sorted(file_list, key=lambda x: int(x.split("_")[-1].split(".")[0]))

    heads_n_4_n_n = []

    for fname in file_roots:
        arr1 = np.load(os.path.join(folder_path, fname))
        arr1 = arr1.view(dtype=np.float16)
        arr1 = jax.nn.softmax(arr1)
        if arr1.shape == (n, 4, n, n):
            heads_n_4_n_n.append(arr1)

    attention_spectrum = []
    for arr in heads_n_4_n_n:
        attn = np.sum(arr, axis=(0, 1, 2))
        attention_spectrum.append(attn)
    attention_spectrum = np.array(attention_spectrum)

    return attention_spectrum


def average(attention_spectrum: np.ndarray) -> np.ndarray:
    return np.mean(attention_spectrum, axis=0)


def min_max(data: np.ndarray) -> np.ndarray:
    data = np.asarray(data)
    min_val = data.min()
    max_val = data.max()
    if max_val == min_val:
        return np.zeros_like(data)
    return (data - min_val) / (max_val - min_val)


def read_sequence_file(sequence_file: str) -> str:
    with open(sequence_file, "r") as f:
        content = f.readlines()
        start_end = []
        found = 0
        sequence = ""
        for linenum, line in enumerate(content):
            if line.startswith(">"):
                start_end.append(linenum)
                found += 1
            if found == 2:
                break

        if found == 1:
            start_line = start_end[0] + 2
            # Keep going until end of file in case the sequence got broken up into multiple lines
            while linecache.getline(sequence_file, start_line):
                sequence += linecache.getline(sequence_file, start_line).strip()
                start_line += 1

        if found == 2:
            start_line = start_end[0]
            end_line = start_end[1]
            # Appending all lines to sequence after stripping
            for line in content[start_line + 1 : end_line]:
                sequence += line.strip()
    return sequence


def read_alignment(protein1: str, protein2: str, alignment: str) -> Tuple[str, str]:
    aligned_sequence1 = ""
    aligned_sequence2 = ""
    with open(alignment, "r") as f:
        content = f.readlines()
        for linenum, line in enumerate(content):
            if protein1.casefold() in line.casefold():
                aligned_sequence1 = content[linenum + 1].strip()
            if protein2.casefold() in line.casefold():
                aligned_sequence2 = content[linenum + 1].strip()
    return aligned_sequence1, aligned_sequence2


def align_attention(
    attention1: np.ndarray, attention2: np.ndarray, sequence1: str, sequence2: str
) -> Tuple[np.ndarray, np.ndarray, List[int], List[int]]:
    og_seq1_len = len(attention1)
    og_seq2_len = len(attention2)

    seq_ratio = og_seq1_len / og_seq2_len

    aligned_attention1 = np.zeros([len(sequence1)])
    aligned_attention2 = np.zeros([len(sequence2)])

    gaps1 = []
    gaps2 = []
    i1_vals = 0
    i2_vals = 0
    for (i1, resi1), (i2, resi2) in zip(enumerate(sequence1), enumerate(sequence2)):
        if resi1 == "-":
            aligned_attention1[i1] = 0
            gaps1.append(i1)
        elif resi1.isalpha():
            aligned_attention1[i1] = attention1[i1_vals]
            i1_vals += 1

        if resi2 == "-":
            aligned_attention2[i2] = 0
            gaps2.append(i2)
        elif resi2.isalpha():
            aligned_attention2[i2] = attention2[i2_vals]
            i2_vals += 1

    # Because I'm always doing seq 1/seq 2
    # If seq 1 is longer then the ratio will be > 1
    # Need to multiple seq 1 by ratio as its longer and needs more attention
    # If seq 1 is shorter the ratio is < 1
    # Multiple seq 1 by the ratio as its shorter and needs less attention
    aligned_attention1 = aligned_attention1 * seq_ratio
    return aligned_attention1, aligned_attention2, gaps1, gaps2
