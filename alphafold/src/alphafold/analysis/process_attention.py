import os
import sys
import jax
import h5py
import logging
import linecache
import collections
import numpy as np

from typing import List, Tuple

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
logger = logging.getLogger(__name__)


def get_n(file_path: str) -> int:
    """Determine the sequence length `n` from an HDF5 archive.

    Peeks at all datasets within the HDF5 file to record shapes and returns
    the first dimension of the most frequent shape (the inferred n).

    Args:
        file_path: Path to the .h5 attention archive.

    Returns:
        int: the inferred `n`.
    """
    shape_counts = collections.Counter()

    try:
        with h5py.File(file_path, "r") as f:

            def visitor_func(name, obj):
                if isinstance(obj, h5py.Dataset):
                    shape_counts[obj.shape] += 1

            f.visititems(visitor_func)
    except Exception as e:
        logger.error("Could not open HDF5 file %s: %s", file_path, e)
        sys.exit(1)

    if not shape_counts:
        logger.error("No datasets found inside %s", file_path)
        sys.exit(1)

    sorted_shapes = shape_counts.most_common()

    logger.info("Unique shapes found in HDF5 (sorted by frequency):")
    for shape, count in sorted_shapes:
        logger.info("%s: %d datasets", shape, count)

    most_frequent_shape = sorted_shapes[0][0]

    return most_frequent_shape[0]


def get_attention(file_path: str, n: int) -> np.ndarray:
    """Load and convert attention data from an HDF5 archive into a per-head spectrum.

    The function performs the following steps:
      - Traverses the HDF5 hierarchy to locate all datasets.
      - Sorts datasets by their 'global_index' attribute to maintain model sequence.
      - Views the raw 2-byte dataset elements as float16.
      - Applies softmax (via jax.nn.softmax) to the restored values.
      - Filters for tensors matching the expected shape (n, 4, n, n).
      - Collapses axes 0, 1, and 2 by summation to produce a length-n vector per head.

    Args:
        file_path: Path to the .h5 attention archive generated during inference.
        n: Expected sequence length (inferred by get_n) used for shape validation.

    Returns:
        np.ndarray: A 2D array of shape (num_heads, n), where each row is the
            per-residue attention vector derived from a single attention head.
            This represents the 'attention spectrum' for the entire model run.
    """
    heads_n_4_n_n = []

    with h5py.File(file_path, "r") as f:
        index_map = []

        def find_indices(name, obj):
            if isinstance(obj, h5py.Dataset):
                idx = obj.attrs.get("global_index")
                if idx is not None:
                    index_map.append((idx, name))

        f.visititems(find_indices)

        index_map.sort(key=lambda x: x[0])

        for _, name in index_map:
            try:
                arr1 = f[name][:]

                arr1 = arr1.view(dtype=np.float16)
                arr1 = jax.nn.softmax(arr1)

                if arr1.shape == (n, 4, n, n):
                    heads_n_4_n_n.append(arr1)

            except Exception as e:
                logger.warning("Error processing head %s: %s", name, e)

    attention_spectrum = []
    for arr in heads_n_4_n_n:
        attn = np.sum(arr, axis=(0, 1, 2))
        attention_spectrum.append(attn)

    return np.array(attention_spectrum)


def average(attention_spectrum: np.ndarray) -> np.ndarray:
    """Compute the mean attention vector across an attention spectrum.

    Args:
        attention_spectrum: 2D array where each row is a per-file per-residue vector.

    Returns:
        np.ndarray: 1D array (length n) containing the mean attention per residue.
    """
    return np.mean(attention_spectrum, axis=0)


def min_max(data: np.ndarray) -> np.ndarray:
    """Scale a numeric array to the [0, 1] range using min-max normalization.

    If the array is constant (max == min) returns a zero array of the same shape.

    Args:
        data: 1D numeric array.

    Returns:
        np.ndarray: normalized array with values in [0, 1], or zeros if constant.
    """
    min_max = np.zeros([len(data)])
    data2 = [x for x in data if x != 0]
    min = np.min(data2)
    max = data.max()

    if min == max:
        for i, value in enumerate(data):
            if value != 0:
                min_max[i] = 1.0
        return min_max

    for i, value in enumerate(data):
        if value == 0:
            min_max[i] = 0
        else:
            min_max[i] = (value - min) / (max - min)
    return min_max


def read_sequence_file(sequence_file: str) -> str:
    """Read a FASTA-like sequence from a file.

    The function handles files with one or multiple FASTA records. If a single
    record is present it will read all lines after the header until EOF. If
    multiple records are present it reads the sequence lines between the first
    and second header lines.

    Args:
        sequence_file: path to the sequence file.

    Returns:
        str: concatenated sequence (header lines removed, whitespace stripped).
    """
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
    """Extract aligned sequences for two protein identifiers from an alignment file.

    The alignment file is expected to contain protein identifiers on a line
    followed by the aligned sequence on the next line. Matching is case-insensitive.

    Args:
        protein1: identifier for the first protein to search in the file.
        protein2: identifier for the second protein to search in the file.
        alignment: path to the alignment file.

    Returns:
        Tuple[str, str]: (aligned_sequence1, aligned_sequence2). Empty strings are
        returned for sequences not found.
    """
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
    """Map per-residue attention arrays onto aligned sequences (introducing gaps).

    For each aligned sequence position:
      - if the residue is '-', the aligned attention at that position is set to 0
        and the index is recorded in the corresponding gaps list;
      - if the residue is alphabetic, the next value from the original attention
        array is consumed and assigned.

    The function also rescales the aligned attention for sequence1 by the ratio
    of original sequence lengths to partially account for length differences.

    Args:
        attention1: 1D attention array for sequence1 (length = original seq1 length).
        attention2: 1D attention array for sequence2 (length = original seq2 length).
        sequence1: aligned sequence string for seq1 (may contain '-').
        sequence2: aligned sequence string for seq2 (may contain '-').

    Returns:
        Tuple containing:
          - aligned_attention1: 1D array mapped to sequence1 (length = len(sequence1)),
          - aligned_attention2: 1D array mapped to sequence2 (length = len(sequence2)),
          - gaps1: list of indices in sequence1 that are gaps,
          - gaps2: list of indices in sequence2 that are gaps.
    """
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

    # If seq 1 is longer then the ratio will be > 1
    # Multiply seq 1 by ratio as its longer and needs more attention
    # If seq 1 is shorter the ratio is < 1
    # Multiply seq 1 by the ratio as its shorter and needs less attention
    aligned_attention1 = aligned_attention1 * seq_ratio
    return aligned_attention1, aligned_attention2, gaps1, gaps2
