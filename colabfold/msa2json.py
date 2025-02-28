"""
This code has been adapted from the original code in `cddlab/alphafold3_tools`

Find the original code at:
https://github.com/cddlab/alphafold3_tools/blob/main/alphafold3tools/msatojson.py
"""

import concurrent.futures
import json
import os
from argparse import ArgumentDefaultsHelpFormatter, ArgumentParser
from dataclasses import dataclass
from pathlib import Path

# from loguru import logger
import logging

# from alphafold3tools.log import log_setup


@dataclass
class Seq:
    name: str
    sequence: str


def get_residuelens_stoichiometries(lines) -> tuple[list[int], list[int]]:
    """Get residue lengths and stoichiometries from msa file.
    Args:
        lines: list[str]
            Lines of input msa file
    Returns:
        residue_lens: list[int]
            Residue lengths of each polypeptide chain
        stoichiometries: list[int]
            Stoichiomerties of each polypeptide chain
    """
    residue_lens_, stoichiometries_ = lines[0].split("\t")
    residue_lens = list(map(int, residue_lens_.lstrip("#").split(",")))
    stoichiometries = list(map(int, stoichiometries_.split(",")))
    return residue_lens, stoichiometries


def split_a3msequences(residue_lens, line) -> list[str]:
    """Split a3m sequences into a list of a3m sequences.
    Note: The a3m-format MSA file represents inserted residues with lowercase.
    The first line (starting with '#') of the MSA file contains residue lengths
    and stoichiometries of each polypeptide chain.
    From the second line, the first sequence is the query.
    After this, the paired MSA blocks are followed by the unpaired MSA.
    Args:
        residue_lens: list[int]
            Residue lengths of each polypeptide chain
        line: str
            A3M sequences
    Returns:
        a3msequences: list[str]
            A3M sequences, len(a3msequences) should be the same as len(residue_lens).
    """
    a3msequences = [""] * len(residue_lens)
    i = 0
    count = 0
    current_residue = []

    for char in line:
        current_residue.append(char)
        if char == "-" or char.isupper():
            count += 1
        if count == residue_lens[i]:
            a3msequences[i] = "".join(current_residue)
            current_residue = []
            count = 0
            i += 1
            if i == len(residue_lens):
                break

    if current_residue and i < len(residue_lens):
        a3msequences[i] = "".join(current_residue)

    return a3msequences


def get_paired_and_unpaired_msa(
    lines: list[str], residue_lens: list[int], cardinality: int
) -> tuple[list[list[Seq]], list[list[Seq]]]:
    """Get paired and unpaired MSAs from input MSA file.
    Args:
        lines: list[str]
            Lines of input MSA file
        residue_lens: list[int]
            Residue lengths of each polypeptide chain
        cardinality: int
            Number of polypeptide chains
        query_seqnames: list[int]
            Query sequence names
    Returns:
        pairedmsas: list[list[Seq]]
            Paired MSAs, len(pairedmsa) should be the cardinality.
            If cardinality is 1, pairedmsas returns [[Seq("", "")]].
        unpairedmsas: list[list[Seq]]
            Unpaired MSAs, len(unpairedmsa) should be the cardinality.
    """
    pairedmsas: list[list[Seq]] = [[] for _ in range(cardinality)]
    unpairedmsas: list[list[Seq]] = [[] for _ in range(cardinality)]
    pairedflag = False
    unpairedflag = False
    seen = False
    seqnames_seen = []
    query_seqnames = [int(101 + i) for i in range(cardinality)]
    chain = -1
    for line in lines[1:]:
        if line.startswith(">"):
            if line not in seqnames_seen:
                seqnames_seen.append(line)
            else:
                seen = True
                continue
            if cardinality > 1 and line.startswith(
                ">" + "\t".join(map(str, query_seqnames)) + "\n"
            ):
                pairedflag = True
                unpairedflag = False
            elif any(line.startswith(f">{seq}\n") for seq in query_seqnames):
                pairedflag = False
                unpairedflag = True
                chain += 1
            seqname = line
        else:
            if seen:
                seen = False
                continue
            if pairedflag:
                a3mseqs = split_a3msequences(residue_lens, line)
                for i in range(cardinality):
                    pairedmsas[i].append(Seq(seqname, a3mseqs[i]))

            elif unpairedflag:
                a3mseqs = split_a3msequences(residue_lens, line)
                for i in range(cardinality):
                    # Remove all-gapped sequences
                    if a3mseqs[i] == "-" * residue_lens[i]:
                        continue
                    unpairedmsas[i].append(Seq(seqname, a3mseqs[i]))
            else:
                raise ValueError("Flag must be either paired or unpaired.")
    return pairedmsas, unpairedmsas


def convert_msas_to_str(msas):
    """convert MSAs to str format for AlphaFold3 input JSON file."""
    if msas == []:
        return ""
    else:
        return "\n".join(f"{seq.name}{seq.sequence}" for seq in msas) + "\n"


def int_id_to_str_id(num: int) -> str:
    """Encodes a number as a string, using reverse spreadsheet style naming.
    This block is cited from
    https://github.com/google-deepmind/alphafold3/blob/main/src/alphafold3/structure/mmcif.py#L40

    Args:
      num: A positive integer.

    Returns:
      A string that encodes the positive integer using reverse spreadsheet style,
      naming e.g. 1 = A, 2 = B, ..., 27 = AA, 28 = BA, 29 = CA, ... This is the
      usual way to encode chain IDs in mmCIF files.
    """
    if num <= 0:
        raise ValueError(f"Only positive integers allowed, got {num}.")

    num = num - 1  # 1-based indexing.
    output = []
    while num >= 0:
        output.append(chr(num % 26 + ord("A")))
        num = num // 26 - 1
    return "".join(output)


def generate_input_json_content(
    name: str,
    cardinality: int,
    stoichiometries: list[int],
    pairedmsas: list[list[Seq]],
    unpairedmsas: list[list[Seq]],
) -> str:
    """generate AlphaFold3 input JSON file.

    Args:
        name (str): Name of the protein complex.
                    Used for the name field in the JSON file.
        cardinality (int): The number of distinct polypeptide chains.
        stoichiometries (list[int]): Stoichiometries of each polypeptide chain.
        pairedmsas (list[list[Seq]]): Paired MSAs.
        unpairedmsas (list[list[Seq]]): Unpaired MSAs.
    Returns:
        str: JSON string for AlphaFold3 input file.
    """
    sequences: list[dict] = []
    chain_id_count = 0
    null = None
    for i in range(cardinality):
        # unpairedmsa[i][0] is more appropriate than pairedmsa[i][0].
        query_seq = unpairedmsas[i][0].sequence
        chain_ids = [
            int_id_to_str_id(chain_id_count + j + 1) for j in range(stoichiometries[i])
        ]
        chain_id_count += stoichiometries[i]
        sequences.append(
            {
                "protein": {
                    "id": chain_ids,
                    "sequence": query_seq,
                    "modifications": [],
                    "unpairedMsa": convert_msas_to_str(unpairedmsas[i]),
                    "pairedMsa": convert_msas_to_str(pairedmsas[i]),
                    "templates": [],
                }
            }
        )
    content = json.dumps(
        {
            "dialect": "alphafold3",
            "version": 2,
            "name": f"{name}",
            "sequences": sequences,
            "modelSeeds": [1],
            "bondedAtomPairs": null,
            "userCCD": null,
        },
        indent=4,
    )
    return content


def write_input_json_file(
    inputmsafile: str | Path,
    name: str,
    outputjsonfile: str | Path,
) -> None:
    """Write AlphaFold3 input JSON file from a3m-format MSA file.

    Args:
        inputmsafile (str): Input MSA file path.
        name (str): Name of the protein complex.
                    Used for the name field in the JSON file.
        cardinality (int): The number of distinct polypeptide chains.
        stoichiometries (list[int]): Stoichiometries of each polypeptide chain.
        pairedmsas (list[list[Seq]]): Paired MSAs.
        unpairedmsas (list[list[Seq]]): Unpaired MSAs.
        outputfile (str): Output file path.
    """
    with open(inputmsafile, "r") as f:
        lines = f.readlines()
    residue_lens, stoichiometries = get_residuelens_stoichiometries(lines)
    if len(residue_lens) != len(stoichiometries):
        raise ValueError("Length of residue_lens and stoichiometries must be the same.")
    cardinality = len(residue_lens)
    # TODO: Change loguru logger to logging
    # logger.info(
    #     f"The input MSA file contains {cardinality} distinct polypeptide chains."
    # )
    # logger.info(f"Residue lengths: {residue_lens}")
    # logger.info(f"Stoichiometries: {stoichiometries}")
    pairedmsas, unpairedmsas = get_paired_and_unpaired_msa(
        lines, residue_lens, cardinality
    )
    content = generate_input_json_content(
        name=f"{name}",
        cardinality=cardinality,
        stoichiometries=stoichiometries,
        pairedmsas=pairedmsas,
        unpairedmsas=unpairedmsas,
    )
    with open(outputjsonfile, "w") as f:
        f.write(content)


def process_a3m_file(a3m_file, output_dir):
    name = Path(a3m_file).stem
    output_file = os.path.join(output_dir, f"{name}.json")
    write_input_json_file(a3m_file, name, output_file)


def main():
    parser = ArgumentParser(
        formatter_class=ArgumentDefaultsHelpFormatter,
        description="Converts a3m-format MSA file to AlphaFold3 input JSON file.",
    )
    parser.add_argument(
        "-i",
        "--input",
        help="Input A3M file or directory containing A3M files. e.g. 1bjp.a3m",
        type=str,
        required=True,
    )
    parser.add_argument(
        "-n",
        "--name",
        help="Name of the protein complex.",
        type=str,
        default="",
    )
    parser.add_argument(
        "-o", "--out", help="Output directory or JSON file.", type=str, required=True
    )
    parser.add_argument(
        "-d",
        "--debug",
        help="Print lots of debugging statements",
        dest="loglevel",
        action="store_const",
        const="DEBUG",
        default="SUCCESS",
    )
    args = parser.parse_args()
    # log_setup(args.loglevel)
    # Default name is the input file name without extension
    if args.name == "":
        args.name = os.path.splitext(os.path.basename(args.input))[0]
    # logger.info(f"Name of the protein complex: {args.name}") # TODO
    input_path = Path(args.input)
    if not input_path.exists():
        raise FileNotFoundError(f"{input_path} does not exist.")
    out_path = Path(args.out)
    if input_path.is_dir():
        # logger.info(f"Input directory: {input_path}") # TODO
        if out_path.suffix == ".json":
            raise ValueError(
                "Now the input is directory, so output name must be a directory."
            )
        # logger.info(f"Output directory: {out_path}") # TODO
        out_path.mkdir(parents=True, exist_ok=True)
        a3m_files = list(input_path.glob("*.a3m"))
        with concurrent.futures.ThreadPoolExecutor() as executor:
            futures = [
                executor.submit(process_a3m_file, a3m_file, out_path)
                for a3m_file in a3m_files
            ]
            concurrent.futures.wait(futures)
    else:
        name = input_path.stem
        if input_path.suffix != ".a3m":
            raise ValueError("Input file must have .a3m extension.")
        # logger.info(f"Input A3M file: {input_path}") # TODO
        if out_path.suffix != ".json":
            raise ValueError("Output file must have .json extension.")
        # logger.info(f"Output JSON file: {out_path}") # TODO
        write_input_json_file(args.input, name, out_path)


if __name__ == "__main__":
    main()