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
                    "unpairedMsa": null,
                    "pairedMsa": null,
                    "templates": null,
                }
            }
        )
    content = json.dumps(
        {
            "dialect": "alphafold3",
            "version": 1,
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
    inputfastafile: str | Path,
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
    with open(inputfastafile, "r") as f:
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
    content = generate_input_json_content(
        name=f"{name}",
        cardinality=cardinality,
        stoichiometries=stoichiometries,
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
        description="Generates AlphaFold3 default input JSON file from FASTA file.",
    )
    parser.add_argument(
        "-i",
        "--input",
        help="Input fasta file. e.g. 1bjp.fasta",
        type=str,
        required=True,
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
    # logger.info(f"Name of the protein complex: {args.name}") # TODO
    input_path = Path(args.input)
    if not input_path.exists():
        raise FileNotFoundError(f"{input_path} does not exist.")
    out_path = Path(args.out)

    model_dict = {} # TODO
    with open(args.input, "r") as f:
        lines = f.readlines()


    # if input_path.is_dir():
    #     # logger.info(f"Input directory: {input_path}") # TODO
    #     if out_path.suffix == ".json":
    #         raise ValueError(
    #             "Now the input is directory, so output name must be a directory."
    #         )
    #     # logger.info(f"Output directory: {out_path}") # TODO
    #     out_path.mkdir(parents=True, exist_ok=True)
    #     a3m_files = list(input_path.glob("*.a3m"))
    #     with concurrent.futures.ThreadPoolExecutor() as executor:
    #         futures = [
    #             executor.submit(process_a3m_file, a3m_file, out_path)
    #             for a3m_file in a3m_files
    #         ]
    #         concurrent.futures.wait(futures)
    # else:
    #     name = input_path.stem
    #     if input_path.suffix != ".a3m":
    #         raise ValueError("Input file must have .a3m extension.")
    #     # logger.info(f"Input A3M file: {input_path}") # TODO
    #     if out_path.suffix != ".json":
    #         raise ValueError("Output file must have .json extension.")
    #     # logger.info(f"Output JSON file: {out_path}") # TODO
    #     write_input_json_file(args.input, name, out_path)


if __name__ == "__main__":
    main()