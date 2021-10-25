"""
colabdfold_search produces two a3m files with null separated msa in them.
We merge the two searches and then split into one a3m file per msa.
"""
import logging
from argparse import ArgumentParser
from pathlib import Path
from subprocess import check_call

from tqdm import tqdm

logger = logging.getLogger(__name__)


def merge_msa(mmseqs: str = "mmseqs", cwd: Path = Path(".")):
    check_call(
        [
            mmseqs,
            "mergedbs",
            "bfd.mgnify30.metaeuk30.smag30.a3m",
            "merged.a3m",
            "uniref.a3m",
            "bfd.mgnify30.metaeuk30.smag30.a3m",
        ],
        cwd=cwd,
    )
    return Path(cwd).joinpath("merged.a3m")


def split_msa(merged_msa: Path, output_folder: Path):
    for msa in tqdm(merged_msa.read_text().split("\0")):
        if not msa.strip():
            continue
        filename = msa.split("\n", 1)[0][1:].split(" ")[0] + ".a3m"
        output_folder.joinpath(filename).write_text(msa)


def main():
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s")

    parser = ArgumentParser(
        description="Take an a3m database from the colabdb search and turn it into a folder of a3m files"
    )
    parser.add_argument(
        "search_folder",
        help="The search folder in which you ran colabfold_search, "
        "which should contain uniref.a3m and bfd.mgnify30.metaeuk30.smag30.a3m",
    )
    parser.add_argument("output_folder", help="Will contain all the a3m files")
    parser.add_argument("--mmseqs", help="Path to the mmseqs2 binary", default="mmseqs")
    args = parser.parse_args()
    output_folder = Path(args.output_folder)
    output_folder.mkdir(exist_ok=True)

    logger.info("Merging MSAs")
    merged_msa = merge_msa(args.mmseqs, Path(args.search_folder))
    logger.info("Splitting MSAs")
    split_msa(merged_msa, output_folder)
    logger.info("Done")


if __name__ == "__main__":
    main()
