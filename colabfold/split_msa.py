from argparse import ArgumentParser
from pathlib import Path

from tqdm import tqdm


def main():
    parser = ArgumentParser(
        description="Take an a3m database from the colabdb search and turn it into a folder of a3m files"
    )
    parser.add_argument("mmseqs_msa")
    parser.add_argument("output_folder")
    args = parser.parse_args()
    output_folder = Path(args.output_folder)
    output_folder.mkdir(exist_ok=True)
    for msa in tqdm(Path(args.mmseqs_msa).read_text().split("\0")):
        filename = msa.split("\n", 1)[0][1:].split(" ")[0] + ".a3m"
        output_folder.joinpath(filename).write_text(msa)


if __name__ == "__main__":
    main()
