import logging
import tarfile
from pathlib import Path

import requests
import tqdm

logger = logging.getLogger(__name__)


def download_alphafold_params(target=Path("params")):
    success_marker = target.joinpath("download_finished.txt")

    if success_marker.is_file():
        logger.info("")

    target.mkdir(exist_ok=True)
    url = "https://storage.googleapis.com/alphafold/alphafold_params_2021-07-14.tar"
    response = requests.get(url, stream=True)
    file_size = int(response.headers.get("Content-Length", 0))
    with tqdm.tqdm.wrapattr(
        response.raw,
        "read",
        total=file_size,
        desc=f"Downloading alphafold2 weights to {target}",
    ) as response_raw:
        # Open in stream mode ("r|"), as our requests response doesn't support seeking)
        file = tarfile.open(fileobj=response_raw, mode="r|")
        file.extractall(path=target)
    success_marker.touch()


if __name__ == "__main__":
    download_alphafold_params()
