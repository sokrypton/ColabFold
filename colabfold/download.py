import logging
import tarfile
from pathlib import Path

import appdirs
import requests
import tqdm

logger = logging.getLogger(__name__)

# The data dir location logic switches between a version with and one without "params" because alphafold
# always internally joins "params". (We should probably patch alphafold)
default_data_dir = Path(appdirs.user_cache_dir(__package__ or "colabfold"))


def download_alphafold_params(data_dir: Path = default_data_dir):
    params_dir = data_dir.joinpath("params")
    success_marker = params_dir.joinpath("download_finished.txt")
    if success_marker.is_file():
        return

    params_dir.mkdir(parents=True, exist_ok=True)
    url = "https://storage.googleapis.com/alphafold/alphafold_params_2021-07-14.tar"
    response = requests.get(url, stream=True)
    file_size = int(response.headers.get("Content-Length", 0))
    with tqdm.tqdm.wrapattr(
        response.raw,
        "read",
        total=file_size,
        desc=f"Downloading alphafold2 weights to {data_dir}",
    ) as response_raw:
        # Open in stream mode ("r|"), as our requests response doesn't support seeking)
        file = tarfile.open(fileobj=response_raw, mode="r|")
        file.extractall(path=params_dir)
    success_marker.touch()


if __name__ == "__main__":
    download_alphafold_params()
