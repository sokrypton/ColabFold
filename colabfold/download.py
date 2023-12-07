import time
import tarfile
from pathlib import Path
from multiprocessing import Queue, Process

import appdirs
import tqdm
import requests

# The data dir location logic switches between a version with and one without "params" because alphafold
# always internally joins "params". (We should probably patch alphafold)
default_data_dir = Path(appdirs.user_cache_dir(__package__ or "colabfold"))

def download(url, params_dir, size_queue, progress_queue):
    try:
        filename = url.split("/")[-1]
        response = requests.get(url, stream=True, timeout=6.02)
        response.raise_for_status()

        file_size = int(response.headers.get('Content-Length', 0))
        size_queue.put(file_size)

        output = params_dir.joinpath(filename)
        if output.is_file() and output.stat().st_size == file_size:
            progress_queue.put(file_size)
            return

        with open(output, "wb") as output:
            for chunk in response.iter_content(chunk_size=8192):
                if chunk:
                    output.write(chunk)
                    progress_queue.put(len(chunk))
    except Exception as e:
        size_queue.put("error")
        progress_queue.put("error")


def download_alphafold_params(model_type: str, data_dir: Path = default_data_dir):
    params_dir = data_dir.joinpath("params")
    if model_type == "alphafold2_multimer_v3":
        url = "https://storage.googleapis.com/alphafold/alphafold_params_colab_2022-12-06.tar"
        success_marker = params_dir.joinpath(
            "download_complexes_multimer_v3_finished.txt"
        )
    elif model_type == "alphafold2_multimer_v2":
        url = "https://storage.googleapis.com/alphafold/alphafold_params_colab_2022-03-02.tar"
        success_marker = params_dir.joinpath(
            "download_complexes_multimer_v2_finished.txt"
        )
    elif model_type == "alphafold2_multimer_v1":
        url = "https://storage.googleapis.com/alphafold/alphafold_params_colab_2021-10-27.tar"
        success_marker = params_dir.joinpath(
            "download_complexes_multimer_v1_finished.txt"
        )
    elif model_type == "AlphaFold2-ptm" or model_type == "alphafold2_ptm" or model_type == "alphafold2":
        url = "https://storage.googleapis.com/alphafold/alphafold_params_2021-07-14.tar"
        success_marker = params_dir.joinpath("download_finished.txt")
    elif model_type == "deepfold_v1":
        url = [
            "https://colabfold.steineggerlab.workers.dev/deepfold_v1/deepfold_model_1.npz",
            "https://colabfold.steineggerlab.workers.dev/deepfold_v1/deepfold_model_2.npz",
            "https://colabfold.steineggerlab.workers.dev/deepfold_v1/deepfold_model_3.npz",
            "https://colabfold.steineggerlab.workers.dev/deepfold_v1/deepfold_model_4.npz",
            "https://colabfold.steineggerlab.workers.dev/deepfold_v1/deepfold_model_5.npz",
        ]
        success_marker = params_dir.joinpath("download_deepfold-v1_finished.txt")
    else:
        raise ValueError(f"Unknown model {model_type}")

    if success_marker.is_file():
        return

    params_dir.mkdir(parents=True, exist_ok=True)
    if isinstance(url, str):
        response = requests.get(url, stream=True)
        file_size = int(response.headers.get("Content-Length", 0))
        with tqdm.tqdm.wrapattr(
            response.raw,
            "read",
            total=file_size,
            desc=f"Downloading {model_type} weights to {data_dir}",
        ) as response_raw:
            # Open in stream mode ("r|"), as our requests response doesn't support seeking)
            file = tarfile.open(fileobj=response_raw, mode="r|")
            file.extractall(path=params_dir)

    elif isinstance(url, list):
        # download files in parallel to increase download speed
        size_queue = Queue()
        progress_queue = Queue()
        processes = []

        for u in url:
            process = Process(target=download, args=(u, params_dir, size_queue, progress_queue))
            process.start()
            processes.append(process)

        total_size = 0
        sizes_collected = 0
        error = 0
        while sizes_collected < len(url):
            size = size_queue.get()
            if size == "error":
                error = 1
                break
            total_size += size
            sizes_collected += 1

        with tqdm.tqdm(
                total=total_size,
                desc=f"Downloading {model_type} weights to {data_dir}",
                unit="B",
                unit_scale=True,
                unit_divisor=1024
            ) as pbar:
            downloaded = 0
            while downloaded < total_size:
                progress = progress_queue.get()
                if progress == "error":
                    error = 1
                    break
                downloaded += progress
                pbar.update(downloaded - pbar.n)

        for process in processes:
            process.join()

        if error:
            raise RuntimeError("Error downloading files")

    success_marker.touch()


if __name__ == "__main__":
    from sys import argv
    if len(argv) == 2:
        download_alphafold_params(argv[1])
    else:
        download_alphafold_params("alphafold2_multimer_v3")
        download_alphafold_params("AlphaFold2-ptm")
