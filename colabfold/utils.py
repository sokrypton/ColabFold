import logging
import warnings
from pathlib import Path

from absl import logging as absl_logging
from tqdm import TqdmExperimentalWarning

NO_GPU_FOUND = """ERROR: Jax could not find GPU. This can be either because your machine doesn't have a GPU
or because jax can't find it. You might need to run

pip install --upgrade "jax[cuda111]" -f https://storage.googleapis.com/jax-releases/jax_releases.html  # Note: wheels only available on linux.

See https://github.com/google/jax/#pip-installation-gpu-cuda for more details.

If you're sure you want to run without a GPU, pass `--cpu`"""

DEFAULT_API_SERVER = "https://a3m.mmseqs.com"

ACCEPT_DEFAULT_TERMS = """WARNING: You are welcome to use the default MSA server, however keep in mind that it's a limited shared resource only capable of processing a few thousand MSAs per day. Please submit jobs only from a single IP address. We reserve the right to limit access to the server case-by-case when usage exceeds fair use.

If you require more MSAs, please host your own API and pass it to `--host-url`"""


class TqdmHandler(logging.StreamHandler):
    """https://stackoverflow.com/a/38895482/3549270"""

    def __init__(self):
        logging.StreamHandler.__init__(self)

    def emit(self, record):
        # We need the native tqdm here
        from tqdm import tqdm

        msg = self.format(record)
        tqdm.write(msg)


def setup_logging(log_file: Path):
    log_file.parent.mkdir(exist_ok=True, parents=True)
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(message)s",
        handlers=[TqdmHandler(), logging.FileHandler(log_file)],
    )
    # otherwise jax will tell us about its search for devices
    absl_logging.set_verbosity("error")
    warnings.simplefilter(action="ignore", category=TqdmExperimentalWarning)


def safe_filename(file: str) -> str:
    return "".join([c if c.isalnum() or c in ["_", ".", "-"] else "_" for c in file])
