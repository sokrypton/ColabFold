import logging

from absl import logging as absl_logging

NO_GPU_FOUND = """ERROR: Jax could not find GPU. This can be either because your machine doesn't have a GPU
or because jax can't find it. You might need to run

pip install --upgrade "jax[cuda111]" -f https://storage.googleapis.com/jax-releases/jax_releases.html  # Note: wheels only available on linux.

See https://github.com/google/jax/#pip-installation-gpu-cuda for more details.

If you're sure you want to run without a GPU, pass `--cpu`"""

class TqdmHandler(logging.StreamHandler):
    """https://stackoverflow.com/a/38895482/3549270"""

    def __init__(self):
        logging.StreamHandler.__init__(self)

    def emit(self, record):
        # We need the native tqdm here
        from tqdm import tqdm

        msg = self.format(record)
        tqdm.write(msg)


def setup_logging():
    logging.basicConfig(
        level=logging.INFO, format="%(asctime)s %(message)s", handlers=[TqdmHandler()]
    )
    # otherwise jax will tell us about its search for devices
    absl_logging.set_verbosity("error")


