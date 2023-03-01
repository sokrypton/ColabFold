import json
import logging
import warnings
from pathlib import Path
from typing import Any, Callable, Dict, List, Optional, Tuple, Union, TYPE_CHECKING

from absl import logging as absl_logging
from importlib_metadata import distribution
from tqdm import TqdmExperimentalWarning


NO_GPU_FOUND = \
"""
WARNING: Jax could not find GPU. This can be either because your machine doesn't have a 
GPU or because jax can't find it. For instructions how to fix this, see: 
https://github.com/YoshitakaMo/localcolabfold
"""
DEFAULT_API_SERVER = "https://api.colabfold.com"
ACCEPT_DEFAULT_TERMS = \
"""
WARNING: You are welcome to use the default MSA server, however keep in mind that it's a
limited shared resource only capable of processing a few thousand MSAs per day. Please
submit jobs only from a single IP address. We reserve the right to limit access to the
server case-by-case when usage exceeds fair use. If you require more MSAs: You can 
precompute all MSAs with `colabfold_search` or host your own API and pass it to `--host-url`
"""

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
  root = logging.getLogger()
  if root.handlers:
    for handler in root.handlers:
      root.removeHandler(handler)
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

def get_commit() -> Optional[str]:
  text = distribution("colabfold").read_text("direct_url.json")
  if not text:
    return None
  direct_url = json.loads(text)
  if "vcs_info" not in direct_url:
    return None
  if "commit_id" not in direct_url["vcs_info"]:
    return None
  return direct_url["vcs_info"]["commit_id"]

class file_manager:
  def __init__(self, prefix: str, result_dir: Path):
    self.prefix = prefix
    self.result_dir = result_dir
    self.tag = None
    self.files = {}
  
  def get(self, x: str, ext:str) -> Path:
    if self.tag not in self.files:
      self.files[self.tag] = []
    file = self.result_dir.joinpath(f"{self.prefix}_{x}_{self.tag}.{ext}")
    self.files[self.tag].append([x,ext,file])
    return file

  def set_tag(self, tag):
    self.tag = tag