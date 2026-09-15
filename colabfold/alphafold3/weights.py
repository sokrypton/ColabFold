"""
Fetch alphafold3-open's converted weights through ColabFold's own download path.

They come from Hugging Face, where alphafold3-open itself fetches them. Once the
open data bucket carries them, put MIRROR back in front in ``urls_for``.
"""
import logging
from pathlib import Path
from typing import List, Optional, Tuple

from colabfold.download import default_data_dir, fetch_file

logger = logging.getLogger(__name__)

MIRROR = "https://opendata.mmseqs.org/colabfold/af3/models"
HUGGINGFACE = "https://huggingface.co/{repo}/resolve/main/{path}"

# Google DeepMind publishes AlphaFold 3's own parameters. They are not ours to
# redistribute, so they are neither mirrored nor fetched without --accept-alphafold3-terms.
OFFICIAL_URLS = {"alphafold3": "https://storage.googleapis.com/alphafold3/af3.bin.zst"}
TERMS = "https://github.com/google-deepmind/alphafold3/blob/main/WEIGHTS_TERMS_OF_USE.md"


def urls_for(path: str, repo: str) -> List[str]:
    """Where a blob is fetched from, tried in order."""
    return [
        # the bucket is not filled yet, and an empty mirror is a 404 per file
        # f"{MIRROR}/{path}",
        HUGGINGFACE.format(repo=repo, path=path),
    ]


def model_dir_for(model_name: str, data_dir: Path, precision: str = "fp32") -> Path:
    """Where ``--data`` keeps a model's weights, beside AlphaFold2's params."""
    suffix = "" if precision == "fp32" else f"-{precision}"
    return Path(data_dir).joinpath("params", "af3", model_name + suffix)


def _wanted(spec, precision: str) -> List[Tuple[str, str]]:
    """(path in the repo, local filename) for the blob and anything beside it."""
    from alphafold3.model import model_config

    out = [(spec.weights_path_for(precision), spec.weights_file_for(precision))]
    if spec.name in model_config.ESMFOLD2_FAMILY:
        # the LM shim, without which the ESM-C path raises on the first fold
        companion = f"{spec.name}.lm.npz"
        out.append((spec.companion_path(companion), companion))
    return out


def ensure_weights(model_name: str, data_dir: Optional[Path] = None,
                   model_dir: Optional[Path] = None, download: bool = True,
                   precision: str = "fp32", accept_terms: bool = False) -> Path:
    """Make ``model_name``'s weights exist on disk; return their directory."""
    from alphafold3.model import model_registry, weights

    spec = model_registry.get(model_name)
    official = OFFICIAL_URLS.get(spec.name)
    if official is not None and precision != "fp32":
        logger.info(f"{spec.name} is published as float32 only, ignoring --weights-precision {precision}")
        precision = "fp32"
    target = (Path(model_dir).expanduser() if model_dir is not None
              else model_dir_for(spec.name, data_dir or default_data_dir, precision))
    success_marker = target.joinpath(f"download_{spec.name}_{precision}_finished.txt")
    if success_marker.is_file():
        return target

    try:
        # weights already in place, converted by hand or pointed at with --model-dir
        return Path(weights.ensure_weights(model_name, model_dir=target,
                                           download=False, precision=precision))
    except FileNotFoundError:
        if not download or (spec.weights_repo is None and official is None):
            raise
    if official is not None and not accept_terms:
        raise RuntimeError(
            f"{spec.name} weights are Google DeepMind's and carry their own terms, which "
            f"do not allow commercial use. Read {TERMS} and pass --accept-alphafold3-terms "
            f"to download them, or point --model-dir at a copy you already have."
        )

    target.mkdir(parents=True, exist_ok=True)
    for path, filename in _wanted(spec, precision):
        dest = target.joinpath(filename)
        if not dest.is_file():
            urls = [official] if official is not None else urls_for(path, spec.weights_repo)
            fetch_file(urls, dest, f"Downloading {spec.name} weights to {target}")
    success_marker.touch()
    return target
