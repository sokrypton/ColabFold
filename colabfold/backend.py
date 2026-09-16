"""
Folding backend interface.

``colabfold.batch`` drives prediction through a :class:`FoldingBackend` chosen
once per run by :func:`get_backend` from the ``model_type`` string.

- ``RunOptions`` carries the run knobs plus an opaque ``backend_opts`` dict for
  backend-private settings (e.g. use_pallas / compile_mode).
- ``backend.predict(...)`` returns ``{"rank", "metric", "result_files"}``:
  ranked tags, scalar score dicts, and file paths.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional, Protocol, runtime_checkable


# the metrics a prediction may report, in the order they are printed
_METRIC_LABELS = (("mean_plddt", "pLDDT"), ("ptm", "pTM"), ("iptm", "ipTM"),
                  ("actifptm", "actifpTM"), ("ranking_score", "ranking_score"))


def metrics_line(scores) -> str:
    """The ` pLDDT=95.4 pTM=0.826` tail both backends put after a tag."""
    return "".join(f" {label}={float(scores[key]):.3g}"
                   for key, label in _METRIC_LABELS if key in scores)


@dataclass
class RunOptions:
    model_type: str
    num_models: int = 1
    num_seeds: int = 1
    num_recycles: Optional[int] = None
    recycle_early_stop_tolerance: Optional[float] = None
    use_templates: bool = False
    max_template_date: str = "2100-01-01"
    max_template_hits: int = 20
    rank_by: str = "auto"
    stop_at_score: float = 100.0
    random_seed: int = 0
    initial_guess: Optional[str] = None
    # relaxation
    num_relax: int = 0
    relax_max_iterations: int = 0
    relax_tolerance: float = 2.39
    relax_stiffness: float = 10.0
    relax_max_outer_iterations: int = 3
    use_gpu_relax: bool = False
    # outputs
    save_all: bool = False
    save_single_representations: bool = False
    save_pair_representations: bool = False
    save_recycles: bool = False
    # resources / sizing
    data_dir: Path = Path(".")
    max_seq: Optional[int] = None
    max_extra_seq: Optional[int] = None
    # opaque, backend-private
    backend_opts: Dict[str, Any] = field(default_factory=dict)

    def opt(self, key: str, default: Any = None) -> Any:
        return self.backend_opts.get(key, default)


@runtime_checkable
class FoldingBackend(Protocol):
    def configure(self, opts: RunOptions, *, max_len, max_num, num_queries,
                  msa_mode, is_complex, use_templates) -> None:
        ...

    def featurize(
        self,
        query_seqs_unique: List[str],
        query_seqs_cardinality: List[int],
        unpaired_msa,
        paired_msa,
        template_results,
        is_complex: bool,
        opts: RunOptions,
        extras=None,
    ):
        ...

    def predict(
        self,
        prefix: str,
        result_dir: Path,
        model_input: Dict[str, Any],
        is_complex: bool,
        sequences_lengths: List[int],
        opts: RunOptions,
        prediction_callback=None,
    ) -> Dict[str, Any]:
        ...

    def config_dict(self, opts: RunOptions) -> Dict[str, Any]:
        ...

    def plot_msa(self, model_input, dpi: int = 200):
        ...

    def plot_extra_metrics(self, scores, fig_path) -> None:
        ...


_current_backend = None

AF3_MODELS = (
    "alphafold3", "af3", "openfold3", "of3", "openbind", "openbind0",
    "protenix", "protenix1", "protenix2", "boltz2", "chai1", "chai",
    "intellifold2", "if2", "intellifold", "opendde", "rosettafold3", "rf3",
    "esmfold2", "esmfold2_fast", "esmfold2_lm300m", "esmfold2_lm600m",
)


def is_af3_model(model_type: str) -> bool:
    """True for anything alphafold3-open's registry runs on its AF3 graph."""
    if model_type.startswith("alphafold2") or model_type.startswith("deepfold"):
        return False
    return model_type.startswith("alphafold3") or model_type in AF3_MODELS


def get_backend(model_type: str, data_dir=None) -> FoldingBackend:
    global _current_backend
    if is_af3_model(model_type):
        from colabfold.alphafold3.backend import AF3Backend
        _current_backend = AF3Backend(model_type, data_dir)
        return _current_backend
    if model_type.startswith("alphafold2") or model_type.startswith("deepfold"):
        from colabfold.alphafold.backend import AF2Backend
        _current_backend = AF2Backend(model_type)
        return _current_backend
    raise NotImplementedError(
        f"model_type {model_type!r} is not implemented yet"
    )


def get_current_backend() -> FoldingBackend | None:
    return _current_backend
