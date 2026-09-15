"""
AlphaFold3 backend: implements :class:`colabfold.backend.FoldingBackend`.
"""
import logging
from datetime import date
from pathlib import Path
from typing import Any, Dict, List, Optional

from colabfold.backend import RunOptions

logger = logging.getLogger(__name__)

_CONFIG_KEYS = ("num_diffusion_samples", "use_dropout", "buckets", "download_weights",
                "af3_pairing", "use_fast_kernels")

_DEFAULTS = {
    "num_diffusion_samples": 5,
    "use_dropout": False,
    "buckets": None,
    "download_weights": True,
    "af3_pairing": "colabfold",
    "use_fast_kernels": False,
}


class AF3Backend:
    """AlphaFold3 and the other models alphafold3-open's registry can run."""

    def __init__(self, model_type: str):
        from colabfold.alphafold3 import require

        require()
        self.model_type = model_type
        self.model_runner = None
        self._fold_input = None
        self._num_seeds = 1
        self._random_seed = 0
        self._use_templates = False
        self._max_template_hits = 20

    def _opt(self, opts: RunOptions, key: str):
        return opts.opt(key, _DEFAULTS[key])

    def configure(self, opts, *, max_len, max_num, num_queries, msa_mode,
                  is_complex, use_templates) -> None:
        self._num_seeds = opts.num_seeds
        self._random_seed = opts.random_seed
        self._use_templates = use_templates
        self._max_template_hits = opts.max_template_hits
        self._warn_about_ignored(opts)
        samples = self._opt(opts, "num_diffusion_samples")
        logger.info(f"{self.model_type}: {opts.num_seeds} seed(s) x {samples} diffusion "
                    f"sample(s) = {opts.num_seeds * samples} structure(s)")

    def _warn_about_ignored(self, opts: RunOptions) -> None:
        """Say which AlphaFold2 options this backend does not honour."""
        ignored = []
        if opts.num_relax:
            ignored.append(("--amber", "alphafold3 results are not relaxed"))
        if opts.initial_guess:
            ignored.append(("--initial-guess", "diffusion takes no starting structure"))
        if opts.max_extra_seq is not None:
            ignored.append(("--max-msa", "only the first number is used"))
        if opts.save_recycles:
            ignored.append(("--save-recycles", "there is no per-recycle output"))
        for flag, reason in ignored:
            logger.warning(f"{self.model_type} ignores {flag}: {reason}")

    def featurize(self, query_seqs_unique, query_seqs_cardinality, unpaired_msa, paired_msa,
                  template_results, is_complex: bool, opts: RunOptions, extras=None):
        import dataclasses

        from colabfold.alphafold3.input import build_fold_input, with_msas

        seeds = [self._random_seed + i for i in range(self._num_seeds)]
        pairing = self._opt(opts, "af3_pairing")
        templates = self._templates_for(query_seqs_unique, template_results)
        given = getattr(extras, "fold_input", None)
        if given is not None:
            fold_input = with_msas(given, unpaired_msa, paired_msa, pairing)
            fold_input = dataclasses.replace(fold_input, rng_seeds=seeds)
        else:
            fold_input = build_fold_input(
                name="colabfold",
                query_seqs_unique=query_seqs_unique,
                query_seqs_cardinality=query_seqs_cardinality,
                unpaired_msa=unpaired_msa,
                paired_msa=paired_msa,
                molecules=extras if isinstance(extras, (list, tuple)) else None,
                seeds=seeds,
                pairing=pairing,
                templates=templates,
            )
        self._fold_input = fold_input
        return fold_input, {}

    def _templates_for(self, query_seqs_unique, template_results):
        """One alphafold3 Template list per unique sequence, from ColabFold's hits."""
        if not self._use_templates or not template_results:
            return None
        from colabfold.alphafold.features import search_templates
        from colabfold.alphafold3.templates import build_templates

        out = []
        for sequence, result in zip(query_seqs_unique, template_results):
            if result is None:
                out.append([])
                continue
            a3m_lines, template_path = result
            try:
                hits = search_templates(a3m_lines, template_path)
                out.append(build_templates(hits, sequence, template_path,
                                           max_templates=self._max_template_hits))
            except Exception as e:
                logger.warning(f"no templates for this chain: {e}")
                out.append([])
        logger.info(f"templates per chain: {[len(t) for t in out]}")
        return out

    def _ensure_loaded(self, opts: RunOptions) -> None:
        if self.model_runner is not None:
            return
        from colabfold.alphafold3.models import load_model

        if self._opt(opts, "use_fast_kernels"):
            from colabfold.alphafold3.attention import install

            install()
        model_dir = opts.opt("model_dir")
        self.model_runner = load_model(
            self.model_type,
            num_recycles=opts.num_recycles,
            num_diffusion_samples=self._opt(opts, "num_diffusion_samples"),
            model_dir=Path(model_dir) if model_dir else None,
            use_dropout=self._opt(opts, "use_dropout"),
            download=self._opt(opts, "download_weights"),
            num_msa=opts.max_seq,
            return_embeddings=opts.save_single_representations or opts.save_pair_representations,
        )

    def predict(self, prefix: str, result_dir: Path, model_input, is_complex: bool,
                sequences_lengths: List[int], opts: RunOptions,
                prediction_callback=None) -> Dict[str, Any]:
        from colabfold.alphafold3.models import featurise
        from colabfold.alphafold3.predict import predict_structure

        self._ensure_loaded(opts)
        import dataclasses

        fold_input = dataclasses.replace(model_input, name=prefix)
        examples = featurise(fold_input, self.model_runner.model_name,
                             self.model_runner.model_dir, buckets=self._opt(opts, "buckets"),
                             ref_max_modified_date=date.fromisoformat(opts.max_template_date))
        return predict_structure(
            prefix=prefix,
            result_dir=result_dir,
            fold_input=fold_input,
            model_runner=self.model_runner,
            model_type=self.model_type,
            featurised_examples=examples,
            save_all=opts.save_all,
            rank_by=opts.rank_by,
            stop_at_score=opts.stop_at_score,
            prediction_callback=prediction_callback,
        )

    def config_dict(self, opts: RunOptions) -> Dict[str, Any]:
        cfg = {k: self._opt(opts, k) for k in _CONFIG_KEYS}
        if opts.num_recycles is None:
            cfg["num_recycles"] = 10
        return cfg

    def plot_msa(self, model_input, dpi: int = 200):
        from colabfold.alphafold3.plot import plot_msa_coverage

        return plot_msa_coverage(model_input, dpi=dpi)

    def plot_extra_metrics(self, scores, fig_path) -> None:
        return None
