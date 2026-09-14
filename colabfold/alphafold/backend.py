"""
AlphaFold2 backend: implements :class:`FoldingBackend`.
"""
import logging
import math
from pathlib import Path
from typing import Any, Dict, List

from colabfold.backend import RunOptions

logger = logging.getLogger(__name__)

_CONFIG_KEYS = (
    "model_order",
    "num_ensemble",
    "use_dropout",
    "use_cluster_profile",
    "use_fuse",
    "use_bfloat16",
    "use_fast_kernels",
    "kernel_backend",
    "compile_mode",
    "recompile_padding",
    "calc_extra_ptm",
    "use_probs_extra",
)

# Defaults for backend_opts parameters
_DEFAULTS = {
    "model_order": [1, 2, 3, 4, 5],
    "num_ensemble": 1,
    "use_dropout": False,
    "use_cluster_profile": True,
    "use_fuse": True,
    "use_bfloat16": True,
    "use_fast_kernels": False,
    "kernel_backend": "auto",
    "compile_mode": "tuned",
    "recompile_padding": 10,
    "calc_extra_ptm": False,
    "use_probs_extra": True,
}


class AF2Backend:
    """AlphaFold2 / AlphaFold2-multimer backend."""

    def __init__(self, model_type: str):
        self.model_type = model_type
        self.model_runner_and_params = None
        # MSA-cluster sizing
        self.max_seq = None
        self.max_extra_seq = None
        # recompilation-avoidance state
        self._max_len = None
        self._pad_len = 0
        self._msa_pad_depth = 0
        self._num_queries = None
        self._msa_mode = None
        self._use_templates = False

    def _opt(self, opts: RunOptions, key: str):
        return opts.opt(key, _DEFAULTS[key])

    def configure(self, opts, *, max_len, max_num, num_queries, msa_mode,
                  is_complex, use_templates) -> None:
        self._max_len = max_len
        self._num_queries = num_queries
        self._msa_mode = msa_mode
        self._use_templates = use_templates

        # MSA cluster sizes per model variant:
        #   512 5120 = alphafold2_ptm (models 1,3,4) / 512 1024 (models 2,5)
        #   508 2048 = alphafold2_multimer_v3 (models 1,2,3) / 508 1152 (models 4,5)
        #   252 1152 = alphafold2_multimer_v[1,2]
        set_if = lambda x, y: y if x is None else x
        max_seq, max_extra_seq = opts.max_seq, opts.max_extra_seq
        if self.model_type in ("alphafold2_multimer_v1", "alphafold2_multimer_v2"):
            max_seq, max_extra_seq = set_if(max_seq, 252), set_if(max_extra_seq, 1152)
        elif self.model_type == "alphafold2_multimer_v3":
            max_seq, max_extra_seq = set_if(max_seq, 508), set_if(max_extra_seq, 2048)
        else:
            max_seq, max_extra_seq = set_if(max_seq, 512), set_if(max_extra_seq, 5120)

        if msa_mode == "single_sequence":
            num_seqs = 1
            if is_complex and "multimer" not in self.model_type:
                num_seqs += max_num
            if use_templates:
                num_seqs += 4
            max_seq = min(num_seqs, max_seq)
            max_extra_seq = max(min(num_seqs - max_seq, max_extra_seq), 1)

        self.max_seq = max_seq
        self.max_extra_seq = max_extra_seq

    def _ensure_loaded(self, opts: RunOptions, model_input) -> None:
        if self.model_runner_and_params is not None:
            return
        # For a single query with a real MSA, shrink max_seq to the actual depth.
        if self._num_queries == 1 and self._msa_mode != "single_sequence":
            if "msa_mask" in model_input:
                num_seqs = int(sum(model_input["msa_mask"].max(-1) == 1))
            else:
                num_seqs = int(len(model_input["msa"]))
            if self._use_templates:
                num_seqs += 4
            self.max_seq = min(num_seqs, self.max_seq)
            self.max_extra_seq = max(min(num_seqs - self.max_seq, self.max_extra_seq), 1)
            logger.info(f"Setting max_seq={self.max_seq}, max_extra_seq={self.max_extra_seq}")
        self.load(opts)

    def load(self, opts: RunOptions) -> None:
        from colabfold.alphafold.models import load_models_and_params

        self.model_runner_and_params = load_models_and_params(
            num_models=opts.num_models,
            use_templates=opts.use_templates,
            num_recycles=opts.num_recycles,
            num_ensemble=self._opt(opts, "num_ensemble"),
            model_order=self._opt(opts, "model_order"),
            model_type=self.model_type,
            data_dir=opts.data_dir,
            stop_at_score=opts.stop_at_score,
            rank_by=opts.rank_by,
            use_dropout=self._opt(opts, "use_dropout"),
            max_seq=self.max_seq,
            max_extra_seq=self.max_extra_seq,
            use_cluster_profile=self._opt(opts, "use_cluster_profile"),
            recycle_early_stop_tolerance=opts.recycle_early_stop_tolerance,
            use_fuse=self._opt(opts, "use_fuse"),
            use_bfloat16=self._opt(opts, "use_bfloat16"),
            save_all=opts.save_all,
            calc_extra_ptm=self._opt(opts, "calc_extra_ptm"),
            use_probs_extra=self._opt(opts, "use_probs_extra"),
            use_fast_kernels=self._opt(opts, "use_fast_kernels"),
            kernel_backend=self._opt(opts, "kernel_backend"),
            compile_mode=self._opt(opts, "compile_mode"),
        )

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
        from colabfold.alphafold.features import (
            build_template_features,
            generate_input_feature,
        )

        template_features = build_template_features(
            query_seqs_unique,
            template_results,
            max_template_date=opts.max_template_date,
            max_template_hits=opts.max_template_hits,
        )
        return generate_input_feature(
            query_seqs_unique,
            query_seqs_cardinality,
            unpaired_msa,
            paired_msa,
            template_features,
            is_complex,
            self.model_type,
            max_seq=self.max_seq,
        )

    def plot_msa(self, model_input, dpi: int = 200):
        from colabfold.plot import plot_msa_v2

        return plot_msa_v2(model_input, dpi=dpi)

    def plot_extra_metrics(self, scores, fig_path):
        from colabfold.alphafold import extra_ptm
        extra_ptm.plot_chain_pairwise_analysis(scores, fig_path=fig_path)

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
        from colabfold.alphafold.predict import predict_structure

        # Load lazily on the first prediction (sizing depends on the first MSA).
        self._ensure_loaded(opts, model_input)

        # Grow the padded length (cap at the longest query) so inputs of similar
        # length share one compiled shape instead of recompiling per length.
        seq_len = sum(sequences_lengths)
        recompile_padding = self._opt(opts, "recompile_padding")
        if seq_len > self._pad_len:
            if isinstance(recompile_padding, float):
                self._pad_len = math.ceil(seq_len * recompile_padding)
            else:
                self._pad_len = seq_len + recompile_padding
            self._pad_len = min(self._pad_len, self._max_len)

        # Track the deepest multimer MSA seen so far so all queries share a shape.
        if "multimer" in self.model_type and "msa" in model_input:
            self._msa_pad_depth = max(self._msa_pad_depth, len(model_input["msa"]))

        return predict_structure(
            prefix=prefix,
            result_dir=result_dir,
            feature_dict=model_input,
            is_complex=is_complex,
            use_templates=opts.use_templates,
            sequences_lengths=sequences_lengths,
            pad_len=self._pad_len,
            model_type=self.model_type,
            model_runner_and_params=self.model_runner_and_params,
            initial_guess=opts.initial_guess,
            msa_pad_depth=self._msa_pad_depth,
            num_relax=opts.num_relax,
            relax_max_iterations=opts.relax_max_iterations,
            relax_tolerance=opts.relax_tolerance,
            relax_stiffness=opts.relax_stiffness,
            relax_max_outer_iterations=opts.relax_max_outer_iterations,
            rank_by=opts.rank_by,
            stop_at_score=opts.stop_at_score,
            prediction_callback=prediction_callback,
            use_gpu_relax=opts.use_gpu_relax,
            random_seed=opts.random_seed,
            num_seeds=opts.num_seeds,
            save_all=opts.save_all,
            save_single_representations=opts.save_single_representations,
            save_pair_representations=opts.save_pair_representations,
            save_recycles=opts.save_recycles,
            calc_extra_ptm=self._opt(opts, "calc_extra_ptm"),
            use_probs_extra=self._opt(opts, "use_probs_extra"),
        )

    def config_dict(self, opts: RunOptions) -> Dict[str, Any]:
        cfg = {k: self._opt(opts, k) for k in _CONFIG_KEYS}
        cfg["max_seq"] = self.max_seq
        cfg["max_extra_seq"] = self.max_extra_seq
        return cfg
