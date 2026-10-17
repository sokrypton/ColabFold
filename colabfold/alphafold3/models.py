"""
Load an alphafold3-open model and run its forward pass.

alphafold3-open ships run_alphafold.py outside the wheel, so the few steps
ColabFold needs from it are here instead.
"""
import functools
import logging
import re
import string
from pathlib import Path
from typing import Any, List, Optional, Sequence

import numpy as np

logger = logging.getLogger(__name__)


def mark_absl_flags_parsed() -> None:
    """Stop alphafold3 reading ColabFold's argv as absl flags."""
    from absl import flags

    flags.FLAGS([""])  # absl wants the program name, and nothing after it


def resolve_model_name(model_type: str) -> str:
    """ColabFold's ``--model-type`` to an alphafold3-open registry name."""
    from alphafold3.model import model_registry

    name = model_type[len("alphafold3_"):] if model_type.startswith("alphafold3_") else model_type
    model_registry.get(name)
    return model_registry.ALIASES.get(name, name)


def known_models() -> List[str]:
    from alphafold3.model import model_registry

    return sorted(set(model_registry.MODEL_SPECS) | set(model_registry.ALIASES))


def ensure_weights(model_name: str, download: bool = True, data_dir: Optional[Path] = None,
                   precision: str = "fp32", accept_terms: bool = False) -> Path:
    from colabfold.alphafold3.weights import ensure_weights as fetch

    return fetch(model_name, data_dir=data_dir, download=download,
                 precision=precision, accept_terms=accept_terms)


def default_num_recycles(model_name: str) -> int:
    """The model's own default, in AF's convention (N recycles = N + 1 trunk passes)."""
    from alphafold3.model import model_registry

    return getattr(model_registry.get(model_name), "default_num_recycles", 10)


def make_config(model_name: str, num_recycles: Optional[int], num_diffusion_samples: int,
                num_msa: Optional[int] = None, return_embeddings: bool = False,
                return_logits: bool = False):
    from alphafold3.model import model, model_registry

    config = model.Model.Config()
    config.heads.diffusion.eval.num_samples = num_diffusion_samples
    config.num_recycles = num_recycles if num_recycles is not None else default_num_recycles(model_name)
    if num_msa is not None:
        config.evoformer.num_msa = num_msa
    config.return_embeddings = return_embeddings
    config.return_confidence_logits = config.return_distogram = return_logits
    model_registry.get(model_name).configure(config)
    from colabfold_kernels import compute_capability

    cc = compute_capability()
    if cc is not None and cc < 80:
        # the sm_70 and sm_75 kernels are float16, and those cards have no bfloat16 units
        config.global_config.half_dtype = "float16"
        # a T4 samples slower in half precision than in float32, and less precisely
        config.global_config.bfloat16 = "intermediate"
    elif cc is not None and config.global_config.bfloat16 != "none":
        # flash triangle attention holds no [rows, heads, N, N] logits, so row chunks only add a scan
        config.global_config.pair_attention_chunk_size = ((2048, None), (None, 256))
    return config


class ModelRunner:
    """The parts of alphafold3-open's run_alphafold.ModelRunner that a fold needs."""

    def __init__(self, config, model_dir: Path, use_dropout: bool = False, device=None,
                 fused_layer_norm: bool = False):
        self._config = config
        self._fused_layer_norm = fused_layer_norm
        self._model_dir = Path(model_dir)
        self._use_dropout = use_dropout
        self._device = device
        self._apply = None

    @property
    def model_dir(self) -> Path:
        return self._model_dir

    @property
    def model_name(self) -> str:
        return self._config.global_config.model

    @functools.cached_property
    def model_params(self):
        from alphafold3.model import params

        return params.get_model_haiku_params(model_dir=str(self._model_dir))

    def _model(self):
        if self._apply is None:
            import haiku as hk
            import jax
            from alphafold3.model import model

            @hk.transform
            def forward_fn(batch):
                if self._fused_layer_norm:
                    from colabfold.alphafold3.layer_norm import interceptor

                    with hk.intercept_methods(interceptor):
                        return model.Model(self._config)(batch, use_dropout=self._use_dropout)
                return model.Model(self._config)(batch, use_dropout=self._use_dropout)

            apply_fn = jax.jit(forward_fn.apply, device=self._device)
            self._preinit_tokamax_context()
            self._apply = functools.partial(apply_fn, self.model_params)
        return self._apply

    @staticmethod
    def _preinit_tokamax_context() -> None:
        # tokamax builds its autotune context on first use; if that happens inside
        # the first trace the whole model is retraced on the second call
        try:
            from tokamax._src.ops import op as tokamax_op

            tokamax_op.get_autotuning_cache_overlay_state()
        except Exception:
            pass

    def run_inference(self, featurised_example, rng_key):
        import jax
        import jax.numpy as jnp
        from alphafold3.model.components import utils

        featurised_example = jax.device_put(
            jax.tree_util.tree_map(jnp.asarray, utils.remove_invalidly_typed_feats(featurised_example)),
            self._device,
        )
        result = self._model()(rng_key, featurised_example)
        result = jax.tree.map(np.asarray, result)
        result = jax.tree.map(
            lambda x: x.astype(jnp.float32) if x.dtype == jnp.bfloat16 else x, result)
        result = dict(result)
        # alphafold3 writes this into the output mmCIF and refuses to build a result
        # without it; converted weights carry one, a hand-made blob may not
        meta = self.model_params.get("__meta__", {}).get("__identifier__")
        result["__identifier__"] = (np.asarray(meta).tobytes() if meta is not None
                                    else self.model_name.encode())
        return result

    def extract_embeddings(self, result, num_tokens: int):
        out = {}
        for key in ("single_embeddings", "pair_embeddings"):
            if key in result:
                value = result[key]
                out[key] = (value[:num_tokens] if key == "single_embeddings"
                            else value[:num_tokens, :num_tokens]).astype(np.float16)
        return out or None

    def extract_logits(self, result, sample: int, num_tokens: int):
        """One sample's confidence logits and its seed's distogram logits, trimmed to the tokens"""
        out = {}
        for key in ("predicted_lddt_logits", "experimentally_resolved_logits"):
            if key in result:
                out[key] = result[key][sample, :num_tokens]
        for key in ("predicted_aligned_error_logits", "predicted_distance_error_logits"):
            if key in result:
                out[key] = result[key][sample, :num_tokens, :num_tokens]
        if "distogram" in result["distogram"]:
            out["distogram_logits"] = result["distogram"]["distogram"][:num_tokens, :num_tokens]
        return {key: value.astype(np.float16) for key, value in out.items()} or None

    def extract_inference_results(self, batch, result, target_name: str):
        from alphafold3.model import model

        return list(model.Model.get_inference_result(batch=batch, result=result,
                                                     target_name=target_name))


def _esm_dir(name: str, model_dir: Path) -> str:
    """A language model lives beside the models, and is fetched on first use."""
    return str(Path(model_dir).expanduser().parent / name)


def _refuse_esm_only_model(model_name: str) -> None:
    """A variant with no MSA encoder has nothing to fold from without its language model."""
    from alphafold3.model import model_registry

    variant = model_registry.ESMFOLD2_VARIANTS.get(model_name)
    if variant is not None and not variant.get("msa"):
        raise ValueError(f"{model_name} does not use MSAs and needs --use-esm")


def _resolve_esm(use_esm: bool, fold_input, model_name: str, model_dir: Path):
    """``use_esm`` -> ``(esm2 rows, (esmc rows or None, shim params))``, either of which may be None."""
    if not use_esm:
        _refuse_esm_only_model(model_name)
        return None, None
    from alphafold3.model import esm, model_registry

    sequences = [chain.sequence for chain in fold_input.chains
                 if type(chain).__name__ == "ProteinChain"]
    if not sequences and model_name not in model_registry.ESMFOLD2_VARIANTS:
        logger.warning(f"--use-esm needs a protein chain, and {fold_input.name} has none")
        return None, None
    if model_name == "chai1":
        logger.info("embedding the sequences with ESM2 for chai1")
        rows = esm.embed(sequences, _esm_dir("esm2", model_dir), "esm2")
        return rows, None
    if model_name not in model_registry.ESMFOLD2_VARIANTS:
        logger.warning(f"--use-esm does nothing for {model_name}, which folds from its MSA")
        return None, None
    # each release trains its own shim; another's reads as noise
    shim = esm.load_shim_params(str(model_dir), model_name)
    if not sequences:
        # native runs the shim over every token all the same, on zero states off protein
        logger.info(f"{fold_input.name} has no protein chain: {model_name}'s language-model pair from zero states")
        return None, (None, shim)
    variant = model_registry.ESMFOLD2_VARIANTS[model_name]["esmc"]
    logger.info(f"embedding {len(sequences)} protein chain(s) with {variant} for {model_name}")
    return None, (esm.embed(sequences, _esm_dir(variant, model_dir), "esmc", variant), shim)


_A3M_ROW = re.compile(r"^[^>\s].*", re.M)
_INSERTIONS = str.maketrans("", "", string.ascii_lowercase)


def _aligned(a3m: Optional[str]) -> bool:
    """True if an a3m holds a row other than its query; a query-only MSA is no alignment at all."""
    rows = (m.group().rstrip() for m in _A3M_ROW.finditer(a3m or ""))
    query = next(rows, None)
    return any(row.translate(_INSERTIONS) != query for row in rows)


def featurise(fold_input, model_name: str, model_dir: Path, buckets: Optional[Sequence[int]] = None,
              use_esm: bool = False, ref_max_modified_date=None) -> List[Any]:
    """Featurise one fold input, applying the model family's own conventions."""
    from alphafold3.constants import decoded_ccd
    from alphafold3.data import featurisation
    from alphafold3.model import model_registry
    from alphafold3.model.pipeline import model_features

    spec = model_registry.get(model_name)
    fold_input = model_features.prepare_input(fold_input, spec)
    ccd = decoded_ccd.get_ccd(user_ccd=fold_input.user_ccd)
    build = functools.partial(
        featurisation.featurise_input,
        fold_input=fold_input,
        buckets=buckets,
        ccd=ccd,
        # a CCD entry without ideal coordinates falls back to its model ones, and
        # alphafold3 raises on the date comparison that decides it if this is None
        ref_max_modified_date=ref_max_modified_date,
        flatten_non_standard_residues=not spec.featurise.get("modified_as_one_token", False),
        ccd_leaving_atoms=bool(spec.featurise.get("ccd_leaving_atoms", False)),
    )
    examples = build(verbose=False)
    if spec.featurise:
        has_msa = any(_aligned(getattr(c, "unpaired_msa", None)) or _aligned(getattr(c, "paired_msa", None))
                      for c in fold_input.chains)
        esm_rows, lm_hidden = _resolve_esm(use_esm, fold_input, model_name, model_dir)
        examples = [
            model_features.apply(example, spec, refeaturise=lambda: build(verbose=False),
                                 model_dir=str(model_dir), esm=esm_rows, lm_hidden=lm_hidden,
                                 has_msa=has_msa, fold_input=fold_input, cyclic=False, seed=seed)
            for seed, example in zip(fold_input.rng_seeds, examples)
        ]
    return examples


def load_model(model_type: str, *, num_recycles: Optional[int], num_diffusion_samples: int,
               use_dropout: bool = False, download: bool = True, num_msa: Optional[int] = None,
               return_embeddings: bool = False, return_logits: bool = False,
               data_dir: Optional[Path] = None, precision: str = "fp32",
               accept_terms: bool = False, fused_layer_norm: bool = False) -> "ModelRunner":
    mark_absl_flags_parsed()
    model_name = resolve_model_name(model_type)
    weights_dir = ensure_weights(model_name, download=download, data_dir=data_dir,
                                 precision=precision, accept_terms=accept_terms)
    config = make_config(model_name, num_recycles, num_diffusion_samples,
                         num_msa=num_msa, return_embeddings=return_embeddings,
                         return_logits=return_logits)
    logger.info(f"Running {model_name} from {weights_dir}")
    return ModelRunner(config, weights_dir, use_dropout=use_dropout,
                       fused_layer_norm=fused_layer_norm)
