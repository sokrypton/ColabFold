"""
Load an alphafold3-open model and run its forward pass.

alphafold3-open ships run_alphafold.py outside the wheel, so the few steps
ColabFold needs from it are here instead.
"""
import functools
import logging
from pathlib import Path
from typing import Any, List, Optional, Sequence

import numpy as np

logger = logging.getLogger(__name__)


def resolve_model_name(model_type: str) -> str:
    """ColabFold's ``--model-type`` to an alphafold3-open registry name."""
    from alphafold3.model import model_registry

    name = model_type[len("alphafold3_"):] if model_type.startswith("alphafold3_") else model_type
    model_registry.get(name)
    return model_registry.ALIASES.get(name, name)


def known_models() -> List[str]:
    from alphafold3.model import model_registry

    return sorted(set(model_registry.MODEL_SPECS) | set(model_registry.ALIASES))


def ensure_weights(model_name: str, model_dir: Optional[Path] = None, download: bool = True) -> Path:
    from alphafold3.model import weights

    return Path(weights.ensure_weights(model_name, model_dir=model_dir, download=download,
                                       log=logger.info))


def make_config(model_name: str, num_recycles: Optional[int], num_diffusion_samples: int):
    from alphafold3.model import model, model_registry

    config = model.Model.Config()
    config.heads.diffusion.eval.num_samples = num_diffusion_samples
    if num_recycles is not None:
        config.num_recycles = num_recycles
    model_registry.get(model_name).configure(config)
    return config


class ModelRunner:
    """The parts of alphafold3-open's run_alphafold.ModelRunner that a fold needs."""

    def __init__(self, config, model_dir: Path, use_dropout: bool = False, device=None):
        self._config = config
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
        from alphafold3.model import features  # noqa: F401
        from alphafold3.jax import utils

        featurised_example = jax.device_put(
            jax.tree_util.tree_map(jnp.asarray, utils.remove_invalidly_typed_feats(featurised_example)),
            self._device,
        )
        result = self._model()(rng_key, featurised_example)
        result = jax.tree.map(np.asarray, result)
        return result

    def extract_inference_results(self, batch, result, target_name: str):
        from alphafold3.model import model

        return list(model.Model.get_inference_result(batch=batch, result=result,
                                                     target_name=target_name))


def featurise(fold_input, model_name: str, model_dir: Path, buckets: Optional[Sequence[int]] = None,
              use_esm: bool = False, ref_max_modified_date=None) -> List[Any]:
    """Featurise one fold input, applying the model family's own conventions."""
    from alphafold3.constants import decoded_ccd
    from alphafold3.data import featurisation
    from alphafold3.model import model_registry
    from alphafold3.model.pipeline import model_features

    spec = model_registry.get(model_name)
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
    )
    examples = build(verbose=False)
    if spec.featurise:
        has_msa = any(getattr(c, "unpaired_msa", None) or getattr(c, "paired_msa", None)
                      for c in fold_input.chains)
        examples = [
            model_features.apply(example, spec, refeaturise=lambda: build(verbose=False),
                                 model_dir=str(model_dir), esm=None, lm_pair=None,
                                 has_msa=has_msa, fold_input=fold_input, cyclic=False)
            for example in examples
        ]
    return examples


def load_model(model_type: str, *, num_recycles: Optional[int], num_diffusion_samples: int,
               model_dir: Optional[Path] = None, use_dropout: bool = False,
               download: bool = True) -> "ModelRunner":
    model_name = resolve_model_name(model_type)
    weights_dir = ensure_weights(model_name, model_dir, download=download)
    config = make_config(model_name, num_recycles, num_diffusion_samples)
    logger.info(f"Running {model_name} from {weights_dir}")
    return ModelRunner(config, weights_dir, use_dropout=use_dropout)
