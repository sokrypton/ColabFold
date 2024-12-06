from pathlib import Path
from functools import wraps, partialmethod
from typing import Tuple, List, Optional
import haiku
from alphafold.model import model, config, data
from alphafold.model.modules import AlphaFold
from alphafold.model.modules_multimer import AlphaFold as AlphaFoldMultimer

def get_model_haiku_params(
    data_dir: str,
    model_type: str,
    model_number: str,
    use_fuse: bool = True,
    to_jnp: bool = True,
) -> haiku.Params:
    """Get the Haiku parameters from a model type and number."""
    import os
    import numpy as np
    from alphafold.model import utils
    import jax.numpy as jnp

    is_deepfold = False
    if model_type == "alphafold2_multimer_v1":
        file = f"params_model_{model_number}_multimer.npz"
    elif model_type == "alphafold2_multimer_v2":
        file = f"params_model_{model_number}_multimer_v2.npz"
    elif model_type == "alphafold2_multimer_v3":
        file = f"params_model_{model_number}_multimer_v3.npz"
    elif model_type == "alphafold2_ptm":
        file = f"params_model_{model_number}_ptm.npz"
    elif model_type == "alphafold2":
        file = f"params_model_{model_number}.npz"
    elif model_type == "deepfold_v1":
        file = f"deepfold_model_{model_number}.npz"
        is_deepfold = True
    else:
        raise ValueError(f"Unknown model_type {model_type}")

    path = os.path.join(data_dir, "params", file)
    params = np.load(path, allow_pickle=False)
    return utils.flat_params_to_haiku(params, fuse=use_fuse, to_jnp=to_jnp)


def model_to_config_name(model_type: str, model_number: str) -> str:
    if model_type == "alphafold2_multimer_v1":
        return f"model_{model_number}_multimer"
    elif model_type == "alphafold2_multimer_v2":
        return f"model_{model_number}_multimer_v2"
    elif model_type == "alphafold2_multimer_v3":
        return f"model_{model_number}_multimer_v3"
    elif model_type == "alphafold2_ptm":
        return f"model_{model_number}_ptm"
    elif model_type == "alphafold2":
        return f"model_{model_number}"
    elif model_type == "deepfold_v1":
        return f"model_{model_number}"
    else:
        raise ValueError(f"Unknown model_type {model_type}")


def load_models_and_params(
    num_models: int,
    use_templates: bool,
    num_recycles: Optional[int] = None,
    recycle_early_stop_tolerance: Optional[float] = None,
    num_ensemble: int = 1,
    model_order: Optional[List[int]] = None,
    model_type: str = "",
    data_dir: Path = Path("."),
    stop_at_score: float = 100,
    rank_by: str = "auto",
    max_seq: Optional[int] = None,
    max_extra_seq: Optional[int] = None,
    use_cluster_profile: bool = True,
    use_fuse: bool = True,
    use_bfloat16: bool = True,
    use_dropout: bool = False,
    save_all: bool = False,
    calc_extra_ptm: bool = False,
    use_probs_extra: bool = True
) -> List[Tuple[str, model.RunModel, haiku.Params]]:
    """We use only two actual models and swap the parameters to avoid recompiling.

    Note that models 1 and 2 have a different number of parameters compared to models 3, 4 and 5,
    so we load model 1 and model 3.
    """

    # Use only two model and later swap params to avoid recompiling
    model_runner_and_params: [Tuple[str, model.RunModel, haiku.Params]] = []

    if model_order is None:
        model_order = [1, 2, 3, 4, 5]
    else:
        model_order.sort()

    model_build_order = [3, 4, 5, 1, 2]
    if "multimer" in model_type:
        models_need_compilation = [3]
    else:
        # only models 1,2 use templates
        models_need_compilation = [1, 3] if use_templates else [3]
    
    model_runner_and_params_build_order: [Tuple[str, model.RunModel, haiku.Params]] = []
    model_runner = None
    for model_number in model_build_order:
        if model_number in models_need_compilation:
            # get configurations
            config_name = model_to_config_name(model_type, model_number)
            model_config = config.model_config(config_name)
            model_config.model.stop_at_score = float(stop_at_score)
            model_config.model.rank_by = rank_by

            # set dropouts
            model_config.model.global_config.eval_dropout = use_dropout

            # set bfloat options
            model_config.model.global_config.bfloat16 = use_bfloat16
            
            # set fuse options
            model_config.model.embeddings_and_evoformer.evoformer.triangle_multiplication_incoming.fuse_projection_weights = use_fuse
            model_config.model.embeddings_and_evoformer.evoformer.triangle_multiplication_outgoing.fuse_projection_weights = use_fuse
            if "multimer" in config_name or model_number in [1,2]:
                model_config.model.embeddings_and_evoformer.template.template_pair_stack.triangle_multiplication_incoming.fuse_projection_weights = use_fuse
                model_config.model.embeddings_and_evoformer.template.template_pair_stack.triangle_multiplication_outgoing.fuse_projection_weights = use_fuse
                        
            # set number of sequences options
            if max_seq is not None:
                if "multimer" in config_name:
                    model_config.model.embeddings_and_evoformer.num_msa = max_seq
                else:
                    model_config.data.eval.max_msa_clusters = max_seq
            
            if max_extra_seq is not None:
                if "multimer" in config_name:
                    model_config.model.embeddings_and_evoformer.num_extra_msa = max_extra_seq
                else:
                    model_config.data.common.max_extra_msa = max_extra_seq

            # disable some outputs if not being saved
            if not save_all:
                if not calc_extra_ptm:
                    model_config.model.heads.distogram.weight = 0.0
                model_config.model.heads.masked_msa.weight = 0.0
                model_config.model.heads.experimentally_resolved.weight = 0.0

            # set number of recycles and ensembles            
            if "multimer" in config_name:
                if num_recycles is not None:
                    model_config.model.num_recycle = num_recycles
                model_config.model.embeddings_and_evoformer.use_cluster_profile = use_cluster_profile
                model_config.model.num_ensemble_eval = num_ensemble
            else:
                if num_recycles is not None:
                    model_config.data.common.num_recycle = num_recycles
                    model_config.model.num_recycle = num_recycles
                model_config.data.eval.num_ensemble = num_ensemble

            if recycle_early_stop_tolerance is not None:
                model_config.model.recycle_early_stop_tolerance = recycle_early_stop_tolerance
            
            # get model runner
            params = get_model_haiku_params(
                model_type=model_type,
                model_number=model_number,
                data_dir=str(data_dir),
                use_fuse=use_fuse
            )

            model_runner = model.RunModel(
                model_config,
                params,
                extended_ptm_config={'calc_extended_ptm': calc_extra_ptm,
                                     'use_probs_extended': use_probs_extra}
            )
        
        params = get_model_haiku_params(
            model_type=model_type,
            model_number=model_number,
            data_dir=str(data_dir),
            use_fuse=use_fuse,
        )
        # keep only parameters of compiled model
        params_subset = {}
        for k in model_runner.params.keys():
            params_subset[k] = params[k]

        model_name = f"model_{model_number}"
        model_runner_and_params_build_order.append(
            (model_name, model_runner, params_subset)
        )
    # reorder model
    for n, model_number in enumerate(model_order):
        if n == num_models:
            break
        model_name = f"model_{model_number}"
        for m in model_runner_and_params_build_order:
            if model_name == m[0]:
                model_runner_and_params.append(m)
                break
    return model_runner_and_params
