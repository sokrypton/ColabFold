from pathlib import Path
from functools import wraps, partialmethod
from typing import Tuple, List, Optional

import haiku

from alphafold.model import model, config, data
from alphafold.model.modules import AlphaFold
from alphafold.model.modules_multimer import AlphaFold as AlphaFoldMultimer


def load_models_and_params(
    num_models: int,
    use_templates: bool,
    num_recycle: int = 3,
    num_ensemble: int = 1,
    model_order: Optional[List[int]] = None,
    model_suffix: str = "_ptm",
    data_dir: Path = Path("."),
    recompile_all_models: bool = False,
    stop_at_score: float = 100,
    rank_by: str = "plddt",
    return_representations: bool = False,
    training: bool = False,
    max_msa: str = None,
) -> List[Tuple[str, model.RunModel, haiku.Params]]:
    """We use only two actual models and swap the parameters to avoid recompiling.

    Note that models 1 and 2 have a different number of parameters compared to models 3, 4 and 5,
    so we load model 1 and model 3.
    """

    if return_representations:
        # this forces the AlphaFold to always return representations
        AlphaFold.__call__ = partialmethod(
            AlphaFold.__call__, return_representations=True
        )

        AlphaFoldMultimer.__call__ = partialmethod(
            AlphaFoldMultimer.__call__, return_representations=True
        )

    if not model_order:
        model_order = [3, 4, 5, 1, 2]

    # Use only two model and later swap params to avoid recompiling
    model_runner_and_params: [Tuple[str, model.RunModel, haiku.Params]] = []

    if recompile_all_models:
        for n, model_number in enumerate(model_order):
            if n == num_models:
                break
            model_name = f"model_{model_number}"
            params = data.get_model_haiku_params(
                model_name=model_name + model_suffix, data_dir=str(data_dir)
            )
            model_config = config.model_config(model_name + model_suffix)
            model_config.model.stop_at_score = float(stop_at_score)
            model_config.model.stop_at_score_ranker = rank_by
            if max_msa != None:
                max_msa_clusters, max_extra_msa = [int(x) for x in max_msa.split(":")]
                model_config.data.eval.max_msa_clusters = max_msa_clusters
                model_config.data.common.max_extra_msa = max_extra_msa
            if model_suffix == "_ptm":
                model_config.data.common.num_recycle = num_recycle
                model_config.model.num_recycle = num_recycle
                model_config.data.eval.num_ensemble = num_ensemble
            elif model_suffix.startswith("_multimer"):
                model_config.model.num_recycle = num_recycle
                if training:
                    model_config.model.num_ensemble_train = num_ensemble
                else:
                    model_config.model.num_ensemble_eval = num_ensemble
            model_runner_and_params.append(
                (
                    model_name,
                    model.RunModel(model_config, params, is_training=training),
                    params,
                )
            )
    else:
        models_need_compilation = [1, 3] if use_templates else [3]
        model_build_order = [3, 4, 5, 1, 2]
        model_runner_and_params_build_order: [
            Tuple[str, model.RunModel, haiku.Params]
        ] = []
        model_runner = None
        for model_number in model_build_order:
            if model_number in models_need_compilation:
                model_config = config.model_config(
                    "model_" + str(model_number) + model_suffix
                )
                model_config.model.stop_at_score = float(stop_at_score)
                model_config.model.stop_at_score_ranker = rank_by
                if max_msa != None:
                    max_msa_clusters, max_extra_msa = [
                        int(x) for x in max_msa.split(":")
                    ]
                    model_config.data.eval.max_msa_clusters = max_msa_clusters
                    model_config.data.common.max_extra_msa = max_extra_msa
                if model_suffix == "_ptm":
                    model_config.data.common.num_recycle = num_recycle
                    model_config.model.num_recycle = num_recycle
                    model_config.data.eval.num_ensemble = num_ensemble
                elif model_suffix.startswith("_multimer"):
                    model_config.model.num_recycle = num_recycle
                    if training:
                        model_config.model.num_ensemble_train = num_ensemble
                    else:
                        model_config.model.num_ensemble_eval = num_ensemble
                model_runner = model.RunModel(
                    model_config,
                    data.get_model_haiku_params(
                        model_name="model_" + str(model_number) + model_suffix,
                        data_dir=str(data_dir),
                    ),
                    is_training=training,
                )
            model_name = f"model_{model_number}"
            params = data.get_model_haiku_params(
                model_name=model_name + model_suffix, data_dir=str(data_dir)
            )
            # keep only parameters of compiled model
            params_subset = {}
            for k in model_runner.params.keys():
                params_subset[k] = params[k]

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
