from typing import Dict, Tuple

import haiku
from alphafold.model import model, config, data


def load_models_and_params(
    num_models: int,
) -> Dict[str, Tuple[model.RunModel, haiku.Params]]:
    # Use only two model and later swap params to avoid recompiling
    # note: models 1,2 have diff number of params compared to models 3,4,5
    model_runner_and_params: Dict[str, Tuple[model.RunModel, haiku.Params]] = dict()
    model_runner_1 = None
    model_runner_3 = None
    for model_number in range(1, num_models + 1):
        if model_number in [1, 2]:
            if not model_runner_1:
                model_config = config.model_config("model_1_ptm")
                model_config.data.eval.num_ensemble = 1
                model_runner_1 = model.RunModel(
                    model_config,
                    data.get_model_haiku_params(model_name="model_1_ptm", data_dir="."),
                )
            model_runner = model_runner_1
        else:
            assert model_number in [3, 4, 5]

            if not model_runner_3:
                model_config = config.model_config("model_3_ptm")
                model_config.data.eval.num_ensemble = 1
                model_runner_3 = model.RunModel(
                    model_config,
                    data.get_model_haiku_params(model_name="model_3_ptm", data_dir="."),
                )
            model_runner = model_runner_3

        model_name = f"model_{model_number}"
        params = data.get_model_haiku_params(
            model_name=model_name + "_ptm", data_dir="."
        )
        model_runner_and_params[model_name] = (model_runner, params)
    return model_runner_and_params
