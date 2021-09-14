import logging
import math
import os
import pickle
import warnings
from argparse import ArgumentParser
from pathlib import Path
from string import ascii_uppercase
from typing import Any, Dict, Tuple, List, Union, Mapping

import haiku
import numpy
import numpy as np
import tensorflow as tf
from absl import logging as absl_logging
from alphafold.common import protein
from alphafold.data import pipeline, templates
from alphafold.data.tools import hhsearch
from alphafold.model import model

from local2fold.citations import write_bibtex
from local2fold.colabfold import run_mmseqs2
from local2fold.download import download_alphafold_params
from local2fold.models import load_models_and_params
from local2fold.msa import make_fixed_size
from local2fold.pdb import set_bfactor
from local2fold.plot import plot_predicted_alignment_error, plot_lddt
from local2fold.utils import TqdmHandler

logger = logging.getLogger(__name__)


def mk_mock_template(query_sequence: str, num_temp: int = 1) -> Mapping[str, Any]:
    ln = len(query_sequence)
    output_templates_sequence = "A" * ln
    output_confidence_scores = np.full(ln, 1.0)
    templates_all_atom_positions = np.zeros(
        (ln, templates.residue_constants.atom_type_num, 3)
    )
    templates_all_atom_masks = np.zeros((ln, templates.residue_constants.atom_type_num))
    templates_aatype = templates.residue_constants.sequence_to_onehot(
        output_templates_sequence, templates.residue_constants.HHBLITS_AA_TO_ID
    )
    template_features = {
        "template_all_atom_positions": np.tile(
            templates_all_atom_positions[None], [num_temp, 1, 1, 1]
        ),
        "template_all_atom_masks": np.tile(
            templates_all_atom_masks[None], [num_temp, 1, 1]
        ),
        "template_sequence": [f"none".encode()] * num_temp,
        "template_aatype": np.tile(np.array(templates_aatype)[None], [num_temp, 1, 1]),
        "template_confidence_scores": np.tile(
            output_confidence_scores[None], [num_temp, 1]
        ),
        "template_domain_names": [f"none".encode()] * num_temp,
        "template_release_date": [f"none".encode()] * num_temp,
    }
    return template_features


def mk_template(
    a3m_lines: str, template_paths: str, query_sequence: str
) -> Mapping[str, Any]:
    template_featurizer = templates.TemplateHitFeaturizer(
        mmcif_dir=template_paths,
        max_template_date="2100-01-01",
        max_hits=20,
        kalign_binary_path="kalign",
        release_dates_path=None,
        obsolete_pdbs_path=None,
    )

    hhsearch_pdb70_runner = hhsearch.HHSearch(
        binary_path="hhsearch", databases=[f"{template_paths}/pdb70"]
    )

    hhsearch_result = hhsearch_pdb70_runner.query(a3m_lines)
    hhsearch_hits = pipeline.parsers.parse_hhr(hhsearch_result)
    templates_result = template_featurizer.get_templates(
        query_sequence=query_sequence,
        query_pdb_code=None,
        query_release_date=None,
        hits=hhsearch_hits,
    )
    return templates_result.features


def predict_structure(
    prefix: str,
    result_dir: Path,
    feature_dict: Dict[str, Any],
    sequences_lengths: List[int],
    crop_len: int,
    model_runner_and_params: Dict[str, Tuple[model.RunModel, haiku.Params]],
    do_relax: bool = False,
    random_seed: int = 0,
    cache: bool = False,
):
    """Predicts structure using AlphaFold for the given sequence."""
    # Run the models.
    plddts, paes = [], []
    unrelaxed_pdb_lines = []
    relaxed_pdb_lines = []
    seq_len = feature_dict["seq_length"][0]
    for model_name, (model_runner, params) in model_runner_and_params.items():
        logger.info(f"running {model_name}")
        # swap params to avoid recompiling
        # note: models 1,2 have diff number of params compared to models 3,4,5 (this was handled on construction)
        model_runner.params = params

        processed_feature_dict = model_runner.process_features(
            feature_dict, random_seed=random_seed
        )
        model_config = model_runner.config
        eval_cfg = model_config.data.eval
        crop_feats = {k: [None] + v for k, v in dict(eval_cfg.feat).items()}

        # templates models
        if model_name == "model_1" or model_name == "model_2":
            pad_msa_clusters = eval_cfg.max_msa_clusters - eval_cfg.max_templates
        else:
            pad_msa_clusters = eval_cfg.max_msa_clusters

        max_msa_clusters = pad_msa_clusters

        # let's try pad (num_res + X)
        input_fix = make_fixed_size(
            processed_feature_dict,
            crop_feats,
            msa_cluster_size=max_msa_clusters,  # true_msa (4, 512, 68)
            extra_msa_size=5120,  # extra_msa (4, 5120, 68)
            num_res=crop_len,  # aatype (4, 68)
            num_templates=4,
        )  # template_mask (4, 4) second value

        if cache:
            prediction_result = run_model_cached(
                input_fix, model_name, model_runner, prefix
            )
        else:
            prediction_result = model_runner.predict(input_fix)

        unrelaxed_protein = protein.from_prediction(input_fix, prediction_result)
        unrelaxed_pdb_lines.append(protein.to_pdb(unrelaxed_protein))
        plddts.append(prediction_result["plddt"][:seq_len])
        paes_res = []
        for i in range(seq_len):
            paes_res.append(prediction_result["predicted_aligned_error"][i][:seq_len])
        paes.append(paes_res)
        if do_relax:
            from alphafold.relax import relax
            from alphafold.common import residue_constants

            residue_constants.stereo_chemical_props_path = "stereo_chemical_props.txt"

            # Relax the prediction.
            amber_relaxer = relax.AmberRelaxation(
                max_iterations=0,
                tolerance=2.39,
                stiffness=10.0,
                exclude_residues=[],
                max_outer_iterations=20,
            )
            relaxed_pdb_str, _, _ = amber_relaxer.process(prot=unrelaxed_protein)
            relaxed_pdb_lines.append(relaxed_pdb_str)

    idx_res = feature_dict["residue_index"]
    chains = list(
        "".join([ascii_uppercase[n] * L for n, L in enumerate(sequences_lengths)])
    )

    # rerank models based on predicted lddt
    lddt_rank = np.mean(plddts, -1).argsort()[::-1]
    out = {}
    logger.info("reranking models based on avg. predicted lDDT")
    for n, r in enumerate(lddt_rank):
        logger.info(f"model_{n + 1} {np.mean(plddts[r])}")

        unrelaxed_pdb_path = result_dir.joinpath(
            f"{prefix}_unrelaxed_model_{n + 1}.pdb"
        )
        unrelaxed_pdb_path.write_text(unrelaxed_pdb_lines[r])
        set_bfactor(unrelaxed_pdb_path, plddts[r], idx_res, chains)

        if do_relax:
            relaxed_pdb_path = result_dir.joinpath(
                f"{prefix}_relaxed_model_{n + 1}.pdb"
            )
            relaxed_pdb_path.write_text(unrelaxed_pdb_lines[r])
            set_bfactor(relaxed_pdb_path, plddts[r], idx_res, chains)

        out[f"model_{n + 1}"] = {"plddt": plddts[r], "pae": paes[r]}
    return out


def run_model_cached(
    input_fix: dict, model_name: str, model_runner: model.RunModel, prefix: str
):
    """Caching the expensive compilation + prediction step - for development only

    We store both input and output to ensure that the input is actually the same that we cached
    """
    pickle_path = Path("pickle").joinpath(prefix)
    if pickle_path.joinpath(f"{model_name}_input_fix.pkl").is_file():
        logger.info("Using cached computation")
        with pickle_path.joinpath(f"{model_name}_input_fix.pkl").open("rb") as fp:
            input_fix2 = pickle.load(fp)
            # Make sure we're actually predicting the same input again
            numpy.testing.assert_equal(input_fix, input_fix2)
        with pickle_path.joinpath(f"{model_name}_prediction_result.pkl").open(
            "rb"
        ) as fp:
            prediction_result = pickle.load(fp)
    else:
        # The actual operation that we cache
        prediction_result = model_runner.predict(input_fix)

        pickle_path.mkdir(parents=True, exist_ok=True)
        with pickle_path.joinpath(f"{model_name}_input_fix.pkl").open("wb") as fp:
            pickle.dump(input_fix, fp)
        with pickle_path.joinpath(f"{model_name}_prediction_result.pkl").open(
            "wb"
        ) as fp:
            pickle.dump(prediction_result, fp)
    return prediction_result


def run(
    input_dir: Union[str, Path],
    result_dir: Union[str, Path],
    use_templates: bool,
    use_amber: bool,
    use_env: bool,
    num_models: int,
    homooligomer: int,
    do_not_overwrite_results: bool,
    cache: bool = False,
):
    # hiding warning messages
    warnings.filterwarnings("ignore")
    absl_logging.set_verbosity("error")
    os.environ["TF_CPP_MIN_LOG_LEVEL"] = "3"
    tf.get_logger().setLevel("ERROR")

    input_dir = Path(input_dir)
    result_dir = Path(result_dir)
    result_dir.mkdir(exist_ok=True)
    # TODO: What's going on with MSA mode?
    write_bibtex(True, use_env, use_templates, use_amber, result_dir)

    queries = []
    for file in input_dir.iterdir():
        if not file.is_file():
            continue
        (seqs, header) = pipeline.parsers.parse_fasta(file.read_text())
        query_sequence = seqs[0]
        if len(seqs) > 1:
            logger.warning(
                f"More than one sequence in {file}, ignoring all but the first sequence"
            )
        queries.append((file.stem, query_sequence, file.suffix))
    # sort by seq. len
    queries.sort(key=lambda t: len(t[1]))

    crop_len = math.ceil(len(queries[0][1]) * 1.1)

    model_runner_and_params = load_models_and_params(num_models)

    for jobname, query_sequence, extension in queries:
        a3m_file = f"{jobname}.a3m"
        if len(query_sequence) > crop_len:
            crop_len = math.ceil(len(query_sequence) * 1.1)
        logger.info(f"Running: {jobname}")
        if (
            do_not_overwrite_results
            and result_dir.joinpath(jobname).with_suffix(".result.zip").is_file()
        ):
            continue
        if use_templates:
            try:
                a3m_lines, template_paths = run_mmseqs2(
                    query_sequence,
                    str(result_dir.joinpath(jobname)),
                    use_env,
                    use_templates=True,
                )
            except Exception as e:
                logger.exception(f"{jobname} could not be processed: {e}")
                continue
            if template_paths is None:
                template_features = mk_mock_template(query_sequence, 100)
            else:
                template_features = mk_template(
                    a3m_lines, template_paths, query_sequence
                )
            if extension.lower() == ".a3m":
                a3m_lines = "".join(input_dir.joinpath(a3m_file).read_text())
        else:
            if extension.lower() == ".a3m":
                a3m_lines = "".join(input_dir.joinpath(a3m_file).read_text())
            else:
                try:
                    a3m_lines = run_mmseqs2(
                        query_sequence, str(result_dir.joinpath(jobname)), use_env
                    )
                except Exception as e:
                    logger.exception(f"{jobname} could not be processed: {e}")
                    continue
            template_features = mk_mock_template(query_sequence, 100)

        result_dir.joinpath(a3m_file).write_text(a3m_lines)
        # parse MSA
        msa, deletion_matrix = pipeline.parsers.parse_a3m(a3m_lines)

        # Gather input features, predict structure
        msas = [msa]
        deletion_matrices = [deletion_matrix]
        try:
            # gather features
            feature_dict = {
                **pipeline.make_sequence_features(
                    sequence=query_sequence,
                    description="none",
                    num_res=len(query_sequence),
                ),
                **pipeline.make_msa_features(
                    msas=msas, deletion_matrices=deletion_matrices
                ),
                **template_features,
            }
        except Exception as e:
            logger.exception(f"{jobname} could not be processed: {e}")
            continue

        outs = predict_structure(
            jobname,
            result_dir,
            feature_dict,
            sequences_lengths=[len(query_sequence)],
            crop_len=crop_len,
            model_runner_and_params=model_runner_and_params,
            do_relax=use_amber,
            cache=cache,
        )

        plot_lddt(homooligomer, jobname, msa, outs, query_sequence, result_dir)
        plot_predicted_alignment_error(jobname, num_models, outs, result_dir)


def main():
    logging.basicConfig(
        level=logging.INFO, format="%(asctime)s %(message)s", handlers=[TqdmHandler()]
    )

    download_alphafold_params()

    parser = ArgumentParser()
    parser.add_argument("--input-dir", default="input")
    parser.add_argument("--result-dir", default="results")
    # TODO: This currently isn't actually used
    parser.add_argument(
        "--msa-mode",
        default="MMseqs2 (UniRef+Environmental)",
        choices=[
            "MMseqs2 (UniRef+Environmental)",
            "MMseqs2 (UniRef only)",
            "single_sequence",
            "custom",
        ],
    )
    parser.add_argument("--use-amber", default=False, action="store_true")
    parser.add_argument("--use-templates", default=False, action="store_true")
    parser.add_argument("--use-env", default=False, action="store_true")
    parser.add_argument(
        "--cache",
        default=False,
        action="store_true",
        help="Caches the model output. For development only",
    )
    parser.add_argument("--num-models", type=int, default=5, choices=[1, 2, 3, 4, 5])
    parser.add_argument("--homooligomer", type=int, default=1)
    parser.add_argument(
        "--do-not-overwrite-results", default=True, action="store_false"
    )
    args = parser.parse_args()

    assert args.msa_mode == "MMseqs2 (UniRef+Environmental)", "Unsupported"
    assert args.homooligomer == 1, "Unsupported"

    run(
        args.input_dir,
        args.result_dir,
        args.use_templates,
        args.use_amber,
        args.use_env,
        args.num_models,
        args.homooligomer,
        args.do_not_overwrite_results,
    )


if __name__ == "__main__":
    main()
