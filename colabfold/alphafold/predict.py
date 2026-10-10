"""
AlphaFold2 inference
"""
import logging
import pickle
import time
from pathlib import Path
from typing import Any, Callable, Dict, List, Tuple, TYPE_CHECKING

import numpy as np

from alphafold.common import protein
from alphafold.data import pipeline_multimer

from colabfold.alphafold import extra_ptm, ipsae
from colabfold.input import pdb_to_string
from colabfold.backend import metrics_line
from colabfold.relax import relax_me

if TYPE_CHECKING:
    import haiku
    from alphafold.model import model

logger = logging.getLogger(__name__)

try:
    import orjson
    hasOrjson = True
except ImportError:
    import json
    hasOrjson = False


def pad_input(
    input_features: "model.features.FeatureDict",
    model_runner: "model.RunModel",
    model_name: str,
    pad_len: int,
    use_templates: bool,
) -> "model.features.FeatureDict":
    from colabfold.alphafold.msa import make_fixed_size

    model_config = model_runner.config
    eval_cfg = model_config.data.eval
    crop_feats = {k: [None] + v for k, v in dict(eval_cfg.feat).items()}

    max_msa_clusters = eval_cfg.max_msa_clusters
    max_extra_msa = model_config.data.common.max_extra_msa
    # templates models
    if (model_name == "model_1" or model_name == "model_2") and use_templates:
        pad_msa_clusters = max_msa_clusters - eval_cfg.max_templates
    else:
        pad_msa_clusters = max_msa_clusters

    max_msa_clusters = pad_msa_clusters

    # let's try pad (num_res + X)
    input_fix = make_fixed_size(
        input_features,
        crop_feats,
        msa_cluster_size=max_msa_clusters,  # true_msa (4, 512, 68)
        extra_msa_size=max_extra_msa,  # extra_msa (4, 5120, 68)
        num_res=pad_len,  # aatype (4, 68)
        num_templates=4,
    )  # template_mask (4, 4) second value
    return input_fix


# Residue axis of each AF2-multimer input feature (None: has none); padding is masked by seq_mask and msa_mask.
_MULTIMER_RESIDUE_AXIS = {
    "aatype": 0, "all_atom_mask": 0, "all_atom_positions": 0, "asym_id": 0, "deletion_mean": 0,
    "entity_id": 0, "entity_mask": 0, "residue_index": 0, "seq_mask": 0, "sym_id": 0,
    "bert_mask": 1, "deletion_matrix": 1, "msa": 1, "msa_mask": 1,
    "template_aatype": 1, "template_all_atom_mask": 1, "template_all_atom_positions": 1,
    "assembly_num_chains": None, "cluster_bias_mask": None, "num_alignments": None, "num_templates": None,
    "seq_length": None,
}


def pad_input_multimer(input_features: Dict[str, Any], pad_len: int) -> Dict[str, Any]:
    """Zero-pad the residue axes of AF2-multimer features to pad_len; an unknown feature raises KeyError."""
    def pad(v, axis):
        if axis is None:
            return v
        widths = [(0, 0)] * np.ndim(v)
        widths[axis] = (0, pad_len - np.shape(v)[axis])
        return np.pad(v, widths)
    return {k: pad(v, _MULTIMER_RESIDUE_AXIS[k]) for k, v in input_features.items()}


# Outputs with a second residue axis; all others have one, axis 0 or after the layer/MSA axis.
_PAIR_OUTPUTS = {"predicted_aligned_error", "aligned_confidence_probs", "pae_matrix_with_logits", "pair",
                 "distogram"}


def crop_residues(tree: Any, seq_len: int, path: Tuple[str, ...] = ()) -> Any:
    """Crop the padded residue axes of model outputs back to seq_len."""
    if isinstance(tree, dict):
        return {k: crop_residues(v, seq_len, path + (k,)) for k, v in tree.items()}
    if np.ndim(tree) == 0 or path[-1] in ("bin_edges", "breaks"):
        return tree
    first = 1 if {"traj", "sidechains", "masked_msa"} & set(path) else 0
    axes = (first, first + 1) if _PAIR_OUTPUTS & set(path) else (first,)
    return tree[tuple(slice(0, seq_len) if a in axes else slice(None) for a in range(np.ndim(tree)))]


class file_manager:
    def __init__(self, prefix: str, result_dir: Path):
        self.prefix = prefix
        self.result_dir = result_dir
        self.tag = None
        self.files = {}

    def get(self, x: str, ext:str) -> Path:
        if self.tag not in self.files:
            self.files[self.tag] = []
        file = self.result_dir.joinpath(f"{self.prefix}_{x}_{self.tag}.{ext}")
        self.files[self.tag].append([x,ext,file])
        return file

    def set_tag(self, tag):
        self.tag = tag


def predict_structure(
    prefix: str,
    result_dir: Path,
    feature_dict: Dict[str, Any],
    is_complex: bool,
    use_templates: bool,
    sequences_lengths: List[int],
    pad_len: int,
    model_type: str,
    model_runner_and_params: List[Tuple[str, "model.RunModel", "haiku.Params"]],
    initial_guess: str = None,
    msa_pad_depth: int = 0,
    num_relax: int = 0,
    relax_max_iterations: int = 0,
    relax_tolerance: float = 2.39,
    relax_stiffness: float = 10.0,
    relax_max_outer_iterations: int = 3,
    rank_by: str = "auto",
    random_seed: int = 0,
    num_seeds: int = 1,
    stop_at_score: float = 100,
    prediction_callback: Callable[[Any, Any, Any, Any, Any], Any] = None,
    use_gpu_relax: bool = False,
    save_all: bool = False,
    save_single_representations: bool = False,
    save_pair_representations: bool = False,
    save_recycles: bool = False,
    calc_extra_ptm: bool = False,
    use_probs_extra: bool = True,
):
    """Predicts structure using AlphaFold for the given sequence."""
    mean_scores = []
    conf = []
    unrelaxed_pdb_lines = []
    prediction_times = []
    model_names = []
    files = file_manager(prefix, result_dir)
    seq_len = sum(sequences_lengths)

    # iterate through random seeds
    for seed_num, seed in enumerate(range(random_seed, random_seed+num_seeds)):

        # iterate through models
        for model_num, (model_name, model_runner, params) in enumerate(model_runner_and_params):

            # swap params to avoid recompiling
            model_runner.params = params

            #########################
            # process input features
            #########################
            if "multimer" in model_type:
                if model_num == 0 and seed_num == 0:
                    input_features = feature_dict
                    input_features["asym_id"] = input_features["asym_id"] - input_features["asym_id"][...,0]
                    # Pad MSA depth up to msa_pad_depth to share one compiled shape
                    # Padded rows are masked (msa_mask=0), so the prediction is unchanged
                    if msa_pad_depth > input_features["msa"].shape[0]:
                        input_features = pipeline_multimer.pad_msa(input_features, min_num_seq=msa_pad_depth)
                    if seq_len < pad_len:
                        logger.info(f"Padding length to {pad_len}")
            else:
                if model_num == 0:
                    input_features = model_runner.process_features(feature_dict, random_seed=seed)
                    r = input_features["aatype"].shape[0]
                    input_features["asym_id"] = np.tile(feature_dict["asym_id"],r).reshape(r,-1)
                    if seq_len < pad_len:
                        input_features = pad_input(input_features, model_runner,
                            model_name, pad_len, use_templates)
                        logger.info(f"Padding length to {pad_len}")


            tag = f"{model_type}_{model_name}_seed_{seed:03d}"
            model_names.append(tag)
            files.set_tag(tag)

            # initial guess
            if initial_guess:
                input_guess = Path(initial_guess)
                if input_guess.suffix == ".pdb":
                    positions = protein.from_pdb_string(pdb_to_string(initial_guess)).atom_positions
                elif input_guess.suffix == ".cif":
                    positions = protein.from_mmcif_string(input_guess.read_text()).atom_positions
                else:
                    raise ValueError(f"Unsupported initial guess file format: {initial_guess}")
                # a padded monomer's features are already pad_len long
                pad = input_features["aatype"].shape[-1] - len(positions)
                input_features["all_atom_positions"] = np.pad(positions, [(0, pad), (0, 0), (0, 0)])

            # Padded residues are masked; multimer outputs are cropped back right after prediction
            crop_output = "multimer" in model_type and seq_len < pad_len
            model_features = pad_input_multimer(input_features, pad_len) if crop_output else input_features


            ########################
            # predict
            ########################
            start = time.time()

            # monitor intermediate results
            def callback(result, recycles):
                if crop_output: result = crop_residues(result, seq_len)
                if recycles == 0: result.pop("tol",None)
                if not is_complex: result.pop("iptm",None)
                print_line = ""
                for x,y in [["mean_plddt","pLDDT"],["ptm","pTM"],["iptm","ipTM"],["tol","tol"]]:
                  if x in result:
                    print_line += f" {y}={result[x]:.3g}"
                logger.info(f"{tag} recycle={recycles}{print_line}")

                if save_recycles:
                    final_atom_mask = result["structure_module"]["final_atom_mask"]
                    b_factors = result["plddt"][:, None] * final_atom_mask
                    unrelaxed_protein = protein.from_prediction(
                        features=input_features,
                        result=result, b_factors=b_factors,
                        remove_leading_feature_dimension=("multimer" not in model_type))
                    files.get("unrelaxed",f"r{recycles}.pdb").write_text(protein.to_pdb(unrelaxed_protein))

                    if save_all:
                        with files.get("all",f"r{recycles}.pickle").open("wb") as handle:
                            pickle.dump(crop_residues(result, seq_len), handle)
                    del unrelaxed_protein

            return_representations = save_all or save_single_representations or save_pair_representations

            # predict
            result, recycles = \
            model_runner.predict(model_features,
                random_seed=seed,
                return_representations=return_representations,
                callback=callback)
            if crop_output:
                result = crop_residues(result, seq_len)

            if calc_extra_ptm and 'predicted_aligned_error' in result.keys():
                extra_ptm_output = extra_ptm.get_chain_and_interface_metrics(result, input_features['asym_id'],
                    use_probs_extra=use_probs_extra,
                    use_jnp=False)
                result.pop('pae_matrix_with_logits', None)
                result['actifptm'] = extra_ptm_output['actifptm']
            else:
                calc_extra_ptm = False
            prediction_times.append(time.time() - start)

            ########################
            # parse results
            ########################

            # summary metrics
            mean_scores.append(result["ranking_confidence"])
            if recycles == 0: result.pop("tol",None)
            if not is_complex: result.pop("iptm",None)
            conf.append({x: float(result[x])
                         for x in ("mean_plddt", "ptm", "iptm", "actifptm") if x in result})
            conf[-1]["print_line"] = metrics_line(result)
            logger.info(f"{tag} took {prediction_times[-1]:.1f}s ({recycles} recycles)")

            # create protein object
            final_atom_mask = result["structure_module"]["final_atom_mask"]
            b_factors = result["plddt"][:, None] * final_atom_mask
            unrelaxed_protein = protein.from_prediction(
                features=input_features,
                result=result,
                b_factors=b_factors,
                remove_leading_feature_dimension=("multimer" not in model_type))

            # callback for visualization
            if prediction_callback is not None:
                prediction_callback(unrelaxed_protein, sequences_lengths,
                                    result, input_features, (tag, False))

            #########################
            # save results
            #########################

            # save pdb
            protein_lines = protein.to_pdb(unrelaxed_protein)
            files.get("unrelaxed","pdb").write_text(protein_lines)
            unrelaxed_pdb_lines.append(protein_lines)

            # save raw outputs
            if save_all:
                with files.get("all","pickle").open("wb") as handle:
                    pickle.dump(crop_residues(result, seq_len), handle)
            if save_single_representations or save_pair_representations:
                # a padded monomer's outputs keep the padding until here
                representations = crop_residues(result["representations"], seq_len, ("representations",))
            if save_single_representations:
                np.save(files.get("single_repr","npy"),representations["single"])
            if save_pair_representations:
                np.save(files.get("pair_repr","npy"),representations["pair"])

            # write an easy-to-use format (pAE and pLDDT)
            plddt = result["plddt"][:seq_len]
            scores = {"plddt": np.around(plddt.astype(float), 2).tolist()}
            if "predicted_aligned_error" in result:
                pae = result["predicted_aligned_error"][:seq_len,:seq_len]
                scores.update({"max_pae": pae.max().astype(float).item(),
                                "pae": np.around(pae.astype(float), 2).tolist()})
                if calc_extra_ptm:
                    scores.update(extra_ptm_output)
                for k in ["ptm", "iptm"]:
                    if k in conf[-1]:
                        scores[k] = np.around(conf[-1][k], 2).item()
                if is_complex:
                    try:
                        asym_id = input_features["asym_id"]
                        if asym_id.ndim > 1: asym_id = asym_id[0]
                        interface_scores = ipsae.get_interface_scores(
                            pae=pae,
                            plddt=plddt,
                            asym_id=asym_id[:seq_len],
                            atom_positions=result["structure_module"]["final_atom_positions"][:seq_len],
                            atom_mask=result["structure_module"]["final_atom_mask"][:seq_len])
                        scores.update(interface_scores)
                        if interface_scores:
                            conf[-1]["print_line"] += (
                                f" ipSAE={ipsae.format_ipsae(interface_scores['ipsae'])}"
                                f" pDockQ2={ipsae.format_ipsae(interface_scores['pdockq2'])}"
                            )
                    except Exception as e:
                        logger.warning(f"Could not compute ipSAE/pDockQ interface scores: {e}")
                del pae
            del plddt
            file = files.get("scores", "json")
            if hasOrjson:
                file.write_bytes(orjson.dumps(scores))
            else:
                file.write_text(json.dumps(scores))

            del result, unrelaxed_protein

            # early stop criteria fulfilled
            if mean_scores[-1] > stop_at_score: break

        # early stop criteria fulfilled
        if mean_scores[-1] > stop_at_score: break

        # cleanup
        if "multimer" not in model_type: del input_features
    if "multimer" in model_type: del input_features

    ###################################################
    # rerank models based on predicted confidence
    ###################################################

    rank, metric = [],[]
    result_files = []
    logger.info(f"reranking models by '{rank_by}' metric")
    model_rank = np.array(mean_scores).argsort()[::-1]
    for n, key in enumerate(model_rank):
        metric.append(conf[key])
        tag = model_names[key]
        files.set_tag(tag)
        # save relaxed pdb
        if n < num_relax:
            start = time.time()
            pdb_lines = relax_me(
                pdb_lines=unrelaxed_pdb_lines[key],
                max_iterations=relax_max_iterations,
                tolerance=relax_tolerance,
                stiffness=relax_stiffness,
                max_outer_iterations=relax_max_outer_iterations,
                use_gpu=use_gpu_relax)
            files.get("relaxed","pdb").write_text(pdb_lines)
            logger.info(f"Relaxation took {(time.time() - start):.1f}s")

        # rename files to include rank
        new_tag = f"rank_{(n+1):03d}_{tag}"
        rank.append(new_tag)
        logger.info(f"{new_tag}{metric[-1]['print_line']}")
        for x, ext, file in files.files[tag]:
            new_file = result_dir.joinpath(f"{prefix}_{x}_{new_tag}.{ext}")
            file.rename(new_file)
            result_files.append(new_file)

    return {"rank":rank,
            "metric":metric,
            "result_files":result_files}
