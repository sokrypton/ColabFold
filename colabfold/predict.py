import time
import pickle
from pathlib import Path
import numpy as np
import json

import logging
logger = logging.getLogger(__name__)
from typing import Any, Callable, Dict, List, Optional, Tuple, Union

from colabfold.inputs import (
  pad_input, pad_input_multimer,
  haiku, model, protein
)
from colabfold.utils import file_manager

################################
# main prediction function
################################

def predict_structure(
  prefix: str,
  result_dir: Path,
  feature_dict: Dict[str, Any],
  is_complex: bool,
  use_templates: bool,
  sequences_lengths: List[int],
  pad_len: int,
  model_type: str,
  model_runner_and_params: List[Tuple[str, model.RunModel, haiku.Params]],
  num_relax: int = 0,
  rank_by: str = "auto",
  random_seed: int = 0,
  num_seeds: int = 1,
  stop_at_score: float = 100,
  outputs_callback: Callable[[Any], Any] = None,
  use_gpu_relax: bool = False,
  save_all: bool = False,
  save_single_representations: bool = False,
  save_pair_representations: bool = False,
  save_recycles: bool = False,
):
  """Predicts structure using AlphaFold for the given sequence."""

  if num_relax > 0:
    from colabfold.alphafold.relax import run_relax

  mean_scores = []
  conf = []
  unrelaxed_pdb_lines = []
  prediction_times = []
  seq_len = sum(sequences_lengths)
  model_names = []
  files = file_manager(prefix, result_dir)
  
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
          # TODO
          if seq_len < pad_len:
            input_features = pad_input_multimer(input_features, model_runner, model_name, pad_len, use_templates)
            logger.info(f"Padding length to {pad_len}")
      else:
        if model_num == 0:
          input_features = model_runner.process_features(feature_dict, random_seed=seed)            
          r = input_features["aatype"].shape[0]
          input_features["asym_id"] = np.tile(feature_dict["asym_id"],r).reshape(r,-1)
          if seq_len < pad_len:
            input_features = pad_input(input_features, model_runner, model_name, pad_len, use_templates)
            logger.info(f"Padding length to {pad_len}")

      tag = f"{model_type}_{model_name}_seed_{seed:03d}"
      model_names.append(tag)
      files.set_tag(tag)

      # monitor intermediate results
      def callback(result, recycles):
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
          unrelaxed_protein = protein.from_prediction(features=input_features,
            result=result, b_factors=b_factors,
            remove_leading_feature_dimension=("multimer" not in model_type))
          
          unrelaxed_pdb_lines = protein.to_pdb(unrelaxed_protein)
          files.get("unrelaxed",f"r{recycles}.pdb").write_text(unrelaxed_pdb_lines)
        
          if save_all:
            with files.get("all",f"r{recycles}.pickle").open("wb") as handle:
              pickle.dump(result, handle)

      ########################
      # predict
      ########################
      start = time.time()
      result, recycles = \
      model_runner.predict(input_features,
        random_seed=seed,
        return_representations=(save_all or save_single_representations or save_pair_representations),
        callback=callback)
      prediction_times.append(time.time() - start)

      ########################
      # parse results
      ########################
      # summary metrics
      mean_scores.append(result["ranking_confidence"])     
      if recycles == 0: result.pop("tol",None)
      if not is_complex: result.pop("iptm",None)
      print_line = ""
      conf.append({})
      for x,y in [["mean_plddt","pLDDT"],["ptm","pTM"],["iptm","ipTM"]]:
        if x in result:
          print_line += f" {y}={result[x]:.3g}"
          conf[-1][x] = float(result[x])
      conf[-1]["print_line"] = print_line
      logger.info(f"{tag} took {prediction_times[-1]:.1f}s ({recycles} recycles)")

      # create protein object
      final_atom_mask = result["structure_module"]["final_atom_mask"]
      b_factors = result["plddt"][:, None] * final_atom_mask
      unrelaxed_protein = protein.from_prediction(
        features=input_features,
        result=result,
        b_factors=b_factors,
        remove_leading_feature_dimension=("multimer" not in model_type))

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
          pickle.dump(result, handle)      
      if save_single_representations:
        np.save(files.get("single_repr","npy"),result["representations"]["single"])
      if save_pair_representations:
        np.save(files.get("pair_repr","npy"),result["representations"]["pair"])

      # write an easy-to-use format (pAE and pLDDT)
      with files.get("scores","json").open("w") as handle:
        plddt = result["plddt"][:seq_len]
        scores = {"plddt": np.around(plddt.astype(float), 2).tolist()}
        if "predicted_aligned_error" in result:
          pae = result["predicted_aligned_error"][:seq_len,:seq_len]
          scores.update({"max_pae": pae.max().astype(float).item(),
                         "pae": np.around(pae.astype(float), 2).tolist()})
          for k in ["ptm","iptm"]:
            if k in conf[-1]: scores[k] = np.around(conf[-1][k], 2).item()
          del pae
        del plddt
        json.dump(scores, handle)

      ###############################
      # callback for visualization
      ###############################
      if outputs_callback is not None:
        outputs_callback({
          "unrelaxed_protein":unrelaxed_protein,
          "sequences_lengths":sequences_lengths,
          "result":result,
          "input_features":input_features,
          "tag":tag,
          "files":files.files[tag]})

      del result, unrelaxed_protein

      # early stop criteria fulfilled
      if mean_scores[-1] > stop_at_score: break

    # early stop criteria fulfilled
    if mean_scores[-1] > stop_at_score: break
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
      pdb_lines = run_relax(pdb_lines=unrelaxed_pdb_lines[key], use_gpu=use_gpu_relax)
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
    
  return {"rank":rank, "metric":metric, "result_files":result_files}