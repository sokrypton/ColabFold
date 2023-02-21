import time
import pickle
from pathlib import Path
import numpy as np
import json

import logging
logger = logging.getLogger(__name__)
from typing import Any, Callable, Dict, List, Optional, Tuple, Union

import haiku, jax

from colabfold.inputs import (
  pad_input, pad_input_multimer,
  model, protein
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
  save_best: bool = False,
  cyclic: bool = False
):
  """Predicts structure using AlphaFold for the given sequence."""

  if num_relax > 0:
    from colabfold.alphafold.relax import run_relax

  mean_scores = []
  conf = []
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
      def cyclic_offset(L):
        i = np.arange(L)
        ij = np.stack([i,i+L],-1)
        offset = i[:,None] - i[None,:]
        c_offset = np.abs(ij[:,None,:,None] - ij[None,:,None,:]).min((2,3))
        return np.sign(offset) * c_offset
      if "multimer" in model_type:
        if model_num == 0 and seed_num == 0:
          input_features = feature_dict
          input_features["asym_id"] = input_features["asym_id"] - input_features["asym_id"][...,0]
          if cyclic:
            input_features["offset"] = cyclic_offset(seq_len)
          
          # TODO
          # if seq_len < pad_len:
          #   input_features = pad_input_multimer(input_features, model_runner, model_name, pad_len, use_templates)
          #   logger.info(f"Padding length to {pad_len}")
      else:
        if model_num == 0:
          input_features = model_runner.process_features(feature_dict, random_seed=seed)            
          r = input_features["aatype"].shape[0]
          input_features["asym_id"] = np.tile(feature_dict["asym_id"],r).reshape(r,-1)
          
          if cyclic:
            B = input_features["aatype"].shape[0]
            input_features["offset"] = np.tile(cyclic_offset(seq_len)[None],(B,1,1))
          
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
          
          files.get("unrelaxed",f"r{recycles}.pdb").write_text(protein.to_pdb(unrelaxed_protein))
        
          if save_all:
            with files.get("all",f"r{recycles}.pickle").open("wb") as handle:
              pickle.dump(result, handle)

      ########################
      # predict
      ########################
      scores_path = files.get("scores","json")
      if scores_path.is_file(): 
        # if score file already exists, log scores and continue        
        
        with scores_path.open("r") as handle:
          scores = json.load(handle)
                  
        mean_scores.append(float(scores["ranking_confidence"]))
        
        print_line = ""
        conf.append({})
        for x,y in [["plddt","pLDDT"],["ptm","pTM"],["iptm","ipTM"]]:
          if x in scores:
            conf[-1][x] = float(np.mean(scores[x]) if x == "plddt" else scores[x])
            print_line += f" {y}={conf[-1][x]:.3g}"
        conf[-1]["print_line"] = print_line
        logger.info(f"{tag}{print_line}")

        files.get("unrelaxed","pdb")
        if save_all: files.get("all","pickle")
        if save_single_representations: files.get("single_repr","npy")
        if save_pair_representations: files.get("pair_repr","npy")

      else:
        ###########################################################       
        start = time.time()
        result, recycles = \
        _predict(
          self=model_runner,
          feat=input_features,
          random_seed=seed,
          return_representations=(save_all or save_single_representations or save_pair_representations),
          callback=callback,
          return_best=save_best,
        )
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
        files.get("unrelaxed","pdb").write_text(protein.to_pdb(unrelaxed_protein))

        # save raw outputs
        if save_all:
          with files.get("all","pickle").open("wb") as handle:
            pickle.dump(result, handle)      
        if save_single_representations:
          np.save(files.get("single_repr","npy"),result["representations"]["single"])
        if save_pair_representations:
          np.save(files.get("pair_repr","npy"),result["representations"]["pair"])

        # write an easy-to-use format (pAE and pLDDT)
        with scores_path.open("w") as handle:
          plddt = result["plddt"][:seq_len]
          scores = {
            "plddt": np.around(plddt.astype(float), 2).tolist(),
            "ranking_confidence": np.around(result["ranking_confidence"], 2).tolist()
          }
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
        ###########################################################       

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
      pdb_filename = result_dir.joinpath(f"{prefix}_unrelaxed_{tag}.pdb")
      pdb_lines = run_relax(pdb_filename=pdb_filename, use_gpu=use_gpu_relax)
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

def _predict(self,
  feat,
  random_seed=0,
  return_representations=False,
  callback=None,
  return_best=False):
  '''
  This is a copy of the predict function from model.py (alphafold).
  Moving here to make easier to modify without making edits to alphafold branch.
  '''
  
  self.init_params(feat)

  # get shapes
  aatype = feat["aatype"]
  num_iters = self.config.model.num_recycle + 1
  if self.multimer_mode:
    L = aatype.shape[0]
  else:
    num_ensemble = self.config.data.eval.num_ensemble
    L = aatype.shape[1]
  
  # initialize
  zeros = lambda shape: np.zeros(shape, dtype=np.float16)
  prev = {'prev_msa_first_row': zeros([L,256]),
          'prev_pair':          zeros([L,L,128]),
          'prev_pos':           zeros([L,37,3])}
  
  def run(key, feat, prev):
    def _jnp_to_np(x):
      for k, v in x.items():
        if isinstance(v, dict):
          x[k] = _jnp_to_np(v)
        else:
          x[k] = np.asarray(v,np.float16)
      return x
    result = _jnp_to_np(self.apply(self.params, key, {**feat, "prev":prev}))
    prev = result.pop("prev")
    return result, prev

  # initialize random key
  key = jax.random.PRNGKey(random_seed)
  
  # iterate through recycles
  best_result, best_score, best_r = {}, 0, 0
  for r in range(num_iters):      
    # grab subset of features
    if self.multimer_mode:
      sub_feat = feat
    else:
      s = r * num_ensemble
      e = (r+1) * num_ensemble
      sub_feat = jax.tree_map(lambda x:x[s:e], feat)
        
    # run
    key, sub_key = jax.random.split(key)
    result, prev = run(sub_key, sub_feat, prev)
    
    if return_representations:
      result["representations"] = {"pair":   prev["prev_pair"],
                                   "single": prev["prev_msa_first_row"]}                                     
    # callback
    if callback is not None: callback(result, r)

    if return_best and result["ranking_confidence"] > best_score:
      best_score = result["ranking_confidence"]
      best_result = jax.tree_map(lambda x:x, result)
      best_r = r

    # decide when to stop
    if result["ranking_confidence"] > self.config.model.stop_at_score:
      break
    if r > 0 and result["tol"] < self.config.model.recycle_early_stop_tolerance:
      break
  
  if return_best:
    return best_result, best_r
  else:
    return result, r