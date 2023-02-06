from __future__ import annotations
import os
os.environ["TF_FORCE_UNIFIED_MEMORY"] = "1"
os.environ["XLA_PYTHON_CLIENT_MEM_FRACTION"] = "2.0"

# test if alphafold installed
try:
  import alphafold
except ModuleNotFoundError:
  raise RuntimeError("alphafold is not installed. Please run `pip install colabfold[alphafold]`")

import json
import math
import sys
import time
import zipfile
import pickle
import argparse
import importlib_metadata
import numpy as np
from typing import Any, Callable, Dict, List, Optional, Tuple, Union, TYPE_CHECKING
from pathlib import Path
import random

import logging
logger = logging.getLogger(__name__)

# import from colabfold
from colabfold.inputs import (
  mk_hhsearch_db, generate_input_feature, msa_to_str, parse_fasta,
  pad_input, get_queries, unserialize_msa, 
  pad_input_multimer, get_msa_and_templates, unpack_a3ms, get_queries_pairwise,
  jnp, haiku, 

  # alphafold specific
  model, protein, residue_constants
)

from colabfold.download import (
  default_data_dir,
  download_alphafold_params)

from colabfold.utils import (
  DEFAULT_API_SERVER,
  ACCEPT_DEFAULT_TERMS,
  get_commit,
  setup_logging,
  safe_filename,
  file_manager, jnp_to_np, class_to_np)

from colabfold.citations import write_bibtex

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
  model_order: List[int],
  num_relax: int = 0,
  rank_by: str = "auto",
  random_seed: int = 0,
  num_seeds: int = 1,
  stop_at_score: float = 100,
  outputs_callback: Callable[Any] = None,
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

  # process input features
  if "multimer" in model_type:
    input_features = feature_dict
    if pad_len > seq_len: # pad inputs
      input_features = pad_input_multimer(feature_dict, model_runner, model_name, pad_len, use_templates)
    input_features["asym_id"] = input_features["asym_id"] - input_features["asym_id"][...,0]
  
  # iterate through random seeds
  for seed in range(random_seed, random_seed+num_seeds):
    
    # process input features
    if "ptm" in model_type:
      (model_name, model_runner, params) = model_runner_and_params[0]
      input_features = model_runner.process_features(feature_dict, random_seed=seed)      
      r = input_features["aatype"].shape[0]
      input_features["asym_id"] = np.tile(feature_dict["asym_id"],r).reshape(r,-1)
      if pad_len > seq_len: # pad inputs
        input_features = pad_input(input_features, model_runner, model_name, pad_len, use_templates)

    # iterate through models
    for m in model_order:
      (model_name, model_runner, params) = model_runner_and_params[m-1]
      model_runner.params = params

      tag = f"{model_type}_{model_name}_seed_{seed:03d}"
      model_names.append(tag)
      files.set_tag(tag)

      # monitor intermediate results
      def prediction_callback(prediction_result, recycles):
        print_line = ""
        for x,y in [["mean_plddt","pLDDT"],["ptm","pTM"],["iptm","ipTM"],["diff","tol"]]:
          if x in prediction_result:
            print_line += f" {y}={prediction_result[x]:.3g}"
        logger.info(f"{tag} recycle={recycles}{print_line}")
      
        if save_recycles:
          prediction_result = jnp_to_np(prediction_result)
          final_atom_mask = prediction_result["structure_module"]["final_atom_mask"]
          b_factors = prediction_result["plddt"][:, None] * final_atom_mask
          unrelaxed_protein = protein.from_prediction(features=input_features,
            result=prediction_result, b_factors=b_factors,
            remove_leading_feature_dimension=("ptm" in model_type))
          
          unrelaxed_pdb_lines = protein.to_pdb(class_to_np(unrelaxed_protein))
          files.get("unrelaxed",f"r{recycles}.pdb").write_text(unrelaxed_pdb_lines)
        
          if save_all:
            with files.get("all",f"r{recycles}.pickle").open("wb") as handle:
              pickle.dump(prediction_result, handle)

      ########################
      # predict
      ########################
      start = time.time()
      prediction_result, recycles = \
      model_runner.predict(input_features, random_seed=seed, prediction_callback=prediction_callback)
      prediction_result = jnp_to_np(prediction_result)
      prediction_times.append(time.time() - start)

      ########################
      # parse results
      ########################
      # summary metrics
      mean_scores.append(prediction_result["ranking_confidence"])     
      print_line = ""
      conf.append({})
      for x,y in [["mean_plddt","pLDDT"],["ptm","pTM"],["iptm","ipTM"]]:
        if x in prediction_result:
          print_line += f" {y}={prediction_result[x]:.3g}"
          conf[-1][x] = float(prediction_result[x])
      conf[-1]["print_line"] = print_line
      logger.info(f"{tag} took {prediction_times[-1]:.1f}s ({recycles} recycles)")

      # create protein object
      final_atom_mask = prediction_result["structure_module"]["final_atom_mask"]
      b_factors = prediction_result["plddt"][:, None] * final_atom_mask
      unrelaxed_protein = protein.from_prediction(
        features=input_features,
        result=prediction_result,
        b_factors=b_factors,
        remove_leading_feature_dimension=("ptm" in model_type))
      unrelaxed_protein = class_to_np(unrelaxed_protein)

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
          pickle.dump(prediction_result, handle)
      
      if save_single_representations or save_pair_representations:
        rep = prediction_result["prev"]
        if save_single_representations:
          np.save(files.get("single_repr","npy"), rep["prev_msa_first_row"])
        if save_pair_representations:
          np.save(files.get("pair_repr","npy"), rep["prev_pair"])

      # write an easy-to-use format (pAE and pLDDT)
      with files.get("scores","json").open("w") as handle:
        pae = prediction_result["predicted_aligned_error"][:seq_len,:seq_len]
        plddt = prediction_result["plddt"][:seq_len]      
        scores = {
          "max_pae": pae.max().astype(float).item(),
          "pae":   np.around(pae.astype(float), 2).tolist(),
          "plddt": np.around(plddt.astype(float), 2).tolist(),
        }
        for k in ["ptm","iptm"]:
          if k in conf[-1]: scores[k] = np.around(conf[-1][k], 2).item()
        json.dump(scores, handle)

      ###############################
      # callback for visualization
      ###############################
      if outputs_callback is not None:
        outputs_callback({
          "unrelaxed_protein":unrelaxed_protein,
          "sequences_lengths":sequences_lengths,
          "prediction_result":prediction_result,
          "input_features":input_features,
          "tag":tag, "files":files.files[tag]})

      # early stop criteria fulfilled
      if mean_scores[-1] > stop_at_score: break

    # early stop criteria fulfilled
    if mean_scores[-1] > stop_at_score: break

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

def run(
  queries: List[Tuple[str, Union[str, List[str]], Optional[List[str]]]],
  result_dir: Union[str, Path],
  num_models: int,
  is_complex: bool,
  num_recycles: Optional[int] = None,
  recycle_early_stop_tolerance: Optional[float] = None,
  model_order: List[int] = [1,2,3,4,5],
  num_ensemble: int = 1,
  model_type: str = "auto",
  msa_mode: str = "mmseqs2_uniref_env",
  use_templates: bool = False,
  custom_template_path: str = None,
  num_relax: int = 0,
  keep_existing_results: bool = True,
  rank_by: str = "auto",
  pair_mode: str = "unpaired_paired",
  data_dir: Union[str, Path] = default_data_dir,
  host_url: str = DEFAULT_API_SERVER,
  random_seed: int = 0,
  num_seeds: int = 1,
  recompile_padding: Union[int, float] = 10,
  zip_results: bool = False,
  save_single_representations: bool = False,
  save_pair_representations: bool = False,
  save_all: bool = False,
  save_recycles: bool = False,
  use_dropout: bool = False,
  use_gpu_relax: bool = False,
  stop_at_score: float = 100,
  dpi: int = 200,
  max_seq: Optional[int] = None,
  max_extra_seq: Optional[int] = None,
  inputs_callback: Callable[[Any], Any] = None,
  outputs_callback: Callable[[Any], Any] = None,
  **kwargs
):
  # check what device is available
  try:
    # check if TPU is available
    import jax.tools.colab_tpu
    jax.tools.colab_tpu.setup_tpu()
    logger.info('Running on TPU')
    DEVICE = "tpu"
    use_gpu_relax = False
  except:
    if jax.local_devices()[0].platform == 'cpu':
      logger.info("WARNING: no GPU detected, will be using CPU")
      DEVICE = "cpu"
      use_gpu_relax = False
    else:
      import tensorflow as tf
      logger.info('Running on GPU')
      DEVICE = "gpu"
      # disable GPU on tensorflow
      tf.config.set_visible_devices([], 'GPU')

  from alphafold.notebooks.notebook_utils import get_pae_json
  from colabfold.alphafold.models import load_models_and_params
  from colabfold.plot import (plot_plddts, plot_paes, plot_msa)

  data_dir = Path(data_dir)
  result_dir = Path(result_dir)
  result_dir.mkdir(exist_ok=True)
  model_type = set_model_type(is_complex, model_type)

  # determine model extension
  if   model_type == "alphafold2_multimer_v1": model_suffix = "_multimer"
  elif model_type == "alphafold2_multimer_v2": model_suffix = "_multimer_v2"
  elif model_type == "alphafold2_multimer_v3": model_suffix = "_multimer_v3"
  elif model_type == "alphafold2_ptm": model_suffix = "_ptm"
  else: raise ValueError(f"Unknown model_type {model_type}")

  # backward-compatibility with old options
  old_names = {"MMseqs2 (UniRef+Environmental)":"mmseqs2_uniref_env",
               "MMseqs2 (UniRef only)":"mmseqs2_uniref",
               "unpaired+paired":"unpaired_paired"}

  msa_mode  = old_names.get(msa_mode,msa_mode)
  pair_mode = old_names.get(pair_mode,pair_mode)
  use_dropout           = kwargs.pop("training", use_dropout)
  use_cluster_profile   = kwargs.pop("use_cluster_profile", None)
  use_fuse              = kwargs.pop("use_fuse", True)
  use_bfloat16          = kwargs.pop("use_bfloat16", True)
  max_msa               = kwargs.pop("max_msa",None)

  if max_msa is not None:
    max_seq, max_extra_seq = [int(x) for x in max_msa.split(":")]

  if len(kwargs) > 0:
    logger.warning(f"WARNING: the following options are not being used: {kwargs}")

  if use_templates and is_complex and "ptm" in model_type:
    logger.warning("WARNING: templates are not supported for alphafold2_ptm + complex.")
    use_templates = False

  # decide how to rank outputs
  if rank_by == "auto":
    rank_by = "multimer" if is_complex else "plddt"

  # Record the parameters of this run
  config = {
    "num_queries": len(queries),
    "use_templates": use_templates,
    "num_relax": num_relax,
    "msa_mode": msa_mode,
    "model_type": model_type,
    "num_models": num_models,
    "num_recycles": num_recycles,
    "recycle_early_stop_tolerance": recycle_early_stop_tolerance,
    "num_ensemble": num_ensemble,
    "model_order": model_order,
    "keep_existing_results": keep_existing_results,
    "rank_by": rank_by,
    "max_seq": max_seq,
    "max_extra_seq": max_extra_seq,
    "pair_mode": pair_mode,
    "host_url": host_url,
    "stop_at_score": stop_at_score,
    "random_seed": random_seed,
    "num_seeds": num_seeds,
    "recompile_padding": recompile_padding,
    "commit": get_commit(),
    "use_dropout": use_dropout,
    "use_cluster_profile": use_cluster_profile,
    "use_fuse": use_fuse,
    "use_bfloat16":use_bfloat16,
    "version": importlib_metadata.version("colabfold"),
  }
  config_out_file = result_dir.joinpath("config.json")
  config_out_file.write_text(json.dumps(config, indent=4))

  bibtex_file = write_bibtex(
    model=model_type,
    use_msa="mmseqs2" in msa_mode, 
    use_env="env" in msa_mode, 
    use_templates=use_templates,
    use_amber=num_relax > 0,
    result_dir=result_dir,
  )

  if custom_template_path is not None:
    mk_hhsearch_db(custom_template_path)

  # get max length (for padding purposes)
  max_len = 0
  for _, query_sequence, _ in queries:
    L = len("".join(query_sequence))
    if L > max_len: max_len = L

  pad_len = 0
  ranks, metrics = [],[]
  first_job = True
  for job_number, (raw_jobname, query_sequence, a3m_lines) in enumerate(queries):
    jobname = safe_filename(raw_jobname)
    
    #######################################
    # check if job has already finished
    #######################################
    # In the colab version and with --zip we know we're done when a zip file has been written
    result_zip = result_dir.joinpath(jobname).with_suffix(".result.zip")
    if keep_existing_results and result_zip.is_file():
      logger.info(f"Skipping {jobname} (result.zip)")
      continue
    # In the local version we use a marker file
    is_done_marker = result_dir.joinpath(jobname + ".done.txt")
    if keep_existing_results and is_done_marker.is_file():
      logger.info(f"Skipping {jobname} (already done)")
      continue

    total_len = len("".join(query_sequence))
    logger.info(f"Query {job_number + 1}/{len(queries)}: {jobname} (length {total_len})")

    ###########################################
    # generate MSA (a3m_lines) and templates
    ###########################################
    try:
      if use_templates or a3m_lines is None:
        (unpaired_msa, paired_msa, query_seqs_unique, query_seqs_cardinality, template_features) \
        = get_msa_and_templates(jobname, query_sequence, result_dir, msa_mode, use_templates, 
          custom_template_path, pair_mode, host_url)
      if a3m_lines is not None:
        (unpaired_msa, paired_msa, query_seqs_unique, query_seqs_cardinality, template_features_) \
        = unserialize_msa(a3m_lines, query_sequence)
        if not use_templates: template_features = template_features_

      # save a3m
      msa = msa_to_str(unpaired_msa, paired_msa, query_seqs_unique, query_seqs_cardinality)
      result_dir.joinpath(f"{jobname}.a3m").write_text(msa)
        
    except Exception as e:
      logger.exception(f"Could not get MSA/templates for {jobname}: {e}")
      continue
    
    #######################
    # generate features
    #######################
    try:
      (feature_dict, domain_names) \
      = generate_input_feature(query_seqs_unique, query_seqs_cardinality, unpaired_msa, paired_msa,
        template_features, is_complex, model_type)
      
      # to allow display of MSA info during colab/chimera run (thanks tomgoddard)
      if inputs_callback is not None:
        inputs_callback({"feature_dict":feature_dict})
    
    except Exception as e:
      logger.exception(f"Could not generate input features {jobname}: {e}")
      continue
    
    ######################
    # predict structures
    ######################
    try:
      # get list of lengths
      query_sequence_len_array = sum([[len(x)] * y 
        for x,y in zip(query_seqs_unique,query_seqs_cardinality)],[])
      
      # decide how much to pad (to avoid recompiling)
      if total_len > pad_len:
        if isinstance(recompile_padding, float):
          pad_len = math.ceil(total_len * recompile_padding)
        else:
          pad_len = total_len + recompile_padding
        pad_len = min(pad_len, max_len)
        logger.info(f"Padding length to {pad_len}")

      # prep model and params
      if first_job:
        # if one job input adjust max settings
        if len(queries) == 1 or msa_mode == "single_sequence":
          # get number of sequences
          if msa_mode == "single_sequence":
            num_seqs = 1
            if "ptm" in model_type and is_complex:
              num_seqs += len(query_sequence_len_array)
          else:
            if "msa_mask" in feature_dict:
              num_seqs = int(sum(feature_dict["msa_mask"].max(-1) == 1))
            else:
              num_seqs = int(len(feature_dict["msa"]))

          # get max settings
          # 512 5120 = alphafold (models 1,3,4)
          # 512 1024 = alphafold (models 2,5)
          # 508 2048 = alphafold-multimer (v3, models 1,2,3)
          # 508 1152 = alphafold-multimer (v3, models 4,5)
          # 252 1152 = alphafold-multimer (v1, v2)
          set_if = lambda x,y: y if x is None else x
          if model_type in ["alphafold2_multimer_v1","alphafold2_multimer_v2"]:
            (max_seq, max_extra_seq) = (set_if(max_seq,252), set_if(max_extra_seq,1152))
          elif model_type == "alphafold2_multimer_v3":
            (max_seq, max_extra_seq) = (set_if(max_seq,508), set_if(max_extra_seq,2048))
          else:
            (max_seq, max_extra_seq) = (set_if(max_seq,512), set_if(max_extra_seq,5120))
            if use_templates: num_seqs = num_seqs + 4
          
          # adjust max settings
          max_seq = min(num_seqs, max_seq)
          max_extra_seq = max(min(num_seqs - max_seq, max_extra_seq), 1)
          logger.info(f"Setting max_seq={max_seq}, max_extra_seq={max_extra_seq}")

        model_runner_and_params = load_models_and_params(
          num_models=num_models,
          use_templates=use_templates,
          num_recycles=num_recycles,
          num_ensemble=num_ensemble,
          model_suffix=model_suffix,
          data_dir=data_dir,
          stop_at_score=stop_at_score,
          rank_by=rank_by,
          use_dropout=use_dropout,
          max_seq=max_seq,
          max_extra_seq=max_extra_seq,
          use_cluster_profile=use_cluster_profile,
          recycle_early_stop_tolerance=recycle_early_stop_tolerance,
          use_fuse=use_fuse,
          use_bfloat16=use_bfloat16,
        )
        first_job = False

      results = predict_structure(
        prefix=jobname,
        result_dir=result_dir,
        feature_dict=feature_dict,
        is_complex=is_complex,
        use_templates=use_templates,
        sequences_lengths=query_sequence_len_array,
        pad_len=pad_len,
        model_type=model_type,
        model_runner_and_params=model_runner_and_params,
        model_order=model_order,
        num_relax=num_relax,
        rank_by=rank_by,
        stop_at_score=stop_at_score,
        outputs_callback=outputs_callback,
        use_gpu_relax=use_gpu_relax,
        random_seed=random_seed,
        num_seeds=num_seeds,
        save_all=save_all,
        save_single_representations=save_single_representations,
        save_pair_representations=save_pair_representations,
        save_recycles=save_recycles,
      )
      result_files = results["result_files"]
      ranks.append(results["rank"])
      metrics.append(results["metric"])

    except RuntimeError as e:
      # This normally happens on OOM. TODO: Filter for the specific OOM error message
      logger.error(f"Could not predict {jobname}. Not Enough GPU memory? {e}")
      continue

    ###############
    # save plots
    ###############

    # make msa plot
    msa_plot = plot_msa(feature_dict, dpi=dpi)
    coverage_png = result_dir.joinpath(f"{jobname}_coverage.png")
    msa_plot.savefig(str(coverage_png), bbox_inches='tight')
    msa_plot.close()
    result_files.append(coverage_png)

    # load the scores
    scores = []
    for r in results["rank"][:5]:
      scores_file = result_dir.joinpath(f"{jobname}_scores_{r}.json")
      with scores_file.open("r") as handle:
        scores.append(json.load(handle))
    
    # write alphafold-db format (pAE)
    af_pae_file = result_dir.joinpath(f"{jobname}_predicted_aligned_error_v1.json")
    af_pae_file.write_text(json.dumps({
      "predicted_aligned_error":scores[0]["pae"],
      "max_predicted_aligned_error":scores[0]["max_pae"]}))
    result_files.append(af_pae_file)
    
    # make pAE plots
    paes_plot = plot_paes([np.asarray(x["pae"]) for x in scores],
      Ls=query_sequence_len_array, dpi=dpi)
    pae_png = result_dir.joinpath(f"{jobname}_pae.png")
    paes_plot.savefig(str(pae_png), bbox_inches='tight')
    paes_plot.close()
    result_files.append(pae_png)

    # make pLDDT plot
    plddt_plot = plot_plddts([np.asarray(x["plddt"]) for x in scores],
      Ls=query_sequence_len_array, dpi=dpi)
    plddt_png = result_dir.joinpath(f"{jobname}_plddt.png")
    plddt_plot.savefig(str(plddt_png), bbox_inches='tight')
    plddt_plot.close()
    result_files.append(plddt_png)

    if use_templates:
      templates_file = result_dir.joinpath(f"{jobname}_template_domain_names.json")
      templates_file.write_text(json.dumps(domain_names))
      result_files.append(templates_file)

    result_files.append(result_dir.joinpath(jobname + ".a3m"))
    result_files += [bibtex_file, config_out_file]

    if zip_results:
      with zipfile.ZipFile(result_zip, "w") as result_zip:
        for file in result_files:
          result_zip.write(file, arcname=file.name)
      
      # Delete only after the zip was successful, and also not the bibtex and config because we need those again
      for file in result_files[:-2]:
        file.unlink()
    else:
      is_done_marker.touch()

  logger.info("Done")
  return {"rank":ranks,"metric":metrics}

def set_model_type(is_complex: bool, model_type: str) -> str:
  # backward-compatibility with old options
  old_names = {"AlphaFold2-multimer-v1":"alphafold2_multimer_v1",
               "AlphaFold2-multimer-v2":"alphafold2_multimer_v2",
               "AlphaFold2-multimer-v3":"alphafold2_multimer_v3",
               "AlphaFold2-ptm":"alphafold2_ptm"}
  model_type = old_names.get(model_type, model_type)
  if model_type == "auto" and is_complex:
    model_type = "alphafold2_multimer_v3"
  elif model_type == "auto" and not is_complex:
    model_type = "alphafold2_ptm"
  return model_type

def main():
  parser = argparse.ArgumentParser()
  parser.add_argument("input",
    default="input",
    help="Can be one of the following: "
    "Directory with fasta/a3m files, a csv/tsv file, a fasta file or an a3m file",
  )
  parser.add_argument("results", help="Directory to write the results to")

  # Main performance parameter
  parser.add_argument("--stop-at-score",
    help="Compute models until plddt (single chain) or ptmscore (complex) > threshold is reached. "
    "This can make colabfold much faster by only running the first model for easy queries.",
    type=float,
    default=100,
  )

  parser.add_argument("--num-recycles",
    help="Number of prediction recycles."
    "Increasing recycles can improve the quality but slows down the prediction.",
    type=int,
    default=None,
  )
  parser.add_argument("--recycle-early-stop-tolerance",
    help="Specify convergence criteria."
    "Run until the distance between recycles is within specified value.",
    type=float,
    default=None,
  )

  parser.add_argument("--num-ensemble",
    help="Number of ensembles."
    "The trunk of the network is run multiple times with different random choices for the MSA cluster centers.",
    type=int,
    default=1,
  )
  parser.add_argument("--num-seeds",
    help="Number of seeds to try. Will iterate from range(random_seed, random_seed+num_seeds)."
    ".",
    type=int,
    default=1,
  )
  parser.add_argument("--random-seed",
    help="Changing the seed for the random number generator can result in different structure predictions.",
    type=int,
    default=0,
  )
  parser.add_argument("--num-models", type=int, default=5, choices=[1, 2, 3, 4, 5])
  parser.add_argument("--recompile-padding",
    type=int,
    default=10,
    help="Whenever the input length changes, the model needs to be recompiled."
    "We pad sequences by specified length, so we can e.g. compute sequence from length 100 to 110 without recompiling."
    "The prediction will become marginally slower for the longer input, "
    "but overall performance increases due to not recompiling. "
    "Set to 0 to disable.",
  )
  parser.add_argument("--model-order", default="1,2,3,4,5", type=str)
  parser.add_argument("--host-url", default=DEFAULT_API_SERVER)
  parser.add_argument("--data")
  parser.add_argument("--msa-mode",
    default="mmseqs2_uniref_env",
    choices=[
      "mmseqs2_uniref_env",
      "mmseqs2_uniref",
      "single_sequence",
    ],
    help="Using an a3m file as input overwrites this option",
  )
  parser.add_argument("--model-type",
    help="predict strucutre/complex using the following model."
    'Auto will pick "alphafold2_ptm" for structure predictions and "alphafold2_multimer_v3" for complexes.',
    type=str,
    default="auto",
    choices=[
      "auto",
      "alphafold2_ptm",
      "alphafold2_multimer_v1",
      "alphafold2_multimer_v2",
      "alphafold2_multimer_v3",
    ],
  )
  parser.add_argument("--num-relax",
    help="specify how many of the top ranked structures to relax using amber.",
    type=int,
    default=0,
  )
  parser.add_argument("--use-templates",
    default=False,
    action="store_true",
    help="Use templates from pdb"
  )
  parser.add_argument("--custom-template-path",
    type=str,
    default=None,
    help="Directory with pdb files to be used as input",
  )
  parser.add_argument("--rank_by",
    help="rank models by auto, plddt or ptmscore",
    type=str,
    default="auto",
    choices=["auto", "plddt", "ptm", "iptm", "multimer"],
  )
  parser.add_argument("--pair-mode",
    help="rank models by auto, unpaired, paired, unpaired_paired",
    type=str,
    default="unpaired_paired",
    choices=["unpaired", "paired", "unpaired_paired"],
  )
  parser.add_argument("--sort-queries-by",
    help="sort queries by: none, length, random",
    type=str,
    default="length",
    choices=["none", "length", "random"],
  )
  parser.add_argument("--save-single-representations",
    default=False,
    action="store_true",
    help="saves the single representation embeddings of all models",
  )
  parser.add_argument("--save-pair-representations",
    default=False,
    action="store_true",
    help="saves the pair representation embeddings of all models",
  )
  parser.add_argument(
    "--use-dropout",
    default=False,
    action="store_true",
    help="activate dropouts during inference to sample from uncertainity of the models",
  )
  parser.add_argument("--max-seq",
    help="number of sequence clusters to use",
    type=int,
    default=None,
  )
  parser.add_argument("--max-extra-seq",
    help="number of extra sequences to use",
    type=int,
    default=None,
  )
  parser.add_argument("--zip-results",
    default=False,
    action="store_true",
    help="zip all results into one <jobname>.result.zip and delete the original files",
  )
  parser.add_argument("--use-gpu-relax",
    default=False,
    action="store_true",
    help="run amber on GPU instead of CPU",
  )
  parser.add_argument("--save-all",
    default=False,
    action="store_true",
    help="save ALL raw outputs from model to a pickle file",
  )
  parser.add_argument("--save-recycles",
    default=False,
    action="store_true",
    help="save all intermediate predictions at each recycle",
  )
  
  # undocumented arguements
  parser.add_argument("--overwrite-existing-results", default=False, action="store_true")
  parser.add_argument("--interaction-scan",           default=False, action="store_true")
  parser.add_argument("--disable-cluster-profile",    default=False, action="store_true")

  # backward compatability
  parser.add_argument('--templates',   default=False, action="store_true", help=argparse.SUPPRESS)
  parser.add_argument('--zip',         default=False, action="store_true", help=argparse.SUPPRESS)
  parser.add_argument('--amber',       default=False, action="store_true", help=argparse.SUPPRESS)
  parser.add_argument('--num-recycle', default=None,  type=int,            help=argparse.SUPPRESS)
  parser.add_argument("--max-msa",     default=None,  type=str,            help=argparse.SUPPRESS)

  # parse arguments
  args = parser.parse_args()

  # backward compatability
  if args.templates: args.use_templates = True
  if args.zip:       args.zip_results = True
  if args.amber and args.num_relax == 0: args.num_relax = args.num_models * args.num_seeds
  if args.num_recycle is not None: args.num_recycles = args.num_recycle
  if args.max_msa is not None: (args.max_seq, args.max_extra_seq) = (int(x) for x in args.max_msa.split(":"))
  
  # setup logging
  setup_logging(Path(args.results).joinpath("log.txt"))
  version = importlib_metadata.version("colabfold")
  commit = get_commit()
  if commit: version += f" ({commit})"

  logger.info(f"Running colabfold {version}")
  data_dir = Path(args.data or default_data_dir)
  model_order = [int(i) for i in args.model_order.split(",")]
  assert args.recompile_padding >= 0, "Can't apply negative padding"

  # parse queries
  if args.interaction_scan:
    # protocol from @Dohyun-s
    batch_size = 10
    queries, is_complex, headers = get_queries_pairwise(args.input, batch_size)
  else:
    queries, is_complex = get_queries(args.input, args.sort_queries_by)

  # download params
  model_type = set_model_type(is_complex, args.model_type)
  download_alphafold_params(model_type, data_dir)
  
  # warning about api
  if "mmseqs2" in args.msa_mode or args.use_templates:
    # TODO: check if server used in the case of templates
    uses_api = any((query[2] is None for query in queries))
    if uses_api and args.host_url == DEFAULT_API_SERVER:
      print(ACCEPT_DEFAULT_TERMS, file=sys.stderr)

  run_params = dict(
    result_dir=args.results,
    use_templates=args.use_templates,
    custom_template_path=args.custom_template_path,
    num_relax=args.num_relax,
    msa_mode=args.msa_mode,
    model_type=model_type,
    num_models=args.num_models,
    num_recycles=args.num_recycles,
    recycle_early_stop_tolerance=args.recycle_early_stop_tolerance,
    num_ensemble=args.num_ensemble,
    model_order=model_order,
    is_complex=is_complex,
    keep_existing_results=not args.overwrite_existing_results,
    rank_by=args.rank_by,
    pair_mode=args.pair_mode,
    data_dir=data_dir,
    host_url=args.host_url,
    random_seed=args.random_seed,
    num_seeds=args.num_seeds,
    stop_at_score=args.stop_at_score,
    recompile_padding=args.recompile_padding,
    zip_results=args.zip_results,
    save_single_representations=args.save_single_representations,
    save_pair_representations=args.save_pair_representations,
    use_dropout=args.use_dropout,
    max_seq=args.max_seq,
    max_extra_seq=args.max_extra_seq,
    use_cluster_profile=not args.disable_cluster_profile,
    use_gpu_relax = args.use_gpu_relax,
    save_all=args.save_all,
    save_recycles=args.save_recycles,
  )

  if args.interaction_scan:

    # protocol from @Dohyun-s
    from colabfold.mmseqs.api import run_mmseqs2
    output = [queries[i:i + batch_size] for i in range(0, len(queries), batch_size)]
    headers_list = [headers[i:i + batch_size] for i in range(0, len(headers), batch_size)]
    headers_list[0].remove(headers_list[0][0])
    header_first = headers[0]
    
    for jobname, batch in enumerate(output):
      query_seqs_unique = []
      for x in batch:
        if x not in query_seqs_unique:
          query_seqs_unique.append(x)
      use_env = "env" in args.msa_mode or "Environmental" in args.msa_mode
      paired_a3m_lines = run_mmseqs2(
        query_seqs_unique,
        str(Path(args.results).joinpath(str(jobname))),
        use_env=use_env,
        use_pairwise=True,
        use_pairing=True,
        host_url=args.host_url,
      )
      
      path_o = Path(args.results).joinpath(f"{jobname}_pairwise")      
      for filenum in path_o.iterdir():
        queries_new = [] 
        if Path(filenum).suffix.lower() == ".a3m":
          outdir = path_o.joinpath("tmp")
          unpack_a3ms(filenum, outdir)
          for i, file in enumerate(sorted(outdir.iterdir())):
            if outdir.joinpath(file).suffix.lower() == ".a3m":
              (seqs, header) = parse_fasta(Path(file).read_text())
              query_sequence = seqs[0]
              a3m_lines = [Path(file).read_text()]
              val = int(header[0].split('\t')[1][1:]) - 102
              queries_new.append((header_first + '_' + headers_list[jobname][val], query_sequence, a3m_lines))

          if args.sort_queries_by == "length":
            queries_new.sort(key=lambda t: len(''.join(t[1])),reverse=True)
          elif args.sort_queries_by == "random":
            random.shuffle(queries_new)

          run(queries=queries_new, **run_params)
  
  else:
    run(queries=queries, **run_params)
