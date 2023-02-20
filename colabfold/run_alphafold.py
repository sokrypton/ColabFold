from __future__ import annotations

import json
import math
import zipfile
import importlib_metadata
import numpy as np
from typing import Any, Callable, Dict, List, Optional, Tuple, Union, TYPE_CHECKING
from pathlib import Path

import logging
logger = logging.getLogger(__name__)

import jax
import jax.numpy as jnp
logging.getLogger('jax._src.lib.xla_bridge').addFilter(lambda _: False)

# import from colabfold
from colabfold.inputs import (
  mk_hhsearch_db, generate_input_feature, msa_to_str,
  unserialize_msa, 
  get_msa_and_templates,
)

from colabfold.predict import predict_structure
from colabfold.download import default_data_dir

from colabfold.utils import (
  DEFAULT_API_SERVER,
  get_commit,
  safe_filename)

from colabfold.citations import write_bibtex

def set_model_type(is_complex: bool, model_type: str) -> str:
  # backward-compatibility with old options
  old_names = {"AlphaFold2-multimer-v1":"alphafold2_multimer_v1",
               "AlphaFold2-multimer-v2":"alphafold2_multimer_v2",
               "AlphaFold2-multimer-v3":"alphafold2_multimer_v3",
               "AlphaFold2-ptm"        :"alphafold2_ptm",
               "AlphaFold2"            :"alphafold2"}
  model_type = old_names.get(model_type, model_type)
  if model_type == "auto" and is_complex:
    model_type = "alphafold2_multimer_v3"
  elif model_type == "auto" and not is_complex:
    model_type = "alphafold2_ptm"
  return model_type

def run(
  queries: List[Tuple[str, Union[str, List[str]], Optional[List[str]]]],
  result_dir: Union[str, Path],
  is_complex: bool,
  num_recycles: Optional[int] = None,
  recycle_early_stop_tolerance: Optional[float] = None,
  model_order: List[int] = [1,2,3,4,5],
  num_models: int = 5,
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
      tf.get_logger().setLevel(logging.ERROR)
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
  elif model_type == "alphafold2_ptm":         model_suffix = "_ptm"
  elif model_type == "alphafold2":             model_suffix = ""
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

  if use_templates and is_complex and "multimer" not in model_type:
    logger.warning("WARNING: templates are not supported for alphafold2_ptm + complex.")
    use_templates = False

  # decide how to rank outputs
  if rank_by == "auto":
    rank_by = "multimer" if is_complex else "plddt"
  if "ptm" not in model_type and "multimer" not in model_type:
    rank_by = "plddt"

  # get max length
  max_len, max_num = 0, 0
  for _, query_sequence, _ in queries:
    N = 1 if isinstance(query_sequence,str) else len(query_sequence)
    L = len("".join(query_sequence))
    if L > max_len: max_len = L
    if N > max_num: max_num = N

  # get max sequences
  # 512 5120 = alphafold_ptm (models 1,3,4)
  # 512 1024 = alphafold_ptm (models 2,5)
  # 508 2048 = alphafold-multimer_v3 (models 1,2,3)
  # 508 1152 = alphafold-multimer_v3 (models 4,5)
  # 252 1152 = alphafold-multimer_v[1,2]

  set_if = lambda x,y: y if x is None else x
  if model_type in ["alphafold2_multimer_v1","alphafold2_multimer_v2"]:
    (max_seq, max_extra_seq) = (set_if(max_seq,252), set_if(max_extra_seq,1152))
  elif model_type == "alphafold2_multimer_v3":
    (max_seq, max_extra_seq) = (set_if(max_seq,508), set_if(max_extra_seq,2048))
  else:
    (max_seq, max_extra_seq) = (set_if(max_seq,512), set_if(max_extra_seq,5120))

  if msa_mode == "single_sequence":
    num_seqs = 1
    if is_complex and "multimer" not in model_type: num_seqs += max_num
    if use_templates: num_seqs += 4
    max_seq = min(num_seqs, max_seq)
    max_extra_seq = max(min(num_seqs - max_seq, max_extra_seq), 1)

  # Record the parameters of this run
  config = {
    "num_queries": len(queries),
    "use_templates": use_templates,
    "num_relax": num_relax,
    "msa_mode": msa_mode,
    "model_type": model_type,
    "num_recycles": num_recycles,
    "recycle_early_stop_tolerance": recycle_early_stop_tolerance,
    "num_ensemble": num_ensemble,
    "model_order": model_order,
    "num_models": num_models,
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

    seq_len = len("".join(query_sequence))
    logger.info(f"Query {job_number + 1}/{len(queries)}: {jobname} (length {seq_len})")

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

      if seq_len != sum(query_sequence_len_array):
        logger.info(f"ERROR: please report bug, something went wrong...")
        seq_len = sum(query_sequence_len_array)
      
      # decide how much to pad (to avoid recompiling)
      if seq_len > pad_len:
        if isinstance(recompile_padding, float):
          pad_len = math.ceil(seq_len * recompile_padding)
        else:
          pad_len = seq_len + recompile_padding
        pad_len = min(pad_len, max_len)

      # prep model and params
      if first_job:
        # if one job input adjust max settings
        if len(queries) == 1 and msa_mode != "single_sequence":
          if "msa_mask" in feature_dict:
            num_seqs = int(sum(feature_dict["msa_mask"].max(-1) == 1))
          else:
            num_seqs = int(len(feature_dict["msa"]))
          if use_templates: num_seqs += 4
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
          save_all=save_all,
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
    if "pae" in scores[0]:
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