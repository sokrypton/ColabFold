#########################################################
# code for generating input features to alphafold
#########################################################
from typing import Any, Callable, Dict, List, Optional, Tuple, Union, TYPE_CHECKING
from pathlib import Path
import pandas
import random
from io import StringIO
import shutil
import numpy as np
import os

import logging
logger = logging.getLogger(__name__)

import jax
import jax.numpy as jnp
import haiku
logging.getLogger('jax._src.lib.xla_bridge').addFilter(lambda _: False)

# TODO: move alphafold specific functions to alphafold/inputs.py

# alphafold imports
from alphafold.model import model
from alphafold.common import protein, residue_constants
from alphafold.data import (
  feature_processing,
  msa_pairing,
  pipeline,
  pipeline_multimer,
  templates,
)
from alphafold.data.tools import hhsearch

# colabfold imports
from colabfold.alphafold.msa import make_fixed_size, make_fixed_size_multimer
from colabfold.utils import (
  DEFAULT_API_SERVER,
  NO_GPU_FOUND,
  CIF_REVISION_DATE,
  CFMMCIFIO,
)
from Bio.PDB import MMCIFParser, PDBParser, MMCIF2Dict

###############################
# INPUTS
###############################
def pad_input(
  input_features: model.features.FeatureDict,
  model_runner: model.RunModel,
  model_name: str,
  pad_len: int,
  use_templates: bool,
) -> model.features.FeatureDict:

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

def pad_input_multimer(
  input_features: model.features.FeatureDict,
  model_runner: model.RunModel,
  model_name: str,
  pad_len: int,
  msa_cluster_size: Optional[int],
  use_templates: bool,
) -> model.features.FeatureDict:
  shape_schema = {
    "aatype": ["num residues placeholder"],
    "residue_index": ["num residues placeholder"],
    "msa": ["msa placeholder", "num residues placeholder"],
    "template_all_atom_positions": [
      "num templates placeholder",
      "num residues placeholder",
      None,
      None,
    ],
    "template_all_atom_mask": [
      "num templates placeholder",
      "num residues placeholder",
      None,
    ],
    "template_aatype": ["num templates placeholder", "num residues placeholder"],
    "asym_id": ["num residues placeholder"],
    "sym_id": ["num residues placeholder"],
    "entity_id": ["num residues placeholder"],
    "deletion_matrix": ["msa placeholder", "num residues placeholder"],
    "deletion_mean": ["num residues placeholder"],
    "all_atom_mask": ["num residues placeholder", None],
    "all_atom_positions": ["num residues placeholder", None, None],
    "entity_mask": ["num residues placeholder"],
    "cluster_bias_mask": ["msa placeholder"],
    "bert_mask": ["msa placeholder", "num residues placeholder"],
    "seq_mask": ["num residues placeholder"],
    "msa_mask": ["msa placeholder", "num residues placeholder"],
    "seq_length": [None],
    "num_alignments": [None],
    "assembly_num_chains": [None],
    "num_templates": [None],
  }
  input_fix = make_fixed_size_multimer(
    input_features,
    shape_schema,
    num_res=pad_len,
    msa_cluster_size=msa_cluster_size,
    num_templates=4,
  )
  return input_fix

def parse_fasta(fasta_string: str) -> Tuple[List[str], List[str]]:
  """Parses FASTA string and returns list of strings with amino-acid sequences.

  Arguments:
    fasta_string: The string contents of a FASTA file.

  Returns:
    A tuple of two lists:
    * A list of sequences.
    * A list of sequence descriptions taken from the comment lines. In the
    same order as the sequences.
  """
  sequences = []
  descriptions = []
  index = -1
  for line in fasta_string.splitlines():
    line = line.strip()
    if line.startswith("#"):
      continue
    if line.startswith(">"):
      index += 1
      descriptions.append(line[1:])  # Remove the '>' at the beginning.
      sequences.append("")
      continue
    elif not line:
      continue  # Skip blank lines.
    sequences[index] += line

  return sequences, descriptions

def get_queries(
  input_path: Union[str, Path], sort_queries_by: str = "length"
) -> Tuple[List[Tuple[str, str, Optional[List[str]]]], bool]:
  """Reads a directory of fasta files, a single fasta file or a csv file and returns a tuple
  of job name, sequence and the optional a3m lines"""

  input_path = Path(input_path)
  if not input_path.exists():
    raise OSError(f"{input_path} could not be found")

  if input_path.is_file():
    if input_path.suffix == ".csv" or input_path.suffix == ".tsv":
      sep = "\t" if input_path.suffix == ".tsv" else ","
      df = pandas.read_csv(input_path, sep=sep)
      assert "id" in df.columns and "sequence" in df.columns
      queries = [
        (seq_id, sequence.upper().split(":"), None)
        for seq_id, sequence in df[["id", "sequence"]].itertuples(index=False)
      ]
      for i in range(len(queries)):
        if len(queries[i][1]) == 1:
          queries[i] = (queries[i][0], queries[i][1][0], None)
    elif input_path.suffix == ".a3m":
      (seqs, header) = parse_fasta(input_path.read_text())
      if len(seqs) == 0:
        raise ValueError(f"{input_path} is empty")
      query_sequence = seqs[0]
      # Use a list so we can easily extend this to multiple msas later
      a3m_lines = [input_path.read_text()]
      queries = [(input_path.stem, query_sequence, a3m_lines)]
    elif input_path.suffix in [".fasta", ".faa", ".fa"]:
      (sequences, headers) = parse_fasta(input_path.read_text())
      queries = []
      for sequence, header in zip(sequences, headers):
        sequence = sequence.upper()
        if sequence.count(":") == 0:
          # Single sequence
          queries.append((header, sequence, None))
        else:
          # Complex mode
          queries.append((header, sequence.upper().split(":"), None))
    else:
      raise ValueError(f"Unknown file format {input_path.suffix}")
  else:
    assert input_path.is_dir(), "Expected either an input file or a input directory"
    queries = []
    for file in sorted(input_path.iterdir()):
      if not file.is_file():
        continue
      if file.suffix.lower() not in [".a3m", ".fasta", ".faa"]:
        logger.warning(f"non-fasta/a3m file in input directory: {file}")
        continue
      (seqs, header) = parse_fasta(file.read_text())
      if len(seqs) == 0:
        logger.error(f"{file} is empty")
        continue
      query_sequence = seqs[0]
      if len(seqs) > 1 and file.suffix in [".fasta", ".faa", ".fa"]:
        logger.warning(
          f"More than one sequence in {file}, ignoring all but the first sequence"
        )

      if file.suffix.lower() == ".a3m":
        a3m_lines = [file.read_text()]
        queries.append((file.stem, query_sequence.upper(), a3m_lines))
      else:
        if query_sequence.count(":") == 0:
          # Single sequence
          queries.append((file.stem, query_sequence, None))
        else:
          # Complex mode
          queries.append((file.stem, query_sequence.upper().split(":"), None))

  # sort by seq. len
  if sort_queries_by == "length":
    queries.sort(key=lambda t: len("".join(t[1])))
  elif sort_queries_by == "random":
    random.shuffle(queries)
  is_complex = False
  for job_number, (raw_jobname, query_sequence, a3m_lines) in enumerate(queries):
    if isinstance(query_sequence, list):
      is_complex = True
      break
    if a3m_lines is not None and a3m_lines[0].startswith("#"):
      a3m_line = a3m_lines[0].splitlines()[0]
      tab_sep_entries = a3m_line[1:].split("\t")
      if len(tab_sep_entries) == 2:
        query_seq_len = tab_sep_entries[0].split(",")
        query_seq_len = list(map(int, query_seq_len))
        query_seqs_cardinality = tab_sep_entries[1].split(",")
        query_seqs_cardinality = list(map(int, query_seqs_cardinality))
        is_single_protein = (
          True
          if len(query_seq_len) == 1 and query_seqs_cardinality[0] == 1
          else False
        )
        if not is_single_protein:
          is_complex = True
          break
  return queries, is_complex

def get_queries_pairwise(
  input_path: Union[str, Path], batch_size: int = 10,
) -> Tuple[List[Tuple[str, str, Optional[List[str]]]], bool]:
  """Reads a directory of fasta files, a single fasta file or a csv file and returns a tuple
  of job name, sequence and the optional a3m lines"""
  input_path = Path(input_path)
  if not input_path.exists():
    raise OSError(f"{input_path} could not be found")
  if input_path.is_file():
    if input_path.suffix == ".csv" or input_path.suffix == ".tsv":
      sep = "\t" if input_path.suffix == ".tsv" else ","
      df = pandas.read_csv(input_path, sep=sep)
      assert "id" in df.columns and "sequence" in df.columns
      queries = []
      seq_id_list = []
      for i, (seq_id, sequence) in enumerate(df[["id", "sequence"]].itertuples(index=False)):
        if i>0 and i % 10 == 0:
          queries.append(queries[0].upper())
        queries.append(sequence.upper())
        seq_id_list.append(seq_id)
      return queries, True, seq_id_list
    elif input_path.suffix == ".a3m":
      raise NotImplementedError()
    elif input_path.suffix in [".fasta", ".faa", ".fa"]:
      (sequences, headers) = parse_fasta(input_path.read_text())
      queries = []
      for i, (sequence, header) in enumerate(zip(sequences, headers)):
        sequence = sequence.upper()
        if sequence.count(":") == 0:
          if i>0 and i % 10 == 0:
            queries.append(sequences[0])
          queries.append(sequence)
        else:
          # Complex mode
          queries.append((header, sequence.upper().split(":"), None))
      return queries, True, headers
    else:
      raise ValueError(f"Unknown file format {input_path.suffix}")
  else:
    raise NotImplementedError()

def pair_sequences(
  a3m_lines: List[str], query_sequences: List[str], query_cardinality: List[int]
) -> str:
  a3m_line_paired = [""] * len(a3m_lines[0].splitlines())
  for n, seq in enumerate(query_sequences):
    lines = a3m_lines[n].splitlines()
    for i, line in enumerate(lines):
      if line.startswith(">"):
        if n != 0:
          line = line.replace(">", "\t", 1)
        a3m_line_paired[i] = a3m_line_paired[i] + line
      else:
        a3m_line_paired[i] = a3m_line_paired[i] + line * query_cardinality[n]
  return "\n".join(a3m_line_paired)

def pad_sequences(
  a3m_lines: List[str], query_sequences: List[str], query_cardinality: List[int]
) -> str:
  _blank_seq = [
    ("-" * len(seq))
    for n, seq in enumerate(query_sequences)
    for _ in range(query_cardinality[n])
  ]
  a3m_lines_combined = []
  pos = 0
  for n, seq in enumerate(query_sequences):
    for j in range(0, query_cardinality[n]):
      lines = a3m_lines[n].split("\n")
      for a3m_line in lines:
        if len(a3m_line) == 0:
          continue
        if a3m_line.startswith(">"):
          a3m_lines_combined.append(a3m_line)
        else:
          a3m_lines_combined.append(
            "".join(_blank_seq[:pos] + [a3m_line] + _blank_seq[pos + 1 :])
          )
      pos += 1
  return "\n".join(a3m_lines_combined)

def get_msa_and_templates(
  jobname: str,
  query_sequences: Union[str, List[str]],
  result_dir: Path,
  msa_mode: str,
  use_templates: bool,
  custom_template_path: str,
  pair_mode: str,
  host_url: str = DEFAULT_API_SERVER,
  use_pairwise: bool = False,
) -> Tuple[
  Optional[List[str]], Optional[List[str]], List[str], List[int], List[Dict[str, Any]]
]:
  from colabfold.mmseqs.api import run_mmseqs2

  use_env = msa_mode == "mmseqs2_uniref_env"
  if isinstance(query_sequences, str): query_sequences = [query_sequences]

  # remove duplicates before searching
  query_seqs_unique = []
  for x in query_sequences:
    if x not in query_seqs_unique:
      query_seqs_unique.append(x)

  # determine how many times is each sequence is used
  query_seqs_cardinality = [0] * len(query_seqs_unique)
  for seq in query_sequences:
    seq_idx = query_seqs_unique.index(seq)
    query_seqs_cardinality[seq_idx] += 1

  # get template features
  template_features = []
  if use_templates:
    a3m_lines_mmseqs2, template_paths = run_mmseqs2(
      query_seqs_unique,
      str(result_dir.joinpath(jobname)),
      use_env=use_env,
      use_templates=True,
      host_url=host_url,
    )
    if custom_template_path is not None:
      template_paths = {}
      for index in range(0, len(query_seqs_unique)):
        template_paths[index] = custom_template_path
    if template_paths is None:
      logger.info("No template detected")
      for index in range(0, len(query_seqs_unique)):
        template_feature = mk_mock_template(query_seqs_unique[index])
        template_features.append(template_feature)
    else:
      for index in range(0, len(query_seqs_unique)):
        if template_paths[index] is not None:
          template_feature = mk_template(
            a3m_lines_mmseqs2[index],
            template_paths[index],
            query_seqs_unique[index],
          )
          if len(template_feature["template_domain_names"]) == 0:
            template_feature = mk_mock_template(query_seqs_unique[index])
            logger.info(f"Sequence {index} found no templates")
          else:
            logger.info(
              f"Sequence {index} found templates: {template_feature['template_domain_names'].astype(str).tolist()}"
            )
        else:
          template_feature = mk_mock_template(query_seqs_unique[index])
          logger.info(f"Sequence {index} found no templates")

        template_features.append(template_feature)
  else:
    for index in range(0, len(query_seqs_unique)):
      template_feature = mk_mock_template(query_seqs_unique[index])
      template_features.append(template_feature)

  if len(query_sequences) == 1:
    pair_mode = "none"

  if pair_mode == "none" or pair_mode == "unpaired" or pair_mode == "unpaired_paired":
    if msa_mode == "single_sequence":
      a3m_lines = []
      num = 101
      for i, seq in enumerate(query_seqs_unique):
        a3m_lines.append(f">{num + i}\n{seq}")
    else:
      # find normal a3ms
      a3m_lines = run_mmseqs2(
        query_seqs_unique,
        str(result_dir.joinpath(jobname)),
        use_env=use_env,
        use_pairwise=use_pairwise,
        use_pairing=False,
        host_url=host_url,
      )
  else:
    a3m_lines = None

  if msa_mode != "single_sequence" and (
    pair_mode == "paired" or pair_mode == "unpaired_paired"
  ):
    # find paired a3m if not a homooligomers
    if len(query_seqs_unique) > 1:
      paired_a3m_lines = run_mmseqs2(
        query_seqs_unique,
        str(result_dir.joinpath(jobname)),
        use_env,
        use_pairing=True,
        host_url=host_url,
      )
    else:
      # homooligomers
      num = 101
      paired_a3m_lines = []
      for i in range(0, query_seqs_cardinality[0]):
        paired_a3m_lines.append(f">{num+i}\n{query_seqs_unique[0]}\n")
  else:
    paired_a3m_lines = None

  return (
    a3m_lines,
    paired_a3m_lines,
    query_seqs_unique,
    query_seqs_cardinality,
    template_features,
  )

def build_monomer_feature(
  sequence: str, unpaired_msa: str, template_features: Dict[str, Any]
):
  msa = pipeline.parsers.parse_a3m(unpaired_msa)
  # gather features
  return {
    **pipeline.make_sequence_features(
      sequence=sequence, description="none", num_res=len(sequence)
    ),
    **pipeline.make_msa_features([msa]),
    **template_features,
  }

def build_multimer_feature(paired_msa: str) -> Dict[str, np.ndarray]:
  parsed_paired_msa = pipeline.parsers.parse_a3m(paired_msa)
  return {
    f"{k}_all_seq": v
    for k, v in pipeline.make_msa_features([parsed_paired_msa]).items()
  }

def process_multimer_features(
  features_for_chain: Dict[str, Dict[str, np.ndarray]]
) -> Dict[str, np.ndarray]:
  all_chain_features = {}
  for chain_id, chain_features in features_for_chain.items():
    all_chain_features[chain_id] = pipeline_multimer.convert_monomer_features(
      chain_features, chain_id
    )

  all_chain_features = pipeline_multimer.add_assembly_features(all_chain_features)
  # np_example = feature_processing.pair_and_merge(
  #  all_chain_features=all_chain_features, is_prokaryote=is_prokaryote)
  feature_processing.process_unmerged_features(all_chain_features)
  np_chains_list = list(all_chain_features.values())
  # noinspection PyProtectedMember
  pair_msa_sequences = not feature_processing._is_homomer_or_monomer(np_chains_list)
  chains = list(np_chains_list)
  chain_keys = chains[0].keys()
  updated_chains = []
  for chain_num, chain in enumerate(chains):
    new_chain = {k: v for k, v in chain.items() if "_all_seq" not in k}
    for feature_name in chain_keys:
      if feature_name.endswith("_all_seq"):
        feats_padded = msa_pairing.pad_features(
          chain[feature_name], feature_name
        )
        new_chain[feature_name] = feats_padded
    new_chain["num_alignments_all_seq"] = np.asarray(
      len(np_chains_list[chain_num]["msa_all_seq"])
    )
    updated_chains.append(new_chain)
  np_chains_list = updated_chains
  np_chains_list = feature_processing.crop_chains(
    np_chains_list,
    msa_crop_size=feature_processing.MSA_CROP_SIZE,
    pair_msa_sequences=pair_msa_sequences,
    max_templates=feature_processing.MAX_TEMPLATES,
  )
  # merge_chain_features crashes if there are additional features only present in one chain
  # remove all features that are not present in all chains
  common_features = set([*np_chains_list[0]]).intersection(*np_chains_list)
  np_chains_list = [
    {key: value for (key, value) in chain.items() if key in common_features}
    for chain in np_chains_list
  ]
  np_example = feature_processing.msa_pairing.merge_chain_features(
    np_chains_list=np_chains_list,
    pair_msa_sequences=pair_msa_sequences,
    max_templates=feature_processing.MAX_TEMPLATES,
  )
  np_example = feature_processing.process_final(np_example)

  # Pad MSA to avoid zero-sized extra_msa.
  np_example = pipeline_multimer.pad_msa(np_example, min_num_seq=512)
  return np_example

def pair_msa(
  query_seqs_unique: List[str],
  query_seqs_cardinality: List[int],
  paired_msa: Optional[List[str]],
  unpaired_msa: Optional[List[str]],
) -> str:
  if paired_msa is None and unpaired_msa is not None:
    a3m_lines = pad_sequences(
      unpaired_msa, query_seqs_unique, query_seqs_cardinality
    )
  elif paired_msa is not None and unpaired_msa is not None:
    a3m_lines = (
      pair_sequences(paired_msa, query_seqs_unique, query_seqs_cardinality)
      + "\n"
      + pad_sequences(unpaired_msa, query_seqs_unique, query_seqs_cardinality)
    )
  elif paired_msa is not None and unpaired_msa is None:
    a3m_lines = pair_sequences(
      paired_msa, query_seqs_unique, query_seqs_cardinality
    )
  else:
    raise ValueError(f"Invalid pairing")
  return a3m_lines

def generate_input_feature(
  query_seqs_unique: List[str],
  query_seqs_cardinality: List[int],
  unpaired_msa: List[str],
  paired_msa: List[str],
  template_features: List[Dict[str, Any]],
  is_complex: bool,
  model_type: str,
) -> Tuple[Dict[str, Any], Dict[str, str]]:

  input_feature = {}
  domain_names = {}
  if is_complex and "multimer" not in model_type:

    full_sequence = ""
    Ls = []
    for sequence_index, sequence in enumerate(query_seqs_unique):
      for cardinality in range(0, query_seqs_cardinality[sequence_index]):
        full_sequence += sequence
        Ls.append(len(sequence))

    # bugfix
    a3m_lines = f">0\n{full_sequence}\n"
    a3m_lines += pair_msa(query_seqs_unique, query_seqs_cardinality, paired_msa, unpaired_msa)    

    input_feature = build_monomer_feature(full_sequence, a3m_lines, mk_mock_template(full_sequence))
    input_feature["residue_index"] = np.concatenate([np.arange(L) for L in Ls])
    input_feature["asym_id"] = np.concatenate([np.full(L,n) for n,L in enumerate(Ls)])
    if any(
      [
        template != b"none"
        for i in template_features
        for template in i["template_domain_names"]
      ]
    ):
      logger.warning(
        "alphafold2_ptm complex does not consider templates. Chose multimer model-type for template support."
      )

  else:
    features_for_chain = {}
    chain_cnt = 0
    # for each unique sequence
    for sequence_index, sequence in enumerate(query_seqs_unique):
      
      # get unpaired msa
      if unpaired_msa is None:
        input_msa = f">{101 + sequence_index}\n{sequence}"
      else:
        input_msa = unpaired_msa[sequence_index]

      feature_dict = build_monomer_feature(
        sequence, input_msa, template_features[sequence_index])

      if "multimer" in model_type:
        # get paired msa
        if paired_msa is None:
          input_msa = f">{101 + sequence_index}\n{sequence}"
        else:
          input_msa = paired_msa[sequence_index]
        feature_dict.update(build_multimer_feature(input_msa))

      # for each copy
      for cardinality in range(0, query_seqs_cardinality[sequence_index]):
        features_for_chain[protein.PDB_CHAIN_IDS[chain_cnt]] = feature_dict
        chain_cnt += 1

    if "multimer" not in model_type:
      input_feature = features_for_chain[protein.PDB_CHAIN_IDS[0]]
      input_feature["asym_id"] = np.zeros(input_feature["aatype"].shape[0],dtype=int)
      domain_names = {
        protein.PDB_CHAIN_IDS[0]: [
          name.decode("UTF-8")
          for name in input_feature["template_domain_names"]
          if name != b"none"
        ]
      }
    else:
      # combine features across all chains
      input_feature = process_multimer_features(features_for_chain)
      domain_names = {
        chain: [
          name.decode("UTF-8")
          for name in feature["template_domain_names"]
          if name != b"none"
        ]
        for (chain, feature) in features_for_chain.items()
      }
  return (input_feature, domain_names)

def generate_msa_size(inputs, query_seqs_unique, use_templates, is_complex, model_type):
  template_features_ = []
  from colabfold.inputs import mk_mock_template
  from colabfold.inputs import generate_input_feature
  for query_seq in query_seqs_unique:
    template_feature = mk_mock_template(query_seq)
    template_features_.append(template_feature)
  if not use_templates: template_features = template_features_
  else: raise NotImplementedError

  (feature_dict, _) \
  = generate_input_feature(*inputs, template_features, is_complex, model_type)
  return feature_dict["bert_mask"].shape[0]

def unserialize_msa(
  a3m_lines: List[str], query_sequence: Union[List[str], str]
) -> Tuple[
  Optional[List[str]],
  Optional[List[str]],
  List[str],
  List[int],
  List[Dict[str, Any]],
]:
  a3m_lines = a3m_lines[0].replace("\x00", "").splitlines()
  if not a3m_lines[0].startswith("#") or len(a3m_lines[0][1:].split("\t")) != 2:
    assert isinstance(query_sequence, str)
    return (
      ["\n".join(a3m_lines)],
      None,
      [query_sequence],
      [1],
      [mk_mock_template(query_sequence)],
    )

  if len(a3m_lines) < 3:
    raise ValueError(f"Unknown file format a3m")
  tab_sep_entries = a3m_lines[0][1:].split("\t")
  query_seq_len = tab_sep_entries[0].split(",")
  query_seq_len = list(map(int, query_seq_len))
  query_seqs_cardinality = tab_sep_entries[1].split(",")
  query_seqs_cardinality = list(map(int, query_seqs_cardinality))
  is_homooligomer = (
    True if len(query_seq_len) == 1 and query_seqs_cardinality[0] > 1 else False
  )
  is_single_protein = (
    True if len(query_seq_len) == 1 and query_seqs_cardinality[0] == 1 else False
  )
  query_seqs_unique = []
  prev_query_start = 0
  # we store the a3m with cardinality of 1
  for n, query_len in enumerate(query_seq_len):
    query_seqs_unique.append(
      a3m_lines[2][prev_query_start : prev_query_start + query_len]
    )
    prev_query_start += query_len
  paired_msa = [""] * len(query_seq_len)
  unpaired_msa = [""] * len(query_seq_len)
  already_in = dict()
  for i in range(1, len(a3m_lines), 2):
    header = a3m_lines[i]
    seq = a3m_lines[i + 1]
    if (header, seq) in already_in:
      continue
    already_in[(header, seq)] = 1
    has_amino_acid = [False] * len(query_seq_len)
    seqs_line = []
    prev_pos = 0
    for n, query_len in enumerate(query_seq_len):
      paired_seq = ""
      curr_seq_len = 0
      for pos in range(prev_pos, len(seq)):
        if curr_seq_len == query_len:
          prev_pos = pos
          break
        paired_seq += seq[pos]
        if seq[pos].islower():
          continue
        if seq[pos] != "-":
          has_amino_acid[n] = True
        curr_seq_len += 1
      seqs_line.append(paired_seq)

    # is sequence is paired add them to output
    if (
      not is_single_protein
      and not is_homooligomer
      and sum(has_amino_acid) == len(query_seq_len)
    ):
      header_no_faster = header.replace(">", "")
      header_no_faster_split = header_no_faster.split("\t")
      for j in range(0, len(seqs_line)):
        paired_msa[j] += ">" + header_no_faster_split[j] + "\n"
        paired_msa[j] += seqs_line[j] + "\n"
    else:
      for j, seq in enumerate(seqs_line):
        if has_amino_acid[j]:
          unpaired_msa[j] += header + "\n"
          unpaired_msa[j] += seq + "\n"
  if is_homooligomer:
    # homooligomers
    num = 101
    paired_msa = [""] * query_seqs_cardinality[0]
    for i in range(0, query_seqs_cardinality[0]):
      paired_msa[i] = ">" + str(num + i) + "\n" + query_seqs_unique[0] + "\n"
  if is_single_protein:
    paired_msa = None
  template_features = []
  for query_seq in query_seqs_unique:
    template_feature = mk_mock_template(query_seq)
    template_features.append(template_feature)

  if unpaired_msa == [""] * len(query_seq_len):
    unpaired_msa = None
  return (
    unpaired_msa,
    paired_msa,
    query_seqs_unique,
    query_seqs_cardinality,
    template_features,
  )

def unpack_a3ms(input, outdir):
  if not os.path.isdir(outdir):
     os.mkdir(outdir)
  i = 1
  count = 0
  with open(input, "r") as f:
    lines = f.readlines()
    a3m = ""
    hasLen1 = False
    len1 = 0
    inMsa2 = False
    hasLen2 = False
    switch = False
    len2 = 0
    left=[]
    right=[]
    descp_l=[]
    descp_r=[]
    for idx, line in enumerate(lines):
      if line.startswith("\x00"):
        i += 1
        line = line[1:]
        inMsa2 = True
        if i % 2 == 1:
          a3m="\n".join([l[:-1]+'\t'+r[:-1]+'\n'+a_[:-1]+b_[:-1]
           for l, r, a_, b_ in zip(descp_l, descp_r, left, right)])

          complex_description = "#" + str(len1) + "," + str(len2) + "\t1,1\n"
          a3m = complex_description + a3m
          out = Path(outdir) / f"job{count}.a3m"
          with open(out, "w") as w:
            w.write(a3m)
          a3m = ""
          count += 1
          hasLen1 = False
          hasLen2 = False
          inMsa2 = False
          left = []
          right = []
          descp_l = []
          descp_r = []
          switch = False
        else:
          switch = True

      if line.startswith(">") and hasLen2 is False:
        descp_l.append(line)
      if not line.startswith(">") and not hasLen1:
        len1 = len(line)-1
        hasLen1 = True
      if not line.startswith(">") and hasLen1 is True and hasLen2 is False:
        left.append(line)
      if line.startswith(">") and switch:
        descp_r.append(line)
      if not line.startswith(">") and not hasLen2 and inMsa2:
        len2 = len(line)-1
        hasLen2 = True
      if not line.startswith(">") and hasLen1 is True and hasLen2 is True:
        right.append(line)

def msa_to_str(
  unpaired_msa: List[str],
  paired_msa: List[str],
  query_seqs_unique: List[str],
  query_seqs_cardinality: List[int],
) -> str:
  msa = "#" + ",".join(map(str, map(len, query_seqs_unique))) + "\t"
  msa += ",".join(map(str, query_seqs_cardinality)) + "\n"
  # build msa with cardinality of 1, it makes it easier to parse and manipulate
  query_seqs_cardinality = [1 for _ in query_seqs_cardinality]
  msa += pair_msa(query_seqs_unique, query_seqs_cardinality, paired_msa, unpaired_msa)
  return msa

###############################
# TEMPLATES
###############################
def mk_mock_template(
  query_sequence: Union[List[str], str], num_temp: int = 1
) -> Dict[str, Any]:
  ln = (
    len(query_sequence)
    if isinstance(query_sequence, str)
    else sum(len(s) for s in query_sequence)
  )
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
    "template_sum_probs": np.zeros([num_temp], dtype=np.float32),
  }
  return template_features


def mk_template(
  a3m_lines: str, template_path: str, query_sequence: str
) -> Dict[str, Any]:
  template_featurizer = templates.HhsearchHitFeaturizer(
    mmcif_dir=template_path,
    max_template_date="2100-01-01",
    max_hits=20,
    kalign_binary_path="kalign",
    release_dates_path=None,
    obsolete_pdbs_path=None,
  )

  hhsearch_pdb70_runner = hhsearch.HHSearch(
    binary_path="hhsearch", databases=[f"{template_path}/pdb70"]
  )

  hhsearch_result = hhsearch_pdb70_runner.query(a3m_lines)
  hhsearch_hits = pipeline.parsers.parse_hhr(hhsearch_result)
  templates_result = template_featurizer.get_templates(
    query_sequence=query_sequence, hits=hhsearch_hits
  )
  return dict(templates_result.features)

def validate_and_fix_mmcif(cif_file: Path):
  """validate presence of _entity_poly_seq in cif file and add revision_date if missing"""
  # check that required poly_seq and revision_date fields are present
  cif_dict = MMCIF2Dict.MMCIF2Dict(cif_file)
  required = [
    "_chem_comp.id",
    "_chem_comp.type",
    "_struct_asym.id",
    "_struct_asym.entity_id",
    "_entity_poly_seq.mon_id",
  ]
  for r in required:
    if r not in cif_dict:
      raise ValueError(f"mmCIF file {cif_file} is missing required field {r}.")
  if "_pdbx_audit_revision_history.revision_date" not in cif_dict:
    logger.info(
      f"Adding missing field revision_date to {cif_file}. Backing up original file to {cif_file}.bak."
    )
    shutil.copy2(cif_file, str(cif_file) + ".bak")
    with open(cif_file, "a") as f:
      f.write(CIF_REVISION_DATE)

def convert_pdb_to_mmcif(pdb_file: Path):
  """convert existing pdb files into mmcif with the required poly_seq and revision_date"""
  i = pdb_file.stem
  cif_file = pdb_file.parent.joinpath(f"{i}.cif")
  if cif_file.is_file():
    return
  parser = PDBParser(QUIET=True)
  structure = parser.get_structure(i, pdb_file)
  cif_io = CFMMCIFIO()
  cif_io.set_structure(structure)
  cif_io.save(str(cif_file))

def mk_hhsearch_db(template_dir: str):
  template_path = Path(template_dir)

  cif_files = template_path.glob("*.cif")
  for cif_file in cif_files:
    validate_and_fix_mmcif(cif_file)

  pdb_files = template_path.glob("*.pdb")
  for pdb_file in pdb_files:
    convert_pdb_to_mmcif(pdb_file)

  pdb70_db_files = template_path.glob("pdb70*")
  for f in pdb70_db_files:
    os.remove(f)

  with open(template_path.joinpath("pdb70_a3m.ffdata"), "w") as a3m, open(
    template_path.joinpath("pdb70_cs219.ffindex"), "w"
  ) as cs219_index, open(
    template_path.joinpath("pdb70_a3m.ffindex"), "w"
  ) as a3m_index, open(
    template_path.joinpath("pdb70_cs219.ffdata"), "w"
  ) as cs219:
    n = 1000000
    index_offset = 0
    cif_files = template_path.glob("*.cif")
    for cif_file in cif_files:
      with open(cif_file) as f:
        cif_string = f.read()
      cif_fh = StringIO(cif_string)
      parser = MMCIFParser(QUIET=True)
      structure = parser.get_structure("none", cif_fh)
      models = list(structure.get_models())
      if len(models) != 1:
        raise ValueError(
          f"Only single model PDBs are supported. Found {len(models)} models."
        )
      model = models[0]
      for chain in model:
        amino_acid_res = []
        for res in chain:
          if res.id[2] != " ":
            raise ValueError(
              f"PDB contains an insertion code at chain {chain.id} and residue "
              f"index {res.id[1]}. These are not supported."
            )
          amino_acid_res.append(
            residue_constants.restype_3to1.get(res.resname, "X")
          )

        protein_str = "".join(amino_acid_res)
        a3m_str = f">{cif_file.stem}_{chain.id}\n{protein_str}\n\0"
        a3m_str_len = len(a3m_str)
        a3m_index.write(f"{n}\t{index_offset}\t{a3m_str_len}\n")
        cs219_index.write(f"{n}\t{index_offset}\t{len(protein_str)}\n")
        index_offset += a3m_str_len
        a3m.write(a3m_str)
        cs219.write("\n\0")
        n += 1