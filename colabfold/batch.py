from __future__ import annotations

import os
ENV = {"TF_FORCE_UNIFIED_MEMORY":"1", "XLA_PYTHON_CLIENT_MEM_FRACTION":"4.0"}
for k,v in ENV.items():
    if k not in os.environ: os.environ[k] = v

import warnings
from Bio import BiopythonDeprecationWarning # what can possibly go wrong...
warnings.simplefilter(action='ignore', category=BiopythonDeprecationWarning)

import json
import logging
import math
import random
import sys
import time
import zipfile
import shutil
import pickle
import gzip

from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
from pathlib import Path
from typing import Any, Callable, Dict, List, Optional, Tuple, Union, TYPE_CHECKING
from io import StringIO

import importlib_metadata
import numpy as np
import pandas

try:
    import alphafold
except ModuleNotFoundError:
    raise RuntimeError(
        "\n\nalphafold is not installed. Please run `pip install colabfold[alphafold]`\n"
    )

from alphafold.common import protein, residue_constants

# delay imports of tensorflow, jax and numpy
# loading these for type checking only can take around 10 seconds just to show a CLI usage message
if TYPE_CHECKING:
    import haiku
    from alphafold.model import model
    from numpy import ndarray

from alphafold.common.protein import Protein
from alphafold.data import (
    feature_processing,
    msa_pairing,
    pipeline,
    pipeline_multimer,
    templates,
)
from alphafold.data.tools import hhsearch
from colabfold.citations import write_bibtex
from colabfold.download import default_data_dir, download_alphafold_params
from colabfold.utils import (
    ACCEPT_DEFAULT_TERMS,
    DEFAULT_API_SERVER,
    NO_GPU_FOUND,
    CIF_REVISION_DATE,
    get_commit,
    safe_filename,
    setup_logging,
    CFMMCIFIO,
)
from colabfold.relax import relax_me

from Bio.PDB import MMCIFParser, PDBParser, MMCIF2Dict
from Bio.PDB.PDBIO import Select

# logging settings
logger = logging.getLogger(__name__)
import jax
import jax.numpy as jnp

# from jax 0.4.6, jax._src.lib.xla_bridge moved to jax._src.xla_bridge
# suppress warnings: Unable to initialize backend 'rocm' or 'tpu'
logging.getLogger('jax._src.xla_bridge').addFilter(lambda _: False) # jax >=0.4.6
logging.getLogger('jax._src.lib.xla_bridge').addFilter(lambda _: False) # jax < 0.4.5

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

modified_mapping = {
  "MSE" : "MET", "MLY" : "LYS", "FME" : "MET", "HYP" : "PRO",
  "TPO" : "THR", "CSO" : "CYS", "SEP" : "SER", "M3L" : "LYS",
  "HSK" : "HIS", "SAC" : "SER", "PCA" : "GLU", "DAL" : "ALA",
  "CME" : "CYS", "CSD" : "CYS", "OCS" : "CYS", "DPR" : "PRO",
  "B3K" : "LYS", "ALY" : "LYS", "YCM" : "CYS", "MLZ" : "LYS",
  "4BF" : "TYR", "KCX" : "LYS", "B3E" : "GLU", "B3D" : "ASP",
  "HZP" : "PRO", "CSX" : "CYS", "BAL" : "ALA", "HIC" : "HIS",
  "DBZ" : "ALA", "DCY" : "CYS", "DVA" : "VAL", "NLE" : "LEU",
  "SMC" : "CYS", "AGM" : "ARG", "B3A" : "ALA", "DAS" : "ASP",
  "DLY" : "LYS", "DSN" : "SER", "DTH" : "THR", "GL3" : "GLY",
  "HY3" : "PRO", "LLP" : "LYS", "MGN" : "GLN", "MHS" : "HIS",
  "TRQ" : "TRP", "B3Y" : "TYR", "PHI" : "PHE", "PTR" : "TYR",
  "TYS" : "TYR", "IAS" : "ASP", "GPL" : "LYS", "KYN" : "TRP",
  "CSD" : "CYS", "SEC" : "CYS"
}

order_to_restype = {v: k for k, v in residue_constants.restype_order_with_x.items()}

class ReplaceOrRemoveHetatmSelect(Select):
  def accept_residue(self, residue):
    hetfield, _, _ = residue.get_id()
    if hetfield != " ":
      if residue.resname in modified_mapping:
        # set unmodified resname
        residue.resname = modified_mapping[residue.resname]
        # clear hetatm flag
        residue._id = (" ", residue._id[1], " ")
        t = residue.full_id
        residue.full_id = (t[0], t[1], t[2], residue._id)
        return 1
      return 0
    else:
      return 1

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
    cif_io.save(str(cif_file), ReplaceOrRemoveHetatmSelect())

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
                logger.warning(f"WARNING: Found {len(models)} models in {cif_file}. The first model will be used as a template.", )
                # raise ValueError(
                #     f"Only single model PDBs are supported. Found {len(models)} models in {cif_file}."
                # )
            model = models[0]
            for chain in model:
                amino_acid_res = []
                for res in chain:
                    if res.id[2] != " ":
                        logger.warning(f"WARNING: Found insertion code at chain {chain.id} and residue index {res.id[1]} of {cif_file}. "
                                       "This file cannot be used as a template.")
                        continue
                        # raise ValueError(
                        #     f"PDB {cif_file} contains an insertion code at chain {chain.id} and residue "
                        #     f"index {res.id[1]}. These are not supported."
                        # )
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

def pad_input(
    input_features: model.features.FeatureDict,
    model_runner: model.RunModel,
    model_name: str,
    pad_len: int,
    use_templates: bool,
) -> model.features.FeatureDict:
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

MODRES = {'MSE':'MET','MLY':'LYS','FME':'MET','HYP':'PRO',
          'TPO':'THR','CSO':'CYS','SEP':'SER','M3L':'LYS',
          'HSK':'HIS','SAC':'SER','PCA':'GLU','DAL':'ALA',
          'CME':'CYS','CSD':'CYS','OCS':'CYS','DPR':'PRO',
          'B3K':'LYS','ALY':'LYS','YCM':'CYS','MLZ':'LYS',
          '4BF':'TYR','KCX':'LYS','B3E':'GLU','B3D':'ASP',
          'HZP':'PRO','CSX':'CYS','BAL':'ALA','HIC':'HIS',
          'DBZ':'ALA','DCY':'CYS','DVA':'VAL','NLE':'LEU',
          'SMC':'CYS','AGM':'ARG','B3A':'ALA','DAS':'ASP',
          'DLY':'LYS','DSN':'SER','DTH':'THR','GL3':'GLY',
          'HY3':'PRO','LLP':'LYS','MGN':'GLN','MHS':'HIS',
          'TRQ':'TRP','B3Y':'TYR','PHI':'PHE','PTR':'TYR',
          'TYS':'TYR','IAS':'ASP','GPL':'LYS','KYN':'TRP',
          'CSD':'CYS','SEC':'CYS'}

def pdb_to_string(
        pdb_file: str, 
        chains: Optional[str] = None, 
        models: Optional[list] = None
    ) -> str:
  '''read pdb file and return as string'''

  if chains is not None:
    if "," in chains: chains = chains.split(",")
    if not isinstance(chains,list): chains = [chains]
  if models is not None:
    if not isinstance(models,list): models = [models]

  modres = {**MODRES}
  lines = []
  seen = []
  model = 1
  
  if "\n" in pdb_file:
    old_lines = pdb_file.split("\n")
  else:
    with open(pdb_file,"rb") as f:
      old_lines = [line.decode("utf-8","ignore").rstrip() for line in f]  
  for line in old_lines:
    if line[:5] == "MODEL":
      model = int(line[5:])
    if models is None or model in models:
      if line[:6] == "MODRES":
        k = line[12:15]
        v = line[24:27]
        if k not in modres and v in residue_constants.restype_3to1:
          modres[k] = v
      if line[:6] == "HETATM":
        k = line[17:20]
        if k in modres:
          line = "ATOM  "+line[6:17]+modres[k]+line[20:]
      if line[:4] == "ATOM":
        chain = line[21:22]
        if chains is None or chain in chains:
          atom = line[12:12+4].strip()
          resi = line[17:17+3]
          resn = line[22:22+5].strip()
          if resn[-1].isalpha(): # alternative atom
            resn = resn[:-1]
            line = line[:26]+" "+line[27:]
          key = f"{model}_{chain}_{resn}_{resi}_{atom}"
          if key not in seen: # skip alternative placements
            lines.append(line)
            seen.append(key)
      if line[:5] == "MODEL" or line[:3] == "TER" or line[:6] == "ENDMDL":
        lines.append(line)
  return "\n".join(lines)

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
    model_runner_and_params: List[Tuple[str, model.RunModel, haiku.Params]],
    initial_guess: str = None,
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
                    # TODO: add pad_input_mulitmer()
                    input_features = feature_dict
                    input_features["asym_id"] = input_features["asym_id"] - input_features["asym_id"][...,0]
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
                    pdb_string = pdb_to_string(initial_guess)
                    input_features["all_atom_positions"] = protein.from_pdb_string(pdb_string).atom_positions
                elif input_guess.suffix == ".cif":
                    input_features["all_atom_positions"] = protein.from_mmcif_string(input_guess.read_text()).atom_positions
                else:
                    raise ValueError(f"Unsupported initial guess file format: {initial_guess}")
                

            ########################
            # predict
            ########################
            start = time.time()

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
                    unrelaxed_protein = protein.from_prediction(
                        features=input_features,
                        result=result, b_factors=b_factors,
                        remove_leading_feature_dimension=("multimer" not in model_type))
                    files.get("unrelaxed",f"r{recycles}.pdb").write_text(protein.to_pdb(unrelaxed_protein))

                    if save_all:
                        with files.get("all",f"r{recycles}.pickle").open("wb") as handle:
                            pickle.dump(result, handle)
                    del unrelaxed_protein

            return_representations = save_all or save_single_representations or save_pair_representations

            # predict
            result, recycles = \
            model_runner.predict(input_features,
                random_seed=seed,
                return_representations=return_representations,
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
                  pae   = result["predicted_aligned_error"][:seq_len,:seq_len]
                  scores.update({"max_pae": pae.max().astype(float).item(),
                                 "pae": np.around(pae.astype(float), 2).tolist()})
                  for k in ["ptm","iptm"]:
                    if k in conf[-1]: scores[k] = np.around(conf[-1][k], 2).item()
                  del pae
                del plddt
                json.dump(scores, handle)

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

def decode_structure_sequences(
    aatype_array: List[int], 
    chain_index_array: List[int], 
    order_dict: Dict[int, str] = order_to_restype
) -> List[str]:
    decoded_sequences = []
    current_sequence = []

    for i in range(len(aatype_array)):
        amino_acid = order_dict[aatype_array[i]]
        if i == 0 or chain_index_array[i] == chain_index_array[i - 1]:
            current_sequence.append(amino_acid)
        else:
            decoded_sequences.append("".join(current_sequence))
            current_sequence = [amino_acid]

    # Append the last sequence
    decoded_sequences.append("".join(current_sequence))

    return decoded_sequences

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
        elif input_path.suffix in [".pdb", ".cif"]:
            if input_path.suffix == ".pdb":
                pdb_string = pdb_to_string(input_path.read_text())
                prot = protein.from_pdb_string(pdb_string)
            elif input_path.suffix == ".cif":
                prot = protein.from_mmcif_string(input_path.read_text())
            header = input_path.stem
            sequences = decode_structure_sequences(prot.aatype, prot.chain_index)

            if len(sequences) == 0:
                raise ValueError(f"{input_path} is empty")
                
            queries = [(header, sequences, None)]

        else:
            raise ValueError(f"Unknown file format {input_path.suffix}")
    else:
        assert input_path.is_dir(), "Expected either an input file or a input directory"
        queries = []
        for file in sorted(input_path.iterdir()):
            if not file.is_file():
                continue
            if file.suffix.lower() not in [".a3m", ".fasta", ".faa", ".pdb", ".cif"]:
                logger.warning(f"non-fasta/a3m/pdb/cif file in input directory: {file}")
                continue
            if file.suffix.lower() in [".pdb", ".cif"]:
                header = file.stem
                if file.suffix.lower() == ".pdb":
                    pdb_string = pdb_to_string(file.read_text())
                    prot = protein.from_pdb_string(pdb_string)
                else:  # file.suffix.lower() == ".cif"
                    prot = protein.from_mmcif_string(file.read_text())
                sequences = decode_structure_sequences(prot.aatype, prot.chain_index)

                if len(sequences) == 0:
                    logger.error(f"{file} is empty")
                    continue

                queries.append((header, sequences, None))
            else:  # file.suffix.lower() in [".a3m", ".fasta", ".faa"]
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
    for job_number, (_, query_sequence, a3m_lines) in enumerate(queries):
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
    a3m_lines: Optional[List[str]],
    result_dir: Path,
    msa_mode: str,
    use_templates: bool,
    custom_template_path: str,
    pair_mode: str,
    pairing_strategy: str = "greedy",
    host_url: str = DEFAULT_API_SERVER,
    user_agent: str = "",
) -> Tuple[
    Optional[List[str]], Optional[List[str]], List[str], List[int], List[Dict[str, Any]]
]:
    from colabfold.colabfold import run_mmseqs2

    use_env = msa_mode == "mmseqs2_uniref_env" or msa_mode == "mmseqs2_uniref_env_envpair"
    use_envpair = msa_mode == "mmseqs2_uniref_env_envpair"
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
        # Skip template search when custom_template_path is provided
        if custom_template_path is not None:
            if msa_mode == "single_sequence":
                a3m_lines = []
                num = 101
                for i, seq in enumerate(query_seqs_unique):
                    a3m_lines.append(f">{num + i}\n{seq}")

            if a3m_lines is None:
                a3m_lines_mmseqs2 = run_mmseqs2(
                    query_seqs_unique,
                    str(result_dir.joinpath(jobname)),
                    use_env,
                    use_templates=False,
                    host_url=host_url,
                    user_agent=user_agent,
                )
            else:
                a3m_lines_mmseqs2 = a3m_lines
            template_paths = {}
            for index in range(0, len(query_seqs_unique)):
                template_paths[index] = custom_template_path
        else:
            a3m_lines_mmseqs2, template_paths = run_mmseqs2(
                query_seqs_unique,
                str(result_dir.joinpath(jobname)),
                use_env,
                use_templates=True,
                host_url=host_url,
                user_agent=user_agent,
            )
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
                use_env,
                use_pairing=False,
                host_url=host_url,
                user_agent=user_agent,
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
                use_envpair,
                use_pairing=True,
                pairing_strategy=pairing_strategy,
                host_url=host_url,
                user_agent=user_agent,
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

def build_multimer_feature(paired_msa: str) -> Dict[str, ndarray]:
    parsed_paired_msa = pipeline.parsers.parse_a3m(paired_msa)
    return {
        f"{k}_all_seq": v
        for k, v in pipeline.make_msa_features([parsed_paired_msa]).items()
    }

def process_multimer_features(
    features_for_chain: Dict[str, Dict[str, ndarray]],
    min_num_seq: int = 512,
) -> Dict[str, ndarray]:
    all_chain_features = {}
    for chain_id, chain_features in features_for_chain.items():
        all_chain_features[chain_id] = pipeline_multimer.convert_monomer_features(
            chain_features, chain_id
        )

    all_chain_features = pipeline_multimer.add_assembly_features(all_chain_features)
    # np_example = feature_processing.pair_and_merge(
    #    all_chain_features=all_chain_features, is_prokaryote=is_prokaryote)
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
    np_example = pipeline_multimer.pad_msa(np_example, min_num_seq=min_num_seq)
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
    max_seq: int,
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
                f"{model_type} complex does not consider templates. Chose multimer model-type for template support."
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

        if "multimer" in model_type:
            # combine features across all chains
            input_feature = process_multimer_features(features_for_chain, min_num_seq=max_seq + 4)
            domain_names = {
                chain: [
                    name.decode("UTF-8")
                    for name in feature["template_domain_names"]
                    if name != b"none"
                ]
                for (chain, feature) in features_for_chain.items()
            }
        else:
            input_feature = features_for_chain[protein.PDB_CHAIN_IDS[0]]
            input_feature["asym_id"] = np.zeros(input_feature["aatype"].shape[0],dtype=int)
            domain_names = {
                protein.PDB_CHAIN_IDS[0]: [
                    name.decode("UTF-8")
                    for name in input_feature["template_domain_names"]
                    if name != b"none"
                ]
            }
    return (input_feature, domain_names)

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

        # if sequence is paired add them to output
        if (
            not is_single_protein
            and not is_homooligomer
            and sum(has_amino_acid) > 1 # at least 2 sequences are paired
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

    return (
        unpaired_msa,
        paired_msa,
        query_seqs_unique,
        query_seqs_cardinality,
        template_features,
    )

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


def put_mmciffiles_into_resultdir(
    pdb_hit_file: Path,
    local_pdb_path: Path,
    result_dir: Path,
    max_num_templates: int = 20,
):
    """Put mmcif files from local_pdb_path into result_dir and unzip them.
    max_num_templates is the maximum number of templates to use (default: 20).
    Args:
        pdb_hit_file (Path): Path to pdb_hit_file
        local_pdb_path (Path): Path to local_pdb_path
        result_dir (Path): Path to result_dir
        max_num_templates (int): Maximum number of templates to use
    """
    pdb_hit_file = Path(pdb_hit_file)
    local_pdb_path = Path(local_pdb_path)
    result_dir = Path(result_dir)
    result_dir.mkdir(parents=True, exist_ok=True)

    query_ids = []
    with open(pdb_hit_file, "r") as f:
        for line in f:
            query_id = line.split("\t")[0]
            query_ids.append(query_id)
            if query_ids.count(query_id) > max_num_templates:
                continue
            else:
                pdb_id = line.split("\t")[1][0:4]
                divided_pdb_id = pdb_id[1:3]
                gzipped_divided_mmcif_file = local_pdb_path / divided_pdb_id / (pdb_id + ".cif.gz")
                gzipped_mmcif_file = local_pdb_path / (pdb_id + ".cif.gz")
                unzipped_mmcif_file = local_pdb_path / (pdb_id + ".cif")
                result_file = result_dir / (pdb_id + ".cif")
                possible_files = [gzipped_divided_mmcif_file, gzipped_mmcif_file, unzipped_mmcif_file]
                for file in possible_files:
                    if file == gzipped_divided_mmcif_file or file == gzipped_mmcif_file:
                        if file.exists():
                            with gzip.open(file, "rb") as f_in:
                                with open(result_file, "wb") as f_out:
                                    shutil.copyfileobj(f_in, f_out)
                                    break
                    else:
                        # unzipped_mmcif_file
                        if file.exists():
                            shutil.copyfile(file, result_file)
                            break
                if not result_file.exists():
                    print(f"WARNING: {pdb_id} does not exist in {local_pdb_path}.")


def run(
    queries: List[Tuple[str, Union[str, List[str]], Optional[List[str]]]],
    result_dir: Union[str, Path],
    num_models: int,
    is_complex: bool,
    num_recycles: Optional[int] = None,
    recycle_early_stop_tolerance: Optional[float] = None,
    model_order: List[int] = [1,2,3,4,5],
    initial_guess: str = None,
    num_ensemble: int = 1,
    model_type: str = "auto",
    msa_mode: str = "mmseqs2_uniref_env",
    use_templates: bool = False,
    custom_template_path: str = None,
    num_relax: int = 0,
    relax_max_iterations: int = 0,
    relax_tolerance: float = 2.39,
    relax_stiffness: float = 10.0,
    relax_max_outer_iterations: int = 3,
    keep_existing_results: bool = True,
    rank_by: str = "auto",
    pair_mode: str = "unpaired_paired",
    pairing_strategy: str = "greedy",
    data_dir: Union[str, Path] = default_data_dir,
    host_url: str = DEFAULT_API_SERVER,
    user_agent: str = "",
    random_seed: int = 0,
    num_seeds: int = 1,
    recompile_padding: Union[int, float] = 10,
    zip_results: bool = False,
    prediction_callback: Callable[[Any, Any, Any, Any, Any], Any] = None,
    save_single_representations: bool = False,
    save_pair_representations: bool = False,
    jobname_prefix: Optional[str] = None,
    save_all: bool = False,
    save_recycles: bool = False,
    use_dropout: bool = False,
    use_gpu_relax: bool = False,
    stop_at_score: float = 100,
    dpi: int = 200,
    max_seq: Optional[int] = None,
    max_extra_seq: Optional[int] = None,
    pdb_hit_file: Optional[Path] = None,
    local_pdb_path: Optional[Path] = None,
    use_cluster_profile: bool = True,
    feature_dict_callback: Callable[[Any], Any] = None,
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
    from colabfold.colabfold import plot_paes, plot_plddts
    from colabfold.plot import plot_msa_v2

    data_dir = Path(data_dir)
    result_dir = Path(result_dir)
    result_dir.mkdir(exist_ok=True)
    model_type = set_model_type(is_complex, model_type)

    # backward-compatibility with old options
    old_names = {"MMseqs2 (UniRef+Environmental)":"mmseqs2_uniref_env",
                 "MMseqs2 (UniRef+Environmental+Env. Pairing)":"mmseqs2_uniref_env_envpair",
                 "MMseqs2 (UniRef only)":"mmseqs2_uniref",
                 "unpaired+paired":"unpaired_paired"}
    msa_mode   = old_names.get(msa_mode,msa_mode)
    pair_mode  = old_names.get(pair_mode,pair_mode)
    feature_dict_callback = kwargs.pop("input_features_callback", feature_dict_callback)
    use_dropout           = kwargs.pop("training", use_dropout)
    use_fuse              = kwargs.pop("use_fuse", True)
    use_bfloat16          = kwargs.pop("use_bfloat16", True)
    max_msa               = kwargs.pop("max_msa",None)
    if max_msa is not None:
        max_seq, max_extra_seq = [int(x) for x in max_msa.split(":")]

    if kwargs.pop("use_amber", False) and num_relax == 0:
        num_relax = num_models * num_seeds

    if len(kwargs) > 0:
        print(f"WARNING: the following options are not being used: {kwargs}")

    # decide how to rank outputs
    if rank_by == "auto":
        rank_by = "multimer" if is_complex else "plddt"
    if "ptm" not in model_type and "multimer" not in model_type:
        rank_by = "plddt"

    # get max length
    max_len = 0
    max_num = 0
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

    # sort model order
    model_order.sort()

    # initial guess
    if initial_guess is not None:
        logger.info(f'Using initial guess: {initial_guess}')

    # Record the parameters of this run
    config = {
        "num_queries": len(queries),
        "use_templates": use_templates,
        "num_relax": num_relax,
        "relax_max_iterations": relax_max_iterations,
        "relax_tolerance": relax_tolerance,
        "relax_stiffness": relax_stiffness,
        "relax_max_outer_iterations": relax_max_outer_iterations,
        "msa_mode": msa_mode,
        "model_type": model_type,
        "num_models": num_models,
        "num_recycles": num_recycles,
        "recycle_early_stop_tolerance": recycle_early_stop_tolerance,
        "num_ensemble": num_ensemble,
        "model_order": model_order,
        "initial_guess": initial_guess,
        "keep_existing_results": keep_existing_results,
        "rank_by": rank_by,
        "max_seq": max_seq,
        "max_extra_seq": max_extra_seq,
        "pair_mode": pair_mode,
        "pairing_strategy": pairing_strategy,
        "host_url": host_url,
        "user_agent": user_agent,
        "stop_at_score": stop_at_score,
        "random_seed": random_seed,
        "num_seeds": num_seeds,
        "recompile_padding": recompile_padding,
        "commit": get_commit(),
        "use_dropout": use_dropout,
        "use_cluster_profile": use_cluster_profile,
        "use_fuse": use_fuse,
        "use_bfloat16": use_bfloat16,
        "version": importlib_metadata.version("colabfold"),
    }
    config_out_file = result_dir.joinpath("config.json")
    config_out_file.write_text(json.dumps(config, indent=4))
    use_env = "env" in msa_mode
    use_msa = "mmseqs2" in msa_mode
    use_amber = num_models > 0 and num_relax > 0

    bibtex_file = write_bibtex(
        model_type if num_models > 0 else "", use_msa, use_env, use_templates, use_amber, result_dir
    )

    if pdb_hit_file is not None:
        if local_pdb_path is None:
            raise ValueError("local_pdb_path is not specified.")
        else:
            custom_template_path = result_dir / "templates"
            put_mmciffiles_into_resultdir(pdb_hit_file, local_pdb_path, custom_template_path)

    if custom_template_path is not None:
        mk_hhsearch_db(custom_template_path)

    pad_len = 0
    ranks, metrics = [],[]
    first_job = True
    job_number = 0
    for job_number, (raw_jobname, query_sequence, a3m_lines) in enumerate(queries):
        if jobname_prefix is not None:
            # pad job number based on number of queries
            fill = len(str(len(queries)))
            jobname = safe_filename(jobname_prefix) + "_" + str(job_number).zfill(fill)
            job_number += 1
        else:
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
            pickled_msa_and_templates = result_dir.joinpath(f"{jobname}.pickle")
            if pickled_msa_and_templates.is_file():
                with open(pickled_msa_and_templates, 'rb') as f:
                    (unpaired_msa, paired_msa, query_seqs_unique, query_seqs_cardinality, template_features) = pickle.load(f)
                logger.info(f"Loaded {pickled_msa_and_templates}")

            else:
                if a3m_lines is None:
                    (unpaired_msa, paired_msa, query_seqs_unique, query_seqs_cardinality, template_features) \
                    = get_msa_and_templates(jobname, query_sequence, a3m_lines, result_dir, msa_mode, use_templates,
                        custom_template_path, pair_mode, pairing_strategy, host_url, user_agent)

                elif a3m_lines is not None:
                    (unpaired_msa, paired_msa, query_seqs_unique, query_seqs_cardinality, template_features) \
                    = unserialize_msa(a3m_lines, query_sequence)
                    if use_templates:
                        (_, _, _, _, template_features) \
                            = get_msa_and_templates(jobname, query_seqs_unique, unpaired_msa, result_dir, 'single_sequence', use_templates,
                                custom_template_path, pair_mode, pairing_strategy, host_url, user_agent)

                if num_models == 0:
                    with open(pickled_msa_and_templates, 'wb') as f:
                        pickle.dump((unpaired_msa, paired_msa, query_seqs_unique, query_seqs_cardinality, template_features), f)
                    logger.info(f"Saved {pickled_msa_and_templates}")

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
                                     template_features, is_complex, model_type, max_seq=max_seq)

            # to allow display of MSA info during colab/chimera run (thanks tomgoddard)
            if feature_dict_callback is not None:
                feature_dict_callback(feature_dict)

        except Exception as e:
            logger.exception(f"Could not generate input features {jobname}: {e}")
            continue

        ###############
        # save plots not requiring prediction
        ###############

        result_files = []

        # make msa plot
        msa_plot = plot_msa_v2(feature_dict, dpi=dpi)
        coverage_png = result_dir.joinpath(f"{jobname}_coverage.png")
        msa_plot.savefig(str(coverage_png), bbox_inches='tight')
        msa_plot.close()
        result_files.append(coverage_png)

        if use_templates:
            templates_file = result_dir.joinpath(f"{jobname}_template_domain_names.json")
            templates_file.write_text(json.dumps(domain_names))
            result_files.append(templates_file)

        result_files.append(result_dir.joinpath(jobname + ".a3m"))
        result_files += [bibtex_file, config_out_file]

        ######################
        # predict structures
        ######################
        if num_models > 0:
            try:
                # get list of lengths
                query_sequence_len_array = sum([[len(x)] * y
                    for x,y in zip(query_seqs_unique, query_seqs_cardinality)],[])

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
                        # get number of sequences
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
                        model_order=model_order,
                        model_type=model_type,
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
                    initial_guess=initial_guess,
                    model_type=model_type,
                    model_runner_and_params=model_runner_and_params,
                    num_relax=num_relax,
                    relax_max_iterations=relax_max_iterations,
                    relax_tolerance=relax_tolerance,
                    relax_stiffness=relax_stiffness,
                    relax_max_outer_iterations=relax_max_outer_iterations,
                    rank_by=rank_by,
                    stop_at_score=stop_at_score,
                    prediction_callback=prediction_callback,
                    use_gpu_relax=use_gpu_relax,
                    random_seed=random_seed,
                    num_seeds=num_seeds,
                    save_all=save_all,
                    save_single_representations=save_single_representations,
                    save_pair_representations=save_pair_representations,
                    save_recycles=save_recycles,
                )
                result_files += results["result_files"]
                ranks.append(results["rank"])
                metrics.append(results["metric"])

            except RuntimeError as e:
                # This normally happens on OOM. TODO: Filter for the specific OOM error message
                logger.error(f"Could not predict {jobname}. Not Enough GPU memory? {e}")
                continue

            ###############
            # save prediction plots
            ###############

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

        if zip_results:
            with zipfile.ZipFile(result_zip, "w") as result_zip:
                for file in result_files:
                    result_zip.write(file, arcname=file.name)

            # Delete only after the zip was successful, and also not the bibtex and config because we need those again
            for file in result_files:
                if file != bibtex_file and file != config_out_file:
                    file.unlink()
        else:
            if num_models > 0:
                is_done_marker.touch()

    logger.info("Done")
    return {"rank":ranks,"metric":metrics}

def set_model_type(is_complex: bool, model_type: str) -> str:
    # backward-compatibility with old options
    old_names = {
        "AlphaFold2-multimer-v1":"alphafold2_multimer_v1",
        "AlphaFold2-multimer-v2":"alphafold2_multimer_v2",
        "AlphaFold2-multimer-v3":"alphafold2_multimer_v3",
        "AlphaFold2-ptm":        "alphafold2_ptm",
        "AlphaFold2":            "alphafold2",
        "DeepFold":              "deepfold_v1",
    }
    model_type = old_names.get(model_type, model_type)
    if model_type == "auto":
        if is_complex:
            model_type = "alphafold2_multimer_v3"
        else:
            model_type = "alphafold2_ptm"
    return model_type

def main():
    parser = ArgumentParser(formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument(
        "input",
        default="input",
        help="One of: 1) directory with FASTA/A3M files, 2) CSV/TSV file, 3) FASTA file or 4) A3M file.",
    )
    parser.add_argument("results", help="Results output directory.")

    msa_group = parser.add_argument_group("MSA arguments", "")
    msa_group.add_argument(
        "--msa-only",
        action="store_true",
        help="Query and store MSAs from the MSA server without structure prediction",
    )
    msa_group.add_argument(
        "--msa-mode",
        default="mmseqs2_uniref_env",
        choices=[
            "mmseqs2_uniref_env",
            "mmseqs2_uniref_env_envpair",
            "mmseqs2_uniref",
            "single_sequence",
        ],
        help="Databases to use to create the MSA: UniRef30+Environmental (default), UniRef30 only or None. "
        "Using an A3M file as input overwrites this option.",
    )
    msa_group.add_argument(
        "--pair-mode",
        help="Multimer MSA pairing mode for complex prediction: unpaired MSA only, paired MSA only, both (default).",
        type=str,
        default="unpaired_paired",
        choices=["unpaired", "paired", "unpaired_paired"],
    )
    msa_group.add_argument(
        "--pair-strategy",
        help="How sequences are paired during MSA pairing for complex prediction. "
        "complete: MSA sequences should only be paired if the same species exists in all MSAs. "
        "greedy: MSA sequences should only be paired if the same species exists in at least two MSAs. "
        "Typically, greedy produces better predictions as it results in more paired sequences. "
        "However, in some cases complete pairing might help, especially if MSAs are already large and can be well paired. ",
        type=str,
        default="greedy",
        choices=["complete", "greedy"],
    )
    msa_group.add_argument(
        "--templates",
        default=False,
        action="store_true",
        help="Query PDB templates from the MSA server. "
        'If this parameter is not set, "--custom-template-path" and "--pdb-hit-file" will not be used. '
        "Warning: This can result in the MSA server being queried with A3M input. "
    )
    msa_group.add_argument(
        "--custom-template-path",
        type=str,
        default=None,
        help="Directory with PDB files to provide as custom templates to the predictor. "
        "No templates will be queried from the MSA server. "
        "'--templates' argument is also required to enable this.",
    )
    msa_group.add_argument(
        "--pdb-hit-file",
        default=None,
        help="Path to a BLAST-m8 formatted PDB hit file corresponding to the input A3M file (e.g. pdb70.m8). "
        "Typically, this parameter should be used for a MSA generated by 'colabfold_search'. "
        "'--templates' argument is also required to enable this.",
    )
    msa_group.add_argument(
        "--local-pdb-path",
        default=None,
        help="Directory of a local mirror of the PDB mmCIF database (e.g. /path/to/pdb/divided). "
        "If provided, PDB files from the directory are used for templates specified by '--pdb-hit-file'. ",
    )

    pred_group = parser.add_argument_group("Prediction arguments", "")
    pred_group.add_argument(
        "--num-recycle",
        help="Number of prediction recycles. "
        "Increasing recycles can improve the prediction quality but slows down the prediction.",
        type=int,
        default=None,
    )
    pred_group.add_argument(
        "--recycle-early-stop-tolerance",
        help="Specify convergence criteria. "
        "Run recycles until the distance between recycles is within the given tolerance value.",
        type=float,
        default=None,
    )
    pred_group.add_argument(
        "--num-ensemble",
        help="Number of ensembles. "
        "The trunk of the network is run multiple times with different random choices for the MSA cluster centers. "
        "This can result in a better prediction at the cost of longer runtime. ",
        type=int,
        default=1,
    )
    pred_group.add_argument(
        "--num-seeds",
        help="Number of seeds to try. Will iterate from range(random_seed, random_seed+num_seeds). "
        "This can result in a better/different prediction at the cost of longer runtime. ",
        type=int,
        default=1,
    )
    pred_group.add_argument(
        "--random-seed",
        help="Changing the seed for the random number generator can result in better/different structure predictions.",
        type=int,
        default=0,
    )
    pred_group.add_argument(
        "--num-models",
        help="Number of models to use for structure prediction. "
        "Reducing the number of models speeds up the prediction but results in lower quality.",
        type=int,
        default=5,
        choices=[1, 2, 3, 4, 5],
    )
    pred_group.add_argument(
        "--model-type",
        help="Predict structure/complex using the given model. "
        'Auto will pick "alphafold2_ptm" for structure predictions and "alphafold2_multimer_v3" for complexes. '
        "Older versions of the AF2 models are generally worse, however they can sometimes result in better predictions. "
        "If the model is not already downloaded, it will be automatically downloaded. ",
        type=str,
        default="auto",
        choices=[
            "auto",
            "alphafold2",
            "alphafold2_ptm",
            "alphafold2_multimer_v1",
            "alphafold2_multimer_v2",
            "alphafold2_multimer_v3",
            "deepfold_v1",
        ],
    )
    pred_group.add_argument("--model-order", default="1,2,3,4,5", type=str)
    pred_group.add_argument(
        "--initial-guess",
        nargs="?",
        const=True,
        help="Specify a starting model for the prediction. If the main input file is a PDB format, "
        "it will be used as the initial guess. Otherwise, you can provide an input file with this flag, "
        "which will override the main input."
    )
    pred_group.add_argument(
        "--use-dropout",
        default=False,
        action="store_true",
        help="Activate dropouts during inference to sample from uncertainty of the models. "
        "This can result in different predictions and can be (carefully!) used for conformations sampling.",
    )
    pred_group.add_argument(
        "--max-seq",
        help="Number of sequence clusters to use. "
        "This can result in different predictions and can be (carefully!) used for conformations sampling.",
        type=int,
        default=None,
    )
    pred_group.add_argument(
        "--max-extra-seq",
        help="Number of extra sequences to use. "
        "This can result in different predictions and can be (carefully!) used for conformations sampling.",
        type=int,
        default=None,
    )
    pred_group.add_argument(
        "--max-msa",
        help="Defines: `max-seq:max-extra-seq` number of sequences to use in one go. "
        '"--max-seq" and "--max-extra-seq" are ignored if this parameter is set.',
        type=str,
        default=None,
    )
    pred_group.add_argument(
        "--disable-cluster-profile",
        default=False,
        action="store_true",
        help="Experimental: For multimer models, disable cluster profiles.",
    )
    pred_group.add_argument("--data", help="Path to AlphaFold2 weights directory.")

    relax_group = parser.add_argument_group("Relaxation arguments", "")
    relax_group.add_argument(
        "--amber",
        default=False,
        action="store_true",
        help="Enable OpenMM/Amber for structure relaxation. "
        "Can improve the quality of side-chains at a cost of longer runtime. "
    )
    relax_group.add_argument(
        "--num-relax",
        help="Specify how many of the top ranked structures to relax using OpenMM/Amber. "
        "Typically, relaxing the top-ranked prediction is enough and speeds up the runtime. ",
        type=int,
        default=0,
    )
    relax_group.add_argument(
        "--relax-max-iterations",
        type=int,
        default=2000,
        help="Maximum number of iterations for the relaxation process. "
        "AlphaFold2 sets this to unlimited (0), however, we found that this can lead to very long relaxation times for some inputs.",
    )
    relax_group.add_argument(
        "--relax-tolerance",
        type=float,
        default=2.39,
        help="Tolerance threshold for relaxation convergence.",
    )
    relax_group.add_argument(
        "--relax-stiffness",
        type=float,
        default=10.0,
        help="Stiffness parameter for relaxation.",
    )
    relax_group.add_argument(
        "--relax-max-outer-iterations",
        type=int,
        default=3,
        help="Maximum number of outer iterations for the relaxation process.",
    )
    relax_group.add_argument(
        "--use-gpu-relax",
        default=False,
        action="store_true",
        help="Run OpenMM/Amber on GPU instead of CPU. "
        "This can significantly speed up the relaxation runtime, however, might lead to compatibility issues with CUDA. "
        "Unsupported on AMD/ROCM and Apple Silicon.",
    )

    output_group = parser.add_argument_group("Output arguments", "")
    output_group.add_argument(
        "--rank",
        help='Choose metric to rank the "--num-models" predicted models.',
        type=str,
        default="auto",
        choices=["auto", "plddt", "ptm", "iptm", "multimer"],
    )
    output_group.add_argument(
        "--stop-at-score",
        help="Compute models until pLDDT (single chain) or pTM-score (multimer) > threshold is reached. "
        "This speeds up prediction by running less models for easier queries.",
        type=float,
        default=100,
    )
    output_group.add_argument(
        "--jobname-prefix",
        help="If set, the jobname will be prefixed with the given string and a running number, instead of the input headers/accession.",
        type=str,
        default=None,
    )
    output_group.add_argument(
        "--save-all",
        default=False,
        action="store_true",
        help="Save all raw outputs from model to a pickle file. "
        "Useful for downstream use in other models."
    )
    output_group.add_argument(
        "--save-recycles",
        default=False,
        action="store_true",
        help="Save all intermediate predictions at each recycle iteration.",
    )
    output_group.add_argument(
        "--save-single-representations",
        default=False,
        action="store_true",
        help="Save the single representation embeddings of all models.",
    )
    output_group.add_argument(
        "--save-pair-representations",
        default=False,
        action="store_true",
        help="Save the pair representation embeddings of all models.",
    )
    output_group.add_argument(
        "--overwrite-existing-results",
        default=False,
        action="store_true",
        help="Do not recompute results, if a query has already been predicted.",
    )
    output_group.add_argument(
        "--zip",
        default=False,
        action="store_true",
        help="Zip all results into one <jobname>.result.zip and delete the original files.",
    )
    output_group.add_argument(
        "--sort-queries-by",
        help="Sort input queries by: none, length, random. "
        "Sorting by length speeds up prediction as models are recompiled less often.",
        type=str,
        default="length",
        choices=["none", "length", "random"],
    )

    adv_group = parser.add_argument_group(
        "Advanced arguments", ""
    )
    adv_group.add_argument(
        "--host-url",
        default=DEFAULT_API_SERVER,
        help="Which MSA server should be queried. By default, the free public MSA server hosted by the ColabFold team is queried. "
    )
    adv_group.add_argument(
        "--disable-unified-memory",
        default=False,
        action="store_true",
        help="If you are getting TensorFlow/Jax errors, it might help to disable this.",
    )
    adv_group.add_argument(
        "--recompile-padding",
        type=int,
        default=10,
        help="Whenever the input length changes, the model needs to be recompiled. "
        "We pad sequences by the specified length, so we can e.g., compute sequences from length 100 to 110 without recompiling. "
        "Individual predictions will become marginally slower due to longer input, "
        "but overall performance increases due to not recompiling. "
        "Set to 0 to disable.",
    )

    args = parser.parse_args()

    if (args.custom_template_path is not None) and (args.pdb_hit_file is not None):
        raise RuntimeError("Arguments --pdb-hit-file and --custom-template-path cannot be used simultaneously.")
    # disable unified memory
    if args.disable_unified_memory:
        for k in ENV.keys():
            if k in os.environ: del os.environ[k]

    setup_logging(Path(args.results).joinpath("log.txt"))

    version = importlib_metadata.version("colabfold")
    commit = get_commit()
    if commit:
        version += f" ({commit})"

    logger.info(f"Running colabfold {version}")

    data_dir = Path(args.data or default_data_dir)

    queries, is_complex = get_queries(args.input, args.sort_queries_by)
    model_type = set_model_type(is_complex, args.model_type)

    # use pdb or cif input as initial guess
    if args.initial_guess is not None:
        if isinstance(args.initial_guess, str) and Path(args.initial_guess).suffix in (".pdb", ".cif"):
            initial_guess = args.initial_guess
        elif Path(args.input).suffix in (".pdb", ".cif"):
            initial_guess = args.input
        else:
            raise ValueError("Provide PDB or CIF file for initial guess.")
    else:
        initial_guess = None

    if args.msa_only:
        args.num_models = 0

    if args.num_models > 0:
        download_alphafold_params(model_type, data_dir)

    if args.msa_mode != "single_sequence" and not args.templates:
        uses_api = any((query[2] is None for query in queries))
        if uses_api and args.host_url == DEFAULT_API_SERVER:
            print(ACCEPT_DEFAULT_TERMS, file=sys.stderr)

    model_order = [int(i) for i in args.model_order.split(",")]

    assert args.recompile_padding >= 0, "Can't apply negative padding"

    # backward compatibility
    if args.amber and args.num_relax == 0:
        args.num_relax = args.num_models * args.num_seeds

    user_agent = f"colabfold/{version}"

    run(
        queries=queries,
        result_dir=args.results,
        use_templates=args.templates,
        custom_template_path=args.custom_template_path,
        num_relax=args.num_relax,
        relax_max_iterations=args.relax_max_iterations,
        relax_tolerance=args.relax_tolerance,
        relax_stiffness=args.relax_stiffness,
        relax_max_outer_iterations=args.relax_max_outer_iterations,
        msa_mode=args.msa_mode,
        model_type=model_type,
        num_models=args.num_models,
        num_recycles=args.num_recycle,
        recycle_early_stop_tolerance=args.recycle_early_stop_tolerance,
        num_ensemble=args.num_ensemble,
        model_order=model_order,
        initial_guess=initial_guess,
        is_complex=is_complex,
        keep_existing_results=not args.overwrite_existing_results,
        rank_by=args.rank,
        pair_mode=args.pair_mode,
        pairing_strategy=args.pair_strategy,
        data_dir=data_dir,
        host_url=args.host_url,
        user_agent=user_agent,
        random_seed=args.random_seed,
        num_seeds=args.num_seeds,
        stop_at_score=args.stop_at_score,
        recompile_padding=args.recompile_padding,
        zip_results=args.zip,
        save_single_representations=args.save_single_representations,
        save_pair_representations=args.save_pair_representations,
        use_dropout=args.use_dropout,
        max_seq=args.max_seq,
        max_extra_seq=args.max_extra_seq,
        max_msa=args.max_msa,
        pdb_hit_file=args.pdb_hit_file,
        local_pdb_path=args.local_pdb_path,
        use_cluster_profile=not args.disable_cluster_profile,
        use_gpu_relax = args.use_gpu_relax,
        jobname_prefix=args.jobname_prefix,
        save_all=args.save_all,
        save_recycles=args.save_recycles,
    )

if __name__ == "__main__":
    main()
