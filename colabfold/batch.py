from __future__ import annotations

import os
import importlib.util
# the jax ROCm plugin is how a ROCm install shows itself before jax is imported
_IS_ROCM = any(importlib.util.find_spec(m) is not None
               for m in ("jax_rocm7_plugin", "jax_rocm60_plugin", "jax_rocm_plugin"))
if _IS_ROCM:
    # ROCm 7.2 segfaults on graph segment scheduling and on the plugin's profiler hook
    ENV = {"DEBUG_HIP_GRAPH_SEGMENT_SCHEDULING":"0", "ROCPROFILER_REGISTER_ENABLED":"0"}
else:
    # spilling past the GPU takes managed memory
    ENV = {"TF_FORCE_UNIFIED_MEMORY":"1", "XLA_PYTHON_CLIENT_MEM_FRACTION":"4.0"}
for k,v in ENV.items():
    if k not in os.environ: os.environ[k] = v

import warnings
from Bio import BiopythonDeprecationWarning # what can possibly go wrong...
warnings.simplefilter(action='ignore', category=BiopythonDeprecationWarning)

import functools
import json
hasOrjson = False
try:
    import orjson
    hasOrjson = True
except ImportError:
    pass
import logging
import sys
import time
import zipfile
import shutil
import pickle
import gzip

from argparse import (ArgumentParser, ArgumentDefaultsHelpFormatter,
                      BooleanOptionalAction, SUPPRESS)
from pathlib import Path
from typing import Any, Callable, Dict, List, Optional, Tuple, Union, TYPE_CHECKING
from io import StringIO

import importlib_metadata
import numpy as np

# delay imports of jax and numpy
# loading these for type checking only can take around 10 seconds just to show a CLI usage message
if TYPE_CHECKING:
    import haiku
    from numpy import ndarray

from colabfold.backend import AF3_MODELS, is_af3_model
from colabfold.citations import write_bibtex
from colabfold.download import default_data_dir, download_alphafold_params
from colabfold.utils import (
    ACCEPT_DEFAULT_TERMS,
    DEFAULT_API_SERVER,
    NO_GPU_FOUND,
    CIF_REVISION_DATE,
    get_commit,
    setup_logging,
    CFMMCIFIO,
    AF3Utils,
)
from colabfold.input import (
    dropped_entities,
    pair_msa,
    refuse_mixed_msas,
    warn_about_dropped_entities,
    msa_to_str,
    get_queries,
    safe_filename,
    modified_mapping,
    pdb_to_string,
    restype_3to1,
)
from colabfold.relax import relax_me

from Bio.PDB import MMCIFParser, PDBParser, MMCIF2Dict
from Bio.PDB.PDBIO import Select

# logging settings
logger = logging.getLogger(__name__)

def _str2bool(v):
    """argparse type for an optional-value bool flag (--flag / --flag true/false)."""
    if isinstance(v, bool):
        return v
    s = str(v).lower()
    if s in ("1", "true", "yes", "on", "t"):
        return True
    if s in ("0", "false", "no", "off", "f"):
        return False
    from argparse import ArgumentTypeError
    raise ArgumentTypeError(f"expected true/false, got {v!r}")

# from jax 0.4.6, jax._src.lib.xla_bridge moved to jax._src.xla_bridge
# suppress warnings: Unable to initialize backend 'rocm' or 'tpu'
logging.getLogger('jax._src.xla_bridge').addFilter(lambda _: False) # jax >=0.4.6
logging.getLogger('jax._src.lib.xla_bridge').addFilter(lambda _: False) # jax < 0.4.5


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

def mk_hhsearch_single_entry_db(cif_file: Path, dbdir_cache_path: str):
    dbdir = Path(dbdir_cache_path)
    dbdir.mkdir(parents=True, exist_ok=True)
    tmp_cif_path = str(dbdir_cache_path) + "/1dmy.cif"
    shutil.copy2(cif_file, tmp_cif_path)
    # clear internal AF2 cache of mmCIF files
    from colabfold.backend import get_current_backend
    if type(get_current_backend()).__name__ == "AF2Backend":
        from alphafold.data import templates
        templates._read_file.cache_clear()
    cif_file = Path(tmp_cif_path)
    pdb70_db_files = dbdir.glob("pdb70*")
    for f in pdb70_db_files:
        os.remove(f)

    with open(dbdir.joinpath("pdb70_a3m.ffdata"), "w") as a3m, open(
        dbdir.joinpath("pdb70_cs219.ffindex"), "w"
    ) as cs219_index, open(
        dbdir.joinpath("pdb70_a3m.ffindex"), "w"
    ) as a3m_index, open(
        dbdir.joinpath("pdb70_cs219.ffdata"), "w"
    ) as cs219:
        n = 1000000
        index_offset = 0
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
                    restype_3to1.get(res.resname, "X")
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
                        restype_3to1.get(res.resname, "X")
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
    max_template_date="2100-01-01",
    max_template_hits=20,
) -> Tuple[
    Optional[List[str]],
    Optional[List[str]],
    List[str],
    List[int],
    List[Optional[Tuple[str, Any]]],
]:
    from colabfold.colabfold import run_mmseqs2

    use_env = msa_mode == "mmseqs2_uniref_env" or msa_mode == "mmseqs2_uniref_env_envpair"
    use_envpair = msa_mode == "mmseqs2_uniref_env_envpair"
    if isinstance(query_sequences, str): query_sequences = [query_sequences]

    # an alphafold3 job of nucleic acids and ligands only has nothing to search for
    if not query_sequences:
        return None, None, [], [], []

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

    # get templates
    template_results = []
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
            template_results = [None] * len(query_seqs_unique)
        else:
            for index in range(0, len(query_seqs_unique)):
                if template_paths[index] is not None:
                    template_results.append(
                        (a3m_lines_mmseqs2[index], template_paths[index])
                    )
                else:
                    template_results.append(None)
                    logger.info(f"Sequence {index} found no templates")
    else:
        template_results = [None] * len(query_seqs_unique)

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
        template_results,
    )


def normalize_a3m(lines: list[str]) -> list[str]:
    out = []
    i = 0

    # keep meta header
    if lines and lines[0].startswith("#"):
        out.append(lines[0].rstrip("\n"))
        i = 1

    header = None
    seq_chunks = []
    while i < len(lines):
        line = lines[i].strip()
        i += 1
        if not line:
            continue
        if line.startswith(">"):
            if header is not None:
                out.append(header)
                out.append("".join(seq_chunks))
            header = line
            seq_chunks = []
        else:
            # remove all whitespace inside sequence lines
            seq_chunks.append("".join(line.split()))
    if header is not None:
        out.append(header)
        out.append("".join(seq_chunks))

    return out

def unserialize_msa(
    a3m_lines: List[str], query_sequence: Union[List[str], str]
) -> Tuple[
    Optional[List[str]],
    Optional[List[str]],
    List[str],
    List[int],
    List[Dict[str, Any]],
]:
    if not a3m_lines:
        raise ValueError("no a3m to unserialize")
    a3m_lines = a3m_lines[0].replace("\x00", "").splitlines()
    a3m_lines = normalize_a3m(a3m_lines)
    if not a3m_lines:
        raise ValueError("this a3m has no lines; it needs at least a query and its sequence")
    if not a3m_lines[0].startswith("#") or len(a3m_lines[0][1:].split("\t")) != 2:
        assert isinstance(query_sequence, str)
        return (
            ["\n".join(a3m_lines)],
            None,
            [query_sequence],
            [1],
            [None],
        )

    if len(a3m_lines) < 3:
        raise ValueError(f"Unknown file format a3m")
    tab_sep_entries = a3m_lines[0][1:].split("\t")
    query_seq_len = tab_sep_entries[0].split(",")
    query_seq_len = list(map(int, query_seq_len))
    query_seqs_cardinality = tab_sep_entries[1].split(",")
    query_seqs_cardinality = list(map(int, query_seqs_cardinality))
    num_chains = len(query_seq_len)
    is_homooligomer = (
        True if num_chains == 1 and query_seqs_cardinality[0] > 1 else False
    )
    is_single_protein = (
        True if num_chains == 1 and query_seqs_cardinality[0] == 1 else False
    )

    query_seqs_unique = []
    qcat = a3m_lines[2]
    prev = 0
    for qlen in query_seq_len:
        nxt = prev + max(int(qlen), 0)
        query_seqs_unique.append(qcat[prev:nxt])
        prev = nxt

    paired_chunks = [[] for _ in range(num_chains)]
    unpaired_chunks = [[] for _ in range(num_chains)]
    already_in = set()

    qlens_local = query_seq_len
    def _split_by_aln(s):
        segments = [""] * num_chains
        has_aa = [False] * num_chains
        seg_idx = 0
        aln_count = 0
        start = 0
        curr_has_aa = False

        dash = '-'
        A, Z = 'A', 'Z'
        a, z = 'a', 'z'
        n = len(s)
        i = 0
        while i < n and seg_idx < num_chains:
            c = s[i]
            # non-lowercase are aligned columns
            if not (a <= c <= z):
                if c != dash and (A <= c <= Z):
                    curr_has_aa = True
                aln_count += 1
                if aln_count == qlens_local[seg_idx]:
                    # close current segment
                    segments[seg_idx] = s[start:i+1]
                    has_aa[seg_idx] = curr_has_aa
                    seg_idx += 1
                    aln_count = 0
                    start = i + 1
                    curr_has_aa = False
            i += 1
        return segments, has_aa

    i = 1
    end = len(a3m_lines)
    while i + 1 < end:
        header = a3m_lines[i]
        seq = a3m_lines[i + 1]
        i += 2

        key = (header, seq)
        if key in already_in:
            continue
        already_in.add(key)

        segments, has_aa = _split_by_aln(seq)

        # Paired if multi-chain (not single protein), not homo-oligomer, >=2 segments have AA
        if (not is_single_protein) and (not is_homooligomer) and (sum(has_aa) > 1):
            header_no_gt = header.replace(">", "")
            header_fields = header_no_gt.split("\t")
            for j, seg in enumerate(segments):
                label = header_fields[j] if j < len(header_fields) else (header_fields[-1] if header_fields else "")
                pc = paired_chunks[j]
                pc.append(">")
                pc.append(label)
                pc.append("\n")
                pc.append(seg)
                pc.append("\n")
        else:
            for j, seg in enumerate(segments):
                if has_aa[j]:
                    uc = unpaired_chunks[j]
                    uc.append(header)
                    uc.append("\n")
                    uc.append(seg)
                    uc.append("\n")

    if is_homooligomer:
        num = 101
        count = max(query_seqs_cardinality[0], 0)
        q = query_seqs_unique[0] if query_seqs_unique else ""
        paired_msa = [f">{num + k}\n{q}\n" for k in range(count)]
    else:
        if is_single_protein:
            paired_msa = None
        else:
            paired_msa = ["".join(ch) for ch in paired_chunks]


    unpaired_msa = ["".join(ch) for ch in unpaired_chunks]
    template_results = [None] * len(query_seqs_unique)

    return (
        unpaired_msa,
        paired_msa,
        query_seqs_unique,
        query_seqs_cardinality,
        template_results,
    )

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
    queries: List[Tuple[str, Union[str, List[str]], Optional[List[str]], Optional[List[Tuple[str, str,int]]]]],
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
    custom_template_cache_path: str = None,
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
    skip_output: List[str] = [],
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
    calc_extra_ptm: bool = False,
    use_probs_extra: bool = True,
    max_template_date: str = "2100-01-01",
    max_template_hits: int = 20,
    **kwargs
):
    # check what device is available
    try:
        # check if TPU is available
        from tpu_info import device
        if len(device.get_local_chips()) > 0:
            import jax.tools.colab_tpu
            jax.tools.colab_tpu.setup_tpu()
            logger.info('Running on TPU')
            use_gpu_relax = False
    except:
        from jax import local_devices
        if local_devices()[0].platform == 'cpu':
            logger.info("WARNING: no GPU detected, will be using CPU")
            use_gpu_relax = False
        else:
            logger.info('Running on GPU')

    from colabfold.backend import RunOptions, get_backend

    data_dir = Path(data_dir)
    result_dir = Path(result_dir)
    result_dir.mkdir(exist_ok=True)
    model_type = set_model_type(is_complex, model_type)
    backend = get_backend(model_type, data_dir)

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
    use_pallas            = kwargs.pop("use_pallas", False)  # old name
    use_fast_kernels      = kwargs.pop("use_fast_kernels", use_pallas)
    use_esm               = kwargs.pop("use_esm", False)
    kernel_backend        = kwargs.pop("kernel_backend", "auto")
    compile_mode          = kwargs.pop("compile_mode", "tuned")
    num_diffusion_samples = kwargs.pop("num_diffusion_samples", 5)
    af3_pairing           = kwargs.pop("af3_pairing", "colabfold")
    weights_precision     = kwargs.pop("weights_precision", "int8")
    buckets               = kwargs.pop("buckets", None)
    accept_alphafold3_terms = kwargs.pop("accept_alphafold3_terms", False)
    if use_fast_kernels and not use_bfloat16:
        raise ValueError("--use-fast-kernels needs half precision, not use_bfloat16=False")
    max_msa               = kwargs.pop("max_msa",None)
    if max_msa is not None:
        max_seq, max_extra_seq = [int(x) for x in max_msa.split(":")]

    if kwargs.pop("use_amber", False) and num_relax == 0:
        num_relax = num_models * num_seeds

    if len(kwargs) > 0:
        print(f"WARNING: the following options are not being used: {kwargs}")

    # decide how to rank outputs
    if is_af3_model(model_type):
        if rank_by == "auto":
            rank_by = "ranking_score"
    else:
        if rank_by == "auto":
            rank_by = "multimer" if is_complex else "plddt"
        if "ptm" not in model_type and "multimer" not in model_type:
            rank_by = "plddt"

    # added for actifptm calculation
    if not is_complex and calc_extra_ptm:
        logger.info("Calculating extra pTM is not supported for single chain prediction, skipping it.")
        calc_extra_ptm = False

    # get max length
    max_len = 0
    max_num = 0
    for _, query_sequence, _, _ in queries:
        N = 1 if isinstance(query_sequence,str) else len(query_sequence)
        L = len("".join(query_sequence))
        if L > max_len: max_len = L
        if N > max_num: max_num = N

    # sort model order
    model_order.sort()

    # initial guess
    if initial_guess is not None:
        logger.info(f'Using initial guess: {initial_guess}')

    # Run options handed to the backend
    opts = RunOptions(
        model_type=model_type,
        num_models=num_models,
        num_seeds=num_seeds,
        num_recycles=num_recycles,
        recycle_early_stop_tolerance=recycle_early_stop_tolerance,
        use_templates=use_templates,
        max_template_date=max_template_date,
        max_template_hits=max_template_hits,
        rank_by=rank_by,
        stop_at_score=stop_at_score,
        random_seed=random_seed,
        initial_guess=initial_guess,
        num_relax=num_relax,
        relax_max_iterations=relax_max_iterations,
        relax_tolerance=relax_tolerance,
        relax_stiffness=relax_stiffness,
        relax_max_outer_iterations=relax_max_outer_iterations,
        use_gpu_relax=use_gpu_relax,
        save_all=save_all,
        save_single_representations=save_single_representations,
        save_pair_representations=save_pair_representations,
        save_recycles=save_recycles,
        data_dir=data_dir,
        max_seq=max_seq,
        max_extra_seq=max_extra_seq,
        backend_opts={
            "model_order": model_order,
            "num_ensemble": num_ensemble,
            "use_dropout": use_dropout,
            "use_cluster_profile": use_cluster_profile,
            "use_fuse": use_fuse,
            "use_bfloat16": use_bfloat16,
            "use_fast_kernels": use_fast_kernels,
            "use_esm": use_esm,
            "kernel_backend": kernel_backend,
            "compile_mode": compile_mode,
            "recompile_padding": recompile_padding,
            "calc_extra_ptm": calc_extra_ptm,
            "use_probs_extra": use_probs_extra,
            "num_diffusion_samples": num_diffusion_samples,
            "af3_pairing": af3_pairing,
            "weights_precision": weights_precision,
            "buckets": buckets,
            "accept_alphafold3_terms": accept_alphafold3_terms,
        },
    )
    backend.configure(
        opts, max_len=max_len, max_num=max_num, num_queries=len(queries),
        msa_mode=msa_mode, is_complex=is_complex, use_templates=use_templates,
    )

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
        "initial_guess": initial_guess,
        "keep_existing_results": keep_existing_results,
        "rank_by": rank_by,
        "pair_mode": pair_mode,
        "pairing_strategy": pairing_strategy,
        "host_url": host_url,
        "user_agent": user_agent,
        "stop_at_score": stop_at_score,
        "random_seed": random_seed,
        "num_seeds": num_seeds,
        "commit": get_commit(),
        "version": importlib_metadata.version("colabfold"),
        "max_template_date": max_template_date,
        "max_template_hits": max_template_hits,
    }
    refuse_mixed_msas(queries, model_type)
    dropped = dropped_entities(queries, model_type)
    warn_about_dropped_entities(dropped, model_type)
    if dropped:
        config["dropped_entities"] = dropped
    config.update(backend.config_dict(opts))
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

    ranks, metrics = [],[]
    job_number = 0
    for job_number, (raw_jobname, query_sequence, a3m_lines, custom_template_path_per_entry) in enumerate(queries):

        if use_templates and custom_template_path_per_entry is not None and isinstance(custom_template_path_per_entry, Path):
            mk_hhsearch_single_entry_db(custom_template_path_per_entry, custom_template_cache_path)
            custom_template_path = custom_template_cache_path

        if jobname_prefix is not None:
            # pad job number based on number of queries
            fill = len(str(len(queries)))
            jobname = safe_filename(jobname_prefix) + "_" + str(job_number).zfill(fill)
            job_number += 1
        else:
            jobname = safe_filename(raw_jobname)

        # nothing survives the drop: skip this job rather than the whole batch
        if not query_sequence and not is_af3_model(model_type):
            logger.error(f"{jobname}: {model_type} folds protein chains only")
            continue

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

        # an alphafold3 JSON knows its own chains, and its nucleic acids count too
        fold_input = getattr(custom_template_path_per_entry, "fold_input", None)
        chain_lengths = None
        own_msa = False
        if fold_input is not None and is_af3_model(model_type):
            from colabfold.alphafold3.input import msa_state, polymer_lengths
            chain_lengths = polymer_lengths(fold_input)
            own_msa = msa_state(fold_input) == "provided"
        seq_len = sum(chain_lengths) if chain_lengths is not None else len("".join(query_sequence))
        logger.info(f"Query {job_number + 1}/{len(queries)}: {jobname} (length {seq_len})")

        ###########################################
        # generate MSA (a3m_lines) and templates
        ###########################################
        if own_msa:
            logger.info(f"{jobname}: folding with the MSAs in the input file")
        try:
            pickled_msa_and_templates = result_dir.joinpath(f"{jobname}.pickle")
            if pickled_msa_and_templates.is_file():
                with open(pickled_msa_and_templates, 'rb') as f:
                    (unpaired_msa, paired_msa, query_seqs_unique, query_seqs_cardinality, template_results) = pickle.load(f)
                logger.info(f"Loaded {pickled_msa_and_templates}")

            else:
                # the file brought its own MSAs, so only templates can still need a search
                if own_msa:
                    from colabfold.alphafold3.input import msas_of, sequences_of

                    query_seqs_unique, query_seqs_cardinality = sequences_of(fold_input)
                    unpaired_msa, paired_msa = msas_of(fold_input)
                    template_results = []
                    if use_templates:
                        (_, _, _, _, template_results) \
                        = get_msa_and_templates(
                            jobname, query_seqs_unique, unpaired_msa, result_dir, 'single_sequence',
                            use_templates, custom_template_path, pair_mode, pairing_strategy,
                            host_url, user_agent, max_template_date=max_template_date,
                            max_template_hits=max_template_hits,
                        )

                elif a3m_lines is None:
                    (unpaired_msa, paired_msa, query_seqs_unique, query_seqs_cardinality, template_results) \
                    = get_msa_and_templates(
                        jobname, query_sequence, a3m_lines, result_dir, msa_mode, use_templates,
                        custom_template_path, pair_mode, pairing_strategy, host_url, user_agent,
                        max_template_date=max_template_date, max_template_hits=max_template_hits,
                    )

                elif a3m_lines is not None:
                    if(isinstance(a3m_lines, Path)):
                        a3m_lines = [a3m_lines.read_text()]
                    (unpaired_msa, paired_msa, query_seqs_unique, query_seqs_cardinality, template_results) \
                    = unserialize_msa(a3m_lines, query_sequence)
                    if use_templates:
                        (_, _, _, _, template_results) \
                            = get_msa_and_templates(
                                jobname, query_seqs_unique, unpaired_msa, result_dir, 'single_sequence', use_templates,
                                custom_template_path, pair_mode, pairing_strategy, host_url, user_agent,
                                max_template_date=max_template_date, max_template_hits=max_template_hits,
                            )

                if num_models == 0:
                    with open(pickled_msa_and_templates, 'wb') as f:
                        pickle.dump((unpaired_msa, paired_msa, query_seqs_unique, query_seqs_cardinality, template_results), f)
                    logger.info(f"Saved {pickled_msa_and_templates}")

            # a3m input and pickles bypass get_msa_and_templates, and num_extra_msa=1 would then sample a random homolog per seed
            if msa_mode == "single_sequence" and unpaired_msa is not None:
                unpaired_msa = [f">{101 + i}\n{seq}" for i, seq in enumerate(query_seqs_unique)]
                paired_msa = None

            # save a3m
            if not 'msa' in skip_output and query_seqs_unique:
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
            = backend.featurize(query_seqs_unique, query_seqs_cardinality, unpaired_msa, paired_msa,
                                template_results, sum(query_seqs_cardinality) > 1, opts,
                                extras=custom_template_path_per_entry)

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
        if not 'plots' in skip_output:
            msa_plot = backend.plot_msa(feature_dict, dpi=dpi)
            if msa_plot is not None:
                coverage_png = result_dir.joinpath(f"{jobname}_coverage.png")
                msa_plot.savefig(str(coverage_png), bbox_inches='tight')
                msa_plot.close()
                result_files.append(coverage_png)

        if use_templates:
            templates_file = result_dir.joinpath(f"{jobname}_template_domain_names.json")
            templates_file.write_text(json.dumps(domain_names))
            result_files.append(templates_file)

        if query_seqs_unique:
            result_files.append(result_dir.joinpath(jobname + ".a3m"))
        result_files += [bibtex_file, config_out_file]

        ######################
        # predict structures
        ######################
        if num_models > 0:
            try:
                # get list of lengths
                query_sequence_len_array = chain_lengths if chain_lengths is not None else sum(
                    [[len(x)] * y for x,y in zip(query_seqs_unique, query_seqs_cardinality)],[])

                results = backend.predict(
                    prefix=jobname,
                    result_dir=result_dir,
                    model_input=feature_dict,
                    is_complex=is_complex,
                    sequences_lengths=query_sequence_len_array,
                    opts=opts,
                    prediction_callback=prediction_callback,
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
            if not 'pae_json' in skip_output:
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
                    if not 'plots' in skip_output:
                        from colabfold.colabfold import plot_paes
                        paes_plot = plot_paes([np.asarray(x["pae"]) for x in scores],
                            Ls=query_sequence_len_array, dpi=dpi)
                        pae_png = result_dir.joinpath(f"{jobname}_pae.png")
                        paes_plot.savefig(str(pae_png), bbox_inches='tight')
                        paes_plot.close()
                        result_files.append(pae_png)

                    # make pairwise interface metric plots and chainwise ptm plot
                    if calc_extra_ptm:
                        ext_metric_png = result_dir.joinpath(f"{jobname}_ext_metrics.png")
                        backend.plot_extra_metrics(scores, ext_metric_png)

                # make pLDDT plot
                if not 'plots' in skip_output:
                    from colabfold.colabfold import plot_plddts
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

def generate_af3_input(
    queries: List[Tuple[str, Union[str, List[str]], Optional[List[str]], Optional[List[Tuple[str, str, int]]]]],
    result_dir: Union[str, Path],
    msa_mode: str = "mmseqs2_uniref_env",
    pair_mode: str = "unpaired_paired",
    pairing_strategy: str = "greedy",
    use_templates: bool = False,
    custom_template_path: str = None,
    jobname_prefix: Optional[str] = None,
    host_url: str = DEFAULT_API_SERVER,
    user_agent: str = "",
    max_template_date: str = "2100-01-01",
    max_template_hits: int = 20,
    # is_complex: bool,
    # model_type: str = "auto",
):
    result_dir = Path(result_dir) 
    result_dir.mkdir(exist_ok=True)

    job_number = 0

    for job_number, (raw_jobname, query_sequences, a3m_lines, other_molecules) in enumerate(queries):
        if jobname_prefix is not None:
            # pad job number based on number of queries
            fill = len(str(len(queries)))
            jobname = safe_filename(jobname_prefix) + "_" + str(job_number).zfill(fill)
            # job_number += 1 # Why add?
        else:
            jobname = safe_filename(raw_jobname)

        ###########################################
        # generate MSA (a3m_lines) and templates
        ###########################################
        try:
            if a3m_lines is None:
                (unpaired_msa, paired_msa, query_seqs_unique, query_seqs_cardinality, template_results) \
                = get_msa_and_templates(
                    jobname, query_sequences, a3m_lines, result_dir, msa_mode, use_templates,
                    custom_template_path, pair_mode, pairing_strategy, host_url, user_agent,
                    max_template_date=max_template_date, max_template_hits=max_template_hits,
                )

            elif a3m_lines is not None:
                (unpaired_msa, paired_msa, query_seqs_unique, query_seqs_cardinality, template_results) \
                = unserialize_msa(a3m_lines, query_sequences)
                # if use_templates:
                #     (_, _, _, _, template_results) \
                #         = get_msa_and_templates(
                #               jobname, query_seqs_unique, unpaired_msa, result_dir, 'single_sequence', use_templates,
                #               custom_template_path, pair_mode, pairing_strategy, host_url, user_agent,
                #               max_template_date=max_template_date, max_template_hits=max_template_hits,
                #           )

            # save json
            af3 = AF3Utils(jobname, query_seqs_unique, query_seqs_cardinality, unpaired_msa, paired_msa, other_molecules)
            with open(result_dir.joinpath(f"{jobname}.json"), "w") as f:
                f.write(json.dumps(af3.content, indent = 4))
                
            # save a3m
            msa = msa_to_str(unpaired_msa, paired_msa, query_seqs_unique, query_seqs_cardinality)
            result_dir.joinpath(f"{jobname}.a3m").write_text(msa)

        except Exception as e:
            logger.exception(f"Failed to generate AF3 input json for {jobname}: Could not get MSA/templates. {e}")
            continue

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
        "--custom-template-cache-path",
        type=str,
        default=None,
        help="Directory to generate temporary HHsearch databases for custom templates.",
    )
    msa_group.add_argument(
        "--max-template-date",
        type=str,
        default="2100-01-01",
        help="Maximum release date (YYYY-MM-DD) for templates to be considered."
    )
    msa_group.add_argument(
        "--max-template-hits",
        type=int,
        default=20,
        help="Maximum number of template hits to consider."
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
            *AF3_MODELS,
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
    pred_group.add_argument(
        "--calc-extra-ptm",
        default=False,
        action="store_true",
        help="Experimental: calculate pairwise metrics (ipTM and actifpTM), and also chain-wise pTM",
    )
    pred_group.add_argument(
        "--no-use-probs-extra",
        default=False,
        action="store_true",
        help="Experimental: instead of contact probabilities form use binary contacts for extra metrics calculation",
    )
    pred_group.add_argument("--data", help="Path to AlphaFold2 weights directory.")

    af3_group = parser.add_argument_group("AlphaFold3 arguments", "")
    af3_group.add_argument(
        "--num-diffusion-samples",
        help="Structures per seed, for the alphafold3 models only.",
        type=int,
        default=5,
    )
    af3_group.add_argument(
        "--use-esm",
        default=False,
        action="store_true",
        help="Fold with the model's protein language model. esmfold2 is a different model "
        "without it and chai1 takes it as an extra feature, while models that fold from an "
        "MSA ignore it. Fetched on first use.",
    )
    af3_group.add_argument(
        "--af3-pairing",
        help="How the alphafold3 models pair a complex MSA: keep ColabFold's row pairing, "
        "or let them pair by UniProt species id.",
        choices=["colabfold", "uniprot"],
        default="colabfold",
    )
    af3_group.add_argument(
        "--accept-alphafold3-terms",
        help="Download AlphaFold3's own parameters from Google DeepMind, under their terms of "
        "use, which do not allow commercial use. Every other model is fetched without this.",
        action="store_true",
    )
    af3_group.add_argument(
        "--buckets",
        help="Token counts the alphafold3 models pad to, so queries of mixed length compile "
        "once per bucket instead of once per length. Comma separated, increasing.",
        default=None,
    )
    af3_group.add_argument(
        "--weights-precision",
        help="Which published form of the alphafold3 weights to fetch. int8 is the same "
        "weights stored smaller and expanded on load, which leaves inference unchanged.",
        choices=["fp32", "int8"],
        default="int8",
    )
    af3_group.add_argument(
        "--af3-json",
        help="Generate input JSON for AlphaFold3 from the provided FASTA/A3M file.",
        action="store_true",
    )

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
    def comma_separated_list(arg_string):
        return [item.strip() for item in arg_string.split(',') if item.strip() in ['msa', 'plots', 'pae_json']]
    output_group.add_argument(
        "--skip-output",
        help="Comma-separated list of output types to skip: msa, plots, pae_json.",
        type=comma_separated_list,
        default="",
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
        help="Sort input queries by: none, length, msa_depth, random. "
        "Sorting by length or depth speeds up monomer or multimer prediction as models are recompiled less often.",
        type=str,
        default="length",
        choices=["none", "length", "msa_depth", "random"],
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
    adv_group.add_argument(
        "--use-fast-kernels",
        dest="use_fast_kernels",
        default=False,
        action=BooleanOptionalAction,
        help="Use fused kernels for faster prediction",
    )
    adv_group.add_argument(
        # old name: takes an optional value, so it must not precede the paths
        "--use-pallas",
        dest="use_fast_kernels",
        nargs="?",
        const=True,
        default=False,
        type=_str2bool,
        metavar="BOOL",
        help=SUPPRESS,
    )
    adv_group.add_argument(
        "--kernel-backend",
        choices=["auto", "pallas", "cuda_legacy"],
        default="auto",
        help="Fused kernels: auto (pallas on sm_80+, cuda_legacy below), pallas, "
             "or cuda_legacy (float16).",
    )
    adv_group.add_argument(
        "--compile-mode",
        choices=["fast", "tuned", "full"],
        default="tuned",
        help="Kernel autotuning effort, trading compile time for inference speed: "
        "'fast': cuBLAS only, fastest compile, ~4%% slower inference, single/one-off predictions. "
        "'tuned': fixed, tuned kernel shape compiled with cuBLAS and Pallas, fast compile and ~1%% off optimal. "
        "'full': unbounded Kernel autotune, can take tens of minutes to compile for longsequences  but optimal inference, can be worth it for large batches."
    )
    adv_group.add_argument(
        "--debug-logging",
        default=False,
        action="store_true",
        help="Enable debug message logging.",
    )

    args = parser.parse_args()

    if (args.custom_template_path is not None) and (args.pdb_hit_file is not None):
        raise RuntimeError("Arguments --pdb-hit-file and --custom-template-path cannot be used simultaneously.")
    # disable unified memory
    if args.disable_unified_memory:
        for k in ENV.keys():
            if k in os.environ: del os.environ[k]

    setup_logging(Path(args.results).joinpath("log.txt"), verbose=args.debug_logging)

    version = importlib_metadata.version("colabfold")
    commit = get_commit()
    if commit:
        version += f" ({commit})"

    logger.info(f"Running colabfold {version}")

    data_dir = Path(args.data or default_data_dir)

    if args.msa_only:
        args.num_models = 0

    if args.num_models > 0 and is_af3_model(args.model_type):
        from colabfold.alphafold3.tokamax_shim import install as install_shim

        # first: whatever imports alphafold3 binds tokamax, and the shim cannot reach it after
        install_shim(force=args.use_fast_kernels)

        from colabfold.alphafold3.models import resolve_model_name
        from colabfold.alphafold3.weights import ensure_ccd, ensure_weights

        ensure_ccd(data_dir)
        ensure_weights(resolve_model_name(args.model_type), data_dir=data_dir,
                       precision=args.weights_precision,
                       accept_terms=args.accept_alphafold3_terms)

    queries, is_complex = get_queries(args.input, args.sort_queries_by)

    has_per_entry_templates = any(isinstance(q[3], Path) for q in queries)
    if has_per_entry_templates and args.custom_template_cache_path is None:
        raise ValueError("--custom-template-cache-path must be set when using per-entry template paths in CSV input")
    if has_per_entry_templates and args.custom_template_path is not None:
        raise ValueError("--custom-template-path and per-entry template paths in CSV input cannot be used simultaneously")

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

    if args.num_models > 0 and not is_af3_model(model_type):
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

    # added for actifptm calculation
    use_probs_extra = False if args.no_use_probs_extra else True

    user_agent = f"colabfold/{version}"

    if args.af3_json:
        generate_af3_input(
            queries=queries,
            result_dir=args.results,
            msa_mode=args.msa_mode,
            pair_mode=args.pair_mode,
            pairing_strategy=args.pair_strategy,
            use_templates=args.templates,
            custom_template_path=args.custom_template_path,
            jobname_prefix=args.jobname_prefix,
            host_url=args.host_url,
            user_agent=user_agent,
            max_template_date=args.max_template_date,
            max_template_hits=args.max_template_hits,
            # extra_molecules=extra_molecules,
        )
        return
    
    run(
        queries=queries,
        result_dir=args.results,
        use_templates=args.templates,
        custom_template_path=args.custom_template_path,
        custom_template_cache_path=args.custom_template_cache_path,
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
        skip_output=args.skip_output,
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
        calc_extra_ptm=args.calc_extra_ptm,
        use_probs_extra=use_probs_extra,
        max_template_date=args.max_template_date,
        max_template_hits=args.max_template_hits,
        use_fast_kernels=args.use_fast_kernels,
        kernel_backend=args.kernel_backend,
        compile_mode=args.compile_mode,
        num_diffusion_samples=args.num_diffusion_samples,
        use_esm=args.use_esm,
        af3_pairing=args.af3_pairing,
        weights_precision=args.weights_precision,
        buckets=[int(b) for b in args.buckets.split(",")] if args.buckets else None,
        accept_alphafold3_terms=args.accept_alphafold3_terms,
    )

if __name__ == "__main__":
    main()
