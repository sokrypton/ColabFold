import json
import logging
import math
import random
import sys
import time
import zipfile
from argparse import ArgumentParser
from pathlib import Path
from typing import Any, Callable, Dict, List, Optional, Tuple, Union

import haiku
import importlib_metadata
import numpy as np
import pandas
from alphafold.notebooks.notebook_utils import get_pae_json
from jax.lib import xla_bridge
from numpy import ndarray

try:
    import alphafold
except ModuleNotFoundError:
    raise RuntimeError(
        "\n\nalphafold is not installed. Please run `pip install colabfold[alphafold]`\n"
    )

from alphafold.common import protein
from alphafold.common import confidence

from alphafold.common.protein import Protein
from alphafold.data import (
    feature_processing,
    msa_pairing,
    pipeline,
    pipeline_multimer,
    templates,
)
from alphafold.data.tools import hhsearch
from alphafold.model import model

from colabfold.alphafold.models import load_models_and_params
from colabfold.alphafold.msa import make_fixed_size
from colabfold.citations import write_bibtex
from colabfold.colabfold import chain_break, plot_paes, plot_plddts, run_mmseqs2
from colabfold.download import default_data_dir, download_alphafold_params
from colabfold.plot import plot_msa
from colabfold.utils import (
    ACCEPT_DEFAULT_TERMS,
    DEFAULT_API_SERVER,
    NO_GPU_FOUND,
    get_commit,
    safe_filename,
    setup_logging,
)

logger = logging.getLogger(__name__)


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


def batch_input(
    input_features: model.features.FeatureDict,
    model_runner: model.RunModel,
    model_name: str,
    crop_len: int,
    use_templates: bool,
) -> model.features.FeatureDict:
    model_config = model_runner.config
    eval_cfg = model_config.data.eval
    crop_feats = {k: [None] + v for k, v in dict(eval_cfg.feat).items()}

    # templates models
    if (model_name == "model_1" or model_name == "model_2") and use_templates:
        pad_msa_clusters = eval_cfg.max_msa_clusters - eval_cfg.max_templates
    else:
        pad_msa_clusters = eval_cfg.max_msa_clusters

    max_msa_clusters = pad_msa_clusters

    # let's try pad (num_res + X)
    input_fix = make_fixed_size(
        input_features,
        crop_feats,
        msa_cluster_size=max_msa_clusters,  # true_msa (4, 512, 68)
        extra_msa_size=5120,  # extra_msa (4, 5120, 68)
        num_res=crop_len,  # aatype (4, 68)
        num_templates=4,
    )  # template_mask (4, 4) second value
    return input_fix


def predict_structure(
    prefix: str,
    result_dir: Path,
    feature_dict: Dict[str, Any],
    is_complex: bool,
    use_templates: bool,
    sequences_lengths: List[int],
    crop_len: int,
    model_type: str,
    model_runner_and_params: List[Tuple[str, model.RunModel, haiku.Params]],
    do_relax: bool = False,
    rank_by: str = "auto",
    random_seed: int = 0,
    stop_at_score: float = 100,
    prediction_callback: Callable[[Any, Any, Any, Any], Any] = None,
):
    """Predicts structure using AlphaFold for the given sequence."""
    if rank_by == "auto":
        # score complexes by ptmscore and sequences by plddt
        rank_by = "plddt" if not is_complex else "ptmscore"
        rank_by = (
            "multimer"
            if is_complex and model_type == "AlphaFold2-multimer"
            else rank_by
        )

    plddts, paes, ptmscore, iptmscore = [], [], [], []
    max_paes = []
    unrelaxed_pdb_lines = []
    relaxed_pdb_lines = []
    prediction_times = []
    seq_len = sum(sequences_lengths)

    model_names = []
    for (model_name, model_runner, params) in model_runner_and_params:
        logger.info(f"Running {model_name}")
        model_names.append(model_name)
        # swap params to avoid recompiling
        # note: models 1,2 have diff number of params compared to models 3,4,5 (this was handled on construction)
        model_runner.params = params

        processed_feature_dict = model_runner.process_features(
            feature_dict, random_seed=random_seed
        )
        if not is_complex:
            input_features = batch_input(
                processed_feature_dict,
                model_runner,
                model_name,
                crop_len,
                use_templates,
            )
        else:
            input_features = processed_feature_dict

        start = time.time()

        # The original alphafold only returns the prediction_result,
        # but our patched alphafold also returns a tuple (recycles,tol)
        prediction_result, recycles = model_runner.predict(input_features)

        prediction_time = time.time() - start
        prediction_times.append(prediction_time)

        mean_plddt = np.mean(prediction_result["plddt"][:seq_len])
        mean_ptm = prediction_result["ptm"]
        if rank_by == "plddt":
            mean_score = mean_plddt
        else:
            mean_score = mean_ptm

        if is_complex:
            logger.info(
                f"{model_name} took {prediction_time:.1f}s ({recycles[0]} recycles) "
                f"with pLDDT {mean_plddt:.3g} and ptmscore {mean_ptm:.3g}"
            )
        else:
            logger.info(
                f"{model_name} took {prediction_time:.1f}s ({recycles[0]} recycles) "
                f"with pLDDT {mean_plddt:.3g}"
            )
        final_atom_mask = prediction_result["structure_module"]["final_atom_mask"]
        b_factors = prediction_result["plddt"][:, None] * final_atom_mask
        if is_complex and model_type == "AlphaFold2-ptm":
            input_features["asym_id"] = feature_dict["asym_id"]
            input_features["aatype"] = input_features["aatype"][0]
            input_features["residue_index"] = input_features["residue_index"][0]
            curr_residue_index = 1
            res_index_array = input_features["residue_index"].copy()
            res_index_array[0] = 0
            for i in range(1, input_features["aatype"].shape[0]):
                if (
                    input_features["residue_index"][i]
                    - input_features["residue_index"][i - 1]
                ) > 1:
                    curr_residue_index = 0
                res_index_array[i] = curr_residue_index
                curr_residue_index += 1
            input_features["residue_index"] = res_index_array

        unrelaxed_protein = protein.from_prediction(
            features=input_features,
            result=prediction_result,
            b_factors=b_factors,
            remove_leading_feature_dimension=not is_complex,
        )

        if prediction_callback is not None:
            prediction_callback(
                unrelaxed_protein, sequences_lengths, prediction_result, input_features
            )

        unrelaxed_pdb_lines.append(protein.to_pdb(unrelaxed_protein))
        plddts.append(prediction_result["plddt"][:seq_len])
        ptmscore.append(prediction_result["ptm"])
        if model_type == "AlphaFold2-multimer":
            iptmscore.append(prediction_result["iptm"])
        max_paes.append(prediction_result["max_predicted_aligned_error"].item())
        paes_res = []

        for i in range(seq_len):
            paes_res.append(prediction_result["predicted_aligned_error"][i][:seq_len])
        paes.append(paes_res)
        if do_relax:
            from alphafold.common import residue_constants
            from alphafold.relax import relax

            # Hack so that we don't need to download into the alphafold package itself
            residue_constants.stereo_chemical_props_path = "stereo_chemical_props.txt"

            # Remove the padding because unlike to_pdb() amber doesn't handle that
            remove_padding_mask = np.array(unrelaxed_protein.atom_mask.sum(axis=-1) > 0)
            unrelaxed_protein = Protein(
                atom_mask=unrelaxed_protein.atom_mask[remove_padding_mask],
                atom_positions=unrelaxed_protein.atom_positions[remove_padding_mask],
                aatype=unrelaxed_protein.aatype[remove_padding_mask],
                residue_index=unrelaxed_protein.residue_index[remove_padding_mask],
                b_factors=unrelaxed_protein.b_factors[remove_padding_mask],
                chain_index=unrelaxed_protein.chain_index[remove_padding_mask],
            )

            # Relax the prediction.
            amber_relaxer = relax.AmberRelaxation(
                max_iterations=0,
                tolerance=2.39,
                stiffness=10.0,
                exclude_residues=[],
                max_outer_iterations=20,
            )
            relaxed_pdb_str, _, _ = amber_relaxer.process(prot=unrelaxed_protein)
            # TODO: Those aren't actually used in batch
            relaxed_pdb_lines.append(relaxed_pdb_str)
        # early stop criteria fulfilled
        if mean_score > stop_at_score:
            break
    # rerank models based on predicted lddt
    if rank_by == "ptmscore":
        model_rank = np.array(ptmscore).argsort()[::-1]
    elif rank_by == "multimer":
        rank_array = np.array(iptmscore) * 0.8 + np.array(ptmscore) * 0.2
        model_rank = rank_array.argsort()[::-1]
    else:
        model_rank = np.mean(plddts, -1).argsort()[::-1]
    out = {}
    logger.info(f"reranking models based on average {rank_by}")
    for n, key in enumerate(model_rank):
        unrelaxed_pdb_path = result_dir.joinpath(
            f"{prefix}_unrelaxed_rank_{n + 1}_{model_names[key]}.pdb"
        )
        unrelaxed_pdb_path.write_text(unrelaxed_pdb_lines[key])

        if do_relax:
            relaxed_pdb_path = result_dir.joinpath(
                f"{prefix}_relaxed_rank_{n + 1}_{model_names[key]}.pdb"
            )
            relaxed_pdb_path.write_text(relaxed_pdb_lines[key])

        # Write an easy-to-use format (PAE and plDDT)
        scores_file = result_dir.joinpath(
            f"{prefix}_unrelaxed_rank_{n + 1}_{model_names[key]}_scores.json"
        )
        with scores_file.open("w") as fp:
            # We use astype(np.float64) to prevent very long stringified floats from float imprecision
            scores = {
                "max_pae": max_paes[key],
                "pae": np.around(np.asarray(paes[key]).astype(np.float64), 2).tolist(),
                "plddt": np.around(np.asarray(plddts[key]), 2).tolist(),
                "ptm": np.around(ptmscore[key], 2).item(),
            }
            json.dump(scores, fp)

        out[key] = {
            "plddt": np.asarray(plddts[key]),
            "pae": np.asarray(paes[key]),
            "max_pae": max_paes[key],
            "pTMscore": ptmscore[key],
            "model_name": model_names[key],
        }
    return out, model_rank


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
        elif input_path.suffix == ".fasta":
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
            if len(seqs) > 1 and file.suffix == ".fasta":
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
        queries.sort(key=lambda t: len(t[1]))
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
    pair_mode: str,
    host_url: str = DEFAULT_API_SERVER,
) -> Tuple[
    Optional[List[str]], Optional[List[str]], List[str], List[int], List[Dict[str, Any]]
]:
    use_env = msa_mode == "MMseqs2 (UniRef+Environmental)"
    # remove duplicates before searching
    query_sequences = (
        [query_sequences] if isinstance(query_sequences, str) else query_sequences
    )
    query_seqs_unique = []
    for x in query_sequences:
        if x not in query_seqs_unique:
            query_seqs_unique.append(x)
    query_seqs_cardinality = [0] * len(query_seqs_unique)
    for seq in query_sequences:
        seq_idx = query_seqs_unique.index(seq)
        query_seqs_cardinality[seq_idx] += 1

    template_features = []
    if use_templates:
        a3m_lines_mmseqs2, template_paths = run_mmseqs2(
            query_seqs_unique,
            str(result_dir.joinpath(jobname)),
            use_env,
            use_templates=True,
            host_url=host_url,
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
                    logger.info(
                        f"Sequence {index} found templates: {template_feature['template_domain_names']}"
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

    if pair_mode == "none" or pair_mode == "unpaired" or pair_mode == "unpaired+paired":
        if msa_mode == "single_sequence":
            a3m_lines = []
            num = 101
            for i, seq in enumerate(query_seqs_unique):
                a3m_lines.append(">" + str(num + i) + "\n" + seq)
        else:
            # find normal a3ms
            a3m_lines = run_mmseqs2(
                query_seqs_unique,
                str(result_dir.joinpath(jobname)),
                use_env,
                use_pairing=False,
                host_url=host_url,
            )
    else:
        a3m_lines = None

    if pair_mode == "paired" or pair_mode == "unpaired+paired":
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
                paired_a3m_lines.append(
                    ">" + str(num + i) + "\n" + query_seqs_unique[0] + "\n"
                )
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
    features_for_chain: Dict[str, Dict[str, ndarray]]
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
) -> Dict[str, Any]:
    input_feature = {}
    if is_complex and model_type == "AlphaFold2-ptm":
        a3m_lines = pair_msa(
            query_seqs_unique, query_seqs_cardinality, paired_msa, unpaired_msa
        )
        total_sequence = ""
        Ls = []
        for sequence_index, sequence in enumerate(query_seqs_unique):
            for cardinality in range(0, query_seqs_cardinality[sequence_index]):
                total_sequence += sequence
                Ls.append(len(sequence))

        input_feature = build_monomer_feature(
            total_sequence, a3m_lines, mk_mock_template(total_sequence)
        )
        input_feature["residue_index"] = chain_break(input_feature["residue_index"], Ls)
        input_feature["asym_id"] = np.array(
            [int(n) for n, l in enumerate(Ls) for _ in range(0, l)]
        )
    else:
        features_for_chain = {}
        chain_cnt = 0
        for sequence_index, sequence in enumerate(query_seqs_unique):
            for cardinality in range(0, query_seqs_cardinality[sequence_index]):
                if unpaired_msa is None:
                    input_msa = ">" + str(101 + sequence_index) + "\n" + sequence
                else:
                    input_msa = unpaired_msa[sequence_index]

                feature_dict = build_monomer_feature(
                    sequence, input_msa, template_features[sequence_index]
                )
                if is_complex:
                    all_seq_features = build_multimer_feature(
                        paired_msa[sequence_index]
                    )
                    feature_dict.update(all_seq_features)
                features_for_chain[protein.PDB_CHAIN_IDS[chain_cnt]] = feature_dict
                chain_cnt += 1

        # Do further feature post-processing depending on the model type.
        if not is_complex:
            input_feature = features_for_chain[protein.PDB_CHAIN_IDS[0]]
        elif model_type == "AlphaFold2-multimer":
            input_feature = process_multimer_features(features_for_chain)
    return input_feature


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


def run(
    queries: List[Tuple[str, Union[str, List[str]], Optional[List[str]]]],
    result_dir: Union[str, Path],
    num_models: int,
    num_recycles: int,
    model_order: List[int],
    is_complex: bool,
    model_type: str = "auto",
    msa_mode: str = "MMseqs2 (UniRef+Environmental)",
    use_templates: bool = False,
    use_amber: bool = False,
    keep_existing_results: bool = True,
    rank_by: str = "auto",
    pair_mode: str = "unpaired+paired",
    data_dir: Union[str, Path] = default_data_dir,
    host_url: str = DEFAULT_API_SERVER,
    stop_at_score: float = 100,
    recompile_padding: float = 1.1,
    recompile_all_models: bool = False,
    zip_results: bool = False,
    prediction_callback: Callable[[Any, Any, Any, Any], Any] = None,
):
    version = importlib_metadata.version("colabfold")
    commit = get_commit()
    if commit:
        version += f" ({commit})"

    logger.info(f"Running colabfold {version}")

    data_dir = Path(data_dir)
    result_dir = Path(result_dir)
    result_dir.mkdir(exist_ok=True)
    model_type = set_model_type(is_complex, model_type)
    if model_type == "AlphaFold2-multimer":
        model_extension = "_multimer"
    elif model_type == "AlphaFold2-ptm":
        model_extension = "_ptm"
    else:
        raise ValueError(f"Unknown model_type {model_type}")

    if rank_by == "auto":
        # score complexes by ptmscore and sequences by plddt
        rank_by = "plddt" if not is_complex else "ptmscore"

    # Record the parameters of this run
    config = {
        "num_queries": len(queries),
        "use_templates": use_templates,
        "use_amber": use_amber,
        "msa_mode": msa_mode,
        "model_type": model_type,
        "num_models": num_models,
        "num_recycles": num_recycles,
        "model_order": model_order,
        "keep_existing_results": keep_existing_results,
        "rank_by": rank_by,
        "pair_mode": pair_mode,
        "host_url": host_url,
        "stop_at_score": stop_at_score,
        "recompile_padding": recompile_padding,
        "recompile_all_models": recompile_all_models,
        "commit": get_commit(),
        "version": importlib_metadata.version("colabfold"),
    }
    config_out_file = result_dir.joinpath("config.json")
    config_out_file.write_text(json.dumps(config, indent=4))
    use_env = msa_mode == "MMseqs2 (UniRef+Environmental)"
    use_msa = (
        msa_mode == "MMseqs2 (UniRef only)"
        or msa_mode == "MMseqs2 (UniRef+Environmental)"
    )

    bibtex_file = write_bibtex(
        model_type, use_msa, use_env, use_templates, use_amber, result_dir
    )

    model_runner_and_params = load_models_and_params(
        num_models,
        use_templates,
        num_recycles,
        model_order,
        model_extension,
        data_dir,
        recompile_all_models,
        stop_at_score=stop_at_score,
        rank_by=rank_by,
    )

    crop_len = 0
    for job_number, (raw_jobname, query_sequence, a3m_lines) in enumerate(queries):
        jobname = safe_filename(raw_jobname)
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

        query_sequence_len = (
            len(query_sequence)
            if isinstance(query_sequence, str)
            else sum(len(s) for s in query_sequence)
        )
        logger.info(
            f"Query {job_number + 1}/{len(queries)}: {jobname} (length {query_sequence_len})"
        )

        try:
            if a3m_lines is not None:
                (
                    unpaired_msa,
                    paired_msa,
                    query_seqs_unique,
                    query_seqs_cardinality,
                    template_features,
                ) = unserialize_msa(a3m_lines, query_sequence)
            else:
                (
                    unpaired_msa,
                    paired_msa,
                    query_seqs_unique,
                    query_seqs_cardinality,
                    template_features,
                ) = get_msa_and_templates(
                    jobname,
                    query_sequence,
                    result_dir,
                    msa_mode,
                    use_templates,
                    pair_mode,
                    host_url,
                )
            msa = msa_to_str(
                unpaired_msa, paired_msa, query_seqs_unique, query_seqs_cardinality
            )
            result_dir.joinpath(jobname + ".a3m").write_text(msa)
        except Exception as e:
            logger.exception(f"Could not get MSA/templates for {jobname}: {e}")
            continue
        try:
            input_features = generate_input_feature(
                query_seqs_unique,
                query_seqs_cardinality,
                unpaired_msa,
                paired_msa,
                template_features,
                is_complex,
                model_type,
            )
        except Exception as e:
            logger.exception(f"Could not generate input features {jobname}: {e}")
            continue
        try:
            query_sequence_len_array = [
                len(query_seqs_unique[i])
                for i, cardinality in enumerate(query_seqs_cardinality)
                for _ in range(0, cardinality)
            ]

            if sum(query_sequence_len_array) > crop_len:
                crop_len = math.ceil(sum(query_sequence_len_array) * recompile_padding)

            outs, model_rank = predict_structure(
                jobname,
                result_dir,
                input_features,
                is_complex,
                use_templates,
                sequences_lengths=query_sequence_len_array,
                crop_len=crop_len,
                model_type=model_type,
                model_runner_and_params=model_runner_and_params,
                do_relax=use_amber,
                rank_by=rank_by,
                stop_at_score=stop_at_score,
                prediction_callback=prediction_callback,
            )
        except RuntimeError as e:
            # This normally happens on OOM. TODO: Filter for the specific OOM error message
            logger.error(f"Could not predict {jobname}. Not Enough GPU memory? {e}")
            continue

        # Write alphafold-db format (PAE)
        alphafold_pae_file = result_dir.joinpath(
            jobname + "_predicted_aligned_error_v1.json"
        )
        alphafold_pae_file.write_text(get_pae_json(outs[0]["pae"], outs[0]["max_pae"]))
        num_alignment = (
            int(input_features["num_alignments"])
            if model_type == "AlphaFold2-multimer"
            else input_features["num_alignments"][0]
        )
        msa_plot = plot_msa(
            input_features["msa"][0:num_alignment],
            input_features["msa"][0],
            query_sequence_len_array,
            query_sequence_len,
        )
        coverage_png = result_dir.joinpath(jobname + "_coverage.png")
        msa_plot.savefig(str(coverage_png))
        msa_plot.close()
        paes_plot = plot_paes(
            [outs[k]["pae"] for k in model_rank], Ls=query_sequence_len_array, dpi=200
        )
        pae_png = result_dir.joinpath(jobname + "_PAE.png")
        paes_plot.savefig(str(pae_png))
        paes_plot.close()
        plddt_plot = plot_plddts(
            [outs[k]["plddt"] for k in model_rank], Ls=query_sequence_len_array, dpi=200
        )
        plddt_png = result_dir.joinpath(jobname + "_plddt.png")
        plddt_plot.savefig(str(plddt_png))
        plddt_plot.close()
        result_files = [
            bibtex_file,
            config_out_file,
            alphafold_pae_file,
            result_dir.joinpath(jobname + ".a3m"),
            pae_png,
            coverage_png,
            plddt_png,
        ]
        for i in model_rank:
            result_files.append(
                result_dir.joinpath(
                    f"{jobname}_unrelaxed_rank_{i + 1}_{outs[i]['model_name']}.pdb"
                )
            )
            result_files.append(
                result_dir.joinpath(
                    f"{jobname}_unrelaxed_rank_{i + 1}_{outs[i]['model_name']}_scores.json"
                )
            )
            if use_amber:
                result_files.append(
                    result_dir.joinpath(
                        f"{jobname}_relaxed_rank_{i + 1}_{outs[i]['model_name']}.pdb"
                    )
                )

        if zip_results:
            with zipfile.ZipFile(result_zip, "w") as result_zip:
                for file in result_files:
                    result_zip.write(file, arcname=file.name)
            # Delete only after the zip was successful, and also not the bibtex and config because we need those again
            for file in result_files[2:]:
                file.unlink()
        else:
            is_done_marker.touch()

    logger.info("Done")


def set_model_type(is_complex: bool, model_type: str) -> str:
    if model_type == "auto" and is_complex:
        model_type = "AlphaFold2-multimer"
    elif model_type == "auto" and not is_complex:
        model_type = "AlphaFold2-ptm"
    return model_type


def main():
    parser = ArgumentParser()
    parser.add_argument(
        "input",
        default="input",
        help="Can be one of the following: "
        "Directory with fasta/a3m files, a csv/tsv file, a fasta file or an a3m file",
    )
    parser.add_argument("results", help="Directory to write the results to")

    # Main performance parameter
    parser.add_argument(
        "--stop-at-score",
        help="Compute models until plddt or ptmscore > threshold is reached. "
        "This can make colabfold much faster by only running the first model for easy queries.",
        type=float,
        default=100,
    )

    parser.add_argument(
        "--num-recycle",
        help="Number of prediction cycles."
        "Increasing recycles can improve the quality but slows down the prediction.",
        type=int,
        default=3,
    )

    parser.add_argument("--num-models", type=int, default=5, choices=[1, 2, 3, 4, 5])
    parser.add_argument(
        "--recompile-padding",
        type=float,
        default=1.1,
        help="Whenever the input length changes, the model needs to be recompiled, which is slow. "
        "We pad sequences by this factor, so we can e.g. compute sequence from length 100 to 110 without recompiling. "
        "The prediction will become marginally slower for the longer input, "
        "but overall performance increases due to not recompiling. "
        "Set to 1 to disable.",
    )

    parser.add_argument("--model-order", default="3,4,5,1,2", type=str)
    parser.add_argument("--host-url", default=DEFAULT_API_SERVER)
    parser.add_argument("--data")

    # TODO: This currently isn't actually used
    parser.add_argument(
        "--msa-mode",
        default="MMseqs2 (UniRef+Environmental)",
        choices=[
            "MMseqs2 (UniRef+Environmental)",
            "MMseqs2 (UniRef only)",
            "single_sequence",
        ],
        help="Using an a3m file as input overwrites this option",
    )

    parser.add_argument(
        "--model-type",
        help="predict strucutre/complex using the following model."
        'Auto will pick "AlphaFold2" (ptm) for structure predictions and "AlphaFold2-multimer" for complexes.',
        type=str,
        default="auto",
        choices=["auto", "AlphaFold2-ptm", "AlphaFold2-multimer"],
    )

    parser.add_argument(
        "--amber",
        default=False,
        action="store_true",
        help="Use amber for structure refinement",
    )
    parser.add_argument(
        "--templates", default=False, action="store_true", help="Use templates from pdb"
    )
    parser.add_argument("--env", default=False, action="store_true")
    parser.add_argument(
        "--cpu",
        default=False,
        action="store_true",
        help="Allow running on the cpu, which is very slow",
    )
    parser.add_argument(
        "--rank",
        help="rank models by auto, plddt or ptmscore",
        type=str,
        default="auto",
        choices=["auto", "plddt", "ptmscore"],
    )
    parser.add_argument(
        "--pair-mode",
        help="rank models by auto, unpaired, paired, unpaired+paired",
        type=str,
        default="unpaired+paired",
        choices=["unpaired", "paired", "unpaired+paired"],
    )
    parser.add_argument(
        "--recompile-all-models",
        help="recompile all models instead of just model 1 ane 3",
        default=False,
        action="store_true",
    )
    parser.add_argument(
        "--sort-queries-by",
        help="sort queries by: none, length, random",
        type=str,
        default="length",
        choices=["none", "length", "random"],
    )
    parser.add_argument(
        "--zip",
        default=False,
        action="store_true",
        help="zip all results into one <jobname>.result.zip and delete the original files",
    )
    parser.add_argument(
        "--overwrite-existing-results", default=False, action="store_true"
    )

    args = parser.parse_args()

    setup_logging(Path(args.results).joinpath("log.txt"))

    data_dir = Path(args.data or default_data_dir)

    # Prevent people from accidentally running on the cpu, which is really slow
    if not args.cpu and xla_bridge.get_backend().platform == "cpu":
        print(NO_GPU_FOUND, file=sys.stderr)
        sys.exit(1)

    queries, is_complex = get_queries(args.input, args.sort_queries_by)
    model_type = set_model_type(is_complex, args.model_type)
    download_alphafold_params(model_type, data_dir)
    uses_api = any((query[2] is None for query in queries))
    if uses_api and args.host_url == DEFAULT_API_SERVER:
        print(ACCEPT_DEFAULT_TERMS, file=sys.stderr)

    model_order = [int(i) for i in args.model_order.split(",")]

    assert 1 <= args.recompile_padding, "Can't apply negative padding"

    run(
        queries=queries,
        result_dir=args.results,
        use_templates=args.templates,
        use_amber=args.amber,
        msa_mode=args.msa_mode,
        model_type=model_type,
        num_models=args.num_models,
        num_recycles=args.num_recycle,
        model_order=model_order,
        is_complex=is_complex,
        keep_existing_results=not args.overwrite_existing_results,
        rank_by=args.rank,
        pair_mode=args.pair_mode,
        data_dir=data_dir,
        host_url=args.host_url,
        stop_at_score=args.stop_at_score,
        recompile_padding=args.recompile_padding,
        recompile_all_models=args.recompile_all_models,
        zip_results=args.zip,
    )


if __name__ == "__main__":
    main()
