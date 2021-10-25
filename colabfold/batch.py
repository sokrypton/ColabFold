import argparse
import json
import logging
import math
import pickle
import random
import sys
import time
from argparse import ArgumentParser
from pathlib import Path
from string import ascii_uppercase
from typing import Any, Dict, Tuple, List, Union, Mapping, Optional

import haiku
import numpy
import numpy as np
import pandas
from alphafold.common import protein
from alphafold.common.protein import Protein
from alphafold.data import pipeline, templates
from alphafold.data.tools import hhsearch
from alphafold.model import model
from alphafold.model.features import FeatureDict
from jax.lib import xla_bridge

from colabfold.citations import write_bibtex
from colabfold.colabfold import run_mmseqs2
from colabfold.download import download_alphafold_params, default_data_dir
from colabfold.models import load_models_and_params
from colabfold.msa import make_fixed_size
from colabfold.pdb import set_bfactor
from colabfold.plot import plot_predicted_alignment_error, plot_lddt
from colabfold.utils import (
    setup_logging,
    safe_filename,
    NO_GPU_FOUND,
    DEFAULT_API_SERVER,
    ACCEPT_DEFAULT_TERMS,
)

logger = logging.getLogger(__name__)


def mk_mock_template(query_sequence, num_temp: int = 1) -> Mapping[str, Any]:
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
    model_runner_and_params: List[Tuple[str, model.RunModel, haiku.Params]],
    do_relax: bool = False,
    rank_by: str = "auto",
    random_seed: int = 0,
    stop_at_score: float = 100,
    cache: Optional[str] = None,
):
    """Predicts structure using AlphaFold for the given sequence."""
    # Run the models.
    if rank_by == "auto":
        # score complexes by ptmscore and sequences by plddt
        rank_by = "plddt" if len(sequences_lengths) == 1 else "ptmscore"

    plddts, paes, ptmscore = [], [], []
    unrelaxed_pdb_lines = []
    relaxed_pdb_lines = []
    prediction_times = []
    seq_len = feature_dict["seq_length"][0]

    # Minkyung's code
    # add big enough number to residue index to indicate chain breaks
    idx_res = feature_dict["residue_index"]
    L_prev = 0
    # Ls: number of residues in each chain
    for L_i in sequences_lengths[:-1]:
        idx_res[L_prev + L_i :] += 200
        L_prev += L_i

    chains = list(
        "".join([ascii_uppercase[n] * L for n, L in enumerate(sequences_lengths)])
    )
    feature_dict["residue_index"] = idx_res
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
                input_fix, model_name, model_runner, prefix, cache
            )
        else:
            start = time.time()
            prediction_result = model_runner.predict(input_fix)
            prediction_time = time.time() - start
            prediction_times.append(prediction_time)
            mean_plddt = np.mean(prediction_result["plddt"][:seq_len])
            logger.info(
                f"{model_name} took {prediction_time:.1f}s with pLDDT {mean_plddt :.1f}"
            )

        unrelaxed_protein = protein.from_prediction(input_fix, prediction_result)
        unrelaxed_pdb_lines.append(protein.to_pdb(unrelaxed_protein))
        plddts.append(prediction_result["plddt"][:seq_len])
        ptmscore.append(prediction_result["ptm"])
        paes_res = []

        for i in range(seq_len):
            paes_res.append(prediction_result["predicted_aligned_error"][i][:seq_len])
        paes.append(paes_res)
        if do_relax:
            from alphafold.relax import relax
            from alphafold.common import residue_constants

            # Hack so that we don't need to download into the alphafold package itself
            residue_constants.stereo_chemical_props_path = "stereo_chemical_props.txt"

            # Remove the padding because unlike to_pdb() amber doesn't handle that
            remove_padding_mask = unrelaxed_protein.atom_mask.sum(axis=-1) > 0
            unrelaxed_protein = Protein(
                atom_mask=unrelaxed_protein.atom_mask[remove_padding_mask],
                atom_positions=unrelaxed_protein.atom_positions[remove_padding_mask],
                aatype=unrelaxed_protein.aatype[remove_padding_mask],
                residue_index=unrelaxed_protein.residue_index[remove_padding_mask],
                b_factors=unrelaxed_protein.b_factors[remove_padding_mask],
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
        if np.mean(prediction_result["plddt"][:seq_len]) > stop_at_score:
            break
    # rerank models based on predicted lddt
    if rank_by == "ptmscore":
        model_rank = np.array(ptmscore).argsort()[::-1]
    else:
        model_rank = np.mean(plddts, -1).argsort()[::-1]
    out = {}
    logger.info("reranking models based on avg. predicted lDDT")
    for n, r in enumerate(model_rank):
        unrelaxed_pdb_path = result_dir.joinpath(
            f"{prefix}_unrelaxed_{model_names[r]}_rank_{n + 1}.pdb"
        )
        unrelaxed_pdb_path.write_text(unrelaxed_pdb_lines[r])
        set_bfactor(unrelaxed_pdb_path, plddts[r], idx_res, chains)

        if do_relax:
            relaxed_pdb_path = result_dir.joinpath(
                f"{prefix}_relaxed_{model_names[r]}_rank_{n + 1}.pdb"
            )
            relaxed_pdb_path.write_text(unrelaxed_pdb_lines[r])
            set_bfactor(relaxed_pdb_path, plddts[r], idx_res, chains)

        out[f"model_{n + 1}"] = {
            "plddt": plddts[r],
            "pae": paes[r],
            "pTMscore": ptmscore,
        }
    return out


def run_model_cached(
    input_fix: FeatureDict,
    model_name: str,
    model_runner: model.RunModel,
    prefix: str,
    cache: str,
):
    """Caching the expensive compilation + prediction step - for development only

    We store both input and output to ensure that the input is actually the same that we cached
    """
    pickle_path = Path(cache).joinpath(prefix)
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


def get_queries(
    input_path: Union[str, Path], sort_queries_by: str = "length"
) -> List[Tuple[str, str, Optional[str]]]:
    """Reads a directory of fasta files, a single fasta file or a csv file and returns a tuple
    of job name, sequence and the optional a3m lines"""

    input_path = Path(input_path)
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
            (seqs, header) = pipeline.parsers.parse_fasta(input_path.read_text())
            query_sequence = seqs[0]
            a3m_lines = input_path.read_text()
            queries = [(input_path.stem, query_sequence, a3m_lines)]
        elif input_path.suffix == ".fasta":
            (sequences, headers) = pipeline.parsers.parse_fasta(input_path.read_text())
            queries = [
                (header, sequence.upper(), None)
                for sequence, header in zip(sequences, headers)
            ]
        else:
            raise ValueError(f"Unknown file format {input_path.suffix}")
    else:
        assert input_path.is_dir(), "Expected either an input file or a input directory"
        queries = []
        for file in sorted(input_path.iterdir()):
            if not file.is_file():
                continue
            (seqs, header) = pipeline.parsers.parse_fasta(file.read_text())
            query_sequence = seqs[0]
            if len(seqs) > 1 and file.suffix == ".fasta":
                logger.warning(
                    f"More than one sequence in {file}, ignoring all but the first sequence"
                )

            if file.suffix.lower() == ".a3m":
                a3m_lines = file.read_text()
            else:
                a3m_lines = None
            queries.append((file.stem, query_sequence.upper(), a3m_lines))

    # sort by seq. len
    if sort_queries_by == "length":
        queries.sort(key=lambda t: len(t[1]))
    elif sort_queries_by == "random":
        random.shuffle(queries)
    return queries


def pair_sequences(
    a3m_lines: List[str], query_sequences: List[str], query_cardinality: List[int]
) -> str:
    a3m_line_paired = [""] * len(a3m_lines[0].splitlines())
    for n, seq in enumerate(query_sequences):
        lines = a3m_lines[n].splitlines()
        for i, line in enumerate(lines):
            if line.startswith(">"):
                if n != 0:
                    line = line.replace(">", "_", 1)
                a3m_line_paired[i] = a3m_line_paired[i] + line
            else:
                a3m_line_paired[i] = a3m_line_paired[i] + line * query_cardinality[n]
    return "\n".join(a3m_line_paired)


def pad_sequences(
    a3m_lines: List[str], query_sequences: List[str], query_cardinality: List[int]
) -> str:
    _blank_seq = [
        ("-" * len(seq)) * query_cardinality[n] for n, seq in enumerate(query_sequences)
    ]
    a3m_lines_combined = []
    for n, seq in enumerate(query_sequences):
        lines = a3m_lines[n].split("\n")
        for a3m_line in lines:
            if len(a3m_line) == 0:
                continue
            if a3m_line.startswith(">"):
                a3m_lines_combined.append(a3m_line)
            else:
                a3m_lines_combined.append(
                    "".join(
                        _blank_seq[:n]
                        + [a3m_line] * query_cardinality[n]
                        + _blank_seq[n + 1 :]
                    )
                )
    return "\n".join(a3m_lines_combined)


def get_msa_and_templates(
    a3m_lines: Optional[str],
    jobname: str,
    query_sequences,
    result_dir: Path,
    use_env: bool,
    use_templates: bool,
    pair_mode: str,
    host_url: str = DEFAULT_API_SERVER,
) -> Tuple[str, Mapping[str, Any]]:
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

    if use_templates:
        a3m_lines_mmseqs2, template_paths = run_mmseqs2(
            query_seqs_unique,
            str(result_dir.joinpath(jobname)),
            use_env,
            use_templates=True,
            host_url=host_url,
        )
        if template_paths is None:
            template_features = mk_mock_template(query_sequences, 100)
        else:
            if len(a3m_lines_mmseqs2) != 1 or len(query_seqs_unique) != 1:
                raise ValueError(
                    "Templates are currently only supported for one sequences"
                )
            template_features = mk_template(
                a3m_lines_mmseqs2[0], template_paths, query_seqs_unique[0]
            )
        if not a3m_lines:
            a3m_lines = a3m_lines_mmseqs2
    else:
        template_features = mk_mock_template(query_sequences, 100)

    if not a3m_lines:
        if len(query_sequences) == 1:
            pair_mode = "none"
        if (
            pair_mode == "none"
            or pair_mode == "unpaired"
            or pair_mode == "unpaired+paired"
        ):
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
            # find paired a3m
            paired_a3m_lines = run_mmseqs2(
                query_seqs_unique,
                str(result_dir.joinpath(jobname)),
                use_env,
                use_pairing=True,
                host_url=host_url,
            )
        else:
            paired_a3m_lines = None

        if pair_mode == "none":
            assert isinstance(a3m_lines, list) and len(a3m_lines) == 1
            [a3m_lines] = a3m_lines
        elif pair_mode == "unpaired":
            a3m_lines = pad_sequences(
                a3m_lines, query_seqs_unique, query_seqs_cardinality
            )
        elif pair_mode == "unpaired+paired":
            a3m_lines = (
                pair_sequences(
                    paired_a3m_lines, query_seqs_unique, query_seqs_cardinality
                )
                + "\n"
                + pad_sequences(a3m_lines, query_seqs_unique, query_seqs_cardinality)
            )
        elif pair_mode == "paired":
            a3m_lines = pair_sequences(
                paired_a3m_lines, query_seqs_unique, query_seqs_cardinality
            )
        else:
            raise ValueError(f"Invalid pair_mod: {pair_mode}")

    return a3m_lines, template_features


def run(
    queries: List[Tuple[str, Union[str, List[str]], Optional[str]]],
    result_dir: Union[str, Path],
    use_templates: bool,
    use_amber: bool,
    msa_mode: str,
    num_models: int,
    model_order: List[int],
    keep_existing_results: bool,
    rank_mode: str,
    pair_mode: str,
    data_dir: Union[str, Path] = default_data_dir,
    host_url: str = DEFAULT_API_SERVER,
    cache: Optional[str] = None,
    stop_at_score: float = 100,
    recompile_padding: float = 1.1,
    recompile_all_models: bool = False,
):
    data_dir = Path(data_dir)
    result_dir = Path(result_dir)
    result_dir.mkdir(exist_ok=True)

    # Record the parameters of this run
    result_dir.joinpath("config.json").write_text(
        json.dumps(
            {
                "num_queries": len(queries),
                "use_templates": use_templates,
                "use_amber": use_amber,
                "msa_mode": msa_mode,
                "num_models": num_models,
                "model_order": model_order,
                "keep_existing_results": keep_existing_results,
                "rank_mode": rank_mode,
                "pair_mode": pair_mode,
                "host_url": host_url,
                "stop_at_score": stop_at_score,
                "recompile_padding": recompile_padding,
                "recompile_all_models": recompile_all_models,
            }
        )
    )

    use_env = msa_mode == "MMseqs2 (UniRef+Environmental)"

    # TODO: What's going on with MSA mode?
    write_bibtex(True, use_env, use_templates, use_amber, result_dir)

    model_runner_and_params = load_models_and_params(
        num_models, model_order, data_dir, recompile_all_models
    )

    crop_len = 0
    for job_number, (raw_jobname, query_sequence, a3m_lines) in enumerate(queries):
        jobname = safe_filename(raw_jobname)
        # In the colab version we know we're done when a zip file has been written
        if (
            keep_existing_results
            and result_dir.joinpath(jobname).with_suffix(".result.zip").is_file()
        ):
            logger.info(f"Skipping {jobname} (result.zip)")
            continue
        # In the local version we don't zip the files, so assume we're done if the last unrelaxed pdb file exists
        last_pdb_file = f"{jobname}_unrelaxed_model_{num_models}.pdb"
        if keep_existing_results and result_dir.joinpath(last_pdb_file).is_file():
            logger.info(f"Skipping {jobname} (pdb)")
            continue
        query_sequence_len_array = (
            [len(query_sequence)]
            if isinstance(query_sequence, str)
            else [len(q) for q in query_sequence]
        )

        logger.info(
            f"Query {job_number + 1}/{len(queries)}: {jobname} (length {sum(query_sequence_len_array)})"
        )

        a3m_file = f"{jobname}.a3m"
        query_sequence_len = (
            len(query_sequence)
            if isinstance(query_sequence, str)
            else sum(len(s) for s in query_sequence)
        )

        if query_sequence_len > crop_len:
            crop_len = math.ceil(query_sequence_len * recompile_padding)
        try:
            a3m_lines, template_features = get_msa_and_templates(
                a3m_lines,
                jobname,
                query_sequence,
                result_dir,
                use_env,
                use_templates,
                pair_mode,
                host_url,
            )
        except Exception as e:
            logger.exception(f"Could not get MSA/templates for {jobname}: {e}")
            continue

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
                    sequence=query_sequence
                    if isinstance(query_sequence, str)
                    else "".join(query_sequence),
                    description="none",
                    num_res=query_sequence_len,
                ),
                **pipeline.make_msa_features(
                    msas=msas, deletion_matrices=deletion_matrices
                ),
                **template_features,
            }
        except Exception as e:
            logger.exception(f"Could not predict {jobname}: {e}")
            continue

        try:
            outs = predict_structure(
                jobname,
                result_dir,
                feature_dict,
                sequences_lengths=query_sequence_len_array,
                crop_len=crop_len,
                model_runner_and_params=model_runner_and_params,
                do_relax=use_amber,
                rank_by=rank_mode,
                stop_at_score=stop_at_score,
                cache=cache,
            )
        except RuntimeError as e:
            # This normally happens on OOM. TODO: Filter for the specific OOM error message
            logger.error(f"Could not predict {jobname}. Not Enough GPU memory? {e}")
            continue

        plot_lddt(jobname, msa, outs, query_sequence, result_dir)
        plot_predicted_alignment_error(jobname, num_models, outs, result_dir)
    logger.info("Done")


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
            "custom",
        ],
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
        "--overwrite-existing-results", default=False, action="store_true"
    )

    # Caches the model output. For development only. Remove eventually
    parser.add_argument("--cache", help=argparse.SUPPRESS)

    args = parser.parse_args()

    setup_logging(Path(args.results).joinpath("log.txt"))

    data_dir = Path(args.data or default_data_dir)
    download_alphafold_params(data_dir)

    assert args.msa_mode == "MMseqs2 (UniRef+Environmental)", "Unsupported"

    # Prevent people from accidentally running on the cpu, which is really slow
    if not args.cpu and xla_bridge.get_backend().platform == "cpu":
        print(NO_GPU_FOUND, file=sys.stderr)
        sys.exit(1)

    queries = get_queries(args.input, args.sort_queries_by)
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
        num_models=args.num_models,
        model_order=model_order,
        keep_existing_results=not args.overwrite_existing_results,
        rank_mode=args.rank,
        pair_mode=args.pair_mode,
        data_dir=data_dir,
        host_url=args.host_url,
        cache=args.cache,
        stop_at_score=args.stop_at_score,
        recompile_padding=args.recompile_padding,
        recompile_all_models=args.recompile_all_models,
    )


if __name__ == "__main__":
    main()
