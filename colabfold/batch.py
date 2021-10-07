import argparse
import logging
import math
import pickle
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
    model_runner_and_params: Dict[str, Tuple[model.RunModel, haiku.Params]],
    do_relax: bool = False,
    rank_by: str = "auto",
    random_seed: int = 0,
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

    for model_name, (model_runner, params) in model_runner_and_params.items():
        logger.info(f"Running {model_name}")
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
            logger.info(
                f"{model_name} took {prediction_time:.1f}s with pLDDT {np.mean(prediction_result['plddt'][:seq_len]):.1f}"
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
            relaxed_pdb_lines.append(relaxed_pdb_str)

    # rerank models based on predicted lddt
    if rank_by == "ptmscore":
        model_rank = np.array(ptmscore).argsort()[::-1]
    else:
        model_rank = np.mean(plddts, -1).argsort()[::-1]
    out = {}
    logger.info("reranking models based on avg. predicted lDDT")
    for n, r in enumerate(model_rank):
        unrelaxed_pdb_path = result_dir.joinpath(
            f"{prefix}_unrelaxed_model_{n + 1}.pdb"
        )
        unrelaxed_pdb_path.write_text(unrelaxed_pdb_lines[r])
        set_bfactor(unrelaxed_pdb_path, plddts[r], idx_res, chains)

        if do_relax:
            relaxed_pdb_path = result_dir.joinpath(
                f"{prefix}_relaxed_model_{n + 1}.pdb"
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


def get_queries(input_path: Union[str, Path]) -> List[Tuple[str, str, Optional[str]]]:
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
            a3m_lines = input_path.read_text().upper()
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
        for file in input_path.iterdir():
            if not file.is_file():
                continue
            (seqs, header) = pipeline.parsers.parse_fasta(file.read_text())
            query_sequence = seqs[0]
            if len(seqs) > 1 and file.suffix == ".fasta":
                logger.warning(
                    f"More than one sequence in {file}, ignoring all but the first sequence"
                )

            if file.suffix.lower() == ".a3m":
                a3m_lines = file.read_text().upper()
            else:
                a3m_lines = None
            queries.append((file.stem, query_sequence.upper(), a3m_lines))

    # sort by seq. len
    queries.sort(key=lambda t: len(t[1]))
    return queries


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
    if use_templates:
        a3m_lines_mmseqs2, template_paths = run_mmseqs2(
            query_sequences,
            str(result_dir.joinpath(jobname)),
            use_env,
            use_templates=True,
            host_url=host_url,
        )
        if template_paths is None:
            template_features = mk_mock_template(query_sequences, 100)
        else:
            template_features = mk_template(
                a3m_lines_mmseqs2, template_paths, query_sequences
            )
        if not a3m_lines:
            a3m_lines = a3m_lines_mmseqs2
    else:
        if not a3m_lines:
            use_pairing = (
                not isinstance(query_sequences, str) and len(query_sequences) > 1
            )
            # find normal a3ms
            if (
                not use_pairing
                or pair_mode == "unpaired"
                or pair_mode == "unpaired+paired"
            ):
                a3m_lines = run_mmseqs2(
                    query_sequences,
                    str(result_dir.joinpath(jobname)),
                    use_env,
                    use_pairing=False,
                    host_url=host_url,
                )

            if use_pairing:
                if pair_mode == "paired" or pair_mode == "unpaired+paired":
                    # find paired a3m
                    paired_a3m_lines = run_mmseqs2(
                        query_sequences,
                        str(result_dir.joinpath(jobname)),
                        use_env,
                        use_pairing=True,
                        host_url=host_url,
                    )

                if pair_mode == "unpaired" or pair_mode == "unpaired+paired":
                    # pad sequences
                    _blank_seq = ["-" * len(seq) for seq in query_sequences]
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
                                        + [a3m_line]
                                        + _blank_seq[n + 1 :]
                                    )
                                )
                    if pair_mode == "unpaired":
                        a3m_lines = "\n".join(a3m_lines_combined)
                    else:
                        a3m_lines = (
                            paired_a3m_lines + "\n" + "\n".join(a3m_lines_combined)
                        )
                else:
                    a3m_lines = paired_a3m_lines

        template_features = mk_mock_template(query_sequences, 100)
    return a3m_lines, template_features


def run(
    queries: List[Tuple[str, Union[str, List[str]], Optional[str]]],
    result_dir: Union[str, Path],
    use_templates: bool,
    use_amber: bool,
    msa_mode: str,
    num_models: int,
    homooligomer: int,
    do_not_overwrite_results: bool,
    rank_mode: str,
    pair_mode: str,
    data_dir: Union[str, Path] = default_data_dir,
    host_url: str = DEFAULT_API_SERVER,
    cache: Optional[str] = None,
):
    result_dir = Path(result_dir)
    result_dir.mkdir(exist_ok=True)
    data_dir = Path(data_dir)

    use_env = msa_mode == "MMseqs2 (UniRef+Environmental)"

    # TODO: What's going on with MSA mode?
    write_bibtex(True, use_env, use_templates, use_amber, result_dir)

    model_runner_and_params = load_models_and_params(num_models, data_dir)

    crop_len = 0
    for job_number, (raw_jobname, query_sequence, a3m_lines) in enumerate(queries):
        jobname = safe_filename(raw_jobname)
        if (
            do_not_overwrite_results
            and result_dir.joinpath(jobname).with_suffix(".result.zip").is_file()
        ):
            logger.info(f"Skipping {jobname}")
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
            crop_len = math.ceil(query_sequence_len * 1.1)
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

        outs = predict_structure(
            jobname,
            result_dir,
            feature_dict,
            sequences_lengths=query_sequence_len_array,
            crop_len=crop_len,
            model_runner_and_params=model_runner_and_params,
            do_relax=use_amber,
            rank_by=rank_mode,
            cache=cache,
        )

        plot_lddt(homooligomer, jobname, msa, outs, query_sequence, result_dir)
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
    # Caches the model output. For development only
    parser.add_argument("--cache", help=argparse.SUPPRESS)
    parser.add_argument("--num-models", type=int, default=5, choices=[1, 2, 3, 4, 5])
    parser.add_argument(
        "--rank",
        help="rank models by auto, plddt or ptmscore",
        type=str,
        default="auto",
        choices=["auto", "plddt", "ptmscore"],
    )
    parser.add_argument("--homooligomer", type=int, default=1)
    parser.add_argument(
        "--pair-mode",
        help="rank models by auto, unpaired, paired, unpaired+paired",
        type=str,
        default="unpaired+paired",
        choices=["unpaired", "paired", "unpaired+paired"],
    )

    parser.add_argument("--data")
    parser.add_argument(
        "--do-not-overwrite-results", default=True, action="store_false"
    )

    parser.add_argument("--host-url", default=DEFAULT_API_SERVER)
    args = parser.parse_args()

    setup_logging(Path(args.results).joinpath("log.txt"))

    data_dir = Path(args.data or default_data_dir)
    download_alphafold_params(data_dir)

    assert args.msa_mode == "MMseqs2 (UniRef+Environmental)", "Unsupported"
    assert args.homooligomer == 1, "Unsupported"

    # Prevent people from accidentally running on the cpu, which is really slow
    if not args.cpu and xla_bridge.get_backend().platform == "cpu":
        print(NO_GPU_FOUND, file=sys.stderr)
        sys.exit(1)

    queries = get_queries(args.input)
    uses_api = any((query[2] is None for query in queries))
    if uses_api and args.host_url == DEFAULT_API_SERVER:
        print(ACCEPT_DEFAULT_TERMS, file=sys.stderr)

    run(
        queries=queries,
        result_dir=args.results,
        use_templates=args.templates,
        use_amber=args.amber,
        msa_mode=args.msa_mode,
        num_models=args.num_models,
        homooligomer=args.homooligomer,
        do_not_overwrite_results=args.do_not_overwrite_results,
        rank_mode=args.rank,
        pair_mode=args.pair_mode,
        data_dir=data_dir,
        host_url=args.host_url,
        cache=args.cache,
    )


if __name__ == "__main__":
    main()
