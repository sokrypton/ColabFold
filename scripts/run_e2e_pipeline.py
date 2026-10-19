import sys
import logging
import argparse

from pathlib import Path
from colabfold.download import download_alphafold_params
from colabfold.batch import run, get_queries
from colabfold.utils import setup_logging
from colabfold.caat_options import add_prediction_options, prediction_kwargs

from alphafold.analysis.utils import _parse_indices
from alphafold.analysis.attention_pipeline import run_pipeline
from alphafold.model.modules import reset_attention_state


logger = logging.getLogger(__name__)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="End-to-End AlphaAttention: Prediction + Attention Analysis",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    prediction = parser.add_argument_group("prediction settings")
    prediction.add_argument(
        "--query-seq-path",
        type=str,
        help="Query MSA/Fasta file.",
    )
    prediction.add_argument("--model-type", type=str, default="alphafold2")
    prediction.add_argument("--num-models", type=int, default=5)
    prediction.add_argument(
        "--random-seed",
        type=int,
        default=0,
        help="Starting random seed for prediction.",
    )
    prediction.add_argument(
        "--num-seeds",
        type=int,
        default=1,
        help="Number of consecutive seeds to run, starting at --random-seed.",
    )
    prediction.add_argument("--result-dir", type=str, default="results")
    prediction.add_argument(
        "--attention-output-dir",
        type=str,
        default="attention_outputs",
        help="Directory to save raw attention head NumPy files.\nDefault: attention_outputs",
    )
    prediction.add_argument(
        "--save-attention-compressed",
        action="store_true",
        help="If set, exports compressed attention weights in H5 format to local disk.",
        default=False,
    )
    prediction.add_argument(
        "--save-intermediate-structures",
        default=None,
        type=str,
        help="Directory to save intermediate structures from each evoformer loop.",
    )
    prediction.add_argument(
        "--distogram-output-dir",
        type=str,
        default=None,
        help="Directory to save distogram .npz files.",
    )

    analysis = parser.add_argument_group("analysis settings")
    analysis.add_argument(
        "--query-name",
        required=True,
        help="ID for query protein.",
    )
    analysis.add_argument(
        "--vis-output-dir",
        type=str,
        default="attention_visualizations",
    )
    analysis.add_argument(
        "--query-highlight-indices",
        default=None,
        help="Comma-separated 1-based indices to highlight in query (e.g. 1,5,10).",
    )
    analysis.add_argument(
        "--target-highlight-indices",
        default=None,
        help="Comma-separated 1-based indices to highlight in target (e.g. 1,5,10).",
    )
    analysis.add_argument(
        "--query-highlight-color",
        default="#AE0639",
        help="Hex color for query sequence highlights.",
    )
    analysis.add_argument(
        "--target-highlight-color",
        default="#1f77b4",
        help="Hex color for target sequence highlights.",
    )
    analysis.add_argument(
        "--save-attention-npy",
        action="store_true",
        help="Retain the raw .npy tensors selected for attention analysis.",
        default=False,
    )
    analysis.add_argument(
        "--generate-alignment",
        action="store_true",
        help="If set, generates a pairwise alignment file from the query and target sequences.",
    )
    analysis.add_argument(
        "--include-extra-msa",
        action="store_true",
        help="Emit and analyze extra-MSA Evoformer attention blocks.",
    )

    comparison = parser.add_argument_group("comparison settings (optional)")
    comparison.add_argument(
        "--target-name", default=None, help="ID for target protein."
    )
    comparison.add_argument(
        "--target-seq-path", default=None, help="Path to target sequence."
    )
    comparison.add_argument(
        "--alignment-path", default=None, help="Path to MSA alignment."
    )

    prediction_option_names = add_prediction_options(parser)
    args: argparse.Namespace = parser.parse_args()
    colabfold_options = prediction_kwargs(args, prediction_option_names)

    download_alphafold_params(args.model_type, Path("."))
    res_dir = Path(args.result_dir)
    res_dir.mkdir(parents=True, exist_ok=True)
    setup_logging(res_dir / "log.txt")

    logging.getLogger().setLevel(logging.INFO)

    base_attn_dir = Path(args.attention_output_dir)
    query_attn_dir = base_attn_dir / args.query_name
    target_attn_dir = base_attn_dir / args.target_name if args.target_name else None

    query_data, is_complex_query = get_queries(args.query_seq_path)
    logger.info("Generating Query Attention: %s", args.query_name)

    query_intermediate_path = (
        f"{args.save_intermediate_structures}/{args.query_name}"
        if args.save_intermediate_structures
        else None
    )
    logger.info(f"Intermediate structures path: {query_intermediate_path}")

    query_distogram_output_dir = (
        f"{args.distogram_output_dir}/{args.query_name}"
        if args.distogram_output_dir
        else None
    )
    logger.info(f"Distogram output path: {query_distogram_output_dir}")

    run(
        **colabfold_options,
        queries=query_data,
        result_dir=args.result_dir,
        num_models=args.num_models,
        random_seed=args.random_seed,
        num_seeds=args.num_seeds,
        attention_output_dir=str(query_attn_dir),
        model_type=args.model_type,
        is_complex=is_complex_query,
        keep_existing_results=False,
        save_attention_compressed=args.save_attention_compressed,
        include_extra_msa_attention=args.include_extra_msa,
        save_intermediate_structures=query_intermediate_path,
        distogram_output_dir=query_distogram_output_dir,
    )

    logging.getLogger().setLevel(logging.INFO)

    if args.target_seq_path and args.target_name:
        reset_attention_state()

        if args.query_name == args.target_name:
            logger.error("Query and Target names must be different.")
            sys.exit(1)

        target_data, is_complex_target = get_queries(args.target_seq_path)

        logger.info("Generating Target Attention: %s", args.target_name)

        target_intermediate_path = (
            f"{args.save_intermediate_structures}/{args.target_name}"
            if args.save_intermediate_structures
            else None
        )

        target_distogram_output_dir = (
            f"{args.distogram_output_dir}/{args.target_name}"
            if args.distogram_output_dir
            else None
        )

        run(
            **colabfold_options,
            queries=target_data,
            result_dir=args.result_dir,
            num_models=args.num_models,
            random_seed=args.random_seed,
            num_seeds=args.num_seeds,
            attention_output_dir=str(target_attn_dir),
            model_type=args.model_type,
            is_complex=is_complex_target,
            keep_existing_results=False,
            save_attention_compressed=args.save_attention_compressed,
            include_extra_msa_attention=args.include_extra_msa,
            save_intermediate_structures=target_intermediate_path,
            distogram_output_dir=target_distogram_output_dir,
        )
        logging.getLogger().setLevel(logging.INFO)

    logger.info("Starting Attention Analysis: %s", args.query_name)

    run_pipeline(
        query_seq_path=args.query_seq_path,
        query_attn_dir=str(query_attn_dir),
        query_name=args.query_name,
        save_path=args.vis_output_dir,
        target_seq_path=args.target_seq_path,
        target_attn_dir=str(target_attn_dir) if target_attn_dir else None,
        target_name=args.target_name,
        alignment_path=args.alignment_path,
        query_highlight_indices=_parse_indices(args.query_highlight_indices),
        target_highlight_indices=_parse_indices(args.target_highlight_indices),
        query_highlight_color=args.query_highlight_color,
        target_highlight_color=args.target_highlight_color,
        save_attention_npy=args.save_attention_npy,
        generate_alignment=args.generate_alignment,
        include_extra_msa=args.include_extra_msa,
    )

    logger.info("End-to-end pipeline complete. Results in %s", args.vis_output_dir)


if __name__ == "__main__":
    main()
