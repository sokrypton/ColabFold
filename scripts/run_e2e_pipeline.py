import sys
import logging
import argparse

from pathlib import Path

from colabfold.download import download_alphafold_params
from colabfold.batch import run, get_queries
from colabfold.utils import setup_logging

from alphafold.analysis.utils import _parse_indices
from alphafold.analysis.attention_pipeline import run_pipeline


logger = logging.getLogger(__name__)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="End-to-End AlphaAttention: Prediction + Attention Analysis",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    prediction = parser.add_argument_group("prediction settings")
    prediction.add_argument("input_fasta", type=str, help="Query MSA/Fasta file.")
    prediction.add_argument("--model-type", type=str, default="alphafold2_ptm")
    prediction.add_argument("--num-models", type=int, default=5)
    prediction.add_argument("--result-dir", type=str, default="results")

    analysis = parser.add_argument_group("analysis settings")
    analysis.add_argument("--query_name", required=True, help="ID for query protein.")
    analysis.add_argument(
        "--vis_output_dir", type=str, default="attention_visualizations"
    )
    analysis.add_argument(
        "--highlight_indices_query",
        default=None,
        help="Comma-separated 1-based indices to highlight in query (e.g. 1,5,10).",
    )
    analysis.add_argument(
        "--highlight_indices_target",
        default=None,
        help="Comma-separated 1-based indices to highlight in target (e.g. 1,5,10).",
    )
    analysis.add_argument(
        "--highlight_color_query",
        default="#AE0639",
        help="Hex color for query sequence highlights.",
    )
    analysis.add_argument(
        "--highlight_color_target",
        default="#1f77b4",
        help="Hex color for target sequence highlights.",
    )

    comparison = parser.add_argument_group("comparison settings (optional)")
    comparison.add_argument(
        "--target_name", default=None, help="ID for target protein."
    )
    comparison.add_argument(
        "--target_seq_path", default=None, help="Path to target sequence."
    )
    comparison.add_argument(
        "--alignment_path", default=None, help="Path to MSA alignment."
    )

    args: argparse.Namespace = parser.parse_args()

    base_attn_dir = Path("attention_outputs")
    query_attn_dir = base_attn_dir / args.query_name
    target_attn_dir = base_attn_dir / args.target_name if args.target_name else None

    download_alphafold_params(args.model_type, Path("."))
    res_dir = Path(args.result_dir)
    res_dir.mkdir(parents=True, exist_ok=True)
    setup_logging(res_dir / "log.txt")

    logging.getLogger().setLevel(logging.INFO)

    query_data, is_complex_query = get_queries(args.input_fasta)
    logger.info("Generating Query Attention: %s", args.query_name)

    run(
        queries=query_data,
        result_dir=args.result_dir,
        num_models=args.num_models,
        attention_output_dir=str(query_attn_dir),
        model_type=args.model_type,
        is_complex=is_complex_query,
    )

    logging.getLogger().setLevel(logging.INFO)

    if args.target_seq_path and args.target_name:
        if args.query_name == args.target_name:
            logger.error("Query and Target names must be different.")
            sys.exit(1)

        target_data, is_complex_target = get_queries(args.target_seq_path)

        logger.info("Generating Target Attention: %s", args.target_name)
        run(
            queries=target_data,
            result_dir=args.result_dir,
            num_models=args.num_models,
            attention_output_dir=str(target_attn_dir),
            model_type=args.model_type,
            is_complex=is_complex_target,
        )
        logging.getLogger().setLevel(logging.INFO)

    logger.info("Starting Attention Analysis: %s", args.query_name)
    run_pipeline(
        query_seq_path=args.input_fasta,
        query_attn_dir=str(query_attn_dir),
        query_name=args.query_name,
        save_path=args.vis_output_dir,
        target_seq_path=args.target_seq_path,
        target_attn_dir=str(target_attn_dir) if target_attn_dir else None,
        target_name=args.target_name,
        alignment_path=args.alignment_path,
        highlight_indices_query=_parse_indices(args.highlight_indices_query),
        highlight_indices_target=_parse_indices(args.highlight_indices_target),
        highlight_color_query=args.highlight_color_query,
        highlight_color_target=args.highlight_color_target,
    )

    logger.info("End-to-end pipeline complete. Results in %s", args.vis_output_dir)


if __name__ == "__main__":
    main()
