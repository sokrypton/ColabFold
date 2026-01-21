import sys
import logging
import argparse

from pathlib import Path

from colabfold.download import download_alphafold_params
from colabfold.batch import run, get_queries
from colabfold.utils import setup_logging

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
    prediction.add_argument("--model-type", type=str, default="alphafold2_ptm")
    prediction.add_argument("--num-models", type=int, default=5)
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
    )
    prediction.add_argument(
        "--save-intermediate-structures",
        default=None,
        type=str,
        help="Directory to save intermediate structures from each evoformer loop.",
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
        help="If set, exports individual uncompressed attention heads (.npy) to local disk.",
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

    args: argparse.Namespace = parser.parse_args()

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

    run(
        queries=query_data,
        result_dir=args.result_dir,
        num_models=args.num_models,
        attention_output_dir=str(query_attn_dir),
        model_type=args.model_type,
        is_complex=is_complex_query,
        save_attention_compressed=args.save_attention_compressed,
        save_intermediate_structures=args.save_intermediate_structures,
    )

    logging.getLogger().setLevel(logging.INFO)

    if args.target_seq_path and args.target_name:
        reset_attention_state()

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
            save_attention_compressed=args.save_attention_compressed,
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
    )

    logger.info("End-to-end pipeline complete. Results in %s", args.vis_output_dir)


if __name__ == "__main__":
    main()
