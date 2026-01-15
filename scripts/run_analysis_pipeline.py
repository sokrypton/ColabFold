import argparse

from alphafold.analysis.utils import _parse_indices
from alphafold.analysis.attention_pipeline import run_pipeline


def main():
    parser = argparse.ArgumentParser(
        description="Visualize and compare AlphaFold attention heads.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument(
        "--query-attn-dir",
        required=True,
        help="Path to the directory containing .npy attention files for the query protein.",
    )
    parser.add_argument(
        "--query-name",
        required=True,
        help="Identifier for the query protein (e.g., 'FUBP1').",
    )
    parser.add_argument(
        "--query-seq-path",
        required=True,
        help="Path to the .a3m or .fasta file for the query protein.",
    )
    parser.add_argument(
        "--output-dir",
        default="attention_visualizations",
        help="Directory where output visualizations will be saved.",
    )
    parser.add_argument(
        "--target-attn-dir",
        default=None,
        help="Path to attention directory for the target protein.",
    )
    parser.add_argument(
        "--target-name", default=None, help="Identifier for the target protein."
    )
    parser.add_argument(
        "--target-seq-path",
        default=None,
        help="Path to sequence file for the target protein.",
    )
    parser.add_argument(
        "--alignment-path",
        default=None,
        help="Path to the alignment file (e.g., .ali) mapping query to target.",
    )
    parser.add_argument(
        "--query-highlight-indices",
        default=None,
        help="Comma-separated 1-based indices to highlight in query (e.g. 1,5,10).",
    )
    parser.add_argument(
        "--target-highlight-indices",
        default=None,
        help="Comma-separated 1-based indices to highlight in target (e.g. 1,5,10).",
    )
    parser.add_argument(
        "--query-highlight-color",
        default="#AE0639",
        help="Hex color for query sequence highlights.",
    )
    parser.add_argument(
        "--target-highlight-color",
        default="#1f77b4",
        help="Hex color for target sequence highlights.",
    )
    parser.add_argument(
        "--save-attention-h5",
        action="store_true",
        help="If set, exports attention weights in H5 format to local disk.",
    )

    args = parser.parse_args()

    if args.alignment_path:
        missing = [
            arg
            for arg, val in {
                "--target_attn_dir": args.target_attn_dir,
                "--target_name": args.target_name,
                "--target_seq_path": args.target_seq_path,
            }.items()
            if val is None
        ]
        if missing:
            parser.error(
                f"--alignment_path requires the following: {', '.join(missing)}"
            )

    run_pipeline(
        query_seq_path=args.query_seq_path,
        target_seq_path=args.target_seq_path,
        query_attn_dir=args.query_attn_dir,
        target_attn_dir=args.target_attn_dir,
        query_name=args.query_name,
        target_name=args.target_name,
        alignment_path=args.alignment_path,
        save_path=args.output_dir,
        query_highlight_indices=_parse_indices(args.query_highlight_indices),
        target_highlight_indices=_parse_indices(args.target_highlight_indices),
        query_highlight_color=args.query_highlight_color,
        target_highlight_color=args.target_highlight_color,
        save_attention_h5=args.save_attention_h5,
    )


if __name__ == "__main__":
    main()
