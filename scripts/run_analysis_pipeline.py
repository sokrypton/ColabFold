import argparse

from alphafold.analysis.attention_pipeline import run_pipeline


def main():
    parser = argparse.ArgumentParser(
        description="Visualize and compare AlphaFold attention heads.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument(
        "--query_attn_dir",
        required=True,
        help="Path to the directory containing .npy attention files for the query protein.",
    )
    parser.add_argument(
        "--query_name",
        required=True,
        help="Identifier for the query protein (e.g., 'FUBP1').",
    )
    parser.add_argument(
        "--query_seq_path",
        required=True,
        help="Path to the .a3m or .fasta file for the query protein.",
    )
    parser.add_argument(
        "--output_dir",
        default="attention_visualizations",
        help="Directory where output visualizations will be saved.",
    )

    # Target Protein (Optional for Comparison)
    parser.add_argument(
        "--target_attn_dir",
        default=None,
        help="Path to attention directory for the target protein.",
    )
    parser.add_argument(
        "--target_name", default=None, help="Identifier for the target protein."
    )
    parser.add_argument(
        "--target_seq_path",
        default=None,
        help="Path to sequence file for the target protein.",
    )
    parser.add_argument(
        "--alignment_path",
        default=None,
        help="Path to the alignment file (e.g., .ali) mapping query to target.",
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
    )


if __name__ == "__main__":
    main()
