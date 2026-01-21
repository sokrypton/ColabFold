import os
import argparse

from pathlib import Path
from alphafold.model.modules import reset_attention_state
from colabfold.download import download_alphafold_params
from colabfold.batch import run, get_queries
from colabfold.utils import setup_logging


def main():
    """
    Parses command-line arguments and executes the ColabFold prediction
    with custom attention head collection.
    """
    parser = argparse.ArgumentParser(
        description="Run custom ColabFold prediction and collect attention heads.",
        formatter_class=argparse.RawTextHelpFormatter,
    )

    parser.add_argument(
        "--query-seq-path",
        type=str,
        help="Path to the input MSA file (.a3m or .fasta).",
    )
    parser.add_argument(
        "--model-type",
        type=str,
        default="alphafold2_ptm",
        help="AlphaFold model type (e.g., alphafold2_ptm, alphafold2_multimer_v3).\nDefault: alphafold2_ptm",
    )
    parser.add_argument(
        "--attention-output-dir",
        type=str,
        default="attention_outputs",
        help="Directory to save raw attention head NumPy files.\nDefault: attention_outputs",
    )
    parser.add_argument(
        "--query-name",
        required=True,
        help="ID for query protein.",
    )
    parser.add_argument(
        "--result-dir",
        type=str,
        default="results",
        help="Directory to save final PDB structures and metadata.\nDefault: results",
    )
    parser.add_argument(
        "--num-models",
        type=int,
        default=5,
        help="Number of models to run for each query.\nDefault: 5",
    )
    parser.add_argument(
        "--save-attention-compressed",
        action="store_true",
        help="If set, exports compressed attention weights in H5 format to local disk.",
    )
    parser.add_argument(
        "--save-intermediate-structures",
        default=None,
        type=str,
        help="Directory to save intermediate structures from each evoformer loop.",
    )

    args = parser.parse_args()

    print(f"Starting ColabFold Prediction for: {args.query_seq_path}")
    print(f"Model Type: {args.model_type}")
    print(f"Attention Output: {args.attention_output_dir}")
    print("-" * 30)

    if not os.path.exists(args.query_seq_path):
        raise FileNotFoundError(f"Input file not found at: {args.query_seq_path}")

    queries, is_complex = get_queries(args.query_seq_path)

    download_alphafold_params(args.model_type, Path("."))

    main_results_dir = Path(args.result_dir)

    log_file_path = main_results_dir.joinpath("log.txt")
    main_results_dir.mkdir(parents=True, exist_ok=True)

    setup_logging(log_file_path, verbose=False)

    base_attn_dir = Path(args.attention_output_dir)
    query_attn_dir = base_attn_dir / args.query_name

    results = run(
        queries=queries,
        result_dir=args.result_dir,
        num_models=args.num_models,
        attention_output_dir=str(query_attn_dir),
        model_type=args.model_type,
        is_complex=is_complex,
        save_attention_compressed=args.save_attention_compressed,
        save_intermediate_structures=args.save_intermediate_structures,
    )

    reset_attention_state()

    print("-" * 30)
    print("Prediction complete!")
    print(f"Structures saved to: {args.result_dir}")
    print(f"Attention heads saved to: {args.attention_output_dir}")


if __name__ == "__main__":
    main()
