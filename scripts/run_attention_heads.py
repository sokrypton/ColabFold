import os
import argparse

from pathlib import Path
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

    # --- Required Arguments ---
    parser.add_argument(
        "input_fasta", type=str, help="Path to the input MSA file (.a3m or .fasta)."
    )

    # --- Optional Arguments ---
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

    args = parser.parse_args()

    print(f"Starting ColabFold Prediction for: {args.input_fasta}")
    print(f"Model Type: {args.model_type}")
    print(f"Attention Output: {args.attention_output_dir}")
    print("-" * 30)

    if not os.path.exists(args.input_fasta):
        raise FileNotFoundError(f"Input file not found at: {args.input_fasta}")

    queries, is_complex = get_queries(args.input_fasta)

    download_alphafold_params(args.model_type, Path("."))

    main_results_dir = Path(args.result_dir)

    log_file_path = main_results_dir.joinpath("log.txt")
    main_results_dir.mkdir(parents=True, exist_ok=True)

    setup_logging(log_file_path, verbose=False)

    results = run(
        queries=queries,
        result_dir=args.result_dir,
        num_models=args.num_models,
        attention_output_dir=args.attention_output_dir,
        model_type=args.model_type,
        is_complex=is_complex,
    )

    print("-" * 30)
    print("Prediction complete!")
    print(f"Structures saved to: {args.result_dir}")
    print(f"Attention heads saved to: {args.attention_output_dir}")


if __name__ == "__main__":
    main()
