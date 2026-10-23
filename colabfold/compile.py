"""colabfold_compile: AlphaFold2 compiled for every length bucket into the compile cache, once per machine.

Takes colabfold_batch's flags, so the compiled models are those its folds load. One synthetic query per bucket, its MSA
only itself: with the fast kernels every MSA runs at the model's full depth, so it has the shape of any real query of
that length. Workers are this tool again on a share of the buckets; XLA compiles mostly on one thread, so they run side
by side.
"""
import os
import subprocess
import sys
import tempfile
from argparse import SUPPRESS
from pathlib import Path

from colabfold.alphafold.backend import BUCKETS
from colabfold.backend import is_af3_model
from colabfold.batch import build_parser, compile_cache_dir, parse_buckets, run_args
from colabfold.download import default_data_dir, download_alphafold_params


def build_compile_parser():
    parser = build_parser(positional=False)
    parser.description = "Compile AlphaFold2 for every length bucket into the compile cache that colabfold_batch loads."
    group = parser.add_argument_group("Compile arguments", "")
    group.add_argument("--max-length", default=2048, type=int, help="Longest bucket compiled.")
    group.add_argument("--workers", default=0, type=int,
                       help="Processes compiling side by side (0: up to 4, by CPU count).")
    group.add_argument("--clean", action="store_true",
                       help="Remove the compiled models other installations (jax, driver, GPU or ColabFold versions) "
                       "left in the compile cache, and exit.")
    group.add_argument("--clean-all", action="store_true",
                       help="Remove all your compiled models from the compile cache, this installation's too, and exit.")
    group.add_argument("--worker", nargs=2, metavar=("INPUTS", "RESULTS"), help=SUPPRESS)
    return parser


def main(argv=None):
    argv = sys.argv[1:] if argv is None else argv
    parser = build_compile_parser()
    args = parser.parse_args(argv)
    if args.worker:
        args.input, args.results = args.worker
        run_args(args, compile_only=True)
        return

    # only the fast kernels run at bucketed lengths: without them nothing here would be reused
    args.use_fast_kernels = True
    model_type = "alphafold2_ptm" if args.model_type == "auto" else args.model_type
    if is_af3_model(model_type):
        parser.error("compiles the alphafold2 models")
    if "multimer" in model_type:
        parser.error("compiles the monomer models: a multimer compiles at its MSA's depth as well as its length")
    cache = compile_cache_dir(args)
    if cache is None:
        parser.error("needs a --compile-cache to compile into")
    if args.clean or args.clean_all:
        from colabfold.alphafold.compile_cache import clean, install_hash
        removed = clean(cache, keep=None if args.clean_all else install_hash())
        print(f"Removed {removed} compiled models from {cache}")
        return
    buckets = parse_buckets(args.buckets)
    lengths = sorted((b for b in buckets or BUCKETS if b <= args.max_length), reverse=True)
    if not lengths:
        parser.error(f"no bucket up to --max-length {args.max_length}")

    download_alphafold_params(model_type, Path(args.data or default_data_dir))  # once, not by every worker
    workers = min(args.workers or max(1, min(4, (os.cpu_count() or 4) // 4)), len(lengths))
    env = {**os.environ, "XLA_PYTHON_CLIENT_PREALLOCATE": "false"}  # only compiling: they share the GPU
    flags = [a for a in argv if a != "--templates"] + ["--use-fast-kernels"]  # template folds compile on first use
    with tempfile.TemporaryDirectory(prefix="colabfold_compile_") as tmp:
        procs = []
        for w in range(workers):
            inputs = Path(tmp, f"worker_{w}")
            inputs.mkdir()
            for n in lengths[w::workers]:  # longest first, dealt round-robin: the workers finish together
                seq = ("ACDEFGHIKLMNPQRSTVWY" * n)[:n]
                (inputs / f"bucket_{n}.a3m").write_text(f"#{n}\t1\n>101\n{seq}\n")
            procs.append(subprocess.Popen(
                [sys.executable, "-m", "colabfold.compile", *flags, "--model-type", model_type,
                 "--compile-cache", str(cache), "--worker", str(inputs), str(inputs) + "_out"], env=env))
        if any(p.wait() for p in procs):
            raise RuntimeError(f"compiling {model_type} failed")
    print(f"Compiled {model_type} for {len(lengths)} lengths up to {lengths[0]} into {cache}")


if __name__ == "__main__":
    main()
