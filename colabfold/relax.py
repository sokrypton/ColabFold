from pathlib import Path

def relax_me(
    pdb_filename=None,
    pdb_lines=None,
    pdb_obj=None,
    use_gpu=False,
    max_iterations=0,
    tolerance=2.39,
    stiffness=10.0,
    max_outer_iterations=3
):
    from alphafold.common import protein
    from alphafold.relax import relax

    if pdb_obj is None:
        if pdb_lines is None:
            pdb_lines = Path(pdb_filename).read_text()
        pdb_obj = protein.from_pdb_string(pdb_lines)

    amber_relaxer = relax.AmberRelaxation(
        max_iterations=max_iterations,
        tolerance=tolerance,
        stiffness=stiffness,
        exclude_residues=[],
        max_outer_iterations=max_outer_iterations,
        use_gpu=use_gpu
    )

    relaxed_pdb_lines, _, _ = amber_relaxer.process(prot=pdb_obj)
    return relaxed_pdb_lines

def main():
    from argparse import ArgumentParser
    import os
    import glob
    from tqdm import tqdm

    parser = ArgumentParser()
    parser.add_argument("input",
        default="input",
        help="Can be one of the following: "
        "Directory with PDB files or a single PDB file",
    )
    parser.add_argument("results", help="Directory to write the results to or single output PDB file")
    parser.add_argument(
        "--max-iterations",
        type=int,
        default=2000,
        help="Maximum number of iterations for the relaxation process. AlphaFold2 sets this to unlimited (0), however, we found that this can lead to very long relaxation times for some inputs."
    )
    parser.add_argument(
        "--tolerance",
        type=float,
        default=2.39,
        help="tolerance level for the relaxation convergence"
    )
    parser.add_argument(
        "--stiffness",
        type=float,
        default=10.0,
        help="stiffness parameter for the relaxation"
    )
    parser.add_argument(
        "--max-outer-iterations",
        type=int,
        default=3,
        help="maximum number of outer iterations for the relaxation process"
    )
    parser.add_argument("--use-gpu",
        default=False,
        action="store_true",
        help="run amber on GPU instead of CPU",
    )
    args = parser.parse_args()

    input_path = Path(args.input)
    output_path = Path(args.results)
    if output_path.is_dir():
        output_path.mkdir(parents=True, exist_ok=True)

    if input_path.is_dir():
        pdb_files = glob.glob(str(input_path / "*.pdb"))
    else:
        pdb_files = [str(input_path)]

    if len(pdb_files) > 1:
        pdb_files = tqdm(pdb_files, desc="Processing PDB files")

    for pdb_file in pdb_files:
        relaxed_pdb = relax_me(
            pdb_filename=pdb_file,
            use_gpu=args.use_gpu,
            max_iterations=args.max_iterations,
            tolerance=args.tolerance,
            stiffness=args.stiffness,
            max_outer_iterations=args.max_outer_iterations
        )
        if output_path.is_dir():
            output_file = output_path / Path(pdb_file).name
        else:
            output_file = output_path

        with open(output_file, 'w') as file:
            file.write(relaxed_pdb)

if __name__ == "__main__":
    main()