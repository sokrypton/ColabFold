import os
ENV = {"TF_FORCE_UNIFIED_MEMORY":"1", "XLA_PYTHON_CLIENT_MEM_FRACTION":"4.0"}
for k,v in ENV.items():
    if k not in os.environ: os.environ[k] = v

# test if alphafold installed
try:
  import alphafold
except ModuleNotFoundError:
  raise RuntimeError("alphafold is not installed. Please run `pip install colabfold[alphafold]`")

import importlib_metadata
from pathlib import Path
import random

from colabfold.run_alphafold import run
from colabfold.utils import (
  DEFAULT_API_SERVER, ACCEPT_DEFAULT_TERMS,
  get_commit, setup_logging
)

from colabfold.inputs import (
  get_queries_pairwise, unpack_a3ms,
  parse_fasta, get_queries,
)

from colabfold.download import default_data_dir, download_alphafold_params

import logging
logger = logging.getLogger(__name__)

import argparse
def main():
  parser = argparse.ArgumentParser()
  parser.add_argument("input",
    default="input",
    help="Can be one of the following: "
    "Directory with fasta/a3m files, a csv/tsv file, a fasta file or an a3m file",
  )
  parser.add_argument("results", help="Directory to write the results to")

  # Main performance parameter
  parser.add_argument("--stop-at-score",
    help="Compute models until plddt (single chain) or ptmscore (complex) > threshold is reached. "
    "This can make colabfold much faster by only running the first model for easy queries.",
    type=float,
    default=100,
  )

  parser.add_argument("--num-recycles",
    help="Number of prediction recycles."
    "Increasing recycles can improve the quality but slows down the prediction.",
    type=int,
    default=None,
  )
  parser.add_argument("--recycle-early-stop-tolerance",
    help="Specify convergence criteria."
    "Run until the distance between recycles is within specified value.",
    type=float,
    default=None,
  )

  parser.add_argument("--num-ensemble",
    help="Number of ensembles."
    "The trunk of the network is run multiple times with different random choices for the MSA cluster centers.",
    type=int,
    default=1,
  )
  parser.add_argument("--num-seeds",
    help="Number of seeds to try. Will iterate from range(random_seed, random_seed+num_seeds)."
    ".",
    type=int,
    default=1,
  )
  parser.add_argument("--random-seed",
    help="Changing the seed for the random number generator can result in different structure predictions.",
    type=int,
    default=0,
  )
  parser.add_argument("--num-models", type=int, default=5, choices=[1, 2, 3, 4, 5])
  parser.add_argument("--recompile-padding",
    type=int,
    default=10,
    help="Whenever the input length changes, the model needs to be recompiled."
    "We pad sequences by specified length, so we can e.g. compute sequence from length 100 to 110 without recompiling."
    "The prediction will become marginally slower for the longer input, "
    "but overall performance increases due to not recompiling. "
    "Set to 0 to disable.",
  )
  parser.add_argument("--model-order", default="1,2,3,4,5", type=str)
  parser.add_argument("--host-url", default=DEFAULT_API_SERVER)
  parser.add_argument("--data")
  parser.add_argument("--msa-mode",
    default="mmseqs2_uniref_env",
    choices=[
      "mmseqs2_uniref_env",
      "mmseqs2_uniref",
      "single_sequence",
    ],
    help="Using an a3m file as input overwrites this option",
  )
  parser.add_argument("--model-type",
    help="predict strucutre/complex using the following model."
    'Auto will pick "alphafold2_ptm" for structure predictions and "alphafold2_multimer_v3" for complexes.',
    type=str,
    default="auto",
    choices=[
      "auto",
      "alphafold2",
      "alphafold2_ptm",
      "alphafold2_multimer_v1",
      "alphafold2_multimer_v2",
      "alphafold2_multimer_v3",
    ],
  )
  parser.add_argument("--num-relax",
    help="specify how many of the top ranked structures to relax using amber.",
    type=int,
    default=0,
  )
  parser.add_argument("--use-templates",
    default=False,
    action="store_true",
    help="Use templates from pdb"
  )
  parser.add_argument("--custom-template-path",
    type=str,
    default=None,
    help="Directory with pdb files to be used as input",
  )
  parser.add_argument("--rank_by",
    help="rank models by auto, plddt or ptmscore",
    type=str,
    default="auto",
    choices=["auto", "plddt", "ptm", "iptm", "multimer"],
  )
  parser.add_argument("--pair-mode",
    help="how to generate MSA for multimeric inputs: unpaired, paired, unpaired_paired",
    type=str,
    default="unpaired_paired",
    choices=["unpaired", "paired", "unpaired_paired"],
  )
  parser.add_argument("--sort-queries-by",
    help="sort queries by: none, length, random",
    type=str,
    default="length",
    choices=["none", "length", "random"],
  )
  parser.add_argument("--save-single-representations",
    default=False,
    action="store_true",
    help="saves the single representation embeddings of all models",
  )
  parser.add_argument("--save-pair-representations",
    default=False,
    action="store_true",
    help="saves the pair representation embeddings of all models",
  )
  parser.add_argument("--use-dropout",
    default=False,
    action="store_true",
    help="activate dropouts during inference to sample from uncertainity of the models",
  )
  parser.add_argument("--disable-masking",
    default=False,
    action="store_true",
    help='by default, 15% of the input MSA is randomly masked, set this flag to disable this',
  )
  parser.add_argument("--max-seq",
    help="number of sequence clusters to use",
    type=int,
    default=None,
  )
  parser.add_argument("--max-extra-seq",
    help="number of extra sequences to use",
    type=int,
    default=None,
  )
  parser.add_argument("--zip-results",
    default=False,
    action="store_true",
    help="zip all results into one <jobname>.result.zip and delete the original files",
  )
  parser.add_argument("--use-gpu-relax",
    default=False,
    action="store_true",
    help="run amber on GPU instead of CPU",
  )
  parser.add_argument("--save-all",
    default=False,
    action="store_true",
    help="save ALL raw outputs from model to a pickle file",
  )
  parser.add_argument("--save-recycles",
    default=False,
    action="store_true",
    help="save all intermediate predictions at each recycle",
  )
  parser.add_argument("--disable-unified-memory",
    default=False,
    action="store_true",
    help="if you are getting tensorflow/jax errors it might help to disable this",
  )  
  # undocumented arguements
  parser.add_argument("--overwrite-existing-results", default=False, action="store_true")
  parser.add_argument("--interaction-scan",           default=False, action="store_true")
  parser.add_argument("--disable-cluster-profile",    default=False, action="store_true")

  # backward compatability
  parser.add_argument('--training',    default=False, action="store_true", help=argparse.SUPPRESS)
  parser.add_argument('--templates',   default=False, action="store_true", help=argparse.SUPPRESS)
  parser.add_argument('--zip',         default=False, action="store_true", help=argparse.SUPPRESS)
  parser.add_argument('--amber',       default=False, action="store_true", help=argparse.SUPPRESS)
  parser.add_argument('--num-recycle', default=None,  type=int,            help=argparse.SUPPRESS)
  parser.add_argument("--max-msa",     default=None,  type=str,            help=argparse.SUPPRESS)

  # parse arguments
  args = parser.parse_args()

  # disable unified memory
  if args.disable_unified_memory:
    for k in ENV.keys():
      if k in os.environ: del os.environ[k]

  # backward compatability
  if args.training:  args.use_dropout   = True
  if args.templates: args.use_templates = True
  if args.zip:       args.zip_results   = True
  if args.amber and args.num_relax == 0: args.num_relax = args.num_models * args.num_seeds
  if args.num_recycle is not None: args.num_recycles = args.num_recycle
  if args.max_msa     is not None: (args.max_seq, args.max_extra_seq) = (int(x) for x in args.max_msa.split(":"))
  
  # setup logging
  setup_logging(Path(args.results).joinpath("log.txt"))
  version = importlib_metadata.version("colabfold")
  commit = get_commit()
  if commit: version += f" ({commit})"

  logger.info(f"Running colabfold {version}")
  data_dir = Path(args.data or default_data_dir)
  model_order = [int(i) for i in args.model_order.split(",")]
  assert args.recompile_padding >= 0, "Can't apply negative padding"

  # parse queries
  if args.interaction_scan:
    # protocol from @Dohyun-s
    batch_size = 10
    queries, is_complex, headers = get_queries_pairwise(args.input, batch_size)
  else:
    queries, is_complex = get_queries(args.input, args.sort_queries_by)

  # download params
  model_type = set_model_type(is_complex, args.model_type)
  download_alphafold_params(model_type, data_dir)
  
  # warning about api
  if "mmseqs2" in args.msa_mode or args.use_templates:
    # TODO: check if server used in the case of templates
    uses_api = any((query[2] is None for query in queries))
    if uses_api and args.host_url == DEFAULT_API_SERVER:
      print(ACCEPT_DEFAULT_TERMS, file=sys.stderr)

  run_params = dict(
    result_dir=args.results,
    use_templates=args.use_templates,
    custom_template_path=args.custom_template_path,
    num_relax=args.num_relax,
    msa_mode=args.msa_mode,
    model_type=model_type,
    num_models=args.num_models,
    num_recycles=args.num_recycles,
    recycle_early_stop_tolerance=args.recycle_early_stop_tolerance,
    num_ensemble=args.num_ensemble,
    model_order=model_order,
    is_complex=is_complex,
    keep_existing_results=not args.overwrite_existing_results,
    rank_by=args.rank_by,
    pair_mode=args.pair_mode,
    data_dir=data_dir,
    host_url=args.host_url,
    random_seed=args.random_seed,
    num_seeds=args.num_seeds,
    stop_at_score=args.stop_at_score,
    recompile_padding=args.recompile_padding,
    zip_results=args.zip_results,
    save_single_representations=args.save_single_representations,
    save_pair_representations=args.save_pair_representations,
    use_dropout=args.use_dropout,
    use_masking=not args.disable_masking,
    max_seq=args.max_seq,
    max_extra_seq=args.max_extra_seq,
    use_cluster_profile=not args.disable_cluster_profile,
    use_gpu_relax=args.use_gpu_relax,
    save_all=args.save_all,
    save_recycles=args.save_recycles,
  )

  if args.interaction_scan:

    # protocol from @Dohyun-s
    from colabfold.mmseqs.api import run_mmseqs2
    output = [queries[i:i + batch_size] for i in range(0, len(queries), batch_size)]
    headers_list = [headers[i:i + batch_size] for i in range(0, len(headers), batch_size)]
    headers_list[0].remove(headers_list[0][0])
    header_first = headers[0]
    
    for jobname, batch in enumerate(output):
      query_seqs_unique = []
      for x in batch:
        if x not in query_seqs_unique:
          query_seqs_unique.append(x)
      use_env = "env" in args.msa_mode or "Environmental" in args.msa_mode
      paired_a3m_lines = run_mmseqs2(
        query_seqs_unique,
        str(Path(args.results).joinpath(str(jobname))),
        use_env=use_env,
        use_pairwise=True,
        use_pairing=True,
        host_url=args.host_url,
      )
      
      path_o = Path(args.results).joinpath(f"{jobname}_pairwise")      
      for filenum in path_o.iterdir():
        queries_new = [] 
        if Path(filenum).suffix.lower() == ".a3m":
          outdir = path_o.joinpath("tmp")
          unpack_a3ms(filenum, outdir)
          for i, file in enumerate(sorted(outdir.iterdir())):
            if outdir.joinpath(file).suffix.lower() == ".a3m":
              (seqs, header) = parse_fasta(Path(file).read_text())
              query_sequence = seqs[0]
              a3m_lines = [Path(file).read_text()]
              val = int(header[0].split('\t')[1][1:]) - 102
              queries_new.append((header_first + '_' + headers_list[jobname][val], query_sequence, a3m_lines))

          if args.sort_queries_by == "length":
            queries_new.sort(key=lambda t: len(''.join(t[1])),reverse=True)
          elif args.sort_queries_by == "random":
            random.shuffle(queries_new)

          run(queries=queries_new, **run_params)
  
  else:
    run(queries=queries, **run_params)

if __name__ == "__main__":
  main()