"""
Functionality for running mmseqs locally. Takes in a fasta file, outputs final.a3m

Note: Currently needs mmseqs compiled from source
"""

import logging
import math
import os
import shutil
import subprocess
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
from pathlib import Path
from typing import List, Union

from colabfold.batch import get_queries, msa_to_str
from colabfold.utils import safe_filename

logger = logging.getLogger(__name__)


def run_mmseqs(mmseqs: Path, params: List[Union[str, Path]]):
    params_log = " ".join(str(i) for i in params)
    logger.info(f"Running {mmseqs} {params_log}")
    # hide MMseqs2 verbose paramters list that clogs up the log
    os.environ["MMSEQS_CALL_DEPTH"] = "1"
    subprocess.check_call([mmseqs] + params)


def mmseqs_search_monomer(
    dbbase: Path,
    base: Path,
    uniref_db: Path = Path("uniref30_2302_db"),
    template_db: Path = Path(""),  # Unused by default
    metagenomic_db: Path = Path("colabfold_envdb_202108_db"),
    mmseqs: Path = Path("mmseqs"),
    use_env: bool = True,
    use_templates: bool = False,
    filter: bool = True,
    expand_eval: float = math.inf,
    align_eval: int = 10,
    diff: int = 3000,
    qsc: float = -20.0,
    max_accept: int = 1000000,
    prefilter_mode: int = 0,
    s: float = 8,
    db_load_mode: int = 2,
    threads: int = 32,
):
    """Run mmseqs with a local colabfold database set

    db1: uniprot db (UniRef30)
    db2: Template (unused by default)
    db3: metagenomic db (colabfold_envdb_202108 or bfd_mgy_colabfold, the former is preferred)
    """
    if filter:
        # 0.1 was not used in benchmarks due to POSIX shell bug in line above
        #  EXPAND_EVAL=0.1
        align_eval = 10
        qsc = 0.8
        max_accept = 100000

    used_dbs = [uniref_db]
    if use_templates:
        used_dbs.append(template_db)
    if use_env:
        used_dbs.append(metagenomic_db)

    for db in used_dbs:
        if not dbbase.joinpath(f"{db}.dbtype").is_file():
            raise FileNotFoundError(f"Database {db} does not exist")
        if (
            (
                not dbbase.joinpath(f"{db}.idx").is_file()
                and not dbbase.joinpath(f"{db}.idx.index").is_file()
            )
            or os.environ.get("MMSEQS_IGNORE_INDEX", False)
        ):
            logger.info("Search does not use index")
            db_load_mode = 0
            dbSuffix1 = "_seq"
            dbSuffix2 = "_aln"
            dbSuffix3 = ""
        else:
            dbSuffix1 = ".idx"
            dbSuffix2 = ".idx"
            dbSuffix3 = ".idx"

    # fmt: off
    # @formatter:off
    search_param = ["--num-iterations", "3", "--db-load-mode", str(db_load_mode), "-a", "-e", "0.1", "--max-seqs", "10000"]
    search_param += ["--prefilter-mode", str(prefilter_mode)]
    if s is not None:
        search_param += ["-s", "{:.1f}".format(s)]
    else:
        search_param += ["--k-score", "'seq:96,prof:80'"]

    filter_param = ["--filter-msa", str(filter), "--filter-min-enable", "1000", "--diff", str(diff), "--qid", "0.0,0.2,0.4,0.6,0.8,1.0", "--qsc", "0", "--max-seq-id", "0.95",]
    expand_param = ["--expansion-mode", "0", "-e", str(expand_eval), "--expand-filter-clusters", str(filter), "--max-seq-id", "0.95",]

    run_mmseqs(mmseqs, ["search", base.joinpath("qdb"), dbbase.joinpath(uniref_db), base.joinpath("res"), base.joinpath("tmp"), "--threads", str(threads)] + search_param)
    run_mmseqs(mmseqs, ["mvdb", base.joinpath("tmp/latest/profile_1"), base.joinpath("prof_res")])
    run_mmseqs(mmseqs, ["lndb", base.joinpath("qdb_h"), base.joinpath("prof_res_h")])
    run_mmseqs(mmseqs, ["expandaln", base.joinpath("qdb"), dbbase.joinpath(f"{uniref_db}{dbSuffix1}"), base.joinpath("res"), dbbase.joinpath(f"{uniref_db}{dbSuffix2}"), base.joinpath("res_exp"), "--db-load-mode", str(db_load_mode), "--threads", str(threads)] + expand_param)
    run_mmseqs(mmseqs, ["align", base.joinpath("prof_res"), dbbase.joinpath(f"{uniref_db}{dbSuffix1}"), base.joinpath("res_exp"), base.joinpath("res_exp_realign"), "--db-load-mode", str(db_load_mode), "-e", str(align_eval), "--max-accept", str(max_accept), "--threads", str(threads), "--alt-ali", "10", "-a"])
    run_mmseqs(mmseqs, ["filterresult", base.joinpath("qdb"), dbbase.joinpath(f"{uniref_db}{dbSuffix1}"),
                        base.joinpath("res_exp_realign"), base.joinpath("res_exp_realign_filter"), "--db-load-mode",
                        str(db_load_mode), "--qid", "0", "--qsc", str(qsc), "--diff", "0", "--threads",
                        str(threads), "--max-seq-id", "1.0", "--filter-min-enable", "100"])
    run_mmseqs(mmseqs, ["result2msa", base.joinpath("qdb"), dbbase.joinpath(f"{uniref_db}{dbSuffix1}"),
                        base.joinpath("res_exp_realign_filter"), base.joinpath("uniref.a3m"), "--msa-format-mode",
                        "6", "--db-load-mode", str(db_load_mode), "--threads", str(threads)] + filter_param)
    subprocess.run([mmseqs] + ["rmdb", base.joinpath("res_exp_realign")])
    subprocess.run([mmseqs] + ["rmdb", base.joinpath("res_exp")])
    subprocess.run([mmseqs] + ["rmdb", base.joinpath("res")])
    subprocess.run([mmseqs] + ["rmdb", base.joinpath("res_exp_realign_filter")])

    if use_env:
        run_mmseqs(mmseqs, ["search", base.joinpath("prof_res"), dbbase.joinpath(metagenomic_db), base.joinpath("res_env"),
                            base.joinpath("tmp3"), "--threads", str(threads)] + search_param)
        run_mmseqs(mmseqs, ["expandaln", base.joinpath("prof_res"), dbbase.joinpath(f"{metagenomic_db}{dbSuffix1}"), base.joinpath("res_env"),
                            dbbase.joinpath(f"{metagenomic_db}{dbSuffix2}"), base.joinpath("res_env_exp"), "-e", str(expand_eval),
                            "--expansion-mode", "0", "--db-load-mode", str(db_load_mode), "--threads", str(threads)])
        run_mmseqs(mmseqs, ["align", base.joinpath("tmp3/latest/profile_1"), dbbase.joinpath(f"{metagenomic_db}{dbSuffix1}"),
                            base.joinpath("res_env_exp"), base.joinpath("res_env_exp_realign"), "--db-load-mode",
                            str(db_load_mode), "-e", str(align_eval), "--max-accept", str(max_accept), "--threads",
                            str(threads), "--alt-ali", "10", "-a"])
        run_mmseqs(mmseqs, ["filterresult", base.joinpath("qdb"), dbbase.joinpath(f"{metagenomic_db}{dbSuffix1}"),
                            base.joinpath("res_env_exp_realign"), base.joinpath("res_env_exp_realign_filter"),
                            "--db-load-mode", str(db_load_mode), "--qid", "0", "--qsc", str(qsc), "--diff", "0",
                            "--max-seq-id", "1.0", "--threads", str(threads), "--filter-min-enable", "100"])
        run_mmseqs(mmseqs, ["result2msa", base.joinpath("qdb"), dbbase.joinpath(f"{metagenomic_db}{dbSuffix1}"),
                            base.joinpath("res_env_exp_realign_filter"),
                            base.joinpath("bfd.mgnify30.metaeuk30.smag30.a3m"), "--msa-format-mode", "6",
                            "--db-load-mode", str(db_load_mode), "--threads", str(threads)] + filter_param)

        run_mmseqs(mmseqs, ["rmdb", base.joinpath("res_env_exp_realign_filter")])
        run_mmseqs(mmseqs, ["rmdb", base.joinpath("res_env_exp_realign")])
        run_mmseqs(mmseqs, ["rmdb", base.joinpath("res_env_exp")])
        run_mmseqs(mmseqs, ["rmdb", base.joinpath("res_env")])

        run_mmseqs(mmseqs, ["mergedbs", base.joinpath("qdb"), base.joinpath("final.a3m"), base.joinpath("uniref.a3m"), base.joinpath("bfd.mgnify30.metaeuk30.smag30.a3m")])
        run_mmseqs(mmseqs, ["rmdb", base.joinpath("bfd.mgnify30.metaeuk30.smag30.a3m")])
    else:
        run_mmseqs(mmseqs, ["mvdb", base.joinpath("uniref.a3m"), base.joinpath("final.a3m")])

    if use_templates:
        run_mmseqs(mmseqs, ["search", base.joinpath("prof_res"), dbbase.joinpath(template_db), base.joinpath("res_pdb"),
                            base.joinpath("tmp2"), "--db-load-mode", str(db_load_mode), "--threads", str(threads), "-s", "7.5", "-a", "-e", "0.1", "--prefilter-mode", str(prefilter_mode)])
        run_mmseqs(mmseqs, ["convertalis", base.joinpath("prof_res"), dbbase.joinpath(f"{template_db}{dbSuffix3}"), base.joinpath("res_pdb"),
                            base.joinpath(f"{template_db}"), "--format-output",
                            "query,target,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,cigar",
                            "--db-output", "1",
                            "--db-load-mode", str(db_load_mode), "--threads", str(threads)])
        run_mmseqs(mmseqs, ["unpackdb", base.joinpath(f"{template_db}"), base.joinpath("."), "--unpack-name-mode", "0", "--unpack-suffix", ".m8"])
        run_mmseqs(mmseqs, ["rmdb", base.joinpath("res_pdb")])
        run_mmseqs(mmseqs, ["rmdb", base.joinpath(f"{template_db}")])

    run_mmseqs(mmseqs, ["unpackdb", base.joinpath("final.a3m"), base.joinpath("."), "--unpack-name-mode", "0", "--unpack-suffix", ".a3m"])
    run_mmseqs(mmseqs, ["rmdb", base.joinpath("final.a3m")])
    run_mmseqs(mmseqs, ["rmdb", base.joinpath("uniref.a3m")])
    run_mmseqs(mmseqs, ["rmdb", base.joinpath("res")])
    # @formatter:on
    # fmt: on

    for file in base.glob("prof_res*"):
        file.unlink()
    shutil.rmtree(base.joinpath("tmp"))
    if use_templates:
        shutil.rmtree(base.joinpath("tmp2"))
    if use_env:
        shutil.rmtree(base.joinpath("tmp3"))


def mmseqs_search_pair(
    dbbase: Path,
    base: Path,
    uniref_db: Path = Path("uniref30_2302_db"),
    mmseqs: Path = Path("mmseqs"),
    prefilter_mode: int = 0,
    s: float = 8,
    threads: int = 64,
    db_load_mode: int = 2,
    pairing_strategy: int = 0,
):
    if not dbbase.joinpath(f"{uniref_db}.dbtype").is_file():
        raise FileNotFoundError(f"Database {uniref_db} does not exist")
    if (
        (
            not dbbase.joinpath(f"{uniref_db}.idx").is_file()
            and not dbbase.joinpath(f"{uniref_db}.idx.index").is_file()
        )
        or os.environ.get("MMSEQS_IGNORE_INDEX", False)
    ):
        logger.info("Search does not use index")
        db_load_mode = 0
        dbSuffix1 = "_seq"
        dbSuffix2 = "_aln"
    else:
        dbSuffix1 = ".idx"
        dbSuffix2 = ".idx"

    # fmt: off
    # @formatter:off
    search_param = ["--num-iterations", "3", "--db-load-mode", str(db_load_mode), "-a", "-e", "0.1", "--max-seqs", "10000",]
    search_param += ["--prefilter-mode", str(prefilter_mode)]
    if s is not None:
        search_param += ["-s", "{:.1f}".format(s)]
    else:
        search_param += ["--k-score", "'seq:96,prof:80'"]
    expand_param = ["--expansion-mode", "0", "-e", "inf", "--expand-filter-clusters", "0", "--max-seq-id", "0.95",]
    run_mmseqs(mmseqs, ["search", base.joinpath("qdb"), dbbase.joinpath(uniref_db), base.joinpath("res"), base.joinpath("tmp"), "--threads", str(threads),] + search_param,)
    run_mmseqs(mmseqs, ["expandaln", base.joinpath("qdb"), dbbase.joinpath(f"{uniref_db}{dbSuffix1}"), base.joinpath("res"), dbbase.joinpath(f"{uniref_db}{dbSuffix2}"), base.joinpath("res_exp"), "--db-load-mode", str(db_load_mode), "--threads", str(threads),] + expand_param,)
    run_mmseqs(mmseqs, ["align", base.joinpath("qdb"), dbbase.joinpath(f"{uniref_db}{dbSuffix1}"), base.joinpath("res_exp"), base.joinpath("res_exp_realign"), "--db-load-mode", str(db_load_mode), "-e", "0.001", "--max-accept", "1000000", "--threads", str(threads), "-c", "0.5", "--cov-mode", "1",],)
    run_mmseqs(mmseqs, ["pairaln", base.joinpath("qdb"), dbbase.joinpath(f"{uniref_db}"), base.joinpath("res_exp_realign"), base.joinpath("res_exp_realign_pair"), "--db-load-mode", str(db_load_mode), "--pairing-mode", str(pairing_strategy), "--pairing-dummy-mode", "0", "--threads", str(threads), ],)
    run_mmseqs(mmseqs, ["align", base.joinpath("qdb"), dbbase.joinpath(f"{uniref_db}{dbSuffix1}"), base.joinpath("res_exp_realign_pair"), base.joinpath("res_exp_realign_pair_bt"), "--db-load-mode", str(db_load_mode), "-e", "inf", "-a", "--threads", str(threads), ],)
    run_mmseqs(mmseqs, ["pairaln", base.joinpath("qdb"), dbbase.joinpath(f"{uniref_db}"), base.joinpath("res_exp_realign_pair_bt"), base.joinpath("res_final"), "--db-load-mode", str(db_load_mode), "--pairing-mode", str(pairing_strategy), "--pairing-dummy-mode", "1", "--threads", str(threads),],)
    run_mmseqs(mmseqs, ["result2msa", base.joinpath("qdb"), dbbase.joinpath(f"{uniref_db}{dbSuffix1}"), base.joinpath("res_final"), base.joinpath("pair.a3m"), "--db-load-mode", str(db_load_mode), "--msa-format-mode", "5", "--threads", str(threads),],)
    run_mmseqs(mmseqs, ["unpackdb", base.joinpath("pair.a3m"), base.joinpath("."), "--unpack-name-mode", "0", "--unpack-suffix", ".paired.a3m",],)
    run_mmseqs(mmseqs, ["rmdb", base.joinpath("qdb")])
    run_mmseqs(mmseqs, ["rmdb", base.joinpath("qdb_h")])
    run_mmseqs(mmseqs, ["rmdb", base.joinpath("res")])
    run_mmseqs(mmseqs, ["rmdb", base.joinpath("res_exp")])
    run_mmseqs(mmseqs, ["rmdb", base.joinpath("res_exp_realign")])
    run_mmseqs(mmseqs, ["rmdb", base.joinpath("res_exp_realign_pair")])
    run_mmseqs(mmseqs, ["rmdb", base.joinpath("res_exp_realign_pair_bt")])
    run_mmseqs(mmseqs, ["rmdb", base.joinpath("res_final")])
    run_mmseqs(mmseqs, ["rmdb", base.joinpath("pair.a3m")])
    shutil.rmtree(base.joinpath("tmp"))
    # @formatter:on
    # fmt: on


def main():
    parser = ArgumentParser(formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument(
        "query",
        type=Path,
        help="fasta files with the queries.",
    )
    parser.add_argument(
        "dbbase",
        type=Path,
        help="The path to the database and indices you downloaded and created with setup_databases.sh",
    )
    parser.add_argument(
        "base", type=Path, help="Directory for the results (and intermediate files)"
    )
    parser.add_argument(
        "--prefilter-mode",
        type=int,
        default=0,
        choices=[0, 1, 2],
        help="Prefiltering algorithm to use: 0: k-mer (high-mem), 1: ungapped (high-cpu), 2: exhaustive (no prefilter, very slow). See wiki for more details: https://github.com/sokrypton/ColabFold/wiki#colabfold_search",
    )
    parser.add_argument(
        "-s",
        type=float,
        default=None,
        help="MMseqs2 sensitivity. Lowering this will result in a much faster search but possibly sparser MSAs. By default, the k-mer threshold is directly set to the same one of the server, which corresponds to a sensitivity of ~8.",
    )
    # dbs are uniref, templates and environmental
    # We normally don't use templates
    parser.add_argument(
        "--db1", type=Path, default=Path("uniref30_2302_db"), help="UniRef database"
    )
    parser.add_argument("--db2", type=Path, default=Path(""), help="Templates database")
    parser.add_argument(
        "--db3",
        type=Path,
        default=Path("colabfold_envdb_202108_db"),
        help="Environmental database",
    )

    # poor man's boolean arguments
    parser.add_argument(
        "--use-env", type=int, default=1, choices=[0, 1], help="Use --db3"
    )
    parser.add_argument(
        "--use-templates", type=int, default=0, choices=[0, 1], help="Use --db2"
    )
    parser.add_argument(
        "--filter",
        type=int,
        default=1,
        choices=[0, 1],
        help="Filter the MSA by pre-defined align_eval, qsc, max_accept",
    )

    # mmseqs params
    parser.add_argument(
        "--mmseqs",
        type=Path,
        default=Path("mmseqs"),
        help="Location of the mmseqs binary.",
    )
    parser.add_argument(
        "--expand-eval",
        type=float,
        default=math.inf,
        help="e-val threshold for 'expandaln'.",
    )
    parser.add_argument(
        "--align-eval", type=int, default=10, help="e-val threshold for 'align'."
    )
    parser.add_argument(
        "--diff",
        type=int,
        default=3000,
        help="filterresult - Keep at least this many seqs in each MSA block.",
    )
    parser.add_argument(
        "--qsc",
        type=float,
        default=-20.0,
        help="filterresult - reduce diversity of output MSAs using min score thresh.",
    )
    parser.add_argument(
        "--max-accept",
        type=int,
        default=1000000,
        help="align - Maximum accepted alignments before alignment calculation for a query is stopped.",
    )
    parser.add_argument(
        "--pairing_strategy", type=int, default=0, help="pairaln - Pairing strategy."
    )
    parser.add_argument(
        "--db-load-mode",
        type=int,
        default=0,
        help="Database preload mode 0: auto, 1: fread, 2: mmap, 3: mmap+touch",
    )
    parser.add_argument(
        "--threads", type=int, default=64, help="Number of threads to use."
    )
    args = parser.parse_args()

    logging.basicConfig(level = logging.INFO)

    queries, is_complex = get_queries(args.query, None)

    queries_unique = []
    for job_number, (raw_jobname, query_sequences, a3m_lines) in enumerate(queries):
        # remove duplicates before searching
        query_sequences = (
            [query_sequences] if isinstance(query_sequences, str) else query_sequences
        )
        query_seqs_unique = []
        for x in query_sequences:
            if x not in query_seqs_unique:
                query_seqs_unique.append(x)
        query_seqs_cardinality = [0] * len(query_seqs_unique)
        for seq in query_sequences:
            seq_idx = query_seqs_unique.index(seq)
            query_seqs_cardinality[seq_idx] += 1

        queries_unique.append([raw_jobname, query_seqs_unique, query_seqs_cardinality])

    args.base.mkdir(exist_ok=True, parents=True)
    query_file = args.base.joinpath("query.fas")
    with query_file.open("w") as f:
        for job_number, (
            raw_jobname,
            query_sequences,
            query_seqs_cardinality,
        ) in enumerate(queries_unique):
            for j, seq in enumerate(query_sequences):
                # The header of first sequence set as 101
                query_seq_headername = 101 + j
                f.write(f">{query_seq_headername}\n{seq}\n")

    run_mmseqs(
        args.mmseqs,
        ["createdb", query_file, args.base.joinpath("qdb"), "--shuffle", "0"],
    )
    with args.base.joinpath("qdb.lookup").open("w") as f:
        id = 0
        file_number = 0
        for job_number, (
            raw_jobname,
            query_sequences,
            query_seqs_cardinality,
        ) in enumerate(queries_unique):
            for seq in query_sequences:
                raw_jobname_first = raw_jobname.split()[0]
                f.write(f"{id}\t{raw_jobname_first}\t{file_number}\n")
                id += 1
            file_number += 1

    mmseqs_search_monomer(
        mmseqs=args.mmseqs,
        dbbase=args.dbbase,
        base=args.base,
        uniref_db=args.db1,
        template_db=args.db2,
        metagenomic_db=args.db3,
        use_env=args.use_env,
        use_templates=args.use_templates,
        filter=args.filter,
        expand_eval=args.expand_eval,
        align_eval=args.align_eval,
        diff=args.diff,
        qsc=args.qsc,
        max_accept=args.max_accept,
        prefilter_mode=args.prefilter_mode,
        s=args.s,
        db_load_mode=args.db_load_mode,
        threads=args.threads,
    )
    if is_complex is True:
        mmseqs_search_pair(
            mmseqs=args.mmseqs,
            dbbase=args.dbbase,
            base=args.base,
            uniref_db=args.db1,
            prefilter_mode=args.prefilter_mode,
            s=args.s,
            db_load_mode=args.db_load_mode,
            threads=args.threads,
            pairing_strategy=args.pairing_strategy,
        )

        id = 0
        for job_number, (
            raw_jobname,
            query_sequences,
            query_seqs_cardinality,
        ) in enumerate(queries_unique):
            unpaired_msa = []
            paired_msa = None
            if len(query_seqs_cardinality) > 1:
                paired_msa = []
            for seq in query_sequences:
                with args.base.joinpath(f"{id}.a3m").open("r") as f:
                    unpaired_msa.append(f.read())
                args.base.joinpath(f"{id}.a3m").unlink()
                if len(query_seqs_cardinality) > 1:
                    with args.base.joinpath(f"{id}.paired.a3m").open("r") as f:
                        paired_msa.append(f.read())
                args.base.joinpath(f"{id}.paired.a3m").unlink()
                id += 1
            msa = msa_to_str(
                unpaired_msa, paired_msa, query_sequences, query_seqs_cardinality
            )
            args.base.joinpath(f"{job_number}.a3m").write_text(msa)

    # rename a3m files
    for job_number, (raw_jobname, query_sequences, query_seqs_cardinality) in enumerate(queries_unique):
        os.rename(
            args.base.joinpath(f"{job_number}.a3m"),
            args.base.joinpath(f"{safe_filename(raw_jobname)}.a3m"),
        )

    # rename m8 files
    if args.use_templates:
        id = 0
        for raw_jobname, query_sequences, query_seqs_cardinality in queries_unique:
            with args.base.joinpath(f"{safe_filename(raw_jobname)}_{args.db2}.m8").open(
                "w"
            ) as f:
                for _ in range(len(query_seqs_cardinality)):
                    with args.base.joinpath(f"{id}.m8").open("r") as g:
                        f.write(g.read())
                    os.remove(args.base.joinpath(f"{id}.m8"))
                    id += 1

    query_file.unlink()
    run_mmseqs(args.mmseqs, ["rmdb", args.base.joinpath("qdb")])
    run_mmseqs(args.mmseqs, ["rmdb", args.base.joinpath("qdb_h")])


if __name__ == "__main__":
    main()
