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
from colabfold.utils import safe_filename, log_function_call

logger = logging.getLogger(__name__)


def run_mmseqs(mmseqs: Path, params: List[Union[str, Path]]):
    """
    Run mmseqs with the given parameters
    :param mmseqs: Path to the mmseqs binary
    :param params: List of parameters to pass to mmseqs
    :return:
    """
    params_log = " ".join(str(i) for i in params)
    logger.info(f"Running {mmseqs} {params_log}")
    # hide MMseqs2 verbose parameters list that clogs up the log
    os.environ["MMSEQS_CALL_DEPTH"] = "1"
    # Open a subprocess and direct stdout and stderr to subprocess.PIPE
    with subprocess.Popen([mmseqs] + params, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True,
                          bufsize=1) as proc:
        # Read the output line by line as it becomes available
        for line in proc.stdout:
            logger.info(line.strip())  # Log each line from the output
        # Wait for the subprocess to finish and get the exit code
        proc.wait()
        if proc.returncode != 0:
            raise subprocess.CalledProcessError(proc.returncode, [mmseqs] + params)


def safe_rename_files(base, old_suffix, new_suffix):
    # Iterate through every item in base directory
    for item in base.iterdir():
        # Check if the item is a file and has the correct extension
        if item.is_file() and item.suffix == old_suffix:
            # Construct new filename using safe_filename function
            new_filename = safe_filename(item.stem) + new_suffix
            # Rename the old file to the new one
            item.rename(base.joinpath(new_filename))


def rename_m8_files(queries_unique, base, template_db):
    id = 0
    for raw_jobname, query_sequences, query_seqs_cardinality in queries_unique:
        with base.joinpath(f"{safe_filename(raw_jobname)}_{template_db}.m8").open(
                "w"
        ) as f:
            for _ in range(len(query_seqs_cardinality)):
                with base.joinpath(f"{id}.m8").open("r") as g:
                    f.write(g.read())
                os.remove(base.joinpath(f"{id}.m8"))
                id += 1


def rename_a3m_files(queries_unique, base):
    logging.info("Renaming a3m files with args %s", queries_unique)
    for job_number, (raw_jobname, query_sequences, query_seqs_cardinality) in enumerate(queries_unique):
        logging.info(f"Renaming {base.joinpath(f'{job_number}.a3m')} to "
                     f"{base.joinpath(f'{safe_filename(raw_jobname)}.a3m')}")
        os.rename(
            base.joinpath(f"{job_number}.a3m"),
            base.joinpath(f"{safe_filename(raw_jobname)}.a3m"),
        )


def get_db_suffix(db, db_load_mode, dbtype="env"):
    if not Path(f"{db}.dbtype").is_file():
        raise FileNotFoundError(f"Database {db} does not exist")
    if (
        (
            not Path(f"{db}.idx").is_file()
            and not Path(f"{db}.idx.index").is_file()
        )
        or os.environ.get("MMSEQS_IGNORE_INDEX", False)
    ):
        logger.info("Search does not use index")
        db_load_mode = 0
        dbSuffix1 = "_seq" if dbtype == "env" else ""
        dbSuffix2 = "_aln" if dbtype == "env" else ""
    else:
        dbSuffix1 = ".idx"
        dbSuffix2 = ".idx"

    return dbSuffix1, dbSuffix2, db_load_mode


@log_function_call
def search_uniref_db(base, uniref_db, filter: int = 1, mmseqs: Path = Path("mmseqs"),align_eval: int = 10,
                     qsc: float = -20.0, max_accept: int = 1000000, db_load_mode: int = 0, threads: int = 32,
                     s: float = 8, diff: int = 3000, expand_eval: float = math.inf, prefilter_mode: int = 0):
    if filter:
        # 0.1 was not used in benchmarks due to POSIX shell bug in line above
        #  EXPAND_EVAL=0.1
        align_eval = 10
        qsc = 0.8
        max_accept = 100000

    dbSuffix1, dbSuffix2, db_load_mode = get_db_suffix(uniref_db, db_load_mode)

    search_param = ["--num-iterations", "3", "--db-load-mode", str(db_load_mode), "-a", "-e", "0.1", "--max-seqs",
                    "10000"]

    search_param += ["--prefilter-mode", str(prefilter_mode)]
    if s is not None:
        search_param += ["-s", "{:.1f}".format(s)]
    else:
        search_param += ["--k-score", "'seq:96,prof:80'"]
    filter_param = ["--filter-msa", str(filter), "--filter-min-enable", "1000", "--diff", str(diff), "--qid",
                    "0.0,0.2,0.4,0.6,0.8,1.0", "--qsc", "0", "--max-seq-id", "0.95",]
    expand_param = ["--expansion-mode", "0", "-e", str(expand_eval), "--expand-filter-clusters", str(filter),
                    "--max-seq-id", "0.95",]

    run_mmseqs(mmseqs, ["search", base.joinpath("qdb"), uniref_db, base.joinpath("res"), base.joinpath("tmp"),
                        "--threads", str(threads)] + search_param)
    run_mmseqs(mmseqs, ["mvdb", base.joinpath("tmp/latest/profile_1"), base.joinpath("prof_res")])
    run_mmseqs(mmseqs, ["lndb", base.joinpath("qdb_h"), base.joinpath("prof_res_h")])
    run_mmseqs(mmseqs, ["expandaln", base.joinpath("qdb"), f"{uniref_db}{dbSuffix1}", base.joinpath("res"),
                        f"{uniref_db}{dbSuffix2}", base.joinpath("res_exp"), "--db-load-mode", str(db_load_mode),
                        "--threads", str(threads)] + expand_param)
    run_mmseqs(mmseqs, ["align", base.joinpath("prof_res"), f"{uniref_db}{dbSuffix1}", base.joinpath("res_exp"),
                        base.joinpath("res_exp_realign"), "--db-load-mode", str(db_load_mode), "-e", str(align_eval),
                        "--max-accept", str(max_accept), "--threads", str(threads), "--alt-ali", "10", "-a"])
    run_mmseqs(mmseqs, ["filterresult", base.joinpath("qdb"), f"{uniref_db}{dbSuffix1}",
                        base.joinpath("res_exp_realign"), base.joinpath("res_exp_realign_filter"), "--db-load-mode",
                        str(db_load_mode), "--qid", "0", "--qsc", str(qsc), "--diff", "0", "--threads",
                        str(threads), "--max-seq-id", "1.0", "--filter-min-enable", "100"])
    run_mmseqs(mmseqs, ["result2msa", base.joinpath("qdb"), f"{uniref_db}{dbSuffix1}",
                        base.joinpath("res_exp_realign_filter"), base.joinpath(f"{uniref_db.name}.a3m"),
                        "--msa-format-mode", "6", "--db-load-mode", str(db_load_mode), "--threads",
                        str(threads)] + filter_param)
    subprocess.run([mmseqs] + ["rmdb", base.joinpath("res_exp_realign")])
    subprocess.run([mmseqs] + ["rmdb", base.joinpath("res_exp")])
    subprocess.run([mmseqs] + ["rmdb", base.joinpath("res")])
    subprocess.run([mmseqs] + ["rmdb", base.joinpath("res_exp_realign_filter")])
    shutil.rmtree(base.joinpath("tmp"))


@log_function_call
def search_env_db(base, metagenomic_db, mmseqs: Path = Path("mmseqs"), align_eval: int = 10, qsc: float = -20.0,
                  max_accept: int = 1000000, db_load_mode: int = 0, threads: int = 32, expand_eval: float = math.inf,
                  diff: int = 3000, filter: int = 1, s: float = 8, split_memory_limit: str = "100G"):

    if filter:
        # 0.1 was not used in benchmarks due to POSIX shell bug in line above
        #  EXPAND_EVAL=0.1
        align_eval = 10
        qsc = 0.8
        max_accept = 100000

    dbSuffix1, dbSuffix2, db_load_mode = get_db_suffix(metagenomic_db, db_load_mode)

    filter_param = ["--filter-msa", str(filter), "--filter-min-enable", "1000", "--diff", str(diff), "--qid",
                    "0.0,0.2,0.4,0.6,0.8,1.0", "--qsc", "0", "--max-seq-id", "0.95", ]
    search_param = ["--num-iterations", "3", "--db-load-mode", str(db_load_mode), "-a", "-e", "0.1", "--max-seqs",
                    "10000", "--split-memory-limit", split_memory_limit]
    if s:
        search_param += ["-s", "{:.1f}".format(s)]
    else:
        search_param += ["--k-score", "'seq:96,prof:80'"]

    run_mmseqs(mmseqs,
               ["search", base.joinpath("prof_res"), metagenomic_db, base.joinpath(f"res_env_{metagenomic_db.name}"),
                base.joinpath(f"tmp_{metagenomic_db.name}"), "--threads", str(threads)] + search_param)
    run_mmseqs(mmseqs,
               ["expandaln", base.joinpath("prof_res"), f"{metagenomic_db}{dbSuffix1}",
                base.joinpath(f"res_env_{metagenomic_db.name}"), f"{metagenomic_db}{dbSuffix2}",
                base.joinpath(f"res_env_{metagenomic_db.name}_exp"), "-e", str(expand_eval), "--expansion-mode", "0",
                "--db-load-mode", str(db_load_mode), "--threads", str(threads)])
    run_mmseqs(mmseqs, ["align", base.joinpath(f"tmp_{metagenomic_db.name}/latest/profile_1"),
                        f"{metagenomic_db}{dbSuffix1}", base.joinpath(f"res_env_{metagenomic_db.name}_exp"),
                        base.joinpath(f"res_env_{metagenomic_db.name}_exp_realign"), "--db-load-mode",
                        str(db_load_mode), "-e", str(align_eval), "--max-accept", str(max_accept), "--threads",
                        str(threads), "--alt-ali", "10", "-a"])
    run_mmseqs(mmseqs, ["filterresult", base.joinpath("qdb"), f"{metagenomic_db}{dbSuffix1}",
                        base.joinpath(f"res_env_{metagenomic_db.name}_exp_realign"),
                        base.joinpath(f"res_env_{metagenomic_db.name}_exp_realign_filter"),
                        "--db-load-mode", str(db_load_mode), "--qid", "0", "--qsc", str(qsc), "--diff", "0",
                        "--max-seq-id", "1.0", "--threads", str(threads), "--filter-min-enable", "100"])
    run_mmseqs(mmseqs, ["result2msa", base.joinpath("qdb"), f"{metagenomic_db}{dbSuffix1}",
                        base.joinpath(f"res_env_{metagenomic_db.name}_exp_realign_filter"),
                        base.joinpath(f"{metagenomic_db.name}.a3m"), "--msa-format-mode", "6",
                        "--db-load-mode", str(db_load_mode), "--threads", str(threads)] + filter_param)

    run_mmseqs(mmseqs, ["rmdb", base.joinpath(f"res_env_{metagenomic_db.name}_exp_realign_filter")])
    run_mmseqs(mmseqs, ["rmdb", base.joinpath(f"res_env_{metagenomic_db.name}_exp_realign")])
    run_mmseqs(mmseqs, ["rmdb", base.joinpath(f"res_env_{metagenomic_db.name}_exp")])
    run_mmseqs(mmseqs, ["rmdb", base.joinpath(f"res_env_{metagenomic_db.name}")])

    shutil.rmtree(base.joinpath(f"tmp_{metagenomic_db.name}"))


@log_function_call
def search_template_db(base, template_db, prefilter_mode: int = 0, mmseqs: Path = Path("mmseqs"), db_load_mode: int = 0,
                       threads: int = 32):

    dbSuffix, _, db_load_mode = get_db_suffix(template_db, db_load_mode, "template")

    run_mmseqs(mmseqs, ["search", base.joinpath("prof_res"), template_db, base.joinpath("res_pdb"),
                        base.joinpath("tmp2"), "--db-load-mode", str(db_load_mode), "--threads", str(threads), "-s",
                        "7.5", "-a", "-e", "0.1", "--prefilter-mode", str(prefilter_mode)])
    run_mmseqs(mmseqs, ["convertalis", base.joinpath("prof_res"), f"{template_db}{dbSuffix}",
                        base.joinpath("res_pdb"), base.joinpath(f"{template_db.name}"), "--format-output",
                        "query,target,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,cigar",
                        "--db-output", "1",
                        "--db-load-mode", str(db_load_mode), "--threads", str(threads)])


@log_function_call
def search_pair(
    base: Path,
    uniref_db: Path = Path("uniref30_2302_db"),
    mmseqs: Path = Path("mmseqs"),
    s: float = 8,
    threads: int = 64,
    db_load_mode: int = 2,
    pairing_strategy: int = 0,
):
    dbSuffix1, dbSuffix2, db_load_mode = get_db_suffix(uniref_db, db_load_mode)

    # fmt: off
    # @formatter:off
    search_param = ["--num-iterations", "3", "--db-load-mode", str(db_load_mode), "-a", "-e", "0.1", "--max-seqs",
                    "10000",]
    if s is not None:
        search_param += ["-s", "{:.1f}".format(s)]
    else:
        search_param += ["--k-score", "'seq:96,prof:80'"]
    expand_param = ["--expansion-mode", "0", "-e", "inf", "--expand-filter-clusters", "0", "--max-seq-id", "0.95",]
    run_mmseqs(mmseqs, ["search", base.joinpath("qdb"), uniref_db, base.joinpath("res"), base.joinpath("tmp"),
                        "--threads", str(threads),] + search_param,)
    run_mmseqs(mmseqs, ["expandaln", base.joinpath("qdb"), f"{uniref_db}{dbSuffix1}", base.joinpath("res"),
                        f"{uniref_db}{dbSuffix2}", base.joinpath("res_exp"), "--db-load-mode", str(db_load_mode),
                        "--threads", str(threads),] + expand_param,)
    run_mmseqs(mmseqs, ["align", base.joinpath("qdb"), f"{uniref_db}{dbSuffix1}", base.joinpath("res_exp"),
                        base.joinpath("res_exp_realign"), "--db-load-mode", str(db_load_mode), "-e", "0.001",
                        "--max-accept", "1000000", "--threads", str(threads), "-c", "0.5", "--cov-mode", "1",],)
    run_mmseqs(mmseqs, ["pairaln", base.joinpath("qdb"), f"{uniref_db}", base.joinpath("res_exp_realign"),
                        base.joinpath("res_exp_realign_pair"), "--db-load-mode", str(db_load_mode), "--pairing-mode",
                        str(pairing_strategy), "--pairing-dummy-mode", "0", "--threads", str(threads), ],)
    run_mmseqs(mmseqs, ["align", base.joinpath("qdb"), f"{uniref_db}{dbSuffix1}",
                        base.joinpath("res_exp_realign_pair"), base.joinpath("res_exp_realign_pair_bt"),
                        "--db-load-mode", str(db_load_mode), "-e", "inf", "-a", "--threads", str(threads), ],)
    run_mmseqs(mmseqs, ["pairaln", base.joinpath("qdb"), f"{uniref_db}",
                        base.joinpath("res_exp_realign_pair_bt"), base.joinpath("res_final"), "--db-load-mode",
                        str(db_load_mode), "--pairing-mode", str(pairing_strategy), "--pairing-dummy-mode", "1",
                        "--threads", str(threads),],)
    run_mmseqs(mmseqs, ["result2msa", base.joinpath("qdb"), f"{uniref_db}{dbSuffix1}",
                        base.joinpath("res_final"), base.joinpath("pair.a3m"), "--db-load-mode", str(db_load_mode),
                        "--msa-format-mode", "5", "--threads", str(threads),],)
    run_mmseqs(mmseqs, ["unpackdb", base.joinpath("pair.a3m"), base.joinpath("."), "--unpack-name-mode", "0",
                        "--unpack-suffix", ".paired.a3m",],)
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


@log_function_call
def merge_dbs(base, mmseqs: Path = Path("mmseqs")):
    msa_files = list(base.glob("*.a3m"))
    run_mmseqs(mmseqs, ["mergedbs", base.joinpath("qdb"), base.joinpath("final.a3m")] + msa_files)
    for msa_file in msa_files:
        run_mmseqs(mmseqs, ["rmdb", msa_file])


@log_function_call
def unpack_msa_files(base, mmseqs: Path = Path("mmseqs"), template_db_name: str = None, use_lookup: bool = False):
    if use_lookup:
        shutil.copyfile(base.joinpath("qdb.lookup"), base.joinpath("final.a3m.lookup"))
        run_mmseqs(mmseqs, ["unpackdb", base.joinpath("final.a3m"), base.joinpath("."), "--unpack-suffix", ".a3m"])
        if template_db_name:
            shutil.copyfile(base.joinpath("qdb.lookup"), base.joinpath(f"{template_db_name}.lookup"))
            run_mmseqs(mmseqs, ["unpackdb", base.joinpath(f"{template_db_name}"), base.joinpath("."),
                                "--unpack-suffix", ".m8"])
    else:
        run_mmseqs(mmseqs, ["unpackdb", base.joinpath("final.a3m"), base.joinpath("."), "--unpack-name-mode",
                            "0", "--unpack-suffix", ".a3m"])
        if template_db_name:
            run_mmseqs(mmseqs, ["unpackdb", base.joinpath(f"{template_db_name}"), base.joinpath("."),
                                "--unpack-name-mode", "0", "--unpack-suffix", ".m8"])

    run_mmseqs(mmseqs, ["rmdb", base.joinpath("final.a3m")])
    run_mmseqs(mmseqs, ["rmdb", base.joinpath("res")])
    run_mmseqs(mmseqs, ["rmdb", base.joinpath("res_pdb")])
    run_mmseqs(mmseqs, ["rmdb", base.joinpath(f"{template_db_name}")])
    shutil.rmtree(base.joinpath("tmp2"))
    for file in base.glob("prof_res*"):
        file.unlink()


@log_function_call
def rename_msa_files(base, queries_unique, template_db_name: str = None, use_lookup=False):
    if use_lookup:
        safe_rename_files(base, ".a3m", ".a3m")
        if template_db_name:
            safe_rename_files(base, ".m8", f"_{template_db_name}.m8")
    else:
        rename_a3m_files(queries_unique, base)
        if template_db_name:
            rename_m8_files(queries_unique, base, template_db_name)


def mmseqs_search_monomer(
    base: Path,
    uniref_db: Path,
    template_db: Path = None,  # Unused by default
    env_dbs: [Path] = None,
    mmseqs: Path = Path("mmseqs"),
    filter: int = 1,
    expand_eval: float = math.inf,
    align_eval: int = 10,
    diff: int = 3000,
    qsc: float = -20.0,
    max_accept: int = 1000000,
    s: float = 8,
    db_load_mode: int = 0,
    threads: int = 32,
    prefilter_mode: int = 0,
):
    """Run mmseqs with a local colabfold database set

    db1: uniprot db (UniRef30)
    db2: Template (unused by default)
    db3: metagenomic db (colabfold_envdb_202108 or bfd_mgy_colabfold, the former is preferred)
    """

    search_uniref_db(base=base, uniref_db=uniref_db, filter=filter, mmseqs=mmseqs, align_eval=align_eval, qsc=qsc,
                     max_accept=max_accept, db_load_mode=db_load_mode, threads=threads, s=s, diff=diff,
                     prefilter_mode=prefilter_mode)
    if env_dbs:
        for metagenomic_db in env_dbs:
            search_env_db(base=base, metagenomic_db=metagenomic_db, mmseqs=mmseqs, align_eval=align_eval, qsc=qsc,
                          max_accept=max_accept, db_load_mode=db_load_mode, threads=threads, expand_eval=expand_eval,
                          diff=diff, filter=filter)
        merge_dbs(base,  mmseqs)
    else:
        run_mmseqs(mmseqs, ["mvdb", base.joinpath(f"{uniref_db.name}.a3m"), base.joinpath("final.a3m")])

    if template_db:
        search_template_db(mmseqs=mmseqs, base=base, template_db=template_db, db_load_mode=db_load_mode,
                           threads=threads, prefilter_mode=prefilter_mode)


@log_function_call
def mmseqs_search_multimer(
    base: Path,
    uniref_db: Path,
    mmseqs: Path = Path("mmseqs"),
    s: float = 8,
    db_load_mode: int = 0,
    threads: int = 32,
    pairing_strategy: int = 0,
    queries_unique: list[list[Union[str, list[str], list[int]]]] = None,
):
    search_pair(
        mmseqs=mmseqs,
        base=base,
        uniref_db=uniref_db,
        s=s,
        db_load_mode=db_load_mode,
        threads=threads,
        pairing_strategy=pairing_strategy,
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
        for _ in query_sequences:
            with base.joinpath(f"{id}.a3m").open("r") as f:
                unpaired_msa.append(f.read())
            base.joinpath(f"{id}.a3m").unlink()
            if len(query_seqs_cardinality) > 1:
                with base.joinpath(f"{id}.paired.a3m").open("r") as f:
                    paired_msa.append(f.read())
            base.joinpath(f"{id}.paired.a3m").unlink()
            id += 1
        msa = msa_to_str(
            unpaired_msa, paired_msa, query_sequences, query_seqs_cardinality
        )
        base.joinpath(f"{job_number}.a3m").write_text(msa)


def create_query_db(query, base, mmseqs: Path = Path("mmseqs")):
    queries, is_complex = get_queries(query)

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

    base.mkdir(exist_ok=True, parents=True)
    query_file = base.joinpath("query.fas")
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
        mmseqs, ["createdb", query_file, base.joinpath("qdb"), "--shuffle", "0"],
    )
    with base.joinpath("qdb.lookup").open("w") as f:
        id = 0
        file_number = 0
        for job_number, (
                raw_jobname,
                query_sequences,
                query_seqs_cardinality,
        ) in enumerate(queries_unique):
            for _ in query_sequences:
                raw_jobname_first = raw_jobname.split()[0]
                f.write(f"{id}\t{raw_jobname_first}\t{file_number}\n")
                id += 1
            file_number += 1

    return is_complex, queries_unique, query_file


@log_function_call
def compute_msas(
    query: Path,
    base: Path,
    dbbase: Path = None,
    s: float = None,
    env_dbs: [Path] = None,
    uniref_db: Path = None,
    template_db: Path = None,
    filter: int = 1,
    mmseqs: Path = Path("mmseqs"),
    expand_eval: float = math.inf,
    align_eval: int = 10,
    diff: int = 3000,
    qsc: float = -20.0,
    max_accept: int = 1000000,
    pairing_strategy: int = 0,
    db_load_mode: int = 0,
    threads: int = 64,
    use_lookup: bool = False,
):

    is_complex, queries_unique, query_file = create_query_db(query, base, mmseqs)
    mmseqs_search_monomer(
        mmseqs=mmseqs,
        base=base,
        uniref_db=dbbase/uniref_db,
        template_db=dbbase/template_db if template_db else None,
        env_dbs=[dbbase/env_db for env_db in env_dbs] if env_dbs else None,
        filter=filter,
        expand_eval=expand_eval,
        align_eval=align_eval,
        diff=diff,
        qsc=qsc,
        max_accept=max_accept,
        s=s,
        db_load_mode=db_load_mode,
        threads=threads,
    )

    unpack_msa_files(base=base, mmseqs=mmseqs, template_db_name=(dbbase / template_db).name if template_db else None,
                     use_lookup=use_lookup)

    if is_complex is True:
        mmseqs_search_multimer(
            base=base,
            uniref_db=dbbase/uniref_db,
            mmseqs=mmseqs,
            s=s,
            db_load_mode=db_load_mode,
            threads=threads,
            pairing_strategy=pairing_strategy,
            queries_unique=queries_unique,
        )

    rename_msa_files(base=base, queries_unique=queries_unique,
                     template_db_name=(dbbase / template_db).name if template_db else None,
                     use_lookup=use_lookup)

    query_file.unlink()
    run_mmseqs(mmseqs, ["rmdb", base.joinpath("qdb")])
    run_mmseqs(mmseqs, ["rmdb", base.joinpath("qdb_h")])


def parse_arguments():
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
        help="Prefiltering algorithm to use: 0: k-mer (high-mem), 1: ungapped (high-cpu), "
             "2: exhaustive (no prefilter, very slow). See wiki for more details: "
             "https://github.com/sokrypton/ColabFold/wiki#colabfold_search",
    )
    parser.add_argument(
        "-s",
        type=float,
        default=None,
        help="MMseqs2 sensitivity. Lowering this will result in a much faster search but possibly sparser MSAs. "
             "By default, the k-mer threshold is directly set to the same one of the server, "
             "which corresponds to a sensitivity of ~8.",
    )

    parser.add_argument("--uniref-db", type=Path, default=Path("uniref30_2302_db"),
                        help="UniRef database")
    parser.add_argument("--template-db", type=Path, default=Path(""),
                        help="Templates database")
    parser.add_argument("--env-dbs", type=Path, nargs='+', default=[Path("colabfold_envdb_202108_db")],
                        help="Environmental databases")
    # Backwards compatibility
    # dbs are uniref, templates and environmental
    # We normally don't use templates
    parser.add_argument("--db1", type=Path, help="UniRef database")
    parser.add_argument("--db2", type=Path, help="Templates database")
    parser.add_argument("--db3", type=Path, help="Environmental database")
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

    # Backwards compatibility
    if args.db1:
        args.uniref_db = args.db1

    if args.use_templates:
        if args.db2:
            args.template_db = args.db2
    else:
        args.template_db = None

    if args.use_env:
        if args.db3:
            args.env_dbs = [args.db3]
    else:
        args.env_dbs = None

    return args


def main():
    args = parse_arguments()
    compute_msas(query=args.query,
                 base=args.base,
                 dbbase=args.dbbase,
                 s=args.s,
                 uniref_db=args.uniref_db,
                 template_db=args.template_db,
                 env_dbs=args.env_dbs,
                 filter=args.filter,
                 mmseqs=args.mmseqs,
                 expand_eval=args.expand_eval,
                 align_eval=args.align_eval,
                 diff=args.diff,
                 qsc=args.qsc,
                 max_accept=args.max_accept,
                 pairing_strategy=args.pairing_strategy,
                 db_load_mode=args.db_load_mode,
                 threads=args.threads,
                 )


if __name__ == "__main__":
    main()
