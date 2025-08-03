from typing import List, Optional, Tuple, Union, Dict
from pathlib import Path
import random
import logging
from colabfold.utils import MolType
logger = logging.getLogger(__name__)

def safe_filename(file: str) -> str:
    return "".join([c if c.isalnum() or c in ["_", ".", "-"] else "_" for c in file])

def pair_sequences(
    a3m_lines: List[str], query_sequences: List[str], query_cardinality: List[int]
) -> str:
    a3m_line_paired = [""] * len(a3m_lines[0].splitlines())
    for n, seq in enumerate(query_sequences):
        lines = a3m_lines[n].splitlines()
        for i, line in enumerate(lines):
            if line.startswith(">"):
                if n != 0:
                    line = line.replace(">", "\t", 1)
                a3m_line_paired[i] = a3m_line_paired[i] + line
            else:
                a3m_line_paired[i] = a3m_line_paired[i] + line * query_cardinality[n]
    return "\n".join(a3m_line_paired)

def pad_sequences(
    a3m_lines: List[str], query_sequences: List[str], query_cardinality: List[int]
) -> str:
    _blank_seq = [
        ("-" * len(seq))
        for n, seq in enumerate(query_sequences)
        for _ in range(query_cardinality[n])
    ]
    a3m_lines_combined = []
    pos = 0
    for n, seq in enumerate(query_sequences):
        for j in range(0, query_cardinality[n]):
            lines = a3m_lines[n].split("\n")
            for a3m_line in lines:
                if len(a3m_line) == 0:
                    continue
                if a3m_line.startswith(">"):
                    a3m_lines_combined.append(a3m_line)
                else:
                    a3m_lines_combined.append(
                        "".join(_blank_seq[:pos] + [a3m_line] + _blank_seq[pos + 1 :])
                    )
            pos += 1
    return "\n".join(a3m_lines_combined)

def pair_msa(
    query_seqs_unique: List[str],
    query_seqs_cardinality: List[int],
    paired_msa: Optional[List[str]],
    unpaired_msa: Optional[List[str]],
) -> str:
    if paired_msa is None and unpaired_msa is not None:
        a3m_lines = pad_sequences(
            unpaired_msa, query_seqs_unique, query_seqs_cardinality
        )
    elif paired_msa is not None and unpaired_msa is not None:
        a3m_lines = (
            pair_sequences(paired_msa, query_seqs_unique, query_seqs_cardinality)
            + "\n"
            + pad_sequences(unpaired_msa, query_seqs_unique, query_seqs_cardinality)
        )
    elif paired_msa is not None and unpaired_msa is None:
        a3m_lines = pair_sequences(
            paired_msa, query_seqs_unique, query_seqs_cardinality
        )
    else:
        raise ValueError(f"Invalid pairing")
    return a3m_lines

def msa_to_str(
    unpaired_msa: List[str],
    paired_msa: List[str],
    query_seqs_unique: List[str],
    query_seqs_cardinality: List[int],
) -> str:
    msa = "#" + ",".join(map(str, map(len, query_seqs_unique))) + "\t"
    msa += ",".join(map(str, query_seqs_cardinality)) + "\n"
    # build msa with cardinality of 1, it makes it easier to parse and manipulate
    query_seqs_cardinality = [1 for _ in query_seqs_cardinality]
    msa += pair_msa(query_seqs_unique, query_seqs_cardinality, paired_msa, unpaired_msa)
    return msa

def parse_fasta(fasta_string: str) -> Tuple[List[str], List[str]]:
    """Parses FASTA string and returns list of strings with amino-acid sequences.

    Arguments:
      fasta_string: The string contents of a FASTA file.

    Returns:
      A tuple of two lists:
      * A list of sequences.
      * A list of sequence descriptions taken from the comment lines. In the
        same order as the sequences.
    """
    sequences = []
    descriptions = []
    index = -1
    for line in fasta_string.splitlines():
        line = line.strip()
        if line.startswith("#"):
            continue
        if line.startswith(">"):
            index += 1
            descriptions.append(line[1:])  # Remove the '>' at the beginning.
            sequences.append("")
            continue
        elif not line:
            continue  # Skip blank lines.
        sequences[index] += line

    return sequences, descriptions

def classify_molecules(query_sequence: str) -> Tuple[List[str], Optional[List[Tuple[MolType, str, int]]]]:
    """Classifies the sequences in the query sequence string into protein and non-protein sequences.

    Returns a tuple of two lists:
    * A list of protein sequences.
    * A list of tuples, each containing a molecule type, a sequence, and number of copies.
    """
    sequences = query_sequence.upper().split(":")
    protein_queries = []
    other_queries = []
    for seq in sequences:
        if seq.count("|") == 0:
            protein_queries.append(seq)
        else:
            parts = seq.split("|")
            moltype, sequence, *rest = parts
            moltype = MolType.get_moltype(moltype)
            if moltype == MolType.SMILES:
                sequence = sequence.replace(";", ":")
            copies = int(rest[0]) if rest else 1
            other_queries.append((moltype, sequence, copies))  # (molecule type, sequence, copies)

    if len(other_queries) == 0:
        other_queries = None

    return protein_queries, other_queries

modified_mapping = {
    "MSE" : "MET", "MLY" : "LYS", "FME" : "MET", "HYP" : "PRO",
    "TPO" : "THR", "CSO" : "CYS", "SEP" : "SER", "M3L" : "LYS",
    "HSK" : "HIS", "SAC" : "SER", "PCA" : "GLU", "DAL" : "ALA",
    "CME" : "CYS", "CSD" : "CYS", "OCS" : "CYS", "DPR" : "PRO",
    "B3K" : "LYS", "ALY" : "LYS", "YCM" : "CYS", "MLZ" : "LYS",
    "4BF" : "TYR", "KCX" : "LYS", "B3E" : "GLU", "B3D" : "ASP",
    "HZP" : "PRO", "CSX" : "CYS", "BAL" : "ALA", "HIC" : "HIS",
    "DBZ" : "ALA", "DCY" : "CYS", "DVA" : "VAL", "NLE" : "LEU",
    "SMC" : "CYS", "AGM" : "ARG", "B3A" : "ALA", "DAS" : "ASP",
    "DLY" : "LYS", "DSN" : "SER", "DTH" : "THR", "GL3" : "GLY",
    "HY3" : "PRO", "LLP" : "LYS", "MGN" : "GLN", "MHS" : "HIS",
    "TRQ" : "TRP", "B3Y" : "TYR", "PHI" : "PHE", "PTR" : "TYR",
    "TYS" : "TYR", "IAS" : "ASP", "GPL" : "LYS", "KYN" : "TRP",
    "CSD" : "CYS", "SEC" : "CYS"
}

restype_1to3 = {
    'A': 'ALA',
    'R': 'ARG',
    'N': 'ASN',
    'D': 'ASP',
    'C': 'CYS',
    'Q': 'GLN',
    'E': 'GLU',
    'G': 'GLY',
    'H': 'HIS',
    'I': 'ILE',
    'L': 'LEU',
    'K': 'LYS',
    'M': 'MET',
    'F': 'PHE',
    'P': 'PRO',
    'S': 'SER',
    'T': 'THR',
    'W': 'TRP',
    'Y': 'TYR',
    'V': 'VAL',
}
restype_3to1 = {v: k for k, v in restype_1to3.items()}

def pdb_to_string(
        pdb_file: str,
        chains: Optional[str] = None,
        models: Optional[list] = None,
    ) -> str:
    '''read pdb file and return as string'''

    if chains is not None:
        if "," in chains: chains = chains.split(",")
        if not isinstance(chains,list): chains = [chains]
    if models is not None:
        if not isinstance(models,list): models = [models]

    modres = {**modified_mapping}
    lines = []
    seen = []
    model = 1

    if "\n" in pdb_file:
        old_lines = pdb_file.split("\n")
    else:
        with open(pdb_file,"rb") as f:
          old_lines = [line.decode("utf-8","ignore").rstrip() for line in f]
    for line in old_lines:
        if line[:5] == "MODEL":
            model = int(line[5:])
        if models is None or model in models:
            if line[:6] == "MODRES":
                k = line[12:15]
                v = line[24:27]
                if k not in modres and v in restype_3to1:
                    modres[k] = v
            if line[:6] == "HETATM":
                k = line[17:20]
                if k in modres:
                    line = "ATOM  "+line[6:17]+modres[k]+line[20:]
            if line[:4] == "ATOM":
                chain = line[21:22]
                if chains is None or chain in chains:
                    atom = line[12:12+4].strip()
                    resi = line[17:17+3]
                    resn = line[22:22+5].strip()
                    if resn[-1].isalpha(): # alternative atom
                        resn = resn[:-1]
                        line = line[:26]+" "+line[27:]
                    key = f"{model}_{chain}_{resn}_{resi}_{atom}"
                    if key not in seen: # skip alternative placements
                        lines.append(line)
                        seen.append(key)
            if line[:5] == "MODEL" or line[:3] == "TER" or line[:6] == "ENDMDL":
                lines.append(line)
    return "\n".join(lines)

restypes = [
    'A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P',
    'S', 'T', 'W', 'Y', 'V'
]
restypes_with_x = restypes + ['X']
restype_order_with_x = {restype: i for i, restype in enumerate(restypes_with_x)}
order_to_restype = {v: k for k, v in restype_order_with_x.items()}
def decode_structure_sequences(
    aatype_array: List[int],
    chain_index_array: List[int],
    order_dict: Dict[int, str] = order_to_restype
) -> List[str]:
    decoded_sequences = []
    current_sequence = []

    for i in range(len(aatype_array)):
        amino_acid = order_dict[aatype_array[i]]
        if i == 0 or chain_index_array[i] == chain_index_array[i - 1]:
            current_sequence.append(amino_acid)
        else:
            decoded_sequences.append("".join(current_sequence))
            current_sequence = [amino_acid]

    # Append the last sequence
    decoded_sequences.append("".join(current_sequence))

    return decoded_sequences

def get_queries(
    input_path: Union[str, Path], sort_queries_by: str = "length"
) -> Tuple[List[Tuple[str, str, Optional[List[str]], Optional[List[Tuple[MolType, str, int]]]]], bool]:
    """Reads a directory of fasta files, a single fasta file or a csv file and returns a tuple
    of job name, sequence, optional a3m lines, and the optional non-protein sequences."""

    input_path = Path(input_path)
    if not input_path.exists():
        raise OSError(f"{input_path} could not be found")

    if input_path.is_file():
        if input_path.suffix == ".csv" or input_path.suffix == ".tsv":
            sep = "\t" if input_path.suffix == ".tsv" else ","
            import pandas
            df = pandas.read_csv(input_path, sep=sep, dtype=str)
            assert "id" in df.columns and "sequence" in df.columns
            queries = [
                (seq_id, sequence.upper().split(":"), None, None)
                for seq_id, sequence in df[["id", "sequence"]].itertuples(index=False)
            ]
            for i in range(len(queries)):
                if len(queries[i][1]) == 1:
                    queries[i] = (queries[i][0], queries[i][1][0], None, None)
        elif input_path.suffix == ".a3m":
            (seqs, header) = parse_fasta(input_path.read_text())
            if len(seqs) == 0:
                raise ValueError(f"{input_path} is empty")
            query_sequence = seqs[0]
            # Use a list so we can easily extend this to multiple msas later
            a3m_lines = [input_path.read_text()]
            queries = [(input_path.stem, query_sequence, a3m_lines, None)]
        elif input_path.suffix in [".fasta", ".faa", ".fa"]:
            (sequences, headers) = parse_fasta(input_path.read_text())
            queries = []
            for sequence, header in zip(sequences, headers):
                sequence = sequence.upper()
                if sequence.count(":") == 0:
                    # Single sequence
                    queries.append((header, sequence, None, None))
                else:
                    # Complex mode
                    protein_queries, other_queries = classify_molecules(sequence)
                    queries.append((header, protein_queries, None, other_queries))
        elif input_path.suffix in [".pdb", ".cif"]:
            from alphafold.common import protein
            if input_path.suffix == ".pdb":
                pdb_string = pdb_to_string(input_path.read_text())
                prot = protein.from_pdb_string(pdb_string)
            elif input_path.suffix == ".cif":
                prot = protein.from_mmcif_string(input_path.read_text())
            header = input_path.stem
            sequences = decode_structure_sequences(prot.aatype, prot.chain_index)

            if len(sequences) == 0:
                raise ValueError(f"{input_path} is empty")

            queries = [(header, sequences, None, None)]

        else:
            raise ValueError(f"Unknown file format {input_path.suffix}")
    else:
        assert input_path.is_dir(), "Expected either an input file or a input directory"
        queries = []
        for file in sorted(input_path.iterdir()):
            if not file.is_file():
                continue
            if file.suffix.lower() not in [".a3m", ".fasta", ".faa", ".fa", ".pdb", ".cif"]:
                logger.warning(f"non-fasta/a3m/pdb/cif file in input directory: {file}")
                continue
            if file.suffix.lower() in [".pdb", ".cif"]:
                header = file.stem
                if file.suffix.lower() == ".pdb":
                    pdb_string = pdb_to_string(file.read_text())
                    prot = protein.from_pdb_string(pdb_string)
                else:  # file.suffix.lower() == ".cif"
                    prot = protein.from_mmcif_string(file.read_text())
                sequences = decode_structure_sequences(prot.aatype, prot.chain_index)

                if len(sequences) == 0:
                    logger.error(f"{file} is empty")
                    continue

                queries.append((header, sequences, None))
            else:  # file.suffix.lower() in [".a3m", ".fasta", ".faa"]
                (seqs, header) = parse_fasta(file.read_text())
            if len(seqs) == 0:
                logger.error(f"{file} is empty")
                continue
            query_sequence = seqs[0]
            if len(seqs) > 1 and file.suffix in [".fasta", ".faa", ".fa"]:
                logger.warning(
                    f"More than one sequence in {file}, ignoring all but the first sequence"
                )

            if file.suffix.lower() == ".a3m":
                a3m_lines = [file.read_text()]
                queries.append((file.stem, query_sequence.upper(), a3m_lines, None))
            else:
                if query_sequence.count(":") == 0:
                    # Single sequence
                    queries.append((file.stem, query_sequence, None, None))
                else:
                    # Complex mode
                    protein_queries, other_queries = classify_molecules(query_sequence)
                    queries.append((file.stem, protein_queries, None, other_queries))

    # sort by seq. len
    if sort_queries_by == "length":
        queries.sort(key=lambda t: len("".join(t[1])))

    elif sort_queries_by == "random":
        random.shuffle(queries)

    is_complex = False
    for job_number, (_, query_sequence, a3m_lines, _) in enumerate(queries):
        if isinstance(query_sequence, list):
            is_complex = True
            break
        if a3m_lines is not None and a3m_lines[0].startswith("#"):
            a3m_line = a3m_lines[0].splitlines()[0]
            tab_sep_entries = a3m_line[1:].split("\t")
            if len(tab_sep_entries) == 2:
                query_seq_len = tab_sep_entries[0].split(",")
                query_seq_len = list(map(int, query_seq_len))
                query_seqs_cardinality = tab_sep_entries[1].split(",")
                query_seqs_cardinality = list(map(int, query_seqs_cardinality))
                is_single_protein = (
                    True
                    if len(query_seq_len) == 1 and query_seqs_cardinality[0] == 1
                    else False
                )
                if not is_single_protein:
                    is_complex = True
                    break
    return queries, is_complex
