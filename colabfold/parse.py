from typing import Any, Callable, Dict, List, Optional, Tuple, Union, TYPE_CHECKING
import numpy as np
import pandas
from pathlib import Path
from Bio.PDB import MMCIFParser, PDBParser, MMCIF2Dict

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

def a3m_queries(a3m_path):
  '''
  parse a3ms to get queries
  
  for monomer:
    #10 1
    >name
    AAAAAAAAAA
    >seq0
    AAAAAAAAAA

  for homooligomer
    #10 2
    >name
    AAAAAAAAAA
    >seq0
    AAAAAAAAAA

  for heterooligomer
    #5,5 1,1
    >name
    AAAAABBBBB
    >seq0
    AAAAABBBBB
    >seq1
    AAAAA-----
    >seq2
    -----BBBBB

  '''
  def get_cardinality(lines):  
    if lines[0].startswith("#"):
      query_seq_len, query_seqs_cardinality = \
      [[int(y) for y in x.split(",")] for x in lines[0][1:].split()]
      header, sequence = lines[1][1:], lines[2]
      query_seqs_unique = []
      n = 0
      for query_len in query_seq_len:
        query_seqs_unique.append(sequence[n:n+query_len])
        n += query_len
      return query_seqs_unique, query_seqs_cardinality, header
    
    else:
      header, sequence = lines[0][1:], lines[1]
      return [sequence], [1], header
    
  # split the a3m file
  queries, a3ms = [], []
  for a3m_line in a3m_path.read_text().splitlines():
    if a3m_line.startswith("#") or len(a3ms) == 0:
      a3ms.append([])
    a3ms[-1].append(a3m_line)      

  # for each a3m
  is_complex = False
  for a3m in a3ms:
    query_seqs_unique, query_seqs_cardinality, query_name = get_cardinality(a3m)
    query_sequence = []
    for seq, copies in zip(query_seqs_unique, query_seqs_cardinality):
      query_sequence += [seq] * copies
    if len(query_sequence) == 1:
      query_sequence = query_sequence[0]
    else:
      is_complex = True

    a3m_str = "\n".join(a3m)
    queries.append([query_name, query_sequence, [a3m_str]])
  
  return queries, is_complex

def fas_queries(fas_path):
  '''
  parse fasta to get queries
  
  for monomer:
    >name
    AAAAA

  for homooligomer
    >name
    AAAAA:AAAAA

  for heterooligomer
    >name
    AAAAA:BBBBB

  '''
  (sequences, headers) = parse_fasta(fas_path.read_text())
  queries, is_complex = [], False
  for sequence, header in zip(sequences, headers):
    sequence = sequence.upper()
    if sequence.count(":") == 0:
      queries.append([header, sequence, None])
    else:
      is_complex = True
      queries.append([header, sequence.upper().split(":"), None])
  return queries, is_complex

def csv_queries(csv_path):
  '''
  parse fasta to get queries
  
  for monomer:
    id,sequence
    name,AAAAA

  for homooligomer
    id,sequence
    name,AAAAA:AAAAA

  for heterooligomer
    id,sequence
    name,AAAAA:BBBBB

  '''
  sep = "\t" if csv_path.suffix == ".tsv" else ","
  df = pandas.read_csv(input_path, sep=sep)
  assert "id" in df.columns and "sequence" in df.columns
  queries, is_complex = [], False
  for header, sequence in df[["id", "sequence"]].itertuples(index=False):
    sequence = sequence.upper()
    if sequence.count(":") == 0:
      queries.append([header, sequence, None])
    else:
      is_complex = True
      queries.append([header, sequence.upper().split(":"), None])
  return queries, is_complex

def get_queries(
  input_path: Union[str, Path], sort_queries_by: str = "length"
) -> Tuple[List[Tuple[str, str, Optional[List[str]]]], bool]:
  """
  Reads a directory of fasta files, a single fasta file or a csv file and returns a tuple
  of job name, sequence and the optional a3m lines
  """
  input_path = Path(input_path)
  if not input_path.exists():
    raise OSError(f"{input_path} could not be found")

  def to_queries(file_path):
    if file_path.suffix in [".csv",".tsv"]:
      queries, is_complex = csv_queries(file_path)
    elif file_path.suffix == ".a3m":
      queries, is_complex = a3m_queries(file_path)
    elif file_path.suffix in [".fasta", ".faa", ".fa"]:
      queries, is_complex = fas_queries(file_path)
    else:
      raise ValueError(f"Unknown file format {file_path.suffix}")
    if len(queries) == 1: queries[0][0] = file_path.stem
    return queries, is_complex
    
  if input_path.is_file():
    queries, is_complex = to_queries(input_path)

  else:
    assert input_path.is_dir(), "Expected either an input file or a input directory"
    queries, is_complex = [], False
    for file in sorted(input_path.iterdir()):
      if not file.is_file():
        continue
      if file.suffix.lower() not in [".a3m", ".fasta", ".faa",".fa",".csv",".tsv"]:
        logger.warning(f"non-parsable file in input directory: {file}")
        continue
      queries_, is_complex_ = to_queries(file)
      queries += queries_
      if is_complex_: is_complex = True

  # sort by seq. len
  if sort_queries_by == "length":
    queries.sort(key=lambda t: len("".join(t[1])))

  elif sort_queries_by == "random":
    random.shuffle(queries)
  
  return queries, is_complex

def get_queries_pairwise(
  input_path: Union[str, Path], batch_size: int = 10,
) -> Tuple[List[Tuple[str, str, Optional[List[str]]]], bool]:
  """Reads a directory of fasta files, a single fasta file or a csv file and returns a tuple
  of job name, sequence and the optional a3m lines"""
  input_path = Path(input_path)
  if not input_path.exists():
    raise OSError(f"{input_path} could not be found")
  if input_path.is_file():
    if input_path.suffix == ".csv" or input_path.suffix == ".tsv":
      sep = "\t" if input_path.suffix == ".tsv" else ","
      df = pandas.read_csv(input_path, sep=sep)
      assert "id" in df.columns and "sequence" in df.columns
      queries = []
      seq_id_list = []
      for i, (seq_id, sequence) in enumerate(df[["id", "sequence"]].itertuples(index=False)):
        if i>0 and i % 10 == 0:
          queries.append(queries[0].upper())
        queries.append(sequence.upper())
        seq_id_list.append(seq_id)
      return queries, True, seq_id_list
    elif input_path.suffix == ".a3m":
      raise NotImplementedError()
    elif input_path.suffix in [".fasta", ".faa", ".fa"]:
      (sequences, headers) = parse_fasta(input_path.read_text())
      queries = []
      for i, (sequence, header) in enumerate(zip(sequences, headers)):
        sequence = sequence.upper()
        if sequence.count(":") == 0:
          if i>0 and i % 10 == 0:
            queries.append(sequences[0])
          queries.append(sequence)
        else:
          # Complex mode
          queries.append((header, sequence.upper().split(":"), None))
      return queries, True, headers
    else:
      raise ValueError(f"Unknown file format {input_path.suffix}")
  else:
    raise NotImplementedError()

def unserialize_msa(
  a3m_lines: List[str], query_sequence: Union[List[str], str]
) -> Tuple[
  Optional[List[str]],
  Optional[List[str]],
  List[str],
  List[int],
  List[Dict[str, Any]],
]:
  a3m_lines = a3m_lines[0].replace("\x00", "").splitlines()
  if not a3m_lines[0].startswith("#") or len(a3m_lines[0][1:].split("\t")) != 2:
    assert isinstance(query_sequence, str)
    return (
      ["\n".join(a3m_lines)],
      None,
      [query_sequence],
      [1],
      [mk_mock_template(query_sequence)],
    )

  if len(a3m_lines) < 3:
    raise ValueError(f"Unknown file format a3m")
  tab_sep_entries = a3m_lines[0][1:].split("\t")
  query_seq_len = tab_sep_entries[0].split(",")
  query_seq_len = list(map(int, query_seq_len))
  query_seqs_cardinality = tab_sep_entries[1].split(",")
  query_seqs_cardinality = list(map(int, query_seqs_cardinality))
  is_homooligomer = (
    True if len(query_seq_len) == 1 and query_seqs_cardinality[0] > 1 else False
  )
  is_single_protein = (
    True if len(query_seq_len) == 1 and query_seqs_cardinality[0] == 1 else False
  )
  query_seqs_unique = []
  prev_query_start = 0
  # we store the a3m with cardinality of 1
  for n, query_len in enumerate(query_seq_len):
    query_seqs_unique.append(
      a3m_lines[2][prev_query_start : prev_query_start + query_len]
    )
    prev_query_start += query_len
  paired_msa = [""] * len(query_seq_len)
  unpaired_msa = None
  already_in = dict()
  for i in range(1, len(a3m_lines), 2):
    header = a3m_lines[i]
    seq = a3m_lines[i + 1]
    if (header, seq) in already_in:
      continue
    already_in[(header, seq)] = 1
    has_amino_acid = [False] * len(query_seq_len)
    seqs_line = []
    prev_pos = 0
    for n, query_len in enumerate(query_seq_len):
      paired_seq = ""
      curr_seq_len = 0
      for pos in range(prev_pos, len(seq)):
        if curr_seq_len == query_len:
          prev_pos = pos
          break
        paired_seq += seq[pos]
        if seq[pos].islower():
          continue
        if seq[pos] != "-":
          has_amino_acid[n] = True
        curr_seq_len += 1
      seqs_line.append(paired_seq)

    # is sequence is paired add them to output
    if (
      not is_single_protein
      and not is_homooligomer
      and sum(has_amino_acid) == len(query_seq_len)
    ):
      header_no_faster = header.replace(">", "")
      header_no_faster_split = header_no_faster.split("\t")
      for j in range(0, len(seqs_line)):
        paired_msa[j] += ">" + header_no_faster_split[j] + "\n"
        paired_msa[j] += seqs_line[j] + "\n"
    else:
      unpaired_msa = [""] * len(query_seq_len)
      for j, seq in enumerate(seqs_line):
        if has_amino_acid[j]:
          unpaired_msa[j] += header + "\n"
          unpaired_msa[j] += seq + "\n"
  if is_homooligomer:
    # homooligomers
    num = 101
    paired_msa = [""] * query_seqs_cardinality[0]
    for i in range(0, query_seqs_cardinality[0]):
      paired_msa[i] = ">" + str(num + i) + "\n" + query_seqs_unique[0] + "\n"
  if is_single_protein:
    paired_msa = None
  template_features = []
  for query_seq in query_seqs_unique:
    template_feature = mk_mock_template(query_seq)
    template_features.append(template_feature)

  return (
    unpaired_msa,
    paired_msa,
    query_seqs_unique,
    query_seqs_cardinality,
    template_features,
  )

def unpack_a3ms(input, outdir):
  if not os.path.isdir(outdir):
     os.mkdir(outdir)
  i = 1
  count = 0
  with open(input, "r") as f:
    lines = f.readlines()
    a3m = ""
    hasLen1 = False
    len1 = 0
    inMsa2 = False
    hasLen2 = False
    switch = False
    len2 = 0
    left=[]
    right=[]
    descp_l=[]
    descp_r=[]
    for idx, line in enumerate(lines):
      if line.startswith("\x00"):
        i += 1
        line = line[1:]
        inMsa2 = True
        if i % 2 == 1:
          a3m="\n".join([l[:-1]+'\t'+r[:-1]+'\n'+a_[:-1]+b_[:-1]
           for l, r, a_, b_ in zip(descp_l, descp_r, left, right)])

          complex_description = "#" + str(len1) + "," + str(len2) + "\t1,1\n"
          a3m = complex_description + a3m
          out = Path(outdir) / f"job{count}.a3m"
          with open(out, "w") as w:
            w.write(a3m)
          a3m = ""
          count += 1
          hasLen1 = False
          hasLen2 = False
          inMsa2 = False
          left = []
          right = []
          descp_l = []
          descp_r = []
          switch = False
        else:
          switch = True

      if line.startswith(">") and hasLen2 is False:
        descp_l.append(line)
      if not line.startswith(">") and not hasLen1:
        len1 = len(line)-1
        hasLen1 = True
      if not line.startswith(">") and hasLen1 is True and hasLen2 is False:
        left.append(line)
      if line.startswith(">") and switch:
        descp_r.append(line)
      if not line.startswith(">") and not hasLen2 and inMsa2:
        len2 = len(line)-1
        hasLen2 = True
      if not line.startswith(">") and hasLen1 is True and hasLen2 is True:
        right.append(line)

#############
# TEMPLATES
#############
def validate_and_fix_mmcif(cif_file: Path):
  """validate presence of _entity_poly_seq in cif file and add revision_date if missing"""
  # check that required poly_seq and revision_date fields are present
  cif_dict = MMCIF2Dict.MMCIF2Dict(cif_file)
  required = [
    "_chem_comp.id",
    "_chem_comp.type",
    "_struct_asym.id",
    "_struct_asym.entity_id",
    "_entity_poly_seq.mon_id",
  ]
  for r in required:
    if r not in cif_dict:
      raise ValueError(f"mmCIF file {cif_file} is missing required field {r}.")
  if "_pdbx_audit_revision_history.revision_date" not in cif_dict:
    logger.info(
      f"Adding missing field revision_date to {cif_file}. Backing up original file to {cif_file}.bak."
    )
    shutil.copy2(cif_file, str(cif_file) + ".bak")
    with open(cif_file, "a") as f:
      f.write(CIF_REVISION_DATE)

def convert_pdb_to_mmcif(pdb_file: Path):
  """convert existing pdb files into mmcif with the required poly_seq and revision_date"""
  i = pdb_file.stem
  cif_file = pdb_file.parent.joinpath(f"{i}.cif")
  if cif_file.is_file():
    return
  parser = PDBParser(QUIET=True)
  structure = parser.get_structure(i, pdb_file)
  cif_io = CFMMCIFIO()
  cif_io.set_structure(structure)
  cif_io.save(str(cif_file))

def mk_hhsearch_db(template_dir: str):
  template_path = Path(template_dir)

  cif_files = template_path.glob("*.cif")
  for cif_file in cif_files:
    validate_and_fix_mmcif(cif_file)

  pdb_files = template_path.glob("*.pdb")
  for pdb_file in pdb_files:
    convert_pdb_to_mmcif(pdb_file)

  pdb70_db_files = template_path.glob("pdb70*")
  for f in pdb70_db_files:
    os.remove(f)

  with open(template_path.joinpath("pdb70_a3m.ffdata"), "w") as a3m, open(
    template_path.joinpath("pdb70_cs219.ffindex"), "w"
  ) as cs219_index, open(
    template_path.joinpath("pdb70_a3m.ffindex"), "w"
  ) as a3m_index, open(
    template_path.joinpath("pdb70_cs219.ffdata"), "w"
  ) as cs219:
    n = 1000000
    index_offset = 0
    cif_files = template_path.glob("*.cif")
    for cif_file in cif_files:
      with open(cif_file) as f:
        cif_string = f.read()
      cif_fh = StringIO(cif_string)
      parser = MMCIFParser(QUIET=True)
      structure = parser.get_structure("none", cif_fh)
      models = list(structure.get_models())
      if len(models) != 1:
        raise ValueError(
          f"Only single model PDBs are supported. Found {len(models)} models."
        )
      model = models[0]
      for chain in model:
        amino_acid_res = []
        for res in chain:
          if res.id[2] != " ":
            raise ValueError(
              f"PDB contains an insertion code at chain {chain.id} and residue "
              f"index {res.id[1]}. These are not supported."
            )
          amino_acid_res.append(
            residue_constants.restype_3to1.get(res.resname, "X")
          )

        protein_str = "".join(amino_acid_res)
        a3m_str = f">{cif_file.stem}_{chain.id}\n{protein_str}\n\0"
        a3m_str_len = len(a3m_str)
        a3m_index.write(f"{n}\t{index_offset}\t{a3m_str_len}\n")
        cs219_index.write(f"{n}\t{index_offset}\t{len(protein_str)}\n")
        index_offset += a3m_str_len
        a3m.write(a3m_str)
        cs219.write("\n\0")
        n += 1