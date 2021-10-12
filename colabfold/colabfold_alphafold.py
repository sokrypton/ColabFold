# fmt: off
# @formatter:off

import os

from urllib import request
from concurrent import futures
import pickle

import jax
from alphafold.data.tools import jackhmmer
from alphafold.data import parsers
from alphafold.data import pipeline
from alphafold.common import protein
from alphafold.model import config
from alphafold.model import model
from alphafold.model import data
from alphafold.model.tf import shape_placeholders

import tensorflow as tf

from string import ascii_uppercase

import numpy as np
import matplotlib.pyplot as plt

import re
import colabfold as cf
try:
  import pairmsa
except:
  pairmsa=None

try:
  from google.colab import files
  IN_COLAB = True
except:
  IN_COLAB = False

import tqdm.notebook
TQDM_BAR_FORMAT = '{l_bar}{bar}| {n_fmt}/{total_fmt} [elapsed: {elapsed} remaining: {remaining}]'

#######################################################################################################################################
# prep_inputs
#######################################################################################################################################

def prep_inputs(sequence, jobname="test", homooligomer="1", output_dir=None, clean=False, verbose=True):
  # process inputs
  sequence = str(sequence)
  sequence = re.sub("[^A-Z:/]", "", sequence.upper())
  sequence = re.sub(":+",":",sequence)
  sequence = re.sub("/+","/",sequence)
  sequence = re.sub("^[:/]+","",sequence)
  sequence = re.sub("[:/]+$","",sequence)
  jobname = re.sub(r'\W+', '', jobname)
  homooligomer = str(homooligomer)
  homooligomer = re.sub("[:/]+",":",homooligomer)
  homooligomer = re.sub("^[:/]+","",homooligomer)
  homooligomer = re.sub("[:/]+$","",homooligomer)

  if len(homooligomer) == 0: homooligomer = "1"
  homooligomer = re.sub("[^0-9:]", "", homooligomer)

  # define inputs
  I = {"ori_sequence":sequence,
       "sequence":sequence.replace("/","").replace(":",""),
       "seqs":sequence.replace("/","").split(":"),
       "homooligomer":homooligomer,
       "homooligomers":[int(h) for h in homooligomer.split(":")],
       "msas":[], "deletion_matrices":[]}

  # adjust homooligomer option
  if len(I["seqs"]) != len(I["homooligomers"]):
    if len(I["homooligomers"]) == 1:
      I["homooligomers"] = [I["homooligomers"][0]] * len(I["seqs"])
    else:
      if verbose:
        print("WARNING: Mismatch between number of breaks ':' in 'sequence' and 'homooligomer' definition")
      while len(I["seqs"]) > len(I["homooligomers"]):
        I["homooligomers"].append(1)
      I["homooligomers"] = I["homooligomers"][:len(I["seqs"])]
    I["homooligomer"] = ":".join([str(h) for h in I["homooligomers"]])

  # define full sequence being modelled
  I["full_sequence"] = ''.join([s*h for s,h in zip(I["seqs"],I["homooligomers"])])
  I["lengths"] = [len(seq) for seq in I["seqs"]]

  # prediction directory
  if output_dir is None:
    I["output_dir"] = 'prediction_' + jobname + '_' + cf.get_hash(I["full_sequence"])[:5]
  else:
    I["output_dir"] = output_dir
  os.makedirs(I["output_dir"], exist_ok=True)
    
  # delete existing files in working directory
  if clean:
    for f in os.listdir(I["output_dir"]):
      os.remove(os.path.join(I["output_dir"], f))

  if verbose and len(I["full_sequence"]) > 1400:
    print(f"WARNING: For a typical Google-Colab-GPU (16G) session, the max total length is ~1400 residues. You are at {len(I['full_sequence'])}!")
    print(f"Run Alphafold may crash, unless you trim to the protein(s) to a short length. (See trim options below).")

  if verbose:
    print(f"homooligomer: {I['homooligomer']}")
    print(f"total_length: {len(I['full_sequence'])}")
    print(f"output_dir: {I['output_dir']}")
    
  return I

#######################################################################################################################################
# prep_msa
#######################################################################################################################################

def run_jackhmmer(sequence, prefix, jackhmmer_binary_path='jackhmmer', verbose=True):

  fasta_path = f"{prefix}.fasta"
  with open(fasta_path, 'wt') as f:
    f.write(f'>query\n{sequence}')

  pickled_msa_path = f"{prefix}.jackhmmer.pickle"
  if os.path.isfile(pickled_msa_path):
    msas_dict = pickle.load(open(pickled_msa_path,"rb"))
    msas, deletion_matrices, names = (msas_dict[k] for k in ['msas', 'deletion_matrices', 'names'])
    full_msa = []
    for msa in msas:
      full_msa += msa
  else:
    # --- Find the closest source ---
    test_url_pattern = 'https://storage.googleapis.com/alphafold-colab{:s}/latest/uniref90_2021_03.fasta.1'
    ex = futures.ThreadPoolExecutor(3)
    def fetch(source):
      request.urlretrieve(test_url_pattern.format(source))
      return source
    fs = [ex.submit(fetch, source) for source in ['', '-europe', '-asia']]
    source = None
    for f in futures.as_completed(fs):
      source = f.result()
      ex.shutdown()
      break
      
    dbs = []

    num_jackhmmer_chunks = {'uniref90': 59, 'smallbfd': 17, 'mgnify': 71}
    total_jackhmmer_chunks = sum(num_jackhmmer_chunks.values())
    disable_tqdm = not verbose
    with tqdm.notebook.tqdm(total=total_jackhmmer_chunks, bar_format=TQDM_BAR_FORMAT, disable=disable_tqdm) as pbar:
      def jackhmmer_chunk_callback(i):
        pbar.update(n=1)

      pbar.set_description('Searching uniref90')
      jackhmmer_uniref90_runner = jackhmmer.Jackhmmer(
          binary_path=jackhmmer_binary_path,
          database_path=f'https://storage.googleapis.com/alphafold-colab{source}/latest/uniref90_2021_03.fasta',
          get_tblout=True,
          num_streamed_chunks=num_jackhmmer_chunks['uniref90'],
          streaming_callback=jackhmmer_chunk_callback,
          z_value=135301051)
      dbs.append(('uniref90', jackhmmer_uniref90_runner.query(fasta_path)))

      pbar.set_description('Searching smallbfd')
      jackhmmer_smallbfd_runner = jackhmmer.Jackhmmer(
          binary_path=jackhmmer_binary_path,
          database_path=f'https://storage.googleapis.com/alphafold-colab{source}/latest/bfd-first_non_consensus_sequences.fasta',
          get_tblout=True,
          num_streamed_chunks=num_jackhmmer_chunks['smallbfd'],
          streaming_callback=jackhmmer_chunk_callback,
          z_value=65984053)
      dbs.append(('smallbfd', jackhmmer_smallbfd_runner.query(fasta_path)))

      pbar.set_description('Searching mgnify')
      jackhmmer_mgnify_runner = jackhmmer.Jackhmmer(
          binary_path=jackhmmer_binary_path,
          database_path=f'https://storage.googleapis.com/alphafold-colab{source}/latest/mgy_clusters_2019_05.fasta',
          get_tblout=True,
          num_streamed_chunks=num_jackhmmer_chunks['mgnify'],
          streaming_callback=jackhmmer_chunk_callback,
          z_value=304820129)
      dbs.append(('mgnify', jackhmmer_mgnify_runner.query(fasta_path)))

    # --- Extract the MSAs and visualize ---
    # Extract the MSAs from the Stockholm files.
    # NB: deduplication happens later in pipeline.make_msa_features.

    mgnify_max_hits = 501
    msas = []
    deletion_matrices = []
    names = []
    for db_name, db_results in dbs:
      unsorted_results = []
      for i, result in enumerate(db_results):
        msa, deletion_matrix, target_names = parsers.parse_stockholm(result['sto'])
        e_values_dict = parsers.parse_e_values_from_tblout(result['tbl'])
        e_values = [e_values_dict[t.split('/')[0]] for t in target_names]
        zipped_results = zip(msa, deletion_matrix, target_names, e_values)
        if i != 0:
          # Only take query from the first chunk
          zipped_results = [x for x in zipped_results if x[2] != 'query']
        unsorted_results.extend(zipped_results)
      sorted_by_evalue = sorted(unsorted_results, key=lambda x: x[3])
      db_msas, db_deletion_matrices, db_names, _ = zip(*sorted_by_evalue)
      if db_msas:
        if db_name == 'mgnify':
          db_msas = db_msas[:mgnify_max_hits]
          db_deletion_matrices = db_deletion_matrices[:mgnify_max_hits]
          db_names = db_names[:mgnify_max_hits]
        msas.append(db_msas)
        deletion_matrices.append(db_deletion_matrices)
        names.append(db_names)
        msa_size = len(set(db_msas))
        print(f'{msa_size} Sequences Found in {db_name}')

      pickle.dump({"msas":msas,
                   "deletion_matrices":deletion_matrices,
                   "names":names}, open(pickled_msa_path,"wb"))
  return msas, deletion_matrices, names

def prep_msa(I, msa_method="mmseqs2", add_custom_msa=False, msa_format="fas",
             pair_mode="unpaired", pair_cov=50, pair_qid=20,
             hhfilter_loc="hhfilter", reformat_loc="reformat.pl", TMP_DIR="tmp",
             custom_msa=None, precomputed=None,
             mmseqs_host_url="https://a3m.mmseqs.com",
             verbose=True):
  
  # make temp directory
  os.makedirs(TMP_DIR, exist_ok=True)

  # clear previous inputs
  I["msas"] = []
  I["deletion_matrices"] = []

  if add_custom_msa:
    if IN_COLAB:
      print(f"upload custom msa in '{msa_format}' format")
      msa_dict = files.upload()      
      lines = msa_dict[list(msa_dict.keys())[0]].decode()
      input_file = os.path.join(I["output_dir"],f"upload.{msa_format}")
      with open(input_file,"w") as tmp_upload:
        tmp_upload.write(lines)
    else:
      input_file = custom_msa
    if input_file is None or not os.path.isfile(input_file):
      raise ValueError("ERROR: `custom_msa` undefined")      
    else:
      # convert to a3m
      output_file = os.path.join(I["output_dir"],f"upload.a3m")
      os.system(f"{reformat_loc} {msa_format} a3m {input_file} {output_file}")

      # parse
      msa, mtx = parsers.parse_a3m(open(output_file,"r").read())
      I["msas"].append(msa)
      I["deletion_matrices"].append(mtx)
      if len(I["msas"][0][0]) != len(I["sequence"]):
        raise ValueError("ERROR: the length of msa does not match input sequence")

  if msa_method == "precomputed":
    if IN_COLAB:
      print("upload precomputed pickled msa from previous run")
      uploaded_dict = files.upload()
      uploaded_filename = list(uploaded_dict.keys())[0]
      I.update(pickle.loads(uploaded_dict[uploaded_filename]))
    elif precomputed is None:
      raise ValueError("ERROR: `precomputed` undefined")
    else:
      I.update(pickle.load(open(precomputed,"rb")))

  elif msa_method == "single_sequence":
    if len(I["msas"]) == 0:
      I["msas"].append([I["sequence"]])
      I["deletion_matrices"].append([[0]*len(I["sequence"])])

  else:
    _blank_seq = ["-" * L for L in I["lengths"]]
    _blank_mtx = [[0] * L for L in I["lengths"]]
    def _pad(ns,vals,mode):
      if mode == "seq": _blank = _blank_seq.copy()
      if mode == "mtx": _blank = _blank_mtx.copy()
      if isinstance(ns, list):
        for n,val in zip(ns,vals): _blank[n] = val
      else: _blank[ns] = vals
      if mode == "seq": return "".join(_blank)
      if mode == "mtx": return sum(_blank,[])

    if len(I["seqs"]) == 1 or "unpaired" in pair_mode:
      # gather msas
      if msa_method == "mmseqs2":
        prefix = cf.get_hash(I["sequence"])
        prefix = os.path.join(TMP_DIR,prefix)
        print(f"running mmseqs2")
        A3M_LINES = cf.run_mmseqs2(I["seqs"], prefix, use_filter=True, host_url=mmseqs_host_url)

      for n, seq in enumerate(I["seqs"]):
        # tmp directory
        prefix = cf.get_hash(seq)
        prefix = os.path.join(TMP_DIR,prefix)

        if msa_method == "mmseqs2":
          # run mmseqs2
          a3m_lines = A3M_LINES[n]
          msa, mtx = parsers.parse_a3m(a3m_lines)
          msas_, mtxs_ = [msa],[mtx]

        elif msa_method == "jackhmmer":
          print(f"running jackhmmer on seq_{n}")
          # run jackhmmer
          msas_, mtxs_, names_ = ([sum(x,())] for x in run_jackhmmer(seq, prefix))
        
        # pad sequences
        for msa_,mtx_ in zip(msas_,mtxs_):
          msa,mtx = [I["sequence"]],[[0]*len(I["sequence"])]      
          for s,m in zip(msa_,mtx_):
            msa.append(_pad(n,s,"seq"))
            mtx.append(_pad(n,m,"mtx"))

          I["msas"].append(msa)
          I["deletion_matrices"].append(mtx)

    # PAIR_MSA
    if len(I["seqs"]) > 1 and (pair_mode == "paired" or pair_mode == "unpaired+paired"):
      print("attempting to pair some sequences...")

      if msa_method == "mmseqs2":
        prefix = cf.get_hash(I["sequence"])
        prefix = os.path.join(TMP_DIR,prefix)
        print(f"running mmseqs2_noenv_nofilter on all seqs")
        A3M_LINES = cf.run_mmseqs2(I["seqs"], prefix, use_env=False, use_filter=False, host_url=mmseqs_host_url)

      _data = []
      for a in range(len(I["seqs"])):
        print(f"prepping seq_{a}")
        _seq = I["seqs"][a]
        _prefix = os.path.join(TMP_DIR,cf.get_hash(_seq))

        if msa_method == "mmseqs2":
          a3m_lines = A3M_LINES[a]
          _msa, _mtx, _lab = pairmsa.parse_a3m(a3m_lines,
                                               filter_qid=pair_qid/100,
                                               filter_cov=pair_cov/100)

        elif msa_method == "jackhmmer":
          _msas, _mtxs, _names = run_jackhmmer(_seq, _prefix)
          _msa, _mtx, _lab = pairmsa.get_uni_jackhmmer(_msas[0], _mtxs[0], _names[0],
                                                       filter_qid=pair_qid/100,
                                                       filter_cov=pair_cov/100)

        if len(_msa) > 1:
          _data.append(pairmsa.hash_it(_msa, _lab, _mtx, call_uniprot=False))
        else:
          _data.append(None)
      
      Ln = len(I["seqs"])
      O = [[None for _ in I["seqs"]] for _ in I["seqs"]]
      for a in range(Ln):
        if _data[a] is not None:
          for b in range(a+1,Ln):
            if _data[b] is not None:
              print(f"attempting pairwise stitch for {a} {b}")            
              O[a][b] = pairmsa._stitch(_data[a],_data[b])
              _seq_a, _seq_b, _mtx_a, _mtx_b = (*O[a][b]["seq"],*O[a][b]["mtx"])

              # filter to remove redundant sequences
              ok = []
              with open(f"{TMP_DIR}/tmp.fas","w") as fas_file:
                fas_file.writelines([f">{n}\n{a+b}\n" for n,(a,b) in enumerate(zip(_seq_a,_seq_b))])
              os.system(f"{hhfilter_loc} -maxseq 1000000 -i {TMP_DIR}/tmp.fas -o {TMP_DIR}/tmp.id90.fas -id 90")
              for line in open(f"{TMP_DIR}/tmp.id90.fas","r"):
                if line.startswith(">"): ok.append(int(line[1:]))
                
              if verbose:      
                print(f"found {len(_seq_a)} pairs ({len(ok)} after filtering)")

              if len(_seq_a) > 0:
                msa,mtx = [I["sequence"]],[[0]*len(I["sequence"])]
                for s_a,s_b,m_a,m_b in zip(_seq_a, _seq_b, _mtx_a, _mtx_b):
                  msa.append(_pad([a,b],[s_a,s_b],"seq"))
                  mtx.append(_pad([a,b],[m_a,m_b],"mtx"))
                I["msas"].append(msa)
                I["deletion_matrices"].append(mtx)
      
  # save MSA as pickle
  pickle.dump({"msas":I["msas"],"deletion_matrices":I["deletion_matrices"]},
              open(os.path.join(I["output_dir"],"msa.pickle"),"wb"))
  return I

#######################################################################################################################################
# prep_filter
#######################################################################################################################################

def trim_inputs(trim, msas, deletion_matrices, ori_seq=None, inverse=False):
  '''
  input: trim, msas, deletion_matrices, ori_seq
  output: msas, deletion_matrices, ori_seq
  '''
  if ori_seq is None: ori_seq = msas[0][0]
  seqs = ori_seq.replace("/","").split(":")
  L_ini = 0
  chain_idx = {}
  idx_chain = []
  for chain,seq in zip(ascii_uppercase,seqs):
    L = len(seq)
    chain_idx[chain] = dict(zip(range(L),range(L_ini,L_ini+L)))
    idx_chain += [f"{chain}{i+1}" for i in range(L)]
    L_ini += L
  global_idx = dict(zip(range(L_ini),range(L_ini)))

  mode = "keeping" if inverse else "trimming"
  trim_set = []
  for idx in trim.split(","):

    i,j = idx.split("-") if "-" in idx else (idx,"")

    # set index reference frame
    trim_idx_i = trim_idx_j = global_idx    
    if i != "" and i[0] in ascii_uppercase:
      trim_idx_i,i = chain_idx[i[0]], i[1:]
    if j != "" and j[0] in ascii_uppercase:
      trim_idx_j,j = chain_idx[j[0]], j[1:]

    # set which positions to trim
    if "-" in idx:
      i = trim_idx_i[int(i)-1] if i != "" else trim_idx_i[0]
      j = trim_idx_j[int(j)-1] if j != "" else trim_idx_j[len(trim_idx_j) - 1]
      trim_set += list(range(i,j+1))
      print(f"{mode} positions: {idx_chain[i]}-{idx_chain[j]}")
    else:
      i = trim_idx_i[int(i)-1]
      trim_set.append(i)
      print(f"{mode} position: {idx_chain[i]}")

  # deduplicate list
  trim_set = set(trim_set)
  if inverse:
    trim_set = set(range(L_ini)) ^ trim_set

  trim_set = sorted(list(trim_set))

  # trim MSA
  mod_msas, mod_mtxs = [],[]
  for msa, mtx in zip(msas, deletion_matrices):
    mod_msa = np.delete([list(s) for s in msa], trim_set, 1)
    ok = (mod_msa != "-").sum(-1) > 0
    mod_msas.append(["".join(s) for s in mod_msa[ok]])
    mod_mtx = np.asarray(mtx)[ok]
    mod_mtxs.append(np.delete(mod_mtx, trim_set, 1).tolist())

  # trim original sequence
  mod_idx = []
  mod_chain = []
  mod_ori_seq = []
  for n,a in enumerate(ori_seq.replace("/","").replace(":","")):
    if n not in trim_set:
      mod_ori_seq.append(a)
      mod_idx.append(n)
      mod_chain.append(idx_chain[n][0])
      if len(mod_idx) > 1:
        if mod_chain[-1] != mod_chain[-2]:
          mod_ori_seq[-1] = ":"
          mod_ori_seq.append(a)
        elif (mod_idx[-1] - mod_idx[-2]) > 1:
          mod_ori_seq[-1] = "/"
          mod_ori_seq.append(a)

  mod_ori_seq = "".join(mod_ori_seq)
  chains = sorted([ascii_uppercase.index(a) for a in set(mod_chain)])
  return {"msas":mod_msas, "deletion_matrices":mod_mtxs,
          "ori_sequence":mod_ori_seq, "chains":chains}

def cov_qid_filter(msas, deletion_matrices, ori_seq=None, cov=0, qid=0):
  if ori_seq is None: ori_seq = msas[0][0]
  seqs = ori_seq.replace("/","").split(":")
  ref_seq_ = np.array(list("".join(seqs)))

  new_msas,new_mtxs = [],[]
  L = np.asarray([len(seq) for seq in seqs])
  Ln = np.cumsum(np.append(0,L))
  for msa, mtx in zip(msas, deletion_matrices):
    msa_ = np.asarray([list(seq) for seq in msa])

    # coverage (non-gap characters)
    cov_ = msa_ != "-"
    # sequence identity to query
    qid_ = msa_ == ref_seq_

    # split by protein (for protein complexes)
    cov__ = np.stack([cov_[:,Ln[i]:Ln[i+1]].sum(-1) for i in range(len(seqs))],-1)
    qid__ = np.stack([qid_[:,Ln[i]:Ln[i+1]].sum(-1) for i in range(len(seqs))],-1)

    not_empty__ = cov__ > 0
    ok = []
    for n in range(len(msa)):
      m = not_empty__[n]
      if m.sum() > 0:
        q = qid__[n][m].sum() / cov__[n][m].sum()
        c = cov__[n][m].sum() / L[m].sum()
        if q > qid and c > cov:
         ok.append(n)

    new_msas.append([msa[n] for n in ok])
    new_mtxs.append([mtx[n] for n in ok])
  return {"msas":new_msas, "deletion_matrices":new_mtxs}

def prep_filter(I, trim="", trim_inverse=False, cov=0, qid=0, verbose=True):
  trim = re.sub("[^0-9A-Z,-]", "", trim.upper())
  trim = re.sub(",+",",",trim)
  trim = re.sub("^[,]+","",trim)
  trim = re.sub("[,]+$","",trim)
  if trim != "" or cov > 0 or qid > 0:
    mod_I = dict(I)
    
    if trim != "":
      mod_I.update(trim_inputs(trim, mod_I["msas"], mod_I["deletion_matrices"],
                               mod_I["ori_sequence"], inverse=trim_inverse))
      
      mod_I["homooligomers"] = [mod_I["homooligomers"][c] for c in mod_I["chains"]]
      mod_I["sequence"] = mod_I["ori_sequence"].replace("/","").replace(":","")
      mod_I["seqs"] = mod_I["ori_sequence"].replace("/","").split(":") 
      mod_I["full_sequence"] = "".join([s*h for s,h in zip(mod_I["seqs"], mod_I["homooligomers"])])
      new_length = len(mod_I["full_sequence"])
      if verbose:
        print(f"total_length: '{new_length}' after trimming")

    if cov > 0 or qid > 0:
      mod_I.update(cov_qid_filter(mod_I["msas"], mod_I["deletion_matrices"],
                                     mod_I["ori_sequence"], cov=cov/100, qid=qid/100))
    return mod_I
  else:
    return I  
  
#######################################################################################################################################
# prep features
#######################################################################################################################################

def prep_feats(I, clean=False):
  def _placeholder_template_feats(num_templates_, num_res_):
    return {
        'template_aatype': np.zeros([num_templates_, num_res_, 22], np.float32),
        'template_all_atom_masks': np.zeros([num_templates_, num_res_, 37, 3], np.float32),
        'template_all_atom_positions': np.zeros([num_templates_, num_res_, 37], np.float32),
        'template_domain_names': np.zeros([num_templates_], np.float32),
        'template_sum_probs': np.zeros([num_templates_], np.float32),
    }
  # delete old files
  if clean:
    for f in os.listdir(I["output_dir"]):
      if "rank_" in f: os.remove(os.path.join(I["output_dir"], f))
        
  if len(I["msas"]) == 0:
    print("WARNING: no MSA found, switching to 'single_sequence' mode")
    I["msas"].append([I["sequence"]])
    I["deletion_matrices"].append([[0]*len(I["sequence"])])

  # homooligomerize
  lengths = [len(seq) for seq in I["seqs"]]
  msas_mod, deletion_matrices_mod = cf.homooligomerize_heterooligomer(I["msas"], I["deletion_matrices"],
                                                                      lengths, I["homooligomers"])
  # define input features
  num_res = len(I["full_sequence"])
  feature_dict = {}
  feature_dict.update(pipeline.make_sequence_features(I["full_sequence"], 'test', num_res))
  feature_dict.update(pipeline.make_msa_features(msas_mod, deletion_matrices=deletion_matrices_mod))
  feature_dict.update(_placeholder_template_feats(0, num_res))

  # set chainbreaks
  Ls = []
  for seq,h in zip(I["ori_sequence"].split(":"), I["homooligomers"]):
    Ls += [len(s) for s in seq.split("/")] * h
  Ls_plot = []
  for seq,h in zip(I["seqs"], I["homooligomers"]):
    Ls_plot += [len(seq)] * h

  feature_dict['residue_index'] = cf.chain_break(feature_dict['residue_index'], Ls)
  feature_dict['Ls'] = Ls_plot
  feature_dict['output_dir'] = I["output_dir"]
  return feature_dict

def make_fixed_size(feat, runner):
  '''pad input features'''
  opt = runner["opt"]
  cfg = runner["model"].config
  shape_schema = {k:[None]+v for k,v in dict(cfg.data.eval.feat).items()}   
  pad_size_map = {
      shape_placeholders.NUM_RES: opt["L"],
      shape_placeholders.NUM_MSA_SEQ: cfg.data.eval.max_msa_clusters,
      shape_placeholders.NUM_EXTRA_SEQ: cfg.data.common.max_extra_msa,
      shape_placeholders.NUM_TEMPLATES: 0,
  }
  for k, v in feat.items():
    # Don't transfer this to the accelerator.
    if k == 'extra_cluster_assignment':
      continue
    shape = list(v.shape)
    schema = shape_schema[k]
    assert len(shape) == len(schema), (
        f'Rank mismatch between shape and shape schema for {k}: '
        f'{shape} vs {schema}')
    pad_size = [pad_size_map.get(s2, None) or s1 for (s1, s2) in zip(shape, schema)]
    padding = [(0, p - tf.shape(v)[i]) for i, p in enumerate(pad_size)]
    if padding:
      feat[k] = tf.pad(v, padding, name=f'pad_to_fixed_{k}')
      feat[k].set_shape(pad_size)
  return {k:np.asarray(v) for k,v in feat.items()}

#######################################################################################################################################
# run alphafold
#######################################################################################################################################

def clear_mem(device=None):
  '''remove all data from device'''
  backend = jax.lib.xla_bridge.get_backend(device)
  if hasattr(backend,'live_buffers'):
    for buf in backend.live_buffers():
      buf.delete()

OPT_DEFAULT = {"N":None, "L":None,
               "use_ptm":True, "use_turbo":True,
               "max_recycles":3, "tol":0, "num_ensemble":1,
               "max_msa_clusters":512, "max_extra_msa":1024,
               "is_training":False}

def prep_model_runner(opt=None, model_name="model_5", old_runner=None, params_loc='./alphafold/data'):
  
  # setup the [opt]ions
  if opt is None:
    opt = OPT_DEFAULT.copy()
  else:
    for k in OPT_DEFAULT:
      if k not in opt: opt[k] = OPT_DEFAULT[k]

  # if old_runner not defined or [opt]ions changed, start new runner
  if old_runner is None or old_runner["opt"] != opt:
    clear_mem()
    name = f"{model_name}_ptm" if opt["use_ptm"] else model_name
    cfg = config.model_config(name)

    if opt["use_turbo"]:
      if opt["N"] is None:
        cfg.data.eval.max_msa_clusters = opt["max_msa_clusters"]
        cfg.data.common.max_extra_msa = opt["max_extra_msa"]
      else:
        msa_clusters = min(opt["N"], opt["max_msa_clusters"])
        cfg.data.eval.max_msa_clusters = msa_clusters
        cfg.data.common.max_extra_msa = max(min(opt["N"] - msa_clusters, opt["max_extra_msa"]),1)

    cfg.data.common.num_recycle = opt["max_recycles"]
    cfg.model.num_recycle = opt["max_recycles"]
    cfg.model.recycle_tol = opt["tol"]
    cfg.data.eval.num_ensemble = opt["num_ensemble"]

    params = data.get_model_haiku_params(name, params_loc)
    return {"model":model.RunModel(cfg, params, is_training=opt["is_training"]), "opt":opt}
  else:
    return old_runner
  
def run_alphafold(feature_dict, opt=None, runner=None, num_models=5, num_samples=1, subsample_msa=True,
                  pad_feats=False, rank_by="pLDDT", show_images=True, params_loc='./alphafold/data', verbose=True):
  
  def do_subsample_msa(F, random_seed=0):
    '''subsample msa to avoid running out of memory'''
    N = len(F["msa"])
    L = len(F["residue_index"])
    N_ = int(3E7/L)
    if N > N_:
      if verbose:
        print(f"whhhaaa... too many sequences ({N}) subsampling to {N_}")
      np.random.seed(random_seed)
      idx = np.append(0,np.random.permutation(np.arange(1,N)))[:N_]
      F_ = {}
      F_["msa"] = F["msa"][idx]
      F_["deletion_matrix_int"] = F["deletion_matrix_int"][idx]
      F_["num_alignments"] = np.full_like(F["num_alignments"],N_)
      for k in F.keys():
        if k not in F_: F_[k] = F[k]      
      return F_
    else:
      return F

  def parse_results(prediction_result, processed_feature_dict, r, t, num_res):
    '''parse results and convert to numpy arrays'''
    
    to_np = lambda a: np.asarray(a)
    def class_to_np(c):
      class dict2obj():
        def __init__(self, d):
          for k,v in d.items(): setattr(self, k, to_np(v))
      return dict2obj(c.__dict__)

    dist_bins = jax.numpy.append(0,prediction_result["distogram"]["bin_edges"])
    dist_logits = prediction_result["distogram"]["logits"][:num_res,:][:,:num_res]
    dist_mtx = dist_bins[dist_logits.argmax(-1)]
    contact_mtx = jax.nn.softmax(dist_logits)[:,:,dist_bins < 8].sum(-1)

    b_factors = prediction_result['plddt'][:,None] * prediction_result['structure_module']['final_atom_mask']
    p = protein.from_prediction(processed_feature_dict, prediction_result, b_factors=b_factors)  
    plddt = prediction_result['plddt'][:num_res]
    out = {"unrelaxed_protein": class_to_np(p),
           "plddt": to_np(plddt),
           "pLDDT": to_np(plddt.mean()),
           "dists": to_np(dist_mtx),
           "adj": to_np(contact_mtx),
           "recycles":to_np(r),
           "tol":to_np(t)}
    if "ptm" in prediction_result:
      out["pae"] = to_np(prediction_result['predicted_aligned_error'][:num_res,:][:,:num_res])
      out["pTMscore"] = to_np(prediction_result['ptm'])      
    return out

  num_res = len(feature_dict["residue_index"])

  # if [opt]ions not defined
  if opt is None:
    opt = OPT_DEFAULT.copy()
    opt["N"] = len(feature_dict["msa"])
    opt["L"] = num_res
  else:
    for k in OPT_DEFAULT.keys():
      if k not in opt: opt[k] = OPT_DEFAULT[k]
          
  model_names = ['model_1', 'model_2', 'model_3', 'model_4', 'model_5'][:num_models]
  total = len(model_names) * num_samples
  outs = {}

  def do_report(key):
    o = outs[key]
    if verbose:
      line = f"{key} recycles:{o['recycles']} tol:{o['tol']:.2f} pLDDT:{o['pLDDT']:.2f}"
      if 'pTMscore' in o:
        line += f" pTMscore:{o['pTMscore']:.2f}"
      print(line)
    if show_images:
      fig = cf.plot_protein(o['unrelaxed_protein'], Ls=feature_dict["Ls"], dpi=100)
      plt.show()  
    tmp_pdb_path = os.path.join(feature_dict["output_dir"],f'unranked_{key}_unrelaxed.pdb')
    pdb_lines = protein.to_pdb(o['unrelaxed_protein'])
    with open(tmp_pdb_path, 'w') as f: f.write(pdb_lines)
  
  disable_tqdm = not verbose
  with tqdm.notebook.tqdm(total=total, bar_format=TQDM_BAR_FORMAT, disable=disable_tqdm) as pbar:
    if opt["use_turbo"]:
      if runner is None:
        runner = prep_model_runner(opt,params_loc=params_loc)
      
      # go through each random_seed
      for seed in range(num_samples):
        # prep input features
        feat = do_subsample_msa(feature_dict, random_seed=seed) if subsample_msa else feature_dict
        processed_feature_dict = runner["model"].process_features(feat, random_seed=seed)
        if pad_feats:
          processed_feature_dict = make_fixed_size(processed_feature_dict, runner)

        # go through each model
        for num, model_name in enumerate(model_names):
          name = model_name+"_ptm" if opt["use_ptm"] else model_name
          key = f"{name}_seed_{seed}"
          pbar.set_description(f'Running {key}')

          # replace model parameters
          params = data.get_model_haiku_params(name, params_loc)
          for k in runner["model"].params.keys():
            runner["model"].params[k] = params[k]

          # predict
          prediction_result, (r, t) = runner["model"].predict(processed_feature_dict, random_seed=seed)
          outs[key] = parse_results(prediction_result, processed_feature_dict, r=r, t=t, num_res=num_res)            

          # cleanup
          del prediction_result, params, r, t

          # report
          do_report(key)
          pbar.update(n=1)

        # cleanup
        del processed_feature_dict
        if subsample_msa: del feat

    else:  
      # go through each model
      for num, model_name in enumerate(model_names):
        name = model_name+"_ptm" if opt["use_ptm"] else model_name
        model_runner = prep_model_runner(opt, model_name=model_name, use_turbo=False, params_loc=params_loc)["model"]

        # go through each random_seed
        for seed in range(num_samples):
          key = f"{name}_seed_{seed}"
          pbar.set_description(f'Running {key}')
          processed_feature_dict = model_runner.process_features(feature_dict, random_seed=seed)
          
          # predict
          prediction_result, (r, t) = model_runner.predict(processed_feature_dict, random_seed=seed)
          outs[key] = parse_results(prediction_result, processed_feature_dict, r=r, t=t, num_res=num_res)

          # cleanup
          del processed_feature_dict, prediction_result, r, t

          # report          
          do_report(key)
          pbar.update(n=1)
        
        # cleanup
        del model_runner
  
  # Find the best model according to the mean pLDDT.
  model_rank = list(outs.keys())
  model_rank = [model_rank[i] for i in np.argsort([outs[x][rank_by] for x in model_rank])[::-1]]

  # Write out the prediction
  for n,key in enumerate(model_rank):
    prefix = f"rank_{n+1}_{key}" 
    pred_output_path = os.path.join(feature_dict["output_dir"],f'{prefix}_unrelaxed.pdb')
    fig = cf.plot_protein(outs[key]["unrelaxed_protein"], Ls=feature_dict["Ls"], dpi=200)
    plt.savefig(os.path.join(feature_dict["output_dir"],f'{prefix}.png'), bbox_inches = 'tight')
    plt.close(fig)
    pdb_lines = protein.to_pdb(outs[key]["unrelaxed_protein"])
    with open(pred_output_path, 'w') as f:
      f.write(pdb_lines)
    
    tmp_pdb_path = os.path.join(feature_dict["output_dir"],f'unranked_{key}_unrelaxed.pdb')
    if os.path.isfile(tmp_pdb_path):
      os.remove(tmp_pdb_path)

  ############################################################
  if verbose:
    print(f"model rank based on {rank_by}")
    for n,key in enumerate(model_rank):
      print(f"rank_{n+1}_{key} {rank_by}:{outs[key][rank_by]:.2f}")

  return outs, model_rank
