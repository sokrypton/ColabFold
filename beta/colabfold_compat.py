############################################
# imports
############################################
import random
import requests
import tarfile
import time
import os
import numpy as np
import tqdm.notebook
    
##########################################
# call mmseqs2
##########################################
TQDM_BAR_FORMAT = '{l_bar}{bar}| {n_fmt}/{total_fmt} [elapsed: {elapsed} remaining: {remaining}]'

def a3m_files_to_lines(a3m_files):
  # gather a3m lines  
  a3m_lines = {}
  for a3m_file in a3m_files:
    update_M,M = True,None
    for line in open(a3m_file,"r"):
      if len(line) > 0:
        if "\x00" in line:
          line = line.replace("\x00","")
          update_M = True
        if line.startswith(">") and update_M:
          M = int(line[1:].rstrip())
          update_M = False
          if M not in a3m_lines: a3m_lines[M] = []
        a3m_lines[M].append(line)
  return ["".join(lines) for lines in a3m_lines]

def run_mmseqs2(x, prefix, use_env=True, filter=True, return_a3m_lines=True):
  
  def submit(seqs, mode, N=101):
    
    n,query = N,""
    for seq in seqs:
      query += f">{n}\n{seq}\n"
      n += 1
      
    res = requests.post('https://a3m.mmseqs.com/ticket/msa', data={'q':query,'mode': mode})
    try: out = res.json()
    except ValueError: out = {"status":"UNKNOWN"}
    return out

  def status(ID):
    res = requests.get(f'https://a3m.mmseqs.com/ticket/{ID}')
    try: out = res.json()
    except ValueError: out = {"status":"UNKNOWN"}
    return out

  def download(ID, path):
    res = requests.get(f'https://a3m.mmseqs.com/result/download/{ID}')
    with open(path,"wb") as out: out.write(res.content)
  
  # process input x
  seqs = [x] if isinstance(x, str) else x
    
  # setup mode
  if filter:
    mode = "env" if use_env else "all"
  else:
    mode = "env-nofilter" if use_env else "nofilter"
  
  # define path
  path = f"{prefix}_{mode}"
  if not os.path.isdir(path): os.mkdir(path)

  # call mmseqs2 api
  tar_gz_file = f'{path}/out.tar.gz'
  if not os.path.isfile(tar_gz_file):
    TIME_ESTIMATE = 150 * len(seqs)
    with tqdm.notebook.tqdm(total=TIME_ESTIMATE, bar_format=TQDM_BAR_FORMAT) as pbar:
      N,REDO = 101,True
      while REDO:
        pbar.set_description("SUBMIT")
        
        # Resubmit job until it goes through
        out = submit(seqs, mode, N)
        while out["status"] in ["UNKNOWN","RATELIMIT"]:
          # resubmit
          time.sleep(5 + random.randint(0,5))
          out = submit(seqs, mode, N)

        # wait for job to finish
        ID,TIME = out["id"],0
        pbar.set_description(out["status"])
        while out["status"] in ["UNKNOWN","RUNNING","PENDING"]:
          t = 5 + random.randint(0,5)
          time.sleep(t)
          out = status(ID)    
          pbar.set_description(out["status"])
          if out["status"] == "RUNNING":
            TIME += t
            pbar.update(n=t)
          #if TIME > 900 and out["status"] != "COMPLETE":
          #  # something failed on the server side, need to resubmit
          #  N += 1
          #  break
        
        if out["status"] == "COMPLETE":
          if TIME < TIME_ESTIMATE:
            pbar.update(n=(TIME_ESTIMATE-TIME))
          REDO = False
          
      # Download results
      download(ID, tar_gz_file)


  # prep list of a3m files
  a3m_files = [f"{path}/uniref.a3m"]
  if use_env:
    a3m_files.append(f"{path}/bfd.mgnify30.metaeuk30.smag30.a3m")

  # extract a3m files
  if not os.path.isfile(a3m_files[0]):
    with tarfile.open(tar_gz_file) as tar_gz:
      tar_gz.extractall(path)  

  if return_a3m_lines:
    a3m_lines = a3m_files_to_lines(a3m_files)  
    
    # return results
    Ms = sorted(list(a3m_lines.keys()))
    if isinstance(x, str):
      return a3m_lines[Ms[0]]
    else:
      return [a3m_lines[n] for n in Ms]
  else:
    return a3m_files

#########################################################################
# run old mmseqs2 filtering schemes
#########################################################################
def mmseqs2_bucket_filter(A3M_LINES, prefix, hhfilter_bin="tmp/bin/hhfilter"):

  def hhfilter(a3m_in, a3m_out, hhf_id=90, hhf_qid=0, hhf_cov=0, hhf_diff=None):
    if not os.path.isfile(a3m_out):
      hhf_diff_ = 0 if hhf_diff is None else hhf_diff
      os.system(f"{hhfilter_bin} -i {a3m_in} -o {a3m_out} -qid {hhf_qid} -id {hhf_id} -cov {hhf_cov} -diff {hhf_diff_}")

  for m, a3m_lines in A3M_LINES.items():
    a3m_in = f"{prefix}.{m}.a3m"
    if not os.path.isfile(a3m_in):
      with open(a3m_in,"w") as a3m_file:
        a3m_file.write(a3m_lines)
            
    a3m_outs = []
    for hhf_qid in [0.5,0.3,0.15]:
      a3m_out = f"{prefix}.{m}.qid{hhf_qid}.a3m"
      a3m_outs.append(a3m_out)
      hhfilter(a3m_in, a3m_out, hhf_qid=hhf_qid)
    
    a3m_all_in = f"{prefix}.{m}.filt.a3m"
    a3m_all_out = f"{prefix}.{m}.filt.id90.a3m"
    if not os.path.isfile(a3m_all_in):
      a3m_outs = " ".join(a3m_outs)
      os.system(f"cat {a3m_outs} > {a3m_all_in}")
      
    hhfilter(a3m_all_in, a3m_all_out)
    A3M_LINES[m] = open(a3m_all_out,"r").read()

def run_mmseqs2_compat(x, prefix, hhfilter_bin="tmp/bin/hhfilter", filter_scheme="18Aug2021"):

  # 22Jul2021 diff_1000 filtering ONLY applied to uniref
  # 16Aug2021 diff_1000 filtering applied to both uniref and env
  # 18Aug2021 three seperate diff_1000 with qid_50, qid_30 and qid_15, merged and -id_90 filtered

  a3m_files = run_mmseqs2(x, prefix, use_env=True, filter=False, return_a3m_lines=False)

  if filter_scheme == "22Jul2021" or filter_scheme == "17Aug2021":
    os.system(f"{hhfilter_bin} -i {a3m_files[0]} -o {a3m_files[0]}.filt -diff 1000")
    a3m_files[0] += ".filt"

  elif filter_scheme == "17Aug2021":
    os.system(f"{hhfilter_bin} -i {a3m_files[1]} -o {a3m_files[1]}.filt -diff 1000")
    a3m_files[1] += ".filt"
  
  a3m_lines = a3m_files_to_lines(a3m_files)

  if filter_scheme == "18Aug2021":
    a3m_lines = mmseqs2_bucket_filter(a3m_lines, f"{prefix}_env-nofilter", hhfilter_bin=hhfilter_bin)

  # return results
  Ms = sorted(list(a3m_lines.keys()))
  if isinstance(x, str):
    return a3m_lines[Ms[0]]
  else:
    return [a3m_lines[n] for n in Ms]
