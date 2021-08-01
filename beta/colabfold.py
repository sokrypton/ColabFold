############################################
# imports
############################################
import jax
import requests
import hashlib
import tarfile
import time
import pickle
import os
import numpy as np

###########################################
# control gpu/cpu memory usage
###########################################
def rm(x):
  '''remove data from device'''
  jax.tree_util.tree_map(lambda y: y.device_buffer.delete(), x)

def to(x,device="cpu"):
  '''move data to device'''
  d = jax.devices(device)[0]
  return jax.tree_util.tree_map(lambda y:jax.device_put(y,d), x)

def clear_mem(device="gpu"):
  '''remove all data from device'''
  backend = jax.lib.xla_bridge.get_backend(device)
  for buf in backend.live_buffers(): buf.delete()
    
##########################################
# call mmseqs2
##########################################

def run_mmseqs2(query_sequence, prefix, use_env=True, filter=False):
    def submit(query_sequence, mode):
      res = requests.post('https://a3m.mmseqs.com/ticket/msa', data={'q':f">1\n{query_sequence}", 'mode': mode})
      return res.json()
    def status(ID):
      res = requests.get(f'https://a3m.mmseqs.com/ticket/{ID}')
      return res.json()
    def download(ID, path):
      res = requests.get(f'https://a3m.mmseqs.com/result/download/{ID}')
      with open(path,"wb") as out: out.write(res.content)
      
    if filter:
      mode = "env" if use_env else "all"
    else:
      mode = "env-nofilter" if use_env else "nofilter"
    
    path = f"{prefix}_{mode}"
    if not os.path.isdir(path): os.mkdir(path)

    # call mmseqs2 api
    tar_gz_file = f'{path}/out.tar.gz'
    if not os.path.isfile(tar_gz_file):
      out = submit(query_sequence, mode)
      while out["status"] in ["RUNNING","PENDING"]:
        time.sleep(1)
        out = status(out["id"])    
      download(out["id"], tar_gz_file)
    
    # parse a3m files
    a3m_lines = []
    a3m = f"{prefix}_{mode}.a3m"
    if not os.path.isfile(a3m):
      with tarfile.open(tar_gz_file) as tar_gz: tar_gz.extractall(path)
      a3m_files = [f"{path}/uniref.a3m"]
      if use_env: a3m_files.append(f"{path}/bfd.mgnify30.metaeuk30.smag30.a3m")
      a3m_out = open(a3m,"w")
      for a3m_file in a3m_files:
        for line in open(a3m_file,"r"):
          line = line.replace("\x00","")
          if len(line) > 0:
            a3m_lines.append(line)
            a3m_out.write(line)
    else:
      a3m_lines = open(a3m).readlines()
    return "".join(a3m_lines)

#########################################################################
# utils
#########################################################################
def get_hash(x):
  return hashlib.sha1(x.encode()).hexdigest()

def cov_filter(msas, deletion_matrices, cov=0):
  if cov > 0:
    filtered = 0
    new_msas = []
    new_mtxs = []
    for msa,mtx in zip(msas,deletion_matrices):
      new_msa = []
      new_mtx = []
      for s,m in zip(msa,mtx):
        c = (np.array(list(s)) != "-").mean()
        if c >= cov/100:
          new_msa.append(s)
          new_mtx.append(m)
        else:
          filtered += 1
      new_msas.append(new_msa)
      new_mtxs.append(new_mtx)
    print(f"Filtered {filtered} number of sequences that don't cover at least {cov}% of query")
    return new_msas, new_mtxs
  else:
    return msas, deletion_matrices

def homooliomerize(msas, deletion_matrices, homooligomer):
 if homooligomer == 1:
  return msas, deletion_matrices
 else:
  new_msas = []
  new_mtxs = []
  for o in range(homooligomer):
    L = num_res * o
    R = num_res * (homooligomer-(o+1))
    for msa,mtx in zip(msas,deletion_matrices):
      new_msas.append(["-"*L+seq+"-"*R for seq in msa])
      new_mtxs.append([[0]*L+mtx+[0]*R for mtx in deletion_matrix])
  return new_msas, new_mtxs

def chain_break(idx_res, Ls, length=200):
  # Minkyung's code
  # add big enough number to residue index to indicate chain breaks
  L_prev = 0
  for L_i in Ls[:-1]:
    idx_res[L_prev+L_i:] += lengths
    L_prev += L_i      
  return idx_res
