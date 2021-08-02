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
import matplotlib.pyplot as plt
import py3Dmol

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

def homooliomerize(msas, deletion_matrices, homooligomer=1):
 if homooligomer == 1:
  return msas, deletion_matrices
 else:
  new_msas = []
  new_mtxs = []
  for o in range(homooligomer):
    for msa,mtx in zip(msas, deletion_matrices):
      num_res = len(msa[0])
      L = num_res * o
      R = num_res * (homooligomer-(o+1))
      new_msas.append(["-"*L+s+"-"*R for s in msa])
      new_mtxs.append([[0]*L+m+[0]*R for m in mtx])
  return new_msas, new_mtxs

def chain_break(idx_res, Ls, length=200):
  # Minkyung's code
  # add big enough number to residue index to indicate chain breaks
  L_prev = 0
  for L_i in Ls[:-1]:
    idx_res[L_prev+L_i:] += length
    L_prev += L_i      
  return idx_res

##################################################
# parsers
##################################################
from alphafold.common import protein
def parse_results(prediction_result, processed_feature_dict):
  b_factors = prediction_result['plddt'][:,None] * prediction_result['structure_module']['final_atom_mask']
  out = {"unrelaxed_protein": protein.from_prediction(processed_feature_dict, prediction_result, b_factors=b_factors),
         "plddt": prediction_result['plddt'],
         "sco": prediction_result['plddt'].mean(),
         "mtx": prediction_result["distogram"]["bin_edges"][prediction_result["distogram"]["logits"].argmax(-1)],
         "adj": jax.nn.softmax(prediction_result["distogram"]["logits"])[:,:,prediction_result["distogram"]["bin_edges"] < 8].sum(-1)}
  if "ptm" in prediction_result:
    out.update({"pae": prediction_result['predicted_aligned_error'],
                "ptm": prediction_result['ptm']})
  return out

##################################################
# plotting
##################################################

def plot_plddt_legend(dpi=100):
  thresh = ['plDDT:','Very low (<50)','Low (60)','OK (70)','Confident (80)','Very high (>90)']
  plt.figure(figsize=(1,0.1),dpi=dpi)
  ########################################
  for c in ["#FFFFFF","#FF0000","#FFFF00","#00FF00","#00FFFF","#0000FF"]:
    plt.bar(0, 0, color=c)
  plt.legend(thresh, frameon=False,
             loc='center', ncol=6,
             handletextpad=1,
             columnspacing=1,
             markerscale=0.5,)
  plt.axis(False)
  return plt

def plot_confidence(plddt, pae=None, Ls=None, dpi=100):
  use_ptm = False if pae is None else True
  if use_ptm:
    plt.figure(figsize=(10,3), dpi=dpi)
    plt.subplot(1,2,1);
  else:
    plt.figure(figsize=(5,3), dpi=dpi)
  plt.title('Predicted lDDT')
  plt.plot(plddt)
  if Ls is not None:
    L_prev = 0
    for L_i in Ls[:-1]:
      L = L_prev + L_i
      L_prev += L_i
      plt.plot([L,L],[0,100],color="black")
  plt.ylabel('plDDT')
  plt.xlabel('position')
  if use_ptm:
    plt.subplot(1,2,2);plt.title('Predicted Aligned Error')
    plt.imshow(pae, cmap="bwr",vmin=0,vmax=30)
    plt.colorbar()
    plt.xlabel('Scored residue')
    plt.ylabel('Aligned residue')
  return plt

def show_pdb(pred_output_path, show_sidechains=False, show_mainchains=False, color="lDDT", chains=1):
  view = py3Dmol.view(js='https://3dmol.org/build/3Dmol.js',)
  view.addModel(open(pred_output_path,'r').read(),'pdb')
  if color == "lDDT":
    view.setStyle({'cartoon': {'colorscheme': {'prop':'b','gradient': 'roygb','min':50,'max':90}}})
  elif color == "rainbow":
    view.setStyle({'cartoon': {'color':'spectrum'}})
  elif color == "chain":
    for n,chain,color in zip(range(chains),list("ABCDEFGH"),
                     ["lime","cyan","magenta","yellow","salmon","white","blue","orange"]):
       view.setStyle({'chain':chain},{'cartoon': {'color':color}})
  if show_sidechains:
    BB = ['C','O','N']
    view.addStyle({'and':[{'resn':["GLY","PRO"],'invert':True},{'atom':BB,'invert':True}]},
                        {'stick':{'colorscheme':f"WhiteCarbon",'radius':0.3}})
    view.addStyle({'and':[{'resn':"GLY"},{'atom':'CA'}]},
                        {'sphere':{'colorscheme':f"WhiteCarbon",'radius':0.3}})
    view.addStyle({'and':[{'resn':"PRO"},{'atom':['C','O'],'invert':True}]},
                        {'stick':{'colorscheme':f"WhiteCarbon",'radius':0.3}})  
  if show_mainchains:
    BB = ['C','O','N','CA']
    view.addStyle({'atom':BB},{'stick':{'colorscheme':f"WhiteCarbon",'radius':0.3}})
  view.zoomTo()
  return view

def plot_plddts(plddts, Ls=None, dpi=100, fig=True):
  if fig: plt.figure(figsize=(8,5),dpi=100)
  plt.title("Predicted lDDT per position")
  for n,plddt in enumerate(plddts):
    plt.plot(plddt,label=f"model_{n+1}")
  if Ls is not None:
    L_prev = 0
    for L_i in Ls[:-1]:
      L = L_prev + L_i
      L_prev += L_i
      plt.plot([L,L],[0,100],color="black")
  plt.legend()
  plt.ylim(0,100)
  plt.ylabel("Predicted lDDT")
  plt.xlabel("Positions")
  return plt

def plot_paes(paes, dpi=100, fig=True):
  num_models = len(paes)
  if fig: plt.figure(figsize=(3*num_models,2), dpi=dpi)
  for n,pae in enumerate(paes):
    plt.subplot(1,num_models,n+1)
    plt.title(f"model_{n+1}")
    plt.imshow(pae,cmap="bwr",vmin=0,vmax=30)
    plt.colorbar()
  return plt

def plot_adjs(adjs, dpi=100, fig=True):
  num_models = len(adjs)
  if fig: plt.figure(figsize=(3*num_models,2), dpi=dpi)
  for n,adj in enumerate(adjs):
    plt.subplot(1,num_models,n+1)
    plt.title(f"model_{n+1}")
    plt.imshow(adj,cmap="binary",vmin=0,vmax=1)
    plt.colorbar()
  return plt
