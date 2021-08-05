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
import matplotlib
from matplotlib import collections as mcoll

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
  plt.ylim(0,100)
  plt.ylabel('plDDT')
  plt.xlabel('position')
  if use_ptm:
    plt.subplot(1,2,2);plt.title('Predicted Aligned Error')
    plt.imshow(pae, cmap="bwr",vmin=0,vmax=30)
    plt.colorbar()
    plt.xlabel('Scored residue')
    plt.ylabel('Aligned residue')
  return plt

def show_pdb(pred_output_path, show_sidechains=False, show_mainchains=False,
             color="lDDT", chains=1, vmin=50, vmax=90):
  view = py3Dmol.view(js='https://3dmol.org/build/3Dmol.js',)
  view.addModel(open(pred_output_path,'r').read(),'pdb')
  if color == "lDDT":
    view.setStyle({'cartoon': {'colorscheme': {'prop':'b','gradient': 'roygb','min':vmin,'max':vmax}}})
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

def plot_dists(dists, dpi=100, fig=True):
  num_models = len(dists)
  if fig: plt.figure(figsize=(3*num_models,2), dpi=dpi)
  for n,dist in enumerate(dists):
    plt.subplot(1,num_models,n+1)
    plt.title(f"model_{n+1}")
    plt.imshow(dist)
    plt.colorbar()
  return plt

def plot_protein(protein=None, pos=None, plddt=None, dpi=100):
  
  if protein is not None:
    pos = protein.atom_positions[:,1,:]
    plddt = protein.b_factors[:,0]

  # get best view
  pos = pos - pos.mean(0,keepdims=True)
  eigen_val, eigen_vec = np.linalg.eigh(np.cov(pos.T))
  pos = (pos @ eigen_vec)[...,::-1]

  def make_segments(x, y):
    points = np.array([x, y]).T.reshape(-1, 1, 2)
    segments = np.concatenate([points[:-1], points[1:]], axis=1)
    return segments
  
  def get_color(x, alpha=None, tint=None, shade=None, vmin=50, vmax=90):
    cmap = matplotlib.cm.get_cmap('gist_rainbow')
    if x < vmin: x = vmin
    if x > vmax: x = vmax
    x = (x - vmin)/(vmax - vmin)
    color = np.array(cmap(x * 0.8, alpha=alpha))
    if tint is not None:
      color[:3] = color[:3] + (1 - color[:3]) * tint
    if shade is not None:
      color[:3] = color[:3] * shade
    return color

  if plddt is not None:
    fig, (ax1, ax3) = plt.subplots(1,2)
    fig.set_figwidth(6)
    fig.set_figheight(3)
  else:
    fig, ax1 = plt.subplots(1,1)
    fig.set_figwidth(3)
    fig.set_figheight(3)
    
  fig.set_dpi(dpi)
  fig.subplots_adjust(top = 0.9, bottom = 0.1, right = 1, left = 0, hspace = 0.1, wspace = 0.1)

  range_xy = (pos[:,:2].min(),pos[:,:2].max())
  range_z = (pos[:,-1].min(),pos[:,-1].max())
  z = (pos[:,-1] - range_z[0]) / (range_z[1] - range_z[0])

  srt = (z[:-1]+z[1:]).argsort()
  seg = make_segments(pos[:,0],pos[:,1])

  # color by rainbow
  colors1 = np.array([get_color(p, vmin=0, vmax=len(pos)) for p in np.arange(len(pos))[::-1]])
  ttl3 = plt.text(0.5, 1.01, f"colored by N->C", horizontalalignment='center', verticalalignment='bottom', transform=ax1.transAxes)
  ax1.axis('scaled')
  ax1.set_xlim(*range_xy)
  ax1.set_ylim(*range_xy)
  ax1.axis(False)
  ax1.add_collection(mcoll.LineCollection(seg[srt], colors=colors1[srt], linewidths=5))
   
  if plddt is not None:
    # color by plddt
    colors3 = np.array([get_color(p) for a,p in zip(z, plddt)])
    ttl3 = plt.text(0.5, 1.01, f"colored by pLDDT", horizontalalignment='center', verticalalignment='bottom', transform=ax3.transAxes)
    ax3.axis('scaled')
    ax3.set_xlim(*range_xy)
    ax3.set_ylim(*range_xy)
    ax3.axis(False)
    ax3.add_collection(mcoll.LineCollection(seg[srt], colors=colors3[srt], linewidths=5))
    plt.show()
