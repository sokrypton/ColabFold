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
import matplotlib.patheffects
from matplotlib import collections as mcoll

import py3Dmol
from string import ascii_uppercase,ascii_lowercase

pymol_color_list = ["#33ff33","#00ffff","#ff33cc","#ffff00","#ff9999","#e5e5e5","#7f7fff","#ff7f00",
                    "#7fff7f","#199999","#ff007f","#ffdd5e","#8c3f99","#b2b2b2","#007fff","#c4b200",
                    "#8cb266","#00bfbf","#b27f7f","#fcd1a5","#ff7f7f","#ffbfdd","#7fffff","#ffff7f",
                    "#00ff7f","#337fcc","#d8337f","#bfff3f","#ff7fff","#d8d8ff","#3fffbf","#b78c4c",
                    "#339933","#66b2b2","#ba8c84","#84bf00","#b24c66","#7f7f7f","#3f3fa5","#a5512b"]

pymol_cmap = matplotlib.colors.ListedColormap(pymol_color_list)
alphabet_list = list(ascii_uppercase+ascii_lowercase)


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
    for n,chain,color in zip(range(chains),alphabet_list,pymol_color_list):
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
    plt.plot(plddt,label=f"rank_{n+1}")
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
    plt.title(f"rank_{n+1}")
    plt.imshow(pae,cmap="bwr",vmin=0,vmax=30)
    plt.colorbar()
  return plt

def plot_adjs(adjs, dpi=100, fig=True):
  num_models = len(adjs)
  if fig: plt.figure(figsize=(3*num_models,2), dpi=dpi)
  for n,adj in enumerate(adjs):
    plt.subplot(1,num_models,n+1)
    plt.title(f"rank_{n+1}")
    plt.imshow(adj,cmap="binary",vmin=0,vmax=1)
    plt.colorbar()
  return plt

def plot_dists(dists, dpi=100, fig=True):
  num_models = len(dists)
  if fig: plt.figure(figsize=(3*num_models,2), dpi=dpi)
  for n,dist in enumerate(dists):
    plt.subplot(1,num_models,n+1)
    plt.title(f"rank_{n+1}")
    plt.imshow(dist)
    plt.colorbar()
  return plt

##########################################################################
##########################################################################

def kabsch(a, b, weights=None, return_v=False):
  a = np.asarray(a)
  b = np.asarray(b)
  if weights is None: weights = np.ones(len(b))
  else: weights = np.asarray(weights)
  B = np.einsum('ji,jk->ik', weights[:, None] * a, b)
  u, s, vh = np.linalg.svd(B)
  if np.linalg.det(u @ vh) < 0: u[:, -1] = -u[:, -1]
  if return_v: return u
  else: return u @ vh

def plot_pseudo_3D(xyz, c=None, ax=None, chainbreak=5,
                   cmap="gist_rainbow", line_w=2.0,
                   cmin=None, cmax=None, zmin=None, zmax=None):

  def rescale(a,amin=None,amax=None):
    a = np.copy(a)
    if amin is None: amin = a.min()
    if amax is None: amax = a.max()
    a[a < amin] = amin
    a[a > amax] = amax
    return (a - amin)/(amax - amin)

  # make segments
  xyz = np.asarray(xyz)
  seg = np.concatenate([xyz[:-1,None,:],xyz[1:,None,:]],axis=-2)
  seg_xy = seg[...,:2]
  seg_z = seg[...,2].mean(-1)
  ord = seg_z.argsort()

  # set colors
  if c is None: c = np.arange(len(seg))[::-1]
  else: c = (c[1:] + c[:-1])/2
  c = rescale(c,cmin,cmax)  

  if isinstance(cmap, str):
    if cmap == "gist_rainbow": c *= 0.8
    colors = matplotlib.cm.get_cmap(cmap)(c)
  else:
    colors = cmap(c)
  
  if chainbreak is not None:
    dist = np.linalg.norm(xyz[:-1] - xyz[1:], axis=-1)
    colors[...,3] = (dist < chainbreak).astype(np.float)

  # add shade/tint based on z-dimension
  z = rescale(seg_z,zmin,zmax)[:,None]
  tint, shade = z/3, (z+2)/3
  colors[:,:3] = colors[:,:3] + (1 - colors[:,:3]) * tint
  colors[:,:3] = colors[:,:3] * shade

  set_lim = False
  if ax is None:
    fig, ax = plt.subplots()
    fig.set_figwidth(5)
    fig.set_figheight(5)
    set_lim = True
  else:
    fig = ax.get_figure()
    if ax.get_xlim() == (0,1):
      set_lim = True
      
  if set_lim:
    xy_min = xyz[:,:2].min() - line_w
    xy_max = xyz[:,:2].max() + line_w
    ax.set_xlim(xy_min,xy_max)
    ax.set_ylim(xy_min,xy_max)

  ax.set_aspect('equal')
    
  # determine linewidths
  width = fig.bbox_inches.width * ax.get_position().width
  linewidths = line_w * 72 * width / np.diff(ax.get_xlim())

  lines = mcoll.LineCollection(seg_xy[ord], colors=colors[ord], linewidths=linewidths,
                               path_effects=[matplotlib.patheffects.Stroke(capstyle="round")])
  
  return ax.add_collection(lines)

def add_text(text, ax):
  return plt.text(0.5, 1.01, text, horizontalalignment='center',
                  verticalalignment='bottom', transform=ax.transAxes)

def plot_protein(protein=None, pos=None, plddt=None, Ls=None, dpi=100, best_view=True, line_w=2.0):
  
  if protein is not None:
    pos = np.asarray(protein.atom_positions[:,1,:])
    plddt = np.asarray(protein.b_factors[:,0])

  # get best view
  if best_view:
    if plddt is not None:
      weights = plddt/100
      pos = pos - (pos * weights[:,None]).sum(0,keepdims=True) / weights.sum()
      pos = pos @ kabsch(pos, pos, weights, return_v=True)
    else:
      pos = pos - pos.mean(0,keepdims=True)
      pos = pos @ kabsch(pos, pos, return_v=True)

  if plddt is not None:
    fig, (ax1, ax2) = plt.subplots(1,2)
    fig.set_figwidth(6); fig.set_figheight(3)
    ax = [ax1, ax2]
  else:
    fig, ax1 = plt.subplots(1,1)
    fig.set_figwidth(3); fig.set_figheight(3)
    ax = [ax1]
    
  fig.set_dpi(dpi)
  fig.subplots_adjust(top = 0.9, bottom = 0.1, right = 1, left = 0, hspace = 0, wspace = 0)

  xy_min = pos[...,:2].min() - line_w
  xy_max = pos[...,:2].max() + line_w
  for a in ax:
    a.set_xlim(xy_min, xy_max)
    a.set_ylim(xy_min, xy_max)
    a.axis(False)

  if Ls is None or len(Ls) == 1:
    # color N->C
    c = np.arange(len(pos))[::-1]
    plot_pseudo_3D(pos,  line_w=line_w, ax=ax1)
    add_text("colored by N→C", ax1)
  else:
    # color by chain
    c = np.concatenate([[n]*L for n,L in enumerate(Ls)])
    if len(Ls) > 40:   plot_pseudo_3D(pos, c=c, line_w=line_w, ax=ax1)
    else:              plot_pseudo_3D(pos, c=c, cmap=pymol_cmap, cmin=0, cmax=39, line_w=line_w, ax=ax1)
    add_text("colored by chain", ax1)
    
  if plddt is not None:
    # color by pLDDT
    plot_pseudo_3D(pos, c=plddt, cmin=50, cmax=90, line_w=line_w, ax=ax2)
    add_text("colored by pLDDT", ax2)

  return fig
