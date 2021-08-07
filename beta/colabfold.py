import matplotlib.pyplot as plt
import matplotlib
import matplotlib.patheffects
from matplotlib import collections as mcoll
def kabsch(P, Q, return_v=False):
  # code borrowed from: https://github.com/charnley/rmsd/blob/master/rmsd/calculate_rmsd.py
  V, S, W = np.linalg.svd(P.T @ Q)
  if (np.linalg.det(V) * np.linalg.det(W)) < 0.0:
    S[-1] = -S[-1]
    V[:,-1] = -V[:,-1]
  if return_v: return V
  else: return np.dot(V, W)

def plot_protein(protein=None, pos=None, plddt=None, Ls=None, dpi=100):
  
  if protein is not None:
    pos = protein.atom_positions[:,1,:]
    plddt = protein.b_factors[:,0]

  # get best view
  pos = pos - pos.mean(0,keepdims=True)
  pos = pos @ kabsch(pos,pos,return_v=True)

  def make_segments(x, y):
    points = np.array([x, y]).T.reshape(-1, 1, 2)
    segments = np.concatenate([points[:-1], points[1:]], axis=1)
    return segments
  
  def get_color(x, alpha=None, tint=None, shade=None, vmin=50, vmax=90, cmap="gist_rainbow"):
    color_map = matplotlib.cm.get_cmap(cmap)
    if x < vmin: x = vmin
    if x > vmax: x = vmax
    x = (x - vmin)/(vmax - vmin)
    if cmap == "gist_rainbow": x *= 0.8
    color = np.array(color_map(x, alpha=alpha))
    if tint is not None:
      color[:3] = color[:3] + (1 - color[:3]) * tint
    if shade is not None:
      color[:3] = color[:3] * shade
    return color

  if plddt is not None:
    fig, (ax1, ax2) = plt.subplots(1,2)
    fig.set_figwidth(6); fig.set_figheight(3)
  else:
    fig, ax1 = plt.subplots(1,1)
    fig.set_figwidth(3); fig.set_figheight(3)
    
  fig.set_dpi(dpi)
  fig.subplots_adjust(top = 0.9, bottom = 0.1, right = 1, left = 0, hspace = 0, wspace = 0)

  range_xy = (pos[:,:2].min(),pos[:,:2].max())
  range_z = (pos[:,-1].min(),pos[:,-1].max())
  z = (pos[:,-1] - range_z[0]) / (range_z[1] - range_z[0])

  srt = (z[:-1]+z[1:]).argsort()
  seg = make_segments(pos[:,0],pos[:,1])
  alpha = (np.sqrt(np.square(pos[:-1] - pos[1:]).sum(-1)) < 4).astype(float)

  def plot_lines(ax, vmin, vmax, values, title=None, cmap="gist_rainbow"):
    shade = (z + 2)/3
    tint = z/3
    colors = np.array([get_color(p, vmin=vmin, vmax=vmax,
                                 alpha=a, shade=s, tint=t, cmap=cmap) for p,a,s,t in zip(values,alpha,shade,tint)])
    plt.text(0.5, 1.01, title, horizontalalignment='center', verticalalignment='bottom', transform=ax.transAxes)
    ax.axis('scaled')
    ax.set_xlim(*range_xy); ax.set_ylim(*range_xy); ax.axis(False)
    ax.add_collection(mcoll.LineCollection(seg[srt], colors=colors[srt], linewidths=5,
                                           path_effects=[matplotlib.patheffects.Stroke(capstyle="round")]))
  if Ls is None or len(Ls) == 1:
    c = np.arange(len(pos))[::-1]
    plot_lines(ax1, c.min(), c.max(), c, "colored by N->C")
  else:
    c = np.concatenate([[n]*L for n,L in enumerate(Ls)])
    plot_lines(ax1, 0, 9, c, "colored by chain", cmap="Set1")

  if plddt is not None: plot_lines(ax2, 50, 90, plddt, "colored by pLDDT")
  plt.show()
