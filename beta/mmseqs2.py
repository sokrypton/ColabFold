import tqdm.notebook
import tarfile
import time
import requests
import random

def run_mmseqs2(x, prefix, use_env=True, filter=False):
  
  def submit(seqs, mode, N=1):
    
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
    with tqdm.notebook.tqdm(total=300, bar_format=TQDM_BAR_FORMAT) as pbar:
      N,REDO = 1,True
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
          if TIME > 900 and out["status"] != "COMPLETE":
            # something failed on the server side, need to resubmit
            N += 1
            break
        
        if out["status"] == "COMPLETE":
          if TIME < 300:
            pbar.update(n=(300-TIME))
          REDO = False
          
      # Download results
      download(ID, tar_gz_file)

  # prep list of a3m files
  a3m_files = [f"{path}/uniref.a3m"]
  if use_env: a3m_files.append(f"{path}/bfd.mgnify30.metaeuk30.smag30.a3m")
  
  # extract a3m files
  if not os.path.isfile(a3m_files[0])
    with tarfile.open(tar_gz_file) as tar_gz:
      tar_gz.extractall(path)  
    
  a3m_lines = {n:[] for n in range(N,N+len(seqs))}
  for a3m_file in a3m_files:
    update_M,M = True,None
    for line in open(a3m_file,"r"):
      if len(line) > 0:
        if "\x00" in line:
          update_M = True
        else:
          if line.startswith(">") and update_M:
            M = int(line[1:].rstrip())
            update_M = False
          a3m_lines[M].append(line)
  if isinstance(x, str):
    return "".join(a3m_lines[N])
  else:
    return ["".join(a3m_lines[n]) for n in range(N,N+len(x))]
