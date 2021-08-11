import tqdm.notebook
import tarfile
import time
import requests
import random

def run(query_sequence, prefix, use_env=True, filter=False):
  
  def submit(query_sequence, mode, N=1):
    res = requests.post('https://a3m.mmseqs.com/ticket/msa', data={'q':f">{N}\n{query_sequence}", 'mode': mode})
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
    
  if filter:
    mode = "env" if use_env else "all"
  else:
    mode = "env-nofilter" if use_env else "nofilter"
  
  path = f"{prefix}_{mode}"
  if not os.path.isdir(path): os.mkdir(path)

  # call mmseqs2 api
  tar_gz_file = f'{path}/out.tar.gz'
  if not os.path.isfile(tar_gz_file):
    with tqdm.notebook.tqdm(bar_format='{l_bar}{bar}') as pbar:
      N,REDO = 1,True
      while REDO:
        pbar.set_description("SUBMIT")
        
        # Resubmit job until it goes through
        out = submit(query_sequence, mode, N)
        while out["status"] in ["UNKNOWN","RATELIMIT"]:
          # resubmit
          time.sleep(5 + random.randint(0,5))
          out = submit(query_sequence, mode, N)

        # wait for job to finish
        ID = out["id"]
        TIME = 0
        pbar.set_description(out["status"])
        while out["status"] in ["UNKNOWN","RUNNING","PENDING"]:
          t = 5 + random.randint(0,5)
          time.sleep(t)
          TIME += t
          out = status(ID)    
          pbar.set_description(out["status"])
          if TIME > 600 and out["status"] != "COMPLETE":
            # something failed on the server side, need to resubmit
            N += 1
            break
        
        if out["status"] == "COMPLETE":
          REDO = False

    download(ID, tar_gz_file)

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
