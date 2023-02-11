import json
import lzma
import os
import pickle
from pathlib import Path
from typing import List, Tuple, Mapping, Any, Dict


import jax
import jax.numpy as jnp
import numpy as np
from alphafold.model.features import FeatureDict
from alphafold.model.model import RunModel
from colabfold.colabfold import run_mmseqs2

def jnp_to_np(output: Dict[str, Any]) -> Dict[str, Any]:
  """Recursively changes jax arrays to numpy arrays."""
  for k, v in output.items():
    if isinstance(v, dict):
      output[k] = jnp_to_np(v)
    elif isinstance(v, jnp.ndarray):
      output[k] = np.array(v)
  return output

# Copy the original method before mocking
original_run_model = RunModel.predict

class MockRunModel:
  """Mocks FeatureDict -> prediction
  The class is stateful, i.e. predictions need to be done in the given order
  msa_feat is a) large and b) has some variance between machines, so we ignore it
  """

  fixture_dir: Path
  predictions: List[str]
  pos: int

  def __init__(self, fixture_dir: Path, predictions: List[str]):
    self.fixture_dir = fixture_dir
    self.predictions = predictions
    self.pos = 0

  def predict(
    self,
    model_runner: RunModel,
    feat: FeatureDict,
    random_seed: int,
    return_representations: bool = False,
    callback: Any = None
  ) -> Mapping[str, Any]:
    """feat["msa"] or feat["msa_feat"] for normal/complexes is non-deterministic, so we remove it before storing,
    but we keep it for predicting or returning, where we need it for plotting"""

    feat_file = self.fixture_dir.joinpath(self.predictions[self.pos]).joinpath("model_feat.pkl.xz")
    pred_file = self.fixture_dir.joinpath(self.predictions[self.pos]).joinpath("model_pred.pkl.xz")

    if os.environ.get("PRED_TEST") or not pred_file.is_file():
      pred, recycles = original_run_model(model_runner, feat)
      pred = jnp_to_np(pred)
        
    if not feat_file.is_file() or not pred_file.is_file():
      print("updating snapshots...")
      prev_feat = feat
      prev_pred = pred
      with lzma.open(feat_file,"wb") as fp:
        pickle.dump(prev_feat, fp)
      with lzma.open(pred_file,"wb") as fp:
        pickle.dump(prev_pred, fp)
    
    else:
      with lzma.open(feat_file) as handle:
        prev_feat = pickle.load(handle)
      with lzma.open(pred_file) as handle:
        prev_pred = pickle.load(handle)
    
    def cmp_dict(x,y):
      ''' check if two dictionaries are "allclose" '''
      
      def chk(a,b):
        test = []
        for k,v in a.items():
          if k in b: # TODO
            if isinstance(v, dict):
              test.append(chk(v,b[k]))
            else:
              if not np.allclose(v,b[k]):
                print("--------------------")
                print(k)
                print(v)
                print(b[k])
                print("--------------------")
              test.append(np.allclose(v,b[k]))
        return test
      
      return all(jax.tree_util.tree_flatten(chk(x,y))[0])
    
    # test input features match
    assert cmp_dict(prev_feat, feat)

    # test output predictions match
    if os.environ.get("PRED_TEST"):
      assert cmp_dict(prev_pred, pred)

    self.pos += 1
    return prev_pred, 3

class MMseqs2Mock:
  """Mocks out the call to the mmseqs2 api

  Each test has its own json file which contains the run_mmseqs2 input data in the
  config field and the saved response. To update responses or to add new tests,
  set the UPDATE_SNAPSHOTS env var (e.g. `UPDATE_SNAPSHOTS=1 pytest`
  """

  data_file: Path
  saved_responses: List[Dict[str, Any]]

  def __init__(self, rootpath: Path, name: str):
    self.data_file = (
      rootpath.joinpath("test-data/mmseqs-api-reponses")
      .joinpath(name)
      .with_suffix(".json")
    )
    if os.environ.get("UPDATE_SNAPSHOTS") and not self.data_file.is_file():
      # TODO: call mmseqs2 server??
      self.data_file.write_text("[]")
    with self.data_file.open() as fp:
      self.saved_responses = []
      for saved_response in json.load(fp):
        # Join lines we've split before
        response = join_lines(saved_response["response"])
        self.saved_responses.append(
          {"config": saved_response["config"], "response": response}
        )

  def mock_run_mmseqs2(
    self,
    query,
    prefix,
    use_env=True,
    use_filter=True,
    use_templates=False,
    filter=None,
    use_pairing=False,
    host_url="https://a3m.mmseqs.com",
  ):
    assert prefix
    config = {
      "query": query,
      "use_env": use_env,
      "use_filter": use_filter,
      "use_templates": use_templates,
      "filter": filter,
      "use_pairing": use_pairing,
    }

    for saved_response in self.saved_responses:
      if saved_response["config"] == config:
        return saved_response["response"]

    if os.environ.get("UPDATE_SNAPSHOTS"):
      print(f"\nrun_mmseqs2 with {config}")
      response = run_mmseqs2(
        x=config["query"],
        prefix=prefix,
        use_env=config["use_env"],
        use_filter=config["use_filter"],
        use_templates=config["use_templates"],
        filter=config["filter"],
        use_pairing=config["use_pairing"],
        host_url=host_url,
      )
      # Split lines so we get a readable json file
      response = split_lines(response)
      self.saved_responses.append({"config": config, "response": response})
      self.data_file.write_text(json.dumps(self.saved_responses, indent=2))
    else:
      assert False, config


def split_lines(x):
  """Split each files into a list of lines"""
  if isinstance(x, list):
    return [split_lines(i) for i in x]
  elif isinstance(x, str):
    return x.splitlines()
  else:
    raise TypeError(f"{type(x)} {str(x)[:20]}")


def join_lines(x):
  """Inverse of split_lines"""
  if all(isinstance(i, str) for i in x):
    return "\n".join(x)
  elif all(isinstance(i, list) for i in x):
    return [join_lines(i) for i in x]
  else:
    raise TypeError(f"{[type(i) for i in x]} {str(x)[:20]}")
