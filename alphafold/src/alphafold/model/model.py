# Copyright 2021 DeepMind Technologies Limited
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#      http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""Code for constructing the model."""
from typing import Any, Mapping, Optional, Union

from absl import logging
from alphafold.common import confidence
from alphafold.model import features
from alphafold.model import modules
from alphafold.model import modules_multimer
from alphafold.common import residue_constants
import haiku as hk
import jax
import jax.numpy as jnp
import ml_collections
import numpy as np
import tensorflow.compat.v1 as tf
import tree

class RunModel:
  """Container for JAX model."""

  def __init__(self,
               config: ml_collections.ConfigDict,
               params: Optional[Mapping[str, Mapping[str, np.ndarray]]] = None,
               attention_output_dir: Optional[str] = None,
               is_training = False, extended_ptm_config=None):
    
    self.config = config
    self.params = params
    self.multimer_mode = config.model.global_config.multimer_mode
    self.config.model.calc_extended_ptm = extended_ptm_config['calc_extended_ptm'] if extended_ptm_config else False
    self.config.model.use_probs_extended = extended_ptm_config['use_probs_extended'] if extended_ptm_config else False
    self.attention_output_dir = attention_output_dir

    if self.multimer_mode:
      def _forward_fn(batch):
        model = modules_multimer.AlphaFold(self.config.model)
        return model(batch, is_training=is_training)
    else:
      def _forward_fn(batch):
        if self.config.data.eval.num_ensemble == 1:
          model = modules.AlphaFold_noE(self.config.model, attention_output_dir=self.attention_output_dir)
          return model(batch, is_training=is_training)
        else:
          model = modules.AlphaFold(self.config.model, attention_output_dir=self.attention_output_dir)
          return model(
              batch,
              is_training=is_training,
              compute_loss=False,
              ensemble_representations=True)

    self.apply = jax.jit(hk.transform(_forward_fn).apply)
    self.init = jax.jit(hk.transform(_forward_fn).init)

  def init_params(self, feat: features.FeatureDict, random_seed: int = 0):
    """Initializes the model parameters.

    If none were provided when this class was instantiated then the parameters
    are randomly initialized.

    Args:
      feat: A dictionary of NumPy feature arrays as output by
        RunModel.process_features.
      random_seed: A random seed to use to initialize the parameters if none
        were set when this class was initialized.
    """
    if not self.params:
      # Init params randomly.
      rng = jax.random.PRNGKey(random_seed)
      self.params = hk.data_structures.to_mutable_dict(
          self.init(rng, feat))
      logging.warning('Initialized parameters randomly')

  def process_features(
      self,
      raw_features: Union[tf.train.Example, features.FeatureDict],
      random_seed: int) -> features.FeatureDict:
    """Processes features to prepare for feeding them into the model.

    Args:
      raw_features: The output of the data pipeline either as a dict of NumPy
        arrays or as a tf.train.Example.
      random_seed: The random seed to use when processing the features.

    Returns:
      A dict of NumPy feature arrays suitable for feeding into the model.
    """

    if self.multimer_mode:
      return raw_features

    # Single-chain mode.
    if isinstance(raw_features, dict):
      return features.np_example_to_features(
          np_example=raw_features,
          config=self.config,
          random_seed=random_seed)
    else:
      return features.tf_example_to_features(
          tf_example=raw_features,
          config=self.config,
          random_seed=random_seed)

  def eval_shape(self, feat: features.FeatureDict) -> jax.ShapeDtypeStruct:
    self.init_params(feat)
    logging.debug('Running eval_shape with shape(feat) = %s',
                 tree.map_structure(lambda x: x.shape, feat))
    shape = jax.eval_shape(self.apply, self.params, jax.random.PRNGKey(0), feat)
    logging.info('Output shape was %s', shape)
    return shape

  def predict(self,
              feat: features.FeatureDict,
              random_seed: int = 0,
              return_representations: bool = False,
              callback: Any = None) -> Mapping[str, Any]:
    """Makes a prediction by inferencing the model on the provided features.

    Args:
      feat: A dictionary of NumPy feature arrays as output by
        RunModel.process_features.
      random_seed: The random seed to use when running the model. In the
        multimer model this controls the MSA sampling.

    Returns:
      A dictionary of model outputs.
    """
    self.init_params(feat)
    logging.info('Running predict with shape(feat) = %s',
                 tree.map_structure(lambda x: x.shape, feat))
    
    # get shapes
    aatype = feat["aatype"]
    num_iters = self.config.model.num_recycle + 1
    if self.multimer_mode:
      L = aatype.shape[0]
    else:
      num_ensemble = self.config.data.eval.num_ensemble
      L = aatype.shape[1]
    
    # initialize
    zeros = lambda shape: np.zeros(shape, dtype=np.float16)
    prev = {'prev_msa_first_row': zeros([L,256]),
            'prev_pair':          zeros([L,L,128])}
    
    # initial guess
    if "all_atom_positions" in feat:
      logging.info("INFO: using provided all_atom_positions as initial guess")
      prev["prev_pos"] = feat["all_atom_positions"]
    else:
      prev["prev_pos"] = np.zeros([L,37,3])
    
    def run(key, feat, prev):
      def _jnp_to_np(x):
        for k, v in x.items():
          if isinstance(v, dict):
            x[k] = _jnp_to_np(v)
          else:
            x[k] = np.asarray(v,np.float16)
        return x
      result = _jnp_to_np(self.apply(self.params, key, {**feat, "prev":prev}))
      prev = result.pop("prev")
      return result, prev


    # initialize random key
    key = jax.random.PRNGKey(random_seed)
    
    # iterate through recycles
    for r in range(num_iters):      
        # set recycle number for this iteration
        modules.set_recycle_number(r)

        # grab subset of features
        if self.multimer_mode:
          sub_feat = feat
        else:
          s = r * num_ensemble
          e = (r+1) * num_ensemble
          sub_feat = jax.tree_util.tree_map(lambda x:x[s:e], feat)
            
        # run
        key, sub_key = jax.random.split(key)
        result, prev = run(sub_key, sub_feat, prev)
        
        if return_representations:

          result["representations"] = {"pair":   prev["prev_pair"],
                                       "single": prev["prev_msa_first_row"]}
                                       
        # callback
        if callback is not None: callback(result, r)

        # decide when to stop
        if result["ranking_confidence"] > self.config.model.stop_at_score:
          break
        if r > 0 and result["tol"] < self.config.model.recycle_early_stop_tolerance:
          break
    logging.info('Output shape was %s', tree.map_structure(lambda x: x.shape, result))
    return result, r
