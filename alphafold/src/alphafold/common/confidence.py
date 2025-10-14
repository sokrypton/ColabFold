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

"""Functions for processing confidence metrics."""

import jax.numpy as jnp
import jax
import numpy as np
from alphafold.common import residue_constants
import scipy.special


def compute_tol(prev_pos, current_pos, mask, use_jnp=False):
    # Early stopping criteria based on criteria used in
    # AF2Complex: https://www.nature.com/articles/s41467-022-29394-2    
    _np = jnp if use_jnp else np
    dist = lambda x:_np.sqrt(((x[:,None] - x[None,:])**2).sum(-1))
    ca_idx = residue_constants.atom_order['CA']
    sq_diff = _np.square(dist(prev_pos[:,ca_idx])-dist(current_pos[:,ca_idx]))
    mask_2d = mask[:,None] * mask[None,:]
    return _np.sqrt((sq_diff * mask_2d).sum()/mask_2d.sum() + 1e-8)


def compute_plddt(logits, use_jnp=False):
  """Computes per-residue pLDDT from logits.
  Args:
    logits: [num_res, num_bins] output from the PredictedLDDTHead.
  Returns:
    plddt: [num_res] per-residue pLDDT.
  """
  if use_jnp:
    _np, _softmax = jnp, jax.nn.softmax
  else:
    _np, _softmax = np, scipy.special.softmax
  
  num_bins = logits.shape[-1]
  bin_width = 1.0 / num_bins
  bin_centers = _np.arange(start=0.5 * bin_width, stop=1.0, step=bin_width)
  probs = _softmax(logits, axis=-1)
  predicted_lddt_ca = (probs * bin_centers[None, :]).sum(-1)
  return predicted_lddt_ca * 100

def _calculate_bin_centers(breaks, use_jnp=False):
  """Gets the bin centers from the bin edges.
  Args:
    breaks: [num_bins - 1] the error bin edges.
  Returns:
    bin_centers: [num_bins] the error bin centers.
  """
  _np = jnp if use_jnp else np
  step = breaks[1] - breaks[0]

  # Add half-step to get the center
  bin_centers = breaks + step / 2

  # Add a catch-all bin at the end.
  return _np.append(bin_centers, bin_centers[-1] + step)

def _calculate_expected_aligned_error(
  alignment_confidence_breaks,
  aligned_distance_error_probs,
  use_jnp=False):
  """Calculates expected aligned distance errors for every pair of residues.
  Args:
    alignment_confidence_breaks: [num_bins - 1] the error bin edges.
    aligned_distance_error_probs: [num_res, num_res, num_bins] the predicted
      probs for each error bin, for each pair of residues.
  Returns:
    predicted_aligned_error: [num_res, num_res] the expected aligned distance
      error for each pair of residues.
    max_predicted_aligned_error: The maximum predicted error possible.
  """
  bin_centers = _calculate_bin_centers(alignment_confidence_breaks, use_jnp=use_jnp)
  # Tuple of expected aligned distance error and max possible error.
  pae = (aligned_distance_error_probs * bin_centers).sum(-1)
  return (pae, bin_centers[-1])

def compute_predicted_aligned_error(logits, breaks, use_jnp=False):
  """Computes aligned confidence metrics from logits.
  Args:
    logits: [num_res, num_res, num_bins] the logits output from
      PredictedAlignedErrorHead.
    breaks: [num_bins - 1] the error bin edges.

  Returns:
    aligned_confidence_probs: [num_res, num_res, num_bins] the predicted
      aligned error probabilities over bins for each residue pair.
    predicted_aligned_error: [num_res, num_res] the expected aligned distance
      error for each pair of residues.
    max_predicted_aligned_error: The maximum predicted error possible.
  """
  _softmax = jax.nn.softmax if use_jnp else scipy.special.softmax
  aligned_confidence_probs = _softmax(logits,axis=-1)
  predicted_aligned_error, max_predicted_aligned_error = \
  _calculate_expected_aligned_error(breaks, aligned_confidence_probs, use_jnp=use_jnp)

  return {
      'aligned_confidence_probs': aligned_confidence_probs,
      'predicted_aligned_error': predicted_aligned_error,
      'max_predicted_aligned_error': max_predicted_aligned_error,
  }

def predicted_tm_score(logits, breaks, residue_weights = None,
    asym_id = None, use_jnp=False):
  """Computes predicted TM alignment or predicted interface TM alignment score.

  Args:
    logits: [num_res, num_res, num_bins] the logits output from
      PredictedAlignedErrorHead.
    breaks: [num_bins] the error bins.
    residue_weights: [num_res] the per residue weights to use for the
      expectation.
    asym_id: [num_res] the asymmetric unit ID - the chain ID. Only needed for
      ipTM calculation.

  Returns:
    ptm_score: The predicted TM alignment or the predicted iTM score.
  """
  if use_jnp:
    _np, _softmax = jnp, jax.nn.softmax
  else:
    _np, _softmax = np, scipy.special.softmax

  # residue_weights has to be in [0, 1], but can be floating-point, i.e. the
  # exp. resolved head's probability.
  if residue_weights is None:
    residue_weights = _np.ones(logits.shape[0])

  bin_centers = _calculate_bin_centers(breaks, use_jnp=use_jnp)
  num_res = residue_weights.shape[0]

  # Clip num_res to avoid negative/undefined d0.
  clipped_num_res = _np.maximum(residue_weights.sum(), 19)

  # Compute d_0(num_res) as defined by TM-score, eqn. (5) in Yang & Skolnick
  # "Scoring function for automated assessment of protein structure template
  # quality", 2004: http://zhanglab.ccmb.med.umich.edu/papers/2004_3.pdf
  d0 = 1.24 * (clipped_num_res - 15) ** (1./3) - 1.8

  # Convert logits to probs.
  probs = _softmax(logits, axis=-1)

  # TM-Score term for every bin.
  tm_per_bin = 1. / (1 + _np.square(bin_centers) / _np.square(d0))
  # E_distances tm(distance).
  predicted_tm_term = (probs * tm_per_bin).sum(-1)

  if asym_id is None:
    pair_mask = _np.full((num_res,num_res),True)
  else:
    pair_mask = asym_id[:, None] != asym_id[None, :]

  predicted_tm_term *= pair_mask

  pair_residue_weights = pair_mask * (residue_weights[None, :] * residue_weights[:, None])
  normed_residue_mask = pair_residue_weights / (1e-8 + pair_residue_weights.sum(-1, keepdims=True))
  per_alignment = (predicted_tm_term * normed_residue_mask).sum(-1)

  return (per_alignment * residue_weights).max()

def get_confidence_metrics(prediction_result, mask, rank_by = "plddt", use_jnp=False, keep_pae=False):
  """Post processes prediction_result to get confidence metrics."""  
  confidence_metrics = {}
  plddt = compute_plddt(prediction_result['predicted_lddt']['logits'], use_jnp=use_jnp)
  confidence_metrics['plddt'] = plddt  
  confidence_metrics["mean_plddt"] = (plddt * mask).sum()/mask.sum()

  if 'predicted_aligned_error' in prediction_result:
    if keep_pae:
        prediction_result['pae_matrix_with_logits'] = prediction_result['predicted_aligned_error']

    confidence_metrics.update(compute_predicted_aligned_error(
        logits=prediction_result['predicted_aligned_error']['logits'],
        breaks=prediction_result['predicted_aligned_error']['breaks'],
        use_jnp=use_jnp))

    confidence_metrics['ptm'] = predicted_tm_score(
        logits=prediction_result['predicted_aligned_error']['logits'],
        breaks=prediction_result['predicted_aligned_error']['breaks'],
        residue_weights=mask,
        use_jnp=use_jnp)

    if "asym_id" in prediction_result["predicted_aligned_error"]:
      # Compute the ipTM only for the multimer model.
      confidence_metrics['iptm'] = predicted_tm_score(
          logits=prediction_result['predicted_aligned_error']['logits'],
          breaks=prediction_result['predicted_aligned_error']['breaks'],
          residue_weights=mask,
          asym_id=prediction_result['predicted_aligned_error']['asym_id'],
          use_jnp=use_jnp)

  # compute mean_score
  if rank_by == "multimer":
    mean_score = 80 * confidence_metrics["iptm"] + 20 * confidence_metrics["ptm"]
  elif rank_by == "iptm":
    mean_score = 100 * confidence_metrics["iptm"]
  elif rank_by == "ptm":
    mean_score = 100 * confidence_metrics["ptm"]
  else:
    mean_score = confidence_metrics["mean_plddt"]
  confidence_metrics["ranking_confidence"] = mean_score
  return confidence_metrics