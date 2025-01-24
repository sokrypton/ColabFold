import numpy as np
import string
from copy import deepcopy
import jax.numpy as jnp
import jax
import scipy
import pandas as pd
import matplotlib.pyplot as plt

import os
from alphafold.common import confidence

""" Functions to calculate interface metrics on actual interfaces"""
def get_dgram_bins(result):
    """From colabdesign, calculate bin boundaries of distogram"""
    dgram = result["distogram"]["logits"]
    if dgram.shape[-1] == 64:
        dgram_bins = jnp.append(0,jnp.linspace(2.3125,21.6875,63))
    if dgram.shape[-1] == 39:
        dgram_bins = jnp.linspace(3.25,50.75,39) + 1.25
    return dgram_bins


def get_contact_map(result, dist=8.0):
    """From colabdesign, get contact map from distogram"""
    dist_logits = result["distogram"]["logits"]
    dist_bins = get_dgram_bins(result)
    return (jax.nn.softmax(dist_logits) * (dist_bins < dist)).sum(-1)


def get_chain_indices(asym_id, use_jnp=True):
    """Returns a list of tuples indicating the start and end indices for each chain."""

    chain_starts_ends = []
    unique_chains = np.unique(asym_id) # chains are numbered 0, 1, 2, ...

    for chain in unique_chains:
        positions = np.where(asym_id == chain)[0]
        chain_starts_ends.append((positions[0], positions[-1]))

    return chain_starts_ends


def predicted_tm_score_modified(logits, breaks, residue_weights=None,
                                asym_id=None, pair_residue_weights=None, use_jnp=False, ):
    """Computes predicted TM alignment or predicted interface TM alignment score.

    Args:
      logits: [num_res, num_res, num_bins] the logits output from
        PredictedAlignedErrorHead.
      breaks: [num_bins] the error bins.
      residue_weights: [num_res] the per residue weights to use for the
        expectation.
      asym_id: [num_res] the asymmetric unit ID - the chain ID. Only needed for
        ipTM calculation.
      pair_residue_weights: [num_res, num_res] unnormalized weights for actifptm calculation

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
    bin_centers = confidence._calculate_bin_centers(breaks, use_jnp=use_jnp)
    num_res = residue_weights.shape[0]

    # Clip num_res to avoid negative/undefined d0.
    clipped_num_res = _np.maximum(num_res, 19)

    # Compute d_0(num_res) as defined by TM-score, eqn. (5) in Yang & Skolnick
    # "Scoring function for automated assessment of protein structure template
    # quality", 2004: http://zhanglab.ccmb.med.umich.edu/papers/2004_3.pdf
    d0 = 1.24 * (clipped_num_res - 15) ** (1. / 3) - 1.8

    # Convert logits to probs.
    probs = _softmax(logits, axis=-1)

    # TM-Score term for every bin.
    tm_per_bin = 1. / (1 + _np.square(bin_centers) / _np.square(d0))

    # E_distances tm(distance).
    predicted_tm_term = (probs * tm_per_bin).sum(-1)

    if asym_id is None:
        pair_mask = _np.full((num_res, num_res), True)
    else:
        pair_mask = asym_id[:, None] != asym_id[None, :]

    predicted_tm_term *= pair_mask

    # If pair_residue_weights is provided (e.g. for if_ptm with contact probabilities),
    # it should not be overwritten
    if pair_residue_weights is None:
        pair_residue_weights = pair_mask * (residue_weights[None, :] * residue_weights[:, None])
    normed_residue_mask = pair_residue_weights / (1e-8 + pair_residue_weights.sum(-1, keepdims=True))

    per_alignment = (predicted_tm_term * normed_residue_mask).sum(-1)
    residuewise_iptm = per_alignment * residue_weights

    return residuewise_iptm


def get_ptm_modified(inputs, outputs, interface=False):
    """ This function is the same as in the original AF2, just calls a modified TM-score calculation function."""

    pae = {"residue_weights": inputs["seq_mask"], **outputs["predicted_aligned_error"]}

    if interface:
        if "asym_id" not in pae:
            pae["asym_id"] = inputs["asym_id"]
    else:
        if "asym_id" in pae:
            pae.pop("asym_id")

    return predicted_tm_score_modified(**pae, use_jnp=True)


def get_actifptm_probs(result, asym_id, cmap, start_i, end_i, start_j, end_j):
    """
    This function calculates the interface PTM score for a given interface, taking into account the contact probabilities, not binary contacts.

    Args:
        af: AlphaFold object
        cmap: Contact map
        start_i, end_i: Start and end indices of the first chain
        start_j, end_j: Start and end indices of the second chain

    Returns:
        actifptm: Interface pTM score
    """

    total_length = len(asym_id)
    outputs = deepcopy(result)

    # Create a new matrix, which contains only the contacts between the two chains
    cmap_copy = np.zeros((total_length, total_length))
    cmap_copy[start_i:end_i + 1, start_j:end_j + 1] = cmap[start_i:end_i + 1, start_j:end_j + 1]
    cmap_copy[start_j:end_j + 1, start_i:end_i + 1] = cmap[start_j:end_j + 1, start_i:end_i + 1]

    # this is for the full-length actifptm
    if end_i == end_j == total_length - 1 and start_i == start_j == 0:
        pair_mask = asym_id[:, None] != asym_id[None, :]
        cmap_copy *= pair_mask

    # Update seq_mask for these positions to True within inputs
    seq_mask = np.full(total_length, 0, dtype=float)
    seq_mask[np.concatenate((np.arange(start_i, end_i + 1),
                             np.arange(start_j, end_j + 1)))] = 1

    # Call get_ptm with updated inputs and outputs
    pae = {"residue_weights": seq_mask,
           **outputs["predicted_aligned_error"],
           "asym_id": asym_id}
    residuewise_actifptm = predicted_tm_score_modified(**pae, use_jnp=True, pair_residue_weights=cmap_copy)

    return residuewise_actifptm


def get_actifptm_contacts(result, asym_id, cmap, start_i, end_i, start_j, end_j):
    """
    This function calculates the interface PTM score for a given interface, taking into account binary contacts.
    In case of no confident contacts, the interface PTM score is set to 0.

    Args:
        af: AlphaFold object
        cmap: Contact map
        start_i, end_i: Start and end indices of the first chain
        start_j, end_j: Start and end indices of the second chain

    Returns:
        actifptm: Interface pTM score
    """

    # Prepare a dictionary to collect the inputs for calculation
    inputs_actifptm = {}
    contacts = np.where(cmap[start_i:end_i + 1, start_j:end_j + 1] >= 0.6)
    total_length = len(asym_id)
    outputs = deepcopy(result)

    if contacts[0].size > 0:  # If there are contacts
        # Convert local chain positions back to global positions using JAX
        global_i_positions = contacts[0] + start_i
        global_j_positions = contacts[1] + start_j
        global_positions = list(set(np.concatenate((global_i_positions, global_j_positions))))
        global_positions = np.array(global_positions, dtype=int)
        global_positions.sort()

        # Initialize new input dictionary
        inputs_actifptm['asym_id'] = asym_id
        inputs_actifptm['seq_mask'] = np.full(total_length, 0, dtype=float)
        inputs_actifptm['seq_mask'][global_positions] = 1

        # Call get_ptm with updated inputs and outputs
        residuewise_actifptm = get_ptm_modified(inputs_actifptm, outputs, interface=True)
    else:
        residuewise_actifptm = np.array([0.0])
        inputs_actifptm['seq_mask'] = np.full(total_length, 0, dtype=float)

    return residuewise_actifptm, inputs_actifptm['seq_mask']


def get_pairwise_iptm(result, asym_id, start_i, end_i, start_j, end_j):
    """This will calculate ipTM as usual, just between given chains"""

    input_pairwise_iptm = {}

    # Prepare inputs
    outputs = deepcopy(result)
    input_pairwise_iptm['seq_mask'] = np.full(len(asym_id), 0, dtype=float)
    input_pairwise_iptm['asym_id'] = asym_id
    input_pairwise_iptm['seq_mask'][np.concatenate((np.arange(start_i, end_i + 1),
                                                    np.arange(start_j, end_j + 1)))] = 1

    return get_ptm_modified(input_pairwise_iptm, outputs, interface=True)


def get_per_chain_ptm(result, cmap, start, end):
    """
    Calculates the chain PTM score for a specified interface region, using contact probabilities.

    Args:
        af: AlphaFold object
        cmap: Contact map
        start: Start index of the interface region
        end: End index of the interface region

    Returns:
        cpTM: Chain pTM score
    """
    # Extract only the relevant subset of the contact map
    cmap_subset = cmap[start:end + 1, start:end + 1]

    # Extract only the relevant subset of the logits map
    pae_copy = deepcopy(result)["predicted_aligned_error"]
    pae_copy['logits'] = pae_copy['logits'][start:end + 1, start:end + 1, :]

    # Prepare inputs for the modified predicted_tm_score function
    pae_copy['residue_weights'] = np.ones(end - start + 1, dtype=float)
    pae_copy.pop('asym_id', None)

    # Calculate and return chain PTM score
    cptm = round(float(predicted_tm_score_modified(**pae_copy, use_jnp=True).max()), 3)

    return cptm


def get_chain_and_interface_metrics(result, asym_id, use_probs_extra=False, use_jnp=True):
    """
    This function iterates over all pairs of chains and calculates the interface and interchain PTM score for each pair.

    Args:
        result: The result from AlphaFold.
        asym_id: Array indicating chain boundaries.
        use_probs_extra: If True, calculate interface pTM score based on contact probabilities. Default is False.
        use_jnp: If True, use JAX numpy. Default is True.
    Returns:
        a dictionary with the pairwise interface pTM-s, and the chain-wise pTM.
        returns None for each, if there was an error finding the logits for the pae matrix
    """
    # this is to deal with the ptm models (af2 monomer)
    if len(asym_id.shape) > 1:
      asym_id = asym_id[0]

    full_length = len(asym_id)
    # Prepare dictionaries to collect results
    output = {'pairwise_actifptm': {}, 'pairwise_iptm': {}, 'per_chain_ptm': {}}
    #residuewise_output = {'residuewise_actifptm': {}, 'residuewise_iptm': {}}
    chain_starts_ends = get_chain_indices(asym_id, use_jnp=use_jnp)
    pair_residue_weights_no_probs = np.zeros((full_length, full_length), dtype=float)

    # Generate chain labels (A, B, C, ...)
    chain_labels = list(string.ascii_uppercase)

    # Define interface with 8A between Cb-s
    cmap = get_contact_map(result, 8)

    # This is for compatibility between colabdesign and colabfold
    results = {}
    if isinstance(result['predicted_aligned_error'], (np.ndarray, list)):
        if 'pae_matrix_with_logits' in result.keys():
            results['predicted_aligned_error'] = deepcopy(result['pae_matrix_with_logits'])
        else:
            print('There was an error retrieving the predicted aligned error matrix.')
            return {"pairwise_actifptm": None, "pairwise_iptm": None, "per_chain_ptm": None, 'actifptm': None}
    else:
        results['predicted_aligned_error'] = deepcopy(result['predicted_aligned_error'])

    for i, (start_i, end_i) in enumerate(chain_starts_ends):
        chain_label_i = chain_labels[i % len(chain_labels)]  # Wrap around if more than 26 chains
        for j, (start_j, end_j) in enumerate(chain_starts_ends):
            chain_label_j = chain_labels[j % len(chain_labels)]  # Wrap around if more than 26 chains
            if i < j:  # Avoid self-comparison and duplicate comparisons
                key = f"{chain_label_i}-{chain_label_j}"
                if not use_probs_extra:
                    residuewise_actifptm, seq_mask = get_actifptm_contacts(results, asym_id, cmap, start_i, end_i, start_j, end_j)
                    pair_residue_weights_no_probs += seq_mask[None, :] * seq_mask[:, None]
                    output['pairwise_actifptm'][key] = round(float(residuewise_actifptm.max()), 3)
                    #residuewise_output['residuewise_actifptm'][key] = residuewise_actifptm
                else:
                    residuewise_actifptm = get_actifptm_probs(results, asym_id, cmap, start_i, end_i, start_j, end_j)
                    output['pairwise_actifptm'][key] = round(float(residuewise_actifptm.max()), 3)
                    #residuewise_output['residuewise_actifptm'][key] = residuewise_actifptm


                # Also add regular i_ptm (interchain), pairwise
                residuewise_iptm = get_pairwise_iptm(results, asym_id, start_i, end_i, start_j, end_j)
                output['pairwise_iptm'][key] = round(float(residuewise_iptm.max()), 3)
                #residuewise_output['residuewise_iptm'][key] = output['pairwise_iptm'][key]

        # Also calculate pTM score for single chain
        output['per_chain_ptm'][chain_label_i] = get_per_chain_ptm(results, cmap, start_i, end_i)

    if not use_probs_extra:
        # we need to recreate the full matrix from the previously calculated contacts
        pair_mask = asym_id[:, None] != asym_id[None, :]
        pair_residue_weights_no_probs *= pair_mask
        pae = {"residue_weights": np.full(full_length, 1, dtype=float),
               **results["predicted_aligned_error"],
               "asym_id": asym_id}
        output['actifptm'] = round(float(predicted_tm_score_modified(**pae, pair_residue_weights=pair_residue_weights_no_probs).max()), 3)
    else:
        output['actifptm'] = round(float(get_actifptm_probs(results, asym_id, cmap, 0, full_length - 1, 0, full_length - 1).max()), 3)

    return output


def plot_matrix(actifptm_dict, iptm_dict, cptm_dict, prefix='rank', ax_in=None, fig_path=None):
    """This function plots the metrics in a matrix. The diagonal will be chain-wise pTM-s,
    the lower triangle displays actifptm and the upper triangle the ipTM (calculated in the original way)."""
    if not ax_in: # In case, we are not plotting multiple models next to each other
        fig, ax = plt.subplots(1, 1, figsize=(5, 5), squeeze=False)

    # get chain ids
    letters = sorted(set([key.split('-')[0] for key in actifptm_dict.keys()] + [key.split('-')[1] for key in actifptm_dict.keys()]))

    # initialize dataframe
    data = pd.DataFrame(np.zeros((len(letters), len(letters))), index=letters, columns=letters)

    # populate dataframe with ifptm and iptm values
    for key, value in actifptm_dict.items():
        i, j = key.split('-')
        data.loc[j, i] = value
        data.loc[i, j] = iptm_dict[f'{i}-{j}']

    # add cptm values to the dataframe
    for chain, value in cptm_dict.items():
        if chain in data.index:
            data.loc[chain, chain] = value

    # create masks for lower, upper triangles, and diagonal
    mask_upper = np.triu(np.ones(data.shape), k=1)
    mask_lower = np.tril(np.ones(data.shape), k=-1)
    mask_diagonal = np.eye(data.shape[0])

    dyn_size_ch = max(- 1.5 * len(letters) + 18, 3)   # resize the font with differently sized figures
    # Plot lower triangle (actifpTM)
    ax_in.imshow(np.ma.masked_where(mask_upper + mask_diagonal, data), cmap='Blues', vmax=1, vmin=0)

    # Plot upper triangle (ipTM)
    ax_in.imshow(np.ma.masked_where(mask_lower + mask_diagonal, data), cmap='Reds', vmax=1, vmin=0)

    # Plot diagonal (cpTM)
    diagonal_data = np.diag(np.diag(data))
    im = ax_in.imshow(np.ma.masked_where(~mask_diagonal.astype(bool), diagonal_data), cmap='Greys', vmax=1, vmin=0)

    # Add colorbar for cpTM (diagonal)
    cbar = plt.colorbar(im, ax=ax_in, fraction=0.046, pad=0.04)
    cbar.ax.tick_params(labelsize=dyn_size_ch)  # Set fontsize for colorbar labels
    cbar.outline.set_edgecolor('grey')
    cbar.outline.set_linewidth(0.5)

    # Add text annotations
    for i in range(data.shape[0]):
        for j in range(data.shape[1]):
            value = data.iloc[i, j]
            if not mask_upper[i, j] and not mask_diagonal[i, j]:
                text_color = 'white' if value > 0.8 else 'black'
                ax_in.text(j, i, f"{value:.2f}", ha='center', va='center', color=text_color, fontsize=dyn_size_ch*1.2)
            elif not mask_lower[i, j] and not mask_diagonal[i, j]:
                text_color = 'white' if value > 0.8 else 'black'
                ax_in.text(j, i, f"{value:.2f}", ha='center', va='center', color=text_color, fontsize=dyn_size_ch*1.2)
            elif mask_diagonal[i, j]:
                text_color = 'white' if value > 0.5 else 'black'
                ax_in.text(j, i, f"{value:.2f}", ha='center', va='center', color=text_color, fontsize=dyn_size_ch*1.2)

    # Custom colored legend (ifpTM, cpTM, ipTM)
    x_start = 0.35
    x_offset = 0.125
    dyn_size = 16
    ax_in.text(0.1, 1.05, prefix, fontsize=dyn_size, fontweight='bold', color='black', ha='center', transform=ax_in.transAxes)
    ax_in.text(x_start + x_offset - 0.06, 1.05, 'actifpTM', fontsize=dyn_size, fontweight='bold', color='darkblue', ha='center', transform=ax_in.transAxes)
    ax_in.text(x_start + 2 * x_offset, 1.05, ' - ', fontsize=dyn_size, fontweight='bold', color='black', ha='center', transform=ax_in.transAxes)
    ax_in.text(x_start + 3 * x_offset, 1.05, 'cpTM', fontsize=dyn_size, fontweight='bold', color='dimgrey', ha='center', transform=ax_in.transAxes)
    ax_in.text(x_start + 4 * x_offset, 1.05, ' - ', fontsize=dyn_size, fontweight='bold', color='black', ha='center', transform=ax_in.transAxes)
    ax_in.text(x_start + 5 * x_offset, 1.05, 'ipTM', fontsize=dyn_size, fontweight='bold', color='firebrick', ha='center', transform=ax_in.transAxes)

    # Format labels
    ax_in.set_yticks(np.arange(len(letters)))
    ax_in.set_yticklabels(letters, rotation=0, fontsize=dyn_size_ch*1.5)
    ax_in.set_xticks(np.arange(len(letters)))
    ax_in.set_xticklabels(letters, fontsize=dyn_size_ch*1.5)

    # If this was only one plot, display and save it.
    # If multiple plots have been appended, this needs to be done from the calling function
    if not ax_in:
        plt.tight_layout()
        plt.savefig(fig_path, dpi=200, bbox_inches='tight')

def plot_chain_pairwise_analysis(info, prefix='rank_', fig_path="chain_pairwise_ptm.png"):
    num_elements = len(info)
    fig, axes = plt.subplots(1, num_elements, figsize=(num_elements * 5, 5), squeeze=False)
    
    # We do this so that the same function can be used both for colabfold and colabdesign
    for idx, ax in enumerate(axes[0]):
        prefix_plot = prefix + str(idx+1)
        if isinstance(info[idx], dict):
            actifptm_dict = info[idx].get("pairwise_actifptm", {})
            iptm_dict = info[idx].get("pairwise_iptm", {})
            cptm_dict = info[idx].get("per_chain_ptm", {})
        else: 
            actifptm_dict = info[idx][4]
            iptm_dict = info[idx][3]
            cptm_dict = info[idx][5]

        plot_matrix(actifptm_dict, iptm_dict, cptm_dict, prefix=prefix_plot, ax_in=ax, fig_path=fig_path)

    plt.tight_layout()
    plt.savefig(fig_path, dpi=200, bbox_inches='tight')
    plt.close(fig)
