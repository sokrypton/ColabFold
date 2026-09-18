"""Interface confidence scores."""

# Reimplemented from https://github.com/DunbrackLab/IPSAE (MIT).

import numpy as np

from alphafold.common import residue_constants
from alphafold.common.protein import PDB_CHAIN_IDS

# Angstrom.
DEFAULT_PAE_CUTOFF = 15.0
PDOCKQ_DIST_CUTOFF = 8.0


def _ptm_score(pae, d0):
    """TM-score kernel."""
    return 1.0 / (1.0 + (pae / d0) ** 2.0)


def _calc_d0_array(n_res):
    """TM-score d0."""
    length = np.maximum(26.0, np.asarray(n_res, dtype=np.float64))
    return np.maximum(1.0, 1.24 * (length - 15.0) ** (1.0 / 3.0) - 1.8)


def get_interface_scores(pae, plddt, asym_id, atom_positions, atom_mask,
                         pae_cutoff=DEFAULT_PAE_CUTOFF):
    """Return ipSAE, pDockQ and pDockQ2 per chain pair."""
    asym_id = np.asarray(asym_id)
    unique_chains = np.unique(asym_id)
    if len(unique_chains) < 2:
        return {}

    pae = np.asarray(pae, dtype=np.float64)
    plddt = np.asarray(plddt, dtype=np.float64)
    atom_positions = np.asarray(atom_positions, dtype=np.float64)

    # Use CA for glycine.
    atom_mask = np.asarray(atom_mask)
    ca_idx = residue_constants.atom_order["CA"]
    cb_idx = residue_constants.atom_order["CB"]
    has_cb = atom_mask[:, cb_idx] > 0.5
    cb_coords = np.where(has_cb[:, None],
                         atom_positions[:, cb_idx], atom_positions[:, ca_idx])
    distances = np.sqrt(((cb_coords[:, None] - cb_coords[None, :]) ** 2).sum(-1))
    # Exclude unplaced residues.
    unplaced = atom_mask[:, ca_idx] <= 0.5
    distances[unplaced, :] = np.inf
    distances[:, unplaced] = np.inf

    chain_label = {chain: PDB_CHAIN_IDS[int(chain)] for chain in unique_chains}
    ipsae, pdockq, pdockq2 = {}, {}, {}
    for chain_1 in unique_chains:
        for chain_2 in unique_chains:
            if chain_1 == chain_2:
                continue
            key = f"{chain_label[chain_1]}-{chain_label[chain_2]}"
            pair_mask = np.outer(asym_id == chain_1, asym_id == chain_2)

            # ipSAE.
            valid_pairs = pair_mask & (pae < pae_cutoff)
            n0_res = valid_pairs.sum(axis=1)
            d0_res = _calc_d0_array(n0_res)
            ptm_sums = (_ptm_score(pae, d0_res[:, None]) * valid_pairs).sum(axis=1)
            ipsae_by_res = np.divide(ptm_sums, n0_res,
                                     out=np.zeros_like(ptm_sums), where=n0_res > 0)
            ipsae[key] = round(float(ipsae_by_res.max()), 6)

            # pDockQ and pDockQ2.
            contacts = pair_mask & (distances <= PDOCKQ_DIST_CUTOFF)
            n_contacts = contacts.sum()
            if n_contacts > 0:
                interface_res = contacts.any(axis=1) | contacts.any(axis=0)
                mean_plddt = plddt[interface_res].mean()
                pdockq_val = 0.724 / (1.0 + np.exp(
                    -0.052 * (mean_plddt * np.log10(n_contacts) - 152.611))) + 0.018
                mean_ptm = _ptm_score(pae[contacts], 10.0).mean()
                pdockq2_val = 1.31 / (1.0 + np.exp(
                    -0.075 * (mean_plddt * mean_ptm - 84.733))) + 0.005
            else:
                pdockq_val, pdockq2_val = 0.0, 0.0
            if chain_1 < chain_2:  # Symmetric.
                pdockq[key] = round(float(pdockq_val), 4)
            pdockq2[key] = round(float(pdockq2_val), 4)

    return {"ipsae": ipsae, "pdockq": pdockq, "pdockq2": pdockq2}


def format_ipsae(ipsae_scores):
    """Format per-interface maxima."""
    pair_max = {}
    for key, value in ipsae_scores.items():
        chain_1, chain_2 = key.split("-")
        pair = key if chain_1 < chain_2 else f"{chain_2}-{chain_1}"
        pair_max[pair] = max(pair_max.get(pair, 0.0), value)
    if len(pair_max) == 1:
        return f"{next(iter(pair_max.values())):.3g}"
    return ",".join(f"{pair}:{value:.3g}" for pair, value in sorted(pair_max.items()))
