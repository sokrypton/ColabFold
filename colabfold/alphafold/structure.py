"""
Re-exports of ``alphafold.common`` protein I/O and residue constants.
"""
from alphafold.common import protein, residue_constants

from_prediction = protein.from_prediction
to_pdb = protein.to_pdb
from_pdb_string = protein.from_pdb_string
from_mmcif_string = protein.from_mmcif_string
PDB_CHAIN_IDS = protein.PDB_CHAIN_IDS


def amber_relaxation(**kwargs):
    from alphafold.relax import relax
    return relax.AmberRelaxation(**kwargs)
