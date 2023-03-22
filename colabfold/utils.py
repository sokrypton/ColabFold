import json
import logging
import warnings
from pathlib import Path
from typing import Optional

from absl import logging as absl_logging
from importlib_metadata import distribution
from tqdm import TqdmExperimentalWarning

NO_GPU_FOUND = """ERROR: Jax could not find GPU. This can be either because your machine doesn't have a GPU
or because jax can't find it. You might need to run

pip install --upgrade "jax[cuda]" -f https://storage.googleapis.com/jax-releases/jax_releases.html  # Note: wheels only available on linux.

See https://github.com/google/jax/#pip-installation-gpu-cuda for more details.

If you're sure you want to run without a GPU, pass `--cpu`"""

DEFAULT_API_SERVER = "https://api.colabfold.com"

ACCEPT_DEFAULT_TERMS = \
"""
WARNING: You are welcome to use the default MSA server, however keep in mind that it's a
limited shared resource only capable of processing a few thousand MSAs per day. Please
submit jobs only from a single IP address. We reserve the right to limit access to the
server case-by-case when usage exceeds fair use. If you require more MSAs: You can 
precompute all MSAs with `colabfold_search` or host your own API and pass it to `--host-url`
"""

class TqdmHandler(logging.StreamHandler):
    """https://stackoverflow.com/a/38895482/3549270"""

    def __init__(self):
        logging.StreamHandler.__init__(self)

    def emit(self, record):
        # We need the native tqdm here
        from tqdm import tqdm

        msg = self.format(record)
        tqdm.write(msg)


def setup_logging(log_file: Path):
    log_file.parent.mkdir(exist_ok=True, parents=True)
    root = logging.getLogger()
    if root.handlers:
        for handler in root.handlers:
            root.removeHandler(handler)
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(message)s",
        handlers=[TqdmHandler(), logging.FileHandler(log_file)],
    )
    # otherwise jax will tell us about its search for devices
    absl_logging.set_verbosity("error")
    warnings.simplefilter(action="ignore", category=TqdmExperimentalWarning)


def safe_filename(file: str) -> str:
    return "".join([c if c.isalnum() or c in ["_", ".", "-"] else "_" for c in file])


def get_commit() -> Optional[str]:
    text = distribution("colabfold").read_text("direct_url.json")
    if not text:
        return None
    direct_url = json.loads(text)
    if "vcs_info" not in direct_url:
        return None
    if "commit_id" not in direct_url["vcs_info"]:
        return None
    return direct_url["vcs_info"]["commit_id"]


# Copied from Bio.PDB to override _save_dict method
# https://github.com/biopython/biopython/blob/biopython-179/Bio/PDB/mmcifio.py
# We add poly_seq and revision_date so that AF2 can read these cif files
# Original license BSD 3-clause

import re
from Bio.PDB import MMCIFIO
from Bio.PDB.Polypeptide import standard_aa_names

CIF_REVISION_DATE = """loop_
_pdbx_audit_revision_history.ordinal
_pdbx_audit_revision_history.data_content_type
_pdbx_audit_revision_history.major_revision
_pdbx_audit_revision_history.minor_revision
_pdbx_audit_revision_history.revision_date
1 'Structure model' 1 0 1971-01-01
#\n"""

### begin section copied from Bio.PDB
mmcif_order = {
    "_atom_site": [
        "group_PDB",
        "id",
        "type_symbol",
        "label_atom_id",
        "label_alt_id",
        "label_comp_id",
        "label_asym_id",
        "label_entity_id",
        "label_seq_id",
        "pdbx_PDB_ins_code",
        "Cartn_x",
        "Cartn_y",
        "Cartn_z",
        "occupancy",
        "B_iso_or_equiv",
        "pdbx_formal_charge",
        "auth_seq_id",
        "auth_comp_id",
        "auth_asym_id",
        "auth_atom_id",
        "pdbx_PDB_model_num",
    ]
}


class CFMMCIFIO(MMCIFIO):
    def _save_dict(self, out_file):
        # Form dictionary where key is first part of mmCIF key and value is list
        # of corresponding second parts
        key_lists = {}
        for key in self.dic:
            if key == "data_":
                data_val = self.dic[key]
            else:
                s = re.split(r"\.", key)
                if len(s) == 2:
                    if s[0] in key_lists:
                        key_lists[s[0]].append(s[1])
                    else:
                        key_lists[s[0]] = [s[1]]
                else:
                    raise ValueError("Invalid key in mmCIF dictionary: " + key)

        # Re-order lists if an order has been specified
        # Not all elements from the specified order are necessarily present
        for key, key_list in key_lists.items():
            if key in mmcif_order:
                inds = []
                for i in key_list:
                    try:
                        inds.append(mmcif_order[key].index(i))
                    # Unrecognised key - add at end
                    except ValueError:
                        inds.append(len(mmcif_order[key]))
                key_lists[key] = [k for _, k in sorted(zip(inds, key_list))]

        # Write out top data_ line
        if data_val:
            out_file.write("data_" + data_val + "\n#\n")
            ### end section copied from Bio.PDB
            # Add poly_seq as default MMCIFIO doesn't handle this
            out_file.write(
                """loop_
_entity_poly_seq.entity_id
_entity_poly_seq.num
_entity_poly_seq.mon_id
_entity_poly_seq.hetero
#\n"""
            )
            poly_seq = []
            chain_idx = 1
            for model in self.structure:
                for chain in model:
                    res_idx = 1
                    for residue in chain:
                        hetatm, _, _ = residue.get_id()
                        if hetatm != " ":
                            continue
                        poly_seq.append(
                            (chain_idx, res_idx, residue.get_resname(), "n")
                        )
                        res_idx += 1
                    chain_idx += 1
            for seq in poly_seq:
                out_file.write(f"{seq[0]} {seq[1]} {seq[2]}  {seq[3]}\n")
            out_file.write("#\n")
            out_file.write(
                """loop_
_chem_comp.id
_chem_comp.type
#\n"""
            )
            for three in standard_aa_names:
                out_file.write(f'{three} "peptide linking"\n')
            out_file.write("#\n")
            out_file.write(
                """loop_
_struct_asym.id
_struct_asym.entity_id
#\n"""
            )
            chain_idx = 1
            for model in self.structure:
                for chain in model:
                    out_file.write(f"{chain.get_id()} {chain_idx}\n")
                    chain_idx += 1
            out_file.write("#\n")

        ### begin section copied from Bio.PDB
        for key, key_list in key_lists.items():
            # Pick a sample mmCIF value, which can be a list or a single value
            sample_val = self.dic[key + "." + key_list[0]]
            n_vals = len(sample_val)
            # Check the mmCIF dictionary has consistent list sizes
            for i in key_list:
                val = self.dic[key + "." + i]
                if (
                    isinstance(sample_val, list)
                    and (isinstance(val, str) or len(val) != n_vals)
                ) or (isinstance(sample_val, str) and isinstance(val, list)):
                    raise ValueError(
                        "Inconsistent list sizes in mmCIF dictionary: " + key + "." + i
                    )
            # If the value is a single value, write as key-value pairs
            if isinstance(sample_val, str) or (
                isinstance(sample_val, list) and len(sample_val) == 1
            ):
                m = 0
                # Find the maximum key length
                for i in key_list:
                    if len(i) > m:
                        m = len(i)
                for i in key_list:
                    # If the value is a single item list, just take the value
                    if isinstance(sample_val, str):
                        value_no_list = self.dic[key + "." + i]
                    else:
                        value_no_list = self.dic[key + "." + i][0]
                    out_file.write(
                        "{k: <{width}}".format(k=key + "." + i, width=len(key) + m + 4)
                        + self._format_mmcif_col(value_no_list, len(value_no_list))
                        + "\n"
                    )
            # If the value is more than one value, write as keys then a value table
            elif isinstance(sample_val, list):
                out_file.write("loop_\n")
                col_widths = {}
                # Write keys and find max widths for each set of values
                for i in key_list:
                    out_file.write(key + "." + i + "\n")
                    col_widths[i] = 0
                    for val in self.dic[key + "." + i]:
                        len_val = len(val)
                        # If the value requires quoting it will add 2 characters
                        if self._requires_quote(val) and not self._requires_newline(
                            val
                        ):
                            len_val += 2
                        if len_val > col_widths[i]:
                            col_widths[i] = len_val
                # Technically the max of the sum of the column widths is 2048

                # Write the values as rows
                for i in range(n_vals):
                    for col in key_list:
                        out_file.write(
                            self._format_mmcif_col(
                                self.dic[key + "." + col][i], col_widths[col] + 1
                            )
                        )
                    out_file.write("\n")
            else:
                raise ValueError(
                    "Invalid type in mmCIF dictionary: " + str(type(sample_val))
                )
            out_file.write("#\n")
            ### end section copied from Bio.PDB
            out_file.write(CIF_REVISION_DATE)
