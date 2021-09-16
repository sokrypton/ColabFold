from pathlib import Path
from typing import Union

import numpy as np


def set_bfactor(pdb_filename: Union[str, Path], bfac, idx_res, chains):
    in_file = open(pdb_filename, "r").readlines()
    out_file = open(pdb_filename, "w")
    for line in in_file:
        if line[0:6] == "ATOM  ":
            seq_id = int(line[22:26].strip()) - 1
            seq_id = np.where(idx_res == seq_id)[0][0]
            out_file.write(
                f"{line[:21]}{chains[seq_id]}{line[22:60]}{bfac[seq_id]:6.2f}{line[66:]}"
            )
    out_file.close()


def show_pdb(
    use_amber: bool,
    jobname: str,
    homooligomer,
    model_num=1,
    show_sidechains=False,
    show_mainchains=False,
    color="lDDT",
):
    import py3Dmol

    model_name = f"model_{model_num}"
    if use_amber:
        pdb_filename = f"{jobname}_relaxed_{model_name}.pdb"
    else:
        pdb_filename = f"{jobname}_unrelaxed_{model_name}.pdb"

    view = py3Dmol.view(js="https://3dmol.org/build/3Dmol.js")
    view.addModel(open(pdb_filename, "r").read(), "pdb")

    if color == "lDDT":
        view.setStyle(
            {
                "cartoon": {
                    "colorscheme": {
                        "prop": "b",
                        "gradient": "roygb",
                        "min": 50,
                        "max": 90,
                    }
                }
            }
        )
    elif color == "rainbow":
        view.setStyle({"cartoon": {"color": "spectrum"}})
    elif color == "chain":
        for n, chain, color in zip(
            range(homooligomer),
            list("ABCDEFGH"),
            ["lime", "cyan", "magenta", "yellow", "salmon", "white", "blue", "orange"],
        ):
            view.setStyle({"chain": chain}, {"cartoon": {"color": color}})
    if show_sidechains:
        BB = ["C", "O", "N"]
        view.addStyle(
            {
                "and": [
                    {"resn": ["GLY", "PRO"], "invert": True},
                    {"atom": BB, "invert": True},
                ]
            },
            {"stick": {"colorscheme": f"WhiteCarbon", "radius": 0.3}},
        )
        view.addStyle(
            {"and": [{"resn": "GLY"}, {"atom": "CA"}]},
            {"sphere": {"colorscheme": f"WhiteCarbon", "radius": 0.3}},
        )
        view.addStyle(
            {"and": [{"resn": "PRO"}, {"atom": ["C", "O"], "invert": True}]},
            {"stick": {"colorscheme": f"WhiteCarbon", "radius": 0.3}},
        )
    if show_mainchains:
        BB = ["C", "O", "N", "CA"]
        view.addStyle(
            {"atom": BB}, {"stick": {"colorscheme": f"WhiteCarbon", "radius": 0.3}}
        )

    view.zoomTo()
    return view
