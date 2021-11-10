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
