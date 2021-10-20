from string import ascii_uppercase
from typing import Tuple, Any

import numpy as np
from alphafold.common import protein
from alphafold.data import pipeline
from alphafold.data import templates
from alphafold.data.tools import hhsearch
from matplotlib import pyplot as plt

from colabfold.colabfold import run_mmseqs2
from colabfold.pdb import set_bfactor
from colabfold.utils import DEFAULT_API_SERVER


def mk_mock_template(query_sequence):
    # since alphafold's model requires a template input
    # we create a blank example w/ zero input, confidence -1
    ln = len(query_sequence)
    output_templates_sequence = "-" * ln
    output_confidence_scores = np.full(ln, -1)
    templates_all_atom_positions = np.zeros(
        (ln, templates.residue_constants.atom_type_num, 3)
    )
    templates_all_atom_masks = np.zeros((ln, templates.residue_constants.atom_type_num))
    templates_aatype = templates.residue_constants.sequence_to_onehot(
        output_templates_sequence, templates.residue_constants.HHBLITS_AA_TO_ID
    )
    template_features = {
        "template_all_atom_positions": templates_all_atom_positions[None],
        "template_all_atom_masks": templates_all_atom_masks[None],
        "template_sequence": [f"none".encode()],
        "template_aatype": np.array(templates_aatype)[None],
        "template_confidence_scores": output_confidence_scores[None],
        "template_domain_names": [f"none".encode()],
        "template_release_date": [f"none".encode()],
    }
    return template_features


def mk_template(query_sequence, a3m_lines, template_paths):
    template_featurizer = templates.TemplateHitFeaturizer(
        mmcif_dir=template_paths,
        max_template_date="2100-01-01",
        max_hits=20,
        kalign_binary_path="kalign",
        release_dates_path=None,
        obsolete_pdbs_path=None,
    )

    hhsearch_pdb70_runner = hhsearch.HHSearch(
        binary_path="hhsearch", databases=[f"{template_paths}/pdb70"]
    )

    hhsearch_result = hhsearch_pdb70_runner.query(a3m_lines)
    hhsearch_hits = pipeline.parsers.parse_hhr(hhsearch_result)
    templates_result = template_featurizer.get_templates(
        query_sequence=query_sequence,
        query_pdb_code=None,
        query_release_date=None,
        hits=hhsearch_hits,
    )
    return templates_result.features


def predict_structure(
    prefix, feature_dict, Ls, model_runner_and_params, do_relax=False, random_seed=0
):
    """Predicts structure using AlphaFold for the given sequence."""

    # Minkyung's code
    # add big enough number to residue index to indicate chain breaks
    idx_res = feature_dict["residue_index"]
    L_prev = 0
    # Ls: number of residues in each chain
    for L_i in Ls[:-1]:
        idx_res[L_prev + L_i :] += 200
        L_prev += L_i
    chains = list("".join([ascii_uppercase[n] * L for n, L in enumerate(Ls)]))
    feature_dict["residue_index"] = idx_res

    # Run the models.
    plddts, paes = [], []
    unrelaxed_pdb_lines = []
    relaxed_pdb_lines = []

    for model_name, (model_runner, params) in model_runner_and_params.items():
        print(f"running {model_name}")
        # swap params to avoid recompiling
        # note: models 1,2 have diff number of params compared to models 3,4,5 (this was handled on construction)
        model_runner.params = params

        processed_feature_dict = model_runner.process_features(
            feature_dict, random_seed=random_seed
        )
        prediction_result = model_runner.predict(processed_feature_dict)
        unrelaxed_protein = protein.from_prediction(
            processed_feature_dict, prediction_result
        )
        unrelaxed_pdb_lines.append(protein.to_pdb(unrelaxed_protein))
        plddts.append(prediction_result["plddt"])
        paes.append(prediction_result["predicted_aligned_error"])

        if do_relax:
            from alphafold.relax import relax
            from alphafold.common import residue_constants

            residue_constants.stereo_chemical_props_path = "stereo_chemical_props.txt"

            # Relax the prediction.
            amber_relaxer = relax.AmberRelaxation(
                max_iterations=0,
                tolerance=2.39,
                stiffness=10.0,
                exclude_residues=[],
                max_outer_iterations=20,
            )
            relaxed_pdb_str, _, _ = amber_relaxer.process(prot=unrelaxed_protein)
            relaxed_pdb_lines.append(relaxed_pdb_str)

    # rerank models based on predicted lddt
    lddt_rank = np.mean(plddts, -1).argsort()[::-1]
    out = {}
    print("reranking models based on avg. predicted lDDT")
    for n, r in enumerate(lddt_rank):
        print(f"model_{n + 1} {np.mean(plddts[r])}")

        unrelaxed_pdb_path = f"{prefix}_unrelaxed_model_{n + 1}.pdb"
        with open(unrelaxed_pdb_path, "w") as f:
            f.write(unrelaxed_pdb_lines[r])
        set_bfactor(unrelaxed_pdb_path, plddts[r], idx_res, chains)

        if do_relax:
            relaxed_pdb_path = f"{prefix}_relaxed_model_{n + 1}.pdb"
            with open(relaxed_pdb_path, "w") as f:
                f.write(relaxed_pdb_lines[r])
            set_bfactor(relaxed_pdb_path, plddts[r], idx_res, chains)

        out[f"model_{n + 1}"] = {"plddt": plddts[r], "pae": paes[r]}
    return out


def handle_homooligomer(deletion_matrix, homooligomer, msa, query_sequence):
    if homooligomer == 1:
        msas = [msa]
        deletion_matrices = [deletion_matrix]
    else:
        # make multiple copies of msa for each copy
        # AAA------
        # ---AAA---
        # ------AAA
        #
        # note: if you concat the sequences (as below), it does NOT work
        # AAAAAAAAA
        msas = []
        deletion_matrices = []
        Ln = len(query_sequence)
        for o in range(homooligomer):
            L = Ln * o
            R = Ln * (homooligomer - (o + 1))
            msas.append(["-" * L + seq + "-" * R for seq in msa])
            deletion_matrices.append(
                [[0] * L + mtx + [0] * R for mtx in deletion_matrix]
            )
    return deletion_matrices, msas


def get_msa(
    query_sequence: str,
    jobname: str,
    homooligomer: int,
    use_templates: bool,
    use_msa: bool,
    use_env: bool,
    a3m_file: str,
    host_url: str = DEFAULT_API_SERVER,
) -> Tuple[Any, Any, Any]:
    if use_templates:
        a3m_lines, template_paths = run_mmseqs2(
            query_sequence, jobname, use_env, use_templates=True, host_url=host_url
        )
        if template_paths is None:
            template_features = mk_mock_template(query_sequence * homooligomer)
        else:
            template_features = mk_template(query_sequence, a3m_lines, template_paths)
    elif use_msa:
        a3m_lines = run_mmseqs2(query_sequence, jobname, use_env, host_url=host_url)
        with open(a3m_file, "w") as text_file:
            text_file.write(a3m_lines)
        template_features = mk_mock_template(query_sequence * homooligomer)
    else:
        template_features = mk_mock_template(query_sequence * homooligomer)
        a3m_lines = "".join(open(a3m_file, "r").read())

    # parse MSA
    msa, deletion_matrix = pipeline.parsers.parse_a3m(a3m_lines)

    return msa, deletion_matrix, template_features


def plot_plddt_legend():
    thresh = [
        "plDDT:",
        "Very low (<50)",
        "Low (60)",
        "OK (70)",
        "Confident (80)",
        "Very high (>90)",
    ]
    plt.figure(figsize=(1, 0.1), dpi=100)
    for c in ["#FFFFFF", "#FF0000", "#FFFF00", "#00FF00", "#00FFFF", "#0000FF"]:
        plt.bar(0, 0, color=c)
    plt.legend(
        thresh,
        frameon=False,
        loc="center",
        ncol=6,
        handletextpad=1,
        columnspacing=1,
        markerscale=0.5,
    )
    plt.axis(False)
    return plt


def plot_confidence(outs, homooligomer: int, query_sequence: str, model_num=1):
    model_name = f"model_{model_num}"
    plt.figure(figsize=(10, 3), dpi=100)
    """Plots the legend for plDDT."""

    plt.subplot(1, 2, 1)
    plt.title("Predicted lDDT")
    plt.plot(outs[model_name]["plddt"])
    for n in range(homooligomer + 1):
        x = n * (len(query_sequence))
        plt.plot([x, x], [0, 100], color="black")
    plt.ylabel("plDDT")
    plt.xlabel("position")

    plt.subplot(1, 2, 2)
    plt.title("Predicted Aligned Error")
    plt.imshow(outs[model_name]["pae"], cmap="bwr", vmin=0, vmax=30)
    plt.colorbar()
    plt.xlabel("Scored residue")
    plt.ylabel("Aligned residue")

    return plt
