# RGI in the standard Colab form

Open [ColabFold2 preview](https://colab.research.google.com/github/sokrypton/ColabFold/blob/main/ColabFold2_preview.ipynb#scrollTo=rgi-restraints).
**RGI: Making restraint-guided protein folding inference accessible to all! (optional)**
is a standard Colab form cell after **Input sequences**.
Its fields are visible and editable before execution.

1. Choose an RGI-supported model such as **openbind0** in **Install dependencies**.
2. Fill in **Input sequences**.
3. In the **RGI** form, turn on **use_rgi** and enter the selections and targets.
4. Use **Runtime → Run all**.

For a 25 Å distance, set `distance_atom_selection1` to `chain A and resid 1 to 10`,
`distance_atom_selection2` to `chain A and resid 40 to 50`, and `target_distance` to
`25`. Use chains and residues from your own input. Leave unused restraint fields empty.

The same form provides conformer, angle, custom and RMSD. Lists in the fields create
multiple restraints; `restraints_config` accepts additional native toolkit entries
or `{"config_path": "restraints.yaml"}`. See the
[full guide and examples](https://github.com/cddlab/rgi_toolkit/blob/main/docs/colabfold.md).

OpenBind-0 uses the same RGI sampler hook as OpenFold3 in the shared JAX runner,
while retaining its own model configuration and checkpoint. It is available in
both ColabFold2 preview and AlphaFold3 / OpenFold3.

Leave **use_rgi** off for vanilla. AF2 and IntelliFold2 are vanilla only.
Vanilla runs use the upstream predictor without the RGI wrapper or RGI dependency.
Any restraint configuration and chain opt-ins left from an earlier run are removed.
The notebooks retain their default models and sampling settings.
After edits, rerun the RGI cell and prediction, or use **Run all**. In Boltz-1 the form
precedes installation; rerun installation too when changing **use_rgi**.
RGI results include the RGI-Toolkit paper in `cite.bibtex`, which is also included
in the downloaded results ZIP. Existing citations are preserved without duplicating
the RGI entry on reruns.
