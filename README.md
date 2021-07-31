# ColabFold
Making Protein folding accessible to all via Google Colab!
- [AlphaFold2_mmseqs2](https://colab.research.google.com/github/sokrypton/ColabFold/blob/main/AlphaFold2.ipynb) - `monomers=Yes, homoligomers=Yes, mmseqs2=Yes, jackhmmer=No, templates=Yes`
- [AlphaFold2_complexes](https://colab.research.google.com/github/sokrypton/ColabFold/blob/main/AlphaFold2_complexes.ipynb) - `monomers=No, heterodimers=Yes, mmseqs2=Yes, jackhmmer=No, templates=No`
- [AlphaFold2_jackhmmer](https://colab.research.google.com/github/sokrypton/ColabFold/blob/main/beta/AlphaFold_wJackhmmer.ipynb) - `monomers=Yes, homoligomers=Yes, mmseqs2=Yes, jackhmmer=Yes, templates=No`

Official Notebook from Deepmind:
- [AlphaFold2](https://colab.research.google.com/github/deepmind/alphafold/blob/main/notebooks/AlphaFold.ipynb) - `monomers=Yes, homoligomers=No, mmseqs2=No, jackhmmer=Yes, templates=No`

Experimental notebooks (WIP):
- [RoseTTAFold](https://colab.research.google.com/github/sokrypton/ColabFold/blob/main/RoseTTAFold.ipynb) - `Work in Progress`
- [AlphaFold2_noTemplates_noMD](https://colab.research.google.com/github/sokrypton/ColabFold/blob/main/verbose/alphafold_noTemplates_noMD.ipynb)
- [AlphaFold2_noTemplates_yesMD](https://colab.research.google.com/github/sokrypton/ColabFold/blob/main/verbose/alphafold_noTemplates_yesMD.ipynb)

Others
- [TrDesign](https://colab.research.google.com/github/gjoni/trDesign/blob/beta/02-GD/notebooks/TrDesign_GD_demo.ipynb) - TrRosetta_v1 for Protein Structure Prediction and Design

Maximum length limits depends on free GPU provided by Google-Colab `fingers-crossed`
- For GPU: `Tesla T4` or `Tesla P100` with ~16G the max length is ~1400
- For GPU: `Tesla K80` with ~12G the max length is ~1000
- To check what GPU you got, open a new code cell and type `!nvidia-smi`

Acknowledgments
- We would like to thank the RoseTTAFold and AlphaFold team for doing an excellent job open sourcing the software. 
- Also credit to [David Koes](https://github.com/dkoes) for his awesome [py3Dmol](https://3dmol.csb.pitt.edu/) plugin, without whom these notebooks would be quite boring!
- A colab by Sergey Ovchinnikov (@sokrypton), Milot Mirdita (@milot_mirdita) and Martin Steinegger (@thesteinegger).


How do I reference this work?

[![DOI](https://zenodo.org/badge/387617756.svg)](https://zenodo.org/badge/latestdoi/387617756)
