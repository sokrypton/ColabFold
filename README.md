# ColabFold
<img src="https://user-images.githubusercontent.com/4187522/128616692-a5f8ba4a-4f08-44ff-9bfd-a931bd8329c2.png" width="256" height="215">

Making Protein folding accessible to all via Google Colab!
- [AlphaFold2_mmseqs2](https://colab.research.google.com/github/sokrypton/ColabFold/blob/main/AlphaFold2.ipynb) - `monomers=Yes, homoligomers=Yes, mmseqs2=Yes, jackhmmer=No, templates=Yes`
- [AlphaFold2_advanced](https://colab.research.google.com/github/sokrypton/ColabFold/blob/main/beta/AlphaFold2_advanced.ipynb) - `monomers=Yes, complexes=Yes, homoligomers=Yes, mmseqs2=Yes, jackhmmer=Yes, templates=No`
- [RoseTTAFold](https://colab.research.google.com/github/sokrypton/ColabFold/blob/main/RoseTTAFold.ipynb) - `monomers=Yes, homoligomers=No, mmseqs2=Yes, jackhmmer=No, templates=No`

Official Notebook from Deepmind:
- [AlphaFold2](https://colab.research.google.com/github/deepmind/alphafold/blob/main/notebooks/AlphaFold.ipynb) - `monomers=Yes, homoligomers=No, mmseqs2=No, jackhmmer=Yes, templates=No`

OLD Experimental notebooks:
- [AlphaFold2_complexes](https://colab.research.google.com/github/sokrypton/ColabFold/blob/main/AlphaFold2_complexes.ipynb) - `monomers=No, complexes=Yes, mmseqs2=Yes, jackhmmer=No, templates=No`
- [AlphaFold2_jackhmmer](https://colab.research.google.com/github/sokrypton/ColabFold/blob/main/beta/AlphaFold_wJackhmmer.ipynb) - `monomers=Yes, homoligomers=Yes, mmseqs2=Yes, jackhmmer=Yes, templates=No`
- [AlphaFold2_noTemplates_noMD](https://colab.research.google.com/github/sokrypton/ColabFold/blob/main/verbose/alphafold_noTemplates_noMD.ipynb)
- [AlphaFold2_noTemplates_yesMD](https://colab.research.google.com/github/sokrypton/ColabFold/blob/main/verbose/alphafold_noTemplates_yesMD.ipynb)

FAQ
- Can I use the models for **Molecular Replacement**?
  - Yes, but be **CAREFUL**, the bfactor column is populated with pLDDT confidence values (higher = better). Phenix.phaser expects a "real" bfactor, where (lower = better). See [post](https://twitter.com/cheshireminima/status/1423929241675120643) from Claudia Millán.
- What is the maximum length?
  - Limits depends on free GPU provided by Google-Colab `fingers-crossed`
  - For GPU: `Tesla T4` or `Tesla P100` with ~16G the max length is ~1400
  - For GPU: `Tesla K80` with ~12G the max length is ~1000
  - To check what GPU you got, open a new code cell and type `!nvidia-smi`

Tutorials & Presentations
- ColabFold Tutorial presented at the Boston Protein Design and Modeling Club. [[video]](https://www.youtube.com/watch?v=Rfw7thgGTwI) [[slides]](https://docs.google.com/presentation/d/1mnffk23ev2QMDzGZ5w1skXEadTe54l8-Uei6ACce8eI/edit?usp=sharing). 

Acknowledgments
- We would like to thank the RoseTTAFold and AlphaFold team for doing an excellent job open sourcing the software. 
- Also credit to [David Koes](https://github.com/dkoes) for his awesome [py3Dmol](https://3dmol.csb.pitt.edu/) plugin, without whom these notebooks would be quite boring!
- A colab by Sergey Ovchinnikov (@sokrypton), Milot Mirdita (@milot_mirdita) and Martin Steinegger (@thesteinegger).


How do I reference this work?

[![DOI](https://zenodo.org/badge/387617756.svg)](https://zenodo.org/badge/latestdoi/387617756)
