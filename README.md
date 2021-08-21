# ColabFold
-----------------
**MMseqs2 API status**
```diff
  16Aug2021: WARNING - MMseqs2 API is undergoing upgrade, you may see error messages.
  17Aug2021: If you see any errors, please report them.
  17Aug2021: We are still debugging the MSA generation procedure...
  20Aug2021: WARNING - MMseqs2 API is undergoing upgrade, you may see error messages.
             To avoid Google Colab from crashing, for large MSA we did -diff 1000 to get 
             1K most diverse sequences. This caused some large MSA to degrade in quality,
             as sequences close to query were being merged to single representive.
             We are working on updating the server (today) to fix this, by making sure
             that both diverse and sequences close to query are included in the final MSA.
             We'll post update here when update is complete.
+ 21Aug2021  The MSA issues should now be resolved! Please report any errors you see.
+            In short, to reduce MSA size we filter (qsc > 0.8, id > 0.95) and take 3K
+            most diverse sequences at different qid (sequence identity to query) intervals 
+            and merge them. More specifically 3K sequences at qid at (0→0.2),(0.2→0.4),
+            (0.4→0.6),(0.6→0.8) and (0.8→1). If you submitted your sequence between
+            16Aug2021 and 20Aug2021, we recommend submitting again for best results!
- 21Aug2021  The use_templates option in AlphaFold2_mmseqs2 is not properly working. We are
-            working on fixing this. If you are not using templates, this does not affect the
-            the results. Other notebooks that do not use_templates are unaffected.
```
-----------------

<p align="center"><img src="https://github.com/sokrypton/ColabFold/raw/main/.github/ColabFold_Marv_Logo.png" height="250"/></p>

### Making Protein folding accessible to all via Google Colab!

| Notebooks | monomers | complexes | mmseqs2 | jackhmmer | templates   |
| :-------- | -------  | --------- | ------- | --------- | ----------- |
| [AlphaFold2_mmseqs2](https://colab.research.google.com/github/sokrypton/ColabFold/blob/main/AlphaFold2.ipynb) | Yes | No | Yes | No | Yes |
| [AlphaFold2_advanced](https://colab.research.google.com/github/sokrypton/ColabFold/blob/main/beta/AlphaFold2_advanced.ipynb) | Yes | Yes | Yes | Yes | No |
| [RoseTTAFold](https://colab.research.google.com/github/sokrypton/ColabFold/blob/main/RoseTTAFold.ipynb) | Yes | No | Yes | No | No |
| [AlphaFold2](https://colab.research.google.com/github/deepmind/alphafold/blob/main/notebooks/AlphaFold.ipynb) (from Deepmind) | Yes | No | No | Yes | No |
||
| **OLD retired notebooks** | **monomers** | **complexes** | **mmseqs2** | **jackhmmer** | **templates** |
| [AlphaFold2_complexes](https://colab.research.google.com/github/sokrypton/ColabFold/blob/main/AlphaFold2_complexes.ipynb) | No | Yes | No | No | No |
| [AlphaFold2_jackhmmer](https://colab.research.google.com/github/sokrypton/ColabFold/blob/main/beta/AlphaFold_wJackhmmer.ipynb) | Yes | No | Yes | Yes | No |
| [AlphaFold2_noTemplates_noMD](https://colab.research.google.com/github/sokrypton/ColabFold/blob/main/verbose/alphafold_noTemplates_noMD.ipynb) |
| [AlphaFold2_noTemplates_yesMD](https://colab.research.google.com/github/sokrypton/ColabFold/blob/main/verbose/alphafold_noTemplates_yesMD.ipynb) |

### FAQ
- Can I use the models for **Molecular Replacement**?
  - Yes, but be **CAREFUL**, the bfactor column is populated with pLDDT confidence values (higher = better). Phenix.phaser expects a "real" bfactor, where (lower = better). See [post](https://twitter.com/cheshireminima/status/1423929241675120643) from Claudia Millán.
- What is the maximum length?
  - Limits depends on free GPU provided by Google-Colab `fingers-crossed`
  - For GPU: `Tesla T4` or `Tesla P100` with ~16G the max length is ~1400
  - For GPU: `Tesla K80` with ~12G the max length is ~1000
  - To check what GPU you got, open a new code cell and type `!nvidia-smi`

### Tutorials & Presentations
- ColabFold Tutorial presented at the Boston Protein Design and Modeling Club. [[video]](https://www.youtube.com/watch?v=Rfw7thgGTwI) [[slides]](https://docs.google.com/presentation/d/1mnffk23ev2QMDzGZ5w1skXEadTe54l8-Uei6ACce8eI). 

### Acknowledgments
- We would like to thank the RoseTTAFold and AlphaFold team for doing an excellent job open sourcing the software. 
- Also credit to [David Koes](https://github.com/dkoes) for his awesome [py3Dmol](https://3dmol.csb.pitt.edu/) plugin, without whom these notebooks would be quite boring!
- A colab by Sergey Ovchinnikov (@sokrypton), Milot Mirdita (@milot_mirdita) and Martin Steinegger (@thesteinegger).


### How do I reference this work?

Mirdita M, Ovchinnikov S and Steinegger M. ColabFold - Making protein folding accessible to all. 
<br />
bioRxiv, doi: [10.1101/2021.08.15.456425](https://www.biorxiv.org/content/10.1101/2021.08.15.456425v1) (2021)

[![DOI](https://zenodo.org/badge/doi/10.5281/zenodo.5123296.svg)](https://doi.org/10.5281/zenodo.5123296)
