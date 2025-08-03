# ColabFold - v1.5.5

For details of what was changed in v1.5, see [change log](https://github.com/sokrypton/ColabFold/wiki/v1.5.0)!

<p align="center"><img src="https://github.com/sokrypton/ColabFold/raw/main/.github/ColabFold_Marv_Logo.png" height="250"/></p>

> [!NOTE]
> 04Aug2025: We changed the taxonomy/pairing files for the UniRef100 database. This might affect multimer predictions. Check [the wiki entry](https://github.com/sokrypton/ColabFold/wiki/MSA-Server-Database-History) for details. 

### Making Protein folding accessible to all via Google Colab!

| Notebooks                                                                                                                                        | monomers | complexes | mmseqs2 | jackhmmer | templates |
| :----------------------------------------------------------------------------------------------------------------------------------------------- | -------- | --------- | ------- | --------- | --------- |
| [AlphaFold2_mmseqs2](https://colab.research.google.com/github/sokrypton/ColabFold/blob/main/AlphaFold2.ipynb)                                    | Yes      | Yes       | Yes     | No        | Yes       |
| [AlphaFold2_batch](https://colab.research.google.com/github/sokrypton/ColabFold/blob/main/batch/AlphaFold2_batch.ipynb)                          | Yes      | Yes       | Yes     | No        | Yes       |
| [AlphaFold2](https://colab.research.google.com/github/deepmind/alphafold/blob/main/notebooks/AlphaFold.ipynb) (from Deepmind)                    | Yes      | Yes       | No      | Yes       | No        |
| [relax_amber](https://colab.research.google.com/github/sokrypton/ColabFold/blob/main/beta/relax_amber.ipynb) (relax input structure)             |          |           |         |           |           |
| [ESMFold](https://colab.research.google.com/github/sokrypton/ColabFold/blob/main/ESMFold.ipynb)                                                  | Yes      | Maybe     | No      | No        | No        |
|                                                                                                                                                  |
| **BETA (in development) notebooks**                                                                                                              |          |           |         |           |           |
| [RoseTTAFold2](https://colab.research.google.com/github/sokrypton/ColabFold/blob/main/RoseTTAFold2.ipynb)                                        | Yes      | Yes       | Yes     | No        | WIP       |
| [Boltz](https://colab.research.google.com/github/sokrypton/ColabFold/blob/main/Boltz1.ipynb)                                        | Yes      | Yes       | Yes     | No        | No       |
| [BioEmu](https://colab.research.google.com/github/sokrypton/ColabFold/blob/main/BioEmu.ipynb)                                        | Yes      | No       | Yes     | No        | No       |
| [OmegaFold](https://colab.research.google.com/github/sokrypton/ColabFold/blob/main/beta/omegafold.ipynb)                                         | Yes      | Maybe     | No      | No        | No        |
| [AlphaFold2_advanced_v2](https://colab.research.google.com/github/sokrypton/ColabDesign/blob/gamma/af/examples/predict.ipynb) (new experimental notebook)                  | Yes      | Yes       | Yes     | No        | Yes       |

Check the wiki page [old retired notebooks](https://github.com/sokrypton/ColabFold/wiki/Old-retired-notebooks) for unsupported notebooks.

### FAQ
- Where can I chat with other ColabFold users?
  - See our [Discord](https://discord.gg/gna8maru7d) channel!
- Can I use the models for **Molecular Replacement**?
  - Yes, but be **CAREFUL**, the bfactor column is populated with pLDDT confidence values (higher = better). Phenix.phaser expects a "real" bfactor, where (lower = better). See [post](https://twitter.com/cheshireminima/status/1423929241675120643) from Claudia Millán.
- What is the maximum length?
  - Limits depends on free GPU provided by Google-Colab `fingers-crossed`
  - For GPU: `Tesla T4` with ~16G the max length is ~2000
  - To check what GPU you got, open a new code cell and type `!nvidia-smi`
- Is it okay to use the MMseqs2 MSA server on a local computer?
  - You can access the server from a local computer if you queries are serial from a single IP. Please do not use multiple computers to query the server.
- Where can I download the databases used by ColabFold?
  - The databases are available at [colabfold.mmseqs.com](https://colabfold.mmseqs.com)
- I want to render my own images of the predicted structures, how do I color by pLDDT?
  - In pymol for AlphaFold structures: `spectrum b, red_yellow_green_cyan_blue, minimum=50, maximum=90`
  - If you want to use AlphaFold Colours (credit: Konstantin Korotkov)
    ```python
    set_color n0, [0.051, 0.341, 0.827]
    set_color n1, [0.416, 0.796, 0.945]
    set_color n2, [0.996, 0.851, 0.212]
    set_color n3, [0.992, 0.490, 0.302]
    color n0, b < 100; color n1, b < 90
    color n2, b < 70;  color n3, b < 50
    ```
  - In pymol for RoseTTAFold structures: `spectrum b, red_yellow_green_cyan_blue, minimum=0.5, maximum=0.9`
- What is the difference between the AlphaFold2_advanced and AlphaFold2_mmseqs2 (_batch) notebook for complex prediction?
  - We currently have two different ways to predict protein complexes: (1) using the AlphaFold2 model with residue index jump and (2) using the AlphaFold2-multimer model. AlphaFold2_advanced supports (1) and AlphaFold2_mmseqs2 (_batch) (2).
- What is the difference between localcolabfold and the pip installable colabfold_batch?
  -  [LocalColabFold](https://github.com/YoshitakaMo/localcolabfold) is an installer script designed to make ColabFold functionality available on local users' machines. It supports wide range of operating systems, such as Windows 10 or later (using Windows Subsystem for Linux 2), macOS, and Linux.
- Is there a way to amber-relax structures without having to rerun alphafold/colabfold from scratch?
  - Yes, see this [notebook](https://colab.research.google.com/github/sokrypton/ColabFold/blob/main/beta/relax_amber.ipynb).
- Where can I find the old notebooks that were previously developed and are now retired?
  - You can find the list of retired notebooks in the [old retired notebooks](https://github.com/sokrypton/ColabFold/wiki/Old-retired-notebooks) wiki page.
- Where can I find the history of MSA Server Databases used in ColabFold?
  - You can view the database version history on the [MSA Server Database History](https://github.com/sokrypton/ColabFold/wiki/MSA-Server-Database-History) wiki page.

### Running locally
For instructions on how to install ColabFold locally refer to [localcolabfold](https://github.com/YoshitakaMo/localcolabfold) or see our [wiki](https://github.com/sokrypton/ColabFold/wiki/Running-ColabFold-in-Docker) on how to run ColabFold within Docker.

### Generating MSAs for small scale local structure/complex predictions using the MSA server

When you pass a FASTA or CSV file containing your sequences to `colabfold_batch` it will automatically query the public MSA server to generate MSAs. You might want to split this into two steps for better GPU resource utilization:

```
# Query the MSA server and predict the structure on local GPU in one go:
colabfold_batch input_sequences.fasta out_dir

# Split querying MSA server and GPU predictions into two steps
colabfold_batch input_sequences.fasta out_dir --msa-only
colabfold_batch input_sequences.fasta out_dir
```

### Generating MSAs for large scale structure/complex predictions

First create a directory for the databases on a disk with sufficient storage (940GB (!)). Depending on where you are, this will take a couple of hours:

Note: [MMseqs2 Release 18](https://github.com/soedinglab/MMseqs2/releases/tag/18-8cc5c) is used to create the databases and perform sequece search in the ColabFold MSA server. Please use this version if you want to obtain the same MSAs as the server.

```shell
MMSEQS_NO_INDEX=1 ./setup_databases.sh /path/to/db_folder
```

If MMseqs2 is not installed in your `PATH`, add `--mmseqs <path to mmseqs>` to your `mmseqs` in `colabfold_search`:

```shell
# This needs a lot of CPU
colabfold_search --mmseqs /path/to/bin/mmseqs input_sequences.fasta /path/to/db_folder msas
# This needs a GPU
colabfold_batch msas predictions
```

This will create intermediate folder `msas` that contains all input multiple sequence alignments formated as a3m files and a `predictions` folder with all predicted pdb,json and png files.

The procedure above disables MMseqs2 preindexing of the various ColabFold databases by setting the `MMSEQS_NO_INDEX=1` environment variable before calling the database setup script. For most use-cases of `colabfold_search` precomputing the index is not required and might hurt search speed. The precomputed index is necessary for fast response times of the ColabFold server, where the whole database is permamently kept in memory. In any case the batch searches will require a machine with about 128GB RAM or, if the databases are to be kept permamently in RAM, with over 1TB RAM.

In some cases using precomputed database can still be useful. For the following cases, call the `setup_databases.sh` script without the `MMSEQS_NO_INDEX` environment variable:

(0) As mentioned above, if you want to set-up a server.

(1) If the precomputed index is stored on a very fast storage system (e.g., NVMe-SSDs) it might be faster to read the index from disk than computing in on the fly.  In this case, the search should be performed on the same machine that called `setup_databases.sh` since the precomputed index is created to fit within the given main memory size. Additionaly, pass the `--db-load-mode 0` option to make sure the database is read once from the storage system before use.

(2) Fast single query searches require the full index (the `.idx` files) to be kept in memory. This can be done with e.g. by using [vmtouch](https://github.com/hoytech/vmtouch). Thus, this type of search requires a machine with at least 768GB to 1TB RAM for the ColabfoldDB. If the index is present in memory, use the `--db-load-mode 2` parameter in `colabfold_search` to avoid index loading overhead.

If no index was created (`MMSEQS_NO_INDEX=1` was set), then `--db-load-mode` does not do anything and can be ignored.

### Saving MSAs in AlphaFold3-compatible JSON format
You can export MSAs into a json format compatible with AlphaFold3 input using the `--af3-json` option. 

**With colabfold_search:**

If you are using the local database setup with colabfold_search, you can add the `--af3-json` option to save the MSAs as AlphaFold3 input json:
```shell
colabfold_search --mmseqs /path/to/bin/mmseqs input_sequences.fasta /path/to/db_folder msas --af3-json
```
This will create a json file in the `msas` folder, using the same name as the a3m file.

**With colabfold_batch:**

If you are using the MSA server via colabfold_batch, you can also use the `--af3-json` option. However, structure prediction will be skipped, and only the json file will be generated. 
```shell
colabfold_batch input_sequences.fasta out_dir --af3-json
```

#### Including non-protein molecules in FASTA
AlphaFold3 supports non-protein components such as ligands and nucleic acids in input complexes. To include these in the generated json file, you can specify them directly in your FASTA input using the following format, `molecule type|sequence|(copies)`. As molecue types, dna, rna, ccd, smiles are allowed.

 > :exclamation: **Substitute aromatic bonds in SMILES**
 > If your SMILES string contains aromatic bonds (`:`), please replace them with semicolons (`;`) to avoid internal parsing issues.

- Examples
  - For DNA: `dna|ATCG`
  - For RNA: `rna|AUGC`
  - For ligands: 
    - SMILES string: `smiles|C1=NC(=C2C(=N1)N(C=N2)[C@H]3[C@@H]([C@@H]([C@H](O3)COP(=O)(O)OP(=O)(O)OP(=O)(O)O)O)O)N`
    - CCD code: `ccd|ATP`
  - To specify multiple copies of a molecule, you can add a number after the sequence, e.g. `ccd|ATP|2` or `dna|ATCG|2`.

Here is an example of biological complex with 2 proteins and 2 ATP ligands:
```
>Complex1|Prot1:Prot2:Lig
FIRSTPROTEIN:SECONDPROTEIN:ccd|ATP|2
>Complex2|Prot1:Prot2:Lig
FIRSTPROTEIN:SECONDPROTEIN:ccd|ATP:ccd|ATP
```
As the `copies` is optional, the `Complex1` and `Complex2` will result in identical json input.

 Note that MMseqs2-based MSAs are only generated for the protein sequences. RNA entries will not have unpaired MSAs in the json file. However, the field is marked as null so that AlphaFold3 can generate MSAs for them. 

### GPU-accelerated search with ⁠`colabfold_search` ⁠
ColabFold supports GPU-accelerated MSA searches through [MMseqs2-GPU](https://www.biorxiv.org/content/10.1101/2024.11.13.623350v1).

#### GPU database setup
To setup the GPU databases, you will need to run the ⁠`setup_databases.sh`⁠ command with ⁠`GPU=1`⁠ as an environment variable:

```
GPU=1 ./setup_databases.sh /path/to/db_folder
```

This will download and setup the GPU databases in the specified folder. Note that here we do not pass ⁠`MMSEQS_NO_INDEX=1`⁠ as an argument since the indices are useful in the GPU search since we will keep them in the GPU memory.

#### GPU search
By default, running `colabfold_search` with the `--gpu 1` option uses all available GPUs for its search.

```
colabfold_search /path/to/bin/mmseqs input_sequences.fasta /path/to/db_folder msas --gpu 1 
```

To select specific GPUs, set the `CUDA_VISIBLE_DEVICES` environment variable:
```
CUDA_VISIBLE_DEVICES=0,1 colabfold_search --mmseqs /path/to/bin/mmseqs input_sequences.fasta /path/to/db_folder msas --gpu 1
```

#### Optional GPU server for enhanced performance:
For frequent searches or to achieve minimal latency, you can run a dedicated GPU server. This server holds databases permanently in GPU memory, largely eliminating search overhead:

Start the GPU server(s):
```
mmseqs gpuserver /path/to/db_folder/colabfold_envdb_202108_db --max-seqs 10000 --db-load-mode 0 --prefilter-mode 1 &
PID1=$!
mmseqs gpuserver /path/to/db_folder/uniref30_2302_db --max-seqs 10000 --db-load-mode 0 --prefilter-mode 1 &
PID2=$!
```

By default, the GPU server distributes the database evenly across all visible GPUs. You can limit GPU usage by setting the CUDA_VISIBLE_DEVICES environment variable (e.g., `CUDA_VISIBLE_DEVICES=0,1`). 
Important: Ensure that the `CUDA_VISIBLE_DEVICES` environment variable is set consistently for both `gpuserver` and `colabfold_search`, otherwise `colabfold_search` will try wait for the `gpuserver` to appear until a set timeout (by default 5 minutes). If your database exceeds GPU memory capacity, the GPU server efficiently streams data between host and GPU memory using asynchronous CUDA streams.

Run searches using the GPU server:
```
colabfold_search --mmseqs /path/to/bin/mmseqs input_sequences.fasta /path/to/db_folder msas --gpu 1 --gpu-server 1
```
To stop the server(s) when done:
```
kill $PID1
kill $PID2
```
For more details, see [GPU-accelerated search](https://github.com/soedinglab/MMseqs2/wiki#gpu-accelerated-search).

### Tutorials & Presentations
- ColabFold Tutorial presented at the Boston Protein Design and Modeling Club. [[video]](https://www.youtube.com/watch?v=Rfw7thgGTwI) [[slides]](https://docs.google.com/presentation/d/1mnffk23ev2QMDzGZ5w1skXEadTe54l8-Uei6ACce8eI).

### Projects based on ColabFold or helpers

- [Run ColabFold on your local computer](https://github.com/YoshitakaMo/localcolabfold) by Yoshitaka Moriwaki
- [ColabFold/AlphaFold2 for protein structure predictions for Discoba species](https://github.com/zephyris/discoba_alphafold) by Richard John Wheeler
- [Cloud-based molecular simulations for everyone](https://github.com/pablo-arantes/Making-it-rain) by Pablo R. Arantes, Marcelo D. Polêto, Conrado Pedebos and Rodrigo Ligabue-Braun
- [getmoonbear is a webserver to predict protein structures](https://www.getmoonbear.com/AlphaFold2) by Stephanie Zhang and Neil Deshmukh
- [ColabFold/AlphaFold2 IDR complex prediction](https://github.com/normandavey/AlphaFold2-IDR-complex-prediction) by Balint Meszaros
- [ColabFold/AlphaFold2 (Phenix version) for macromolecular structure determination](https://colab.research.google.com/github/phenix-project/Colabs/blob/main/alphafold2/AlphaFold2.ipynb) by Tom Terwilliger
- [AlphaPickle: making AlphaFold2/ColabFold outputs interpretable](https://colab.research.google.com/github/mattarnoldbio/alphapickle/blob/main/AlphaPickle.ipynb) by Matt Arnold

### Acknowledgments
- We would like to thank the [RoseTTAFold](https://github.com/RosettaCommons/RoseTTAFold) and [AlphaFold](https://github.com/deepmind/alphafold) team for doing an excellent job open sourcing the software.
- Also credit to [David Koes](https://github.com/dkoes) for his awesome [py3Dmol](https://3dmol.csb.pitt.edu/) plugin, without whom these notebooks would be quite boring!
- A colab by Sergey Ovchinnikov (@sokrypton), Milot Mirdita (@milot_mirdita) and Martin Steinegger (@thesteinegger).

### How do I reference this work?

- Mirdita M, Schütze K, Moriwaki Y, Heo L, Ovchinnikov S and Steinegger M. ColabFold: Making protein folding accessible to all. <br />
  Nature Methods (2022) doi: [10.1038/s41592-022-01488-1](https://www.nature.com/articles/s41592-022-01488-1)
- If you’re using **AlphaFold**, please also cite: <br />
  Jumper et al. "Highly accurate protein structure prediction with AlphaFold." <br />
  Nature (2021) doi: [10.1038/s41586-021-03819-2](https://doi.org/10.1038/s41586-021-03819-2)
- If you’re using **AlphaFold-multimer**, please also cite: <br />
  Evans et al. "Protein complex prediction with AlphaFold-Multimer." <br />
  biorxiv (2021) doi: [10.1101/2021.10.04.463034v1](https://www.biorxiv.org/content/10.1101/2021.10.04.463034v1)
- If you are using **RoseTTAFold**, please also cite: <br />
  Minkyung et al. "Accurate prediction of protein structures and interactions using a three-track neural network." <br />
  Science (2021) doi: [10.1126/science.abj8754](https://doi.org/10.1126/science.abj8754)

[![DOI](https://zenodo.org/badge/doi/10.5281/zenodo.5123296.svg)](https://doi.org/10.5281/zenodo.5123296)
