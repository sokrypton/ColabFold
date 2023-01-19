import logging
from pathlib import Path

logger = logging.getLogger(__name__)

citations = {
    "Mirdita2021": """@article{Mirdita2022,
author= {Mirdita, Milot and Schütze, Konstantin and Moriwaki, Yoshitaka and Heo, Lim and Ovchinnikov, Sergey and Steinegger, Martin },
doi = {10.1038/s41592-022-01488-1},
journal = {Nature Methods},
title = {{ColabFold: Making Protein folding accessible to all}},
year = {2022},
comment = {ColabFold including MMseqs2 MSA server}
}""",
    "Mitchell2019": """@article{Mitchell2019,
author = {Mitchell, Alex L and Almeida, Alexandre and Beracochea, Martin and Boland, Miguel and Burgin, Josephine and Cochrane, Guy and Crusoe, Michael R and Kale, Varsha and Potter, Simon C and Richardson, Lorna J and Sakharova, Ekaterina and Scheremetjew, Maxim and Korobeynikov, Anton and Shlemov, Alex and Kunyavskaya, Olga and Lapidus, Alla and Finn, Robert D},
doi = {10.1093/nar/gkz1035},
journal = {Nucleic Acids Res.},
title = {{MGnify: the microbiome analysis resource in 2020}},
year = {2019},
comment = {MGnify database}
}""",
    "Eastman2017": """@article{Eastman2017,
author = {Eastman, Peter and Swails, Jason and Chodera, John D. and McGibbon, Robert T. and Zhao, Yutong and Beauchamp, Kyle A. and Wang, Lee-Ping and Simmonett, Andrew C. and Harrigan, Matthew P. and Stern, Chaya D. and Wiewiora, Rafal P. and Brooks, Bernard R. and Pande, Vijay S.},
doi = {10.1371/journal.pcbi.1005659},
journal = {PLOS Comput. Biol.},
number = {7},
title = {{OpenMM 7: Rapid development of high performance algorithms for molecular dynamics}},
volume = {13},
year = {2017},
comment = {Amber relaxation}
}""",
    "Jumper2021": """@article{Jumper2021,
author = {Jumper, John and Evans, Richard and Pritzel, Alexander and Green, Tim and Figurnov, Michael and Ronneberger, Olaf and Tunyasuvunakool, Kathryn and Bates, Russ and {\v{Z}}{\'{i}}dek, Augustin and Potapenko, Anna and Bridgland, Alex and Meyer, Clemens and Kohl, Simon A. A. and Ballard, Andrew J. and Cowie, Andrew and Romera-Paredes, Bernardino and Nikolov, Stanislav and Jain, Rishub and Adler, Jonas and Back, Trevor and Petersen, Stig and Reiman, David and Clancy, Ellen and Zielinski, Michal and Steinegger, Martin and Pacholska, Michalina and Berghammer, Tamas and Bodenstein, Sebastian and Silver, David and Vinyals, Oriol and Senior, Andrew W. and Kavukcuoglu, Koray and Kohli, Pushmeet and Hassabis, Demis},
doi = {10.1038/s41586-021-03819-2},
journal = {Nature},
pmid = {34265844},
title = {{Highly accurate protein structure prediction with AlphaFold.}},
year = {2021},
comment = {AlphaFold2 + BFD Database}
}""",
    "Evans2021": """@article{Evans2021,
  author   = {Evans, Richard and O'Neill, Michael and Pritzel, Alexander and Antropova, Natasha and Senior, Andrew and Green, Tim and  Zidek, Augustin and Bates, Russ and Blackwell, Sam and Yim, Jason and Ronneberger, Olaf and Bodenstein, Sebastian and Zielinski, Michal and Bridgland, Alex and Potapenko, Anna and Cowie, Andrew and Tunyasuvunakool, Kathryn and Jain, Rishub and Clancy, Ellen and Kohli, Pushmeet and Jumper, John and Hassabis, Demis},
  doi    = {10.1101/2021.10.04.463034v1},
  journal  = {bioRxiv},
  title    = {{Protein complex prediction with AlphaFold-Multimer}},
  year     =  {2021},
  comment = {AlphaFold2-multimer}
}""",
    "Mirdita2019": """@article{Mirdita2019,
author = {Mirdita, Milot and Steinegger, Martin and S{\"{o}}ding, Johannes},
doi = {10.1093/bioinformatics/bty1057},
journal = {Bioinformatics},
number = {16},
pages = {2856--2858},
pmid = {30615063},
title = {{MMseqs2 desktop and local web server app for fast, interactive sequence searches}},
volume = {35},
year = {2019},
comment = {MMseqs2 search server}
}""",
    "Steinegger2019": """@article{Steinegger2019,
author = {Steinegger, Martin and Meier, Markus and Mirdita, Milot and V{\"{o}}hringer, Harald and Haunsberger, Stephan J. and S{\"{o}}ding, Johannes},
doi = {10.1186/s12859-019-3019-7},
journal = {BMC Bioinform.},
number = {1},
pages = {473},
pmid = {31521110},
title = {{HH-suite3 for fast remote homology detection and deep protein annotation}},
volume = {20},
year = {2019},
comment = {PDB70 database}
}""",
    "Mirdita2017": """@article{Mirdita2017,
author = {Mirdita, Milot and von den Driesch, Lars and Galiez, Clovis and Martin, Maria J. and S{\"{o}}ding, Johannes and Steinegger, Martin},
doi = {10.1093/nar/gkw1081},
journal = {Nucleic Acids Res.},
number = {D1},
pages = {D170--D176},
pmid = {27899574},
title = {{Uniclust databases of clustered and deeply annotated protein sequences and alignments}},
volume = {45},
year = {2017},
comment = {Uniclust30/UniRef30 database},
}""",
    "Berman2003": """@misc{Berman2003,
author = {Berman, Helen and Henrick, Kim and Nakamura, Haruki},
booktitle = {Nat. Struct. Biol.},
doi = {10.1038/nsb1203-980},
number = {12},
pages = {980},
pmid = {14634627},
title = {{Announcing the worldwide Protein Data Bank}},
volume = {10},
year = {2003},
comment = {templates downloaded from wwPDB server}
}""",
}


def write_bibtex(
    model: str,
    use_msa: bool,
    use_env: bool,
    use_templates: bool,
    use_amber: bool,
    result_dir: Path,
    bibtex_file: str = "cite.bibtex",
) -> Path:
    to_cite = ["Mirdita2021"]
    if model == "AlphaFold2-ptm":
        to_cite += ["Jumper2021"]
    if model.startswith("AlphaFold2-multimer"):
        to_cite += ["Evans2021"]
    if use_msa:
        to_cite += ["Mirdita2019"]
    if use_msa:
        to_cite += ["Mirdita2017"]
    if use_env:
        to_cite += ["Mitchell2019"]
    if use_templates:
        to_cite += ["Steinegger2019"]
    if use_templates:
        to_cite += ["Berman2003"]
    if use_amber:
        to_cite += ["Eastman2017"]

    bibtex_file = result_dir.joinpath(bibtex_file)
    with bibtex_file.open("w", encoding="utf-8") as writer:
        for i in to_cite:
            writer.write(citations[i])
            writer.write("\n")

    logger.info(f"Found {len(to_cite)} citations for tools or databases")
    return bibtex_file