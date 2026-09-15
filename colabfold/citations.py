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
    "Abramson2024": """@article{Abramson2024,
author = {Abramson, Josh and Adler, Jonas and Dunger, Jack and Evans, Richard and Green, Tim and Pritzel, Alexander and Ronneberger, Olaf and Willmore, Lindsay and Ballard, Andrew J. and Bambrick, Joshua and Bodenstein, Sebastian W. and Evans, David A. and Hung, Chia-Chun and O'Neill, Michael and Reiman, David and Tunyasuvunakool, Kathryn and Wu, Zachary and {\\v{Z}}emgulyt{\\.{e}}, Akvil{\\.{e}} and Arvaniti, Eirini and Beattie, Charles and Bertolli, Ottavia and Bridgland, Alex and Cherepanov, Alexey and Congreve, Miles and Cowen-Rivers, Alexander I. and Cowie, Andrew and Figurnov, Michael and Fuchs, Fabian B. and Gladman, Hannah and Jain, Rishub and Khan, Yousuf A. and Low, Caroline M. R. and Perlin, Kuba and Potapenko, Anna and Savy, Pascal and Singh, Sukhdeep and Stecula, Adrian and Thillaisundaram, Ashok and Tong, Catherine and Yakneen, Sergei and Zhong, Ellen D. and Zielinski, Michal and {\\v{Z}}{\\'{i}}dek, Augustin and Bapst, Victor and Kohli, Pushmeet and Jaderberg, Max and Hassabis, Demis and Jumper, John M.},
doi = {10.1038/s41586-024-07487-w},
journal = {Nature},
title = {{Accurate structure prediction of biomolecular interactions with AlphaFold 3}},
year = {2024},
comment = {AlphaFold3 architecture}
}""",
    "OpenFold3": """@software{OpenFold3,
author = {{The OpenFold3 Team}},
doi = {10.5281/zenodo.19001000},
title = {{OpenFold3-preview}},
version = {0.4.2},
year = {2025},
comment = {OpenFold3 and OpenBind weights}
}""",
    "OpenBind": """@article{OpenBind,
author = {Nelen, Jochem and Khan, Omeir and Adams, Etowah and Aschenbrenner, Jasmin C. and Thompson, Warren and Ebrahim, Ali and {\\c C}apkin, Eda and Vall{\\'e}e, C{\\'e}dric and OpenBind and Shotton, Elizabeth J. and Griffen, Ed J. and Chodera, John D. and Deane, Charlotte M. and von Delft, Frank and AlQuraishi, Mohammed and Imrie, Fergus},
doi = {10.64898/2026.08.27.747600},
journal = {bioRxiv},
title = {The first OpenBind release: An open experimental structure{\\textendash}affinity dataset and benchmark for structure-based AI},
year = {2026},
comment = {OpenFold3 and OpenBind}
}""",
    "ProtenixV1": """@article{ProtenixV1,
author = {Protenix Team and Zhang, Yuxuan and Gong, Chengyue and Zhang, Hanyu and Ma, Wenzhi and Liu, Zhenyu and Chen, Xinshi and Guan, Jiaqi and Wang, Lan and Yang, Yanping and Xia, Yu and Xiao, Wenzhi},
doi = {10.64898/2026.02.05.703733},
journal = {bioRxiv},
title = {Protenix-v1: Toward High-Accuracy Open-Source Biomolecular Structure Prediction},
year = {2026},
comment = {Protenix-v1 weights}
}""",
    "ProtenixV2": """@article{ProtenixV2,
author = {Zhang, Yuxuan and Gong, Chengyue and Sun, Jinyuan and Guan, Jiaqi and Ren, Milong and Xue, Song and Zhang, Hanyu and Ma, Wenzhi and Liu, Zhenyu and Chen, Xinshi and Xiao, Wenzhi},
doi = {10.64898/2026.04.10.717613},
journal = {bioRxiv},
title = {Protenix-v2: Broadening the Reach of Structure Prediction and Biomolecular Design},
year = {2026},
comment = {Protenix-v2 weights}
}""",
    "Boltz2": """@article{Boltz2,
author = {Passaro, Saro and Corso, Gabriele and Wohlwend, Jeremy and Reveiz, Mateo and Thaler, Stephan and Somnath, Vignesh Ram and Getz, Noah and Portnoi, Tally and Roy, Julien and Stark, Hannes and Kwabi-Addo, David and Beaini, Dominique and Jaakkola, Tommi and Barzilay, Regina},
doi = {10.1101/2025.06.14.659707},
journal = {bioRxiv},
title = {Boltz-2: Towards Accurate and Efficient Binding Affinity Prediction},
year = {2025},
comment = {Boltz-2 weights}
}""",
    "Chai1": """@article{Chai1,
author = {Chai Discovery and Boitreaud, Jacques and Dent, Jack and McPartlon, Matthew and Meier, Joshua and Reis, Vinicius and Rogozhnikov, Alex and Wu, Kevin},
doi = {10.1101/2024.10.10.615955},
journal = {bioRxiv},
title = {Chai-1: Decoding the molecular interactions of life},
year = {2024},
comment = {Chai-1 weights}
}""",
    "IntelliFold2": """@article{IntelliFold2,
author = {Qiao, Lifeng and Yan, He and Liu, Gary and Guo, Gaoxing and Sun, Siqi},
doi = {10.64898/2026.02.09.704787},
journal = {bioRxiv},
title = {IntelliFold-2: Surpassing AlphaFold 3 via Architectural Refinement and Structural Consistency},
year = {2026},
comment = {IntelliFold-2 weights}
}""",
    "RoseTTAFold3": """@article{RoseTTAFold3,
author = {Corley, Nathaniel and Mathis, Simon and Krishna, Rohith and Bauer, Magnus S. and Thompson, Tuscan R. and Ahern, Woody and Kazman, Maxwell W. and Brent, Rafael I. and Didi, Kieran and Kubaney, Andrew and McHugh, Lilian and Nagle, Arnav and Favor, Andrew and Kshirsagar, Meghana and Sturmfels, Pascal and Li, Yanjing and Butcher, Jasper and Qiang, Bo and Schaaf, Lars L. and Mitra, Raktim and Campbell, Katelyn and Zhang, Odin and Weissman, Roni and Humphreys, Ian R. and Cong, Qian and Funk, Jonathan and Sonthalia, Shreyash and Li{\\`o}, Pietro and Baker, David and DiMaio, Frank},
doi = {10.1101/2025.08.14.670328},
journal = {bioRxiv},
title = {Accelerating Biomolecular Modeling with AtomWorks and RF3},
year = {2025},
comment = {RoseTTAFold3 weights}
}""",
    "OpenDDE": """@misc{OpenDDE,
author = {Aureka AI OpenDDE project},
title = {{Folding, Reasoning, and Scaling with Open-source Drug Discovery Engine}},
year = {2026},
eprint = {2607.03787},
archivePrefix = {arXiv},
primaryClass = {cs.AI},
comment = {OpenDDE weights}
}""",
    "ESMFold2": """@article{ESMFold2,
author = {Candido, Salvatore and Hayes, Thomas and Derry, Alexander and Rao, Roshan and Lin, Zeming and Verkuil, Robert and Wu, Bryan Z. and Lee, Jin Sub and Bruguera, Elise S. and Keval, Jehan A. and Kopylov, Mykhailo and Pak, John E. and Wu, Wesley and Thomas, Neil and Mataraso, Samson and Hsu, Alvin and Trotman-Grant, Ashton C. and Fatras, Kilian and dos Santos Costa, Allan and Badkundri, Rohil and Ak{\\i}n, Halil and Oktay, Deniz and Deaton, Jonathan and Montabana, Elizabeth and Sitwala, Hrishita and Yu, Yue and Wiggert, Marius and Carlin, Dylan Alexander and Goering, Anthony W. and Blazejewski, Tomasz and Sandora, McCullen and Hla, Michael and Jia, Tina Z. and Kloker, Leon H. and Sofroniew, Nicholas J. and Uehara, Masatoshi and Pannu, Jassi and Bachas, Sharrol and Liu, Daniel S. and Sercu, Tom and Rives, Alexander},
doi = {10.64898/2026.06.03.729735},
journal = {bioRxiv},
title = {Language Modeling Materializes a World Model of Protein Biology},
year = {2026},
comment = {ESMFold2 weights}
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
    "VanKempen2023": """@article{VanKempen2023,
author = {van Kempen, Michel and Kim, Stephanie S and Tumescheit, Charlotte and Mirdita, Milot and Lee, Jeongjae and Gilchrist, Cameron L M and S{\"{o}}ding, Johannes and Steinegger, Martin},
doi = {10.1038/s41587-023-01773-0},
journal = {Nature Biotechnology},
title = {{Fast and accurate protein structure search with Foldseek}},
year = {2023},
comment = {PDB100 database}
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
comment = {Uniclust30/UniRef30 database}
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
    "Lee2023": """@article{Lee2023,
author = {Lee, Jae-Won and Won, Jong-Hyun and Jeon, Seonggwang and Choo, Yujin and Yeon, Yubin and Oh, Jin-Seon and Kim, Minsoo and Kim, SeonHwa and Joung, InSuk and Jang, Cheongjae and Lee, Sung Jong and Kim, Tae Hyun and Jin, Kyong Hwan and Song, Giltae and Kim, Eun-Sol and Yoo, Jejoong and Paek, Eunok and Noh, Yung-Kyun and Joo, Keehyoung},
title = "{DeepFold: enhancing protein structure prediction through optimized loss functions, improved template features, and re-optimized energy function}",
journal = {Bioinformatics},
volume = {39},
number = {12},
pages = {btad712},
year = {2023},
month = {11},
doi = {10.1093/bioinformatics/btad712},
comment = {DeepFold-v1 Model}
}
""",
}


# whose parameters a model in the alphafold3 family actually runs on
_AF3_WEIGHTS = {"openfold3": ("OpenFold3", "OpenBind"), "of3": ("OpenFold3", "OpenBind"),
                "openbind": ("OpenFold3", "OpenBind"), "openbind0": ("OpenFold3", "OpenBind"),
                "protenix1": "ProtenixV1",
                "protenix": "ProtenixV2", "protenix2": "ProtenixV2", "boltz2": "Boltz2",
                "chai": "Chai1", "chai1": "Chai1", "intellifold": "IntelliFold2",
                "if2": "IntelliFold2", "intellifold2": "IntelliFold2",
                "rf3": "RoseTTAFold3", "rosettafold3": "RoseTTAFold3", "opendde": "OpenDDE",
                "esmfold2": "ESMFold2", "esmfold2_fast": "ESMFold2",
                "esmfold2_lm300m": "ESMFold2", "esmfold2_lm600m": "ESMFold2"}


def af3_citations(model: str) -> list:
    from colabfold.backend import is_af3_model

    if not model or not is_af3_model(model):
        return []
    name = model[len("alphafold3_"):] if model.startswith("alphafold3_") else model
    weights = _AF3_WEIGHTS.get(name, ())
    if isinstance(weights, str):
        weights = (weights,)
    return ["Abramson2024"] + list(weights)


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
    if model == "alphafold2_ptm" or model == "alphafold2":
        to_cite += ["Jumper2021"]
    if model == "deepfold_v1":
        to_cite += ["Lee2023"]
    if model.startswith("alphafold2_multimer"):
        to_cite += ["Evans2021"]
    to_cite += af3_citations(model)
    if use_msa:
        to_cite += ["Mirdita2019"]
    if use_msa:
        to_cite += ["Mirdita2017"]
    if use_env:
        to_cite += ["Mitchell2019"]
    if use_templates:
        to_cite += ["VanKempen2023"]
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