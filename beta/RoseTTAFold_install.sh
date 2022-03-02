#!/bin/sh -e
conda create -n rosettafold -c conda-forge cudnn cudatoolkit==11.3.1 cudatoolkit-dev==11.3.1 python==3.7
conda activate rosettafold
pip install jupyter papermill matplotlib tensorflow-gpu==1.15
pip install torch==1.10.0+cu113 -f https://download.pytorch.org/whl/cu113/torch_stable.html
