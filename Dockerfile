ARG CUDA_VERSION=11.8.0
ARG COLABFOLD_VERSION=1.5.5
FROM nvidia/cuda:${CUDA_VERSION}-base-ubuntu22.04

RUN apt-get update && apt-get install -y wget cuda-nvcc-$(echo $CUDA_VERSION | cut -d'.' -f1,2 | tr '.' '-') --no-install-recommends --no-install-suggests && rm -rf /var/lib/apt/lists/* && \
    wget -qnc https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-Linux-x86_64.sh && \
    bash Mambaforge-Linux-x86_64.sh -bfp /usr/local && \
    mamba config --set auto_update_conda false && \
    rm -f Mambaforge-Linux-x86_64.sh && \
    CONDA_OVERRIDE_CUDA=$(echo $CUDA_VERSION | cut -d'.' -f1,2) mamba create -y -n colabfold -c conda-forge -c bioconda colabfold=$COLABFOLD_VERSION jaxlib==*=cuda* && \
    mamba clean -afy

ENV PATH /usr/local/envs/colabfold/bin:$PATH
ENV MPLBACKEND Agg

VOLUME cache
ENV MPLCONFIGDIR /cache
ENV XDG_CACHE_HOME /cache
