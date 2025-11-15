ARG MMSEQS_HASH=archive/8cc5ce367b5638c4306c2d7cfc652dd099a4643f
ARG downloader=${TARGETARCH}_downloader

FROM scratch AS amd64_downloader
ARG MMSEQS_HASH
WORKDIR /opt/build
ONBUILD ADD https://mmseqs.com/${MMSEQS_HASH}/mmseqs-linux-gpu.tar.gz  .

FROM scratch AS arm64_downloader
ARG MMSEQS_HASH
WORKDIR /opt/build
ONBUILD ADD https://mmseqs.com/${MMSEQS_HASH}/mmseqs-linux-gpu-arm64.tar.gz  .

FROM $downloader AS downloader

FROM debian:trixie-slim AS builder
WORKDIR /opt/build
COPY --from=downloader /opt/build/* .
RUN mkdir binaries; \
    for i in *.tar.gz; do \
        if [ -e ${i} ]; then \
            tar -xzvf ${i}; \
            mv -f -- */bin/* binaries/; \
        fi; \
    done; \
    chmod -R +x binaries;

FROM debian:trixie-slim

VOLUME cache
ENV MPLBACKEND=Agg
ENV MPLCONFIGDIR=/cache
ENV XDG_CACHE_HOME=/cache
ENV CONDA_VERSION=25.9.1-0

RUN apt-get update; \
    apt-get install -y wget git --no-install-suggests; \
    rm -rf /var/lib/apt/lists/*;

SHELL ["/bin/bash", "--login", "-x", "-c"]
RUN wget -qnc https://github.com/conda-forge/miniforge/releases/download/${CONDA_VERSION}/Miniforge3-${CONDA_VERSION}-Linux-x86_64.sh; \
    bash Miniforge3-${CONDA_VERSION}-Linux-x86_64.sh -bfp /usr/local; \
    conda config --set auto_update_conda false; \
    rm -f Miniforge3-${CONDA_VERSION}-Linux-x86_64.sh; \
    conda install -y -c conda-forge -c bioconda kalign2=2.04 hhsuite=3.3.0; \
    conda clean -afy; \
    conda shell.bash hook;
COPY --from=builder /opt/build/binaries/* /usr/local/bin/

WORKDIR /app
COPY . /app
RUN pip install --no-cache-dir \
        ".[alphafold-minus-jax]" \
        "jax[cuda]<0.8" \
        "dm-haiku[flax]==0.0.14" \
        "OpenMM[cuda12]==8.2.0" \
        'pdbfixer @ git+https://github.com/openmm/pdbfixer#94cfa4c0ca551cdc5f13320f9a658efd59f2b881'; \
    rm -rf /root/.cache/pip