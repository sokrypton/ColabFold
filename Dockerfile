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

ARG CUDA=cuda12

VOLUME cache
ENV MPLBACKEND=Agg
ENV MPLCONFIGDIR=/cache
ENV XDG_CACHE_HOME=/cache

RUN apt-get update; \
    apt-get install -y wget git python3 python3-venv --no-install-suggests; \
    rm -rf /var/lib/apt/lists/*; \
    python3 -m venv /usr/local;
COPY --from=builder /opt/build/binaries/* /usr/local/bin/

WORKDIR /app
COPY . /app
RUN pip install --no-cache-dir \
        ".[alphafold,openmm]" \
        "jax[${CUDA}]<0.12" \
        "openmm[${CUDA}]"; \
    rm -rf /root/.cache/pip
