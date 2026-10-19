FROM python:3.11-slim-bookworm

WORKDIR /opt/CAAT
COPY pyproject.toml poetry.lock README.md LICENSE ./
COPY colabfold ./colabfold
COPY scripts ./scripts
COPY dashboard ./dashboard
COPY alphafold ./alphafold
COPY tests ./tests
COPY test-data ./test-data
COPY images ./images
RUN sh /opt/CAAT/images/install.sh

ENV PATH="/opt/CAAT/.venv/bin:/opt/CAAT:${PATH}" \
    MPLBACKEND=Agg \
    PYTHONNOUSERSITE=1
RUN JAX_PLATFORMS=cpu python /opt/CAAT/images/smoke_test.py
WORKDIR /work
ENTRYPOINT ["/opt/CAAT/images/entrypoint.sh"]
CMD ["--help"]
