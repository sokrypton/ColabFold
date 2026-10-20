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
# Keep installed code root-owned; only runtime data belongs to the app user.
RUN groupadd --gid 10001 caat \
    && useradd --uid 10001 --gid caat --create-home caat \
    && mkdir -p /home/caat/.cache /home/caat/.config /work \
    && chown -R caat:caat /home/caat /work
ENV HOME=/home/caat
WORKDIR /work
USER 10001:10001
RUN test "$(id -u)" -ne 0 \
    && touch /work/.write-check "$HOME/.write-check" "$HOME/.cache/.write-check" \
    && rm /work/.write-check "$HOME/.write-check" "$HOME/.cache/.write-check" \
    && JAX_PLATFORMS=cpu python /opt/CAAT/images/smoke_test.py
ENTRYPOINT ["/opt/CAAT/images/entrypoint.sh"]
CMD ["--help"]
