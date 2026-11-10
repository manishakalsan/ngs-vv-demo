# Reproducible runtime: the exact package set from the explicit lock file.
#   docker build -t ngs-vv-demo:1.0.0 .
#   docker run --rm ngs-vv-demo:1.0.0                      # runs the full test suite
#   docker run --rm -v "$PWD/out:/app/results" ngs-vv-demo:1.0.0 snakemake --cores 2
FROM mambaorg/micromamba:2.0.5

COPY --chown=$MAMBA_USER:$MAMBA_USER envs/conda-linux-64.lock /tmp/conda.lock
RUN micromamba install -y -n base -f /tmp/conda.lock && micromamba clean --all --yes

WORKDIR /app
COPY --chown=$MAMBA_USER:$MAMBA_USER . /app
ARG MAMBA_DOCKERFILE_ACTIVATE=1

LABEL org.opencontainers.image.title="ngs-vv-demo" \
      org.opencontainers.image.version="1.0.0" \
      org.opencontainers.image.description="Verified mini variant-calling pipeline (demo)"

CMD ["pytest", "-q"]
