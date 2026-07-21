# ampliresist-figures — ampliresist-owned runtime for stage 18 publication figures.
#
# Replaces the borrowed `haema-figures` image (which belongs to the HÆMA pipeline
# and bundles seaborn/geopandas/shapely/pyproj that stage 18 never imports).
#
# bin/07_build_roadmap_figures.py imports ONLY numpy, pandas and matplotlib:
#   * GeoJSON is parsed with the standard-library `json` module (no geopandas);
#   * the site-adjusted penalised logistic model is a hand-rolled numpy IRLS
#     (no scipy / scikit-learn / statsmodels).
# Package versions are pinned to those the previous image shipped, so the
# deterministic figure source CSVs are reproduced exactly.
#
# Build:
#   docker build -t ghcr.io/owner/ampliresist-figures:1.0.0 \
#       -f containers/figures.Dockerfile containers/
#
# (Replace `owner` with your GitHub org/user before `docker push`.)
FROM python:3.11.15-slim-bookworm

LABEL org.opencontainers.image.title="ampliresist-figures"
LABEL org.opencontainers.image.version="1.0.0"
LABEL org.opencontainers.image.description="Matplotlib publication-figure runtime for ampliresist stage 18"
LABEL org.opencontainers.image.source="https://github.com/cyril-kuma/ampliresist"
LABEL org.opencontainers.image.licenses="MIT"

# procps provides `ps`, which Nextflow requires in every task container to
# collect run metrics; matplotlib needs a writable config dir under an arbitrary
# UID (Nextflow's docker.runOptions uses -u $(id -u):$(id -g)).
RUN apt-get update \
    && apt-get install -y --no-install-recommends procps \
    && rm -rf /var/lib/apt/lists/*
ENV MPLCONFIGDIR=/tmp/mplconfig \
    HOME=/tmp \
    PYTHONDONTWRITEBYTECODE=1
RUN mkdir -p /tmp/mplconfig && chmod 1777 /tmp/mplconfig

# Pinned to match the previous (haema-figures:0.4.0) runtime for output parity.
# DejaVu Sans (the figure font) ships inside the matplotlib wheel — no system fonts required.
RUN pip install --no-cache-dir \
        numpy==2.1.3 \
        pandas==2.2.3 \
        matplotlib==3.9.2

# Fail the build, not the pipeline, if the runtime is incomplete.
RUN python -c "import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot, numpy, pandas; \
print('ampliresist-figures OK numpy', numpy.__version__, 'pandas', pandas.__version__, 'matplotlib', matplotlib.__version__)"

CMD ["python3"]
