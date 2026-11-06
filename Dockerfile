# syntax=docker/dockerfile:1.7

###############################################################################
# Stage 1 — builder: compiles R packages with dev headers + toolchain
###############################################################################
FROM rocker/r-ver:4.6.0 AS builder

ENV DEBIAN_FRONTEND=noninteractive

# Dev headers + toolchain needed ONLY to compile R packages from source.
RUN apt-get update && apt-get install -y --no-install-recommends \
      curl wget \
      cmake pkg-config git \
      libuv1-dev \
      libcurl4-openssl-dev libssl-dev libxml2-dev \
      libfontconfig1-dev libfreetype6-dev \
      libpng-dev libtiff-dev libjpeg-dev \
      libharfbuzz-dev libfribidi-dev \
      libhdf5-dev libglpk-dev libgmp-dev libgit2-dev \
      libmagick++-dev libx11-dev \
 && rm -rf /var/lib/apt/lists/*

WORKDIR /app

# Cache MUST be disabled here: with the cache enabled, renv::restore installs
# packages into RENV_PATHS_CACHE and only symlinks them into /app/renv/library.
# Those symlinks then dangle in the runtime image (the cache directory is not
# carried by `COPY --from=builder /app /app`), so renv tries to re-bootstrap
# at container startup and fails with no network.
ENV RENV_CONFIG_SANDBOX_ENABLED=FALSE \
    RENV_CONFIG_CACHE_ENABLED=FALSE

COPY renv.lock renv.lock
COPY renv/ renv/
COPY .Rprofile .Rprofile

# Install renv + BiocManager up front so Bioc repos resolve cleanly,
# then materialise the full library directly under /app/renv/library.
RUN R -e "install.packages(c('renv', 'BiocManager'))" \
 && R -e "renv::restore(confirm=FALSE)"

# Slim the library: drop help/html/doc/demo/tests and strip .so debug symbols.
# This alone removes ~800 MB – 1.2 GB from the R library.
RUN set -eux; \
    for LIB in $(R --no-save --slave -e "cat(.libPaths(), sep='\n')"); do \
      [ -d "$LIB" ] || continue; \
      find "$LIB" -mindepth 2 -maxdepth 2 -type d \
        \( -name help -o -name html -o -name doc -o -name demo -o -name tests \) \
        -exec rm -rf {} +; \
      find "$LIB" -name "*.so" -exec strip --strip-unneeded {} + 2>/dev/null || true; \
    done


###############################################################################
# Stage 2 — runtime: only the shared libs + Chromium + Python
###############################################################################
FROM rocker/r-ver:4.6.0

ENV DEBIAN_FRONTEND=noninteractive

# Runtime shared libraries only (no -dev, no cmake, no build-essential).
# Package names below are for Ubuntu 24.04 (noble), which is what
# rocker/r-ver:4.6.0 is based on. If any name is wrong for your base image
# the build will fail with a clear "unable to locate package" error.
RUN apt-get update && apt-get install -y --no-install-recommends \
      software-properties-common \
 && add-apt-repository -y ppa:xtradeb/apps \
 && apt-get update && apt-get install -y --no-install-recommends \
      # runtime shared libs (match the -dev list from stage 1)
      libcurl4t64 libssl3t64 libxml2 \
      libuv1t64 \
      libfontconfig1 libfreetype6 \
      libpng16-16t64 libtiff6 libjpeg-turbo8 \
      libharfbuzz0b libfribidi0 \
      libhdf5-103-1t64 libhdf5-cpp-103-1t64 \
      libglpk40 libgmp10 libgit2-1.7 \
      libmagick++-6.q16-9t64 libmagickcore-6.q16-7t64 \
      libx11-6 \
      pandoc \
      # Chromium (chromote/webshot2)
      chromium fonts-liberation \
      libnss3 libatk-bridge2.0-0 libxkbcommon0 libxcomposite1 \
      libxdamage1 libxrandr2 libgbm1 libgtk-3-0 libasound2t64 \
      libpangocairo-1.0-0 libxshmfence1 \
      # Python for reticulate
      python3 python3-pip libpython3-dev\
 && pip3 install --no-cache-dir --break-system-packages \
      numpy umap-learn scikit-learn \
 && rm -rf /var/lib/apt/lists/* /root/.cache/pip \
 && apt-get clean

# Chromium wrapper (required so chromote can run as root inside Docker)
RUN printf '#!/bin/bash\nexec /usr/bin/chromium --no-sandbox --disable-dev-shm-usage "$@"\n' \
      > /usr/local/bin/chromium-wrapper \
 && chmod +x /usr/local/bin/chromium-wrapper

ENV CHROMOTE_CHROME=/usr/local/bin/chromium-wrapper \
    RETICULATE_PYTHON=/usr/bin/python3 \
    RENV_CONFIG_SANDBOX_ENABLED=FALSE \
    RENV_CONFIG_CACHE_ENABLED=FALSE \
    LANG=C.UTF-8 \
    LC_ALL=C.UTF-8

WORKDIR /app

# Bring the fully-built R project (library + renv metadata) from the builder.
COPY --from=builder /app /app

# Pipeline code (changes often, kept last to maximise layer cache hits).
COPY R/ /app/R/
COPY _targets.R /app/
COPY run_pipeline.R /app/
COPY config_builder.R /app/
COPY config_template.yml /app/

RUN mkdir -p /data
ENV DRUGSEQ_CONFIG=/data/config.yml

ARG VERSION=1.0
LABEL org.opencontainers.image.version=$VERSION \
      org.opencontainers.image.title="ADEPT" \
      org.opencontainers.image.description="ADEPT: Automated Differential Expression Pipeline for high-throughput Transcriptomics screens"

ENTRYPOINT ["Rscript", "/app/run_pipeline.R"]
