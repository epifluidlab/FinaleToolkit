# syntax=docker/dockerfile:1

# Two-stage build. The first stage has the tooling needed to produce a wheel;
# the second installs that wheel into a clean image, so the compilers and
# build-time packages never ship in the final artifact.
#
# Build with docker/build.sh, which derives the version from git tags and
# passes it in. See docker/README.md.

ARG PYTHON_VERSION=3.12


# --- Stage 1: build the wheel ------------------------------------------------
FROM python:${PYTHON_VERSION}-slim AS builder

# setuptools-scm normally derives the version from git history, but .git is
# excluded from the build context (it is ~200MB, and shipping it to the daemon
# on every build is slow). The version is passed in instead. Without this,
# setuptools-scm silently falls back to "0.0.0.dev0+unknown" -- the runtime
# stage asserts against exactly that.
ARG SETUPTOOLS_SCM_PRETEND_VERSION
ENV SETUPTOOLS_SCM_PRETEND_VERSION=${SETUPTOOLS_SCM_PRETEND_VERSION}

WORKDIR /src

RUN --mount=type=cache,target=/root/.cache/pip \
    pip install build

# MANIFEST.in and README.md are referenced by the build (readme is declared in
# pyproject.toml), so the build fails without them even though only src/ ends
# up in the wheel.
COPY pyproject.toml MANIFEST.in README.md LICENSE ./
COPY src/ ./src/

RUN --mount=type=cache,target=/root/.cache/pip \
    python -m build --wheel --outdir /dist


# --- Stage 2: runtime --------------------------------------------------------
FROM python:${PYTHON_VERSION}-slim AS runtime

LABEL org.opencontainers.image.title="FinaleToolkit" \
      org.opencontainers.image.description="Extract fragmentation features of cell-free DNA from paired-end sequencing data." \
      org.opencontainers.image.source="https://github.com/epifluidlab/FinaleToolkit" \
      org.opencontainers.image.documentation="https://epifluidlab.github.io/FinaleToolkit/" \
      org.opencontainers.image.licenses="MIT"

# numba and matplotlib both want a writable cache/config directory and fall
# back to $HOME, which the non-root user below cannot always write to
# (e.g. when the image is run with an arbitrary --user). Point them at /tmp,
# which is writable regardless.
#
# PYTHONDONTWRITEBYTECODE is deliberately *not* set: pip precompiles .pyc
# files at install time, and keeping them means the CLI does not re-parse
# scipy/numba/matplotlib source on every invocation. This is a short-lived
# process started many times, so startup cost matters more than image size.
ENV NUMBA_CACHE_DIR=/tmp/numba-cache \
    MPLCONFIGDIR=/tmp/matplotlib \
    PYTHONUNBUFFERED=1

# Bind-mounting the wheel rather than COPYing it keeps it out of the image's
# layers -- the installed package is all we need.
RUN --mount=type=bind,from=builder,source=/dist,target=/dist \
    --mount=type=cache,target=/root/.cache/pip \
    pip install /dist/*.whl

# Two smoke checks that catch the failure modes that would otherwise ship
# silently:
#
#  1. Version. Without SETUPTOOLS_SCM_PRETEND_VERSION the package falls back to
#     "0.0.0.dev0+unknown" and the image is untraceable.
#  2. Package data. The .tsv/.gz files under frag/data and genome/data are
#     pulled in by [tool.setuptools.package-data]; if that ever stops working,
#     the CLI installs cleanly and then fails at run time on the subcommands
#     that need them.
#  3. Imports. Loading the CLI module confirms the compiled wheels (pysam,
#     numba, pyBigWig, py2bit) actually load on this base image.
#
# Bind-mounted rather than COPYed so the script leaves no layer behind.
RUN --mount=type=bind,source=docker/verify.py,target=/tmp/verify.py \
    python /tmp/verify.py

RUN useradd --create-home --uid 1000 finaletoolkit \
    && mkdir -p /data \
    && chown finaletoolkit:finaletoolkit /data

USER finaletoolkit

# Conventional mount point for input BAM/CRAM files and outputs:
#   docker run -v /host/path:/data <image> <subcommand> /data/sample.bam ...
WORKDIR /data

# ENTRYPOINT makes the image behave as the CLI itself, so `docker run <image>
# frag-length-bins ...` works the same way the installed command does.
ENTRYPOINT ["finaletoolkit"]
CMD ["--help"]


# --- Stage 3: test -----------------------------------------------------------
# A strict superset of the runtime image: same environment, plus pytest. The
# 91MB of fixtures under tests/ are mounted at run time rather than baked in.
#
#   docker build --target test -t finaletoolkit:test .
#   docker run --rm -v "$PWD/tests:/app/tests:ro" --entrypoint pytest \
#       finaletoolkit:test
#
# The ENTRYPOINT is deliberately left as-is so this image still works as a
# plain CLI if it is built by accident (it is the last stage, so a bare
# `docker build .` lands here).
FROM runtime AS test

USER root
RUN --mount=type=cache,target=/root/.cache/pip \
    pip install pytest pytest-cov

# pytest reads [tool.pytest.ini_options] from pyproject.toml, which sets
# testpaths = ["tests"]; this also makes /app the rootdir. It has to be
# writable by the non-root user, since pytest writes .pytest_cache there and
# --cov writes .coverage.
RUN mkdir -p /app && chown finaletoolkit:finaletoolkit /app
WORKDIR /app
COPY --chown=finaletoolkit:finaletoolkit pyproject.toml ./

USER finaletoolkit
