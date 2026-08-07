#!/usr/bin/env bash
#
# Build the FinaleToolkit container image.
#
# .git is excluded from the build context, so setuptools-scm cannot derive the
# version inside the container. This script resolves it here and passes it in.
#
# Usage:
#   docker/build.sh                  # build the runtime image
#   docker/build.sh --target test    # build the test image instead
#   FTK_VERSION=1.1.0 docker/build.sh
#
# Any extra arguments are forwarded to `docker build`.

set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$repo_root"

image="${FTK_IMAGE:-finaletoolkit}"
platform="${FTK_PLATFORM:-linux/amd64}"

# Mirrors setuptools-scm's default version_scheme ("guess-next-dev") combined
# with the local_scheme = "no-local-version" set in pyproject.toml:
#
#   exactly on tag v1.1.0, clean   -> 1.1.0
#   4 commits past v1.1.0          -> 1.1.1.dev4
#
# Only used when setuptools-scm is not importable. If the version_scheme in
# pyproject.toml ever changes, this needs to change with it.
derive_version_from_git() {
    local describe tag distance dirty base
    describe="$(git describe --tags --long --dirty 2>/dev/null)" || return 1

    dirty=""
    if [[ "$describe" == *-dirty ]]; then
        dirty="yes"
        describe="${describe%-dirty}"
    fi

    # v1.1.0-4-g3ca99fd -> tag=v1.1.0, distance=4
    tag="${describe%-*-g*}"
    distance="${describe%-g*}"
    distance="${distance##*-}"

    base="${tag#v}"

    if [[ "$distance" == "0" && -z "$dirty" ]]; then
        printf '%s\n' "$base"
        return 0
    fi

    # Bump the last numeric component: 1.1.0 -> 1.1.1
    local prefix last
    last="${base##*.}"
    prefix="${base%.*}"
    if [[ "$last" =~ ^[0-9]+$ ]]; then
        printf '%s.%s.dev%s\n' "$prefix" "$((last + 1))" "$distance"
    else
        return 1
    fi
}

# setuptools-scm is authoritative when available; the git fallback keeps this
# script usable in an environment that only has the runtime dependencies.
version="${FTK_VERSION:-}"
if [[ -z "$version" ]]; then
    version="$(python -m setuptools_scm 2>/dev/null || true)"
fi
if [[ -z "$version" ]]; then
    version="$(derive_version_from_git || true)"
    [[ -n "$version" ]] && echo "note: setuptools-scm unavailable; derived ${version} from git" >&2
fi
if [[ -z "$version" ]]; then
    cat >&2 <<'EOF'
error: could not determine the package version.

setuptools-scm is not importable and the version could not be derived from git
(no tags reachable, or not a git checkout). Either:

    pip install setuptools-scm

or set it explicitly:

    FTK_VERSION=1.1.0 docker/build.sh
EOF
    exit 1
fi

# Default to the runtime stage. Without this a bare `docker build` would build
# the last stage in the Dockerfile, which is the larger test image.
target="runtime"
args=()
while [[ $# -gt 0 ]]; do
    case "$1" in
        --target)
            target="$2"
            shift 2
            ;;
        --target=*)
            target="${1#*=}"
            shift
            ;;
        *)
            args+=("$1")
            shift
            ;;
    esac
done

# Non-runtime stages get a suffixed tag so that building the test image does
# not overwrite the runtime image's tags.
tags=()
if [[ "$target" == "runtime" ]]; then
    tags+=(--tag "${image}:${version}" --tag "${image}:latest")
else
    tags+=(--tag "${image}:${version}-${target}" --tag "${image}:${target}")
fi

echo "Building ${image} ${version} (target=${target}, platform=${platform})"

set -x
docker build \
    --platform "$platform" \
    --target "$target" \
    --build-arg "SETUPTOOLS_SCM_PRETEND_VERSION=${version}" \
    "${tags[@]}" \
    "${args[@]+"${args[@]}"}" \
    .
