# Container image

A container image bundles FinaleToolkit with a pinned Python and all of its
compiled dependencies (numpy, pysam, numba, pyBigWig, py2bit), so it runs
identically on a laptop, in CI, and on a cloud instance without anyone having
to reproduce the environment by hand.

## Build

```sh
docker/build.sh
```

This tags `finaletoolkit:<version>` and `finaletoolkit:latest`, where
`<version>` comes from git tags via setuptools-scm — the same version the wheel
would get.

Why a script rather than a bare `docker build`: `.git` is excluded from the
build context (it is ~200MB), so setuptools-scm cannot derive the version
inside the container. The script resolves it and passes it in as a build arg.
Building by hand works too, as long as you supply it:

```sh
docker build --build-arg SETUPTOOLS_SCM_PRETEND_VERSION=1.1.0 \
    --target runtime -t finaletoolkit:1.1.0 .
```

The build fails rather than producing an image that reports itself as
`0.0.0.dev0+unknown`.

Overrides, all optional:

| Variable | Default | Purpose |
| --- | --- | --- |
| `FTK_VERSION` | derived from git | Set the version explicitly |
| `FTK_IMAGE` | `finaletoolkit` | Image name |
| `FTK_PLATFORM` | `linux/amd64` | Target architecture |

`linux/amd64` is the default because py2bit and pyBigWig do not reliably
publish prebuilt ARM Linux wheels; on ARM they may have to compile from source.
On Apple Silicon the resulting image runs under emulation — correct, but slow.

## Run

The image behaves as the `finaletoolkit` command itself:

```sh
docker run --rm finaletoolkit:latest --help
docker run --rm finaletoolkit:latest --version
```

A container has its own filesystem that is discarded on exit, so input and
output have to be bind-mounted from the host. `/data` is the conventional
mount point and the image's working directory:

```sh
docker run --rm -v "/path/to/data:/data" finaletoolkit:latest \
    frag-length-bins /data/sample.bam -o /data/output.tsv
```

Paths in the command are container paths, not host paths — `/data/sample.bam`
refers to `/path/to/data/sample.bam` on the host.

The container runs as a non-root user (uid 1000). On Linux, output files land
on the host owned by that uid; add `--user "$(id -u):$(id -g)"` if you need
them owned by you. Docker Desktop on macOS handles this mapping already.

## Test

The `test` stage is the runtime image plus pytest. The ~91MB of fixtures under
`tests/` are mounted at run time rather than baked into the image:

```sh
docker/build.sh --target test
docker run --rm -v "$PWD/tests:/app/tests:ro" \
    --entrypoint pytest finaletoolkit:test
```

Tests write through pytest's `tmp_path` fixture, so the mount can stay
read-only.

Passing tests here is a stronger signal than passing on the host: it confirms
the packaged environment is correct, not just that the code imports.

## What is in the image

Two stages. The first installs the build tooling and produces a wheel; the
second installs that wheel into a clean `python:3.12-slim` base, so compilers
and build-time packages never ship in the final artifact.

The Python version is pinned to one version (3.12), unlike CI, which tests a
3.10–3.13 matrix. That is intentional: the matrix proves the library is
portable, while the image is a single reproducible deployment target.

`NUMBA_CACHE_DIR` and `MPLCONFIGDIR` are pointed at `/tmp` because numba and
matplotlib both want a writable cache directory and default to `$HOME`, which
is not writable under an arbitrary `--user`.
