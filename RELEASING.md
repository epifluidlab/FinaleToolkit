# Releasing FinaleToolkit

This describes the CI automation around `CHANGELOG.md` and cutting a release.
The workflows referenced here live in `.github/workflows/`; the CHANGELOG
helper script lives in `scripts/prepare_release.py`.

## While working on a PR

`CHANGELOG.md` follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/).
If your PR changes anything under `src/finaletoolkit/`, add a bullet point
under the `## [Unreleased]` heading at the top of the file, in the relevant
`### Added` / `### Changed` / `### Fixed` / `### Removed` subsection (create
the subsection if it doesn't exist yet). Write it for a user of the package,
not as a commit-log summary -- what changed and why it matters to them.

**[changelog-check.yml](.github/workflows/changelog-check.yml)** enforces
this: a PR that touches `src/finaletoolkit/` without also touching
`CHANGELOG.md` fails CI. If your change genuinely has no user-visible effect
(test-only changes, internal refactors, CI/infra work), add the
`skip-changelog` label to the PR instead of writing an entry -- adding or
removing the label re-runs the check.

## Cutting a release

Releases are ad hoc: whenever a maintainer decides enough has landed in
`[Unreleased]` to justify one, not on a fixed schedule. There's no
conventional-commit history to infer a version bump from automatically, so a
human picks the version number.

1. **Go to Actions -> "Prepare release" -> Run workflow**, and enter the new
   version (e.g. `1.1.0`, no leading `v`).
   ([release-prepare.yml](.github/workflows/release-prepare.yml))

   This moves everything currently under `[Unreleased]` into a new
   `## [1.1.0] - <today>` heading, leaves a fresh empty `[Unreleased]` section
   above it, and opens a PR (branch `release-v1.1.0`, labeled `release`) with
   that change. It fails if `[Unreleased]` is empty -- nothing to release.

2. **Review and merge that PR** like any other. This is the only manual
   review gate in the whole process -- check the CHANGELOG reads well before
   merging, since everything after this step is automatic.

3. **Merging triggers the rest automatically**
   ([release-tag.yml](.github/workflows/release-tag.yml)):
   - Tags the merge commit `v1.1.0`.
   - Creates a GitHub Release from that version's CHANGELOG section.
   - That Release's `published` event is what triggers
     [python-publish.yml](.github/workflows/python-publish.yml), which builds
     and publishes to PyPI.

No git tag is created by hand, no version string is edited anywhere (the
package version is derived from the tag via `setuptools-scm`), and there's no
manual "Draft a new release" click on GitHub.

## One-time repository setup

Both `release-prepare.yml` and `release-tag.yml` push/create things with a
**`RELEASE_TOKEN` secret** (a Personal Access Token), not the default
`GITHUB_TOKEN`. This isn't a style choice -- it's required for the automation
to actually work, because of two GitHub Actions behaviors:

1. **GitHub does not let a `GITHUB_TOKEN`-authored event trigger further
   workflow runs** (anti-recursion). A GitHub Release created with the
   default token would never fire `python-publish.yml`'s `release: published`
   trigger -- the publish step would just silently never run. (We initially
   worked around this by having `release-tag.yml` directly invoke
   `python-publish.yml` as a reusable workflow instead of relying on that
   event. That produced a *different* failure: PyPI's Trusted Publisher is
   registered for the exact file `python-publish.yml`, but a workflow invoked
   this way reports the *calling* workflow's identity to PyPI's OIDC
   verification, not the reusable one -- so PyPI rejected the upload as an
   unrecognized publisher. A real token sidesteps both problems at once: the
   release is a normal, non-recursion-limited event, and `python-publish.yml`
   runs as itself.)
2. **GitHub does not run `pull_request`-triggered workflows** (the test
   matrix, `changelog-check.yml`) **on a PR opened by `GITHUB_TOKEN`**. Using
   a real token for `release-prepare.yml` means the release PR actually gets
   CI checks, like any other PR.

To set this up:

1. Create a [Personal Access Token](https://github.com/settings/tokens) --
   either a classic token with the `repo` scope, or a fine-grained token
   scoped to this repository with **Contents: Read and write** and **Pull
   requests: Read and write** permissions.
2. Add it as a repository secret named `RELEASE_TOKEN`: **Settings > Secrets
   and variables > Actions > New repository secret**.

Without this secret, both workflows fail fast with a clear error rather than
silently doing the wrong thing.

## Retrying a failed PyPI publish

If the tag and GitHub Release already exist but the PyPI upload failed (check
the failed run's logs under the repo's Actions tab), nothing needs to be
re-tagged or re-released -- **Actions -> "Upload Python Package" -> Run
workflow**, entering the existing tag (e.g. `v1.1.0`), rebuilds from that
exact tag and retries the publish.
([python-publish.yml](.github/workflows/python-publish.yml)'s
`workflow_dispatch` trigger.)

## Bioconda

Bioconda's `bioconda-recipes` repo has its own bot that watches PyPI for new
stable releases and opens a version-bump PR there automatically -- nothing in
this repo needs to do anything further once the PyPI publish above succeeds.
