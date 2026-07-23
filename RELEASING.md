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
   - Publishes to PyPI, by directly invoking
     [python-publish.yml](.github/workflows/python-publish.yml) as a reusable
     workflow (**not** by relying on the `release: published` event it also
     listens for -- see the note in that file for why that event would never
     fire here).

No git tag is created by hand, no version string is edited anywhere (the
package version is derived from the tag via `setuptools-scm`), and there's no
manual "Draft a new release" click on GitHub.

## One-time repository setup

`release-prepare.yml` pushes a branch and opens a PR using the default
`GITHUB_TOKEN`. This requires, under **Settings > Actions > General >
Workflow permissions**:

- "Read and write permissions"
- "Allow GitHub Actions to create and approve pull requests"

Without these, step 1 above will fail to push the branch or open the PR.

## Bioconda

Bioconda's `bioconda-recipes` repo has its own bot that watches PyPI for new
stable releases and opens a version-bump PR there automatically -- nothing in
this repo needs to do anything further once the PyPI publish above succeeds.
