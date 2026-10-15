# CI and automated Freyja updates

This page describes how the CI pieces fit together, what has to be configured on GitHub, and how to run every check
locally. Design background: `docs/PLAN.md` section 9 (Freyja auto-update, D1) and section 10 (digest lock, F13).

## How the pieces fit

```
                   daily / manual
 freyja-update.yml ───────────────► bin/bump_freyja.py --check --json
        │                                   │ update_available? (quay.io + bioconda + Galaxy depot all have it)
        │                                   ▼
        │             bin/bump_freyja.py --apply --tag T   rewrites every Freyja pin + CHANGELOG
        │             bin/lock_containers.py --update      refreshes conf/containers.lock.config
        ▼
 peter-evans/create-pull-request ──► PR  bot/freyja-<version>   (opened with BOT_TOKEN, so CI runs)
        │
        ├─ same major  ─► gh pr merge --auto --squash ─► merges once the required checks pass
        └─ new MAJOR   ─► label `needs-review`, left open for a human

 ci.yml (every PR, push to main, weekly)  ──►  ci-required  (the single required check)
```

Everything that matters is also a plain script, so it runs identically on a laptop (see "Run it locally").

### Workflows

| File                    | Trigger                                  | Purpose                                                                                                                            |
| ----------------------- | ---------------------------------------- | ---------------------------------------------------------------------------------------------------------------------------------- |
| `ci.yml`                | PR, push to `main`, weekly, manual       | The pipeline CI (jobs below).                                                                                                      |
| `freyja-update.yml`     | daily 05:17 UTC, manual                  | Detects a new Freyja release and opens the bump PR.                                                                                |
| `linting.yml`           | PR, release                              | nf-core template: `pre-commit` (prek: prettier, whitespace, `nextflow lint`) and `nf-core` jobs.                                   |
| `nf-test.yml`           | release, manual (never on PRs)           | nf-core template: the extended, sharded nf-test matrix (docker + singularity + conda; Nextflow `26.04.0` and `latest-everything`). |
| `download_pipeline.yml` | manual                                   | nf-core template: `nf-core pipelines download` plus an offline stub run. Run it before a release.                                  |
| `branch.yml`            | PR to `main`                             | nf-core template: only `dev`/`patch` may PR into `main` (see "Known template friction").                                           |
| others                  | `pr-comment.yml`, `fix_linting.yml`, ... | nf-core template helpers, unchanged.                                                                                               |

### Jobs in `ci.yml`

| Job (check name)  | What                                                                                                                                                                                 | Blocking                                                  |
| ----------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ | --------------------------------------------------------- |
| `ruff`            | `ruff check` and `ruff format --check` on `bin` and `tests`                                                                                                                          | yes (via `ci-required`)                                   |
| `pytest`          | `pytest tests/python`, including the invariant "every Freyja pin agrees"                                                                                                             | yes (via `ci-required`)                                   |
| `nf-core-lint`    | `nf-core pipelines lint` (never `--release`), on Nextflow `26.04.6`                                                                                                                  | yes (via `ci-required`)                                   |
| `nf-test-stub`    | `nf-test test --tag stub --profile test_stub`, on Nextflow `26.04.0` and `26.04.6`; fails if no test carries the tag                                                                 | yes (via `ci-required`)                                   |
| `nf-test-docker`  | `nf-test test --tag docker --profile=+docker`, on Nextflow `26.04.0` and `26.04.6`; only when a PR touches a Freyja pin, on `main`, weekly, manual. A skipped job counts as passing. | yes (via `ci-required`)                                   |
| `upstream-canary` | `nf-test test --tag upstream_canary` (is plain `freyja plot` still broken upstream, PLAN D6?)                                                                                        | **no** (`continue-on-error`); same trigger rule as docker |
| `container-lock`  | `bin/lock_containers.py --check`; weekly/manual also `--verify-digests`                                                                                                              | **no** (`continue-on-error`)                              |
| `ci-required`     | aggregates the blocking jobs above; fails if any failed or was cancelled                                                                                                             | **the required check**                                    |

"Touches a Freyja pin" means the PR diff adds or removes a line matching `biocontainers/freyja:`, `singularity/freyja:`
or `bioconda::freyja=`, or the head branch starts with `bot/freyja-`.

### nf-test tag conventions (needed by `ci.yml`)

`nf-test` can only _include_ tags (no `not`), so the jobs select by positive tags. Set them at test level
(`test("...") { tag "stub" ... }`), as `tests/default.nf.test` already does for `docker`:

| Tag               | Put it on                                     |
| ----------------- | --------------------------------------------- |
| `stub`            | every test that runs with `options "-stub"`   |
| `docker`          | every real test that needs a container engine |
| `upstream_canary` | `tests/canary/*` (already set)                |

Both required nf-test jobs run a `--dry-run` first and fail with an explicit error if their tag selects nothing, so a
missing tag can never turn a required check into a silent no-op.

## What you must configure on GitHub

These are repository settings and secrets. The workflows cannot set them for you.

### 1. Secret `BOT_TOKEN` (Settings > Secrets and variables > Actions > New repository secret)

Pull requests opened with the default `GITHUB_TOKEN` do **not** trigger other workflows (a GitHub safeguard against
recursion), so the bump PR would sit there with no CI and could never auto-merge. `freyja-update.yml` therefore opens
the PR, creates the labels and enables auto-merge with `BOT_TOKEN`. Use one of:

- **Fine-grained personal access token** (simplest). Resource owner: you; repository access: only
  `charlesfoster/wasteline`; repository permissions: _Contents_ read and write, _Pull requests_ read and write,
  _Issues_ read and write (labels), _Metadata_ read. Set an expiry and put a reminder in your calendar: when it
  expires the daily job fails loudly with "BOT_TOKEN secret is not set" or a 401.
- **GitHub App** (better for an organisation; no expiry, not tied to a person). Create an App with the same four
  permissions, install it on the repository, and mint a token in the workflow with `actions/create-github-app-token`
  instead of reading `secrets.BOT_TOKEN` (replace the three `secrets.BOT_TOKEN` uses in `freyja-update.yml`).

The job fails (rather than skipping quietly) when a newer Freyja exists but `BOT_TOKEN` is unset.

### 2. Settings > General > Pull Requests: tick "Allow auto-merge"

Without it `gh pr merge --auto` fails and the bump PR stays open.

### 3. Branch protection (or a ruleset) on `main`: this is what makes auto-merge safe

Auto-merge merges the PR **only when branch protection says it may**. With no required checks, GitHub would consider the
PR immediately mergeable. Settings > Branches > Add branch protection rule for `main`:

- Require a pull request before merging (required approvals: 0 is fine; the bot cannot approve its own PR, and a
  required approval would block every auto-merge).
- Require status checks to pass before merging, and add exactly these **required checks**:
  - `ci-required` (from `ci.yml`: ruff, pytest, nf-core lint, nf-test stub, nf-test docker)
  - `pre-commit` (from `linting.yml`: prettier, whitespace, `nextflow lint`)
- Do **not** require `nf-core` (the `linting.yml` job): on PRs to `main` it runs `nf-core pipelines lint --release`,
  which rejects the `0.1.0dev` version, so it would fail every PR. `ci.yml`'s `nf-core-lint` runs the same linter
  without `--release` and is already part of `ci-required`.
- Do not require `branch`, `Run nf-test` jobs, `upstream-canary` or `container-lock`.
- Optional but recommended: "Do not allow bypassing the above settings", so that not even an admin token can merge past
  failing checks (and the bot's `BOT_TOKEN` owner is not exempt).
- A `merge queue` is not required.

A new **major** Freyja version never gets auto-merge; its PR carries the `needs-review` label and waits for you.

### 4. Optional: Actions secrets/variables

None other. The third-party actions are pinned to full commit SHAs (table below); no other secret is read by any
workflow, and the workflows contain no secrets.

## Run it locally

```bash
pixi install

# lint and unit tests (what ruff + pytest + part of pre-commit run in CI)
pixi run ruff check bin tests
pixi run ruff format --check bin tests
pixi run pytest tests/python
pixi run pre-commit run --all-files
pixi run nf-core pipelines lint              # needs a git repo with a commit; lint a committed copy (it rewrites some files), see docs/development.md

# nf-test, same selections as CI (Docker Desktop must be running for the docker ones)
pixi run nf-test test --tag stub --profile test_stub --ci --verbose
pixi run nf-test test --tag docker --profile=+docker --ci --verbose
pixi run nf-test test --tag upstream_canary --profile=+docker --verbose

# Freyja updates (network; stdlib Python only, no pixi needed)
python bin/bump_freyja.py --require-consistent        # offline invariant: all pins agree
python bin/bump_freyja.py --check                     # pinned vs latest
python bin/bump_freyja.py --check --json              # what CI consumes
python bin/bump_freyja.py --apply --dry-run           # what a bump would rewrite
python bin/bump_freyja.py --apply                     # rewrite pins + CHANGELOG (to latest)
python bin/bump_freyja.py --apply --tag 2.0.4--pyhdfd78af_0   # explicit tag, no lookups: also the way to roll back

# digest lock
python bin/lock_containers.py --check                 # offline: is conf/containers.lock.config stale?
python bin/lock_containers.py --check --verify-digests  # also re-resolve every digest against the registries
python bin/lock_containers.py --update                # (re)write the lock; unchanged pins are not re-resolved
python bin/lock_containers.py --update --refresh      # re-resolve everything
```

To exercise the whole bump flow without touching the repository: `cp -R` the tree to a scratch directory, point
`--root` at it, and run `--apply --tag <something>`; the unit tests do exactly that against fixtures.

### `bin/bump_freyja.py`

- `--check` exit codes: `0` up to date (or a newer release exists but is **blocked**, e.g. bioconda or the Galaxy
  singularity depot has not got it yet; it is retried tomorrow), `10` update available, `2` error.
- `--check --json` keys: `status` (`up_to_date`, `update_available`, `blocked`), `update_available`, `pinned`,
  `latest`, `pinned_version`, `latest_version`, `major_bump`, `blocked_reasons`, `checks`.
- "Latest" means the highest numeric version among quay.io `biocontainers/freyja` release tags
  (`<version>--<build>_<n>`), then the highest build number; rc/dev/beta/arch-suffixed/`latest` tags are ignored.
  The same version must exist on bioconda (anaconda.org API) and as
  `https://depot.galaxyproject.org/singularity/freyja:<tag>` (HEAD request; an unreachable depot is reported as
  `skipped`, not as a blocker; `--no-depot-check` skips it).
- `--require-consistent` (exit `1` if not) is the invariant also asserted by
  `tests/python/test_bump_freyja.py::test_real_repo_pins_are_consistent`. It checks that all `quay.io/biocontainers/freyja:`,
  `depot.galaxyproject.org/singularity/freyja:` and `bioconda::freyja=` strings agree, that the conda version equals the container
  version, that every module that pins Freyja pins it in all three places, and that version stamps in tests agree.
- `--apply` finds the pins by grep over the whole tree (no module list is hard-coded). Besides the three pin
  spellings it also rewrites **version stamps in tests**: `"freyja": "x.y.z"` in nf-test snapshots (`*.snap`) and
  `freyja_version == 'x.y.z'` / `freyja: x.y.z` in `*.nf.test`, because those assert the running Freyja version and
  would fail every bump PR. It appends one bullet under the unreleased section's `Dependencies` heading of
  `CHANGELOG.md` (creating the heading or section when absent) and is idempotent.
- Deliberately **not scanned**: `CHANGELOG.md`, `pixi.lock`, `conf/containers.lock.config` (generated),
  `docs/PLAN.md`, `docs/test_data.md` (a record of how the test data was made), `docs/ci.md`, the two scripts, their
  tests and `tests/fixtures/ci/`. Put `bump-freyja: ignore` on any other line that must keep an old pin.

### `bin/lock_containers.py` and `-profile locked`

`conf/containers.lock.config` maps each process to `registry/repo@sha256:<digest>`, resolved read-only from the
registry (`Docker-Content-Digest` of a manifest HEAD; a multi-arch tag records the manifest-list digest and is marked
`(index)`; the trailing `// pin:` comment is the exact tag it was resolved from). The file is generated and
deterministic (no timestamps). `--check` compares the pins found in `modules/**/main.nf` and `conf/*.config`
overrides against the lock without any network; it exits `1` when stale.

Wiring (for `nextflow.config`, inside the existing `profiles { ... }` block):

```groovy
    locked {
        // Immutable image digests for every process (bin/lock_containers.py). Use with docker/singularity:
        //   nextflow run . -profile docker,locked ...
        // Process container overrides only, not params; has no effect with conda/mamba.
        includeConfig 'conf/containers.lock.config'
    }
```

It sets `process.withName` container selectors, which take precedence over the `container` directive in the modules,
so a locked run uses the docker image for every engine (Singularity should then pull the OCI image from the registry
instead of the Galaxy depot `.img`; not tested). It has been checked with Docker: `nextflow config -profile locked` parses and
`docker run` with the digest works. Provenance can read the lock with `lock_containers.parse_lock()` (PLAN 5.13: "all
container tags/digests (if lockfile)").

Pins that cannot be locked make `--update` fail with a clear message: a floating tag (`latest`), a missing tag, or a
module whose `container` expression contains no non-URL image.

## Operating the Freyja bot

1. Every day the job runs `--check`. Nothing newer, or newer but **blocked** (written to the job summary): no PR.
2. Newer and available everywhere: it applies the bump and refreshes the lock, then opens or updates the PR on
   `bot/freyja-<version>` (so re-runs are idempotent). The PR title ends in `(MAJOR)` for a major bump.
3. Required checks run on that PR (they are triggered because the PR is opened with `BOT_TOKEN`). Because the PR
   touches a Freyja pin, `nf-test-docker` runs the real container tests.
4. Same-major: auto-merge squashes it once `ci-required` and `pre-commit` are green. Major: stays open, labelled
   `needs-review`; read the Freyja release notes, then merge manually.
5. If CI is red, auto-merge simply never fires; the PR stays open and the next daily run does not disturb it (it
   re-applies the same bump to the same branch). Fix forward on the branch or close the PR.

Two upstream-sensitive checks to watch after a bump: the `upstream-canary` job turns red when a Freyja release fixes
plain `freyja plot` (remove the D6 workaround then), and `tests/default.nf.test.snap` stamps the Freyja version, which
`bin/bump_freyja.py --apply` rewrites; any other snapshot content that changes with Freyja output needs a human to
regenerate it.

## Third-party actions (looked up from the GitHub API on 2026-10-15, pinned to the commit SHA of the release tag)

| Action                            | Version | Commit SHA                                 |
| --------------------------------- | ------- | ------------------------------------------ |
| `actions/checkout`                | v7.0.1  | `3d3c42e5aac5ba805825da76410c181273ba90b1` |
| `actions/setup-python`            | v7.0.0  | `5fda3b95a4ea91299a34e894583c3862153e4b97` |
| `actions/upload-artifact`         | v7.0.2  | `cf430e030ddbb5b0abf93d22962f4752f3646cd9` |
| `prefix-dev/setup-pixi`           | v0.10.2 | `d3f436a425481402e6a95a1d1fc10331c708cd9e` |
| `peter-evans/create-pull-request` | v8.1.1  | `5f6978faf089d4d20b00c7766989d076bb2fc7f1` |
| `nf-core/setup-nextflow`          | v3.0.1  | `893c28b667aedeba26e37f296d260ccc5bc4d914` |
| `nf-core/setup-nf-test`           | v2.1.0  | `4069fbbaabe94c08faba4ad261bfa88225ba133f` |

`prefix-dev/setup-pixi` v0.10.2 is deliberately not the newest release (v0.11.2 was a day old when this was written).
Tool versions used by CI (`ci.yml`): pixi `v0.67.0` and nf-test `0.9.5` (both as locked in `pixi.lock`). Nextflow: the
minimum is `26.04.0` (`manifest.nextflowVersion = '!>=26.04.0'`); the `nf-test-stub` and `nf-test-docker` jobs run on a matrix
of `26.04.0` (the minimum) and `26.04.6` (latest stable used), and `nf-core-lint` and `upstream-canary` run on `26.04.6`.
The template's extended `nf-test.yml` (release/manual only) uses `26.04.0` and `latest-everything` with nf-test `0.9.4`.

## Known template friction

- `branch.yml` (nf-core) fails any PR to `main` whose head is not `dev` or `patch`. The repository now uses the
  `dev` -> `main` flow (see docs/development.md), so PRs from `dev` pass; the Freyja bot's `bot/freyja-<version>` PRs
  into `main` still show a red, **non-required** `test` job from it. `nf-core pipelines lint` compares `branch.yml` and
  `linting.yml` against the template (`files_unchanged`), which is why they were not edited here; add them to
  `lint: files_unchanged:` in `.nf-core.yml` before changing them.
- `linting.yml`'s `nf-core` job lints with `--release` on PRs to `main` (see branch protection above), for the same
  reason.
- `nf-test.yml` is the template's heavy matrix. It is kept for releases/manual use; its `pull_request:` trigger only
  exists because `nf-core pipelines lint` requires one (it uses `paths-ignore: ["**"]`, so it never fires).
