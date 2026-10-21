# CLAUDE.md

wasteline is a Nextflow DSL2 (nf-core style) wastewater SARS-CoV-2 pipeline for Illumina amplicon data, ported from
`charlesfoster/covid-illumina-snakemake` (**wastewater workflow only**; the Snakemake "routine"/clinical workflow is
out of scope). Authoritative design, decisions and phase gates: `docs/PLAN.md`. Read it before changing anything.

## Delegation

- Opus plans, makes judgement calls, reviews, and does the high-stakes work itself.
- Hand implementation of an agreed plan to `coder` (Sonnet; see `.claude/agents/coder.md`).
- Fan out `scout` (Haiku; see `.claude/agents/scout.md`) in parallel for bulk read-only work (one per file, sample or shard).
- Review everything `coder` and `scout` return before relying on it.

## Dev environment

```bash
pixi install                  # nf-core, nf-test, pre-commit, pytest, ruff, jinja2, pandas
pixi run nf-test test --profile docker --verbose
pixi run pytest tests/python
pixi run pre-commit run --all-files
pixi run nf-core pipelines lint
```

Develop on macOS + Docker; production is Linux + Docker. Keep `docker`, `singularity`, `conda`/`mamba` profiles working.
Docker Desktop must be running for container tests.

## Conventions

- Process names: nf-core modules `UPPERCASE`; local modules `UPPERCASE_WITH_UNDERSCORES`. Params: `snake_case`.
- Labels: `process_single`, `process_low`, `process_medium`, `process_high`.
- Pass pipeline params via CLI flags or `-params-file`, not custom `-c` files (nf-schema validates params).
- Every local module: `main.nf`, `environment.yml`, a `stub:` block, version reporting via `topic: versions` tuples (as the
  installed nf-core modules do; collated into `pipeline_info/software_versions.yml`), and an nf-test.
- Every container is pinned to an exact tag. Never `latest`. Freyja is the one tool bumped automatically (see PLAN §9).
- Sample IDs are never rewritten (the Snakemake pipeline turned `_` into `-`; we do not).
- A failed/skipped sample must produce a row with a status + reason; it must not abort the run or yield an empty file
  that looks like a result.
- No hard-coded `NC_045512` / SARS-CoV-2 strings outside `conf/pathogens/sars2.config` and `assets/`.

## Non-obvious gotchas

- Freyja >= 2.0 barcodes are `usher_barcodes.feather` (the Snakemake pipeline used 1.x `usher_barcodes.csv`). Freyja also
  accepts a custom `.csv`/`.feather` via `--barcodes`.
- `freyja variants` runs `samtools mpileup | ivar variants` itself and writes both the variants TSV and depths file.
  ivar may append `.tsv` to the prefix; do not assume the output name, verify it.
- `freyja aggregate` takes a directory **with trailing slash** and labels samples by the variants filename recorded in
  the demix file (e.g. `S1.variants.tsv`); strip that suffix afterwards in our own code (`clean_aggregate_labels.py`).
- `freyja plot` is broken in Freyja 2.0.3-2.0.5 containers (pandas 3 `TypeError`, exits 0 with no PDF). See PLAN D6 for
  the workaround and the "write a failure marker" fallback.
- Containers are `linux/amd64` only; Apple Silicon needs `--platform linux/amd64` (handled in a profile, not by hand).
- The Snakemake pipeline used `freyja variants --varthresh 0.25`, discarding low-frequency variants. Freyja's default is
  0.0 and that is our default. Do not "restore parity" here.
- Nextflow >= 26.04's strict (v2) parser passes CLI booleans as strings (`--x`, `--x=false`) and nf-schema rejects them;
  `-params-file` works, as does `NXF_SYNTAX_PARSER=v1`. Named entries (`-entry`) are rejected too: dispatch is `--mode` only.
- `nf-test test --profile docker,test` (not `docker` alone, until resource defaults are fixed).
- On macOS use `-profile docker`. `-profile mamba` also works there (verified end to end, outputs byte-identical to docker): use mamba, never
  plain conda, with `CONDA_SUBDIR=osx-64` (most bioconda packages have no arm64 build) and `NXF_CONDA_CACHEDIR` in a scratch directory.
- Do not commit or push unless the user asks.
- Never delete directories wholesale or clean up "empty" directories (`find -empty -delete`, `rmdir` sweeps, `git clean`).
  Empty-directory sweeps once deleted `.git/objects` and `.git/refs` in this repo. Remove only specific paths you created.
