# In-pipeline multi-target plan (PLAN D25 phase 2, "intermediate plan"): to pick up later on a new branch

Written 2026-10-22 for a later session. Nothing here is built. Start it on a new branch from `main` (suggested name
`feature/multi-target-in-pipeline`; `main` was `fa3e62f` when this was written; this file and `panel_mode_plan.md` are on `dev`).

## Read this first: the non-negotiables

1. **The pipeline is primarily a SARS-CoV-2 wastewater pipeline, and SARS-CoV-2 stays the default.** A command without `--pathogen` or
   `--panel` is a SARS-CoV-2 run, exactly as today. Other pathogens (RSVa, RSVb, influenza, measles, mpox; later enterovirus, HCV) are an
   extension that must never cost SARS users anything: no extra runtime, no extra required options, no changed results.
2. **Results stay numerically identical for SARS-CoV-2** (abundances, coverage and QC values, consensus sequences, statuses). What this plan
   adds is permission to change the LAYOUT (file names, column order, folder structure) in a controlled, staged way (see "Layout staging").
3. **Tests stay minimal** (standing user rule, docs/development.md "Keep the suite small"): no new test files except one small pytest per
   new script; one stub nf-test per new local module (CLAUDE.md requires it); at most ONE new e2e step for the whole feature.
4. **Test and verification RUNTIME stays minimal.** The e2e takes ~14 min on the user's Apple Silicon Mac (emulated containers). Rules: no
   bootstraps in any new step (`--freyja_boot` off; the two existing bootstrap steps are the only ones; step 1 needs 20 draws, the
   group-intervals step uses 5); new synthetic data tiny (a few hundred read pairs, small genomes); never run the full e2e more than once per
   merged unit of work; run targeted manual checks instead; remove stale `.claude/worktrees` before nf-test (nf-test scans them).
5. Standing repo rules: sample IDs never rewritten; a failed/absent sample/target is a row with a status and reason, never an empty file;
   containers pinned to exact tags; no pathogen names in pipeline code outside `conf/pathogens/` and `assets/`; use mamba (never plain conda)
   with `CONDA_SUBDIR=osx-64` on macOS; do not commit/push unless asked.

## Why this plan exists (what changed after `panel_mode_plan.md`)

`panel_mode_plan.md` recommends architecture B: a typing + partition front end and the UNCHANGED pipeline once per pathogen, because it leaves
SARS byte-identical. Its costs: reads mapped and trimmed twice, one pipeline launch per detected pathogen (fixed overhead each), no single
combined report, and a driver script. The user asked what would be done differently if byte-identity were not required. Answer: build the
in-pipeline design (architecture A). The compromise here ("intermediate plan") keeps the SARS RESULTS identical but allows layout changes,
staged so the SARS default never breaks without notice. Decide between B and this plan before starting (see "Decisions").

## Concepts

- **Target**: one thing measured in a sample: a pathogen, subtype, genotype or segmented genome (one or more contigs; one `freyja_contig`
  carries the Freyja barcodes; may have none).
- **Panel**: a table of targets. **A normal run is a panel of one**: `--pathogen X` is shorthand for "the built-in row called X", and the default
  row is SARS-CoV-2. Columns (CSV): `target,fasta,gff,freyja_pathogen,freyja_contig,marker_contig,family,primer_scheme,barcodes,lineage_groups`.
  Today these facts live in `conf/pathogens/*.config` (11 files) plus hidden "profile facts"; they become rows of one registry
  (`assets/panels/registry.csv` or the profile files generated from it), read per sample from `meta.target` instead of global `params.*`.
- **Reference selection and mixed-sample calls** (non-SARS panels only; a one-row SARS panel never selects): competitive `bwa mem` against the
  concatenated panel; adapt the user's own h2seq `bin/select_reference_from_bam.py` (per-reference stats, ranking, `>= 5 %` secondary =
  mixed, `qc_fail` rows) and `bin/filter_mixed_reference_alignments.awk` (fragment counts for a target only if both mates agree and MAPQ >= 10,
  else ambiguous); group references by an explicit `family` column, never a name prefix; type on `marker_contig` (HA for flu; internal segments
  cross-map between subtypes). Statuses per (sample, target): `detected`, `trace`, `not_detected`; thresholds as params (starting points
  100 fragments, 5 % breadth, 5 % family share; UNCALIBRATED, calibrate on real data). Full statistics and decision rules:
  `docs/wiring/panel_mode_plan.md` (reuse as is). Phase 1 (multi-contig targets, `freyja_contig`, `FREYJA_SELECT_CONTIG`) is done and merged.

## Data model change (the invasive part)

Counts from the planner (grep on `dev` at `fa3e62f`): `params.pathogen` 26 sites, `params.fasta` 21, `params.gff` 4, `params.freyja_contig`
14, primer/BED logic ~78 lines, `FREYJA_DATA`/dataset root 23 lines, publishDir sites 76 in `conf/modules.config`, 57 process calls in
`workflows/wasteline.nf`, 113 `meta.id` uses in modules, six config-time hidden facts.
- `meta.target` (a map with the row's facts) is added next to `meta.id`; per-sample process closures read `meta.target.*`.
- `meta.id` stays the sample ID in files shown to users. Internally work names/joins/collectFile keys use `meta.id` + `meta.target.name`
  (keep a `meta.key`), so aggregate, QC, status files and `carried_low[...]` never mix targets.
- Run-level steps (aggregate, plot, time series, site aggregation, lineage summary, unexplained mutations, report, `demix_status.tsv`) are
  called on `groupTuple(by: target)` channels; one Freyja dataset per target (`FREYJA_DATA` keyed by `freyja_pathogen`).
- Config-time hidden facts (`lineage_groups_bundled`, `primer_schemes_bundled`, `freyja_container_pathogen`, ...) move into the registry row.
- With a single target every channel degenerates to today's shape (this is how SARS stays identical in stage S1).

## Layout staging (keeps SARS safe)

- **S1: internal refactor, outputs byte-identical.** Introduce `meta.target`, the registry and the per-target grouping; a one-row (SARS)
  panel produces EXACTLY today's files (the e2e byte-identity diff against a baseline is the guard, as in D25 phase 1). Nothing user-visible
  changes. This is the largest piece and the one to do first.
- **S2: multi-target runs and the new layout under a flag.** `--panel` with several targets (and `--layout targets`) writes
  `samples/<id>/<target>/...`, `summary/<target>/...`, `summary/targets.tsv`. Single-target runs keep the legacy layout unless `--layout targets`.
  Panel mapping/selection modules and `--mode panel`-style typing-only output (`targets.tsv` without launching analysis) live here.
- **S3: reports.** Per-target report (essentially today's report, per target) PLUS a global overview page (sample x target matrix with status,
  fragments, breadth, top lineage, run-level flags, links), `--report_mode per_target|combined|both` (default `both`). One combined
  self-contained HTML is optional, not required (size grows with targets; per-target files keep the SARS report as is).
- **S4 (only if the user wants): flip the default layout** to the unified one, with a compatibility reader in `reanalyse_runs.py`,
  `runs_store.py`, `make_samplesheet.py --from-results` and `results_input` for old result folders, and a migration note. Not before the user
  has used S2/S3 on real data.

## Change list outline (details to be made by a planning pass at the start of the branch)

New: `bin/panel_select.py` (prepare/call), modules `panel_prepare`, `panel_assign`, `panel_call`, `panel_partition` (or in-pipeline split by
target without writing partitioned FASTQs: since mapping is shared, split the BAM per target with `samtools view` of the target's contigs;
this removes the double mapping/trimming of plan B), `assets/panels/*.csv`, `bin/compare_results.py` only if a layout-insensitive comparison
is needed for S2 (small: compares values, not bytes). Edited: `main.nf`, `nextflow.config`, schema, utils subworkflow, `workflows/wasteline.nf`
(the bulk), `conf/modules.config`, `conf/pathogens/*`, `results_input`, report builders (S3). Unchanged in S1 outputs: everything user-visible.

## Verification plan (minimal)

- S1: full e2e once per merged unit; `diff -r` of every shared output directory against the previous `dev` baseline (excluding timestamps,
  `pipeline_info/`, report provenance keys) = no real differences; zero nf-test snapshot diffs; identical `run_info.json` params.
- S2/S3: ONE new e2e step: a tiny two-target synthetic panel without Freyja barcodes network access if possible (RSVa + RSVb pseudo-genome
  reads from `tests/simulate_pathogen_mixture.py`, ~200 pairs per sample; S1 pure RSVa, S2 pure RSVb, S3 70/30 mix, S4 a trace, one negative
  control; `--protocol capture`; NO bootstraps) asserting `targets.tsv` statuses and mixed fractions, per-target outputs exist, the global
  page lists every sample x target. Plus one pytest per new script. Manual (documented, not committed): a real whole-genome capture check when
  the user has data.
- Per-target report content is covered by the existing report tests; add one assertion for the matrix.

## Risks

RSV-A/B cross-mapping inflating a trace to `detected` (floors, measured with the trace sample); MAPQ ties among same-target candidate
references in candidate-rich panels (HCV/enterovirus; consider the XA tag later); thinly covered flu internal segments in mixed samples (WARN);
join/collectFile keys that forget `meta.target` and silently mix targets (review every `join`/`collectFile`/`groupTuple`); config-time
facts that cannot be per target (move into the registry); demix-mode and `--previous_results` reading old and new layouts; effort: roughly a few
times the work of plan B; the e2e byte-identity guard does not cover multi-target behaviour (only the new step does).

## Decisions pending (ask the user at the start)

1. Architecture: this in-pipeline plan (one run, one DAG, no double mapping, layout staging) versus plan B (`panel_mode_plan.md`: front end +
   per-target runs, SARS untouched). 2. Starting thresholds for detection (100 fragments / 5 % breadth / 5 % family share) and calibration data.
3. Should `trace` targets be analysed or only listed? 4. Flu panel contents (H1N1, H3N2, H5Nx whole genome, B/Vic HA only) and whether flu data is
   capture or amplicon (both possible). 5. Is stage S4 (default layout flip) wanted, and when? 6. Phase 3 (targets without Freyja barcodes:
   `--skip_demix`, best-reference call) after the RSV and flu panels work.

## How to start

1. New branch from `main` (check `git log` first). 2. Run the full e2e once on that commit into a scratch dir as the S1 baseline. 3. Get a
fresh read-only planning pass (Plan agent) for S1 only, using this file and `panel_mode_plan.md`, and review it before any coding. 4. Implement
S1 with strict file ownership, prove byte-identity, merge; only then S2. Keep each coder brief's test limits explicit (hard maxima).
