# charlesfoster/wasteline: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v0.1.0dev - Unreleased

Initial release of charlesfoster/wasteline: a Nextflow DSL2 (nf-core style) port of the wastewater workflow of
[covid-illumina-snakemake](https://github.com/charlesfoster/covid-illumina-snakemake) for Illumina amplicon SARS-CoV-2
data. The Snakemake "routine"/clinical workflow is out of scope.

### `Added`

- Low-confidence estimates (PLAN D18): `--demix_floor_coverage` (default 25). Samples with at least that % of the genome at
  `--freyja_covcut` but below `--min_genome_coverage_for_demix` are demixed and reported as `low_confidence` (a `.demix` exists)
  instead of `skipped_low_coverage`. They appear in `freyja_aggregated.tsv`, the lineage summary (status and reason shown) and the
  report (badge and `demix_low_confidence` info flag), but are left out of the time series, growth rates, dashboard and site
  aggregation unless `--include_low_confidence true`. `depth_stats.py` gains `--floor-genome-coverage` (exit code 11),
  `FILTER_AGGREGATED` / `bin/filter_aggregated.py` filter the aggregate for the time-based analyses and `site_aggregate.py` gains
  `--include-low-confidence`. `--demix_floor_coverage 0` restores the previous behaviour.
- Group-level bootstrap intervals (PLAN D19): with `--freyja_boot` (and the lineage summary) `summary/lineage_group_intervals.tsv` gives a 95 % interval for every lineage group and `other` per sample, computed from the raw Freyja bootstrap draws (`freyja boot --rawboots`, summed into groups draw by draw). `lineage_summary.csv` is unchanged. New module `BOOT_GROUP_INTERVALS`, script `bin/boot_group_intervals.py`; `samples/<s>/freyja/<s>.boot_lineages_boot.csv` holds the raw draws.
- Report privacy: `--report_anonymize none|paths|full` removes file paths, user names (`paths`) and also sample IDs and site names (`full`) from
  `report_data.json` and `report.html`; `full` writes the private `summary/sample_id_map.csv`. In mode `report` the report data are rebuilt
  from the summary files when it is not `none`.
- `summary/lineage_summary.xlsx`: the lineage summary as a formatted Excel workbook (`--lineage_summary_xlsx`, `--lineage_summary_xlsx_decimals`);
  `bin/xlsx_writer.py` converts any other table.
- Primer scheme check: `--primer_check_min_hit_fraction` (0.05) so that a sparse scheme (`eden`) can no longer pass on its enrichment ratio alone;
  `primer_scheme_summary.tsv` gains `declared_hit_fraction` and `best_hit_fraction`; new reason `declared_sparse_inflated`.
- Dropout impact now uses **defining** sites (what separates a lineage from its parent) instead of any informative site, so the flag no longer
  fires on nearly every dropped amplicon; `--dropout_min_affect_abundance` (0.05); new columns `n_defining_sites`, `n_differential_sites`,
  `defining_sites_total`, `defining_sites_lost`, `pct_defining_sites_lost`.
- Unexplained mutations: p-values are computed in log space and floored at `1e-300` (never `0`); new last column `in_dataset`; `--max-rows` rejects
  negative values.
- `bin/reanalyse_runs.py` and `docs/reanalysis.md`: re-analyse past runs with newer Freyja barcodes without re-mapping and see what changed
  (`summary/reanalysis_diff.csv`).
- Lineage groups: `summary/lineage_summary.csv` with one row per sample and one column per user-defined group of lineages
  (patterns `JN.1*`, `NB.1.8.1`, `RECOMBINANT*`, with `exclude` and `in_total`), resolved against the lineage hierarchy of the Freyja
  dataset. On by default with the bundled list `assets/lineage_groups.csv`; `--lineage_groups my.csv` replaces it and
  `--lineage_summary false` turns it off. Also published: `lineage_groups_used.csv` (byte-identical copy of the list used),
  `lineage_groups_resolved.tsv` (audit) and `lineage_ungrouped_summary.tsv` (lineages in no group, with a suggested pattern);
  `run_info.json` records `lineage_groups_file`, `lineage_groups_source` and `lineage_groups_sha256`. Skipped with one warning for
  a custom barcodes file (no `lineages.yml`). `bin/lineage_groups.py --check-groups` and `--print-default-groups`.
- Site aggregation: `--site_aggregation sites.csv` combines per-site Freyja results into user-defined groups of sites by period
  (`--site_aggregation_period week|day`) with a weighted average of the results (`--site_aggregation_weights auto|viral_load|equal|site_weight`,
  `--site_aggregation_min_sites`); reads are never pooled. `summary/site_aggregates_freyja.tsv`, `site_aggregates_meta.tsv` and, with the
  lineage summary, `lineage_summary_site_aggregates.csv`.
- Mistakes in the groups or sites CSV stop the run within minutes (`LINEAGE_GROUPS_CHECK`, `SITE_AGGREGATE_CHECK`, first processes of
  the run) with a message naming the row or line.
- Quiet-error checks: a primer-scheme match check on the read 5' ends (`primer_scheme_check`, `primer_scheme_best`,
  `primer_scheme_enrichment`, `primer_trimmed_pct` in `qc_summary.tsv`, now **27 columns**; WARN reason `primer_scheme_mismatch`,
  optional `low_primer_trimming` via `--qc_min_primer_trimmed_pct`; `--primer_scheme_check`), `summary/dropout_impact[_summary].tsv`
  (which dropped amplicons hide lineage-informative sites, `--dropout_impact`), and `summary/unexplained_mutations.tsv`
  (frequent variants the lineage mixture does not predict, `--unexplained_*`). The ivar trim log is now written in full (nf-core module
  patch).
- Surveillance outputs: `summary/freyja_timeseries.tsv` and `freyja_timeseries_lineages.tsv` (all modes, >= 2 collection dates,
  `--timeseries`), opt-in `--freyja_growth` (`freyja relgrowthrate`: `freyja_growth_rates.csv`, or a `.failed.txt` note, and
  `freyja_growth_inputs.json`) and `--freyja_dash` (`freyja_dashboard.html`); optional numeric `viral_load` samplesheet column for
  load-weighted tables and growth rates.
- `--bam_format cram` (reference-based CRAM output), `--aligner minibwa` (experimental; `docs/aligner_evaluation.md`,
  `tests/compare_aligners.sh`) and `bin/replay_run.py` (rebuilds the command of a recorded run from `run_info.json`).
- Report: sections for the new tables (see `docs/report.md`), run flags `primer_scheme_mismatch`, `dropout_affects_lineages`,
  `consensus_failed`, `consensus_warnings`.
- Probe-capture libraries (e.g. Twist): `--protocol capture` (default `amplicon`) skips primer trimming and the amplicon-dropout
  analysis (no `amplicon_coverage.tsv`, `amplicons_total` / `amplicons_dropped` are `NA`, no `amplicon_dropout` QC reason, the report's
  amplicon section says "not applicable") and marks PCR duplicates (`--duplicates auto|mark|none`; `samtools markdup`, flagged never
  removed, so Freyja's pileup ignores them; `samples/<s>/bam/<s>.markdup_stats.txt`). `--skip_primer_trimming` covers already-trimmed
  amplicon reads. `--primer_scheme`/`--primer_bed` are ignored with a warning in these cases. The protocol, `primer_trimming`
  and the effective `duplicates` are recorded in `run_info.json` and the report header. Synthetic capture test data in
  `tests/data_capture/` (`tests/make_test_data.py --protocol capture`).
- Five run modes selected with `--mode`: `full` (FASTQ to results), `demix` (re-run `freyja demix` on existing variants
  and depths, e.g. with newer barcodes, from a samplesheet or `--previous_results`, including legacy Snakemake results),
  `aggregate` (aggregate existing `.demix` files), `report` (rebuild the HTML report of a finished run) and `update_db`
  (download the latest Freyja data).
- Read processing: fastp (parameters exposed as `--fastp_*`), `bwa mem`, `ivar trim` with nine bundled primer schemes
  (`midnight`, `midnightv2`, `articv3`, `articv4.1`, `articv5.3.2`, `eden`, `swift`, `qiagenv1`, `xgen`; all use the contig
  `NC_045512.2` of the bundled reference) or `--primer_bed` (its contig must match the reference, checked at start-up), and
  primer-trimmed BAM output (`--save_bam`).
- Freyja lineage deconvolution: `freyja variants`, `freyja demix` per sample, `freyja aggregate` with exact sample labels,
  `freyja plot` and the time-aware plot, with all main knobs exposed as `--freyja_*`.
- Freyja data handled as a first-class input: `--barcodes` takes `latest`, `container` (offline), a stored version
  (`YYYY-MM-DD`, `MM_DD_YYYY-HH-MM`, `YYYY-MM-DDTHH-MM`) or a dataset directory / `.csv` / `.feather` file. The dataset is
  resolved once so every sample uses the same data, downloaded datasets are kept in `--db_dir`, and the barcode version
  and checksums are recorded everywhere (`demix_status.tsv`, `barcode_version.json`, `run_info.json`, the report).
- A coverage gate before demixing (`--min_genome_coverage_for_demix`) and a per-sample demix status (`ok`, `failed`,
  `skipped_low_coverage`, `missing`) with reasons. A failed or skipped sample never aborts the run.
- Coverage QC computed from the Freyja depths file, `samtools flagstat` and the fastp JSON (`qc_summary.tsv`, PASS/WARN/FAIL
  with reason tokens; the published table is complete, including the demix, consensus and barcode columns), per-amplicon coverage and dropout (`amplicon_coverage.tsv`), and negative-control assessment
  (`sample_type=negative`, `--ntc_max_mapped_reads`, `--ntc_max_genome_coverage_pct`).
- Per-sample consensus genomes built from the Freyja variants (no second variant caller): depth masking (`N`), IUPAC
  ambiguity codes for mixed sites, `-` for deletions, a status per consensus (`ok`, `empty`, `failed`), and
  `consensus_genomes.fasta` / `consensus_stats.tsv`; tunable with `--consensus_*` and `--indel_freq`, or off with
  `--skip_consensus`.
- Opt-in `--mutations_of_interest` (SNVs, indels and amino-acid changes) giving per-sample allele frequency and depth.
- Opt-in `--compare_barcodes <latest|container|date|stamp|path>` (modes `full` and `demix`): every sample is demixed a second
  time with another Freyja dataset; `summary/barcode_compare.tsv` and `barcode_compare_status.tsv`, a `freyja_data_compare`
  block in `run_info.json`, and a barcode-comparison panel in the report. Ignored with a warning in the other modes.
- Opt-in `--freyja_boot` (with `--freyja_boot_n`, default 100, and `--freyja_boot_seed`, default 0): `freyja boot` bootstrap
  95% intervals (2.5/97.5 percentiles) per lineage in `summary/freyja_boot_summary.tsv`, shown in the report tooltips.
- Opt-in `--annotate_variants`: amino-acid consequences of the variants applied to the consensus (`bcftools csq` with the GFF3),
  `samples/<s>/variants/<s>.annotated.tsv` and `summary/annotated_variants.tsv`; refused at start-up with `--skip_consensus`.
- Negative controls are demixed when they pass the coverage gate but kept out of the aggregate unless
  `--include_negatives_in_aggregate`; their lineages are in `summary/freyja_aggregated_controls.tsv` and a **Negative controls**
  report section, and a control that demixed `ok` stays `ok` in the report.
- The time-series plot and report panel also work in `aggregate` and `demix` modes (collection dates come from the samplesheet).
- `--pathogen` is validated (only `SARS-CoV-2` is accepted) and passed to `freyja update`, `demix`, `boot` and `plot`.
- A self-contained, offline HTML report (`summary/report.html`) built from `summary/report_data.json`: run-level flags,
  QC table, lineage abundance (summarized groups and all lineages), coverage tracks, amplicon dropout heatmap, demix
  status, mutations of interest, time series per site, consensus table and provenance.
- Full provenance in `pipeline_info/run_info.json` (pipeline and tool versions, git commit, parameters, input checksums,
  per-process containers, and the `freyja_data` block) and `-profile locked` with a generated container digest lock
  (`conf/containers.lock.config`, `bin/lock_containers.py`).
- Helper scripts: `bin/make_samplesheet.py` (from FASTQ directories or existing results), `bin/check_samplesheet.py`.
- A deterministic synthetic test dataset (`tests/data`, `tests/make_test_data.py`), pytest unit tests (about 440), nf-test tests
  for every local module and mode (stub and Docker), an upstream canary, and end-to-end and demix-equivalence scripts.
- CI: lint, pytest, nf-test and a daily Freyja release check that opens a bump pull request and auto-merges same-major
  bumps when CI passes (`.github/workflows/ci.yml`, `.github/workflows/freyja-update.yml`, `bin/bump_freyja.py`).
- Documentation: usage, output, parameters, report (with screenshots), development, CI, QC scripts, report data schema, test
  data, migration from Snakemake, troubleshooting, a quick start, and a status page mapping the plan to the implementation
  (`docs/PLAN_STATUS.md`).

### `Changed`

- **QC defaults for wastewater** (calibrated on 28 real QIAseq samples): `--qc_pass_min_pct` 80 -> 60, `--qc_fail_below_pct` 50 -> 30,
  `--min_genome_coverage_for_demix` 60 -> 50 (now described as the high-confidence threshold). The old values remain available as options.
- **The default primer scheme is now `qiagenv1`** (QIAseq Direct SARS-CoV-2) in mode `full`; the Snakemake pipeline defaulted to
  `midnight`. Give `--primer_scheme` explicitly if you use another protocol (the test profiles do). The effective scheme is
  recorded in `run_info.json`.
- `consensus_stats.tsv` has a 7th column, `reason`; `consensus_status.json` has a `warnings` array.
- Compared with the Snakemake pipeline: `freyja variants` keeps Freyja's default `--varthresh 0.0` instead of 0.25;
  sample IDs are never rewritten (`_` is no longer turned into `-`); demix failures are recorded instead of swallowed
  with `|| touch`; Freyja barcodes are a per-run, versioned, checksummed input instead of a mutable shared dataset;
  negative controls are judged from mapped reads and genome coverage rather than kraken2. See
  `docs/migration_from_snakemake.md`.
- Requires Nextflow >= 26.04.0 (strict syntax parser): modes are dispatched with `--mode`, and boolean/numeric
  parameters are typed in `main.nf`.

### `Fixed`

- Consensus genomes of samples with an ivar deletion: ivar prints deleted bases as `N` (`-NNNNNNNNN`), so VCF alleles built from
  them did not match the reference and `bcftools consensus` aborted, giving an all-`N` consensus. Alleles are now built from the
  reference; a variant that still mismatches is skipped (reported in `reason` / `warnings`) instead of failing the sample.

### `Dependencies`

- Freyja `2.0.5` (`quay.io/biocontainers/freyja:2.0.5--pyhdfd78af_0`), the one tool bumped automatically.
- fastp `1.3.6`, bwa `0.7.19`, samtools `1.24`, ivar `1.4.4`, bcftools `1.23`, Python `3.12.12`; the report renderer runs in a Seqera Wave image
  (Python 3.12.12 + Jinja2 3.1.6, `community.wave.seqera.io/library/python_jinja2:b2615084ab562e7a`); every container is pinned to an exact tag.
- Nextflow `>=26.04.0` (CI tests `26.04.0` and `26.04.6`), nf-schema `2.5.1`, generated with nf-core tools `4.1.0`.
- Bump Freyja container `2.0.5--pyhdfd78af_0` -> `2.0.6--pyhdfd78af_0` (automated by `bin/bump_freyja.py`)

### `Deprecated`

- Nothing.

### Known issues

- Only `docker` has been run end to end; the `singularity`/`apptainer` and `conda`/`mamba` profiles are configured but untested,
  and the GitHub workflows have not run yet.
- `--mode demix --previous_results <dir>` lists a `db/` directory inside the results (as left by the `test` profiles) as a
  sample with status `missing` (`docs/troubleshooting.md`).
- Not implemented from the plan: `--kraken2_db` (optional negative-control screening), other pathogens, historical barcode versions (`docs/development.md`, `docs/PLAN_STATUS.md`).
