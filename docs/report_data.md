# `report_data.json` schema (version 1)

Single machine-readable input of the HTML report (`bin/render_report.py`), written by
`bin/build_report_data.py`, consumed by `--mode report`. UTF-8 JSON (written compactly, on one line). Missing values are
`null` (never the string `"NA"`); numbers are JSON numbers. Only `schema_version`, `mode`, `provenance` and `samples`
are always present; the renderer tolerates any other key being absent. See [report.md](report.md) for how the report
uses each field and [output.md](output.md) for the files it is built from.

```
build_report_data.py --run-info run_info.json [--qc qc_summary.tsv] [--amplicons amplicon_coverage.tsv]
    [--aggregated freyja_aggregated.tsv] [--demix-status demix_status.tsv] [--moi mutations_of_interest.tsv]
    [--consensus consensus_stats.tsv] [--depths 'a.depths b.depths' | DIR] [--barcode-compare TSV]
    [--controls-aggregated freyja_aggregated_controls.tsv] [--samplesheet samplesheet.csv]
    [--boot freyja_boot_summary.tsv] [--compare-status barcode_compare_status.tsv]
    [--mode MODE] --out report_data.json [--write-qc-tsv qc_summary.merged.tsv]
render_report.py --data report_data.json --template assets/templates/report.html.j2 --out report.html
```

The scripts need Python >= 3.11 (the builder imports `datetime.UTC`; the system Python 3.9 of macOS fails with an
`ImportError`) and, for the renderer, `jinja2`. `--depths` takes depths files, one quoted string holding a
space-separated list of files, or a directory (searched recursively); the script does **not** expand wildcards, so
`--depths 'x/*.depths'` in quotes silently gives no coverage tracks: let the shell expand an unquoted glob.

Only `--run-info` and `--out` are required. A path argument that does not exist or is named `NO_FILE*` is treated as
absent (the Nextflow process stages `[]` as no file). Options added with the negative-control, `--freyja_boot` and
`--compare_barcodes` work:

| option                  | file                                     | effect                                                                                                                                                                                                                                                        |
| ----------------------- | ---------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `--controls-aggregated` | `summary/freyja_aggregated_controls.tsv` | lineages of the negative controls that demixed `ok` (they are not in the main aggregate); feeds `negative_controls`                                                                                                                                           |
| `--samplesheet`         | `summary/samplesheet.csv` (CSV or TSV)   | `sample_type`, `collection_date`, `site` for samples that have no QC row (mode `aggregate`), so the time-series panel works without a QC table                                                                                                                |
| `--boot`                | `summary/freyja_boot_summary.tsv`        | 95% bootstrap intervals; feeds `demix.<sample>.intervals` and `available.boot_intervals`                                                                                                                                                                      |
| `--compare-status`      | `summary/barcode_compare_status.tsv`     | outcome of the second demix per sample; only used for the `barcode_compare_failed` run flag                                                                                                                                                                   |
| `--write-qc-tsv`        | output                                   | also writes the QC table with `demix_status`, `demix_reason`, `consensus_length`, `consensus_n_pct` and `barcode_version` filled in (cells the QC script already filled are kept byte for byte). The pipeline publishes this file as `summary/qc_summary.tsv` |

`--barcode-compare` is now fed by the pipeline (`summary/barcode_compare.tsv`, see below).

## Top level

| key                     | type         | notes                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 |
| ----------------------- | ------------ | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `schema_version`        | int          | `1`                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   |
| `generated`             | string       | UTC timestamp `YYYY-MM-DDTHH:MM:SSZ`                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  |
| `mode`                  | string       | `full`, `demix`, `aggregate`, `report` (from `--mode`, else `run_info.mode` / `params.mode`, else `unknown`)                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          |
| `header`                | object       | derived summary used by the report header: `pipeline_name`, `pipeline_version`, `git_commit`, `nextflow_version`, `mode`, `freyja_version`, `barcode_version`, `barcode_version_iso`, `barcode_source`, `barcode_sha256`, `pathogen`, `run_started`, `run_info_generated` (all nullable). Looked up from `run_info.json` under several key spellings (see below)                                                                                                                                                                                                                                                                      |
| `provenance`            | object       | verbatim copy of `run_info.json` as it was when the report data were built (the mandatory `freyja_data` block, `freyja_data_compare` with `--compare_barcodes`, `params`, ...; the completion keys `containers`, `run_completed`, `duration` and `success` are added to `run_info.json` only after the run, so they are not here), plus `previous_runs` (list, default `[]`) of `{path, run_info_sha256, barcode_version}` (`barcode_version: "unknown"` for legacy Snakemake results)                                                                                                                                                |
| `params_summary`        | object       | selected scalar params from `provenance.params`, in a fixed order, only those present: `mode`, `primer_scheme`, `primer_bed`, `barcodes`, `freyja_eps`, `freyja_covcut`, `freyja_depthcutoff`, `freyja_varthresh`, `freyja_adapt`, `freyja_autoadapt`, `freyja_relaxedmrca`, `freyja_solver`, `freyja_confirmedonly`, `min_genome_coverage_for_demix`, `amplicon_min_depth`, `consensus_min_depth`, `consensus_freq`, `indel_freq`, `ntc_max_mapped_reads`, `ntc_max_genome_coverage_pct`, `mutations_of_interest`, `compare_barcodes`, `skip_consensus` (23 in a real run). `freyja_boot` and `annotate_variants` are not part of it |
| `thresholds`            | object       | `amplicon_min_depth` (builder default 10 when absent), `consensus_min_depth` (builder default 10 when the param is absent from `run_info.json`; the **pipeline default is 10**, which is what a real run records), `ntc_max_mapped_reads`, `ntc_max_genome_coverage_pct`, `min_genome_coverage_for_demix` (null when not in params)                                                                                                                                                                                                                                                                                                   |
| `available`             | object       | booleans: `qc`, `amplicons`, `demix`, `coverage_tracks`, `mutations_of_interest`, `consensus`, `barcode_comparison`, `negative_controls`, `boot_intervals`, `time_series`                                                                                                                                                                                                                                                                                                                                                                                                                                                             |
| `run_flags`             | list         | `{level: "error"\|"warn"\|"info", code, message, samples?}`, sorted error, warn, info (see codes below)                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               |
| `samples`               | list         | one object per sample, sorted by `sample` (see below)                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 |
| `demix`                 | object       | sample -> demix result (see below); one key per sample in `samples`                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   |
| `amplicons`             | object       | sample -> list of `{amplicon, start, end, mean_depth, dropped}` sorted by start; `start`/`end` 1-based inclusive. `{}` when absent                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    |
| `coverage_tracks`       | object\|null | `{contig, genome_length, bin_size, n_bins, samples: {sample: [mean depth per bin]}}`; positions 1-based, bin `i` covers `[i*bin_size+1, (i+1)*bin_size]`; `n_bins <= 600`; values are ints (>= 10) or 1-dp floats; positions missing from the depths file count as 0; only the first contig of a depths file is used                                                                                                                                                                                                                                                                                                                  |
| `mutations_of_interest` | list\|null   | `{sample, name, mutation, ref_depth, alt_depth, depth, alt_freq}` rows (`alt_freq` 0-1). `null` = not run; `[]` = run with no rows                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    |
| `consensus`             | list         | `{sample, length, n_count, n_pct, ambiguous_count, status}`                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           |
| `barcode_comparison`    | object\|null | `{label_a, label_b, rows: [{sample, lineage, a, b, delta}]}` (fractions 0-1, 5 dp, `delta = b - a`, sorted by `abs(delta)` desc); `null` when there are no rows (no `--compare_barcodes`). `label_a`/`label_b` are the barcode versions of the two datasets (`custom` for a bare file)                                                                                                                                                                                                                                                                                                                                                |
| `negative_controls`     | list         | one entry per negative control that demixed `ok` or has a row in the controls aggregate (see below); `[]` otherwise                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   |
| `time_series`           | object\|null | present only when >= 2 distinct collection dates exist among non-control samples with a successful demix: `{value, groups: [top 6 summarized groups, "Other" excluded], sites: [{site, dates: [ISO], n: [samples per date], series: [{group, values: [mean fraction per date]}]}], n_samples, n_dates}`; `site: null` = samples without a site. Needs sample metadata: the QC table, or `--samplesheet`                                                                                                                                                                                                                               |

`available` in a real `full` run with `--compare_barcodes container --freyja_boot`: `qc`, `amplicons`, `demix`,
`coverage_tracks`, `consensus`, `barcode_comparison`, `negative_controls`, `boot_intervals` and `time_series` are true;
`mutations_of_interest` is false (the option was not given).

### `samples[]`

All 23 `qc_summary.tsv` columns (PLAN section 6a) as keys, typed. `NA` -> `null`.

- integers: `raw_reads`, `trimmed_reads`, `mapped_reads`, `amplicons_total`, `amplicons_dropped`, `consensus_length`
- floats (percentages 0-100): `mapped_pct`, `mean_depth`, `depth_sd`, `median_depth`, `pct_genome_ge_1x`,
  `pct_genome_ge_10x`, `pct_genome_ge_min_depth`, `consensus_n_pct`
- strings: `sample`, `sample_type` (`sample`|`negative`, default `sample`), `collection_date`, `site`, `qc_status`
  (`PASS`|`WARN`|`FAIL`|`null`), `demix_status`, `demix_reason`, `barcode_version`
- `qc_reasons`: **list of strings** (the `;`-joined tokens split), `[]` when empty

`demix_status`, `demix_reason`, `consensus_length`, `consensus_n_pct`, `barcode_version` are filled by the builder
from the dedicated files (`demix_status.tsv`, `consensus_stats.tsv`); the dedicated file wins over a value in
`qc_summary.tsv`. `barcode_version` falls back to the run's `freyja_data.barcode_version` for samples that have a
demix status. `amplicons_total/dropped` are filled from `amplicon_coverage.tsv` when `NA`. The same filled-in values
are written to `summary/qc_summary.tsv` (`--write-qc-tsv`).

If `--qc` is absent (mode `aggregate`) the sample list is the union of all other per-sample inputs and QC fields are
`null`; `sample_type`, `collection_date` and `site` are then taken from `--samplesheet` when the sample is in it. A
sample that has no QC row **and** is not in the samplesheet gets `sample_type: "sample"` and no date or site (for
example an `aggregate` run that discovered a skipped negative control in `--previous_results`: the control is listed as
an ordinary sample flagged `demix_missing`). Samples present in `demix_status.tsv`, `freyja_aggregated.tsv` or the
controls aggregate but not in `qc_summary.tsv` are added; sample IDs in other files that are in no sample table are
ignored (with an `info` flag).

### `demix{sample}`

`{status, reason, barcode_version, groups, lineages, resid, coverage, intervals?}`

- `status`: `ok` | `failed` | `skipped_low_coverage` | `missing` (a demix input was given but this non-control
  sample has no result) | `null` (no demix information at all, e.g. a negative control that is in no demix input)
- `groups`: `[[name, fraction], ...]` Freyja `summarized`, descending, fractions 0-1 (5 dp)
- `lineages`: `[[name, fraction], ...]` Freyja `lineages`/`abundances`, descending
- `resid`, `coverage`: floats or null (`coverage` is Freyja's percent genome coverage; absent in older Freyja)
- `intervals`: only with `--freyja_boot`, and only for samples that have bootstrap rows: `{lineage: [lo, hi]}`, the 2.5
  and 97.5 percentiles as fractions 0-1 (5 dp), keyed by individual **lineage** name (not summarized group); it
  includes lineages with a zero point estimate. The key is absent (not `null`) otherwise
- Row labels in `freyja_aggregated.tsv` are matched to sample IDs; a leftover `.variants.tsv`/`.tsv` suffix is stripped
  only when the label matches no known sample. An `ok` status with no lineages is turned into `failed`, **except for
  negative controls**: they are deliberately not in the main aggregate (decision D8), so a control with status `ok`
  keeps `ok`, with empty `groups`/`lineages` here and its lineages under `negative_controls`. (An earlier version
  downgraded such a control to `failed`; that defect is fixed.)

### `negative_controls[]`

One entry per negative control (`sample_type = negative`) that demixed `ok` or appears in the controls aggregate:

| key                                | type               | notes                                                                                                                                                   |
| ---------------------------------- | ------------------ | ------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `sample`                           | string             | sample ID                                                                                                                                               |
| `status`, `reason`                 | string             | the control's demix status and reason (`ok`, empty)                                                                                                     |
| `qc_status`, `qc_reasons`          | string\|null, list | from the QC row (`FAIL` + `["ntc_reads"]` for a contaminated control)                                                                                   |
| `mapped_reads`, `pct_genome_ge_1x` | number\|null       | from the QC row                                                                                                                                         |
| `groups`, `lineages`               | list               | `[[name, fraction], ...]` from the controls aggregate (or, with `--include_negatives_in_aggregate`, from the main aggregate); `[]` if no file was given |
| `coverage`, `resid`                | number\|null       | Freyja's genome coverage % and fit residual                                                                                                             |
| `excluded_from_aggregate`          | bool               | `true` when the control is not in `freyja_aggregated.tsv` (the default)                                                                                 |

A control that was `skipped_low_coverage` (a clean control) has no entry; it is only in `samples[]` and `demix`. Real
entry (the contaminated control `NEG_CONTAM`, from the full run with `--compare_barcodes`/`--freyja_boot`):

```json
{
  "sample": "NEG_CONTAM",
  "status": "ok",
  "reason": "",
  "qc_status": "FAIL",
  "qc_reasons": ["ntc_reads"],
  "mapped_reads": 5998,
  "pct_genome_ge_1x": 99.37,
  "groups": [["Omicron", 1.0]],
  "lineages": [["BA.2", 1.0]],
  "coverage": 94.55,
  "resid": 0.0428,
  "excluded_from_aggregate": true
}
```

### Parsing `freyja_aggregated.tsv` (and `freyja_aggregated_controls.tsv`)

`summarized` is a Python list of tuples and is parsed with `ast.literal_eval` (regex fallback for numpy-2 reprs such as
`np.float64(0.3)`). `lineages` and `abundances` are either **numpy array reprs** (`['A' 'Q.3']`, no commas, newlines
inside a quoted field) or plain space-separated lists (`BA.2 AY.4`, as Freyja 2.x writes them), which
`ast.literal_eval` would silently mangle (adjacent string literals concatenate), so they are parsed by regex, never
evaluated. The controls file has the same format and goes through the same parser.

### `run_flags` codes

`freyja_data_missing`, `barcode_unknown`, `barcode_custom`, `barcode_container`, `mixed_barcode_versions`,
`previous_run_barcode_unknown`, `previous_run_barcode_differs`, `no_samples`, `qc_not_available`, `no_ntc`,
`ntc_reads`, `ntc_fail` (negative-control problems, level `error`), `qc_fail`, `qc_warn`, `demix_failed` (error),
`demix_skipped`, `demix_missing`, `all_samples_failed` (error), `barcode_compare_failed` (warn: the second demix of
`--compare_barcodes` failed for the listed samples, from `--compare-status`), `input_problems` (unparseable
rows/values, for example an unknown status in `demix_status.tsv`), `input_note`.

The control flags (`ntc_reads`, `ntc_fail`, `no_ntc`) are raised from the QC table. A negative control whose demix
status is `ok` raises no demix flag.

### `header` / `run_info.json` key spellings looked up

`git_commit` | `git_sha` | `commit_id` | `commit` | `git_revision` | `revision`; `nextflow_version` | `nextflow.version`;
`pipeline_version` | `pipeline.version`; `start` | `run_started`; barcode and Freyja fields from `freyja_data`.

### `--barcode-compare` input

Long TSV `sample, lineage, abundance_a, abundance_b` with optional constant columns `barcode_version_a`,
`barcode_version_b` used as labels (A defaults to `run_info.freyja_data.barcode_version`, B to `B`). The pipeline writes
exactly this, with all six columns, as `summary/barcode_compare.tsv` when `--compare_barcodes` is given (see
[output.md](output.md#summarybarcode_comparetsv-and-summarybarcode_compare_statustsv)); rows of samples that are in no
sample table are dropped.

## Report behaviour per mode

Sections whose data are absent show a "not available in mode X" note (QC, lineage, coverage tracks, amplicon
dropout, demix status, consensus) or are omitted entirely (negative controls, mutations of interest, time series,
barcode comparison). In mode `aggregate` the QC-derived sections are absent (`available.qc` is false), but the time
series is built when the samplesheet carries collection dates. Mode `report` re-renders an existing `report_data.json`
unchanged (or rebuilds it from the files in `summary/` with the same builder options, including `--controls-aggregated`,
`--samplesheet`, `--boot` and `--compare-status`).
