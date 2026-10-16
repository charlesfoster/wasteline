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
    [--dropout-impact dropout_impact.tsv] [--dropout-summary dropout_impact_summary.tsv]
    [--unexplained unexplained_mutations.tsv]
    [--timeseries freyja_timeseries.tsv] [--timeseries-lineages freyja_timeseries_lineages.tsv]
    [--growth-rates freyja_growth_rates.csv] [--growth-failed freyja_growth_rates.failed.txt]
    [--growth-inputs freyja_growth_inputs.json] [--dashboard-name freyja_dashboard.html]
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

Options for the optional analyses (each input is optional; absent means the section is simply not in the report):

| option                  | file                                                  | effect                                                                                                                                                                      |
| ----------------------- | ----------------------------------------------------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `--dropout-impact`      | `summary/dropout_impact.tsv`                          | one row per dropped amplicon (`DROPOUT_IMPACT`); feeds `dropout_impact.<sample>.amplicons` and the `dropout_affects_lineages` flag                                          |
| `--dropout-summary`     | `summary/dropout_impact_summary.tsv`                  | one row per sample; feeds `dropout_impact.<sample>.summary`. Either file alone is accepted                                                                                  |
| `--unexplained`         | `summary/unexplained_mutations.tsv` (the concat table) | feeds `unexplained_mutations` (`[]` for a header-only table). A per-sample table (13 columns, no `n_samples`/`recurrent`) loads too                                        |
| `--timeseries`          | `summary/freyja_timeseries.tsv`                       | summarized-group time series (preferred source of the time-series panel); feeds `timeseries.summarized`                                                                     |
| `--timeseries-lineages` | `summary/freyja_timeseries_lineages.tsv`              | full-lineage time series; feeds `timeseries.lineages`, capped at the 30 most abundant lineages                                                                             |
| `--growth-rates`        | `summary/freyja_growth_rates.csv`                     | Freyja's relgrowthrate table (`Lineage`, `Estimated Advantage`, `Bootstrap 95% interval`); feeds `growth_rates.rows`                                                        |
| `--growth-failed`       | `summary/freyja_growth_rates.failed.txt`              | why no growth rates were produced; its text becomes `growth_rates.reason` (only when there is no table)                                                                     |
| `--growth-inputs`       | `summary/freyja_growth_inputs.json`                   | what went into the fit (`n_samples`, `n_dates`, `window_span_days`, ...); feeds `growth_rates.inputs`                                                                       |
| `--dashboard-name`      | a string                                              | file name of the optional Freyja dashboard; stored as `dashboard.file` and shown as text only (never linked)                                                                |

## Top level

| key                     | type         | notes                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 |
| ----------------------- | ------------ | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `schema_version`        | int          | `1`                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   |
| `generated`             | string       | UTC timestamp `YYYY-MM-DDTHH:MM:SSZ`                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  |
| `mode`                  | string       | `full`, `demix`, `aggregate`, `report` (from `--mode`, else `run_info.mode` / `params.mode`, else `unknown`)                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          |
| `header`                | object       | derived summary used by the report header: `pipeline_name`, `pipeline_version`, `git_commit`, `nextflow_version`, `mode`, `freyja_version`, `barcode_version`, `barcode_version_iso`, `barcode_source`, `barcode_sha256`, `pathogen`, `run_started`, `run_info_generated` (all nullable). Looked up from `run_info.json` under several key spellings (see below)                                                                                                                                                                                                                                                                      |
| `provenance`            | object       | verbatim copy of `run_info.json` as it was when the report data were built (the mandatory `freyja_data` block, `freyja_data_compare` with `--compare_barcodes`, `params`, ...; the completion keys `containers`, `run_completed`, `duration` and `success` are added to `run_info.json` only after the run, so they are not here), plus `previous_runs` (list, default `[]`) of `{path, run_info_sha256, barcode_version}` (`barcode_version: "unknown"` for legacy Snakemake results)                                                                                                                                                |
| `params_summary`        | object       | selected scalar params from `provenance.params`, in a fixed order, only those present: `mode`, `primer_scheme`, `primer_bed`, `barcodes`, `freyja_eps`, `freyja_covcut`, `freyja_depthcutoff`, `freyja_varthresh`, `freyja_adapt`, `freyja_autoadapt`, `freyja_relaxedmrca`, `freyja_solver`, `freyja_confirmedonly`, `min_genome_coverage_for_demix`, `amplicon_min_depth`, `consensus_min_depth`, `consensus_freq`, `indel_freq`, `ntc_max_mapped_reads`, `ntc_max_genome_coverage_pct`, `mutations_of_interest`, `compare_barcodes`, `skip_consensus` (23 in a real run). `freyja_boot` and `annotate_variants` are not part of it |
| `thresholds`            | object       | `qc_min_primer_trimmed_pct`, `dropout_min_abundance`, `unexplained_min_freq`, `unexplained_min_excess`, `unexplained_min_depth`, `unexplained_min_alt_depth`, `unexplained_min_recurrence`, `unexplained_max_pvalue` are added **only when the run recorded them** in its params. `amplicon_min_depth` (builder default 10 when absent), `consensus_min_depth` (builder default 10 when the param is absent from `run_info.json`; the **pipeline default is 10**, which is what a real run records), `ntc_max_mapped_reads`, `ntc_max_genome_coverage_pct`, `min_genome_coverage_for_demix` (null when not in params)                                                                                                                                                                                                                                                                                                   |
| `available`             | object       | booleans: `qc`, `amplicons`, `demix`, `coverage_tracks`, `mutations_of_interest`, `consensus`, `barcode_comparison`, `negative_controls`, `boot_intervals`, `time_series`, `dropout_impact`, `unexplained_mutations`, `timeseries`, `growth_rates` (the last four are true when the corresponding input was given and is non-empty; `unexplained_mutations` is true also for a header-only table, `time_series` is true when either the derived or the table-based series has >= 2 dates)                                                                                                                                                                                                                                                                                                                                                                                                                                                             |
| `run_flags`             | list         | `{level: "error"\|"warn"\|"info", code, message, samples?}`, sorted error, warn, info (see codes below)                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               |
| `samples`               | list         | one object per sample, sorted by `sample` (see below)                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 |
| `demix`                 | object       | sample -> demix result (see below); one key per sample in `samples`                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   |
| `amplicons`             | object       | sample -> list of `{amplicon, start, end, mean_depth, dropped}` sorted by start; `start`/`end` 1-based inclusive. `{}` when absent                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    |
| `coverage_tracks`       | object\|null | `{contig, genome_length, bin_size, n_bins, samples: {sample: [mean depth per bin]}}`; positions 1-based, bin `i` covers `[i*bin_size+1, (i+1)*bin_size]`; `n_bins <= 600`; values are ints (>= 10) or 1-dp floats; positions missing from the depths file count as 0; only the first contig of a depths file is used                                                                                                                                                                                                                                                                                                                  |
| `mutations_of_interest` | list\|null   | `{sample, name, mutation, ref_depth, alt_depth, depth, alt_freq}` rows (`alt_freq` 0-1). `null` = not run; `[]` = run with no rows                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    |
| `consensus`             | list         | `{sample, length, n_count, n_pct, ambiguous_count, status, reason}`; `reason` is the 7th column of `consensus_stats.tsv` (text, `null` for `NA` or when an older 6-column table has no such column): why a consensus `failed`/is `empty`, or a warning for an `ok` one such as `skipped variant at NC_045512.2:22942 (overlaps another variant)`                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           |
| `barcode_comparison`    | object\|null | `{label_a, label_b, rows: [{sample, lineage, a, b, delta}]}` (fractions 0-1, 5 dp, `delta = b - a`, sorted by `abs(delta)` desc); `null` when there are no rows (no `--compare_barcodes`). `label_a`/`label_b` are the barcode versions of the two datasets (`custom` for a bare file)                                                                                                                                                                                                                                                                                                                                                |
| `negative_controls`     | list         | one entry per negative control that demixed `ok` or has a row in the controls aggregate (see below); `[]` otherwise                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   |
| `dropout_impact`        | object\|null | sample -> `{summary, amplicons}` (see below); `null` when neither dropout file was given                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               |
| `unexplained_mutations` | list\|null   | rows of the unexplained-mutations table (see below); `null` = not run, `[]` = run with no rows                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         |
| `timeseries`            | object\|null | the freyja_timeseries tables (see below); preferred over `time_series` by the report                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 |
| `growth_rates`          | object\|null | Freyja relgrowthrate result or failure note (see below)                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               |
| `dashboard`             | object\|null | `{file}` when `--dashboard-name` was given                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              |
| `time_series`           | object\|null | present only when >= 2 distinct collection dates exist among non-control samples with a successful demix: `{value, groups: [top 6 summarized groups, "Other" excluded], sites: [{site, dates: [ISO], n: [samples per date], series: [{group, values: [mean fraction per date]}]}], n_samples, n_dates}`; `site: null` = samples without a site. Needs sample metadata: the QC table, or `--samplesheet`                                                                                                                                                                                                                               |

`available` in a real `full` run with `--compare_barcodes container --freyja_boot`: `qc`, `amplicons`, `demix`,
`coverage_tracks`, `consensus`, `barcode_comparison`, `negative_controls`, `boot_intervals` and `time_series` are true;
`mutations_of_interest` is false (the option was not given).

### `samples[]`

All 27 `qc_summary.tsv` columns (PLAN section 6a plus the four appended primer columns) as keys, typed. `NA` -> `null`.
A qc table written before the primer columns existed loads fine: the four keys are `null`.

- integers: `raw_reads`, `trimmed_reads`, `mapped_reads`, `amplicons_total`, `amplicons_dropped`, `consensus_length`
- floats (percentages 0-100 unless noted): `primer_trimmed_pct`, `primer_scheme_enrichment` (a ratio, not a percentage), `mapped_pct`, `mean_depth`, `depth_sd`, `median_depth`, `pct_genome_ge_1x`,
  `pct_genome_ge_10x`, `pct_genome_ge_min_depth`, `consensus_n_pct`
- strings: `sample`, `sample_type` (`sample`|`negative`, default `sample`), `collection_date`, `site`, `qc_status`
  (`PASS`|`WARN`|`FAIL`|`null`), `demix_status`, `demix_reason`, `barcode_version`, `primer_scheme_check`
  (`ok`|`mismatch`|`uninformative`|`null`), `primer_scheme_best` (scheme name or `null`)
- `qc_reasons`: **list of strings** (the `;`-joined tokens split), `[]` when empty. Tokens include `low_coverage`,
  `amplicon_dropout`, `no_reads`, `no_depths`, `ntc_reads`, `primer_scheme_mismatch` (WARN) and `low_primer_trimming`
  (only when the user set `--qc_min_primer_trimmed_pct`)
- `consensus_status`, `consensus_reason`: copied from `consensus_stats.tsv` for samples that have a consensus row
  (`consensus_reason` is `null` when there is none)
- only present when the corresponding input was given: `dropout_pct_informative_lost` (float or `null`),
  `dropout_affects_lineages` (bool, `null` when the sample has no lineage estimate) with `--dropout-impact`/`--dropout-summary`;
  `n_unexplained` (int, `null` when the sample has no successful demix) with `--unexplained`

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

### `dropout_impact{sample}`

From `dropout_impact.tsv` / `dropout_impact_summary.tsv` (definitions in `bin/dropout_impact.py`). Every sample of the
summary appears; samples in neither file are not listed; sample IDs that are in no sample table are ignored (with an
`info` flag).

```json
{
  "summary": {
    "n_dropped_amplicons": 2,
    "informative_sites_total": 8006,
    "informative_sites_in_dropped": 850,
    "pct_informative_sites_lost": 10.62,
    "n_dropped_affecting_detected_lineages": 1
  },
  "amplicons": [
    {
      "amplicon": "nCoV-2019_22",
      "start": 21533,
      "end": 22612,
      "mean_depth": 3.1,
      "n_informative_sites": 532,
      "n_sites_in_detected_lineages": 10,
      "affects_detected_lineage": true,
      "affected_lineages": [{ "lineage": "XBB.1.5", "n_sites": 9 }]
    }
  ]
}
```

`summary` is `null` when only `--dropout-impact` was given. `NA` -> `null`: `n_sites_in_detected_lineages` and
`n_dropped_affecting_detected_lineages` are `null` for a sample without a lineage estimate (skipped, failed, control).
`amplicons` is sorted by start, `[]` when nothing was dropped. `affected_lineages` is `[]` for `NA`.

### `unexplained_mutations[]`

`{sample, position, ref, alt, observed_freq, expected_freq, excess, depth, alt_depth, gene, aa_change, effect,
n_samples, recurrent, p_value}` in file order (recurrent first, as the concat table is sorted). Frequencies are
fractions 0-1 (6 dp), `p_value` a float (`2.39e-42`, not rounded), `recurrent` a boolean; `gene`/`aa_change` are
`null` outside genes. From a per-sample table `n_samples` is `null` and `recurrent` `false`.

### `timeseries`

From the long freyja_timeseries tables (`bin/timeseries.py`): `{enabled, has_viral_load, min_coverage, sites, dates,
summarized, lineages, lineages_total_groups}`. `summarized`/`lineages` are lists of `{site, collection_date,
lineage_group, abundance, abundance_load_weighted, n_samples, mean_viral_load, total_viral_load}` (`NA` -> `null`,
`site: null` for samples without a site); `has_viral_load` is true when any `abundance_load_weighted` is non-null;
`min_coverage` comes from the `freyja_plot_mincov` param; `lineages` keeps the 30 most abundant lineages and
`lineages_total_groups` is then the number of lineages in the table (else `null`). The report draws the top 6 groups.

### `growth_rates`

`{status, reason, inputs, params, rows}`: `status` is `ok` (a table was given) or `failed` (only a failure note);
`reason` is the text of `freyja_growth_rates.failed.txt` (1000 characters at most) when failed; `inputs` the
`freyja_growth_inputs.json` content (or `null`); `params` the `freyja_growth_*` params recorded by the run (`days`,
`nboots`, `seed`, `serial_interval`, `scale_by_viral_load`); `rows` is `[{lineage, advantage_pct, interval_text}]`
(`advantage_pct` parsed from `60.3%`; the interval is kept as Freyja's text, e.g. `[44.74% , 49.89%]`, because it is on
a different scale than the estimate).

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

`consensus_failed` (warn: the consensus status is `failed`, so the sequence is all-N; `samples` lists them and the
consensus table has the reason) and `consensus_warnings` (info: status `ok` with a non-NA reason, e.g. a skipped
overlapping variant; `samples` lists them).

Added for the primer scheme check and the optional analyses: `primer_scheme_mismatch` (warn: samples whose
`primer_scheme_check` is `mismatch`; the message names the best-matching scheme(s) and the declared
`params.primer_scheme`), `dropout_affects_lineages` (warn: any dropped amplicon hides lineage-informative sites of a
detected lineage; `samples` lists them) and `unexplained_mutations` (info: at least one recurrent unexplained
mutation; `samples` lists the samples carrying a recurrent row).

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
barcode comparison, dropout impact, unexplained mutations). An unexplained-mutations table that exists but has no rows
shows "None found with the current thresholds". The primer-trimming/primer-scheme columns of the QC table appear only
when at least one sample has a value. A report built without any of the optional analysis inputs is identical to one
built before they existed. In mode `aggregate` the QC-derived sections are absent (`available.qc` is false), but the time
series is built when the samplesheet carries collection dates. Mode `report` re-renders an existing `report_data.json`
unchanged (or rebuilds it from the files in `summary/` with the same builder options, including `--controls-aggregated`,
`--samplesheet`, `--boot` and `--compare-status`).
