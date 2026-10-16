# Wiring notes: time-series surveillance (TIMESERIES_TABLE, FREYJA_RELGROWTHRATE, FREYJA_DASH)

New files (nothing existing was edited): `bin/timeseries.py`, `tests/make_timeseries_fixture.py`, `tests/data_timeseries/`,
`tests/python/test_timeseries.py`, `modules/local/{timeseries_table,freyja_relgrowthrate,freyja_dash}/`.
All three modules have a `stub:`, a `topic: versions` tuple, `environment.yml` and an nf-test (stub + real Docker, all green).
The Freyja ones use the pinned `freyja:2.0.6--pyhdfd78af_0` (so `bin/bump_freyja.py` picks them up); `TIMESERIES_TABLE` uses
`python:3.12.12` like `QC_CONCAT`.

## What the scripts do

`bin/timeseries.py --aggregated freyja_aggregated.tsv --metadata sheet.csv [--level summarized|lineages] [--min-abundance 0.0]
[--min-coverage 60] -o freyja_timeseries.tsv` (stdlib only, tolerant parser for Python/numpy reprs).

`freyja_timeseries.tsv` (tab-separated, header always written, `NA` for missing):
`site, collection_date, lineage_group, abundance, abundance_load_weighted, n_samples, mean_viral_load, total_viral_load`.

- One row per site, date and lineage group. `abundance` = mean over the samples of that site and date (a group a sample lacks counts as 0).
  `abundance_load_weighted` = sum(a_i * load_i) / sum(load_i) over the samples that have a load (`NA` if none). `n_samples` counts every
  contributing sample; `mean_viral_load`/`total_viral_load` use those with a load. Rows sorted by site, date, abundance (desc), group.
- Excluded, counted on stderr (`timeseries.py: 33 aggregate samples, 30 used; excluded: no_metadata=0, negative_control=1, no_date=1,
  bad_date=0, low_coverage=1, unparsable_rows=0`): sample not in the metadata, `sample_type` negative (case-insensitive), blank or
  non-ISO `collection_date`, aggregate `coverage` < `--min-coverage` (percent; our gate is `>=`, Freyja's own is strictly `>`), rows whose
  groups cannot be parsed. Blank/non-numeric/zero/negative `viral_load` = unknown (counted as `invalid_viral_load` when non-blank).
  Empty `site` -> `NA`. Header-only output is valid (nothing usable); a missing `sample` column is exit 1.
- Metadata columns are matched case-insensitively and every one except `sample` may be absent (a BOM is tolerated). A duplicate sample
  keeps the first row.
- Helper mode used by the two Freyja modules: `--freyja-metadata-out m.csv --freyja-aggregated-out a.tsv [--freyja-summary-out s.json]
  [--require-viral-load] [--window-days N]` writes Freyja's metadata CSV (`Sample,sample_collection_datetime[,viral_load]`; the load column
  only when some kept sample has one), the aggregate cut to exactly those samples in the original order (Freyja raises a KeyError on a
  sample without metadata), and a JSON summary (`n_samples, n_dates, n_dates_window, window_span_days, n_with_viral_load, first_date, last_date,
  span_days, viral_load_column`).

## Includes

```nextflow
include { TIMESERIES_TABLE                              } from '../modules/local/timeseries_table/main'
include { TIMESERIES_TABLE as TIMESERIES_TABLE_LINEAGES } from '../modules/local/timeseries_table/main'  // optional second level
include { FREYJA_RELGROWTHRATE                          } from '../modules/local/freyja_relgrowthrate/main'
include { FREYJA_DASH                                   } from '../modules/local/freyja_dash/main'
```

## Signatures

| process | inputs | outputs (emit) |
|---|---|---|
| `TIMESERIES_TABLE` | `path aggregated` (`freyja_aggregated.tsv`); `path metadata` (samplesheet-like CSV, staged as `metadata.csv`) | `table` = `freyja_timeseries.tsv` |
| `FREYJA_RELGROWTHRATE` | `path aggregated`; `path metadata` (staged as `metadata.csv`); `path dataset` (Freyja dataset dir, same channel as `FREYJA_AGGREGATE_PLOT`) | `growth_rates` = `freyja_growth_rates.csv` (optional); `failed` = `freyja_growth_rates.failed.txt` (optional); `inputs` = `freyja_growth_inputs.json` |
| `FREYJA_DASH` | `path aggregated`; `path metadata`; `path dataset`; `val title`; `val intro` (both inserted as HTML, quotes are handled) | `dashboard` = `freyja_dashboard.html` (optional); `growth_rates` = `freyja_dashboard_rel_growth_rates.csv` (optional); `failed` = `freyja_dashboard.failed.txt` (optional); `inputs` = `freyja_dashboard_inputs.json` |

Exactly one of `growth_rates`/`failed` (resp. `dashboard`/`failed`) is produced; a soft failure never fails the run and the failure note
always says why. Run-level inputs must be single-item (value) channels: use `.first()` / `channel.value` for `ch_samplesheet`,
`FREYJA_AGGREGATE_PLOT.out.aggregated` is already one item and `FREYJA_DATA.out.dataset` is a value.

`metadata` is the samplesheet in effect (`ch_samplesheet`, copied to `summary/samplesheet.csv`): the scripts read `sample`, `sample_type`,
`collection_date`, `site`, `viral_load` and ignore every other column (fastq paths etc.).

`ext` keys handled by the Freyja modules (everything else in `ext.args` goes to `freyja relgrowthrate` / `freyja dash` verbatim):
`ext.days` (56), `ext.mincov` (60), `ext.seed` (1), `ext.min_abundance` (0.01), `ext.min_dates` (4), `ext.min_span_days` (14).
Do not put `--days`/`--mincov` in `ext.args` (the data-sufficiency check uses `ext.days`/`ext.mincov`). `--scale_by_viral_load` in `ext.args`
automatically restricts the fit to samples that have a load (`--require-viral-load`).

## Workflow placement (after `FREYJA_AGGREGATE_PLOT`, `full`/`demix`/`aggregate` modes)

```nextflow
// two or more distinct collection dates among the aggregated non-negative samples
def ch_ndates = ch_aggregate_in
    .filter { meta, _demix -> meta.collection_date && meta.sample_type != 'negative' }
    .map { meta, _demix -> meta.collection_date }
    .unique()
    .count()
def ch_ts_go = ch_ndates.filter { n -> n >= 2 }

if (params.timeseries) {
    def ch_ts_agg = FREYJA_AGGREGATE_PLOT.out.aggregated.combine(ch_ts_go).map { agg, _n -> agg }
    TIMESERIES_TABLE(ch_ts_agg, ch_samplesheet.first())
    TIMESERIES_TABLE_LINEAGES(ch_ts_agg, ch_samplesheet.first())   // report toggle "All lineages"
}
if (params.freyja_growth) {
    FREYJA_RELGROWTHRATE(ch_ts_agg, ch_samplesheet.first(), FREYJA_DATA.out.dataset)
}
if (params.freyja_dash) {
    FREYJA_DASH(ch_ts_agg, ch_samplesheet.first(), FREYJA_DATA.out.dataset, params.freyja_dash_title, params.freyja_dash_intro)
}
```

(`ch_ts_agg` has to be defined once outside the `if`s.) With fewer than two dates the processes do not start, so no empty-looking table
appears; the report shows no panel (existing rule: needs at least two distinct dates). `FREYJA_RELGROWTHRATE` and `FREYJA_DASH` apply their
own, stricter requirement (see below) and write a visible `*.failed.txt` when it is not met.

## conf/modules.config

```groovy
withName: 'TIMESERIES_TABLE' {
    ext.args   = { "--level summarized --min-abundance ${params.timeseries_min_abundance} --min-coverage ${params.freyja_plot_mincov}" }
    publishDir = [ path: { "${params.outdir}/summary" }, mode: params.publish_dir_mode ]          // summary/freyja_timeseries.tsv
}
withName: 'TIMESERIES_TABLE_LINEAGES' {
    ext.args   = { "--level lineages --min-abundance ${params.timeseries_min_abundance} --min-coverage ${params.freyja_plot_mincov}" }
    publishDir = [ path: { "${params.outdir}/summary" }, mode: params.publish_dir_mode,
                   saveAs: { fn -> fn == 'freyja_timeseries.tsv' ? 'freyja_timeseries_lineages.tsv' : null } ]
}
withName: 'FREYJA_RELGROWTHRATE' {
    ext.args           = { "--nboots ${params.freyja_growth_nboots} --serial_interval ${params.freyja_growth_serial_interval} --grthresh ${params.freyja_growth_grthresh} --pathogen ${params.pathogen}" + (params.freyja_growth_scale_by_viral_load ? ' --scale_by_viral_load' : '') }
    ext.days           = { params.freyja_growth_days }
    ext.mincov         = { params.freyja_plot_mincov }
    ext.seed           = { params.freyja_growth_seed }
    ext.min_abundance  = { params.freyja_growth_min_abundance }
    ext.min_dates      = { params.freyja_growth_min_dates }
    ext.min_span_days  = { params.freyja_growth_min_span_days }
    publishDir = [ path: { "${params.outdir}/summary" }, mode: params.publish_dir_mode,
                   pattern: 'freyja_growth_{rates.csv,rates.failed.txt,inputs.json}' ]
}
withName: 'FREYJA_DASH' {
    ext.args           = { "--nboots ${params.freyja_growth_nboots} --serial_interval ${params.freyja_growth_serial_interval} --grthresh ${params.freyja_growth_grthresh} --pathogen ${params.pathogen}" + (params.freyja_growth_scale_by_viral_load ? ' --scale_by_viral_load' : '') }
    ext.days           = { params.freyja_growth_days }
    ext.mincov         = { params.freyja_plot_mincov }
    ext.seed           = { params.freyja_growth_seed }
    ext.min_abundance  = { params.freyja_growth_min_abundance }
    ext.min_dates      = { params.freyja_growth_min_dates }
    ext.min_span_days  = { params.freyja_growth_min_span_days }
    publishDir = [ path: { "${params.outdir}/summary" }, mode: params.publish_dir_mode,
                   pattern: 'freyja_dashboard{.html,.failed.txt}' ]    // not the inputs json or the dashboard's own growth csv
}
```

Notes: `--pathogen` and `--lineageyml` conflict only in that `--lineageyml` (set by the modules when the dataset has `lineages.yml`) wins. The
`withName` selector `TIMESERIES_TABLE` also matches `TIMESERIES_TABLE_LINEAGES` only if written as a regex, so keep them as two exact blocks.
`conf/containers.lock.config` needs three new `withName` entries (regenerate with `bin/lock_containers.py`, needs the registry): `FREYJA_DASH`,
`FREYJA_RELGROWTHRATE` (Freyja digest) and `TIMESERIES_TABLE`, `TIMESERIES_TABLE_LINEAGES` (python digest).

## Params to add (nextflow.config + nextflow_schema.json)

| param | default | meaning |
|---|---|---|
| `timeseries` | `true` | write `summary/freyja_timeseries.tsv` (+ `_lineages.tsv`) when >= 2 distinct collection dates exist |
| `timeseries_min_abundance` | `0.0` | drop site-date groups below this mean abundance (0-1) |
| `freyja_growth` | `false` | run `freyja relgrowthrate` -> `summary/freyja_growth_rates.csv` |
| `freyja_growth_days` | `56` | Freyja `--days`: window (back from the latest date) used for the fit |
| `freyja_growth_nboots` | `1000` | bootstrap replicates (the module fixes the RNG seed) |
| `freyja_growth_seed` | `1` | numpy seed for the bootstrap (Freyja has no seed option; the module seeds numpy's global RNG before `cli()`) |
| `freyja_growth_serial_interval` | `5.5` | days; the table is `exp(serial_interval * r) - 1` |
| `freyja_growth_grthresh` | `0.001` | Freyja `--grthresh`: only lineages with a mean abundance above this over the window get a rate |
| `freyja_growth_min_abundance` | `0.01` | per-sample abundance below which a lineage is lumped into `Other` (see the `--thresh` quirk below) |
| `freyja_growth_min_dates` | `4` | fewer distinct dates in the window -> `freyja_growth_rates.failed.txt` |
| `freyja_growth_min_span_days` | `14` | dates must span at least this many days |
| `freyja_growth_scale_by_viral_load` | `false` | Freyja `--scale_by_viral_load` (needs `viral_load`; samples without one are dropped from the fit) |
| `freyja_dash` | `false` | run `freyja dash` -> `summary/freyja_dashboard.html` (opt-in extra) |
| `freyja_dash_title` | `'wasteline wastewater dashboard'` | HTML title |
| `freyja_dash_intro` | `'Lineage prevalence over time from wasteline.'` | HTML introduction |

The existing `freyja_plot_mincov` (60) is reused for `--min-coverage`/`--mincov`.

## The `viral_load` column (schema and helper changes still to do)

Optional numeric column on every samplesheet flavour (empty = unknown; positive number, e.g. gene copies per litre):

- `assets/schema_input.json`, `schema_input_demix.json`, `schema_input_aggregate.json`: add
  `"viral_load": {"type": "number", "exclusiveMinimum": 0, "errorMessage": "viral_load must be a positive number (empty = unknown)", "meta": ["viral_load"]}`
  (not in `required`; nf-schema accepts an empty cell as absent for `number` only if the property has no `default`; if empty cells are rejected,
  use `"type": ["number", "string"]` with `"pattern"` or `"anyOf"` of a number and an empty-string enum, as done for `collection_date`).
- `bin/make_samplesheet.py`: add `viral_load` to `META_COLUMNS` / the `--metadata` CSV and the written columns (only when any sample has one),
  validate it is a positive number, carry it over from `summary/samplesheet.csv` like `site` (`load_results_sheet`).
- `bin/check_samplesheet.py`: validate (blank or positive finite number).
- `summary/samplesheet.csv` is written from `ch_samplesheet` as-is, so the column travels to `--previous_results` runs for free; the skipped-sample
  sheet (`skipped_samplesheet.csv`, seed `sample,sample_type,collection_date,site`) may stay without it.
- The pipeline's `meta` map does not need `viral_load`: the modules read the CSV directly.

## Report integration (proposal)

`build_report_data.py` gets optional inputs `freyja_timeseries.tsv`, `freyja_timeseries_lineages.tsv`, `freyja_growth_rates.csv`,
`freyja_growth_rates.failed.txt`, `freyja_growth_inputs.json`; `report_data.json` gets:

```json
"timeseries": {
  "enabled": true,
  "has_viral_load": true,
  "min_coverage": 60,
  "sites": ["SiteA", "SiteB", "SiteC"],
  "dates": ["2026-01-05", "..."],
  "summarized": [{"site": "SiteA", "collection_date": "2026-01-05", "lineage_group": "Omicron", "abundance": 0.8006,
                  "abundance_load_weighted": 0.8006, "n_samples": 1, "mean_viral_load": 175675.0, "total_viral_load": 175675.0}],
  "lineages": [ "...same row shape..." ]
},
"growth_rates": {
  "status": "ok",                      // ok | failed | disabled
  "reason": null,                      // text of the .failed.txt first line when failed
  "inputs": {"n_samples": 30, "n_dates": 8, "n_dates_window": 8, "window_span_days": 49, "n_with_viral_load": 26},
  "params": {"days": 56, "nboots": 1000, "serial_interval": 5.5, "seed": 1, "scale_by_viral_load": false},
  "rows": [{"lineage": "XBB.1.5", "advantage_pct": 60.3, "interval_text": "[44.74% , 49.89%]"}]
}
```

Rows are the TSV rows (numbers parsed, `NA` -> `null`); `growth_rates.rows` come from parsing `Estimated Advantage` (`60.3%`); keep the
Freyja interval as text (see the unit warning below). Panel: site selector (all sites + each site; "all sites" = n_samples-weighted mean of the
site-date rows), a toggle "Summarized groups / All lineages" (existing), a toggle "Unweighted / Load-weighted" shown only when `has_viral_load`
and at least one row has a non-null `abundance_load_weighted`, hover text with `n_samples` and `mean_viral_load`, and the growth-rate table under the
chart (or the failure reason in a notice). The existing Mean-abundance time-series panel can switch to this data when present.

## Findings with the real Freyja 2.0.6 container

1. **Same pandas 3 bug as `freyja plot` (PLAN D6).** Plain `freyja relgrowthrate` and `freyja dash` both die with
   `TypeError: Invalid value '[list([...]) ...]' for dtype 'str'` and write nothing. Both modules therefore run Freyja through a small
   `python3` wrapper that sets `pd.options.future.infer_string = False` before `freyja._cli.cli()`: attempt 1 is the wrapper without that
   switch (so the module self-heals when upstream is fixed), attempt 2 with it. Both succeeded on attempt 2 here. If neither yields output, a
   `*.failed.txt` with the last log lines is written and the run continues; exit 0 without output counts as failure.
2. **Seed.** No `--seed`. The bootstrap uses `np.random.randint` (global RNG), so the wrapper calls `np.random.seed(ext.seed)`; two independent
   runs gave byte-identical tables (`XBB.1.5,60.3%,"[44.74% , 49.89%]"`).
3. **Minimum data.** Freyja never refuses: with one date every lineage gets `+73.3%` (the optimiser start value `exp(5.5 * 0.1) - 1`), with two dates the
   numbers are wild (`BA.2 -100%`, interval `[-5499.99% , 399.95%]`), with three they are erratic. The modules therefore require, within the
   last `days` days, at least `min_dates` (4) distinct dates spanning at least `min_span_days` (14) days, else write
   `... insufficient data: 3 distinct collection date(s) spanning 14 day(s) within the last 56 days; need at least 4 dates spanning at least 14 days`.
   Freyja's own fit window is the last `--days` days plus one earlier point; it fits one logistic curve per lineage to the mean abundance of every date,
   pooling all sites. More dates and a wider span give tighter results; 8 weekly points are a comfortable minimum for a stable estimate.
4. **`--thresh` scales with the number of samples (upstream quirk).** Freyja merges a sample's lineage into `Other` when its abundance is below
   `thresh * n_samples`, so the default `--thresh 0.01` with 30 samples means a 30% cut-off: our synthetic XBB.1.5 and JN.1 vanished into `Other`
   and the fit gave XBB.1.5 +102.6% (truth about +73% for XBB.1.5 over BA.2; about +60% in the one-vs-rest model Freyja fits). The modules divide the per-sample cut-off
   (`min_abundance`, 0.01) by the number of kept samples unless `--thresh` is in `ext.args`.
5. **Output semantics.** The table is per lineage against all others (logistic `1 / (1 + b * exp(-r t))` fitted to the lineage's share); `Estimated Advantage` =
   `exp(serial_interval * r) - 1`, the bootstrap interval is `serial_interval * r` percentiles (a *different scale*: the estimate can lie outside its own
   interval, e.g. `60.3%` vs `[44.74% , 49.89%]`). Treat the interval as a rough indication of uncertainty, not as a CI of the estimate.
6. **Synthetic check** (`tests/data_timeseries`, 3 sites x 8 weekly dates, 30 usable samples; truth: XBB.1.5 over BA.2 log-odds slope 0.10/day = +73% per 5.5 days):
   Freyja gave XBB.1.5 +60.3%, BA.2 -35.1%, JN.1 (stable) -1.0%, BA.1 (noise lineage) -4.8%; with `--scale_by_viral_load` (26 samples with a load) +59.5% / -35.1% / -0.3%.
   Sign and order of magnitude are right; the magnitude is lower because of the one-vs-rest fit with a stable third lineage.
7. **Dashboard.** `freyja dash` needs `agg`, `metadata`, and *files* for title and intro (the module writes them), writes `<output minus .html>_rel_growth_rates.csv`
   next to the HTML, and also fails on pandas 3. The HTML is 4.8 MB: Plotly is embedded, but the page loads Bootstrap CSS and JS from
   `cdn.jsdelivr.net` and a logo from `andersen-lab.com`, so **without internet the charts work but the layout is unstyled and the logo is missing**. It is an opt-in extra;
   the pipeline's own `report.html` stays offline. Title/intro are inserted as raw HTML. Dates with several samples at several sites are pooled by Freyja (it has no site notion).
8. Freyja writes `agg_df.csv`, `df_ab_lin_raw.csv` into its cwd; the modules run it in `fg/` so they never reach the outputs. Freyja's progress bar goes to stderr; the failure notes
   filter it out.
9. Nextflow detail: an `exit` inside a module script skips the trailing `eval()` (versions) capture and fails the task with `.command.env` missing, so both Freyja scripts
   run their body in a shell function and `return` early for soft failures.

## Test data and tests

- `python tests/make_timeseries_fixture.py --out tests/data_timeseries --seed 42 [--aggregate]` (deterministic; `--aggregate` runs the real `freyja aggregate` in Docker, exactly like
  `FREYJA_AGGREGATE_PLOT`). 33 samples (30 usable), 152 KB: `samplesheet.csv` (`sample,demix,sample_type,collection_date,site,viral_load`), `demix/*.demix`, `freyja_aggregated.tsv`, `truth.json`.
  Edge cases: `NEG_1` (negative), `NODATE_1` (blank date), `LOWCOV_1` (coverage 42.5), 4+ blank loads, site-dates with two samples.
- `pixi run --frozen pytest tests/python/test_timeseries.py` (36 tests, includes byte-identical regeneration against the committed fixture and hand-computed weighting).
- `NFT_WORKDIR=<scratch> pixi run --frozen nf-test test modules/local/{timeseries_table,freyja_relgrowthrate,freyja_dash}/tests/main.nf.test --profile docker`
  (stub 3/3, real 1 + 3 + 2 tests; about 3 min in total). `bin/*` is an nf-test trigger, so adding `bin/timeseries.py` triggers a full run in CI.
- The relgrowthrate/dash tests use `tests/nextflow.config` (`--nboots 200` / `50`) and `tests/scale.config` inside the module directories.
