# S2 U5a: the modes demix / aggregate / report / update_db on panel results

Branch `feature/multi-target-in-pipeline`. Follows `multi_target_u4b.md`. The temporary "panel results cannot be re-demixed/aggregated/reported
yet" error is gone. Results written in the old layouts are read exactly as before (pipeline snapshots unchanged, a baseline-vs-new stub
comparison of demix, aggregate and report on legacy results differs only in paths and run names).

## Which targets are in play (`targetsForMode(mode)`, utils)

`targetsForModeAll(mode)`, then `--targets` selects (`targetsForMode`; names that are not among the targets are an error):

1. `--panel` given: its rows (registry + overrides, for example a `barcodes` column). For demix/aggregate/report only the rows whose target is in
   the input (an absent one: one warning at start-up; none left: error). The input must then be panel-layout (`checkPanelResultsInput`):
   `--previous_results` dirs with `summary/targets.tsv`, a `--input` sheet with a `target` column, `--results_dir` panel results; `--demix_dir` is refused.
2. no `--panel`, mode demix/aggregate/report, the input is a panel run's: the target names found (`resultsTargetNames(dir)`: `analysed=true` rows of
   `summary/targets.tsv`, else `summary/<t>/` folders beside `samples/<id>/<t>/freyja/`; or the `--input` sheet's `target` column), each resolved
   through the registry (`resultTarget`; a name that is not in `assets/panels/registry.csv` is an error that asks for `--panel <csv>`). Panel and
   single-target results in one `--previous_results` are refused. `--pathogen`/`--fasta` cannot be told apart from the SARS defaults here, so they are not checked.
3. otherwise `[ targetFromParams() ]` (every legacy run).

The target's layout flag follows the input: panel results give `target.panel` true (`sample_sub`/`summary_sub` `/<name>`) also for one selected target.
`targetsFromParams()` keeps its meaning (mode full). `update_db` uses (1) or (3); `--targets` filters it; a panel row's own `barcodes` is dropped
there (it would otherwise replace `latest`). `finaliseRunInfo` uses `targetsForMode(params.mode)`.

## Discovery (`bin/make_samplesheet.py --from-results`, `DISCOVER_RESULTS`)

Panel layout = `summary/targets.tsv`, or a `summary/<t>/` folder beside `samples/<id>/<t>/freyja/` (`is_panel_layout`). Rows are (sample, target) for
every `samples/<id>/<t>/freyja/` with the required files; the sheet gets a `target` column right after `sample`; fastp JSON comes from
`samples/<id>/fastp/`; sample_type/date/site/viral_load/protocol/primer_scheme from `summary/<t>/samplesheet.csv` (then `summary/samplesheet.csv`).
The skipped TSV is `sample, target, reason, sample_type, collection_date, site` (a target whose `.demix` is missing in aggregate mode, an unusable one in demix mode).
New `--targets a,b` (DISCOVER_RESULTS gets it as a third input; `''` = all). Every other layout: no target column, output byte-identical to before.
`check_samplesheet.py` accepts the `target` column in demix/aggregate sheets and checks duplicates on (sample, target).
`RESULTS_INPUT` gives each row the meta of its target; the duplicate check uses `meta.key`; `ch_skipped` items got a sixth element, the target (`''` for the other layouts).

## Workflow

- demix/aggregate: rows carry `meta.target`; everything of U3/U4 is reused. `carriedLowConfidence(targets)` reads `summary/<t>/demix_status.tsv` for panel
  results and is keyed like `meta.key` (`ID` or `ID__target`). A skipped sample goes to the status file of its own target (legacy: every target of the run = the one).
  The summary samplesheet given to the site/lineage summaries has one row per sample (a sheet row beats a skipped one); the report's skipped sheet is made unique.
- `--input` demix/aggregate: `sheetRowTarget` matches the row's `target` to the targets of the run (rows of targets left out by `--targets`/`--panel` are dropped).
- report on panel results: per target, `reportBuildInputs(results, t)` reads `summary/<t>/*` and `samples/*/<t>/freyja/*.depths`, `pipeline_info/<t>/run_info.json`
  (required at start-up); a target with an existing `report_data.json` is only rendered unless `--report_anonymize` is set. One `BUILD_REPORT_DATA` and one `RENDER_REPORT`
  call (a process can only be called once; the two kinds of target are mixed into one channel). Output: `summary/<t>/` of `--outdir`. The legacy path is untouched.
- `previous_runs` of a panel target's `run_info.json` points at `pipeline_info/<t>/run_info.json` of the results (sha256, barcode version); `panel` is `from results` without `--panel`.

## Leftovers

- `--demix_dir` cannot be used with panel results (no target information), and a `--previous_results` mix of panel and single-target results is refused.
- Two result directories with different optional columns (viral_load, protocol, primer_scheme) are concatenated by `collectFile(keepHeader)` in `RESULTS_INPUT`
  under the first header (existing behaviour; harmless for panel results of one run).
- `reanalyse_runs.py`, `runs_store.py`, `replay_run.py` and the overview are U5b/U6.
