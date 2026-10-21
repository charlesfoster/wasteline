# Low-confidence estimates (PLAN D18): wiring notes

Samples whose genome coverage (% of positions at `freyja_covcut` or deeper) lies between `demix_floor_coverage` (default 25)
and `min_genome_coverage_for_demix` (default 50) are demixed and flagged `low_confidence`, instead of `skipped_low_coverage`.
Below the floor nothing changes. Real data: 10 of 28 samples sit at 25-50 %, half of them with mean depth >= 200x (amplicon
dropout, not read shortage), so dropping them threw away usable signal.

Status vocabulary after this change: `ok | low_confidence | failed | skipped_low_coverage | missing`. A `.demix` file exists
for `ok` and `low_confidence`.

## Rules (decisions already made)

- `low_confidence` samples ARE in `freyja_aggregated.tsv`, `lineage_summary.csv/.xlsx` (status shown, values filled, reason
  filled) and the report (badge + one INFO flag `demix_low_confidence` listing them). Freyja's plot mincov (60) already hides them
  from `freyja_plot.pdf`: no change there.
- They are EXCLUDED from the time series tables, growth rates, the dashboard and site aggregation unless
  `--include_low_confidence true` (default false).
- The band is off when `demix_floor_coverage` is 0 or >= `min_genome_coverage_for_demix`: behaviour is then identical to today.
- Negative controls: same rule for the demix status; they are excluded from the aggregate anyway.
- Aggregate mode with given `.demix` files treats them as `ok` (there is no coverage information), except that `--previous_results` carries `low_confidence` over from the earlier run's `summary/demix_status.tsv`.

## Pieces

1. `bin/depth_stats.py`: `--floor-genome-coverage <pct>` (default none). New exit code 11 = coverage >= floor and < the gate
   (only when a floor is given); 10 stays "below the floor / below the gate when no floor". The JSON gains `floor_genome_coverage`
   and `low_confidence` (bool). `passes_gate` keeps meaning "at or above the gate".
2. `bin/demix_status.py`: add `low_confidence` to `STATUSES` (a reason is required, as for any non-ok status).
3. `modules/local/freyja_demix/main.nf`: pass `--floor-genome-coverage ${floor}` when `task.ext.floor_genome_coverage` is set;
   gate 11 -> still run `freyja demix`; success -> `write_status low_confidence "<cov>% of the genome covered at >= <covcut>x:
below the high-confidence threshold of <min>%; treat the estimate with care"`; a demix failure is still `failed`. Update the
   output comment (`.demix` present for ok and low_confidence). Keep the stub as is. FREYJA_DEMIX_COMPARE shares the module and the
   `ext` block (conf/modules.config), so it behaves the same.
4. Params (`nextflow.config`, `nextflow_schema.json`, typed block in `main.nf`, `conf/modules.config`): `demix_floor_coverage`
   (number 0-100, default 25), `include_low_confidence` (boolean, default false). `docs/parameters.md` must match the schema
   exactly (a test checks it). `ext.floor_genome_coverage = { params.demix_floor_coverage }`.
5. Exclusion from the time-based analyses: new tiny local module `FILTER_AGGREGATED` (`bin/filter_aggregated.py`: aggregated TSV +
   demix_status.tsv + `--exclude-status low_confidence` -> same-format TSV without those samples; first column is the sample
   label, keep the header and column order; if nothing is left write the header only). In `workflows/wasteline.nf` use its output
   (instead of `FREYJA_AGGREGATE_PLOT.out.aggregated`) for TIMESERIES_TABLE(x2), FREYJA_RELGROWTHRATE and FREYJA_DASH, only when
   `!params.include_low_confidence` AND the band is enabled; otherwise wire exactly as today. The two-dates test (`ch_ts_go`) must
   not count excluded samples. `freyja_aggregated.tsv` itself stays unfiltered. Module needs main.nf, environment.yml (same
   python/pandas-free stdlib approach as the other local python modules: copy the closest one), a `stub:` block, a versions
   topic entry and meta.yml if the sibling modules have one. Check `conf/modules.config` for a publishDir need (none: internal).
6. `bin/site_aggregate.py`: `--include-low-confidence` flag; with it `low_confidence` counts as ok, without it the sample is
   reported as not ok (reason text mentions low confidence). Pass the flag via `ext.args` in `conf/modules.config` when the param
   is true.
7. `bin/lineage_groups.py` (around the `state == "ok"` branch): `low_confidence` fills the values like `ok` but keeps its own
   status text and its reason in the row. xlsx gets the same rows (check that nothing there assumes `ok`).
8. Report: `bin/build_report_data.py` (`DEMIX_STATUSES`, every `status == "ok"` that decides whether lineages are shown or used:
   low_confidence must behave like ok there; new INFO flag `demix_low_confidence`; the skipped-sample flag text must not mention
   them), `bin/render_report.py` and `assets/templates/report.html.j2` (badge class `warn`, label `low confidence`, the demix counter
   tile, the status legend/tooltips), `docs/report_data.md`. Check `bin/qc_summary.py` copies the status as text (no whitelist).
9. Docs: `docs/usage.md` (a short section "Low-confidence estimates"), `docs/parameters.md`, `docs/output.md` (status vocabulary,
   which files include/exclude them), `docs/troubleshooting.md` (why a sample is `low_confidence`), `docs/report.md`,
   `docs/PLAN_STATUS.md`, `CHANGELOG.md`; the D18 row in `docs/PLAN.md` is already written (adjust only if you deviate, and say why).

## Tests: keep them MINIMAL (the user thinks the suite is far too big)

One focused test per new behaviour, no more: depth_stats band (exit 11), filter_aggregated, lineage_groups with a low_confidence
row, build_report_data flag + status, site_aggregate include/exclude. Module nf-test (stub) for FILTER_AGGREGATED only if every
sibling module has one (it must, by CLAUDE.md). One e2e step, no new synthetic data (see below). Update existing snapshots only
where they really change.

e2e (tests/e2e.sh + tests/check_outputs.py): demix mode on the finished `out_full` run with
`--min_genome_coverage_for_demix 98 --demix_floor_coverage 50` (+ the usual `--previous_results`, `--primer_scheme midnight`,
`--barcodes`, lineage summary off as in the sibling steps): Site_1_week_2 (93.1 %) and NEG_CONTAM (94.5 %) become
`low_confidence`, MIX_A/MIX_B/PURE_C stay ok, LOWCOV_D/NEG_E stay skipped. Assert: demix_status.tsv statuses; Site_1_week_2 present
in freyja_aggregated.tsv; absent from freyja_timeseries.tsv; the report mentions it with the low-confidence badge. Then the same
run with `--include_low_confidence true` (`-resume` is fine): Site_1_week_2 present in freyja_timeseries.tsv.
