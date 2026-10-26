# S2 U6b: the panel overview, report back-links and the panel documentation

Branch `feature/multi-target-in-pipeline`. Follows `multi_target_u5a.md`. Historical note for the module authors; the user documentation is
`docs/usage.md` (Panels) and `docs/other_pathogens.md` (Panels).

## What was added

- `bin/build_overview.py` + `assets/templates/overview.html.j2` + `modules/local/build_overview` (`BUILD_OVERVIEW`, the same python_jinja2 image as
  `RENDER_REPORT`; its jinja2 version is a plain emit, not `topic: versions`, for the same reason as `RENDER_REPORT`: `PROVENANCE` needs every topic
  entry before it writes `run_info.json`, and the overview runs after the reports). Published to `summary/overview.html`
  (`conf/modules.config`; `ext.when = !params.skip_report`).
- Inputs: `summary/targets.tsv` and a manifest (tab separated `target, role, path`, roles `target|aggregated|qc|demix_status|run_info|report`)
  with the per-target files staged as `in/<n>/<name>` (`path inputs, stageAs: 'in/?/*'`). The `target` role lists the columns of the matrix (every
  target of the run, also one with no report); `report` marks a target whose `<target>/report.html` exists (the cells link to it). The workflow
  function `overviewInputs(target_names, items)` builds the manifest (sorted by target, written by `workTextFile`) and the file list.
- Wiring (`workflows/wasteline.nf`): mode full with `--panel`: `ch_report_parts.join(RENDER_REPORT.out.html)` -> `toList()` (the barrier; a run
  with no reported target still builds the overview) -> `combine(wrapItem(...))` with `PANEL_FRONT.out.targets` (combine would spread the list).
  Mode report on panel results: after the per-target `RENDER_REPORT`, from `summary/targets.tsv`, `summary/<t>/*` and `pipeline_info/<t>/run_info.json`
  of `--results_dir` (no overview, with one info line, when the results have no `targets.tsv`). Single-target runs: nothing changes.
- Report: `bin/build_report_data.py` adds `header.panel` (the panel's file name) only when `run_info.panel` is set; `report.html.j2` shows
  `target <name>` in the subtitle and a "&larr; overview" link to `../overview.html` only then. Single-target `report_data.json` and
  `report.html` are byte-identical (checked by rendering the same data with the previous template).

## Anonymisation of the overview

`--report_anonymize paths|full` is passed to `build_overview.py`. The panel is always shown by file name (never a path). With `full`
the sample IDs are replaced by pseudonyms numbered over ALL samples of `targets.tsv` (sorted), with the private
`summary/overview_sample_id_map.csv`. The reports' own maps (`summary/<target>/sample_id_map.csv`) cannot be reused: each report numbers only
its own samples, so the same sample has different pseudonyms in different targets' reports and a sample absent from a target has none.
The overview pseudonyms therefore can differ from those inside a report; this is documented in usage.md.

## Leftovers

- A report written before this version has no `header.panel`, hence no back-link, until it is rebuilt (`--mode report --report_anonymize paths`
  or without its `report_data.json`).
- In `--mode report` the overview columns are the selected targets only (`--targets`); links to reports of targets that were not rebuilt into
  the same `--outdir` would be dead, so they are not shown.
