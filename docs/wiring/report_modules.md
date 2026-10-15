> **Design-era note.** Written before final integration; the implemented behaviour is in the workflow code, `docs/output.md` and `docs/report.md`.

# Wiring notes: report modules (from the report-modules task; wire in the integration round)

Modules: `modules/local/build_report_data`, `modules/local/render_report` (7/7 nf-tests pass, stub + real Docker).

## Includes
```
include { BUILD_REPORT_DATA } from '../modules/local/build_report_data/main'
include { RENDER_REPORT     } from '../modules/local/render_report/main'
```
(adjust the relative path to the including file)

## Signatures
- `BUILD_REPORT_DATA` input[0]: `tuple path(run_info), path(qc_summary), path(amplicon_coverage), path(freyja_aggregated),
  path(demix_status), path(consensus_stats), path(mutations_of_interest), path(barcode_compare), path(depths)`.
  - Only `run_info` is required; pass `[]` for every absent file. `depths` is one file or a list of `<sample>.depths`
    (unique names), not a directory.
  - Keep the standard file names: `qc_summary.tsv`, `amplicon_coverage.tsv`, `freyja_aggregated.tsv`, `demix_status.tsv`,
    `consensus_stats.tsv`, `mutations_of_interest.tsv`, `barcode_compare.tsv`.
  - Every run-level channel must emit exactly ONE item (`.collect()` / `.ifEmpty([])`), otherwise the tuple multiplies
    or vanishes.
  - Outputs: `.out.json` (`report_data.json`) + versions topic.
- `RENDER_REPORT` input[0] `report_data.json`, input[1] `file("${projectDir}/assets/templates/report.html.j2")`.
  Outputs: `.out.html` (`report.html`) + versions topic.

## conf/modules.config
```
withName: 'BUILD_REPORT_DATA' {
    ext.args   = { "--mode ${params.mode}" }
    ext.when   = { !params.skip_report }
    publishDir = [ path: { "${params.outdir}/summary" }, mode: params.publish_dir_mode, pattern: 'report_data.json' ]
}
withName: 'RENDER_REPORT' {
    ext.when   = { !params.skip_report }
    publishDir = [ path: { "${params.outdir}/summary" }, mode: params.publish_dir_mode, pattern: 'report.html' ]
}
```

## Params
- `skip_report` (boolean, default false); ignored/rejected in `report` mode. `mode` enum must include `report`.
- The report reads thresholds from `run_info.json` `params`: `amplicon_min_depth`, `consensus_min_depth`,
  `ntc_max_mapped_reads`, `ntc_max_genome_coverage_pct`, `min_genome_coverage_for_demix`, `freyja_*`. The run_info base
  JSON must carry them.

## Per mode (everything not listed is `[]`)
| mode | BUILD_REPORT_DATA inputs |
|---|---|
| full | run_info, qc_summary, amplicon_coverage, freyja_aggregated, demix_status, consensus_stats, mutations_of_interest (if given), all `*.depths` |
| demix | run_info, demix_status, freyja_aggregated, `*.depths` from the samplesheet |
| aggregate | run_info, demix_status, freyja_aggregated (depths only if supplied) |
| report | see below |

`demix_status.tsv` is produced in the workflow from the `FREYJA_DEMIX.out.status` JSONs (already done in Phase 1; its
columns are a superset of the contract).

## `report` mode (`--results_dir`)
1. If `<results_dir>/summary/report_data.json` exists: skip BUILD_REPORT_DATA, feed it to RENDER_REPORT with the pipeline's
   own template.
2. Else rebuild with BUILD_REPORT_DATA from what exists: `pipeline_info/run_info.json` is REQUIRED (clear error if
   missing, never invent one); `qc_summary.tsv`, `amplicon_coverage.tsv`, `freyja_aggregated.tsv`, `demix_status.tsv`,
   `consensus_stats.tsv`, `mutations_of_interest.tsv` from `summary/`; depths from `samples/*/freyja/*.depths` (legacy:
   `<s>/freyja/<s>.depths`).
3. `--outdir` must differ from `--results_dir`. No samplesheet or Freyja data needed.
4. `consensus_stats.tsv` MUST be published to `summary/` in full mode (see PLAN §6), otherwise a rebuild loses the
   consensus section.

## Container note
RENDER_REPORT runs in a Seqera Wave image built from `modules/local/render_report/environment.yml` (Python 3.12.12 +
Jinja2 3.1.6; see `docs/development.md`, "Rebuilding the report renderer image"). An earlier design used a MultiQC image
for its Jinja2 (about 4 times larger); that was replaced.
