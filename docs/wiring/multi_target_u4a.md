# S2 U4a: per-target run-level backbone (summaries, provenance, report)

Branch `feature/multi-target-in-pipeline`. Follows `multi_target_s2_plan.md` (U4). U4b does the rest (list at the end).

## What is per target now

Run-level steps run once per target, for panels and for normal runs (a normal run is a panel of one: `summary_sub` is `''`, so
layout, names and values are unchanged):

| Step | Input grouping | Output |
| --- | --- | --- |
| `FREYJA_AGGREGATE_PLOT`, `FREYJA_AGGREGATE_CONTROLS` | `.demix` files of the target (`groupSamples`, sorted by sample ID) | `summary<sub>/freyja_aggregated*.tsv`, plots |
| `FILTER_AGGREGATED` | aggregate joined to the target's `demix_status.tsv` on the key | (feeds the U4b steps) |
| `QC_CONCAT`, `CONSENSUS_STATS_CONCAT`, `CONSENSUS_CONCAT`, `MUTATIONS_OF_INTEREST_CONCAT` | per-sample files of the target (`groupByTarget`, sorted by name) | `summary<sub>/qc_summary.tsv`, `amplicon_coverage.tsv`, `consensus_stats.tsv`, `consensus_genomes.fasta`, `mutations_of_interest.tsv` |
| `PROVENANCE` | metas of the target's samples; that target's primer BEDs, FASTA, `freyja_data.json` | `pipeline_info<sub>/run_info.json` |
| `BUILD_REPORT_DATA`, `RENDER_REPORT` | one map (role -> file) per target, joined from every per-target table | `summary<sub>/report_data.json`, `report.html` |

`<sub>` is `meta.target.summary_sub` (`''` or `/<target>`), read in `conf/modules.config`.

The key (the `meta` of every one of these processes) is `tkey(meta) = [ id: target name, target: target map ]`. A target without
any item produces no group, so no process runs for it and it has no summary (an undetected target has none).
Sidecars (Freyja dataset, `freyja_data.json`) are delivered per target in lockstep with `combine(<*_by_target>, by: 0)` on the
target name plus one `multiMap`; shared single values (samplesheet, versions, plot config) are value channels (`.first()`).

Signature change (all other inputs unchanged, outputs became `[ meta, file ]`): `tuple val(meta), ...` is the first input of the 10
modules above (merged into the existing first tuple where there was one); `QC_CONCAT` takes `tuple val(meta), path(qc_rows), path(amplicon_tables)`.

## Files that were `collectFile` outputs

`collectFile` cannot write below a sub-directory (name `RSVa/demix_status.tsv` fails: "No such file or directory", also with `storeDir`)
and its result carries no key, so the per-target files are assembled in Groovy from the `groupTuple` group and returned as
`[ tkey, file ]` (workflow functions at the end of `workflows/wasteline.nf`):

- published ones (`demix_status.tsv`, `barcode_version.json`, `samplesheet.csv`, `contig_coverage.tsv`, `consensus_stats_contigs.tsv`,
  `annotated_variants.tsv`): `summaryFile(outdir, key, name, text)` writes `<outdir>/summary<sub>/<name>` (what `storeDir` did) and
  rewrites only when the text changed;
- internal ones (`times.csv`, the failed-consensus rows, the base run info): `workTextFile(name, text)` under
  `<workDir>/wasteline_run/<digest>/<name>`, again only rewritten when the content changes, so `-resume` is not disturbed (same idea as
  `panelResolvedFile` in `panel_front`).

Failed-consensus rows are files named by `meta.key` (`id` or `id__target`), so two targets of one sample never share a name.

## Provenance per target

`runInfoBase(mode, metas, target)`: metas and facts (FASTA, contigs, `freyja_contig`, default primer scheme) are those of the target;
`target: <name>` is recorded in every `run_info.json`, `panel: <spec>` for panels (the run's `params` stay the run's). Lineage-group
info is not recorded for panel targets until U4b. `finaliseRunInfo` (completion) adds containers/duration to every
`pipeline_info<sub>/run_info.json`. `bin/build_report_data.py` copies `run_info.target` into `header.target`.

## Still gated for panels (U4b), kept working unchanged for one target

One `log.warn` lists them: `UNEXPLAINED_MUTATIONS_CONCAT`, `DROPOUT_IMPACT`, `LINEAGE_GROUPS_CHECK`, `LINEAGE_SUMMARY`,
`SITE_AGGREGATE(_CHECK)`, `BOOT_GROUP_INTERVALS`, `TIMESERIES_TABLE(_LINEAGES)`, `FREYJA_RELGROWTHRATE`, `FREYJA_DASH`,
`--compare_barcodes` (`FREYJA_DEMIX_COMPARE`, `BARCODE_COMPARE`, `barcode_compare*.tsv`), `freyja_boot_summary.tsv`. In `workflows/wasteline.nf`
they sit in the `else` of `if (is_panel)` after "The steps that are not per target yet" and read plain channels derived from the keyed
ones (`.map { _key, f -> f }`). The report gets `[]` for their tables for a panel target (one shared item, see `ch_report_shared`).
`--mode report` still builds one report (meta from `targetFromParams()`); the per-target report modes are U5.
