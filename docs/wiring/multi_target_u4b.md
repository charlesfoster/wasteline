# S2 U4b: the optional analyses per target (panel gate removed)

Branch `feature/multi-target-in-pipeline`. Follows `multi_target_u4a.md`. Every step that U4a left gated now runs once per target, for
panel runs and for normal runs (a normal run is a panel of one: `summary_sub` is `''`, outputs unchanged, guarded by the pipeline snapshots
and a byte-for-byte stub comparison). The panel `log.warn` and the `if (is_panel)` branch are gone.

## What is per target

Same pattern as U4a: leading `tuple val(meta)` (meta = `tkey`: `[ id: target name, target: map ]`) on the module, outputs `[ meta, file ]`,
`publishDir summary${meta.target?.summary_sub ?: ''}`, per-target sidecars in lockstep with `combine(..., by: 0)` on the target name plus one `multiMap`.

| Step | Per-target input | Notes |
| --- | --- | --- |
| `UNEXPLAINED_MUTATIONS_CONCAT` | `groupByTarget(UNEXPLAINED_MUTATIONS.out.tsv)` | |
| `DROPOUT_IMPACT` | the target's `amplicon_coverage.tsv` joined to its aggregate, its dataset | `ext.contig` = `meta.target.freyja_contig` |
| `LINEAGE_GROUPS_CHECK` | `[ key, groups file ]` and the target's dataset (`[]` when it has no hierarchy: format check only) | warnings are prefixed `[<target>]` for panels |
| `SITE_AGGREGATE`, `LINEAGE_SUMMARY` | the target's aggregate and a complete `demix_status.tsv` (header only if the target has no row) | the sites CSV, the summary samplesheet and the groups file are shared or per target as below |
| `BOOT_GROUP_INTERVALS` | aggregate, the target's raw bootstrap draws, groups, dataset | not run for a target none of whose samples was bootstrapped |
| `TIMESERIES_TABLE(_LINEAGES)` | aggregate (low-confidence filtered) joined to `ch_ts_go` | `ch_ts_go`: the target has >= 2 distinct collection dates among its non-negative samples |
| `FREYJA_RELGROWTHRATE`, `FREYJA_DASH` | the same aggregate, the target's dataset | `--pathogen` = `meta.target.freyja_pathogen` |
| `--compare_barcodes` | `FREYJA_DATA_COMPARE` per target, `FREYJA_DEMIX_COMPARE` per (sample, target), `BARCODE_COMPARE` with the two `barcode_version.json` of the sample's target | `summary<sub>/barcode_compare.tsv`, `barcode_compare_status.tsv` written by `summaryFile` |
| bootstrap summary | `groupByTarget(FREYJA_BOOT.out.summary)` | `summary<sub>/freyja_boot_summary.tsv` written by `summaryFile` |
| report | every table is now a role in the per-target map (`taggedByTarget`); only the samplesheets are shared | the `[]` placeholder item for panels is gone |

## Decisions

- Lineage groups: file = `--lineage_groups` (applies to every target of the panel; one info line says so) else the target's `lineage_groups_bundled`.
  A target with neither has its summary off (one info line naming it). The hierarchy check is per target: static reason (`lineageHierarchyUnsupportedReason(target)`,
  e.g. a barcodes CSV in the panel row) or, for a downloaded dataset that may lack it (`hierarchyMayBeMissing(target)`), `lineages_yml_sha256` of that
  target's `freyja_data.json`; the target is then skipped with ONE warning naming it. `mutationTablesUnsupportedReason` was already per target.
- `SITE_AGGREGATE_CHECK` runs once; `SITE_AGGREGATE` and `LINEAGE_SUMMARY` per target on the shared summary samplesheet.
- `--compare_barcodes` is one spec for the panel. Rejected at start-up (`validatePanelParameters`): `container` when a target is not the container's
  pathogen, and a dataset directory or file with more than one target. `latest` and version stamps resolve per target pathogen (every target is
  resolved at start-up). Note `latest` downloads the second dataset of every panel target, undetected ones included.
- `runInfoBase` records lineage-group info for panel targets too.
- `FREYJA_DATASET_INFO` still runs for every target, but `pipeline_info<sub>/freyja_data.json` (and `freyja_data_compare.json`) are now written by the
  workflow (`pipelineInfoFile`) only for targets that have something to summarise; the module's `publishDir` is kept for `--mode update_db` only.
  An undetected panel target no longer gets a `pipeline_info/<t>/` folder.

## Signatures

Leading `tuple val(meta), path(...)` on: `UNEXPLAINED_MUTATIONS_CONCAT`, `DROPOUT_IMPACT` (first input), `LINEAGE_GROUPS_CHECK`, `LINEAGE_SUMMARY`,
`SITE_AGGREGATE`, `BOOT_GROUP_INTERVALS`, `TIMESERIES_TABLE`, `FREYJA_RELGROWTHRATE`, `FREYJA_DASH`; their outputs became `[ meta, file ]`.
Utils: `lineageHierarchyUnsupportedReason(target)`, `hierarchyMayBeMissing(target)`, `lineageGroupsInfo(target)` (all default to `targetFromParams()`);
`validatePanelParameters` checks `--compare_barcodes`. `FREYJA_DATA` no longer emits the unkeyed `dataset`, `freyja_data`, `barcode_version`.
Workflow functions: `pipelineInfoFile`, `targetOutputFile` (`summaryFile` uses it).
