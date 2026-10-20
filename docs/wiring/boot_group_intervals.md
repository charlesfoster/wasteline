# Group-level bootstrap intervals: wiring notes

With `--freyja_boot` the pipeline already writes per-LINEAGE 95 % intervals (`summary/freyja_boot_summary.tsv`). The lineage
summary reports abundances per user-defined GROUP, so this adds the interval of each group:
`summary/lineage_group_intervals.tsv`. The confirmed layout of `lineage_summary.csv/.xlsx` is NOT changed (a separate file).

## What `freyja boot --rawboots` writes (verified in the pinned Freyja 2.0.6 container)

With `--output_base S.boot --rawboots` Freyja writes, besides the quantile tables it writes today (`S.boot_lineages.csv`,
`S.boot_summarized.csv`), the RAW draws `S.boot_lineages_boot.csv` and `S.boot_summarized_boot.csv`: a CSV whose first column is
an unnamed draw index (0..nb-1), one column per lineage, value = abundance in that draw, absent lineages filled with 0. Each
row sums to ~1. We use `*_lineages_boot.csv` (lineage level: group membership is decided by lineage name).

## Method

Per sample and per draw: assign every lineage to the groups exactly as `bin/lineage_groups.py` `summarise()` does
(`hierarchy.canonical(name)` in `group["members"]`; a lineage can count in several groups; `other` = lineages in no `in_total`
group), sum the draw's abundances per group, then take the 2.5 and 97.5 percentiles over the draws (linear interpolation = pandas
`quantile`, the same definition `freyja_boot_summary.py` uses for the lineage level). The point estimate is the group sum of the
sample's own demix (the value in `lineage_summary.csv`), not the bootstrap mean. IMPORT the helpers from `lineage_groups.py`
(`read_groups`, `resolve_groups`, `load_dataset`, `Hierarchy`, `read_aggregate` ...): do not copy the matching logic and do not
change its behaviour.

## Pieces

1. `bin/boot_group_intervals.py`: `--aggregated freyja_aggregated.tsv --groups groups.csv --dataset <dir/file> --boots <*_lineages_boot.csv ...>
--output lineage_group_intervals.tsv [--decimals N]`. Sample ID = file name minus `.boot_lineages_boot.csv` (sample IDs are
   never rewritten). Long TSV, header `sample group abundance lo hi n_boot`, one row per sample x (every group + `other`), sorted by
   sample then the groups' CSV order. A sample with no usable boot file is simply absent (not an error). If the dataset has no
   `lineages.yml` (custom barcode datasets) or the groups CSV is unusable: exit non-zero with the same `UserError` messages as
   `lineage_groups.py`; the workflow only runs it when the lineage summary itself runs.
2. `modules/local/freyja_boot/main.nf`: when `task.ext.args` contains `--rawboots` the raw file `${prefix}.boot_lineages_boot.csv`
   is an extra OPTIONAL output (`emit: raw`); the failure path removes it. Stub writes one when the flag is present.
3. `conf/modules.config`: FREYJA_BOOT adds `--rawboots` when the lineage summary runs (`params.lineage_summary`; the workflow
   decides the exact condition). `main.nf`/schema: no new parameter.
4. New local module `BOOT_GROUP_INTERVALS` (main.nf, environment.yml, stub, versions topic, nf-test like its siblings, container
   pinned like LINEAGE_SUMMARY, `bin/lock_containers.py --update` then `--check`). publishDir `summary/` for
   `lineage_group_intervals.tsv`.
5. `workflows/wasteline.nf`: after FREYJA_BOOT and under the same condition as LINEAGE_SUMMARY (`run_lineage_groups &&
!lineage_reason`) and `params.freyja_boot`: collect all raw boot files; call BOOT_GROUP_INTERVALS with the aggregated table,
   `lineageGroupsFile()` and `FREYJA_DATA.out.dataset` (the same inputs LINEAGE_SUMMARY gets). Never fatal to the run (errorStrategy
   `ignore` like FREYJA_BOOT).
6. Report: ONE sentence in the lineage-group provenance line is enough only if trivially available; otherwise no report change.
   Do not add a report section.
7. Docs: `docs/usage.md` (bootstrap/lineage-group sections), `docs/output.md` (the new file and its columns), `docs/parameters.md`
   (only if text of `--freyja_boot` changes: keep the schema/parameters two-way match), `docs/PLAN.md` decision D19 (below),
   `docs/PLAN_STATUS.md`, `CHANGELOG.md`.

PLAN D19 text: "**Group-level bootstrap intervals**: with `--freyja_boot` the lineage-group summary gets 95 % intervals in
`summary/lineage_group_intervals.tsv` (long format), from the RAW Freyja bootstrap draws (`freyja boot --rawboots`): lineage
abundances are summed into groups draw by draw (same membership rules as the lineage summary) and the 2.5/97.5 percentiles
taken, so group intervals are correct even though lineages within a group are negatively correlated (never add per-lineage
intervals). A separate file keeps the confirmed `lineage_summary.csv` layout unchanged."

## Tests (MINIMAL: the user considers the suite far too big)

One unit test file for the script with: a hand-computed case (two groups, an `other`, a lineage in two groups, 5 draws, known
percentiles), and the failure on a dataset without `lineages.yml`. The module nf-test (stub). One e2e addition only if
`--freyja_boot` is already exercised in tests/e2e.sh (extend that step: the new TSV exists, has the header, one row per
sample x group; for the sample's `other`/groups the point estimates equal the lineage_summary.csv values). Do not add new synthetic data.
