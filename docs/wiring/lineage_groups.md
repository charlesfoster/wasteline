# Lineage groups: wiring notes

New files (all additive; nothing existing was edited):

| File | Role |
| --- | --- |
| `bin/lineage_groups.py` | The script (standard library only; runs in `quay.io/biocontainers/python:3.12.12`). |
| `modules/local/lineage_summary/` | `LINEAGE_SUMMARY` process, `environment.yml`, nf-tests (stub + real). |
| `tests/python/test_lineage_groups.py`, `tests/fixtures/lineage_groups/` | pytest suite and the hand-made hierarchy and samples. |
| `assets/lineage_groups.csv` | Best-effort example built from the user's spreadsheet header (interpretations below). |

## 1. The groups file (for people who maintain it in Excel)

One row per output column. Save as **CSV** (Excel: *Save As -> CSV UTF-8*). Order of the rows is the order of the
columns in the results. A European-locale Excel that saves with `;` between columns also works.

| Column | Required | What to put |
| --- | --- | --- |
| `group` | yes | The column heading, exactly as you want to see it. Spaces, `*`, brackets and commas are all fine. Must be unique. |
| `include` | yes | One or more lineage patterns, separated by `;`. |
| `exclude` | no | Patterns to take out again, separated by `;`. |
| `in_total` | no | `yes` (default) or `no`. Say `no` for a group that only repeats lineages already counted in another group. |

### Patterns

| Pattern | Matches |
| --- | --- |
| `NB.1.8.1*` | `NB.1.8.1` **and every descendant** (`RK.1.2`, `PQ.17.7.2`, ...). |
| `NB.1.8.1` | exactly that lineage, nothing below it. |
| `RECOMBINANT*` | every lineage that descends from a recombinant (`XBB`, `XDV`, `XEC`, `XFG`, and all their sub-lineages such as `NB.1.8.1`, `RK.1.2`, `XFG.1`). |

* Upper or lower case, spaces around names: all ignored (`  jn.1* ; kp.3* ` works).
* Use either the short name or the full Pango name: `NB.1.8` and `XDV.1.5.1.1.8` are the same lineage, and
  `KP.3` is `B.1.1.529.2.86.1.1.11.1.3`. Descendants are found through the parent links of the Freyja dataset, so
  the alias does not matter.
* A lineage belongs to a group when it matches an `include` pattern and **no** `exclude` pattern.
* Groups may overlap on purpose (that is how "includes JN.1*" next to "excludes KP.2*, KP.3*" works): every group
  column is the sum of the sample's lineages that match it.
* Do not write `X*` for recombinants: it means a lineage called `X`, which does not exist. Write `RECOMBINANT*`.
* A name that is not in the dataset (typo, or a lineage newer than the dataset) is reported **once** in the run
  log and in `lineage_groups_resolved.tsv`, and matches nothing. There is no guessing from the spelling of the name.

### Worked example

```
group,include,exclude,in_total
JN.1* (excludes KP.3*),JN.1*,KP.3*,yes
JN.1* (includes KP.3*),JN.1*,,no
KP.3*,KP.3*,,yes
Recombinants (excludes NB.1.8.1*),RECOMBINANT*,NB.1.8.1*,yes
NB.1.8.1*,NB.1.8.1*,,yes
```

A sample with `JN.1.16 0.55`, `KP.3.1 0.25`, `RK.1.2 0.15`, `XDV 0.05` gives `JN.1* (excludes KP.3*)` 0.55,
`JN.1* (includes KP.3*)` 0.80 (informational, not in the total), `KP.3*` 0.25, `NB.1.8.1*` 0.15,
`Recombinants (excludes NB.1.8.1*)` 0.05, `other` blank, `proportion_summary` 1.00.

### Common mistakes

* Two groups both `in_total = yes` that share a lineage (`JN.1*` and `KP.3*` without excluding KP.3 from JN.1): the
  total goes above 1. The run log says which lineages overlap; set one of the groups to `no` or add an `exclude`.
* A comma inside a pattern cell: use `;` between patterns (a `,` also works, but only inside a quoted cell).
* Saving as "CSV (Macintosh)" or "Text (tab delimited)" instead of CSV UTF-8: the first accented character or the
  columns will go wrong. A tab-delimited file is detected, but UTF-8 CSV is the safe choice.
* Two rows with the same `group`, an empty `include`, or `in_total` other than yes/no: the run stops with a message
  naming the row (nothing is silently dropped).
* Group names `sample`, `site`, `other`, `status`, ... collide with the fixed output columns and are refused.

## 2. Outputs

`summary/lineage_summary.csv` (comma separated, UTF-8 without BOM, quoted where needed; opens in Excel). One row per
sample of `demix_status.tsv` (negative controls per the samplesheet are left out; samples missing from
`demix_status.tsv` but present in the aggregate are added as `ok`):

`sample, site, collection_date, status, genome_coverage_pct, <one column per group, file order>, other,
proportion_summary, major_lineages, major_abundances, minor_lineages, minor_abundances, ungrouped_lineages,
status_reason`

* Group values are abundances (0-1) summed over **all** lineages in the sample's Freyja solution that match the
  group; blank when zero (`--zeros zero` writes `0.000000`). `--min-abundance` (0.001) never changes these sums.
* `other` = lineages that match no `in_total` group. `proportion_summary` = sum of the `in_total` groups + `other`.
  It is ~1 when no two `in_total` groups overlap, but note that Freyja's own `lineages` list sums to 0.98-0.997 (it
  drops the tail below its cutoff), so expect 0.98-1.00 on real data, not exactly 1.
* `major_*` = lineages at >= 10% (`--major-threshold`, exactly 0.10 counts as major), `minor_*` = the rest down to
  `--min-abundance`; names space separated, descending abundance (ties alphabetical), abundances in the same
  order. Only the specific lineages (`lineages`/`abundances`) are used, never Freyja's `summarized` groups.
* `ungrouped_lineages` = the lineages behind `other` as `NAME:abundance`; this is where maintainers see which new
  variants to add to the groups file.
* A sample that is `failed`, `skipped_low_coverage` or `missing` still gets a row: `status`/`status_reason`
  filled, all value columns blank. (`status_reason` is an addition to the requested column list; drop it if unwanted.)
* `genome_coverage_pct` comes from `demix_status.tsv` (else the aggregate's `coverage`), rounded to 2 decimals.

`summary/lineage_groups_resolved.tsv` (audit): `group, include, exclude, in_total, n_lineages_matched, warnings,
dataset_version, groups_sha256` (one row per group; the dataset version is `last_barcode_update.txt`, the hash is of the
groups file as supplied).

`summary/lineage_summary_site_aggregates.csv` (only with the site-aggregate inputs): `aggregate_name, period,
period_start, n_sites_included, sites_included, sites_missing, weight_basis, <groups...>, other, proportion_summary,
major_lineages, major_abundances, minor_lineages, minor_abundances, ungrouped_lineages`, one row per row of the
meta TSV (an aggregate with no row in the Freyja-format file stays, with blank values, and a warning).

## 3. How the hierarchy is read (findings from the real Freyja dataset)

`lineages.yml` (6015 entries in `2026-10-06T05-25`) has, per lineage: `name`, `alias` (the full un-aliased Pango
name, e.g. `RK.1.2` -> `XDV.1.5.1.1.8.1.17.7.2.1.2`), `children` (**all** descendants plus itself, already
transitive), `parent` (immediate parent, by *name*, so it crosses aliases correctly: `RK.1.2 -> RK.1 -> ... ->
NB.1.8.1 -> NB.1.8 -> NB.1 -> XDV`) and, on roots only, `recombinant_parents`.

* The roots are `A`, `B` and the 189 recombinants; **every recombinant root (and only those) has
  `recombinant_parents`**, and no entry has both `parent` and `recombinant_parents`. `RECOMBINANT*` therefore means "the
  root of the parent chain has `recombinant_parents`" (2191 lineages), which does not rely on names starting with `X`
  (`NB`, `RK`, `PQ` are recombinant descendants with ordinary two-letter aliases).
* Pango aliasing is resolved by the parent links, not by name prefixes: `PQ.17.7.2` is in `NB.1.8.1*` and **not** in
  `JN.1*`, `JN.1.16.3` is in both `JN.1*` and `BA.2.86*`, `KP.3.1.1` is in `KP.3*`, `JN.1*` and `BA.2.86*`, `XFG.1` is in
  `XFG*`, `BA.3.2.2` in `BA.3.2*`, and the recombinants (`XBB`, `XDV`, `XFG`, `XEC`) are **not** descendants of
  `B.1.1.529`, so `Omicron (B.1.1.529*)` never contains them. A pattern is matched by the name *or* the alias of a
  lineage, case-insensitively.
* `children` and the `parent` closure agreed exactly for every entry; the script unions both.
* `yaml` is not available in the Python container, so `lineage_groups.py` has a small reader for exactly this file
  layout (checked equal to PyYAML on the real file and on the fixture).
* A lineage in the aggregate that is not in the dataset (dataset older than the barcodes, custom barcodes) is warned
  about once and counted as `other` (and listed in `ungrouped_lineages` when above `--min-abundance`). A dataset that
  is only a barcodes `.csv`/`.feather` has no hierarchy: the process stops with a message, so the workflow should only
  call it when `FREYJA_DATA.out.dataset` is a directory containing `lineages.yml` (as for the dropout/unexplained
  tables, `mutationTablesUnsupportedReason()`), and log a warning otherwise.

## 4. `assets/lineage_groups.csv`: interpretations to confirm

Built from the group names in the user's sheet. Rows that are my reading, not something the spreadsheet states:

* `Omicron Lineage (B.1.1.529*)`: written as `B.1.1.529*` excluding `BA.2.86*` and `BA.3.2*`, i.e. the "old Omicron"
  remainder (BA.1, BA.2, BA.4/5, ...). If the user's column is the informational whole-Omicron total, set it to
  `in_total = no` and drop the exclude.
* `Recombinants (excludes XEC*, NB.1.8.1*, XFG*)`: `RECOMBINANT*` minus those three, so `XDV` (not NB.1.8.1), `XBB*`,
  `XGA`, etc.
* `JN.1* (excludes KP.2*, KP.3*, KW.1.1*, LP.8.1*)` is `in_total = yes` and the `(includes ...)` twins (and
  `BA.2.86* (includes JN.1*)`) are `no`.
* `PE.1.4*` is `in_total = no`: PE.1.4 descends from KP.3 in this dataset, so counting it as well as `KP.3*` would
  double count. If it should be a separate total column, add `PE.1.4*` to the `exclude` of `KP.3*` and set it to `yes`.
* `NB.1.8.1* (XDV.1.5.1.1.8*)`: the text in brackets is only the column label; the pattern is `NB.1.8.1*`.
* `BA.2.86*`, `KP.2*`, `KP.3*`, `KW.1.1*`, `LP.8.1*`, `XEC*`, `XFG*`, `BA.3.2*` are literal.

Applied to the user's real run1/run2 aggregates with the real hierarchy: no unknown lineages, no overlapping
`in_total` groups, `proportion_summary` 0.980-0.998 (equal to the sum of Freyja's own `abundances`; the shortfall is
Freyja's trimmed tail, not a grouping loss); the only `other` was `B.52` (0.0010) in `WW-1277E`.

## 5. Wiring (planner)

### Params

```groovy
// nextflow.config params { }
lineage_groups          = null   // path: groups CSV; the summary is skipped when not set
lineage_major_threshold = 0.10   // abundance at/above which a lineage is "major"
lineage_min_abundance   = 0.001  // lineages below this are not listed (group sums still include them)
lineage_summary_zeros   = 'blank' // 'blank' | 'zero': how a group with no abundance is written
```

Schema: `lineage_groups` `format: file-path`, `pattern: ^\S+\.(csv|tsv|txt)$` (exists: true); `lineage_major_threshold` and
`lineage_min_abundance` `number`, `minimum 0`, `maximum 1`; `lineage_summary_zeros` `enum ['blank','zero']`.

### Include / publishDir

```groovy
// workflows/wasteline.nf
include { LINEAGE_SUMMARY } from '../modules/local/lineage_summary/main'

// conf/modules.config
withName: 'LINEAGE_SUMMARY' {
    errorStrategy = 'ignore'   // a typo in the user's CSV must not abort a finished run; the error is in .command.err
    ext.args = { "--major-threshold ${params.lineage_major_threshold} --min-abundance ${params.lineage_min_abundance} --zeros ${params.lineage_summary_zeros}" }
    publishDir = [
        path: { "${params.outdir}/summary" },
        mode: params.publish_dir_mode,
        pattern: 'lineage_*'
    ]
}
```

(`lineage_summary.csv`, `lineage_groups_resolved.tsv` and `lineage_summary_site_aggregates.csv` all start with
`lineage_`. `errorStrategy = 'ignore'` is a judgement call: an invalid groups CSV is a configuration error, so
dropping the line to fail the run loudly is also defensible. The error text names the file and row either way.)

### Signature

```groovy
LINEAGE_SUMMARY(
    FREYJA_AGGREGATE_PLOT.out.aggregated,   // path  freyja_aggregated.tsv
    ch_demix_status,                         // path  summary/demix_status.tsv (value channel, e.g. .first())
    ch_samplesheet.first(),                  // path  samplesheet in effect; [] is allowed
    channel.value(file(params.lineage_groups, checkIfExists: true)),
    FREYJA_DATA.out.dataset,                 // path  Freyja dataset directory (value channel)
    ch_site_aggregates_freyja ?: channel.value([]),   // optional: site_aggregates_freyja.tsv
    ch_site_aggregates_meta   ?: channel.value([])    // optional: site_aggregates_meta.tsv
)
// outputs: .summary, .resolved, .sites (optional, empty channel when no site aggregates)
```

### Placement

After `FREYJA_AGGREGATE_PLOT`, wherever it runs (`full`, `demix`, `aggregate`),
inside `if (params.lineage_groups) { ... }`. The three value inputs (`demix_status`, `samplesheet`, `groups`, `dataset`)
must be value channels so the single aggregate fires the process once. Skip it, with a `log.warn`, when the dataset is
not a directory with `lineages.yml`. Run it after the site-aggregation step when the latter is enabled, passing its
Freyja-format table and meta TSV as the last two inputs (else `channel.value([])` for both; passing only one of them is
an error).

### Report suggestion

A "Lineage groups" section in the HTML report, shown only when `lineage_summary.csv` exists: a table with the same
columns (samples as rows; the group columns as a heat-coloured numeric block; `major_lineages` with their abundances
as `name (xx%)`), the status/reason on non-`ok` rows, and a note under the table listing `ungrouped_lineages` across
samples ("Lineages not in any group: ... add them to the groups file") plus any `warnings` from
`lineage_groups_resolved.tsv` and the dataset version / groups hash in the footer. `report_data.json` can simply
embed the rows of the CSV.
