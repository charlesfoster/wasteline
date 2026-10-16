# DROPOUT_IMPACT: wiring notes

Amplicon dropout biases lineage abundance estimates when the dropped amplicon hides mutations that tell the
lineages apart. `DROPOUT_IMPACT` intersects the dropped amplicons of each sample with the mutations of the Freyja
dataset and the lineages Freyja detected in that sample. Script: `bin/dropout_impact.py` (stdlib only, definitions in
its docstring). Module: `modules/local/dropout_impact/` (container `quay.io/biocontainers/python:3.12.12`, the same
tag as `QC_CONCAT`). Tests: `tests/python/test_dropout_impact.py` (43), `modules/local/dropout_impact/tests/main.nf.test`
(stub + 4 real docker tests), fixture `tests/fixtures/dropout_impact/`. Nothing existing was edited; this file is the
hand-over.

## Definitions (exact)

- **mutation / site**: a single-nucleotide substitution from the dataset (`C241T`). One site = one distinct
  (position, alternate base); a position with two alternate bases is two sites. Deletions (`C11288-`, `-11288-9`,
  `ATG100A`), insertions (`+T123TA`, `A12AT`) and anything else unparseable are skipped, counted, and reported once on
  stderr (`skipped N mutation(s) ... (D deletion, I insertion, O other)`); never fatal. They are invisible to this
  analysis. (The Freyja 2.0.6 `lineage_mutations.json` contains only substitutions: 507,135 of them, 0 skipped.)
- **informative site**: carried by `n` of the `N` lineages of the dataset with `0 < n/N < --informative-max-fraction`
  (default 0.95). Fixed mutations (>= 95% of lineages) cannot separate lineages; a mutation carried by none is no
  lineage marker (only matters for the `.csv` form, where such columns exist).
- **detected lineage**: abundance `>= --min-abundance` (default 0.01) in the `lineages`/`abundances` columns of the
  sample's row of `freyja_aggregated.tsv` (never `summarized`). A detected lineage missing from the dataset carries no
  mutations (counted on stderr).
- **affects_detected_lineage**: a detected lineage carries `>= 1` **informative** site inside the dropped span.
  Informative-only was chosen on purpose: hiding a mutation that nearly every lineage has (or none has) cannot move the
  estimate from one lineage to another, so counting those would flag almost every dropout.
- Per-sample totals: `informative_sites_total` counts informative sites inside the **union of that sample's amplicon
  spans** (what the scheme can observe, not genome-wide; 8006 for midnight + the 2.0.6 dataset). `in_dropped` counts the
  union of the dropped spans, so overlapping dropped amplicons are not double counted (the per-amplicon
  `n_informative_sites` does count the shared sites in both rows).

## Outputs (both always written, header always present)

`dropout_impact.tsv`, one row per **dropped** amplicon, sorted by sample then start:

```
sample  amplicon  start  end  mean_depth  n_informative_sites  n_sites_in_detected_lineages  affects_detected_lineage  affected_lineages
```

`affected_lineages` = detected lineages carrying informative sites in the span, `name:count`, ordered by count
(descending, then name), at most 5, comma-joined (`XBB.1.5:9,BA.2:3,BA.5:3`); `NA` if none.

`dropout_impact_summary.tsv`, one row for **every sample of the amplicon table** (zeros when nothing dropped):

```
sample  n_dropped_amplicons  informative_sites_total  informative_sites_in_dropped  pct_informative_sites_lost  n_dropped_affecting_detected_lineages
```

`pct_informative_sites_lost` has 2 decimals (`NA` when total is 0).

Samples without an aggregate row (skipped, failed, negative controls, so absent from `freyja_aggregated.tsv`) still
appear: `n_sites_in_detected_lineages = NA`, `affects_detected_lineage = false`, `affected_lineages = NA`, and
`n_dropped_affecting_detected_lineages = NA` in the summary. A sample WITH an aggregate row but an empty lineage list
gets `0`, `false`, `NA` and `0`. Aggregate rows for samples absent from the amplicon table are ignored.

## Include line and signature

```nextflow
include { DROPOUT_IMPACT } from '../modules/local/dropout_impact/main'   // path relative to workflows/wasteline.nf
```

| process | inputs | outputs (emit) |
|---|---|---|
| `DROPOUT_IMPACT` | `path amplicons` (`amplicon_coverage.tsv`); `path aggregated` (`freyja_aggregated.tsv`); `path dataset` (a Freyja dataset directory, or a single barcodes file) | `impact` = `dropout_impact.tsv`; `summary` = `dropout_impact_summary.tsv`; `versions_python` (topic `versions`) |

All three inputs are value-like (single file/dir, queue channels of one item or value channels).

Dataset resolution (`--dataset`): a directory with a non-empty `lineage_mutations.json` uses it; a directory without it
falls back to its single `*barcode*.csv` (e.g. `usher_barcodes.csv`); a `.json` or `.csv` file is used as given; a
`.feather` file, or a directory holding only a `.feather` (and a missing/empty json), is an error (exit 2) telling the
user to pass a dataset directory. So `--barcodes <file>.feather` and zero-byte stub datasets make the process fail; see
"Failure handling".

## Workflow sketch

Runs after `QC_CONCAT` (amplicon table) and `FREYJA_AGGREGATE_PLOT` (aggregate), with the same resolved dataset that
`FREYJA_AGGREGATE_PLOT` receives (`FREYJA_DATA.out.dataset`, a value channel). Only in `full` / `demix` modes with the
`amplicon` protocol and a primer scheme (no amplicon table, nothing to do; never for `--protocol capture`).

```nextflow
// in workflows/wasteline.nf, after QC_CONCAT and FREYJA_AGGREGATE_PLOT
def ch_dropout_impact = channel.value([])
def ch_dropout_summary = channel.value([])
if (params.dropout_impact && use_amplicons && (mode == 'full' || mode == 'demix')) {
    DROPOUT_IMPACT(
        QC_CONCAT.out.amplicon_coverage,
        FREYJA_AGGREGATE_PLOT.out.aggregated,
        FREYJA_DATA.out.dataset
    )
    ch_dropout_impact  = DROPOUT_IMPACT.out.impact
    ch_dropout_summary = DROPOUT_IMPACT.out.summary
}
```

(`QC_CONCAT.out.amplicon_coverage` and `FREYJA_AGGREGATE_PLOT.out.aggregated` are optional outputs: when either is
absent the process simply does not run. Then pass `ch_dropout_impact.ifEmpty([])` / `ch_dropout_summary.ifEmpty([])` to
`BUILD_REPORT_DATA` through `wrapItem(...)` like the other optional report inputs, with new `--dropout-impact` /
`--dropout-summary` options.)

## `conf/modules.config` block

```groovy
    withName: 'DROPOUT_IMPACT' {
        ext.args   = { "--min-abundance ${params.dropout_min_abundance}" }
        errorStrategy = 'ignore'     // an unusable dataset (.feather only) must not abort the run; the report just omits the table
        publishDir = [
            path: { "${params.outdir}/summary" },
            mode: params.publish_dir_mode
        ]
    }
```

Both files are published to `summary/` (default `saveAs`; no `pattern` needed, there are no intermediates).
`--informative-max-fraction` is deliberately not a pipeline param (an analysis constant, 0.95); add it to `ext.args`
only if a user needs it.

## Params (suggested)

| param | default | meaning |
|---|---|---|
| `dropout_impact` | `true` | run the analysis whenever an amplicon table exists (it is simply skipped otherwise); `false` disables it |
| `dropout_min_abundance` | `0.01` | a lineage counts as detected at abundance `>=` this (same value as the script default) |

`nextflow_schema.json`: both under the QC/amplicon group (`dropout_impact` boolean, `dropout_min_abundance` number,
min 0, max 1). Document in `docs/parameters.md`, and the two outputs in `docs/output.md` (`summary/dropout_impact.tsv`,
`summary/dropout_impact_summary.tsv`, only with a primer scheme and the `amplicon` protocol).

## Failure handling

The script exits 2 with `error: ...` for an unusable dataset (`.feather` only, missing, malformed JSON, no lineages), an
aggregate without `lineages`/`abundances` columns, or an amplicon table missing columns. With `errorStrategy 'ignore'`
the run continues and the outputs are absent, so the report must treat both inputs as optional (`[]`), as for the other
optional tables. Alternatively the workflow can skip the process up front when `params.barcodes` ends in `.feather`
and log why. Warnings (skipped indel notations, unknown detected lineages, label/length mismatches) go to stderr only.

## How the report should use it

1. **Table "Amplicon dropout impact"** directly under the amplicon dropout heatmap (visible only when
   `dropout_impact.tsv` exists). One row per dropped amplicon with a detected-lineage effect first; columns: sample,
   amplicon, span, mean depth, informative sites, sites in detected lineages, affected lineages. A per-sample caption
   from the summary: `pct_informative_sites_lost` % of lineage-informative sites lost in N dropped amplicons. Rows
   with `affects_detected_lineage = true` get a warn highlight; samples with `NA` (no aggregate row) say "no lineage
   estimate".
2. **Run-level flag** `dropout_affects_lineages`, level `warn`: raised when any sample has
   `affects_detected_lineage == true`. `samples` = those sample IDs (sorted), message e.g. *"Amplicon dropout hides
   lineage-informative sites for detected lineages in N sample(s); abundance estimates may be biased."* Add the code to
   the `run_flags` codes list in `docs/report_data.md`.
3. **`report_data.json` fields** (new, all optional, `{}`/absent when the table is absent):
   - `available.dropout_impact`: boolean.
   - `dropout_impact`: object `sample -> {summary: {n_dropped_amplicons, informative_sites_total,
     informative_sites_in_dropped, pct_informative_sites_lost, n_dropped_affecting_detected_lineages}, amplicons:
     [{amplicon, start, end, mean_depth, n_informative_sites, n_sites_in_detected_lineages,
     affects_detected_lineage, affected_lineages}]}`, `NA` -> `null`, `affects_detected_lineage` -> boolean,
     `affected_lineages` -> list of `{lineage, n_sites}`, amplicons sorted by start. Every sample of the summary appears
     (empty `amplicons` list when nothing dropped).
   - `thresholds.dropout_min_abundance` (from the run params).
   - optionally per-sample in `samples[*]`: `dropout_pct_informative_lost`, `dropout_affects_lineages` (boolean or null).

## What the real data showed (Freyja 2.0.6 dataset: 5,977 lineages; tests/data synthetic run)

Using the container's bundled `lineage_mutations.json` (10_06_2026-05-25) with the full e2e run's tables
(`summary/amplicon_coverage.tsv`, `freyja_aggregated.tsv`); the e2e run itself used the 14-lineage custom
`barcodes_test.csv`, which also works (`.csv` path, 244 informative sites). The whole analysis takes 0.6 s.

`Site_1_week_2` (truth: dropped `nCoV-2019_7` and `nCoV-2019_22`; detected `BA.5` 0.607, `XBB.1.5` 0.359, `BA.2` 0.034):

| amplicon | span | informative sites | in detected lineages | affected lineages |
|---|---|---|---|---|
| `nCoV-2019_7` | 6284-7401 | 318 | 0 | none (`NA`) |
| `nCoV-2019_22` | 21533-22612 | 532 | 10 | `XBB.1.5:9,BA.2:3,BA.5:3` (affects: true) |

Summary: 2 dropped amplicons, 850 of 8006 informative sites lost (10.62%), 1 dropped amplicon affecting detected
lineages. Hand check against `lineage_mutations.json` (counting carriers over all 5,977 lineages):

- span 6284-7401: `BA.5`, `XBB.1.5`, `BA.2` carry no mutation at all there, so 0 is right even though 318 informative
  sites (of other lineages) are hidden.
- span 21533-22612: `BA.5` and `BA.2` each carry `C21618T` (71.0% of lineages), `T22200G` (54.1%), `G22578A`
  (73.1%); `XBB.1.5` carries `C21618T`, `T21810C` (16.4%), `C22000A` (16.9%), `C22109G` (16.4%), `T22200A` (16.5%),
  `G22317T` (14.2%), `G22577C` (57.3%), `G22578A`, `G22599C` (46.8%). That gives `XBB.1.5:9`, `BA.2:3`, `BA.5:3`, and
  10 distinct sites in the union (the three shared `C21618T`/`G22578A` and the BA.2/BA.5 `T22200G` counted once).
  All are informative (< 95%): this is the RBD region, so the spike-defining sites of the detected lineages are hidden.

Other samples: `MIX_A`, `MIX_B`, `PURE_C` have no dropped amplicons (summary zeros, no impact rows); `LOWCOV_D` and
`NEG_E` (all 29 amplicons dropped, no aggregate row) get `NA` detected-lineage fields and 100% informative sites lost;
`NEG_CONTAM` demixed but is excluded from the aggregate, so it is `NA` too.

## Unverified / caveats

- Not run through the whole pipeline (wiring is the planner's); the process was tested alone (stub + docker).
- Deletion/insertion notations are only handled by skipping; no real dataset available here contains any, so the
  skip-and-count path is covered by the fixture only.
- Lineage names in the aggregate are assumed to match the dataset's names (both come from the same Freyja dataset);
  mismatches show up in the `not in the dataset` warning.
