# SITE_AGGREGATE: wiring notes (aggregate sites into groups at the result level)

Combine the per-site Freyja results into user-defined groups of sites (for example `Site1`, `Site2`, `Site3` into
`Metro_Sites`) **without pooling reads or variant tables**. Every site-sample is demixed on its own, exactly as today;
the finished abundance solutions of the member sites are then averaged with weights. Script: `bin/site_aggregate.py`
(stdlib only). Module: `modules/local/site_aggregate/` (container `quay.io/biocontainers/python:3.12.12`, the same tag as
`TIMESERIES_TABLE`). Tests: `tests/python/test_site_aggregate.py` (79), `modules/local/site_aggregate/tests/main.nf.test`
(stub + 3 real docker tests), fixtures `tests/fixtures/site_aggregate/` (read-only use of `tests/data_timeseries/`).
Example input: `assets/site_aggregation_example.csv`. Nothing existing was edited; this file is the hand-over.

## 1. Why result-level, not pooled (the evidence)

Pooling reads (or variant tables) from several sites into one library and demixing once has two problems: it
**implicitly weights each site by its sequencing depth** (a well-sequenced site dominates, whatever its share of the
population or of the virus), and the single fit on a mixture of mixtures can report lineages that **no** site contains.
Averaging the per-site solutions cannot invent a lineage (each output lineage is a weighted mean of values that are
non-negative and zero for sites that lack it), and the weights are explicit.

Experiment (real Freyja 2.0.6 container, `linux/amd64`, `tests/data/barcodes_test.csv`, the same `freyja demix` options as
the pipeline: `--eps 0.001 --covcut 10 --depthcutoff 0 --adapt 0.0 --solver CLARABEL`). Four finished synthetic samples of the
e2e run were treated as four sites: `MIX_A` (truth BA.2 0.7 / AY.4 0.3), `MIX_B` (BA.1 0.5 / XBB.1.5 0.3 / KP.3 0.2),
`PURE_C` (JN.1 1.0) and `Site_1_week_2` (BA.5 0.6 / XBB.1.5 0.4). Unequal depth was made by multiplying every depth and
read count of a site by an integer (frequencies are unchanged, so the site's own result is unchanged; mean depths 55x,
64x, 73x, 60x become 110x, 385x, 73x, 179x in the `unequal` rows). For each set:

* (i) **aggregate**: each site demixed alone, combined by `bin/site_aggregate.py` (equal weights, and "depth-proportional"
  weights only to mimic what pooling does implicitly; the latter is **not** a recommended basis, see section 3);
* (ii) **pooled variants**: ALT_DP / ALT_RV summed per (position, ref, alt) over the sites' ivar TSVs, TOTAL_DP and the depths
  file summed over the sites' depths, ALT_FREQ recomputed, one `freyja demix` of the pooled files;
* (iii) **truth**: the design fractions of `tests/data/truth.json` mixed with the same weights.

"max |diff|" is the largest absolute difference in abundance over lineages. (Sites' own Freyja errors, for example 2% spurious
BA.2 in `MIX_B`, are inside every number.)

| sites | depths | aggregate (equal w) vs truth | pooled vs truth (equal w) | aggregate (depth w) vs truth (depth w) | pooled vs truth (depth w) | pooled vs aggregate (depth w) | lineage in pooled result, in no site |
|---|---|---|---|---|---|---|---|
| MIX_A + PURE_C | equal | 0.006 | 0.050 | 0.005 | 0.022 | 0.021 | B.1 1.8% |
| MIX_A + PURE_C | unequal | 0.006 | 0.120 | 0.007 | 0.034 | 0.027 | B.1 0.7% |
| MIX_A + MIX_B + PURE_C | equal | 0.011 | 0.045 | 0.011 | 0.041 | 0.030 | B.1 0.8% |
| MIX_A + MIX_B + PURE_C | unequal | 0.011 | 0.201 | 0.017 | 0.051 | 0.033 | B.1 1.9% |
| MIX_A + MIX_B + Site_1_week_2 | equal | 0.023 | 0.076 | 0.023 | 0.089 | 0.066 | B.1 0.8% |
| MIX_A + MIX_B + Site_1_week_2 | unequal | 0.023 | 0.079 | 0.024 | 0.080 | 0.056 | B.1 1.2% |
| all four | equal | 0.017 | 0.062 | 0.016 | 0.058 | 0.043 | B.1 0.9% |
| all four | unequal | 0.017 | 0.152 | 0.021 | 0.068 | 0.047 | B.1 0.7% |

Worked rows (abundances; `all four`, unequal depths): aggregate with depth weights BA.1 0.258, XBB.1.5 0.229, BA.5 0.145,
BA.2 0.124, KP.3 0.103, JN.1 0.098, AY.4 0.042; pooled BA.1 0.232, XBB.1.5 0.204, BA.2 0.171, BA.5 0.145, KP.3 0.102, JN.1 0.098,
AY.4 0.041, **B.1 0.007**; truth (depth weights) BA.1 0.258, XBB.1.5 0.250, BA.5 0.144, KP.3 0.103, BA.2 0.103, JN.1 0.098, AY.4 0.044.

What this shows, honestly:

* The result-level aggregate stays within 0.025 (absolute) of the weighted truth in all eight settings; that residual is
  the sites' own Freyja error carried through, as it should be. With equal weights it is unaffected by depth (0.006-0.023).
* The pooled fit was 0.02-0.09 away from the depth-weighted truth and 0.05-0.20 away from the equal-weight truth when depths
  differ: pooling cannot be told to weight sites equally (or by load), it weights by depth. This is the larger, systematic
  effect.
* **In all eight pooled fits Freyja reported a lineage that no site contains** (`B.1`, 0.7-1.9%); the aggregate never does
  (asserted in the tests over random inputs). This is the concern that motivated the feature, and it is confirmed here on
  synthetic data. The cause was not investigated further (the likely reason is a sparse non-negative fit of a mixture of
  mixtures, but that is a hypothesis).
* Limits of the evidence: synthetic data with exactly known, linear lineage signatures and no real wastewater noise; four
  sites; one barcode set. It shows direction and rough size, not what real wastewater would give. The two approaches
  agree more closely than the truth-based columns suggest when depths are equal; the difference then is mostly small.
  The experiment driver is scratch work and is not committed.

## 2. The sites file (for non-bioinformaticians)

A plain table, one row per site per group. Keep it in Excel and **save as CSV**.

```
aggregate_name,site_name,weight
Metro_Sites,Site1,250000
Metro_Sites,Site2,120000
Metro_Sites,Site3,80000
```

| column | required | meaning |
|---|---|---|
| `aggregate_name` | yes | the name of the group; becomes part of the output row label (`Metro_Sites_2026-W41`). Use letters, digits, `_` and `-`; no tabs. |
| `site_name` | yes | a site as written in the `site` column of your samplesheet. |
| `weight` | no | a positive number, usually the population served (`250000`); plain digits, no units and no thousands separators. Only used when the weights basis says so (section 3). |

Rules and what the program does:

* A site can be in several groups: repeat it on another row (`Metro_Sites,Site1` and `North,Site1`).
* Site names are matched to the samplesheet `site` column **ignoring upper/lower case and extra spaces** (`site 1`, ` SITE 1 `
  and `Site 1` are the same site). `Site1` and `Site 1` are **not** the same.
* A site name that no sample uses is reported once on stderr as a warning (probably a typo); the group still works with the
  remaining sites.
* A group name written with different capitals on different rows (`Metro`, `metro`) is treated as one group (a warning names the
  line).
* A UTF-8 byte-order mark, Windows line endings, empty rows at the end, a semicolon- or tab-separated file (Excel in many
  countries writes `;`), header names such as `Aggregate Name` and an extra empty column are all accepted.
* Errors name the line and stop the run, for example:
  `sites.csv: line 7: site_name is blank (aggregate 'Metro_Sites')`,
  `sites.csv: line 5: site 'site1' is listed twice for aggregate 'Metro_Sites' (first on line 3)`,
  `sites.csv: line 1: unknown column 'site_names' (did you mean 'site_name'?)`,
  `sites.csv: line 4: weight '1 200' must be a positive number (plain digits such as 125000, no units or thousands separators) or blank`,
  and for an empty file or a file with a header only.

Common mistakes:

1. **Spelling** (`Site_1` in the samplesheet, `Site1` in the sites file): the warning "site 'Site1' in the sites file matches no
   sample's 'site'" is the symptom.
2. **Saving as the wrong CSV.** Excel offers "CSV (Comma delimited)" (Windows-1252; accented letters break) and "CSV UTF-8". Use
   **CSV UTF-8**. A non-UTF-8 file is still read (as Windows-1252) with a warning.
3. Typing a unit or separator in `weight` (`250 000`, `250,000`, `1.2M`, `250000 people`).
4. Leaving the header out, or renaming a column.
5. A site's samples with **no or a wrongly formatted `collection_date`** (use `2026-10-12`, not `12/10/2026`): they cannot be put
   in any week and are counted on stderr (`no_date=N`).
6. Expecting negative controls or failed/low-coverage samples to count: they never do.

A second group in the same file (commented here because CSV has no comments; add the rows to use it):

```
aggregate_name,site_name,weight
Metro_Sites,Site1,250000
Metro_Sites,Site2,120000
Metro_Sites,Site3,80000
Coastal,Site3,80000
Coastal,Site4,45000
```

Early check of the file alone (no samplesheet needed; with `--samplesheet` it also warns about unknown sites):
`site_aggregate.py --check-sites --sites sites.csv [--samplesheet samplesheet.csv]`.

## 3. Algorithm (exact)

Inputs: `freyja_aggregated.tsv` (row label = cleaned sample ID), `demix_status.tsv`, the samplesheet in effect
(`sample, sample_type, collection_date, site, viral_load`; other columns ignored, names matched case-insensitively), the sites file.

**Which samples contribute.** A sample contributes if it is in the aggregate table, its `demix_status` is `ok`,
`sample_type` is not `negative`, it has a `site` and a valid ISO `collection_date` (`YYYY-MM-DD`; a time part is dropped).
Samples without a site or a date cannot contribute and are counted on stderr (`no_site=`, `no_date=`); negative controls
(`negative_controls=`), non-ok samples (`not_ok=`) and aggregate rows without a samplesheet row are counted too. A
`viral_load` that is blank, not a number, or <= 0 means unknown.

**Period** (`--period`): `week` (default) is the ISO week, Monday to Sunday, labelled `2026-W41` with `period_start` = that
Monday (`2025-12-29` to `2026-01-04` is `2026-W01`; `2027-01-01` is `2026-W53`); `day` is the date itself (`2026-10-12`).

For each aggregate and each period in which at least one member site has a contributing sample:

1. **Replicates first.** For each member site with >= 1 contributing sample in the period, its samples are combined into one
   site solution. Weights are the samples' `viral_load` if **every** one of them has a load, else equal. A site with two
   samples therefore never counts double. The site's load is the **mean** of its samples' loads (only if all have one,
   else unknown); a mean, not a sum, so that the number of replicates does not change the site's weight.
2. **Sites.** `--weights` chooses the basis (written to `weight_basis`):
   * `viral_load`: the site's mean load. `site_weight`: the `weight` of the sites file. `equal`: 1 for each site.
   * `auto` (default): `viral_load` if every contributing site has a load; else `site_weight` if every contributing site has a
     weight; else `equal`. In `auto` the label says when a better basis was only partly available: `viral_load`, `site_weight`,
     `site_weight (viral_load missing)`, `equal`, `equal (viral_load missing)`, `equal (site_weight missing)`.
   * An explicit choice that cannot be met does not fail: it falls back to `equal` and says why
     (`equal (viral_load missing)`, `equal (site_weight missing)`).
   * The basis is decided per row, so it can change between periods (a missing load in one week only affects that week).
3. **Average.** For every lineage, `abundance = sum_i(w_i * a_i) / sum_i(w_i)` over the contributing sites (a lineage a site lacks
   counts as 0). The same weights apply to Freyja's `summarized` groups, `resid` and `coverage`. Lineages are matched by name, so
   an aggregate can never contain a lineage that no member site detected (asserted in the tests over random inputs, also for
   `summarized`). Nothing is renormalised (the weighted sum of the inputs' sums is kept); values are written with 8 decimals,
   lineages sorted by abundance (descending, then name).
4. **`--min-sites`** (default 1): an aggregate-period with fewer contributing sites is omitted (the number omitted is on stderr).

**Why not sequencing depth.** Weighting by read depth is deliberately **not** offered. Depth reflects how much was sequenced
(library size, run balancing, lane luck), not how much virus a site contributes; weighting by it would bias the group toward
well-sequenced sites, and it is exactly the implicit weighting that read pooling cannot avoid. The right weights describe the
population being summarised: **viral load** (the concentration measured for that sample, for example gene copies per litre;
it follows the day's actual shedding) or **population served** (`weight`; stable, known in advance, but blind to the day's
signal). If you want a group value that is proportional to the total virus entering all the sites, supply `viral_load` as a
load (concentration x flow, copies per day); a plain concentration is a reasonable but different weighting. Equal weights are
the honest fallback ("average site").

**Rows record their own completeness** (meta table): `n_sites_included`, `sites_included`, `sites_missing` and a parallel
`missing_reasons` (`no_sample`: the site has no sample in that period; `not_ok:<status>`: it has a dated sample but its demix
status is not `ok` (`failed`, `skipped_low_coverage`, `missing`, `no_status`, `not_in_aggregate`); `no_date`: it has no
dated sample in that period but has samples without a valid date), `weight_basis`, `n_samples` and
`weight_fraction_included`: included weight / total weight of the member sites, using the chosen basis (`viral_load`: a missing
site's load is the mean load of its non-ok samples in that period, `NA` unless every missing site has one; `site_weight`: `NA`
unless every missing site has a weight; `equal`: included sites / member sites; 1 when nothing is missing).

### Hand-computed week (`tests/data_timeseries`, aggregate `Metro` = SiteA + SiteB + SiteC)

ISO week `2026-W03` (Monday 2026-01-12). XBB.1.5 and loads of the six samples:

| site | sample XBB.1.5 (load) | replicates by load | site load (mean) |
|---|---|---|---|
| SiteA | 0.16710542 (851825), 0.23093447 (301960) | 0.183810 | 576892.5 |
| SiteB | 0.08106117 (1169960), 0.11608544 (938723) | 0.096653 | 1054341.5 |
| SiteC | 0.05187679 (110143), 0.04421386 (964206) | 0.044999 | 537174.5 |

`Metro_2026-W03` XBB.1.5 = (0.183810 * 576892.5 + 0.096653 * 1054341.5 + 0.044999 * 537174.5) / 2168408.5 = **0.107045**
(basis `viral_load`, 6 samples; equal site weights would give 0.108488). In `2026-W02` SiteC's only sample has no viral load, so
that row's basis is `site_weight (viral_load missing)`. `Pair` (SiteA + `siteb `, matched ignoring case and spaces)
shows the matching; the negative control `NEG_1` (SiteA, 2026-01-12) and the skipped low-coverage `LOWCOV_1` never contribute.

## 4. Outputs (contract with the lineage-group work; do not change without the planner)

`site_aggregates_freyja.tsv`: Freyja aggregate format, tab-separated, the first column unnamed (the row label), columns
`summarized, lineages, abundances, resid, coverage` (no `pathogen`), list strings formatted like the real file
(`[('Omicron', 0.8383), ...]`, `BA.2 JN.1`, `0.83213495 0.11711788`), row label `<aggregate_name>_<period>`. Header always written.
Readable by `timeseries.py` and `build_report_data.py` parsers (checked).

`site_aggregates_meta.tsv`: `row_label, aggregate_name, period, period_start, n_sites_included, sites_included, sites_missing,
missing_reasons, weight_basis, weight_fraction_included, n_samples`; `;`-joined lists in the file's spelling (an empty list is an
empty cell), `NA` for an unknown `weight_fraction_included`. Header always written. Both files are sorted by aggregate (order of
first appearance in the sites file), then period.

```
row_label	aggregate_name	period	period_start	n_sites_included	sites_included	sites_missing	missing_reasons	weight_basis	weight_fraction_included	n_samples
Metro_2026-W03	Metro	2026-W03	2026-01-12	3	SiteA;SiteB;SiteC			viral_load	1.00000000	6
Metro_2026-W04	Metro	2026-W04	2026-01-19	2	SiteA;SiteC	SiteB	not_ok:failed	viral_load	NA	2
```

CLI: `site_aggregate.py --aggregated freyja_aggregated.tsv --demix-status demix_status.tsv --samplesheet samplesheet.csv --sites sites.csv
[--period week|day] [--weights auto|viral_load|equal|site_weight] [--min-sites 1] --out site_aggregates_freyja.tsv --meta-out site_aggregates_meta.tsv`.
Exit 1 with `site_aggregate.py: error: <file>: line N: ...` for input problems; never a traceback. Empty inputs write
header-only outputs (exit 0). stderr summary: `site_aggregate.py: 2 aggregate(s), 16 aggregate-period row(s) written (week); 30 of 33
aggregated samples contributed; not contributing: negative_controls=1, not_ok=1, no_site=0, no_date=1, ...`.

## 5. Module and wiring

```nextflow
include { SITE_AGGREGATE } from '../modules/local/site_aggregate/main'
```

| process | inputs | outputs (emit) |
|---|---|---|
| `SITE_AGGREGATE` | `path aggregated` (`freyja_aggregated.tsv`); `path demix_status` (staged as `demix_status.tsv`); `path samplesheet` (staged as `samplesheet.csv`); `path sites` (staged as `sites.csv`) | `aggregated` = `site_aggregates_freyja.tsv`; `meta` = `site_aggregates_meta.tsv`; `versions_python` (`topic: versions`) |

All four inputs are single-item (value) channels. `ext.args` takes `--period`, `--weights`, `--min-sites`.

Placement in `workflows/wasteline.nf` (`full`, `demix` and `aggregate` modes), after `FREYJA_AGGREGATE_PLOT`, next to the time-series block.
It needs the **summary** samplesheet in effect (with `site`, `collection_date` and, optionally, `viral_load`) and the
`demix_status.tsv` channel that is published to `summary/`:

```nextflow
if (params.site_aggregation) {
    SITE_AGGREGATE(
        FREYJA_AGGREGATE_PLOT.out.aggregated,
        ch_demix_status,                         // the collectFile channel behind summary/demix_status.tsv (one file)
        ch_samplesheet.first(),
        channel.value(file(params.site_aggregation, checkIfExists: true))
    )
    ch_site_aggregates_freyja = SITE_AGGREGATE.out.aggregated
    ch_site_aggregates_meta   = SITE_AGGREGATE.out.meta
}
```

`conf/modules.config`:

```groovy
withName: 'SITE_AGGREGATE' {
    ext.args   = { "--period ${params.site_aggregation_period} --weights ${params.site_aggregation_weights} --min-sites ${params.site_aggregation_min_sites}" }
    publishDir = [ path: { "${params.outdir}/summary" }, mode: params.publish_dir_mode,
                   pattern: 'site_aggregates_{freyja,meta}.tsv' ]          // summary/site_aggregates_freyja.tsv, summary/site_aggregates_meta.tsv
}
```

`conf/containers.lock.config` needs a `SITE_AGGREGATE` entry (python digest; `bin/lock_containers.py`, needs the registry).
A bad sites file fails this process after the long steps; because everything before it is cached, `-resume` after fixing the file
only reruns this step. For an early failure, call `site_aggregate.py --check-sites --sites <file>` in the parameter validation
(or reproduce its few rules there); the validator itself is `read_sites()`.

Params (`nextflow.config` + `nextflow_schema.json`; the report/summary group):

| param | default | meaning |
|---|---|---|
| `site_aggregation` | `null` | path to the sites file (CSV); **giving it enables the feature** |
| `site_aggregation_period` | `'week'` | `week` (ISO week, Monday-Sunday) or `day` |
| `site_aggregation_weights` | `'auto'` | `auto`, `viral_load`, `equal` or `site_weight` |
| `site_aggregation_min_sites` | `1` | omit aggregate-periods with fewer contributing sites |

```groovy
site_aggregation            = null     // sites file (CSV: aggregate_name,site_name[,weight]); enables result-level site aggregation
site_aggregation_period     = 'week'
site_aggregation_weights    = 'auto'
site_aggregation_min_sites  = 1
```

```json
"site_aggregation": { "type": "string", "format": "file-path", "exists": true, "mimetype": "text/csv", "pattern": "^\\S+\\.csv$",
  "description": "Sites file (CSV) that groups sites into aggregates; enables `summary/site_aggregates_*.tsv`.",
  "help_text": "Columns `aggregate_name`, `site_name` and optional `weight`. Per-site Freyja results are combined by a weighted average; reads are never pooled. See docs/usage.md." },
"site_aggregation_period": { "type": "string", "default": "week", "enum": ["week", "day"] },
"site_aggregation_weights": { "type": "string", "default": "auto", "enum": ["auto", "viral_load", "equal", "site_weight"] },
"site_aggregation_min_sites": { "type": "integer", "default": 1, "minimum": 1 }
```

The `^\S+\.csv$` pattern rejects paths with spaces; relax it if that matters.

Published files (add to `docs/output.md`, `summary/` table, modes `full`/`demix`/`aggregate`, only with `--site_aggregation`):
`summary/site_aggregates_freyja.tsv`, `summary/site_aggregates_meta.tsv`.

## 6. Interaction with the lineage-group summary and the report

* The parallel lineage-group work (`LINEAGE_SUMMARY`) consumes both files as its optional site-aggregate inputs
  (`ch_site_aggregates_freyja`, `ch_site_aggregates_meta`) and writes `lineage_summary_site_aggregates.csv`. It reads the
  aggregate in exactly the per-sample Freyja format, so each aggregate-period is treated like a sample; the meta table
  supplies `aggregate_name`, `period` and the completeness columns it needs.
* **Report** (`build_report_data.py` / template; not implemented here): a "Site aggregates" section, shown only when the
  meta table exists. For every aggregate-period one stacked bar built like the per-sample bars from
  `site_aggregates_freyja.tsv` (Summarized groups / All lineages toggle), x axis = `period_start` ordered by time, one panel (or
  selector) per aggregate. Hover/table text from the meta table: `n_sites_included` of the member sites, `sites_missing` with
  `missing_reasons`, `weight_basis`, `weight_fraction_included`, `n_samples`. An **INFO** flag (a run-level flag, never an
  error) when any aggregate-period is incomplete (`sites_missing` non-empty) or its `weight_basis` starts with `equal (` or
  `site_weight (` (a weaker basis than requested was used). `report_data.json` can carry `site_aggregates: {enabled, period,
  rows: [meta + {summarized, lineages}]}`.

## 7. Tests

* `pixi run --frozen pytest tests/python/test_site_aggregate.py` (79 tests): replicates (load vs equal), `viral_load` / `site_weight` /
  `equal` / `auto` incl. soft fallbacks and per-period basis, missing viral load, missing sites and every reason code,
  Sunday/Monday and year-boundary weeks (incl. `2026-W53`), `day` period, `--min-sites`, several aggregates sharing a site,
  case/space matching, negatives and non-ok exclusion, `summarized`/`resid`/`coverage` weighting, no renormalisation, output
  format, numpy/repr input forms, the sites-file quirks and every friendly error, empty inputs, `--check-sites`, the hand-computed
  fixture week, and the no-new-lineage / convex-combination property over random inputs (12 seeds x 4 weight modes).
* `NFT_WORKDIR=<scratch> pixi run --frozen nf-test test modules/local/site_aggregate/tests/main.nf.test --profile docker`
  (stub; fixture week by hand with weights and the pair aggregate; `ext.args --weights equal --min-sites 3`; a bad sites file fails
  naming the line). `bin/*` is an nf-test trigger, so adding `bin/site_aggregate.py` triggers a full run in CI.
