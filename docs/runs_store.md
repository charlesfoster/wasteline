# A store of past runs, with trend and emergence alerts

Every wasteline run answers "what is in this batch?". Wastewater surveillance asks a different question: "what is **changing**
across batches?". `bin/runs_store.py` keeps the results of many runs in one SQLite file and applies a few simple rules to
flag lineages, mutations and watchlist entries worth a second look.

**The alerts are screening hints, not epidemiological conclusions.** They come from a handful of noisy samples, depend on genome
coverage and on the barcode version used, and a lineage call can change with the next dataset. Look at the samples, the report
and the bootstrap intervals before acting on one.

It is a host-side helper: it uses the Python standard library only (`sqlite3`), runs on your host Python (3.9 or newer), is not
part of the pipeline and never touches the results it reads or the network.

## Quick start

```bash
python bin/runs_store.py add    ~/wasteline/store.sqlite /data/results/run1 /data/results/run2   # idempotent
python bin/runs_store.py runs   ~/wasteline/store.sqlite
python bin/runs_store.py alerts ~/wasteline/store.sqlite --out alerts.tsv --markdown alerts.md
```

```
  run                      date        mode  barcodes          samples  ok  low_conf  other  dated
  run1                     2026-10-16  full  2026-10-06T05-25  10       8   0         2      0
  2026-10-19               2026-10-19  full  2026-10-06T05-25  12       5   0         7      0

2 run(s), 22 sample row(s), 22 distinct sample(s) (latest run each)
```

`alerts` with fewer runs than a rule needs, or no alert, prints `no alerts`; it never fails because the store is small.

## `add`: ingest results

`add <store.sqlite> <results_dir>...` creates the store if needed and reads, for each results directory (or each results
directory one level below a folder you give), these files:

| File                                                | Goes into                                                                                          |
| --------------------------------------------------- | -------------------------------------------------------------------------------------------------- |
| `pipeline_info/run_info.json`                       | `runs`: run identity (`run_name@start`), mode, start date, Freyja and barcode version              |
| `summary/demix_status.tsv` (required)               | `samples`: status, reason, genome coverage, barcode version                                        |
| `summary/samplesheet.csv`, `summary/qc_summary.tsv` | `samples`: site, collection date, sample type (samplesheet first), QC status                       |
| `summary/freyja_aggregated.tsv`                     | `lineage_abundance` (with `lo`/`hi` from `freyja_boot_summary.tsv` when present)                   |
| `summary/lineage_summary.csv`                       | `group_abundance` (a blank cell is stored as 0; with `lo`/`hi` from `lineage_group_intervals.tsv`) |
| `summary/unexplained_mutations.tsv`                 | `unexplained` (including `in_dataset`)                                                             |
| `summary/mutations_of_interest.tsv`                 | `moi`                                                                                              |

Files that a run did not produce (no bootstrap, no lineage groups, no mutations of interest) are skipped. A directory without
`summary/demix_status.tsv` is skipped with a warning. Sample IDs are stored exactly as they are.

- **Idempotent.** Adding a directory again replaces its rows. A run is recognised by its identity (`run_name@start`) or by its
  path, so re-running `add` over a whole results folder after every batch is safe.
- **The same sample in several runs** (for example a re-analysis with newer barcodes, [reanalysis.md](reanalysis.md)): all rows
  are kept, each tagged with the run's barcode version. Queries use the **latest run per sample**: the newest barcode version,
  then the newest run start. The latest run is chosen first and only then filtered by status, so a later failed result is not
  silently replaced by an older one. A sample is identified by its ID together with its collection date.
- A directory without `run_info.json` is accepted (its path is the identity), but has no run date, so its samples count only
  when the samplesheet gives a collection date.

- **Panel results.** The results of a `--panel` run (`summary/<target>/demix_status.tsv`, `pipeline_info/<target>/run_info.json`)
  are added as **one run per target** (identity `run_name@start#target`; idempotent per directory and target). The `runs` and
  `samples` tables have a `target` column, empty for single-target results; a store made before this column existed is
  migrated by `ALTER TABLE ... ADD COLUMN target TEXT DEFAULT ''` the first time it is opened. Targets are never mixed: the latest
  run per sample is chosen per target, and the four alert rules group by target and site. `runs` shows a `target` column and
  `alerts.tsv` (and the Markdown and terminal output) a `target` column after `rule`, **only when some run in the store has a
  target**, so the output for single-target stores is unchanged.

The tables are plain SQLite: `runs`, `samples`, `lineage_abundance`, `group_abundance`, `unexplained`, `moi`. Query them with
`sqlite3` if you want something the alerts do not give.

## `runs`

`runs <store.sqlite>` lists the runs with their date, mode, barcode version and sample counts (`ok`, `low_confidence`, other
statuses, and how many have a collection date).

## `alerts`

```
alerts <store.sqlite> [--watchlist mutations.csv] [--since YYYY-MM-DD] [--out alerts.tsv] [--markdown alerts.md]
                      [--include-low-confidence] [--levels group,lineage] [thresholds]
```

**Which samples count.** Status `ok` or `low_confidence`, `sample_type` other than `negative`, and a date: the collection date
if the samplesheet has one, else the date the run started (the Markdown report says how many samples used the run date; a run date
is only a rough stand-in for the sampling date). Periods are ISO weeks. `low_confidence` samples
(`--min_genome_coverage_for_demix`, [PLAN D18](PLAN.md)) are **excluded from rules 1 and 2** unless you pass
`--include-low-confidence`; rules 3 and 4 include them and say so in `detail`. Sites are kept apart in every rule; samples with
no site form the site `(none)`. The catch-all group `other` is never alerted on.

| Rule                    | Fires when                                                                                                                                                 | Options (default)                                                                                                  |
| ----------------------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------ |
| `emerging`              | At a site, the mean abundance of a lineage group (and of a lineage) in the **latest** ISO week is >= X while **no earlier sample** of the site had >= 1 %. | `--emerging-min-pct` (5), `--emerging-max-earlier-pct` (1), `--min-history-periods` (2 earlier weeks needed)       |
| `sustained_growth`      | At a site, the abundance (averaged per date) is still rising over the last >= 3 consecutive dates, and the total rise is >= 10 percentage points.          | `--growth-min-samples` (3), `--growth-min-rise-pts` (10)                                                           |
| `recurring_unexplained` | The same unexplained mutation (`unexplained_mutations.tsv`) is in >= 3 samples on >= 2 distinct dates of a site within the last N weeks.                   | `--unexplained-min-samples` (3), `--unexplained-weeks` (8; the window ends at the latest sample date in the store) |
| `watchlist`             | A mutation of your `--watchlist` CSV is at >= 5 % frequency in a sample.                                                                                   | `--watchlist-min-pct` (5)                                                                                          |

Notes on the rules:

- Rules 1 and 2 work on both the **group** level (`lineage_summary.csv`) and the **lineage** level (`freyja_aggregated.tsv`);
  `--levels group` keeps only the groups. Lineage-level alerts are noisier: Freyja often shuffles abundance between sister
  lineages, and overlapping groups ("includes JN.1*" and "excludes JN.1*") can alert together.
- `emerging` needs `--min-history-periods` earlier weeks at the site, so a store with a week or two of history gives no
  emerging alerts. "Earlier" is judged on every single earlier sample, not the weekly mean.
- `recurring_unexplained` reports `in_dataset` (see [output.md](output.md#summaryunexplained_mutationstsv)): `in_dataset=yes`
  means some lineage of the barcode dataset carries the mutation but the fitted mixture does not predict it; `in_dataset=no`
  means **no** lineage of the dataset carries it: a blind spot of the dataset or a genuinely new change (the table cannot tell
  which). Read the second kind first.
- `watchlist` takes the same CSV as `--mutations_of_interest` (`name,mutation`; [usage.md](usage.md#10-mutations-of-interest)).
  It matches the `mutation` text exactly against `summary/mutations_of_interest.tsv` of the stored runs, and also against the
  unexplained table (nucleotide `C6501A` and amino-acid `ORF1ab:P2079Q` forms). A watchlist mutation that no stored sample
  measured is listed in a **note** (absence of an alert then says nothing): run the pipeline with `--mutations_of_interest`
  and the same list.
- `--since YYYY-MM-DD` keeps only alerts whose latest evidence is on or after that date; the history behind them is still used.

**Output.** `--out` writes a TSV (a header only when there is no alert) with the columns `rule`, `level` (`group`, `lineage`
or `mutation`), `site`, `group_or_mutation`, `detail`, `first_date`, `last_date`, `samples` (comma-separated sample IDs, at most
eight shown). `--markdown` writes a short report for people with the same table per rule, the thresholds used and the notes.
The same alerts are printed to the terminal.

## Known limits

- Only as good as the dates: without collection dates the run date is used, and several batches run on one day make a single
  "week".
- No cross-site pooling: a mutation seen once at each of three sites does not fire rule 3.
- Sample identity is the sample ID plus its collection date; the same ID reused for a different sample on the same date is not
  distinguished.
- Abundances come from Freyja at the barcode version of each run. Mixing barcode versions in one store is allowed, but a lineage
  that appears only because the newer barcodes can name it looks like an emerging lineage; re-analyse the older runs
  ([reanalysis.md](reanalysis.md)) and `add` them again to compare like with like.
- A store is a convenience, not a record: the results directories stay the source of truth, and `add` can rebuild it.
