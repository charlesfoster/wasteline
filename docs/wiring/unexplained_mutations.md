# UNEXPLAINED_MUTATIONS: wiring notes

"Unexplained mutations" are high-frequency SNVs in a sample that the demixed lineage mixture does not account for: an
early signal for a novel or cryptic lineage (or a lineage missing from the barcodes). It is deliberately **not a variant
caller**. Script: `bin/unexplained_mutations.py` (stdlib only; it imports the GFF/codon translator from
`bin/mutations_of_interest.py`, which sits in the same `bin/` directory, so both must stay together). Modules:
`modules/local/unexplained_mutations/` and `modules/local/unexplained_mutations_concat/` (container
`quay.io/biocontainers/python:3.12.12`, the same tag as `MUTATIONS_OF_INTEREST`). Tests:
`tests/python/test_unexplained_mutations.py` (62), the two modules' `tests/main.nf.test` (5 + 5, stub + real docker),
fixtures in `tests/fixtures/unexplained/`. Nothing existing was edited; this file is the hand-over.

## Method (exact)

For each sample, per SNV `m` (position, ref, alt) seen in the sample's ivar variants TSV:

```
expected_freq(m) = min(1, sum over demixed lineages l of  a_l * barcode(l, m))
excess(m)        = observed_freq(m) - expected_freq(m)
```

- `a_l`: abundance from the `lineages` / `abundances` rows of `<s>.demix` (never `summarized`). Both rows are
  whitespace-separated in the real file (`lineages  BA.2 AY.4`, `abundances  0.71212100 0.28787900`); a Python/numpy repr
  (`['BA.2', 'AY.4']`, `[0.71 0.29]`) is accepted too.
- `barcode(l, m)`: 1 if the nucleotide SNV (`C241T`) is in lineage `l`'s barcode (a fractional value in a custom CSV is
  used as given, clipped to [0, 1]). Source, in order: dataset directory with `usher_barcodes.feather` -> its
  `lineage_mutations.json` (verified identical to the feather for JN.1, BA.2, AY.4 in Freyja 2.0.6; same 5,977 lineages);
  directory with only `usher_barcodes.csv` -> that CSV; directory with only a JSON -> the JSON; a `.csv` file (custom
  barcodes, e.g. `tests/data/barcodes_test.csv`) or `.json` file as given. `.feather` only / a bare `.feather` -> exit 2
  with an actionable message (so `--barcodes x.feather` makes the process fail: use `errorStrategy 'ignore'`, below).
  The CSV is streamed and only the demixed lineages' rows are parsed.
- Observed: ivar `ALT_FREQ` (`ALT_DP / TOTAL_DP` if it is not numeric), `TOTAL_DP`, `ALT_DP`, `ALT_QUAL`. The `PASS`
  column is ignored (as in `MUTATIONS_OF_INTEREST`).
- v1 considers SNVs only. Indel rows in the TSV (`ALT` starting `+`/`-`) and non-SNV barcode tokens are skipped and
  counted in one stderr note per sample.

A mutation is reported when ALL of (all inclusive; a 1e-9 tolerance so values exactly on a threshold are kept):

| rule                                                                               | option                                        | default                      |
| ---------------------------------------------------------------------------------- | --------------------------------------------- | ---------------------------- |
| `TOTAL_DP >=`                                                                      | `--min-depth`                                 | 20                           |
| `ALT_DP >=`                                                                        | `--min-alt-depth`                             | 10                           |
| `ALT_QUAL >=` (if numeric)                                                         | `--min-alt-qual`                              | 20                           |
| observed `>=`                                                                      | `--min-freq`                                  | 0.05                         |
| observed - expected `>=`                                                           | `--min-excess`                                | 0.05                         |
| `P(X >= ALT_DP) <=` for `X ~ Binomial(TOTAL_DP, max(expected, --background-rate))` | `--max-pvalue` (1 = off), `--background-rate` | 0.05 / 29903 = 1.67e-6, 0.01 |

**The binomial guard is an addition to the requested rules** and is on by default. Reason (measured, see below): at the
20-130x depths of the synthetic data a perfectly explained mixture still produces many positions whose observed
frequency differs from the expectation by more than 0.05 by chance (MIX_A: 18 rows, MIX_B: 15, Site_1_week_2: 9 without
the guard; 0 with it). The default is a **Bonferroni-style** level, `0.05 / 29903` (= 1.67e-6; in the code
`DEFAULT_MAX_PVALUE = 0.05 / GENOME_LENGTH`): about 30,000 positions are tested per sample, so a per-site 0.001 would
let ~30 chance calls through and did let the Site_1_week_2 chance call (p = 3.3e-4) pass. It is a fixed constant of a
genome-sized test count, not derived from the pathogen; override it with `unexplained_max_pvalue`. The three spiked SNVs
(30 %, 18 %, 12 %) pass with a very large margin (p = 2.4e-42, 1.4e-17, 2.3e-14, i.e. 8 to 36 orders of magnitude below
the cut-off), so no compromise default was needed. At real wastewater depths (hundreds to thousands of x) the guard
is irrelevant for genuine signals. Set `unexplained_max_pvalue = 1` to get exactly the requested plain rules. The
p-value of every reported row is the last output column (`p_value`, `%.3g`, e.g. `2.39e-42`), also when the guard is
off, so borderline calls can be judged by eye.

## Outputs

Per sample `<sample>.unexplained_mutations.tsv`:

```
sample  position  ref  alt  observed_freq  expected_freq  excess  depth  alt_depth  gene  aa_change  effect  p_value
```

(frequencies with 6 significant digits; `depth` = TOTAL_DP, `alt_depth` = ALT_DP; `p_value` = the binomial tail
probability of the Method table, 3 significant digits.)

Run level `unexplained_mutations.tsv` (concat) = the same, with two columns inserted before `p_value` (which stays last)

```
...  effect  n_samples  recurrent  p_value
```

`n_samples` = number of distinct samples in the input having the same position + alt unexplained; `recurrent` =
`true` when `n_samples >= --min-recurrence` (default 2). Sorted recurrent first, then excess descending, then sample,
position. Header-only inputs are fine (zero rows -> header-only output); a table with another header is exit 2.
Negative controls and samples without a demix are simply not callers (header-only or absent).

## Include lines and signatures

```nextflow
include { UNEXPLAINED_MUTATIONS        } from '../modules/local/unexplained_mutations/main'
include { UNEXPLAINED_MUTATIONS_CONCAT } from '../modules/local/unexplained_mutations_concat/main'
```

| process                        | inputs                                                                                                                                                                                             | outputs (emit)                                                                           |
| ------------------------------ | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | ---------------------------------------------------------------------------------------- |
| `UNEXPLAINED_MUTATIONS`        | `tuple val(meta), path(variants), path(demix)` (demix may be `[]`); `path dataset` (Freyja dataset dir or custom barcodes `.csv`); `path gff`; `path reference` (`[]` for both to skip annotation) | `tsv` = `[ meta, <id>.unexplained_mutations.tsv ]`; `versions_python` (topic `versions`) |
| `UNEXPLAINED_MUTATIONS_CONCAT` | `path tables` (collected, staged as `unexplained/*`)                                                                                                                                               | `tsv` = `unexplained_mutations.tsv`; `versions_python` (topic `versions`)                |

Both honour `ext.args` (thresholds), `ext.prefix` (per-sample only) and `ext.when`. The three value-like inputs
(`dataset`, `gff`, `reference`) are single items / value channels.

## Workflow placement

In `workflows/wasteline.nf`, inside the `mode == 'full' || mode == 'demix'` block next to `MUTATIONS_OF_INTEREST`
(needs `ch_variants`, `ch_demixed` = `FREYJA_DEMIX.out.demix`, `FREYJA_DATA.out.dataset`, `ch_gff`, `ch_fasta`).
Aggregate and report modes have no variants file, so nothing runs there. Negative controls are excluded.

```nextflow
def ch_unexplained = channel.value([])
if (params.unexplained_mutations) {
    // every non-negative sample; a sample without a demix (skipped / failed) gets a header-only table
    def ch_unexplained_in = ch_variants
        .filter { meta, _variants -> meta.sample_type != 'negative' }
        .join(ch_demixed, remainder: true)
        .filter { meta, variants, _demix -> meta && variants }
        .map { meta, variants, demix -> [ meta, variants, demix ?: [] ] }
    UNEXPLAINED_MUTATIONS(ch_unexplained_in, FREYJA_DATA.out.dataset, ch_gff, ch_fasta)
    UNEXPLAINED_MUTATIONS_CONCAT(UNEXPLAINED_MUTATIONS.out.tsv.map { _meta, tsv -> tsv }.collect())
    ch_unexplained = UNEXPLAINED_MUTATIONS_CONCAT.out.tsv.ifEmpty([])
}
```

Notes: `ch_variants.join(ch_demixed, remainder: true)` can emit demix-only tuples (negatives are in `ch_demixed`) which
the `filter` drops; if the planner prefers, replace it by an inner join of the demixed samples only (skipped samples
then simply have no per-sample file; the concat is unaffected). `FREYJA_DATA.out.dataset` is the same resolved dataset
the demix used (`--barcodes latest|container|<date>|<csv>`), which is what makes the expectation consistent with the
abundances. `ch_gff`/`ch_fasta` are the same value channels `MUTATIONS_OF_INTEREST` gets. An empty `collect()` (every
sample negative) emits nothing, so the concat does not run: hence `.ifEmpty([])` and the report input `[]`.

## `conf/modules.config` blocks

```groovy
    // ---- unexplained mutations ------------------------------------------------------------
    withName: 'UNEXPLAINED_MUTATIONS' {
        ext.args = {
            [
                "--min-depth ${params.unexplained_min_depth}",
                "--min-alt-depth ${params.unexplained_min_alt_depth}",
                "--min-freq ${params.unexplained_min_freq}",
                "--min-excess ${params.unexplained_min_excess}",
                "--max-pvalue ${params.unexplained_max_pvalue}",
            ].join(' ')
        }
        errorStrategy = 'ignore'     // e.g. a .feather-only dataset: the run continues, the table is just absent
        publishDir = [
            path: { "${params.outdir}/samples/${meta.id}/variants" },
            mode: params.publish_dir_mode
        ]
    }

    withName: 'UNEXPLAINED_MUTATIONS_CONCAT' {
        ext.args = { "--min-recurrence ${params.unexplained_min_recurrence}" }
        errorStrategy = 'ignore'
        publishDir = [
            path: { "${params.outdir}/summary" },
            mode: params.publish_dir_mode
        ]
    }
```

(Alternatively add `UNEXPLAINED_MUTATIONS_CONCAT` to the existing `summary/` group
`'CONSENSUS_STATS_CONCAT|CONSENSUS_CONCAT|MUTATIONS_OF_INTEREST_CONCAT'` and keep only `ext.args` in its own block.)

## Params (suggested)

| param                        | default                    | meaning                                                                                            |
| ---------------------------- | -------------------------- | -------------------------------------------------------------------------------------------------- |
| `unexplained_mutations`      | `true`                     | run the analysis (only meaningful in `full` / `demix`); `false` skips it                           |
| `unexplained_min_freq`       | `0.05`                     | minimum observed frequency (0-1)                                                                   |
| `unexplained_min_excess`     | `0.05`                     | minimum observed minus expected frequency (0-1)                                                    |
| `unexplained_min_depth`      | `20`                       | minimum total depth at the position                                                                |
| `unexplained_min_alt_depth`  | `10`                       | minimum alt-supporting reads                                                                       |
| `unexplained_min_recurrence` | `2`                        | samples with the same mutation unexplained for `recurrent = true`                                  |
| `unexplained_max_pvalue`     | `1.67e-6` (= 0.05 / 29903) | binomial guard against depth noise, Bonferroni-style over the genome; `1` disables it (see Method) |

`nextflow_schema.json`: a group "Unexplained mutations" (booleans/numbers with min/max: frequencies 0-1, depths
integer >= 0, `unexplained_min_recurrence` integer >= 1, `unexplained_max_pvalue` 0-1). Declare the defaults in the
typed `params {}` block of `nextflow.config` (D12). Document in `docs/parameters.md`, and in `docs/output.md`:
`samples/<s>/variants/<s>.unexplained_mutations.tsv` and `summary/unexplained_mutations.tsv` (`full` / `demix` only),
and in `docs/qc_scripts.md` (the script follows the `mutations_of_interest.py` pattern).

## Report integration

1. **Section "Unexplained mutations"** (visible only when `unexplained_mutations.tsv` exists and has rows; with zero rows
   a one-line "No unexplained mutations" note). One table: sample, position, ref>alt, gene, aa_change, effect,
   observed %, expected %, excess, depth, alt depth, p-value, n_samples. Rows with `recurrent = true` get a highlight (and sort
   first, as in the file). A caption explains the method in one sentence ("high-frequency mutations that the demixed
   lineages do not predict; a hint of a novel lineage, not a variant call") and states the thresholds from
   `thresholds`.
2. **Run-level flag** `unexplained_mutations`, level `info`: raised when any row has `recurrent == true`; `samples` =
   the sorted IDs of the samples carrying a recurrent row; message e.g. _"N mutation(s) unexplained by the demixed
   lineages recur in several samples (first: S:E180D in 3 samples); possible novel lineage."_ Add the code to the
   `run_flags` list in `docs/report_data.md`.
3. **`report_data.json` fields** (all new and optional):
   - `available.unexplained_mutations`: boolean.
   - `unexplained_mutations`: `null` = not run; `[]` = run with no rows; else the concat rows as
     `{sample, position, ref, alt, observed_freq, expected_freq, excess, depth, alt_depth, gene, aa_change, effect,
n_samples, recurrent, p_value}` (numbers as numbers, `NA` -> `null`, `recurrent` -> boolean), file order.
   - `thresholds.unexplained_min_freq`, `unexplained_min_excess`, `unexplained_min_depth`, `unexplained_min_alt_depth`,
     `unexplained_min_recurrence`, `unexplained_max_pvalue` (from the run params). `p_value` is a float (`2.39e-42` -> 2.39e-42), or null.
   - optionally per sample in `samples[*]`: `n_unexplained` (integer, or null when the sample has no demix).
     `BUILD_REPORT_DATA` takes a new `--unexplained-mutations` option, fed through `wrapItem(...)` like the other
     optional tables; `report`-mode `optionalFile(results, 'unexplained_mutations.tsv')` picks it up from old results.
4. Optional: `qc_summary.tsv` is unchanged (no new column is needed).

## What the real data showed (synthetic `tests/data`, Freyja 2.0.6 e2e run `out_full`)

Dataset used: **`tests/data/barcodes_test.csv`** (custom CSV, 14 lineages), the one the samples were demixed with; and,
as the "messier" case, the **container-bundled 10_06_2026 dataset** (`lineage_mutations.json`, 5,977 lineages) after
re-running `freyja demix` of the same variants/depths with it in `quay.io/biocontainers/freyja:2.0.6--pyhdfd78af_0`
(MIX_A then gets 9 lineages, e.g. `BA.2 0.692, AY.4 0.176, AY.1 0.071, ...`; Site_1_week_2 gets 9; PURE_C stays
`JN.1 1.0`).

| sample        | custom CSV, default | bundled dataset, default | custom CSV, `--max-pvalue 1`            | custom CSV, `--max-pvalue 0.001` (old default)                                                                                                                   |
| ------------- | ------------------- | ------------------------ | --------------------------------------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| PURE_C        | 0 rows              | 0 rows                   | 0                                       | 0                                                                                                                                                                |
| MIX_A         | 0 rows              | 0 rows                   | 18 (all depth 35-86x, excess 0.05-0.15) | 0                                                                                                                                                                |
| MIX_B         | 0 rows              | 0 rows                   | 15                                      | 0                                                                                                                                                                |
| NEG_CONTAM    | 0 rows              | n/a                      | 0                                       | 0                                                                                                                                                                |
| Site_1_week_2 | 0 rows              | 0 rows                   | 9                                       | 1 row: `9866 C>T` ORF1ab L3201F, 17/22 reads (0.77) vs expected 0.39 (XBB.1.5 + a spurious 3 % BA.2 both carry it), p = 3.3e-4 (8.0e-4 with the bundled dataset) |

All five synthetic samples therefore give 0 rows with both datasets at the new default; every row that appears with
`--max-pvalue 1` has `p_value` above the cut-off (test `test_real_samples_without_the_guard_show_the_noise`). The old
0.001 default is the only setting here that produced a false call. With 20-50x data this signal is weak; it is meant for
deep wastewater data.

**Controlled positive test** (`tests/fixtures/unexplained/PURE_C_spiked.variants.tsv.gz` = PURE_C's real variants plus
three made-up SNVs at positions in no barcode of the 5,977-lineage dataset). Reported, and nothing else:

```
PURE_C  22102  A  T  0.3   0  0.3   120  36  S   E180D  missense     2.39e-42
PURE_C  28802  A  T  0.18  0  0.18  100  18  N   R177*  stop_gained  1.41e-17
PURE_C  29800  G  A  0.12  0  0.12  150  18  NA  NA     intergenic   2.28e-14
```

(hand check: S starts at 21563, `22102 - 21563 = 539`, codon 180 third base, GAA/GAG -> GAT = Asp; N starts at 28274,
`28802 - 28274 = 528`, codon 177 first base, AGA -> TGA = stop; 29800 is past ORF10 which ends at 29674.) JN.1's own
mutations (C241T, A23063T = S:N501Y, ...) are in the file at ~100 % and are not reported. Feeding three copies of that
file as samples A, B, C through the concat gives 9 rows, all `n_samples = 3`, `recurrent = true`, sorted with the three
`22102` rows first. When the spiked mutation is declared part of a hypothetical lineage in the mixture at 30 %, it
disappears (test `test_spike_explained_by_a_lineage_in_the_dataset_is_not_reported`).

## Tests

- `pixi run pytest tests/python/test_unexplained_mutations.py`: 62 passing. Hand-checked mini fixture (expected = sum of
  abundances: 0.8, 0.65, 0.5, 0.3 and 0.15 for a minor-lineage-only mutation), exact thresholds for min-freq, min-excess,
  min-depth, min-alt-depth, ALT_QUAL, indel skipping, missing-lineage tolerance (counted on stderr, also when all are
  missing), demix parsing (plain and repr forms, malformed), dataset resolution and errors, annotation (S D614G,
  ORF1ab P4715L, intergenic, start_lost, overlapping ORF7a/ORF7b), header-only for skipped samples, `--max-rows`, binomial
  tail values, the real-data checks above, concat (recurrence, order, `--min-recurrence`, zero rows, bad header).
- Module nf-tests (docker): `nf-test test modules/local/unexplained_mutations/tests/main.nf.test --profile docker` (stub
  - spiked PURE_C + MIX_A header-only + dataset-directory/no-GFF + no-demix) and
    `.../unexplained_mutations_concat/tests/main.nf.test` (stub + 3-sample recurrence with a skipped sample +
    header-only + `ext.args --min-recurrence 4` + wrong header rejected).
- Quality: `ruff check` and `ruff format --check` clean on the script and test file; `nextflow lint` on both modules
  has only the `when:` deprecation warning that every module in the repository has.
