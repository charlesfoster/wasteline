# Primer scheme check: hit-fraction floor (wiring notes)

Edited by this change: `bin/primer_scheme_check.py`, `modules/local/primer_scheme_check/` (main.nf stub, tests, snapshot),
`tests/python/test_primer_scheme_check.py`, `tests/fixtures/qc/QIAGEN_LIKE_EDEN.read_starts.tsv` (new, synthetic, 6000 reads).
Nothing else. `bin/qc_summary.py` needs **no change** (checked: `parse_primer_scheme_summary` reads columns by name, so
the two appended columns are ignored; run on the new 9-column summary it returns `('mismatch', 'qiagenv1', 6.8)`;
`tests/python/test_qc_summary.py` 76 passed).

## Why

First real QIAseq data (11 samples, declared `qiagenv1`): the right scheme is clearly enriched, but the sparse `eden`
scheme (28 sites, window coverage 0.23 %) also reached enrichment 6.0-7.8 (above the bar of 5, p at the 1e-300 floor)
with only 1.4-1.8 % of reads at its sites, because its tiny expected count inflates the ratio. Re-scoring the real read
starts of all 11 samples against every possible declared scheme (99 combinations) with the old rule gave **3 false `ok`**:
`eden` declared on WW-1287E, WW-1289E, WW-1290E (the best, qiagenv1, was only 2.2-2.7x eden's enrichment, so the old
`>= 3x` rule did not even call a mismatch there). The hit fraction (hits / reads) is an absolute signal that is not
inflated by a small window set.

## New rule (the output columns/formats of the old rule are kept; two columns are appended)

- New option `--min-hit-fraction` (default **0.05**, 0 disables = the old behaviour).
- A scheme **reaches the bar** iff `enrichment >= --min-enrichment` (5) AND `p_value <= --max-p-value` (0.001) AND
  `hit_fraction >= --min-hit-fraction` (0.05).
- `ok`: the declared scheme reaches the bar (`declared_enriched`). Unchanged.
- `mismatch` (`best_other_scheme`): declared does not reach the bar; the best scheme (highest enrichment among those that
  reach the bar, so its hit fraction is necessarily >= the floor) is not equivalent to the declared one AND
  (`best_enrichment >= 3 x declared_enrichment` OR `best_hit_fraction >= 3 x declared_hit_fraction`; an NA declared
  enrichment or hit fraction counts as satisfied). The second alternative handles a sparse declared scheme with an
  inflated ratio.
- `uninformative`: `no_scheme_enriched` (nothing reaches the bar), `best_not_decisive` (a different best exists but
  neither 3x test is met), and the **new reason `declared_sparse_inflated`**: nothing else reaches the bar (or the best
  is equivalent to the declared scheme) and the declared scheme passes the enrichment and p-value tests but not the
  hit-fraction floor, i.e. its ratio is inflated by a tiny window set. Not an `ok`, not a warning.
- `NA`: `too_few_reads`, `missing_input`, `declared_not_found`. Unchanged. All other reason strings are unchanged;
  `qc_summary.py` does not read `reason`.
- Summary TSV: the seven existing columns keep their positions; **`declared_hit_fraction`, `best_hit_fraction`** (4 dp,
  `NA` when not computable; `best_hit_fraction` is the top scheme's, `NA` when `best_scheme` is `NA`) are appended last.
  The per-scheme detail table is unchanged. The pipeline always writes the extended summary (no flag needed: the only
  consumer, `qc_summary.py`, is name-based).

## Why 0.05 (all numbers)

| data                                                                       | scheme(s)                         | enrichment       | hit_fraction                                      |
| -------------------------------------------------------------------------- | --------------------------------- | ---------------- | ------------------------------------------------- |
| real QIAseq (11 samples), TRUE scheme                                      | qiagenv1                          | 12.7-17.4        | **0.584-0.799**                                   |
| real QIAseq, wrong sparse scheme                                           | eden (28 sites)                   | 0.0-7.8          | **0.000-0.018**                                   |
| real QIAseq, wrong dense schemes                                           | articv3 / articv4.1 / articv5.3.2 | 0.5-3.4          | 0.009-0.060                                       |
| real QIAseq, wrong dense schemes                                           | swift = xgen                      | 1.0-2.0          | 0.056-0.116                                       |
| real QIAseq, wrong sparse schemes                                          | midnight / midnightv2             | 0.3-1.0          | 0.0015-0.005                                      |
| synthetic Midnight-shaped (PURE_C fixture; MIX_A, Site_1_week_2 from docs) | midnight                          | 41.7 (41.4-42.8) | **0.20** (about 0.20-0.21 = enrichment x 0.00485) |
| synthetic, Midnight reads, wrong schemes (PURE_C)                          | eden / articv3 / qiagenv1 / swift | 1.4-2.6          | 0.006-0.091                                       |
| random fragments (UNANCHORED fixture, 3000 reads)                          | all                               | 0.9-1.3          | 0.002-0.051                                       |
| new fixture QIAGEN_LIKE_EDEN (72 % at qiagenv1 ends + 48 reads on eden)    | eden / qiagenv1                   | 6.8 / 15.6       | 0.016 / 0.717                                     |

A dense wrong scheme (swift/xgen, articv3) can have a hit fraction above 0.05 but never an enrichment near 5 unless its
reads are genuinely anchored there; the floor exists for sparse schemes. The floor must lie above the largest wrong sparse
scheme (0.018) and below the lowest true scheme (0.20 for the synthetic Midnight data, 0.58 real QIAseq): any value in
**[0.019, 0.19]** gives the same 99-combination matrix on the real data (verified at 0.019, 0.05 and 0.19; 0.0181 still
leaves one false `ok`). 0.05 is the geometric-mean region (3x the largest wrong sparse value, 4x below the synthetic
Midnight data) and keeps every synthetic fixture `ok`.

Trade-off to document: for a very sparse scheme such as Midnight (window coverage 0.49 %), enrichment 5 alone needs only
about 2.4 % of reads at primer ends; with the floor, Midnight-like libraries need >= 5 % (enrichment about 10). A weakly
anchored Midnight library with 2.4-5 % anchored reads therefore becomes `uninformative` (`declared_sparse_inflated`)
instead of `ok`. Never a false warning: only `mismatch` warns.

## Validation: 11 real samples x 9 declared schemes = 99 combinations

(read starts regenerated for the run's pre-trim BAMs, in scratch only; the rescored qiagenv1 enrichments equal the published
ones to 0.1 for all 11 samples.)

| declared             | verdicts, new rule                                                            | old rule (no floor)                                            |
| -------------------- | ----------------------------------------------------------------------------- | -------------------------------------------------------------- |
| qiagenv1 (right)     | 11 / 11 `ok` (`declared_enriched`)                                            | 11 / 11 `ok`                                                   |
| 8 wrong schemes (88) | 88 / 88 `mismatch` (`best_other_scheme`, best `qiagenv1`, best enr 12.7-17.4) | 85 `mismatch`, **3 false `ok`** (eden: WW-1287E, 1289E, 1290E) |

The 3 eden cases are caught only by the new `best_hit_fraction >= 3x declared_hit_fraction` alternative (their enrichment
ratio is only 2.2-2.7x), 85 by either; no case is `uninformative`. Reproduction: `tests/python/test_primer_scheme_check.py`
(synthetic) and `tests/fixtures/qc/QIAGEN_LIKE_EDEN.read_starts.tsv` (a real-data-shaped stand-in; the real starts are not in
the repo).

## Wiring to apply (not done here)

1. `nextflow.config`, next to `primer_check_min_enrichment` (line ~78):
   `primer_check_min_hit_fraction = 0.05   // at least this fraction of reads must start at the declared scheme's primers (guards sparse schemes)`
2. `conf/modules.config` `PRIMER_SCHEME_CHECK` `ext.args`:
   `{ "--window ${params.primer_check_window} --min-enrichment ${params.primer_check_min_enrichment} --min-hit-fraction ${params.primer_check_min_hit_fraction}" }`
3. `nextflow_schema.json` (next to `primer_check_min_enrichment`):
   `"primer_check_min_hit_fraction": {"type": "number", "default": 0.05, "description": "Fraction of reads that must start within the window of a primer 5' end for a scheme to count as matching (0 disables). Stops sparse schemes whose enrichment ratio is inflated by a tiny window set.", "fa_icon": "fas fa-ruler", "minimum": 0, "maximum": 1}`
4. `docs/parameters.md`: add the row `--primer_check_min_hit_fraction | number (0-1) | 0.05 | ...` (same text).
5. `docs/qc_scripts.md` (primer_scheme_check section):
   - usage line: add `[--min-hit-fraction 0.05]`;
   - "A scheme reaches the bar" bullet: add `and hit_fraction >= --min-hit-fraction (0.05)`, with the reason (sparse schemes
     inflate the ratio);
   - summary-row sentence: columns `... reason, declared_hit_fraction, best_hit_fraction`;
   - `mismatch` row of the table: "... has >= 3x its enrichment or >= 3x its hit fraction"; `uninformative` row: add
     reason `declared_sparse_inflated`;
   - "Measured" section: add the real-QIAseq numbers table above (no read data is in the repo: only the summary numbers);
     replace the QIAseq limit bullet ("actual enrichment is unknown until a real QIAseq run") with: real QIAseq enrichment
     12.7-17.4 and hit fraction 0.58-0.80 for the declared scheme; sparse eden reached 6-7.8 at 1.4-1.8 %;
   - Limits: add the Midnight sensitivity trade-off above.
6. `docs/usage.md` / `docs/output.md`: wherever `scheme_check` reasons or `primer_check_*` params are listed, add the new
   param and note that a declared sparse scheme (e.g. `eden`) is no longer `ok` on ratio alone. `qc_summary.tsv` itself is
   unchanged (27 columns); the two new columns exist only in `coverage/<s>.primer_scheme_summary.tsv`.
7. Pipeline-level snapshots that hash `primer_scheme_summary.tsv` (if any) change with the two extra columns; I found none
   (`tests/default.nf.test*` only read the `qc_summary` columns).
