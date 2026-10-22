# QC helper scripts

Standalone, stdlib-only Python scripts in `bin/` (PLAN section 5 steps 9-12, section 6a, decision D7: coverage QC comes
from the Freyja depths file, not mosdepth). Each is unit-tested in `tests/python/` (fixtures in `tests/fixtures/qc/`).
Exit codes: 0 ok, 2 bad input/usage. Warnings go to stderr.

Inputs shared by the scripts:

- **depths**: Freyja `*.depths`, four tab-separated columns, no header: `contig, 1-based position, ref base, depth`
  (Q>=20 bases), one row per genome position.
- **variants**: Freyja/ivar variants TSV (`REGION POS REF ALT REF_DP REF_RV REF_QUAL ALT_DP ALT_RV ALT_QUAL ALT_FREQ
TOTAL_DP PVAL PASS GFF_FEATURE REF_CODON REF_AA ALT_CODON ALT_AA POS_AA`). Indels: `ALT` starts with `+` (insertion,
  anchored at `POS`) or `-` (deletion of the bases that follow `POS`, i.e. first deleted base is `POS + 1`). The AA columns
  are `NA` unless ivar ran with `--annot`.

## amplicon_coverage.py

```
amplicon_coverage.py --depths S.depths --bed primers.bed --sample S [--min-depth 10] -o amplicon_coverage.tsv
```

Columns: `sample, amplicon, start, end, mean_depth, min_depth, pct_ge_min_depth, dropped`.

- Amplicons come from the primer BED names. The direction is parsed from `_F/_R`, `_LEFT/_RIGHT` or a bare `F/R`
  suffix (xgen `..._88R`); `_altN` and QIAseq `-N` alternates are folded into their amplicon (`nCoV-2019_7_alt0_F` ->
  `nCoV-2019_7`, `QIAseq_221-2_LEFT` -> `QIAseq_221`). Amplicon = leftmost left-primer start to rightmost right-primer
  end. The amplicon label is the shared primer key (e.g. `nCoV-2019_12`).
- Primers without a partner of the same key are resolved: first by pairing a lone left and lone right orphan with the
  same trailing number (midnightv2 `SARSCoV_1200_22_F` / `nCoV_1200_22_R`), then by joining the amplicon with the
  nearest left start / right end (swift and xgen `M1`/`M2` start primers, midnightv2 `ARTIC_V4_1200_71_R`).
- `start`/`end` are **1-based inclusive** genome positions (BED start + 1, BED end), matching the depths file.
- `mean_depth` and `pct_ge_min_depth` are over the whole span (primers included); positions absent from the depths file
  count as 0. `dropped` = `true` when `mean_depth < --min-depth`.
- The header is written even for zero amplicons. A BED contig that differs from the depths contig is accepted when
  both files have one contig each (MN908947.3 vs NC_045512.2 are the same sequence), with a warning.

## qc_summary.py

```
qc_summary.py --sample S --sample-type sample|negative [--collection-date D] [--site X] --depths S.depths
              [--flagstat S.flagstat] [--fastp-json S.fastp.json] [--amplicons amplicon_coverage.tsv]
              [--ivar-trim-stats S.ivar_trim_stats.tsv] [--min-primer-trimmed-pct 0]
              [--primer-scheme-summary S.primer_scheme_summary.tsv]
              [--min-depth 10] [--consensus-min-depth 10] [thresholds] -o S.qc.tsv
qc_summary.py --concat a.tsv b.tsv ... -o qc_summary.tsv      # sorted by sample; headers validated; duplicates rejected
```

Columns are those of PLAN 6a plus four **appended last, in this order**: `primer_trimmed_pct`, `primer_scheme_check`,
`primer_scheme_best`, `primer_scheme_enrichment` (27 columns; every older column keeps its position).
Semantics:

- `mean_depth`, `median_depth`, `depth_sd` (population SD) are over every position in the depths file, zeros included.
- `pct_genome_ge_1x`, `pct_genome_ge_10x` are fixed thresholds; `pct_genome_ge_min_depth` uses `--consensus-min-depth`.
- `raw_reads`/`trimmed_reads`: fastp `summary.before_filtering.total_reads` / `after_filtering.total_reads`. These count
  **reads, not pairs** (a pair is 2).
- `mapped_reads`: samtools flagstat `primary mapped` (falls back to `mapped` on old samtools). Also records, so a pair
  is 2. `mapped_pct` = `mapped_reads / trimmed_reads * 100`, falling back to raw reads, then the flagstat primary total
  (on a BAM already filtered with `-F 4` the flagstat total equals mapped, so prefer supplying fastp).
- `amplicons_total`/`amplicons_dropped` are counted from `--amplicons` (rows of `--sample`).
- `primer_trimmed_pct` (0-100, 2 dp): % of mapped reads `ivar trim` trimmed primers from, taken from
  `--ivar-trim-stats` (see [ivar_trim_stats.py](#ivar_trim_statspy)). `NA` without that input (capture protocol, demix
  mode, skipped trimming) or when the value is unavailable (missing/empty/garbled log, zero reads).
- `primer_scheme_check` (`ok|mismatch|uninformative|NA`), `primer_scheme_best` (name(s) joined by `/`, or `NA`) and
  `primer_scheme_enrichment` (the **declared** scheme's enrichment, 1 dp) come from `--primer-scheme-summary` (see
  [primer_scheme_check.py](#primer_scheme_checkpy)); all `NA` without it (capture, demix, skipped trimming) or when it
  is unreadable.
- Missing/unreadable optional inputs give `NA` and a warning; the run continues. `demix_*`, `consensus_*` and
  `barcode_version` are always `NA` in the script's own output; the report-data builder fills them in
  (`build_report_data.py --write-qc-tsv`), and the pipeline publishes that completed table as `summary/qc_summary.tsv`
  (see [output.md](output.md#summaryqc_summarytsv)).

QC status (`--min-depth` is the depth used for the coverage test and is independent of the reported `ge_10x` column):

| sample_type | status | rule                                                                                                                                                     | reason tokens            |
| ----------- | ------ | -------------------------------------------------------------------------------------------------------------------------------------------------------- | ------------------------ |
| sample      | FAIL   | all depths 0 or `mapped_reads == 0`                                                                                                                      | `no_reads`               |
| sample      | FAIL   | % genome >= min-depth is `< --fail-below-pct` (30; script and pipeline default since D17, it was 50)                                                     | `low_coverage`           |
| sample      | WARN   | below `--pass-min-pct` (60; it was 80)                                                                                                                   | `low_coverage`           |
| sample      | WARN   | `amplicons_dropped > --max-dropped` (script default 10; the pipeline passes `qc_max_dropped_amplicons`, default **1**), if PASS otherwise                | `amplicon_dropout`       |
| sample      | WARN   | `primer_trimmed_pct < --min-primer-trimmed-pct` (script default **0 = rule off**), if PASS otherwise; never FAIL by itself, only when the value is known | `low_primer_trimming`    |
| sample      | WARN   | `primer_scheme_check == mismatch` (always on, no threshold; the check's own `--min-enrichment`/`--window` are module `ext.args`); never FAIL by itself   | `primer_scheme_mismatch` |
| sample      | FAIL   | depths unreadable/empty                                                                                                                                  | `no_depths`              |
| negative    | FAIL   | `mapped_reads > --ntc-max-mapped-reads` (1000) or `pct_genome_ge_1x > --ntc-max-genome-coverage-pct` (10)                                                | `ntc_reads`              |
| negative    | PASS   | otherwise (never `low_coverage`/`no_reads`)                                                                                                              |                          |
| negative    | NA     | no usable depths and no flagstat                                                                                                                         | `no_depths`              |

Reasons are `;`-joined (`low_coverage;amplicon_dropout;low_primer_trimming;primer_scheme_mismatch`); the last two are appended in this order after the coverage reasons, never
hides another reason and never turns a FAIL into a WARN. Negatives have no such rule (the value is still reported).
`--concat` rejects (exit 2) any table whose header is not exactly the current column list, so tables written before
`primer_trimmed_pct` (23 columns) or before the `primer_scheme_*` columns (24) existed are refused rather than padded.

### Primer-trimming QC (`primer_trimmed_pct`) and its limits

Purpose: catch a **wrong primer scheme** (the Snakemake pipeline could silently run with one). This ivar-based metric
does not separate right from wrong (below); the statistically calibrated check is
[primer_scheme_check.py](#primer_scheme_checkpy). The rule is a plain
"too few reads were primer-trimmed" check, and it is **weak and protocol dependent**; the default is therefore off
(`qc_min_primer_trimmed_pct` = 0). Measured with ivar 1.4.4 on the synthetic test samples (whose reads are anchored at
midnight primer starts for ~40 % of pairs, `docs/test_data.md`), same reads, different `-b` scheme:

| scheme given to ivar (data are midnight-shaped) | MIX_A            | PURE_C           | Site_1_week_2   |
| ----------------------------------------------- | ---------------- | ---------------- | --------------- |
| midnight (correct)                              | 22.33            | 21.97            | 21.70           |
| midnightv2 (same coordinates, correct)          | 22.33            | 21.97            | 21.70           |
| articv3 / articv4.1 / articv5.3.2 (wrong)       | 11.0 / 8.7 / 9.4 | 10.7 / 8.7 / 8.9 | 9.5 / 7.8 / 8.7 |
| eden (wrong)                                    | 1.59             | 1.49             | 1.29            |
| qiagenv1 (wrong for these data)                 | 20.46            | 21.73            | 21.33           |
| swift (wrong, `-x 5`)                           | 40.26            | 39.93            | 40.70           |
| xgen (wrong)                                    | 52.94            | 52.71            | 50.55           |

Wrong schemes with many, dense primers (swift, xgen, qiagenv1) trim **as much or more** than the correct sparse scheme,
because any read that happens to start inside a primer window is trimmed. So the percentage does not separate right from
wrong across schemes; it only flags schemes much sparser than the data's real one, and the expected value depends on the
scheme and, much more, on the library chemistry:

- Chemistries whose reads begin at the primer sites (e.g. QIAseq single-primer extension): a large fraction is expected
  with the right scheme, so a low value is informative. (Not measurable on the synthetic data, which are midnight-shaped;
  the expectation is unverified until a real QIAseq run is seen.)
- Tagmentation/fragmentation libraries made from long amplicons (e.g. Midnight with Nextera-style prep): fragment starts
  are random, so even the correct scheme trims only a few % of reads (rough estimate: primer windows are ~5 % of the
  genome); a fixed threshold would warn on healthy runs.

Choose a threshold only after looking at `primer_trimmed_pct` across several known-good runs of your own protocol (for
example half of the lowest good value). A strong wrong-scheme check is a different test (pooled amplicon-dropout pattern
or trying candidate schemes); this metric is only a cheap tripwire.

**Where the number comes from.** `ivar trim` writes its summary to **stderr**, but the nf-core `IVAR_TRIM` module
redirects only stdout (`> <prefix>.ivar.log`), so the published `*.ivar.log` is **empty** unless the pipeline merges stderr
into it (see [wiring/primer_trimming_qc.md](wiring/primer_trimming_qc.md)). An empty log gives `NA` with a note.

## ivar_trim_stats.py

```
ivar_trim_stats.py --log S.ivar.log --sample S -o S.ivar_trim_stats.tsv
```

Parses the summary lines of `ivar trim` (checked against ivar 1.4.4 output and binary strings; wording of 1.2-1.4.x):
`Trimmed primers from X% (N) of reads.`, `X% (N) of reads were quality trimmed below the minimum length of L bp ...` and
`X% (N) of reads [that] started outside of primer regions. ...` (with `-e`, without `-e`, with `-k`, or "no primers found
in BED file"). The per-primer table and BED warnings are ignored. One row:
`sample, reads_total, primer_trimmed_reads, primer_trimmed_pct, no_primer_reads, no_primer_pct, quality_trimmed_reads,
quality_trimmed_pct, status`. `reads_total` = trimmed + no-primer reads (ivar's denominator: mapped reads); percentages
are recomputed from the counts when both are known (ivar prints 2 dp); quality-trimmed reads overlap the other two
groups. `status`: `ok`, `missing` (no/unreadable log), `empty`, `unparsed`. Anything not found is `NA`; the script never
fails on bad input (exit 0 with a note on stderr; 2 only for bad usage). `nan` (zero reads) gives `NA`.

## primer_scheme_check.py

```
primer_scheme_check.py --starts S.read_starts.tsv --declared SCHEME_OR_BED --sample S [--beds-dir assets/primer_schemes]
                       [--window 2] [--min-enrichment 5] [--min-hit-fraction 0.05] [--min-reads 500] [--genome-length 29903]
                       [--max-p-value 0.001] [--equivalence-jaccard 0.9]
                       -o S.primer_scheme_check.tsv --summary-out S.primer_scheme_summary.tsv
```

A scheme-match check that, unlike `primer_trimmed_pct`, **can** tell a wrong scheme. In anchored libraries
(ligation/ARTIC-style, QIAseq single-primer extension) a read's 5' end coincides with a primer's 5' end on the same
strand. For every candidate scheme (all `*.bed` in `--beds-dir`, plus the declared one if a BED path) it counts the mapped
primary reads whose 5' alignment end lies within +-`--window` bp of a primer 5' end **of the same strand**, and compares
with what random read starts would give.

- `read_starts.tsv` (module `READ_STARTS`, from the mapped BAM **before** `ivar trim`, which moves 5' ends): header
  `strand, five_prime_pos`; `five_prime_pos` is the 1-based reference coordinate of the 5' aligned base (forward: `POS`;
  reverse: `POS + reference length - 1`, reference length = M/D/N/=/X CIGAR operations). Soft clips are ignored.
  `samtools view -F 2308` drops unmapped, secondary and supplementary records. Single-contig references only.
- Primer 5' ends (1-based): `+` primer = BED start + 1; `-` primer = BED end. The strand is BED column 6, or the name
  suffix (`F`/`LEFT` = `+`, `R`/`RIGHT` = `-`, a trailing `_alt<N>` ignored) when column 6 is absent (xgen). Swift's
  column 7 is ignored.
- Per strand: window set = all positions within +-W of a primer 5' end of that strand (overlapping windows of dense
  schemes are counted once); `expected = n_reads_on_strand * |window set| / genome_length`; `hits` = reads whose 5'
  position is in the set. Totals over both strands: `enrichment = hits / expected` (`NA` when expected is 0),
  `hit_fraction = hits / n_reads`, `p_value` = one-sided binomial `P(X >= hits)`, `X ~ Bin(n_reads, expected / n_reads)`
  (pooled over the strands, which is slightly conservative), computed exactly in log space (no normal approximation;
  floored at 1e-300).
- A scheme **reaches the bar** when `enrichment >= --min-enrichment` (5) **and** `p_value <= --max-p-value` (0.001) **and**
  `hit_fraction >= --min-hit-fraction` (0.05). The p-value guard is an addition: with few expected hits a handful of chance hits
  gives a large ratio. The hit-fraction floor is another: a sparse scheme (`eden`, 28 sites, window coverage 0.23 %) can reach an
  enrichment of 6-8 with only 1.4-1.8 % of the reads at its sites, because its tiny expected count inflates the ratio; the fraction of
  reads is an absolute signal that a small window set cannot inflate (`0` switches the floor off).
- Schemes are **equivalent** when the Jaccard index of their (strand, 5' position) site sets is >= 0.9 (midnight vs
  midnightv2: 58/59; swift and xgen have identical sites). Reads cannot separate them; they are reported together
  (`best_scheme = midnight/midnightv2`).
- Detail table, one row per scheme: `scheme, n_primer_5prime_sites, n_reads, hits, expected_hits, hit_fraction,
enrichment, p_value, equivalent_to_declared`. Summary row: `sample, declared_scheme, declared_enrichment, best_scheme,
best_enrichment, scheme_check, reason, declared_hit_fraction, best_hit_fraction` (the last two appended, 4 decimals, `NA` when
  not computable; `qc_summary.tsv` reads columns by name and is unchanged).

| `scheme_check`  | meaning                                                                                                                                                                                                                                       | `reason`                                                              |
| --------------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | --------------------------------------------------------------------- |
| `ok`            | the declared scheme reaches the bar                                                                                                                                                                                                           | `declared_enriched`                                                   |
| `mismatch`      | declared does not reach the bar; the best scheme does, is not equivalent to the declared one and has >= 3x its enrichment **or** >= 3x its hit fraction                                                                                       | `best_other_scheme`                                                   |
| `uninformative` | no scheme reaches the bar (e.g. fragmentation/tagmentation libraries), the best other scheme is not decisively better, or the declared scheme passes enrichment and p-value but not the hit-fraction floor (a sparse scheme's inflated ratio) | `no_scheme_enriched`, `best_not_decisive`, `declared_sparse_inflated` |
| `NA`            | fewer than `--min-reads` (500) reads, unreadable starts, or unknown declared scheme                                                                                                                                                           | `too_few_reads`, `missing_input`, `declared_not_found`                |

Exit 0 always (2 only for unusable arguments); `uninformative` is not an error and does not warn.

**Measured (bwa-mapped synthetic data, window 2).** Midnight-shaped samples (reads anchored at midnight primer starts for
part of the pairs) against all nine bundled schemes (enrichment; every p-value for midnight is at the 1e-300 floor):

| scheme                | MIX_A       | PURE_C      | Site_1_week_2 | random fragments (8000 pairs, throwaway) | capture data CAP_PURE |
| --------------------- | ----------- | ----------- | ------------- | ---------------------------------------- | --------------------- |
| midnight / midnightv2 | 42.8 / 42.1 | 41.6 / 40.9 | 41.4 / 40.7   | 1.11 / 1.09                              | 0.81 / 0.81           |
| articv3               | 2.40        | 2.19        | 1.95          | 0.94                                     | 0.89                  |
| articv4.1             | 1.16        | 1.14        | 1.06          | 1.05                                     | 1.06                  |
| articv5.3.2           | 1.26        | 1.08        | 1.06          | 0.95                                     | 1.06                  |
| eden                  | 3.31        | 2.67        | 2.43          | 1.07                                     | 0.67                  |
| qiagenv1              | 1.60        | 1.81        | 1.80          | 1.00                                     | 0.98                  |
| swift = xgen          | 1.36        | 1.36        | 1.38          | 0.97                                     | 0.99                  |

Wrong schemes stay between ~1 and 3.3. The elevated ones are real, not noise (p-values far below 0.001 for articv3,
eden, qiagenv1, swift/xgen): they share primer positions with midnight (articv3 6 exact sites of midnight's 58, eden 1,
swift/xgen 4 exact and 10 within +-2, qiagenv1 3 and 11) and the anchored reads pile onto those; eden is high only because
its tiny window set gives a tiny expected count. All stay under the bar of 5. Declared `articv4.1`, `swift`, `qiagenv1`,
`eden`, `xgen`, `articv3`, `articv5.3.2` on midnight reads are all flagged `mismatch` with best `midnight/midnightv2`;
declared `midnight` or `midnightv2` is `ok`. Randomly fragmented reads (the null) and the capture data: every scheme at
0.7-1.1, `uninformative`. Fewer than 500 reads: `NA`.

**Limits (read before trusting a verdict).**

- _Needs anchored reads._ Tagmentation/fragmentation libraries made from long amplicons (e.g. Midnight with Nextera-style
  preparation) have random read starts, so no scheme is enriched: the verdict is `uninformative`, never `ok`. The check
  can confirm a scheme or flag a mismatch only for protocols whose reads start at primer sites.
- _Sensitivity depends on the scheme's density._ Chance enrichment is `(f + (1 - f) c) / c` for an anchored fraction f and
  window coverage c (fraction of the genome inside the windows): c = 0.49 % for midnight, 4.5 % for qiagenv1, 5.6 % for
  swift/xgen (measured on random starts, window 2). Reaching 5 needs f of about 2 % for midnight but about 19 % for
  qiagenv1 and 24 % for swift/xgen; the maximum possible enrichment is 1/c (about 22 for qiagenv1).
- _Near-identical schemes cannot be told apart_ (Jaccard >= 0.9, e.g. midnight vs midnightv2; and a scheme that differs
  from the declared one by a few primers yields `ok`).
- A scheme that shares most positions with the true one (a superset, a previous version) can score close to the true one.
  `mismatch` needs the declared scheme below the bar AND the best at least 3x higher.
- Genome length defaults to 29903 (SARS-CoV-2); pass `--genome-length` for another reference. Duplicate-marked reads are
  counted (PCR duplicates of anchored reads raise the fraction).
- _QIAseq (the default scheme)_, measured on 11 real samples declared `qiagenv1`: the right scheme has enrichment 12.7-17.4 and
  **0.58-0.80** of the reads at its primers; the sparse `eden` reached enrichment 0.0-7.8 at 0.000-0.018 of the reads (the old rule gave
  3 false `ok` for `eden` declared on QIAseq reads, the hit-fraction floor and the 3x hit-fraction alternative remove them); `articv*` 0.5-3.4
  at 0.009-0.060, `swift`/`xgen` 1.0-2.0 at 0.056-0.116, `midnight*` 0.3-1.0 at 0.0015-0.005. Re-scoring all 11 samples against all
  nine declared schemes (99 combinations) gives 11/11 `ok` for the right scheme and 88/88 `mismatch` (best `qiagenv1`) for the wrong
  ones. Any floor between 0.019 and 0.19 gives the same matrix; 0.05 is the default.
- _The hit-fraction floor costs Midnight-like libraries some sensitivity_: for a very sparse scheme (window coverage 0.49 %), enrichment 5
  alone needs about 2.4 % of the reads at primer ends; with the floor it needs 5 % (enrichment about 10). A weakly anchored Midnight
  library with 2.4-5 % anchored reads is therefore `uninformative` (`declared_sparse_inflated`) instead of `ok`. Never a false warning:
  only `mismatch` warns. The synthetic Midnight-shaped samples have a hit fraction of about 0.20.

## dropout_impact.py

```
dropout_impact.py --amplicons amplicon_coverage.tsv --aggregated freyja_aggregated.tsv --dataset <dataset dir | barcodes .csv>
                  [--min-abundance 0.01] [--min-affect-abundance 0.05] [--informative-max-fraction 0.95]
                  [--lineages-yml PATH] [--strict-roots] --out dropout_impact.tsv --summary-out dropout_impact_summary.tsv
```

Which dropped amplicons hide mutations that matter for the lineages detected in a sample. Output columns:
[output.md](output.md#summarydropout_impacttsv-and-summarydropout_impact_summarytsv). Definitions:

- an **informative** site is a single-nucleotide substitution carried by some but not nearly all lineages of the dataset (< 95 %);
- a **defining** site of lineage `l` is an informative site that `l` carries and its **parent** does not
  (`barcode(l) - barcode(parent(l))`). The parent is the `parent:` of `lineages.yml` (by name or alias, case-insensitive), the nearest
  ancestor with a barcode if the parent has none; a root has all its informative sites as defining. **Recombinant roots**
  (entries with `recombinant_parents:`, e.g. `XFG`, `XBB`) inherit the sites of their recombinant parents: those are not defining
  (`--strict-roots` counts them as defining; it is a script option, not exposed as a pipeline parameter, and is only useful to compare
  with the literal definition). Without `lineages.yml` (a bare barcodes `.csv`) the parent is guessed by dropping trailing `.N` name
  components, with a warning: aliased names (`BA.5`) then act as roots and defining sites are over-counted;
- a lineage is **detected** at abundance `>= --min-abundance` (0.01); `affects_detected_lineage` needs one at
  `>= --min-affect-abundance` (0.05) with a defining site in the dropped span. `n_differential_sites` (carried by some but not all
  detected lineages) is reported but is not part of the flag: with many detected lineages nearly every amplicon has one.

Why defining and not informative: a Freyja barcode lists hundreds of mutations inherited from ancestors, so under the "any informative
site" rule 24 of 29, 6 of 10, 17 of 22, 26 of 36 and 19 of 26 dropped amplicons of five real samples were flagged and the run-level
warning fired on every sample; with the defining rule 12, 2, 5, 2 and 0 are. The 0.05 floor only stops lineages at 1-5 % from
raising the flag; they are still listed. Runs in 0.2 s per run; the new columns are appended, the old ones keep name and position.

## unexplained_mutations.py

```
unexplained_mutations.py --variants S.variants.tsv --dataset <dataset> --sample S [--demix S.demix] [--gff g.gff3 --reference r.fa]
                         [--min-depth 20] [--min-alt-depth 10] [--min-freq 0.05] [--min-excess 0.05] [--max-pvalue 1.67e-6]
                         [--max-rows 100] --output S.unexplained_mutations.tsv
unexplained_mutations.py --concat S1.tsv S2.tsv ... [--min-recurrence 2] --output unexplained_mutations.tsv
```

Frequent variants that the demixed lineage mixture does not predict (method and columns: [usage.md](usage.md#unexplained-mutations),
[output.md](output.md#summaryunexplained_mutationstsv)). Notes on the numbers: the binomial tail is summed in **log space**, so deep calls
do not underflow to `0`; `p_value` is written with 3 significant digits and **floored at `1e-300`** (it means "at most 1e-300"; always a
plain positive float), and `--max-pvalue` is compared in log space (`--max-pvalue 0` rejects everything). `in_dataset` (last column) says
whether any lineage of the dataset carries the exact mutation. `--max-rows` (default 100; negative values are refused, `0` writes a header
only) keeps the rows with the largest excess per sample and warns; the cut is made before `--concat`, so `n_samples` and `recurrent` count only
the rows that survived it, and ties on excess (many `excess = 1`) are cut by position, not by evidence. `--concat` refuses per-sample tables
with a different header (for example from before `in_dataset`).

## low_depth_mask.py

```
low_depth_mask.py --depths S.depths [--min-depth 10] [--keep-bed variants.bed] [--length 29903] -o mask.bed
```

BED (0-based, half-open, merged) of positions with `depth < --min-depth`; contig from the depths file. `--keep-bed`
positions (BED3+; exact contig match) are subtracted, replacing the old `bedtools subtract` of confidently called
variants. With `--length`, positions not listed in the depths file (gaps, tail) have depth 0 and are masked (single
contig only). An empty mask is a legitimately empty file; an empty/unreadable depths file is an error.

## consensus_stats.py

```
consensus_stats.py --fasta S.consensus.fa --sample S [--reference-length 29903] -o S.stats.tsv
consensus_stats.py --concat a.tsv b.tsv ... -o consensus_stats.tsv
```

Columns `sample, length, n_count, consensus_coverage_pct, ambiguous_count, status, reason`. `n_count` counts `N`/`n`/`-`; `ambiguous_count` counts
other non-ACGT characters; `consensus_coverage_pct` (100 minus the % `N`; `0.00` when `empty`, `NA` when `failed`) uses `--reference-length` as the denominator when given (a length mismatch only
warns). Status: `ok`; `empty` (all N/`-`); `failed` (missing/unreadable/not FASTA/not exactly one record/zero length;
numeric fields are `NA`). Multi-line FASTA is supported.

## mutations_of_interest.py

```
mutations_of_interest.py --variants S.variants.tsv --moi moi.csv --sample S [--depths S.depths]
                         [--gff NC_045512.gff3 --reference NC_045512.fasta] -o S.moi.tsv
```

`moi.csv` has columns `name,mutation`. Mutation syntax:

| form       | example                               | notes                                                                                                    |
| ---------- | ------------------------------------- | -------------------------------------------------------------------------------------------------------- |
| SNV        | `A23063T`                             | ref base, 1-based position, alt base; ref must equal the TSV `REF`                                       |
| deletion   | `C75-`, `AGTTCA22029-`                | deleted bases + position of the **first deleted base** (ivar row: `POS` = that - 1, `ALT` = `-` + bases) |
| insertion  | `A14+C`                               | ref base at the anchor position + inserted bases (ivar row: `POS` = anchor)                              |
| amino acid | `S:N501Y`, `ORF1ab:P4715L`, `N:R203*` | CDS/transcript name, ref residue, protein position, alt residue                                          |

Amino-acid resolution: if the TSV row has `GFF_FEATURE` equal to the gene and `POS_AA`/`ALT_AA`, those are used;
otherwise the SNV is translated with the GFF3 CDS model (exon-spliced, ORF1ab slippage handled: ORF1ab positions beyond
the 13468 slippage site use the repeated base, so `P4715L` is C14408T) and the reference. Names are the transcript
`Name` (`ORF1a` and `ORF1ab` are distinct). SNV rows are translated individually; several SNVs in one codon are not
combined; rows giving the same amino-acid change are summed (alt depth, frequency capped at 1). Only `+` strand CDS
and SNVs are supported for amino-acid mutations. AA mutations need `--gff` and `--reference` when the TSV is not annotated.

Output columns `sample, name, mutation, ref_depth, alt_depth, depth, alt_freq` (one row per MOI line, input order).
Present: `REF_DP, ALT_DP, TOTAL_DP, ALT_FREQ` of the matching row (the `PASS` column is not applied). Absent from the
file: `alt_depth 0`, `alt_freq 0`; with `--depths`, `depth` = depth at the position (amino acid: minimum over the codon)
and `ref_depth` = `depth`, else both `NA`. Variants file missing/unreadable: all values `NA` plus a warning. Malformed
mutation, unknown gene, or a protein position beyond the CDS: exit 2. A reference-residue mismatch warns.
