# Consensus fix: competing ALT alleles at one position are merged into one IUPAC call

Hand-over for the wiring/docs owner. No workflow, `conf/modules.config` or schema change: the same `ext.*` thresholds are
used. Edited: `modules/local/consensus_filter_vcf/main.nf` (+ tests), `modules/local/consensus_bcftools/main.nf` (+ tests),
`tests/python/test_consensus_multiallelic.py` (new), `tests/python/test_ivar_variants_to_vcf.py` (one test),
`tests/fixtures/consensus_multiallelic/` (new), `docs/usage.md` and `docs/output.md` (IUPAC paragraphs only).
`bin/ivar_variants_to_vcf.py` is unchanged (it still writes one record per ALT, frequencies intact).

## Problem

ivar writes one row per ALT. `bcftools consensus` applies one record per position and skips the others
(`The site X overlaps with another variant, skipping...`), so the base was wrong for the mixture: WW-1275E 22942, ref T,
`T>G 0.6849` + `T>A 0.2941` gave `K` (G/T, only the first ALT applied together with the reference) instead of `R` (A/G).

## Rule as implemented (`CONSENSUS_FILTER_VCF`)

1. `bcftools norm -a -m -` then `bcftools view -i 'FORMAT/ALT_QUAL >= q && INFO/DP >= d'` keep the candidate records
   (per-record thresholds only; the AF rules moved into step 2).
2. awk, per position (consecutive records with the same CHROM/POS; SNVs with the same REF are one site):
   - **Indel** (REF and ALT of different length): unchanged. AF >= `indel_freq` -> GT `1/1`, else dropped. Never merged,
     never IUPAC. Two indels at the same anchor stay two records.
   - **SNV site**: f_i = ALT_FREQ of each candidate ALT, f_ref = max(0, 1 - sum f_i) (rounded to 6 decimals).
     - (a) exactly one ALT with f >= `consensus_freq` -> that ALT alone, GT `1/1` (any other ALT at the site is dropped).
       More than one ALT >= `consensus_freq` (only possible when `consensus_freq` <= 0.5): with IUPAC on they go to (b),
       with IUPAC off the highest-frequency one is applied.
     - (b) IUPAC on only: S = ALTs with f > `iupac_min_af` and ALT_DP >= `iupac_min_alt_reads`, plus the reference if
       f_ref > `iupac_min_af`. |S| >= 2 -> **one record**, `ALT` = the S alts ordered by decreasing frequency, GT lists every
       member (ref = 0): `0/1`, `1/2`, `0/1/2`, `1/2/3`, `0/1/2/3`. FORMAT values are those of the highest-frequency member.
     - (c) |S| == 1 and it is an ALT (so f_ref <= `iupac_min_af`: the reference is not a contributor): a site with exactly
       one candidate row -> GT `0/1` (that ALT with the reference: the single-ALT rule of the old module, unchanged, e.g.
       `consensus_freq` 0.9 with one ALT at 0.8); a site with two or more candidate rows -> that ALT alone, GT `1/1`
       (plain base: a reference that is not in the mixture is not paired with it). No new threshold: f_ref and
       `iupac_min_af` are compared directly.
     - otherwise (S empty or only the reference, or IUPAC off without (a)) -> no record (reference base).
   - Comparisons with `iupac_min_af` are strict (`>`), as in the old single-ALT expression (an ALT at exactly 0.25 is
     still ignored; the spec text said ">= 0.25", the existing test "AF == .25 -> reference" pins `>`). `consensus_freq`
     stays `>=`.
   - ALTs failing DP/ALT_QUAL are not candidates and do not count in f_ref (the mass is treated as reference); an ALT
     with ALT_DP < `iupac_min_alt_reads` is not an IUPAC member but does count in the sum (so f_ref is not inflated).
3. Single-ALT sites and all indels produce byte-identical VCF lines to the previous module
   (`test_single_alt_output_is_byte_identical_to_the_previous_filter`: 200 random sites x {defaults, IUPAC off,
   consensus_freq 0.5} compared line by line with the old `view -i` + awk filter).

Worked examples (all synthetic, run through the real chain): `T>G 0.6849 + T>A 0.2941` -> `R`; `T>G 0.576 + T>C 0.424` -> `S`;
`C>G 0.61 + C>A 0.39` -> `R`; three ALTs 0.40/0.30/0.30 -> the three-base code (`B`/`V`/...); ref 0.30 + 0.40 + 0.30 -> the
three-base code; `0.80 + 0.15` -> the plain 0.80 base; `0.20 + 0.20` (ref 0.60) -> reference; an allele with ALT_DP < 5 is
left out; `iupac false` -> only an allele >= 0.75.

Example for (c): `A 0.58 + C 0.24 + G 0.18`, reference ~0 (WW-1278E 22896, ref T): only A is a member (C is below the strict
0.25 floor) and f_ref = 0, so the site is the plain base `A`, not `W` (A/T). With `iupac_min_af` 0.1 the same site is the
three-base code (A, C, G), with `consensus_freq` 0.5 plain `A` by rule (a). Other effects of the same rule on run1: a
second row that is not a member now leaves the ALT plain when the reference is also at or below `iupac_min_af` (e.g.
`0.6 + 0.15`, f_ref 0.25, gives the plain ALT, where the old module gave ALT/reference).

## What bcftools does with polyploid genotypes (bcftools 1.23, `consensus -H I`)

Tested in `quay.io/biocontainers/bcftools:1.23--h3a4d415_0` on a hand-made VCF over `ACGTACGT...`: polyploid GTs ARE
honoured and give the union of the listed alleles, independent of allele order.

| REF | ALT | GT | consensus |
| --- | --- | --- | --- |
| A | C,G | `1/2` | `S` |
| C | G,T | `1/2` | `K` |
| G | A,T | `0/1/2` | `D` (G,A,T) |
| T | A,C,G | `0/1/2/3` | `N` (all four) |
| A | C,G,T | `1/2/3` | `B` |
| C | A,G | `0/1/2` | `V` |
| G | C,T | `0/1/2` | `B` |
| T | A,C | `0/1/2` | `H` |
| A | G | `0/1` | `R` |
| T | A,G | `1/2` | `R` |
| G | A,T | `1/2` | `W` |

So no fallback to "the two highest-frequency members" is needed. If a future bcftools stopped honouring them, the pytest
`test_cases_default` (three/four-member sites) and the nf-test "multi-allelic sites: one IUPAC code ..." fail.
`--mark-del '-'` and `-m mask` are unaffected (masks only change positions outside the called variants).

## `CONSENSUS_BCFTOOLS`

Only the wording of one warning changed. A site that bcftools skips because it lies inside a deletion that was applied
(site in `(POS, END]` of a deletion record that was not itself skipped) is listed as
`skipped variant at <contig>:<pos> (inside a deleted region)`; every other overlap keeps
`skipped variant at <contig>:<pos> (overlaps another variant)`. It stays in `warnings` and in the `n variant(s) not applied`
reason, so `consensus_stats.tsv.reason` and the report INFO flag `consensus_warnings` still fire for it.
Suggestion for the report owner: do not raise the run-level INFO `consensus_warnings` for warnings that all end in
`(inside a deleted region)` (expected, e.g. an SNV at 28370 inside the 9-base deletion 28362-28370), or show them in a
separate "expected" tone. Existing test `SNV inside a deleted span` was updated to the new wording.

## Real-data check (read-only on run1, scratch outputs only; before = HEAD modules, after = new modules, same inputs)

| Sample | Position | Old | New | Why |
| --- | --- | --- | --- | --- |
| WW-1274E | 23604 | S | R | C>G 0.61 + C>A 0.39: {G,A} |
| WW-1275E | 22942 | K | R | T>G 0.68 + T>A 0.29: {G,A} |
| WW-1275E | 23485 | W | A | G 0.056 + C 0.001 + A 0.728 (f_ref 0.21, one member): plain A (rule c, absent reference) |
| WW-1277E | 22896 | K | S | T>G 0.58 + T>C 0.42: {G,C} |
| WW-1277E | 22942 | W | A | G 0.27 (3 reads, not a member) + A 0.55 (f_ref 0.18): plain A (rule c) |
| WW-1278E | 22896 | W | A | A 0.58 + C 0.24 + G 0.18, reference absent: plain A (rule c) |
| WW-1280E | 22896 | K | S | T>C 0.56 + T>G 0.43: {C,G} |
| WW-1276E, 1279E, 1281E, 1282E | - | - | - | no consensus difference |

All 9 samples: status `ok`, length 29903 before and after, every `-` is identical (WW-1275E: 20 in total, the nine at
11288-11296 among them). Warnings before: 9 in 7 samples (22942, 23604, 22896 x2, 28370 x4, 29759); after: 5, all
`(inside a deleted region)` notes (WW-1276E 28370; WW-1278E 28370 and 29759; WW-1280E 28370; WW-1282E 28370); no
`(overlaps another variant)` warnings remain in run1. (The published run1 consensus files predate the consensus fix, e.g.
WW-1275E is `failed` there, so they are not a usable "before".)

## Proposed text for `docs/PLAN.md` section 5.11 (replace the rule sentence; I did not edit PLAN)

> Rules (parity with the old pipeline's `setGT` chain, plus a read-count guard borrowed from h2seq): DP >=
> `consensus_min_depth_variants` (10) and ALT_QUAL >= 20 for any call. **SNVs are decided per position, not per record**:
> ivar writes one row per ALT but `bcftools consensus` applies one record per position, so `CONSENSUS_FILTER_VCF` merges the
> SNV rows of a position into one record. With f_i the ALT frequencies and f_ref = max(0, 1 - sum f_i): one ALT with f >=
> `consensus_freq` (0.75) -> that alt allele; otherwise (IUPAC on) the set S of alleles with f > `consensus_iupac_min_af`
> (0.25) (an ALT also needs ALT_DP >= `consensus_iupac_min_alt_reads` (5); the reference counts with f_ref) -> if |S| >= 2 a
> single record with ALT = the S alts and a polyploid genotype listing all members (`0/1`, `1/2`, `0/1/2`, `1/2/3`,
> `0/1/2/3`), emitted as one IUPAC code by `bcftools consensus -H I` (T>G 0.68 + T>A 0.29 gives `R`, not `K`); if S is a
> single ALT below `consensus_freq`, the site is ALT/reference (`0/1`) when it has one candidate row (as before) and the plain
> ALT (`1/1`) when several rows compete and the reference is not a contributor (f_ref <= `iupac_min_af`). `--consensus_iupac false` applies
> only an allele >= `consensus_freq`. Indels (AF >= `indel_freq`, 0.75) are never merged and never IUPAC. A SNV inside an
> applied deletion is skipped by bcftools (the deletion wins) and noted in the status JSON `warnings` as
> `(inside a deleted region)`.

## Other doc follow-ups (outside my edit scope)

- `docs/output.md` row for `variants/<s>.consensus_variants.vcf.gz`: updated (done).
- `docs/PLAN.md` section 5.11 as above.
