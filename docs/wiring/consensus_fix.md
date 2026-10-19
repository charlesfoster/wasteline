# Consensus fix: indel alleles from the reference, bcftools retry, status/reason columns

Hand-over for the wiring agent (workflow, `conf/modules.config`, report). Nothing outside the consensus modules,
`bin/ivar_variants_to_vcf.py`, `bin/consensus_stats.py`, their tests and fixtures was edited.

## Root cause

ivar (as run by `freyja variants`) prints the deleted bases of a deletion as `N`: the row
`NC_045512.2 11287 G -NNNNNNNNN ... ALT_FREQ 0.988636` is the common 9-base Omicron deletion at 11288-11296. The
old converter built the VCF record from the printed bases (`REF=GNNNNNNNNN ALT=G`), but the reference at 11287..11296 is
`GTCTGGTTTT`. `CONSENSUS_FILTER_VCF` kept the record (AF 0.99), `bcftools consensus` aborted with
`The fasta sequence does not match the REF allele at NC_045512.2:11287`, and the whole genome became N with status
`failed`. Every `-N...` deletion in the TSV was affected (293 of 11050 records of the failing sample had a REF that did
not match the reference; only the AF >= 0.75 ones reach the consensus, but one is enough).

## 1. `IVAR_VARIANTS_TO_VCF` (signature change)

```nextflow
IVAR_VARIANTS_TO_VCF(ch_consensus_samples, ch_fasta)        // was: IVAR_VARIANTS_TO_VCF(ch_consensus_samples)
```

- New second input `path fasta` (the reference). `[]` is accepted and keeps the old behaviour with a loud stderr warning
  (`WARNING: no --reference given ...`); the pipeline must pass `ch_fasta`.
- The module adds `--reference <fasta>` to `bin/ivar_variants_to_vcf.py`. No `ext.args` / `params` change.
- Alleles (ivar POS `p` is the anchor base, 1-based): deletion `-X`: `L = len(X)`, `REF = ref[p-1 : p+L]`,
  `ALT = ref[p-1]` (the printed bases are ignored, only their number is used); insertion `+X`: `REF = ref[p-1]`,
  `ALT = ref[p-1] + X`; SNV unchanged.
- Skipped (counted on stderr, no row, no failure): insertions containing anything but A/C/G/T (e.g. `+ANT`: an unresolved
  inserted base must not enter a consensus), indels on a contig not in the reference, deletions with nothing left to
  delete. A deletion overrunning the contig end is truncated to the contig end. A printed REF differing from the reference
  is counted (`ref_mismatch`); for SNVs the record is left as is and handled by `CONSENSUS_BCFTOOLS` (see 2).
- Also fixed: a 14-column TSV (no trailing annotation columns) made every row `FAIL`.
- Existing module test `modules/local/annotate_variants/tests/main.nf.test` calls `IVAR_VARIANTS_TO_VCF` in its `setup`
  with one input: add `input[1] = file("\$projectDir/assets/references/NC_045512.fasta")` there (I may not edit it).

## 2. `CONSENSUS_BCFTOOLS` (no signature change)

- When `bcftools consensus` fails with `... does not match the REF allele at <contig>:<pos>`, the records at exactly that
  position are dropped from the VCF and bcftools is re-run (max 10 retries, one offending position per retry). The sample
  ends `ok` with the dropped variants in `warnings`. More than 10 offending positions, or any other error, gives
  `failed` + all-N as before (`warnings` is then `[]`; the retries stay in `<sample>.consensus.log`).
- Variants that bcftools itself skips (`The site X overlaps with another variant, skipping...`: a SNV inside a deleted
  span, or two ALT alleles of one site both applied) are also listed in `warnings`; this is not an error and needs no retry.
- `consensus_status.json` gains `warnings` (always present, array of strings, max 20 + a `(N more)` entry):

  ```json
  {
    "sample": "S",
    "status": "ok",
    "reason": "1 variant(s) not applied: skipped variant at NC_045512.2:11287 (REF mismatch)",
    "warnings": ["skipped variant at NC_045512.2:11287 (REF mismatch)"]
  }
  ```

  Warning texts: `skipped variant at <contig>:<pos> (REF mismatch)` and
  `skipped variant at <contig>:<pos> (overlaps another variant)`. `reason` stays `""` for a clean `ok`, is the old
  "no variant passed ..." text when nothing was applied, and is `"<n> variant(s) not applied: <w1>; <w2>"` with warnings.
  `status` values are unchanged (`ok|empty|failed`).

- Deletions are written as `-` (`--mark-del '-'`), genome length is kept: the 9-base deletion at 11287 gives
  `G` + 9 x `-` at 11287-11296.

## 3. `CONSENSUS_STATS` (signature change) and `CONSENSUS_STATS_CONCAT`

```nextflow
CONSENSUS_STATS(CONSENSUS_BCFTOOLS.out.fasta.join(CONSENSUS_BCFTOOLS.out.status), ch_fasta)
// was: CONSENSUS_STATS(CONSENSUS_BCFTOOLS.out.fasta, ch_fasta)
```

- The input tuple is now `[ meta, fasta, status_json ]` (`status_json` may be `[]`: then status/reason are recomputed
  from the sequence, as before). The module adds `--status-json <file>` to `bin/consensus_stats.py`.
- Precedence: JSON `failed` -> `failed`, JSON `empty` -> `empty` (a failed sample is stored as an all-N genome, which alone
  read as `empty`); JSON `ok` (or no JSON) -> status recomputed from the sequence (`ok`, or `empty` if all N/`-`, or
  `failed` if the FASTA is unreadable). Numeric columns are computed from the FASTA whenever it is readable.
- `consensus_stats.tsv` has a new last column: `sample, length, n_count, n_pct, ambiguous_count, status, reason`.
  `reason` is `NA` for a clean `ok`; otherwise the one-line text (JSON reason incl. the warnings text; or the parse error
  for an unreadable FASTA; `consensus contains only N or - characters` for a computed `empty`). All existing columns and
  their order are unchanged.
- `--concat` accepts the new 7-column header; a legacy 6-column file is still accepted (reason `NA`), so the workflow's
  synthesised failed row does not break the run, but update it:

  ```nextflow
  [ "${meta.id}.consensus_stats.tsv", "sample\tlength\tn_count\tn_pct\tambiguous_count\tstatus\treason\n${meta.id}\tNA\tNA\tNA\tNA\tfailed\tconsensus steps did not complete\n" ]
  ```

- Also update `docs/PLAN.md` (the `consensus_stats.tsv` contract line) and `tests/check_outputs.py` if it pins the header.

## 4. Report implications (for the report builder / template)

- `bin/build_report_data.py::load_consensus` reads columns by name; add `"reason": clean_str(row.get("reason"))` and show it
  in the consensus table (tooltip or a column), most visibly for non-`ok` rows.
- Run-level flags:
  - **WARN `consensus_failed`**: any sample with consensus `status == failed`; list the samples (and their reason). The
    consensus of such a sample is all-N, the sample has no usable genome.
  - **INFO `consensus_warnings`**: any sample with `status == ok` and a non-`NA` reason, i.e. variants were skipped; list
    sample + reason. Expected occasionally (e.g. WW-1275E: two ALT alleles at 22942 both applied, one skipped by bcftools).
  - `empty` keeps its existing handling.
