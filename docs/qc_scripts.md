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
              [--min-depth 10] [--consensus-min-depth 10] [thresholds] -o S.qc.tsv
qc_summary.py --concat a.tsv b.tsv ... -o qc_summary.tsv      # sorted by sample; headers validated; duplicates rejected
```

Columns are exactly those of PLAN 6a. Semantics:

- `mean_depth`, `median_depth`, `depth_sd` (population SD) are over every position in the depths file, zeros included.
- `pct_genome_ge_1x`, `pct_genome_ge_10x` are fixed thresholds; `pct_genome_ge_min_depth` uses `--consensus-min-depth`.
- `raw_reads`/`trimmed_reads`: fastp `summary.before_filtering.total_reads` / `after_filtering.total_reads`. These count
  **reads, not pairs** (a pair is 2).
- `mapped_reads`: samtools flagstat `primary mapped` (falls back to `mapped` on old samtools). Also records, so a pair
  is 2. `mapped_pct` = `mapped_reads / trimmed_reads * 100`, falling back to raw reads, then the flagstat primary total
  (on a BAM already filtered with `-F 4` the flagstat total equals mapped, so prefer supplying fastp).
- `amplicons_total`/`amplicons_dropped` are counted from `--amplicons` (rows of `--sample`).
- Missing/unreadable optional inputs give `NA` and a warning; the run continues. `demix_*`, `consensus_*` and
  `barcode_version` are always `NA` in the script's own output; the report-data builder fills them in
  (`build_report_data.py --write-qc-tsv`), and the pipeline publishes that completed table as `summary/qc_summary.tsv`
  (see [output.md](output.md#summaryqc_summarytsv)).

QC status (`--min-depth` is the depth used for the coverage test and is independent of the reported `ge_10x` column):

| sample_type | status | rule                                                                                                                                      | reason tokens      |
| ----------- | ------ | ----------------------------------------------------------------------------------------------------------------------------------------- | ------------------ |
| sample      | FAIL   | all depths 0 or `mapped_reads == 0`                                                                                                       | `no_reads`         |
| sample      | FAIL   | % genome >= min-depth is `< --fail-below-pct` (50)                                                                                        | `low_coverage`     |
| sample      | WARN   | below `--pass-min-pct` (80)                                                                                                               | `low_coverage`     |
| sample      | WARN   | `amplicons_dropped > --max-dropped` (script default 10; the pipeline passes `qc_max_dropped_amplicons`, default **1**), if PASS otherwise | `amplicon_dropout` |
| sample      | FAIL   | depths unreadable/empty                                                                                                                   | `no_depths`        |
| negative    | FAIL   | `mapped_reads > --ntc-max-mapped-reads` (1000) or `pct_genome_ge_1x > --ntc-max-genome-coverage-pct` (10)                                 | `ntc_reads`        |
| negative    | PASS   | otherwise (never `low_coverage`/`no_reads`)                                                                                               |                    |
| negative    | NA     | no usable depths and no flagstat                                                                                                          | `no_depths`        |

Reasons are `;`-joined (`low_coverage;amplicon_dropout`).

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

Columns `sample, length, n_count, n_pct, ambiguous_count, status`. `n_count` counts `N`/`n`/`-`; `ambiguous_count` counts
other non-ACGT characters; `n_pct` uses `--reference-length` as the denominator when given (a length mismatch only
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
