> **Design-era note.** Written before final integration; the implemented defaults/params are in `nextflow_schema.json` and `docs/parameters.md`.

# Phase 3 modules: wiring notes

Local modules for coverage/amplicon QC, consensus genomes and mutations of interest (PLAN section 5 steps 9-12,
section 6a, D7, D8). Script semantics are in `docs/qc_scripts.md`. Every module has a `stub:` block, a `topic: versions`
tuple and an nf-test (`modules/local/<name>/tests/main.nf.test`). The shared wiring (workflow, `conf/modules.config`,
`nextflow.config`, `nextflow_schema.json`) is deliberately not touched by the module author; this file is the hand-over.

Containers (all exact tags, verified on quay.io, the Galaxy depot singularity mirror, and run under linux/amd64 Docker):

| tool                                        | docker                                                                                                         | singularity                                                             | used by                                  |
| ------------------------------------------- | -------------------------------------------------------------------------------------------------------------- | ----------------------------------------------------------------------- | ---------------------------------------- |
| Python 3.12.12 (stdlib scripts from `bin/`) | `quay.io/biocontainers/python:3.12.12`                                                                         | `https://depot.galaxyproject.org/singularity/python:3.12.12`            | all Python modules                       |
| bcftools 1.23 (+ htslib 1.23: bgzip/tabix)  | `quay.io/biocontainers/bcftools:1.23--h3a4d415_0`                                                              | `https://depot.galaxyproject.org/singularity/bcftools:1.23--h3a4d415_0` | CONSENSUS_FILTER_VCF, CONSENSUS_BCFTOOLS |
| samtools 1.24                               | nf-core module `samtools/flagstat` (`community.wave.seqera.io/library/htslib_samtools:1.24--d697cfb9dce007cd`) | (module default)                                                        | SAMTOOLS_FLAGSTAT                        |

## Includes

```nextflow
include { SAMTOOLS_FLAGSTAT           } from '../modules/nf-core/samtools/flagstat/main'
include { AMPLICON_COVERAGE           } from '../modules/local/amplicon_coverage/main'
include { QC_SAMPLE                   } from '../modules/local/qc_sample/main'
include { QC_CONCAT                   } from '../modules/local/qc_concat/main'
include { IVAR_VARIANTS_TO_VCF        } from '../modules/local/ivar_variants_to_vcf/main'
include { CONSENSUS_FILTER_VCF        } from '../modules/local/consensus_filter_vcf/main'
include { LOW_DEPTH_MASK              } from '../modules/local/low_depth_mask/main'
include { CONSENSUS_BCFTOOLS          } from '../modules/local/consensus_bcftools/main'
include { CONSENSUS_STATS             } from '../modules/local/consensus_stats/main'
include { CONSENSUS_STATS_CONCAT      } from '../modules/local/consensus_stats_concat/main'
include { CONSENSUS_CONCAT            } from '../modules/local/consensus_concat/main'
include { MUTATIONS_OF_INTEREST       } from '../modules/local/mutations_of_interest/main'
include { MUTATIONS_OF_INTEREST_CONCAT} from '../modules/local/mutations_of_interest_concat/main'
```

(Paths relative to `workflows/wasteline.nf`.) `modules.json` already lists `samtools/flagstat`.

## Channel signatures

`meta` is the samplesheet map (`id`, and for QC_SAMPLE also `sample_type`, `collection_date`, `site`; the last three
may be null/absent). Optional inputs are given as `[]`.

| process                        | inputs                                                                                                      | outputs (emit)                                                                                                               |
| ------------------------------ | ----------------------------------------------------------------------------------------------------------- | ---------------------------------------------------------------------------------------------------------------------------- |
| `SAMTOOLS_FLAGSTAT` (nf-core)  | `tuple val(meta), path(bam), path(bai)`                                                                     | `flagstat` = `[meta, <id>.flagstat]`                                                                                         |
| `AMPLICON_COVERAGE`            | `tuple val(meta), path(depths)`; `path bed`                                                                 | `tsv` = `[meta, <id>.amplicon_coverage.tsv]`                                                                                 |
| `QC_SAMPLE`                    | `tuple val(meta), path(depths), path(flagstat) \| [], path(fastp_json) \| [], path(amplicons) \| []`        | `tsv` = `[meta, <id>.qc_row.tsv]`                                                                                            |
| `QC_CONCAT`                    | `path qc_rows` (collected files); `path amplicon_tables` (collected files)                                  | `qc_summary` = `qc_summary.tsv`; `amplicon_coverage` = `amplicon_coverage.tsv`                                               |
| `IVAR_VARIANTS_TO_VCF`         | `tuple val(meta), path(variants)` (Freyja `<id>.variants.tsv`)                                              | `vcf` = `[meta, <id>.ivar.vcf]`                                                                                              |
| `CONSENSUS_FILTER_VCF`         | `tuple val(meta), path(vcf)`                                                                                | `vcf` = `[meta, <id>.consensus_variants.vcf.gz, ...vcf.gz.tbi]`; `bed` = `[meta, <id>.called_variants.bed]`                  |
| `LOW_DEPTH_MASK`               | `tuple val(meta), path(depths), path(keep_bed) \| []`; `path fasta`                                         | `bed` = `[meta, <id>.mask.bed]`                                                                                              |
| `CONSENSUS_BCFTOOLS`           | `tuple val(meta), path(vcf), path(tbi), path(mask)`; `path fasta`                                           | `fasta` = `[meta, <id>.consensus.fa]`; `status` = `[meta, <id>.consensus_status.json]`; `log` = `[meta, <id>.consensus.log]` |
| `CONSENSUS_STATS`              | `tuple val(meta), path(fasta)`; `path reference`                                                            | `tsv` = `[meta, <id>.consensus_stats.tsv]`                                                                                   |
| `CONSENSUS_STATS_CONCAT`       | `path stats` (collected files)                                                                              | `tsv` = `consensus_stats.tsv`                                                                                                |
| `CONSENSUS_CONCAT`             | `path fastas` (collected files)                                                                             | `fasta` = `consensus_genomes.fasta`                                                                                          |
| `MUTATIONS_OF_INTEREST`        | `tuple val(meta), path(variants), path(depths) \| []`; `path moi`; `path gff \| []`; `path reference \| []` | `tsv` = `[meta, <id>.mutations_of_interest.tsv]`                                                                             |
| `MUTATIONS_OF_INTEREST_CONCAT` | `path tables` (collected files)                                                                             | `tsv` = `mutations_of_interest.tsv`                                                                                          |

`CONSENSUS_FILTER_VCF` is an addition to the brief: PLAN 5.11 needs the called-variant positions _before_ the low-depth
mask can be built (they are subtracted from it, replacing the old `bedtools subtract`), and the mask must exist before
`bcftools consensus`, so the filter step cannot live inside `CONSENSUS_BCFTOOLS`. Chain:
`IVAR_VARIANTS_TO_VCF -> CONSENSUS_FILTER_VCF -> LOW_DEPTH_MASK -> CONSENSUS_BCFTOOLS`.

All concat modules reject mismatched headers and (QC_CONCAT/CONSENSUS_STATS_CONCAT) duplicate samples; all sort by
sample ID in byte order (`LC_ALL=C`), matching the Python scripts. Never call them with zero files (guard with
`.collect().filter { it }` or the module-level `if` of the caller; with no samples there is nothing to summarise).

## Workflow sketch

```nextflow
// ch_depths    : [meta, depths]            FREYJA_VARIANTS.out.depths
// ch_variants  : [meta, variants]          FREYJA_VARIANTS.out.variants
// ch_bam_bai   : [meta, bam, bai]          primer-trimmed BAM + index
// ch_fastp     : [meta, fastp_json]        FASTP.out.json
// ch_primer_bed, ch_fasta, ch_gff          value channels / files

// -- coverage QC (every sample incl. negatives) --
SAMTOOLS_FLAGSTAT(ch_bam_bai)
AMPLICON_COVERAGE(ch_depths, ch_primer_bed)

ch_qc_in = ch_depths
    .join(SAMTOOLS_FLAGSTAT.out.flagstat, remainder: true)
    .join(ch_fastp, remainder: true)
    .join(AMPLICON_COVERAGE.out.tsv, remainder: true)
    .map { meta, depths, flagstat, fastp, amplicons -> [ meta, depths, flagstat ?: [], fastp ?: [], amplicons ?: [] ] }
QC_SAMPLE(ch_qc_in)

QC_CONCAT(
    QC_SAMPLE.out.tsv.map { _meta, f -> f }.collect(),
    AMPLICON_COVERAGE.out.tsv.map { _meta, f -> f }.collect()
)   // -> summary/qc_summary.tsv, summary/amplicon_coverage.tsv

// -- consensus (skip with params.skip_consensus; samples with sample_type == 'negative' are not built) --
IVAR_VARIANTS_TO_VCF(ch_variants_consensus)
CONSENSUS_FILTER_VCF(IVAR_VARIANTS_TO_VCF.out.vcf)
LOW_DEPTH_MASK(
    ch_depths.join(CONSENSUS_FILTER_VCF.out.bed).map { meta, depths, bed -> [ meta, depths, bed ] },
    ch_fasta
)
CONSENSUS_BCFTOOLS(CONSENSUS_FILTER_VCF.out.vcf.join(LOW_DEPTH_MASK.out.bed), ch_fasta)
CONSENSUS_STATS(CONSENSUS_BCFTOOLS.out.fasta, ch_fasta)
CONSENSUS_STATS_CONCAT(CONSENSUS_STATS.out.tsv.map { _meta, f -> f }.collect())
CONSENSUS_CONCAT(CONSENSUS_BCFTOOLS.out.fasta.map { _meta, f -> f }.collect())

// -- optional: mutations of interest (params.mutations_of_interest) --
MUTATIONS_OF_INTEREST(ch_variants.join(ch_depths), file(params.mutations_of_interest), ch_gff, ch_fasta)
MUTATIONS_OF_INTEREST_CONCAT(MUTATIONS_OF_INTEREST.out.tsv.map { _meta, f -> f }.collect())
```

Notes for the wiring:

- `join` needs identical meta maps on both sides (all of these channels descend from the same samplesheet map).
- **Empty samples must still enter the consensus chain.** PLAN section 8 says a sample with ~0 mapped reads must not abort
  `IVAR_TRIM`/`FREYJA_VARIANTS` and yields empty-but-valid outputs. Feed those (a header-only variants TSV and an empty or
  all-zero depths file) to `IVAR_VARIANTS_TO_VCF`/`LOW_DEPTH_MASK`: the chain turns them into an all-N consensus with
  `status: "empty"` (tested). `QC_SAMPLE` likewise turns an empty depths file into `FAIL` / `no_depths`; `qc_summary.py`
  with empty depths and a flagstat still rates a negative control.
- `IVAR_VARIANTS_TO_VCF`, `CONSENSUS_FILTER_VCF` and `LOW_DEPTH_MASK` fail loudly on malformed (not merely empty) input.
  Only `CONSENSUS_BCFTOOLS` converts its own failure to `status: "failed"` + an all-N FASTA. If the wiring uses
  `errorStrategy 'ignore'` on the upstream steps, a sample with a malformed variants TSV simply has no consensus row;
  the report-data builder should treat a missing `consensus_status.json` as `failed`.
- `consensus_stats.tsv` is the per-sample consensus table (status `ok|empty|failed` from the sequence itself).
  `<id>.consensus_status.json` (`{sample, status, reason}`) carries the _why_; it is published per sample and is not
  concatenated by any module.
- Negatives (D8): do not build consensus genomes for `sample_type == 'negative'` (filter the channel). They still go
  through `AMPLICON_COVERAGE`, `SAMTOOLS_FLAGSTAT` and `QC_SAMPLE`.
- `QC_SAMPLE` is told `--consensus-min-depth` (for `pct_genome_ge_min_depth`) and `LOW_DEPTH_MASK` `--min-depth`; both
  must come from the same param (`consensus_min_depth`).
- `SAMTOOLS_FLAGSTAT` needs BAM + `.bai`. Run it on the primer-trimmed BAM (same file `FREYJA_VARIANTS` sees); the BAM is
  already `-F 4` filtered, so `mapped_pct` should come from fastp reads (supply `fastp_json`).
- Singularity: all modules honour `task.ext.singularity_pull_docker_container` like the existing ones.

## `conf/modules.config` blocks

Add inside the existing `process { ... }` block.

```groovy
    // ---- coverage QC -----------------------------------------------------------------------
    withName: 'SAMTOOLS_FLAGSTAT' {
        ext.prefix = { "${meta.id}.primertrimmed" }
        publishDir = [
            path: { "${params.outdir}/samples/${meta.id}/coverage" },
            mode: params.publish_dir_mode
        ]
    }

    withName: 'AMPLICON_COVERAGE' {
        ext.args   = { "--min-depth ${params.amplicon_min_depth}" }
        publishDir = [
            path: { "${params.outdir}/samples/${meta.id}/coverage" },
            mode: params.publish_dir_mode
        ]
    }

    withName: 'QC_SAMPLE' {
        ext.args = {
            [
                "--min-depth ${params.qc_min_depth}",
                "--consensus-min-depth ${params.consensus_min_depth}",
                "--pass-min-pct ${params.qc_pass_min_pct}",
                "--fail-below-pct ${params.qc_fail_below_pct}",
                "--max-dropped ${params.qc_max_dropped_amplicons}",
                "--ntc-max-mapped-reads ${params.ntc_max_mapped_reads}",
                "--ntc-max-genome-coverage-pct ${params.ntc_max_genome_coverage_pct}",
            ].join(' ')
        }
        publishDir = [
            path: { "${params.outdir}/samples/${meta.id}/coverage" },
            mode: params.publish_dir_mode
        ]
    }

    withName: 'QC_CONCAT' {
        publishDir = [
            path: { "${params.outdir}/summary" },
            mode: params.publish_dir_mode
        ]
    }

    // ---- consensus ---------------------------------------------------------------------------
    withName: 'IVAR_VARIANTS_TO_VCF' {
        // no --pass_only / --allele_freq_thresh: thresholds are applied in CONSENSUS_FILTER_VCF
        publishDir = [
            path: { "${params.outdir}/samples/${meta.id}/variants" },
            mode: params.publish_dir_mode
        ]
    }

    withName: 'CONSENSUS_FILTER_VCF' {
        ext.consensus_freq = { params.consensus_freq }
        ext.indel_freq     = { params.indel_freq }
        ext.min_alt_qual   = { params.consensus_min_alt_qual }
        ext.min_depth      = { params.consensus_min_depth_variants }   // DP of a call; NOT the mask depth (consensus_min_depth)
        ext.iupac               = { params.consensus_iupac }
        ext.iupac_min_af        = { params.consensus_iupac_min_af }
        ext.iupac_min_alt_reads = { params.consensus_iupac_min_alt_reads }
        publishDir = [
            path: { "${params.outdir}/samples/${meta.id}/variants" },
            mode: params.publish_dir_mode,
            pattern: '*.vcf.gz*'                  // the called-variants BED is an intermediate
        ]
    }

    withName: 'LOW_DEPTH_MASK' {
        ext.args   = { "--min-depth ${params.consensus_min_depth}" }
        publishDir = [
            path: { "${params.outdir}/samples/${meta.id}/consensus" },
            mode: params.publish_dir_mode
        ]
    }

    withName: 'CONSENSUS_BCFTOOLS' {
        publishDir = [
            path: { "${params.outdir}/samples/${meta.id}/consensus" },
            mode: params.publish_dir_mode
        ]
    }

    withName: 'CONSENSUS_STATS' {
        publishDir = [
            path: { "${params.outdir}/samples/${meta.id}/consensus" },
            mode: params.publish_dir_mode
        ]
    }

    withName: 'CONSENSUS_STATS_CONCAT|CONSENSUS_CONCAT' {
        publishDir = [
            path: { "${params.outdir}/summary" },
            mode: params.publish_dir_mode
        ]
    }

    // ---- mutations of interest (opt-in) ----------------------------------------------------
    withName: 'MUTATIONS_OF_INTEREST' {
        publishDir = [
            path: { "${params.outdir}/samples/${meta.id}/variants" },
            mode: params.publish_dir_mode
        ]
    }

    withName: 'MUTATIONS_OF_INTEREST_CONCAT' {
        publishDir = [
            path: { "${params.outdir}/summary" },
            mode: params.publish_dir_mode
        ]
    }
```

If the output tree should not carry the intermediate raw `.ivar.vcf`, drop the `IVAR_VARIANTS_TO_VCF` block (nothing is
published without a `publishDir`). `CONSENSUS_FILTER_VCF` reads its seven `ext.*` threshold keys directly (defaults
0.75 / 0.75 / 20 / 10 / true / 0.25 / 5 when unset, so the module also works unconfigured), the
Python modules take their flags through `ext.args`.

## Params

Add to `nextflow.config` `params { }` and to `nextflow_schema.json` (group them under "Coverage QC" and "Consensus").

| param                           | default                                        | schema type                                    | used by                   | meaning                                                                                 |
| ------------------------------- | ---------------------------------------------- | ---------------------------------------------- | ------------------------- | --------------------------------------------------------------------------------------- |
| `amplicon_min_depth`            | `10`                                           | integer, min 1                                 | AMPLICON_COVERAGE         | amplicon dropped when mean depth is below this                                          |
| `qc_min_depth`                  | `10`                                           | integer, min 1                                 | QC_SAMPLE                 | depth for the PASS/WARN/FAIL coverage test                                              |
| `qc_pass_min_pct`               | `80`                                           | number, 0-100                                  | QC_SAMPLE                 | PASS needs at least this % genome at `qc_min_depth`                                     |
| `qc_fail_below_pct`             | `50`                                           | number, 0-100                                  | QC_SAMPLE                 | FAIL below this % genome at `qc_min_depth`                                              |
| `qc_max_dropped_amplicons`      | `1` (was 10 in this design note; see PLAN §6a) | integer, min 0                                 | QC_SAMPLE                 | WARN (`amplicon_dropout`) when more amplicons dropped                                   |
| `ntc_max_mapped_reads`          | `1000`                                         | integer, min 0                                 | QC_SAMPLE                 | negative control FAIL (`ntc_reads`) above this many mapped reads                        |
| `ntc_max_genome_coverage_pct`   | `10`                                           | number, 0-100                                  | QC_SAMPLE                 | negative control FAIL above this % genome at >= 1x                                      |
| `consensus_min_depth`           | `10` (design note said 15; see PLAN)           | integer, min 1                                 | LOW_DEPTH_MASK, QC_SAMPLE | positions below this depth are masked (N); also `pct_genome_ge_min_depth`               |
| `consensus_freq`                | `0.75`                                         | number, 0-1                                    | CONSENSUS_FILTER_VCF      | SNV ALT_FREQ needed to apply the alt allele                                             |
| `indel_freq`                    | `0.75`                                         | number, 0-1                                    | CONSENSUS_FILTER_VCF      | indel ALT_FREQ needed to apply it                                                       |
| `consensus_min_alt_qual`        | `20`                                           | integer, min 0                                 | CONSENSUS_FILTER_VCF      | minimum ALT_QUAL                                                                        |
| `consensus_min_depth_variants`  | `10`                                           | integer, min 1                                 | CONSENSUS_FILTER_VCF      | minimum TOTAL_DP (INFO/DP) at a variant (not the mask depth `consensus_min_depth`)      |
| `consensus_iupac`               | `true`                                         | boolean                                        | CONSENSUS_FILTER_VCF      | emit IUPAC codes for ambiguous SNVs; false leaves the reference base                    |
| `consensus_iupac_min_af`        | `0.25`                                         | number, 0-1                                    | CONSENSUS_FILTER_VCF      | SNV ALT_FREQ must be strictly above this (and below `consensus_freq`) for an IUPAC code |
| `consensus_iupac_min_alt_reads` | `5`                                            | integer, min 1                                 | CONSENSUS_FILTER_VCF      | ALT_DP needed for an IUPAC code                                                         |
| `skip_consensus`                | `false`                                        | boolean                                        | workflow                  | skip the whole consensus chain and the consensus summaries                              |
| `mutations_of_interest`         | `null`                                         | string, format file-path, pattern `^\S+\.csv$` | workflow                  | CSV `name,mutation`; when null MUTATIONS_OF_INTEREST(+CONCAT) are not run               |

The reference FASTA / GFF come from the existing `params.fasta` / `params.gff`; the primer BED from the existing
primer-scheme resolution. No module hard-codes a contig or reference name.

## What the consensus does (PLAN 5.11, ported from `ivar_bcftools_setGT` / `generate_consensus`)

1. `ivar_variants_to_vcf.py` converts the Freyja TSV (all ALTs, `--varthresh 0.0`) to a VCF (no thresholds).
2. `bcftools norm -a -m -` then one `bcftools view -i` keeps a record only if `FORMAT/ALT_QUAL >= 20` and `INFO/DP >=
consensus_min_depth_variants`, and either (SNV, `TYPE!="indel"`) ALT_FREQ `>= consensus_freq`, or (SNV, only when
   `consensus_iupac`) `consensus_iupac_min_af < ALT_FREQ < consensus_freq` with `ALT_DP >= consensus_iupac_min_alt_reads`,
   or (indel) ALT_FREQ `>= indel_freq`. An awk step then sets the genotype: `1/1` (alt allele) for AF `>= consensus_freq`
   and all indels, `0/1` for the IUPAC candidates, which `bcftools consensus -H I` writes as an ambiguity code (A/G = R,
   C/T = Y, ...). Everything else (SNV AF `<= consensus_iupac_min_af`, too few alt reads, indels below `indel_freq`,
   and every minor ALT) is not applied; indels are never IUPAC. No bcftools plugins are used.
3. `low_depth_mask.py --min-depth consensus_min_depth --length <reference length> --keep-bed <called variants>`.
4. `bcftools consensus --mark-del '-' -H I -m <mask>`; header rewritten to the sample ID.
5. If every position outside the called variants is below the minimum depth (or the depths file is empty) the output is an
   all-N sequence of reference length with `status: "empty"`; if `bcftools consensus` fails the same sequence is written
   with `status: "failed"` and the first line of the bcftools error as `reason`. A sample with good coverage whose
   filtered VCF is empty (identical to the reference) is `status: "ok"` with a reason noting that no variant passed.
