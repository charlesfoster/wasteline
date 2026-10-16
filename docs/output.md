# Output

Reference for everything wasteline writes to `--outdir`. All examples below are trimmed from a real `--mode full` run
on the synthetic dataset in `tests/data` (`--barcodes tests/data/barcodes_test.csv --primer_scheme midnight
--compare_barcodes container --freyja_boot --annotate_variants`, Freyja 2.0.6, Nextflow 26.04.6; the surveillance tables come from the same data and the `tests/data_timeseries` series). The one exception is
the `mutations_of_interest.tsv` example, which comes from an earlier real run on the same data with
`--mutations_of_interest moi.csv` (the final run did not use that option). Paths are relative to `--outdir`. Parameter
names link to [parameters.md](parameters.md); how to read the HTML report is in [report.md](report.md).

Contents: [directory tree](#directory-tree) | [what each mode writes](#what-each-mode-writes) |
[summary files](#summary-files) | [per-sample files](#per-sample-files) | [consensus FASTA conventions](#consensus-fasta-conventions) |
[pipeline_info](#pipeline_info) | [run_info.json](#run_infojson) | [how to reproduce a run](#how-to-reproduce-a-run)

Conventions used by all tables:

- TSV files are tab-separated with a header line and UTF-8 text. Missing values are `NA`. Percentages are 0-100.
- Sample IDs are exactly the IDs of the samplesheet (letters, digits, `.`, `_`, `-`); they are never rewritten, so
  `Site_1_week_2` stays `Site_1_week_2`.
- A failed or skipped sample always gets a row with a status and a reason; it never silently disappears and never
  becomes an empty file that looks like a result.
- Results are published with `copy` mode (`--publish_dir_mode`) so that an output directory is self-contained and can be
  fed back to `--mode demix`, `--mode aggregate` or `--mode report`.

## Directory tree

Real output of the run described above (7 samples; `MIX_A` and the two negative controls are expanded, the other
samples have the same shape as `MIX_A`, minus the files a skipped sample does not get). Files marked "only with" exist
only when that option is on.

```
results/
├── summary/
│   ├── report.html                  self-contained HTML report
│   ├── report_data.json             machine-readable input of the report (docs/report_data.md)
│   ├── qc_summary.tsv               one QC row per sample (incl. negative controls), complete (demix/consensus/barcode filled in)
│   ├── amplicon_coverage.tsv        per-amplicon depth and dropout, long format
│   ├── demix_status.tsv             demix outcome per sample (ok / failed / skipped_low_coverage / missing)
│   ├── freyja_aggregated.tsv        Freyja aggregate: lineage abundances per sample (negative controls excluded)
│   ├── freyja_aggregated_controls.tsv  the same table for the negative controls that demixed
│   ├── freyja_plot.pdf              `freyja plot` stacked bars
│   ├── freyja_plot_times.pdf        `freyja plot --times` (samples with a collection_date)
│   ├── freyja_timeseries.tsv        site x date x lineage-group abundance (>= 2 collection dates)
│   ├── freyja_timeseries_lineages.tsv  the same at the level of individual lineages
│   ├── freyja_growth_rates.csv      only with --freyja_growth: relative growth rates (.failed.txt says why not; freyja_growth_inputs.json)
│   ├── freyja_dashboard.html        only with --freyja_dash: Freyja's interactive dashboard
│   ├── dropout_impact.tsv           dropped amplicons and the lineage-informative sites they hide
│   ├── dropout_impact_summary.tsv   per-sample totals of the above
│   ├── unexplained_mutations.tsv    frequent variants the lineage mixture does not predict (usually empty)
│   ├── freyja_boot_summary.tsv      only with --freyja_boot: 95% bootstrap interval per lineage
│   ├── barcode_compare.tsv          only with --compare_barcodes: abundance with barcode set A and B
│   ├── barcode_compare_status.tsv   only with --compare_barcodes: outcome of the second demix per sample
│   ├── annotated_variants.tsv       only with --annotate_variants: amino-acid consequences, all samples
│   ├── consensus_genomes.fasta      all consensus genomes in one file
│   ├── consensus_stats.tsv          length / N / ambiguity / status per consensus
│   ├── mutations_of_interest.tsv    only with --mutations_of_interest
│   ├── barcode_version.json         barcode dataset record (sidecar of freyja_aggregated.tsv)
│   └── samplesheet.csv              the samplesheet in effect (copied or discovered)
├── samples/
│   ├── MIX_A/
│   │   ├── fastp/        MIX_A.fastp.json  MIX_A.fastp.html  MIX_A.fastp.log
│   │   ├── bam/          MIX_A.primertrimmed.bam  MIX_A.primertrimmed.bam.bai  MIX_A.ivar_trimmed.ivar.log
│   │   │                 with --bam_format cram: MIX_A.primertrimmed.cram  MIX_A.primertrimmed.cram.crai (instead of the BAM)
│   │   ├── freyja/       MIX_A.variants.tsv  MIX_A.depths  MIX_A.demix  MIX_A.demix_status.json  MIX_A.depth_stats.json
│   │   │                 with --compare_barcodes: MIX_A.compare.demix  MIX_A.compare.demix_status.json  MIX_A.compare.depth_stats.json
│   │   │                 with --freyja_boot: MIX_A.boot_summary.tsv  MIX_A.boot_lineages.csv  MIX_A.boot_summarized.csv  MIX_A.boot_status.json
│   │   ├── coverage/     MIX_A.primertrimmed.flagstat  MIX_A.amplicon_coverage.tsv  MIX_A.qc_row.tsv
│   │   │                 MIX_A.primer_scheme_check.tsv  MIX_A.primer_scheme_summary.tsv
│   │   ├── variants/     MIX_A.ivar.vcf  MIX_A.consensus_variants.vcf.gz(.tbi)  MIX_A.mask.bed  MIX_A.unexplained_mutations.tsv
│   │   │                 with --annotate_variants: MIX_A.annotated.tsv;  with --mutations_of_interest: MIX_A.mutations_of_interest.tsv
│   │   └── consensus/    MIX_A.consensus.fa  MIX_A.consensus_stats.tsv  MIX_A.consensus_status.json  MIX_A.consensus.log
│   ├── NEG_CONTAM/       (contaminated negative control, demixed ok: has .demix and boot files, but no variants/ or consensus/)
│   │   ├── fastp/  bam/  coverage/
│   │   └── freyja/       NEG_CONTAM.variants.tsv  NEG_CONTAM.depths  NEG_CONTAM.demix  NEG_CONTAM.demix_status.json  NEG_CONTAM.depth_stats.json
│   │                     (+ the .compare.* and .boot_* files)
│   ├── NEG_E/            (clean negative control, skipped_low_coverage: no .demix, no boot files, no variants/ or consensus/)
│   │   ├── fastp/  bam/  coverage/
│   │   └── freyja/       NEG_E.variants.tsv (header only)  NEG_E.depths  NEG_E.demix_status.json  NEG_E.depth_stats.json
│   │                     NEG_E.compare.demix_status.json  NEG_E.compare.depth_stats.json
│   └── ... (LOWCOV_D, MIX_B, PURE_C, Site_1_week_2)
└── pipeline_info/
    ├── run_info.json                full provenance (see below)
    ├── freyja_data.json             the freyja_data block on its own
    ├── freyja_data_compare.json     only with --compare_barcodes: the same for the second dataset
    ├── software_versions.yml        tool versions per process
    ├── params_<timestamp>.json      the parameters of the run (usable as -params-file)
    ├── execution_report_<timestamp>.html
    ├── execution_timeline_<timestamp>.html
    ├── execution_trace_<timestamp>.txt
    └── pipeline_dag_<timestamp>.html
```

`<timestamp>` is `--trace_report_suffix` (default: the launch date and time, `yyyy-MM-dd_HH-mm-ss`). The whole run above
is about 15 MB (13 MB under `samples/`, mostly the per-sample BAM and fastp HTML).

Downloaded Freyja data is **not** written to `--outdir`: it is stored in `--db_dir` (default `~/.wasteline/db`) under
`freyja/<YYYY-MM-DDTHH-MM>/`, see [Freyja data](#freyja-data---db_dir).

## What each mode writes

| Path                                                                                         |               `full`               |                                                            `demix`                                                            |                       `aggregate`                       |                                                       `report`                                                        |   `update_db`    | Controlled by                                                         |
| -------------------------------------------------------------------------------------------- | :--------------------------------: | :---------------------------------------------------------------------------------------------------------------------------: | :-----------------------------------------------------: | :-------------------------------------------------------------------------------------------------------------------: | :--------------: | --------------------------------------------------------------------- |
| `summary/report.html`                                                                        |                yes                 |                                                              yes                                                              |                           yes                           |                                                     yes (rebuilt)                                                     |        no        | `--skip_report`                                                       |
| `summary/report_data.json`                                                                   |                yes                 |                                                              yes                                                              |                           yes                           | no (an existing one in `--results_dir` is reused; if absent it is rebuilt from the other summary files and published) |        no        | `--skip_report`                                                       |
| `summary/qc_summary.tsv`                                                                     |                yes                 |                                                              yes                                                              |                           no                            |                                                          no                                                           |        no        | `--qc_*`, `--ntc_*`, `--amplicon_min_depth`, `--consensus_min_depth`  |
| `summary/amplicon_coverage.tsv`                                                              |                yes                 |                                         only with `--primer_scheme` or `--primer_bed`                                         |                           no                            |                                                          no                                                           |        no        | `--primer_scheme`, `--primer_bed`, `--amplicon_min_depth`             |
| `summary/demix_status.tsv`                                                                   |                yes                 |                                                              yes                                                              | only if some sample is `missing` (then only those rows) |                                                          no                                                           |        no        | `--min_genome_coverage_for_demix`                                     |
| `summary/freyja_aggregated.tsv`, `freyja_plot.pdf`                                           |                yes                 |                                                              yes                                                              |                           yes                           |                                                          no                                                           |        no        | `--freyja_plot_*`, `--include_negatives_in_aggregate`                 |
| `summary/freyja_plot_times.pdf`                                                              |       with collection dates        |                                                     with collection dates                                                     |                  with collection dates                  |                                                          no                                                           |        no        | `--freyja_plot_interval`, `--freyja_plot_windowsize`                  |
| `summary/freyja_aggregated_controls.tsv`                                                     | if a negative control demixed `ok` |                                                             same                                                              |                          same                           |                                                          no                                                           |        no        | not written with `--include_negatives_in_aggregate`                   |
| `summary/freyja_timeseries.tsv`, `freyja_timeseries_lineages.tsv`                            |     with >= 2 collection dates     |                                                             same                                                              |                          same                           |                                                          no                                                           |        no        | `--timeseries`, `--timeseries_min_abundance`                          |
| `summary/freyja_growth_rates.csv` (or `.failed.txt`), `freyja_growth_inputs.json`            |       with `--freyja_growth`       |                                                             same                                                              |                          same                           |                                                          no                                                           |        no        | `--freyja_growth*`                                                    |
| `summary/freyja_dashboard.html` (or `.failed.txt`)                                           |        with `--freyja_dash`        |                                                             same                                                              |                          same                           |                                                          no                                                           |        no        | `--freyja_dash*`                                                      |
| `summary/dropout_impact.tsv`, `dropout_impact_summary.tsv`                                   |      yes (amplicon protocol)       |                                         only with `--primer_scheme` or `--primer_bed`                                         |                           no                            |                                                          no                                                           |        no        | `--dropout_impact`, `--dropout_min_abundance`                         |
| `summary/unexplained_mutations.tsv`                                                          |                yes                 |                                                              yes                                                              |                           no                            |                                                          no                                                           |        no        | `--unexplained_mutations`, `--unexplained_*`                          |
| `summary/freyja_boot_summary.tsv`                                                            |        with `--freyja_boot`        |       with `--freyja_boot` <!-- UNVERIFIED: demix mode with --freyja_boot not run; read from workflows/wasteline.nf -->       |                           no                            |                                                          no                                                           |        no        | `--freyja_boot`, `--freyja_boot_n`, `--freyja_boot_seed`              |
| `summary/barcode_compare.tsv`, `barcode_compare_status.tsv`                                  |     with `--compare_barcodes`      |                                                   with `--compare_barcodes`                                                   |                           no                            |                                                          no                                                           |        no        | `--compare_barcodes`                                                  |
| `summary/annotated_variants.tsv`                                                             |     with `--annotate_variants`     | with `--annotate_variants` <!-- UNVERIFIED: demix mode with --annotate_variants not run; read from workflows/wasteline.nf --> |                           no                            |                                                          no                                                           |        no        | `--annotate_variants` (needs consensus)                               |
| `summary/consensus_genomes.fasta`, `consensus_stats.tsv`                                     |                yes                 |                                                              yes                                                              |                           no                            |                                                          no                                                           |        no        | `--skip_consensus`, `--consensus_*`, `--indel_freq`                   |
| `summary/mutations_of_interest.tsv`                                                          |   with `--mutations_of_interest`   |                                                with `--mutations_of_interest`                                                 |                           no                            |                                                          no                                                           |        no        | `--mutations_of_interest`                                             |
| `summary/barcode_version.json`                                                               |                yes                 |                                                              yes                                                              |                           no                            |                                                          no                                                           |        no        | `--barcodes`                                                          |
| `summary/samplesheet.csv`                                                                    |                yes                 |                                                              yes                                                              |                           yes                           |                                                          no                                                           |        no        |                                                                       |
| `samples/<s>/fastp/`                                                                         |                yes                 |                                                              no                                                               |                           no                            |                                                          no                                                           |        no        | `--fastp_*`; trimmed FASTQ only with `--save_trimmed_reads`           |
| `samples/<s>/bam/`                                                                           |                yes                 |                                                              no                                                               |                           no                            |                                                          no                                                           |        no        | `--save_bam` (the `.ivar_trimmed.ivar.log` is always written)         |
| `samples/<s>/freyja/*.variants.tsv`, `*.depths`                                              |                yes                 |                                                     no (not re-published)                                                     |                           no                            |                                                          no                                                           |        no        |                                                                       |
| `samples/<s>/freyja/*.demix`, `*.demix_status.json`, `*.depth_stats.json`                    |                yes                 |                                                              yes                                                              |                           no                            |                                                          no                                                           |        no        | `--freyja_*`                                                          |
| `samples/<s>/freyja/*.compare.*`                                                             |     with `--compare_barcodes`      |                                                   with `--compare_barcodes`                                                   |                           no                            |                                                          no                                                           |        no        | `--compare_barcodes`                                                  |
| `samples/<s>/freyja/*.boot_*`                                                                |        with `--freyja_boot`        |                                                     with `--freyja_boot`                                                      |                           no                            |                                                          no                                                           |        no        | `--freyja_boot`                                                       |
| `samples/<s>/coverage/*.qc_row.tsv`                                                          |                yes                 |                                                              yes                                                              |                           no                            |                                                          no                                                           |        no        |                                                                       |
| `samples/<s>/coverage/*.flagstat`, `*.amplicon_coverage.tsv`                                 |                yes                 |                                      no (`amplicon_coverage` only with a primer scheme)                                       |                           no                            |                                                          no                                                           |        no        |                                                                       |
| `samples/<s>/variants/`, `consensus/`                                                        |                yes                 |                                                              yes                                                              |                           no                            |                                                          no                                                           |        no        | `--skip_consensus`; `*.annotated.tsv` only with `--annotate_variants` |
| `pipeline_info/run_info.json`, `freyja_data.json`, `software_versions.yml`, Nextflow reports |                yes                 |                                                              yes                                                              |                           yes                           |                                   Nextflow reports and `software_versions.yml` only                                   | Nextflow reports |                                                                       |
| `pipeline_info/freyja_data_compare.json`                                                     |     with `--compare_barcodes`      |                                                   with `--compare_barcodes`                                                   |                           no                            |                                                          no                                                           |        no        | `--compare_barcodes`                                                  |
| `pipeline_info/discovery_<n>.log`                                                            |                 no                 |                                                   with `--previous_results`                                                   |   with `--previous_results` (not with `--demix_dir`)    |                                                          no                                                           |        no        |                                                                       |

Notes verified on real `demix` and `aggregate` runs of the output above:

- In `demix` mode the new run's `samples/<s>/freyja/` holds only `.demix`, `.demix_status.json` and
  `.depth_stats.json` (plus the `.compare.*` files with `--compare_barcodes`). The variants and depths files are read from the earlier results and **not copied**, so the
  output directory of a `demix` run cannot itself be used as `--previous_results` for another `demix` run; point the
  next run at the original `full` results.
- In `aggregate` mode QC and consensus are not available, so the report has no QC, coverage, amplicon or consensus
  panels. The **time-series panel and `freyja_plot_times.pdf` do work** when collection dates are known: they come from
  the samplesheet in effect (an `--input` sheet with `collection_date`, or the `summary/samplesheet.csv` of a
  `--previous_results` directory). `--demix_dir` carries no metadata, so such a run has neither (verified: an
  `--previous_results` run had a 4-date time series, the `--demix_dir` run none).
- Samples that `--previous_results` discovery skips (no usable `.demix`, for example `LOWCOV_D` and the control `NEG_E`
  above) are listed as `missing` in `demix_status.tsv` but are not in the discovered samplesheet, so the report treats
  them as ordinary samples without date or site (`NEG_E` is flagged `demix_missing` like a sample, not shown as a
  control).
- Negative controls are demixed (if they pass the coverage gate) and have a `.demix` file and a `demix_status.tsv` row,
  but are left out of `freyja_aggregated.tsv` and both plots unless `--include_negatives_in_aggregate`. What they contain
  is written to `summary/freyja_aggregated_controls.tsv` instead, and shown in the report's "Negative controls" section.
  An `aggregate` run applies the same rule to the `.demix` files it is given (a control demix from an earlier run is
  left out of the aggregate and listed in its own `freyja_aggregated_controls.tsv`) unless that flag is set.

## Summary files

### `summary/qc_summary.tsv`

One row per sample, including negative controls, sorted by `sample`. Written in `full` and `demix` modes. Produced from
the Freyja depths file (Q>=20 bases, one row per genome position), `samtools flagstat` (mode `full` only) and the
fastp JSON (when known); the last five of the original 23 columns are filled in afterwards from the demix and consensus results (see the
table). The file has **27 columns**: the four primer columns at the end (`primer_trimmed_pct`, `primer_scheme_check`,
`primer_scheme_best`, `primer_scheme_enrichment`) are described below. Real rows (trimmed to the first 18 columns):

```
sample         sample_type collection_date site   raw_reads trimmed_reads mapped_reads mapped_pct mean_depth depth_sd median_depth pct_genome_ge_1x pct_genome_ge_10x pct_genome_ge_min_depth amplicons_total amplicons_dropped qc_status qc_reasons
LOWCOV_D       sample      2026-01-26      SiteB  120       120           119          99.17      0.55       0.77     0.00         39.73            0.00              0.00                    29              29                FAIL      low_coverage;amplicon_dropout
MIX_A          sample      2026-01-05      SiteA  12000     12000         11998        99.98      54.86      24.14    50.00        99.37            99.33             98.62                   29              0                 PASS
NEG_CONTAM     negative    2026-01-19      SiteB  6000      6000          5998         99.97      27.39      14.88    25.00        99.37            94.55             83.96                   29              0                 FAIL      ntc_reads
Site_1_week_2  sample      2026-02-02      Site_1 13000     13000         12998        99.98      59.55      41.86    51.00        93.06            93.06             93.06                   29              2                 WARN      amplicon_dropout
```

and the remaining columns of the same file (the demix, consensus and barcode columns are complete in the published file):

```
sample         qc_status  qc_reasons                      demix_status          demix_reason                                   consensus_length  consensus_n_pct  barcode_version
LOWCOV_D       FAIL       low_coverage;amplicon_dropout   skipped_low_coverage  less than 60% of the genome covered at >= 10x  29903             100.00           custom
MIX_A          PASS                                       ok                    NA                                             29903             1.38             custom
NEG_CONTAM     FAIL       ntc_reads                       ok                    NA                                             NA                NA               custom
NEG_E          PASS                                       skipped_low_coverage  less than 60% of the genome covered at >= 10x  NA                NA               custom
Site_1_week_2  WARN       amplicon_dropout                ok                    NA                                             29903             6.94             custom
```

and the four primer-scheme columns (mode `full`, amplicon protocol with primer trimming; `NA` otherwise):

```
sample         primer_trimmed_pct  primer_scheme_check  primer_scheme_best   primer_scheme_enrichment
LOWCOV_D       22.69               NA                   NA                   43.0
MIX_A          22.33               ok                   midnight/midnightv2  42.8
NEG_CONTAM     22.79               ok                   midnight/midnightv2  43.3
NEG_E          60.00               NA                   NA                   103.1
Site_1_week_2  21.70               ok                   midnight/midnightv2  41.4
```

The merged table is written by the report builder (`BUILD_REPORT_DATA`, `--write-qc-tsv`), which is why it exists only
when the report is built: with `--skip_report` the published `qc_summary.tsv` is the raw QC output and the last five
columns stay `NA` <!-- UNVERIFIED: --skip_report not run; derived from conf/modules.config (QC_CONCAT publishes qc_summary.tsv only when skip_report is set) -->.

| Column                                   | Meaning                                                                                                                                                                                                                                                               |
| ---------------------------------------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `sample`                                 | Sample ID from the samplesheet.                                                                                                                                                                                                                                       |
| `sample_type`                            | `sample` or `negative` (negative control / NTC).                                                                                                                                                                                                                      |
| `collection_date`                        | ISO `YYYY-MM-DD` from the samplesheet, or `NA`.                                                                                                                                                                                                                       |
| `site`                                   | Free-text site label from the samplesheet, or `NA`.                                                                                                                                                                                                                   |
| `raw_reads`, `trimmed_reads`             | Reads (not pairs: a pair counts 2) before and after fastp, from the fastp JSON. `NA` in `demix` mode without a fastp JSON.                                                                                                                                            |
| `mapped_reads`                           | `samtools flagstat` "primary mapped" of the primer-trimmed BAM (a pair counts 2). `NA` in `demix` mode.                                                                                                                                                               |
| `mapped_pct`                             | `mapped_reads / trimmed_reads * 100` (falls back to raw reads, then to the flagstat total).                                                                                                                                                                           |
| `mean_depth`, `depth_sd`, `median_depth` | Mean, population SD and median of the per-position depth over **every** position of the depths file, zeros included.                                                                                                                                                  |
| `pct_genome_ge_1x`, `pct_genome_ge_10x`  | % of the genome with depth >= 1 and >= 10 (fixed thresholds).                                                                                                                                                                                                         |
| `pct_genome_ge_min_depth`                | % of the genome with depth >= `--consensus_min_depth` (default 10), i.e. the part that will not be masked `N` in the consensus.                                                                                                                                       |
| `amplicons_total`, `amplicons_dropped`   | Amplicons of the primer scheme, and how many have mean depth `< --amplicon_min_depth` (see `amplicon_coverage.tsv`). `NA` when no primer scheme is known (`demix` mode without `--primer_scheme`).                                                                    |
| `qc_status`                              | `PASS`, `WARN`, `FAIL`, or `NA` (a negative control with no usable data).                                                                                                                                                                                             |
| `qc_reasons`                             | `;`-joined machine tokens, empty when PASS with no reason (see below).                                                                                                                                                                                                |
| `demix_status`, `demix_reason`           | The sample's row of `demix_status.tsv` (`ok`, `skipped_low_coverage`, ...), with its reason; `demix_reason` is `NA` for `ok`.                                                                                                                                         |
| `consensus_length`, `consensus_n_pct`    | From `consensus_stats.tsv`: length and % `N` of the consensus. `NA` for negative controls (they get no consensus), for a `failed` consensus and with `--skip_consensus`.                                                                                              |
| `barcode_version`                        | Barcode dataset used for the sample's demix (as in `demix_status.tsv`; `custom` in the example).                                                                                                                                                                      |
| `primer_trimmed_pct`                     | % of mapped reads from which `ivar trim` removed primer sequence (from the ivar log). `NA` without primer trimming. Informational unless `--qc_min_primer_trimmed_pct` is set.                                                                                        |
| `primer_scheme_check`                    | `ok`, `mismatch`, `uninformative` or `NA`: does the declared primer scheme match where the reads start ([usage.md](usage.md#primer-scheme-check)). `uninformative` means the reads do not start at primers (for example tagmentation libraries) and is not a warning. |
| `primer_scheme_best`                     | The scheme(s) that fit the reads best, joined by `/` (`midnight/midnightv2` have the same coordinates); `NA` when none stands out.                                                                                                                                    |
| `primer_scheme_enrichment`               | Enrichment of read starts at the declared scheme's primers, as a multiple of chance (1 decimal).                                                                                                                                                                      |

QC rules (thresholds are parameters; defaults shown):

| Sample type | Status                   | Rule                                                                                                               | Token                    |
| ----------- | ------------------------ | ------------------------------------------------------------------------------------------------------------------ | ------------------------ |
| sample      | FAIL                     | no mapped reads or all depths 0                                                                                    | `no_reads`               |
| sample      | FAIL                     | depths unreadable or empty                                                                                         | `no_depths`              |
| sample      | FAIL                     | % genome at `--qc_min_depth` (10x) below `--qc_fail_below_pct` (50)                                                | `low_coverage`           |
| sample      | WARN                     | % genome at `--qc_min_depth` below `--qc_pass_min_pct` (80)                                                        | `low_coverage`           |
| sample      | WARN (if otherwise PASS) | more than `--qc_max_dropped_amplicons` (1) amplicons dropped                                                       | `amplicon_dropout`       |
| sample      | WARN (if otherwise PASS) | the declared primer scheme does not match the reads (`primer_scheme_check` is `mismatch`)                          | `primer_scheme_mismatch` |
| sample      | WARN (if otherwise PASS) | `primer_trimmed_pct` below `--qc_min_primer_trimmed_pct` (off by default)                                          | `low_primer_trimming`    |
| sample      | PASS                     | otherwise                                                                                                          |                          |
| negative    | FAIL                     | mapped reads above `--ntc_max_mapped_reads` (1000) or % genome at >= 1x above `--ntc_max_genome_coverage_pct` (10) | `ntc_reads`              |
| negative    | PASS                     | otherwise (never `low_coverage`/`no_reads`)                                                                        |                          |
| negative    | NA                       | no usable depths and no flagstat                                                                                   | `no_depths`              |

`NEG_E` above (5 mapped reads) is a PASS; `NEG_CONTAM` (a contaminated control with 5998 mapped reads) is a FAIL with
`ntc_reads`, which also raises a run-level error in the report. A negative control that is `PASS` here can still be
`skipped_low_coverage` in `demix_status.tsv` (`NEG_E`): that is the expected outcome for a clean control, not a failure.

### `summary/amplicon_coverage.tsv`

Long format, one row per sample x amplicon (7 samples x 29 midnight amplicons = 203 rows here), negative controls
included. Only written when the primer scheme is known (always in `full` with the `amplicon` protocol; in `demix` only
with `--primer_scheme` / `--primer_bed`). **Never written for `--protocol capture` or `--skip_primer_trimming`**: there
are no amplicons (and then `amplicons_total` / `amplicons_dropped` in `qc_summary.tsv` are `NA` and no `amplicon_dropout`
reason is raised).

```
sample    amplicon      start  end   mean_depth  min_depth  pct_ge_min_depth  dropped
MIX_A     nCoV-2019_1   31     1205  89.02       0          97.96             false
MIX_A     nCoV-2019_2   1101   2266  81.84       43         100.00            false
```

| Column             | Meaning                                                                                                                 |
| ------------------ | ----------------------------------------------------------------------------------------------------------------------- |
| `sample`           | Sample ID.                                                                                                              |
| `amplicon`         | Amplicon label, the shared primer key (`nCoV-2019_12`); alternates (`_alt0`) are folded into their amplicon.            |
| `start`, `end`     | 1-based inclusive genome span from the left primer start to the right primer end (outer coordinates, primers included). |
| `mean_depth`       | Mean depth over the span; positions absent from the depths file count as 0.                                             |
| `min_depth`        | Minimum depth over the span.                                                                                            |
| `pct_ge_min_depth` | % of the span with depth >= `--amplicon_min_depth`.                                                                     |
| `dropped`          | `true` if `mean_depth < --amplicon_min_depth` (default 10), else `false`.                                               |

The same table for one sample is in `samples/<s>/coverage/<s>.amplicon_coverage.tsv`.

### `summary/demix_status.tsv`

One row per sample, **whatever happened to it**. Sorted by `sample`. Real content:

```
sample         status                reason                                          genome_coverage_pct  barcode_version  barcode_version_iso
LOWCOV_D       skipped_low_coverage  less than 60% of the genome covered at >= 10x   0.0                  custom
MIX_A          ok                                                                    99.3345              custom
NEG_CONTAM     ok                                                                    94.5457              custom
NEG_E          skipped_low_coverage  less than 60% of the genome covered at >= 10x   0.0                  custom
Site_1_week_2  ok                                                                    93.0576              custom
```

| Column                | Meaning                                                                                                                                                                                                                                                                                                                                                |
| --------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| `sample`              | Sample ID.                                                                                                                                                                                                                                                                                                                                             |
| `status`              | `ok`: demix ran and a `.demix` file exists. `failed`: `freyja demix` errored or the depths file was unreadable (the reason holds the tail of the error); no `.demix`. `skipped_low_coverage`: the coverage gate stopped the sample (below). `missing`: a `demix`/`aggregate` run found the sample in earlier results but could not use it (see below). |
| `reason`              | Empty for `ok`; otherwise a human-readable reason.                                                                                                                                                                                                                                                                                                     |
| `genome_coverage_pct` | % of the genome at depth >= `--freyja_covcut` (default 10x), the value compared with `--min_genome_coverage_for_demix`. Empty for `missing`.                                                                                                                                                                                                           |
| `barcode_version`     | Barcode dataset used: the verbatim `last_barcode_update.txt` (`MM_DD_YYYY-HH-MM`, e.g. `09_17_2026-00-46`) or `custom` for a bare `.csv`/`.feather` file.                                                                                                                                                                                              |
| `barcode_version_iso` | The same version as `YYYY-MM-DDTHH-MM`; empty for `custom`.                                                                                                                                                                                                                                                                                            |

The coverage gate exists because Freyja reports a confident-looking answer even at very low depth (the 0.6x sample
`LOWCOV_D` demixes to a single lineage with Freyja `coverage` 0.0); only the gate keeps such a result out. It is on by
default (`--min_genome_coverage_for_demix 60`). A sample with essentially no mapped reads (`NEG_E`) is routed to
`skipped_low_coverage` as well, because `freyja demix` would crash on its header-only variants file.

`missing` rows appear in `demix`/`aggregate` runs for samples found in `--previous_results` that cannot be used: a
legacy zero-byte `.demix` file (written by the old Snakemake pipeline when demix failed; reason
`empty demix file (failed in the old pipeline)`), or required files that are absent (for example
`missing demix (new layout)`). They are listed but never aggregated.

### `summary/freyja_aggregated.tsv`

Freyja's own aggregate table for every sample with status `ok` (negative controls excluded unless
`--include_negatives_in_aggregate`; see [`freyja_aggregated_controls.tsv`](#summaryfreyja_aggregated_controlstsv)), with sample labels cleaned by wasteline (exactly the sample IDs; no `_` -> `-`
rewriting) and rows sorted by sample. Real content:

```
	summarized	lineages	abundances	resid	coverage	pathogen
MIX_A	[('Omicron', 0.7121209999278063), ('Delta', 0.287879000132854)]	BA.2 AY.4	0.71212100 0.28787900	3.324458426135655	99.33451493161222	SARS-CoV-2
MIX_B	[('Omicron', 0.5222220001591816), ('XBB.1.5* [Omicron (XBB.1.5.X)]', 0.2777780000352889), ('KP.3* (KP.3.X)', 0.20000000003738969)]	BA.1 XBB.1.5 KP.3 BA.2	0.50000000 0.27777800 0.20000000 0.02222200	4.1167517193691125	99.36795639233522	SARS-CoV-2
PURE_C	[('BA.2.86* (BA.2.86.X)', 0.9999999999939578)]	JN.1	1.00000000	0.09315044335376077	99.36795639233522	SARS-CoV-2
```

The first column (empty header) is the sample label.

| Column       | How to read it                                                                                                                                                                                                                                                                                     |
| ------------ | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `summarized` | Freyja's **lineage groups**: a Python-style list of `(group, fraction)` tuples. Several lineages are collapsed into a broader group (WHO-style names such as `Omicron`, `Delta`, or `XBB.1.5* [Omicron (XBB.1.5.X)]`). Good for an overview and for plots, bad for tracking a specific sublineage. |
| `lineages`   | Space-separated individual lineages (`BA.2 AY.4`), in order of decreasing abundance.                                                                                                                                                                                                               |
| `abundances` | Space-separated fractions (0-1), **parallel to `lineages`**: the n-th abundance belongs to the n-th lineage.                                                                                                                                                                                       |
| `resid`      | Freyja's residual of the fit (lower is a better fit; it is the goodness-of-fit of the deconvolution, not an abundance). Here about 0.09 for the pure sample and 3-4 for mixtures.                                                                                                                  |
| `coverage`   | Freyja's % of the genome covered (equal to `genome_coverage_pct` in `demix_status.tsv`).                                                                                                                                                                                                           |
| `pathogen`   | Freyja `--pathogen` the data belongs to (`SARS-CoV-2`).                                                                                                                                                                                                                                            |

Important when comparing against a truth set or another tool: `summarized` and `lineages` are **different groupings of
the same estimate**. Freyja collapses lineages in `summarized`, for example **JN.1 appears as `BA.2.86* (BA.2.86.X)`**
(sample `PURE_C` above: `lineages = JN.1`, `summarized = BA.2.86*`). Likewise the `BA.2` and `BA.1`/`BA.5` calls fall
into `Omicron`, while `XBB.1.5` and `KP.3` get their own named groups. Compare on `lineages`/`abundances` when you need
exact lineages. Small spurious minor lineages (for instance `BA.2` at 0.02-0.03 in the BA.1/BA.5 mixtures) are normal;
the test suite allows lineages absent from the truth up to 0.05.

The HTML report offers both views with a toggle ("Summarized groups" / "All lineages").

### `summary/freyja_aggregated_controls.tsv`

The negative controls (`sample_type=negative`) that demixed `ok`, in exactly the format of `freyja_aggregated.tsv`.
They are demixed so that contamination is visible, but kept out of the main aggregate and the plots (decision D8), so
this is where to look at what a control contains. Real content (the contaminated control `NEG_CONTAM`; the clean
`NEG_E` was `skipped_low_coverage` and has no row):

```
	summarized	lineages	abundances	resid	coverage	pathogen
NEG_CONTAM	[('Omicron', 0.9999999999495606)]	BA.2	1.00000000	0.0427725076591042	94.54569775607798	SARS-CoV-2
```

Columns are as in [`freyja_aggregated.tsv`](#summaryfreyja_aggregatedtsv). The file is not written when no control
demixed `ok`, and not at all with `--include_negatives_in_aggregate` (the controls are then in the main aggregate). It
feeds the report's "Negative controls" section (`negative_controls` in [report_data.md](report_data.md)).

### `summary/freyja_plot.pdf`, `summary/freyja_plot_times.pdf`

`freyja plot` stacked-bar plots of the aggregated results (`--freyja_plot_mincov` hides samples below that coverage).
`freyja_plot_times.pdf` is the time-aware plot (`freyja plot --times`), drawn for the samples that have a
`collection_date` (`--freyja_plot_interval`, `--freyja_plot_windowsize`). If a plot cannot be produced the run still
succeeds and a note `freyja_plot.failed.txt` / `freyja_plot_times.failed.txt` with the reason is published instead of
a PDF. Background: `freyja plot` is broken in the pinned Freyja container by a pandas 3 incompatibility, so wasteline
retries it through a small workaround (PLAN decision D6); the HTML report draws its own plots and never depends on it.

### `summary/consensus_stats.tsv` and `summary/consensus_genomes.fasta`

One consensus per non-negative sample (negative controls get none, so they have no row). `--skip_consensus` disables
both files.

```
sample         length  n_count  n_pct   ambiguous_count  status  reason
LOWCOV_D       29903   29903    100.00  0                empty
MIX_A          29903   413      1.38    59               ok
PURE_C         29903   189      0.63    0                ok
Site_1_week_2  29903   2076     6.94    27               ok
```

| Column            | Meaning                                                                                                                                                                                                                                                               |
| ----------------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `sample`          | Sample ID (also the FASTA header).                                                                                                                                                                                                                                    |
| `length`          | Sequence length (equals the reference length, 29903, as variants are SNVs; deletions are written as `-`).                                                                                                                                                             |
| `n_count`         | Number of `N` and `-` characters.                                                                                                                                                                                                                                     |
| `n_pct`           | `n_count` as % of the reference length.                                                                                                                                                                                                                               |
| `ambiguous_count` | Number of other non-ACGT characters, i.e. IUPAC ambiguity codes (`R Y K M S W`...).                                                                                                                                                                                   |
| `status`          | `ok`, `empty`, or `failed` (see [consensus conventions](#consensus-fasta-conventions)). For `failed`, the numeric fields are `NA`.                                                                                                                                    |
| `reason`          | Last column. `NA` for a clean `ok`; otherwise a one-line text: why a sample is `empty` or `failed`, or (for an `ok` sample) which variants were not applied, for example `1 variant(s) not applied: skipped variant at NC_045512.2:22942 (overlaps another variant)`. |

### `summary/mutations_of_interest.tsv`

Only with `--mutations_of_interest moi.csv` (a CSV with columns `name,mutation`). One row per sample x line of the input
file, in input order, negative controls included. Real content (`moi.csv`: `S_N501Y,S:N501Y`, `S_L452R,S:L452R`,
`S_K417N,A22813T`, `ORF1b_P314L,ORF1ab:P4715L`):

```
sample  name         mutation       ref_depth  alt_depth  depth  alt_freq
MIX_A   S_N501Y      S:N501Y        11         41         52     0.788462
MIX_A   S_L452R      S:L452R        31         16         47     0.340426
MIX_A   S_K417N      A22813T        70         0          70     0
LOWCOV_D S_N501Y     S:N501Y        0          0          0      0
```

The last column, `compound`, is `NA` for ordinary single-variant changes. When adjacent variants fall in one codon,
`bcftools csq` reports one joint consequence on the first variant and a back-reference on the others; wasteline copies
the joint consequence onto every variant involved and lists the positions in `compound` (for example `23063+23064` for
two SNVs giving `N501S`).

| Column                   | Meaning                                                                                                                                                                                |
| ------------------------ | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `name`, `mutation`       | As in the input file. Syntax: SNV `A23063T`; deletion by its first deleted base `C75-`; insertion by anchor `A14+C`; amino acid `S:N501Y`, `ORF1ab:P4715L` (resolved from the GFF3).   |
| `ref_depth`, `alt_depth` | Read depth supporting the reference and the alternate allele (from the Freyja/ivar variants TSV).                                                                                      |
| `depth`                  | Total depth at the site (for amino-acid mutations, the minimum over the three codon positions).                                                                                        |
| `alt_freq`               | Alternate allele frequency, 0-1. A mutation not in the variants file has `alt_depth 0`, `alt_freq 0` and `depth` taken from the depths file; `depth 0` means no coverage, not absence. |

The report shows these as a heatmap; cells with depth below `--consensus_min_depth` are drawn with a dashed outline.
The per-sample piece is `samples/<s>/variants/<s>.mutations_of_interest.tsv`.

### `summary/freyja_boot_summary.tsv`

Only with `--freyja_boot`. Bootstrap intervals from `freyja boot` (`--freyja_boot_n` replicates, default 100, seed
`--freyja_boot_seed`), one row per sample x lineage, sorted by sample. Real content (trimmed):

```
sample      lineage  abundance  lo        hi
MIX_A       BA.2     0.712121   0.670435  0.727273
MIX_A       AY.4     0.287879   0.269521  0.304565
MIX_A       B.1      0.000000   0.000000  0.019677
MIX_A       BA.5     0.000000   0.000000  0.034250
MIX_B       BA.1     0.500000   0.474049  0.528547
NEG_CONTAM  BA.2     1.000000   1.000000  1.000000
PURE_C      JN.1     1.000000   1.000000  1.000000
```

| Column              | Meaning                                                                                                      |
| ------------------- | ------------------------------------------------------------------------------------------------------------ |
| `sample`, `lineage` | Sample ID and an individual lineage (the `lineages` row of the `.demix`, not the `summarized` groups).       |
| `abundance`         | The point estimate of the sample's own `.demix` file (0 for a lineage that only appears in some bootstraps). |
| `lo`, `hi`          | The 2.5 and 97.5 percentiles of the bootstrap replicates (a 95% interval, fractions 0-1).                    |

A lineage is listed when it has a point estimate or a non-zero upper bound. The intervals come from resampling the
sample's own data; they do not include the effect of choosing another barcode dataset. Every sample that demixed `ok`
is bootstrapped, negative controls included; skipped or failed samples have no rows. A failed bootstrap never fails the
run: the sample gets no rows and a `failed` `<s>.boot_status.json` with the reason. The per-sample files
`samples/<s>/freyja/<s>.boot_summary.tsv` (the sample's rows, no header), `<s>.boot_lineages.csv` and
`<s>.boot_summarized.csv` (Freyja's quantile tables: rows `0.025, 0.05, 0.25, 0.5, 0.75, 0.95, 0.975`, one column per
lineage or summarized group) hold the full distributions. The report shows `lo`-`hi` in the tooltips of the
lineage view.

### `summary/barcode_compare.tsv` and `summary/barcode_compare_status.tsv`

Only with `--compare_barcodes <dataset>`. Every sample that demixes is demixed a second time with the second dataset
(`samples/<s>/freyja/<s>.compare.demix`, `.compare.demix_status.json`, `.compare.depth_stats.json`; same variants and
depths, same demix parameters). `barcode_compare.tsv` puts the two estimates side by side, one row per sample x lineage
that either dataset reports (0 where one does not), sorted by sample. A = `--barcodes`, B = `--compare_barcodes`. Real
content (trimmed; A is the test barcode CSV, B the barcodes bundled in the container):

```
sample  lineage  abundance_a  abundance_b  barcode_version_a  barcode_version_b
MIX_A   BA.2     0.712121     0.692308     custom             09_17_2026-00-46
MIX_A   AY.4     0.287879     0.176471     custom             09_17_2026-00-46
MIX_A   AY.1     0.000000     0.071429     custom             09_17_2026-00-46
MIX_B   KP.3     0.200000     0.024492     custom             09_17_2026-00-46
MIX_B   KP.2     0.000000     0.126662     custom             09_17_2026-00-46
PURE_C  JN.1     1.000000     1.000000     custom             09_17_2026-00-46
```

The abundances are those of the `lineages`/`abundances` rows (individual lineages). Different barcode sets can assign
the same signal to different lineages, so read a change in the context of the whole sample rather than row by row: in
the example `MIX_A`'s 28.8% `AY.4` of dataset A becomes 17.6% `AY.4` plus `AY.1`, `AY.43.8`, ... in dataset B. Negative
controls are included.

`barcode_compare_status.tsv` records the outcome of the second demix per sample (`sample, status, reason,
barcode_version` with the barcode version of dataset B), with the same statuses as `demix_status.tsv` (a sample the
coverage gate skipped is `skipped_low_coverage` again):

```
sample         status                reason                                         barcode_version
LOWCOV_D       skipped_low_coverage  less than 60% of the genome covered at >= 10x  09_17_2026-00-46
MIX_A          ok                                                                   09_17_2026-00-46
NEG_E          skipped_low_coverage  less than 60% of the genome covered at >= 10x  09_17_2026-00-46
```

A second demix that fails is never fatal: it is `failed` here and in `.compare.demix_status.json`, the sample has no
rows in `barcode_compare.tsv`, and the report raises the warning `barcode_compare_failed`. Dataset B is resolved
independently of dataset A and its provenance is recorded as `freyja_data_compare` (see
[run_info.json](#run_infojson) and `pipeline_info/freyja_data_compare.json`).

### `summary/annotated_variants.tsv`

Only with `--annotate_variants` (and the consensus chain: not with `--skip_consensus`). The amino-acid consequences of
the variants that were **applied to the consensus** (the calls that passed the consensus filters, i.e. the
`variants/<s>.consensus_variants.vcf.gz` of each sample), from `bcftools csq` with the reference GFF3; a port of the
old pipeline's `amino_acid_consequences` rule. One row per variant, samples in order, negative controls excluded (no
consensus). The per-sample piece is `samples/<s>/variants/<s>.annotated.tsv`. Real content (trimmed):

```
reference    sample_id  gene    nt_pos  ref  alt  effect       aa_bcsq      aa_standard  vaf       depth
NC_045512.2  MIX_A      ORF1ab  670     T    G    missense     135S>135R    S135R        0.692308  65
NC_045512.2  MIX_A      ORF1ab  3037    C    T    synonymous   924F         F924F       1         82
NC_045512.2  MIX_A      S       21618   C    T    missense     19T>19I      T19I         0.733333  105
NC_045512.2  MIX_A      NA      27383   A    T    intergenic   NA           NA           0.589744  39
```

| Column                   | Meaning                                                                                                                                                                                                                                                                    |
| ------------------------ | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `reference`, `sample_id` | Contig and sample ID.                                                                                                                                                                                                                                                      |
| `gene`                   | Gene name from the GFF3 (`ORF1ab`, `S`, `N`, ...); `NA` when there is no coding consequence.                                                                                                                                                                               |
| `nt_pos`, `ref`, `alt`   | 1-based genome position and alleles.                                                                                                                                                                                                                                       |
| `effect`                 | First `bcftools csq` consequence: `missense`, `synonymous`, `stop_gained`, ... or `intergenic` (no coding consequence, `gene` and the amino-acid columns `NA`). The real run has 271 `missense`, 75 `synonymous`, 24 `intergenic` and 2 `stop_gained` rows (373 in total). |
| `aa_bcsq`, `aa_standard` | The change in bcftools notation (`135S>135R`) and in the usual notation (`S135R`).                                                                                                                                                                                         |
| `vaf`, `depth`           | `ALT_FREQ` and total depth of the call in the variants TSV.                                                                                                                                                                                                                |

Caveat (found in this run): when several adjacent SNVs fall in one codon, `bcftools csq` reports the combined
consequence on the first one and the others only as a continuation, which the module writes as `intergenic` with `NA`
gene. All 24 `intergenic` rows of the real run are of this kind (for example positions 27383 and 27384 next to the
`ORF6` `D61L` call at 27382, and 28882 between the two `N` calls at 28881 and 28883), so an `intergenic` row next to a
coding call in the same codon is not a non-coding variant <!-- UNVERIFIED: inferred from the positions of the rows and from bcftools csq reporting multi-base codon changes once; the raw BCSQ field was not inspected -->.

### `summary/barcode_version.json`

The small record of the Freyja barcode dataset, published next to `freyja_aggregated.tsv` so that the table is never
separated from the dataset that produced it (written in `full` and `demix`; an `aggregate` run demixes nothing and
writes none). The same record is embedded in every `<s>.demix_status.json`.

```json
{
  "barcode_version": "custom",
  "barcode_version_iso": null,
  "barcodes_sha256": "0c3cbd7a1e6a7e6a8ce3e36ca22f4060c483a65e698fff0bec8e3226639509c7",
  "source": "custom",
  "freyja_version": "2.0.5"
}
```

For a `latest` or `container` run `barcode_version` is the stamp (for example `09_17_2026-00-46`) and
`barcode_version_iso` is `2026-09-17T00-46`. Field meanings: see [`freyja_data`](#the-freyja_data-block).

### `summary/dropout_impact.tsv` and `summary/dropout_impact_summary.tsv`

Which dropped amplicons matter for the lineages that were detected ([usage.md](usage.md#dropout-impact)). Written in `full`
mode (amplicon protocol) and in `demix` mode with an explicit primer scheme; not in other modes, and skipped with one warning
when the Freyja dataset has no readable mutation lists (a `.feather`-only custom `--barcodes`). `dropout_impact.tsv` has one row per
**dropped** amplicon of each sample (the header is always written):

```
sample         amplicon      start  end    mean_depth  n_informative_sites  n_sites_in_detected_lineages  affects_detected_lineage  affected_lineages
Site_1_week_2  nCoV-2019_7   6284   7401   4.66        5                    0                             false                     NA
Site_1_week_2  nCoV-2019_22  21533  22612  5.22        28                   10                            true                      XBB.1.5:9,BA.2:3,BA.5:3
```

| Column                                             | Meaning                                                                                                                                               |
| -------------------------------------------------- | ----------------------------------------------------------------------------------------------------------------------------------------------------- |
| `sample`, `amplicon`, `start`, `end`, `mean_depth` | As in `amplicon_coverage.tsv`.                                                                                                                        |
| `n_informative_sites`                              | Single-nucleotide sites inside the amplicon that some, but not nearly all, lineages of the dataset carry (fixed mutations cannot separate lineages).  |
| `n_sites_in_detected_lineages`                     | Of those, the number carried by lineages detected in the sample (abundance >= `--dropout_min_abundance`); `NA` for a sample without an aggregate row. |
| `affects_detected_lineage`                         | `true` when a detected lineage carries at least one informative site in the amplicon.                                                                 |
| `affected_lineages`                                | Up to five detected lineages with their site counts (`name:count`), most affected first; `NA` if none.                                                |

`dropout_impact_summary.tsv` has one row per sample of the amplicon table (zeros when nothing dropped):

```
sample         n_dropped_amplicons  informative_sites_total  informative_sites_in_dropped  pct_informative_sites_lost  n_dropped_affecting_detected_lineages
MIX_A          0                    244                      0                             0.00                        0
NEG_E          29                   244                      244                           100.00                      NA
Site_1_week_2  2                    244                      33                            13.52                       1
```

(The 244 informative sites are those of the small test barcode file; with a real Freyja dataset there are thousands.) Samples
that were not demixed (skipped, negative controls) have `NA` for the lineage columns.

### `summary/unexplained_mutations.tsv`

Frequent single-nucleotide variants that the demixed lineage mixture does not predict
([usage.md](usage.md#unexplained-mutations)), for all non-negative samples of modes `full` and `demix`. The header is always
written; on the synthetic samples at the default thresholds the table is empty. The per-sample tables are in
`samples/<s>/variants/<s>.unexplained_mutations.tsv` (without `n_samples` and `recurrent`, which are added when the tables are
combined). A real row (a spiked variant in the synthetic `PURE_C` test fixture):

```
sample  position  ref  alt  observed_freq  expected_freq  excess  depth  alt_depth  gene  aa_change  effect       n_samples  recurrent  p_value
PURE_C  22102     A    T    0.3            0              0.3     120    36         S     E180D      missense     1          false      2.39e-42
PURE_C  28802     A    T    0.18           0              0.18    100    18         N     R177*      stop_gained  1          false      1.41e-17
```

| Column                                | Meaning                                                                                                                                                 |
| ------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `observed_freq`, `depth`, `alt_depth` | From the ivar variants (`ALT_FREQ`, `TOTAL_DP`, `ALT_DP`).                                                                                              |
| `expected_freq`                       | What the demixed lineages predict: the sum over detected lineages of abundance times "carries the mutation".                                            |
| `excess`                              | `observed_freq - expected_freq`.                                                                                                                        |
| `gene`, `aa_change`, `effect`         | Annotation from the GFF3 (`NA` / `intergenic` outside coding sequence).                                                                                 |
| `n_samples`, `recurrent`              | In how many samples the same position and alternate base are unexplained; `true` at `--unexplained_min_recurrence` samples or more.                     |
| `p_value`                             | Binomial p-value of the alternate-read count against the expected frequency (or a small background rate), see [parameters.md](parameters.md#qc-checks). |

### `summary/freyja_timeseries.tsv` and `summary/freyja_timeseries_lineages.tsv`

Long tables of abundance by site, collection date and lineage group, from the aggregate and the samplesheet
([usage.md](usage.md#time-series-growth-rates-and-the-dashboard)). Written in all modes when there are at least two distinct
collection dates among the non-control samples. `freyja_timeseries.tsv` uses Freyja's summarized groups, the `_lineages` table
the individual lineages. Real rows (the `tests/data_timeseries` series, aggregate mode, with viral loads):

```
site   collection_date  lineage_group  abundance  abundance_load_weighted  n_samples  mean_viral_load  total_viral_load
SiteA  2026-01-05       BA.2           0.800633   0.800633                 1          175675.00        175675.00
SiteA  2026-01-05       JN.1           0.105472   0.105472                 1          175675.00        175675.00
SiteA  2026-01-12       BA.2           0.689777   0.701996                 2          576892.50        1153785.00
SiteA  2026-01-12       XBB.1.5        0.199020   0.183810                 2          576892.50        1153785.00
```

`abundance` is the mean over the samples of that site and date (a group a sample lacks counts as 0);
`abundance_load_weighted` weights each sample by its `viral_load` (`NA` when no sample has one); `n_samples` counts every
contributing sample, `mean_viral_load`/`total_viral_load` only those with a load. Negative controls, undated samples, samples
below `--freyja_plot_mincov` coverage and samples missing from the samplesheet are left out. The header is always written.

### `summary/freyja_growth_rates.csv`, `freyja_growth_rates.failed.txt`, `freyja_growth_inputs.json`

Only with `--freyja_growth`: Freyja's relative growth rates for the last `--freyja_growth_days` days. Exactly one of the CSV and
the `.failed.txt` note is written (too few dates or too short a span, or a Freyja error; never a failed run). Real output
(`tests/data_timeseries`, where XBB.1.5 rises and BA.2 declines):

```
Lineage,Estimated Advantage,Bootstrap 95% interval
XBB.1.5,60.3%,"[44.57% , 50.42%]"
BA.1,-4.8%,"[-30.00% , 5.73%]"
BA.2,-35.1%,"[-50.50% , -36.83%]"
JN.1,-1.0%,"[-2.64% , 1.87%]"
```

The advantage is per serial interval (`--freyja_growth_serial_interval`, 5.5 days) relative to the other lineages. Freyja's bootstrap
interval does not always contain its point estimate (as above), so read the sign and rough size, not the digits.
`freyja_growth_inputs.json` records what went into the fit (`n_samples`, `n_dates`, `n_dates_window`, `first_date`, `last_date`,
`span_days`, `n_with_viral_load`, `window_days`).

### `summary/freyja_dashboard.html`

Only with `--freyja_dash`: Freyja's interactive dashboard (Plotly is embedded; the Bootstrap styling and a logo are fetched from the
internet when the page is opened). If it cannot be built, `summary/freyja_dashboard.failed.txt` says why. The report names the
file but does not embed it.

### `summary/samplesheet.csv`

The samplesheet in effect, copied to the results so that a later `demix`/`aggregate`/`report` run on this directory
finds `sample_type`, `collection_date` and `site`. In `full` mode it is your input sheet verbatim (note it holds the
absolute FASTQ paths you supplied); for runs that discover samples from `--previous_results` it is the discovered
sheet. Its sha256 is recorded in `run_info.json`.

## Per-sample files

Under `samples/<sample>/`. Which files exist depends on the mode and on the sample (see the mode table above).

| File                                                                                  | Content                                                                                                                                                                                                                                                                                                                                                                                  |
| ------------------------------------------------------------------------------------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `fastp/<s>.fastp.json`, `.html`, `.log`                                               | fastp report and log. The JSON feeds `raw_reads`/`trimmed_reads`. Trimmed FASTQs (`<s>.fastp.fastq.gz`) are published only with `--save_trimmed_reads` <!-- UNVERIFIED: file name taken from the nf-core fastp module output pattern; not exercised in this run -->.                                                                                                                     |
| `bam/<s>.primertrimmed.bam`, `.bai`                                                   | Reads aligned with `bwa mem` (or minibwa with `--aligner minibwa`) to the reference, mapped reads only, primers soft-clipped with `ivar trim`, sorted and indexed (`--save_bam`, default true). For a sample with no mapped reads (a clean control) `ivar trim` is skipped and the (empty) BAM is carried on.                                                                            |
| `bam/<s>.mapped.bam`, `.bai`                                                          | Only with `--protocol capture` or `--skip_primer_trimming`: the mapped reads as aligned (nothing trimmed), sorted and indexed. Not published when duplicates are marked (the marked BAM below replaces it).                                                                                                                                                                              |
| `bam/<s>.<primertrimmed\|mapped>.markdup.bam`, `.bai`                                 | Only with duplicate marking (`--duplicates mark`, or `auto` with `--protocol capture`): the BAM Freyja worked on, PCR duplicates flagged (0x400), none removed. Published unless `--save_bam false`.                                                                                                                                                                                     |
| `bam/<s>.markdup_stats.txt`                                                           | Only with duplicate marking: `samtools markdup -s` statistics (`DUPLICATE PAIR`, `EXAMINED`, ...).                                                                                                                                                                                                                                                                                       |
| `bam/<s>.primertrimmed.cram`, `.crai`                                                 | Only with `--bam_format cram`: a reference-based CRAM (much smaller) instead of the BAM; read it with `samtools view -T <reference.fa>`. Freyja and QC always use the BAM inside the run.                                                                                                                                                                                                |
| `bam/<s>.ivar_trimmed.ivar.log`                                                       | ivar trim log (its summary: reads trimmed, with and without primers; BED warnings). `primer_trimmed_pct` in `qc_summary.tsv` is computed from it.                                                                                                                                                                                                                                        |
| `freyja/<s>.variants.tsv`                                                             | `freyja variants` output (ivar TSV, 20 columns: `REGION POS REF ALT REF_DP REF_RV REF_QUAL ALT_DP ALT_RV ALT_QUAL ALT_FREQ TOTAL_DP PVAL PASS GFF_FEATURE REF_CODON REF_AA ALT_CODON ALT_AA POS_AA`). With the default `--freyja_varthresh 0.0` it lists **every** alternate allele seen, including 1-read noise; filter on `ALT_FREQ` yourself. Header-only for a sample without reads. |
| `freyja/<s>.depths`                                                                   | Four tab-separated columns, no header: contig, 1-based position, reference base, depth (Q>=20 bases), one row per genome position. All coverage QC and the consensus mask are computed from it.                                                                                                                                                                                          |
| `freyja/<s>.demix`                                                                    | Freyja demix result for the sample (key/value rows `summarized`, `lineages`, `abundances`, `resid`, `coverage`, `pathogen`; first line holds the label). Present only for status `ok`.                                                                                                                                                                                                   |
| `freyja/<s>.demix_status.json`                                                        | `{sample, status, reason, genome_coverage_pct, covcut, min_genome_coverage, barcode_version, barcode_version_iso}`; the source of the `demix_status.tsv` row.                                                                                                                                                                                                                            |
| `freyja/<s>.depth_stats.json`                                                         | `{sample, positions, positions_ge_covcut, mean_depth, pct_genome_ge_covcut, covcut, min_genome_coverage, passes_gate}`: the coverage-gate computation.                                                                                                                                                                                                                                   |
| `freyja/<s>.compare.demix`, `.compare.demix_status.json`, `.compare.depth_stats.json` | Only with `--compare_barcodes`: the same three files for the second barcode dataset (the status JSON carries dataset B's barcode version). `.compare.demix` only for status `ok`.                                                                                                                                                                                                        |
| `freyja/<s>.boot_summary.tsv`                                                         | Only with `--freyja_boot`: the sample's rows of `freyja_boot_summary.tsv` (data rows, no header).                                                                                                                                                                                                                                                                                        |
| `freyja/<s>.boot_lineages.csv`, `.boot_summarized.csv`                                | Only with `--freyja_boot`: Freyja's bootstrap quantile tables (rows `0.025` ... `0.975`; one column per lineage / summarized group).                                                                                                                                                                                                                                                     |
| `freyja/<s>.boot_status.json`                                                         | Only with `--freyja_boot`: `{sample, status, reason}`; `failed` with the reason when the bootstrap failed (never fatal).                                                                                                                                                                                                                                                                 |
| `coverage/<s>.primertrimmed.flagstat`                                                 | `samtools flagstat` of the final BAM (mode `full`); named `<s>.mapped.flagstat` when nothing was primer-trimmed.                                                                                                                                                                                                                                                                         |
| `coverage/<s>.amplicon_coverage.tsv`                                                  | The sample's rows of `amplicon_coverage.tsv`.                                                                                                                                                                                                                                                                                                                                            |
| `coverage/<s>.primer_scheme_check.tsv`, `.primer_scheme_summary.tsv`                  | Mode `full`, amplicon protocol with trimming: the primer scheme check (one row per bundled scheme with hits, expected hits and enrichment; and the one-row verdict `ok`/`mismatch`/`uninformative`/`NA`).                                                                                                                                                                                |
| `coverage/<s>.qc_row.tsv`                                                             | The sample's row of `qc_summary.tsv`.                                                                                                                                                                                                                                                                                                                                                    |
| `variants/<s>.unexplained_mutations.tsv`                                              | The sample's unexplained mutations (see `summary/unexplained_mutations.tsv`); header only when there are none. Not for negative controls.                                                                                                                                                                                                                                                |
| `variants/<s>.ivar.vcf`                                                               | The variants TSV converted to VCF (all alleles; `FILTER` is ivar's p-value filter).                                                                                                                                                                                                                                                                                                      |
| `variants/<s>.consensus_variants.vcf.gz`, `.tbi`                                      | The calls that passed the consensus filters (depth, ALT quality, frequency) with genotype `1/1` (alternate) or `0/1` (ambiguous, becomes an IUPAC code). Not written for negative controls.                                                                                                                                                                                              |
| `variants/<s>.mask.bed`                                                               | 0-based BED of positions masked `N` (depth below `--consensus_min_depth`, confident variant sites excepted).                                                                                                                                                                                                                                                                             |
| `variants/<s>.mutations_of_interest.tsv`                                              | Only with `--mutations_of_interest`: per-sample rows of `mutations_of_interest.tsv`.                                                                                                                                                                                                                                                                                                     |
| `variants/<s>.annotated.tsv`                                                          | Only with `--annotate_variants`: the sample's rows of `annotated_variants.tsv` (with header).                                                                                                                                                                                                                                                                                            |
| `consensus/<s>.consensus.fa`                                                          | The consensus genome (header `><s>`, 60 columns).                                                                                                                                                                                                                                                                                                                                        |
| `consensus/<s>.consensus_stats.tsv`                                                   | The sample's row of `consensus_stats.tsv`.                                                                                                                                                                                                                                                                                                                                               |
| `consensus/<s>.consensus_status.json`                                                 | `{sample, status, reason, warnings}`; `warnings` is an array (variants that were not applied, at most 20), e.g. `{"status": "empty", "reason": "no position reaches the minimum consensus depth (no usable coverage)", "warnings": []}`.                                                                                                                                                 |
| `consensus/<s>.consensus.log`                                                         | `bcftools consensus` log (applied variants, skipped overlapping sites).                                                                                                                                                                                                                                                                                                                  |

## Consensus FASTA conventions

`summary/consensus_genomes.fasta` (all samples) and `samples/<s>/consensus/<s>.consensus.fa` (one sample):

- **Header** is the bare sample ID: `>MIX_A`.
- **Reference-length output.** SNVs keep the length; the sequence is 29903 bases for NC_045512.2. Every non-negative
  sample has a record; a negative control has none.
- **`N` masking.** A position is masked `N` when its depth (from the depths file, Q>=20 bases) is below
  `--consensus_min_depth` (default 10). Confidently called variant sites are exempt, as in the old pipeline.
- **IUPAC ambiguity codes.** An SNV whose alternate frequency is above `--consensus_iupac_min_af` (0.25) and below
  `--consensus_freq` (0.75), supported by at least `--consensus_iupac_min_alt_reads` (5) reads, is written as the IUPAC
  code of reference and alternate base (`R Y K M S W`...): `K` at position 210 of `MIX_A`, where G (REF) and T (ALT, 31 of
  95 reads) coexist. `--consensus_iupac false` keeps the reference base instead. SNVs with ALT_FREQ at or above
  `--consensus_freq` become the alternate base; below `--consensus_iupac_min_af` they are ignored. Indels are applied only
  if ALT_FREQ >= `--indel_freq` (0.75) and are never written as IUPAC. Every called variant also needs total depth >=
  `--consensus_min_depth_variants` (10) and ALT_QUAL >= `--consensus_min_alt_qual` (20).
- **`-` deletions.** A deletion that is applied is written as `-` characters (bcftools `--mark-del '-'`); `-` is counted
  with `N` in `n_count`. ivar prints the deleted bases of a deletion as `N` (`-NNNNNNNNN`), so the alleles of
  indels are rebuilt from the reference (`--reference` of the converter) and only the _length_ of a printed deletion is used; the
  end-to-end test applies the common 9-base deletion at 11288-11296 and checks that exactly those nine positions are `-`. The
  synthetic dataset itself contains no deletions (the test barcodes carry only SNVs), so no `-` occurs in the example above.
- **Thresholds are applied here, not in Freyja.** Because `--freyja_varthresh` is 0.0, the variants TSV holds everything;
  the consensus filters decide what is called.
- **Status** (`consensus_stats.tsv`, `consensus_status.json`):
  - `ok`: a consensus was built. This includes a clean sample identical to the reference (header-only variants file).
  - `empty`: no usable depth (no position outside called variants reaches `--consensus_min_depth`, or an empty depths
    file). The FASTA record is all `N` at reference length (`LOWCOV_D`).
  - `failed`: a consensus step errored. The record is all `N` where one could be written. `reason` (last column of
    `consensus_stats.tsv`) says why. A missing status file is treated as `failed` by the report builder. A variant whose
    reference allele does not match the reference is not a failure: that variant is left out, the consensus is built from the
    rest (up to 10 such positions), and the sample is `ok` with the skipped variants in `reason` and in the `warnings` array of
    `consensus_status.json`. The report raises the WARN flag `consensus_failed` for `failed` samples and an INFO flag
    `consensus_warnings` for `ok` samples with skipped variants.
  - A consensus is never an empty file.

A low `n_pct` (and `ok`) is what you want; `empty`/high `n_pct` samples should not be used for sequence-based work.

## pipeline_info

| File                                                                                                               | Content                                                                                                                                                                                                          |
| ------------------------------------------------------------------------------------------------------------------ | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `run_info.json`                                                                                                    | Full provenance, see below. Written in `full`, `demix` and `aggregate`.                                                                                                                                          |
| `freyja_data.json`                                                                                                 | The `freyja_data` block alone (also inside `run_info.json`). Written by every mode that resolves Freyja data (`full`, `demix`, `aggregate` verified; `update_db` <!-- UNVERIFIED: not run, needs network -->).   |
| `freyja_data_compare.json`                                                                                         | Only with `--compare_barcodes` (`full`, `demix`): the same block for the second dataset; also `freyja_data_compare` in `run_info.json`.                                                                          |
| `software_versions.yml`                                                                                            | Tool versions per process (`FASTP: fastp: 1.3.6`, `FREYJA_DEMIX: freyja: 2.0.5`, ...).                                                                                                                           |
| `params_<ts>.json`                                                                                                 | All parameters of the run (76 entries incl. defaults and `null`s). Valid as `-params-file`.                                                                                                                      |
| `execution_report_<ts>.html`, `execution_timeline_<ts>.html`, `execution_trace_<ts>.txt`, `pipeline_dag_<ts>.html` | Nextflow reports. The trace also has a `container` column (the image of every task).                                                                                                                             |
| `discovery_<n>.log`                                                                                                | In `demix`/`aggregate` runs from results directories: what was found in the n-th results directory, which samples were skipped and why, and where `sample_type`/`collection_date`/`site` were carried over from. |
| `pipeline_report.{html,txt}`                                                                                       | Only with `--email`/`--email_on_fail` (nf-core completion e-mail). <!-- UNVERIFIED: not exercised -->                                                                                                            |

### Freyja data (`--db_dir`)

`--barcodes latest` downloads a dataset with `freyja update` and stores it in
`<db_dir>/freyja/<YYYY-MM-DDTHH-MM>/` (default `db_dir` is `~/.wasteline/db`): `usher_barcodes.feather`,
`last_barcode_update.txt`, `lineage_mutations.json`, `curated_lineages.json`, `lineages.yml`, `pathogen_config.yml`, and a
small `wasteline_dataset.json` with the download time. The directory name is the ISO form of `last_barcode_update.txt`
(the stamp has a time of day, so two updates on one day are two datasets). Later runs can pin any stored dataset with
`--barcodes 2026-09-17` (newest of that day), `--barcodes 09_17_2026-00-46` or `--barcodes 2026-09-17T00-46`. A run never
falls back silently to a different dataset.

## run_info.json

`pipeline_info/run_info.json` is the provenance record of the run: enough to say what produced the results and to
re-run it. It is assembled in two steps: Nextflow supplies what it knows, a script adds checksums, tool versions and
the Freyja dataset block; a completion hook appends `containers`, `run_completed`, `duration` and `success` when the run
ends (so a run that crashed before completion has no `success` key). Real example (params truncated):

```json
{
  "pipeline_name": "charlesfoster/wasteline",
  "pipeline_version": "0.1.0dev",
  "git_commit": "ea3709bf24f34b1dca3a91b881bc63a8c965ffff",
  "nextflow_version": "26.04.6",
  "mode": "full",
  "profile": "docker",
  "command_line": "nextflow run /path/to/wasteline -profile docker --outdir ... --mode full --input ... --barcodes .../barcodes_test.csv --primer_scheme midnight ...",
  "run_name": "clever_murdock",
  "start": "2026-10-16T02:59:59.631501+11:00",
  "params": { "outdir": "...", "mode": "full", "primer_scheme": "midnight", "...": "78 params in total" },
  "protocol": "amplicon",
  "primer_trimming": true,
  "duplicates": "none",
  "previous_runs": [],
  "tool_versions": {
    "Nextflow": "26.04.6",
    "bcftools": "1.23",
    "bwa": "0.7.19-r1273",
    "charlesfoster/wasteline": "v0.1.0dev",
    "fastp": "1.3.6",
    "freyja": "2.0.5",
    "ivar": "1.4.4",
    "python": "3.12.12",
    "samtools": "1.24"
  },
  "freyja_version": "2.0.5",
  "freyja_data": { "source": "custom", "barcode_version": "custom", "...": "see below" },
  "freyja_data_compare": {
    "source": "container",
    "barcode_version": "09_17_2026-00-46",
    "...": "only with --compare_barcodes"
  },
  "input_checksums": {
    "samplesheet": { "file": "samplesheet.csv", "sha256": "dd5bf2b4...", "size_bytes": 1302 },
    "primer_bed": { "file": "midnight.bed", "sha256": "51af40dc...", "size_bytes": 3014 },
    "reference_fasta": { "file": "NC_045512.fasta", "sha256": "0891c00c...", "size_bytes": 30429 }
  },
  "provenance_generated": "2026-10-15T16:01:15Z",
  "containers": {
    "FASTP": "community.wave.seqera.io/library/fastp:1.3.6--4df8d6c11b471bde",
    "FREYJA_DEMIX": "quay.io/biocontainers/freyja:2.0.6--pyhdfd78af_0",
    "...": "one entry per process"
  },
  "run_completed": "2026-10-16T03:01:17.380280+11:00",
  "duration": "1m 18s",
  "success": true
}
```

`protocol`, `primer_trimming` and `duplicates` are the **effective** library settings (PLAN D13): `protocol` is
`amplicon` or `capture`; `primer_trimming` is `false` for `--protocol capture` and `--skip_primer_trimming`; `duplicates` is
`mark` or `none` after `--duplicates auto` has been resolved. `primer_trimming` and `duplicates` are `null` in modes
without reads (`demix`, `aggregate`). For capture data `input_checksums` has no `primer_bed`.

| Key                                    | Meaning                                                                                                                                                                                                                                                      |
| -------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| `pipeline_name`, `pipeline_version`    | From the manifest.                                                                                                                                                                                                                                           |
| `git_commit`                           | Commit of the pipeline checkout that ran (the revision when launched from GitHub, otherwise `git rev-parse HEAD` of the pipeline directory). It does **not** record uncommitted changes.                                                                     |
| `nextflow_version`                     | Nextflow that ran it (>= 26.04.0 is required).                                                                                                                                                                                                               |
| `mode`                                 | `full`, `demix` or `aggregate`. (`report` and `update_db` write no `run_info.json`.)                                                                                                                                                                         |
| `profile`, `command_line`, `run_name`  | The `-profile`, the launch command and the Nextflow run name.                                                                                                                                                                                                |
| `start`                                | Launch time (with time zone).                                                                                                                                                                                                                                |
| `params`                               | Every pipeline parameter with its effective value (paths as given). `primer_scheme` is filled in as `midnight` even when the default applied. Excludes the nf-core `pipelines_testdata_base_path`.                                                           |
| `previous_runs`                        | `demix`/`aggregate` runs: the earlier runs the results build on, `[{path, run_info_sha256, barcode_version}]`. `run_info_sha256` is `null` and `barcode_version` is `"unknown"` for legacy Snakemake results, which recorded neither. Empty list for `full`. |
| `tool_versions`                        | Flat tool -> version map collected from the processes (`software_versions.yml`).                                                                                                                                                                             |
| `freyja_version`                       | The Freyja package version (also inside `freyja_data`).                                                                                                                                                                                                      |
| `freyja_data`                          | The Freyja dataset block, mandatory for every run (below).                                                                                                                                                                                                   |
| `freyja_data_compare`                  | Only with `--compare_barcodes`: the same block for the second dataset (`source`, `barcode_version`, sha256 values, ...). `freyja_data` always describes dataset A (`--barcodes`). Also written to `pipeline_info/freyja_data_compare.json`.                  |
| `input_checksums`                      | sha256 of the samplesheet, and (modes that map reads) the primer BED and reference FASTA, with file name and size.                                                                                                                                           |
| `provenance_generated`                 | UTC time the record was assembled.                                                                                                                                                                                                                           |
| `containers`                           | Image of every process, exact tag (or `@sha256:` digest with `-profile locked`). Added at completion.                                                                                                                                                        |
| `run_completed`, `duration`, `success` | Added at completion.                                                                                                                                                                                                                                         |
| `freyja_data_note`                     | `aggregate` mode only: explains that `freyja_data` describes the dataset used to group the plot, not the barcodes the `.demix` files were made with (look at `previous_runs`).                                                                               |

### The `freyja_data` block

```json
{
  "source": "custom",
  "barcode_version": "custom",
  "barcode_version_iso": null,
  "barcodes_file": "barcodes_test.csv",
  "barcodes_sha256": "0c3cbd7a1e6a7e6a8ce3e36ca22f4060c483a65e698fff0bec8e3226639509c7",
  "barcodes_format": "csv",
  "lineage_mutations_sha256": null,
  "curated_lineages_sha256": null,
  "lineages_yml_sha256": null,
  "pathogen": "SARS-CoV-2",
  "freyja_version": "2.0.5",
  "freyja_reported_barcode_version": null,
  "freyja_bundled_barcode_version": "09_17_2026-00-46",
  "dataset_path": "/path/to/tests/data/barcodes_test.csv",
  "downloaded_at": null
}
```

| Key                                                                          | Meaning                                                                                                                                                                                                                                                                                                                                                    |
| ---------------------------------------------------------------------------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `source`                                                                     | How `--barcodes` was resolved: `latest` (downloaded now with `freyja update` and stored in `--db_dir`), `container` (the snapshot bundled in the pinned Freyja container, offline), `pinned` (a version/date resolved to a dataset already stored in `--db_dir`), `custom` (a dataset directory or a `.csv`/`.feather` file given by path, used verbatim). |
| `barcode_version`                                                            | The verbatim content of the dataset's `last_barcode_update.txt` (`MM_DD_YYYY-HH-MM`), or `custom` for a bare barcode file.                                                                                                                                                                                                                                 |
| `barcode_version_iso`                                                        | ISO form `YYYY-MM-DDTHH-MM` (names the directory in `--db_dir`); `null` for `custom`.                                                                                                                                                                                                                                                                      |
| `barcodes_file`, `barcodes_format`                                           | Name and `csv`/`feather` format of the barcode file.                                                                                                                                                                                                                                                                                                       |
| `barcodes_sha256`                                                            | sha256 of the barcode file: the strongest proof of which barcodes were used.                                                                                                                                                                                                                                                                               |
| `lineage_mutations_sha256`, `curated_lineages_sha256`, `lineages_yml_sha256` | sha256 of the other dataset files; `null` when absent (always for a bare barcode file).                                                                                                                                                                                                                                                                    |
| `pathogen`                                                                   | Freyja `--pathogen` (`SARS-CoV-2`; set in `conf/pathogens/sars2.config`).                                                                                                                                                                                                                                                                                  |
| `freyja_version`                                                             | Freyja package version that ran.                                                                                                                                                                                                                                                                                                                           |
| `freyja_reported_barcode_version`                                            | The barcode version printed by `freyja demix --version`, **only when `source` is `container`**; otherwise `null`.                                                                                                                                                                                                                                          |
| `freyja_bundled_barcode_version`                                             | What `freyja demix --version` printed, i.e. the version of the data **bundled in the container's package**.                                                                                                                                                                                                                                                |
| `dataset_path`                                                               | Where the dataset came from (for `latest`: the `<db_dir>/freyja` parent; for `container`: `bundled in the Freyja container`; otherwise the path given).                                                                                                                                                                                                    |
| `downloaded_at`                                                              | UTC download time from the dataset's `wasteline_dataset.json` sidecar (datasets fetched with `latest`), else `null`.                                                                                                                                                                                                                                       |

Why two version fields (decision D9): `freyja demix --version` can only ever report the data bundled in the Freyja
package; it cannot confirm an external dataset passed with `--barcodes`. The authoritative record of an external dataset
is its own `last_barcode_update.txt` plus the sha256 values, which is what `barcode_version` and `barcodes_sha256` hold.
So `freyja_reported_barcode_version` is filled only for `source=container` (where the two coincide), and
`freyja_bundled_barcode_version` always shows what is bundled, as information. In the example above a custom barcode CSV
was used, so `barcode_version` is `custom`, the reported version is `null` and the bundled container data happens to be
`09_17_2026-00-46`: the last field says nothing about the barcodes used.

The barcode version is also written to every `<s>.demix_status.json`, `summary/demix_status.tsv`,
`summary/barcode_version.json`, `report_data.json` and the HTML report header (with `--compare_barcodes` also in
`barcode_compare.tsv` and `barcode_compare_status.tsv`, for dataset B). A report flag warns when a barcode is
`custom`, from the `container` snapshot, `unknown` (legacy results) or when samples were demixed with different versions.

## How to reproduce a run

`run_info.json` and `params_<ts>.json` record what is needed. In order:

1. **Pipeline version.** Check out `git_commit` (and note it is the commit only; uncommitted changes are not
   recorded). Or run the released tag that corresponds to `pipeline_version`.

   ```bash
   python3 -c "import json; d=json.load(open('results/pipeline_info/run_info.json')); print(d['git_commit'], d['nextflow_version'], d['freyja_version'])"
   ```

2. **Nextflow.** Use `nextflow_version` (`NXF_VER=26.04.6 nextflow run ...`); the minimum is 26.04.0.
3. **Containers.** The Freyja version and every other tool version are pinned exactly by the pipeline commit; the images
   used are in `containers`. For immutable images add `-profile docker,locked` (digest-pinned, `conf/containers.lock.config`).
4. **Inputs.** Use the same samplesheet and reads; confirm with the checksum:

   ```bash
   sha256sum samplesheet.csv   # must equal run_info.json input_checksums.samplesheet.sha256
   ```

5. **Barcodes** (the part that changes over time; Freyja data is updated almost daily):

   | `freyja_data.source` | What to do                                                                                                                                                                                                                                                                                                                             |
   | -------------------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
   | `latest` or `pinned` | The dataset is still in `<db_dir>/freyja/<barcode_version_iso>/`: use `--barcodes <barcode_version_iso>` (or the stamp). If that directory is gone it **cannot be downloaded again** (old barcode versions are not available from `freyja update`), so archive `--db_dir` datasets you may need to reproduce. Check `barcodes_sha256`. |
   | `container`          | Same pinned Freyja container => same data: `--barcodes container`.                                                                                                                                                                                                                                                                     |
   | `custom`             | Supply the same file or directory: `--barcodes <dataset_path>`; compare `sha256sum` with `barcodes_sha256`.                                                                                                                                                                                                                            |

   If the run used `--compare_barcodes`, dataset B is recorded in `freyja_data_compare` and is reproduced the same way
   (`params_<ts>.json` holds the `--compare_barcodes` value as given, which for `latest` is not enough on its own: pin the
   recorded `barcode_version_iso` instead).

6. **Parameters.** Reuse the recorded parameters; edit `outdir` (and any paths that moved):

   ```bash
   cp results/pipeline_info/params_*.json params.json    # edit "outdir" in params.json
   nextflow run charlesfoster/wasteline -r <git_commit> -profile docker -params-file params.json
   ```

   (Verified: the file written by a real run works as `-params-file`.) The exact original command is in
   `run_info.json` `command_line`.

To redo only the lineage estimate with a different barcode dataset, do not re-run everything: use
`--mode demix --previous_results results/ --barcodes <other>`; the new run records the earlier one in `previous_runs`.
For a byte-for-byte check, `<s>.demix` files and `freyja_aggregated.tsv` of a `demix` run with the same barcodes are
identical to the `full` run (verified on the synthetic dataset).
