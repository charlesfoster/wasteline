# charlesfoster/wasteline: Usage

wasteline estimates SARS-CoV-2 lineage abundances from wastewater Illumina amplicon data with
[Freyja](https://github.com/andersen-lab/Freyja), and also produces per-sample QC, amplicon-dropout tables, consensus
genomes, negative-control checks and a single self-contained HTML report. This page is the user guide. For a five-minute
trial on synthetic data start with [QUICKSTART.md](QUICKSTART.md). Every parameter is listed in
[parameters.md](parameters.md) (generated from `nextflow_schema.json`); every output file in [output.md](output.md);
problems are in [troubleshooting.md](troubleshooting.md).

> **Examples use `nextflow run charlesfoster/wasteline`.** To run a local checkout instead, replace that with
> `nextflow run /path/to/wasteline` (or `nextflow run .` from inside it).
> <!-- UNVERIFIED: the GitHub repository was not reachable (404) when this page was written. Every command here was run from a local checkout with `nextflow run <path>`; `-r <version>`/`nextflow pull` were not tested. -->

Contents: [1 Prerequisites](#1-prerequisites) | [2 Modes](#2-modes-at-a-glance) | [3 Samplesheets](#3-samplesheets) |
[4 Reference and primer schemes](#4-reference-and-primer-schemes) | [5 Running each mode](#5-running-each-mode) |
[6 Freyja data (barcodes)](#6-freyja-data-barcodes) | [7 Freyja settings](#7-freyja-settings) |
[8 Consensus genomes](#8-consensus-genomes) | [9 Negative controls and QC](#9-negative-controls-and-qc-thresholds) |
[10 Mutations of interest](#10-mutations-of-interest) | [11 Resuming](#11-resuming-a-run--resume) |
[12 Profiles](#12-profiles) | [13 Resources](#13-resources) | [14 Parameter files](#14-parameter-files-and-boolean-flags) |
[15 HPC](#15-running-on-an-hpc-cluster) | [16 Updating Freyja](#16-updating-the-freyja-version-of-the-pipeline) |
[17 Not implemented](#17-not-implemented)

Library types: [probe-capture data](#probe-capture-data-eg-twist---protocol-capture) (`--protocol capture`).
Optional extras, each with a worked example: [comparing barcode versions](#comparing-barcode-versions---compare_barcodes) (`--compare_barcodes`),
[bootstrap intervals](#bootstrap-confidence-intervals---freyja_boot) (`--freyja_boot`),
[variant annotation](#annotating-variants---annotate_variants) (`--annotate_variants`),
[mutations of interest](#10-mutations-of-interest), and [time series](#time-series-growth-rates-and-the-dashboard).

## 1. Prerequisites

| Need                         | Details                                                                                                                                                                                                                                                                                                                                                                                              |
| ---------------------------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| Nextflow **>= 26.04.0**      | `manifest.nextflowVersion` is `!>=26.04.0`. Older releases cannot parse the pipeline (see [troubleshooting](troubleshooting.md#script-compilation-error-on-an-old-nextflow)). Check with `nextflow -version`; upgrade with `nextflow self-update`, or run one command with `NXF_VER=26.04.6 nextflow run ...`.                                                                                       |
| A container engine           | **Docker** (the tested one; Docker Desktop must be running on macOS/Windows), or Singularity/Apptainer/Podman. Conda/mamba profiles exist but are not the recommended route.                                                                                                                                                                                                                         |
| Internet on the first run    | Container images are pulled from `quay.io` and `community.wave.seqera.io` (Singularity/Apptainer: `depot.galaxyproject.org` and `community-cr-prod.seqera.io`). `--barcodes latest` (the default) also downloads the current Freyja data. See [offline use](#offline-use).                                                                                                                           |
| Python **>= 3.9** (optional) | Only for running the helper scripts in `bin/` yourself (`make_samplesheet.py`, `check_samplesheet.py`, ...). The pipeline runs its own copies inside containers. `make_samplesheet.py` and `check_samplesheet.py` are verified on Python 3.9 (3.8 syntax is checked statically); the test-only `tests/check_outputs.py` needs Python >= 3.12 (the repository's `pixi install` environment has 3.12). |
| Disk                         | One work directory (`work/`) holding every intermediate file, the results (`--outdir`), and about 16 MB per stored Freyja dataset in `--db_dir`. Delete `work/` when you are done (after which `-resume` no longer works).                                                                                                                                                                           |

Memory and CPUs: the heaviest steps are `BWA_MEM` (label `process_high`: 8 CPUs, 16 GB requested), `FASTP`, `IVAR_TRIM`
and `SAMTOOLS_SORT` (`process_medium`: 4 CPUs, 8 GB). Requests are automatically capped at what the machine has (see
[Resources](#13-resources)). The synthetic dataset runs in about two minutes on a laptop.

macOS note: the Freyja and most bioconda images are `linux/amd64` only. `-profile docker` adds `--platform=linux/amd64`
automatically on Apple Silicon, so everything runs under emulation (slower, but the test dataset still finishes in
about two minutes on an M-series Mac). Use `-profile docker`, not conda, on macOS.

## 2. Modes at a glance

`--mode full|demix|aggregate|report|update_db` is the **only** way to choose what runs (default `full`). Nextflow's
`-entry` is not supported. `--outdir` is required in every mode.

| Mode             | Input                                                                                                                | Does                                                                                              | Does not run                                          |
| ---------------- | -------------------------------------------------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------- | ----------------------------------------------------- |
| `full` (default) | FASTQ samplesheet (`--input`)                                                                                        | everything: trim, map, primer-trim, Freyja variants, demix, aggregate/plot, QC, consensus, report |                                                       |
| `demix`          | `--input` sheet of existing Freyja `variants` + `depths`, **or** `--previous_results <dir>`                          | resolve Freyja data, `freyja demix`, aggregate/plot, QC, consensus, report                        | fastp, bwa, ivar, `freyja variants`                   |
| `aggregate`      | `--input` sheet of existing `.demix` files, **or** `--demix_dir <dir>`, **or** `--previous_results <dir>[,<dir>...]` | `freyja aggregate` + plots, report                                                                | demix and everything upstream (so no QC or consensus) |
| `report`         | `--results_dir <dir>` of an earlier run                                                                              | re-render `summary/report.html`                                                                   | everything else                                       |
| `update_db`      | none                                                                                                                 | download the current Freyja data into `--db_dir` and exit                                         | everything else                                       |

Input rules (checked at start-up, before any work is done): the input options that a mode accepts are `full`: `--input`;
`demix`: `--input` or `--previous_results`; `aggregate`: `--input`, `--demix_dir` or `--previous_results`; `report` and
`update_db`: none (a stray `--input` is ignored with a warning). Give exactly one source. `--results_dir` is only for `report`.

## 3. Samplesheets

A samplesheet is a comma-separated file with a header. Its columns depend on the mode. The schemas are
`assets/schema_input.json` (full), `assets/schema_input_demix.json` and `assets/schema_input_aggregate.json`; they are
enforced with nf-schema.

Rules common to all modes:

- `sample` matches `^[A-Za-z0-9._-]+$` (letters, digits, `.`, `_`, `-`; no spaces) and is **unique**. IDs are used verbatim
  everywhere: **sample IDs are never rewritten** (the old Snakemake pipeline turned `_` into `-`; wasteline does not).
  Duplicate IDs are an error, not a merge: concatenate re-sequenced runs or lanes of one sample yourself first.
- `sample_type` is `sample` (default) or `negative`. Negative controls are assessed differently (see
  [section 9](#9-negative-controls-and-qc-thresholds)).
- `collection_date` is an ISO date `YYYY-MM-DD` or empty. If any samples have dates, a time-aware plot
  (`freyja_plot_times.pdf`) and the report's time-series panel are produced from the dated samples (see
  [time series](#time-series-growth-rates-and-the-dashboard)).
- `site` is free text.
- `viral_load` (optional, any mode) is a positive number: the concentration of virus in the sample, in any unit you use
  consistently (for example gene copies per litre of wastewater). Leave it blank where unknown. It is used only for the
  load-weighted time series and, with `--freyja_growth_scale_by_viral_load`, the growth-rate fit
  ([below](#time-series-growth-rates-and-the-dashboard)). `bin/make_samplesheet.py --metadata` accepts the column and carries it
  over from `summary/samplesheet.csv` of earlier results; it is written to the sheet only when some sample has a value.
- Relative file paths are resolved against the **directory you launch Nextflow from**, not the samplesheet's location.
  Use absolute paths (the default of `make_samplesheet.py`).

### `--mode full`: reads

```csv
sample,fastq_1,fastq_2,sample_type,collection_date,site
MIX_A,/data/run7/MIX_A_R1.fastq.gz,/data/run7/MIX_A_R2.fastq.gz,sample,2026-01-05,SiteA
NEG_E,/data/run7/NEG_E_R1.fastq.gz,/data/run7/NEG_E_R2.fastq.gz,negative,2026-01-05,SiteA
```

Required: `sample,fastq_1,fastq_2`. **Paired-end reads only** (`fastq_2` may not be empty). FASTQ paths must end in
`.fastq.gz` or `.fq.gz` and contain no whitespace; the files must exist.

### `--mode demix`: existing Freyja variants and depths

```csv
sample,variants,depths,fastp_json,sample_type,collection_date,site
MIX_A,results/samples/MIX_A/freyja/MIX_A.variants.tsv,results/samples/MIX_A/freyja/MIX_A.depths,results/samples/MIX_A/fastp/MIX_A.fastp.json,sample,2026-01-05,SiteA
```

Required: `sample,variants,depths`. `variants` is a Freyja/ivar variants table (`.tsv` or `.tsv.gz`), `depths` the Freyja depths
file (plain or gzipped); `fastp_json` is optional and only adds read counts to the QC table. Rather than writing this
by hand, generate it with `make_samplesheet.py --from-results` or give `--previous_results` (below).

### `--mode aggregate`: existing `.demix` files

```csv
sample,demix,sample_type,collection_date,site
MIX_A,results/samples/MIX_A/freyja/MIX_A.demix,sample,2026-01-05,SiteA
```

Required: `sample,demix`.

### Making a samplesheet: `bin/make_samplesheet.py`

From a directory of reads (Illumina names such as `<id>_S<n>_L001_R1_001.fastq.gz`; `_S<n>_L001` is stripped from the ID,
R2 is derived from R1, `Undetermined*` is skipped, and `NEG_`/`NC_` in the file name gives `sample_type=negative`):

```bash
python bin/make_samplesheet.py /data/run7/ -o samplesheet.csv
# other file naming, with dates and sites from a metadata table (keyed by `sample`)
python bin/make_samplesheet.py /data/run7/ --suffix _R1.fastq.gz --metadata meta.csv -o samplesheet.csv
```

```
wrote 7 sample(s) to sheet_full.csv
```

Options: `--suffix` (default `_L001_R1_001.fastq.gz`), `--negative-pattern` (default `NEG_|NC_`), `--exclude-pattern`
(default `Undet`), `--metadata`, `--recursive`, `--relative`, `-o`. The script exits 1 and writes nothing if R2 is missing, if two
files give the same ID, if an ID has characters outside `[A-Za-z0-9._-]`, or if the metadata is invalid. Full reference:
[helper_scripts.md](helper_scripts.md). `--relative` writes paths relative to the **output file's** directory, which only
works with the pipeline if you then launch Nextflow from that same directory.

`meta.csv` example (a sample with no row only produces a warning):

```csv
sample,collection_date,site,sample_type
Site_1_week_2,2026-02-02,Site_1,
```

From the results of an earlier run (wasteline **or** the old Snakemake pipeline; both layouts are detected per sample):

```bash
python bin/make_samplesheet.py --from-results results/ --mode demix     -o demix_sheet.csv
python bin/make_samplesheet.py --from-results results/ --mode aggregate -o aggregate_sheet.csv
```

```
NOTE: carrying over sample_type/collection_date/site from results/summary/samplesheet.csv
NOTE: skipped 'LOWCOV_D' (new layout): missing demix
NOTE: skipped 'NEG_E' (new layout): missing demix
wrote 5 sample(s) to agg_sheet.csv
```

Samples lacking the required files are skipped with a reason; if nothing usable remains the command fails.
`sample_type`, `collection_date` and `site` are carried over from `results/summary/samplesheet.csv` when present.

### Checking a samplesheet: `bin/check_samplesheet.py`

```bash
python bin/check_samplesheet.py samplesheet.csv --mode full --check-files
```

It prints a JSON object to stdout and a table to stderr and exits 0 if the sheet is valid, 1 if not. `--check-files` also
verifies that every file exists, FASTQs are readable gzip files starting with `@`, and other files are not empty.

```
FAIL  bad.csv  mode=full  samples=1  errors=1  warnings=0
LEVEL  MESSAGE
ERROR  line 2 (A): empty fastq_2 (paired-end data is required)
```

The pipeline performs the same checks itself (with nf-schema) at start-up, so this script is optional; it is useful to
validate a sheet before queuing a long run and to check files without starting Nextflow.

## 4. Reference and primer schemes

**Reference.** The default is the bundled SARS-CoV-2 reference `NC_045512.2` (`assets/references/NC_045512.fasta`) and its
GFF3 annotation. Override with `--fasta` and `--gff`. Only SARS-CoV-2 is supported at present (the pathogen defaults live in
`conf/pathogens/sars2.config`).

**Primer scheme** (modes `full`, and `demix` if you want per-amplicon coverage): choose a bundled scheme with
`--primer_scheme`, or give your own BED with `--primer_bed <file>` (which overrides `--primer_scheme`). The default in
`full` mode is **`qiagenv1`** (QIAseq Direct SARS-CoV-2; the old Snakemake pipeline defaulted to `midnight`, see
[migration_from_snakemake.md](migration_from_snakemake.md#the-default-primer-scheme-changed)): if you use another protocol, say
so with `--primer_scheme`. The scheme in effect is recorded in `pipeline_info/run_info.json` (`params.primer_scheme`), also when
you did not give one. In mode `demix` amplicon-level QC happens only when you give a scheme or BED explicitly (the scheme of
old results is not known). The BED's contig name **must** match a contig of the reference FASTA, otherwise the run stops
at start-up (otherwise `ivar trim` would silently trim nothing).

| `--primer_scheme`    | Amplicons | Primer names look like        |
| -------------------- | --------- | ----------------------------- |
| `midnight`           | 29        | `nCoV-2019_1_F`               |
| `midnightv2`         | 29        | `SARSCoV_1200_1_F`            |
| `articv3`            | 98        | `nCoV-2019_1_F`               |
| `articv4.1`          | 99        | `nCoV-2019_1_F`               |
| `articv5.3.2`        | 96        | `SARS-CoV-2_400_1_F`          |
| `eden`               | 14        | `SARS2_A1_F`                  |
| `swift`              | 343       | `covid19genome_0-200_s0_M1_F` |
| `qiagenv1` (default) | 223       | `QIAseq_221_LEFT`             |
| `xgen`               | 343       | `covid19genome_0-200_s0_M1F`  |

All nine bundled BEDs use the contig name `NC_045512.2`, the same as the bundled reference, so every scheme works as-is
(the start-up contig check passes for each of them; verified with `-preview`). Amplicon counts were measured with
`bin/amplicon_coverage.py` on each bundled BED.

**A custom BED must use the contig name of the reference.** A BED copied from elsewhere often names the contig
`MN908947.3` (the same Wuhan-Hu-1 sequence, a different accession). The run then stops at start-up, before any work is done:

```
[ERROR] Primer BED '/path/my_scheme.bed' refers to contig(s) [MN908947.3] that are not in the reference FASTA '.../NC_045512.fasta' (which has [NC_045512.2]). Use a matching --primer_bed/--primer_scheme or --fasta.
```

The check exists on purpose: with a mismatched contig `ivar trim` would silently trim nothing. Fix it by renaming the contig
in the first column of your BED (the coordinates are unchanged when the sequence is the same) or by giving the matching
`--fasta`:

```bash
sed 's/^MN908947\.3/NC_045512.2/' my_scheme.bed > my_scheme.NC_045512.2.bed
nextflow run charlesfoster/wasteline -profile docker --input samplesheet.csv --outdir results --primer_bed my_scheme.NC_045512.2.bed
```

Custom BEDs: tab-separated BED6-like (`contig  start  end  name  pool  strand`), with primer names that end `_F`/`_R`,
`_LEFT`/`_RIGHT` or `F`/`R` so that amplicons can be reconstructed for `amplicon_coverage.tsv` (see
[qc_scripts.md](qc_scripts.md) for the naming rules). `--ivar_offset` sets `ivar trim -x` (default 5 for `swift`,
otherwise 0). If you use a different amplicon length distribution, expect more or fewer `amplicon_dropout` flags; tune
`--amplicon_min_depth` and `--qc_max_dropped_amplicons`.

### Probe-capture data (e.g. Twist): `--protocol capture`

Libraries made by **probe capture** (hybrid capture, e.g. Twist panels) are not PCR amplicons: fragments start anywhere
on the genome and the reads do **not** begin with primer sequences. Trimming primers from them would cut real sequence
and an "amplicon dropout" analysis has no meaning, so for these libraries use `--protocol capture` (default `amplicon`).

```bash
nextflow run charlesfoster/wasteline -profile docker \
    --protocol capture --input samplesheet.csv --outdir results --barcodes latest
```

What changes with `--protocol capture`:

|                                      | `amplicon` (default)                                  | `capture`                                                                                                               |
| ------------------------------------ | ----------------------------------------------------- | ----------------------------------------------------------------------------------------------------------------------- |
| Primer trimming (`ivar trim`)        | yes                                                   | **no**: the mapped BAM goes straight to the sorted, indexed BAM                                                         |
| Primer scheme / BED                  | required (default `qiagenv1`)                         | **ignored** (a warning says so; no start-up contig check)                                                               |
| Amplicon coverage and dropout        | `amplicon_coverage.tsv`, `amplicon_dropout` QC reason | **not produced**; `amplicons_total` / `amplicons_dropped` are `NA`; the report's amplicon section says "not applicable" |
| PCR duplicates (`--duplicates auto`) | left alone                                            | **marked** (`samtools markdup`; flagged, never removed)                                                                 |
| BAM names                            | `<id>.primertrimmed.bam`                              | `<id>.mapped.markdup.bam` (+ `<id>.markdup_stats.txt`)                                                                  |

Everything else is unchanged: Freyja variants/demix/aggregate, coverage-based QC, negative controls, consensus,
mutations of interest, `--annotate_variants`, `--freyja_boot`, `--compare_barcodes` and the report (whose header shows the
protocol). The protocol, whether primers were trimmed and the effective duplicate handling are recorded in
`pipeline_info/run_info.json` (`protocol`, `primer_trimming`, `duplicates`).

**Duplicates.** `--duplicates auto|mark|none` (default `auto`: `mark` for capture, `none` for amplicon). Marking flags
duplicate read pairs (same fragment start and end) in the BAM; `freyja variants` runs `samtools mpileup`, which skips
flagged reads, so duplicates do not count towards depth or allele frequencies. Capture libraries have high duplicate rates
that inflate depth without adding information; amplicon libraries have many natural duplicates (every read of an amplicon
starts at the primer), so marking them is off by default. `--duplicates mark` also works with `amplicon`;
`--duplicates none` leaves capture duplicates in. The marking statistics are published as
`samples/<id>/bam/<id>.markdup_stats.txt`. `qc_summary.tsv` read counts (`mapped_reads`) still include duplicates; the
mean depth does not.

**`--skip_primer_trimming`** is for amplicon reads that were already primer-trimmed elsewhere: no `ivar trim` and, as for
capture, no amplicon-level QC (a primer scheme/BED would be ignored). It is implied by `--protocol capture`.

**Not available for capture:** per-amplicon coverage/dropout. In mode `demix` the protocol only matters for the optional
amplicon coverage (it is never computed for `capture`); demix has no BAM, so duplicates are not touched there.

### Primer scheme check

A wrong primer scheme is a quiet error: `ivar trim` trims the wrong places and the amplicon dropout calls are wrong, but
everything still runs. In mode `full` (amplicon protocol with primer trimming) the pipeline therefore checks the declared
scheme against the reads, on the mapped BAM **before** trimming: it counts the mapped reads whose 5' end lies within
`--primer_check_window` (2) bp of a same-strand primer 5' end, for the declared scheme and for every bundled scheme, and reports
the enrichment over chance (and a binomial p-value). Per sample, `summary/qc_summary.tsv` gets four columns:

| Column                     | Meaning                                                                                                                                                                                                                                                                                                                                                                                                                                       |
| -------------------------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `primer_scheme_check`      | `ok` (the declared scheme is enriched at least `--primer_check_min_enrichment` times, 5, and significantly), `mismatch` (the declared scheme is not, and another, non-equivalent scheme is enriched with at least 3 times its enrichment), `uninformative` (no scheme stands out, for example because the reads do not start at primers; or the best one is not decisively better) or `NA` (fewer than 500 usable reads, no pre-trimming BAM) |
| `primer_scheme_best`       | the best-matching scheme(s), joined by `/` (`midnight/midnightv2` have identical coordinates)                                                                                                                                                                                                                                                                                                                                                 |
| `primer_scheme_enrichment` | enrichment of the declared scheme (1 decimal)                                                                                                                                                                                                                                                                                                                                                                                                 |
| `primer_trimmed_pct`       | percentage of mapped reads from which `ivar trim` removed primer sequence (from the ivar log)                                                                                                                                                                                                                                                                                                                                                 |

A `mismatch` adds the WARN reason `primer_scheme_mismatch` to the sample and a run-level flag to the report; `uninformative` is
not a warning. The per-scheme table is in `samples/<s>/coverage/<s>.primer_scheme_check.tsv` (one row per scheme) and the
one-row summary next to it (`<s>.primer_scheme_summary.tsv`). Real example: midnight-shaped reads with `--primer_scheme midnight`
give enrichment about 41-43 (`ok`, best `midnight/midnightv2`); declared `articv4.1` gives `mismatch`, best
`midnight/midnightv2`.

**Limits.** The check needs libraries whose reads begin at primer sites, such as QIAseq or ligation/PCR libraries with short
amplicons. Tagmentation-based libraries (a Nextera-style prep of long amplicons, for example the Midnight protocol on
Nextera) break amplicons at random positions, so read starts do not mark primers and the result is `uninformative` or
`ok`-with-low-enrichment even when the scheme is right. A scheme that is a subset of or very similar to another cannot be told
apart (`midnight` and `midnightv2`). The method has been calibrated on the synthetic data only (it separates the right scheme
by a factor of about 10 from all others there); on real runs watch the enrichment of known-good runs first and consider
`--primer_check_min_enrichment` if your libraries are less strongly anchored. `--primer_scheme_check false` switches it off;
`--qc_min_primer_trimmed_pct` (default 0, off) is a blunter WARN on the ivar trimming percentage alone. Details of the
method are in [qc_scripts.md](qc_scripts.md#primer_scheme_checkpy).

## 5. Running each mode

All examples use `-profile docker`; add `-resume` to continue an interrupted run. Pipeline parameters take **two** dashes
(`--outdir`), Nextflow's own options one (`-profile`, `-resume`).

### `full`: from reads

```bash
nextflow run charlesfoster/wasteline -profile docker \
    --input samplesheet.csv \
    --outdir results \
    --primer_scheme midnight \
    --barcodes latest
```

(`--primer_scheme` is given here because the reads are Midnight amplicons; without it `qiagenv1` is used.)

Expect (per sample): fastp trimming, `bwa mem`, `ivar trim` with the primer BED, `freyja variants` (variants and depths),
then, once per run, the Freyja data and `freyja demix` per sample, `freyja aggregate` + plots, QC, consensus genomes and the
report. A real run on the synthetic data finishes with:

```
-[charlesfoster/wasteline] Pipeline completed successfully-
[SUCCESS] completed=112 failed=0 cached=0
```

Outputs are in `results/summary/` (run-level tables and the report), `results/samples/<sample>/` and
`results/pipeline_info/`; see [output.md](output.md). **A failed or skipped sample never aborts the run**: it gets a row in
`summary/demix_status.tsv` with a status (`ok`, `failed`, `skipped_low_coverage` or, for unusable inputs found in old
results, `missing`) and a reason:

```
sample         status                reason                                         genome_coverage_pct  barcode_version
LOWCOV_D       skipped_low_coverage  less than 60% of the genome covered at >= 10x  0.0                  custom
MIX_A          ok                                                                   99.3345              custom
NEG_CONTAM     ok                                                                   94.5457              custom
NEG_E          skipped_low_coverage  less than 60% of the genome covered at >= 10x  0.0                  custom
```

(`NEG_CONTAM` is a contaminated negative control: demixed, so you can see what it contains, but not in the aggregate; see
[section 9](#9-negative-controls-and-qc-thresholds).)

### `demix`: re-estimate lineages from existing variants and depths

Use this to re-run the lineage estimate with a **different barcode dataset** (or Freyja settings) without repeating
mapping. It replaces the old pipeline's `--redo_demix`, but **never modifies** the old results; everything is written to
`--outdir`.

```bash
# from a results directory made by wasteline or by the old Snakemake pipeline
nextflow run charlesfoster/wasteline -profile docker \
    --mode demix \
    --previous_results results \
    --barcodes 2026-10-06 \
    --outdir results_redemix

# or from a samplesheet (edit it to choose samples)
python bin/make_samplesheet.py --from-results results/ --mode demix -o demix_sheet.csv
nextflow run charlesfoster/wasteline -profile docker --mode demix --input demix_sheet.csv --barcodes latest --outdir results_redemix
```

What you get: `summary/` with `demix_status.tsv`, `freyja_aggregated.tsv`, plots, `qc_summary.tsv`, `consensus_genomes.fasta`,
`report.html`, and `samples/<s>/freyja/<s>.demix` (+ status), consensus and coverage files. The variants and depths files
are not copied again, so a `demix` output directory can be aggregated or reported on, but not demixed again; point later
`demix` runs at the original results. The optional extras [`--compare_barcodes`](#comparing-barcode-versions---compare_barcodes),
[`--freyja_boot`](#bootstrap-confidence-intervals---freyja_boot) and [`--annotate_variants`](#annotating-variants---annotate_variants)
also work in `demix` mode (that is the cheap way to add them to an existing run).

- With the same barcodes, `demix` reproduces the `.demix` files and `freyja_aggregated.tsv` of the `full` run
  byte-for-byte (verified with `cmp` on the synthetic data).
- Per-amplicon coverage and dropout are only computed in `demix` mode when you also name the primer scheme used
  (`--primer_scheme` or `--primer_bed`): the scheme of old results is not known to the pipeline.
- **Legacy Snakemake results** (`<s>/freyja/<s>.freyja_variants.tsv`, `<s>/freyja/<s>.depths`, gzipped variants/depths
  accepted) are supported; their previous barcode version is recorded as `unknown`. Their variants were made with
  `--varthresh 0.25`, so a re-demix reuses that filtering. See [migration_from_snakemake.md](migration_from_snakemake.md).
  The synthetic legacy fixture (`tests/fixtures/demix_legacy`) demixes with `--mode demix --previous_results`;
  <!-- UNVERIFIED: real Freyja 1.5.3 (Snakemake-era) variants/depths files were not available to test -->
- A negative control with real coverage is demixed too (so contamination is visible) but kept out of the aggregate
  ([section 9](#9-negative-controls-and-qc-thresholds)).
- `--include_negatives_in_aggregate` puts the controls into `freyja_aggregated.tsv` and the plots instead (verified: `NEG_CONTAM`
  then appears in the aggregate and `freyja_aggregated_controls.tsv` is not written).

### `aggregate`: aggregate and plot existing `.demix` files

```bash
# a results directory (carries over sample_type, collection_date and site; several directories: comma-separated)
nextflow run charlesfoster/wasteline -profile docker --mode aggregate \
    --previous_results results --barcodes container --outdir results_aggregate

# or every *.demix file below a directory (the ID is the file name without .demix)
nextflow run charlesfoster/wasteline -profile docker --mode aggregate \
    --demix_dir /data/old_results/run7 --barcodes container --outdir results_aggregate

# or a samplesheet
nextflow run charlesfoster/wasteline -profile docker --mode aggregate --input aggregate_sheet.csv --barcodes container --outdir results_aggregate
```

Outputs: `summary/freyja_aggregated.tsv`, `freyja_aggregated_controls.tsv` (when the inputs contain negative controls),
`freyja_plot.pdf` (and `freyja_plot_times.pdf` if dates are known), `demix_status.tsv`, `report.html`. Only samples that have a non-empty `.demix` file are aggregated; the others appear as
`missing` in `demix_status.tsv` with a reason (for example `missing demix (new layout)`; zero-byte legacy `.demix` files
are reported as `empty demix file (failed in the old pipeline)`). Nothing is demixed, so `--barcodes` here only
supplies the lineage hierarchy used to group the plot: `--barcodes latest` (the default) would download data just for
that, so `container` is the sensible choice offline. With `--demix_dir` there is no metadata, so no time-aware plot; use
`--previous_results` or a samplesheet to keep `collection_date` (then the plot and the report's time-series panel are produced,
see [time series](#time-series-growth-rates-and-the-dashboard)). `--compare_barcodes`, `--freyja_boot` and `--annotate_variants` are
**ignored with a warning** in `aggregate` mode (they need variants and depths):

````
WARN: --compare_barcodes is only used by --mode full and demix and is ignored in --mode aggregate
WARN: --freyja_boot needs variants and depths: it is only used by --mode full and demix and is ignored in --mode aggregate
WARN: --annotate_variants is only used by --mode full and demix and is ignored in --mode aggregate
``` Sample labels come from the sample ID, not from the file
names recorded inside the `.demix` files. Mixing results made with different barcode versions is allowed but the
plot then cannot be traced to a single dataset (the old directories are recorded as `unknown`).

### `report`: re-render the HTML report

```bash
nextflow run charlesfoster/wasteline -profile docker --mode report --results_dir results --outdir results_report
````

Needs `results/pipeline_info/run_info.json`. It reuses `summary/report_data.json` when present and otherwise rebuilds it from
the other files in `results/summary/`. Output: `results_report/summary/report.html` only. `--outdir` **must differ** from
`--results_dir` (report mode never writes into the results it reads). Use it after upgrading the pipeline to refresh the
report layout of an old run.

### `update_db`: download the current Freyja data and stop

```bash
nextflow run charlesfoster/wasteline -profile docker --mode update_db --db_dir /data/wasteline_db --outdir update_run
```

```
Freyja data updated: barcode version 10_06_2026-05-25 (2026-10-06T05-25), stored in /data/wasteline_db/freyja/2026-10-06T05-25
Use it with: --barcodes 2026-10-06   (or --barcodes 2026-10-06T05-25 for this exact update)
```

It runs `freyja update` in the pinned container and stores the dataset in `<db_dir>/freyja/<YYYY-MM-DDTHH-MM>/`
(about 16 MB). It writes no `run_info.json` (only `pipeline_info/freyja_data.json` and Nextflow's own reports). Running it
twice with the same upstream data overwrites the same directory.

## 6. Freyja data (barcodes)

Freyja's results depend on a barcode dataset (lineage definitions) that is updated almost daily. wasteline resolves it
**once per run** into a single dataset used by every sample, records exactly which one it was, and never silently falls
back to a different dataset. Freyja >= 2 barcodes are `usher_barcodes.feather`.

### Choosing a dataset: `--barcodes`

| Value                                    | Meaning                                                                                                    | `freyja_data.source` |
| ---------------------------------------- | ---------------------------------------------------------------------------------------------------------- | -------------------- |
| `latest` (default)                       | download the current dataset with `freyja update` (needs internet), store it under `--db_dir`, use it      | `latest`             |
| `container`                              | the snapshot bundled in the pinned Freyja container: fully offline, same data for a given pipeline release | `container`          |
| `YYYY-MM-DD`                             | the **newest** dataset of that day already stored in `--db_dir`; error if none                             | `pinned`             |
| `MM_DD_YYYY-HH-MM` or `YYYY-MM-DDTHH-MM` | one exact stored dataset (the time of day is part of the version)                                          | `pinned`             |
| a directory                              | a dataset directory (as `freyja update --outdir` writes it); used verbatim and checksummed                 | `custom`             |
| a `.csv` or `.feather` file              | a barcode file used verbatim (a bare file has no `lineages.yml`, so lineage grouping hints are not passed) | `custom`             |

How `latest` works: the process `FREYJA_UPDATE` runs `freyja update --outdir` in the Freyja container, the result is named
after the ISO form of its own `last_barcode_update.txt` (`MM_DD_YYYY-HH-MM`, for example `10_06_2026-05-25`), and is copied to
`<db_dir>/freyja/2026-10-06T05-25/`. Because the stamp has a time of day, two updates on one day can be two datasets.
Only `latest` downloads are stored in `--db_dir`; `container` and path datasets are not.

Gotchas:

- With `-resume`, `--barcodes latest` re-uses the **cached download** of the run being resumed and does not fetch again
  (verified: `update_db` with `-resume` reports `cached=2`). To refresh, run without `-resume` or run `--mode update_db`
  first and pin the result.
- Run-time download cannot change mid-run: the dataset is resolved before any sample is demixed.
- `--barcodes` is validated at start-up: a typo or a missing date stops the run before any work is done:

  ```
  [ERROR] No Freyja dataset from 2020-01-01 found in /data/wasteline_db/freyja. Available: [2026-10-06T05-25]. Download one with --barcodes latest, or point --barcodes at a dataset directory or .csv/.feather file.
  [ERROR] --barcodes 'nope' is not 'latest', 'container', a barcode version (YYYY-MM-DD, MM_DD_YYYY-HH-MM or YYYY-MM-DDTHH-MM) or an existing dataset directory / .csv / .feather file.
  ```

### Where datasets are stored: `--db_dir`

Default `~/.wasteline/db` (so datasets land in `~/.wasteline/db/freyja/<YYYY-MM-DDTHH-MM>/`). Set `--db_dir` to a shared,
**writable** directory to share datasets between users or runs. Each dataset directory holds `usher_barcodes.feather`,
`last_barcode_update.txt`, `lineage_mutations.json`, `curated_lineages.json`, `lineages.yml`, `pathogen_config.yml` and a
small `wasteline_dataset.json` with the download time. If the directory is not writable the run fails (see
[troubleshooting](troubleshooting.md#the-freyja-database-directory-is-not-writable)).

### Updating the barcodes before or during runs

Two equivalent patterns:

```bash
# A. update during the run (the default): every `full`/`demix` run with --barcodes latest fetches and stores the current data
nextflow run charlesfoster/wasteline -profile docker --input samplesheet.csv --outdir results --db_dir /data/wasteline_db

# B. update first (e.g. a nightly cron job on a node with internet), then pin by date, so runs are repeatable and offline
nextflow run charlesfoster/wasteline -profile docker --mode update_db --db_dir /data/wasteline_db --outdir /tmp/update_db
nextflow run charlesfoster/wasteline -profile docker --input samplesheet.csv --outdir results --db_dir /data/wasteline_db --barcodes 2026-10-06
```

### Offline use

`--barcodes container` uses the data bundled in the Freyja container and needs no network for the data (the images
themselves must be available locally or in a mirror: pull them once, or use `nf-core pipelines download`).
<!-- UNVERIFIED: `nf-core pipelines download` and a fully air-gapped run were not tried; a CI workflow (download_pipeline.yml) is meant to cover it. -->

Pinned dates and custom files are also offline. A stored dataset is a plain directory: copy `--db_dir/freyja/<stamp>/` to
the offline machine and use `--barcodes <stamp>` or `--barcodes <that directory>`. With `-profile docker` and no network,
`--barcodes latest` fails (verified by simulating no network inside the container):

```
urllib.error.URLError: <urlopen error [Errno -3] Temporary failure in name resolution>
```

### Reproducing the barcodes of an old run

Every result records its dataset: `pipeline_info/run_info.json` has a `freyja_data` block, and the barcode version is in
every `<sample>.demix_status.json`, `summary/demix_status.tsv`, `summary/barcode_version.json` and the report header.

```bash
python3 -c "import json; d=json.load(open('results/pipeline_info/run_info.json'))['freyja_data']; print(d['source'], d['barcode_version_iso'], d['barcodes_sha256'])"
```

- `latest` or `pinned`: use `--barcodes <barcode_version_iso>` while that directory still exists in `--db_dir`. Old versions
  **cannot be downloaded again** (`freyja update` only gives the current data), so archive the datasets you may need.
  Compare `barcodes_sha256` to confirm.
- `container`: use the same pipeline release (same Freyja container) with `--barcodes container`.
- `custom`: supply the same file or directory (`dataset_path`) and compare the checksum.

Then redo only the lineage estimate with `--mode demix --previous_results results --barcodes <that dataset>`. See
[output.md](output.md#how-to-reproduce-a-run) for the full recipe. The Freyja `--version` only ever reports the data
bundled in the package, so `freyja_reported_barcode_version` is filled only for `container`; the authoritative record is the
dataset's own `last_barcode_update.txt` and the checksums.

### Replaying a run: `bin/replay_run.py`

Step 0 of reproducing a run: `python3 bin/replay_run.py results --pipeline-dir /path/to/wasteline` reads
`results/pipeline_info/run_info.json` and the recorded parameters and prints the whole recipe: the pipeline version
(`git checkout <commit>`), a params file, the barcodes (`--barcodes <iso stamp>` while that dataset is still stored, the recorded
custom file, or `container`), a check of the recorded inputs (samplesheet, primer BED, reference, earlier results) by checksum, the
profile (`locked` is added when the containers were recorded with digests) and the command, with a **new** `--outdir`
(`results_replay` unless you give `--outdir`; the original is never overwritten). `--execute` runs it after printing it, but
only when no hard problem was found (`--force` overrides checksum and commit mismatches, never a missing file). Caveats: it
needs `--pipeline-dir` to leave out parameters that equal the schema defaults; stored Freyja datasets cannot be downloaded
again (it tells you when one has to be restored); uncommitted changes of the original checkout are not recorded;
`--compare_barcodes` and email options are deliberately not repeated. The script uses only the standard library and runs on the
host Python (3.8 or newer). The end-to-end test replays the full synthetic run and checks the demix files are byte-identical.

### Comparing barcode versions: `--compare_barcodes`

`--compare_barcodes <spec>` demixes every sample a **second** time against another dataset and puts the two estimates side by
side, so you can see how much a barcode update changes your lineage calls. `<spec>` takes the same values as `--barcodes`
(`latest`, `container`, a date or stamp already stored in `--db_dir`, a dataset directory, a `.csv`/`.feather` file). It works in
`full` and `demix` modes; the primary result (the `--barcodes` dataset) is unchanged, and the second dataset is resolved and
recorded independently.

```bash
# cheapest: reuse the variants and depths of an earlier run (here the synthetic test data; the primary
# barcodes are the test subset, the comparison dataset is the one bundled in the Freyja container)
nextflow run charlesfoster/wasteline -profile docker \
    --mode demix --previous_results results \
    --barcodes tests/data/barcodes_test.csv --compare_barcodes container \
    --outdir results_compare

# typical real use: how does today's dataset change last month's calls?
nextflow run charlesfoster/wasteline -profile docker \
    --mode demix --previous_results results --barcodes 2026-09-17 --compare_barcodes latest --outdir results_compare
```

(The first command was run and finished with `completed=63 failed=0`. The second form is the same code path with other
dataset specs; it was not run, because it needs a stored dataset and network access.)

What it writes:

- `summary/barcode_compare.tsv`: one row per sample and lineage, with `abundance_a` (the `--barcodes` dataset) and `abundance_b`
  (`--compare_barcodes`) and the barcode version of each. Real output for `MIX_A`:

  ```
  sample  lineage  abundance_a  abundance_b  barcode_version_a  barcode_version_b
  MIX_A   BA.2     0.712121     0.692308     custom             09_17_2026-00-46
  MIX_A   AY.4     0.287879     0.176471     custom             09_17_2026-00-46
  MIX_A   AY.1     0.000000     0.071429     custom             09_17_2026-00-46
  ```

  (The full container dataset splits `AY.4` into relatives that the 14-lineage test subset does not contain; that is the kind
  of change the table makes visible.)

- `summary/barcode_compare_status.tsv`: the outcome of the second demix per sample (`ok`, `skipped_low_coverage` or `failed`
  with a reason). A failed second demix is flagged, never fatal.
- `samples/<s>/freyja/<s>.compare.demix` (+ status and depth-stat files), the raw second demix.
- In `pipeline_info/run_info.json` a `freyja_data_compare` block (same fields as `freyja_data`), and
  `pipeline_info/freyja_data_compare.json`.
- A **Barcode comparison** panel in the report: per-sample turnover (half the summed absolute change in abundance) and the
  largest per-lineage changes ([report.md](report.md)).

In `aggregate`, `report` and `update_db` modes `--compare_barcodes` is ignored with a warning. A bad value
(`--compare_barcodes nope`) stops the run at start-up, like `--barcodes`; note that the message currently names `--barcodes`.
Reproducing a comparison later needs both datasets: archive them ([section 6](#reproducing-the-barcodes-of-an-old-run)).

## 7. Freyja settings

Every `--freyja_*` option maps to the Freyja option of the same name (`--freyja_eps` to `freyja demix --eps`, and so on):
`freyja_minq` and `freyja_varthresh` go to `freyja variants`; `freyja_eps`, `freyja_covcut`, `freyja_depthcutoff`,
`freyja_adapt`, `freyja_autoadapt`, `freyja_relaxedmrca`, `freyja_solver` and `freyja_confirmedonly` go to `freyja demix`;
`freyja_plot_mincov`, `freyja_plot_interval` (`MS` month, `W` week, `D` day) and `freyja_plot_windowsize` go to `freyja plot`;
`freyja_boot_n` and `freyja_boot_seed` go to `freyja boot` ([below](#bootstrap-confidence-intervals---freyja_boot)).
All defaults, ranges and descriptions are in [parameters.md](parameters.md) (the source of truth is `nextflow_schema.json`).

Two defaults deserve attention. `--freyja_varthresh` is **0.0** (every variant is kept; the old Snakemake pipeline discarded
variants below 25%, so abundances can differ from old results, see
[migration_from_snakemake.md](migration_from_snakemake.md)). `--freyja_covcut` is also the depth used by the coverage gate
below.

**Coverage gate.** Freyja will happily report a confident single lineage from a near-empty sample (the synthetic 0.6x
sample demixes to 100% BA.2 with Freyja's `coverage` 0.0), and crashes on a header-only variants file. wasteline therefore
computes the % of the genome covered at `--freyja_covcut`x from the depths file and, below `--min_genome_coverage_for_demix`,
records `skipped_low_coverage` with a reason instead of demixing. Do not set it to 0 unless you want that: with the gate
off the 0.0%-covered `LOWCOV_D` is reported as `ok`, and the near-empty `NEG_E` is recorded as `failed` with the Freyja
traceback as its reason (the run still completes).

**The plot.** `freyja plot` is broken in Freyja 2.0.3 to 2.0.5 containers (a pandas 3 `TypeError`). The pipeline tries the plain
command, retries with a documented workaround, and, if no PDF results, writes `freyja_plot.failed.txt` (and
`freyja_plot_times.failed.txt`) with the reason and **continues**. The report draws its own abundance plots and does not
depend on the PDFs. In a normal run the workaround succeeds; the failure note also appears when no sample exceeds
`--freyja_plot_mincov`:

```
freyja plot produced no PDF (plain attempt and infer_string workaround both failed).
--- attempt 1 (freyja plot), last lines ---
ERROR: No samples matching coverage requirements, so no plot will be generated. Try changing --mincov threshold.
```

The workaround is a stopgap that will be removed once Freyja fixes the upstream bug.

### Bootstrap confidence intervals: `--freyja_boot`

`--freyja_boot` runs `freyja boot` for every sample that was demixed (including a negative control that demixed) and gives a
95% interval for each lineage abundance. It works in `full` and `demix` modes and is off by default (one more task per sample).
`--freyja_boot_n` sets the number of bootstraps (default 100) and `--freyja_boot_seed` the seed (default 0, fixed, so a
re-run gives identical intervals). The bootstrap uses the same barcodes and the same `--freyja_eps`, `--freyja_depthcutoff`,
`--freyja_solver` (and `--freyja_relaxedmrca`/`--freyja_confirmedonly`) as the demix. Example (adds the bootstrap to an
existing run; real output on the synthetic data):

```bash
nextflow run charlesfoster/wasteline -profile docker --mode demix --previous_results results \
    --barcodes tests/data/barcodes_test.csv --freyja_boot --outdir results_boot
```

```
sample  lineage  abundance  lo        hi
MIX_A   BA.2     0.712121   0.670435  0.727273
MIX_A   AY.4     0.287879   0.269521  0.304565
MIX_A   B.1      0.000000   0.000000  0.019677
MIX_B   BA.1     0.500000   0.474049  0.528547
```

This is `summary/freyja_boot_summary.tsv` (`lo` and `hi` are the 2.5 and 97.5 percentiles of the bootstrap replicates; `abundance`
is the demix estimate). The per-sample raw tables (`<s>.boot_lineages.csv`, `<s>.boot_summarized.csv`, `<s>.boot_summary.tsv`)
and a `<s>.boot_status.json` are in `samples/<s>/freyja/`; the report shows the intervals in the abundance tooltips. A
failed bootstrap is recorded in its status file (never fatal) and that sample simply has no interval. In `aggregate` mode the
option is ignored with a warning. The intervals describe sampling noise in the reads given the barcodes; they do not include
barcode-dataset uncertainty (use [`--compare_barcodes`](#comparing-barcode-versions---compare_barcodes) for that).
<!-- UNVERIFIED: only the defaults (--freyja_boot_n 100, --freyja_boot_seed 0) were run; other values were not tried. -->

### Time series, growth rates and the dashboard

When at least some samples have a `collection_date`, a time-aware plot (`summary/freyja_plot_times.pdf`, `freyja plot --times`;
`--freyja_plot_interval MS|W|D` and `--freyja_plot_windowsize`) and the report's time-series panel (needs at least two distinct
dates) are produced from the dated, non-control samples. Dates come from the samplesheet in `full`, `demix` and `aggregate` modes
(`--previous_results` carries them over from the earlier run's `summary/samplesheet.csv`; with `--demix_dir` there is no metadata).
The report panel shows the mean abundance of the top summarized groups per date and per `site`
([report.md](report.md#time-series)).

**Time-series tables** (`--timeseries`, on by default, all modes; written when there are at least two distinct collection dates
among the non-control samples): `summary/freyja_timeseries.tsv` (summarized lineage groups) and
`summary/freyja_timeseries_lineages.tsv` (every lineage Freyja reported), both long tables with the columns
`site collection_date lineage_group abundance abundance_load_weighted n_samples mean_viral_load total_viral_load`. `abundance` is the
mean over the samples of that site and date (a group a sample lacks counts as 0); `abundance_load_weighted` is the mean weighted by
`viral_load` (`NA` when no sample has one). Negative controls, samples without a date, samples below `--freyja_plot_mincov`
genome coverage and samples missing from the samplesheet are left out (the process log counts each reason).

**Relative growth rates** (`--freyja_growth`, off by default) run `freyja relgrowthrate` on the last `--freyja_growth_days` (56) days:
`summary/freyja_growth_rates.csv` (lineage, estimated relative advantage, bootstrap 95 % interval) and
`summary/freyja_growth_inputs.json` (what went into the fit). The fit needs enough data: at least `--freyja_growth_min_dates` (4)
distinct dates spanning at least `--freyja_growth_min_span_days` (14) days in the window. With less, there is no table but a
`summary/freyja_growth_rates.failed.txt` saying why (never a failed run). The bootstrap is seeded (`--freyja_growth_seed`) so a
repeated run gives the same intervals. `--freyja_growth_scale_by_viral_load` weights the fit by `viral_load`; samples without a
load are then dropped from it. On the synthetic series in `tests/data_timeseries` (a lineage rising with a known logistic growth)
the estimate is positive for the rising lineage; treat the numbers as indications of direction, not as precise estimates, unless
your series is long and your samples comparable.

**Dashboard** (`--freyja_dash`, off by default): `summary/freyja_dashboard.html`, Freyja's interactive dashboard, with
`--freyja_dash_title` and `--freyja_dash_intro` (inserted as HTML). Plotly is embedded, but the Bootstrap styling and a logo are
loaded from the internet when the page is opened, so it looks plain offline. The same data requirements as the growth rates apply
(otherwise `summary/freyja_dashboard.failed.txt`). The report links the dashboard file by name; it does not embed it.

## 8. Consensus genomes

One consensus FASTA per **non-negative** sample (`summary/consensus_genomes.fasta`, one record per sample, header = sample ID;
`summary/consensus_stats.tsv`). It is built from the same `freyja variants` table (ivar-based) that Freyja uses: there is no second
variant caller. Because `freyja_varthresh` is 0.0, the table holds every alternate allele and the thresholds are applied
**here**:

| Rule                                                                                                                      | Parameter (default)                                                      |
| ------------------------------------------------------------------------------------------------------------------------- | ------------------------------------------------------------------------ |
| A call needs total depth at the site >=                                                                                   | `--consensus_min_depth_variants` (10)                                    |
| and alt base quality >=                                                                                                   | `--consensus_min_alt_qual` (20)                                          |
| SNV with ALT frequency >= this is written as the alt allele                                                               | `--consensus_freq` (0.75)                                                |
| SNV with frequency above this and below `consensus_freq`, and enough alt reads, is written as an **IUPAC ambiguity code** | `--consensus_iupac_min_af` (0.25), `--consensus_iupac_min_alt_reads` (5) |
| SNV at or below `consensus_iupac_min_af` is ignored (reference base)                                                      |                                                                          |
| indel with frequency >= this is applied; below it ignored, **never** IUPAC                                                | `--indel_freq` (0.75)                                                    |
| positions with depth below this (from the depths file, variant sites exempt) become `N`                                   | `--consensus_min_depth` (10)                                             |

The defaults are listed in [parameters.md](parameters.md). `--consensus_iupac false` switches ambiguity codes off (the reference base is left instead). `--skip_consensus`
disables the whole consensus chain and its summaries. Deletions are written as `-`. `consensus_stats.tsv` gives
`length, n_count, n_pct, ambiguous_count` and a `status`: `ok` (including a clean reference-identical sample), `empty` (no
usable depth; an all-N sequence of reference length) or `failed`. A consensus is never an empty file, and a sample whose
consensus steps failed still gets a `failed` row. Do not use `empty`/high-`n_pct` consensuses for sequence-based work.
The consensus depth mask is different from the QC depth (`--qc_min_depth`, below): the same sample can PASS QC and still have
a few percent `N`.

### Annotating variants: `--annotate_variants`

`--annotate_variants` lists the amino-acid consequences of the variants that were applied to each consensus (a `bcftools csq`
port of the old pipeline's `amino_acid_consequences` rule, using the GFF3 annotation). It works in `full` and `demix` modes. It
runs on the consensus filter's VCF, i.e. exactly the calls that pass the consensus thresholds above (including those written as
IUPAC codes), so it **needs the consensus chain: combining it with `--skip_consensus` stops the run at start-up**:

```
[ERROR] --annotate_variants annotates the variants the consensus is built from: it cannot be combined with --skip_consensus
```

```bash
nextflow run charlesfoster/wasteline -profile docker -resume \
    --input samplesheet.csv --outdir results --annotate_variants        # or add it to an existing run with --mode demix
```

Output: `samples/<s>/variants/<s>.annotated.tsv` per sample and the concatenation `summary/annotated_variants.tsv` (no negative
controls: they have no consensus). Real output (first rows for `MIX_A`, synthetic data):

```
reference    sample_id  gene    nt_pos  ref  alt  effect      aa_bcsq      aa_standard  vaf       depth
NC_045512.2  MIX_A      ORF1ab  670     T    G    missense    135S>135R    S135R        0.692308  65
NC_045512.2  MIX_A      ORF1ab  3037    C    T    synonymous  924F         F924F        1         82
NC_045512.2  MIX_A      ORF1ab  4181    G    T    missense    1306A>1306S  A1306S       0.3       60
```

`aa_bcsq` is bcftools' notation (`135S>135R`), `aa_standard` the usual one (`S135R`); a variant without a coding consequence has
`effect` `intergenic` and `NA` gene/amino acid. Only the first consequence of each variant is reported (as in the old pipeline).
For a chosen list of mutations (including amino-acid changes), [mutations of interest](#10-mutations-of-interest) reports the
frequency and depth in every sample, whether or not the mutation was called.

## 9. Negative controls and QC thresholds

Mark negative controls with `sample_type=negative` (`make_samplesheet.py` does this for file names matching `NEG_|NC_`). Then:

- they are **demixed** if they pass the coverage gate, so contamination is visible, but they are **excluded from the aggregate
  and `freyja plot`** unless `--include_negatives_in_aggregate` (the log says so: `Negative control NEG_E is not included in the
aggregate and plot ...`). What they contain is in their own table, `summary/freyja_aggregated_controls.tsv` (the same format
  as `freyja_aggregated.tsv`), and in the report's **Negative controls** section (QC, reasons, mapped reads, genome coverage,
  lineages detected). A control with no coverage (clean, as it should be) is `skipped_low_coverage` and has no row in that table;
- a control that demixed `ok` stays `ok` in `demix_status.tsv` and in the report: being left out of the aggregate does not make
  it a failure;
- they are judged on reads, not on coverage: `qc_status` is `FAIL` with reason `ntc_reads` if mapped reads exceed
  `--ntc_max_mapped_reads` (1000) **or** more than `--ntc_max_genome_coverage_pct` (10%) of the genome is covered at >= 1x;
  otherwise `PASS`. A failed negative control is a **run-level error flag** in the report;
- they get no consensus genome (and no variant annotation).

Real example (synthetic data): `NEG_E` (6 reads) is `skipped_low_coverage` and QC `PASS`; `NEG_CONTAM` (a negative control that
carries BA.2 reads) is demixed (`ok`), QC `FAIL` (`ntc_reads`), listed in `freyja_aggregated_controls.tsv` with `BA.2 1.0` and absent
from `freyja_aggregated.tsv`, and raises the error flag
`Negative control NEG_CONTAM has reads above the allowed limit (5998 mapped reads, limit 1000): possible contamination.`
Interpretation: a failed control means the run's samples may carry contamination; read the lineages in the controls table to see
what contaminated it.

QC for ordinary samples uses the depths file (bases with Q >= 20), samtools flagstat and the fastp JSON:

- `--qc_min_depth` is the depth used for the coverage test; `PASS` needs at least `--qc_pass_min_pct` of the genome at that
  depth, below `--qc_fail_below_pct` the sample is `FAIL` (`low_coverage`), in between `WARN`.
- An amplicon is _dropped_ when its mean depth is below `--amplicon_min_depth`; more than `--qc_max_dropped_amplicons`
  dropped amplicons gives `WARN` (`amplicon_dropout`). The default is strict for 96/99-amplicon ARTIC schemes: raise it
  (the schema suggests 10 for them).
- `--min_genome_coverage_for_demix` is the coverage gate of [section 7](#7-freyja-settings).

Defaults are in [parameters.md](parameters.md).

Other reason tokens: `no_reads`, `no_depths`. The columns of `qc_summary.tsv` are in [output.md](output.md). A QC `FAIL` does not
stop a sample being demixed (only the coverage gate does).

### Dropout impact

A dropped amplicon matters most when it hides mutations that tell the detected lineages apart. `summary/dropout_impact.tsv` has
one row per dropped amplicon of each sample (modes `full`, and `demix` with an explicit primer scheme; `--dropout_impact false`
turns it off). For every dropped amplicon it counts the _informative_ sites inside it (single-nucleotide substitutions carried
by some but not nearly all lineages of the Freyja dataset: fixed mutations cannot separate lineages) and how many belong to
lineages detected in that sample (abundance of at least `--dropout_min_abundance`, 0.01). `affects_detected_lineage` is `true`
when a detected lineage carries at least one informative site in the dropped span, and `affected_lineages` names them with
the number of sites (`XBB.1.5:9,BA.2:3`). `summary/dropout_impact_summary.tsv` has one row per sample: how many informative
sites the scheme can see at all, how many were lost (`pct_informative_sites_lost`) and in how many dropped amplicons it matters.
Real example (synthetic `Site_1_week_2`, XBB.1.5 + BA.5): amplicon `nCoV-2019_22` (spike) hides 28 informative sites, 10 of them
in the detected lineages (`XBB.1.5:9,BA.2:3,BA.5:3`), so the abundance estimate may be biased; `nCoV-2019_7` hides 5 sites, none
in a detected lineage, so it does not matter for the lineage estimate. The report lists the affected samples and raises the flag
`dropout_affects_lineages`. Limits: only substitutions are considered (the dataset's mutation lists contain no indels);
it needs the dataset's `lineage_mutations.json` or a barcodes `.csv`, so with a custom `.feather`-only `--barcodes` the table is
skipped with one warning; it says which sites are hidden, not by how much the estimate moves.

### Unexplained mutations

`summary/unexplained_mutations.tsv` lists frequent single-nucleotide variants in a sample that the demixed lineage mixture does not
predict: for each variant the expected frequency is the sum over the detected lineages of abundance times "lineage carries the
mutation" (from the dataset); the variant is reported when it has depth of at least `--unexplained_min_depth` (20), at least
`--unexplained_min_alt_depth` (10) supporting reads, an observed frequency of at least `--unexplained_min_freq` (0.05), exceeds
the expectation by at least `--unexplained_min_excess` (0.05), and passes a binomial guard against depth noise
(`--unexplained_max_pvalue`, `0.05 / 29903`; `1` disables it). Columns include the gene and amino-acid change, and
`n_samples`/`recurrent`: the same mutation unexplained in at least `--unexplained_min_recurrence` (2) samples is flagged
(recurrence across samples or sites is the more convincing signal). Per-sample tables are in
`samples/<s>/variants/<s>.unexplained_mutations.tsv`. Modes `full` and `demix`; negative controls are not analysed.

**What it is and is not.** It is a screening aid: an early hint of a lineage missing from the barcodes, a cryptic lineage or a
recurrent artefact. It is **not a variant caller** and says nothing about lineage assignment. Sequencing error, primer
effects, a lineage too rare to be demixed and mixtures that the barcodes represent badly all produce rows; confirm anything
interesting on the reads. On the synthetic samples at the default thresholds the table is empty (a deliberately spiked sample
in the test data is detected, with a p-value many orders of magnitude below the cut-off).

### Alignment format and aligner: `--bam_format`, `--aligner`

`--bam_format cram` publishes a reference-based CRAM (`samples/<s>/bam/<s>.primertrimmed.cram` + `.crai`) instead of the BAM
(needs `--save_bam`, the default). CRAM is much smaller; reading it needs the reference FASTA that made the run (it is recorded
in `run_info.json`), `samtools view -T reference.fa`. Everything else (Freyja, QC, consensus) still uses the BAM inside the run,
so results do not change. `--aligner minibwa` maps with minibwa instead of `bwa mem`; it is **experimental**: `bwa mem` stays the
default and the validated choice, and the evaluation (speed is irrelevant at 30 kb, the tool is a young 0.x release with a preprint
only) is in [aligner_evaluation.md](aligner_evaluation.md). `tests/compare_aligners.sh` runs both on the synthetic data and
compares abundances and consensus genomes.

## 10. Mutations of interest

Opt-in: `--mutations_of_interest mutations.csv` reports the alternate-allele frequency and depth of listed mutations in every
sample (`summary/mutations_of_interest.tsv` and a report heatmap). CSV with a header `name,mutation`:

```csv
name,mutation
S N501Y,S:N501Y
S N501Y (nt),A23063T
S L452R,S:L452R
```

`mutation` is a nucleotide SNV (`A23063T`: reference base, 1-based position, alternate), a deletion by its first deleted base
(`C21995-`, `CTT21991-`), an insertion after a position (`T22204+GAGCCAGAA`), or an amino-acid change `GENE:RefPosAlt` (`S:N501Y`;
`*` = stop; genes from the GFF3, including `ORF1a`/`ORF1ab`). A mutation that is absent is reported with `alt_depth 0`,
`alt_freq 0`; an unparseable entry is an error. Real output (columns `sample name mutation ref_depth alt_depth depth alt_freq`):

```
MIX_A   S N501Y   S:N501Y   11  41  52  0.788462
MIX_B   S L452R   S:L452R    22   6  28  0.214286
PURE_C  S N501Y   S:N501Y    0  40  40  1
```

## 11. Resuming a run: `-resume`

Nextflow caches every task in `work/`. After an interruption, or to add a sample you have not changed anything else
about, repeat the **identical** command with `-resume`:

```bash
nextflow run charlesfoster/wasteline -profile docker --input samplesheet.csv --outdir results -resume
# [SUCCESS] completed=3 failed=0 cached=109     (the provenance and report steps always re-run: they record the run time)
```

- `-resume` continues from the **most recent run in that launch directory** (see `nextflow log`). If you have run other
  commands in between, give the session id of the run you want: `nextflow log -f name,session` then `-resume <session-id>`
  (verified). The run-name form of `-resume` did not pick up that run's cache in a test on 26.04.6.
- A re-run in the same `--outdir` **without** `-resume` recomputes everything.
- Changing a parameter that feeds a process invalidates that process and everything downstream, nothing else.
- `--barcodes latest` + `-resume` re-uses the earlier download (see [section 6](#6-freyja-data-barcodes)).
- Deleting `work/` removes the cache.

## 12. Profiles

`-profile` selects the software and environment; combine with commas (later profiles win): `-profile docker,test`.
Only **docker** was exercised end-to-end for this guide (macOS/Apple Silicon); the others parse (`nextflow config -profile <name>`)
and are set up to work, and the release CI matrix (`.github/workflows/nf-test.yml`) is meant to cover singularity and conda
<!-- UNVERIFIED: singularity, apptainer, podman, conda and mamba profiles were only checked with `nextflow config`, not run; the CI matrix has not been run on GitHub -->.

| Profile                    | Use                                                                                                                                                                                                       |
| -------------------------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `docker`                   | Docker. Adds `--platform=linux/amd64` automatically on arm64 hosts.                                                                                                                                       |
| `singularity`, `apptainer` | Singularity / Apptainer (with auto-mounts).                                                                                                                                                               |
| `podman`                   | Podman.                                                                                                                                                                                                   |
| `conda`, `mamba`           | Conda/mamba environments (`environment.yml` per module). Last resort; not reliable on macOS.                                                                                                              |
| `locked`                   | Use with an engine profile (`-profile docker,locked`): pin every process to an immutable image digest (`conf/containers.lock.config`, generated by `bin/lock_containers.py`). No effect with conda/mamba. |
| `emulate_amd64`            | Force `--platform=linux/amd64` on Docker on any host.                                                                                                                                                     |
| `test`                     | Tiny synthetic dataset from `tests/data` and `--barcodes container`; offline. Launch from the repository root.                                                                                            |
| `test_stub`                | Same inputs with `--barcodes latest` for `-stub-run` wiring checks (no download happens in a stub run).                                                                                                   |
| `test_full`                | Placeholder: the `test` data with `--barcodes latest`.                                                                                                                                                    |
| `debug`                    | Keeps work files and dumps task hashes.                                                                                                                                                                   |

Other generic profiles inherited from the nf-core template (`shifter`, `charliecloud`, `wave`, `arm64`, `gpu`) are present but not used by
the pipeline. Do not put pipeline parameters in a custom `-c` config file (they are validated against the schema); use
`--param` or `-params-file`. Use `-c` only for resources, executors and similar infrastructure settings.

## 13. Resources

Every process requests CPUs, memory and time from its label (`process_single`: 1 CPU, 4 GB; `process_low`: 2 CPUs; `process_medium`:
4 CPUs, 8 GB; `process_high`: 8 CPUs, 16 GB; `FREYJA_DEMIX` and `FREYJA_VARIANTS` are `process_low`), scaled up on
retry. `process.resourceLimits` caps every request: by default at the CPUs and 90% of the memory of the machine, and 240 h. To
change the caps, or to tune single processes, use a small config file passed with `-c` (**not** pipeline parameters):

```groovy
// resources.config
process.resourceLimits = [ cpus: 8, memory: '32.GB', time: '48.h' ]
process {
    withName: 'BWA_MEM'     { cpus = 8; memory = '24.GB' }
    withName: 'FREYJA_DEMIX' { cpus = 4 }
}
```

```bash
nextflow run charlesfoster/wasteline -profile docker -c resources.config --input samplesheet.csv --outdir results
```

(Verified with `nextflow -c resources.config config`: the `resourceLimits` block becomes `cpus = 4`, `memory = '16.GB'`,
`time = '24.h'` for a test file.) The profiles `test` and `test_stub` set smaller limits of their own.

## 14. Parameter files and boolean flags

Put repeated settings in a YAML or JSON file and pass it with `-params-file` (single dash; it is a Nextflow option):

```yaml
# params.yaml
input: /data/run7/samplesheet.csv
outdir: results/run7
primer_scheme: midnight
barcodes: 2026-10-06
db_dir: /data/wasteline_db
save_bam: false
skip_consensus: true
```

```bash
nextflow run charlesfoster/wasteline -profile docker -params-file params.yaml
```

(Verified: with `save_bam: false` and `skip_consensus: true` no BAM and no consensus files were produced.) Command-line parameters
override the file. Each run also writes `pipeline_info/params_<timestamp>.json` with every parameter of the run, which is
valid as a `-params-file` for repeating it (edit `outdir`).

**Boolean flags** work as `--skip_consensus`, `--skip_consensus true`, `--skip_consensus false`, `--skip_consensus=false`,
and in a `-params-file`; this needs Nextflow >= 26.04. A value that is not a boolean is rejected:
`Parameter `skip_consensus` with type Boolean cannot be assigned to maybe [String]`. A misspelt parameter is rejected too:
`Parameter `primer_sheme` was specified on the command line or params file but is not declared in the script or config`.
`nextflow run charlesfoster/wasteline --help` lists every option (use `--help_full` / `--show_hidden` for more).

## 15. Running on an HPC cluster

Generic advice (not cluster-specific, and untested here):

- Use `-profile singularity` (or `apptainer`) and an executor config. Put your site's settings in a `-c` file
  (`process.executor = 'slurm'`, `process.queue = ...`, `singularity.cacheDir = ...`); never parameters.
- Run the Nextflow "head" job where it survives your logout: `nextflow run ... -bg`, `tmux`/`screen`, or a long-running
  scheduler job. Cap its memory with `export NXF_OPTS='-Xms1g -Xmx4g'`.
- Compute nodes often have no internet. Run `--mode update_db --db_dir <shared dir>` on a node that has it, then use
  `--barcodes <date>` (or `container`) in the real runs; make sure the container images are in a shared cache (set
  `NXF_SINGULARITY_CACHEDIR`/`singularity.cacheDir`, and pull them first on a node with access to the registries).
- `--db_dir` and `--outdir` must be on storage that all nodes can reach and that you can write to.
- Reproducible installs: `-r <release>` pins the pipeline version, `-profile <engine>,locked` pins every container by digest.
- Keep the work directory on fast scratch storage and clean it when finished.

## 16. Updating the Freyja version of the pipeline

The Freyja container is pinned to one exact tag per wasteline release (currently `quay.io/biocontainers/freyja:2.0.6--pyhdfd78af_0`, in every
`freyja_*` module and its `environment.yml`); every other container is pinned too, none uses `latest`. Freyja is the one tool
that is bumped automatically:

- A scheduled GitHub workflow (`.github/workflows/freyja-update.yml`, daily) runs `bin/bump_freyja.py --check`; when a newer
  release exists on quay.io, bioconda and the Galaxy singularity depot, it opens a pull request that rewrites all pins and
  refreshes the digest lock. The normal CI (lint, pytest, nf-test on the synthetic data) must pass; same-major bumps then
  auto-merge, a new major version waits for a human. Setup and operation: [ci.md](ci.md).
- By hand, from a checkout (standard library Python, no extra installs):

  ```bash
  python bin/bump_freyja.py --check              # pinned vs latest; exit 10 = update available
  python bin/bump_freyja.py --apply --dry-run    # what would change
  python bin/bump_freyja.py --apply              # rewrite pins + a CHANGELOG line
  python bin/bump_freyja.py --apply --tag 2.0.4--pyhdfd78af_0   # an exact tag, e.g. to roll back
  python bin/bump_freyja.py --require-consistent # fail if any Freyja pin disagrees
  ```

  ```
  Freyja pins consistent (51 pins)
  pinned : 2.0.5--pyhdfd78af_0
  latest : 2.0.5--pyhdfd78af_0
  checks : quay=ok, bioconda=ok, singularity_depot=ok
  Freyja is up to date.
  ```

- Barcode data is **not** part of this: it is data, chosen per run with `--barcodes` ([section 6](#6-freyja-data-barcodes)).
- Users of a released pipeline get the Freyja version of that release: pick the pipeline version with `-r`, update the
  cached copy with `nextflow pull charlesfoster/wasteline`, and read `freyja_version` in `run_info.json`. The version of
  every tool used is also in `pipeline_info/software_versions.yml`.

## 17. Not implemented

These are in [PLAN.md](PLAN.md) or discussed as future work but are **not** available (they are not in
[parameters.md](parameters.md)): the optional `--kraken2_db` screen of negative controls, other pathogens than SARS-CoV-2
(`--pathogen` accepts only `SARS-CoV-2` and rejects anything else at start-up), CRAM output, resolving historical barcode
versions from the Freyja git history, and an alternative aligner. See [development.md](development.md#future-work).

## Where next

[output.md](output.md) (every file), [report.md](report.md) (reading the HTML report), [parameters.md](parameters.md),
[troubleshooting.md](troubleshooting.md), [migration_from_snakemake.md](migration_from_snakemake.md),
[development.md](development.md).
