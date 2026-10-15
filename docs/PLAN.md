# wasteline: design and implementation plan

Status: approved by the user on 2026-10-15 (scope, feature list F1-F20 with the F17 change below, decisions D1-D4).
Owner of this document: the planner (Opus). `coder` implements phase by phase; the planner reviews each phase before the
next begins.

## 1. Scope

Port the **wastewater** workflow of `charlesfoster/covid-illumina-snakemake` (v0.10.x) to Nextflow, nf-core style.
SARS-CoV-2 only for now; structure so other pathogens (Freyja `--pathogen`) can be added later (F15).

In scope: Illumina amplicon reads -> QC -> primer-trimmed BAM -> Freyja variants/depths -> demix -> aggregate/plot ->
HTML report; per-sample consensus genomes; negative-control assessment.
Out of scope (dropped): the routine/clinical workflow, pangolin, nextclade, pango-collapse, `--legacy` QC columns,
qualimap, kraken2 as a default, snpEff.

## 2. Decisions

| # | Decision |
|---|---|
| D1 | Use the latest Freyja release (currently biocontainer `quay.io/biocontainers/freyja:2.0.5--pyhdfd78af_0`; verify at implementation time). A scheduled CI job opens a PR on each new release; it auto-merges only if CI passes (see §9). |
| D2 | Keep `bwa mem` for parity. Future: compare against minibwa (see §11). |
| D3 | No validation dataset available. We generate a synthetic one (see §8). |
| D4 | Production: Linux + Docker. Profiles: `docker`, `singularity`, `conda`/`mamba` (+ `test`, `test_stub`). Dev: macOS + Docker. |
| D6 | `freyja plot` is broken in Freyja 2.0.3-2.0.5 containers (pandas 3.0.x `TypeError` in `prepSummaryDict`; exits 0 with no PDF). Verified by the planner on Freyja's own bundled example data. Decision: keep the single pinned Freyja container (no custom pandas<3 env, which would break the one-container auto-bump). `FREYJA_AGGREGATE_PLOT` tries `freyja plot` normally; if no PDF appears it retries through `python -c` with `pd.options.future.infer_string=False` set before calling `freyja._cli.cli()`; if that also yields no PDF it writes `freyja_plot.failed.txt` (with the reason) and the run continues. The same check applies when `freyja plot` exits 0 without a PDF because no sample exceeds `--mincov`. Our HTML report draws its own abundance plots and never depends on `freyja plot`. The workaround is a documented stopgap: a test (`freyja plot` on the bundled example, run without the workaround) flags when upstream is fixed so it can be removed. Upstream issue to be reported by the user. |
| D7 | Coverage QC is computed from the Freyja **depths file** (4-col, per position, Q>=20 bases) plus samtools flagstat and the fastp JSON, not from mosdepth. Reason: the same QC then works in `demix` mode, where there is no BAM, and it is one fewer tool/container. This supersedes the mosdepth wording in §5.9 and F12's "mosdepth or samtools coverage". |
| D8 | Negative controls (`sample_type=negative`) are demixed if they pass the coverage gate (so contamination is visible in the demix-status table/report) but are **excluded from the aggregate and `freyja plot`** by default; `--include_negatives_in_aggregate` (default false) overrides. A negative that passes the gate (i.e. has real coverage) is also a run-level error flag (`ntc_reads`). |
| D10 | Legacy zero-byte `.demix` files (the old pipeline's `|| touch` on failure) are never aggregated: `make_samplesheet.py --from-results` skips them with a reason, and the workflow records the sample in `demix_status.tsv` with status `missing` and reason `empty demix file (failed in the old pipeline)`. |
| D11 | Boolean CLI params break on Nextflow >= 26.04 with the strict parser (strings, not booleans). Handled in the integration round (typed `params` block or an explicit actionable error); see CLAUDE.md gotchas. |
| D12 | Minimum Nextflow is **26.04.0** (was 25.10.4). Boolean/number CLI params only convert correctly through a typed `params {}` block under Nextflow >= 26.04's strict parser; on 25.10 included scripts and config saw unconverted strings (e.g. `--save_bam false` still published the BAM, silently), so older versions are refused rather than mis-run. h2seq and monty already use 26.04. Made by the integration agent; reversible if production Linux hosts cannot run 26.04. |
| D13 | **Probe-capture data (e.g. Twist) is supported** via `--protocol amplicon|capture` (default `amplicon`). `capture`: no primer trimming (IVAR_TRIM skipped; reads go straight from alignment to sorted/indexed BAM), no amplicon-dropout analysis (amplicon outputs/section omitted, QC dropout rule not applied), duplicate handling per `--duplicates auto|mark|none` (`auto` = mark for capture, none for amplicon; samtools excludes DUP-flagged reads from depth/variant counting by default). `--skip_primer_trimming` also exists for pre-trimmed amplicon data. Primer scheme/BED are ignored with a warning in capture. The protocol is recorded in run_info. Synthetic capture test data live in `tests/data_capture/` (tests/data is unchanged). |
| D9 | `freyja demix --version` only ever reports the data bundled in the package, so it cannot confirm an external dataset. The authoritative barcode record is the dataset's own `last_barcode_update.txt` + sha256s. `freyja_reported_barcode_version` is filled only for `source=container`; otherwise null, with the container's bundled version in `freyja_bundled_barcode_version`. This supersedes the "parsed from demix --version against that exact dataset" wording in §5.13. |
| D5 | Per-sample consensus genomes are kept (F17 changed), built from the ivar-based Freyja variants TSV (no LoFreq). pangolin/nextclade/pango-collapse are dropped. |

## 3. Pipeline modes

`--mode full|demix|aggregate|report|update_db` (schema enum, default `full`). `--mode` is the only dispatch: named
entry workflows (`-entry`) are rejected by Nextflow >= 26.04's strict parser, so they were dropped.

| Mode | Input | Runs | Does NOT run |
|---|---|---|---|
| `full` | FASTQ samplesheet | everything | |
| `demix` | samplesheet of existing `variants` + `depths` (or `--previous_results <dir>`) | resolve Freyja data -> `freyja demix` -> aggregate/plot -> report | fastp, bwa, ivar, variants, consensus |
| `aggregate` | existing `.demix` files (samplesheet, `--demix_dir`, or several results dirs) + optional metadata | aggregate (+ time-series plot) -> report | demix and everything upstream |
| `report` | `--results_dir` of a previous run | re-render `report.html` only | everything else |
| `update_db` | none | download Freyja data into `--db_dir/freyja/<barcode-date>/` and exit | everything else |

`demix` replaces the Snakemake `--redo_demix` (which deleted `*.demix` and re-ran). It must work from a results
directory produced by this pipeline **and** by the old Snakemake pipeline (legacy layout:
`<s>/freyja/<s>.depths`, `<s>/freyja/<s>.freyja_variants.tsv`), independent of any Nextflow work directory. It never
modifies `--previous_results`; output goes to `--outdir`.

## 4. Inputs

nf-schema validated (`assets/schema_input.json` etc.). IDs match `^[A-Za-z0-9._-]+$`, unique, never rewritten.

- `full`: `sample,fastq_1,fastq_2[,sample_type][,collection_date][,site]`. `sample_type` in `sample|negative`
  (default `sample`). `collection_date` is ISO `YYYY-MM-DD`.
- `demix`: `sample,variants,depths[,fastp_json][,sample_type][,collection_date][,site]`.
- `aggregate`: `sample,demix[,sample_type][,collection_date][,site]`.
- Reference/primers: `--genome` default bundled NC_045512.2 (+ GFF) from `assets/references/` (copy from the Snakemake
  repo `CIS/bin/references/`); `--primer_scheme midnight|midnightv2|articv3|articv4.1|articv5.3.2|eden|swift|qiagenv1|xgen`
  (bundled BEDs from `CIS/bin/primer_schemes/`, default `midnight`) or `--primer_bed <path>`. Validate that the BED
  contig name matches the reference contig. `swift` implies ivar `-x 5`.
- Freyja data: `--barcodes latest|<YYYY-MM-DD>|<path to dir or .csv/.feather>` (default `latest`);
  `--db_dir` (default `${HOME}/.wasteline/db` unless set; must be writable).

### Helper scripts (`bin/`, stdlib + pandas only, unit-tested with pytest)

- `make_samplesheet.py <fastq_dir> [--suffix _L001_R1_001.fastq.gz] [--negative-pattern 'NEG_|NC_'] [--metadata meta.csv] -o samplesheet.csv`:
  port the Snakemake discovery (`_S\d+_L001` stripping; R2 derived from R1; `Undetermined` excluded; `NEG_`/`NC_` ->
  `sample_type=negative`). Fails loudly on missing R2 or ID collisions.
- `make_samplesheet.py --from-results <dir> --mode demix|aggregate -o sheet.csv`: discovers new and legacy layouts.
- `check_samplesheet.py`: fail-fast validation, stdout JSON + stderr table (pattern from quasispecies pipeline).
- `bump_freyja.py`: see §9.

## 5. Workflow (mode `full`)

1. `INPUT_CHECK` (nf-schema + custom checks).
2. `FASTP` (nf-core module). Parity params: `--detect_adapter_for_pe --cut_front --cut_tail --trim_poly_x
   --cut_mean_quality 20 --qualified_quality_phred 20 --unqualified_percent_limit 10 --correction --length_required 50`,
   all exposed as `fastp_*` params. Trimmed reads not published unless `--save_trimmed_reads`.
3. `BWA_INDEX` once (cached by reference), `BWA_MEM` -> sort -> `samtools view -F 4` -> index.
4. `IVAR_TRIM` (`-b <bed> -x <offset> -e`, sort, index). Publish the primer-trimmed BAM (`--save_bam`, default true;
   CRAM is a future option).
5. `FREYJA_VARIANTS`: `freyja variants <bam> --ref --variants --depths --minq 20 --varthresh ${params.freyja_varthresh}`
   (default 0.0). One process produces `*.variants.tsv` and `*.depths`; do not run a separate mpileup. Always published.
6. `FREYJA_DATA` (mode-independent): resolves `--barcodes` once per run into a single value channel
   (dir containing `usher_barcodes.feather|csv`, `lineage_mutations.json`, `curated_lineages.json`, `lineages.yml`,
   `last_barcode_update.txt`) so every sample uses the identical dataset. `latest` runs `freyja update --outdir`
   in the pinned Freyja container and publishes to `--db_dir/freyja/<barcode-version>/`. The barcode version is the
   verbatim content of `last_barcode_update.txt` (format `MM_DD_YYYY-HH-MM`, e.g. `09_17_2026-00-46`; it has a time of
   day, so same-day updates differ); directories are named by its ISO form `YYYY-MM-DDTHH-MM`. `--barcodes` accepts
   `latest`, a full stamp (either form), a bare `YYYY-MM-DD` (newest dataset of that day already in `--db_dir`; error
   with a helpful message if none), a dataset directory, or a `.csv`/`.feather` file (used verbatim, checksummed; version
   recorded as `custom`). Also supported: `--barcodes container`, which uses the snapshot bundled in the pinned Freyja
   container (`freyja demix --version` reports its version) for fully offline runs. Run-time update therefore cannot
   change mid-run, and a run never silently falls back to a different dataset.
7. `FREYJA_DEMIX` per sample, with `--barcodes <file>` (+ `--lineageyml` if present) and exposed knobs: `freyja_eps`,
   `freyja_covcut`, `freyja_depthcutoff`, `freyja_adapt`/`freyja_autoadapt`, `freyja_relaxedmrca`, `freyja_solver`,
   `freyja_confirmedonly`. **Never swallow failure**: use `errorStrategy` ignore-with-status; emit
   `<sample>.demix_status.json` `{status: ok|failed|skipped_low_coverage, reason}`. A coverage gate
   (`--min_genome_coverage_for_demix`, % of genome >= `freyja_covcut`, default 60 to match Freyja's plot `--mincov`) skips
   demix with a recorded reason (F2, F3).
8. `FREYJA_AGGREGATE_PLOT` (see decision D6 for the plot workaround): stage ok `.demix` files in a directory (trailing slash), `freyja aggregate`, then our own
   post-step to strip filename extensions from sample labels (no `_` mangling), then `freyja plot` (and `--times` plot
   when `collection_date` is present). Emit `freyja_aggregated.tsv`.
9. Coverage QC (F8, F12): `MOSDEPTH` (per-base + summary) on the trimmed BAM; `bin/amplicon_coverage.py` computes
   per-amplicon mean depth and dropout (< `amplicon_min_depth`, default 10) from the primer BED + depths; QC row per
   sample: reads (fastp), mapped reads, mean depth, SD, % genome >= 10x, % >= `consensus_min_depth`, amplicons dropped.
   Sample QC `PASS|WARN|FAIL` thresholds are params (default PASS if >= 80% genome at >=10x, mirroring the old `>79`).
10. Negative controls (F7): for `sample_type=negative`, QC from mapped reads and genome coverage
    (`ntc_max_mapped_reads`, `ntc_max_genome_coverage_pct`); FAIL flags the whole run in the report. Optional
    `--kraken2_db` (opt-in) is a later extra.
11. Consensus (D5): built from the **Freyja variants TSV** (already `samtools mpileup | ivar variants`), so no second
    variant caller and no LoFreq. Chain: `IVAR_VARIANTS_TO_VCF` (bundled viralrecon converter) -> `CONSENSUS_FILTER_VCF`
    (one `bcftools view -i` expression, no plugins) -> `LOW_DEPTH_MASK` -> `CONSENSUS_BCFTOOLS`
    (`bcftools consensus --mark-del '-' -H I -m mask`), FASTA header = sample ID. Rules (parity with the old pipeline's
    `setGT` chain, plus a read-count guard borrowed from h2seq): DP >= `consensus_min_depth_variants` (10) and ALT_QUAL >= 20
    for any call; **SNV with AF >= `consensus_freq` (0.75) -> alt allele; SNV with `consensus_iupac_min_af` (0.25) < AF <
    `consensus_freq` and ALT_DP >= `consensus_iupac_min_alt_reads` (5) -> heterozygous genotype, emitted as an IUPAC
    ambiguity code (`-H I`)** (`--consensus_iupac false` disables, leaving the reference base); SNV AF <=
    `consensus_iupac_min_af` -> ignored; indel with AF >= `indel_freq` (0.75) -> applied, indel AF < `indel_freq` -> ignored
    (never IUPAC). Positions with depth < `consensus_min_depth` (10, as in the old pipeline, from the depths file; variant sites exempt as in the
    old pipeline) are masked N. Because `freyja_varthresh` defaults to 0.0 the TSV holds every ALT; thresholds are applied
    here, not in Freyja. Status: `ok` (including a clean reference-identical sample with a header-only variants file),
    `empty` (no usable depth: no position outside called variants reaches `consensus_min_depth`, or empty depths file; output
    all-N of reference length), `failed` (e.g. REF mismatch; all-N). Never an empty file. A missing `consensus_status.json`
    is treated as `failed` by the report builder. `--skip_consensus` disables. Concatenated to
    `summary/consensus_genomes.fasta`, stats to `summary/consensus_stats.tsv`.
12. Opt-in extras: `--annotate_variants` (F16; `bcftools csq` port of `amino_acid_consequences`, run on the consensus-filter VCF),
    `--mutations_of_interest <csv>` (F10: columns `name,mutation` e.g. `S:N501Y` or nucleotide; frequency + depth per
    sample from the Freyja variants TSV, output `summary/mutations_of_interest.tsv`),
    `--freyja_boot` (F6; `freyja boot`, `--freyja_boot_n` default 100), `--compare_barcodes <date|path>` (F5; second
    demix with another dataset; report shows per-lineage abundance deltas).
13. `DUMP_SOFTWARE_VERSIONS` + `PROVENANCE` (F1): `pipeline_info/run_info.json` with pipeline version, git
    commit, nextflow version, all container tags/digests (if lockfile), Freyja version, full params, input checksums
    (samplesheet, BED, reference), mode, timestamps, and a **`freyja_data` block that is mandatory for any run that
    demixes**: `{source: latest|container|pinned|custom, barcode_version (verbatim last_barcode_update.txt),
    barcode_version_iso, barcodes_file, barcodes_sha256, barcodes_format, lineage_mutations_sha256,
    curated_lineages_sha256, lineages_yml_sha256, pathogen, freyja_version, freyja_reported_barcode_version (parsed
    from `freyja demix --version --pathogen <p>` run against the exact dataset), dataset_path, downloaded_at}`.
    The barcode version is also written into: every `<sample>.demix_status.json`, `summary/demix_status.tsv`,
    `summary/freyja_aggregated.tsv` provenance sidecar (`summary/barcode_version.json`), `report_data.json`, and the
    HTML report header. In `demix`/`aggregate`/`report` modes the new run records its own `freyja_data` block plus a
    `previous_runs` list (path, previous run_info sha256, and that run's barcode version when discoverable; `unknown`
    for legacy Snakemake results, which recorded none). `--compare_barcodes` (F5) records both datasets. A pytest
    asserts these fields are present in every output that claims them.

**`run_info.json` top-level keys (contract; the report builder reads exactly these):** `pipeline_name, pipeline_version,
git_commit, nextflow_version, mode, start, params, freyja_version, tool_versions, previous_runs, freyja_data,
input_checksums, provenance_generated`. `previous_runs` = list of `{path, run_info_sha256, barcode_version}`. The
Nextflow-supplied base JSON must use these names.

14. `REPORT` (§7).

## 6. Outputs

```
<outdir>/
  samples/<sample>/{fastp,bam,freyja,coverage,consensus,variants}/...
  summary/{freyja_aggregated.tsv,freyja_plot.pdf,qc_summary.tsv,consensus_genomes.fasta,consensus_stats.tsv,
           amplicon_coverage.tsv,mutations_of_interest.tsv,demix_status.tsv,report_data.json,report.html}
  pipeline_info/{run_info.json,software_versions.yml,execution_*.{html,txt}}
```

`summary/report_data.json` is the single machine-readable input of the report and is what `--mode report` consumes.
Legacy-layout results are normalised into the same structure when read.


## 6a. File contracts (parallel work depends on these; change only via the planner)

All TSVs are tab-separated with a header, UTF-8, `NA` for missing values, sorted by `sample` unless stated.

**`qc_summary.tsv`** (one row per sample incl. negatives), columns in this order:
`sample, sample_type, collection_date, site, raw_reads, trimmed_reads, mapped_reads, mapped_pct, mean_depth, depth_sd,
median_depth, pct_genome_ge_1x, pct_genome_ge_10x, pct_genome_ge_min_depth, amplicons_total, amplicons_dropped,
qc_status, qc_reasons, demix_status, demix_reason, consensus_length, consensus_n_pct, barcode_version`.
`qc_status` in `PASS|WARN|FAIL|NA`; `qc_reasons` is `;`-joined machine tokens (e.g. `low_coverage`, `ntc_reads`,
`amplicon_dropout`, `no_reads`) or empty. Counts are integers; percentages are 0-100 floats with 2 dp.
`pct_genome_ge_min_depth` uses `consensus_min_depth`. `demix_*`, `consensus_*`, `barcode_version` may be `NA` when
produced downstream and are filled by the report-data builder, not by the QC script.

**`amplicon_coverage.tsv`** (long): `sample, amplicon, start, end, mean_depth, min_depth, pct_ge_min_depth, dropped`
(`dropped` is `true|false`). Amplicon coordinates are derived from a primer BED: amplicon = `<scheme>_<n>` pools from
the `_LEFT/_RIGHT` or `_F/_R` primers (also handle names like `nCoV-2019_1_F`, `SARS-CoV-2_1_LEFT`, and the swift
`covid19genome_..._F/_R` style) spanning from the left primer start to the right primer end (inner coordinates are
not used).

**`demix_status.tsv`**: `sample, status, reason, genome_coverage_pct, barcode_version, barcode_version_iso` (status
`ok|failed|skipped_low_coverage|missing`; `missing` = discovery-skipped or zero-byte legacy demix, barcode columns empty).

**`freyja_aggregated.tsv`**: Freyja's own aggregate format (tab-separated; index column = sample label; columns
`summarized, lineages, abundances, resid, coverage, pathogen`; list-valued fields stored as Python-style strings). Look at the
bundled example `aggregated_result.tsv` inside the Freyja container to see it.

**`annotated_variants.tsv`** (opt-in, `--annotate_variants`): `reference, sample_id, gene, nt_pos, ref, alt, effect, aa_bcsq, aa_standard, vaf, depth, compound`
(`compound` = positions of a joint multi-SNV codon change, e.g. `23063+23064`, else `NA`).

**`mutations_of_interest.tsv`**: `sample, name, mutation, ref_depth, alt_depth, depth, alt_freq`.

**`consensus_stats.tsv`**: `sample, length, n_count, n_pct, ambiguous_count, status` (status `ok|empty|failed`).


**Contract clarifications (resolved by the planner after the QC-scripts task):**
- `amplicon_coverage.tsv` `start`/`end` are **1-based inclusive** (BED start + 1, BED end), matching the depths file.
  `dropped` = `mean_depth < min_depth` over the amplicon span (primers included); positions absent from the depths file
  count as depth 0. Amplicon counts: articv3 98, articv4.1 99, articv5.3.2 96, midnight 29, midnightv2 29, eden 14,
  qiagenv1 223, swift 343, xgen 343.
- `qc_summary.tsv`: read counts are reads as fastp reports (a pair = 2); `mapped_reads` = flagstat `primary mapped`
  (fallback `mapped`); `mapped_pct` denominator = trimmed reads, else raw, else flagstat primary total, else `NA`.
  Extra `qc_reasons` token `no_depths` (unreadable depths -> FAIL for samples; `NA` for negatives with no data).
  `--concat` exists for `qc_summary.py` and `consensus_stats.py` only; the pipeline joins `amplicon_coverage` and
  `mutations_of_interest` per-sample files by keeping a single header.
- `mutations_of_interest.tsv`: rows in input order; mutation syntax: SNV `A23063T`; deletion by first deleted base
  `C75-` (converted to ivar's anchor convention); insertion by anchor `A14+C`; amino acid `S:N501Y`. Absent mutation =>
  `alt_depth 0`, `alt_freq 0`, `ref_depth = depth` (from `--depths`; amino-acid = min over the 3 codon positions).
  Amino-acid fallback uses the GFF3 (needs `--gff` and `--reference`); ORF1a vs ORF1ab named from the mRNA `Name`.
- The real ivar TSV has a 20th column `POS_AA`.


**Contract changes from the integration round:**
- `qc_max_dropped_amplicons` default is **1** (was 10): >1 dropped amplicons -> WARN `amplicon_dropout`.
- `tests/data/truth.json` gained `expected_qc_status`, `expected_in_aggregate`, `expected_consensus_status`; new sample
  `NEG_CONTAM` (contaminated negative that passes the coverage gate).
- `demix_status.tsv` may contain `missing` rows (discovery-skipped or zero-byte legacy demix files) with an empty barcode version.
- `consensus_stats.tsv` may contain `failed` rows with `NA` values for samples whose consensus steps failed (synthesised by the
  workflow; consensus/QC steps use `errorStrategy 'ignore'`).
- `summary/samplesheet.csv` is published in every mode that has an input sheet. Aggregate mode writes no `barcode_version.json`
  (it demixes nothing). BUILD_REPORT_DATA/RENDER_REPORT do not emit `topic: versions` (it deadlocked PROVENANCE), so their
  python/jinja2 versions are only in `run_info.containers`.
- **Capture protocol (D13) contract:** `run_info.json` gains top-level `protocol`, `primer_trimming`, `duplicates` (effective
  values; the last two `null` in modes without reads). Capture runs have no `amplicon_coverage.tsv`, `amplicons_total` and
  `amplicons_dropped` are `NA` and `amplicon_dropout` is never raised. New per-sample files `samples/<s>/bam/<s>.markdup_stats.txt`
  and `<s>.{mapped|primertrimmed}.markdup.bam(.bai)` (duplicates are flagged, never removed; the unmarked BAM is not published when
  marking is on); `<s>.mapped.bam`/`<s>.mapped.flagstat` appear when nothing was trimmed. Report header gains `protocol`,
  `primer_trimming`, `duplicates`. `tests/data_capture/truth.json` adds `protocol`, `n_duplicate_pairs`,
  `expected_duplicate_fraction_range`. Verified: duplicate marking lowers Freyja depths (CAP_MIX mean depth 54.67 -> 45.27).
- `make_samplesheet.py --from-results --skipped-out` writes `sample, reason, sample_type, collection_date, site`; the helper
  scripts `make_samplesheet.py`/`check_samplesheet.py` run on Python >= 3.9 (3.8 checked statically).
- `primer_scheme` defaults to null (meaning midnight in `full`; in `demix` amplicon coverage runs only when a scheme/BED is
  given explicitly); `run_info` records the effective scheme. `run_info.params` omits `pipelines_testdata_base_path`.
- `process.resourceLimits` is computed from the host (cpus, 90% of memory, 240 h) because nf-core lint forbids `max_*` params.

**`report_data.json`**: documented in `docs/report_data.md` (owned by the report task); assembled by
`bin/build_report_data.py` from the files above plus `run_info.json`.

## 7. HTML report

Self-contained, offline (inline CSS/SVG/vanilla JS, no CDN), Jinja2 template `assets/templates/report.html.j2`, rendered
by `bin/render_report.py` from `report_data.json`; degrades gracefully when inputs are missing (pattern from
quasispecies pipeline `render_sample_report.py`; sentinel statuses become badges). Sections: run header/provenance
(versions, barcode date, key params, git commit); run-level flags (NTC fail, samples failed/skipped); QC table with
PASS/WARN/FAIL; lineage abundance stacked bars + sample x lineage heatmap (summarized lineage groups from Freyja);
per-sample coverage track and amplicon dropout view; demix status with reasons; mutations of interest (if given);
time-series panel when dates exist; barcode comparison (if run); consensus table (N content, length). Must remain
readable at 100+ samples (sortable/filterable table, no per-sample page required).

## 8. Testing (F14)

- `nf-test`: per local module with `-stub` and real tiny-data tests; pipeline test for each mode.
- `pytest tests/python`: helper scripts, report renderer, samplesheet logic, `bump_freyja.py`.
- **Synthetic dataset** (`tests/data/`, generator `tests/make_test_data.py`, deterministic seed): derive two or three
  "lineage" genomes by applying SNVs from a small, committed barcode subset (`tests/data/barcodes_test.csv`, real
  lineage names present in Freyja's lineage hierarchy) onto NC_045512.2; simulate paired-end amplicon reads from the
  chosen primer BED at known mixture ratios (pure-Python simulator, no extra deps); include one negative control and
  one low-coverage sample. Keep total size small (< 5 MB) for CI.
- Output checks use tolerances (abundance within +-0.05 of truth, consensus differences, QC flags), pattern from
  quasispecies `check_outputs.py`. Tests pass `--barcodes tests/data/barcodes_test.csv` so they are deterministic and
  offline; one separate non-blocking job exercises `--barcodes latest`.
- Demix-only equivalence test: run `full`, then `demix --previous_results` with the same barcodes; demix outputs must be
  byte-identical. Also a legacy-layout fixture for `make_samplesheet.py --from-results`.


### Edge cases found with real data (synthetic set, Freyja 2.0.5)

- A sample with ~0 mapped reads (e.g. `NEG_E`, 3 read pairs) makes `freyja variants` write a header-only TSV, and
  `freyja demix` then **crashes** (`TypeError: operation 'radd' not supported for dtype 'str'`). The coverage gate
  (depths file) must route such samples to `skipped_low_coverage`; independently, any `freyja demix` failure is recorded
  as `failed` + reason (never fatal). Also a sample with 0 mapped reads must not abort `IVAR_TRIM`/`FREYJA_VARIANTS`:
  guard on an empty BAM (flagstat) and emit empty-but-valid outputs plus a status flag.
- Freyja demixes `LOWCOV_D` (0.6x mean depth) to a confident single lineage with `coverage` 0.0; **only our gate stops that
  from being reported.** Keep the gate on by default.
- Small spurious minor lineages (BA.2 ~0.03 in the BA.1/BA.5 mixtures) are normal; test tolerance: true lineages within
  +-0.05, and any lineage absent from truth below 0.05.
- Freyja's `summarized` row collapses e.g. JN.1 into `BA.2.86*`; tests compare the `lineages`/`abundances` rows.
- The synthetic set (`tests/data`, owned by `tests/make_test_data.py`) is the only test dataset; `truth.json` holds
  expected demix status, QC reasons and %genome>=10x ranges. Samples: MIX_A, MIX_B, PURE_C, LOWCOV_D, NEG_E, Site_1_week_2.

## 9. Freyja auto-update CI (D1)

- Freyja version is pinned in exactly the local `freyja_*` modules (container + `environment.yml`); an invariant test
  fails if they disagree.
- `bin/bump_freyja.py`: `--check` prints pinned vs latest (queries the quay.io biocontainers tag API, ignores
  pre-releases, verifies the bioconda package exists), `--apply` rewrites all pins deterministically and appends a
  CHANGELOG line. It never edits anything silently outside this explicit command.
- `.github/workflows/freyja-update.yml`: scheduled (daily) + manual; runs `--check`; if newer, creates branch
  `bot/freyja-<ver>` and a PR via a GitHub App token or PAT (PRs from the default `GITHUB_TOKEN` do not trigger CI).
  The normal CI (lint + nf-test docker + pytest + all-mode pipeline tests on the synthetic data) must pass.
  Auto-merge (`gh pr merge --auto`) is enabled only for **same-major** bumps; a new major version leaves the PR open for
  human review. Merging requires branch protection with required checks (a repo setting the user must enable).
- Barcode data is not CI-managed (it is data, resolved per run via `--barcodes`).

## 10. Repo, conventions, branding

- Generated with `nf-core pipelines create` (same style as h2seq: nf-schema, utils subworkflows, `conf/base.config`,
  `resourceLimits`, `.nf-core.yml` lint exceptions). Pipeline name `wasteline`; `manifest.author` Charles Foster.
- Container pins: exact tags now; `conf/containers.lock.config` (digests) generated by `bin/lock_containers.py`
  (F13, optional use via `-profile locked`).
- **Logo** (printed by `nextflow run .` and `--help`): ASCII art of "WASTEline" with toilet paper. Implemented by the
  planner in the nf-core logo helper (colour off with `--monochrome_logs`/`NO_COLOR`).

## 11. Phases and gates

| Phase | Content | Owner |
|---|---|---|
| 0 | Scaffold (nf-core create, CLAUDE.md, agents, pixi, logo, bundled references/primers) | planner |
| 1 | Core `full` path to aggregate/plot, `FREYJA_DATA`, params/schema, stub tests, versions/provenance | coder, planner review |
| 2 | Modes `demix`/`aggregate`/`report`/`update_db`, helper scripts, legacy-layout support | coder |
| 3 | Coverage/amplicon/NTC QC, consensus, status files | coder |
| 4 | HTML report | coder (template), planner reviews design |
| 5 | Extras F5, F6, F10, F11, F13, F16 | coder |
| 6 | Synthetic data, full tests, CI, Freyja auto-update, docs (usage, output, parameters) | coder |

Gate between phases: `nf-core pipelines lint` clean (agreed exceptions only), nf-test green (stub at minimum),
planner review of the diff. Anything needing Docker is verified with the Docker daemon running.

Future (not now): test `minibwa` as an alternative aligner against `bwa mem` on identical inputs; additional pathogens via
`conf/pathogens/<name>.config` + Freyja `--pathogen`; CRAM output; resolving historical barcode versions from the
Freyja git history; optional kraken2 negative-control screening.

## 12. Verified facts about Freyja (main branch source, 2.x) and the old pipeline

- `freyja variants` = `samtools mpileup -aa -A -d 600000 -Q minq -q 0 -B | tee depths | ivar variants -t varthresh`;
  default `--varthresh 0.0`. Old pipeline passed 0.25.
- `freyja demix --barcodes` accepts `.csv` or `.feather`; default SARS-CoV-2 data is `usher_barcodes.feather`.
  Options include `--eps --covcut --depthcutoff --lineageyml --adapt --autoadapt --relaxedmrca --solver --pathogen`.
- `freyja update --outdir` downloads barcodes, `last_barcode_update.txt`, `lineage_mutations.json`,
  `curated_lineages.json`, `lineages.yml`, `pathogen_config.yml`.
- `freyja aggregate <dir/> --ext X --output f`; verified: the sample label is the *variants file name recorded in the demix file header* (e.g. `S1.variants.tsv`), not the demix file name; `clean_aggregate_labels.py` strips `.variants.tsv`/`.tsv` and never touches underscores. `freyja plot` drops samples with
  `coverage <= --mincov` (default 60) and supports `--times`.
- `freyja ampliconstat` exists but needs primer BED naming conventions our bundled BEDs may not follow; we compute
  amplicon coverage ourselves.
- Verified: `freyja variants --variants S1.variants` writes `S1.variants.tsv` (ivar appends `.tsv`); depths file is 4 columns, no header, one row per position. Still glob for both names.
- Verified: the Freyja 2.0.5 container bundles a dataset (`usher_barcodes.feather|csv`, `last_barcode_update.txt` = `09_17_2026-00-46`, ...) under `/usr/local/lib/python3.12/site-packages/freyja/data`; `freyja demix --version --pathogen SARS-CoV-2` reports the same version.
- Verified: containers run only as `linux/amd64`; on Apple Silicon Docker needs `--platform linux/amd64`.
- Old-pipeline defects we deliberately do not reproduce: empty coverage PDFs (`touch`), demix errors swallowed with
  `|| touch`, `_`->`-` sample-ID rewriting, mutable shared dataset updated at run start, `sed` edits of the Snakefile to
  bump versions, VCF fed to `freyja demix` for lofreq/ivar callers.
