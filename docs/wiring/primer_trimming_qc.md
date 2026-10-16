# Primer-trimming effectiveness QC: wiring notes

New per-sample QC column `primer_trimmed_pct` and WARN reason `low_primer_trimming`: the share of mapped reads from which
`ivar trim` removed primer sequence. Aim: catch a wrong primer scheme. **Weak by nature** (see "Calibration"); the rule
ships **off**. Nothing outside the files listed below was edited, so the pipeline is **not yet wired** and the workflow's
current `QC_SAMPLE(...)` call (5-element tuple) no longer matches the module (7 elements).

**Follow-on: the statistically calibrated scheme-match check (section 5) is the real wrong-scheme detector**; the ivar
percentage above is informational (rule off). Both add columns to `qc_summary.tsv`, now **27 columns**.

Files: `bin/primer_scheme_check.py`, `modules/local/read_starts/`, `modules/local/primer_scheme_check/` (new; section 5),
`bin/ivar_trim_stats.py` (new), `bin/qc_summary.py`, `modules/local/qc_sample/` (main.nf, tests, snapshot),
`modules/local/qc_concat/tests/` (re-pointed at new 24-column fixtures), `tests/python/test_ivar_trim_stats.py`,
`tests/python/test_qc_summary.py`, `tests/fixtures/qc/` (real ivar logs, `concat/`), `docs/qc_scripts.md`.

## 0. Blocking finding: the published ivar log is empty

`ivar trim` prints its summary to **stderr**. `modules/nf-core/ivar/trim/main.nf` runs `ivar trim ... > ${prefix}.ivar.log`
(stdout only), so every published `samples/<s>/bam/<s>.ivar_trimmed.ivar.log` is **0 bytes** (verified on the complete
e2e output). The stats script reports `status=empty`/`NA` for such a log. Fix one of:

1. (preferred) patch the nf-core module (`nf-core modules patch ivar/trim`): `> ${prefix}.ivar.log 2>&1`.
2. (no module edit) `ext.args` in `conf/modules.config` for `IVAR_TRIM`, appending `2> ${meta.id}.ivar_trimmed.ivar.log`
   (verified in a container: `ivar trim -x 0 -e 2> X.log -i in.bam -b b.bed -p pre > X.log` leaves the full summary in `X.log`):

```groovy
withName: 'IVAR_TRIM' {
    ext.args   = { "-x ${...} -e 2> ${meta.id}.ivar_trimmed.ivar.log" }   // stderr (ivar's summary) into the log
```

(The log would then also hold ivar's BED warnings and the per-primer read table; the parser ignores them.)

## 1. Parameter

`qc_min_primer_trimmed_pct`: number, 0-100, **default 0** (0 = rule disabled; the column is still reported). Description:
"WARN (low_primer_trimming) when a sample's % of reads that ivar primer-trimmed is below this. 0 disables. Only meaningful
for protocols whose reads start at primer sites; calibrate on your own good runs." Pass it as
`--min-primer-trimmed-pct ${params.qc_min_primer_trimmed_pct}` in the existing `QC_SAMPLE` `ext.args` list in
`conf/modules.config`. Add to `nextflow.config`, `nextflow_schema.json` (minimum 0, maximum 100) and
`docs/parameters.md`.

## 2. Channel into `QC_SAMPLE` (workflows/wasteline.nf)

`QC_SAMPLE` input is now `tuple val(meta), path(depths), path(flagstat), path(fastp_json), path(amplicons), path(ivar_log), path(primer_scheme_summary)`
(7 elements; the last two are `[]` when absent; section 5 wires the last one). In mode `full` with trimming run, join `IVAR_TRIM.out.log` (a `[meta, log]` channel); in
capture protocol, `skip_primer_trimming` and demix mode there is no log, so pass `[]` (`primer_trimmed_pct` = `NA`, rule
inert):

```groovy
def ch_ivar_log = channel.empty()
if (<full mode> && !capture && !params.skip_primer_trimming) { ch_ivar_log = IVAR_TRIM.out.log }
def ch_qc_in = ch_depths
    .join(ch_flagstat, remainder: true)
    .join(ch_fastp, remainder: true)
    .join(ch_amplicons, remainder: true)
    .join(ch_ivar_log, remainder: true)
    .join(ch_scheme_summary, remainder: true)   // section 5
    .map { meta, depths, flagstat, fastp_json, amplicons, ivar_log, scheme_summary ->
        [ meta, depths, flagstat ?: [], fastp_json ?: [], amplicons ?: [], ivar_log ?: [], scheme_summary ?: [] ] }
```

Inside the process, `ivar_trim_stats.py` runs first (`|| true`, so a log problem can never cost the sample its QC row) and
its TSV is passed via `--ivar-trim-stats`. The stats TSV is not an output of the module.

## 3. Contracts and downstream

- `qc_summary.tsv` has **27 columns**: the 23 of PLAN 6a unchanged, then appended in this order `primer_trimmed_pct` (0-100,
  2 dp), `primer_scheme_check` (`ok|mismatch|uninformative|NA`), `primer_scheme_best` (name(s) joined by `/`, or `NA`),
  `primer_scheme_enrichment` (declared scheme, 1 dp); `NA` when unavailable. Update PLAN 6a ("Contract changes"), `docs/output.md` column dictionary and the
  `tests/default.nf.test.snap` snapshot of the pipeline test (qc_summary md5 changes).
- `QC_CONCAT` (`qc_summary.py --concat`) rejects (exit 2) any per-sample table with a different header, so a mixed
  old/new set fails loudly. `bin/build_report_data.py` reads columns by name and copies the header in `--write-qc-tsv`,
  so it tolerates the new columns; to show them, add `primer_trimmed_pct` (and the three `primer_scheme_*`) to `QC_COLUMNS`/`FLOAT_COLUMNS` (report data
  `samples[].primer_trimmed_pct`, float or null) and to the template.
- Reason tokens `low_primer_trimming` and `primer_scheme_mismatch`, appended in this order after `low_coverage`/
  `amplicon_dropout`; status WARN only, never FAIL; negatives: none. Add both to the report's reason labels/legend and
  `primer_scheme_check/best/enrichment` to `samples[]` (strings/float or null). Report wording for `uninformative`:
  "reads are not primer-anchored (e.g. tagmentation), scheme could not be checked" - not a warning.
- Docs to update: `docs/output.md` (column table + reason list), `docs/parameters.md`, PLAN 6a.

## 4. Calibration (ivar 1.4.4, synthetic tests/data, midnight-shaped reads)

Real pipeline runs (`-profile docker`, samples MIX_A, PURE_C, Site_1_week_2; `--primer_scheme` midnight, articv4.1, swift,
midnightv2; logs taken from each task's `.command.err` because the published logs are empty), plus `ivar trim` run
directly on the same mapped BAMs for the other bundled schemes. Primer-trimmed % of mapped reads:

| scheme                   | MIX_A | PURE_C | Site_1_week_2 |
| ------------------------ | ----- | ------ | ------------- |
| midnight (correct)       | 22.33 | 21.97  | 21.70         |
| midnightv2 (same coords) | 22.33 | 21.97  | 21.70         |
| articv3                  | 10.96 | 10.68  | 9.51          |
| articv4.1                | 8.74  | 8.71   | 7.83          |
| articv5.3.2              | 9.38  | 8.88   | 8.69          |
| eden                     | 1.59  | 1.49   | 1.29          |
| qiagenv1                 | 20.46 | 21.73  | 21.33         |
| swift (-x 5)             | 40.26 | 39.93  | 40.70         |
| xgen                     | 52.94 | 52.71  | 50.55         |

The numbers do **not** separate right from wrong: dense schemes (swift, xgen, qiagenv1) trim as much or more than the
correct one on midnight-shaped data. Only sparse wrong schemes (eden, artic*) give lower values. Conclusion: default 0
(disabled). Expectations that could not be measured here: QIAseq (qiagenv1, reads begin at primer sites, planned default
scheme) should trim a large fraction with the right scheme, so a threshold becomes meaningful once real QIAseq runs are
seen; tagmentation-prepared long-amplicon libraries (Midnight) would trim only a few % even when correct. Suggested
practice: look at `primer_trimmed_pct` over known-good runs, then set the param to about half the lowest good value.

## 5. Scheme-match check (READ_STARTS + PRIMER_SCHEME_CHECK)

Method, statistics, limits and measurements: [qc_scripts.md](../qc_scripts.md#primer_scheme_checkpy). Summary: counts
mapped reads whose 5' end is within +-2 bp of a same-strand primer 5' end, for all nine bundled schemes plus a custom
BED, as an enrichment over chance with a binomial p-value. Declared midnight on midnight-shaped reads: enrichment ~41-43;
the other eight schemes 1.1-3.3; random fragments ~1. Declared `articv4.1` on midnight reads: `mismatch`, best
`midnight/midnightv2`.

Modules (both stub + real docker tested):

- `READ_STARTS`: `tuple val(meta), path(bam)` -> `tuple val(meta), path('<id>.read_starts.tsv')` (`emit: tsv`). Input is the
  mapped BAM **before** `IVAR_TRIM` (the `ch_mapped.mapped` BAM; no index needed). samtools container as
  `SAMTOOLS_MARKDUP` (`htslib_samtools:1.24`).
- `PRIMER_SCHEME_CHECK`: `tuple val(meta), path(read_starts)`, `val declared`, `path beds_dir`, `path custom_bed` ->
  `emit: tsv` (per scheme), `emit: summary` (one row). Python container as `QC_CONCAT`.

Wiring (workflows/wasteline.nf), mode `full` with trimming and protocol `amplicon` only; capture, demix and
`skip_primer_trimming` have no primer-anchored reads or no pre-trim BAM, so pass no summary (columns `NA`):

```groovy
READ_STARTS(ch_mapped.mapped.map { meta, bam, _bai, _count -> [ meta, bam ] })          // before IVAR_TRIM
PRIMER_SCHEME_CHECK(
    READ_STARTS.out.tsv,
    params.primer_bed ? 'custom' : (params.primer_scheme ?: 'midnight'),               // declared name (default scheme per D14 later)
    file("${projectDir}/assets/primer_schemes"),
    params.primer_bed ? file(params.primer_bed) : [])                                  // custom BED = declared scheme
ch_scheme_summary = PRIMER_SCHEME_CHECK.out.summary.map { meta, tsv -> [ meta, tsv ] }
```

Notes for the wiring:

- Negatives: run the check on samples only (`meta.sample_type == 'sample'`) or on all; the QC rule only acts on
  `sample_type=sample` and a negative has too few reads (`NA`, harmless). Skipping negatives saves little.
- `conf/modules.config`: `PRIMER_SCHEME_CHECK` `ext.args` optional (`--window`, `--min-enrichment`, `--min-reads`,
  `--genome-length`; defaults 2, 5, 500, 29903); `errorStrategy = 'ignore'` like `QC_SAMPLE` (a missing summary is just
  `NA` columns). No pipeline param is needed (the rule is always on and WARN only); a param to disable the check is
  optional. The outputs need not be published (optionally `samples/<id>/qc/<id>.primer_scheme_check.tsv`).
- Pass `--genome-length` from the reference when it is not SARS-CoV-2 (29903 default).
- When the declared scheme is the pipeline default (`midnight` now, `qiagenv1` planned), `params.primer_scheme` null must
  resolve to that default so the declared name exists in `assets/primer_schemes`.
- `docs/output.md`: document the three columns and the `primer_scheme_mismatch` token; mention `uninformative`.
- Tests touching `qc_summary.tsv` (`tests/default.nf.test.snap`, `tests/check_outputs.py` column checks) must expect 27
  columns.
- QIAseq: reads should be strongly anchored, so `ok` with a high enrichment is expected; the real number is unknown until a
  QIAseq run is checked (dense scheme: the window covers ~4.5 % of the genome, so the maximum enrichment is ~22). That run
  is the best calibration for `--min-enrichment`.
