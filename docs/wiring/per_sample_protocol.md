# Per-sample protocol and primer scheme (PLAN D20): wiring as built

One run can mix amplicon and probe-capture libraries, or several primer schemes, through two optional samplesheet columns
(`protocol`, `primer_scheme`; modes `full` and `demix`). Blank or absent cells mean the run defaults, so a samplesheet without
the columns behaves exactly as before D20 (the default e2e steps prove it).

## Where the settings live

`libraryMeta(sample, protocol_cell, scheme_cell, mode)` in `subworkflows/local/utils_nfcore_wasteline_pipeline/main.nf` is the one
place that computes them. `sampleMeta(meta, mode)` adds the result to the sample meta in `full` and `demix` mode (the three
call sites: the samplesheet rows in `PIPELINE_INITIALISATION`, `resultsRow()` and the `--demix_dir` rows of `RESULTS_INPUT`).

| `meta` key        | Value                                                                                                                                                                                                                        |
| ----------------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `protocol`        | the cell, else `--protocol`                                                                                                                                                                                                  |
| `primer_trimming` | `protocol == amplicon && !--skip_primer_trimming`                                                                                                                                                                            |
| `duplicates`      | `--duplicates`, with `auto` = `mark` for capture and `none` for amplicon                                                                                                                                                     |
| `primer_scheme`   | only when `primer_trimming`: a bundled name, or `custom` (then `primer_bed` is the path). The cell, else `--primer_bed`, else `--primer_scheme`, else (full mode only) `qiagenv1` (D14). `null` = no amplicon-level analysis |
| `primer_bed`      | the BED path for `custom`, else `null`                                                                                                                                                                                       |

`protocolSettings(protocol)` holds the part that depends only on the protocol and is also used for the run defaults
(`validateInputParameters`, `runInfoBase`). `primerBedFor(meta)` resolves a sample's BED. The old `protocolSettings(mode)`
(global, with an `amplicons` flag), `primerBed()` as the BED of every sample and the unused `ivarOffset()` are gone; a
sample has amplicon analysis exactly when `meta.primer_scheme` is set.

A cell value is a BED path when it contains `/` or ends in `.bed` (relative to the launch directory, like the FASTQ paths),
otherwise it must be a scheme name: bundled for the run's pathogen or a `<name>.bed` of `--primer_scheme_dir` (PLAN D26, see [multi_pathogen.md](multi_pathogen.md): the bundled schemes are per pathogen and a user-directory scheme is stored as `custom` with its path). Errors name the sample. The samplesheet rows are built eagerly at start-up, so a bad
cell stops the run before any work; `checkSampleBeds()` checks that each distinct BED is on a contig of `--fasta` (the run
default BED is still checked in `validateInputParameters`).

## Workflow (`workflows/wasteline.nf`)

Nothing is switched on or off for the whole run any more; every step takes its samples from the shared channels.

- **Trimming**: `ch_mapped.mapped` is filtered on `meta.primer_trimming`. `IVAR_TRIM` is an nf-core module whose BED is one
  `path bed` input, which must not be edited, so the BED travels beside the BAM in a `multiMap` (both queue channels are emitted in
  lockstep, so Nextflow pairs the i-th BAM with the i-th BED). The module is unchanged and no local variant exists. `ext.args`
  (the `-x` offset: `meta.primer_scheme == 'swift'` gives 5) is a closure over `meta`.
- **Primer scheme check**: `READ_STARTS` runs for the trimmed samples; `PRIMER_SCHEME_CHECK` gets `declared` (`meta.primer_scheme`,
  `custom` for a BED) and `custom_bed` per sample through a `multiMap` the same way. The module is unchanged.
- **Sorting**: `SAMTOOLS_SORT` takes the trimmed BAMs (from `IVAR_TRIM`) plus the untrimmed mapped and empty BAMs of the other
  samples; empty BAMs of trimmed samples bypass it as before.
- **Duplicates**: `ch_final_bam` is branched on `meta.duplicates == 'mark'`; `SAMTOOLS_MARKDUP` gets the marked ones.
- **Amplicon coverage**: `AMPLICON_COVERAGE` gets the samples with `meta.primer_scheme` and their BED (`multiMap`). The QC dropout
  rule needs the amplicon table, so it is never applied to a sample without one (`amplicons_*` are `NA`). `QC_CONCAT` writes
  `amplicon_coverage.tsv` only when some table exists. `DROPOUT_IMPACT` is called whenever `--dropout_impact` is on and reads
  that table (capture samples are simply not in it); in an all-capture run it never gets an input, so it never runs.
- **Config** (`conf/modules.config`): the prefixes of `SAMTOOLS_SORT/MARKDUP/FLAGSTAT/CRAM` and the publish condition of the
  unmarked BAM read `meta.primer_trimming` and `meta.duplicates` (a `meta` without them, as in a module test, counts as
  trimmed and not marked).
- **Provenance**: `runInfoBase(mode, metas)` is built from the metas of the input rows (`ch_rows.toList()`), so `libraryInfo()` can
  record `protocol` (common value or `mixed`), `primer_trimming`/`duplicates` (common value; `null`/`mixed` for a mixed run) and
  `sample_libraries`. `PROVENANCE` gets every distinct BED (`stageAs: 'primer_beds/?/*'`); `make_run_info.py` writes
  `input_checksums.primer_bed` for one BED (unchanged) and `primer_beds` for several.
- **Report**: `build_report_data.py` adds `header.protocol_counts` (`5 capture, 2 amplicon`) for `protocol: mixed`; the template shows it
  in the "Library protocol" fact and adds a note to the amplicon section. The amplicon section itself only holds the amplicon
  samples because `amplicon_coverage.tsv` only does.

## Demix mode

`make_samplesheet.py --from-results` carries `protocol` and `primer_scheme` over from the previous `summary/samplesheet.csv` (the
columns are written only when some sample has a value). `summary/samplesheet.csv` is a plain copy of the input sheet, so a global
`--protocol capture` of the earlier run is not in it: demix applies the current `--protocol` default to blank cells (state it again;
`bin/reanalyse_runs.py` carries the parameter). In demix mode `meta.primer_scheme` stays `null` unless a cell or
`--primer_scheme`/`--primer_bed` names a scheme (the scheme of old results is unknown).

## Validation

- `bin/check_samplesheet.py` (modes full and demix): unknown `protocol`; `primer_scheme` that is neither a bundled scheme nor a `.bed`
  path; with `--check-files` also a missing BED file; a warning for a `primer_scheme` on a `capture` row.
- `assets/schema_input.json` and `assets/schema_input_demix.json`: `protocol` enum (blank allowed) and `primer_scheme` string, both
  mapped into `meta`.
- `bin/make_samplesheet.py`: `--metadata` tables and discovered results sheets pass the two columns through (and check the protocol).

## Limits

- `--skip_primer_trimming`, `--duplicates` and `--ivar_offset` stay global. There is no per-sample `duplicates` column.
- Several scheme-specific `--primer_scheme_check` thresholds are global too.
- A custom BED in a demix run is not contig-checked against the reference (only samplesheet rows of `full` and `demix` `--input`
  runs are).
- `DROPOUT_IMPACT` and the mutation-table warning now also appear in an all-capture run when the Freyja dataset has no mutation lists.
- `docs/example_report.html` was not regenerated.

## Tests

- pytest: `test_check_samplesheet.py` (validation), `test_make_samplesheet.py` (carry-over), `test_make_run_info.py` (several BEDs).
- e2e step 6 (`tests/e2e.sh`): five capture samples (run default `--protocol capture`) and MIX_A, PURE_C (`amplicon`, `midnight` in
  their cells) in one samplesheet; `tests/check_outputs.py --protocol mixed --amplicon-samples MIX_A,PURE_C` checks `run_info.json`
  (`mixed`, `sample_libraries`), per-sample QC status and amplicon columns, the amplicon table, the BAMs of each protocol, that
  `IVAR_TRIM` ran for the amplicon samples only, and the report header.
- Existing stub nf-tests and snapshots are unchanged and pass (no sample meta key reaches a snapshot).
