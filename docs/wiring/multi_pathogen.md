# Other pathogens (PLAN D22): wiring as built

`--pathogen X` runs the pipeline for the single-sequence Freyja pathogens; SARS-CoV-2 behaviour is unchanged (same params
recorded, same commands, same outputs). User view: [../other_pathogens.md](../other_pathogens.md).

## Profile selection

`nextflow.config` ends with a dynamic include chosen by `params.pathogen` (a command-line or `-params-file` value is visible
while the config is parsed; verified with `--pathogen H1N1` and `--pathogen MEASLESgenome`):

```groovy
includeConfig "conf/pathogens/${['SARS-CoV-2': 'sars2', 'MEASLESgenome': 'measles', 'RSVa': 'rsva', 'RSVb': 'rsvb', 'H1N1': 'h1n1', 'H3N2': 'h3n2', 'FLU-B-VIC': 'flu_b_vic', 'H5Nx': 'h5nx', 'H5Nx-cattle': 'h5nx_cattle'].get(params.pathogen ?: 'SARS-CoV-2', 'generic')}.config"
```

| File                                                           | Sets                                                                                                                         |
| -------------------------------------------------------------- | ---------------------------------------------------------------------------------------------------------------------------- |
| `conf/pathogens/sars2.config`                                  | `pathogen`, `fasta`, `gff` (unchanged) plus the two facts below                                                              |
| `conf/pathogens/measles.config`                                | `pathogen = 'MEASLESgenome'`, `fasta = assets/references/NC_001498.1.fasta`, `gff = null`                                    |
| `conf/pathogens/rsva.config`, `rsvb.config`                    | `pathogen = 'RSVa'` / `'RSVb'`, `fasta` / `gff` = the bundled `PP109421.1` / `OP975389.1` files (PLAN D23)                   |
| `conf/pathogens/{h1n1,h3n2,flu_b_vic,h5nx,h5nx_cattle}.config` | `pathogen`, `fasta` / `gff` = the bundled HA record (`CY121680.1`, `CY163680.1`, `KX058884.1`, `AF144305.1` twice; PLAN D24) |
| `conf/pathogens/generic.config`                                | `pathogen = null` (the strict parser needs the key declared; the command line fills it), `fasta = null`, `gff = null`        |

Every profile also sets two hidden facts (schema entries `hidden: true`, listed in `docs/parameters.md`, filtered out of
`run_info.json` so a SARS-CoV-2 run records exactly the parameters it did before):

- `lineage_groups_bundled`: path of the bundled lineage groups (`assets/lineage_groups.csv` for SARS-CoV-2), else `null`.
  `lineageGroupsFile()` returns `--lineage_groups`, else this, else `null`; `null` switches the lineage summary off with one
  `INFO` line.
- `freyja_container_pathogen`: the pathogen whose data the pinned Freyja container bundles (`SARS-CoV-2`), else `null`. Used by
  `--barcodes container` (refused for any other pathogen) and by `datasetRoot()`.

Values live in the profile files because the strict config parser does not let an included file override a value that
`nextflow.config` already sets, and because the pipeline code may not name a pathogen. `tests/python/test_boolean_params.py`
checks that the map, the files and the schema enum agree.

## Start-up checks (`validateInputParameters`, `subworkflows/local/utils_nfcore_wasteline_pipeline`)

All fail before any work. Modes `full`, `demix`, `aggregate`, `update_db`:

1. `unsupportedPathogenReason()`: `MTB` stops with the reason and what would be needed (scaling to 4.4 Mb). Influenza was refused here
   until PLAN D24 (see the end of this page).
2. Modes `full`/`demix`: `--fasta` is required when the profile has none (message: must be the reference of the Freyja barcodes);
   a FASTA with other than one sequence is refused (any pathogen, a safety net for segmented genomes).
3. `--barcodes container` / `--compare_barcodes container` need `freyja_container_pathogen` to equal `--pathogen`.
4. The default primer scheme (nothing given) whose contig is not in the FASTA gives a dedicated message (the bundled schemes belong
   to one pathogen: pass `--primer_bed` or `--protocol capture`); an explicit BED on the wrong contig keeps the old message.
5. `--annotate_variants` needs `--gff` when the profile has none.

## Pipeline changes

- `ch_fasta` / `ch_gff` are `[]` when the parameter is null (aggregate mode never needs the FASTA; `--gff` is optional:
  `UNEXPLAINED_MUTATIONS` and `MUTATIONS_OF_INTEREST` already accepted `[]`).
- **Datasets**: `datasetRoot(db_dir)` is `<db_dir>/freyja` for the container's pathogen and `<db_dir>/freyja_<pathogen>` otherwise
  (`resolveBarcodes`, `availableDatasets`, `FREYJA_DATA`, and the `FREYJA_UPDATE` publishDir in `conf/modules.config`), so a date
  can never resolve to another pathogen's data. `freyja_dataset_info.py` recognises what `freyja update --pathogen X` writes:
  `<X>_barcodes.csv|feather`, `<X>_lineages.yml` (ignored when it is not a YAML list: Freyja saves a GitHub HTML page for MPX),
  `last_barcode_update_other.txt` (`X:YYYY-MM-DD`; the version is the date, ISO `...T00-00`). `lineage_groups.py`,
  `dropout_impact.py` and `unexplained_mutations.py` find `<X>_lineages.yml` / `<X>_barcodes.csv` the same way (no
  `lineage_mutations.json`: they use the barcodes CSV, as they already could). `make_samplesheet.py` and `reanalyse_runs.py` treat
  `last_barcode_update_other.txt` as a dataset marker; `reanalyse_runs.py` carries `pathogen` over.
- **Lineage summary**: off without bundled groups or `--lineage_groups`. For a dataset downloaded during the run (`latest`) of a
  pathogen that may lack a hierarchy (`hierarchyMayBeMissing()`), `LINEAGE_GROUPS_CHECK`, `LINEAGE_SUMMARY` and
  `BOOT_GROUP_INTERVALS` take a dataset channel that only emits when `freyja_data.json` has a `lineages_yml_sha256`; otherwise they
  do not run and one warning is logged. SARS-CoV-2 keeps the original wiring (`ch_hierarchy_dataset` is `FREYJA_DATA.out.dataset`).
- **Primer scheme check**: `bundledSchemesApply()` (the default bundled scheme's contig is in the FASTA) gates `READ_STARTS` /
  `PRIMER_SCHEME_CHECK`; otherwise the check, which scores against the bundled SARS-CoV-2 schemes, is skipped.
- **Unexplained mutations**: the module passes `--genome-length <reference length>`; `unexplained_mutations.py` rescales its default
  significance level (the pipeline's 1.67e-6 = 0.05 / 29903) to 0.05 / length for another genome length. Equal for SARS-CoV-2; an
  explicit `--unexplained_max_pvalue` other than the default is never rescaled.
- `bin/amplicon_coverage.py` understands primalscheme `<scheme>_<n>_LEFT_<alt>` / `_RIGHT_<alt>` names (the ARTIC measles panel;
  before, all primers collapsed into one amplicon). No bundled SARS-CoV-2 BED has such names.
- The report subtitle reads `Wastewater <pathogen> lineage deconvolution`, with the pathogen taken from the report header data (falling back to SARS-CoV-2).

## Not changed

Genome length 29903 remains only as the SARS-CoV-2 default of `primer_scheme_check.py` (`--genome-length`; the check is skipped for
other references) and `unexplained_mutations.py`; `consensus_stats.py`, `low_depth_mask.py` and the coverage code derive lengths from
the reference or the depths. The `container` snapshot, the bundled primer schemes and lineage groups stay SARS-CoV-2 only.

## Tests

Three new pytest tests (kept to the minimum on purpose): the profile/schema consistency guard (`test_boolean_params.py`), the
other-pathogen dataset names (`test_freyja_dataset_info.py`) and the unexplained-mutations dataset names plus significance scaling
(`test_unexplained_mutations.py`); one extra assertion in `test_amplicon_coverage.py`. No new nf-test and no new e2e step. The measles
validation is a manual run recorded in [../other_pathogens.md](../other_pathogens.md).

## RSV-A and RSV-B (PLAN D23): what changed

Two profiles and four reference files; no pipeline code changed.

- `nextflow.config`: the profile map is now `['SARS-CoV-2': 'sars2', 'MEASLESgenome': 'measles', 'RSVa': 'rsva', 'RSVb': 'rsvb']`.
- `conf/pathogens/rsva.config` / `rsvb.config` (pattern of `measles.config`): `pathogen`, `fasta = assets/references/PP109421.1.fasta` /
  `OP975389.1.fasta`, `gff = ...gff3`, and the two hidden facts as `null` (no lineage groups; the container does not bundle RSV).
- `assets/references/`: the GenBank FASTA records (NCBI efetch, unchanged) and GFF3 files. The GFF3 is NCBI's with one synthetic
  `mRNA` and `exon` line added per gene (and the NCBI `gene` / `CDS` lines rewritten with `gene:` / `transcript:` IDs), because
  `bcftools csq` stops with "No usable transcripts" on a file that has only `gene` and `CDS` lines. `read_cds_models()` in
  `bin/mutations_of_interest.py` reads the converted file as it does the SARS-CoV-2 one (feature name = the transcript's Name).
- Existing D22 behaviour that was exercised for the first time with a pathogen that ships no `lineages.yml`: `latest` downloads
  `<P>_barcodes.csv` only; `freyja_data.json` records `lineages_yml_sha256: null`; the lineage summary is off with one `INFO` line
  (no groups bundled) or one `WARN` line (`--lineage_groups` given); dropout impact and unexplained mutations read the barcodes
  CSV; the default-primer-scheme start-up error fires for `qiagenv1` with the RSV contig name; per-sample `protocol` /
  `primer_scheme` (BED path) cells mix capture and amplicon samples in one run. All behaved as designed; nothing needed changing.
- Tests: the guard `test_pathogen_profiles_are_consistent` (`tests/python/test_boolean_params.py`) now also checks that a profile's
  bundled FASTA exists and has one sequence and that a bundled GFF3 is on that sequence. No other test was added; the validation is
  the manual run in [../other_pathogens.md](../other_pathogens.md), made with `tests/simulate_pathogen_mixture.py` (not run by
  any test).
- Not changed but found: the console banner (`assets/logo.txt`) names SARS-CoV-2 for every pathogen; the influenza refusal text
  (`unsupportedPathogenReason()`, the schema help text) still gave only the "8 segments" reason although Freyja's influenza
  barcodes are HA-only. Resolved by PLAN D24 below.

## Influenza, HA segment only (PLAN D24): what changed

Five profiles, four reference records, one start-up helper; the single-contig pipeline is otherwise untouched.

- `nextflow.config`: the profile map gained `'H1N1': 'h1n1', 'H3N2': 'h3n2', 'FLU-B-VIC': 'flu_b_vic', 'H5Nx': 'h5nx', 'H5Nx-cattle': 'h5nx_cattle'`.
  `conf/pathogens/<profile>.config` (pattern of `rsva.config`): `pathogen`, `fasta`, `gff` (the `h5nx` and `h5nx_cattle` files name the
  same `AF144305.1` record) and the two hidden facts as `null`.
- `assets/references/`: `CY121680.1` (H1N1), `CY163680.1` (H3N2), `KX058884.1` (FLU-B-VIC), `AF144305.1` (H5Nx, H5Nx-cattle): the NCBI
  FASTA (efetch, unchanged; sequence-identical to Freyja's `reference.fasta`) and a GFF3 that is NCBI's gene + CDS lines with an `mRNA` and
  an `exon` line added per gene (`bcftools csq`; coordinates unchanged; HA is one CDS per record). The pytest profile guard now also
  asserts that a bundled GFF3 has gene, mRNA, exon and CDS lines.
- `unsupportedPathogenReason()` keeps `MTB` only; the schema `help_text` / `errorMessage` and `docs/parameters.md` say HA-only.
- **Primer BED with several contigs**: `reconcileBed(bed)` (`utils_nfcore_wasteline_pipeline`) reads the reference's contigs
  (`fastaContigs`) and the BED; lines whose contig is not in the reference are dropped, comment / blank / `track` / `browser` lines are
  kept, and when anything was dropped (and something is left) the filtered text is written to
  `<workflow.workDir>/wasteline_primer_beds/<first 12 hex of its sha256>/<BED file name>` unless it exists (idempotent and deterministic,
  so `-resume` keeps its cache key). It returns `[bed, no_match, dropped_lines, dropped_contigs]`; a BED with no dropped line is
  returned as the very same Path (SARS-CoV-2 and every matching BED: the staged file and its checksum are unchanged). `sourceBedFor(meta)`
  is the BED as given; `primerBedFor(meta)` (used by `IVAR_TRIM`, `AMPLICON_COVERAGE`, `PROVENANCE` and now also the custom-BED input of
  `PRIMER_SCHEME_CHECK`) is `reconcileBed(sourceBedFor(meta)).bed`. Validation: `checkSampleBeds` (samplesheet rows) refuses a BED with
  `no_match` (the previous message) and calls `warnDroppedBedLines` once for all distinct BEDs; `validateInputParameters` does the same
  for the run-default BED and warns itself only when there is no samplesheet (`!params.input`), so there is one warning per run. The
  previous check refused a BED when any of its contigs was missing; now only a BED with no line on the reference is refused.
- The multi-sequence FASTA refusal stays (any pathogen); its message now points to the bundled influenza HA references.
- Nothing in the QC, report or consensus code changed. Reviewed for a 1.7 kb reference: `qc_summary.py` has no rule on `mapped_pct`;
  the coverage rules are percentages of the reference; `build_report_data.py` bins the coverage track by reference length;
  `unexplained_mutations.py` already scales its significance level; `consensus_*` and `low_depth_mask` take lengths from the data.
  What does not scale is the negative-control rule (`--ntc_max_mapped_reads` 1000 reads, `--ntc_max_genome_coverage_pct` 10 %): on a
  1.75 kb reference two read pairs exceed 10 % of the genome. Left as is and documented; a per-length default would be a design change.
- Tests: none added (the existing profile guard was extended by two assertions); the BED helper is Groovy and was checked manually
  (stub run with a BED listing `seg4` / `seg6`: warning, 14-line filtered file, run_info checksum of the filtered BED; BED with only
  `seg4`: the old error; full docker runs in [../other_pathogens.md](../other_pathogens.md)). The simulator
  `tests/simulate_pathogen_mixture.py` is unchanged.

## Mpox, per-pathogen primer schemes and user schemes (PLAN D26): what changed

User view: [../other_pathogens.md](../other_pathogens.md) ("Running mpox", "Primer schemes"), [../usage.md](../usage.md#using-your-own-primer-scheme).

**Mpox.** One profile and two reference files; no pipeline code changed for the pathogen itself.

- `nextflow.config` profile map gained `'MPX': 'mpx'`; `conf/pathogens/mpx.config` sets `pathogen = 'MPX'`, `fasta` / `gff` = `assets/references/NC_063383.1.{fasta,gff3}`.
- The FASTA is the NCBI record `NC_063383.1` (efetch), sequence-identical to Freyja's `MPX/latest/reference.fasta` (header `MPXV-M5312_HM12_Rivers`); the
  GFF3 is NCBI's gene + CDS lines (179 genes, `OPGnnn`) with `mRNA` and `exon` lines added, like the RSV/influenza files. The four genes inside the inverted
  terminal repeats occur twice in the record; ids and names must be unique for `bcftools csq` and `read_cds_models()`, so the copy near the 3' end has the suffix
  `_R` (`OPG003_R`). The NCBI `long_terminal_repeat` and `sequence_feature` lines are not kept.
- **What broke on a 197 kb reference**: `modules/nf-core/bwa/index/main.nf` declares `memory { 7.B * fasta.size() }`. Nextflow passes `--memory 1m` to Docker for a
  1.4 MB value, Docker's minimum is 6 MB, so `BWA_INDEX` died with exit 125 on the mpox FASTA (by the same arithmetic a FASTA of roughly 150-860 kb lands between 1 MB and Docker's 6 MB minimum; the smaller references that were run never did).
  `conf/modules.config` overrides the directive for `BWA_INDEX` (`memory = { 7.B * fasta.size() > 256.MB ? 7.B * fasta.size() : 256.MB }`); a config selector wins over the
  module's own directive (checked), and the nf-core module is untouched. Nothing else in the pipeline depended on the genome length: see the checks in "Running mpox".

**Per-pathogen bundled primer schemes.** Two more hidden profile facts (schema `hidden: true`, filtered out of `run_info.json` in `runInfoBase()` with the other two):

| Fact                     | SARS-CoV-2                                        | measles / mpox                                            | others |
| ------------------------ | ------------------------------------------------- | --------------------------------------------------------- | ------ |
| `primer_schemes_bundled` | `${projectDir}/assets/primer_schemes` (as before) | `${projectDir}/assets/primer_schemes/measles` / `.../mpx` | `null` |
| `primer_scheme_default`  | `'qiagenv1'` (PLAN D14)                           | `null`                                                    | `null` |

`subworkflows/local/utils_nfcore_wasteline_pipeline/main.nf`:

- `defaultPrimerScheme()` returns `params.primer_scheme_default` (was the literal `'qiagenv1'`); `bundledPrimerSchemes()` lists the `*.bed` of `primer_schemes_bundled`;
  `sourceBedFor()` and `bundledSchemesApply()` read the same directory. For SARS-CoV-2 every path resolves to the same files as before, so its `meta`, commands,
  `run_info.json` (the two new facts are filtered) and `PRIMER_SCHEME_CHECK` wiring are unchanged. `PRIMER_SCHEME_CHECK` still stages the top-level `assets/primer_schemes`
  (it globs `*.bed` non-recursively, so the new subdirectories are not scored) and `bundledSchemesApply()` is false for a pathogen without a default scheme.
- `libraryMeta()` resolves a scheme name (cell, `--primer_scheme`) with `chooseScheme(name, who)`; with no default scheme and nothing given in mode `full` it stops with
  `noDefaultSchemeMessage()` naming the sample and the schemes bundled for the pathogen (this replaces the earlier "default scheme is for a different reference"
  error for pathogens without a default; that message remains for a SARS-CoV-2 run with a `--fasta` of another sequence). Because the check moved into
  `libraryMeta()` it is per sample, so a samplesheet that names a scheme for every sample works without a run default. `validateInputParameters()` only checks the BED of a
  run default that exists (`primerBed()` returns `null` otherwise).
- `chooseScheme()` also produces the "unknown scheme" error: it lists the bundled names, the names in `--primer_scheme_dir`, and says when the name is bundled for another
  pathogen (`schemesOfOtherPathogens()`). The schema enum of `--primer_scheme` was removed (names are pathogen dependent); the schema now has a name pattern.
- **`--primer_scheme_dir <dir>`** (`nextflow.config`, schema `format: directory-path`, `exists: true`): `chooseScheme()` looks in the directory first. A user scheme is stored
  in `meta` as `primer_scheme: 'custom'` with `primer_bed` = the absolute path, i.e. it takes the route of a per-sample BED path: `IVAR_TRIM`, `AMPLICON_COVERAGE`, `PROVENANCE`
  (sha256, `custom: <name>.bed` in `sample_libraries`), `PRIMER_SCHEME_CHECK` (custom BED) and the contig check (`checkSampleBeds`, `validateInputParameters`) need no change.
  `validateInputParameters()` warns once about names that exist in both places (user file used), about an empty directory, and about use outside modes `full`/`demix`.
- `bin/check_samplesheet.py` gained `--pathogen` (reads the profile from the map in `nextflow.config`) and `--primer-scheme-dir`; its name check uses the pathogen's directory
  (`assets/primer_schemes/<profile>/`). `bin/replay_run.py` finds a recorded bundled BED in a profile subdirectory (`rglob`); `bin/reanalyse_runs.py` carries `primer_scheme_dir`.
- `assets/primer_schemes/{measles,mpx}/*.bed` are unmodified `primer.bed` files (md5 equal to the scheme's own record) with a `# artic-bed-version v3.0` first line; `ivar trim`
  (0 changes in the measles results), `amplicon_coverage.py` and the contig check skip it. Sources, licences and the schemes that could not be bundled: `assets/primer_schemes/README.md`.
- Tests: no new pytest test. `test_pathogen_profiles_are_consistent` checks the two new facts, that a profile's scheme directory has BEDs on the profile's reference contig and that a
  default scheme is one of them; `test_per_sample_protocol_and_primer_scheme` (samplesheet) also covers per-pathogen names and a user directory. No nf-test, no e2e step.
  `tests/simulate_pathogen_mixture.py` gained `--shear` (fragments inside long amplicons).
