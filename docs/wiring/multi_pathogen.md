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
- The report subtitle is `Wastewater {{ h.pathogen or 'SARS-CoV-2' }} lineage deconvolution`.

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
