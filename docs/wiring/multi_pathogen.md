# Other pathogens (PLAN D22): wiring as built

`--pathogen X` runs the pipeline for the single-sequence Freyja pathogens; SARS-CoV-2 behaviour is unchanged (same params
recorded, same commands, same outputs). User view: [../other_pathogens.md](../other_pathogens.md).

## Profile selection

`nextflow.config` ends with a dynamic include chosen by `params.pathogen` (a command-line or `-params-file` value is visible
while the config is parsed; verified with `--pathogen H1N1` and `--pathogen MEASLESgenome`):

```groovy
includeConfig "conf/pathogens/${['SARS-CoV-2': 'sars2', 'MEASLESgenome': 'measles'].get(params.pathogen ?: 'SARS-CoV-2', 'generic')}.config"
```

| File                            | Sets                                                                                                                  |
| ------------------------------- | --------------------------------------------------------------------------------------------------------------------- |
| `conf/pathogens/sars2.config`   | `pathogen`, `fasta`, `gff` (unchanged) plus the two facts below                                                       |
| `conf/pathogens/measles.config` | `pathogen = 'MEASLESgenome'`, `fasta = assets/references/NC_001498.1.fasta`, `gff = null`                             |
| `conf/pathogens/generic.config` | `pathogen = null` (the strict parser needs the key declared; the command line fills it), `fasta = null`, `gff = null` |

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

1. `unsupportedPathogenReason()`: influenza (`H5Nx`, `H5Nx-cattle`, `H1N1`, `H3N2`, `FLU-B-VIC`) and `MTB` stop with the reason and what
   would be needed (multi-contig handling; scaling to 4.4 Mb).
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
