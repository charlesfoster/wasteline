# Helper scripts

Standard-library Python scripts in `bin/`. Run them with `python bin/<script>.py` (or via `pixi run python ...`).
Sample IDs are never rewritten: underscores and other valid characters are kept exactly as found.

## `make_samplesheet.py`

### From a FASTQ directory (mode `full`)

```bash
python bin/make_samplesheet.py reads/ -o samplesheet.csv
python bin/make_samplesheet.py reads/ --suffix _R1.fq.gz --metadata meta.csv --recursive -o samplesheet.csv
```

Output columns: `sample,fastq_1,fastq_2,sample_type,collection_date,site` plus `viral_load` when any sample has one (sorted by `sample`; absolute paths unless
`--relative`, which writes paths relative to the output file's directory).

| Option               | Default                 | Meaning                                                                                                                                  |
| -------------------- | ----------------------- | ---------------------------------------------------------------------------------------------------------------------------------------- |
| `--suffix`           | `_L001_R1_001.fastq.gz` | R1 file suffix; must end `.fastq.gz`/`.fq.gz`. An `R2` suffix is treated as `R1`.                                                        |
| `--negative-pattern` | `NEG_\|NC_`             | Regex searched in the file name; a match gives `sample_type=negative`.                                                                   |
| `--exclude-pattern`  | `Undet`                 | Regex searched in the file name; matching files are skipped (e.g. `Undetermined_*`).                                                     |
| `--metadata`         | none                    | CSV keyed by `sample`; may supply `collection_date`, `site`, `viral_load` (a positive number; blank allowed) and override `sample_type`. |
| `--recursive`        | off                     | Search sub-directories.                                                                                                                  |
| `--relative`         | off                     | Relative instead of absolute paths.                                                                                                      |
| `-o`                 | stdout                  | Output CSV.                                                                                                                              |

Sample ID rules:

- Default Illumina suffix: `<id>_S<n>_L001_R1_001.fastq.gz` gives `<id>`, e.g.
  `Site_1_week_2_S12_L001_R1_001.fastq.gz` -> `Site_1_week_2`.
- Any other suffix: the ID is the file name with the suffix removed (`alpha_1_R1.fq.gz` -> `alpha_1`).
- R2 is derived from the suffix (`R1` -> `R2`, else `_1` -> `_2`) and must exist.

It fails (exit 1, nothing written, all problems listed) on: missing R2, duplicate IDs after stripping (e.g.
`dup_S1_...` and `dup_S7_...`), IDs not matching `^[A-Za-z0-9._-]+$`, no samples found, bad `--suffix`, and, in the
metadata, dates that are not real ISO `YYYY-MM-DD` dates, invalid `sample_type`, duplicate or empty sample keys.
It warns (stderr) about metadata rows with no matching sample and about samples with no metadata row.

Metadata example:

```csv
sample,collection_date,site,sample_type
Site_1_week_2,2024-03-01,SiteOne,
```

### From existing results (modes `demix` and `aggregate`)

```bash
python bin/make_samplesheet.py --from-results old_run/ --mode demix -o demix_sheet.csv
python bin/make_samplesheet.py --from-results old_run/ --mode aggregate -o aggregate_sheet.csv
```

- `demix` sheet: `sample,variants,depths,fastp_json,sample_type,collection_date,site[,viral_load]`
- `aggregate` sheet: `sample,demix,sample_type,collection_date,site[,viral_load]`

Both layouts are detected automatically (per sample):

|            | New (wasteline)                                             | Legacy (Snakemake)                                                                   |
| ---------- | ----------------------------------------------------------- | ------------------------------------------------------------------------------------ |
| variants   | `samples/<s>/freyja/<s>.variants.tsv`                       | `<s>/freyja/<s>.freyja_variants.tsv`, else `<s>/variants/<s>.freyja_variants.tsv.gz` |
| depths     | `samples/<s>/freyja/<s>.depths`                             | `<s>/freyja/<s>.depths`                                                              |
| demix      | `samples/<s>/freyja/<s>.demix` (`.demix.tsv` also accepted) | `<s>/freyja/<s>.demix`                                                               |
| fastp JSON | `samples/<s>/fastp/<s>.fastp.json`                          | `<s>/fastp/<s>_fastp.json`                                                           |

- If both layouts exist in one directory a warning is printed; a sample present in both uses the new layout.
- Samples missing required files (variants+depths for `demix`, demix for `aggregate`) are skipped with a printed
  `NOTE: skipped ...` reason; if nothing usable remains the command fails.
- gzipped variants (legacy `.tsv.gz`) are emitted as-is with a note; check the pipeline can read them.
- `sample_type=negative` is inferred with `--negative-pattern` on the sample ID. `sample_type`, `collection_date` and
  `site` (and a valid `viral_load`) are carried over from `<dir>/samplesheet.csv` or `<dir>/summary/samplesheet.csv` when present
  (the sheet's `sample_type` wins over inference), and `--metadata` can also be given.
- Paths are absolute (or relative to the output file with `--relative`).
- Results of a panel run (`summary/targets.tsv`, `samples/<s>/<target>/freyja/`) give one row per (sample, target) and a `target` column after
  `sample`; `--targets A,B` keeps only those targets.

#### Several runs, filters, `--list`

`--from-results` is repeatable: the rows of all the directories go into one sheet (one header; the columns are the union, blank where a run has
none). The single-directory output is unchanged. Filters are applied to the discovered rows and combined with AND; site, date and sample type are
the values carried over from the run's `summary/samplesheet.csv` (or `--metadata`).

| Option                        | Meaning                                                                                                                                                                                                                                                                            |
| ----------------------------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `--samples ID[,ID...]\|@file` | Keep only these exact IDs (a file: one ID per line, blank lines and `#` lines ignored). In panel results all targets of a sample. IDs not found give a warning.                                                                                                                    |
| `--include-id REGEX`          | Keep IDs matching the regex (searched, not anchored). Not the FASTQ-mode `--exclude-pattern`.                                                                                                                                                                                      |
| `--exclude-id REGEX`          | Drop IDs matching the regex.                                                                                                                                                                                                                                                       |
| `--site A,B`                  | Keep these sites (exact). Rows without a site are dropped.                                                                                                                                                                                                                         |
| `--date-from`, `--date-to`    | Keep `collection_date` within `YYYY-MM-DD` bounds (inclusive). Rows without a date are dropped when either is given.                                                                                                                                                               |
| `--sample-type`               | `sample` or `negative`.                                                                                                                                                                                                                                                            |
| `--status ok`                 | Keep samples whose `summary/demix_status.tsv` (`summary/<target>/demix_status.tsv` for panel results) says `ok`; a sample with no status is kept, with a note.                                                                                                                     |
| `--on-duplicate`              | The same (sample, target) in two directories, checked after the filters. `error` (default): exit 1 naming the directories; `prefix`: rename to `<dirname>__<id>` in the `sample` column only (directory names must be distinct and use `[A-Za-z0-9._-]`); `first`: keep the first. |
| `--list`                      | Print a table (run, sample, target when present, sample_type, site, collection_date, demix_status, selected, reason) to stdout and write nothing. Unusable samples (missing files) are listed with the reason.                                                                     |

Everything removed by a filter is counted in one `NOTE: N sample(s) removed by filters` line. If nothing is left the command fails.
Mixing panel and single-target results works but gives the single-target rows a blank `target` (a warning says so).

```bash
# Run1 samples 1-3 and Run2 samples 8-10, clean demixes only
python bin/make_samplesheet.py --from-results Run1/ --from-results Run2/ --mode aggregate \
    --include-id '^(S[1-3]|S([89]|10))$' --status ok -o aggregate_sheet.csv
python bin/make_samplesheet.py --from-results Run1/ --from-results Run2/ --mode aggregate --list
```

## `check_samplesheet.py`

```bash
python bin/check_samplesheet.py samplesheet.csv --mode full --check-files
```

Prints one JSON object on stdout and a human-readable table on stderr; exit code 0 if valid, 1 otherwise.

```json
{
  "ok": false,
  "mode": "full",
  "n_samples": 3,
  "errors": ["line 3 (B): empty fastq_2 (paired-end data is required)"],
  "warnings": []
}
```

Rules:

| Mode        | Required columns         | Optional                                                 |
| ----------- | ------------------------ | -------------------------------------------------------- |
| `full`      | `sample,fastq_1,fastq_2` | `sample_type,collection_date,site,viral_load`            |
| `demix`     | `sample,variants,depths` | `fastp_json,sample_type,collection_date,site,viral_load` |
| `aggregate` | `sample,demix`           | `sample_type,collection_date,site,viral_load`            |

- Unknown columns give a warning. Duplicate header names are an error.
- `sample` matches `^[A-Za-z0-9._-]+$` and is unique.
- `sample_type` is `sample` or `negative`; `collection_date` is a real ISO `YYYY-MM-DD` date; `viral_load` is blank or a positive finite number.
- Paired-end is required: an empty `fastq_2` is an error. FASTQ paths must end `.fastq.gz`/`.fq.gz`, contain no
  whitespace, and `fastq_1` and `fastq_2` must differ.
- `--check-files` verifies that every referenced file exists (relative paths resolve against the sheet's directory),
  FASTQs are readable gzip whose first line starts with `@`, and other required files are not empty. Without it the
  file system is not touched.
