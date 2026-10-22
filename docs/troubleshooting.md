# Troubleshooting

Each entry is symptom, cause, fix. The messages quoted below were produced by real runs (Nextflow 26.04.6, Docker 29,
Freyja 2.0.5, macOS on Apple Silicon) or copied from the pipeline's source; the wording of messages that come from
Nextflow, Docker or nf-schema may differ slightly in other versions. First stop for any failure: the **last lines of the
terminal output** (they name the failed process and its `Work dir`), then `.nextflow.log` in the launch directory, then the
`.command.err` / `.command.log` / `.command.sh` files in that work directory.

Contents:
[Start-up and environment](#start-up-and-environment) |
[Samplesheet and parameters](#samplesheet-and-parameters) |
[Freyja data](#freyja-data) |
[Per-sample results](#per-sample-results) |
[Report and modes](#report-and-modes) |
[Resuming and harmless warnings](#resuming-and-harmless-warnings)

## Start-up and environment

### Script compilation error on an old Nextflow

**Symptom** (verified with Nextflow 25.10.4):

```
ERROR ~ Script compilation error
- file : /.../wasteline/main.nf
- cause: Unexpected input: '{' @ line 36, column 8.
   params {
          ^
```

**Cause.** wasteline needs Nextflow **>= 26.04.0**. Its `main.nf` declares typed parameters in a `params { ... }` block so that
boolean and numeric command-line flags work (`--skip_consensus`, `--save_bam false`, `--consensus_min_depth 20`). Nextflow
before 26.04 cannot parse that block, and on Nextflow 26.04's strict parser untyped command-line flags would arrive as
strings and be rejected by the schema.

**Fix.** Upgrade: `nextflow self-update`, or run a single command with a specific version, for example
`NXF_VER=26.04.6 nextflow run charlesfoster/wasteline ...`. `nextflow -version` shows what you have.

### Docker is not running

**Symptom.** The run starts, then every task fails at once and the run ends with `completed with errors`. The task error
contains (Docker 29; older Docker versions say `Cannot connect to the Docker daemon at unix:///var/run/docker.sock. Is the
docker daemon running?`):

```
stderr: failed to connect to the docker API at unix:///tmp/nodocker.sock; check if the path is correct and if the daemon is running: dial unix /tmp/nodocker.sock: connect: no such file or directory
[FAILED] completed=9 failed=9 cached=0
```

(Reproduced by pointing `DOCKER_HOST` at a socket that does not exist.)

**Fix.** Start Docker Desktop (macOS/Windows) or `sudo systemctl start docker` (Linux), check `docker info`, and re-run with
`-resume`. If Docker needs `sudo`, add your user to the `docker` group instead of running Nextflow with `sudo`.

### It is very slow on an Apple Silicon Mac

**Cause.** The Freyja and most bioconda images exist only for `linux/amd64`. On an arm64 host Docker emulates amd64: running
an image by hand prints `WARNING: The requested image's platform (linux/amd64) does not match the detected host platform
(linux/arm64/v8)`. With `-profile docker` the pipeline adds `--platform=linux/amd64` for you, so it works, but under
emulation. The synthetic dataset takes about two minutes this way; real runs are proportionally slower, with `BWA_MEM` and
`FREYJA_VARIANTS` (samtools mpileup) hurting most.

**Fix.** There is nothing to configure for correctness. For speed, give Docker Desktop more CPUs and memory, and (if your
Docker Desktop offers it) enable Rosetta for amd64 emulation
<!-- UNVERIFIED: the Docker Desktop Rosetta setting was not tried here -->.

Run large datasets on Linux/amd64. `-profile emulate_amd64` forces the platform flag on any Docker host. Do not use conda on
macOS.

### Pulling images fails or hangs

**Cause.** The first run pulls several images (Freyja, fastp, bwa, samtools, ivar, bcftools, Python) from `quay.io` and `community.wave.seqera.io` (Singularity/Apptainer also use
`depot.galaxyproject.org` and `community-cr-prod.seqera.io`). A firewall, proxy or no network stops it.
**Fix.** Pull from a machine with access, or mirror the registries and set `docker.registry` (infrastructure setting, in a `-c`
file). The images are listed in `pipeline_info/execution_trace_*.txt` (`container` column) and in
`conf/containers.lock.config`. <!-- UNVERIFIED: proxy/mirror setups were not tested -->

### `TypeError: zip() takes no keyword arguments` or a similar error from a test script

**Cause.** `tests/check_outputs.py` (a developer test script, not a user tool) needs Python >= 3.12; the macOS system `python3` is 3.9. The user-facing helpers `bin/make_samplesheet.py` and `bin/check_samplesheet.py` run on Python >= 3.9.
**Fix.** Use the pixi environment (`pixi install` then `pixi run python tests/check_outputs.py ...`). The pipeline itself is
unaffected: it runs its scripts inside containers.

## Samplesheet and parameters

### Samplesheet: the file or directory does not exist

**Symptom** (the files are there, but the paths in the sheet are relative and you launched from another directory):

```
ERROR ~ Validation of samplesheet failed!
-> Entry 1: Error for field 'fastq_2' (tests/data/fastq/MIX_A_R2.fastq.gz): the file or directory 'tests/data/fastq/MIX_A_R2.fastq.gz' does not exist (FastQ file for reads 2 must be provided, cannot contain spaces and must have extension '.fq.gz' or '.fastq.gz')
```

**Cause.** Relative paths in a samplesheet are resolved against the directory you **launch Nextflow from**, not against the
samplesheet's own directory. The synthetic test sheet `tests/data/samplesheet.csv` uses paths relative to the repository root.
**Fix.** Use absolute paths (`make_samplesheet.py` writes absolute paths unless you pass `--relative`), or launch from the directory
the paths are relative to. For the test profiles, launch from the repository root.

### Samplesheet validation errors

| Message                                                                                       | Cause and fix                                                                                                                                                                                                                                                 |
| --------------------------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `Entry 1: Missing required field(s): fastq_2`                                                 | single-end data, or an empty `fastq_2`. Only paired-end reads are supported.                                                                                                                                                                                  |
| `Entry 1: Missing required field(s): fastq_2, fastq_1` on every row                           | the wrong sheet for the mode: a demix/aggregate sheet given to `--mode full`. Conversely `Missing required field(s): depths, variants` means a FASTQ sheet was given to `--mode demix`. Required columns per mode are in [usage.md](usage.md#3-samplesheets). |
| `Entry 2: Detected duplicate entries: [sample:A]`                                             | two rows with the same `sample`. IDs must be unique; wasteline does not merge repeats. Concatenate the FASTQs of one sample first.                                                                                                                            |
| `Error for field 'sample' (A B): "A B" does not match regular expression [^[A-Za-z0-9._-]+$]` | an ID with a space or other disallowed character. IDs are never rewritten, so rename it (letters, digits, `.`, `_`, `-`).                                                                                                                                     |
| `must have extension '.fq.gz' or '.fastq.gz'`                                                 | FASTQs must be gzipped and named `.fastq.gz` or `.fq.gz`, with no spaces in the path.                                                                                                                                                                         |
| `collection_date must be an ISO date (YYYY-MM-DD)`                                            | fix the date or leave it empty.                                                                                                                                                                                                                               |

`python bin/check_samplesheet.py samplesheet.csv --mode full --check-files` reports all problems in one go.

### Missing `--outdir`, unknown parameters, wrong values

```
ERROR ~ Validation of pipeline parameters failed!
* Missing required parameter(s): outdir
```

`--outdir` is required in **every** mode (including `report` and `update_db`). Other messages:

```
[ERROR] Parameter `primer_sheme` was specified on the command line or params file but is not declared in the script or config
[ERROR] Parameter `skip_consensus` with type Boolean cannot be assigned to maybe [String]
[ERROR] --primer_scheme 'foo' is not a primer scheme available for --pathogen SARS-CoV-2 and not a path to a BED file. Bundled: articv3, articv4.1, ...; give your own with --primer_bed <file> or --primer_scheme_dir <dir of <name>.bed files>
```

A typo in a parameter name, a non-boolean value for a boolean flag, or a value outside an enum. Pipeline parameters take **two**
dashes; Nextflow options (`-profile`, `-resume`, `-params-file`) one. Do not pass parameters in a `-c` config file.
Boolean flags accept `--x`, `--x true`, `--x false`, `--x=false` and `-params-file` entries on Nextflow >= 26.04.

### The mode and the input do not match

```
[ERROR] No input given for --mode demix. --mode demix needs existing Freyja variants and depths: --input samplesheet.csv (...) or --previous_results <results dir> (...)
[ERROR] --previous_results cannot be used with --mode full. --mode full needs a FASTQ samplesheet: --input samplesheet.csv (...)
[ERROR] Give only one of --input, --previous_results with --mode demix. ...
[ERROR] --results_dir is only used by --mode report (this is --mode full). ...
```

`--mode` is the only dispatch (there is no `-entry`). The accepted inputs per mode are in [usage.md](usage.md#2-modes-at-a-glance).
Samples with the same ID in two `--previous_results` directories are also refused: `Sample ID(s) found more than once across
the inputs: ...`.

### A custom primer BED does not match the reference

```
[ERROR] Primer BED '/path/my_scheme.bed' refers to contig(s) [MN908947.3] that are not in the reference FASTA '.../NC_045512.fasta' (which has [NC_045512.2]). Use a matching --primer_bed/--primer_scheme or --fasta.
```

**Cause.** None of the lines of a custom `--primer_bed` is on a contig of the reference (lines on other contigs alone are only
dropped, with one warning: `Primer BED: 10 line(s) on 2 contig(s) that are not in the reference FASTA ... were dropped`, which is
what a whole-genome influenza panel gives against the HA reference). A common case is a
BED made for the `MN908947.3` accession while the bundled reference is `NC_045512.2`. All nine bundled schemes already use
`NC_045512.2`, so this only happens with your own BED (or your own `--fasta`). **Fix.** Rename the contig in the BED to the
reference's name (the coordinates are the same when it is the same sequence), for example
`sed 's/^MN908947\.3/NC_045512.2/' my_scheme.bed > my_scheme.NC_045512.2.bed`, or pass the matching `--fasta`. Details:
[usage.md](usage.md#4-reference-and-primer-schemes). The check exists on purpose: with a mismatched contig `ivar trim` would
silently trim nothing. (Reproduced with a copy of `articv4.1.bed` with the contig renamed.)

### `primer_scheme_check` says `mismatch` (or the default scheme is not mine)

The default primer scheme in mode `full` is `qiagenv1` (QIAseq Direct SARS-CoV-2; the old Snakemake pipeline used `midnight`).
If your library is another protocol and you did not pass `--primer_scheme`, trimming and amplicon-dropout calls use the wrong
scheme; the pipeline notices when the reads begin at primers (`primer_scheme_check` = `mismatch`, `primer_scheme_best` names the
scheme that fits, the report shows the flag `primer_scheme_mismatch`). **Fix.** Re-run with `--primer_scheme <best>` (the
effective scheme is in `pipeline_info/run_info.json`). `uninformative` is not a problem: the reads do not start at primers (for
example tagmentation libraries), so the check cannot tell. Details: [usage.md](usage.md#primer-scheme-check).

### The dropout-impact and unexplained-mutation tables are missing

A single warning at start-up says `Skipping the dropout-impact and unexplained-mutation tables`: the Freyja dataset in effect has
no readable lineage-to-mutation lists. They need `lineage_mutations.json` (every dataset downloaded or bundled with Freyja 2.x
has it) or a barcodes `.csv`; a custom `.feather`-only `--barcodes` has neither. **Fix.** Point `--barcodes` at a dataset directory
(or a `.csv`). `dropout_impact.tsv` also needs an amplicon table (amplicon protocol; in `--mode demix` an explicit `--primer_scheme`).
Errors inside these analyses never fail the run: the table is just absent.

### `freyja_growth_rates.failed.txt`, `freyja_dashboard.failed.txt` or no time series

The time-series tables need at least two distinct collection dates among the non-control samples. The growth rates and the
dashboard need `--freyja_growth_min_dates` (4) distinct dates within the last `--freyja_growth_days` (56) days spanning at least
`--freyja_growth_min_span_days` (14) days; otherwise the `.failed.txt` file says what was missing (the run is not failed). Samples without
a `collection_date`, negative controls and samples below the plot coverage threshold (`--freyja_plot_mincov`, default: the demix threshold) are not used. The dashboard page fetches its
styling from the internet when opened; offline it looks plain.

### A lineage groups or sites file stops the run at the start

The groups file (`--lineage_groups`, or the bundled default) and the sites file (`--site_aggregation`) are checked by two small
processes, `LINEAGE_GROUPS_CHECK` and `SITE_AGGREGATE_CHECK`, **first thing in the run**, in parallel with the read processing. A
format error stops the run within minutes (before mapping finishes) with a message that names the file and the row or line. Look
for the `Command error:` block in the Nextflow output. Typical messages and what to do:

| Message (shortened)                                                                   | What is wrong and how to fix it                                                                                |
| ------------------------------------------------------------------------------------- | -------------------------------------------------------------------------------------------------------------- |
| `groups.csv: the group name 'JN.1*' appears twice (rows 2 and 3)`                     | Two rows with the same `group`. Group names are the column headings and must be unique (case does not matter). |
| `groups.csv, row 4: group 'X' has no 'include' pattern`                               | The `include` cell is empty. Every group needs at least one pattern.                                           |
| `groups.csv, row 2: in_total must be yes or no for group 'X', got 'maybe'`            | `in_total` is `yes`, `no` or empty.                                                                            |
| `groups.csv, row 3: 'other' is a reserved output column name`                         | `sample`, `site`, `collection_date`, `status`, `other`, ... are fixed columns of the output; rename the group. |
| `the groups file groups.csv needs the columns 'group' and 'include' in its first row` | The header row is missing or misspelled (an Excel title row above the header also causes this).                |
| `sites.csv: line 3: site_name is blank (aggregate 'Metro')`                           | A row without a site. Fill it in or delete the row.                                                            |
| `sites.csv: line 5: site 'site1' is listed twice for aggregate 'Metro_Sites'`         | The same site twice in one group (matching ignores case and extra spaces).                                     |
| `sites.csv: line 4: weight '1 200' must be a positive number ...`                     | Write plain digits in `weight` (`1200`): no units, spaces or thousands separators; or leave it blank.          |
| `sites.csv: line 1: unknown column 'site_names' (did you mean 'site_name'?)`          | A misspelled column name; the columns are `aggregate_name`, `site_name` and optional `weight`.                 |

Unknown lineages in the groups file and site names that no sample uses are only **warnings** (shown once in the log): they match
nothing, so a typo silently gives an empty column or a missing site in an aggregate; check
`summary/lineage_groups_resolved.tsv` (column `warnings`) and the log. Other common slips: saving as "CSV (Macintosh)" or "Text (tab
delimited)" instead of **CSV UTF-8** (accents or columns go wrong); writing `X*` for recombinants (use `RECOMBINANT*`); a comma inside a
pattern cell (use `;`); and a sample whose `site` is spelled differently from the sites file (`Site_1` vs `Site1`).

### The lineage summary is missing

With the default settings `summary/lineage_summary.csv` is written in modes `full`, `demix` and `aggregate`. If it is not: the run
log says `Skipping the lineage-group summary` when the Freyja dataset has no lineage hierarchy (a custom `--barcodes` `.csv` or
`.feather` file; use `--barcodes latest`, `container` or a dataset directory that contains `lineages.yml`), or you gave
`--lineage_summary false`. An `ungrouped_lineages` column or `lineage_ungrouped_summary.tsv` full of lineages means the groups list is
out of date for the data: see [usage.md](usage.md#editing-the-list).

### Warning: `--primer_scheme/--primer_bed are ignored`

```
WARN: --primer_scheme/--primer_bed are ignored: --protocol capture reads are not primer-trimmed and there is no amplicon-level analysis
```

**Cause.** You gave a primer scheme or BED together with `--protocol capture` or `--skip_primer_trimming`. Probe-capture
reads have no primers to trim and no amplicons, so the scheme is not used (and the start-up contig check is skipped).
**Fix.** Nothing is wrong; drop the option to silence the warning. If your data really are PCR amplicons, do not use
`--protocol capture`. See [usage.md](usage.md#probe-capture-data-eg-twist---protocol-capture).

### `--pathogen` is rejected

```
* --pathogen (Foo): Expected any of [[SARS-CoV-2, MPX, H5Nx, ...]] (--pathogen must be one of the Freyja pathogens: ...)
[ERROR] --pathogen MTB is not supported: Mycobacterium tuberculosis is a 4.4 Mb bacterial genome ...
[ERROR] --fasta is required with --pathogen RSVa: wasteline bundles a reference only for the pathogens that have a file in conf/pathogens/ ...
[ERROR] The bundled default primer scheme 'qiagenv1' is for a different reference ...
[ERROR] Sample 'S1' has no primer scheme: --pathogen MPX has no default primer scheme (only the pipeline's default pathogen has one). Name one with --primer_scheme ... (bundled for MPX: artic-inrb-mpox-2500-v1.0.0-cladeiib, bccdc-mpox-v2.3.0, yale-mpox-2000-v1.0.0-lifted)
```

**Cause.** `--pathogen` must be a Freyja pathogen name. MTB (4.4 Mb) is accepted by Freyja but not by wasteline (the influenza
pathogens run on their HA segment, a single 1.7-1.9 kb reference); a pathogen without a bundled reference needs `--fasta` (the reference the Freyja barcodes
were built on), and only SARS-CoV-2 has a default primer scheme (bundled schemes are per pathogen: measles one, mpox three, RSV-A and RSV-B two each), so a scheme name, a primer BED or
`--primer_scheme_dir` (or `--protocol capture`) is needed.
**Fix.** See [other_pathogens.md](other_pathogens.md) for what each pathogen needs.

### `--annotate_variants` with `--skip_consensus`

```
[ERROR] --annotate_variants annotates the variants the consensus is built from: it cannot be combined with --skip_consensus
```

**Cause.** The annotation runs on the consensus filter's VCF, which only exists when the consensus chain runs. **Fix.** Drop
`--skip_consensus` (the consensus is cheap), or drop `--annotate_variants`.

### A parameter is ignored in `aggregate`, `report` or `update_db` mode

```
WARN: --compare_barcodes is only used by --mode full and demix and is ignored in --mode aggregate
WARN: --freyja_boot needs variants and depths: it is only used by --mode full and demix and is ignored in --mode aggregate
WARN: --annotate_variants is only used by --mode full and demix and is ignored in --mode aggregate
WARN: --primer_scheme/--primer_bed are not used by --mode aggregate and are ignored
WARN: --input is not used by --mode report and is ignored
```

**Cause.** Those options need variants and depths (or reads), which these modes do not read; the run goes on without them
(reproduced for `aggregate`). **Fix.** Use `--mode demix --previous_results <results>` for `--compare_barcodes`, `--freyja_boot`
and `--annotate_variants`: it re-uses the variants and depths of the earlier run without mapping again.

## Freyja data

### `--barcodes` date or version not found

```
[ERROR] No Freyja dataset from 2020-01-01 found in <db_dir>/freyja. Available: [2026-10-06T05-25]. Download one with --barcodes latest, or point --barcodes at a dataset directory or .csv/.feather file.
[ERROR] Freyja dataset '09_01_2026-00-00' (directory 2026-09-01T00-00) not found in <db_dir>/freyja. Available: none. ...
[ERROR] --barcodes 'nope' is not 'latest', 'container', a barcode version (YYYY-MM-DD, MM_DD_YYYY-HH-MM or YYYY-MM-DDTHH-MM) or an existing dataset directory / .csv / .feather file.
```

(The same check runs on `--compare_barcodes`; its message currently says `--barcodes` too, so look at which of the two options you typed.)

**Cause.** `--barcodes <date>` only selects a dataset **already stored** in `--db_dir`; wasteline cannot download old versions
(`freyja update` only gives the current data). Often the run is also looking in a different `--db_dir` than the one you
filled (the default is `~/.wasteline/db`; the test profiles use `<outdir>/db`).
**Fix.** List what is stored with `ls <db_dir>/freyja`, pass the same `--db_dir` you used when downloading, or fetch the current
data with `--mode update_db`. Use `--barcodes container` or a path for fully offline runs. See
[usage.md](usage.md#6-freyja-data-barcodes).

### `--barcodes latest` fails without network

**Symptom.** The `FREYJA_UPDATE (latest)` task fails with a Python traceback ending

```
urllib.error.URLError: <urlopen error [Errno -3] Temporary failure in name resolution>
```

(reproduced by running the container with `--network none`). **Cause.** `freyja update` downloads the dataset and the machine
has no internet access to it (typical for HPC compute nodes). **Fix.** Use `--barcodes container`, or download on a connected
machine with `--mode update_db --db_dir <shared dir>` and pin with `--barcodes <date>`.

### The Freyja database directory is not writable

**Symptom** (`--db_dir` or `~/.wasteline/db` is read-only or its parent cannot be created):

```
ERROR ~ Unexpected error while finalizing task 'CHARLESFOSTER_WASTELINE:WASTELINE:FREYJA_DATA:FREYJA_UPDATE (latest)' - cause: Failed to create publish directory: /path/to/db/freyja
```

**Cause.** `--barcodes latest` (and `--mode update_db`) store the downloaded dataset in `<db_dir>/freyja/<stamp>/`, which must be
writable by the user running Nextflow. **Fix.** Point `--db_dir` at a directory you can write to. Pinned runs
(`--barcodes <date>`, a path, `container`) only read from it.

### `--barcodes latest` did not give me the newest data on a resumed run

**Cause.** With `-resume`, the cached `FREYJA_UPDATE` task is re-used, so the earlier download is used again (verified:
`cached=2`). **Fix.** Run without `-resume`, or run `--mode update_db` first and pin the result with `--barcodes <date>`.

### Results from two runs differ even with the same reads

The barcode dataset changes almost daily, and a new dataset changes lineage abundances. Compare `freyja_data.barcode_version` in
the two `pipeline_info/run_info.json` files. Also note that the old Snakemake pipeline used `--varthresh 0.25` while wasteline
defaults to `0.0` ([migration_from_snakemake.md](migration_from_snakemake.md)).

## Per-sample results

### A sample is `low_confidence`

**Symptom.** `summary/demix_status.tsv` has, for the sample, `low_confidence   45.0% of the genome covered at >= 10x: below the
high-confidence threshold of 50%; treat the estimate with care`. It has a `.demix` file and a row in `freyja_aggregated.tsv`,
`lineage_summary.csv` and the report (amber `low confidence` badge, `demix_low_confidence` flag), but it is missing from
`freyja_timeseries.tsv`, the growth rates, the dashboard and the site aggregates.

**Cause.** Its genome coverage (% of positions at depth >= `--freyja_covcut`) lies between `--demix_floor_coverage` (default 25)
and `--min_genome_coverage_for_demix` (default 50), so it was demixed but is not trusted as much as an `ok` sample (typically
amplicon dropout: check `amplicons_dropped` in `qc_summary.tsv`). **Fix.** Nothing is wrong with the run. To trend these samples
anyway, set `--include_low_confidence true`; to treat them as skipped (the behaviour before D18), set `--demix_floor_coverage 0`
or to the value of `--min_genome_coverage_for_demix`.

### A sample is `skipped_low_coverage`

**Symptom.** `summary/demix_status.tsv` has, for the sample,
`skipped_low_coverage   less than 50% of the genome covered at >= 10x`; it has no `.demix` file and is absent from
`freyja_aggregated.tsv` and the plots (and `qc_summary.tsv` shows `low_coverage` or `no_reads`). Negative controls with almost
no reads are skipped by design and are fine.

**Cause.** The coverage gate: fewer than `--demix_floor_coverage` (default 25) percent of the genome has depth >=
`--freyja_covcut` (default 10) (samples from there up to `--min_genome_coverage_for_demix`, default 50, are `low_confidence`, see above). Freyja would otherwise return a confident-looking lineage from far too little data.
**Fix.** Look at `mean_depth` and `pct_genome_ge_10x` in `qc_summary.tsv` and the coverage track in the report. Usually the
library needs resequencing, or amplicons dropped out (`amplicons_dropped`). If a lower coverage is acceptable for you, lower
`--demix_floor_coverage` (and/or `--freyja_covcut`), keeping in mind that low-coverage abundances are unreliable;
setting it to 0 turns the gate off, and then a sample with 0% coverage is reported as `ok`.

### A sample is `failed` with a pandas traceback (`freyja demix` crash on near-empty samples)

**Symptom.** `demix_status.tsv`:

```
NEG_E   failed   freyja demix failed: ... TypeError: operation 'radd' not supported for dtype 'str' with dtype 'object'
```

**Cause.** A sample with almost no mapped reads makes `freyja variants` write a header-only table, and `freyja demix` crashes on it.
The coverage gate normally routes such samples to `skipped_low_coverage` before `freyja demix` runs, so you only see this when
the gate has been disabled or lowered to 0. The failure is recorded (status `failed`, reason from Freyja's own output) and the
run **continues**; it never aborts. **Fix.** Restore the gate. For any other `failed` demix read the reason in
`demix_status.tsv`/`<sample>.demix_status.json`, and the full Freyja output in `<work dir>/<sample>.demix.log` of the
`FREYJA_DEMIX (<sample>)` task.

### `freyja_plot.failed.txt` is in `summary/` instead of `freyja_plot.pdf`

**Symptom.** `summary/freyja_plot.failed.txt` (and/or `freyja_plot_times.failed.txt`) exists and there is no PDF. The run still
succeeds and the report still shows abundances. Example content:

```
freyja plot produced no PDF (plain attempt and infer_string workaround both failed).
--- attempt 1 (freyja plot), last lines ---
ERROR: No samples matching coverage requirements, so no plot will be generated. Try changing --mincov threshold.
```

**Cause.** Two possible reasons, both visible in the file. (a) No sample has genome coverage above `--freyja_plot_mincov`
(by default the demix threshold, 50; or an explicit value you passed), as above: lower it, for example `--freyja_plot_mincov 30` (this is the case that was reproduced). (b) Freyja's own
`freyja plot` is broken in Freyja 2.0.3 to 2.0.5 containers (a pandas 3 `TypeError`, exit status 0 with no PDF). The pipeline
retries with a documented workaround, which succeeds in normal runs (the real test run needed it); the note is written only if
that also fails. **Fix.** Use the report's own plots, which do not depend on `freyja plot`; for (a) lower
`--freyja_plot_mincov`; for (b) report the traceback from the note upstream and update the pipeline when its Freyja pin moves
(see [usage.md](usage.md#16-updating-the-freyja-version-of-the-pipeline)).

### A sample has no bootstrap interval or no barcode comparison

With `--freyja_boot` a sample only has rows in `summary/freyja_boot_summary.tsv` if it was demixed and its bootstrap worked: a
skipped sample has none, and a failed bootstrap is written to `samples/<s>/freyja/<s>.boot_status.json` (status `failed`, reason
from Freyja) without stopping the run. With `--compare_barcodes`, `summary/barcode_compare_status.tsv` has one row per sample
with the outcome of the **second** demix (`ok`, `skipped_low_coverage`, `failed` + reason); only samples that are `ok` in both
demixes appear in `barcode_compare.tsv`.
<!-- UNVERIFIED: a failing bootstrap or second demix was not provoked; the behaviour is from reading modules/local/freyja_boot and the workflow. -->

### A sample is missing from the aggregate or the plot

In order of likelihood: it is a **negative control** (excluded by default, listed in `summary/freyja_aggregated_controls.tsv` instead; add
`--include_negatives_in_aggregate` to put it in the aggregate); its demix was
`skipped_low_coverage` or `failed` (see `demix_status.tsv`); in `aggregate` mode its `.demix` file was missing or zero-byte
(status `missing`, reason `missing demix (new layout)` or `empty demix file (failed in the old pipeline)`); or Freyja's plot
drops samples below `--freyja_plot_mincov` (by default the demix threshold, or the low-confidence floor with `--include_low_confidence true`; the table `freyja_aggregated.tsv` keeps them; the PDF does not).

### The consensus is mostly `N`, or `consensus_stats.tsv` says `empty` or `failed`

Positions with depth below `--consensus_min_depth` (default 10) are written as `N`, so a sample with patchy coverage has a high
a low `consensus_coverage_pct`; a sample with no usable depth is `empty` (all `N`, reference length); `failed` means a consensus step failed for that
sample (the run continues, and the row is still present; `reason`, the last column of `consensus_stats.tsv`, says why). A sample
whose variants could not all be applied is `ok` with the skipped variants listed in `reason` and in `consensus_status.json`
(`warnings`); the report shows an INFO flag. `failed` raises the WARN flag `consensus_failed`. Check `amplicons_dropped` and the coverage track, and lower
`--consensus_min_depth` only if you accept less supported bases. Negative controls get no consensus.

### QC says `WARN` (`amplicon_dropout`) for almost everything

`--qc_max_dropped_amplicons` defaults to 1, which is strict for 96/99-amplicon ARTIC schemes where a few dropped amplicons are
routine in wastewater. Raise it (the schema suggests 10) and/or adjust `--amplicon_min_depth` (default 10). Amplicon
dropout is only computed when the scheme is known: in `--mode demix` pass `--primer_scheme`/`--primer_bed`.

### Sample IDs do not match my old results

wasteline never rewrites sample IDs; the old Snakemake pipeline turned `_` into `-` in the aggregate. Normalise the old side when
joining. See [migration_from_snakemake.md](migration_from_snakemake.md).

## Report and modes

### `--mode report` stops with an error

```
[ERROR] --results_dir 'results' has no pipeline_info/run_info.json: report mode needs the provenance record of the run (it never invents one)
[ERROR] --outdir must differ from --results_dir ('results'): report mode never writes into the results it reads
[ERROR] No --results_dir given. --mode report needs --results_dir <results directory of a previous run> (it takes no samplesheet)
```

**Fix.** Point `--results_dir` at the directory of a finished run that has `pipeline_info/run_info.json` (runs in `full`, `demix`
and `aggregate` write it; `report` and `update_db` runs do not, so you cannot re-report a report), and give a **different**
`--outdir`. The new `summary/report.html` is the only file written.

### The report flags `demix_missing` for samples that were not demixed

`No demix result was found for 2 sample(s): LOWCOV_D, NEG_E.` Those samples are listed as `missing` in `demix_status.tsv`
(`aggregate` mode: no `.demix` file in the earlier results, or a zero-byte legacy file). It is a warning about the inputs, not
an error; the other samples are processed normally. If a name you do not recognise is listed, see the next entry.

### A sample called `db` is listed as `missing` after `--previous_results`

**Symptom.** `demix_status.tsv` has `db  missing  missing variants, depths (legacy layout)` and the report warns
`demix_missing ... db`. **Cause.** The previous results directory contains a `db/` directory, which the results discovery
reads as a sample directory of the old Snakemake layout. The `test` profiles set `--db_dir` to `<outdir>/db`, so a run there that
stores or resolves a dataset (an empty `db/freyja` was left by `--compare_barcodes container`; `--barcodes latest` stores real
datasets) leaves one, as does any run where you pointed `--db_dir` inside `--outdir`. **Fix.** The row is harmless (nothing is demixed
for it). To avoid it, keep `--db_dir` outside `--outdir`, or list the samples yourself with `make_samplesheet.py --from-results`
and delete the `db` row. (Reproduced with `--mode demix --previous_results` on a `-profile test` output made with
`--compare_barcodes container`.)

### The anonymised report still shows a real name, or the id map is missing

`--report_anonymize` changes `report_data.json` and `report.html` only; the other files in `summary/` and the per-sample directories keep
the real IDs, so share the report file, not the results folder. `summary/sample_id_map.csv` is written only with `full` (never with `paths`).
A name that is not a sample ID or a site of this run (free text in a flag, a custom file name) is not replaced by `full`; `paths` handles
paths, the user name and e-mail addresses. To anonymise a finished run, use `--mode report --results_dir <run> --outdir <new>
--report_anonymize full`: the report data are rebuilt from the summary files even when `summary/report_data.json` exists. See
[usage.md](usage.md#sharing-a-report).

### The report looks empty in places

Sections are only drawn for what the run produced: no QC and consensus panels in `aggregate` mode; no amplicon panels in `demix`
mode without a primer scheme; no mutation panel without `--mutations_of_interest`; no barcode-comparison panel without
`--compare_barcodes` (and no bootstrap intervals in the tooltips without `--freyja_boot`); no time series without at least two
distinct `collection_date`s (the samplesheet supplies them in all of `full`, `demix` and `aggregate`); the **Negative controls**
section only appears when a control was demixed. The report is a single offline HTML file (no network is used to open it).

## Resuming and harmless warnings

### `-resume` re-ran everything

`-resume` continues from the **most recent run in the launch directory**. If you ran something else in between, find the run
you want with `nextflow log -f name,status,session` and use `-resume <session-id>` (verified). The cache is keyed on the inputs of each task
(file contents, parameters, container), so changing a parameter re-runs the steps that use it, and deleting `work/` removes the cache.
The provenance and report steps always re-run (they contain the run time); 3 of 112 tasks re-running is normal.

### Warnings you can ignore

| Message                                                                                                                 | Why it is harmless                                                                                                      |
| ----------------------------------------------------------------------------------------------------------------------- | ----------------------------------------------------------------------------------------------------------------------- |
| `WARN: nf-core pipelines do not accept positional arguments. The positional argument 'true' has been detected.`         | printed with `-stub-run`/`-stub` on Nextflow 26.04 (the boolean after the flag is seen as an argument). Stub runs work. |
| `Nextflow 26.10.0 is available - Please consider updating`                                                              | informational. Stay on >= 26.04.                                                                                        |
| `Negative control NEG_E is not included in the aggregate and plot (use --include_negatives_in_aggregate to include it)` | the default (decision D8); not an error.                                                                                |
| `WARNING: 5 sample(s) have no metadata row: ...` (from `make_samplesheet.py`)                                           | those samples simply get no date or site.                                                                               |

### The run ended with `Pipeline completed with errors`

A task failed in a way that is not recorded as a per-sample status (for example Docker could not start, a process ran out of
memory, a missing file). Nextflow finishes the tasks in flight, then stops. Read the first `[ERROR]` block: it gives the
process, exit status, the command, the `Work dir` and the container. Fix the cause and re-run the same command with
`-resume`. Out-of-memory exits (137) are retried once automatically with more memory; if they still fail, raise the process memory in
a `-c` config ([usage.md](usage.md#13-resources)). Per-sample QC, consensus and demix problems do **not** end a run like this:
they become rows with a status and a reason.
