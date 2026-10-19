# Re-analysing past runs with new Freyja barcodes

Freyja's barcode dataset (the lineage definitions) is updated almost daily, and the Freyja release itself moves on too. A
lineage call from three weeks ago was made with the barcodes of three weeks ago. `bin/reanalyse_runs.py` re-estimates the
lineages of recent runs with a newer dataset **without mapping the reads again**, and shows what changed.

It is a thin driver around two things that already exist:

- `--mode demix --previous_results <run>` ([usage.md](usage.md#demix-re-estimate-lineages-from-existing-variants-and-depths)),
  which re-runs `freyja demix`, the aggregate, the QC, the lineage groups and the report from the `variants` and `depths` files
  a run kept, and writes everything to a **new** directory;
- `--mode update_db` ([usage.md](usage.md#update_db-download-the-current-freyja-data-and-stop)), used once to fetch the
  current dataset when you ask for `--barcodes latest`.

It uses the Python standard library only and runs on your host Python (3.9 or newer), not in a container. It never touches the
network itself and never modifies or overwrites existing results.

## Quick start

1. **Dry run** (the default). Lists the runs it found and prints the exact Nextflow command for each. Nothing is written.

   ```bash
   python bin/reanalyse_runs.py /data/wasteline/results --last 4 --barcodes latest
   ```

   ```
   Target barcodes: latest (resolved when run)   (abridged)

     run         date        mode  barcodes          freyja  samples  action
     run1        2026-10-16  full  2026-10-06T05-25  2.0.6   10       reanalyse
     run2        2026-10-16  full  2026-10-06T05-25  2.0.6   9        reanalyse
     2026-10-19  2026-10-19  full  2026-10-06T05-25  2.0.6   12       reanalyse

   [1/3] run1  ->  /data/wasteline/results/run1_reanalysis_latest
       params file /data/wasteline/results/run1_reanalysis_latest.params.json (written with --execute)
       carried over: min_genome_coverage_for_demix=40, primer_scheme=qiagenv1

   nextflow run /home/me/wasteline \
       -profile docker \
       -params-file /data/wasteline/results/run1_reanalysis_latest.params.json \
       --mode demix \
       --previous_results /data/wasteline/results/run1 \
       --outdir /data/wasteline/results/run1_reanalysis_latest \
       --barcodes latest
   ```

2. **Execute** when the commands look right. The runs go one after the other; each command is printed first.

   ```bash
   python bin/reanalyse_runs.py /data/wasteline/results --last 4 --barcodes latest --execute \
       --report /data/wasteline/reanalysis_diff.csv --summary-md /data/wasteline/reanalysis_diff.md
   ```

   After each run the old and new results are compared and the comparison is written to
   `<new results>/summary/reanalysis_diff.csv` (and `.md`). `--report` and `--summary-md` additionally collect all runs in one
   file (the CSV gets a leading `run` column).

3. **Look at the differences**: [how to read the diff](#how-to-read-the-diff).

`<root>` can be a results directory itself or a folder that contains results directories (found by
`pipeline_info/run_info.json`, up to two levels down). Several roots can be given. Legacy Snakemake results
(`<sample>/freyja/...`, [migration_from_snakemake.md](migration_from_snakemake.md)) are found too and are listed with barcode
version `unknown`.

## Choosing the runs

| Option               | Meaning                                                                                                             |
| -------------------- | ------------------------------------------------------------------------------------------------------------------- |
| `--since YYYY-MM-DD` | Runs that started on that date or later (the `start` of `run_info.json`; for legacy results the directory's mtime). |
| `--last N`           | The N most recent runs that can be re-analysed.                                                                     |
| `--runs run1,run2`   | Named runs (directory name or path).                                                                                |
| none of them         | Every run found.                                                                                                    |

They combine (names first, then `--since`, then `--last`). Output is always sorted by date, then path.

A run is listed but **skipped** (with the reason in the table) when it cannot be re-analysed:

- its mode is `demix`, `aggregate`, `report` or `update_db`: those directories keep no variants/depths files. This includes
  the `..._reanalysis_...` directories this script makes. Re-analyse the original run instead;
- the recorded run did not succeed (`success: false`), or no sample has both a variants and a depths file (for instance after
  clearing `samples/*/freyja/` to save space);
- the new directory already exists (see [never overwritten](#original-results-are-never-touched)).

## The target barcodes

`--barcodes` takes the same values as the pipeline (default `latest`):

| Value                                    | What happens                                                                                                                                                                                                                                                                                                                                               |
| ---------------------------------------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `latest`                                 | With `--execute` the script first runs `nextflow run ... --mode update_db` **once**, reads the stored dataset's version from `pipeline_info/freyja_data.json` and uses that exact stamp for every run, so all runs share one dataset. Needs the network. In a dry run `latest` is not resolved: nothing is downloaded and no run is skipped as up to date. |
| `YYYY-MM-DD`                             | The newest stored dataset of that day in `--db_dir` (default `~/.wasteline/db`).                                                                                                                                                                                                                                                                           |
| `MM_DD_YYYY-HH-MM` or `YYYY-MM-DDTHH-MM` | One exact stored dataset.                                                                                                                                                                                                                                                                                                                                  |
| `container`                              | The data bundled in the Freyja container of the pipeline checkout. Offline. Never skipped as up to date (the checkout decides what that is).                                                                                                                                                                                                               |
| a directory, `.csv` or `.feather`        | A dataset used verbatim (compared with the recorded run by sha256 of the barcode file).                                                                                                                                                                                                                                                                    |

A run whose recorded barcode version already equals the target is printed as `up to date` and not run again (for a date: the
newest stored stamp of that day in `--db_dir` equals the recorded one, or, when no dataset of that day is stored, the recorded stamp
falls on that day). Runs with `custom` or `unknown` barcodes are
never up to date.

The new directory is named after the resolved version. With a dry run of `latest` the name shows `latest`; with `--execute` it
is the stamp, for example `run1_reanalysis_2026-10-20T04-10`.

## What carries over

The command is built from the run's own `run_info.json`:

| Part                                                     | Source                                                                                                                                     |
| -------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------ |
| pipeline                                                 | `--pipeline-dir` (default: the clone this script lives in). The **current** code and Freyja container, not the ones of the old run.        |
| `-profile`                                               | The recorded profile without `locked` and the `test*` profiles (default `docker`); `--profile` overrides.                                  |
| `--mode`, `--previous_results`, `--outdir`, `--barcodes` | `demix`, the old run's directory, the new directory, the target.                                                                           |
| `--db_dir`                                               | Only when you give it.                                                                                                                     |
| `-params-file`                                           | `<new dir>.params.json`, next to the new directory (written with `--execute`). Holds the recorded values that matter in demix mode, below. |
| extra Nextflow args                                      | `--extra-args '-work-dir /scratch/work -resume'` (one quoted string, split like a shell would).                                            |

Carried into the params file (only when the value is not the default of the pipeline in `--pipeline-dir`):
`min_genome_coverage_for_demix`, `freyja_*`, `consensus_*`, `qc_*`, `ntc_*`, `amplicon_*`, `unexplained_*`, `dropout_*`,
`primer_check_*`, `lineage_*`, `site_*`, `timeseries*`, `mutations_of_interest`, `annotate_variants`,
`include_negatives_in_aggregate`, `skip_consensus`, `skip_report`, `protocol`, `skip_primer_trimming`, and a non-bundled
`fasta`/`gff`. Parameters the pipeline in `--pipeline-dir` does not have are left out with a note, so an older or newer checkout
does not make nf-schema reject the file.

Specific rules:

- **Primer scheme.** In `demix` mode amplicon coverage runs only when a scheme or BED is given explicitly, so the run's
  effective `primer_scheme` is always passed on (even when it is the default). A custom `primer_bed` that still exists replaces
  the scheme. Both are left out for `--protocol capture` and `--skip_primer_trimming` runs (`protocol` itself is carried).
- **Lineage groups.** A custom `--lineage_groups` file is used while it exists. Otherwise the copy published in the results
  (`summary/lineage_groups_used.csv`, byte-identical to the file used; checked against the recorded sha256) is used, with a note.
  If there is neither, the bundled default of the current pipeline applies, with a note. A run that used the bundled groups gets
  the bundled groups of the checkout (which may have changed since).
- **Paths** (`mutations_of_interest`, `primer_bed`, ...) are made absolute when they still exist, and left out with a note when they
  do not.
- **Left out on purpose**: `input`, `outdir`, `barcodes`, `compare_barcodes`, `db_dir`, mapping and read-handling parameters
  (`fastp_*`, `duplicates`, ...), `email` and the other run-specific values.
- **Legacy Snakemake results** record nothing: the pipeline defaults apply and the table says so. The old pipeline discarded
  variants below 25% frequency when it made the variants file, and that filtering stays in the files you re-analyse.

## Original results are never touched

- Every re-analysis goes to `<run>_reanalysis_<barcode version>` next to the run. `--outdir-suffix` changes the suffix
  (`{barcode}` stands for the version; the default is `_reanalysis_{barcode}`).
- If that directory exists the run is **skipped** (`output exists`), never overwritten. `--force-new-suffix` adds a counter
  (`_2`, `_3`, ...) instead.
- The old run is only read. The only files the script writes are the new directory's `summary/reanalysis_diff.*` and the params file beside it.
- Execution stops at the first failed run (exit status 1) unless `--keep-going`. Nextflow's `work/` and `.nextflow/` are created in the
  directory you launch from, like any other run (`--extra-args '-work-dir ...'` moves the work directory).

## How to read the diff

`--diff` does the comparison on its own for any two results directories, for example a run and its re-analysis, or two runs made with
`--mode demix` from the same variants:

```bash
python bin/reanalyse_runs.py --diff results/run1 results/run1_reanalysis_2026-10-20T04-10 \
    --threshold 0.05 --report run1_diff.csv --summary-md run1_diff.md
```

It reads `summary/freyja_aggregated.tsv` (specific lineages from `lineages`/`abundances`, and Freyja's `summarized` groups) and,
when **both** runs have `summary/lineage_summary.csv`, the lineage-group columns (the groups, and `other`). The list columns are
parsed tolerantly: Freyja's space-separated form, python lists, numpy reprs (also wrapped over several lines) and `np.float64(...)`
values. Sample labels such as `S1.variants.tsv` are cleaned like the pipeline does; sample IDs are otherwise never rewritten.

`--level lineages|summarized|groups|all` chooses what is compared (default `all`: every level present in both runs).

The CSV has one row per sample and item:

| Column                           | Meaning                                                                                                                                       |
| -------------------------------- | --------------------------------------------------------------------------------------------------------------------------------------------- |
| `sample`                         | Sample ID.                                                                                                                                    |
| `level`                          | `lineages` (a specific lineage), `summarized` (Freyja's summarized group) or `groups` (a column of `lineage_summary.csv`).                    |
| `item`                           | The lineage or group.                                                                                                                         |
| `old_abundance`, `new_abundance` | Abundance (0-1) in the old and new run; an item missing on one side is 0.                                                                     |
| `change`                         | `new - old`.                                                                                                                                  |
| `flagged`                        | `yes` when \|change\| >= `--threshold` (default 0.05; exactly the threshold counts, float noise is ignored).                                  |
| `note`                           | `sample only in the new run`, `sample missing in the new run`, and `status ok -> skipped_low_coverage` style changes from `demix_status.tsv`. |

Items below 0.001 in both runs are left out. Rows are sorted by sample, level, then the largest change first, so the file is stable
between runs. Filter on `flagged = yes` in a spreadsheet to see only what matters.

The console summary (and the Markdown file) gives: the barcode version and Freyja version **old -> new**, the number of samples
compared and how many have a flagged change, samples that appear or disappear and status changes, the largest changes, and notes
(for example that the group columns of the two runs differ because the groups file changed, or that both runs used the same
barcode and Freyja version so no difference is expected).

Reading it:

- A **sample that appears or disappears** (or changes status) is the first thing to check: a different Freyja version or barcode
  can make a borderline sample pass or fail the demix. Its other rows then show large changes by construction.
- At `level groups` a change that moves between two columns (`other` falling while a new group rises) usually means the groups file
  changed, not the biology: compare `summary/lineage_groups_used.csv` of both runs.
- A change in `lineages` that is absorbed in `summarized` or `groups` (a lineage renamed or split inside the same group) is a
  barcode-update effect, not a different virus.
- The comparison cannot tell a barcode effect from a Freyja version effect when both changed. The versions are printed so you can see
  which it was; `--compare_barcodes` in a normal run isolates the barcode effect.

## Suggested routine

After each Freyja barcode update, re-analyse the last four weeks and review the differences:

```bash
python bin/reanalyse_runs.py /data/wasteline/results \
    --since "$(date -d '28 days ago' +%F)" --barcodes latest --execute --keep-going \
    --report /data/wasteline/reanalysis/$(date +%F)_reanalysis_diff.csv \
    --summary-md /data/wasteline/reanalysis/$(date +%F)_reanalysis_diff.md
```

(On macOS `date -v-28d +%F`; or use `--last 4` when you make one run a week.) Then open `reanalysis_diff.csv`, filter `flagged = yes`,
and look at the samples whose dominant lineage moved, or that appeared or disappeared. Runs already re-analysed with that dataset are
not repeated: their new directory exists, so they are skipped as `output exists`.

Re-analysis directories are full results (report, consensus, lineage summary), so they can be archived or deleted when you have looked.
Keep the stored datasets (`--db_dir`): an old barcode version cannot be downloaded again, and the replay of an old run needs it
(`bin/replay_run.py`, see [usage.md](usage.md)).

### Scheduling (examples only, not installed by anything)

cron, every Monday at 06:00 (the dataset changes almost daily; once a week is usually enough):

```cron
0 6 * * 1  cd /data/wasteline && /usr/bin/python3 /home/me/wasteline/bin/reanalyse_runs.py /data/wasteline/results --since "$(date -d '28 days ago' +\%F)" --barcodes latest --execute --keep-going --report /data/wasteline/reanalysis/reanalysis_diff.csv >> /data/wasteline/reanalysis/cron.log 2>&1
```

launchd on macOS (`~/Library/LaunchAgents/org.wasteline.reanalyse.plist`; load with `launchctl load ~/Library/LaunchAgents/org.wasteline.reanalyse.plist`).
Docker Desktop must be running, and launchd gives jobs a minimal `PATH`, so name `nextflow` explicitly:

```xml
<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE plist PUBLIC "-//Apple//DTD PLIST 1.0//EN" "http://www.apple.com/DTDs/PropertyList-1.0.dtd">
<plist version="1.0">
<dict>
  <key>Label</key><string>org.wasteline.reanalyse</string>
  <key>ProgramArguments</key>
  <array>
    <string>/usr/bin/python3</string>
    <string>/Users/me/wasteline/bin/reanalyse_runs.py</string>
    <string>/Users/me/wasteline_results</string>
    <string>--last</string><string>4</string>
    <string>--barcodes</string><string>latest</string>
    <string>--execute</string>
    <string>--keep-going</string>
    <string>--nextflow</string><string>/usr/local/bin/nextflow</string>
    <string>--report</string><string>/Users/me/wasteline_results/reanalysis_diff.csv</string>
  </array>
  <key>WorkingDirectory</key><string>/Users/me/wasteline_work</string>
  <key>StartCalendarInterval</key><dict><key>Weekday</key><integer>1</integer><key>Hour</key><integer>6</integer></dict>
  <key>StandardOutPath</key><string>/Users/me/wasteline_work/reanalyse.log</string>
  <key>StandardErrorPath</key><string>/Users/me/wasteline_work/reanalyse.err</string>
</dict>
</plist>
```

## Exit status and limits

- Exit 0: the dry run printed, or every executed run succeeded (including "nothing to re-analyse"). Exit 1: unreadable input, no results
  found, `update_db` failed, or a Nextflow run failed.
- Re-analysis is sequential. The pipeline version and the Freyja container are those of `--pipeline-dir` (the re-analysis record in
  `pipeline_info/run_info.json` of the new directory says which, and lists the old run in `previous_runs`).
- The script never decides which of two datasets is "right"; it only shows how the estimates differ.
