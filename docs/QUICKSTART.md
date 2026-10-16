# Quick start: test the pipeline on the synthetic dataset

This takes about five minutes on a laptop (the first run also pulls the container images) and needs only
**Nextflow >= 26.04** and **Docker** (running). No real data and no internet access to Freyja's data servers are needed:
the commands below use the small synthetic dataset and barcode subset that ship in `tests/data/`.

> Run these commands from the **root of the repository checkout**. The sample sheet `tests/data/samplesheet.csv`
> uses paths relative to the repository root, and Nextflow resolves relative sample-sheet paths against the directory
> you launch from (see [troubleshooting](troubleshooting.md#samplesheet-the-file-or-directory-does-not-exist)).

```bash
git clone https://github.com/charlesfoster/wasteline.git && cd wasteline   # or cd into your existing checkout
# 1. Wiring check (stub run: no real analysis, placeholder outputs, ~1 min once images are cached; 113 tasks)
nextflow run . -profile docker,test_stub -stub-run --outdir results_stub
# 2. Real run on the synthetic reads (~1-2 min with the images cached; 112 tasks)
nextflow run . -profile docker,test --primer_scheme midnight --barcodes tests/data/barcodes_test.csv --outdir results
# (the synthetic reads are Midnight amplicons; the pipeline default scheme is qiagenv1, so the scheme is given explicitly)
# 3. Look at the results
ls results/summary                       # tables, plots, consensus genomes, report.html (13 files)
open results/summary/report.html         # macOS; use xdg-open on Linux, or double-click the file
cat results/summary/demix_status.tsv     # one row per sample: ok / skipped_low_coverage / failed, with the reason
column -t -s$'\t' results/summary/qc_summary.tsv | cut -c1-150   # PASS/WARN/FAIL per sample
# 4. (optional) check the results against the known truth; needs Python >= 3.12 (the pixi env has it)
python3 tests/check_outputs.py results --mode full --barcodes tests/data/barcodes_test.csv
```

(The `git clone` line only works once the repository is published at that address
<!-- UNVERIFIED: github.com/charlesfoster/wasteline was not reachable (404) when this page was written; every command below was run from a local checkout -->.)

Both runs were repeated for this page (Nextflow 26.04.6, Docker; stub: `completed=113 failed=0`, real run: `completed=112
failed=0`, about one minute each with the images already pulled), and `tests/check_outputs.py` reported `128/128 checks passed`.

The stub run (step 1) is a syntax and wiring test: every process writes placeholder files, so its `freyja_aggregated.tsv`
and report are not meaningful. It prints a harmless `WARN: nf-core pipelines do not accept positional arguments. The
positional argument 'true' has been detected` line: Nextflow 26.04 passes the boolean of `-stub-run` through as an
argument. Ignore it.

## What you should see

The dataset (`tests/data/`, described in [test_data.md](test_data.md)) has seven samples. The truth is in
`tests/data/truth.json`. With `--barcodes tests/data/barcodes_test.csv` (a 14-lineage subset of Freyja's barcodes) a real
run gives, within the tolerance of 0.05:

| Sample          | Truth                                   | Freyja result (from a real run)                                     | `demix_status`         | `qc_status` (reasons)                                    |
| --------------- | --------------------------------------- | ------------------------------------------------------------------- | ---------------------- | -------------------------------------------------------- |
| `MIX_A`         | BA.2 0.70, AY.4 0.30                    | BA.2 0.71, AY.4 0.29                                                | `ok`                   | PASS                                                     |
| `MIX_B`         | BA.1 0.50, XBB.1.5 0.30, KP.3 0.20      | BA.1 0.50, XBB.1.5 0.28, KP.3 0.20, BA.2 0.02 (small spurious call) | `ok`                   | PASS                                                     |
| `PURE_C`        | JN.1 1.00                               | JN.1 1.00                                                           | `ok`                   | PASS                                                     |
| `Site_1_week_2` | BA.5 0.60, XBB.1.5 0.40                 | BA.5 0.61, XBB.1.5 0.36, BA.2 0.03 (small spurious call)            | `ok`                   | WARN (`amplicon_dropout`: 2 amplicons dropped by design) |
| `LOWCOV_D`      | BA.2 1.00 at about 0.6x depth           | no result                                                           | `skipped_low_coverage` | FAIL (`low_coverage;amplicon_dropout`)                   |
| `NEG_E`         | negative control, 3 read pairs          | no result                                                           | `skipped_low_coverage` | PASS                                                     |
| `NEG_CONTAM`    | negative control contaminated with BA.2 | demixed, but kept out of the aggregate and the plots                | `ok`                   | FAIL (`ntc_reads`)                                       |

Things to notice:

- `summary/freyja_aggregated.tsv` has only **four** rows (`MIX_A`, `MIX_B`, `PURE_C`, `Site_1_week_2`): skipped samples
  have no `.demix` file, and negative controls are demixed but excluded from the aggregate unless you pass
  `--include_negatives_in_aggregate`. What the demixed control contains is in `summary/freyja_aggregated_controls.tsv`
  (one row: `NEG_CONTAM`, `BA.2` 1.0).
- `Site_1_week_2` keeps its underscores: sample IDs are never rewritten.
- `LOWCOV_D` would be demixed by Freyja to a confident single lineage if the coverage gate did not stop it; the gate
  (`--min_genome_coverage_for_demix`, default 60%) is why it is `skipped_low_coverage` with a reason, not a result.
- `NEG_CONTAM` raises a run-level error flag (`ntc_reads`, "5998 mapped reads, limit 1000") in the report: that is the
  negative-control check doing its job. Its demix status stays `ok`, and it has its own **Negative controls** section in the report.
- `qc_summary.tsv` is complete: besides the coverage columns it carries `demix_status`, `demix_reason`, `consensus_length`,
  `consensus_n_pct` and `barcode_version` for every sample.
- With the default `-profile test` barcodes (`--barcodes container`, the full Freyja dataset bundled in the pinned container)
  the mixtures still demix, but Freyja splits lineages into many close relatives (for example `MIX_A` becomes BA.2 0.69,
  AY.4 0.18, AY.1 0.07, ...). The barcode subset gives the clean table above, which is why the command passes it.

Where to look:

| File                                                                                                                  | What it is                                                                                |
| --------------------------------------------------------------------------------------------------------------------- | ----------------------------------------------------------------------------------------- |
| `results/summary/report.html`                                                                                         | the self-contained report (works offline; open it in any browser)                         |
| `results/summary/freyja_aggregated.tsv`, `freyja_aggregated_controls.tsv`, `freyja_plot.pdf`, `freyja_plot_times.pdf` | Freyja's aggregate table (samples), the table of demixed negative controls, and the plots |
| `results/summary/qc_summary.tsv`, `demix_status.tsv`, `amplicon_coverage.tsv`                                         | QC, per-sample demix outcome, per-amplicon coverage                                       |
| `results/summary/consensus_genomes.fasta`, `consensus_stats.tsv`                                                      | consensus genomes (none for negatives) and their N content                                |
| `results/summary/barcode_version.json`, `results/pipeline_info/run_info.json`                                         | which Freyja barcode version and which software produced the results                      |
| `results/samples/<sample>/`                                                                                           | per-sample files: fastp, BAM, Freyja variants/depths/demix, consensus                     |

Next: [usage.md](usage.md) to run it on your own data, [output.md](output.md) for every output file.

## Trying the optional extras

The extras are opt-in. On the synthetic data (the second command re-uses the cached tasks of step 2, so only the new
processes run: 35 completed, 109 cached, under a minute):

```bash
cat > moi.csv <<'EOF'
name,mutation
S N501Y,S:N501Y
S L452R,S:L452R
EOF
nextflow run . -profile docker,test --primer_scheme midnight --barcodes tests/data/barcodes_test.csv \
    --compare_barcodes container --freyja_boot --annotate_variants --mutations_of_interest moi.csv \
    --outdir results_extras -resume
ls results_extras/summary
```

New files in `results_extras/summary/`:

| File                                                | From                           | What it shows                                                                                                                                                                                               |
| --------------------------------------------------- | ------------------------------ | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `barcode_compare.tsv`, `barcode_compare_status.tsv` | `--compare_barcodes container` | each sample demixed again with the Freyja container's dataset: e.g. `MIX_A` `BA.2` 0.712 (test subset) vs 0.692 (container dataset), `AY.4` 0.288 vs 0.176, plus `AY.1` 0.071 only in the container dataset |
| `freyja_boot_summary.tsv`                           | `--freyja_boot`                | 95% bootstrap intervals: `MIX_A` `BA.2` 0.712, 0.670 to 0.727                                                                                                                                               |
| `annotated_variants.tsv`                            | `--annotate_variants`          | amino-acid consequences of the applied variants (e.g. `MIX_A` `ORF1ab` `S135R`)                                                                                                                             |
| `mutations_of_interest.tsv`                         | `--mutations_of_interest`      | `MIX_A` `S:N501Y` alt frequency 0.788 (41 of 52 reads), `PURE_C` 1.0                                                                                                                                        |

The report gains a **Barcode comparison** panel and the interval tooltips; each option is explained in
[usage.md](usage.md#comparing-barcode-versions---compare_barcodes).
(The numbers are real outputs of this run. Note: the `test` profile stores `--db_dir` in
`<outdir>/db`, so `results_extras/` contains an empty `db/freyja` directory. Do not use that directory with
`--previous_results` unless you accept an extra `db ... missing` row, see [troubleshooting](troubleshooting.md#a-sample-called-db-is-listed-as-missing-after---previous_results).)

## Trying the other modes on these results

All of these reuse `results/` read-only (about a minute each):

```bash
# re-demix the existing variants/depths with other barcodes (here the same test barcodes)
nextflow run . -profile docker --mode demix --previous_results results --barcodes tests/data/barcodes_test.csv --outdir results_demix
# re-aggregate existing .demix files (negatives and unusable samples are listed in demix_status.tsv)
nextflow run . -profile docker --mode aggregate --previous_results results --barcodes tests/data/barcodes_test.csv --outdir results_aggregate
# re-render the HTML report only (needs results/pipeline_info/run_info.json; --outdir must differ from --results_dir)
nextflow run . -profile docker --mode report --results_dir results --outdir results_report
```

A repeat of the identical command with `-resume` re-uses the cached tasks (in this test: 109 cached, 3 re-run, the
provenance and report steps, which embed the run time), provided it is the most recent run in the launch directory; otherwise
give `-resume <session-id>` (see [usage.md](usage.md#11-resuming-a-run--resume)). The three modes above were re-run for this
page: `demix` (`completed=49`) reproduces the `.demix` files, `freyja_aggregated.tsv`, `demix_status.tsv` and
`consensus_stats.tsv` byte for byte, `aggregate` (`completed=7`) and `report` (`completed=1`) succeed.
