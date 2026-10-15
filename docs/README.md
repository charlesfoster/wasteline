# charlesfoster/wasteline: Documentation

Start here. New to the pipeline? Read the [quick start](QUICKSTART.md) (five minutes, synthetic data), then the
[usage guide](usage.md).

## For users

| Page | What is in it |
|---|---|
| [QUICKSTART.md](QUICKSTART.md) | A copy-paste test on the bundled synthetic data: stub run, real run, expected results, the optional extras, the other modes. |
| [usage.md](usage.md) | The user guide: prerequisites, samplesheets for every mode, primer schemes, running `full`/`demix`/`aggregate`/`report`/`update_db`, Freyja data (barcode) management, comparing barcode versions, bootstrap intervals, consensus and variant annotation, QC settings, negative controls, mutations of interest, time series, `-resume`, profiles, resources, parameter files, HPC, updating Freyja. |
| [parameters.md](parameters.md) | Every pipeline parameter with its default, grouped as in `nextflow_schema.json`. |
| [output.md](output.md) | Every output file and column, `run_info.json` and how to reproduce a run. |
| [report.md](report.md) | Reading the self-contained HTML report (`summary/report.html`). |
| [example_report.html](example_report.html) | A real report from the synthetic test data (open it in a browser; paths anonymised). |
| [troubleshooting.md](troubleshooting.md) | Failure modes seen in practice: symptom, cause, fix. |
| [migration_from_snakemake.md](migration_from_snakemake.md) | For users of `covid-illumina-snakemake` (`CIS --workflow wastewater`): option mapping, behavioural differences, re-demixing old results. |
| [helper_scripts.md](helper_scripts.md) | `bin/make_samplesheet.py` and `bin/check_samplesheet.py` in detail. |

## For developers and maintainers

| Page | What is in it |
|---|---|
| [development.md](development.md) | Repo layout, dev environment, test layers, how to add things, Freyja bumps, conventions. |
| [CONTRIBUTING.md](CONTRIBUTING.md) | Contribution guidelines. |
| [test_data.md](test_data.md) | The deterministic synthetic dataset in `tests/data` and how it was validated. |
| [qc_scripts.md](qc_scripts.md) | The QC, consensus and mutation-of-interest helper scripts (`amplicon_coverage.py`, `qc_summary.py`, ...). |
| [report_data.md](report_data.md) | The `report_data.json` schema that feeds the HTML report. |
| [ci.md](ci.md) | GitHub Actions, the Freyja auto-update bot, branch protection, the digest lock; every check runnable locally. |
| [PLAN.md](PLAN.md) | The design and decision record (some details are superseded). |
| [PLAN_STATUS.md](PLAN_STATUS.md) | What of the plan is implemented, partial or missing. |
| [wiring/](wiring/) | **Design-era notes**: hand-over notes for the module authors written before the integration round ([phase3_modules.md](wiring/phase3_modules.md), [report_modules.md](wiring/report_modules.md)). Historical: where they differ from the code, the code is right. |
| [images/](images/) | Report screenshots used by the top-level [README.md](../README.md) (`report_overview.png`, `report_lineage_coverage.png`, `report_controls_comparison.png`). |

Other files at the top of the repository: [README.md](../README.md), [CHANGELOG.md](../CHANGELOG.md),
[CITATIONS.md](../CITATIONS.md) (tools to cite), [LICENSE](../LICENSE) (MIT) and [CLAUDE.md](../CLAUDE.md)
(conventions for the AI coding agents that helped write the code).
