# Excel output (`lineage_summary.xlsx`): wiring notes

Besides `lineage_summary.csv`, the lineage summary can be written as a formatted workbook that opens in Excel (and
Numbers / LibreOffice / Google Sheets) with no import step: numbers are real numbers, long lineage lists wrap, the
header row is frozen and filterable. The CSV and TSV outputs are unchanged and stay the machine-readable results.

| File | Role |
| --- | --- |
| `bin/xlsx_writer.py` | Minimal Office Open XML writer (standard library only: `zipfile` + hand-written XML) and a CSV/TSV converter CLI. Written to run on Python 3.8+. |
| `bin/lineage_groups.py` | New optional `--xlsx-out PATH`, `--xlsx-decimals N` (default 3) and `--wasteline-version TEXT`. Without `--xlsx-out` nothing changes (outputs are byte-identical, tested). |
| `modules/local/lineage_summary/` | New optional output `xlsx` (`lineage_summary.xlsx`), produced only when `ext.args` contains `--xlsx-out lineage_summary.xlsx`. No new container: the module's python container is enough. |
| `tests/python/test_xlsx_writer.py`, `tests/python/test_lineage_groups.py`, `tests/fixtures/xlsx/` | pytest suite (every XML part parsed, style ids, escaping, determinism, performance) and hostile-input fixtures. |

## Params and wiring (planner)

```groovy
// nextflow.config params { }
lineage_summary_xlsx          = true   // write summary/lineage_summary.xlsx next to the CSV
lineage_summary_xlsx_decimals = 3      // decimals shown for abundances (display only; stored values are the CSV's)
```

Schema: `lineage_summary_xlsx` `boolean`; `lineage_summary_xlsx_decimals` `integer`, `minimum 0`, `maximum 6`. (The
second one is optional; drop it and the default of 3 applies.)

```groovy
// conf/modules.config: extend the existing LINEAGE_SUMMARY ext.args
withName: 'LINEAGE_SUMMARY' {
    ext.args = {
        def xlsx = params.lineage_summary_xlsx
            ? " --xlsx-out lineage_summary.xlsx --xlsx-decimals ${params.lineage_summary_xlsx_decimals} --wasteline-version ${workflow.manifest.version}"
            : ''
        "--major-threshold ${params.lineage_major_threshold} --min-abundance ${params.lineage_min_abundance} --zeros ${params.lineage_summary_zeros}${xlsx}"
    }
    publishDir = [
        path: { "${params.outdir}/summary" },
        mode: params.publish_dir_mode,
        pattern: 'lineage_*'          // lineage_summary.xlsx starts with lineage_ too: no change needed
    ]
}
```

No workflow change is needed: the new `LINEAGE_SUMMARY.out.xlsx` channel is optional and does not have to be
consumed (the publishDir pattern picks the file up). The `--ungrouped-out` and `--sites-out` tables already passed by the module become sheets when present.

A workbook problem never costs the CSV results: if the workbook cannot be written, `lineage_groups.py` prints a
`WARNING: could not write ...` and exits 0 with all CSV/TSV files in place (`xlsx` is then simply absent from the
outputs).

In a `-stub` run the module leaves an **empty placeholder** `lineage_summary.xlsx` (only when `--xlsx-out` is in
`ext.args`); the real run writes a valid workbook.

## What is in the workbook

| Sheet | Content |
| --- | --- |
| `Lineage summary` | The CSV, same columns in the same order. Group columns, `other`, `proportion_summary` and `genome_coverage_pct` are numeric cells; blank zeros are empty cells; `major_*`/`minor_*`/`ungrouped_lineages`/`status_reason` are text. |
| `Ungrouped lineages` | `lineage_ungrouped_summary.tsv` (lineages falling into `other`, biggest first, with the suggested pattern). |
| `Site aggregates` | `lineage_summary_site_aggregates.csv`, only when site aggregates were produced. |
| `Groups used` | The groups: `group, include, exclude, in_total`. |
| `About` | `wasteline <version>`, Freyja dataset version, groups-file sha256, the thresholds used and short notes on the columns. No timestamp, so the file is reproducible. |

Formatting: bold header with a light grey fill, wrapped, thin bottom border; header row and first column frozen;
auto-filter on the header; column widths from the content (capped at 50 characters; the long lineage lists wrap and
Excel autosizes the row heights on open because no explicit height is stored); numbers right-aligned with format
`0.000` (coverage `0.00`), text left-aligned; sample IDs and dates are always text (never reinterpreted, `007` stays
`007`). Group and `other` abundance cells are shaded by six static blue buckets: `0<x<0.05`, `<0.10`, `<0.20`, `<0.35`,
`<0.50`, `>=0.50` (blank and zero cells are not shaded).

Safety: every string is XML-escaped and written as an inline *string* cell, never as a formula (a sample ID such as
`=cmd|...`, `+1`, `-2` or `@SUM(1)` stays text); control characters XML forbids are dropped; cells are cut at Excel's
32,767-character limit; sheet names are cleaned to Excel's rules (<= 31 characters, none of `[]:*?/\`, unique). Two
runs on the same input give byte-identical files (fixed zip timestamps, order and compression).

## Converting other tables

```bash
bin/xlsx_writer.py --from-csv summary/qc_summary.tsv --sheet-name QC \
                   --from-csv summary/demix_status.tsv --sheet-name "Demix status" \
                   -o summary/qc_summary.xlsx
```

`--from-csv` is repeatable (one sheet each, delimiter `,`, tab or `;` is detected, UTF-8 with or without BOM);
`--sheet-name` is matched by position (default: the file name); `--heat-column HEADER` shades a numeric column;
`--freeze-rows` / `--freeze-cols` default to 1 / 1. In the converter every column except the first whose non-empty cells
are all plain numbers (no leading zeros, no `nan`/`inf`) becomes a numeric column shown with as many decimals as its
data has; everything else stays text. From Python:

```python
from xlsx_writer import write_workbook
write_workbook("out.xlsx", [{"name": "QC", "header": [...], "rows": [[...]], "numeric_columns": [2, 3],
                             "number_formats": {2: "0.0"}, "heat_columns": [3], "freeze": (1, 1)}])
```

## Limitations

* No formulas, conditional formatting, charts, tables objects, merged cells, images, macros or data validation. The
  shading is static (computed when the file is written), so it does not follow later edits.
* No shared strings (inline strings are used): fine for Excel/Numbers/LibreOffice/pandas/openpyxl; Excel converts them
  when the user saves.
* `major_lineages` and `major_abundances` wrap independently, so the names and values do not line up line by line in
  a wrapped cell (same text as the CSV).
* Row heights are left to the spreadsheet application; Excel, Numbers and LibreOffice autosize them, a viewer that
  does not (some web previews) shows wrapped rows at one line.
* Dates are text (`2026-01-12`), as in the CSV, so Excel does not re-interpret them by locale.
* Very large tables: one sheet holds at most 1,048,575 data rows and 16,384 columns (an error is raised beyond that).

## What to document in `docs/usage.md` / `docs/output.md`

* `usage.md`: `--lineage_summary_xlsx` (default `true`) and `--lineage_summary_xlsx_decimals` (default 3, display only);
  that the CSV is unaffected; that `bin/xlsx_writer.py --from-csv` converts any other table.
* `output.md`: `summary/lineage_summary.xlsx` with the five sheets above, the formatting notes (blank = none, shaded
  buckets, wrapped lineage lists, frozen header/first column, filter), "the CSV is the source of truth; the workbook is a
  formatted copy for people", and that no date/time is stored so the file only changes when the results change.
* The HTML report needs no change.
