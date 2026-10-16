#!/usr/bin/env python3
"""Turn a Freyja aggregate plus sample metadata into a long time-series table.

Input: `freyja_aggregated.tsv` (Freyja's aggregate; list-valued columns are Python/numpy reprs) and a
samplesheet-like CSV with `sample, sample_type, collection_date, site` and an optional numeric
`viral_load` (for example gene copies per litre; empty = unknown). Every optional column may be absent.

Output (`freyja_timeseries.tsv`, tab-separated, header always written):
`site, collection_date, lineage_group, abundance, abundance_load_weighted, n_samples,
mean_viral_load, total_viral_load`.

* `abundance` is the mean over the samples of that site and date (a group a sample does not
  contain counts as 0).
* `abundance_load_weighted` is sum(abundance_i * load_i) / sum(load_i) over the samples that have a
  viral load (`NA` when none has).
* `n_samples` counts every contributing sample; `mean_viral_load` / `total_viral_load` use the ones
  with a load.
* Excluded, with a count on stderr: samples without metadata, negative controls
  (`sample_type` negative), samples without a valid ISO `collection_date`, samples below
  `--min-coverage` (percent of genome; Freyja's `coverage` column) and unparsable aggregate rows.
* `--level summarized` (default) uses Freyja's lineage groups, `--level lineages` the full lineages.
* Rows are sorted by site, date, then abundance (descending). An empty site is written as `NA`.

Helper mode for `freyja relgrowthrate` / `freyja dash`: `--freyja-metadata-out` and
`--freyja-aggregated-out` write Freyja's metadata CSV (`Sample,sample_collection_datetime[,viral_load]`)
and the aggregate restricted to the samples that have it (Freyja fails on a sample without metadata),
in the aggregate's row order. See `prepare_freyja`.

Standard library only.
"""

import argparse
import csv
import datetime as dt
import json
import math
import re
import sys

OUT_COLUMNS = [
    "site",
    "collection_date",
    "lineage_group",
    "abundance",
    "abundance_load_weighted",
    "n_samples",
    "mean_viral_load",
    "total_viral_load",
]
NA = "NA"

_QUOTED = re.compile(r"""'([^']*)'|"([^"]*)\"""")
_NUMBER = re.compile(r"[-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?|\bnan\b|\binf\b", re.IGNORECASE)
_PAIR = re.compile(r"""\(\s*(?:'([^']*)'|"([^"]*)")\s*,\s*([^()]*(?:\([^()]*\))?[^()]*?)\s*\)""")
_DATE = re.compile(r"^(\d{4})-(\d{2})-(\d{2})(?:$|[T ])")


# --------------------------------------------------------------------------------------------------
# parsing
# --------------------------------------------------------------------------------------------------
def is_na(value):
    return value is None or str(value).strip() in ("", "NA", "NaN", "nan", "None", "[]")


def to_float(value):
    """Finite float or None."""
    try:
        number = float(str(value).strip())
    except (TypeError, ValueError):
        return None
    return number if math.isfinite(number) else None


def parse_float_list(text):
    """Numbers in a Python list repr or a numpy array repr (whitespace/newline separated)."""
    if is_na(text):
        return []
    values = []
    for token in _NUMBER.findall(str(text)):
        try:
            number = float(token)
        except ValueError:
            continue
        values.append(number if math.isfinite(number) else 0.0)
    return values


def parse_name_list(text):
    """Names of a Python list repr, a numpy array repr (`['A' 'B']`) or plain `A B` text."""
    if is_na(text):
        return []
    quoted = [a or b for a, b in _QUOTED.findall(str(text))]
    if quoted:
        return quoted
    return str(text).strip().strip("[]").replace(",", " ").split()


def parse_summarized(text):
    """`[('Omicron', 0.7), ('Delta', 0.3)]` (also with `np.float64(0.3)` values) -> [(name, frac)]."""
    if is_na(text):
        return []
    pairs = []
    for match in _PAIR.finditer(str(text)):
        name = match.group(1) if match.group(1) is not None else match.group(2)
        numbers = _NUMBER.findall(match.group(3))
        value = to_float(numbers[-1]) if numbers else None
        if value is not None:
            pairs.append((name, value))
    return pairs


def parse_date(text):
    """ISO `YYYY-MM-DD` (a time part is tolerated and dropped) -> `YYYY-MM-DD`, or None."""
    match = _DATE.match(str(text or "").strip())
    if not match:
        return None
    try:
        return dt.date(int(match.group(1)), int(match.group(2)), int(match.group(3))).isoformat()
    except ValueError:
        return None


def read_aggregate(path):
    """Parse a Freyja aggregate into rows.

    Returns `(rows, n_malformed, header)`; a row is a dict `sample, summarized, lineages` (both lists of
    `(name, fraction)`), `coverage` (float or None), `raw` (the original cells) and `index`. Rows whose
    label is empty, or whose groups cannot be parsed, are counted in `n_malformed` and skipped.
    """
    with open(path, newline="") as handle:
        table = list(csv.reader(handle, delimiter="\t"))
    if not table:
        return [], 0, []
    header = [cell.strip() for cell in table[0]]
    col = {name: i for i, name in enumerate(header) if name}

    def cell(row, name):
        i = col.get(name)
        return row[i] if i is not None and i < len(row) else ""

    rows, malformed = [], 0
    for line in table[1:]:
        if not line or not any(c.strip() for c in line):
            continue
        sample = line[0].strip()
        summarized = parse_summarized(cell(line, "summarized"))
        names = parse_name_list(cell(line, "lineages"))
        values = parse_float_list(cell(line, "abundances"))
        lineages = list(zip(names, values)) if len(names) == len(values) else []
        if not sample or (not summarized and not lineages):
            malformed += 1
            continue
        rows.append(
            {
                "sample": sample,
                "summarized": summarized,
                "lineages": lineages,
                "coverage": to_float(cell(line, "coverage")) if "coverage" in col else None,
                "raw": line,
            }
        )
    return rows, malformed, table[0]


def read_metadata(path, counts=None):
    """Read the samplesheet-like CSV -> {sample: dict(sample_type, date, site, viral_load)}.

    Column names are matched case-insensitively; `sample_type`, `site` and `viral_load` may be
    absent (defaults: `sample`, `NA`, unknown). `date` is None for a blank or invalid date.
    """
    counts = counts if counts is not None else {}
    meta = {}
    with open(path, newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        reader.fieldnames = [(name or "").strip().lower() for name in (reader.fieldnames or [])]
        if "sample" not in reader.fieldnames:
            raise ValueError(f"{path}: no 'sample' column (found {reader.fieldnames})")
        for record in reader:
            sample = (record.get("sample") or "").strip()
            if not sample:
                continue
            if sample in meta:
                counts["duplicate_metadata"] = counts.get("duplicate_metadata", 0) + 1
                continue
            raw_date = (record.get("collection_date") or "").strip()
            load_text = (record.get("viral_load") or "").strip()
            load = to_float(load_text) if load_text else None
            if load_text and (load is None or load <= 0):
                counts["invalid_viral_load"] = counts.get("invalid_viral_load", 0) + 1
                load = None
            meta[sample] = {
                "sample_type": (record.get("sample_type") or "sample").strip().lower() or "sample",
                "date": parse_date(raw_date),
                "date_given": bool(raw_date),
                "site": (record.get("site") or "").strip() or NA,
                "viral_load": load,
            }
    return meta


# --------------------------------------------------------------------------------------------------
# selection and aggregation
# --------------------------------------------------------------------------------------------------
def select_samples(rows, meta, min_coverage):
    """Split aggregate rows into usable ones (with their metadata) and exclusion counts."""
    used = []
    excluded = {
        "no_metadata": 0,
        "negative_control": 0,
        "no_date": 0,
        "bad_date": 0,
        "low_coverage": 0,
    }
    for row in rows:
        info = meta.get(row["sample"])
        if info is None:
            excluded["no_metadata"] += 1
        elif info["sample_type"] == "negative":
            excluded["negative_control"] += 1
        elif not info["date_given"]:
            excluded["no_date"] += 1
        elif info["date"] is None:
            excluded["bad_date"] += 1
        elif row["coverage"] is not None and row["coverage"] < min_coverage:
            excluded["low_coverage"] += 1
        else:
            used.append((row, info))
    return used, excluded


def build_series(used, level="summarized", min_abundance=0.0):
    """Long rows (dicts of OUT_COLUMNS, numbers unformatted) for the usable samples."""
    cells = {}
    for row, info in used:
        cells.setdefault((info["site"], info["date"]), []).append((row[level], info["viral_load"]))
    out = []
    for (site, date), samples in cells.items():
        n = len(samples)
        loads = [load for _, load in samples if load is not None]
        total_load = sum(loads) if loads else None
        mean_load = total_load / len(loads) if loads else None
        groups = sorted({name for fractions, _ in samples for name, _ in fractions})
        for group in groups:
            values = []
            for fractions, load in samples:
                value = sum(frac for name, frac in fractions if name == group)
                values.append((value, load))
            abundance = sum(v for v, _ in values) / n
            if abundance <= 0 or abundance < min_abundance:
                continue
            weighted = (
                sum(v * load for v, load in values if load is not None) / total_load
                if loads
                else None
            )
            out.append(
                {
                    "site": site,
                    "collection_date": date,
                    "lineage_group": group,
                    "abundance": abundance,
                    "abundance_load_weighted": weighted,
                    "n_samples": n,
                    "mean_viral_load": mean_load,
                    "total_viral_load": total_load,
                }
            )
    out.sort(key=lambda r: (r["site"], r["collection_date"], -r["abundance"], r["lineage_group"]))
    return out


def fmt(value, places=6):
    return NA if value is None else f"{value:.{places}f}"


def write_series(path, series):
    lines = ["\t".join(OUT_COLUMNS)]
    for r in series:
        cells = [
            r["site"],
            r["collection_date"],
            r["lineage_group"],
            fmt(r["abundance"]),
            fmt(r["abundance_load_weighted"]),
            str(r["n_samples"]),
            fmt(r["mean_viral_load"], 2),
            fmt(r["total_viral_load"], 2),
        ]
        lines.append("\t".join(cells))
    with open(path, "w", newline="") as handle:
        handle.write("\n".join(lines) + "\n")


# --------------------------------------------------------------------------------------------------
# helper mode: Freyja's own metadata format
# --------------------------------------------------------------------------------------------------
def prepare_freyja(
    header,
    used,
    metadata_out,
    aggregated_out,
    require_viral_load=False,
    summary_out=None,
    window_days=None,
):
    """Write Freyja's metadata CSV and the matching aggregate; return the summary dict.

    `freyja relgrowthrate` / `freyja dash` index the metadata by the aggregate's sample label and fail
    on a sample without a row, so the aggregate is cut to the samples that have usable metadata
    (known, not negative, valid date, coverage ok). With `require_viral_load` samples without a load
    are dropped too (needed for `--scale_by_viral_load`: Freyja cannot scale a date without loads).
    The `viral_load` column is written only when at least one kept sample has a load. Rows keep the
    aggregate's order. The summary (also written as JSON to `summary_out`) counts the kept samples and
    distinct dates, and the dates and span inside the last `window_days` days (Freyja fits only that
    window), which is what decides whether a growth-rate fit can be trusted.
    """
    kept = [
        (row, info)
        for row, info in used
        if info["viral_load"] is not None or not require_viral_load
    ]
    has_load = any(info["viral_load"] is not None for _, info in kept)
    with open(metadata_out, "w", newline="") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(
            ["Sample", "sample_collection_datetime"] + (["viral_load"] if has_load else [])
        )
        for row, info in kept:
            extra = []
            if has_load:
                extra = ["" if info["viral_load"] is None else repr(info["viral_load"])]
            writer.writerow([row["sample"], info["date"]] + extra)
    with open(aggregated_out, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(header)
        for row, _ in kept:
            writer.writerow(row["raw"])
    dates = sorted({info["date"] for _, info in kept})
    last = dt.date.fromisoformat(dates[-1]) if dates else None
    in_window = [
        d
        for d in dates
        if window_days is None or (last - dt.date.fromisoformat(d)).days <= window_days
    ]
    summary = {
        "n_samples": len(kept),
        "n_dates": len(dates),
        "first_date": dates[0] if dates else None,
        "last_date": dates[-1] if dates else None,
        "span_days": (dt.date.fromisoformat(dates[-1]) - dt.date.fromisoformat(dates[0])).days
        if dates
        else 0,
        "window_days": window_days,
        "n_dates_window": len(in_window),
        "window_span_days": (last - dt.date.fromisoformat(in_window[0])).days if in_window else 0,
        "n_with_viral_load": sum(1 for _, info in kept if info["viral_load"] is not None),
        "viral_load_column": has_load,
    }
    if summary_out:
        with open(summary_out, "w") as handle:
            json.dump(summary, handle, indent=2, sort_keys=True)
            handle.write("\n")
    return summary


# --------------------------------------------------------------------------------------------------
def parse_args(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--aggregated", required=True, help="freyja_aggregated.tsv")
    parser.add_argument("--metadata", required=True, help="samplesheet-like CSV")
    parser.add_argument("--level", choices=["summarized", "lineages"], default="summarized")
    parser.add_argument(
        "--min-abundance",
        type=float,
        default=0.0,
        help="drop rows whose (unweighted) mean abundance is below this fraction (default 0)",
    )
    parser.add_argument(
        "--min-coverage",
        type=float,
        default=60.0,
        help="exclude samples with a Freyja `coverage` (percent) below this (default 60)",
    )
    parser.add_argument("-o", "--output", help="freyja_timeseries.tsv")
    parser.add_argument("--freyja-metadata-out", help="write Freyja's metadata CSV here")
    parser.add_argument("--freyja-aggregated-out", help="write the matching aggregate here")
    parser.add_argument("--freyja-summary-out", help="write a JSON summary of that selection here")
    parser.add_argument(
        "--require-viral-load",
        action="store_true",
        help="helper mode: keep only samples with a viral load (for --scale_by_viral_load)",
    )
    parser.add_argument(
        "--window-days",
        type=int,
        help="helper mode: Freyja's --days; the summary counts the dates inside that window",
    )
    args = parser.parse_args(argv)
    helper = [args.freyja_metadata_out, args.freyja_aggregated_out]
    if any(helper) and not all(helper):
        parser.error("--freyja-metadata-out and --freyja-aggregated-out go together")
    if not args.output and not any(helper):
        parser.error("give -o/--output and/or --freyja-metadata-out with --freyja-aggregated-out")
    return args


def main(argv=None):
    args = parse_args(argv)
    counts = {}
    try:
        rows, malformed, header = read_aggregate(args.aggregated)
        meta = read_metadata(args.metadata, counts)
    except (OSError, ValueError, csv.Error) as err:
        print(f"timeseries.py: {err}", file=sys.stderr)
        return 1
    used, excluded = select_samples(rows, meta, args.min_coverage)
    if args.output:
        write_series(args.output, build_series(used, args.level, args.min_abundance))
    if args.freyja_metadata_out:
        prepare_freyja(
            header,
            used,
            args.freyja_metadata_out,
            args.freyja_aggregated_out,
            args.require_viral_load,
            args.freyja_summary_out,
            args.window_days,
        )
    dropped = ", ".join(f"{k}={v}" for k, v in {**excluded, "unparsable_rows": malformed}.items())
    notes = ", ".join(f"{k}={v}" for k, v in counts.items())
    print(
        f"timeseries.py: {len(rows)} aggregate samples, {len(used)} used; excluded: {dropped}"
        + (f"; metadata notes: {notes}" if notes else ""),
        file=sys.stderr,
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
