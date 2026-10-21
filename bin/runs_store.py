#!/usr/bin/env python3
"""A persistent cross-run store of wasteline results, with trend and emergence alerts (screening hints).

    runs_store.py add    <store.sqlite> <results_dir> [<results_dir> ...]
    runs_store.py runs   <store.sqlite>
    runs_store.py alerts <store.sqlite> [--watchlist mutations.csv] [--since YYYY-MM-DD]
                         [--out alerts.tsv] [--markdown alerts.md] [thresholds, see --help]

`add` reads the `summary/` tables and `pipeline_info/run_info.json` of finished results directories (or of the
results directories one level below a folder) into one SQLite file. Adding a directory again replaces its rows
(matched by run identity `run_name@start` or by path), so the command can be repeated safely. The same sample may
be in several runs (for example a re-analysis with newer barcodes): all of them are kept, each with its barcode
version, and queries use the latest run per sample (newest barcode version, then newest run start).

`alerts` applies four rules to the samples of the store; they are screening hints, not epidemiological conclusions.
Only samples with status ok or low_confidence, a date and a sample_type other than `negative` count. The date is
the collection date, else the date of the run. Periods are ISO weeks. See docs/runs_store.md for the rules.

Runs on the user's host Python: standard library only (sqlite3), Python >= 3.8 syntax. Exit status: 0 ok,
1 error (unreadable input, nothing added, bad option).
"""

from __future__ import annotations

import argparse
import csv
import datetime
import json
import os
import re
import sqlite3
import sys
from pathlib import Path

SCHEMA_VERSION = 1
COUNTED = ("ok", "low_confidence")
NO_SITE = "(none)"
CATCH_ALL = {"other"}  # the lineage summary's catch-all column: never alerted on
RULES = ("emerging", "sustained_growth", "recurring_unexplained", "watchlist")
ALERT_COLUMNS = [
    "rule",
    "level",
    "site",
    "group_or_mutation",
    "detail",
    "first_date",
    "last_date",
    "samples",
]

SCHEMA = """
CREATE TABLE runs (
    run_id TEXT PRIMARY KEY, name TEXT, path TEXT, mode TEXT, start TEXT, run_date TEXT, success TEXT,
    pipeline_version TEXT, freyja_version TEXT, barcode_version TEXT, barcode_iso TEXT, added_at TEXT);
CREATE TABLE samples (
    run_id TEXT, sample TEXT, site TEXT, collection_date TEXT, sample_type TEXT, status TEXT, reason TEXT,
    genome_coverage_pct REAL, qc_status TEXT, barcode_version TEXT, PRIMARY KEY (run_id, sample));
CREATE TABLE lineage_abundance (
    run_id TEXT, sample TEXT, lineage TEXT, abundance REAL, lo REAL, hi REAL);
CREATE TABLE group_abundance (
    run_id TEXT, sample TEXT, grp TEXT, abundance REAL, lo REAL, hi REAL);
CREATE TABLE unexplained (
    run_id TEXT, sample TEXT, position INTEGER, ref TEXT, alt TEXT, mutation TEXT, aa TEXT,
    observed_freq REAL, expected_freq REAL, depth REAL, alt_depth REAL, n_samples INTEGER, in_dataset TEXT);
CREATE TABLE moi (
    run_id TEXT, sample TEXT, name TEXT, mutation TEXT, depth REAL, alt_depth REAL, alt_freq REAL);
CREATE INDEX idx_lineage ON lineage_abundance (run_id, sample);
CREATE INDEX idx_group ON group_abundance (run_id, sample);
CREATE INDEX idx_unexplained ON unexplained (run_id, sample);
CREATE INDEX idx_moi ON moi (run_id, sample);
"""
CHILD_TABLES = ("samples", "lineage_abundance", "group_abundance", "unexplained", "moi")


class Problem(Exception):
    """A friendly, fatal error (printed without a traceback)."""


# =========================================================================================================
# reading results
# =========================================================================================================


def num(text):
    """float, or None for empty / NA / unparseable."""
    text = (text or "").strip()
    if not text or text.upper() in ("NA", "NAN", "NONE"):
        return None
    try:
        return float(text)
    except ValueError:
        return None


def text_or_empty(text):
    text = (text or "").strip()
    return "" if text.upper() in ("NA", "NAN", "NONE") else text


def iso_date(text):
    """`YYYY-MM-DD` (also the date part of a timestamp) -> that string; anything else -> ''."""
    match = re.match(r"(\d{4}-\d{2}-\d{2})", text or "")
    if not match:
        return ""
    try:
        datetime.date.fromisoformat(match.group(1))
    except ValueError:
        return ""
    return match.group(1)


def read_tsv(path, delimiter="\t"):
    """[dict per row] of a delimited file; {} columns missing in a short row are ''. [] when absent."""
    path = Path(path)
    if not path.is_file():
        return []
    with open(str(path), encoding="utf-8-sig", newline="") as handle:
        rows = [r for r in csv.reader(handle, delimiter=delimiter) if any(c.strip() for c in r)]
    if not rows:
        return []
    header = rows[0]
    return [{h: (r[i] if i < len(r) else "") for i, h in enumerate(header)} for r in rows[1:]]


def read_header(path, delimiter=","):
    with open(str(path), encoding="utf-8-sig", newline="") as handle:
        return next(csv.reader(handle, delimiter=delimiter), [])


def newest(directory, pattern):
    hits = sorted(Path(directory).glob(pattern), key=lambda p: (p.stat().st_mtime, p.name))
    return hits[-1] if hits else None


def parse_aggregated(path):
    """{sample: {lineage: abundance}} from freyja_aggregated.tsv (names and abundances are space separated)."""
    result = {}
    with open(str(path), encoding="utf-8-sig", newline="") as handle:
        rows = [r for r in csv.reader(handle, delimiter="\t") if any(c.strip() for c in r)]
    if not rows:
        return result
    header = rows[0]
    at = {name: i for i, name in enumerate(header) if name}
    li, ai = at.get("lineages"), at.get("abundances")
    if li is None or ai is None:
        return result
    for row in rows[1:]:
        if not row or not row[0].strip() or len(row) <= max(li, ai):
            continue
        names = [t for t in re.split(r"[\s,\[\]()'\"]+", row[li]) if t]
        values = [num(t) for t in re.split(r"[\s,\[\]()'\"]+", row[ai]) if t]
        if len(names) == len(values) and None not in values:
            result[row[0].strip()] = dict(zip(names, values))
    return result


def parse_group_table(path):
    """({sample: {group: abundance}}) from lineage_summary.csv: the columns between genome_coverage_pct (or
    status) and proportion_summary; an empty cell is 0."""
    header = read_header(path)
    if "sample" not in header or "proportion_summary" not in header:
        return {}
    if "genome_coverage_pct" in header:
        start = header.index("genome_coverage_pct") + 1
    else:
        start = header.index("status") + 1 if "status" in header else 1
    columns = header[start : header.index("proportion_summary")]
    result = {}
    for row in read_tsv(path, ","):
        sample = row.get("sample", "").strip()
        if sample:
            result[sample] = {c: (num(row.get(c, "")) or 0.0) for c in columns}
    return result


class RunData:
    """Everything of one results directory that goes into the store."""

    def __init__(self, directory):
        self.dir = Path(os.path.abspath(str(directory)))
        self.path = str(self.dir)
        self.warnings = []
        summary = self.dir / "summary"
        status_file = summary / "demix_status.tsv"
        if not status_file.is_file():
            raise Problem(
                f"{directory}: no summary/demix_status.tsv (not a finished results directory)"
            )
        info = {}
        info_path = self.dir / "pipeline_info" / "run_info.json"
        if info_path.is_file():
            try:
                with open(str(info_path)) as handle:
                    info = json.load(handle)
            except (OSError, ValueError) as err:
                self.warnings.append(
                    f"{info_path}: unreadable ({err}); the run is identified by its path"
                )
        else:
            self.warnings.append(
                f"{directory}: no pipeline_info/run_info.json; the run is identified by its path"
            )
        self.name = self.dir.name
        start = str(info.get("start") or "")
        if info.get("run_name") and start:
            self.run_id = f"{info['run_name']}@{start}"
        else:
            self.run_id = "path:" + self.path
        self.start = start
        self.run_date = iso_date(start)
        block = info.get("freyja_data") or {}
        self.meta = {
            "mode": str(info.get("mode") or ""),
            "success": "" if "success" not in info else str(info["success"]).lower(),
            "pipeline_version": str(info.get("pipeline_version") or ""),
            "freyja_version": str(info.get("freyja_version") or ""),
            "barcode_version": str(block.get("barcode_version") or ""),
            "barcode_iso": str(block.get("barcode_version_iso") or ""),
        }
        self.load_samples(summary, status_file)
        self.lineages = self.load_lineages(summary)
        self.groups = self.load_groups(summary)
        self.unexplained = read_tsv(summary / "unexplained_mutations.tsv")
        self.moi = read_tsv(summary / "mutations_of_interest.tsv")

    def load_samples(self, summary, status_file):
        sheet = {r.get("sample", ""): r for r in read_tsv(summary / "samplesheet.csv", ",")}
        qc_path = summary / "qc_summary.tsv"
        if not qc_path.is_file():
            qc_path = newest(summary, "qc_summary*.tsv") or qc_path
        qc = {r.get("sample", ""): r for r in read_tsv(qc_path)}
        self.samples = []
        for row in read_tsv(status_file):
            sample = row.get("sample", "").strip()
            if not sample:
                continue
            sheet_row, qc_row = sheet.get(sample, {}), qc.get(sample, {})

            def pick(column, sheet_row=sheet_row, qc_row=qc_row):
                return text_or_empty(sheet_row.get(column, "")) or text_or_empty(
                    qc_row.get(column, "")
                )

            raw_date = pick("collection_date")
            date = iso_date(raw_date)
            if raw_date and not date:
                self.warnings.append(
                    f"{sample}: collection_date {raw_date!r} is not YYYY-MM-DD: ignored"
                )
            self.samples.append(
                {
                    "sample": sample,
                    "site": pick("site"),
                    "collection_date": date,
                    "sample_type": pick("sample_type") or "sample",
                    "status": row.get("status", "").strip(),
                    "reason": row.get("reason", "").strip(),
                    "genome_coverage_pct": num(row.get("genome_coverage_pct", "")),
                    "qc_status": text_or_empty(qc_row.get("qc_status", "")),
                    "barcode_version": row.get("barcode_version", "").strip()
                    or self.meta["barcode_version"],
                }
            )
        self.known = {s["sample"] for s in self.samples}

    def load_lineages(self, summary):
        path = summary / "freyja_aggregated.tsv"
        found = parse_aggregated(path) if path.is_file() else {}
        boot = {}
        for row in read_tsv(summary / "freyja_boot_summary.tsv"):
            boot[(row.get("sample", ""), row.get("lineage", ""))] = (
                num(row.get("lo", "")),
                num(row.get("hi", "")),
            )
        rows = []
        for sample, items in found.items():
            if sample not in self.known:
                self.warnings.append(
                    f"freyja_aggregated.tsv: sample {sample!r} is not in demix_status.tsv: left out"
                )
                continue
            for lineage, value in items.items():
                lo, hi = boot.get((sample, lineage), (None, None))
                rows.append((sample, lineage, value, lo, hi))
        return rows

    def load_groups(self, summary):
        path = summary / "lineage_summary.csv"
        found = parse_group_table(path) if path.is_file() else {}
        intervals = {}
        for row in read_tsv(summary / "lineage_group_intervals.tsv"):
            intervals[(row.get("sample", ""), row.get("group", ""))] = (
                num(row.get("lo", "")),
                num(row.get("hi", "")),
            )
        rows = []
        for sample, items in found.items():
            if sample not in self.known:
                continue
            for group, value in items.items():
                lo, hi = intervals.get((sample, group), (None, None))
                rows.append((sample, group, value, lo, hi))
        return rows


# =========================================================================================================
# the store
# =========================================================================================================


def open_store(path, create):
    path = Path(path)
    if not path.is_file():
        if not create:
            raise Problem(f"{path}: no such store (create it with `runs_store.py add`)")
        path.parent.mkdir(parents=True, exist_ok=True)
    conn = sqlite3.connect(str(path))
    conn.row_factory = sqlite3.Row
    try:
        version = conn.execute("PRAGMA user_version").fetchone()[0]
        if version == 0 and not conn.execute("SELECT 1 FROM sqlite_master LIMIT 1").fetchone():
            conn.executescript(SCHEMA)
            conn.execute(f"PRAGMA user_version = {SCHEMA_VERSION}")
            conn.commit()
        elif version != SCHEMA_VERSION:
            conn.close()
            raise Problem(f"{path}: not a runs store of this version (user_version {version})")
    except sqlite3.DatabaseError as err:
        conn.close()
        raise Problem(f"{path}: not an SQLite runs store ({err})")
    return conn


def find_results(directory):
    """The results directories at `directory` or one level below it."""
    base = Path(directory)
    if not base.is_dir():
        raise Problem(f"not a directory: {directory}")
    if (base / "summary" / "demix_status.tsv").is_file():
        return [base]
    kids = sorted(p for p in base.iterdir() if (p / "summary" / "demix_status.tsv").is_file())
    if not kids:
        raise Problem(f"{directory}: no summary/demix_status.tsv here or one level below")
    return kids


def store_run(conn, run):
    """Replace the rows of `run` (matched by run_id or path) with the current ones. One transaction."""
    old = [
        r[0]
        for r in conn.execute(
            "SELECT run_id FROM runs WHERE run_id = ? OR path = ?", (run.run_id, run.path)
        )
    ]
    with conn:
        for run_id in old:
            for table in CHILD_TABLES:
                conn.execute(f"DELETE FROM {table} WHERE run_id = ?", (run_id,))
            conn.execute("DELETE FROM runs WHERE run_id = ?", (run_id,))
        rid = run.run_id
        added = datetime.datetime.now().strftime("%Y-%m-%dT%H:%M:%S")
        conn.execute(
            "INSERT INTO runs VALUES (?,?,?,?,?,?,?,?,?,?,?,?)",
            (
                rid,
                run.name,
                run.path,
                run.meta["mode"],
                run.start,
                run.run_date,
                run.meta["success"],
                run.meta["pipeline_version"],
                run.meta["freyja_version"],
                run.meta["barcode_version"],
                run.meta["barcode_iso"],
                added,
            ),
        )
        conn.executemany(
            "INSERT INTO samples VALUES (?,?,?,?,?,?,?,?,?,?)",
            [
                (
                    rid,
                    s["sample"],
                    s["site"],
                    s["collection_date"],
                    s["sample_type"],
                    s["status"],
                    s["reason"],
                    s["genome_coverage_pct"],
                    s["qc_status"],
                    s["barcode_version"],
                )
                for s in run.samples
            ],
        )
        conn.executemany(
            "INSERT INTO lineage_abundance VALUES (?,?,?,?,?,?)", [(rid,) + r for r in run.lineages]
        )
        conn.executemany(
            "INSERT INTO group_abundance VALUES (?,?,?,?,?,?)", [(rid,) + r for r in run.groups]
        )
        unexplained = []
        for r in run.unexplained:
            sample = r.get("sample", "")
            aa = text_or_empty(r.get("aa_change", ""))
            gene = text_or_empty(r.get("gene", ""))
            unexplained.append(
                (
                    rid,
                    sample,
                    int(num(r.get("position", "")) or 0),
                    r.get("ref", ""),
                    r.get("alt", ""),
                    f"{r.get('ref', '')}{int(num(r.get('position', '')) or 0)}{r.get('alt', '')}",
                    f"{gene}:{aa}" if gene and aa else "",
                    num(r.get("observed_freq", "")),
                    num(r.get("expected_freq", "")),
                    num(r.get("depth", "")),
                    num(r.get("alt_depth", "")),
                    int(num(r.get("n_samples", "")) or 0),
                    r.get("in_dataset", "").strip().lower(),
                )
            )
        conn.executemany(
            "INSERT INTO unexplained VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?)",
            [u for u in unexplained if u[1] in run.known],
        )
        conn.executemany(
            "INSERT INTO moi VALUES (?,?,?,?,?,?,?)",
            [
                (
                    rid,
                    r.get("sample", ""),
                    r.get("name", ""),
                    r.get("mutation", ""),
                    num(r.get("depth", "")),
                    num(r.get("alt_depth", "")),
                    num(r.get("alt_freq", "")),
                )
                for r in run.moi
                if r.get("sample", "") in run.known
            ],
        )
    return len(old)


def cmd_add(args):
    conn = open_store(args.store, create=True)
    added = 0
    try:
        for given in args.results:
            for directory in find_results(given):
                try:
                    run = RunData(directory)
                except Problem as err:
                    print(f"WARNING: {err}", file=sys.stderr)
                    continue
                for warning in run.warnings[:10]:
                    print(f"WARNING: {warning}", file=sys.stderr)
                if len(run.warnings) > 10:
                    print(
                        f"WARNING: ... and {len(run.warnings) - 10} more for {run.name}",
                        file=sys.stderr,
                    )
                replaced = store_run(conn, run)
                added += 1
                print(
                    "{} {}: {} samples, {} lineage rows, {} group rows, {} unexplained, {} moi{}".format(
                        "replaced" if replaced else "added",
                        run.name,
                        len(run.samples),
                        len(run.lineages),
                        len(run.groups),
                        len(run.unexplained),
                        len(run.moi),
                        "" if run.meta["barcode_version"] else " (barcode version unknown)",
                    )
                )
    finally:
        conn.close()
    if not added:
        raise Problem("nothing was added")
    return 0


def print_table(rows, header):
    widths = [max(len(str(r[i])) for r in [header] + rows) for i in range(len(header))]
    for r in [header] + rows:
        print("  " + "  ".join(str(c).ljust(widths[i]) for i, c in enumerate(r)).rstrip())


def cmd_runs(args):
    conn = open_store(args.store, create=False)
    try:
        runs = conn.execute("SELECT * FROM runs ORDER BY run_date, start, name").fetchall()
        if not runs:
            print("The store holds no runs.")
            return 0
        rows = []
        for run in runs:
            counts = {"ok": 0, "low_confidence": 0, "other": 0, "dated": 0}
            for s in conn.execute(
                "SELECT status, collection_date FROM samples WHERE run_id = ?", (run["run_id"],)
            ):
                counts[s["status"] if s["status"] in COUNTED else "other"] += 1
                counts["dated"] += 1 if s["collection_date"] else 0
            rows.append(
                [
                    run["name"],
                    run["run_date"] or "?",
                    run["mode"] or "?",
                    run["barcode_iso"] or run["barcode_version"] or "unknown",
                    str(sum(counts[k] for k in ("ok", "low_confidence", "other"))),
                    str(counts["ok"]),
                    str(counts["low_confidence"]),
                    str(counts["other"]),
                    str(counts["dated"]),
                ]
            )
        print_table(
            rows, ["run", "date", "mode", "barcodes", "samples", "ok", "low_conf", "other", "dated"]
        )
        total = conn.execute("SELECT COUNT(*) FROM samples").fetchone()[0]
        distinct = len(current_samples(conn))
        print(
            f"\n{len(runs)} run(s), {total} sample row(s), {distinct} distinct sample(s) (latest run each)"
        )
    finally:
        conn.close()
    return 0


# =========================================================================================================
# alerts
# =========================================================================================================


def current_samples(conn):
    """The latest run's row per sample: {(sample, collection_date): row dict with run fields}. Latest = newest
    barcode version, then newest run start, then run_id. The row is chosen first and filtered afterwards, so a
    newer run that failed a sample hides its older result."""
    chosen = {}
    query = (
        "SELECT s.*, r.run_date AS run_date, r.start AS start, r.barcode_iso AS barcode_iso "
        "FROM samples s JOIN runs r USING (run_id)"
    )
    for row in conn.execute(query):
        row = dict(row)
        key = (row["sample"], row["collection_date"])
        rank = (row["barcode_iso"] or "", row["start"] or "", row["run_id"])
        if key not in chosen or rank > chosen[key][0]:
            chosen[key] = (rank, row)
    return {k: v[1] for k, v in chosen.items()}


def period_of(date_text):
    iso = datetime.date.fromisoformat(date_text).isocalendar()
    return f"{iso[0]}-W{iso[1]:02d}"


class Obs:
    """The dated, counted samples of the store (latest run per sample) with their per-sample values."""

    def __init__(self, conn, include_low_confidence):
        self.notes = []
        self.rows = []
        undated = fallback = 0
        for row in current_samples(conn).values():
            if row["status"] not in COUNTED or (row["sample_type"] or "").lower() == "negative":
                continue
            date = row["collection_date"] or row["run_date"]
            if not date:
                undated += 1
                continue
            fallback += 0 if row["collection_date"] else 1
            row["date"] = date
            row["site_label"] = row["site"] or NO_SITE
            row["period"] = period_of(date)
            self.rows.append(row)
        self.rows.sort(key=lambda r: (r["date"], r["sample"]))
        if undated:
            self.notes.append(
                f"{undated} counted sample(s) have neither a collection date nor a run date: not used"
            )
        if fallback:
            self.notes.append(
                f"{fallback} of {len(self.rows)} sample(s) have no collection date: the date of their run is used"
            )
        low = sum(1 for r in self.rows if r["status"] == "low_confidence")
        if low and not include_low_confidence:
            self.notes.append(
                f"{low} low_confidence sample(s) are left out of the emerging and growth rules"
            )
        self.by_key = {(r["run_id"], r["sample"]): r for r in self.rows}
        self.include_low_confidence = include_low_confidence
        self.conn = conn

    def table(self, name, column):
        """{(site, level item): {(run_id, sample): value}} from lineage_abundance / group_abundance."""
        data = {}
        for r in self.conn.execute(
            f"SELECT run_id, sample, {column} AS item, abundance FROM {name}"
        ):
            key = (r["run_id"], r["sample"])
            if key in self.by_key:
                data.setdefault(r["item"], {})[key] = r["abundance"] or 0.0
        return data


def sample_list(rows, limit=8):
    names = sorted({r["sample"] for r in rows})
    text = ",".join(names[:limit])
    return text + (f",(+{len(names) - limit} more)" if len(names) > limit else "")


def alert(rule, level, site, item, detail, rows):
    dates = sorted(r["date"] for r in rows)
    return {
        "rule": rule,
        "level": level,
        "site": site,
        "group_or_mutation": item,
        "detail": detail,
        "first_date": dates[0],
        "last_date": dates[-1],
        "samples": sample_list(rows),
    }


def abundance_levels(obs, levels):
    out = []
    for level, table, column in (
        ("group", "group_abundance", "grp"),
        ("lineage", "lineage_abundance", "lineage"),
    ):
        if level in levels:
            out.append((level, obs.table(table, column)))
    return out


def rules_abundance(obs, args):
    """Rules 1 (emerging) and 2 (sustained growth)."""
    alerts = []
    usable = [r for r in obs.rows if obs.include_low_confidence or r["status"] == "ok"]
    sites = {}
    for r in usable:
        sites.setdefault(r["site_label"], []).append(r)
    for level, data in abundance_levels(obs, args.levels):
        # only samples that have rows at this level count: the others have no data (they are not zeros)
        have = {k for values in data.values() for k in values}
        for site, rows in sorted(sites.items()):
            rows = [r for r in rows if (r["run_id"], r["sample"]) in have]
            if not rows:
                continue
            periods = sorted({r["period"] for r in rows})
            latest = periods[-1]
            now = [r for r in rows if r["period"] == latest]
            earlier = [r for r in rows if r["period"] != latest]
            dates = sorted({r["date"] for r in rows})
            for item, values in sorted(data.items()):
                if item in CATCH_ALL:
                    continue

                def value(r, values=values):
                    return values.get((r["run_id"], r["sample"]), 0.0)

                # rule 1
                if len(periods) - 1 >= args.min_history_periods:
                    mean_now = sum(value(r) for r in now) / len(now)
                    seen_before = max([value(r) for r in earlier] or [0.0])
                    if (
                        mean_now * 100 >= args.emerging_min_pct
                        and seen_before * 100 < args.emerging_max_earlier_pct
                    ):
                        hits = [r for r in now if value(r) * 100 >= args.emerging_max_earlier_pct]
                        top = max(value(r) for r in now) * 100
                        detail = (
                            f"mean {mean_now * 100:.1f}% in {latest} (max {top:.1f}%); at most "
                            f"{seen_before * 100:.2f}% in the {len(periods) - 1} earlier period(s) "
                            f"({periods[0]} to {periods[-2]})"
                        )
                        alerts.append(
                            alert(
                                "emerging",
                                level,
                                site,
                                item,
                                detail,
                                hits,
                            )
                        )
                # rule 2: the trailing strictly rising run over distinct dates (replicates averaged)
                by_date = {d: [value(r) for r in rows if r["date"] == d] for d in dates}
                series = [(d, sum(v) / len(v)) for d, v in sorted(by_date.items())]
                end = len(series) - 1
                begin = end
                while begin > 0 and series[begin - 1][1] < series[begin][1]:
                    begin -= 1
                count = end - begin + 1
                rise = (series[end][1] - series[begin][1]) * 100
                if count >= args.growth_min_samples and rise >= args.growth_min_rise_pts:
                    span = {d for d, _ in series[begin:]}
                    alerts.append(
                        alert(
                            "sustained_growth",
                            level,
                            site,
                            item,
                            "rose over {} consecutive dates: {} (+{:.1f} points)".format(
                                count,
                                " -> ".join(f"{v * 100:.1f}%" for _, v in series[begin:]),
                                rise,
                            ),
                            [r for r in rows if r["date"] in span],
                        )
                    )
    return alerts


def low_note(rows):
    low = {r["sample"] for r in rows if r["status"] == "low_confidence"}
    return f"; {len(low)} low_confidence sample(s)" if low else ""


def rule_unexplained(obs, args):
    """Rule 3: the same unexplained mutation in >= N samples on >= 2 dates within the last K weeks."""
    if not obs.rows:
        return []
    end = datetime.date.fromisoformat(max(r["date"] for r in obs.rows))
    start = (end - datetime.timedelta(weeks=args.unexplained_weeks)).isoformat()
    window = [r for r in obs.rows if r["date"] > start]
    sites = {}
    for r in window:
        sites.setdefault(r["site_label"], []).append(r)
    wanted = {(r["run_id"], r["sample"]) for r in window}
    hits = {}
    for u in obs.conn.execute("SELECT * FROM unexplained"):
        key = (u["run_id"], u["sample"])
        if key in wanted:
            hits.setdefault((obs.by_key[key]["site_label"], u["mutation"]), []).append(
                (obs.by_key[key], dict(u))
            )
    alerts = []
    for (site, mutation), items in sorted(hits.items()):
        rows = [r for r, _ in items]
        if len({r["sample"] for r in rows}) < args.unexplained_min_samples:
            continue
        if len({r["date"] for r in rows}) < 2:
            continue
        first = items[0][1]
        in_dataset = "yes" if any(u["in_dataset"] == "true" for _, u in items) else "no"
        mean = sum(u["observed_freq"] or 0.0 for _, u in items) / len(items)
        meaning = (
            "some lineage of the dataset carries it, but the fitted mixture does not predict it"
            if in_dataset == "yes"
            else "no lineage of the dataset carries it: a dataset blind spot or a novel change"
        )
        name = mutation + (f" ({first['aa']})" if first["aa"] else "")
        alerts.append(
            alert(
                "recurring_unexplained",
                "mutation",
                site,
                name,
                "in {} of {} samples on {} dates since {}; mean frequency {:.2f}; in_dataset={} ({}){}".format(
                    len({r["sample"] for r in rows}),
                    len(sites[site]),
                    len({r["date"] for r in rows}),
                    start,
                    mean,
                    in_dataset,
                    meaning,
                    low_note(rows),
                ),
                rows,
            )
        )
    return alerts


MUTATION_FORMS = (
    r"[ACGT]\d+[ACGT]",  # A23063T
    r"[ACGT]+\d+-",  # C21995-, CTT21991-
    r"[ACGT]\d+\+[ACGT]+",  # T22204+GAGCCAGAA
    r"[A-Za-z0-9_.]+:[A-Z*]\d+[A-Z*]",  # S:N501Y
)


def read_watchlist(path):
    """[(name, mutation)] from a `name,mutation` CSV (the --mutations_of_interest format)."""
    try:
        with open(str(path), encoding="utf-8-sig", newline="") as handle:
            reader = csv.DictReader(handle)
            if (
                not reader.fieldnames
                or "name" not in reader.fieldnames
                or "mutation" not in reader.fieldnames
            ):
                raise Problem(f"{path}: expected a header with columns name,mutation")
            entries = []
            for line, rec in enumerate(reader, start=2):
                name = (rec.get("name") or "").strip()
                mutation = (rec.get("mutation") or "").strip()
                if not name and not mutation:
                    continue
                if not any(re.fullmatch(form, mutation) for form in MUTATION_FORMS):
                    raise Problem(
                        f"{path}: line {line}: cannot parse mutation {mutation!r} "
                        "(expected e.g. A23063T, C21995-, T22204+GAGCCAGAA or S:N501Y)"
                    )
                entries.append((name or mutation, mutation))
    except OSError as err:
        raise Problem(f"cannot read {path}: {err}")
    return entries


def rule_watchlist(obs, args, entries):
    """Rule 4: a watchlist mutation at >= N % in a sample (mutations-of-interest table, else unexplained table).
    Returns (alerts, mutations never measured in the store)."""
    by_mutation = {m: n for n, m in entries}
    seen = {}  # mutation -> {(site, run_id, sample): freq}
    measured = set()
    for r in obs.conn.execute("SELECT * FROM moi"):
        key = (r["run_id"], r["sample"])
        if key in obs.by_key and r["mutation"] in by_mutation and r["alt_freq"] is not None:
            measured.add(r["mutation"])
            seen.setdefault(r["mutation"], {})[key] = max(
                r["alt_freq"], seen.get(r["mutation"], {}).get(key, 0.0)
            )
    for r in obs.conn.execute("SELECT * FROM unexplained"):
        key = (r["run_id"], r["sample"])
        if key not in obs.by_key or r["observed_freq"] is None:
            continue
        for form in (r["mutation"], r["aa"]):
            if form and form in by_mutation:
                measured.add(form)
                current = seen.setdefault(form, {}).get(key, 0.0)
                seen[form][key] = max(current, r["observed_freq"])
    alerts = []
    for mutation, per_sample in sorted(seen.items()):
        sites = {}
        for key, freq in per_sample.items():
            if freq * 100 >= args.watchlist_min_pct:
                sites.setdefault(obs.by_key[key]["site_label"], []).append((obs.by_key[key], freq))
        for site, items in sorted(sites.items()):
            rows = [r for r, _ in items]
            top = max(items, key=lambda x: x[1])
            alerts.append(
                alert(
                    "watchlist",
                    "mutation",
                    site,
                    f"{by_mutation[mutation]} ({mutation})"
                    if by_mutation[mutation] != mutation
                    else mutation,
                    "at >= {:g}% in {} sample(s); highest {:.1f}% in {}{}".format(
                        args.watchlist_min_pct,
                        len(rows),
                        top[1] * 100,
                        top[0]["sample"],
                        low_note(rows),
                    ),
                    rows,
                )
            )
    return alerts, sorted(set(by_mutation) - measured)


def compute_alerts(conn, args):
    obs = Obs(conn, args.include_low_confidence)
    alerts = []
    notes = list(obs.notes)
    if not obs.rows:
        notes.append("the store holds no dated sample with status ok or low_confidence")
    else:
        alerts += rules_abundance(obs, args)
        alerts += rule_unexplained(obs, args)
        if args.watchlist:
            found, unmeasured = rule_watchlist(obs, args, read_watchlist(args.watchlist))
            alerts += found
            if unmeasured:
                notes.append(
                    "watchlist mutation(s) not measured in any stored sample (not in a mutations-of-interest "
                    "table or the unexplained table, so absence of an alert says nothing): "
                    + ", ".join(unmeasured)
                )
    if args.since:
        alerts = [a for a in alerts if a["last_date"] >= args.since]
    alerts.sort(
        key=lambda a: (RULES.index(a["rule"]), a["site"], a["level"], a["group_or_mutation"])
    )
    return alerts, notes, obs


DISCLAIMER = (
    "These alerts are screening hints from a few samples, not epidemiological conclusions. Wastewater "
    "abundances are noisy, depend on coverage and on the barcode version, and a lineage call can change with "
    "the next dataset: confirm anything important before acting on it."
)


def write_tsv(path, alerts):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    with open(str(path), "w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(ALERT_COLUMNS)
        for a in alerts:
            writer.writerow([a[c] for c in ALERT_COLUMNS])


def markdown(alerts, notes, obs, args):
    out = ["# Wastewater alerts", "", "> " + DISCLAIMER, ""]
    sample_count = len(obs.rows)
    dates = sorted(r["date"] for r in obs.rows)
    out.append(
        f"- Samples used: {sample_count}"
        + (f" ({dates[0]} to {dates[-1]})" if dates else "")
        + f"; alerts: {len(alerts)}"
    )
    out.append(
        f"- Thresholds: emerging >= {args.emerging_min_pct:g}% now and < {args.emerging_max_earlier_pct:g}% "
        f"before (>= {args.min_history_periods} earlier periods); growth over >= {args.growth_min_samples} "
        f"dates and >= {args.growth_min_rise_pts:g} points; unexplained mutation in >= "
        f"{args.unexplained_min_samples} samples on >= 2 dates in {args.unexplained_weeks} weeks; "
        f"watchlist >= {args.watchlist_min_pct:g}%"
    )
    if args.since:
        out.append(f"- Only alerts with evidence on or after {args.since}")
    for note in notes:
        out.append(f"- Note: {note}")
    out.append("")
    if not alerts:
        out += ["No alerts.", ""]
    for rule in RULES:
        chosen = [a for a in alerts if a["rule"] == rule]
        if not chosen:
            continue
        out += [f"## {rule.replace('_', ' ').capitalize()} ({len(chosen)})", ""]
        out += [
            "| site | level | item | detail | first | last | samples |",
            "| --- | --- | --- | --- | --- | --- | --- |",
        ]
        for a in chosen:
            cells = [
                a["site"],
                a["level"],
                a["group_or_mutation"],
                a["detail"],
                a["first_date"],
                a["last_date"],
                a["samples"],
            ]
            out.append("| " + " | ".join(str(c).replace("|", "/") for c in cells) + " |")
        out.append("")
    return "\n".join(out)


def cmd_alerts(args):
    conn = open_store(args.store, create=False)
    try:
        alerts, notes, obs = compute_alerts(conn, args)
        if args.out:
            write_tsv(args.out, alerts)
        if args.markdown:
            Path(args.markdown).parent.mkdir(parents=True, exist_ok=True)
            Path(args.markdown).write_text(markdown(alerts, notes, obs, args), encoding="utf-8")
    finally:
        conn.close()
    print("Screening hints, not epidemiological conclusions.")
    for note in notes:
        print("note: " + note)
    if not alerts:
        print("no alerts")
    for a in alerts:
        print(
            "  {:<22} {:<8} {:<12} {:<34} {}..{}  {}".format(
                a["rule"],
                a["level"],
                a["site"],
                a["group_or_mutation"][:34],
                a["first_date"],
                a["last_date"],
                a["detail"],
            )
        )
    if args.out:
        print(f"TSV: {args.out}")
    if args.markdown:
        print(f"Markdown: {args.markdown}")
    return 0


# =========================================================================================================
# command line
# =========================================================================================================


def parse_args(argv):
    parser = argparse.ArgumentParser(
        description="A cross-run store of wasteline results (SQLite) with trend and emergence alerts. "
        "Alerts are screening hints, not epidemiological conclusions.",
        epilog="Exit status: 0 ok, 1 error.",
    )
    sub = parser.add_subparsers(dest="command", required=True)

    add = sub.add_parser(
        "add", help="ingest finished results directories (idempotent per directory)"
    )
    add.add_argument("store", metavar="store.sqlite")
    add.add_argument(
        "results", nargs="+", metavar="results_dir", help="a results directory, or a folder of them"
    )

    runs = sub.add_parser("runs", help="list the runs in the store")
    runs.add_argument("store", metavar="store.sqlite")

    al = sub.add_parser("alerts", help="apply the alert rules to the store")
    al.add_argument("store", metavar="store.sqlite")
    al.add_argument(
        "--watchlist",
        metavar="mutations.csv",
        help="name,mutation CSV (as --mutations_of_interest)",
    )
    al.add_argument(
        "--since",
        metavar="YYYY-MM-DD",
        help="only alerts whose latest evidence is on or after this date",
    )
    al.add_argument(
        "--out", metavar="TSV", help="write the alerts as a TSV (header only when there are none)"
    )
    al.add_argument("--markdown", metavar="MD", help="write a short human-readable report")
    al.add_argument(
        "--include-low-confidence",
        action="store_true",
        help="count low_confidence samples in the emerging and growth rules",
    )
    al.add_argument(
        "--levels",
        default="group,lineage",
        help="abundance levels for rules 1-2: group, lineage or both (default group,lineage)",
    )
    al.add_argument(
        "--emerging-min-pct",
        type=float,
        default=5.0,
        help="rule 1: abundance now at least this %% (default 5)",
    )
    al.add_argument(
        "--emerging-max-earlier-pct",
        type=float,
        default=1.0,
        help="rule 1: abundance in every earlier period below this %% (default 1)",
    )
    al.add_argument(
        "--min-history-periods",
        type=int,
        default=2,
        help="rule 1: earlier ISO weeks needed at the site (default 2)",
    )
    al.add_argument(
        "--growth-min-samples",
        type=int,
        default=3,
        help="rule 2: consecutive dates that rise (default 3)",
    )
    al.add_argument(
        "--growth-min-rise-pts",
        type=float,
        default=10.0,
        help="rule 2: total rise in percentage points (default 10)",
    )
    al.add_argument(
        "--unexplained-min-samples",
        type=int,
        default=3,
        help="rule 3: samples with the mutation (default 3; also needs 2 dates)",
    )
    al.add_argument(
        "--unexplained-weeks",
        type=int,
        default=8,
        help="rule 3: window in weeks ending at the latest sample date (default 8)",
    )
    al.add_argument(
        "--watchlist-min-pct",
        type=float,
        default=5.0,
        help="rule 4: frequency at least this %% (default 5)",
    )

    args = parser.parse_args(argv)
    if args.command == "alerts":
        levels = [x.strip() for x in args.levels.split(",") if x.strip()]
        if not levels or any(x not in ("group", "lineage") for x in levels):
            parser.error("--levels takes group, lineage or group,lineage")
        args.levels = levels
        if args.since and not iso_date(args.since):
            parser.error("--since must be YYYY-MM-DD")
        for name in (
            "emerging_min_pct",
            "emerging_max_earlier_pct",
            "growth_min_rise_pts",
            "watchlist_min_pct",
        ):
            if getattr(args, name) < 0:
                parser.error("--" + name.replace("_", "-") + " must be >= 0")
        if (
            args.growth_min_samples < 2
            or args.unexplained_min_samples < 1
            or args.min_history_periods < 1
        ):
            parser.error(
                "--growth-min-samples must be >= 2; --unexplained-min-samples and --min-history-periods >= 1"
            )
    return args


def main(argv=None):
    args = parse_args(sys.argv[1:] if argv is None else argv)
    try:
        return {"add": cmd_add, "runs": cmd_runs, "alerts": cmd_alerts}[args.command](args)
    except (Problem, OSError, sqlite3.Error) as err:
        print(f"error: {err}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())
