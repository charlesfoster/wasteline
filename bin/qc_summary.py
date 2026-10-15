#!/usr/bin/env python3
"""Per-sample coverage QC row (qc_summary.tsv) from a Freyja depths file, samtools flagstat and fastp JSON.

Two modes:

  qc_summary.py --sample S --sample-type sample --depths S.depths [--flagstat S.flagstat]
                [--fastp-json S.fastp.json] [--amplicons S.amplicon_coverage.tsv] -o S.qc.tsv
      writes one row (with header).

  qc_summary.py --concat a.tsv b.tsv ... -o qc_summary.tsv
      merges rows, sorted by sample; headers must all equal the contract header.

Depth statistics are over every genome position in the depths file (zeros included); SD is the population
SD. Missing optional inputs give `NA`. `raw_reads`/`trimmed_reads` are fastp's read counts (reads, not
pairs: a pair counts as 2). `mapped_reads` is `primary mapped` (or `mapped` for old samtools) from flagstat
(also records, so a pair counts as 2). `mapped_pct` = mapped_reads / trimmed_reads * 100 (fallback: raw
reads, then the flagstat primary total); `NA` when no denominator is available.

QC status:
  sample_type=sample
    FAIL  pct genome >= --min-depth is below --fail-below-pct (default 50), or no reads at all
    WARN  below --pass-min-pct (default 80), or amplicons_dropped > --max-dropped (default 10)
    PASS  otherwise
    reasons: low_coverage, amplicon_dropout, no_reads
  sample_type=negative
    FAIL  mapped_reads > --ntc-max-mapped-reads (1000) or pct_genome_ge_1x > --ntc-max-genome-coverage-pct (10)
    PASS  otherwise (never low_coverage)
    reasons: ntc_reads
  Unreadable/empty depths: samples FAIL with `no_depths`; negatives are judged on whatever is available
  (flagstat only) and are NA/`no_depths` when nothing is.
"""

import argparse
import contextlib
import json
import re
import statistics
import sys

COLUMNS = [
    "sample",
    "sample_type",
    "collection_date",
    "site",
    "raw_reads",
    "trimmed_reads",
    "mapped_reads",
    "mapped_pct",
    "mean_depth",
    "depth_sd",
    "median_depth",
    "pct_genome_ge_1x",
    "pct_genome_ge_10x",
    "pct_genome_ge_min_depth",
    "amplicons_total",
    "amplicons_dropped",
    "qc_status",
    "qc_reasons",
    "demix_status",
    "demix_reason",
    "consensus_length",
    "consensus_n_pct",
    "barcode_version",
]
NA = "NA"


def read_depth_values(path):
    """Return the list of depths from a Freyja depths file; ValueError if unreadable or empty."""
    values = []
    try:
        with open(path) as handle:
            for line_number, line in enumerate(handle, start=1):
                line = line.rstrip("\r\n")
                if not line.strip():
                    continue
                fields = line.split("\t")
                if len(fields) < 4:
                    raise ValueError(
                        f"{path}: line {line_number} has {len(fields)} columns, expected at least 4"
                    )
                try:
                    values.append(int(float(fields[3])))
                except ValueError as exc:
                    raise ValueError(f"{path}: line {line_number}: non-numeric depth") from exc
    except OSError as exc:
        raise ValueError(f"cannot read depths {path}: {exc}") from exc
    if not values:
        raise ValueError(f"{path}: no positions")
    return values


def depth_stats(values, min_depth):
    n = len(values)

    def pct(threshold):
        return 100.0 * sum(1 for v in values if v >= threshold) / n

    return {
        "positions": n,
        "mean": statistics.fmean(values),
        "sd": statistics.pstdev(values),
        "median": statistics.median(values),
        "ge1": pct(1),
        "ge10": pct(10),
        "ge_min": pct(min_depth),
        "all_zero": not any(values),
    }


def parse_flagstat(path):
    """Return dict with `mapped` (primary mapped when present) and `total` (primary total when present)."""
    result = {"mapped": None, "total": None}
    pat = re.compile(r"^(\d+)\s*\+\s*(\d+)\s+(.*)$")
    mapped_any = None
    with open(path) as handle:
        for line in handle:
            m = pat.match(line.strip())
            if not m:
                continue
            count = int(m.group(1)) + int(m.group(2))
            label = m.group(3)
            if label.startswith("primary mapped"):
                result["mapped"] = count
            elif label.startswith("mapped"):
                mapped_any = count
            elif label.startswith("primary") and not label.startswith("primary duplicates"):
                result["total"] = count
            elif label.startswith("in total") and result["total"] is None:
                result["total_all"] = count
    if result["mapped"] is None:
        result["mapped"] = mapped_any
    if result["total"] is None:
        result["total"] = result.get("total_all")
    if result["mapped"] is None:
        raise ValueError(f"{path}: no 'mapped' line found")
    return result


def parse_fastp(path):
    with open(path) as handle:
        data = json.load(handle)
    summary = data["summary"]
    return int(summary["before_filtering"]["total_reads"]), int(
        summary["after_filtering"]["total_reads"]
    )


def parse_amplicons(path, sample):
    """Return (total, dropped) for the sample from an amplicon_coverage.tsv."""
    total = dropped = 0
    with open(path) as handle:
        header = handle.readline().rstrip("\r\n").split("\t")
        if "dropped" not in header:
            raise ValueError(f"{path}: no 'dropped' column")
        i_drop = header.index("dropped")
        i_sample = header.index("sample") if "sample" in header else None
        for line in handle:
            fields = line.rstrip("\r\n").split("\t")
            if len(fields) <= i_drop:
                continue
            if i_sample is not None and fields[i_sample] != sample:
                continue
            total += 1
            dropped += fields[i_drop].lower() == "true"
    return total, dropped


def fmt(value, digits=2):
    return NA if value is None else f"{value:.{digits}f}"


def fmt_int(value):
    return NA if value is None else str(int(value))


def assess(sample_type, stats, mapped_reads, amplicons_dropped, args):
    """Return (qc_status, [reasons])."""
    reasons = []
    if sample_type == "negative":
        if stats is None and mapped_reads is None:
            return "NA", ["no_depths"]
        if mapped_reads is not None and mapped_reads > args.ntc_max_mapped_reads:
            reasons.append("ntc_reads")
        if (
            stats is not None
            and stats["ge1"] > args.ntc_max_genome_coverage_pct
            and "ntc_reads" not in reasons
        ):
            reasons.append("ntc_reads")
        return ("FAIL" if reasons else "PASS"), reasons

    if stats is None:
        return "FAIL", ["no_depths"]
    if mapped_reads == 0 or stats["all_zero"]:
        return "FAIL", ["no_reads"]
    status = "PASS"
    cov = stats["ge_status"]
    if cov < args.fail_below_pct:
        status = "FAIL"
        reasons.append("low_coverage")
    elif cov < args.pass_min_pct:
        status = "WARN"
        reasons.append("low_coverage")
    if amplicons_dropped is not None and amplicons_dropped > args.max_dropped:
        reasons.append("amplicon_dropout")
        if status == "PASS":
            status = "WARN"
    return status, reasons


def build_row(args):
    notes = []

    stats = None
    try:
        values = read_depth_values(args.depths)
        stats = depth_stats(values, args.consensus_min_depth)
        stats["ge_status"] = 100.0 * sum(1 for v in values if v >= args.min_depth) / len(values)
    except ValueError as exc:
        notes.append(str(exc))

    mapped_reads = flag_total = None
    if args.flagstat:
        try:
            fs = parse_flagstat(args.flagstat)
            mapped_reads, flag_total = fs["mapped"], fs["total"]
        except (OSError, ValueError) as exc:
            notes.append(f"flagstat: {exc}")

    raw = trimmed = None
    if args.fastp_json:
        try:
            raw, trimmed = parse_fastp(args.fastp_json)
        except (OSError, ValueError, KeyError, TypeError) as exc:
            notes.append(f"fastp: {exc}")

    amp_total = amp_dropped = None
    if args.amplicons:
        try:
            amp_total, amp_dropped = parse_amplicons(args.amplicons, args.sample)
        except (OSError, ValueError) as exc:
            notes.append(f"amplicons: {exc}")

    mapped_pct = None
    if mapped_reads is not None:
        denom = next((d for d in (trimmed, raw, flag_total) if d), None)
        if denom:
            mapped_pct = 100.0 * mapped_reads / denom

    status, reasons = assess(args.sample_type, stats, mapped_reads, amp_dropped, args)

    row = {c: NA for c in COLUMNS}
    row.update(
        {
            "sample": args.sample,
            "sample_type": args.sample_type,
            "collection_date": args.collection_date or NA,
            "site": args.site or NA,
            "raw_reads": fmt_int(raw),
            "trimmed_reads": fmt_int(trimmed),
            "mapped_reads": fmt_int(mapped_reads),
            "mapped_pct": fmt(mapped_pct),
            "amplicons_total": fmt_int(amp_total),
            "amplicons_dropped": fmt_int(amp_dropped),
            "qc_status": status,
            "qc_reasons": ";".join(reasons),
        }
    )
    if stats is not None:
        row.update(
            {
                "mean_depth": fmt(stats["mean"]),
                "depth_sd": fmt(stats["sd"]),
                "median_depth": fmt(stats["median"]),
                "pct_genome_ge_1x": fmt(stats["ge1"]),
                "pct_genome_ge_10x": fmt(stats["ge10"]),
                "pct_genome_ge_min_depth": fmt(stats["ge_min"]),
            }
        )
    for note in notes:
        print(f"warning: {args.sample}: {note}", file=sys.stderr)
    return row


def write_rows(path, rows):
    handle = sys.stdout if path == "-" else open(path, "w")  # noqa: SIM115
    with contextlib.nullcontext(handle) if path == "-" else handle as out:
        out.write("\t".join(COLUMNS) + "\n")
        for row in rows:
            out.write("\t".join(row[c] for c in COLUMNS) + "\n")


def concat(paths):
    rows = []
    for path in paths:
        with open(path) as handle:
            lines = [line.rstrip("\r\n") for line in handle if line.strip()]
        if not lines:
            raise ValueError(f"{path}: empty file")
        header = lines[0].split("\t")
        if header != COLUMNS:
            raise ValueError(f"{path}: header does not match the qc_summary contract")
        for line in lines[1:]:
            fields = line.split("\t")
            if len(fields) != len(COLUMNS):
                raise ValueError(f"{path}: row has {len(fields)} fields, expected {len(COLUMNS)}")
            rows.append(dict(zip(COLUMNS, fields)))
    samples = [r["sample"] for r in rows]
    dupes = sorted({s for s in samples if samples.count(s) > 1})
    if dupes:
        raise ValueError(f"duplicate samples: {', '.join(dupes)}")
    rows.sort(key=lambda r: r["sample"])
    return rows


def main(argv=None):
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    p.add_argument(
        "--concat", nargs="*", metavar="TSV", help="merge per-sample rows instead of computing one"
    )
    p.add_argument("--sample")
    p.add_argument("--sample-type", choices=["sample", "negative"], default="sample")
    p.add_argument("--collection-date", default=None)
    p.add_argument("--site", default=None)
    p.add_argument("--depths")
    p.add_argument("--flagstat")
    p.add_argument("--fastp-json")
    p.add_argument("--amplicons", help="amplicon_coverage.tsv (counts amplicons_total/dropped)")
    p.add_argument(
        "--min-depth",
        type=int,
        default=10,
        help="depth for the coverage PASS/WARN/FAIL test (default 10)",
    )
    p.add_argument(
        "--consensus-min-depth",
        type=int,
        default=10,
        help="depth for pct_genome_ge_min_depth (default 10)",
    )
    p.add_argument(
        "--pass-min-pct",
        type=float,
        default=80.0,
        help="PASS needs >= this %% genome at --min-depth",
    )
    p.add_argument(
        "--fail-below-pct",
        type=float,
        default=50.0,
        help="FAIL below this %% genome at --min-depth",
    )
    p.add_argument(
        "--max-dropped", type=int, default=10, help="WARN when more amplicons than this are dropped"
    )
    p.add_argument("--ntc-max-mapped-reads", type=int, default=1000)
    p.add_argument("--ntc-max-genome-coverage-pct", type=float, default=10.0)
    p.add_argument("-o", "--output", default="-")
    args = p.parse_args(argv)

    if args.concat is not None:
        try:
            rows = concat(args.concat)
        except (OSError, ValueError) as exc:
            print(f"error: {exc}", file=sys.stderr)
            return 2
        write_rows(args.output, rows)
        return 0

    if not args.sample or not args.depths:
        p.error("--sample and --depths are required (or use --concat)")
    write_rows(args.output, [build_row(args)])
    return 0


if __name__ == "__main__":
    sys.exit(main())
