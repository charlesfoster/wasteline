#!/usr/bin/env python3
"""Write the per-sample `<sample>.demix_status.json` record.

status is one of ok | failed | skipped_low_coverage. A non-ok status must carry a reason.
"""

import argparse
import json
import sys

STATUSES = ("ok", "failed", "skipped_low_coverage")
MAX_REASON_CHARS = 500


def tail_reason(path, max_lines=5):
    """Last non-empty lines of a log file, joined into a one-line reason."""
    try:
        with open(path, errors="replace") as handle:
            lines = [line.strip() for line in handle if line.strip()]
    except OSError as err:
        return f"could not read log {path}: {err}"
    return " | ".join(lines[-max_lines:])


def build_record(sample, status, reason="", stats=None, barcode=None):
    """Status record. `barcode` (barcode_version.json) is mandatory: every record, whatever its
    status, states which Freyja barcode version the run used."""
    barcode = barcode or {}
    if not barcode.get("barcode_version"):
        raise ValueError("barcode_version is missing: refusing to write a status without it")
    if status not in STATUSES:
        raise ValueError(f"status must be one of {', '.join(STATUSES)}, got {status!r}")
    reason = " ".join(str(reason).split())
    if status != "ok" and not reason:
        raise ValueError(f"a reason is required when status is {status!r}")
    if len(reason) > MAX_REASON_CHARS:
        reason = reason[: MAX_REASON_CHARS - 3] + "..."
    stats = stats or {}
    return {
        "sample": sample,
        "status": status,
        "reason": reason,
        "genome_coverage_pct": stats.get("pct_genome_ge_covcut"),
        "covcut": stats.get("covcut"),
        "min_genome_coverage": stats.get("min_genome_coverage"),
        "barcode_version": barcode["barcode_version"],
        "barcode_version_iso": barcode.get("barcode_version_iso"),
    }


def parse_args(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--sample", required=True)
    parser.add_argument("--status", required=True, choices=STATUSES)
    parser.add_argument("--reason", default="")
    parser.add_argument(
        "--reason-file", default=None, help="use the tail of this log as the reason"
    )
    parser.add_argument(
        "--stats", default=None, help="depth_stats.py JSON to copy coverage fields from"
    )
    parser.add_argument(
        "--barcode-info", required=True, help="barcode_version.json from freyja_dataset_info.py"
    )
    parser.add_argument("--output", required=True)
    return parser.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    reason = args.reason
    if args.reason_file:
        tail = tail_reason(args.reason_file)
        reason = f"{reason}: {tail}" if reason and tail else (reason or tail)
    stats = None
    if args.stats:
        try:
            with open(args.stats) as handle:
                stats = json.load(handle)
        except (OSError, ValueError):
            stats = None
    try:
        with open(args.barcode_info) as handle:
            barcode = json.load(handle)
        record = build_record(args.sample, args.status, reason, stats, barcode)
    except (OSError, ValueError) as err:
        print(f"demix_status.py: {err}", file=sys.stderr)
        return 2
    with open(args.output, "w") as handle:
        json.dump(record, handle, indent=2)
        handle.write("\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
