#!/usr/bin/env python3
"""Drop the samples of a given demix status from a `freyja aggregate` table.

Used to keep `low_confidence` samples (PLAN D18) out of the time series, growth rates and the
dashboard while `freyja_aggregated.tsv` itself stays complete. The first column of the aggregated
table is the sample label; the header, the column order and every kept row are copied as they are.
Samples that the status table does not list are kept. When nothing is left only the header is
written.
"""

import argparse
import csv
import sys


def excluded_samples(status_path, statuses):
    """Sample IDs whose status (demix_status.tsv) is one of `statuses`."""
    with open(status_path, newline="") as handle:
        return {
            row["sample"]
            for row in csv.DictReader(handle, delimiter="\t")
            if row.get("status") in statuses
        }


def filter_lines(lines, excluded):
    """Yield the header and the rows whose label is not excluded."""
    for number, line in enumerate(lines):
        if number == 0 or line.partition("\t")[0] not in excluded:
            yield line


def parse_args(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("aggregated", help="freyja_aggregated.tsv")
    parser.add_argument("--demix-status", required=True, help="demix_status.tsv")
    parser.add_argument(
        "--exclude-status",
        action="append",
        required=True,
        help="status whose samples are dropped; repeatable",
    )
    parser.add_argument("--output", required=True)
    return parser.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    try:
        excluded = excluded_samples(args.demix_status, set(args.exclude_status))
        with open(args.aggregated, newline="") as handle:
            lines = handle.read().splitlines(keepends=True)
    except OSError as err:
        print(f"filter_aggregated.py: {err}", file=sys.stderr)
        return 2
    if not lines:
        print(f"filter_aggregated.py: {args.aggregated} is empty", file=sys.stderr)
        return 2
    with open(args.output, "w", newline="") as handle:
        handle.writelines(filter_lines(lines, excluded))
    return 0


if __name__ == "__main__":
    sys.exit(main())
