#!/usr/bin/env python3
"""Genome coverage statistics from a `freyja variants` depths file.

The depths file is the tab-separated output of `samtools mpileup -aa` reduced to its first four
columns (contig, position, reference base, depth) with no header. `-aa` reports every position of
the reference, so the number of rows is the genome length.

Exit codes: 0 = coverage gate passed (or no gate requested), 10 = coverage below the gate,
2 = unreadable / empty / malformed input.
"""

import argparse
import json
import sys

EXIT_OK = 0
EXIT_BELOW_GATE = 10
EXIT_ERROR = 2


def compute_stats(path, covcut):
    """Return coverage statistics for a depths file.

    Raises ValueError for malformed rows or a file without any positions.
    """
    positions = 0
    positions_ge_covcut = 0
    depth_sum = 0
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
                depth = int(float(fields[3]))
            except ValueError:
                raise ValueError(
                    f"{path}: line {line_number} has a non-numeric depth {fields[3]!r}"
                ) from None
            positions += 1
            depth_sum += depth
            if depth >= covcut:
                positions_ge_covcut += 1
    if positions == 0:
        raise ValueError(f"{path}: no positions found in depths file")
    return {
        "positions": positions,
        "positions_ge_covcut": positions_ge_covcut,
        "mean_depth": round(depth_sum / positions, 3),
        "pct_genome_ge_covcut": round(100.0 * positions_ge_covcut / positions, 4),
        "covcut": covcut,
    }


def apply_gate(stats, min_genome_coverage):
    """Add the gate fields to `stats` and return it."""
    stats["min_genome_coverage"] = min_genome_coverage
    stats["passes_gate"] = (
        True
        if min_genome_coverage is None
        else stats["pct_genome_ge_covcut"] >= min_genome_coverage
    )
    return stats


def parse_args(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("depths", help="depths file (contig, pos, base, depth; no header)")
    parser.add_argument(
        "--covcut", type=int, default=10, help="minimum depth counted as covered (default: 10)"
    )
    parser.add_argument(
        "--min-genome-coverage",
        type=float,
        default=None,
        help="gate: exit with code 10 if %% genome >= covcut is below this value",
    )
    parser.add_argument("--sample", default=None, help="sample ID to record in the JSON")
    parser.add_argument("--output", default=None, help="write JSON here (default: stdout)")
    return parser.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    try:
        stats = compute_stats(args.depths, args.covcut)
    except (OSError, ValueError) as err:
        print(f"depth_stats.py: {err}", file=sys.stderr)
        return EXIT_ERROR
    apply_gate(stats, args.min_genome_coverage)
    if args.sample is not None:
        stats = {"sample": args.sample, **stats}
    text = json.dumps(stats, indent=2) + "\n"
    if args.output:
        with open(args.output, "w") as handle:
            handle.write(text)
    else:
        sys.stdout.write(text)
    return EXIT_OK if stats["passes_gate"] else EXIT_BELOW_GATE


if __name__ == "__main__":
    sys.exit(main())
