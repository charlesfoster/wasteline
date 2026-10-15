#!/usr/bin/env python3
"""Low-depth mask BED from a Freyja depths file.

Writes merged BED intervals (0-based, half-open) of every position with depth < --min-depth. The contig
name is taken from the depths file (several contigs are handled independently, in file order).

  --keep-bed BED   positions inside these intervals are removed from the mask (e.g. the positions of
                   confidently called variants, so that bcftools consensus can still apply them; this
                   replaces the old pipeline's `bedtools subtract`). Only the first three BED columns are
                   used; contig names must match the depths file exactly.
  --length N       reference length. Positions of the contig not listed in the depths file (including a
                   tail beyond the last listed position) have depth 0 and are masked. Without it the contig
                   is taken to end at its last listed position (unlisted gaps before that are masked too).
                   Only valid for a single-contig depths file.

The depths file has four tab-separated columns, no header: contig, 1-based position, ref base, depth.
An empty or unreadable depths file is an error (exit 2), never an empty mask.
"""

import argparse
import contextlib
import sys
from collections import OrderedDict


def read_depths(path):
    """Return OrderedDict contig -> dict(position -> depth)."""
    contigs = OrderedDict()
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
                pos = int(fields[1])
                depth = int(float(fields[3]))
            except ValueError as exc:
                raise ValueError(f"{path}: line {line_number}: non-numeric position/depth") from exc
            contigs.setdefault(fields[0], {})[pos] = depth
    if not contigs:
        raise ValueError(f"{path}: no positions")
    return contigs


def read_keep_bed(path):
    """Return dict contig -> sorted merged list of (start, end) 0-based half-open."""
    raw = {}
    with open(path) as handle:
        for line_number, line in enumerate(handle, start=1):
            line = line.rstrip("\r\n")
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            fields = line.split("\t")
            if len(fields) < 3:
                raise ValueError(f"{path}: line {line_number} has fewer than 3 columns")
            try:
                start, end = int(fields[1]), int(fields[2])
            except ValueError as exc:
                raise ValueError(f"{path}: line {line_number}: non-integer coordinates") from exc
            raw.setdefault(fields[0], []).append((start, end))
    return {contig: merge(intervals) for contig, intervals in raw.items()}


def merge(intervals):
    """Merge overlapping or adjacent half-open intervals."""
    merged = []
    for start, end in sorted(i for i in intervals if i[1] > i[0]):
        if merged and start <= merged[-1][1]:
            if end > merged[-1][1]:
                merged[-1] = (merged[-1][0], end)
        else:
            merged.append((start, end))
    return merged


def low_depth_intervals(depth_by_pos, min_depth, length=None):
    """Merged 0-based half-open intervals of positions (1..length) with depth < min_depth."""
    last = max(depth_by_pos)
    total = max(last, length) if length else last
    intervals = []
    run_start = None
    for pos in range(1, total + 1):
        if depth_by_pos.get(pos, 0) < min_depth:
            if run_start is None:
                run_start = pos
        elif run_start is not None:
            intervals.append((run_start - 1, pos - 1))
            run_start = None
    if run_start is not None:
        intervals.append((run_start - 1, total))
    return intervals


def subtract(intervals, keep):
    """Remove the (merged, sorted) `keep` intervals from the (merged, sorted) `intervals`."""
    result = []
    for start, end in intervals:
        cursor = start
        for k_start, k_end in keep:
            if k_end <= cursor:
                continue
            if k_start >= end:
                break
            if k_start > cursor:
                result.append((cursor, k_start))
            cursor = max(cursor, k_end)
            if cursor >= end:
                break
        if cursor < end:
            result.append((cursor, end))
    return result


def build_mask(depths_by_contig, min_depth, length=None, keep=None):
    """Return list of (contig, start, end)."""
    if length is not None and len(depths_by_contig) > 1:
        raise ValueError("--length is only valid for a single-contig depths file")
    out = []
    for contig, depth_by_pos in depths_by_contig.items():
        intervals = low_depth_intervals(depth_by_pos, min_depth, length)
        if keep:
            intervals = subtract(intervals, keep.get(contig, []))
        out.extend((contig, s, e) for s, e in intervals)
    return out


def main(argv=None):
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    p.add_argument("--depths", required=True)
    p.add_argument(
        "--min-depth", type=int, default=15, help="mask positions with depth < this (default 15)"
    )
    p.add_argument("--keep-bed", help="BED of positions never masked (confident variant positions)")
    p.add_argument("--length", type=int, help="reference length; unlisted positions are depth 0")
    p.add_argument("-o", "--output", default="-")
    args = p.parse_args(argv)

    try:
        depths = read_depths(args.depths)
        keep = read_keep_bed(args.keep_bed) if args.keep_bed else None
        mask = build_mask(depths, args.min_depth, args.length, keep)
    except (OSError, ValueError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2
    handle = sys.stdout if args.output == "-" else open(args.output, "w")  # noqa: SIM115
    with contextlib.nullcontext(handle) if args.output == "-" else handle as out:
        for contig, start, end in mask:
            out.write(f"{contig}\t{start}\t{end}\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
