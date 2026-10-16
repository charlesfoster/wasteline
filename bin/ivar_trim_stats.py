#!/usr/bin/env python3
"""Primer-trimming effectiveness from the text `ivar trim` prints (iVar 1.2 - 1.4.x).

  ivar_trim_stats.py --log S.ivar.log --sample S -o S.ivar_trim_stats.tsv

`ivar trim` summarises its run in lines such as

  Trimmed primers from 21.97% (3515) of reads.
  0% (0) of reads were quality trimmed below the minimum length of 75 bp and were not written to file.
  78.03% (12483) of reads started outside of primer regions. Since the -e flag was given, these reads were written to file.

(the last line reads "... of reads that started outside of primer regions were not written to file" without
`-e`, and "... Since there were no primers found in BED file ..." for an empty BED). These go to STDERR; the log
passed here must therefore contain stderr (the pipeline merges it into the published `*.ivar.log`). Other lines
(the per-primer table, BED warnings) are ignored.

Output: one row with columns

  sample, reads_total, primer_trimmed_reads, primer_trimmed_pct, no_primer_reads, no_primer_pct,
  quality_trimmed_reads, quality_trimmed_pct, status

`primer_trimmed_pct` is the share of mapped reads that started inside a primer (a wrong primer scheme gives a
low value). `reads_total` = primer_trimmed_reads + no_primer_reads (ivar's own denominator). Quality-trimmed
reads overlap those two groups. `status` is `ok`, `missing` (no/unreadable log), `empty` or `unparsed` (text
present but no summary line; for example garbled or a very different ivar version). Anything not found is `NA`.
Never fails on bad input: exit 0 with a note on stderr (exit 2 only for bad usage).
"""

import argparse
import math
import re
import sys

NA = "NA"
COLUMNS = [
    "sample",
    "reads_total",
    "primer_trimmed_reads",
    "primer_trimmed_pct",
    "no_primer_reads",
    "no_primer_pct",
    "quality_trimmed_reads",
    "quality_trimmed_pct",
    "status",
]

PCT = r"(-?[0-9.]+|-?nan|-?inf)"
TRIMMED = re.compile(r"Trimmed primers from\s+" + PCT + r"\s*%\s*\((\d+)\)\s+of reads")
QUALITY = re.compile(PCT + r"\s*%\s*\((\d+)\)\s+of reads were quality trimmed")
NO_PRIMER = re.compile(
    PCT + r"\s*%\s*\((\d+)\)\s+of reads (?:that )?started outside of primer regions"
)


def to_pct(text):
    """Percent text -> float, or None when not a finite number (e.g. ivar prints nan for 0 reads)."""
    try:
        value = float(text)
    except ValueError:
        return None
    return value if math.isfinite(value) else None


def parse_log(text):
    """Return a dict with keys trimmed, no_primer, quality; each (pct or None, count) or absent."""
    found = {}
    for line in text.splitlines():
        for key, pattern in (
            ("trimmed", TRIMMED),
            ("quality", QUALITY),
            ("no_primer", NO_PRIMER),
        ):
            m = pattern.search(line)
            if m and key not in found:
                found[key] = (to_pct(m.group(1)), int(m.group(2)))
    return found


def fmt_pct(value):
    return NA if value is None else f"{value:.2f}"


def fmt_count(entry):
    return NA if entry is None else str(entry[1])


def build_row(sample, log_path):
    row = {c: NA for c in COLUMNS}
    row["sample"] = sample
    row["status"] = "missing"
    if not log_path:
        print(f"note: {sample}: no ivar log given", file=sys.stderr)
        return row
    try:
        with open(log_path, errors="replace") as handle:
            text = handle.read()
    except OSError as exc:
        print(f"note: {sample}: cannot read ivar log: {exc}", file=sys.stderr)
        return row
    if not text.strip():
        row["status"] = "empty"
        print(
            f"note: {sample}: ivar log {log_path} is empty (ivar prints its summary to stderr)",
            file=sys.stderr,
        )
        return row

    found = parse_log(text)
    if not found:
        row["status"] = "unparsed"
        print(f"note: {sample}: no ivar trim summary line in {log_path}", file=sys.stderr)
        return row

    trimmed, no_primer, quality = found.get("trimmed"), found.get("no_primer"), found.get("quality")
    total = None
    if trimmed is not None and no_primer is not None:
        total = trimmed[1] + no_primer[1]
    pct = {}
    for key, entry in (("trimmed", trimmed), ("no_primer", no_primer), ("quality", quality)):
        if entry is None:
            pct[key] = None
        elif total:
            # recomputed from counts (ivar prints 2 dp); quality-trimmed reads are a subset of the total
            pct[key] = 100.0 * entry[1] / total
        else:
            pct[key] = entry[0]
    row.update(
        {
            "reads_total": NA if total is None else str(total),
            "primer_trimmed_reads": fmt_count(trimmed),
            "primer_trimmed_pct": fmt_pct(pct["trimmed"]),
            "no_primer_reads": fmt_count(no_primer),
            "no_primer_pct": fmt_pct(pct["no_primer"]),
            "quality_trimmed_reads": fmt_count(quality),
            "quality_trimmed_pct": fmt_pct(pct["quality"]),
            "status": "ok" if trimmed is not None else "unparsed",
        }
    )
    if trimmed is None:
        print(f"note: {sample}: no 'Trimmed primers from' line in {log_path}", file=sys.stderr)
    return row


def main(argv=None):
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    p.add_argument("--log", help="ivar trim log containing ivar's stderr (missing/empty -> NA)")
    p.add_argument("--sample", required=True)
    p.add_argument("-o", "--output", default="-")
    args = p.parse_args(argv)

    row = build_row(args.sample, args.log)
    text = "\t".join(COLUMNS) + "\n" + "\t".join(row[c] for c in COLUMNS) + "\n"
    if args.output == "-":
        sys.stdout.write(text)
    else:
        with open(args.output, "w") as out:
            out.write(text)
    return 0


if __name__ == "__main__":
    sys.exit(main())
