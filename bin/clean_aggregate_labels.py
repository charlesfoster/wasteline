#!/usr/bin/env python3
"""Strip filename extensions from the sample labels of a `freyja aggregate` table.

`freyja aggregate` labels each row with the *file name* that was given to `freyja demix` (for
wasteline `<sample>.variants.tsv`). This removes the first matching suffix from the label column
only. Nothing else is rewritten: underscores, dots and the remaining columns are left untouched.
Exits non-zero if two rows would end up with the same label. With `--sort` the rows are ordered by
label (Freyja's own row order follows the directory listing, which is not reproducible).
"""

import argparse
import sys


def clean_label(label, suffixes):
    """Remove the longest matching suffix (never the whole label)."""
    for suffix in sorted(suffixes, key=len, reverse=True):
        if suffix and label.endswith(suffix) and len(label) > len(suffix):
            return label[: -len(suffix)]
    return label


def clean_lines(lines, suffixes):
    """Yield the lines with a cleaned label column; the first line is the header."""
    seen = {}
    for number, line in enumerate(lines):
        if number == 0 or not line.strip():
            yield line
            continue
        label, sep, rest = line.partition("\t")
        cleaned = clean_label(label, suffixes)
        if cleaned in seen:
            raise ValueError(
                f"duplicate sample label {cleaned!r} (rows {seen[cleaned]} and {number + 1})"
            )
        seen[cleaned] = number + 1
        yield f"{cleaned}{sep}{rest}"


def parse_args(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("input")
    parser.add_argument("output")
    parser.add_argument(
        "--suffix",
        action="append",
        default=[],
        help="suffix to strip from labels; repeatable (the longest matching one is removed)",
    )
    parser.add_argument("--sort", action="store_true", help="sort the rows by (cleaned) label")
    return parser.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    with open(args.input, newline="") as handle:
        lines = handle.read().splitlines(keepends=True)
    try:
        out = list(clean_lines([entry.rstrip("\r\n") for entry in lines], args.suffix))
        if args.sort and out:
            out = [
                out[0],
                *sorted(
                    (row for row in out[1:] if row.strip()), key=lambda r: r.partition("\t")[0]
                ),
            ]
    except ValueError as err:
        print(f"clean_aggregate_labels.py: {err}", file=sys.stderr)
        return 1
    with open(args.output, "w") as handle:
        handle.write("\n".join(out) + ("\n" if out else ""))
    return 0


if __name__ == "__main__":
    sys.exit(main())
