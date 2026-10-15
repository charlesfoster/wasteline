#!/usr/bin/env python3
"""Rewrite the sample label of a `freyja demix` result file.

A `.demix` file is a two-column TSV whose first line is `<TAB><label>`; the label is the name of
the variants file that was demixed (`MIX_A.variants.tsv` for wasteline results,
`MIX_A.freyja_variants.tsv` or even a full path for legacy Snakemake results). `freyja aggregate`
uses that label as the row name. `aggregate` mode knows the real sample ID from the samplesheet, so
this sets the label to the sample ID exactly as given (no characters are rewritten) in a COPY of
the file; the original is never modified.
"""

import argparse
import sys


def relabel(text, sample):
    """Return `text` with the label on the first line replaced by `sample`."""
    first, sep, rest = text.partition("\n")
    head = first.rstrip("\r")
    if not head.startswith("\t") or "\t" in head[1:]:
        raise ValueError(
            "not a freyja demix file: expected the first line to be '<TAB><label>', "
            f"got {head[:60]!r}"
        )
    return f"\t{sample}{sep}{rest}"


def parse_args(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("input")
    parser.add_argument("--sample", required=True)
    parser.add_argument("--output", required=True)
    return parser.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    try:
        with open(args.input, newline="") as handle:
            text = handle.read()
        out = relabel(text, args.sample)
    except (OSError, ValueError) as err:
        print(f"relabel_demix.py: {args.input}: {err}", file=sys.stderr)
        return 1
    with open(args.output, "w", newline="") as handle:
        handle.write(out)
    return 0


if __name__ == "__main__":
    sys.exit(main())
