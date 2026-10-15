#!/usr/bin/env python3
"""Per-sample abundance comparison of two `freyja demix` results made with different barcode datasets.

  barcode_compare.py --sample S --a S.demix --b S.compare.demix \
      --version-a barcode_version.json --version-b compare_barcode_version.json -o S.barcode_compare.tsv

Writes data rows (no header) of `summary/barcode_compare.tsv`:
  sample, lineage, abundance_a, abundance_b, barcode_version_a, barcode_version_b
one row per lineage that either result reports (0 where one of them does not). A and B are the
`lineages` / `abundances` rows of the two .demix files (not the `summarized` groups).
"""

import argparse
import json
import sys

HEADER = [
    "sample",
    "lineage",
    "abundance_a",
    "abundance_b",
    "barcode_version_a",
    "barcode_version_b",
]


def read_demix(path):
    """{lineage: abundance} from the `lineages` and `abundances` rows of a .demix file."""
    rows = {}
    with open(path) as handle:
        for line in handle:
            key, _, value = line.rstrip("\n").partition("\t")
            rows[key.strip()] = value.strip()
    names = rows.get("lineages", "").replace("[", " ").replace("]", " ").replace(",", " ").split()
    values = (
        rows.get("abundances", "").replace("[", " ").replace("]", " ").replace(",", " ").split()
    )
    if len(names) != len(values):
        raise ValueError(f"{path}: {len(names)} lineages but {len(values)} abundances")
    return {name: float(value) for name, value in zip(names, values, strict=True)}


def read_version(path):
    with open(path) as handle:
        return json.load(handle).get("barcode_version") or "NA"


def compare_rows(sample, a, b, version_a, version_b):
    lineages = sorted(
        set(a) | set(b), key=lambda name: (-max(a.get(name, 0), b.get(name, 0)), name)
    )
    return [
        [sample, name, f"{a.get(name, 0.0):.6f}", f"{b.get(name, 0.0):.6f}", version_a, version_b]
        for name in lineages
    ]


def parse_args(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--sample", required=True)
    parser.add_argument("--a", required=True, help=".demix made with the primary barcodes")
    parser.add_argument(
        "--b", required=True, help=".demix made with the --compare_barcodes dataset"
    )
    parser.add_argument("--version-a", required=True, help="barcode_version.json of dataset A")
    parser.add_argument("--version-b", required=True, help="barcode_version.json of dataset B")
    parser.add_argument("-o", "--output", required=True)
    return parser.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    try:
        rows = compare_rows(
            args.sample,
            read_demix(args.a),
            read_demix(args.b),
            read_version(args.version_a),
            read_version(args.version_b),
        )
    except (OSError, ValueError) as err:
        print(f"barcode_compare.py: {err}", file=sys.stderr)
        return 1
    with open(args.output, "w") as handle:
        handle.writelines("\t".join(row) + "\n" for row in rows)
    return 0


if __name__ == "__main__":
    sys.exit(main())
