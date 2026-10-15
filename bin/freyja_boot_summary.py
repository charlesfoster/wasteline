#!/usr/bin/env python3
"""Per-lineage bootstrap intervals from `freyja boot`.

  freyja_boot_summary.py --sample S --lineages S.boot_lineages.csv --demix S.demix -o S.boot_summary.tsv

`freyja boot` writes a quantile table (rows 0.025, 0.05, 0.25, 0.5, 0.75, 0.95, 0.975; one column per
lineage). Output rows (no header; `summary/freyja_boot_summary.tsv` has the header
`sample lineage abundance lo hi`):

  abundance  the point estimate of the sample's own .demix (lineages/abundances rows)
  lo, hi     the 2.5 and 97.5 percentiles of the bootstrap replicates (rows 0.025 and 0.975)

Lineages with a point estimate or a non-zero upper bound are listed, ordered by abundance.
"""

import argparse
import csv
import sys


def read_quantiles(path):
    """{lineage: {quantile: value}} from the freyja boot CSV."""
    with open(path, newline="") as handle:
        reader = csv.reader(handle)
        header = next(reader)
        table = {name: {} for name in header[1:]}
        for row in reader:
            quantile = float(row[0])
            for name, value in zip(header[1:], row[1:], strict=True):
                table[name][quantile] = float(value)
    return table


def read_demix(path):
    rows = {}
    with open(path) as handle:
        for line in handle:
            key, _, value = line.rstrip("\n").partition("\t")
            rows[key.strip()] = value.split()
    names, values = rows.get("lineages", []), rows.get("abundances", [])
    if len(names) != len(values):
        raise ValueError(f"{path}: {len(names)} lineages but {len(values)} abundances")
    return {n: float(v) for n, v in zip(names, values, strict=True)}


def summary_rows(sample, table, demix):
    rows = []
    for name, quantiles in table.items():
        lo, hi = quantiles.get(0.025), quantiles.get(0.975)
        if lo is None or hi is None:
            raise ValueError(f"lineage {name}: the 0.025 / 0.975 rows are missing")
        abundance = demix.get(name, 0.0)
        if abundance > 0 or hi > 0:
            rows.append((abundance, name, lo, hi))
    rows.sort(key=lambda r: (-r[0], r[1]))
    return [[sample, name, f"{a:.6f}", f"{lo:.6f}", f"{hi:.6f}"] for a, name, lo, hi in rows]


def parse_args(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--sample", required=True)
    parser.add_argument("--lineages", required=True, help="<base>_lineages.csv from freyja boot")
    parser.add_argument("--demix", required=True, help="the sample's .demix (point estimates)")
    parser.add_argument("-o", "--output", required=True)
    return parser.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    try:
        rows = summary_rows(args.sample, read_quantiles(args.lineages), read_demix(args.demix))
    except (OSError, ValueError, StopIteration) as err:
        print(f"freyja_boot_summary.py: {err}", file=sys.stderr)
        return 1
    with open(args.output, "w") as handle:
        handle.writelines("\t".join(row) + "\n" for row in rows)
    return 0


if __name__ == "__main__":
    sys.exit(main())
