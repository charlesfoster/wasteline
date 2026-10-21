#!/usr/bin/env python3
"""95 % bootstrap intervals for the lineage GROUPS of the lineage summary.

  boot_group_intervals.py --aggregated freyja_aggregated.tsv --groups groups.csv --dataset dataset \\
      --boots S1.boot_lineages_boot.csv S2.boot_lineages_boot.csv ... --output lineage_group_intervals.tsv

Inputs
  --aggregated  freyja_aggregated.tsv (the point estimates; the same table lineage_groups.py reads)
  --groups      the groups CSV (docs/wiring/lineage_groups.md)
  --dataset     the Freyja dataset directory holding `lineages.yml`
  --samplesheet optional summary/samplesheet.csv (columns sample, sample_type): negative controls are left
                out, as in lineage_summary.csv (lineage_groups.py)
  --boots       the RAW bootstrap draws of `freyja boot --rawboots`, one `<sample>.boot_lineages_boot.csv`
                per sample: a CSV whose first column is an unnamed draw index and with one column per
                lineage (abundance in that draw, 0 when absent). The sample ID is the file name minus
                `.boot_lineages_boot.csv`; IDs are never rewritten.

Method: per sample and per draw every lineage is assigned to the groups exactly as `lineage_groups.py`
does (canonical name in the group's members; a lineage counts in every group it matches; `other` = the
lineages in no `in_total` group, unknown lineages included). The draw's abundances are summed per group
and the 2.5 and 97.5 percentiles are taken over the draws (linear interpolation, like pandas
`quantile`). The point estimate is the group sum of the sample's own demix (the value in
lineage_summary.csv), not the bootstrap mean. Intervals are taken on the group sums, never from adding
per-lineage intervals: lineages inside a group are negatively correlated.

Output (TSV): `sample group abundance lo hi n_boot`, one row per sample x (every group + `other`), samples
sorted, groups in the groups file's order. A sample without a usable boot file or aggregate row is absent.

Standard library only (the helpers are imported from lineage_groups.py).
"""

import argparse
import csv
import math
import sys
from pathlib import Path

import lineage_groups as lg

BOOT_SUFFIX = ".boot_lineages_boot.csv"
HEADER = ["sample", "group", "abundance", "lo", "hi", "n_boot"]
LOW, HIGH = 0.025, 0.975


def warn(message):
    print(f"boot_group_intervals.py: WARNING: {message}", file=sys.stderr)


def sample_id(path):
    name = Path(path).name
    if name.endswith(BOOT_SUFFIX) and len(name) > len(BOOT_SUFFIX):
        return name[: -len(BOOT_SUFFIX)]
    return name.removesuffix(".csv")


def percentile(values, q):
    """Linear-interpolation percentile (numpy / pandas default) of an unsorted list; q in 0..1."""
    ordered = sorted(values)
    pos = (len(ordered) - 1) * q
    below = math.floor(pos)
    above = min(below + 1, len(ordered) - 1)
    return ordered[below] + (ordered[above] - ordered[below]) * (pos - below)


def read_boots(path):
    """(lineage names, draws) from a raw boot CSV; draws[i] = list of abundances, one per lineage."""
    with open(path, newline="", encoding="utf-8-sig") as handle:
        reader = csv.reader(handle)
        header = next(reader, None)
        if not header or len(header) < 2:
            raise ValueError("no lineage columns")
        names = [c.strip() for c in header[1:]]
        draws = []
        for row in reader:
            if not row or not any(c.strip() for c in row):
                continue
            if len(row) != len(header):
                raise ValueError(f"a row has {len(row)} cells, the header {len(header)}")
            draws.append([float(c) for c in row[1:]])
    if not draws:
        raise ValueError("no bootstrap draws")
    if any(not math.isfinite(v) for draw in draws for v in draw):
        raise ValueError("non-finite abundance")
    return names, draws


def membership(name, groups, hierarchy):
    """(indices of the groups the lineage counts in, whether it counts in `other`), as summarise()."""
    canon = hierarchy.canonical(name)
    hits, hit_total = [], False
    for i, group in enumerate(groups):
        if canon is not None and canon in group["members"]:
            hits.append(i)
            hit_total = hit_total or group["in_total"]
    return hits, not hit_total


def group_draws(names, draws, groups, hierarchy):
    """Per group (+ `other` last): the group's summed abundance in every draw."""
    plan = [membership(name, groups, hierarchy) for name in names]
    out = []
    for draw in draws:
        sums = [[] for _ in range(len(groups) + 1)]
        for value, (hits, is_other) in zip(draw, plan, strict=True):
            for i in hits:
                sums[i].append(value)
            if is_other:
                sums[-1].append(value)
        out.append([math.fsum(s) for s in sums])
    return [[row[i] for row in out] for i in range(len(groups) + 1)]


def point_estimates(lineages, groups, hierarchy):
    """The sample's own group sums (and `other`), computed like lineage_groups.summarise()."""
    sums = [[] for _ in range(len(groups) + 1)]
    for name, value in lineages:
        hits, is_other = membership(name, groups, hierarchy)
        for i in hits:
            sums[i].append(value)
        if is_other:
            sums[-1].append(value)
    return [math.fsum(s) for s in sums]


def interval_rows(boots, aggregate, groups, hierarchy, decimals, negatives=frozenset()):
    labels = [g["group"].replace("\t", " ") for g in groups] + ["other"]
    rows = []
    for sample in sorted(boots):
        if sample in negatives:
            continue
        agg = aggregate.get(sample)
        if agg is None or agg["error"]:
            warn(f"{sample}: no usable row in the aggregated table; left out")
            continue
        names, draws = boots[sample]
        estimates = point_estimates(agg["lineages"], groups, hierarchy)
        per_group = group_draws(names, draws, groups, hierarchy)
        for label, estimate, values in zip(labels, estimates, per_group, strict=True):
            rows.append(
                [
                    sample,
                    label,
                    f"{estimate:.{decimals}f}",
                    f"{percentile(values, LOW):.{decimals}f}",
                    f"{percentile(values, HIGH):.{decimals}f}",
                    str(len(values)),
                ]
            )
    return rows


def parse_args(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("--aggregated", required=True, help="freyja_aggregated.tsv")
    p.add_argument("--groups", required=True, help="groups CSV (group,include,exclude,in_total)")
    p.add_argument("--dataset", required=True, help="Freyja dataset directory with lineages.yml")
    p.add_argument(
        "--boots", nargs="*", default=[], help="<sample>.boot_lineages_boot.csv files (raw draws)"
    )
    p.add_argument("--samplesheet", help="summary/samplesheet.csv; negative controls are left out")
    p.add_argument("-o", "--output", required=True, help="lineage_group_intervals.tsv")
    p.add_argument("--decimals", type=int, default=6)
    return p.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    try:
        groups, _sha = lg.read_groups(args.groups)
        hierarchy, _version = lg.load_dataset(args.dataset)
        lg.resolve_groups(groups, hierarchy)
        boots = {}
        for path in args.boots:
            sample = sample_id(path)
            try:
                boots[sample] = read_boots(path)
            except (OSError, ValueError) as err:
                warn(f"{sample}: cannot use {Path(path).name} ({err}); left out")
        aggregate = lg.read_aggregate(args.aggregated, set(boots))
        negatives = set()
        if args.samplesheet:
            _, records = lg.read_table(args.samplesheet, ",")
            negatives = {
                r["sample"] for r in records if r.get("sample_type", "").lower() == "negative"
            }
        rows = interval_rows(boots, aggregate, groups, hierarchy, args.decimals, negatives)
    except lg.UserError as err:
        print(f"boot_group_intervals.py: error: {err}", file=sys.stderr)
        return 2
    with open(args.output, "w", newline="", encoding="utf-8") as handle:
        handle.write("\t".join(HEADER) + "\n")
        handle.writelines("\t".join(row) + "\n" for row in rows)
    return 0


if __name__ == "__main__":
    sys.exit(main())
