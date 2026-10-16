#!/usr/bin/env python3
"""Which dropped amplicons remove lineage-informative sites for the lineages detected in a sample.

Amplicon dropout biases lineage abundance estimates when the dropped amplicon hides mutations that
distinguish the lineages present. This script intersects the dropped amplicons of each sample
(`amplicon_coverage.tsv`) with the mutations of the Freyja lineage dataset and the lineages Freyja
detected in that sample (`freyja_aggregated.tsv`).

Inputs
------
--amplicons   summary/amplicon_coverage.tsv: sample, amplicon, start, end (1-based inclusive),
              mean_depth, min_depth, pct_ge_min_depth, dropped (`true`/`false`).
--aggregated  summary/freyja_aggregated.tsv: Freyja aggregate (unnamed index column = sample; the
              `lineages` and `abundances` columns are parallel lists, as plain space-separated text or
              as a Python/numpy repr). `summarized` is NOT used (it is a regrouping of the same estimate).
--dataset     a Freyja dataset directory (preferred: its `lineage_mutations.json`, {lineage: [mutation]}),
              or a custom barcodes `.csv` (index = lineage, columns = mutations, a cell > 0 means the
              lineage carries the mutation). A directory without `lineage_mutations.json` falls back to a
              `usher_barcodes.csv` inside it. A `.feather` barcodes file (alone, or a directory holding only
              `usher_barcodes.feather`) is an error: pass the dataset directory instead.

Definitions
-----------
mutation  A single-nucleotide substitution such as `C241T` -> (position 241, alternate base T), 1-based
          positions in the reference of the dataset. Deletions (`C11288-`, `-11288-9`), insertions
          (`+T123TA`, `A12AT`) and anything else that is not a substitution cannot be placed as one site;
          they are skipped, counted, and reported once on stderr (never fatal). They are therefore invisible
          to this analysis: lineages that differ only by an indel look identical.
site      One distinct mutation (position AND alternate base). A position with two different alternate bases
          is two sites.
informative site
          A mutation carried by more than 0 and fewer than `--informative-max-fraction` (default 0.95) of
          ALL lineages in the dataset (carried by n of N lineages: 0 < n/N < 0.95). A mutation carried by
          (almost) every lineage (fixed) cannot separate lineages; one carried by none is not a lineage
          marker.
detected lineage
          A lineage whose abundance is >= `--min-abundance` (default 0.01) in the `lineages`/`abundances`
          columns of the aggregate row of the sample. Detected lineages missing from the dataset carry no
          mutations (counted and reported on stderr).
dropped amplicon
          A row of the amplicon table with `dropped == true`. Its span is [start, end], 1-based inclusive.

Per dropped amplicon (`dropout_impact.tsv`)
-------------------------------------------
n_informative_sites              informative sites with position in [start, end] (overlapping amplicons
                                 each count the shared sites).
n_sites_in_detected_lineages     those carried by at least one detected lineage. NA when the sample has no
                                 aggregate row (skipped / failed / negative control).
affected_lineages                the detected lineages carrying informative sites in the span, ordered by
                                 count (descending, then name), at most 5, comma-joined as `name:count`
                                 (`BA.2:7,BA.5:6`). NA when there are none or no aggregate row.
affects_detected_lineage         `true` when any detected lineage carries >= 1 INFORMATIVE site in the span,
                                 else `false` (always `false` without an aggregate row). Informative-only is
                                 deliberate: a mutation carried by nearly all lineages (or none) being hidden
                                 cannot shift the estimate between lineages, so counting it would flag almost
                                 every dropout.

Per sample (`dropout_impact_summary.tsv`), one row for every sample of the amplicon table
-----------------------------------------------------------------------------------------
n_dropped_amplicons                    number of dropped amplicons.
informative_sites_total                informative sites inside the union of all amplicon spans of the sample
                                       (what the scheme can observe at all; not genome-wide).
informative_sites_in_dropped           informative sites inside the union of the dropped spans (overlaps
                                       counted once).
pct_informative_sites_lost             100 * in_dropped / total, 2 decimals; NA when total is 0.
n_dropped_affecting_detected_lineages  dropped amplicons with affects_detected_lineage true; NA when the
                                       sample has no aggregate row.

Both outputs always carry their header. Rows are sorted by sample, then start, end, amplicon (byte order).
"""

import argparse
import csv
import json
import math
import re
import sys
from bisect import bisect_left, bisect_right
from pathlib import Path

IMPACT_COLUMNS = [
    "sample",
    "amplicon",
    "start",
    "end",
    "mean_depth",
    "n_informative_sites",
    "n_sites_in_detected_lineages",
    "affects_detected_lineage",
    "affected_lineages",
]
SUMMARY_COLUMNS = [
    "sample",
    "n_dropped_amplicons",
    "informative_sites_total",
    "informative_sites_in_dropped",
    "pct_informative_sites_lost",
    "n_dropped_affecting_detected_lineages",
]
AMPLICON_REQUIRED = ["sample", "amplicon", "start", "end", "mean_depth", "dropped"]
NA = "NA"
NA_TOKENS = {"", "na", "nan", "n/a", "none", "null"}
MAX_AFFECTED = 5
LABEL_SUFFIXES = (".variants.tsv", ".variants", ".demix", ".tsv")

_SNV = re.compile(r"^([ACGTN])(\d+)([ACGTN])$", re.IGNORECASE)
_INSERTION = re.compile(r"^[ACGTN]\d+[ACGTN]{2,}$", re.IGNORECASE)
_DELETION = re.compile(r"^[ACGTN]{2,}\d+[ACGTN]$", re.IGNORECASE)
_QUOTED = re.compile(r"""'([^']*)'|"([^"]*)\"""")
_NUMBER = re.compile(r"[-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?|\bnan\b|\binf\b", re.IGNORECASE)


class DatasetError(Exception):
    """The lineage dataset cannot be used (reported as a clear message, exit status 2)."""


def is_na(value):
    return value is None or str(value).strip().lower() in NA_TOKENS


# ---------------------------------------------------------------------------------------------
# Mutation notation
# ---------------------------------------------------------------------------------------------
class SkipCounter:
    """Counts mutations that cannot be used as a substitution site."""

    def __init__(self):
        self.deletion = 0
        self.insertion = 0
        self.other = 0

    @property
    def total(self):
        return self.deletion + self.insertion + self.other

    def add(self, text):
        if "-" in text or _DELETION.match(text):
            self.deletion += 1
        elif "+" in text or _INSERTION.match(text):
            self.insertion += 1
        else:
            self.other += 1


def parse_mutation(text, skipped):
    """'C241T' -> (241, 'T'); anything else is counted in `skipped` and gives None."""
    match = _SNV.match(str(text).strip())
    if not match:
        skipped.add(str(text).strip())
        return None
    return int(match.group(2)), match.group(3).upper()


# ---------------------------------------------------------------------------------------------
# Lineage dataset
# ---------------------------------------------------------------------------------------------
def resolve_dataset(path):
    """Return ('json'|'csv', file) for a dataset directory or file."""
    path = Path(path)
    if path.is_dir():
        json_file = path / "lineage_mutations.json"
        if json_file.is_file() and json_file.stat().st_size > 0:
            return "json", json_file
        csvs = sorted(p for p in path.glob("*.csv") if "barcode" in p.name.lower())
        if len(csvs) == 1:
            return "csv", csvs[0]
        if len(csvs) > 1:
            raise DatasetError(
                f"{path} has no usable lineage_mutations.json and several barcodes .csv files "
                f"({', '.join(p.name for p in csvs)}); pass one of them as --dataset"
            )
        if any(path.glob("*.feather")):
            raise DatasetError(
                f"{path} has no usable lineage_mutations.json (only a .feather barcodes file); "
                "this analysis needs the mutation lists, so pass a Freyja dataset directory that "
                "contains lineage_mutations.json (e.g. one made by `freyja update`)"
            )
        raise DatasetError(f"{path} contains no lineage_mutations.json and no barcodes .csv")
    if not path.is_file():
        raise DatasetError(f"dataset not found: {path}")
    suffix = path.suffix.lower()
    if suffix == ".json":
        return "json", path
    if suffix == ".csv":
        return "csv", path
    if suffix == ".feather":
        raise DatasetError(
            f"{path} is a .feather barcodes file, which carries no mutation lists we can read here; "
            "pass a Freyja dataset directory containing lineage_mutations.json (or a .csv barcodes file)"
        )
    raise DatasetError(f"unsupported dataset file {path} (need a directory, .json or .csv)")


def load_json_dataset(path, skipped):
    """lineage_mutations.json -> {lineage: set of (position, alt)}."""
    try:
        with open(path, encoding="utf-8") as handle:
            data = json.load(handle)
    except (OSError, ValueError) as exc:
        raise DatasetError(f"cannot read {path}: {exc}") from exc
    if not isinstance(data, dict):
        raise DatasetError(f"{path}: expected a JSON object {{lineage: [mutations]}}")
    lineages = {}
    intern = {}
    for lineage, muts in data.items():
        if isinstance(muts, dict):
            muts = list(muts)
        elif isinstance(muts, str):
            muts = re.split(r"[,\s]+", muts.strip())
        elif not isinstance(muts, list):
            muts = []
        keys = set()
        for text in muts:
            key = parse_mutation(text, skipped)
            if key is not None:
                keys.add(intern.setdefault(key, key))
        lineages[str(lineage)] = keys
    return lineages


def load_csv_dataset(path, skipped):
    """Barcodes .csv (index lineage, columns mutations; cell > 0 = carried) -> {lineage: set}."""
    lineages = {}
    try:
        with open(path, newline="", encoding="utf-8-sig") as handle:
            reader = csv.reader(handle)
            header = next(reader, None)
            if not header:
                raise DatasetError(f"{path}: empty barcodes file")
            columns = []
            for index, name in enumerate(header[1:], start=1):
                key = parse_mutation(name, skipped)
                if key is not None:
                    columns.append((index, key))
            for row in reader:
                if not row or not row[0].strip():
                    continue
                keys = set()
                for index, key in columns:
                    if index >= len(row):
                        break
                    cell = row[index]
                    if cell in ("", "0", "0.0", "0.00"):
                        continue
                    try:
                        if float(cell) > 0:
                            keys.add(key)
                    except ValueError:
                        continue
                lineages[row[0].strip()] = keys
    except OSError as exc:
        raise DatasetError(f"cannot read {path}: {exc}") from exc
    return lineages


class Dataset:
    """Informative sites of a lineage dataset and per-lineage lookups."""

    def __init__(self, lineages, max_fraction):
        self.lineages = lineages
        self.n_lineages = len(lineages)
        counts = {}
        for muts in lineages.values():
            for key in muts:
                counts[key] = counts.get(key, 0) + 1
        n = self.n_lineages
        self.informative = frozenset(
            key for key, carried in counts.items() if 0 < carried and carried / n < max_fraction
        )
        self.sorted_sites = sorted(self.informative)
        self.positions = [pos for pos, _alt in self.sorted_sites]
        self._lineage_cache = {}

    def sites_in(self, start, end):
        """Informative sites with start <= position <= end (sorted list of (pos, alt))."""
        lo = bisect_left(self.positions, start)
        hi = bisect_right(self.positions, end)
        return self.sorted_sites[lo:hi]

    def count_in_intervals(self, intervals):
        """Informative sites inside the union of (start, end) intervals (overlaps counted once)."""
        total = 0
        for start, end in merge_intervals(intervals):
            total += bisect_right(self.positions, end) - bisect_left(self.positions, start)
        return total

    def informative_of(self, lineage):
        """Informative sites carried by `lineage` (empty for lineages not in the dataset)."""
        cached = self._lineage_cache.get(lineage)
        if cached is None:
            cached = self.lineages.get(lineage, frozenset()) & self.informative
            self._lineage_cache[lineage] = cached
        return cached


def merge_intervals(intervals):
    merged = []
    for start, end in sorted(intervals):
        if merged and start <= merged[-1][1] + 1:
            merged[-1][1] = max(merged[-1][1], end)
        else:
            merged.append([start, end])
    return merged


# ---------------------------------------------------------------------------------------------
# Freyja aggregate
# ---------------------------------------------------------------------------------------------
def parse_name_list(text):
    """Names in plain space-separated text, a Python list repr, or a numpy array repr.

    Deliberately not ast.literal_eval: for the numpy form (`['A' 'B']`, no commas) it silently
    concatenates adjacent string literals.
    """
    if is_na(text):
        return []
    quoted = [a or b for a, b in _QUOTED.findall(str(text))]
    if quoted:
        return quoted
    return str(text).strip().strip("[]").replace(",", " ").split()


def parse_float_list(text):
    """Numbers in plain text, a Python list repr or a numpy array repr."""
    if is_na(text):
        return []
    values = []
    for token in _NUMBER.findall(str(text)):
        try:
            number = float(token)
        except ValueError:
            continue
        values.append(number if math.isfinite(number) else 0.0)
    return values


def resolve_label(label, known):
    """Map an aggregate row label to a sample id (labels may still carry file suffixes)."""
    label = label.strip()
    if label in known:
        return label
    for suffix in sorted(LABEL_SUFFIXES, key=len, reverse=True):
        if label.endswith(suffix) and len(label) > len(suffix):
            stripped = label[: -len(suffix)]
            if stripped in known:
                return stripped
    return label


def load_aggregated(path, known, warn):
    """freyja_aggregated.tsv -> {sample: [(lineage, abundance), ...]} (all lineages, file order)."""
    out = {}
    with open(path, newline="", encoding="utf-8-sig") as handle:
        reader = csv.reader(handle, delimiter="\t", quotechar='"')
        header = next(reader, None)
        if not header:
            return out
        names = [h.strip() for h in header]
        if "lineages" not in names or "abundances" not in names:
            raise ValueError(f"{path}: needs `lineages` and `abundances` columns, found {names}")
        i_lin, i_abu = names.index("lineages"), names.index("abundances")
        for row in reader:
            if not row or not row[0].strip():
                continue
            sample = resolve_label(row[0], known)
            lineages = parse_name_list(row[i_lin] if i_lin < len(row) else "")
            abundances = parse_float_list(row[i_abu] if i_abu < len(row) else "")
            if len(lineages) != len(abundances):
                warn(f"{sample}: {len(lineages)} lineages but {len(abundances)} abundances")
            if sample in out:
                warn(f"duplicate aggregate row for {sample!r}; the last one is used")
            out[sample] = list(zip(lineages, abundances, strict=False))
    return out


def detected_lineages(pairs, min_abundance):
    """Names of lineages with abundance >= min_abundance (tolerating float round-off)."""
    return [name for name, ab in pairs if ab >= min_abundance - 1e-12]


# ---------------------------------------------------------------------------------------------
# Amplicon table
# ---------------------------------------------------------------------------------------------
def load_amplicons(path):
    """amplicon_coverage.tsv -> {sample: [row dict]} (every sample, dropped or not)."""
    samples = {}
    with open(path, newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        missing = [c for c in AMPLICON_REQUIRED if c not in (reader.fieldnames or [])]
        if missing:
            raise ValueError(f"{path}: missing column(s) {missing}")
        for line_number, row in enumerate(reader, start=2):
            sample = (row["sample"] or "").strip()
            if not sample:
                continue
            try:
                start, end = int(float(row["start"])), int(float(row["end"]))
            except (TypeError, ValueError):
                raise ValueError(f"{path}: line {line_number}: non-numeric start/end") from None
            samples.setdefault(sample, []).append(
                {
                    "amplicon": row["amplicon"],
                    "start": start,
                    "end": end,
                    "mean_depth": row["mean_depth"],
                    "dropped": str(row["dropped"]).strip().lower() in {"true", "t", "1", "yes"},
                }
            )
    return samples


# ---------------------------------------------------------------------------------------------
# Analysis
# ---------------------------------------------------------------------------------------------
def analyse_sample(rows, dataset, detected):
    """Return (impact_rows, summary_row_fields) for one sample.

    `detected` is None when the sample has no aggregate row, else the list of detected lineages.
    """
    detected_sets = None if detected is None else [(d, dataset.informative_of(d)) for d in detected]
    impact = []
    affecting = 0
    dropped = [r for r in rows if r["dropped"]]
    for row in sorted(dropped, key=lambda r: (r["start"], r["end"], r["amplicon"])):
        sites = dataset.sites_in(row["start"], row["end"])
        if detected_sets is None:
            n_detected_sites, affected, affects = NA, NA, False
        else:
            carried = set()
            per_lineage = []
            for name, lineage_sites in detected_sets:
                hits = [s for s in sites if s in lineage_sites]
                if hits:
                    per_lineage.append((name, len(hits)))
                    carried.update(hits)
            per_lineage.sort(key=lambda p: (-p[1], p[0]))
            n_detected_sites = len(carried)
            affects = bool(per_lineage)
            affected = ",".join(f"{n}:{c}" for n, c in per_lineage[:MAX_AFFECTED]) or NA
        affecting += affects
        impact.append(
            {
                "sample": None,
                "amplicon": row["amplicon"],
                "start": row["start"],
                "end": row["end"],
                "mean_depth": row["mean_depth"],
                "n_informative_sites": len(sites),
                "n_sites_in_detected_lineages": n_detected_sites,
                "affects_detected_lineage": "true" if affects else "false",
                "affected_lineages": affected,
            }
        )
    total = dataset.count_in_intervals([(r["start"], r["end"]) for r in rows])
    lost = dataset.count_in_intervals([(r["start"], r["end"]) for r in dropped])
    summary = {
        "n_dropped_amplicons": len(dropped),
        "informative_sites_total": total,
        "informative_sites_in_dropped": lost,
        "pct_informative_sites_lost": f"{100.0 * lost / total:.2f}" if total else NA,
        "n_dropped_affecting_detected_lineages": NA if detected is None else affecting,
    }
    return impact, summary


def write_table(path, columns, rows):
    with open(path, "w", newline="") as out:
        out.write("\t".join(columns) + "\n")
        out.writelines("\t".join(str(row[c]) for c in columns) + "\n" for row in rows)


def run(args):
    def warn(message):
        print(f"warning: {message}", file=sys.stderr)

    kind, dataset_file = resolve_dataset(args.dataset)
    skipped = SkipCounter()
    loader = load_json_dataset if kind == "json" else load_csv_dataset
    lineages = loader(dataset_file, skipped)
    if skipped.total:
        warn(
            f"{dataset_file.name}: skipped {skipped.total} mutation(s) that are not single-nucleotide "
            f"substitutions ({skipped.deletion} deletion, {skipped.insertion} insertion, "
            f"{skipped.other} other); they are not part of this analysis"
        )
    if not lineages:
        raise DatasetError(f"{dataset_file}: no lineages found")
    dataset = Dataset(lineages, args.informative_max_fraction)

    amplicons = load_amplicons(args.amplicons)
    aggregated = load_aggregated(args.aggregated, set(amplicons), warn)

    unknown = set()
    impact_rows, summary_rows = [], []
    for sample in sorted(amplicons):
        if sample in aggregated:
            detected = detected_lineages(aggregated[sample], args.min_abundance)
            unknown.update(d for d in detected if d not in dataset.lineages)
        else:
            detected = None
        impact, fields = analyse_sample(amplicons[sample], dataset, detected)
        for row in impact:
            row["sample"] = sample
        impact_rows.extend(impact)
        summary_rows.append({"sample": sample, **fields})
    if unknown:
        warn(
            f"{len(unknown)} detected lineage(s) not in the dataset (treated as carrying no "
            f"mutations): {', '.join(sorted(unknown)[:10])}{' ...' if len(unknown) > 10 else ''}"
        )
    write_table(args.out, IMPACT_COLUMNS, impact_rows)
    write_table(args.summary_out, SUMMARY_COLUMNS, summary_rows)
    return 0


def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("--amplicons", required=True, help="summary/amplicon_coverage.tsv")
    parser.add_argument("--aggregated", required=True, help="summary/freyja_aggregated.tsv")
    parser.add_argument(
        "--dataset",
        required=True,
        help="Freyja dataset directory (lineage_mutations.json) or a barcodes .csv",
    )
    parser.add_argument(
        "--min-abundance",
        type=float,
        default=0.01,
        help="a lineage is detected at abundance >= this (default 0.01)",
    )
    parser.add_argument(
        "--informative-max-fraction",
        type=float,
        default=0.95,
        help="a mutation is informative when carried by fewer than this fraction of the dataset's "
        "lineages (default 0.95)",
    )
    parser.add_argument("--out", required=True, help="dropout_impact.tsv")
    parser.add_argument("--summary-out", required=True, help="dropout_impact_summary.tsv")
    args = parser.parse_args(argv)
    if not 0 <= args.min_abundance <= 1:
        parser.error("--min-abundance must be between 0 and 1")
    if not 0 < args.informative_max_fraction <= 1:
        parser.error("--informative-max-fraction must be in (0, 1]")
    return args


def main(argv=None):
    args = parse_args(argv)
    try:
        return run(args)
    except (DatasetError, ValueError, OSError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    sys.exit(main())
