#!/usr/bin/env python3
"""Which dropped amplicons hide the DEFINING sites of the lineages detected in a sample.

Amplicon dropout biases lineage abundance estimates when the dropped amplicon hides mutations that
identify the lineages present. This script intersects the dropped amplicons of each sample
(`amplicon_coverage.tsv`) with the mutations of the Freyja lineage dataset and the lineages Freyja
detected in that sample (`freyja_aggregated.tsv`).

Inputs
------
--amplicons   summary/amplicon_coverage.tsv: sample, amplicon, start, end (1-based inclusive),
              mean_depth, min_depth, pct_ge_min_depth, dropped (`true`/`false`).
--aggregated  summary/freyja_aggregated.tsv: Freyja aggregate (unnamed index column = sample; the
              `lineages` and `abundances` columns are parallel lists, as plain space-separated text or
              as a Python/numpy repr). `summarized` is NOT used (it is a regrouping of the same estimate).
--dataset     a Freyja dataset directory (preferred: its `lineage_mutations.json`, {lineage: [mutation]},
              and its `lineages.yml` for the lineage hierarchy), or a custom barcodes `.csv` (index =
              lineage, columns = mutations, a cell > 0 means the lineage carries the mutation). A directory
              without `lineage_mutations.json` falls back to a `usher_barcodes.csv` inside it. A `.feather`
              barcodes file (alone, or a directory holding only `usher_barcodes.feather`) is an error: pass
              the dataset directory instead.
--lineages-yml  optional explicit `lineages.yml`; default: the one inside the dataset directory (or next to
              a dataset file). Without any, parents are guessed from the lineage names (see "defining site").

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
          marker. Informative sites are the unit of the `n_informative_sites` / `pct_informative_sites_lost`
          columns (how much lineage-discriminating signal the scheme loses, regardless of the sample).
detected lineage
          A lineage whose abundance is >= `--min-abundance` (default 0.01) in the `lineages`/`abundances`
          columns of the aggregate row of the sample. Detected lineages missing from the dataset carry no
          mutations (counted and reported on stderr).
defining site (of a lineage)
          An informative site that the lineage carries and its PARENT lineage does not: barcode(l) minus
          barcode(parent(l)). The parent comes from the dataset's `lineages.yml` (`parent`, by name or
          alias, case-insensitive). A parent that has no barcode in the dataset is skipped: the nearest
          ancestor with a barcode is used. A lineage with no parent / no usable ancestor is a root and has
          ALL its (informative) mutations as defining; a RECOMBINANT root (XBB, XFW, ... carry
          `recombinant_parents` in lineages.yml) is the exception: the sites carried by any of its
          recombinant parents (`BA.1*` = BA.1, resolved by name/alias, only those with a barcode) are
          inherited, not defining, exactly like the sites of an ordinary parent. Without this a detected
          recombinant root counted ~130 "defining" sites that are all inherited (`--strict-roots`
          switches this off).
          Without a usable `lineages.yml` (or for a lineage missing from it) the parent is guessed by
          dropping trailing `.N` components of the name until a lineage with a barcode is found
          (`B.1.1.7` -> `B.1.1` -> `B.1` -> `B`); aliased names (`BA.5`) then become roots, so the result is
          an over-count; a warning says so. Why: a Freyja barcode lists hundreds of mutations inherited from
          ancestors; hiding one of those does not stop the lineage (or its relatives) from being called,
          so counting them flagged practically every dropout. Defining sites are what separates the lineage
          from its parent, i.e. the evidence for calling the lineage at all.
differential site (of a sample)
          A site carried by some but not all of the sample's detected lineages (only when >= 2 are
          detected); the sites the deconvolution can use to split the lineages of this sample. Uses the
          full barcodes (not only informative sites). Sites two detected lineages inherit from a common
          ancestor are NOT differential; sites separating different clades are, so with several clades
          present this count is large: it is reported as context and is not part of the affects rule.
dropped amplicon
          A row of the amplicon table with `dropped == true`. Its span is [start, end], 1-based inclusive.

Per dropped amplicon (`dropout_impact.tsv`)
-------------------------------------------
n_informative_sites              informative sites with position in [start, end] (overlapping amplicons
                                 each count the shared sites).
n_sites_in_detected_lineages     those carried by at least one detected lineage (inherited or defining).
                                 NA when the sample has no aggregate row (skipped / failed / negative
                                 control).
affected_lineages                the detected lineages with >= 1 DEFINING site in the span, ordered by count
                                 (descending, then name), at most 5, comma-joined as `name:count`
                                 (`RK.1.5:3,BA.3.2:1`). NA when there are none or no aggregate row.
affects_detected_lineage         `true` when a detected lineage whose abundance is >=
                                 `--min-affect-abundance` (default 0.05) has >= 1 defining site in the span,
                                 else `false` (always `false` without an aggregate row). Lineages between
                                 `--min-abundance` and `--min-affect-abundance` are listed in
                                 `affected_lineages` and counted in `n_defining_sites` but do not raise the
                                 flag: a lineage at a few per mille is noise-level and its loss does not
                                 bias the estimate in any useful way.
n_defining_sites                 distinct defining sites in the span, union over the detected lineages.
                                 NA without an aggregate row.
n_differential_sites             differential sites in the span (0 when fewer than 2 lineages are detected).
                                 NA without an aggregate row.

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
defining_sites_total                   distinct defining sites of the detected lineages inside the union of
                                       the sample's amplicon spans; NA without an aggregate row.
defining_sites_lost                    those inside the union of the dropped spans (overlaps counted once).
pct_defining_sites_lost                100 * lost / total, 2 decimals; NA when total is 0 or no aggregate row.

New columns are always appended at the end; the older columns keep their position and name. Both outputs
always carry their header. Rows are sorted by sample, then start, end, amplicon (byte order).
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
    "n_defining_sites",
    "n_differential_sites",
]
SUMMARY_COLUMNS = [
    "sample",
    "n_dropped_amplicons",
    "informative_sites_total",
    "informative_sites_in_dropped",
    "pct_informative_sites_lost",
    "n_dropped_affecting_detected_lineages",
    "defining_sites_total",
    "defining_sites_lost",
    "pct_defining_sites_lost",
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


# ---------------------------------------------------------------------------------------------
# Lineage hierarchy (lineages.yml)
# ---------------------------------------------------------------------------------------------
def _unquote(text):
    text = text.strip()
    if len(text) >= 2 and text[0] == text[-1] and text[0] in "'\"":
        return text[1:-1]
    return text


def parse_lineages_yml(path):
    """name -> {"alias": ..., "parent": ...} from Freyja's `lineages.yml`.

    Minimal reader (a list of maps of scalars plus block lists such as `children`, which are skipped),
    so no YAML library is needed. Same shape as the reader in bin/lineage_groups.py (kept as a copy so
    this script stays stand-alone).
    """
    entries, current = {}, None
    with open(path, encoding="utf-8-sig") as handle:
        for raw in handle:
            text = raw.strip()
            if not text or text.startswith("#"):
                continue
            if text.startswith("- ") or text == "-":
                item = text[1:].strip()
                if ":" not in item:  # a list value (children): not needed here
                    continue
                current = {}
                text = item
            elif current is None:
                continue
            key, sep, value = text.partition(":")
            value = value.strip()
            if not sep or value in ("", "[]"):
                continue
            if key.strip() == "name":
                entries[_unquote(value)] = current
            else:
                current[key.strip()] = _unquote(value)
    return entries


class Hierarchy:
    """Parent lookup over a lineage dataset (barcodes) with an optional `lineages.yml`."""

    def __init__(self, barcode_names, entries=None):
        self.barcodes = set(barcode_names)
        self.entries = entries or {}
        self.lookup = {name.casefold(): name for name in self.entries}
        for name, item in self.entries.items():  # aliases never override a real name
            if item.get("alias"):
                self.lookup.setdefault(item["alias"].casefold(), name)
        self.guessed = set()  # lineages whose parent was guessed from the name
        self._cache = {}

    def canonical(self, name):
        """Dataset (barcodes) name for a name or alias, or None."""
        if name in self.barcodes:
            return name
        hit = self.lookup.get(name.strip().casefold())
        return hit if hit in self.barcodes else None

    def _next_up(self, name, via_yml):
        if via_yml:
            parent = self.entries[name].get("parent")
            return self.lookup.get(parent.casefold(), parent) if parent else None
        head, dot, _tail = name.rpartition(".")
        return head if dot else None

    def parent(self, name):
        """Nearest ancestor of `name` that has a barcode, or None (root)."""
        if name in self._cache:
            return self._cache[name]
        via_yml = name in self.entries
        if not via_yml:
            self.guessed.add(name)
        seen, node, found = {name}, name, None
        while True:
            node = self._next_up(node, node in self.entries)
            if not node or node in seen:
                break
            seen.add(node)
            if node in self.barcodes:
                found = node
                break
        self._cache[name] = found
        return found

    def recombinant_parents(self, name):
        """Barcoded lineages in the `recombinant_parents` of a lineage (`BA.1*` = BA.1), not itself."""
        found = []
        for token in (self.entries.get(name, {}).get("recombinant_parents") or "").split(","):
            parent = self.canonical(token.strip().rstrip("*"))
            if parent and parent != name and parent not in found:
                found.append(parent)
        return found


def find_lineages_yml(explicit, dataset_path):
    if explicit:
        path = Path(explicit)
        if not path.is_file():
            raise DatasetError(f"lineages.yml not found: {path}")
        return path
    base = Path(dataset_path)
    base = base if base.is_dir() else base.parent
    candidate = base / "lineages.yml"
    return candidate if candidate.is_file() and candidate.stat().st_size > 0 else None


class Dataset:
    """Informative sites of a lineage dataset and per-lineage lookups."""

    def __init__(self, lineages, max_fraction, hierarchy=None, strict_roots=False):
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
        self.hierarchy = hierarchy or Hierarchy(lineages)
        self._lineage_cache = {}
        self._defining_cache = {}
        self.strict_roots = strict_roots

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

    def resolve(self, name):
        """The dataset name of a detected lineage (aliases accepted), or None if it is not in it."""
        return self.hierarchy.canonical(name)

    def defining_of(self, lineage):
        """Defining sites of `lineage`: informative sites it carries and its parent does not.

        Sorted list of (position, alt); empty for lineages not in the dataset. A root has all its
        informative sites defining, except that a recombinant root loses the sites of its
        `recombinant_parents` (unless `strict_roots`).
        """
        cached = self._defining_cache.get(lineage)
        if cached is None:
            if lineage in self.lineages:
                mine = self.informative_of(lineage)
                parent = self.hierarchy.parent(lineage)
                if parent:
                    cached = sorted(mine - self.informative_of(parent))
                else:
                    inherited = set()
                    if not self.strict_roots:
                        for relative in self.hierarchy.recombinant_parents(lineage):
                            inherited |= self.informative_of(relative)
                    cached = sorted(mine - inherited)
            else:
                cached = []
            self._defining_cache[lineage] = cached
        return cached


def sites_between(sorted_sites, start, end):
    """Sites of a sorted (position, alt) list with start <= position <= end."""
    lo = bisect_left(sorted_sites, (start, ""))
    hi = bisect_right(sorted_sites, (end, "\U0010ffff"))
    return sorted_sites[lo:hi]


def count_in_union(sorted_sites, intervals):
    """Sites of a sorted (position, alt) list inside the union of the intervals (overlaps once)."""
    return sum(len(sites_between(sorted_sites, a, b)) for a, b in merge_intervals(intervals))


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
    """(name, abundance) of lineages with abundance >= min_abundance (tolerating float round-off)."""
    return [(name, ab) for name, ab in pairs if ab >= min_abundance - 1e-12]


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
def analyse_sample(rows, dataset, detected, min_affect_abundance=0.05):
    """Return (impact_rows, summary_row_fields) for one sample.

    `detected` is None when the sample has no aggregate row, else the list of (lineage, abundance) of
    the detected lineages.
    """
    dropped = [r for r in rows if r["dropped"]]
    spans = [(r["start"], r["end"]) for r in rows]
    dropped_spans = [(r["start"], r["end"]) for r in dropped]
    if detected is not None:
        members = []  # (name as reported, abundance, informative sites, defining sites, full barcode)
        for name, abundance in detected:
            key = dataset.resolve(name)
            members.append(
                (
                    name,
                    abundance,
                    dataset.informative_of(key) if key else frozenset(),
                    dataset.defining_of(key) if key else [],
                    dataset.lineages[key] if key else frozenset(),
                )
            )
        all_defining = sorted({site for m in members for site in m[3]})
        differential = []
        if len(members) >= 2:
            counts = {}
            for m in members:
                for site in m[4]:
                    counts[site] = counts.get(site, 0) + 1
            differential = sorted(site for site, n in counts.items() if n < len(members))
    impact = []
    affecting = 0
    for row in sorted(dropped, key=lambda r: (r["start"], r["end"], r["amplicon"])):
        start, end = row["start"], row["end"]
        sites = dataset.sites_in(start, end)
        if detected is None:
            n_detected_sites, affected, affects, n_defining, n_differential = NA, NA, False, NA, NA
        else:
            carried, defining, per_lineage, affects = set(), set(), [], False
            for name, abundance, informative, own, _full in members:
                carried.update(s for s in sites if s in informative)
                hits = sites_between(own, start, end)
                if hits:
                    per_lineage.append((name, len(hits)))
                    defining.update(hits)
                    if abundance >= min_affect_abundance - 1e-12:
                        affects = True
            per_lineage.sort(key=lambda p: (-p[1], p[0]))
            n_detected_sites = len(carried)
            n_defining = len(defining)
            n_differential = len(sites_between(differential, start, end))
            affected = ",".join(f"{n}:{c}" for n, c in per_lineage[:MAX_AFFECTED]) or NA
        affecting += affects
        impact.append(
            {
                "sample": None,
                "amplicon": row["amplicon"],
                "start": start,
                "end": end,
                "mean_depth": row["mean_depth"],
                "n_informative_sites": len(sites),
                "n_sites_in_detected_lineages": n_detected_sites,
                "affects_detected_lineage": "true" if affects else "false",
                "affected_lineages": affected,
                "n_defining_sites": n_defining,
                "n_differential_sites": n_differential,
            }
        )
    total = dataset.count_in_intervals(spans)
    lost = dataset.count_in_intervals(dropped_spans)
    summary = {
        "n_dropped_amplicons": len(dropped),
        "informative_sites_total": total,
        "informative_sites_in_dropped": lost,
        "pct_informative_sites_lost": f"{100.0 * lost / total:.2f}" if total else NA,
        "n_dropped_affecting_detected_lineages": NA if detected is None else affecting,
        "defining_sites_total": NA,
        "defining_sites_lost": NA,
        "pct_defining_sites_lost": NA,
    }
    if detected is not None:
        d_total = count_in_union(all_defining, spans)
        d_lost = count_in_union(all_defining, dropped_spans)
        summary.update(
            defining_sites_total=d_total,
            defining_sites_lost=d_lost,
            pct_defining_sites_lost=f"{100.0 * d_lost / d_total:.2f}" if d_total else NA,
        )
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
    yml = find_lineages_yml(args.lineages_yml, args.dataset)
    entries = parse_lineages_yml(yml) if yml else {}
    if not entries:
        warn(
            "no lineages.yml found: parent lineages are guessed from the names, which treats aliased "
            "lineages (BA.5, XBB.1.5, ...) as roots and over-counts defining sites; use a Freyja "
            "dataset directory (it ships lineages.yml) or --lineages-yml"
        )
    dataset = Dataset(
        lineages, args.informative_max_fraction, Hierarchy(lineages, entries), args.strict_roots
    )

    amplicons = load_amplicons(args.amplicons)
    aggregated = load_aggregated(args.aggregated, set(amplicons), warn)

    unknown = set()
    impact_rows, summary_rows = [], []
    for sample in sorted(amplicons):
        if sample in aggregated:
            detected = detected_lineages(aggregated[sample], args.min_abundance)
            unknown.update(d for d, _ab in detected if dataset.resolve(d) is None)
        else:
            detected = None
        impact, fields = analyse_sample(
            amplicons[sample], dataset, detected, args.min_affect_abundance
        )
        for row in impact:
            row["sample"] = sample
        impact_rows.extend(impact)
        summary_rows.append({"sample": sample, **fields})
    if unknown:
        warn(
            f"{len(unknown)} detected lineage(s) not in the dataset (treated as carrying no "
            f"mutations): {', '.join(sorted(unknown)[:10])}{' ...' if len(unknown) > 10 else ''}"
        )
    if entries and dataset.hierarchy.guessed:
        warn(
            f"{len(dataset.hierarchy.guessed)} lineage(s) of the dataset are not in lineages.yml; their "
            "parents were guessed from the names"
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
        "--min-affect-abundance",
        type=float,
        default=0.05,
        help="a dropped amplicon is flagged (affects_detected_lineage) only when it hides defining "
        "sites of a detected lineage with abundance >= this (default 0.05)",
    )
    parser.add_argument(
        "--strict-roots",
        action="store_true",
        help="count ALL mutations of a recombinant root as defining (default: only those none of its "
        "recombinant_parents carry)",
    )
    parser.add_argument(
        "--lineages-yml",
        default=None,
        help="lineage hierarchy (default: lineages.yml inside/next to --dataset)",
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
    if not 0 <= args.min_affect_abundance <= 1:
        parser.error("--min-affect-abundance must be between 0 and 1")
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
