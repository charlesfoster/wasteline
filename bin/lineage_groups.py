#!/usr/bin/env python3
"""Summarise each sample's Freyja lineage solution into a user-defined list of lineage groups.

Inputs
  --aggregated   freyja_aggregated.tsv (Freyja's aggregate; the `lineages` / `abundances` columns are
                 used, never `summarized`)
  --groups       the groups CSV (see docs/wiring/lineage_groups.md): columns `group, include,
                 [exclude], [in_total]`. Row order = output column order.
  --dataset      the Freyja dataset directory holding `lineages.yml` (and `last_barcode_update.txt`);
                 patterns are resolved against that hierarchy.
  --demix-status summary/demix_status.tsv (every sample listed there gets a row, whatever happened)
  --samplesheet  summary/samplesheet.csv (site, collection_date; negative controls are left out)
  --aggregated-sites / --aggregated-sites-meta / --sites-out
                 optional site-aggregate rows (Freyja aggregate format + a meta TSV) -> a second CSV

Pattern syntax (case-insensitive; one or more per cell separated by `;`)
  NAME*          the lineage and all its descendants in the hierarchy
  NAME           exactly that lineage
  RECOMBINANT*   every lineage that descends from a recombinant root (a root with `recombinant_parents`)
NAME may be a Pango alias (`XDV.1.5.1.1.8`) or the full name; both resolve through `lineages.yml`.
A lineage belongs to a group when it matches an include pattern and no exclude pattern. A pattern naming
a lineage that is not in the dataset is reported once and matches nothing (no prefix guessing).

Output (`lineage_summary.csv`, comma separated, UTF-8 without BOM): `sample, site, collection_date,
status, genome_coverage_pct, <groups...>, other, proportion_summary, major_lineages,
major_abundances, minor_lineages, minor_abundances, ungrouped_lineages, status_reason`.
Group values are sums over every lineage of the sample's solution; `--min-abundance` only limits the
lists. `other` is what no `in_total` group covers, so `proportion_summary` is ~1 unless two `in_total`
groups overlap (reported as a warning).

`--check-groups` only validates the groups file (format errors exit 2 with a message naming the row)
and, when --dataset is given, resolves its patterns against the dataset hierarchy, printing warnings
(unknown lineages, overlapping `in_total` groups); no aggregate or output files are needed. Run it early so a typo stops a run at once.

Standard library only.
"""

import argparse
import csv
import hashlib
import math
import re
import sys
from pathlib import Path

RESERVED = {
    "sample",
    "site",
    "collection_date",
    "status",
    "genome_coverage_pct",
    "other",
    "proportion_summary",
    "major_lineages",
    "major_abundances",
    "minor_lineages",
    "minor_abundances",
    "ungrouped_lineages",
    "status_reason",
}
TAIL_COLUMNS = [
    "other",
    "proportion_summary",
    "major_lineages",
    "major_abundances",
    "minor_lineages",
    "minor_abundances",
    "ungrouped_lineages",
]
SITE_META_COLUMNS = [
    "aggregate_name",
    "period",
    "period_start",
    "n_sites_included",
    "sites_included",
    "sites_missing",
    "weight_basis",
]
RECOMBINANT_TOKEN = "recombinant"
LABEL_SUFFIXES = (".variants.tsv", ".demix", ".variants", ".tsv")
EPS = 1e-9  # tolerance so 0.0999999999 counts as the 0.10 threshold
YES = {"yes", "y", "true", "t", "1"}
NO = {"no", "n", "false", "f", "0"}
_NUMBER = re.compile(r"[-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?|\bnan\b|\binf\b", re.IGNORECASE)
_QUOTED = re.compile(r"""'([^']*)'|"([^"]*)\"""")


class UserError(Exception):
    """A problem with the user's input that deserves a friendly message, not a traceback."""


def warn(message):
    print(f"lineage_groups.py: WARNING: {message}", file=sys.stderr)


# --------------------------------------------------------------------------------------------------
# lineage hierarchy
# --------------------------------------------------------------------------------------------------
def _unquote(text):
    text = text.strip()
    if len(text) >= 2 and text[0] == text[-1] and text[0] in "'\"":
        return text[1:-1]
    return text


def parse_lineages_yml(path):
    """Minimal reader for Freyja's `lineages.yml` (a list of `name/alias/children/parent/...` maps).

    Only what that file uses is supported (scalar values and one block list, `children`), so the
    pipeline needs no YAML library. Returns a list of dicts.
    """
    entries, current, list_key = [], None, None
    with open(path, encoding="utf-8-sig") as handle:
        for number, raw in enumerate(handle, 1):
            line = raw.rstrip("\r\n")
            text = line.strip()
            if not text or text.startswith("#"):
                continue
            if text.startswith("- ") or text == "-":
                item = text[1:].strip()
                if ":" not in item and current is not None and list_key is not None:
                    current[list_key].append(_unquote(item))
                    continue
                current, list_key = {}, None
                entries.append(current)
                text = item
            if current is None:
                raise UserError(f"{path}:{number}: unexpected content before the first list item")
            key, sep, value = text.partition(":")
            if not sep:
                raise UserError(f"{path}:{number}: cannot read line {text!r}")
            key, value = key.strip(), value.strip()
            if value == "":
                current[key] = []
                list_key = key
            elif value == "[]":
                current[key] = []
                list_key = None
            else:
                current[key] = _unquote(value)
                list_key = None
    return entries


class Hierarchy:
    """Names, aliases, ancestors and descendants of the dataset's lineages."""

    def __init__(self, entries):
        self.entry = {}
        self.lookup = {}  # casefolded name or alias -> canonical name
        for item in entries:
            name = item.get("name")
            if not name:
                continue
            self.entry[name] = item
            self.lookup[name.casefold()] = name
        for name, item in self.entry.items():  # aliases never override a real name
            alias = item.get("alias")
            if alias:
                self.lookup.setdefault(alias.casefold(), name)
        self.descendants = {name: {name} for name in self.entry}
        for name, item in self.entry.items():
            seen, node = {name}, item.get("parent")
            while node and node not in seen:
                node = self.lookup.get(node.casefold(), node)
                if node not in self.entry or node in seen:
                    break
                self.descendants[node].add(name)
                seen.add(node)
                node = self.entry[node].get("parent")
            for child in item.get("children") or []:
                canon = self.lookup.get(str(child).casefold())
                if canon:
                    self.descendants[name].add(canon)
        self.recombinants = set()
        for name in self.entry:
            if self.root(name) in self.entry and self.entry[self.root(name)].get(
                "recombinant_parents"
            ):
                self.recombinants.add(name)

    def root(self, name):
        seen = {name}
        while True:
            parent = self.entry[name].get("parent")
            parent = self.lookup.get(parent.casefold()) if parent else None
            if not parent or parent in seen or parent not in self.entry:
                return name
            seen.add(parent)
            name = parent

    def canonical(self, name):
        """Canonical dataset name for a name or alias (any case), or None."""
        return self.lookup.get(name.strip().casefold())


def load_dataset(path):
    """Return (Hierarchy, version) for a dataset directory (or a lineages.yml file)."""
    path = Path(path)
    yml = path if path.is_file() else path / "lineages.yml"
    if not yml.is_file():
        raise UserError(
            f"no lineages.yml found in {path}; a Freyja dataset directory with the lineage "
            "hierarchy is needed to resolve the groups"
        )
    entries = parse_lineages_yml(yml)
    if not entries:
        raise UserError(f"{yml} contains no lineages")
    version_file = yml.parent / "last_barcode_update.txt"
    version = version_file.read_text().strip() if version_file.is_file() else "NA"
    return Hierarchy(entries), version or "NA"


# --------------------------------------------------------------------------------------------------
# groups file
# --------------------------------------------------------------------------------------------------
def split_patterns(cell):
    return [p for p in re.split(r"[;,\s]+", (cell or "").strip()) if p]


def _delimiter(text):
    """Excel in some locales saves `;`-separated CSV: pick the delimiter that splits the header."""
    first = next((ln for ln in text.splitlines() if ln.strip()), "")
    best, best_n = ",", 0
    for delim in (",", ";", "\t"):
        cells = [c.strip().lower() for c in next(csv.reader([first], delimiter=delim), [])]
        if "group" in cells and "include" in cells and len(cells) > best_n:
            best, best_n = delim, len(cells)
    return best


def read_groups(path):
    """Read and validate the groups CSV -> (list of group dicts, sha256 of the file)."""
    raw = Path(path).read_bytes()
    sha = hashlib.sha256(raw).hexdigest()
    text = raw.decode("utf-8-sig", errors="replace")
    if not text.strip():
        raise UserError(
            f"the groups file {path} is empty: it needs a header (group,include) and rows"
        )
    reader = csv.reader(text.splitlines(), delimiter=_delimiter(text))
    rows = [r for r in reader if any(c.strip() for c in r)]
    header = [c.strip().lower() for c in rows[0]]
    for needed in ("group", "include"):
        if needed not in header:
            raise UserError(
                f"the groups file {path} needs the columns 'group' and 'include' in its first "
                f"row (found: {', '.join(c for c in header if c) or 'nothing'})"
            )
    col = {name: i for i, name in reversed(list(enumerate(header)))}

    def cell(row, name):
        i = col.get(name)
        return row[i].strip() if i is not None and i < len(row) else ""

    groups, seen = [], {}
    for number, row in enumerate(rows[1:], 2):
        name = cell(row, "group")
        include = cell(row, "include")
        if not name:
            raise UserError(f"{path}, row {number}: the 'group' cell is empty but the row is not")
        if name.casefold() in seen:
            raise UserError(
                f"{path}: the group name {name!r} appears twice (rows {seen[name.casefold()]} "
                f"and {number}); group names must be unique"
            )
        if name.casefold() in RESERVED:
            raise UserError(f"{path}, row {number}: {name!r} is a reserved output column name")
        seen[name.casefold()] = number
        patterns = split_patterns(include)
        if not patterns:
            raise UserError(f"{path}, row {number}: group {name!r} has no 'include' pattern")
        flag = cell(row, "in_total").casefold()
        if flag == "" or flag in YES:
            in_total = True
        elif flag in NO:
            in_total = False
        else:
            raise UserError(
                f"{path}, row {number}: in_total must be yes or no for group {name!r}, got "
                f"{cell(row, 'in_total')!r}"
            )
        groups.append(
            {
                "group": name,
                "include": patterns,
                "exclude": split_patterns(cell(row, "exclude")),
                "in_total": in_total,
            }
        )
    if not groups:
        raise UserError(f"the groups file {path} has a header but no groups")
    return groups, sha


def resolve_pattern(pattern, hierarchy, problems):
    """Set of canonical lineages matched by one pattern; unknown ones are noted in `problems`."""
    token = pattern.strip()
    star = token.endswith("*")
    base = token[:-1] if star else token
    if base.casefold() == RECOMBINANT_TOKEN:
        return set(hierarchy.recombinants)
    if not base or "*" in base:
        problems.append(f"pattern {pattern!r} is not valid (only one trailing * is allowed)")
        return set()
    canon = hierarchy.canonical(base)
    if canon is None:
        hint = ""
        if re.fullmatch(r"x\*?", token, re.IGNORECASE):
            hint = "; for all recombinants write RECOMBINANT*"
        problems.append(f"lineage {base!r} (pattern {pattern!r}) is not in the dataset{hint}")
        return set()
    return set(hierarchy.descendants[canon]) if star else {canon}


def resolve_groups(groups, hierarchy):
    """Add `members` (frozenset of canonical lineages) and `warnings` to each group; check overlaps."""
    reported = set()
    for group in groups:
        problems = []
        include = set()
        for pattern in group["include"]:
            include |= resolve_pattern(pattern, hierarchy, problems)
        exclude = set()
        for pattern in group["exclude"]:
            exclude |= resolve_pattern(pattern, hierarchy, problems)
        group["members"] = frozenset(include - exclude)
        if not group["members"] and not problems:
            problems.append("matches no lineage in the dataset")
        group["warnings"] = problems
        for problem in problems:
            if problem not in reported:
                reported.add(problem)
                warn(f"group {group['group']!r}: {problem}")
    counted = [g for g in groups if g["in_total"]]
    for i, later in enumerate(counted):
        for earlier in counted[:i]:
            both = sorted(earlier["members"] & later["members"])
            if both:
                shown = ", ".join(both[:8]) + (" ..." if len(both) > 8 else "")
                text = (
                    f"in_total group overlaps in_total group {earlier['group']!r} on {len(both)} "
                    f"lineage(s) ({shown}); proportion_summary can exceed 1 "
                    "(set in_total to no for the informational group)"
                )
                later["warnings"].append(text)
                warn(f"group {later['group']!r}: {text}")


# --------------------------------------------------------------------------------------------------
# aggregate / metadata parsing
# --------------------------------------------------------------------------------------------------
def is_na(value):
    return value is None or str(value).strip() in ("", "NA", "NaN", "nan", "None", "[]")


def parse_float_list(text):
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


def parse_name_list(text):
    """`A B`, `['A', 'B']` and numpy `['A' 'B']` forms."""
    if is_na(text):
        return []
    quoted = [a or b for a, b in _QUOTED.findall(str(text))]
    if quoted:
        return quoted
    return str(text).strip().strip("[]").replace(",", " ").split()


def clean_label(label, known=()):
    label = label.strip()
    if label in known:
        return label
    for suffix in LABEL_SUFFIXES:
        if label.endswith(suffix) and len(label) > len(suffix):
            return label[: -len(suffix)]
    return label


def read_aggregate(path, known=()):
    """Freyja aggregate -> {label: {lineages: [(name, frac)], coverage, error}}."""
    with open(path, newline="", encoding="utf-8-sig") as handle:
        table = list(csv.reader(handle, delimiter="\t"))
    if not table:
        raise UserError(f"{path} is empty")
    header = [c.strip() for c in table[0]]
    col = {name: i for i, name in enumerate(header) if name}
    for needed in ("lineages", "abundances"):
        if needed not in col:
            raise UserError(f"{path}: no '{needed}' column (is this a Freyja aggregate?)")
    rows = {}
    for line in table[1:]:
        if not line or not any(c.strip() for c in line) or not line[0].strip():
            continue
        label = clean_label(line[0], known)

        def cell(name, line=line):
            i = col.get(name)
            return line[i] if i is not None and i < len(line) else ""

        names = parse_name_list(cell("lineages"))
        values = parse_float_list(cell("abundances"))
        error = None
        if not names or len(names) != len(values):
            error = "the lineages / abundances of this row could not be parsed"
        coverage = cell("coverage").strip()
        rows[label] = {
            "lineages": [] if error else list(zip(names, values)),
            "coverage": coverage,
            "error": error,
        }
    return rows


def read_table(path, delimiter, encoding="utf-8-sig"):
    with open(path, newline="", encoding=encoding) as handle:
        reader = csv.DictReader(handle, delimiter=delimiter)
        reader.fieldnames = [(n or "").strip().lower() for n in (reader.fieldnames or [])]
        return reader.fieldnames, [
            {k: (v or "").strip() for k, v in row.items() if k} for row in reader
        ]


# --------------------------------------------------------------------------------------------------
# summarising
# --------------------------------------------------------------------------------------------------
class Formatter:
    def __init__(self, decimals, zeros):
        self.decimals = decimals
        self.zeros = zeros

    def number(self, value):
        text = f"{value:.{self.decimals}f}"
        if self.zeros == "blank" and float(text) == 0:
            return ""
        return text

    def plain(self, value):
        return f"{value:.{self.decimals}f}"


def _by_abundance(item):
    """Descending abundance; equal abundances alphabetically."""
    return (-item[1], item[0].casefold(), item[0])


def summarise(
    lineages, groups, hierarchy, fmt, major_threshold, min_abundance, unknown, collect=None
):
    """One row of group values for a list of (name, abundance)."""
    sums = [[] for _ in groups]
    other, ungrouped, listed = [], [], []
    for name, value in lineages:
        canon = hierarchy.canonical(name)
        if canon is None and name not in unknown:
            unknown.add(name)
            warn(
                f"lineage {name!r} is in the aggregate but not in the dataset hierarchy; counted as 'other'"
            )
        hit = False
        for i, group in enumerate(groups):
            if canon is not None and canon in group["members"]:
                sums[i].append(value)
                hit = hit or group["in_total"]
        if not hit:
            other.append(value)
            if value >= min_abundance:
                ungrouped.append((name, value))
                if collect is not None:
                    collect.setdefault(name, []).append(value)
        if value >= min_abundance:
            listed.append((name, value))
    listed.sort(key=_by_abundance)
    ungrouped.sort(key=_by_abundance)
    major = [item for item in listed if item[1] >= major_threshold - EPS]
    minor = [item for item in listed if item[1] < major_threshold - EPS]
    group_totals = [math.fsum(s) for s in sums]
    other_total = math.fsum(other)
    total = (
        math.fsum(t for t, g in zip(group_totals, groups, strict=True) if g["in_total"])
        + other_total
    )
    out = [fmt.number(t) for t in group_totals]
    out += [fmt.number(other_total), fmt.number(total)]
    out += [
        " ".join(n for n, _ in major),
        " ".join(fmt.plain(v) for _, v in major),
        " ".join(n for n, _ in minor),
        " ".join(fmt.plain(v) for _, v in minor),
        " ".join(f"{n}:{fmt.plain(v)}" for n, v in ungrouped),
    ]
    return out


def sample_rows(args, groups, hierarchy, fmt, ungrouped=None):
    sheet = {}
    if args.samplesheet:
        _, records = read_table(args.samplesheet, ",")
        sheet = {r.get("sample", ""): r for r in records if r.get("sample")}
    status = {}
    if args.demix_status:
        _, records = read_table(args.demix_status, "\t")
        status = {r.get("sample", ""): r for r in records if r.get("sample")}
    aggregate = read_aggregate(args.aggregated, set(status) | set(sheet))
    names = sorted(set(status) | set(aggregate))
    width = len(groups) + len(TAIL_COLUMNS)
    unknown, rows = set(), []
    for sample in names:
        info = sheet.get(sample, {})
        if info.get("sample_type", "").lower() == "negative":
            continue
        st = status.get(sample)
        agg = aggregate.get(sample)
        state, reason = (st.get("status") or "ok", st.get("reason", "")) if st else ("ok", "")
        coverage = (st.get("genome_coverage_pct", "") if st else "") or (
            agg["coverage"] if agg else ""
        )
        if coverage:
            try:
                coverage = f"{float(coverage):.2f}"
            except ValueError:
                pass
        values = [""] * width
        if state == "ok":
            if agg is None:
                state, reason = "missing", "no row for this sample in the aggregated table"
            elif agg["error"]:
                state, reason = "failed", agg["error"]
            else:
                values = summarise(
                    agg["lineages"],
                    groups,
                    hierarchy,
                    fmt,
                    args.major_threshold,
                    args.min_abundance,
                    unknown,
                    ungrouped,
                )
        rows.append(
            [sample, info.get("site", ""), info.get("collection_date", ""), state, coverage]
            + values
            + [reason]
        )
    return rows


def site_rows(args, groups, hierarchy, fmt):
    meta_fields, meta = read_table(args.aggregated_sites_meta, "\t")
    for needed in ("row_label", "aggregate_name", "period"):
        if needed not in meta_fields:
            raise UserError(f"{args.aggregated_sites_meta}: no '{needed}' column")
    aggregate = read_aggregate(args.aggregated_sites, {m["row_label"] for m in meta})
    used, unknown, rows = set(), set(), []
    width = len(groups) + len(TAIL_COLUMNS)
    for record in meta:
        label = record["row_label"]
        agg = aggregate.get(label)
        values = [""] * width
        if agg is None or agg["error"]:
            warn(f"site aggregate {label!r}: no usable row in {args.aggregated_sites}")
        else:
            used.add(label)
            values = summarise(
                agg["lineages"],
                groups,
                hierarchy,
                fmt,
                args.major_threshold,
                args.min_abundance,
                unknown,
            )
        rows.append([record.get(c, "") for c in SITE_META_COLUMNS] + values)
    for label in sorted(set(aggregate) - used):
        warn(f"site aggregate row {label!r} has no entry in {args.aggregated_sites_meta}; skipped")
    return rows


# --------------------------------------------------------------------------------------------------
# output
# --------------------------------------------------------------------------------------------------
def write_csv(path, header, rows):
    with open(path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(header)
        writer.writerows(rows)


UNGROUPED_COLUMNS = ["lineage", "n_samples", "mean_abundance", "max_abundance", "suggested_pattern"]


def write_ungrouped(path, ungrouped, fmt):
    """Lineages that fell into `other` in at least one ok sample, biggest first (max abundance)."""
    items = sorted(ungrouped.items(), key=lambda kv: (-max(kv[1]), kv[0].casefold(), kv[0]))
    with open(path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(UNGROUPED_COLUMNS)
        for name, values in items:
            writer.writerow(
                [
                    name,
                    len(values),
                    fmt.plain(math.fsum(values) / len(values)),
                    fmt.plain(max(values)),
                    f"{name}*",
                ]
            )


def default_groups_path():
    return Path(__file__).resolve().parent.parent / "assets" / "lineage_groups.csv"


def write_resolved(path, groups, version, sha):
    with open(path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(
            [
                "group",
                "include",
                "exclude",
                "in_total",
                "n_lineages_matched",
                "warnings",
                "dataset_version",
                "groups_sha256",
            ]
        )
        for g in groups:
            writer.writerow(
                [
                    g["group"],
                    ";".join(g["include"]),
                    ";".join(g["exclude"]),
                    "yes" if g["in_total"] else "no",
                    len(g["members"]),
                    " | ".join(g["warnings"]),
                    version,
                    sha,
                ]
            )


def parse_args(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("--aggregated", help="freyja_aggregated.tsv")
    p.add_argument("--groups", help="groups CSV (group,include,exclude,in_total)")
    p.add_argument("--dataset", help="Freyja dataset directory with lineages.yml")
    p.add_argument("--demix-status", help="summary/demix_status.tsv")
    p.add_argument("--samplesheet", help="summary/samplesheet.csv")
    p.add_argument("--major-threshold", type=float, default=0.10)
    p.add_argument("--min-abundance", type=float, default=0.001)
    p.add_argument("--decimals", type=int, default=6)
    p.add_argument("--zeros", choices=["blank", "zero"], default="blank")
    p.add_argument("--label-column", default="sample", help="header of the first output column")
    p.add_argument("-o", "--output", help="lineage_summary.csv")
    p.add_argument("--resolved-out", help="lineage_groups_resolved.tsv")
    p.add_argument("--aggregated-sites", help="site aggregates in Freyja aggregate format")
    p.add_argument("--aggregated-sites-meta", help="site_aggregates_meta.tsv")
    p.add_argument("--sites-out", help="lineage_summary_site_aggregates.csv")
    p.add_argument(
        "--ungrouped-out",
        help="lineage_ungrouped_summary.tsv: the lineages that fall into `other`, biggest first",
    )
    p.add_argument(
        "--print-default-groups",
        action="store_true",
        help="print the bundled default groups CSV (assets/lineage_groups.csv next to bin/) and exit",
    )
    p.add_argument(
        "--check-groups",
        action="store_true",
        help="only validate the groups file and resolve it against --dataset (warnings, exit 2 on errors)",
    )
    args = p.parse_args(argv)
    if args.print_default_groups:
        return args
    if not args.groups:
        p.error("the following arguments are required: --groups")
    if not args.check_groups:
        missing = [
            f"--{n.replace('_', '-')}"
            for n in ("dataset", "aggregated", "output", "resolved_out")
            if not getattr(args, n)
        ]
        if missing:
            p.error("the following arguments are required: " + ", ".join(missing))
    if bool(args.aggregated_sites) != bool(args.aggregated_sites_meta):
        p.error("--aggregated-sites and --aggregated-sites-meta go together")
    if args.aggregated_sites and not args.sites_out:
        p.error("--aggregated-sites needs --sites-out")
    return args


def main(argv=None):
    args = parse_args(argv)
    if args.print_default_groups:
        path = default_groups_path()
        if not path.is_file():
            print(
                f"lineage_groups.py: error: the bundled default {path} was not found; "
                "--print-default-groups only works from a clone of the pipeline "
                "(copy assets/lineage_groups.csv instead)",
                file=sys.stderr,
            )
            return 2
        sys.stdout.write(path.read_text(encoding="utf-8"))
        return 0
    try:
        groups, sha = read_groups(args.groups)
        if args.check_groups and not args.dataset:
            # format only: the dataset has no lineage hierarchy (custom barcodes) or is not known yet
            print(
                f"lineage_groups.py: groups file format OK, {len(groups)} group(s); "
                "no dataset given, so the patterns were not resolved"
            )
            return 0
        hierarchy, version = load_dataset(args.dataset)
        resolve_groups(groups, hierarchy)
        if args.check_groups:
            print(
                f"lineage_groups.py: groups file OK, {len(groups)} group(s) "
                f"({sum(1 for g in groups if g['in_total'])} counted in the total); "
                f"dataset {version}; "
                f"{sum(1 for g in groups if g['warnings'])} group(s) with warnings"
            )
            return 0
        fmt = Formatter(args.decimals, args.zeros)
        names = [g["group"] for g in groups]
        header = [args.label_column, "site", "collection_date", "status", "genome_coverage_pct"]
        header += names + TAIL_COLUMNS + ["status_reason"]
        ungrouped = {}
        write_csv(args.output, header, sample_rows(args, groups, hierarchy, fmt, ungrouped))
        if args.ungrouped_out:
            write_ungrouped(args.ungrouped_out, ungrouped, fmt)
        write_resolved(args.resolved_out, groups, version, sha)
        if args.aggregated_sites:
            write_csv(
                args.sites_out,
                SITE_META_COLUMNS + names + TAIL_COLUMNS,
                site_rows(args, groups, hierarchy, fmt),
            )
    except UserError as err:
        print(f"lineage_groups.py: error: {err}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    sys.exit(main())
