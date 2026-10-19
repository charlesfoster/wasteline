#!/usr/bin/env python3
"""Combine per-site Freyja results into user-defined groups of sites, without pooling reads.

Every site-sample is demixed on its own, exactly as without this step. This script then combines the
finished abundance solutions of the member sites of an aggregate by a weighted average, per period
(ISO week or day). Nothing is re-fitted, so an aggregate can never contain a lineage that none of its
member sites reported.

Inputs: `--aggregated` (Freyja aggregate `freyja_aggregated.tsv`, row label = sample ID),
`--demix-status` (`demix_status.tsv`; only status `ok` samples contribute), `--samplesheet`
(`sample, sample_type, collection_date, site, viral_load`, extra columns ignored) and `--sites`, a
long-format CSV with columns `aggregate_name, site_name[, weight]`.

Algorithm, per aggregate and period:

1. A sample contributes if its demix status is `ok`, its `sample_type` is not negative, and it has a
   site and a valid ISO `collection_date`. Sites are matched to the sites file case-insensitively
   after trimming (and collapsing) spaces.
2. For each member site with >= 1 contributing sample in the period, its samples (replicates) are
   combined first: weighted by `viral_load` if every one of them has a load, else equally. A site with
   two samples therefore never counts double.
3. The sites are combined with `--weights`: `viral_load` (the site's mean load), `site_weight` (the
   `weight` column of the sites file), `equal`, or `auto` (viral_load if every contributing site has a
   load, else site_weight if every contributing site has a weight, else equal). An explicit choice that
   cannot be met falls back to `equal` and says so in `weight_basis`.
4. abundance = sum(w_i * a_i) / sum(w_i) over the sites, for the lineages, Freyja's `summarized`
   groups, `resid` and `coverage` alike (a lineage a site lacks counts as 0). Nothing is renormalised.
5. Periods with fewer than `--min-sites` contributing sites are omitted.

Outputs (both always written with a header): `--out` in Freyja's aggregate format (row label
`<aggregate_name>_<period>`; columns `summarized, lineages, abundances, resid, coverage`) and
`--meta-out` with `row_label, aggregate_name, period, period_start, n_sites_included, sites_included,
sites_missing, missing_reasons, weight_basis, weight_fraction_included, n_samples`.

Standard library only.
"""

import argparse
import csv
import datetime as dt
import difflib
import math
import re
import sys

NA = "NA"
SITES_COLUMNS = ["aggregate_name", "site_name", "weight"]
META_COLUMNS = [
    "row_label",
    "aggregate_name",
    "period",
    "period_start",
    "n_sites_included",
    "sites_included",
    "sites_missing",
    "missing_reasons",
    "weight_basis",
    "weight_fraction_included",
    "n_samples",
]
AGG_HEADER = "\tsummarized\tlineages\tabundances\tresid\tcoverage"
PLACES = 8

_QUOTED = re.compile(r"""'([^']*)'|"([^"]*)\"""")
_NUMBER = re.compile(r"[-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?|\bnan\b|\binf\b", re.IGNORECASE)
_PAIR = re.compile(r"""\(\s*(?:'([^']*)'|"([^"]*)")\s*,\s*([^()]*(?:\([^()]*\))?[^()]*?)\s*\)""")
_DATE = re.compile(r"^(\d{4})-(\d{2})-(\d{2})(?:$|[T ])")


class InputError(Exception):
    """A problem with an input file that the user can fix; the message says where."""


# --------------------------------------------------------------------------------------------------
# small parsers (same conventions as timeseries.py)
# --------------------------------------------------------------------------------------------------
def is_na(value):
    return value is None or str(value).strip() in ("", "NA", "NaN", "nan", "None", "[]")


def to_float(value):
    """Finite float or None."""
    try:
        number = float(str(value).strip())
    except (TypeError, ValueError):
        return None
    return number if math.isfinite(number) else None


def parse_float_list(text):
    """Numbers in a Python list repr or a numpy array repr (whitespace/newline separated)."""
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
    """Names of a Python list repr, a numpy array repr (`['A' 'B']`) or plain `A B` text."""
    if is_na(text):
        return []
    quoted = [a or b for a, b in _QUOTED.findall(str(text))]
    if quoted:
        return quoted
    return str(text).strip().strip("[]").replace(",", " ").split()


def parse_summarized(text):
    """`[('Omicron', 0.7), ('Delta', 0.3)]` (also with `np.float64(0.3)` values) -> [(name, frac)]."""
    if is_na(text):
        return []
    pairs = []
    for match in _PAIR.finditer(str(text)):
        name = match.group(1) if match.group(1) is not None else match.group(2)
        numbers = _NUMBER.findall(match.group(3))
        value = to_float(numbers[-1]) if numbers else None
        if value is not None:
            pairs.append((name, value))
    return pairs


def parse_date(text):
    """ISO `YYYY-MM-DD` (a time part is tolerated and dropped) -> `datetime.date`, or None."""
    match = _DATE.match(str(text or "").strip())
    if not match:
        return None
    try:
        return dt.date(int(match.group(1)), int(match.group(2)), int(match.group(3)))
    except ValueError:
        return None


def site_key(name):
    """Matching key of a site name: trimmed, inner whitespace collapsed, case-folded."""
    return " ".join(str(name or "").split()).casefold()


def _add(mapping, name, value):
    mapping[name] = mapping.get(name, 0.0) + value


# --------------------------------------------------------------------------------------------------
# reading the inputs
# --------------------------------------------------------------------------------------------------
def read_text(path, warnings):
    """File text; UTF-8 with or without BOM, falling back to Windows-1252 (old Excel) with a warning."""
    try:
        with open(path, "rb") as handle:
            raw = handle.read()
    except OSError as err:
        raise InputError(f"cannot read {path}: {err.strerror or err}") from err
    try:
        return raw.decode("utf-8-sig")
    except UnicodeDecodeError:
        warnings.append(
            f"{path}: not valid UTF-8, read as Windows-1252; save it as 'CSV UTF-8' in Excel to be safe"
        )
        return raw.decode("cp1252", errors="replace")


def sniff_delimiter(first_line):
    """Comma unless the header clearly uses semicolons or tabs (Excel in some locales)."""
    counts = {d: first_line.count(d) for d in (",", ";", "\t")}
    best = max(counts, key=lambda d: counts[d])
    return best if counts[best] > 0 else ","


def normalise_column(name):
    return "_".join(str(name or "").strip().lower().split())


def read_sites(path, warnings=None):
    """Read the sites file -> list of aggregates `{name, sites: [{name, key, weight, line}]}`.

    Raises `InputError` (with the line number) for an empty file, unknown or missing columns, blank
    names, bad weights and duplicate (aggregate, site) rows.
    """
    warnings = warnings if warnings is not None else []
    text = read_text(path, warnings)
    lines = text.splitlines()
    first = next((ln for ln in lines if ln.strip().strip(",;\t")), None)
    if first is None:
        raise InputError(
            f"{path}: the sites file is empty; it needs a header line "
            f"'aggregate_name,site_name,weight' and one row per site"
        )
    delimiter = sniff_delimiter(first)
    reader = csv.reader(text.splitlines(), delimiter=delimiter)
    header = None
    header_line = 0
    rows = []
    for cells in reader:
        line_no = reader.line_num
        if not any(c.strip() for c in cells):
            continue
        if header is None:
            header = [normalise_column(c) for c in cells]
            header_line = line_no
            continue
        rows.append((line_no, cells))
    unknown = [c for c in header if c and c not in SITES_COLUMNS]
    if unknown:
        hints = []
        for col in unknown:
            close = difflib.get_close_matches(col, SITES_COLUMNS, n=1, cutoff=0.5)
            hints.append(f"'{col}'" + (f" (did you mean '{close[0]}'?)" if close else ""))
        raise InputError(
            f"{path}: line {header_line}: unknown column {', '.join(hints)}; "
            f"the columns are aggregate_name, site_name and (optional) weight"
        )
    for needed in ("aggregate_name", "site_name"):
        if needed not in header:
            raise InputError(
                f"{path}: line {header_line}: no '{needed}' column (found: "
                f"{', '.join(c for c in header if c) or 'none'}); "
                f"the columns are aggregate_name, site_name and (optional) weight"
            )
    if len({c for c in header if c}) != len([c for c in header if c]):
        raise InputError(f"{path}: line {header_line}: a column name appears twice")
    index = {c: i for i, c in enumerate(header) if c}
    aggregates = {}  # normalised aggregate name -> aggregate
    seen = {}  # (aggregate key, site key) -> line
    for line_no, cells in rows:
        if len(cells) > len(header) and any(c.strip() for c in cells[len(header) :]):
            raise InputError(f"{path}: line {line_no}: more cells than the header has columns")

        def cell(name, cells=cells):
            i = index.get(name)
            return cells[i].strip() if i is not None and i < len(cells) else ""

        agg_name = " ".join(cell("aggregate_name").split())
        site_name = " ".join(cell("site_name").split())
        if not agg_name:
            raise InputError(f"{path}: line {line_no}: aggregate_name is blank")
        if not site_name:
            raise InputError(f"{path}: line {line_no}: site_name is blank (aggregate '{agg_name}')")
        weight_text = cell("weight")
        weight = None
        if weight_text:
            weight = to_float(weight_text)
            if weight is None or weight <= 0:
                raise InputError(
                    f"{path}: line {line_no}: weight '{weight_text}' must be a positive number "
                    f"(plain digits such as 125000, no units or thousands separators) or blank"
                )
        akey = agg_name.casefold()
        skey = site_key(site_name)
        if (akey, skey) in seen:
            raise InputError(
                f"{path}: line {line_no}: site '{site_name}' is listed twice for aggregate "
                f"'{agg_name}' (first on line {seen[(akey, skey)]})"
            )
        seen[(akey, skey)] = line_no
        aggregate = aggregates.get(akey)
        if aggregate is None:
            aggregate = aggregates[akey] = {"name": agg_name, "sites": []}
        elif aggregate["name"] != agg_name:
            warnings.append(
                f"{path}: line {line_no}: aggregate '{agg_name}' is spelled differently from "
                f"'{aggregate['name']}' used earlier; treated as the same aggregate"
            )
        aggregate["sites"].append(
            {"name": site_name, "key": skey, "weight": weight, "line": line_no}
        )
    if not aggregates:
        raise InputError(f"{path}: the sites file has a header but no rows")
    return list(aggregates.values())


def read_aggregate(path):
    """Freyja aggregate -> ({sample: {summ, lin, resid, coverage}}, n_unparsable)."""
    with open(path, newline="", encoding="utf-8-sig") as handle:
        table = list(csv.reader(handle, delimiter="\t"))
    if not table:
        return {}, 0
    header = [cell.strip() for cell in table[0]]
    col = {name: i for i, name in enumerate(header) if name}

    def cell(row, name):
        i = col.get(name)
        return row[i] if i is not None and i < len(row) else ""

    out, malformed = {}, 0
    for line in table[1:]:
        if not line or not any(c.strip() for c in line):
            continue
        sample = line[0].strip()
        summ, lin = {}, {}
        for name, value in parse_summarized(cell(line, "summarized")):
            _add(summ, name, value)
        names = parse_name_list(cell(line, "lineages"))
        values = parse_float_list(cell(line, "abundances"))
        if len(names) == len(values):
            for name, value in zip(names, values):
                _add(lin, name, value)
        if not sample or (not summ and not lin):
            malformed += 1
            continue
        if sample in out:
            continue
        out[sample] = {
            "summ": summ,
            "lin": lin,
            "resid": to_float(cell(line, "resid")),
            "coverage": to_float(cell(line, "coverage")),
        }
    return out, malformed


def read_demix_status(path):
    """demix_status.tsv -> {sample: status}."""
    with open(path, newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = [(f or "").strip().lower() for f in (reader.fieldnames or [])]
        if not fields:
            return {}
        if "sample" not in fields or "status" not in fields:
            raise InputError(f"{path}: needs 'sample' and 'status' columns (found {fields})")
        reader.fieldnames = fields
        return {
            (r.get("sample") or "").strip(): (r.get("status") or "").strip().lower() or "unknown"
            for r in reader
            if (r.get("sample") or "").strip()
        }


def read_samplesheet(path, counts=None):
    """Samplesheet CSV -> {sample: dict(sample_type, date, date_given, site, viral_load)}."""
    counts = counts if counts is not None else {}
    meta = {}
    with open(path, newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        reader.fieldnames = [(name or "").strip().lower() for name in (reader.fieldnames or [])]
        if "sample" not in reader.fieldnames:
            raise InputError(f"{path}: no 'sample' column (found {reader.fieldnames})")
        for record in reader:
            sample = (record.get("sample") or "").strip()
            if not sample:
                continue
            if sample in meta:
                counts["duplicate_samplesheet_rows"] = (
                    counts.get("duplicate_samplesheet_rows", 0) + 1
                )
                continue
            raw_date = (record.get("collection_date") or "").strip()
            load_text = (record.get("viral_load") or "").strip()
            load = to_float(load_text) if load_text else None
            if load_text and (load is None or load <= 0):
                counts["invalid_viral_load"] = counts.get("invalid_viral_load", 0) + 1
                load = None
            meta[sample] = {
                "sample_type": (record.get("sample_type") or "sample").strip().lower() or "sample",
                "date": parse_date(raw_date),
                "date_given": bool(raw_date),
                "site": " ".join((record.get("site") or "").split()),
                "viral_load": load,
            }
    return meta


# --------------------------------------------------------------------------------------------------
# periods
# --------------------------------------------------------------------------------------------------
def period_of(date, period):
    """`(label, start_date)`: ISO week `2026-W41` starting on its Monday, or the day itself."""
    if period == "day":
        return date.isoformat(), date
    year, week, _ = date.isocalendar()
    return f"{year}-W{week:02d}", date - dt.timedelta(days=date.weekday())


# --------------------------------------------------------------------------------------------------
# selection
# --------------------------------------------------------------------------------------------------
def select_samples(agg, status, meta, period):
    """Split the samplesheet samples into contributing ones and everything else.

    Returns `(contrib, others, undated, counts)`:
    `contrib[(site_key, (label, start))]` -> list of sample records (contributing);
    `others[(site_key, period)]` -> list of `(status, viral_load)` of dated, non-negative samples that
    cannot contribute (demix not ok, or no aggregate row); `undated[site_key]` -> number of
    non-negative samples with a site but no valid date.
    """
    contrib, others, undated = {}, {}, {}
    counts = {
        "in_aggregate": len(agg),
        "negative_control": 0,
        "not_ok": 0,
        "no_site": 0,
        "no_date": 0,
        "used": 0,
        "aggregate_rows_without_metadata": sum(1 for s in agg if s not in meta),
    }
    for sample in sorted(meta):
        info = meta[sample]
        if info["sample_type"] == "negative":
            counts["negative_control"] += 1
            continue
        state = status.get(sample, "no_status")
        if state == "ok" and sample not in agg:
            state = "not_in_aggregate"
        ok = state == "ok"
        if not ok:
            counts["not_ok"] += 1
        key = site_key(info["site"])
        if not key:
            if ok:
                counts["no_site"] += 1
            continue
        if info["date"] is None:
            undated[key] = undated.get(key, 0) + 1
            if ok:
                counts["no_date"] += 1
            continue
        per = period_of(info["date"], period)
        if ok:
            record = dict(agg[sample], sample=sample, load=info["viral_load"])
            contrib.setdefault((key, per), []).append(record)
            counts["used"] += 1
        else:
            others.setdefault((key, per), []).append((state, info["viral_load"]))
    return contrib, others, undated, counts


# --------------------------------------------------------------------------------------------------
# combining
# --------------------------------------------------------------------------------------------------
def weighted_vector(items):
    """`[(weight, {name: value})]` -> {name: sum(w * v) / sum(w)} (a missing name counts as 0)."""
    total = sum(w for w, _ in items)
    out = {}
    for w, vec in items:
        for name, value in vec.items():
            out[name] = out.get(name, 0.0) + w * value
    return {name: value / total for name, value in out.items()}


def weighted_scalar(items):
    """`[(weight, value or None)]` -> weighted mean over the values that exist, else None."""
    known = [(w, v) for w, v in items if v is not None]
    total = sum(w for w, _ in known)
    return sum(w * v for w, v in known) / total if known else None


def combine(units):
    """Weighted combination of units `{w, summ, lin, resid, coverage}` -> same keys (no `w`)."""
    return {
        "summ": weighted_vector([(u["w"], u["summ"]) for u in units]),
        "lin": weighted_vector([(u["w"], u["lin"]) for u in units]),
        "resid": weighted_scalar([(u["w"], u["resid"]) for u in units]),
        "coverage": weighted_scalar([(u["w"], u["coverage"]) for u in units]),
    }


def combine_site(samples):
    """Replicates of one site and period -> one unit; also the site's mean viral load (or None).

    Samples are weighted by viral load if every one has a load, else equally.
    """
    loads = [s["load"] for s in samples]
    all_loads = all(load is not None for load in loads)
    units = [dict(s, w=s["load"] if all_loads else 1.0) for s in samples]
    result = combine(units)
    result["n"] = len(samples)
    result["load"] = sum(loads) / len(loads) if all_loads else None
    return result


def choose_basis(mode, units):
    """`(kind, label)` for the site weights. `units` carry `load` and `site_weight` (or None)."""
    have_load = [u["load"] is not None for u in units]
    have_weight = [u["site_weight"] is not None for u in units]
    all_load, all_weight = all(have_load), all(have_weight)
    if mode == "equal":
        return "equal", "equal"
    if mode == "viral_load":
        return ("viral_load", "viral_load") if all_load else ("equal", "equal (viral_load missing)")
    if mode == "site_weight":
        return (
            ("site_weight", "site_weight")
            if all_weight
            else ("equal", "equal (site_weight missing)")
        )
    # auto
    if all_load:
        return "viral_load", "viral_load"
    note = " (viral_load missing)" if any(have_load) else ""
    if all_weight:
        return "site_weight", f"site_weight{note}"
    if any(have_weight) and not note:
        note = " (site_weight missing)"
    return "equal", f"equal{note}"


def weight_fraction(kind, included, missing, members_weight):
    """Included weight / total weight, or None when a missing site's weight is unknown.

    `included` and `missing` are lists of site keys. For `viral_load` a missing site's load is the mean
    load of its samples in that period that could not contribute (all must have one).
    """
    if not missing:
        return 1.0
    if kind == "equal":
        return len(included) / (len(included) + len(missing))
    if kind == "site_weight":
        total_inc = sum(members_weight[k] for k in included)
        weights = [members_weight.get(k) for k in missing]
        if any(w is None for w in weights):
            return None
        return total_inc / (total_inc + sum(weights))
    return None  # viral_load: handled by the caller, which has the loads


def missing_reason(key, per, others, undated):
    """Short reason a member site has no contributing sample in the period."""
    bad = others.get((key, per), [])
    if bad:
        return "not_ok:" + min({state for state, _ in bad})
    if undated.get(key):
        return "no_date"
    return "no_sample"


def aggregate_one(aggregate, periods, contrib, others, undated, weights_mode, min_sites):
    """All rows (dicts) of one aggregate; also the number of periods dropped by `min_sites`."""
    members = aggregate["sites"]
    by_key = {m["key"]: m for m in members}
    rows, dropped = [], 0
    for per in sorted(periods, key=lambda p: p[1]):
        included = [m for m in members if (m["key"], per) in contrib]
        if not included:
            continue
        if len(included) < min_sites:
            dropped += 1
            continue
        states = {m["key"]: combine_site(contrib[(m["key"], per)]) for m in included}
        for m in included:
            states[m["key"]]["site_weight"] = m["weight"]
        kind, label = choose_basis(weights_mode, list(states.values()))
        for m in included:
            st = states[m["key"]]
            st["w"] = {"viral_load": st["load"], "site_weight": m["weight"], "equal": 1.0}[kind]
        result = combine([states[m["key"]] for m in included])
        missing = [m for m in members if (m["key"], per) not in contrib]
        reasons = [missing_reason(m["key"], per, others, undated) for m in missing]
        inc_keys = [m["key"] for m in included]
        miss_keys = [m["key"] for m in missing]
        if kind == "viral_load":
            fraction = 1.0 if not missing else None
            if missing:
                loads = []
                for k in miss_keys:
                    ls = [load for _, load in others.get((k, per), [])]
                    loads.append(
                        sum(ls) / len(ls) if ls and all(x is not None for x in ls) else None
                    )
                if all(x is not None for x in loads):
                    total_inc = sum(states[k]["w"] for k in inc_keys)
                    fraction = total_inc / (total_inc + sum(loads))
        else:
            weights = {k: by_key[k]["weight"] for k in by_key}
            fraction = weight_fraction(kind, inc_keys, miss_keys, weights)
        rows.append(
            {
                "aggregate": aggregate["name"],
                "period": per[0],
                "period_start": per[1].isoformat(),
                "result": result,
                "included": [m["name"] for m in included],
                "missing": [m["name"] for m in missing],
                "reasons": reasons,
                "basis": label,
                "fraction": fraction,
                "n_samples": sum(st["n"] for st in states.values()),
            }
        )
    return rows, dropped


def build_aggregates(aggregates, agg, status, meta, period="week", weights="auto", min_sites=1):
    """Run the algorithm; returns `(rows, notes)` where notes carries the counts for stderr."""
    contrib, others, undated, counts = select_samples(agg, status, meta, period)
    all_rows, dropped = [], {}
    for aggregate in aggregates:
        keys = {m["key"] for m in aggregate["sites"]}
        periods = {per for (key, per) in contrib if key in keys}
        rows, n_dropped = aggregate_one(
            aggregate, periods, contrib, others, undated, weights, min_sites
        )
        all_rows.extend(rows)
        dropped[aggregate["name"]] = n_dropped
    return all_rows, counts, dropped


# --------------------------------------------------------------------------------------------------
# writing
# --------------------------------------------------------------------------------------------------
def fmt(value):
    """Eight decimals; `NA` for None."""
    return NA if value is None else f"{value:.{PLACES}f}"


def sorted_items(vec):
    return sorted(((n, v) for n, v in vec.items() if v > 0), key=lambda kv: (-kv[1], kv[0]))


def format_summarized(vec):
    pairs = ", ".join(
        f"({name!r}, {float(f'{value:.{PLACES}f}')!r})" for name, value in sorted_items(vec)
    )
    return f"[{pairs}]"


def write_outputs(rows, out_path, meta_path):
    with open(out_path, "w", newline="") as handle:
        handle.write(AGG_HEADER + "\n")
        for r in rows:
            res = r["result"]
            lin = sorted_items(res["lin"])
            cells = [
                f"{r['aggregate']}_{r['period']}",
                format_summarized(res["summ"]),
                " ".join(name for name, _ in lin),
                " ".join(fmt(value) for _, value in lin),
                fmt(res["resid"]),
                fmt(res["coverage"]),
            ]
            handle.write("\t".join(cells) + "\n")
    with open(meta_path, "w", newline="") as handle:
        handle.write("\t".join(META_COLUMNS) + "\n")
        for r in rows:
            cells = [
                f"{r['aggregate']}_{r['period']}",
                r["aggregate"],
                r["period"],
                r["period_start"],
                str(len(r["included"])),
                ";".join(r["included"]),
                ";".join(r["missing"]),
                ";".join(r["reasons"]),
                r["basis"],
                fmt(r["fraction"]),
                str(r["n_samples"]),
            ]
            handle.write("\t".join(cells) + "\n")


# --------------------------------------------------------------------------------------------------
def parse_args(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--aggregated", help="freyja_aggregated.tsv")
    parser.add_argument("--demix-status", help="demix_status.tsv")
    parser.add_argument("--samplesheet", help="summary/samplesheet.csv")
    parser.add_argument(
        "--sites", required=True, help="sites file: aggregate_name,site_name[,weight]"
    )
    parser.add_argument("--period", choices=["week", "day"], default="week")
    parser.add_argument(
        "--weights", choices=["auto", "viral_load", "equal", "site_weight"], default="auto"
    )
    parser.add_argument("--min-sites", type=int, default=1, help="omit periods with fewer sites")
    parser.add_argument("--out", help="site_aggregates_freyja.tsv")
    parser.add_argument("--meta-out", help="site_aggregates_meta.tsv")
    parser.add_argument(
        "--check-sites",
        action="store_true",
        help="only validate the sites file (and, with --samplesheet, warn about unknown sites)",
    )
    args = parser.parse_args(argv)
    if args.min_sites < 1:
        parser.error("--min-sites must be >= 1")
    if not args.check_sites:
        needed = ["aggregated", "demix_status", "samplesheet", "out", "meta_out"]
        missing = [f"--{n.replace('_', '-')}" for n in needed if not getattr(args, n)]
        if missing:
            parser.error("the following arguments are required: " + ", ".join(missing))
    return args


def unknown_sites(aggregates, meta):
    """Names (sites-file spelling, once each) that no samplesheet row uses as its site."""
    used = {site_key(m["site"]) for m in meta.values()}
    unknown = []
    for aggregate in aggregates:
        for site in aggregate["sites"]:
            if site["key"] not in used and site["name"] not in unknown:
                unknown.append(site["name"])
    return unknown


def check_sites(aggregates, warnings, samplesheet=None):
    """`--check-sites`: report what was read; unknown sites (given a samplesheet) are warnings."""
    for warning in warnings:
        print(f"site_aggregate.py: warning: {warning}", file=sys.stderr)
    if samplesheet:
        for name in unknown_sites(aggregates, read_samplesheet(samplesheet)):
            print(
                f"site_aggregate.py: warning: site '{name}' in the sites file matches no sample's "
                f"'site' (check the spelling)",
                file=sys.stderr,
            )
    summary = "; ".join(
        f"{a['name']}: {', '.join(s['name'] for s in a['sites'])}" for a in aggregates
    )
    print(f"site_aggregate.py: sites file OK, {len(aggregates)} aggregate(s) ({summary})")
    return 0


def main(argv=None):
    args = parse_args(argv)
    warnings = []
    counts = {}
    try:
        aggregates = read_sites(args.sites, warnings)
        if args.check_sites:
            return check_sites(aggregates, warnings, args.samplesheet)
        meta = read_samplesheet(args.samplesheet, counts)
        status = read_demix_status(args.demix_status)
        agg, malformed = read_aggregate(args.aggregated)
    except (OSError, InputError, csv.Error) as err:
        print(f"site_aggregate.py: error: {err}", file=sys.stderr)
        return 1
    rows, sel, dropped = build_aggregates(
        aggregates, agg, status, meta, args.period, args.weights, args.min_sites
    )
    write_outputs(rows, args.out, args.meta_out)
    unknown = unknown_sites(aggregates, meta)
    for warning in warnings:
        print(f"site_aggregate.py: warning: {warning}", file=sys.stderr)
    for name in unknown:
        print(
            f"site_aggregate.py: warning: site '{name}' in the sites file matches no sample's "
            f"'site' (check the spelling)",
            file=sys.stderr,
        )
    for aggregate in aggregates:
        if not any(r["aggregate"] == aggregate["name"] for r in rows):
            print(
                f"site_aggregate.py: warning: aggregate '{aggregate['name']}' has no usable "
                f"period (no contributing samples"
                + (f", or fewer than --min-sites {args.min_sites}" if args.min_sites > 1 else "")
                + ")",
                file=sys.stderr,
            )
    n_dropped = sum(dropped.values())
    notes = ", ".join(f"{k}={v}" for k, v in counts.items())
    print(
        f"site_aggregate.py: {len(aggregates)} aggregate(s), {len(rows)} aggregate-period row(s) "
        f"written ({args.period}); {sel['used']} of {sel['in_aggregate']} aggregated samples "
        f"contributed; not contributing: negative_controls={sel['negative_control']}, "
        f"not_ok={sel['not_ok']}, no_site={sel['no_site']}, no_date={sel['no_date']}, "
        f"aggregate_rows_without_samplesheet_row={sel['aggregate_rows_without_metadata']}, "
        f"unparsable_aggregate_rows={malformed}; periods omitted by --min-sites "
        f"{args.min_sites}: {n_dropped}" + (f"; samplesheet notes: {notes}" if notes else ""),
        file=sys.stderr,
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
