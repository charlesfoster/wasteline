#!/usr/bin/env python3
"""Assemble `report_data.json`, the single machine-readable input of the HTML report.

Inputs are the per-file contracts of docs/PLAN.md section 6a plus `run_info.json` (provenance) and
the per-sample Freyja `*.depths` files. Only `--run-info` is required; every other input is
optional and its absence degrades the report instead of failing the build (the report is also
built in `demix` and `aggregate` modes, where many of these files do not exist). The schema is
documented in docs/report_data.md.

Standard library only.
"""

import argparse
import ast
import csv
import json
import math
import re
import sys
from collections import Counter, defaultdict
from datetime import UTC, date, datetime
from pathlib import Path

SCHEMA_VERSION = 1
TRACK_POINTS = 600
NA_TOKENS = {"", "na", "nan", "n/a", "none", "null"}
TOP_FLAG_IDS = 8  # sample ids listed in a run flag message

QC_COLUMNS = (
    "sample",
    "sample_type",
    "collection_date",
    "site",
    "raw_reads",
    "trimmed_reads",
    "mapped_reads",
    "mapped_pct",
    "mean_depth",
    "depth_sd",
    "median_depth",
    "pct_genome_ge_1x",
    "pct_genome_ge_10x",
    "pct_genome_ge_min_depth",
    "amplicons_total",
    "amplicons_dropped",
    "qc_status",
    "qc_reasons",
    "demix_status",
    "demix_reason",
    "consensus_length",
    "consensus_n_pct",
    "barcode_version",
    # appended later (older qc_summary.tsv files lack them: they load as null)
    "primer_trimmed_pct",
    "primer_scheme_check",
    "primer_scheme_best",
    "primer_scheme_enrichment",
)
PRIMER_SCHEME_CHECKS = ("ok", "mismatch", "uninformative")
INT_COLUMNS = {
    "raw_reads",
    "trimmed_reads",
    "mapped_reads",
    "amplicons_total",
    "amplicons_dropped",
    "consensus_length",
}
FLOAT_COLUMNS = {
    "mapped_pct",
    "mean_depth",
    "depth_sd",
    "median_depth",
    "pct_genome_ge_1x",
    "pct_genome_ge_10x",
    "pct_genome_ge_min_depth",
    "consensus_n_pct",
    "primer_trimmed_pct",
    "primer_scheme_enrichment",
}
QC_STATUSES = ("PASS", "WARN", "FAIL")
DEMIX_STATUSES = ("ok", "failed", "skipped_low_coverage", "missing")
NEGATIVE_ID_RE = re.compile(r"NEG_|NC_")
LABEL_SUFFIXES = (".variants.tsv", ".variants", ".demix", ".tsv")

# Params copied into `params_summary` (in this order) when present in run_info.params.
SUMMARY_PARAMS = (
    "mode",
    "protocol",
    "skip_primer_trimming",
    "duplicates",
    "primer_scheme",
    "primer_bed",
    "barcodes",
    "freyja_eps",
    "freyja_covcut",
    "freyja_depthcutoff",
    "freyja_varthresh",
    "freyja_adapt",
    "freyja_autoadapt",
    "freyja_relaxedmrca",
    "freyja_solver",
    "freyja_confirmedonly",
    "min_genome_coverage_for_demix",
    "amplicon_min_depth",
    "consensus_min_depth",
    "consensus_freq",
    "indel_freq",
    "ntc_max_mapped_reads",
    "ntc_max_genome_coverage_pct",
    "qc_pass_min_pct",
    "qc_fail_below_pct",
    "mutations_of_interest",
    "compare_barcodes",
    "freyja_boot",
    "annotate_variants",
    "skip_consensus",
    "qc_min_primer_trimmed_pct",
    "dropout_impact",
    "dropout_min_abundance",
    "unexplained_mutations",
    "unexplained_min_freq",
    "unexplained_min_excess",
    "unexplained_min_depth",
    "unexplained_min_alt_depth",
    "unexplained_min_recurrence",
    "unexplained_max_pvalue",
    "timeseries",
    "timeseries_min_abundance",
    "freyja_growth",
    "freyja_growth_days",
    "freyja_growth_nboots",
    "freyja_growth_seed",
    "freyja_growth_serial_interval",
    "freyja_growth_grthresh",
    "freyja_growth_min_abundance",
    "freyja_growth_min_dates",
    "freyja_growth_min_span_days",
    "freyja_growth_scale_by_viral_load",
    "freyja_dash",
)
# Thresholds that only appear in report_data.thresholds when the run recorded them
OPTIONAL_THRESHOLDS = (
    "qc_min_primer_trimmed_pct",
    "dropout_min_abundance",
    "unexplained_min_freq",
    "unexplained_min_excess",
    "unexplained_min_depth",
    "unexplained_min_alt_depth",
    "unexplained_min_recurrence",
    "unexplained_max_pvalue",
)
MAX_TIMESERIES_LINEAGES = 30  # full-lineage time series keep the 30 most abundant lineages
MAX_GROWTH_TEXT = 1000
GIT_KEYS = ("git_commit", "git_sha", "commit_id", "commit", "git_revision", "revision")


# ---------------------------------------------------------------------------------------------
# Small value helpers
# ---------------------------------------------------------------------------------------------
def is_na(value):
    return value is None or str(value).strip().lower() in NA_TOKENS


def to_int(value):
    if is_na(value):
        return None
    try:
        return int(float(str(value).strip().replace(",", "")))
    except (ValueError, OverflowError):
        raise ValueError(f"not an integer: {value!r}") from None


def to_float(value):
    if is_na(value):
        return None
    try:
        number = float(str(value).strip())
    except ValueError:
        raise ValueError(f"not a number: {value!r}") from None
    return number if math.isfinite(number) else None


def to_bool(value):
    return str(value).strip().lower() in {"true", "t", "1", "yes", "y"}


def clean_str(value):
    return None if is_na(value) else str(value).strip()


def round_or_none(value, digits=2):
    return None if value is None else round(value, digits)


def read_tsv(path):
    """Rows of a tab-separated file as dicts (UTF-8, header required). Empty file -> []."""
    with open(path, newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return [{(k or "").strip(): v for k, v in row.items() if k is not None} for row in reader]


class Problems:
    """Collects non-fatal parse problems so they can be surfaced as one run flag per file."""

    def __init__(self):
        self.items = defaultdict(list)

    def add(self, source, message):
        self.items[source].append(message)


# ---------------------------------------------------------------------------------------------
# Per-file loaders
# ---------------------------------------------------------------------------------------------
def load_qc(path, problems):
    samples = {}
    for number, row in enumerate(read_tsv(path), start=2):
        sample = clean_str(row.get("sample"))
        if sample is None:
            problems.add("qc_summary.tsv", f"line {number}: empty sample id")
            continue
        record = {}
        for column in QC_COLUMNS:
            raw = row.get(column)
            try:
                if column in INT_COLUMNS:
                    record[column] = to_int(raw)
                elif column in FLOAT_COLUMNS:
                    record[column] = round_or_none(to_float(raw), 4)
                else:
                    record[column] = clean_str(raw)
            except ValueError as err:
                problems.add("qc_summary.tsv", f"{sample}.{column}: {err}")
                record[column] = None
        record["sample"] = sample
        if sample in samples:
            problems.add("qc_summary.tsv", f"duplicate sample {sample!r} (last row wins)")
        status = (record["qc_status"] or "").upper()
        record["qc_status"] = status if status in QC_STATUSES else None
        record["sample_type"] = (record["sample_type"] or "sample").lower()
        check = (record["primer_scheme_check"] or "").lower()
        record["primer_scheme_check"] = check if check in PRIMER_SCHEME_CHECKS else None
        reasons = record["qc_reasons"]
        record["qc_reasons"] = (
            [r.strip() for r in reasons.split(";") if r.strip()] if reasons else []
        )
        samples[sample] = record
    return samples


def load_amplicons(path, problems):
    by_sample = defaultdict(list)
    for number, row in enumerate(read_tsv(path), start=2):
        sample = clean_str(row.get("sample"))
        if sample is None:
            continue
        try:
            by_sample[sample].append(
                {
                    "amplicon": clean_str(row.get("amplicon")),
                    "start": to_int(row.get("start")),
                    "end": to_int(row.get("end")),
                    "mean_depth": round_or_none(to_float(row.get("mean_depth")), 2),
                    "dropped": to_bool(row.get("dropped")),
                }
            )
        except ValueError as err:
            problems.add("amplicon_coverage.tsv", f"line {number}: {err}")
    for rows in by_sample.values():
        rows.sort(key=lambda r: (r["start"] if r["start"] is not None else 0, r["amplicon"] or ""))
    return dict(by_sample)


def load_demix_status(path, problems):
    out = {}
    for row in read_tsv(path):
        sample = clean_str(row.get("sample"))
        if sample is None:
            continue
        status = clean_str(row.get("status"))
        if status is not None and status not in DEMIX_STATUSES:
            problems.add("demix_status.tsv", f"{sample}: unknown status {status!r}")
        out[sample] = {
            "status": status,
            "reason": clean_str(row.get("reason")) or "",
            "barcode_version": clean_str(row.get("barcode_version")),
        }
    return out


_QUOTED = re.compile(r"""'([^']*)'|"([^"]*)\"""")
_NUMBER = re.compile(r"[-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?|\bnan\b|\binf\b", re.IGNORECASE)
_PAIR = re.compile(r"""\(\s*(?:'([^']*)'|"([^"]*)")\s*,\s*([^()]*(?:\([^()]*\))?[^()]*?)\s*\)""")


def parse_name_list(text):
    """Names in a Python list repr *or* a numpy array repr (`['A' 'B']`, no commas).

    Deliberately not ast.literal_eval: for the numpy form it silently concatenates adjacent
    string literals (`['A' 'Q.3']` -> `['AQ.3']`).
    """
    if is_na(text):
        return []
    quoted = [a or b for a, b in _QUOTED.findall(str(text))]
    if quoted:
        return quoted
    # Freyja 2.x writes the `lineages` field of a .demix file as plain space-separated names
    # (`BA.2 AY.4`), which `freyja aggregate` passes through unchanged
    return str(text).strip().strip("[]").replace(",", " ").split()


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


def _finite(value):
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return number if math.isfinite(number) else None


def parse_summarized(text):
    """Freyja `summarized`: a Python list of (group, fraction) tuples, as a string.

    Parsed with ast.literal_eval (never eval); if that fails (e.g. numpy 2 reprs such as
    `np.float64(0.3)`), a conservative regex fallback extracts the pairs.
    """
    if is_na(text):
        return []
    text = str(text).strip()
    pairs = None
    try:
        parsed = ast.literal_eval(text)
        if isinstance(parsed, dict):
            parsed = list(parsed.items())
        pairs = [(str(name), _finite(frac)) for name, frac in parsed]
    except (ValueError, SyntaxError, TypeError, MemoryError, RecursionError):
        pairs = None
    if pairs is None:
        pairs = []
        for match in _PAIR.finditer(text):
            name = match.group(1) if match.group(1) is not None else match.group(2)
            numbers = _NUMBER.findall(match.group(3))
            pairs.append((name, _finite(numbers[-1]) if numbers else None))
    return [[name, round(frac, 5)] for name, frac in pairs if frac is not None]


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
    for suffix in sorted(LABEL_SUFFIXES, key=len, reverse=True):
        if label.endswith(suffix) and len(label) > len(suffix):
            return label[: -len(suffix)]
    return label


def load_aggregated(path, known, problems):
    """freyja_aggregated.tsv -> {sample: {groups, lineages, resid, coverage}}."""
    out = {}
    with open(path, newline="", encoding="utf-8-sig") as handle:
        reader = csv.reader(handle, delimiter="\t", quotechar='"')
        header = next(reader, None)
        if not header:
            return out
        names = [h.strip() for h in header]
        # Freyja writes an unnamed index column; treat the first column as the label whatever
        # its header says.
        index = {name: i for i, name in enumerate(names) if name}
        for row in reader:
            if not row or not row[0].strip():
                continue

            def cell(name, row=row, index=index):
                i = index.get(name)
                return row[i] if i is not None and i < len(row) else None

            sample = resolve_label(row[0], known)
            lineages = parse_name_list(cell("lineages"))
            abundances = parse_float_list(cell("abundances"))
            if len(lineages) != len(abundances):
                problems.add(
                    "freyja_aggregated.tsv",
                    f"{sample}: {len(lineages)} lineages but {len(abundances)} abundances",
                )
            pairs = [[n, round(a, 5)] for n, a in zip(lineages, abundances, strict=False)]
            pairs.sort(key=lambda p: -p[1])
            groups = parse_summarized(cell("summarized"))
            groups.sort(key=lambda p: -p[1])
            if sample in out:
                problems.add("freyja_aggregated.tsv", f"duplicate sample {sample!r}")
            try:
                resid = to_float(cell("resid"))
                coverage = to_float(cell("coverage"))
            except ValueError as err:
                problems.add("freyja_aggregated.tsv", f"{sample}: {err}")
                resid = coverage = None
            out[sample] = {
                "groups": groups,
                "lineages": pairs,
                "resid": round_or_none(resid, 4),
                "coverage": round_or_none(coverage, 2),
            }
    return out


def load_moi(path, problems):
    rows = []
    for number, row in enumerate(read_tsv(path), start=2):
        sample = clean_str(row.get("sample"))
        if sample is None:
            continue
        try:
            rows.append(
                {
                    "sample": sample,
                    "name": clean_str(row.get("name")) or clean_str(row.get("mutation")),
                    "mutation": clean_str(row.get("mutation")),
                    "ref_depth": to_int(row.get("ref_depth")),
                    "alt_depth": to_int(row.get("alt_depth")),
                    "depth": to_int(row.get("depth")),
                    "alt_freq": round_or_none(to_float(row.get("alt_freq")), 4),
                }
            )
        except ValueError as err:
            problems.add("mutations_of_interest.tsv", f"line {number}: {err}")
    return rows


def load_consensus(path, problems):
    rows = []
    for number, row in enumerate(read_tsv(path), start=2):
        sample = clean_str(row.get("sample"))
        if sample is None:
            continue
        try:
            rows.append(
                {
                    "sample": sample,
                    "length": to_int(row.get("length")),
                    "n_count": to_int(row.get("n_count")),
                    "n_pct": round_or_none(to_float(row.get("n_pct")), 3),
                    "ambiguous_count": to_int(row.get("ambiguous_count")),
                    "status": clean_str(row.get("status")),
                    "reason": clean_str(row.get("reason")),  # 7th column; older tables lack it
                }
            )
        except ValueError as err:
            problems.add("consensus_stats.tsv", f"line {number}: {err}")
    return rows


def load_samplesheet(path, problems):
    """Sample metadata (sample_type, collection_date, site) from the samplesheet in effect (CSV or TSV)."""
    with open(path, newline="", encoding="utf-8-sig") as handle:
        text = handle.read()
    first = text.splitlines()[0] if text.strip() else ""
    delimiter = "\t" if first.count("\t") > first.count(",") else ","
    out = {}
    for row in csv.DictReader(text.splitlines(), delimiter=delimiter):
        row = {(k or "").strip(): v for k, v in row.items() if k is not None}
        sample = clean_str(row.get("sample"))
        if sample is None:
            continue
        sample_type = (clean_str(row.get("sample_type")) or "sample").lower()
        out[sample] = {
            "sample_type": sample_type if sample_type in ("sample", "negative") else "sample",
            "collection_date": clean_str(row.get("collection_date")),
            "site": clean_str(row.get("site")),
        }
    return out


def load_boot(path, problems):
    """freyja_boot_summary.tsv: sample, lineage, abundance, lo, hi -> {sample: {lineage: [lo, hi]}}."""
    out = defaultdict(dict)
    for number, row in enumerate(read_tsv(path), start=2):
        sample, lineage = clean_str(row.get("sample")), clean_str(row.get("lineage"))
        if sample is None or lineage is None:
            continue
        try:
            lo, hi = to_float(row.get("lo")), to_float(row.get("hi"))
        except ValueError as err:
            problems.add("freyja_boot_summary.tsv", f"line {number}: {err}")
            continue
        out[sample][lineage] = [round_or_none(lo, 5), round_or_none(hi, 5)]
    return dict(out)


def load_compare_status(path, problems):
    """barcode_compare_status.tsv: sample, status, reason of the second demix."""
    return {
        sample: {
            "status": clean_str(row.get("status")),
            "reason": clean_str(row.get("reason")) or "",
        }
        for row in read_tsv(path)
        if (sample := clean_str(row.get("sample")))
    }


def load_barcode_compare(path, label_a, problems):
    """Long TSV `sample, lineage, abundance_a, abundance_b` (+ optional barcode_version_a/b)."""
    rows = []
    label_b = None
    for number, row in enumerate(read_tsv(path), start=2):
        sample = clean_str(row.get("sample"))
        lineage = clean_str(row.get("lineage"))
        if sample is None or lineage is None:
            continue
        try:
            a = to_float(row.get("abundance_a"))
            b = to_float(row.get("abundance_b"))
        except ValueError as err:
            problems.add("barcode_compare", f"line {number}: {err}")
            continue
        a = a or 0.0
        b = b or 0.0
        label_a = clean_str(row.get("barcode_version_a")) or label_a
        label_b = clean_str(row.get("barcode_version_b")) or label_b
        rows.append(
            {
                "sample": sample,
                "lineage": lineage,
                "a": round(a, 5),
                "b": round(b, 5),
                "delta": round(b - a, 5),
            }
        )
    rows.sort(key=lambda r: (-abs(r["delta"]), r["sample"], r["lineage"]))
    return {"label_a": label_a or "A", "label_b": label_b or "B", "rows": rows}


# ---------------------------------------------------------------------------------------------
# Dropout impact, unexplained mutations, time series, growth rates
# ---------------------------------------------------------------------------------------------
def _parse_affected(text):
    """`XBB.1.5:9,BA.2:3` -> [{lineage, n_sites}] (NA -> [])."""
    if is_na(text):
        return []
    out = []
    for part in str(text).split(","):
        name, sep, count = part.strip().rpartition(":")
        if not sep or not name:
            continue
        try:
            out.append({"lineage": name, "n_sites": int(count)})
        except ValueError:
            continue
    return out


def load_dropout_impact(impact_path, summary_path, problems):
    """dropout_impact.tsv (+ summary) -> {sample: {summary, amplicons}} or None when neither exists."""
    if impact_path is None and summary_path is None:
        return None
    out = {}
    if summary_path is not None:
        for number, row in enumerate(read_tsv(summary_path), start=2):
            sample = clean_str(row.get("sample"))
            if sample is None:
                continue
            try:
                out[sample] = {
                    "summary": {
                        "n_dropped_amplicons": to_int(row.get("n_dropped_amplicons")),
                        "informative_sites_total": to_int(row.get("informative_sites_total")),
                        "informative_sites_in_dropped": to_int(
                            row.get("informative_sites_in_dropped")
                        ),
                        "pct_informative_sites_lost": round_or_none(
                            to_float(row.get("pct_informative_sites_lost")), 2
                        ),
                        "n_dropped_affecting_detected_lineages": to_int(
                            row.get("n_dropped_affecting_detected_lineages")
                        ),
                    },
                    "amplicons": [],
                }
            except ValueError as err:
                problems.add("dropout_impact_summary.tsv", f"line {number}: {err}")
    if impact_path is not None:
        for number, row in enumerate(read_tsv(impact_path), start=2):
            sample = clean_str(row.get("sample"))
            if sample is None:
                continue
            try:
                amplicon = {
                    "amplicon": clean_str(row.get("amplicon")),
                    "start": to_int(row.get("start")),
                    "end": to_int(row.get("end")),
                    "mean_depth": round_or_none(to_float(row.get("mean_depth")), 2),
                    "n_informative_sites": to_int(row.get("n_informative_sites")),
                    "n_sites_in_detected_lineages": to_int(row.get("n_sites_in_detected_lineages")),
                    "affects_detected_lineage": to_bool(row.get("affects_detected_lineage")),
                    "affected_lineages": _parse_affected(row.get("affected_lineages")),
                }
            except ValueError as err:
                problems.add("dropout_impact.tsv", f"line {number}: {err}")
                continue
            out.setdefault(sample, {"summary": None, "amplicons": []})["amplicons"].append(amplicon)
    for entry in out.values():
        entry["amplicons"].sort(key=lambda a: a["start"] if a["start"] is not None else 0)
    return dict(sorted(out.items()))


def load_unexplained(path, problems):
    """unexplained_mutations.tsv (concat table; a per-sample table also loads) -> list of rows."""
    rows = []
    for number, row in enumerate(read_tsv(path), start=2):
        sample = clean_str(row.get("sample"))
        if sample is None:
            continue
        try:
            rows.append(
                {
                    "sample": sample,
                    "position": to_int(row.get("position")),
                    "ref": clean_str(row.get("ref")),
                    "alt": clean_str(row.get("alt")),
                    "observed_freq": round_or_none(to_float(row.get("observed_freq")), 6),
                    "expected_freq": round_or_none(to_float(row.get("expected_freq")), 6),
                    "excess": round_or_none(to_float(row.get("excess")), 6),
                    "depth": to_int(row.get("depth")),
                    "alt_depth": to_int(row.get("alt_depth")),
                    "gene": clean_str(row.get("gene")),
                    "aa_change": clean_str(row.get("aa_change")),
                    "effect": clean_str(row.get("effect")),
                    "n_samples": to_int(row.get("n_samples")),
                    "recurrent": to_bool(row.get("recurrent")),
                    "p_value": to_float(row.get("p_value")),
                }
            )
        except ValueError as err:
            problems.add("unexplained_mutations.tsv", f"line {number}: {err}")
    return rows


def load_timeseries(summarized_path, lineages_path, problems):
    """freyja_timeseries.tsv (+ the _lineages variant) -> the `timeseries` object, or None."""
    if summarized_path is None and lineages_path is None:
        return None

    def rows_of(path):
        out = []
        if path is None:
            return out
        for number, row in enumerate(read_tsv(path), start=2):
            group = clean_str(row.get("lineage_group"))
            when = clean_str(row.get("collection_date"))
            if group is None or when is None:
                continue
            try:
                out.append(
                    {
                        "site": clean_str(row.get("site")),
                        "collection_date": when[:10],
                        "lineage_group": group,
                        "abundance": round_or_none(to_float(row.get("abundance")), 5),
                        "abundance_load_weighted": round_or_none(
                            to_float(row.get("abundance_load_weighted")), 5
                        ),
                        "n_samples": to_int(row.get("n_samples")),
                        "mean_viral_load": round_or_none(to_float(row.get("mean_viral_load")), 2),
                        "total_viral_load": round_or_none(to_float(row.get("total_viral_load")), 2),
                    }
                )
            except ValueError as err:
                problems.add("freyja_timeseries.tsv", f"line {number}: {err}")
        return out

    summarized = rows_of(summarized_path)
    lineages = rows_of(lineages_path)
    truncated = None
    totals = Counter()
    for r in lineages:
        totals[r["lineage_group"]] += r["abundance"] or 0.0
    if len(totals) > MAX_TIMESERIES_LINEAGES:
        keep = {name for name, _ in totals.most_common(MAX_TIMESERIES_LINEAGES)}
        truncated = len(totals)
        lineages = [r for r in lineages if r["lineage_group"] in keep]
    allrows = summarized + lineages
    return {
        "enabled": True,
        "has_viral_load": any(r["abundance_load_weighted"] is not None for r in allrows),
        "min_coverage": None,
        "sites": sorted({r["site"] for r in allrows if r["site"] is not None}),
        "dates": sorted({r["collection_date"] for r in allrows}),
        "summarized": summarized,
        "lineages": lineages,
        "lineages_total_groups": truncated,
    }


_PCT = re.compile(r"[-+]?\d+(?:\.\d+)?(?:[eE][-+]?\d+)?")


def load_growth_rates(csv_path, failed_path, inputs_path, problems):
    """Freyja's relgrowthrate table (+ failure note, + inputs json) -> the `growth_rates` object."""
    if csv_path is None and failed_path is None:
        return None
    reason = None
    if failed_path is not None:
        text = " ".join(Path(failed_path).read_text(encoding="utf-8", errors="replace").split())
        reason = text[:MAX_GROWTH_TEXT] or "no growth rates were produced (no reason recorded)"
    rows = []
    if csv_path is not None:
        with open(csv_path, newline="", encoding="utf-8-sig") as handle:
            reader = csv.reader(handle)
            header = next(reader, None) or []
            low = [h.strip().lower() for h in header]

            def col(word, default):
                for i, name in enumerate(low):
                    if word in name:
                        return i
                return default

            i_lin, i_adv, i_int = col("lineage", 0), col("advantage", 1), col("interval", 2)
            for number, row in enumerate(reader, start=2):
                if not row or not any(c.strip() for c in row):
                    continue
                lineage = clean_str(row[i_lin]) if i_lin < len(row) else None
                if lineage is None:
                    continue
                adv = None
                if i_adv < len(row):
                    found = _PCT.search(row[i_adv])
                    adv = float(found.group()) if found else None
                    if adv is None and not is_na(row[i_adv]):
                        problems.add("freyja_growth_rates.csv", f"line {number}: {row[i_adv]!r}")
                rows.append(
                    {
                        "lineage": lineage,
                        "advantage_pct": adv,
                        "interval_text": clean_str(row[i_int]) if i_int < len(row) else None,
                    }
                )
    inputs = None
    if inputs_path is not None:
        try:
            loaded = json.loads(Path(inputs_path).read_text(encoding="utf-8"))
            inputs = loaded if isinstance(loaded, dict) and loaded else None
        except (OSError, ValueError):
            inputs = None
    status = "ok" if csv_path is not None else "failed"
    return {
        "status": status,
        "reason": reason if status == "failed" else None,
        "inputs": inputs,
        "params": {},
        "rows": rows,
    }


# ---------------------------------------------------------------------------------------------
# Depth files -> binned coverage tracks
# ---------------------------------------------------------------------------------------------
def collect_depth_files(items):
    """Accepts files, directories (searched recursively) or whitespace-separated lists."""
    files = []
    for item in items or []:
        path = Path(item)
        if path.is_dir():
            files.extend(sorted(path.rglob("*.depths")))
        elif path.exists():
            files.append(path)
        else:
            files.extend(Path(part) for part in str(item).split() if Path(part).exists())
    seen, unique = set(), []
    for path in files:
        if str(path) not in seen:
            seen.add(str(path))
            unique.append(path)
    return unique


def sample_from_depth_name(path):
    name = Path(path).name
    for suffix in (".depths", ".depth"):
        if name.endswith(suffix) and len(name) > len(suffix):
            return name[: -len(suffix)]
    return name


def read_depths(path):
    """Per-position depth list (index = position - 1) for the first contig of a depths file."""
    contig = None
    depths = []
    with open(path, encoding="utf-8", errors="replace") as handle:
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 4:
                parts = line.split()
                if len(parts) < 4:
                    continue
            if contig is None:
                contig = parts[0]
            elif parts[0] != contig:
                continue
            try:
                pos = int(parts[1])
                depth = float(parts[3])
            except ValueError:
                continue
            if pos < 1:
                continue
            if pos > len(depths):
                depths.extend([0.0] * (pos - len(depths)))
            depths[pos - 1] = depth
    return contig, depths


def _fmt_depth(value):
    return round(value) if value >= 10 else round(value, 1)


def build_coverage_tracks(depth_files, problems, points=TRACK_POINTS):
    raw = {}
    contigs = Counter()
    for path in depth_files:
        sample = sample_from_depth_name(path)
        try:
            contig, depths = read_depths(path)
        except OSError as err:
            problems.add("depths", f"{path}: {err}")
            continue
        if not depths:
            problems.add("depths", f"{path}: no usable rows")
            continue
        raw[sample] = depths
        contigs[contig] += 1
    if not raw:
        return None
    genome_length = max(len(d) for d in raw.values())
    bin_size = max(1, math.ceil(genome_length / points))
    n_bins = math.ceil(genome_length / bin_size)
    samples = {}
    for sample, depths in sorted(raw.items()):
        track = []
        for i in range(n_bins):
            lo = i * bin_size
            hi = min(lo + bin_size, genome_length)
            chunk = depths[lo:hi]
            track.append(_fmt_depth(sum(chunk) / (hi - lo)) if hi > lo else 0)
        samples[sample] = track
    return {
        "contig": contigs.most_common(1)[0][0],
        "genome_length": genome_length,
        "bin_size": bin_size,
        "n_bins": n_bins,
        "samples": samples,
    }


# ---------------------------------------------------------------------------------------------
# Provenance, flags, time series
# ---------------------------------------------------------------------------------------------
def dig(obj, dotted):
    for part in dotted.split("."):
        if not isinstance(obj, dict) or part not in obj:
            return None
        obj = obj[part]
    return obj


def first_of(obj, keys):
    for key in keys:
        value = dig(obj, key)
        if value not in (None, ""):
            return value
    return None


def build_header(run_info, mode):
    freyja_data = run_info.get("freyja_data") if isinstance(run_info, dict) else None
    freyja_data = freyja_data if isinstance(freyja_data, dict) else {}
    nextflow = first_of(run_info, ("nextflow_version", "nextflow.version", "nextflow"))
    return {
        "pipeline_name": first_of(run_info, ("pipeline_name", "pipeline.name")) or "wasteline",
        "pipeline_version": first_of(run_info, ("pipeline_version", "pipeline.version")),
        "git_commit": first_of(run_info, GIT_KEYS),
        "nextflow_version": nextflow if isinstance(nextflow, str) else None,
        "mode": mode,
        "freyja_version": freyja_data.get("freyja_version") or run_info.get("freyja_version"),
        "barcode_version": freyja_data.get("barcode_version"),
        "barcode_version_iso": freyja_data.get("barcode_version_iso"),
        "barcode_source": freyja_data.get("source"),
        "barcode_sha256": freyja_data.get("barcodes_sha256"),
        "pathogen": freyja_data.get("pathogen"),
        # library protocol (PLAN D13): older run_info files have none (= amplicon)
        "protocol": first_of(run_info, ("protocol", "params.protocol")),
        "primer_trimming": run_info.get("primer_trimming"),
        "duplicates": run_info.get("duplicates"),
        "run_started": first_of(run_info, ("start", "run_started", "timestamp.start")),
        "run_info_generated": first_of(run_info, ("provenance_generated", "generated")),
    }


def build_params_summary(run_info):
    params = run_info.get("params") if isinstance(run_info.get("params"), dict) else {}
    summary = {}
    for key in SUMMARY_PARAMS:
        if key in params and isinstance(params[key], (str, int, float, bool, type(None))):
            summary[key] = params[key]
    if "mode" not in summary and run_info.get("mode"):
        summary["mode"] = run_info["mode"]
    return summary


def num_param(params, key, default=None):
    value = params.get(key) if isinstance(params, dict) else None
    try:
        number = float(value)
    except (TypeError, ValueError):
        return default
    return int(number) if number == int(number) else number


def build_thresholds(run_info):
    params = run_info.get("params") if isinstance(run_info.get("params"), dict) else {}
    thresholds = {
        "amplicon_min_depth": num_param(params, "amplicon_min_depth", 10),
        "consensus_min_depth": num_param(params, "consensus_min_depth", 10),
        "ntc_max_mapped_reads": num_param(params, "ntc_max_mapped_reads"),
        "ntc_max_genome_coverage_pct": num_param(params, "ntc_max_genome_coverage_pct"),
        "min_genome_coverage_for_demix": num_param(params, "min_genome_coverage_for_demix"),
    }
    for key in OPTIONAL_THRESHOLDS:  # only the ones this run recorded
        value = num_param(params, key)
        if value is not None:
            thresholds[key] = value
    return thresholds


def build_time_series(samples, demix, top_n=6):
    """Per site: mean abundance per summarized group per collection date (ok samples only)."""
    points = defaultdict(lambda: defaultdict(list))  # site -> date -> [group dicts]
    totals = Counter()
    count = 0
    for sample in samples:
        info = demix.get(sample["sample"])
        if sample["sample_type"] == "negative" or not info or info["status"] != "ok":
            continue
        if not sample["collection_date"] or not info["groups"]:
            continue
        try:
            when = date.fromisoformat(sample["collection_date"][:10])
        except ValueError:
            continue
        groups = {name: frac for name, frac in info["groups"]}
        points[sample["site"] or ""][when].append(groups)
        count += 1
        for name, frac in groups.items():
            totals[name] += frac
    distinct = {when for dates in points.values() for when in dates}
    if len(distinct) < 2:
        return None
    top = [name for name, _ in totals.most_common() if name != "Other"][:top_n]
    sites = []
    for site in sorted(points):
        dates = sorted(points[site])
        series = []
        for name in top:
            values = []
            for when in dates:
                fracs = [g.get(name, 0.0) for g in points[site][when]]
                values.append(round(sum(fracs) / len(fracs), 5))
            series.append({"group": name, "values": values})
        sites.append(
            {
                "site": site or None,
                "dates": [d.isoformat() for d in dates],
                "n": [len(points[site][d]) for d in dates],
                "series": series,
            }
        )
    return {
        "value": "summarized_group_mean_fraction",
        "groups": top,
        "sites": sites,
        "n_samples": count,
        "n_dates": len(distinct),
    }


def _ids(items):
    items = sorted(items)
    shown = ", ".join(items[:TOP_FLAG_IDS])
    return shown + (f" and {len(items) - TOP_FLAG_IDS} more" if len(items) > TOP_FLAG_IDS else "")


def build_run_flags(ctx):
    """Run-level flags. ctx: dict with samples, demix, provenance, thresholds, availability."""
    flags = []

    def flag(level, code, message, samples=None):
        item = {"level": level, "code": code, "message": message}
        if samples:
            item["samples"] = sorted(samples)
        flags.append(item)

    samples = ctx["samples"]
    demix = ctx["demix"]
    freyja_data = (ctx["provenance"] or {}).get("freyja_data")
    thr = ctx["thresholds"]
    available = ctx["available"]

    if not samples:
        flag("warn", "no_samples", "No samples were found in the inputs; the report is empty.")
    if not available["qc"]:
        flag("info", "qc_not_available", "No QC summary was provided: QC statistics are omitted.")

    # Barcode / Freyja dataset provenance.
    barcode = freyja_data.get("barcode_version") if isinstance(freyja_data, dict) else None
    any_demix = any(info.get("status") is not None for info in demix.values())
    if any_demix and not barcode:
        flag(
            "error",
            "freyja_data_missing",
            "run_info.json has no freyja_data.barcode_version although lineage results are "
            "present: the Freyja barcode dataset used for these results is not recorded.",
        )
    elif barcode and str(barcode).lower() == "unknown":
        flag(
            "warn",
            "barcode_unknown",
            "The Freyja barcode version used for these results is unknown (legacy results).",
        )
    elif barcode and str(barcode).lower() == "custom":
        flag(
            "info",
            "barcode_custom",
            "A custom barcode file was used; its version is not tied to a Freyja release "
            "(see the checksum in the run provenance).",
        )
    if isinstance(freyja_data, dict) and freyja_data.get("source") == "container":
        flag(
            "info",
            "barcode_container",
            "The barcode snapshot bundled in the Freyja container was used; it may be older "
            "than the latest Freyja release.",
        )
    versions = {
        info["barcode_version"]
        for info in demix.values()
        if info.get("barcode_version") and info.get("status") is not None
    }
    if len(versions) > 1:
        flag(
            "error",
            "mixed_barcode_versions",
            "Samples were demixed with different barcode versions: " + ", ".join(sorted(versions)),
        )
    previous = (ctx["provenance"] or {}).get("previous_runs") or []
    unknown_prev = [
        p for p in previous if str(p.get("barcode_version") or "unknown").lower() == "unknown"
    ]
    if unknown_prev:
        flag(
            "warn",
            "previous_run_barcode_unknown",
            f"{len(unknown_prev)} previous run(s) used for these results did not record a barcode "
            "version (legacy Snakemake results).",
        )
    differing = [
        p
        for p in previous
        if p.get("barcode_version")
        and str(p["barcode_version"]).lower() != "unknown"
        and barcode
        and p["barcode_version"] != barcode
    ]
    if differing:
        flag(
            "info",
            "previous_run_barcode_differs",
            f"{len(differing)} previous run(s) used a different barcode version than this run.",
        )

    # Negative controls.
    negatives = [s for s in samples if s["sample_type"] == "negative"]
    if available["qc"] and samples and not negatives:
        flag("info", "no_ntc", "No negative-control sample is present in this run.")
    for s in negatives:
        reasons = s["qc_reasons"]
        too_many = (
            thr["ntc_max_mapped_reads"] is not None
            and s["mapped_reads"] is not None
            and s["mapped_reads"] > thr["ntc_max_mapped_reads"]
        )
        if too_many or "ntc_reads" in reasons:
            detail = ""
            if s["mapped_reads"] is not None:
                limit = thr["ntc_max_mapped_reads"]
                detail = f" ({s['mapped_reads']} mapped reads" + (
                    f", limit {limit})" if limit is not None else ")"
                )
            flag(
                "error",
                "ntc_reads",
                f"Negative control {s['sample']} has reads above the allowed limit{detail}: "
                "possible contamination.",
                [s["sample"]],
            )
        elif s["qc_status"] == "FAIL":
            flag(
                "error",
                "ntc_fail",
                f"Negative control {s['sample']} failed QC ({'; '.join(reasons) or 'no reason given'}).",
                [s["sample"]],
            )

    # Sample QC and demix outcomes.
    real = [s for s in samples if s["sample_type"] != "negative"]
    failed_qc = [s["sample"] for s in real if s["qc_status"] == "FAIL"]
    warn_qc = [s["sample"] for s in real if s["qc_status"] == "WARN"]
    if failed_qc:
        flag(
            "warn",
            "qc_fail",
            f"{len(failed_qc)} sample(s) failed QC: {_ids(failed_qc)}.",
            failed_qc,
        )
    if warn_qc:
        flag(
            "info",
            "qc_warn",
            f"{len(warn_qc)} sample(s) have QC warnings: {_ids(warn_qc)}.",
            warn_qc,
        )
    by_status = defaultdict(list)
    for s in real:
        info = demix.get(s["sample"])
        if info and info["status"]:
            by_status[info["status"]].append(s["sample"])
    if by_status.get("failed"):
        ids = by_status["failed"]
        flag(
            "error", "demix_failed", f"Demixing failed for {len(ids)} sample(s): {_ids(ids)}.", ids
        )
    if by_status.get("skipped_low_coverage"):
        ids = by_status["skipped_low_coverage"]
        flag(
            "warn",
            "demix_skipped",
            f"Demixing was skipped (low genome coverage) for {len(ids)} sample(s): {_ids(ids)}.",
            ids,
        )
    if by_status.get("missing"):
        ids = by_status["missing"]
        flag(
            "warn",
            "demix_missing",
            f"No demix result was found for {len(ids)} sample(s): {_ids(ids)}.",
            ids,
        )
    if real:
        bad = {
            s["sample"]
            for s in real
            if s["qc_status"] == "FAIL"
            or (
                demix.get(s["sample"], {}).get("status")
                in ("failed", "skipped_low_coverage", "missing")
            )
        }
        if len(bad) == len(real):
            flag(
                "error",
                "all_samples_failed",
                "Every non-control sample failed QC or has no usable lineage result.",
            )

    # Consensus genomes: failed (all-N) and built-with-warnings
    consensus_rows = ctx.get("consensus") or []
    failed_cons = sorted(r["sample"] for r in consensus_rows if r["status"] == "failed")
    if failed_cons:
        flag(
            "warn",
            "consensus_failed",
            f"consensus genome failed for {len(failed_cons)} sample(s): {_ids(failed_cons)}; "
            "the consensus is all-N: see the consensus table for the reason.",
            failed_cons,
        )
    warned_cons = sorted(
        r["sample"] for r in consensus_rows if r["status"] == "ok" and r.get("reason")
    )
    if warned_cons:
        flag(
            "info",
            "consensus_warnings",
            f"consensus built with warnings for {len(warned_cons)} sample(s): {_ids(warned_cons)}; "
            "see the consensus table.",
            warned_cons,
        )

    # Primer scheme check (qc_summary.tsv primer_scheme_check): reads line up with another scheme
    mismatched = [s for s in real if s.get("primer_scheme_check") == "mismatch"]
    if mismatched:
        declared = dig(ctx["provenance"] or {}, "params.primer_scheme")
        declared_text = f"the declared scheme {declared}" if declared else "the declared scheme"
        by_best = defaultdict(list)
        for s in mismatched:
            by_best[s.get("primer_scheme_best") or "another scheme"].append(s["sample"])
        ids = [s["sample"] for s in mismatched]
        if len(by_best) == 1:
            ((best, _),) = by_best.items()
            message = (
                f"{len(ids)} sample(s) ({_ids(ids)}): reads line up much better with {best} "
                f"than with {declared_text}: check --primer_scheme."
            )
        else:
            detail = "; ".join(f"{b}: {_ids(v)}" for b, v in sorted(by_best.items()))
            message = (
                f"{len(ids)} sample(s): reads line up much better with a different scheme than "
                f"with {declared_text} ({detail}): check --primer_scheme."
            )
        flag("warn", "primer_scheme_mismatch", message, ids)

    # Amplicon dropout that hides lineage-informative sites of the detected lineages
    affected = sorted(
        sid
        for sid, entry in (ctx.get("dropout") or {}).items()
        if any(a["affects_detected_lineage"] for a in entry["amplicons"])
    )
    if affected:
        flag(
            "warn",
            "dropout_affects_lineages",
            f"Amplicon dropout hides lineage-informative sites for detected lineages in "
            f"{len(affected)} sample(s): {_ids(affected)}; abundance estimates may be biased.",
            affected,
        )

    # Mutations unexplained by the demixed lineages that recur across samples
    recurrent = [r for r in (ctx.get("unexplained") or []) if r["recurrent"]]
    if recurrent:
        keys = {(r["position"], r["alt"]) for r in recurrent}
        first = recurrent[0]
        label = (
            f"{first['gene']}:{first['aa_change']}"
            if first["gene"] and first["aa_change"]
            else f"{first['position']}{first['ref'] or ''}>{first['alt'] or ''}"
        )
        n_first = first["n_samples"] if first["n_samples"] is not None else "several"
        flag(
            "info",
            "unexplained_mutations",
            f"{len(keys)} mutation(s) unexplained by the demixed lineages recur in several samples "
            f"(first: {label} in {n_first} samples); possible novel lineage.",
            {r["sample"] for r in recurrent},
        )

    for source, messages in ctx["problems"].items():
        shown = "; ".join(messages[:3]) + ("; ..." if len(messages) > 3 else "")
        flag(
            "warn", "input_problems", f"{source}: {len(messages)} problem(s) while reading: {shown}"
        )
    for note in ctx["notes"]:
        flag("info", "input_note", note)

    order = {"error": 0, "warn": 1, "info": 2}
    flags.sort(key=lambda f: order[f["level"]])
    return flags


# ---------------------------------------------------------------------------------------------
# Assembly
# ---------------------------------------------------------------------------------------------
def blank_sample(sample_id):
    record = dict.fromkeys(QC_COLUMNS)
    record["sample"] = sample_id
    record["sample_type"] = "sample"
    record["qc_reasons"] = []
    return record


def build(
    run_info,
    qc=None,
    amplicons=None,
    aggregated=None,
    demix_status=None,
    moi=None,
    consensus=None,
    depths=None,
    barcode_compare=None,
    controls=None,
    samplesheet=None,
    boot=None,
    compare_status=None,
    dropout_impact=None,
    unexplained=None,
    timeseries=None,
    growth=None,
    dashboard=None,
    mode=None,
    now=None,
    problems=None,
):
    """Assemble report_data from already-loaded pieces. See `main` for the file-based entry."""
    problems = problems or Problems()
    notes = []
    mode = mode or run_info.get("mode") or dig(run_info, "params.mode") or "unknown"
    freyja_data = (
        run_info.get("freyja_data") if isinstance(run_info.get("freyja_data"), dict) else None
    )
    run_barcode = freyja_data.get("barcode_version") if freyja_data else None
    provenance = dict(run_info)
    provenance.setdefault("previous_runs", [])
    thresholds = build_thresholds(run_info)

    # Sample table.
    if qc is not None:
        samples = {sid: dict(rec) for sid, rec in qc.items()}
    else:
        samples = {}
    for source in (demix_status, aggregated, controls):
        for sid in source or {}:
            samples.setdefault(sid, blank_sample(sid))
    if qc is None:
        # depths is the coverage-track dict: its sample IDs are in depths["samples"], the other keys
        # (contig, genome_length, bin_size, n_bins) are not samples
        for source in (
            amplicons,
            consensus_by_sample(consensus),
            moi_by_sample(moi),
            depths["samples"] if depths else None,
        ):
            for sid in source or {}:
                samples.setdefault(sid, blank_sample(sid))
    # Sample metadata from the samplesheet for samples without a QC row (aggregate mode: no QC exists)
    for sid, rec in samples.items():
        meta = (samplesheet or {}).get(sid)
        if meta and (qc is None or sid not in qc):
            rec["sample_type"] = meta["sample_type"]
            rec["collection_date"] = meta["collection_date"]
            rec["site"] = meta["site"]
        elif (
            (qc is None or sid not in qc)
            and (demix_status or {}).get(sid, {}).get("status") == "missing"
            and NEGATIVE_ID_RE.search(sid)
        ):
            # a control without a demix result is only known from discovery: NEG_/NC_ in the ID
            rec["sample_type"] = "negative"
    known = set(samples)
    for label, source in (
        ("amplicons", amplicons),
        ("consensus", consensus_by_sample(consensus)),
        ("depths", depths and depths["samples"]),
    ):
        stray = sorted(set(source or {}) - known)
        if stray:
            notes.append(
                f"{label}: ignored data for sample(s) not in the sample table: {_ids(stray)}."
            )

    # Demix section.
    demix = {}
    have_demix_inputs = demix_status is not None or aggregated is not None
    for sid, rec in samples.items():
        agg = (aggregated or {}).get(sid)
        status_row = (demix_status or {}).get(sid)
        entry = {
            "status": None,
            "reason": "",
            "barcode_version": None,
            "groups": [],
            "lineages": [],
            "resid": None,
            "coverage": None,
        }
        if agg:
            entry.update(agg)
        if status_row:
            entry["status"] = status_row["status"]
            entry["reason"] = status_row["reason"]
            entry["barcode_version"] = status_row["barcode_version"]
        elif agg or sid in (controls or {}):
            entry["status"] = "ok"
        elif rec["demix_status"]:
            entry["status"] = rec["demix_status"]
            entry["reason"] = rec["demix_reason"] or ""
        elif have_demix_inputs and rec["sample_type"] != "negative":
            entry["status"] = "missing"
            entry["reason"] = "no demix result was found for this sample"
        # A negative control is demixed but not aggregated (PLAN D8): "ok" without lineages in the
        # main aggregate is expected for it, never a failure.
        if (
            entry["status"] == "ok"
            and not entry["groups"]
            and not entry["lineages"]
            and rec["sample_type"] != "negative"
        ):
            entry["status"] = "failed"
            entry["reason"] = entry["reason"] or "demix reported ok but no lineages were returned"
        if boot is not None and sid in boot:
            entry["intervals"] = boot[sid]
        if entry["status"] == "ok" and status_row is None:
            entry["barcode_version"] = run_barcode
        demix[sid] = entry

    # Fill the downstream qc_summary columns from the dedicated files.
    cons_rows = consensus_by_sample(consensus)
    for sid, rec in samples.items():
        info = demix[sid]
        if info["status"] is not None:
            rec["demix_status"] = info["status"]
            rec["demix_reason"] = info["reason"] or None
        if info["barcode_version"]:
            rec["barcode_version"] = info["barcode_version"]
        elif rec["barcode_version"] is None and info["status"] is not None:
            rec["barcode_version"] = run_barcode
        row = cons_rows.get(sid)
        if row:
            if row["length"] is not None:
                rec["consensus_length"] = row["length"]
            if row["n_pct"] is not None:
                rec["consensus_n_pct"] = row["n_pct"]
            rec["consensus_status"] = row.get("status")
            rec["consensus_reason"] = row.get("reason")

    sample_list = sorted(samples.values(), key=lambda r: r["sample"])

    # Optional analyses (each absent unless its input was given)
    params = run_info.get("params") if isinstance(run_info.get("params"), dict) else {}
    dropout_out = None
    if dropout_impact is not None:
        dropout_out = {sid: v for sid, v in dropout_impact.items() if sid in known}
        stray = sorted(set(dropout_impact) - known)
        if stray:
            notes.append(
                f"dropout impact: ignored data for sample(s) not in the sample table: {_ids(stray)}."
            )
        for rec in sample_list:
            entry = dropout_out.get(rec["sample"])
            summ = entry["summary"] if entry else None
            rec["dropout_pct_informative_lost"] = (
                summ["pct_informative_sites_lost"] if summ else None
            )
            n_affecting = summ["n_dropped_affecting_detected_lineages"] if summ else None
            if (
                n_affecting is None
                and entry
                and any(a["affects_detected_lineage"] for a in entry["amplicons"])
            ):
                n_affecting = 1
            rec["dropout_affects_lineages"] = None if n_affecting is None else n_affecting > 0
    unexplained_rows = None
    if unexplained is not None:
        unexplained_rows = [r for r in unexplained if r["sample"] in known]
        counts = Counter(r["sample"] for r in unexplained_rows)
        for rec in sample_list:
            has_demix = demix[rec["sample"]]["status"] == "ok"
            rec["n_unexplained"] = counts.get(rec["sample"], 0) if has_demix else None
    timeseries_out = None
    if timeseries is not None:
        timeseries_out = dict(timeseries)
        timeseries_out["min_coverage"] = num_param(params, "freyja_plot_mincov")
    growth_out = None
    if growth is not None:
        growth_out = dict(growth)
        growth_out["params"] = {
            key: params[f"freyja_growth_{key}"]
            for key in ("days", "nboots", "serial_interval", "seed", "scale_by_viral_load")
            if f"freyja_growth_{key}" in params
        }
    negative_controls = []
    for rec in sample_list:
        sid = rec["sample"]
        if rec["sample_type"] != "negative":
            continue
        info = demix[sid]
        found = (controls or {}).get(sid) or (aggregated or {}).get(sid)
        if info["status"] not in ("ok", "missing") and not found:
            continue
        found = found or {}
        negative_controls.append(
            {
                "sample": sid,
                "status": info["status"],
                "reason": info["reason"],
                "qc_status": rec["qc_status"],
                "qc_reasons": rec["qc_reasons"],
                "mapped_reads": rec["mapped_reads"],
                "pct_genome_ge_1x": rec["pct_genome_ge_1x"],
                "groups": found.get("groups", []),
                "lineages": found.get("lineages", []),
                "coverage": found.get("coverage"),
                "resid": found.get("resid"),
                "excluded_from_aggregate": sid not in (aggregated or {}),
            }
        )
    amplicon_out = {sid: rows for sid, rows in sorted((amplicons or {}).items()) if sid in known}
    for sid, rows in amplicon_out.items():
        rec = samples[sid]
        if rec["amplicons_total"] is None:
            rec["amplicons_total"] = len(rows)
        if rec["amplicons_dropped"] is None:
            rec["amplicons_dropped"] = sum(1 for r in rows if r["dropped"])
    tracks = None
    if depths:
        tracks = dict(depths)
        tracks["samples"] = {k: v for k, v in depths["samples"].items() if k in known}
    moi_rows = [r for r in (moi or []) if r["sample"] in known]
    cons_out = sorted(
        (r for r in (consensus or []) if r["sample"] in known), key=lambda r: r["sample"]
    )

    comparison = None
    if barcode_compare is not None:
        if barcode_compare["label_a"] == "A" and run_barcode:
            barcode_compare["label_a"] = run_barcode
        comparison = barcode_compare
        comparison["rows"] = [r for r in comparison["rows"] if r["sample"] in known]

    available = {
        "qc": qc is not None,
        "amplicons": bool(amplicon_out),
        "demix": any(d["status"] is not None for d in demix.values()),
        "coverage_tracks": bool(tracks and tracks["samples"]),
        "mutations_of_interest": moi is not None,
        "consensus": bool(cons_out),
        "barcode_comparison": bool(comparison and comparison["rows"]),
        "negative_controls": bool(negative_controls),
        "boot_intervals": bool(boot),
        "dropout_impact": bool(dropout_out),
        "unexplained_mutations": unexplained_rows is not None,
        "timeseries": bool(timeseries_out),
        "growth_rates": growth_out is not None,
    }
    time_series = build_time_series(sample_list, demix)
    available["time_series"] = time_series is not None or bool(
        timeseries_out and len(timeseries_out["dates"]) >= 2
    )

    flags = build_run_flags(
        {
            "samples": sample_list,
            "demix": demix,
            "provenance": provenance,
            "thresholds": thresholds,
            "available": available,
            "problems": problems.items,
            "notes": notes,
            "dropout": dropout_out,
            "unexplained": unexplained_rows,
            "consensus": cons_out,
        }
    )
    failed_compare = sorted(
        sid for sid, st in (compare_status or {}).items() if st["status"] == "failed"
    )
    if failed_compare:
        flags.append(
            {
                "level": "warn",
                "code": "barcode_compare_failed",
                "message": f"The second barcode dataset could not be demixed for {len(failed_compare)} sample(s): "
                + ", ".join(failed_compare)
                + ".",
                "samples": failed_compare,
            }
        )
    now = now or datetime.now(UTC)
    return {
        "schema_version": SCHEMA_VERSION,
        "generated": now.strftime("%Y-%m-%dT%H:%M:%SZ"),
        "mode": mode,
        "header": build_header(run_info, mode),
        "provenance": provenance,
        "params_summary": build_params_summary(run_info),
        "thresholds": thresholds,
        "available": available,
        "run_flags": flags,
        "samples": sample_list,
        "demix": {sid: demix[sid] for sid in sorted(demix)},
        "amplicons": amplicon_out,
        "coverage_tracks": tracks,
        "mutations_of_interest": moi_rows if moi is not None else None,
        "consensus": cons_out,
        "barcode_comparison": comparison if available["barcode_comparison"] else None,
        "negative_controls": negative_controls,
        "time_series": time_series,
        "timeseries": timeseries_out,
        "growth_rates": growth_out,
        "dashboard": {"file": dashboard} if dashboard else None,
        "dropout_impact": dropout_out,
        "unexplained_mutations": unexplained_rows,
    }


def consensus_by_sample(rows):
    return {r["sample"]: r for r in rows or []}


def moi_by_sample(rows):
    return {r["sample"]: r for r in rows or []}


def parse_args(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--run-info", required=True, help="pipeline_info/run_info.json")
    parser.add_argument("--qc", help="qc_summary.tsv")
    parser.add_argument("--amplicons", help="amplicon_coverage.tsv")
    parser.add_argument("--aggregated", help="freyja_aggregated.tsv")
    parser.add_argument("--demix-status", help="demix_status.tsv")
    parser.add_argument("--moi", help="mutations_of_interest.tsv")
    parser.add_argument("--consensus", help="consensus_stats.tsv")
    parser.add_argument(
        "--depths",
        nargs="+",
        help="*.depths files, a whitespace-separated list of them, or a directory",
    )
    parser.add_argument(
        "--barcode-compare", help="long TSV: sample, lineage, abundance_a, abundance_b"
    )
    parser.add_argument(
        "--controls-aggregated",
        help="freyja_aggregated_controls.tsv: aggregate of the ok negative controls (excluded from the main aggregate)",
    )
    parser.add_argument(
        "--samplesheet",
        nargs="+",
        help="samplesheet(s) in effect: sample_type, collection_date, site",
    )
    parser.add_argument(
        "--boot", help="freyja_boot_summary.tsv: sample, lineage, abundance, lo, hi"
    )
    parser.add_argument(
        "--compare-status", help="barcode_compare_status.tsv: sample, status, reason"
    )
    parser.add_argument(
        "--dropout-impact", help="dropout_impact.tsv: one row per dropped amplicon (DROPOUT_IMPACT)"
    )
    parser.add_argument(
        "--dropout-summary", help="dropout_impact_summary.tsv: one row per sample (DROPOUT_IMPACT)"
    )
    parser.add_argument(
        "--unexplained",
        help="unexplained_mutations.tsv: the concatenated table (UNEXPLAINED_MUTATIONS_CONCAT)",
    )
    parser.add_argument(
        "--timeseries",
        help="freyja_timeseries.tsv: site, collection_date, lineage_group, abundance, ...",
    )
    parser.add_argument(
        "--timeseries-lineages", help="freyja_timeseries_lineages.tsv (same columns, full lineages)"
    )
    parser.add_argument(
        "--growth-rates", help="freyja_growth_rates.csv (Freyja relgrowthrate table)"
    )
    parser.add_argument(
        "--growth-failed",
        help="freyja_growth_rates.failed.txt: why growth rates could not be computed",
    )
    parser.add_argument("--growth-inputs", help="freyja_growth_inputs.json")
    parser.add_argument(
        "--dashboard-name",
        help="file name of the optional Freyja dashboard (shown as text only; never linked)",
    )
    parser.add_argument("--mode", help="override the mode recorded in run_info.json")
    parser.add_argument("--out", required=True, help="report_data.json")
    parser.add_argument(
        "--write-qc-tsv",
        help="also write the QC table with the demix_*, consensus_* and barcode_version columns filled in",
    )
    return parser.parse_args(argv)


MERGED_QC_COLUMNS = {
    "demix_status": "demix_status",
    "demix_reason": "demix_reason",
    "consensus_length": "consensus_length",
    "consensus_n_pct": "consensus_n_pct",
    "barcode_version": "barcode_version",
}


def write_merged_qc(qc_path, samples, out_path):
    """Copy qc_summary.tsv, filling the columns that only exist downstream (demix status, consensus,
    barcode version) from the assembled sample table. Cells the QC script already filled, and the
    formatting of every other cell, are kept byte for byte."""
    by_sample = {rec["sample"]: rec for rec in samples}
    with open(qc_path, newline="", encoding="utf-8-sig") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader)
        rows = list(reader)
    index = {name: i for i, name in enumerate(header)}
    for row in rows:
        rec = by_sample.get(row[index["sample"]])
        if rec is None:
            continue
        for column, key in MERGED_QC_COLUMNS.items():
            i = index.get(column)
            value = rec.get(key)
            if i is not None and row[i] in ("", "NA") and value not in (None, ""):
                row[i] = f"{value:.2f}" if isinstance(value, float) else str(value)
    with open(out_path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(header)
        writer.writerows(rows)


def optional_file(path):
    """A path argument that Nextflow may stage as an empty placeholder (NO_FILE) -> None."""
    if not path:
        return None
    candidate = Path(path)
    if not candidate.is_file() or candidate.name.startswith("NO_FILE"):
        return None
    return candidate


def main(argv=None):
    args = parse_args(argv)
    try:
        run_info = json.loads(Path(args.run_info).read_text(encoding="utf-8"))
        if not isinstance(run_info, dict):
            raise TypeError("run_info.json must contain a JSON object")
    except (OSError, ValueError, TypeError) as err:
        print(f"build_report_data.py: cannot read run_info.json: {err}", file=sys.stderr)
        return 2
    problems = Problems()
    inputs = {}
    try:
        qc_path = optional_file(args.qc)
        inputs["qc"] = load_qc(qc_path, problems) if qc_path else None
        path = optional_file(args.amplicons)
        inputs["amplicons"] = load_amplicons(path, problems) if path else None
        path = optional_file(args.demix_status)
        inputs["demix_status"] = load_demix_status(path, problems) if path else None
        path = optional_file(args.aggregated)
        known = set(inputs["qc"] or {}) | set(inputs["demix_status"] or {})
        inputs["aggregated"] = load_aggregated(path, known, problems) if path else None
        path = optional_file(args.moi)
        inputs["moi"] = load_moi(path, problems) if path else None
        path = optional_file(args.consensus)
        inputs["consensus"] = load_consensus(path, problems) if path else None
        files = collect_depth_files(args.depths)
        inputs["depths"] = build_coverage_tracks(files, problems) if files else None
        path = optional_file(args.barcode_compare)
        label_a = dig(run_info, "freyja_data.barcode_version")
        inputs["barcode_compare"] = load_barcode_compare(path, label_a, problems) if path else None
        path = optional_file(args.controls_aggregated)
        inputs["controls"] = load_aggregated(path, known, problems) if path else None
        sheets = [optional_file(item) for item in args.samplesheet or []]
        merged = {}
        for sheet in sheets:
            if sheet:
                merged.update(load_samplesheet(sheet, problems))
        inputs["samplesheet"] = merged or None
        path = optional_file(args.boot)
        inputs["boot"] = load_boot(path, problems) if path else None
        path = optional_file(args.compare_status)
        inputs["compare_status"] = load_compare_status(path, problems) if path else None
        impact, summary = optional_file(args.dropout_impact), optional_file(args.dropout_summary)
        inputs["dropout_impact"] = load_dropout_impact(impact, summary, problems)
        path = optional_file(args.unexplained)
        inputs["unexplained"] = load_unexplained(path, problems) if path else None
        inputs["timeseries"] = load_timeseries(
            optional_file(args.timeseries), optional_file(args.timeseries_lineages), problems
        )
        inputs["growth"] = load_growth_rates(
            optional_file(args.growth_rates),
            optional_file(args.growth_failed),
            optional_file(args.growth_inputs),
            problems,
        )
        inputs["dashboard"] = (args.dashboard_name or "").strip() or None
    except (OSError, UnicodeDecodeError, csv.Error) as err:
        print(f"build_report_data.py: {err}", file=sys.stderr)
        return 2
    data = build(run_info, mode=args.mode, problems=problems, **inputs)
    Path(args.out).write_text(
        json.dumps(data, indent=None, separators=(",", ":")) + "\n", encoding="utf-8"
    )
    if args.write_qc_tsv and qc_path:
        write_merged_qc(qc_path, data["samples"], args.write_qc_tsv)
    return 0


if __name__ == "__main__":
    sys.exit(main())
