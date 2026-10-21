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
DEMIX_STATUSES = ("ok", "low_confidence", "failed", "skipped_low_coverage", "missing")
# statuses with a demix result: lineages are shown and used (low_confidence is flagged, PLAN D18)
DEMIXED_STATUSES = ("ok", "low_confidence")
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
    "dropout_min_affect_abundance",
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
    "lineage_groups",
    "lineage_major_threshold",
    "lineage_min_abundance",
    "lineage_summary_zeros",
    "site_aggregation",
    "site_aggregation_period",
    "site_aggregation_weights",
    "site_aggregation_min_sites",
)
# Thresholds that only appear in report_data.thresholds when the run recorded them
OPTIONAL_THRESHOLDS = (
    "qc_min_primer_trimmed_pct",
    "dropout_min_abundance",
    "dropout_min_affect_abundance",
    "unexplained_min_freq",
    "unexplained_min_excess",
    "unexplained_min_depth",
    "unexplained_min_alt_depth",
    "unexplained_min_recurrence",
    "unexplained_max_pvalue",
    "lineage_major_threshold",
    "lineage_min_abundance",
)
UNGROUPED_FLAG_THRESHOLD = 0.02  # an ok sample with `other` >= this raises lineage_groups_ungrouped
SITE_PERIOD_RE = re.compile(r"^(?P<name>.+)_(?P<period>\d{4}-W\d{2}|\d{4}-\d{2}-\d{2})$")
LG_LEFT_ANCHORS = ("genome_coverage_pct", "weight_basis")
LG_RIGHT_ANCHORS = (
    "other",
    "proportion_summary",
    "major_lineages",
    "major_abundances",
    "minor_lineages",
    "minor_abundances",
    "ungrouped_lineages",
    "status_reason",
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
                        # appended columns: absent in older tables -> null
                        "defining_sites_total": to_int(row.get("defining_sites_total")),
                        "defining_sites_lost": to_int(row.get("defining_sites_lost")),
                        "pct_defining_sites_lost": round_or_none(
                            to_float(row.get("pct_defining_sites_lost")), 2
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
                    # appended columns: absent in older tables -> null
                    "n_defining_sites": to_int(row.get("n_defining_sites")),
                    "n_differential_sites": to_int(row.get("n_differential_sites")),
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
                    # last column of newer tables: does any lineage of the dataset model this site?
                    "in_dataset": None
                    if is_na(row.get("in_dataset"))
                    else to_bool(row.get("in_dataset")),
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
# Lineage groups (lineage_summary.csv) and site aggregates
# ---------------------------------------------------------------------------------------------
def read_delimited(path):
    """(header, rows) of a CSV/TSV whose delimiter is sniffed (`,` `;` or tab); UTF-8 with or without BOM."""
    with open(path, newline="", encoding="utf-8-sig") as handle:
        text = handle.read()
    first = text.splitlines()[0] if text.strip() else ""
    delimiter = max(",;\t", key=first.count) if first else ","
    reader = csv.reader(text.splitlines(), delimiter=delimiter)
    header = [h.strip() for h in next(reader, [])]
    rows = []
    for raw in reader:
        if not raw or not any(c.strip() for c in raw):
            continue
        raw = raw + [""] * (len(header) - len(raw))
        rows.append(dict(zip(header, raw, strict=False)))
    return header, rows


def _split_list(text, sep=";"):
    return [] if is_na(text) else [part.strip() for part in str(text).split(sep) if part.strip()]


def _ungrouped_pairs(text):
    """`B.52:0.001 A.1:0.03` -> [[B.52, 0.001], [A.1, 0.03]]."""
    out = []
    for token in str(text or "").split():
        name, sep, value = token.rpartition(":")
        if not sep or not name:
            continue
        try:
            out.append([name, round(float(value), 6)])
        except ValueError:
            continue
    return out


def load_lineage_audit(path, problems):
    """lineage_groups_resolved.tsv -> {groups: [...], dataset_version, groups_sha256}."""
    _, rows = read_delimited(path)
    groups = []
    for row in rows:
        name = (row.get("group") or "").strip()
        if not name:
            continue
        try:
            matched = to_int(row.get("n_lineages_matched"))
        except ValueError as err:
            problems.add("lineage_groups_resolved.tsv", f"{name}: {err}")
            matched = None
        groups.append(
            {
                "name": name,
                "include": _split_list(row.get("include")),
                "exclude": _split_list(row.get("exclude")),
                "in_total": str(row.get("in_total", "yes")).strip().lower() != "no",
                "n_matched": matched,
                "warnings": [
                    w.strip() for w in str(row.get("warnings") or "").split(" | ") if w.strip()
                ],
            }
        )
    first = rows[0] if rows else {}
    return {
        "groups": groups,
        "dataset_version": clean_str(first.get("dataset_version")),
        "groups_sha256": clean_str(first.get("groups_sha256")),
    }


def load_lineage_ungrouped(path):
    """Only what the report needs of lineage_summary.csv: per ok sample, `other` and the lineages behind it.

    The per-sample group table is deliberately not read: it lives in the CSV, not in the report.
    Returns a list of {sample, status, other, ungrouped}."""
    _, rows = read_delimited(path)
    out = []
    for row in rows:
        sample = clean_str(row.get("sample"))
        if sample is None:
            continue
        try:
            other = round_or_none(to_float(row.get("other")), 6)
        except ValueError:
            other = None
        out.append(
            {
                "sample": sample,
                "status": clean_str(row.get("status")),
                "other": other,
                "ungrouped": _ungrouped_pairs(row.get("ungrouped_lineages")),
            }
        )
    return out


def build_lineage_groups(audit, summary_rows, params, threshold=UNGROUPED_FLAG_THRESHOLD):
    """The (small) `lineage_groups` object: audit facts, warnings and the inputs of the ungrouped flag."""
    if audit is None and summary_rows is None:
        return None
    audit = audit or {}
    configured = params.get("lineage_groups") if isinstance(params, dict) else None
    ok_rows = [r for r in (summary_rows or []) if r["status"] in DEMIXED_STATUSES]
    with_ungrouped = [r for r in ok_rows if r["ungrouped"]]
    top = {}
    for r in with_ungrouped:
        for name, value in r["ungrouped"]:
            entry = top.setdefault(name, [name, 0.0, 0])
            entry[1] = max(entry[1], value)
            entry[2] += 1
    top_ungrouped = sorted(top.values(), key=lambda e: (-e[1], e[0]))[:10]
    return {
        "n_groups": len(audit.get("groups", [])),
        "groups_sha256": audit.get("groups_sha256"),
        "dataset_version": audit.get("dataset_version"),
        "file": Path(str(configured)).name if configured else "bundled default",
        "warnings": [
            {"group": g["name"], "text": w} for g in audit.get("groups", []) for w in g["warnings"]
        ],
        "n_samples": len(ok_rows) if summary_rows is not None else None,
        "n_samples_with_ungrouped": len(with_ungrouped),
        "top_ungrouped": top_ungrouped,
        "samples_above_threshold": sorted(
            r["sample"] for r in ok_rows if r["other"] is not None and r["other"] >= threshold
        ),
        "ungrouped_flag_threshold": threshold,
    }


def load_site_aggregates(freyja_path, meta_path, problems):
    """site_aggregates_freyja.tsv + site_aggregates_meta.tsv -> list of aggregate-periods."""
    if freyja_path is None and meta_path is None:
        return None
    meta_rows = []
    if meta_path is not None:
        _, raw = read_delimited(meta_path)
        for number, row in enumerate(raw, start=2):
            label = clean_str(row.get("row_label"))
            if label is None:
                continue
            try:
                included = _split_list(row.get("sites_included"))
                missing = _split_list(row.get("sites_missing"))
                reasons = _split_list(row.get("missing_reasons"))
                meta_rows.append(
                    {
                        "row_label": label,
                        "aggregate_name": clean_str(row.get("aggregate_name")),
                        "period": clean_str(row.get("period")),
                        "period_start": clean_str(row.get("period_start")),
                        "n_sites_included": to_int(row.get("n_sites_included"))
                        if not is_na(row.get("n_sites_included"))
                        else len(included),
                        "sites_included": included,
                        "sites_missing": missing,
                        "missing_reasons": dict(zip(missing, reasons, strict=False)),
                        "weight_basis": clean_str(row.get("weight_basis")),
                        "weight_fraction_included": round_or_none(
                            to_float(row.get("weight_fraction_included")), 4
                        ),
                        "n_samples": to_int(row.get("n_samples")),
                    }
                )
            except ValueError as err:
                problems.add("site_aggregates_meta.tsv", f"line {number}: {err}")
    solutions = {}
    if freyja_path is not None:
        solutions = load_aggregated(freyja_path, {m["row_label"] for m in meta_rows}, problems)
    if meta_path is None:
        for label in solutions:
            match = SITE_PERIOD_RE.match(label)
            meta_rows.append(
                {
                    "row_label": label,
                    "aggregate_name": match["name"] if match else label,
                    "period": match["period"] if match else None,
                    "period_start": None,
                    "n_sites_included": None,
                    "sites_included": [],
                    "sites_missing": [],
                    "missing_reasons": {},
                    "weight_basis": None,
                    "weight_fraction_included": None,
                    "n_samples": None,
                }
            )
    out = []
    for meta in meta_rows:
        sol = solutions.get(meta["row_label"]) or {}
        out.append(
            dict(
                meta,
                groups=sol.get("groups", []),
                lineages=sol.get("lineages", []),
                resid=sol.get("resid"),
                coverage=sol.get("coverage"),
            )
        )
    return out


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


def protocol_counts(run_info):
    """'5 capture, 2 amplicon' from run_info.sample_libraries when the run mixed protocols, else None."""
    libraries = run_info.get("sample_libraries") if isinstance(run_info, dict) else None
    if not isinstance(libraries, dict) or run_info.get("protocol") != "mixed":
        return None
    counts = Counter(
        entry.get("protocol") for entry in libraries.values() if isinstance(entry, dict)
    )
    return ", ".join(f"{n} {name}" for name, n in counts.most_common())


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
        # mixed run (PLAN D20): how many samples per protocol, e.g. "5 capture, 2 amplicon"
        "protocol_counts": protocol_counts(run_info),
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


def build_time_series(samples, demix, top_n=6, include_low_confidence=False):
    """Per site: mean abundance per summarized group per collection date (ok samples only; low_confidence
    ones only with `include_low_confidence`, like the time series tables)."""
    used = DEMIXED_STATUSES if include_low_confidence else ("ok",)
    points = defaultdict(lambda: defaultdict(list))  # site -> date -> [group dicts]
    totals = Counter()
    count = 0
    for sample in samples:
        info = demix.get(sample["sample"])
        if sample["sample_type"] == "negative" or not info or info["status"] not in used:
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
    if by_status.get("low_confidence"):
        ids = by_status["low_confidence"]
        threshold = thr.get("min_genome_coverage_for_demix")
        below = f" below the {threshold:g} % high-confidence threshold" if threshold else ""
        flag(
            "info",
            "demix_low_confidence",
            f"{len(ids)} sample(s) were demixed with low genome coverage{below}: treat their "
            f"lineage estimates with care: {_ids(ids)}.",
            ids,
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

    # Lineage groups: audit warnings and lineages that no group covers
    lg = ctx.get("lineage_groups")
    if lg:
        warnings = lg["warnings"]
        if warnings:
            shown = "; ".join(f"{w['group']}: {w['text']}" for w in warnings[:3])
            more = f" (and {len(warnings) - 3} more)" if len(warnings) > 3 else ""
            flag(
                "warn",
                "lineage_groups_warnings",
                f"{len(warnings)} warning(s) in the lineage groups audit: {shown}{more}.",
            )
        high = lg["samples_above_threshold"]
        if high:
            limit = lg["ungrouped_flag_threshold"]
            top = ", ".join(f"{n} ({v * 100:.1f}%)" for n, v, _ in lg["top_ungrouped"][:3])
            flag(
                "info",
                "lineage_groups_ungrouped",
                f"{len(high)} sample(s) have {limit * 100:g}% or more of their lineages outside every "
                f"lineage group (top ungrouped: {top or 'none listed'}); add them to your lineage "
                "groups file; see summary/lineage_ungrouped_summary.tsv.",
                high,
            )

    # Site aggregates that miss member sites or used a weaker weighting than asked for
    sa = ctx.get("site_aggregates") or []
    incomplete = [
        a
        for a in sa
        if a["sites_missing"]
        or str(a.get("weight_basis") or "").startswith(("equal (", "site_weight ("))
    ]
    if incomplete:
        parts = [
            f"{a['aggregate_name'] or a['row_label']} {a['period'] or ''}".strip()
            + (
                f" (missing {', '.join(a['sites_missing'])})"
                if a["sites_missing"]
                else " (weaker weighting)"
            )
            for a in incomplete[:4]
        ]
        more = f" and {len(incomplete) - 4} more" if len(incomplete) > 4 else ""
        flag(
            "info",
            "site_aggregate_incomplete",
            f"{len(incomplete)} of {len(sa)} site aggregate period(s) are incomplete or use a weaker "
            f"weighting: {'; '.join(parts)}{more}.",
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
        floor = thr.get("dropout_min_affect_abundance")
        floor = 0.05 if floor is None else floor
        flag(
            "warn",
            "dropout_affects_lineages",
            "Amplicon dropout hides defining mutations (those that separate a lineage from its "
            f"parent) of lineages detected at >= {floor * 100:g}% in {len(affected)} sample(s): "
            f"{_ids(affected)}; abundance estimates for those lineages may be biased.",
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
        blind = {(r["position"], r["alt"]) for r in recurrent if r.get("in_dataset") is False}
        blind_text = (
            f" ({len(blind)} of them at sites no lineage in the dataset models)" if blind else ""
        )
        flag(
            "info",
            "unexplained_mutations",
            f"{len(keys)} mutation(s) unexplained by the demixed lineages recur in several samples"
            f"{blind_text} (first: {label} in {n_first} samples); possible novel lineage.",
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
# Anonymisation for sharing (--anonymize paths|full)
# ---------------------------------------------------------------------------------------------
ANONYMIZE_MODES = ("none", "paths", "full")
_HOME_ROOTS = r"(?:/Users|/home|/Volumes|/mnt|/media|/root|~)"
# a path token: quoted (may hold spaces), a home-rooted path (spaces allowed inside directories),
# a Windows drive/UNC path, or a plain absolute / ./ / ../ path without spaces
_PATH_RES = (
    re.compile(r"""(?P<q>['"])(?P<body>(?:(?!(?P=q)).)*[\\/](?:(?!(?P=q)).)*)(?P=q)"""),
    re.compile(rf"""(?<![\w.:~/-]){_HOME_ROOTS}(?:/(?:(?!\s-)[^/'"\n;])+)*/[^\s/'"\n;]+"""),
    re.compile(r"""(?<![\w.:-])[A-Za-z]:[\\/](?:(?:(?!\s-)[^\\/'"\n;])+[\\/])*[^\s\\/'"\n;]*"""),
    re.compile(r"""(?<![\w.:-])\\\\[^\s\\/'"\n;]+(?:\\[^\\/'"\n;]+)*"""),
    re.compile(r"""(?<![\w.:~/<>-])(?:\.{1,2})?/(?=[^\s/])[^\s'"\n;,)\]]*"""),
)
_USER_RES = (
    re.compile(r"(?:/Users|/home)/([^/\\\s'\"]+)"),
    re.compile(r"[A-Za-z]:[\\/]Users[\\/]([^\\/'\"\n;]+)", re.IGNORECASE),
    re.compile(r"(?<![\w.])~([A-Za-z][\w.-]*)/"),
)
_PATH_KEY_RE = re.compile(
    r"(dir|path|file|input|bed|fasta|genome|sheet|results|groups|sites|gff|bam|reads|command|outdir)",
    re.IGNORECASE,
)
_EMAIL_RE = re.compile(r"[\w.+-]+@[\w-]+(?:\.[\w-]+)+")


def _basename(text):
    return re.split(r"[\\/]", text.rstrip("\\/"))[-1] or text


class Anonymizer:
    """Rewrites an assembled report_data dict in place (paths; and sample / site names for `full`)."""

    def __init__(self, mode, data, extra_ids=()):
        self.mode = mode
        self.users = self._find_users(data)
        self.ids, self.sites, self.aggs = {}, {}, {}
        self._text_re = self._loose_re = None
        self._lookup = {}
        if mode == "full":
            ids, sites, aggs = self._collect(data)
            ids |= {str(i) for i in extra_ids}
            self.ids = self._pseudonyms(ids, "S", 3)
            self.sites = self._pseudonyms(sites, "SITE", 1)
            self.aggs = self._pseudonyms(aggs, "AGG", 1)
            names = sorted({**self.ids, **self.sites, **self.aggs}, key=lambda n: (-len(n), n))
            if names:
                alt = "|".join(re.escape(n) for n in names)
                self._text_re = re.compile(
                    rf"(?<![A-Za-z0-9_-])(?:{alt})(?![A-Za-z0-9_-]|\.[A-Za-z0-9])"
                )
            self._lookup = {**self.aggs, **self.sites, **self.ids}
            self._loose_re = re.compile("|".join(re.escape(n) for n in names)) if names else None

    # ---- discovery ----------------------------------------------------------------------
    @staticmethod
    def _iter_strings(obj):
        if isinstance(obj, str):
            yield obj
        elif isinstance(obj, dict):
            for key, value in obj.items():
                yield from Anonymizer._iter_strings(str(key))
                yield from Anonymizer._iter_strings(value)
        elif isinstance(obj, list):
            for item in obj:
                yield from Anonymizer._iter_strings(item)

    def _find_users(self, data):
        users = set()
        for text in self._iter_strings(data):
            for rx in _USER_RES:
                users.update(m.group(1).strip() for m in rx.finditer(text))
        users.discard("")
        return sorted(users, key=lambda u: (-len(u), u))

    @staticmethod
    def _collect(data):
        ids, sites, aggs = set(), set(), set()

        def rows(key):
            return [r for r in (data.get(key) or []) if isinstance(r, dict)]

        for key in ("samples", "mutations_of_interest", "consensus", "negative_controls",
                    "unexplained_mutations"):  # fmt: skip
            ids |= {str(r["sample"]) for r in rows(key) if r.get("sample") is not None}
        comparison = data.get("barcode_comparison") or {}
        ids |= {str(r["sample"]) for r in comparison.get("rows", []) if isinstance(r, dict)}
        for key in ("demix", "amplicons", "dropout_impact"):
            ids |= {str(k) for k in (data.get(key) or {})}
        ids |= {str(k) for k in ((data.get("coverage_tracks") or {}).get("samples") or {})}
        for flag in rows("run_flags"):
            ids |= {str(x) for x in flag.get("samples", [])}
        sites |= {str(r["site"]) for r in rows("samples") if r.get("site")}
        for r in rows("negative_controls"):
            pass
        ts = data.get("time_series") or {}
        sites |= {
            str(x["site"]) for x in ts.get("sites", []) if isinstance(x, dict) and x.get("site")
        }
        t2 = data.get("timeseries") or {}
        sites |= {str(x) for x in t2.get("sites", [])}
        for level in ("summarized", "lineages"):
            sites |= {
                str(r["site"]) for r in t2.get(level, []) if isinstance(r, dict) and r.get("site")
            }
        for agg in data.get("site_aggregates") or []:
            if not isinstance(agg, dict):
                continue
            if agg.get("aggregate_name"):
                aggs.add(str(agg["aggregate_name"]))
            for key in ("sites_included", "sites_missing"):
                sites |= {str(x) for x in agg.get(key, [])}
        return ids, sites, aggs

    @staticmethod
    def _pseudonyms(names, prefix, width):
        names = sorted(names)
        width = max(width, len(str(len(names))))
        return {n: f"{prefix}{i:0{width}d}" for i, n in enumerate(names, start=1)}

    # ---- strings ------------------------------------------------------------------------
    def scrub_text(self, text, substring=False):
        """Paths -> `.../basename`; home directories -> `~`; the user name -> `<user>`."""

        def path(match):
            if re.fullmatch(
                r"(?:/Users|/home)(?:/[^/]+)?|/root|[A-Za-z]:[\\/]Users(?:[\\/][^\\/]+)?",
                match.group(0),
                re.IGNORECASE,
            ):
                return "~"  # a bare home directory
            token = (
                match.group("body")
                if "body" in match.groupdict() and match.group("body")
                else match.group(0)
            )
            quote = match.group("q") if "q" in match.groupdict() and match.group("q") else ""
            if quote and not self._looks_like_path(token):
                return match.group(0)
            return f"{quote}.../{self._safe_base(_basename(token))}{quote}"

        for rx in _PATH_RES:
            text = rx.sub(path, text)
        text = _EMAIL_RE.sub("<email>", text)
        for user in self.users:
            if len(user) >= 1:
                text = re.sub(rf"(?<![\w.-]){re.escape(user)}(?![\w-])", "<user>", text)
            if substring and len(user) >= 4:
                text = re.sub(re.escape(user), "<user>", text, flags=re.IGNORECASE)
        text = re.sub(r"(?<![\w.:-])(?:/Users|/home)(?=$|[\s'\";,)\]])", "~", text)
        text = re.sub(
            r"(?<![\w.:-])[A-Za-z]:[\\/]Users(?=$|[\s'\";,)\]])", "~", text, flags=re.IGNORECASE
        )
        return text

    @staticmethod
    def _looks_like_path(token):
        return bool(
            re.match(r"^(?:~|/|\.{1,2}/|[A-Za-z]:[\\/]|\\\\)", token)
            or ("/" in token and " " not in token and not re.search(r"://", token))
        )

    def _safe_base(self, base):
        """A basename may still contain the user name (`jdoe_samples.csv`): mask long names there too."""
        for user in self.users:
            if len(user) >= 4:
                base = re.sub(re.escape(user), "<user>", base, flags=re.IGNORECASE)
        if self._loose_re:  # file names may embed a sample / site name (`sites_Metro_Region.csv`)
            base = self._loose_re.sub(lambda m: self._lookup[m.group(0)], base)
        return base

    def path_value(self, key, value):
        """A whole config / provenance value: absolute paths (and relative ones under path-like keys)
        are reduced to the base name, everything else goes through `scrub_text`."""
        stripped = value.strip()
        if "://" in stripped and " " not in stripped:
            return value  # URL (container depot, ...)
        absolute = bool(re.match(r"^(?:~|/|\.{1,2}/|[A-Za-z]:[\\/]|\\\\)", stripped))
        relative = (
            ("/" in stripped or "\\" in stripped)
            and key is not None
            and bool(_PATH_KEY_RE.search(str(key)))
        )
        if (absolute or relative) and not (key and "command" in str(key).lower()):
            if (
                " " in stripped
                and absolute
                and not stripped.startswith(
                    ("/Users/", "/home/", "~", "/Volumes/", "/mnt/", "C:", "c:", "\\\\")
                )
            ):
                return self.scrub_text(value)
            return f".../{self._safe_base(_basename(stripped))}"
        return self.scrub_text(value, substring=True)

    def names(self, text):
        """Replace known sample / site / aggregate names inside free text (boundary aware)."""
        if not self._text_re or not isinstance(text, str):
            return text
        return self._text_re.sub(lambda m: self._lookup[m.group(0)], text)

    # ---- walk ---------------------------------------------------------------------------
    def _walk_paths(self, obj, key=None, config=False):
        """Scrub every string leaf; `config` subtrees get the whole-value path rule."""
        if isinstance(obj, str):
            return self.path_value(key, obj) if config else self.scrub_text(obj)
        if isinstance(obj, dict):
            return {
                self.scrub_text(str(k)) if isinstance(k, str) else k: self._walk_paths(v, k, config)
                for k, v in obj.items()
            }
        if isinstance(obj, list):
            return [self._walk_paths(v, key, config) for v in obj]
        return obj

    def apply(self, data):
        out = dict(data)
        for key in data:
            out[key] = self._walk_paths(
                data[key], key, config=key in ("provenance", "params_summary", "header")
            )
        if self.mode == "full":
            self._apply_names(out)
        out["anonymized"] = self.mode
        return out

    def _apply_names(self, d):
        i, st, ag, nm = self.ids, self.sites, self.aggs, self.names
        for r in d.get("samples") or []:
            r["sample"] = i.get(r["sample"], r["sample"])
            if r.get("site"):
                r["site"] = st.get(r["site"], r["site"])
            for k in ("demix_reason", "consensus_reason"):
                if k in r:
                    r[k] = nm(r[k])
        for key in ("demix", "amplicons", "dropout_impact"):
            if d.get(key):
                d[key] = {i.get(k, k): v for k, v in d[key].items()}
        for v in (d.get("demix") or {}).values():
            v["reason"] = nm(v.get("reason"))
        if d.get("coverage_tracks"):
            d["coverage_tracks"]["samples"] = {
                i.get(k, k): v for k, v in d["coverage_tracks"]["samples"].items()
            }
        for key in (
            "mutations_of_interest",
            "consensus",
            "negative_controls",
            "unexplained_mutations",
        ):
            for r in d.get(key) or []:
                if isinstance(r, dict) and r.get("sample") is not None:
                    r["sample"] = i.get(r["sample"], r["sample"])
                if isinstance(r, dict) and "reason" in r:
                    r["reason"] = nm(r["reason"])
        for r in (d.get("barcode_comparison") or {}).get("rows", []):
            r["sample"] = i.get(r["sample"], r["sample"])
        for f in d.get("run_flags") or []:
            f["message"] = nm(f.get("message"))
            if "samples" in f:
                f["samples"] = sorted(i.get(x, x) for x in f["samples"])
        lg = d.get("lineage_groups")
        if lg and lg.get("samples_above_threshold"):
            lg["samples_above_threshold"] = sorted(
                i.get(x, x) for x in lg["samples_above_threshold"]
            )
        for x in (d.get("time_series") or {}).get("sites", []):
            if x.get("site"):
                x["site"] = st.get(x["site"], x["site"])
        t2 = d.get("timeseries")
        if t2:
            t2["sites"] = [st.get(x, x) for x in t2.get("sites", [])]
            for level in ("summarized", "lineages"):
                for r in t2.get(level, []):
                    if r.get("site"):
                        r["site"] = st.get(r["site"], r["site"])
        for a in d.get("site_aggregates") or []:
            old = a.get("aggregate_name")
            new = ag.get(old, old)
            if old and a.get("row_label", "").startswith(old):
                a["row_label"] = new + a["row_label"][len(old) :]
            a["aggregate_name"] = new
            for k in ("sites_included", "sites_missing"):
                a[k] = [st.get(x, x) for x in a.get(k, [])]
            a["missing_reasons"] = {
                st.get(k, k): nm(v) for k, v in a.get("missing_reasons", {}).items()
            }
        for key in ("growth_rates",):
            if d.get(key) and d[key].get("reason"):
                d[key]["reason"] = nm(d[key]["reason"])
        for w in (lg or {}).get("warnings", []):
            w["text"] = nm(w["text"])
        # provenance: only the command line and path-derived strings can carry names
        prov = d.get("provenance")
        if isinstance(prov, dict):
            self._names_in_provenance(prov)
        for k, v in list((d.get("params_summary") or {}).items()):
            if isinstance(v, str) and v.startswith(".../"):
                d["params_summary"][k] = nm(v)

    def _names_in_provenance(self, prov):
        nm = self.names

        def walk(obj, key=None):
            if isinstance(obj, str):
                if key == "command_line" or obj.startswith(".../") or " .../" in obj:
                    return nm(obj)
                return obj
            if isinstance(obj, dict):
                return {k: walk(v, k) for k, v in obj.items()}
            if isinstance(obj, list):
                return [walk(v, key) for v in obj]
            return obj

        for k in list(prov):
            prov[k] = walk(prov[k], k)

    def mapping_rows(self):
        rows = [("sample", k, v) for k, v in self.ids.items()]
        rows += [("site", k, v) for k, v in self.sites.items()]
        rows += [("aggregate", k, v) for k, v in self.aggs.items()]
        return rows


def write_anonymize_map(path, rows):
    with open(path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(["type", "original", "pseudonym"])
        writer.writerows(rows)


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
    lineage_groups=None,
    site_aggregates=None,
    mode=None,
    now=None,
    problems=None,
    anonymize="none",
    anonymize_map=None,
):
    """Assemble report_data from already-loaded pieces. See `main` for the file-based entry."""
    problems = problems or Problems()
    notes = []
    stray_ids = set()
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
        stray_ids.update(stray)
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
            entry["status"] in DEMIXED_STATUSES
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
        stray_ids.update(stray)
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
            rec["dropout_pct_defining_lost"] = summ["pct_defining_sites_lost"] if summ else None
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
            has_demix = demix[rec["sample"]]["status"] in DEMIXED_STATUSES
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
        if info["status"] not in (*DEMIXED_STATUSES, "missing") and not found:
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
        "lineage_groups": lineage_groups is not None,
        "site_aggregates": bool(site_aggregates),
    }
    time_series = build_time_series(
        sample_list, demix, include_low_confidence=bool(params.get("include_low_confidence"))
    )
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
            "lineage_groups": lineage_groups,
            "site_aggregates": site_aggregates,
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
    data = {
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
        "lineage_groups": lineage_groups,
        "site_aggregates": site_aggregates,
    }
    if anonymize != "none":
        anonymizer = Anonymizer(anonymize, data, extra_ids=stray_ids)
        data = anonymizer.apply(data)
        if anonymize_map is not None:
            anonymize_map.extend(anonymizer.mapping_rows())
    return data


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
    parser.add_argument("--lineage-summary", help="lineage_summary.csv (LINEAGE_SUMMARY)")
    parser.add_argument(
        "--lineage-groups-resolved", help="lineage_groups_resolved.tsv: the groups audit table"
    )
    parser.add_argument("--site-aggregates", help="site_aggregates_freyja.tsv (SITE_AGGREGATE)")
    parser.add_argument("--site-aggregates-meta", help="site_aggregates_meta.tsv (SITE_AGGREGATE)")
    parser.add_argument(
        "--ungrouped-flag-threshold",
        type=float,
        default=UNGROUPED_FLAG_THRESHOLD,
        help="an ok sample whose `other` share is at least this raises lineage_groups_ungrouped",
    )
    parser.add_argument(
        "--anonymize",
        choices=ANONYMIZE_MODES,
        default="none",
        help="none (default): as is. paths: local paths, home directories and the user name are removed. "
        "full: additionally sample IDs and site names become stable pseudonyms (S001, SITE1, ...)",
    )
    parser.add_argument(
        "--anonymize-map-out",
        help="with --anonymize full: write the original -> pseudonym table here (never into the report)",
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
        audit_path = optional_file(args.lineage_groups_resolved)
        audit = load_lineage_audit(audit_path, problems) if audit_path else None
        path = optional_file(args.lineage_summary)
        summary_rows = load_lineage_ungrouped(path) if path else None
        params = run_info.get("params") if isinstance(run_info.get("params"), dict) else {}
        inputs["lineage_groups"] = build_lineage_groups(
            audit, summary_rows, params, args.ungrouped_flag_threshold
        )
        inputs["site_aggregates"] = load_site_aggregates(
            optional_file(args.site_aggregates), optional_file(args.site_aggregates_meta), problems
        )
    except (OSError, UnicodeDecodeError, csv.Error) as err:
        print(f"build_report_data.py: {err}", file=sys.stderr)
        return 2
    mapping = []
    data = build(
        run_info,
        mode=args.mode,
        problems=problems,
        anonymize=args.anonymize,
        anonymize_map=mapping,
        **inputs,
    )
    if args.anonymize_map_out:
        if args.anonymize == "full":
            write_anonymize_map(args.anonymize_map_out, mapping)
        else:
            print(
                "build_report_data.py: --anonymize-map-out is only used with --anonymize full "
                "(no map written)",
                file=sys.stderr,
            )
    Path(args.out).write_text(
        json.dumps(data, indent=None, separators=(",", ":")) + "\n", encoding="utf-8"
    )
    if args.write_qc_tsv and qc_path:
        # the published QC table is a pipeline result, not part of the shared report: it keeps the real IDs
        merged = data["samples"]
        if args.anonymize != "none":
            merged = build(run_info, mode=args.mode, problems=problems, **inputs)["samples"]
        write_merged_qc(qc_path, merged, args.write_qc_tsv)
    return 0


if __name__ == "__main__":
    sys.exit(main())
