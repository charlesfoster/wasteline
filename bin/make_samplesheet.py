#!/usr/bin/env python3
"""Build wasteline samplesheets.

Two uses:

1. FASTQ directory -> `full` mode samplesheet
   make_samplesheet.py <fastq_dir> [--suffix S] [--metadata meta.csv] -o samplesheet.csv

2. Existing results directory -> `demix` or `aggregate` mode samplesheet
   make_samplesheet.py --from-results <dir> --mode demix|aggregate [--targets A,B] -o sheet.csv
   Results of a panel run (summary/targets.tsv, samples/<id>/<target>/freyja/) give one row per (sample, target) and a
   `target` column after `sample`; --targets keeps only those targets. Other layouts give no target column.
   --from-results is repeatable (rows of all the runs in one sheet) and the rows can be filtered (all filters are
   combined with AND): --samples, --include-id, --exclude-id, --site, --date-from, --date-to, --sample-type,
   --status ok; --on-duplicate says what to do when a sample ID is in several runs; --list prints the table of
   candidate rows (selected or not, with the reason) and writes nothing.

Sample IDs are NEVER rewritten (no underscore-to-dash); the only transformation is the removal of
the Illumina `_S<n>_L001_R1_001.fastq.gz` tail for the default suffix, or of the user suffix.
Every problem is collected and reported; on any error nothing is written and the exit code is 1.
"""

from __future__ import annotations

import argparse
import contextlib
import csv
import datetime
import io
import math
import os
import re
import sys
from pathlib import Path

DEFAULT_SUFFIX = "_L001_R1_001.fastq.gz"
DEFAULT_NEGATIVE_PATTERN = "NEG_|NC_"
DEFAULT_EXCLUDE_PATTERN = "Undet"
ID_RE = re.compile(r"^[A-Za-z0-9._-]+$")
DATE_RE = re.compile(r"^\d{4}-\d{2}-\d{2}$")
SAMPLE_TYPES = ("sample", "negative")
FASTQ_EXTS = (".fastq.gz", ".fq.gz")
META_COLUMNS = ("sample_type", "collection_date", "site")
# optional numeric column (e.g. gene copies per litre); written to a sheet only when some sample has a value
VIRAL_LOAD = "viral_load"
# optional per-sample library columns (PLAN D20): written to a sheet only when some sample has a value
LIBRARY_COLUMNS = ("protocol", "primer_scheme")
PROTOCOLS = ("amplicon", "capture")
META_INPUT_COLUMNS = (*META_COLUMNS, VIRAL_LOAD, *LIBRARY_COLUMNS)

FULL_COLUMNS = ["sample", "fastq_1", "fastq_2", *META_COLUMNS]
DEMIX_COLUMNS = ["sample", "variants", "depths", "fastp_json", *META_COLUMNS]
AGGREGATE_COLUMNS = ["sample", "demix", *META_COLUMNS]

# Directories under a results dir that are never legacy sample directories.
# `db` is --db_dir when it sits inside an output directory (the test profiles do that)
RESERVED_DIRS = {"samples", "summary", "pipeline_info", "db"}
# files that mark a directory as a Freyja dataset (<db_dir>/freyja/...), not as a sample
DATASET_MARKERS = (
    "last_barcode_update.txt",
    "last_barcode_update_other.txt",
    "usher_barcodes.feather",
    "usher_barcodes.csv",
    "lineages.yml",
)


class SamplesheetError(Exception):
    """Collects one or more fatal problems."""

    def __init__(self, errors):
        if isinstance(errors, str):
            errors = [errors]
        self.errors = list(errors)
        super().__init__("; ".join(self.errors))


# ----------------------------------------------------------------------------------------------
# shared helpers
# ----------------------------------------------------------------------------------------------
def valid_viral_load(value: str) -> bool:
    """A positive, finite number."""
    try:
        number = float(value)
    except ValueError:
        return False
    return math.isfinite(number) and number > 0


def valid_date(value: str) -> bool:
    if not DATE_RE.match(value):
        return False
    try:
        datetime.date.fromisoformat(value)
    except ValueError:
        return False
    return True


def read_csv_rows(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with open(path, newline="", encoding="utf-8-sig") as fh:
        text = fh.read()
    # allow tab-separated input by sniffing the header line
    first = text.splitlines()[0] if text.strip() else ""
    delim = "\t" if first.count("\t") > first.count(",") else ","
    reader = csv.DictReader(io.StringIO(text), delimiter=delim)
    header = [h.strip() for h in (reader.fieldnames or [])]
    reader.fieldnames = header
    rows = []
    for row in reader:
        rows.append({k: (v or "").strip() for k, v in row.items() if k is not None})
    return header, rows


def load_metadata(path: Path) -> tuple[dict[str, dict[str, str]], list[str], list[str]]:
    """Read a metadata CSV keyed by `sample`. Returns (meta, errors, warnings)."""
    errors: list[str] = []
    warnings: list[str] = []
    meta: dict[str, dict[str, str]] = {}
    try:
        header, rows = read_csv_rows(path)
    except OSError as exc:
        return {}, [f"cannot read metadata file {path}: {exc}"], []
    if "sample" not in header:
        return {}, [f"metadata file {path} has no 'sample' column (found: {header})"], []
    extra = [h for h in header if h != "sample" and h not in META_INPUT_COLUMNS]
    if extra:
        warnings.append(
            f"metadata columns ignored (not one of {', '.join(META_INPUT_COLUMNS)}): {extra}"
        )
    for i, row in enumerate(rows, start=2):
        sample = row.get("sample", "")
        if not sample:
            errors.append(f"metadata line {i}: empty sample")
            continue
        if sample in meta:
            errors.append(f"metadata line {i}: duplicate sample '{sample}'")
            continue
        entry = {c: row.get(c, "") for c in META_INPUT_COLUMNS}
        if entry[VIRAL_LOAD] and not valid_viral_load(entry[VIRAL_LOAD]):
            errors.append(
                f"metadata sample '{sample}': viral_load '{entry[VIRAL_LOAD]}' "
                "must be a positive number (or blank)"
            )
        if entry["collection_date"] and not valid_date(entry["collection_date"]):
            errors.append(
                f"metadata sample '{sample}': collection_date '{entry['collection_date']}' "
                "is not ISO YYYY-MM-DD"
            )
        if entry["sample_type"] and entry["sample_type"] not in SAMPLE_TYPES:
            errors.append(
                f"metadata sample '{sample}': sample_type '{entry['sample_type']}' "
                f"must be one of {'|'.join(SAMPLE_TYPES)}"
            )
        if entry["protocol"] and entry["protocol"] not in PROTOCOLS:
            errors.append(
                f"metadata sample '{sample}': protocol '{entry['protocol']}' "
                f"must be one of {'|'.join(PROTOCOLS)}"
            )
        meta[sample] = entry
    return meta, errors, warnings


def apply_metadata(records: list[dict[str, str]], meta: dict[str, dict[str, str]]) -> list[str]:
    """Merge metadata into records in place; return warnings."""
    warnings = []
    known = {r["sample"] for r in records}
    unmatched = sorted(set(meta) - known)
    if unmatched:
        warnings.append(f"{len(unmatched)} metadata row(s) match no sample: {', '.join(unmatched)}")
    missing = []
    for rec in records:
        entry = meta.get(rec["sample"])
        if entry is None:
            missing.append(rec["sample"])
            continue
        if entry["collection_date"]:
            rec["collection_date"] = entry["collection_date"]
        if entry["site"]:
            rec["site"] = entry["site"]
        if entry["sample_type"]:
            rec["sample_type"] = entry["sample_type"]
        if entry[VIRAL_LOAD]:
            rec[VIRAL_LOAD] = entry[VIRAL_LOAD]
        for col in LIBRARY_COLUMNS:
            if entry[col]:
                rec[col] = entry[col]
    if missing:
        warnings.append(
            f"{len(missing)} sample(s) have no metadata row: {', '.join(sorted(missing))}"
        )
    return warnings


def make_path(path: str, relative_to: Path | None) -> str:
    """Absolute by default (no symlink resolution); relative to a directory if given."""
    if not path:
        return ""
    absolute = os.path.abspath(path)
    if relative_to is None:
        return absolute
    return os.path.relpath(absolute, start=os.path.abspath(relative_to))


def write_sheet(records: list[dict[str, str]], columns: list[str], out: str) -> None:
    with contextlib.ExitStack() as stack:
        fh = (
            sys.stdout
            if out == "-"
            else stack.enter_context(open(out, "w", newline="", encoding="utf-8"))
        )
        writer = csv.DictWriter(fh, fieldnames=columns, lineterminator="\n", extrasaction="ignore")
        writer.writeheader()
        for rec in records:
            writer.writerow({c: rec.get(c, "") for c in columns})


# ----------------------------------------------------------------------------------------------
# FASTQ discovery
# ----------------------------------------------------------------------------------------------
def r2_suffix(suffix: str) -> str | None:
    """Derive the R2 suffix from the R1 suffix (last `R1` -> `R2`, else last `_1` -> `_2`)."""
    for token, repl in (("R1", "R2"), ("_1", "_2")):
        idx = suffix.rfind(token)
        if idx != -1:
            return suffix[:idx] + repl + suffix[idx + len(token) :]
    return None


def discover_fastqs(
    fastq_dir: Path,
    suffix: str = DEFAULT_SUFFIX,
    negative_pattern: str = DEFAULT_NEGATIVE_PATTERN,
    exclude_pattern: str = DEFAULT_EXCLUDE_PATTERN,
    recursive: bool = False,
) -> tuple[list[dict[str, str]], list[str], list[str]]:
    """Return (records, errors, warnings). Records carry absolute-path fastq_1/fastq_2."""
    errors: list[str] = []
    warnings: list[str] = []
    if not fastq_dir.is_dir():
        return [], [f"FASTQ directory not found: {fastq_dir}"], []
    # Mirror the Snakemake behaviour: a user may hand over the R2 suffix.
    if "R2" in suffix and "R1" not in suffix:
        suffix = suffix.replace("R2", "R1")
    if not suffix.endswith(FASTQ_EXTS):
        return [], [f"--suffix '{suffix}' must end with .fastq.gz or .fq.gz"], []
    suffix2 = r2_suffix(suffix)
    if suffix2 is None:
        return [], [f"cannot derive the R2 suffix from '{suffix}' (no 'R1' or '_1' in it)"], []
    try:
        neg_re = re.compile(negative_pattern)
        excl_re = re.compile(exclude_pattern) if exclude_pattern else None
    except re.error as exc:
        return [], [f"invalid regular expression: {exc}"], []

    is_illumina = suffix == DEFAULT_SUFFIX
    illumina_re = re.compile(r"^(?P<id>.+)_S\d+" + re.escape(suffix) + "$")

    pattern = "*" + suffix
    found = sorted(fastq_dir.rglob(pattern) if recursive else fastq_dir.glob(pattern))
    found = [p for p in found if p.is_file()]

    by_id: dict[str, list[Path]] = {}
    for r1 in found:
        name = r1.name
        if excl_re and excl_re.search(name):
            continue
        if is_illumina:
            m = illumina_re.match(name)
            if not m:
                errors.append(
                    f"{r1}: cannot extract a sample ID (expected '<id>_S<n>{suffix}'); "
                    "use --suffix to describe a different naming scheme"
                )
                continue
            sample = m.group("id")
        else:
            sample = name[: -len(suffix)]
        by_id.setdefault(sample, []).append(r1)

    records = []
    for sample in sorted(by_id):
        files = by_id[sample]
        if len(files) > 1:
            errors.append(
                f"duplicate sample ID '{sample}' from {len(files)} files: "
                + ", ".join(str(f) for f in files)
            )
            continue
        r1 = files[0]
        if not ID_RE.match(sample):
            errors.append(
                f"sample ID '{sample}' (from {r1.name}) does not match ^[A-Za-z0-9._-]+$ "
                "(IDs are never rewritten; rename the file or fix the ID in the sheet)"
            )
            continue
        r2 = r1.with_name(r1.name[: -len(suffix)] + suffix2)
        if not r2.is_file():
            errors.append(f"missing R2 for '{sample}': expected {r2}")
            continue
        records.append(
            {
                "sample": sample,
                "fastq_1": str(r1),
                "fastq_2": str(r2),
                "sample_type": "negative" if neg_re.search(r1.name) else "sample",
                "collection_date": "",
                "site": "",
            }
        )
    if not records and not errors:
        errors.append(f"no samples found matching '*{suffix}' in {fastq_dir}")
    return records, errors, warnings


# ----------------------------------------------------------------------------------------------
# results discovery
# ----------------------------------------------------------------------------------------------
EMPTY_DEMIX_REASON = "empty demix file (failed in the old pipeline)"


def first_existing(*paths: Path) -> Path | None:
    for p in paths:
        if p.is_file():
            return p
    return None


def first_nonempty(*paths: Path) -> Path | None:
    """First existing, non-zero-byte file (a zero-byte one is skipped, never returned)."""
    for p in paths:
        if p.is_file() and p.stat().st_size > 0:
            return p
    return None


def find_fastp_json(directory: Path, sample: str) -> Path | None:
    hit = first_existing(
        directory / f"{sample}.fastp.json",
        directory / f"{sample}_fastp.json",
        directory / f"{sample}.json",
    )
    if hit:
        return hit
    if directory.is_dir():
        jsons = sorted(directory.glob("*.json"))
        if len(jsons) == 1:
            return jsons[0]
    return None


def locate_new(results: Path, sample: str, target: str | None = None) -> dict[str, Path | None]:
    """Files of a sample in samples/<sample>/freyja/, or (panel layout, `target` given) of a
    (sample, target) pair in samples/<sample>/<target>/freyja/; fastp stays per sample."""
    base = results / "samples" / sample
    fj = (base / target if target else base) / "freyja"
    return {
        "variants": first_existing(fj / f"{sample}.variants.tsv", fj / f"{sample}.variants.tsv.gz"),
        "depths": first_existing(fj / f"{sample}.depths", fj / f"{sample}.depths.gz"),
        "demix": first_nonempty(fj / f"{sample}.demix", fj / f"{sample}.demix.tsv"),
        "fastp_json": find_fastp_json(base / "fastp", sample),
        "demix_any": first_existing(fj / f"{sample}.demix", fj / f"{sample}.demix.tsv"),
    }


def locate_legacy(results: Path, sample: str) -> dict[str, Path | None]:
    base = results / sample
    fj = base / "freyja"
    return {
        "variants": first_existing(
            fj / f"{sample}.freyja_variants.tsv",
            base / "variants" / f"{sample}.freyja_variants.tsv",
            base / "variants" / f"{sample}.freyja_variants.tsv.gz",
        ),
        "depths": first_existing(fj / f"{sample}.depths", fj / f"{sample}.depths.gz"),
        "demix": first_nonempty(fj / f"{sample}.demix", fj / f"{sample}.demix.tsv"),
        "fastp_json": first_existing(base / "fastp" / f"{sample}_fastp.json")
        or find_fastp_json(base / "fastp", sample),
        "demix_any": first_existing(fj / f"{sample}.demix", fj / f"{sample}.demix.tsv"),
    }


def load_results_sheet(
    results: Path, target: str | None = None
) -> tuple[dict[str, dict[str, str]], str | None]:
    """The samplesheet the run kept. Panel layout (`target` given): summary/<target>/samplesheet.csv first,
    and only the rows of that target (or without a target) of a sheet that has a target column."""
    rels = ["samplesheet.csv", "summary/samplesheet.csv"]
    if target:
        rels = [f"summary/{target}/samplesheet.csv", "summary/samplesheet.csv", "samplesheet.csv"]
    for rel in rels:
        p = results / rel
        if p.is_file():
            try:
                header, rows = read_csv_rows(p)
            except (OSError, csv.Error):
                continue
            if "sample" not in header:
                continue
            return (
                {
                    r["sample"]: {c: r.get(c, "") for c in META_INPUT_COLUMNS}
                    for r in rows
                    if r["sample"] and (not target or r.get("target", "") in ("", target))
                },
                str(p),
            )
    return {}, None


def carried_viral_load(meta):
    """The viral_load of a results samplesheet row when it is a valid positive number, else blank."""
    value = meta.get(VIRAL_LOAD, "")
    return value if value and valid_viral_load(value) else ""


def skipped_meta(sample, sheet_meta, neg_re):
    """(sample_type, collection_date, site) of a sample that is skipped: from the results' own
    samplesheet, else inferred from the ID (NEG_/NC_ = negative control)."""
    meta = sheet_meta.get(sample, {})
    sample_type = meta.get("sample_type") or ("negative" if neg_re.search(sample) else "sample")
    return (sample_type, meta.get("collection_date", ""), meta.get("site", ""))


def is_dataset_dir(path: Path) -> bool:
    """A `freyja/` dir that holds Freyja data (any dataset file, directly or one level down) is a
    downloaded dataset, not a sample's results."""
    freyja = path / "freyja"
    if not freyja.is_dir():
        return False
    for child in [freyja, *(c for c in freyja.iterdir() if c.is_dir())]:
        if any((child / marker).exists() for marker in DATASET_MARKERS):
            return True
    return False


def is_panel_layout(results: Path) -> bool:
    """Results of a panel run: summary/targets.tsv, or a summary/<target>/ folder beside samples/<id>/<target>/freyja/."""
    summary = results / "summary"
    if (summary / "targets.tsv").is_file():
        return True
    samples = results / "samples"
    if not (summary.is_dir() and samples.is_dir()):
        return False
    names = {d.name for d in summary.iterdir() if d.is_dir()}
    return any(
        (target / "freyja").is_dir()
        for sample in samples.iterdir()
        if sample.is_dir()
        for target in sample.iterdir()
        if target.is_dir() and target.name in names
    )


def discover_panel_results(
    results: Path,
    mode: str,
    neg_re,
    skipped: list[tuple[str, ...]] | None,
    targets: set[str] | None,
) -> tuple[list[dict[str, str]], list[str], list[str], list[str]]:
    """Panel layout (samples/<sample>/<target>/freyja/): one record per (sample, target) with a `target` key.
    Skipped entries are (sample, target, reason, sample_type, collection_date, site)."""
    errors: list[str] = []
    warnings: list[str] = []
    notes: list[str] = [
        f"panel layout: one row per (sample, target) in {results}/samples/<sample>/<target>/"
    ]
    required = ("variants", "depths") if mode == "demix" else ("demix",)
    pairs = []
    for sample_dir in sorted(p for p in (results / "samples").iterdir() if p.is_dir()):
        for target_dir in sorted(p for p in sample_dir.iterdir() if p.is_dir()):
            if (target_dir / "freyja").is_dir():
                pairs.append((sample_dir.name, target_dir.name))
    found = {target for _, target in pairs}
    if targets is not None:
        missing_targets = sorted(targets - found)
        if missing_targets:
            warnings.append(
                f"--targets {', '.join(missing_targets)}: no results for them in {results} "
                f"(targets found: {', '.join(sorted(found)) or 'none'})"
            )
        pairs = [(sample, target) for sample, target in pairs if target in targets]
    sheets: dict[str, dict[str, dict[str, str]]] = {}
    for target in sorted({target for _, target in pairs}):
        sheets[target], sheet_path = load_results_sheet(results, target)
        if sheet_path:
            notes.append(
                f"target {target}: carrying over sample_type/collection_date/site/viral_load/protocol/primer_scheme "
                f"from {sheet_path}"
            )

    records = []
    for sample, target in pairs:
        sheet_meta = sheets[target]
        files = locate_new(results, sample, target)
        if not ID_RE.match(sample):
            reason = "sample ID does not match ^[A-Za-z0-9._-]+$"
            notes.append(f"skipped '{sample}' ({target}): ID does not match ^[A-Za-z0-9._-]+$")
            if skipped is not None:
                skipped.append((sample, target, reason, *skipped_meta(sample, sheet_meta, neg_re)))
            continue
        missing = [k for k in required if files[k] is None]
        if missing:
            if mode == "aggregate" and files.get("demix_any") is not None:
                reason = EMPTY_DEMIX_REASON
                notes.append(
                    f"skipped '{sample}' ({target}, panel layout): {reason}: {files['demix_any']}"
                )
            else:
                reason = f"missing {', '.join(missing)} (panel layout)"
                notes.append(
                    f"skipped '{sample}' ({target}, panel layout): missing {', '.join(missing)}"
                )
            if skipped is not None:
                skipped.append((sample, target, reason, *skipped_meta(sample, sheet_meta, neg_re)))
            continue
        for key in required:
            f = files[key]
            if f is not None and f.name.endswith(".gz"):
                notes.append(
                    f"'{sample}' ({target}): {key} is gzipped ({f}); make sure the pipeline can read .gz"
                )
        meta = sheet_meta.get(sample, {})
        sample_type = meta.get("sample_type") or ("negative" if neg_re.search(sample) else "sample")
        rec = {
            "sample": sample,
            "target": target,
            "sample_type": sample_type,
            "collection_date": meta.get("collection_date", ""),
            "site": meta.get("site", ""),
            VIRAL_LOAD: carried_viral_load(meta),
            **{col: meta.get(col, "") for col in LIBRARY_COLUMNS},
        }
        for key in ("variants", "depths", "demix", "fastp_json"):
            rec[key] = str(files[key]) if files[key] is not None else ""
        records.append(rec)
    if not records:
        errors.append(
            f"no usable samples found in {results} for mode '{mode}'"
            + (f" and target(s) {', '.join(sorted(targets))}" if targets is not None else "")
            + " (looked for samples/<s>/<target>/freyja/)"
        )
    return records, errors, warnings, notes


def discover_results(
    results: Path,
    mode: str,
    negative_pattern: str = DEFAULT_NEGATIVE_PATTERN,
    skipped: list[tuple[str, ...]] | None = None,
    targets: set[str] | None = None,
) -> tuple[list[dict[str, str]], list[str], list[str], list[str]]:
    """Return (records, errors, warnings, notes). `notes` are informational lines.

    Samples that are found but unusable are appended to `skipped` as (sample, reason, sample_type,
    collection_date, site) when a list is given (the type is inferred like for usable samples): missing required files, or (aggregate mode) a zero-byte .demix file, which the old
    pipeline left behind when `freyja demix` failed (`|| touch`). Such a file is never aggregated.

    Results of a panel run (is_panel_layout) give one record per (sample, target), with a `target` key, and
    skipped entries (sample, target, reason, sample_type, collection_date, site); `targets` (a set of names)
    keeps only those targets. For every other layout `targets` is ignored.
    """
    errors: list[str] = []
    warnings: list[str] = []
    notes: list[str] = []
    if not results.is_dir():
        return [], [f"results directory not found: {results}"], [], []
    try:
        neg_re = re.compile(negative_pattern)
    except re.error as exc:
        return [], [f"invalid regular expression: {exc}"], [], []
    if is_panel_layout(results):
        return discover_panel_results(results, mode, neg_re, skipped, targets)
    if targets:
        notes.append(
            f"--targets {', '.join(sorted(targets))} ignored: {results} is not a panel run's results"
        )

    new_dir = results / "samples"
    new_ids = (
        {p.name for p in new_dir.iterdir() if p.is_dir() and (p / "freyja").is_dir()}
        if new_dir.is_dir()
        else set()
    )
    legacy_ids = {
        p.name
        for p in results.iterdir()
        if p.is_dir()
        and p.name not in RESERVED_DIRS
        and not p.name.startswith(".")
        and ((p / "freyja").is_dir() or (p / "variants").is_dir())
        and not is_dataset_dir(p)
    }
    if new_ids and legacy_ids:
        warnings.append(
            f"both result layouts present in {results}: {len(new_ids)} sample(s) in "
            f"samples/<sample>/ (new) and {len(legacy_ids)} in <sample>/ (legacy); "
            "preferring the new layout where a sample exists in both"
        )
    overlap = sorted(new_ids & legacy_ids)
    if overlap:
        warnings.append(
            f"{len(overlap)} sample(s) present in both layouts, using the new one: "
            + ", ".join(overlap)
        )

    sheet_meta, sheet_path = load_results_sheet(results)
    if sheet_path:
        notes.append(
            f"carrying over sample_type/collection_date/site/viral_load/protocol/primer_scheme from {sheet_path}"
        )

    required = ("variants", "depths") if mode == "demix" else ("demix",)
    records = []
    for sample in sorted(new_ids | legacy_ids):
        layout = "new" if sample in new_ids else "legacy"
        files = locate_new(results, sample) if layout == "new" else locate_legacy(results, sample)
        if not ID_RE.match(sample):
            reason = "sample ID does not match ^[A-Za-z0-9._-]+$"
            notes.append(f"skipped '{sample}': ID does not match ^[A-Za-z0-9._-]+$")
            if skipped is not None:
                skipped.append((sample, reason, *skipped_meta(sample, sheet_meta, neg_re)))
            continue
        missing = [k for k in required if files[k] is None]
        if missing:
            if mode == "aggregate" and files.get("demix_any") is not None:
                reason = EMPTY_DEMIX_REASON
                notes.append(
                    f"skipped '{sample}' ({layout} layout): {reason}: {files['demix_any']}"
                )
            else:
                reason = f"missing {', '.join(missing)} ({layout} layout)"
                notes.append(f"skipped '{sample}' ({layout} layout): missing {', '.join(missing)}")
            if skipped is not None:
                skipped.append((sample, reason, *skipped_meta(sample, sheet_meta, neg_re)))
            continue
        for key in required + ("fastp_json",):
            f = files[key]
            if f is not None and f.name.endswith(".gz") and key != "fastp_json":
                notes.append(
                    f"'{sample}': {key} is gzipped ({f}); make sure the pipeline can read .gz"
                )
        meta = sheet_meta.get(sample, {})
        sample_type = meta.get("sample_type") or ("negative" if neg_re.search(sample) else "sample")
        rec = {
            "sample": sample,
            "sample_type": sample_type,
            "collection_date": meta.get("collection_date", ""),
            "site": meta.get("site", ""),
            VIRAL_LOAD: carried_viral_load(meta),
            **{col: meta.get(col, "") for col in LIBRARY_COLUMNS},
        }
        for key in ("variants", "depths", "demix", "fastp_json"):
            rec[key] = str(files[key]) if files[key] is not None else ""
        records.append(rec)
    if not records:
        errors.append(
            f"no usable samples found in {results} for mode '{mode}' "
            "(looked for samples/<s>/freyja/ and <s>/freyja/)"
        )
    return records, errors, warnings, notes


# ----------------------------------------------------------------------------------------------
# row selection for --from-results (filters, duplicate handling, --list)
# ----------------------------------------------------------------------------------------------
ON_DUPLICATE = ("error", "prefix", "first")
FILTER_OPTIONS = (
    "samples",
    "include_id",
    "exclude_id",
    "site",
    "date_from",
    "date_to",
    "sample_type",
    "status",
)
STATUS_FILE = "demix_status.tsv"


def parse_sample_list(spec: str) -> tuple[set[str], str | None]:
    """`ID[,ID...]` or `@file` (one ID per line; blank lines and # comments ignored)."""
    if spec.startswith("@"):
        try:
            text = Path(spec[1:]).read_text(encoding="utf-8")
        except OSError as exc:
            return set(), f"cannot read --samples file {spec[1:]}: {exc}"
        ids = [ln.strip() for ln in text.splitlines()]
        return {i for i in ids if i and not i.startswith("#")}, None
    return {i.strip() for i in spec.split(",") if i.strip()}, None


def build_filters(args: argparse.Namespace) -> tuple[dict, list[str]]:
    """The filters of the command line as a dict (only those given) and the errors in them."""
    flt: dict = {}
    errors: list[str] = []
    if args.samples:
        ids, err = parse_sample_list(args.samples)
        if err:
            errors.append(err)
        flt["samples"] = ids
    for key, option in (("include_id", "--include-id"), ("exclude_id", "--exclude-id")):
        value = getattr(args, key)
        if value:
            try:
                flt[key] = re.compile(value)
            except re.error as exc:
                errors.append(f"invalid regular expression for {option}: {exc}")
    if args.site:
        flt["site"] = {s.strip() for s in args.site.split(",") if s.strip()}
    for key, option in (("date_from", "--date-from"), ("date_to", "--date-to")):
        value = getattr(args, key)
        if value:
            if not valid_date(value):
                errors.append(f"{option} '{value}' is not an ISO YYYY-MM-DD date")
            flt[key] = value
    if args.sample_type:
        flt["sample_type"] = args.sample_type
    if args.status:
        flt["status"] = args.status
    return flt, errors


def filter_reason(rec: dict, flt: dict) -> str | None:
    """Why a row is removed by the filters (the first that fails), or None."""
    sample = rec["sample"]
    if "samples" in flt and sample not in flt["samples"]:
        return "not in --samples"
    if "include_id" in flt and not flt["include_id"].search(sample):
        return "ID does not match --include-id"
    if "exclude_id" in flt and flt["exclude_id"].search(sample):
        return "ID matches --exclude-id"
    if "site" in flt and rec.get("site", "") not in flt["site"]:
        return f"site '{rec.get('site', '')}' not in --site" if rec.get("site") else "no site"
    date = rec.get("collection_date", "")
    if ("date_from" in flt or "date_to" in flt) and not date:
        return "no collection_date"
    if "date_from" in flt and date < flt["date_from"]:
        return f"date {date} before --date-from"
    if "date_to" in flt and date > flt["date_to"]:
        return f"date {date} after --date-to"
    if "sample_type" in flt and rec.get("sample_type", "") != flt["sample_type"]:
        return f"sample_type {rec.get('sample_type', '')}"
    if "status" in flt and rec.get("_status") not in (None, flt["status"]):
        return f"demix status {rec['_status']}"
    return None


class StatusLookup:
    """demix_status.tsv of a results directory (summary/<target>/ for panel results), read once."""

    def __init__(self) -> None:
        self.cache: dict[tuple[str, str], dict[str, str] | None] = {}
        self.missing: dict[str, int] = {}

    def get(self, results: Path, target: str, sample: str) -> str | None:
        key = (str(results), target)
        if key not in self.cache:
            path = results / "summary" / target / STATUS_FILE if target else None
            if path is None:
                path = results / "summary" / STATUS_FILE
            table = None
            if path.is_file():
                try:
                    header, rows = read_csv_rows(path)
                except (OSError, csv.Error):
                    header, rows = [], []
                if "sample" in header and "status" in header:
                    table = {r["sample"]: r["status"] for r in rows if r.get("sample")}
            self.cache[key] = table
        table = self.cache[key]
        status = None if table is None else table.get(sample)
        if status is None:
            self.missing[str(results)] = self.missing.get(str(results), 0) + 1
        return status


def run_labels(dirs: list[str]) -> dict[str, str]:
    """A short label per results directory: its name, or its full path when two share a name."""
    names = [os.path.basename(os.path.abspath(d)) for d in dirs]
    return {d: (n if names.count(n) == 1 else os.path.abspath(d)) for d, n in zip(dirs, names)}


def resolve_duplicates(
    records: list[dict], policy: str
) -> tuple[list[dict], list[tuple[dict, str]], list[str]]:
    """The same (sample, target) in several results directories. Returns (kept, dropped with reason, errors)."""
    groups: dict[tuple[str, str], list[dict]] = {}
    for rec in records:
        groups.setdefault((rec["sample"], rec.get("target", "")), []).append(rec)
    kept: list[dict] = []
    dropped: list[tuple[dict, str]] = []
    errors: list[str] = []
    for rec in records:
        group = groups[(rec["sample"], rec.get("target", ""))]
        if len(group) == 1:
            kept.append(rec)
            continue
        runs = ", ".join(r["_run"] for r in group)
        where = f" (target {rec['target']})" if rec.get("target") else ""
        if policy == "error":
            dropped.append((rec, f"duplicate ID, also in {runs}; set --on-duplicate prefix|first"))
            if rec is group[0]:
                errors.append(
                    f"sample '{rec['sample']}'{where} is in several results directories: {runs} "
                    "(use --on-duplicate prefix|first)"
                )
        elif policy == "first":
            if rec is group[0]:
                kept.append(rec)
            else:
                dropped.append((rec, f"duplicate of the one in {group[0]['_run']}"))
        else:  # prefix
            names = [os.path.basename(os.path.abspath(r["_dir"])) for r in group]
            name = os.path.basename(os.path.abspath(rec["_dir"]))
            if rec is group[0]:
                if len(set(names)) < len(names):
                    errors.append(
                        f"--on-duplicate prefix: the directories of '{rec['sample']}' share a name "
                        f"({runs}); cannot build distinct prefixes"
                    )
                if not ID_RE.match(name):
                    errors.append(
                        f"--on-duplicate prefix: directory name '{name}' has characters outside "
                        "[A-Za-z0-9._-]"
                    )
            rec["_renamed_from"] = rec["sample"]
            rec["sample"] = f"{name}__{rec['sample']}"
            kept.append(rec)
    seen = set()
    for rec in kept:
        key = (rec["sample"], rec.get("target", ""))
        if key in seen:
            errors.append(f"sample ID '{rec['sample']}' is not unique after --on-duplicate prefix")
        seen.add(key)
    return kept, dropped, errors


def select_rows(
    records: list[dict], unusable: list[dict], flt: dict, policy: str
) -> tuple[list[dict], list[dict], list[str], list[str]]:
    """Apply the filters, then the duplicate policy. Returns (selected, table, errors, notes).
    `table` has every candidate row (selected or not, with the reason) for --list."""
    notes: list[str] = []
    after_filters, table = [], []
    n_removed = 0
    for rec in records:
        reason = filter_reason(rec, flt)
        if reason:
            n_removed += 1
            table.append({**rec, "_selected": "no", "_reason": reason})
        else:
            after_filters.append(rec)
    if flt:
        notes.append(f"{n_removed} sample(s) removed by filters")
    selected, dropped, errors = resolve_duplicates(after_filters, policy)
    for rec, reason in dropped:
        table.append({**rec, "_selected": "no", "_reason": reason})
    if dropped and policy == "first":
        notes.append(f"{len(dropped)} duplicate row(s) dropped (--on-duplicate first)")
    for rec in selected:
        reason = f"renamed from {rec['_renamed_from']}" if "_renamed_from" in rec else ""
        table.append({**rec, "_selected": "yes", "_reason": reason})
    for rec in unusable:
        table.append({**rec, "_selected": "no", "_reason": f"unusable: {rec['_reason']}"})
    return selected, table, errors, notes


def print_table(table: list[dict], run_order: dict[str, int], with_target: bool) -> None:
    columns = [("run", "_run"), ("sample", "_shown")]
    if with_target:
        columns.append(("target", "target"))
    columns += [
        ("sample_type", "sample_type"),
        ("site", "site"),
        ("collection_date", "collection_date"),
        ("demix_status", "_status"),
        ("selected", "_selected"),
        ("reason", "_reason"),
    ]
    table = sorted(
        table, key=lambda r: (run_order[r["_run"]], r.get("_shown", ""), r.get("target", ""))
    )
    cells = [[h for h, _ in columns]]
    for rec in table:
        cells.append([str(rec.get(key) or "-") for _, key in columns])
    widths = [max(len(row[i]) for row in cells) for i in range(len(columns))]
    for i, row in enumerate(cells):
        line = "  ".join(c.ljust(w) for c, w in zip(row, widths)).rstrip()
        print(line)
        if i == 0:
            print("  ".join("-" * w for w in widths))


# ----------------------------------------------------------------------------------------------
# CLI
# ----------------------------------------------------------------------------------------------
def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="Create a wasteline samplesheet from a FASTQ directory or existing results.",
    )
    p.add_argument("fastq_dir", nargs="?", help="directory containing R1/R2 FASTQ files")
    p.add_argument("--suffix", default=DEFAULT_SUFFIX, help="R1 suffix (default: %(default)s)")
    p.add_argument(
        "--negative-pattern",
        default=DEFAULT_NEGATIVE_PATTERN,
        help="regex; matching file names are sample_type=negative (default: %(default)s)",
    )
    p.add_argument(
        "--exclude-pattern",
        default=DEFAULT_EXCLUDE_PATTERN,
        help="regex; matching file names are skipped entirely (default: %(default)s)",
    )
    p.add_argument(
        "--metadata",
        help="CSV keyed by 'sample' with collection_date, site, sample_type, viral_load",
    )
    p.add_argument(
        "--recursive", action="store_true", help="search the FASTQ directory recursively"
    )
    p.add_argument(
        "--from-results",
        metavar="DIR",
        action="append",
        help="build a demix/aggregate sheet from an existing results directory "
        "(repeatable: the rows of all of them go into one sheet)",
    )
    p.add_argument("--mode", choices=["demix", "aggregate"], help="sheet type for --from-results")
    p.add_argument(
        "--relative",
        action="store_true",
        help="write paths relative to the output file's directory (default: absolute)",
    )
    p.add_argument("-o", "--output", default="-", help="output CSV (default: stdout)")
    p.add_argument(
        "--targets",
        metavar="A,B",
        help="with --from-results on the results of a panel run: only these targets (comma-separated)",
    )
    p.add_argument(
        "--skipped-out",
        metavar="TSV",
        help="with --from-results: write (sample, reason) of every sample that was found but "
        "skipped (missing files, zero-byte legacy .demix) to this TSV",
    )
    g = p.add_argument_group(
        "row filters for --from-results (combined with AND)",
        "Site, date and sample type are the values carried over from the results' samplesheet.",
    )
    g.add_argument(
        "--samples",
        metavar="ID[,ID...]|@FILE",
        help="keep only these sample IDs (exact; @FILE has one ID per line); all targets of a sample",
    )
    g.add_argument("--include-id", metavar="REGEX", help="keep only IDs matching this regex")
    g.add_argument("--exclude-id", metavar="REGEX", help="drop IDs matching this regex")
    g.add_argument("--site", metavar="A,B", help="keep only these sites (exact, comma-separated)")
    g.add_argument("--date-from", metavar="YYYY-MM-DD", help="keep collection_date >= this date")
    g.add_argument("--date-to", metavar="YYYY-MM-DD", help="keep collection_date <= this date")
    g.add_argument("--sample-type", choices=list(SAMPLE_TYPES), help="keep only this sample_type")
    g.add_argument(
        "--status",
        choices=["ok"],
        help="keep only samples whose demix_status.tsv says ok (drops low_confidence, skipped_*, ...); "
        "samples with no status file are kept",
    )
    p.add_argument(
        "--on-duplicate",
        choices=list(ON_DUPLICATE),
        default="error",
        help="the same sample ID in two --from-results directories: error (default, exit 1), "
        "prefix (rename to <dirname>__<id> in the sample column only) or first (keep the first)",
    )
    p.add_argument(
        "--list",
        action="store_true",
        help="with --from-results: print the table of candidate rows (selected or not, with the "
        "reason) to stdout and write nothing",
    )
    return p


def discover_all(args, dirs, wanted, status_lookup):
    """Discover every --from-results directory. Returns (records, unusable, skipped_rows, panel_flags,
    errors, warnings, notes); the rows carry _run/_dir/_shown/_status for the selection."""
    labels = run_labels(dirs)
    records: list[dict] = []
    unusable: list[dict] = []
    skipped_rows: list[tuple[str, ...]] = []
    errors: list[str] = []
    warnings: list[str] = []
    notes: list[str] = []
    panel_flags = []
    for d in dirs:
        results = Path(d)
        skipped: list[tuple[str, ...]] = []
        recs, errs, warns, nts = discover_results(
            results, args.mode, args.negative_pattern, skipped, wanted
        )
        errors += errs
        warnings += warns
        notes += nts
        panel = is_panel_layout(results) if results.is_dir() else False
        panel_flags.append(panel)
        for rec in recs:
            rec["_run"] = labels[d]
            rec["_dir"] = d
            rec["_shown"] = rec["sample"]
            rec["_status"] = status_lookup.get(results, rec.get("target", ""), rec["sample"])
            records.append(rec)
        for row in skipped:
            skipped_rows.append(row if panel else (row[0], "", *row[1:]))
            sample, target, reason, sample_type, date, site = skipped_rows[-1]
            unusable.append(
                {
                    "_run": labels[d],
                    "_shown": sample,
                    "sample": sample,
                    "target": target,
                    "sample_type": sample_type,
                    "collection_date": date,
                    "site": site,
                    "_status": status_lookup.get(results, target, sample),
                    "_reason": reason,
                }
            )
    return records, unusable, skipped_rows, panel_flags, errors, warnings, notes


def run(args: argparse.Namespace) -> int:
    errors: list[str] = []
    warnings: list[str] = []
    notes: list[str] = []

    dirs = args.from_results or []
    if bool(args.fastq_dir) == bool(dirs):
        print("ERROR: give exactly one of <fastq_dir> or --from-results", file=sys.stderr)
        return 2
    if dirs and not args.mode:
        print("ERROR: --from-results requires --mode demix|aggregate", file=sys.stderr)
        return 2
    if args.mode and not dirs:
        print("ERROR: --mode is only valid with --from-results", file=sys.stderr)
        return 2

    if args.targets and not dirs:
        print("ERROR: --targets is only valid with --from-results", file=sys.stderr)
        return 2
    given = ["--" + k.replace("_", "-") for k in (*FILTER_OPTIONS, "list") if getattr(args, k)]
    if args.on_duplicate != "error":
        given.append("--on-duplicate")
    if given and not dirs:
        print(f"ERROR: {', '.join(given)} only valid with --from-results", file=sys.stderr)
        return 2
    flt, flt_errors = build_filters(args)
    if flt_errors:
        for e in flt_errors:
            print(f"ERROR: {e}", file=sys.stderr)
        return 2

    unusable: list[dict] = []
    status_lookup = StatusLookup()
    panel_flags: list[bool] = []
    if dirs:
        wanted = {t.strip() for t in args.targets.split(",") if t.strip()} if args.targets else None
        records, unusable, skipped_rows, panel_flags, errors, warnings, notes = discover_all(
            args, dirs, wanted, status_lookup
        )
        panel = any(panel_flags)
        if panel and not all(panel_flags):
            warnings.append(
                "panel and single-target results are mixed: the rows of the single-target "
                "results have a blank target"
            )
        if args.skipped_out and not errors and not args.list:
            with open(args.skipped_out, "w", newline="", encoding="utf-8") as fh:
                if panel:
                    fh.write("sample\ttarget\treason\tsample_type\tcollection_date\tsite\n")
                    rows = skipped_rows
                else:
                    fh.write("sample\treason\tsample_type\tcollection_date\tsite\n")
                    rows = [(r[0], *r[2:]) for r in skipped_rows]
                fh.writelines("\t".join(row) + "\n" for row in sorted(rows))
        columns = DEMIX_COLUMNS if args.mode == "demix" else AGGREGATE_COLUMNS
        if panel:
            columns = [columns[0], "target", *columns[1:]]
        path_cols = ("variants", "depths", "fastp_json") if args.mode == "demix" else ("demix",)
    else:
        records, errors, warnings = discover_fastqs(
            Path(args.fastq_dir),
            args.suffix,
            args.negative_pattern,
            args.exclude_pattern,
            args.recursive,
        )
        columns = FULL_COLUMNS
        path_cols = ("fastq_1", "fastq_2")

    if args.metadata:
        meta, m_err, m_warn = load_metadata(Path(args.metadata))
        errors += m_err
        warnings += m_warn
        if not m_err:
            warnings += apply_metadata(records, meta)

    table: list[dict] = []
    if dirs and not errors:
        if "samples" in flt:
            known = {r["sample"] for r in records} | {r["sample"] for r in unusable}
            absent = sorted(flt["samples"] - known)
            if absent:
                warnings.append(f"--samples: no such sample in the results: {', '.join(absent)}")
        if "status" in flt and status_lookup.missing:
            for d, n in sorted(status_lookup.missing.items()):
                notes.append(
                    f"--status ok: no demix status for {n} sample(s) of {d} "
                    f"({STATUS_FILE} missing or without them); kept"
                )
        records, table, sel_errors, sel_notes = select_rows(
            records, unusable, flt, args.on_duplicate
        )
        notes += sel_notes
        if not args.list:
            errors += sel_errors
            if not records and not sel_errors:
                errors.append("no sample is left after the filters")

    for n in notes:
        print(f"NOTE: {n}", file=sys.stderr)
    for w in warnings:
        print(f"WARNING: {w}", file=sys.stderr)
    if errors:
        for e in errors:
            print(f"ERROR: {e}", file=sys.stderr)
        print(f"ERROR: no samplesheet written ({len(errors)} problem(s))", file=sys.stderr)
        return 1

    if args.list:
        labels = run_labels(dirs)
        order = {labels[d]: i for i, d in enumerate(dirs)}
        print_table(table, order, any(panel_flags))
        print("NOTE: --list: no samplesheet written", file=sys.stderr)
        return 0

    rel_to = None
    if args.relative:
        rel_to = Path(args.output).parent if args.output != "-" else Path.cwd()
    for rec in records:
        for col in path_cols:
            rec[col] = make_path(rec.get(col, ""), rel_to)
    for optional in (VIRAL_LOAD, *LIBRARY_COLUMNS):
        if any(rec.get(optional) for rec in records):
            columns = [*columns, optional]
    write_sheet(records, columns, args.output)
    print(
        f"wrote {len(records)} sample(s) to {'stdout' if args.output == '-' else args.output}",
        file=sys.stderr,
    )
    return 0


def main(argv: list[str] | None = None) -> int:
    return run(build_parser().parse_args(argv))


if __name__ == "__main__":
    sys.exit(main())
