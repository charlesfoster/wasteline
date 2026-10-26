#!/usr/bin/env python3
"""Fail-fast wasteline samplesheet validator.

    check_samplesheet.py sheet.csv --mode full|demix|aggregate [--check-files]

stdout: one JSON object {"ok", "mode", "n_samples", "errors", "warnings"}.
stderr: human-readable table of problems.
Exit code: 0 if valid, 1 if any error.

Relative paths in the sheet are resolved against the sheet's own directory (nf-schema semantics).
Without --check-files the filesystem is not touched.
"""

from __future__ import annotations

import argparse
import csv
import datetime
import gzip
import json
import math
import re
import sys
from pathlib import Path

ID_RE = re.compile(r"^[A-Za-z0-9._-]+$")
DATE_RE = re.compile(r"^\d{4}-\d{2}-\d{2}$")
SAMPLE_TYPES = ("sample", "negative")
FASTQ_EXTS = (".fastq.gz", ".fq.gz")
PROTOCOLS = ("amplicon", "capture")
OPTIONAL_META = ["sample_type", "collection_date", "site", "viral_load"]
# demix / aggregate sheets of panel results name the target of each row: one row per (sample, target)
TARGET_COLUMN = "target"
# per-sample library settings (PLAN D20); empty = the run's --protocol / --primer_scheme / --primer_bed
LIBRARY_COLUMNS = ["protocol", "primer_scheme"]
PIPELINE_DIR = Path(__file__).resolve().parent.parent
PRIMER_SCHEME_DIR = PIPELINE_DIR / "assets" / "primer_schemes"

REQUIRED = {
    "full": ["sample", "fastq_1", "fastq_2"],
    "demix": ["sample", "variants", "depths"],
    "aggregate": ["sample", "demix"],
}
OPTIONAL = {
    "full": [*OPTIONAL_META, *LIBRARY_COLUMNS],
    "demix": ["fastp_json", *OPTIONAL_META, *LIBRARY_COLUMNS, TARGET_COLUMN],
    "aggregate": [*OPTIONAL_META, TARGET_COLUMN],
}
# Columns holding file paths, per mode (required ones must be non-empty; optional may be empty).
PATH_COLUMNS = {
    "full": ["fastq_1", "fastq_2"],
    "demix": ["variants", "depths", "fastp_json"],
    "aggregate": ["demix"],
}


def valid_date(value: str) -> bool:
    if not DATE_RE.match(value):
        return False
    try:
        datetime.date.fromisoformat(value)
    except ValueError:
        return False
    return True


def valid_viral_load(value: str) -> bool:
    """A positive, finite number (blank is checked by the caller)."""
    try:
        number = float(value)
    except ValueError:
        return False
    return math.isfinite(number) and number > 0


def is_bed_path(value: str) -> bool:
    """A primer_scheme cell is a BED file path when it has a '/' or ends in .bed, else a bundled scheme name."""
    return "/" in value or value.endswith(".bed")


def pathogen_profile(pathogen: str) -> str:
    """The conf/pathogens profile of a --pathogen name, read from the map in nextflow.config ('generic' if unlisted)."""
    try:
        text = (PIPELINE_DIR / "nextflow.config").read_text()
        line = re.search(r'includeConfig "conf/pathogens/.*', text).group(0)
        return dict(re.findall(r"'([^']+)': '(\w+)'", line)).get(pathogen, "generic")
    except (OSError, AttributeError):
        return "generic"


def bundled_schemes(pathogen: str | None = None, scheme_dirs: tuple = ()) -> list[str]:
    """Names usable in a primer_scheme cell: the SARS-CoV-2 schemes (no pathogen given, or SARS-CoV-2), or the schemes
    bundled for the pathogen (assets/primer_schemes/<profile>/), plus the <name>.bed files of the user's directories."""
    profile = pathogen_profile(pathogen) if pathogen else "sars2"
    bundled = PRIMER_SCHEME_DIR if profile == "sars2" else PRIMER_SCHEME_DIR / profile
    dirs = [bundled, *[Path(d) for d in scheme_dirs]]
    return sorted({p.stem for d in dirs for p in d.glob("*.bed")})


def resolve(value: str, base: Path) -> Path:
    p = Path(value)
    return p if p.is_absolute() else base / p


def check_fastq_header(path: Path) -> str | None:
    """Return a problem description, or None if the file looks like a FASTQ."""
    try:
        opener = gzip.open if path.name.endswith(".gz") else open
        with opener(path, "rt", encoding="ascii", errors="replace") as fh:
            first = fh.readline()
    except (OSError, EOFError) as exc:
        return f"cannot read ({exc})"
    if not first:
        return "file is empty"
    if not first.startswith("@"):
        return "first line does not start with '@' (not FASTQ?)"
    return None


def check_file(path: Path, column: str) -> str | None:
    if not path.exists():
        return "file does not exist"
    if not path.is_file():
        return "not a regular file"
    if column in ("fastq_1", "fastq_2"):
        return check_fastq_header(path)
    try:
        if path.stat().st_size == 0:
            return "file is empty"
    except OSError as exc:
        return f"cannot stat ({exc})"
    return None


def validate(
    sheet: Path,
    mode: str,
    check_files: bool = False,
    pathogen: str | None = None,
    scheme_dirs: tuple = (),
) -> dict:
    errors: list[str] = []
    warnings: list[str] = []
    result = {"ok": False, "mode": mode, "n_samples": 0, "errors": errors, "warnings": warnings}

    try:
        with open(sheet, newline="", encoding="utf-8-sig") as fh:
            reader = csv.reader(fh)
            rows = list(reader)
    except (OSError, UnicodeDecodeError, csv.Error) as exc:
        errors.append(f"cannot read samplesheet {sheet}: {exc}")
        return result

    rows = [r for r in rows if any(c.strip() for c in r)]  # drop blank lines
    if not rows:
        errors.append("samplesheet is empty")
        return result

    header = rows[0]
    if any(h != h.strip() for h in header):
        warnings.append("header names have leading/trailing whitespace")
        header = [h.strip() for h in header]
    dups = sorted({h for h in header if header.count(h) > 1})
    if dups:
        errors.append(f"duplicate column(s) in header: {', '.join(dups)}")
    missing = [c for c in REQUIRED[mode] if c not in header]
    if missing:
        errors.append(
            f"missing required column(s) for mode '{mode}': {', '.join(missing)} "
            f"(header: {', '.join(header)})"
        )
    known = set(REQUIRED[mode]) | set(OPTIONAL[mode])
    for col in header:
        if col not in known:
            warnings.append(f"unknown column '{col}' (ignored by mode '{mode}')")
    if errors:
        return result

    data = rows[1:]
    if not data:
        errors.append("samplesheet has a header but no samples")
        return result

    base = sheet.resolve().parent
    seen: dict[str, int] = {}
    seen_pairs: dict[tuple, int] = {}
    path_cols = [c for c in PATH_COLUMNS[mode] if c in header]
    for lineno, raw in enumerate(data, start=2):
        if len(raw) > len(header):
            errors.append(f"line {lineno}: {len(raw)} fields but header has {len(header)}")
            continue
        # no zip(strict=): this script runs on the host Python (>= 3.8)
        row = dict(zip(header, raw + [""] * (len(header) - len(raw))))
        sample = row["sample"]
        tag = f"line {lineno} ({sample or 'no sample'})"

        if not sample:
            errors.append(f"line {lineno}: empty sample ID")
            continue
        if not ID_RE.match(sample):
            errors.append(f"{tag}: sample ID must match ^[A-Za-z0-9._-]+$")
        target = row.get(TARGET_COLUMN, "").strip()
        if TARGET_COLUMN in header:
            # with a target column the same sample may appear once per target
            if (sample, target) in seen_pairs:
                what = f"sample ID and target '{target}'" if target else "sample ID"
                errors.append(
                    f"{tag}: duplicate {what} (first seen on line {seen_pairs[(sample, target)]})"
                )
            else:
                seen_pairs[(sample, target)] = lineno
            seen.setdefault(sample, lineno)
        elif sample in seen:
            errors.append(f"{tag}: duplicate sample ID (first seen on line {seen[sample]})")
        else:
            seen[sample] = lineno

        stype = row.get("sample_type", "")
        if stype and stype not in SAMPLE_TYPES:
            errors.append(f"{tag}: sample_type '{stype}' must be one of {'|'.join(SAMPLE_TYPES)}")
        date = row.get("collection_date", "")
        if date and not valid_date(date):
            errors.append(f"{tag}: collection_date '{date}' is not ISO YYYY-MM-DD")
        load = row.get("viral_load", "").strip()
        if load and not valid_viral_load(load):
            errors.append(f"{tag}: viral_load '{load}' must be a positive number (or blank)")

        protocol = row.get("protocol", "").strip()
        if protocol and protocol not in PROTOCOLS:
            errors.append(
                f"{tag}: protocol '{protocol}' must be one of {'|'.join(PROTOCOLS)} (or blank)"
            )
        scheme = row.get("primer_scheme", "").strip()
        if scheme and protocol == "capture":
            warnings.append(f"{tag}: primer_scheme '{scheme}' is ignored for protocol capture")
        elif scheme and is_bed_path(scheme):
            if check_files and not resolve(scheme, base).is_file():
                errors.append(f"{tag}: primer_scheme BED file '{scheme}' does not exist")
        elif scheme and scheme not in bundled_schemes(pathogen, scheme_dirs):
            available = ", ".join(bundled_schemes(pathogen, scheme_dirs)) or "none"
            errors.append(
                f"{tag}: primer_scheme '{scheme}' is not a primer scheme available for "
                f"{pathogen or 'the default pathogen'} (bundled or in --primer-scheme-dir: {available}) "
                f"and not a path to a BED file (.bed)"
            )

        for col in path_cols:
            val = row.get(col, "")
            required = col in REQUIRED[mode]
            if not val:
                if required:
                    what = (
                        "empty fastq_2 (paired-end data is required)"
                        if col == "fastq_2"
                        else f"empty {col}"
                    )
                    errors.append(f"{tag}: {what}")
                continue
            if re.search(r"\s", val):
                errors.append(f"{tag}: {col} path contains whitespace: '{val}'")
                continue
            if col in ("fastq_1", "fastq_2") and not val.endswith(FASTQ_EXTS):
                errors.append(f"{tag}: {col} must end with .fastq.gz or .fq.gz: '{val}'")
                continue
            if check_files:
                problem = check_file(resolve(val, base), col)
                if problem:
                    errors.append(f"{tag}: {col} '{val}': {problem}")

        if mode == "full" and row.get("fastq_1") and row.get("fastq_1") == row.get("fastq_2"):
            errors.append(f"{tag}: fastq_1 and fastq_2 are the same file")

    result["n_samples"] = len(seen)
    result["ok"] = not errors
    return result


def print_table(result: dict, sheet: str) -> None:
    status = "PASS" if result["ok"] else "FAIL"
    print(
        f"{status}  {sheet}  mode={result['mode']}  samples={result['n_samples']}  "
        f"errors={len(result['errors'])}  warnings={len(result['warnings'])}",
        file=sys.stderr,
    )
    rows = [("ERROR", m) for m in result["errors"]] + [("WARN", m) for m in result["warnings"]]
    if rows:
        print(f"{'LEVEL':<6} MESSAGE", file=sys.stderr)
        for level, msg in rows:
            print(f"{level:<6} {msg}", file=sys.stderr)


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description="Validate a wasteline samplesheet.")
    p.add_argument("samplesheet")
    p.add_argument("--mode", required=True, choices=sorted(REQUIRED))
    p.add_argument(
        "--check-files",
        action="store_true",
        help="verify referenced files exist (FASTQs: readable, first line starts '@')",
    )
    p.add_argument(
        "--pathogen",
        help="the run's --pathogen: primer_scheme names are checked against the schemes bundled for it "
        "(default: the SARS-CoV-2 schemes)",
    )
    p.add_argument(
        "--primer-scheme-dir",
        action="append",
        default=[],
        help="directory of your own <name>.bed schemes (the run's --primer_scheme_dir); repeatable",
    )
    args = p.parse_args(argv)
    result = validate(
        Path(args.samplesheet),
        args.mode,
        args.check_files,
        args.pathogen,
        tuple(args.primer_scheme_dir),
    )
    print(json.dumps(result, indent=2))
    print_table(result, args.samplesheet)
    return 0 if result["ok"] else 1


if __name__ == "__main__":
    sys.exit(main())
