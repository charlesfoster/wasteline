#!/usr/bin/env python3
"""Build wasteline samplesheets.

Two uses:

1. FASTQ directory -> `full` mode samplesheet
   make_samplesheet.py <fastq_dir> [--suffix S] [--metadata meta.csv] -o samplesheet.csv

2. Existing results directory -> `demix` or `aggregate` mode samplesheet
   make_samplesheet.py --from-results <dir> --mode demix|aggregate -o sheet.csv

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


def locate_new(results: Path, sample: str) -> dict[str, Path | None]:
    base = results / "samples" / sample
    fj = base / "freyja"
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


def load_results_sheet(results: Path) -> tuple[dict[str, dict[str, str]], str | None]:
    for rel in ("samplesheet.csv", "summary/samplesheet.csv"):
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
                    if r["sample"]
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


def discover_results(
    results: Path,
    mode: str,
    negative_pattern: str = DEFAULT_NEGATIVE_PATTERN,
    skipped: list[tuple[str, ...]] | None = None,
) -> tuple[list[dict[str, str]], list[str], list[str], list[str]]:
    """Return (records, errors, warnings, notes). `notes` are informational lines.

    Samples that are found but unusable are appended to `skipped` as (sample, reason, sample_type,
    collection_date, site) when a list is given (the type is inferred like for usable samples): missing required files, or (aggregate mode) a zero-byte .demix file, which the old
    pipeline left behind when `freyja demix` failed (`|| touch`). Such a file is never aggregated.
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
        help="build a demix/aggregate sheet from an existing results directory",
    )
    p.add_argument("--mode", choices=["demix", "aggregate"], help="sheet type for --from-results")
    p.add_argument(
        "--relative",
        action="store_true",
        help="write paths relative to the output file's directory (default: absolute)",
    )
    p.add_argument("-o", "--output", default="-", help="output CSV (default: stdout)")
    p.add_argument(
        "--skipped-out",
        metavar="TSV",
        help="with --from-results: write (sample, reason) of every sample that was found but "
        "skipped (missing files, zero-byte legacy .demix) to this TSV",
    )
    return p


def run(args: argparse.Namespace) -> int:
    errors: list[str] = []
    warnings: list[str] = []
    notes: list[str] = []

    if bool(args.fastq_dir) == bool(args.from_results):
        print("ERROR: give exactly one of <fastq_dir> or --from-results", file=sys.stderr)
        return 2
    if args.from_results and not args.mode:
        print("ERROR: --from-results requires --mode demix|aggregate", file=sys.stderr)
        return 2
    if args.mode and not args.from_results:
        print("ERROR: --mode is only valid with --from-results", file=sys.stderr)
        return 2

    if args.from_results:
        skipped: list[tuple[str, str]] = []
        records, errors, warnings, notes = discover_results(
            Path(args.from_results), args.mode, args.negative_pattern, skipped
        )
        if args.skipped_out and not errors:
            with open(args.skipped_out, "w", newline="", encoding="utf-8") as fh:
                fh.write("sample\treason\tsample_type\tcollection_date\tsite\n")
                fh.writelines("\t".join(row) + "\n" for row in sorted(skipped))
        columns = DEMIX_COLUMNS if args.mode == "demix" else AGGREGATE_COLUMNS
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

    for n in notes:
        print(f"NOTE: {n}", file=sys.stderr)
    for w in warnings:
        print(f"WARNING: {w}", file=sys.stderr)
    if errors:
        for e in errors:
            print(f"ERROR: {e}", file=sys.stderr)
        print(f"ERROR: no samplesheet written ({len(errors)} problem(s))", file=sys.stderr)
        return 1

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
