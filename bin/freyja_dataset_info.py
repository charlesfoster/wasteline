#!/usr/bin/env python3
"""Describe a Freyja dataset (directory or single barcode file) for provenance.

A dataset directory is what `freyja update --outdir` writes: `usher_barcodes.feather` (or `.csv`),
`last_barcode_update.txt`, `lineage_mutations.json`, `curated_lineages.json`, `lineages.yml`,
`pathogen_config.yml`. A single `.csv` / `.feather` barcode file is also accepted.

For the pathogens other than the default one, `freyja update --pathogen X` writes `X_barcodes.csv` (`.feather`
for MTB), `X_lineages.yml` (only for the few pathogens that have a hierarchy), `pathogen_config.yml` and
`last_barcode_update_other.txt` (one `X:YYYY-MM-DD` line per pathogen) instead; no `lineage_mutations.json`. These
names are recognised too. Their barcode version is the date (verbatim), its ISO form is `YYYY-MM-DDT00-00`.

The barcode version is the verbatim content of `last_barcode_update.txt` (`MM_DD_YYYY-HH-MM`, e.g.
`09_17_2026-00-46`); its ISO form `YYYY-MM-DDTHH-MM` names directories in the local database. A
single barcode file has version `custom`.

Outputs: the full `freyja_data` block (`--output`) and a small `barcode_version.json`
(`--version-output`) that is attached to every per-sample status file and published in `summary/`.
"""

import argparse
import hashlib
import json
import re
import sys
from pathlib import Path

BARCODE_NAMES = ("usher_barcodes.feather", "usher_barcodes.csv")
DATE_FILE = "last_barcode_update.txt"
DATE_FILE_OTHER = "last_barcode_update_other.txt"  # `PATHOGEN:YYYY-MM-DD` lines, written for the non-default pathogens
SIDECAR = "wasteline_dataset.json"
BARCODE_SUFFIXES = (".csv", ".feather")
CUSTOM = "custom"
SOURCES = ("latest", "container", "pinned", "custom")

# Keys that every `freyja_data` block must carry (PLAN section 5.13).
FREYJA_DATA_KEYS = (
    "source",
    "barcode_version",
    "barcode_version_iso",
    "barcodes_file",
    "barcodes_sha256",
    "barcodes_format",
    "lineage_mutations_sha256",
    "curated_lineages_sha256",
    "lineages_yml_sha256",
    "pathogen",
    "freyja_version",
    "freyja_reported_barcode_version",
    "dataset_path",
    "downloaded_at",
)
# Keys of the small record copied into status files and summary/barcode_version.json.
BARCODE_VERSION_KEYS = (
    "barcode_version",
    "barcode_version_iso",
    "barcodes_sha256",
    "source",
    "freyja_version",
)

_STAMP = re.compile(r"^(\d{2})_(\d{2})_(\d{4})-(\d{2})-(\d{2})$")
_ISO = re.compile(r"^(\d{4})-(\d{2})-(\d{2})T(\d{2})[-:](\d{2})$")
_DATE = re.compile(r"^\d{4}-\d{2}-\d{2}$")


def _valid(month, day, hour, minute):
    return 1 <= int(month) <= 12 and 1 <= int(day) <= 31 and int(hour) < 24 and int(minute) < 60


def stamp_to_iso(text):
    """`MM_DD_YYYY-HH-MM` or `YYYY-MM-DDTHH-MM` (`:` allowed) -> `YYYY-MM-DDTHH-MM`, else None."""
    text = (text or "").strip()
    match = _STAMP.match(text)
    if match:
        month, day, year, hour, minute = match.groups()
    else:
        match = _ISO.match(text)
        if not match:
            return None
        year, month, day, hour, minute = match.groups()
    if not _valid(month, day, hour, minute):
        return None
    return f"{year}-{month}-{day}T{hour}-{minute}"


def iso_to_stamp(iso):
    """Inverse of stamp_to_iso for a valid ISO form."""
    year, month, day, hour, minute = _ISO.match(iso).groups()
    return f"{month}_{day}_{year}-{hour}-{minute}"


def is_bare_date(text):
    return bool(_DATE.match((text or "").strip()))


def sha256_file(path, chunk=1 << 20):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(chunk), b""):
            digest.update(block)
    return digest.hexdigest()


def _sha_if_present(path):
    return sha256_file(path) if path.is_file() else None


def find_barcode_file(directory):
    for name in BARCODE_NAMES:
        candidate = directory / name
        if candidate.is_file():
            return candidate
    # other pathogens: <PATHOGEN>_barcodes.csv / .feather (exactly one)
    others = sorted(
        p for p in directory.glob("*_barcodes.*") if p.is_file() and p.suffix in BARCODE_SUFFIXES
    )
    if len(others) == 1:
        return others[0]
    if others:
        raise ValueError(
            f"{directory}: several barcode files ({', '.join(p.name for p in others)}); "
            "use a dataset directory made for one pathogen"
        )
    raise ValueError(
        f"{directory}: no {' or '.join(BARCODE_NAMES)} or <PATHOGEN>_barcodes.csv found; "
        "not a Freyja dataset directory"
    )


def find_lineages_yml(directory):
    """`lineages.yml`, or the `<PATHOGEN>_lineages.yml` of another pathogen; None when the dataset has none."""
    candidate = directory / "lineages.yml"
    if candidate.is_file():
        return candidate
    others = sorted(directory.glob("*_lineages.yml"))
    if len(others) == 1 and _is_yaml_list(others[0]):
        return others[0]
    return None


def _is_yaml_list(path):
    """Freyja's pathogen_config.yml points MPX at a GitHub `blob` page, so `freyja update` saves HTML as the
    hierarchy: a real lineages.yml is a YAML list, i.e. its first line starts with `-`."""
    try:
        with open(path, encoding="utf-8", errors="replace") as handle:
            for line in handle:
                if line.strip():
                    return line.lstrip().startswith("-")
    except OSError:
        pass
    return False


def other_pathogen_date(directory, pathogen=""):
    """Date (YYYY-MM-DD) of the `PATHOGEN:DATE` line of last_barcode_update_other.txt, or None.

    The line of `pathogen` when there is one, otherwise the only line the file has.
    """
    path = directory / DATE_FILE_OTHER
    if not path.is_file():
        return None
    entries = {}
    for line in path.read_text().splitlines():
        name, _, date = line.strip().partition(":")
        if name and date:
            entries[name.strip()] = date.strip()
    if pathogen in entries:
        return entries[pathogen]
    return next(iter(entries.values())) if len(entries) == 1 else None


def dataset_paths(path):
    """(barcode file, lineages.yml or None) of a dataset directory or single barcode file."""
    path = Path(path)
    if path.is_dir():
        return find_barcode_file(path), find_lineages_yml(path)
    if path.is_file() and path.suffix in BARCODE_SUFFIXES:
        return path, None
    raise ValueError(f"{path}: not a Freyja dataset directory or .csv/.feather barcode file")


def parse_reported_version(text):
    """Barcode version from `freyja demix --version` output (None if it printed a path only)."""
    match = re.search(r"Barcode version:\s*\n\s*(\S+)", text or "")
    return match.group(1) if match else None


def describe(
    path,
    source,
    pathogen,
    freyja_version=None,
    dataset_path=None,
    reported_text=None,
):
    """Return the `freyja_data` block. Raises ValueError for unusable datasets or sources."""
    if source not in SOURCES:
        raise ValueError(f"source must be one of {', '.join(SOURCES)}, got {source!r}")
    path = Path(path)
    downloaded_at = None
    if path.is_dir():
        barcode = find_barcode_file(path)
        date_file = path / DATE_FILE
        version = date_file.read_text().strip() if date_file.is_file() else ""
        version = version.splitlines()[0].strip() if version else CUSTOM
        iso = stamp_to_iso(version)
        other_date = None if date_file.is_file() else other_pathogen_date(path, pathogen)
        if other_date and is_bare_date(other_date):
            # another pathogen: the version is the date of its barcodes; the ISO form (db directory name) has no time
            version, iso = other_date, f"{other_date}T00-00"
        sidecar = path / SIDECAR
        if sidecar.is_file():
            downloaded_at = json.loads(sidecar.read_text()).get("downloaded_at")
        lineage_mutations = _sha_if_present(path / "lineage_mutations.json")
        curated = _sha_if_present(path / "curated_lineages.json")
        lineages_yml = _sha_if_present(find_lineages_yml(path) or path / "lineages.yml")
    elif path.is_file():
        if path.suffix not in BARCODE_SUFFIXES:
            raise ValueError(f"{path}: barcode file must end in {' or '.join(BARCODE_SUFFIXES)}")
        barcode = path
        version, iso = CUSTOM, None
        lineage_mutations = curated = lineages_yml = None
    else:
        raise ValueError(f"{path}: does not exist")
    if version != CUSTOM and iso is None:
        raise ValueError(f"{path}: unrecognised {DATE_FILE} content {version!r}")
    # `freyja demix --version` always reports the data bundled in the Freyja package, never an
    # --outdir dataset, so it is only meaningful as "the version used" for source == container.
    reported = parse_reported_version(reported_text)
    return {
        "source": source,
        "barcode_version": version,
        "barcode_version_iso": iso,
        "barcodes_file": barcode.name,
        "barcodes_sha256": sha256_file(barcode),
        "barcodes_format": barcode.suffix.lstrip("."),
        "lineage_mutations_sha256": lineage_mutations,
        "curated_lineages_sha256": curated,
        "lineages_yml_sha256": lineages_yml,
        "pathogen": pathogen,
        "freyja_version": freyja_version,
        "freyja_reported_barcode_version": reported if source == "container" else None,
        "freyja_bundled_barcode_version": reported,
        "dataset_path": dataset_path if dataset_path is not None else str(path),
        "downloaded_at": downloaded_at,
    }


def barcode_version_record(freyja_data):
    """The small record carried by status files and `summary/barcode_version.json`."""
    return {key: freyja_data.get(key) for key in BARCODE_VERSION_KEYS}


def _read(path):
    return Path(path).read_text() if path else None


def parse_args(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("dataset", help="dataset directory or .csv/.feather barcode file")
    parser.add_argument("--source", choices=SOURCES, default="custom")
    parser.add_argument("--pathogen", default="")
    parser.add_argument("--freyja-version", default=None)
    parser.add_argument("--dataset-path", default=None, help="original location (recorded)")
    parser.add_argument("--reported-file", default=None, help="stdout of `freyja demix --version`")
    parser.add_argument("--output", default=None, help="write the freyja_data JSON here")
    parser.add_argument("--version-output", default=None, help="write barcode_version.json here")
    parser.add_argument(
        "--print-barcodes-path",
        action="store_true",
        help="print only the path of the barcode file to pass to `freyja demix --barcodes`",
    )
    parser.add_argument(
        "--print-lineageyml-path",
        action="store_true",
        help="print only the path of lineages.yml (empty line if the dataset has none)",
    )
    parser.add_argument(
        "--print-iso",
        action="store_true",
        help="print only the ISO barcode version (YYYY-MM-DDTHH-MM); fail if undeterminable",
    )
    return parser.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    if args.print_barcodes_path or args.print_lineageyml_path:
        try:
            barcode, lineage_yml = dataset_paths(args.dataset)
        except (OSError, ValueError) as err:
            print(f"freyja_dataset_info.py: {err}", file=sys.stderr)
            return 1
        print(barcode if args.print_barcodes_path else (lineage_yml or ""))
        return 0
    try:
        info = describe(
            args.dataset,
            args.source,
            args.pathogen,
            freyja_version=args.freyja_version,
            dataset_path=args.dataset_path,
            reported_text=_read(args.reported_file),
        )
    except (OSError, ValueError) as err:
        print(f"freyja_dataset_info.py: {err}", file=sys.stderr)
        return 1
    if args.print_iso:
        if not info["barcode_version_iso"]:
            print(
                f"freyja_dataset_info.py: {args.dataset} has no usable {DATE_FILE} "
                f"(version {info['barcode_version']!r})",
                file=sys.stderr,
            )
            return 1
        print(info["barcode_version_iso"])
        return 0
    text = json.dumps(info, indent=2) + "\n"
    if args.output:
        Path(args.output).write_text(text)
    else:
        sys.stdout.write(text)
    if args.version_output:
        Path(args.version_output).write_text(
            json.dumps(barcode_version_record(info), indent=2) + "\n"
        )
    return 0


if __name__ == "__main__":
    sys.exit(main())
