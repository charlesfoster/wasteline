#!/usr/bin/env python3
"""Assemble `pipeline_info/run_info.json` (run provenance).

The Nextflow side supplies a base JSON (pipeline version, git commit, Nextflow version, mode,
container tags, parameters, start time). This script adds sha256 checksums of the key inputs, the
Freyja barcode dataset description, the tool versions collected during the run (and the Freyja
version in particular) and the time the record was generated.
"""

import argparse
import hashlib
import json
import re
import sys
from datetime import UTC, datetime
from pathlib import Path

from freyja_dataset_info import FREYJA_DATA_KEYS


def sha256_file(path, chunk=1 << 20):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(chunk), b""):
            digest.update(block)
    return digest.hexdigest()


def checksum_entry(path):
    path = Path(path)
    return {"file": path.name, "sha256": sha256_file(path), "size_bytes": path.stat().st_size}


def parse_versions_yaml(text):
    """Parse the flat `PROCESS:\\n  tool: version` YAML written by the pipeline (no PyYAML needed)."""
    versions = {}
    section = None
    for raw in text.splitlines():
        if not raw.strip() or raw.lstrip().startswith("#"):
            continue
        if not raw.startswith((" ", "\t")):
            key = raw.strip().rstrip(":").strip().strip("\"'")
            section = key
            versions.setdefault(section, {})
            continue
        if section is None:
            continue
        match = re.match(r"^\s+([^:]+):\s*(.*)$", raw)
        if match:
            tool = match.group(1).strip().strip("\"'")
            versions[section][tool] = match.group(2).strip().strip("\"'")
    return versions


def tool_versions(versions_by_process):
    """Collapse {process: {tool: version}} to {tool: version} (first one wins)."""
    tools = {}
    for entries in versions_by_process.values():
        for tool, version in entries.items():
            tools.setdefault(tool, version)
    return dict(sorted(tools.items()))


def check_freyja_data(freyja_data):
    """The freyja_data block is mandatory and must carry every key of PLAN section 5.13."""
    missing = [key for key in FREYJA_DATA_KEYS if key not in freyja_data]
    if missing:
        raise ValueError(f"freyja_data block is missing keys: {', '.join(missing)}")
    if not freyja_data["barcode_version"]:
        raise ValueError("freyja_data.barcode_version is empty")


def build_run_info(
    base, samplesheet, bed, fasta, dataset_info, versions_text, now=None, dataset_info_compare=None
):
    now = now or datetime.now(UTC)
    by_process = parse_versions_yaml(versions_text) if versions_text else {}
    tools = tool_versions(by_process)
    info = dict(base)
    info["tool_versions"] = tools
    info["freyja_version"] = tools.get("freyja")
    info["freyja_data"] = dataset_info
    if dataset_info_compare is not None:
        # --compare_barcodes: the second dataset is recorded with the same full block
        check_freyja_data(dataset_info_compare)
        info["freyja_data_compare"] = dataset_info_compare
    # Only the inputs the mode actually uses are checksummed (demix/aggregate need no BED/FASTA)
    info["input_checksums"] = {"samplesheet": checksum_entry(samplesheet)}
    # one primer BED: primer_bed; samples with different primer schemes (PLAN D20): primer_beds, one entry each
    beds = [bed] if isinstance(bed, (str, Path)) else list(bed or [])
    if len(beds) == 1:
        info["input_checksums"]["primer_bed"] = checksum_entry(beds[0])
    elif beds:
        info["input_checksums"]["primer_beds"] = [checksum_entry(b) for b in beds]
    if fasta:
        info["input_checksums"]["reference_fasta"] = checksum_entry(fasta)
    info.setdefault("previous_runs", [])
    info["provenance_generated"] = now.strftime("%Y-%m-%dT%H:%M:%SZ")
    return info


def parse_args(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--base", required=True, help="JSON with fields known to Nextflow")
    parser.add_argument("--samplesheet", required=True)
    parser.add_argument(
        "--bed", nargs="+", default=None, help="primer BED(s) of the samples (modes that map reads)"
    )
    parser.add_argument("--fasta", default=None, help="reference FASTA (modes that map reads)")
    parser.add_argument("--dataset-info", required=True, help="freyja_dataset_info.py JSON")
    parser.add_argument(
        "--dataset-info-compare",
        default=None,
        help="freyja_dataset_info.py JSON of --compare_barcodes",
    )
    parser.add_argument("--versions", default=None, help="software_versions.yml")
    parser.add_argument("--output", required=True)
    return parser.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    try:
        base = json.loads(Path(args.base).read_text())
        dataset_info = json.loads(Path(args.dataset_info).read_text())
        versions_text = Path(args.versions).read_text() if args.versions else ""
        compare = (
            json.loads(Path(args.dataset_info_compare).read_text())
            if args.dataset_info_compare
            else None
        )
        info = build_run_info(
            base,
            args.samplesheet,
            args.bed,
            args.fasta,
            dataset_info,
            versions_text,
            dataset_info_compare=compare,
        )
    except (OSError, ValueError) as err:
        print(f"make_run_info.py: {err}", file=sys.stderr)
        return 1
    Path(args.output).write_text(json.dumps(info, indent=2) + "\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
