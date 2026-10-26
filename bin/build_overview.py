#!/usr/bin/env python3
"""Build the OVERVIEW page of a panel run (`summary/overview.html`): a sample x target matrix.

Inputs: `summary/targets.tsv` (the typing of every sample against every target) and, per target, the
tables the run wrote for it: `freyja_aggregated.tsv` (top lineage), `qc_summary.tsv` (QC status),
`demix_status.tsv` (why a sample has no lineage) and `run_info.json` (header). The per-target files
are listed in a manifest (`--inputs`, tab separated: target, role, path; roles `aggregated`, `qc`,
`demix_status`, `run_info`, `report`); any may be missing for a target and the cell then shows `-`.
`target` lists the targets to show (without it, every target of targets.tsv); `report` marks a target whose
`<target>/report.html` exists (the cells link to it).

The page is one self-contained offline HTML file (same look and light/dark tokens as the report).
Needs jinja2; everything else is the standard library. Docs: docs/usage.md (Panels).
"""

import argparse
import csv
import json
import sys
from datetime import date
from pathlib import Path
from urllib.parse import quote

import jinja2

DEFAULT_TEMPLATE = (
    Path(__file__).resolve().parent.parent / "assets" / "templates" / "overview.html.j2"
)
STATUSES = ("detected", "trace", "not_detected")
ROLES = ("target", "aggregated", "qc", "demix_status", "run_info", "report")
SAMPLE_FLAGS = ("qc_fail", "no_target_detected")
PSEUDONYM_PREFIX = "S"


def read_tsv(path):
    """Rows of a tab separated file as lists of dicts keyed by header; [] when missing or empty."""
    if not path or not Path(path).is_file():
        return []
    with open(path, newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def read_manifest(path):
    """{target: {role: path}} from the manifest; unknown roles and blank lines are ignored."""
    out = {}
    if not path or not Path(path).is_file():
        return out
    for line in Path(path).read_text(encoding="utf-8").splitlines():
        cells = line.split("\t")
        if len(cells) >= 3 and cells[1] in ROLES:
            out.setdefault(cells[0], {})[cells[1]] = cells[2]
    return out


def to_float(value):
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def top_lineage(row):
    """(lineage, abundance 0-1, number of lineages) from a freyja_aggregated row, or None."""
    lineages = (row.get("lineages") or "").split()
    abundances = [to_float(x) for x in (row.get("abundances") or "").split()]
    pairs = [(a, name) for name, a in zip(lineages, abundances, strict=False) if a is not None]
    if not pairs:
        return None
    abundance, name = max(pairs, key=lambda p: p[0])
    return name, abundance, len(pairs)


def load_aggregated(path):
    """{sample: (lineage, abundance, n)} ; the sample is the first column (its header is empty)."""
    out = {}
    if not path or not Path(path).is_file():
        return out
    with open(path, newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        first = reader.fieldnames[0] if reader.fieldnames else None
        for row in reader:
            top = top_lineage(row)
            if first is not None and row.get(first) and top:
                out[row[first]] = top
    return out


def load_by_sample(path, key="sample"):
    return {r[key]: r for r in read_tsv(path) if r.get(key)}


def pseudonyms(samples):
    names = sorted(samples)
    width = max(3, len(str(len(names))))
    return {n: f"{PSEUDONYM_PREFIX}{i:0{width}d}" for i, n in enumerate(names, start=1)}


def scrub(text, mapping):
    """Replace known sample IDs (longest first) inside a free-text reason."""
    for name in sorted(mapping, key=lambda n: (-len(n), n)):
        text = text.replace(name, mapping[name])
    return text


def truthy(value):
    return str(value).strip().lower() == "true"


def basename(text):
    return str(text).replace("\\", "/").rstrip("/").split("/")[-1] if text else ""


def build_context(targets_rows, manifest, args):
    """Everything the template needs: header, target strip, matrix rows."""
    target_names = []
    sample_names = []
    for row in targets_rows:
        if row.get("target") and row["target"] not in target_names:
            target_names.append(row["target"])
        if row.get("sample") and row["sample"] not in sample_names:
            sample_names.append(row["sample"])
    only = {t for t, files in manifest.items() if "target" in files}
    if only:
        target_names = [t for t in target_names if t in only]
    targets_rows = [r for r in targets_rows if r.get("target") in target_names]

    id_map = pseudonyms(sample_names) if args.anonymize == "full" else {}

    def label(sample):
        return id_map.get(sample, sample)

    per_target = {}
    for name in target_names:
        files = manifest.get(name, {})
        per_target[name] = {
            "aggregated": load_aggregated(files.get("aggregated")),
            "qc": load_by_sample(files.get("qc")),
            "demix": load_by_sample(files.get("demix_status")),
            "has_report": "report" in files,
        }

    cells = {}
    sample_info = {}
    for row in targets_rows:
        sample, name = row["sample"], row["target"]
        data = per_target[name]
        analysed = truthy(row.get("analysed"))
        status = row.get("status") if row.get("status") in STATUSES else "not_detected"
        cell = {
            "status": status,
            "analysed": analysed,
            "fragments": row.get("fragments") or "0",
            "breadth": to_float(row.get("breadth_pct")),
            "mixed": truthy(row.get("mixed")),
            "secondary": [x for x in (row.get("secondary_targets") or "").split(";") if x],
            "reason": scrub(row.get("reason") or "", id_map),
            "lineage": None,
            "qc": None,
            "qc_reasons": "",
            "demix": None,
            "demix_reason": "",
            "link": f"{quote(name)}/report.html" if analysed and data["has_report"] else None,
        }
        if analysed:
            top = data["aggregated"].get(sample)
            if top:
                cell["lineage"] = {
                    "name": top[0],
                    "pct": round(top[1] * 100, 1),
                    "more": top[2] - 1,
                }
            qc = data["qc"].get(sample)
            if qc:
                cell["qc"] = (qc.get("qc_status") or "").upper() or None
                cell["qc_reasons"] = qc.get("qc_reasons") or ""
            dm = data["demix"].get(sample)
            if dm:
                cell["demix"] = dm.get("status") or None
                cell["demix_reason"] = dm.get("reason") or ""
        cells[(sample, name)] = cell
        info = sample_info.setdefault(
            sample,
            {
                "type": row.get("sample_type") or "sample",
                "status": row.get("sample_status") or "",
                "reasons": [],
            },
        )
        if cell["reason"] and cell["reason"] not in info["reasons"]:
            info["reasons"].append(cell["reason"])

    samples = []
    for sample in sample_names:
        info = sample_info.get(sample)
        if info is None:
            continue
        flag = info["status"] if info["status"] in SAMPLE_FLAGS else None
        samples.append(
            {
                "name": label(sample),
                "type": info["type"],
                "negative": info["type"] == "negative",
                "flag": flag,
                "reasons": "; ".join(info["reasons"]),
                "cells": [cells.get((sample, t)) for t in target_names],
            }
        )

    strip = []
    for name in target_names:
        mine = [c for (s, t), c in cells.items() if t == name]
        strip.append(
            {
                "name": name,
                "n_detected": sum(c["status"] == "detected" for c in mine),
                "n_trace": sum(c["status"] == "trace" for c in mine),
                "n_analysed": sum(c["analysed"] for c in mine),
                "n_lineage": sum(c["lineage"] is not None for c in mine),
                "n_samples": len(mine),
                "link": f"{quote(name)}/report.html" if per_target[name]["has_report"] else None,
            }
        )

    run_info = {}
    for name in target_names:
        path = manifest.get(name, {}).get("run_info")
        if path and Path(path).is_file():
            try:
                run_info = json.loads(Path(path).read_text(encoding="utf-8"))
                break
            except (OSError, ValueError):
                run_info = {}
    started = run_info.get("start") if isinstance(run_info.get("start"), str) else ""
    panel_spec = run_info.get("panel") if isinstance(run_info.get("panel"), str) else ""
    header = {
        "pipeline_version": args.pipeline_version or run_info.get("pipeline_version") or "",
        "mode": args.mode or run_info.get("mode") or "",
        # a path would leak the user's directories: only the file name (or registry name) is shown
        "panel": basename(args.panel or panel_spec),
        "date": started[:10] or date.today().isoformat(),
        "n_samples": len(samples),
        "n_targets": len(target_names),
        "anonymize": args.anonymize,
    }
    ctx = {
        "h": header,
        "targets": strip,
        "samples": samples,
        "n_unanalysed_targets": sum(1 for t in strip if t["n_analysed"] == 0),
        "n_flagged": sum(1 for s in samples if s["flag"]),
    }
    return ctx, id_map


def render(ctx, template_path=DEFAULT_TEMPLATE):
    env = jinja2.Environment(
        loader=jinja2.FileSystemLoader(str(Path(template_path).parent)),
        autoescape=True,
        undefined=jinja2.StrictUndefined,
    )
    return env.get_template(Path(template_path).name).render(**ctx)


def parse_args(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--targets", required=True, help="summary/targets.tsv")
    parser.add_argument("--inputs", help="manifest of the per-target files (target, role, path)")
    parser.add_argument("--template", default=str(DEFAULT_TEMPLATE))
    parser.add_argument("--pipeline-version", default="")
    parser.add_argument("--mode", default="")
    parser.add_argument("--panel", default="", help="panel spec (only its file name is shown)")
    parser.add_argument("--anonymize", choices=("none", "paths", "full"), default="none")
    parser.add_argument(
        "--anonymize-map-out", help="with --anonymize full: original -> pseudonym table"
    )
    parser.add_argument("--out", required=True)
    return parser.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    targets_rows = read_tsv(args.targets)
    ctx, id_map = build_context(targets_rows, read_manifest(args.inputs), args)
    try:
        html = render(ctx, args.template)
    except jinja2.TemplateError as err:
        print(f"build_overview.py: template error: {err}", file=sys.stderr)
        return 2
    Path(args.out).write_text(html, encoding="utf-8")
    if args.anonymize == "full" and args.anonymize_map_out:
        with open(args.anonymize_map_out, "w", newline="", encoding="utf-8") as handle:
            writer = csv.writer(handle)
            writer.writerow(["type", "original", "pseudonym"])
            writer.writerows(("sample", k, v) for k, v in sorted(id_map.items()))
    return 0


if __name__ == "__main__":
    sys.exit(main())
