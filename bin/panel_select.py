#!/usr/bin/env python3
"""Panel typing: which targets of a panel does a sample contain? (standard library only)

  prepare   a resolved panel table (target,fasta,marker_contig,family) -> panel.fasta + contig_map.tsv
  call      <id>.panel_stats.tsv (PANEL_ASSIGN) + contig_map.tsv + thresholds -> <id>.targets.tsv

`prepare` concatenates the FASTA of every target into one competitive-mapping reference. With more than one
target each contig is renamed `<target>|<contig>`; with one target the names are kept. `call` ranks the targets
of each family and gives every target of the panel exactly one row (detected | trace | not_detected), also for a
sample with no panel reads: the file is never empty. The thresholds are uncalibrated starting values.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import os
import sys

MAP_COLUMNS = ["panel_contig", "target", "contig", "is_marker", "family"]
STATS_COLUMNS = [
    "sample",
    "target",
    "fragments",
    "ambiguous_fragments",
    "aligned_bases",
    "mean_mapq",
    "breadth_pct",
    "mean_depth",
    "panel_fragments",
]
TARGET_COLUMNS = [
    "sample",
    "sample_type",
    "target",
    "family",
    "status",
    "reason",
    "sample_status",
    "mixed",
    "secondary_targets",
    "fragments",
    "ambiguous_fragments",
    "fraction_of_family",
    "breadth_pct",
    "mean_depth",
    "mean_mapq",
    "second_target",
    "margin",
    "rescue",
    "analysed",
]


class PanelError(Exception):
    """A problem with the panel or its inputs; main() prints it and exits 2."""


# --- prepare ---------------------------------------------------------------------------------------------------


def open_text(path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8")
    return open(path, encoding="utf-8")


def locate_fasta(cell, fasta_dir):
    """The FASTA of a table row: the path as written, else its file name in fasta_dir (staged by Nextflow)."""
    if os.path.isfile(cell):
        return cell
    staged = os.path.join(fasta_dir, os.path.basename(cell))
    if os.path.isfile(staged):
        return staged
    raise PanelError(f"FASTA not found: {cell} (also looked for {staged})")


def read_fasta(path):
    """[(contig, [sequence lines])] with the contig name = first word of the header."""
    records = []
    with open_text(path) as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                words = line[1:].split()
                if not words:
                    raise PanelError(f"{path}: a record has an empty header")
                records.append((words[0], []))
            elif records:
                records[-1][1].append(line)
            else:
                raise PanelError(f"{path}: sequence before the first '>' header")
    if not records:
        raise PanelError(f"{path}: no FASTA records")
    return records


def read_panel_table(path):
    with open(path, newline="", encoding="utf-8") as handle:
        rows = [
            {k.strip(): (v or "").strip() for k, v in row.items() if k}
            for row in csv.DictReader(handle)
        ]
    if not rows:
        raise PanelError(f"{path}: the panel table has no rows")
    for column in ("target", "fasta"):
        if column not in rows[0]:
            raise PanelError(f"{path}: required column '{column}' is missing")
    return rows


def prepare(panel, fasta_dir, out_fasta, out_map):
    rows = read_panel_table(panel)
    targets = [row["target"] for row in rows]
    for target in targets:
        if not target or "|" in target or any(c.isspace() for c in target):
            raise PanelError(
                f"invalid target name {target!r}: it must be non-empty, without whitespace or '|'"
            )
    duplicated = sorted({t for t in targets if targets.count(t) > 1})
    if duplicated:
        raise PanelError(f"target listed more than once in the panel: {', '.join(duplicated)}")
    keep_names = len(rows) == 1  # one target: the contig names stay as they are (legacy naming)

    map_rows = []
    seen = set()
    with open(out_fasta, "w", encoding="utf-8") as fasta_out:
        for row in rows:
            target = row["target"]
            records = read_fasta(locate_fasta(row["fasta"], fasta_dir))
            names = [name for name, _ in records]
            repeated = sorted({n for n in names if names.count(n) > 1})
            if repeated:
                raise PanelError(
                    f"target {target}: contig name {repeated[0]!r} appears twice in {row['fasta']}"
                )
            markers = [
                m.strip()
                for m in row.get("marker_contig", "").replace(";", ",").split(",")
                if m.strip()
            ]
            for marker in markers:
                if marker not in names:
                    raise PanelError(
                        f"target {target}: marker_contig {marker!r} is not a contig of {row['fasta']}"
                    )
            family = row.get("family", "") or target
            for name, lines in records:
                panel_contig = name if keep_names else f"{target}|{name}"
                if (
                    panel_contig in seen
                ):  # cannot happen with distinct targets, but a clash would corrupt the BAM
                    raise PanelError(f"panel contig name {panel_contig!r} is not unique")
                seen.add(panel_contig)
                is_marker = (name in markers) if markers else True
                map_rows.append(
                    [panel_contig, target, name, "true" if is_marker else "false", family]
                )
                fasta_out.write(f">{panel_contig}\n")
                fasta_out.write("\n".join(lines) + "\n")
    with open(out_map, "w", encoding="utf-8") as handle:
        handle.write("\t".join(MAP_COLUMNS) + "\n")
        handle.writelines("\t".join(map_row) + "\n" for map_row in map_rows)
    return map_rows


# --- call ------------------------------------------------------------------------------------------------------


def read_tsv(path):
    with open(path, newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def read_contig_map(path):
    """([targets in panel order], {target: family})."""
    rows = read_tsv(path)
    if not rows or any(c not in rows[0] for c in MAP_COLUMNS):
        raise PanelError(f"{path}: not a contig map (columns {', '.join(MAP_COLUMNS)})")
    families = {}
    for row in rows:
        families.setdefault(row["target"], row["family"] or row["target"])
    return list(families), families


def read_stats(path, targets):
    """{target: stats dict}; a target absent from the file is all zero. An empty/missing file is all zero too."""
    stats = {}
    panel_fragments = 0
    if os.path.isfile(path) and os.path.getsize(path) > 0:
        for row in read_tsv(path):
            if row["target"] not in targets:
                raise PanelError(f"{path}: target {row['target']!r} is not in the contig map")
            stats[row["target"]] = {
                "fragments": int(row["fragments"]),
                "ambiguous_fragments": int(row["ambiguous_fragments"]),
                "aligned_bases": int(row["aligned_bases"]),
                "mean_mapq": float(row["mean_mapq"]),
                "breadth_pct": float(row["breadth_pct"]),
                "mean_depth": float(row["mean_depth"]),
            }
            panel_fragments = max(panel_fragments, int(row["panel_fragments"]))
    zero = {
        "fragments": 0,
        "ambiguous_fragments": 0,
        "aligned_bases": 0,
        "mean_mapq": 0.0,
        "breadth_pct": 0.0,
        "mean_depth": 0.0,
    }
    return {t: stats.get(t, dict(zero)) for t in targets}, panel_fragments


def call(
    sample,
    sample_type,
    stats,
    panel_fragments,
    families,
    min_fragments,
    min_breadth_pct,
    family_share_pct,
):
    """One row (dict of TARGET_COLUMNS) per target, in panel order."""
    rows = {}
    members = {}
    for target, family in families.items():
        members.setdefault(family, []).append(target)
    for family, names in members.items():
        total = sum(stats[t]["fragments"] for t in names)
        ranked = sorted(
            names,
            key=lambda t: (
                -stats[t]["fragments"],
                -stats[t]["aligned_bases"],
                -stats[t]["breadth_pct"],
                t,
            ),
        )
        for rank, target in enumerate(ranked):
            s = stats[target]
            share = s["fragments"] / total if total else 0.0
            other = ranked[1] if rank == 0 else ranked[0]
            has_second = len(ranked) > 1 and stats[other]["fragments"] > 0 and s["fragments"] > 0
            failed = []
            if s["fragments"] < min_fragments:
                failed.append(f"low_fragments({s['fragments']}<{min_fragments:g})")
            if s["breadth_pct"] < min_breadth_pct:
                failed.append(f"low_breadth({s['breadth_pct']:.2f}<{min_breadth_pct:g})")
            if rank > 0 and share * 100 < family_share_pct:
                failed.append(f"low_family_share({share * 100:.2f}<{family_share_pct:g})")
            if s["fragments"] == 0:
                status, reason = (
                    "not_detected",
                    "no_panel_reads" if panel_fragments == 0 else "no_fragments",
                )
            elif failed:
                status, reason = "trace", ";".join(failed)
            else:
                status, reason = "detected", ""
            second_share = stats[other]["fragments"] / total if total else 0.0
            rows[target] = {
                "sample": sample,
                "sample_type": sample_type,
                "target": target,
                "family": family,
                "status": status,
                "reason": reason,
                "mixed": "false",
                "secondary_targets": "",
                "fragments": s["fragments"],
                "ambiguous_fragments": s["ambiguous_fragments"],
                "fraction_of_family": f"{share:.4f}",
                "breadth_pct": f"{s['breadth_pct']:.2f}",
                "mean_depth": f"{s['mean_depth']:.2f}",
                "mean_mapq": f"{s['mean_mapq']:.1f}" if s["fragments"] else "NA",
                "second_target": other if has_second else "",
                "margin": f"{share - second_share:.4f}" if has_second else "",
                "rescue": "",
                "analysed": "true" if status == "detected" else "false",
            }
    for names in members.values():
        detected = [t for t in names if rows[t]["status"] == "detected"]
        if len(detected) >= 2:
            for t in detected:
                rows[t]["mixed"] = "true"
                rows[t]["secondary_targets"] = ";".join(d for d in detected if d != t)
    if panel_fragments == 0:
        sample_status = "qc_fail"
    elif any(r["status"] == "detected" for r in rows.values()):
        sample_status = "ok"
    else:
        sample_status = "no_target_detected"
    out = []
    for target in families:
        rows[target]["sample_status"] = sample_status
        out.append(rows[target])
    return out


def write_targets(path, rows):
    with open(path, "w", encoding="utf-8") as handle:
        handle.write("\t".join(TARGET_COLUMNS) + "\n")
        handle.writelines("\t".join(str(row[c]) for c in TARGET_COLUMNS) + "\n" for row in rows)


# --- CLI -------------------------------------------------------------------------------------------------------


def build_parser():
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    sub = parser.add_subparsers(dest="command", required=True)
    p = sub.add_parser("prepare", help="resolved panel table -> panel.fasta + contig_map.tsv")
    p.add_argument(
        "--panel", required=True, help="resolved panel CSV (target,fasta,marker_contig,family,...)"
    )
    p.add_argument(
        "--fasta-dir",
        default=".",
        help="where a FASTA is looked for by file name when its path does not exist",
    )
    p.add_argument("--out-fasta", required=True)
    p.add_argument("--out-map", required=True)
    c = sub.add_parser("call", help="panel_stats.tsv -> targets.tsv")
    c.add_argument("--sample", required=True)
    c.add_argument("--sample-type", default="sample")
    c.add_argument("--stats", required=True)
    c.add_argument("--contig-map", required=True)
    c.add_argument("--min-fragments", type=float, default=100)
    c.add_argument("--min-breadth-pct", type=float, default=5)
    c.add_argument("--family-share-pct", type=float, default=5)
    c.add_argument("--output", required=True)
    return parser


def main(argv=None):
    args = build_parser().parse_args(argv)
    try:
        if args.command == "prepare":
            prepare(args.panel, args.fasta_dir, args.out_fasta, args.out_map)
        else:
            targets, families = read_contig_map(args.contig_map)
            stats, panel_fragments = read_stats(args.stats, set(targets))
            rows = call(
                args.sample,
                args.sample_type,
                stats,
                panel_fragments,
                families,
                args.min_fragments,
                args.min_breadth_pct,
                args.family_share_pct,
            )
            write_targets(args.output, rows)
    except (PanelError, OSError, KeyError, ValueError) as err:
        print(f"ERROR: {err}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    sys.exit(main())
