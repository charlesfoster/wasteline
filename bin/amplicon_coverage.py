#!/usr/bin/env python3
"""Per-amplicon coverage and dropout from a Freyja depths file and a primer BED.

Amplicons are derived from the primer BED (no amplicon BED is needed). Each primer name is split into an
amplicon key and a direction (left/forward or right/reverse). Primers sharing a key form one amplicon and
alternate primers (`_alt1_F`, `QIAseq_221-2_LEFT`, several swift/xgen left primers on one right primer) are
merged: the amplicon spans the leftmost left-primer start to the rightmost right-primer end (outer
coordinates, primers included).

Output columns: sample, amplicon, start, end, mean_depth, min_depth, pct_ge_min_depth, dropped.
`start`/`end` are 1-based inclusive genome positions (the same coordinates as the depths file), i.e.
BED start + 1 and BED end. `dropped` is `true` when mean_depth < --min-depth.

The depths file has four tab-separated columns and no header: contig, 1-based position, ref base, depth.
Positions missing from the file are treated as depth 0.
"""

import argparse
import contextlib
import re
import sys
from collections import OrderedDict

COLUMNS = [
    "sample",
    "amplicon",
    "start",
    "end",
    "mean_depth",
    "min_depth",
    "pct_ge_min_depth",
    "dropped",
]

_DIRECTION_RE = re.compile(r"^(.*?)_?(LEFT|RIGHT|F|R)$")
_LEFT = {"LEFT", "F"}


def read_depths(path):
    """Return OrderedDict contig -> list of depths (index = position - 1; gaps are 0)."""
    contigs = OrderedDict()
    with open(path) as handle:
        for line_number, line in enumerate(handle, start=1):
            line = line.rstrip("\r\n")
            if not line.strip():
                continue
            fields = line.split("\t")
            if len(fields) < 4:
                raise ValueError(
                    f"{path}: line {line_number} has {len(fields)} columns, expected at least 4"
                )
            try:
                pos = int(fields[1])
                depth = int(float(fields[3]))
            except ValueError as exc:
                raise ValueError(f"{path}: line {line_number}: non-numeric position/depth") from exc
            if pos < 1:
                raise ValueError(f"{path}: line {line_number}: position must be >= 1")
            arr = contigs.setdefault(fields[0], [])
            if len(arr) < pos:
                arr.extend([0] * (pos - len(arr)))
            arr[pos - 1] = depth
    return contigs


def split_primer_name(name, strand=None):
    """Return (amplicon_key, 'L'|'R') for a primer name.

    Handles `X_LEFT/X_RIGHT`, `X_F/X_R`, `X_altN_F`, `QIAseq_221-2_LEFT`, and the swift/xgen
    `covid19genome_..._88_F` / `..._88F` styles. Falls back to the BED strand column when the name has no
    direction suffix.
    """
    match = _DIRECTION_RE.match(name)
    if match:
        key, direction = match.group(1), "L" if match.group(2) in _LEFT else "R"
    elif strand in ("+", "-"):
        key, direction = name, "L" if strand == "+" else "R"
    else:
        raise ValueError(f"cannot determine primer direction from name {name!r}")
    key = re.sub(r"_alt\d+$", "", key)
    key = re.sub(r"(?<=\d)-\d+$", "", key)
    return key, direction


def read_primers(path):
    """Parse a primer BED into a list of (contig, start0, end, key, direction)."""
    primers = []
    with open(path) as handle:
        for line_number, line in enumerate(handle, start=1):
            line = line.rstrip("\r\n")
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            fields = line.split("\t")
            if len(fields) < 4:
                raise ValueError(f"{path}: line {line_number} has fewer than 4 columns")
            try:
                start, end = int(fields[1]), int(fields[2])
            except ValueError as exc:
                raise ValueError(f"{path}: line {line_number}: non-integer coordinates") from exc
            strand = fields[5] if len(fields) > 5 else None
            key, direction = split_primer_name(fields[3], strand)
            primers.append((fields[0], start, end, key, direction))
    return primers


def _trailing_number(key):
    match = re.search(r"(\d+)_*$", key)
    return match.group(1) if match else None


def derive_amplicons(primers):
    """Return a list of (contig, amplicon, start0, end) sorted by position.

    Resolution order for primers whose key lacks an opposite-strand partner ("orphans"):
    1. a lone left orphan and a lone right orphan on the same contig that share the trailing amplicon
       number are paired (e.g. midnightv2 `SARSCoV_1200_22_F` with `nCoV_1200_22_R`);
    2. remaining orphans join the paired amplicon whose left start (left orphan) or right end (right
       orphan) is closest, extending its span only outwards (swift/xgen `M1`/`M2` start primers; the
       midnightv2 `ARTIC_V4_1200_71` alternate right primer);
    3. otherwise they are dropped (with a warning).
    """
    groups = OrderedDict()
    for contig, start, end, key, direction in primers:
        entry = groups.setdefault((contig, key), {"L": [], "R": []})
        entry[direction].append((start, end))

    paired = OrderedDict()
    orphans = OrderedDict()
    for gk, entry in groups.items():
        (paired if entry["L"] and entry["R"] else orphans)[gk] = entry

    # 1. pair orphans by trailing number
    by_num = OrderedDict()
    for (contig, key), entry in orphans.items():
        num = _trailing_number(key)
        if num is None:
            continue
        side = "L" if entry["L"] else "R"
        by_num.setdefault((contig, num), {"L": [], "R": []})[side].append((contig, key))
    resolved = set()
    for (contig, num), sides in by_num.items():
        if len(sides["L"]) == 1 and len(sides["R"]) == 1:
            lk, rk = sides["L"][0], sides["R"][0]
            paired[lk] = {
                "L": orphans[lk]["L"],
                "R": orphans[rk]["R"],
            }
            resolved.update([lk, rk])

    spans = []
    for (contig, key), entry in paired.items():
        spans.append(
            [
                contig,
                key,
                min(s for s, _ in entry["L"]),
                max(e for _, e in entry["R"]),
            ]
        )

    # 2. nearest-neighbour fallback for remaining orphans
    for gk, entry in orphans.items():
        if gk in resolved:
            continue
        contig, key = gk
        candidates = [s for s in spans if s[0] == contig]
        target = None
        if entry["L"]:
            lstart = min(s for s, _ in entry["L"])
            ahead = [s for s in candidates if s[3] > lstart]
            if ahead:
                target = min(ahead, key=lambda s: abs(s[2] - lstart))
                target[2] = min(target[2], lstart)
        else:
            rend = max(e for _, e in entry["R"])
            behind = [s for s in candidates if s[2] < rend]
            if behind:
                target = min(behind, key=lambda s: abs(s[3] - rend))
                target[3] = max(target[3], rend)
        if target is None:
            print(
                f"warning: primer group {key!r} on {contig} has no partner; ignored",
                file=sys.stderr,
            )

    spans = [tuple(s) for s in spans if s[3] > s[2]]
    spans.sort(key=lambda s: (s[0], s[2], s[3]))
    return spans


def amplicon_rows(sample, depths_by_contig, amplicons, min_depth):
    """Compute output rows (lists of strings) for each amplicon."""
    single_depth_contig = next(iter(depths_by_contig)) if len(depths_by_contig) == 1 else None
    bed_contigs = {a[0] for a in amplicons}
    rows = []
    for contig, name, start0, end in amplicons:
        if contig in depths_by_contig:
            arr = depths_by_contig[contig]
        elif single_depth_contig is not None and len(bed_contigs) == 1:
            arr = depths_by_contig[single_depth_contig]
        else:
            arr = []
        window = [arr[i] if i < len(arr) else 0 for i in range(start0, end)]
        n = len(window)
        mean = sum(window) / n
        rows.append(
            [
                sample,
                name,
                str(start0 + 1),
                str(end),
                f"{mean:.2f}",
                str(min(window)),
                f"{100.0 * sum(1 for d in window if d >= min_depth) / n:.2f}",
                "true" if mean < min_depth else "false",
            ]
        )
    return rows


def write_table(path, rows):
    handle = sys.stdout if path == "-" else open(path, "w")  # noqa: SIM115
    with contextlib.nullcontext(handle) if path == "-" else handle as out:
        out.write("\t".join(COLUMNS) + "\n")
        for row in rows:
            out.write("\t".join(row) + "\n")


def main(argv=None):
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("--depths", required=True, help="Freyja depths file (4 columns, no header)")
    parser.add_argument("--bed", required=True, help="primer BED")
    parser.add_argument("--sample", required=True)
    parser.add_argument(
        "--min-depth", type=int, default=10, help="dropout/coverage threshold (default 10)"
    )
    parser.add_argument("-o", "--output", default="-")
    args = parser.parse_args(argv)

    try:
        depths = read_depths(args.depths)
        amplicons = derive_amplicons(read_primers(args.bed))
    except (OSError, ValueError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2
    if not amplicons:
        print(f"warning: no amplicons could be derived from {args.bed}", file=sys.stderr)
    bed_contigs = {a[0] for a in amplicons}
    if (
        depths
        and bed_contigs
        and not (bed_contigs & set(depths))
        and not (len(bed_contigs) == 1 and len(depths) == 1)
    ):
        print(
            f"error: BED contigs {sorted(bed_contigs)} not found in depths contigs {sorted(depths)}",
            file=sys.stderr,
        )
        return 2
    if bed_contigs and not (bed_contigs & set(depths)):
        print(
            f"warning: BED contig {sorted(bed_contigs)} differs from depths contig {sorted(depths)}; assuming same reference",
            file=sys.stderr,
        )
    write_table(args.output, amplicon_rows(args.sample, depths, amplicons, args.min_depth))
    return 0


if __name__ == "__main__":
    sys.exit(main())
