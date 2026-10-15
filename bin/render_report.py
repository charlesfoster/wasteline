#!/usr/bin/env python3
"""Render the self-contained wasteline HTML report from `report_data.json`.

The page is fully offline: inline CSS, inline SVG built here (only numbers and class names are
computed in Python; every string that comes from the data is emitted by the Jinja2 template with
autoescape on), and a few lines of vanilla JS for table sort/filter/toggles. No external URLs.

Requires jinja2 (standard library otherwise). Schema: docs/report_data.md.
"""

import argparse
import json
import math
import sys
from collections import Counter
from datetime import date
from pathlib import Path

import jinja2

N_COLORS = 10  # categorical classes c0..c9; everything else is the grey `co` class
HEAT_BUCKETS = 9
MAX_HEAT_COLS = {"groups": 20, "lineages": 25}
BAR_UNITS = 1000.0
MIN_SEG_UNITS = 2.0  # segments narrower than this are merged into "Other" for drawing
MAX_BC_ROWS = 100
MAX_BC_SAMPLES = 50
DEFAULT_TEMPLATE = (
    Path(__file__).resolve().parent.parent / "assets" / "templates" / "report.html.j2"
)


# ---------------------------------------------------------------------------------------------
# Normalisation: the renderer must cope with minimal / partial report_data.json files
# ---------------------------------------------------------------------------------------------
def _as_dict(value):
    return value if isinstance(value, dict) else {}


def _as_list(value):
    return value if isinstance(value, list) else []


def normalise(data):
    data = _as_dict(data)
    out = dict(data)
    out["schema_version"] = data.get("schema_version", 0)
    out["mode"] = data.get("mode") or "unknown"
    out["provenance"] = _as_dict(data.get("provenance"))
    out["header"] = _as_dict(data.get("header"))
    out["params_summary"] = _as_dict(data.get("params_summary"))
    out["thresholds"] = _as_dict(data.get("thresholds"))
    out["run_flags"] = [f for f in _as_list(data.get("run_flags")) if isinstance(f, dict)]
    samples = []
    for rec in _as_list(data.get("samples")):
        if isinstance(rec, dict) and rec.get("sample") is not None:
            rec = dict(rec)
            rec["sample"] = str(rec["sample"])
            rec["sample_type"] = rec.get("sample_type") or "sample"
            rec["qc_reasons"] = [str(r) for r in _as_list(rec.get("qc_reasons"))]
            samples.append(rec)
    out["samples"] = samples
    out["demix"] = _as_dict(data.get("demix"))
    out["amplicons"] = _as_dict(data.get("amplicons"))
    out["coverage_tracks"] = _as_dict(data.get("coverage_tracks")) or None
    out["mutations_of_interest"] = data.get("mutations_of_interest")
    out["consensus"] = _as_list(data.get("consensus"))
    out["barcode_comparison"] = _as_dict(data.get("barcode_comparison")) or None
    out["time_series"] = _as_dict(data.get("time_series")) or None
    out["negative_controls"] = [
        c for c in _as_list(data.get("negative_controls")) if isinstance(c, dict)
    ]
    return out


def _num(value):
    return value if isinstance(value, (int, float)) and not isinstance(value, bool) else None


def fmt_int(value):
    value = _num(value)
    return "–" if value is None else f"{round(value):,}"


def fmt_dec(value, digits=1):
    value = _num(value)
    return "–" if value is None else f"{value:,.{digits}f}"


def fmt_pct(value, digits=1):
    """A 0-100 percentage."""
    value = _num(value)
    return "–" if value is None else f"{value:.{digits}f}"


def fmt_frac(value, digits=1):
    """A 0-1 fraction shown as a percentage."""
    value = _num(value)
    return "–" if value is None else f"{value * 100:.{digits}f}"


def sortval(value):
    value = _num(value)
    return "" if value is None else repr(value)


def short(value, length=26):
    value = "" if value is None else str(value)
    return value if len(value) <= length else value[: length - 1] + "…"


def dash(value):
    return "–" if value in (None, "") else value


# ---------------------------------------------------------------------------------------------
# View models
# ---------------------------------------------------------------------------------------------
def _palette(totals):
    ranked = [name for name, _ in totals.most_common() if name != "Other"][:N_COLORS]
    return {name: i for i, name in enumerate(ranked)}


def _heat_bucket(value):
    if value is None or value < 0.0005:
        return 0
    return 1 + min(HEAT_BUCKETS - 1, int(math.sqrt(min(value, 1.0)) * HEAT_BUCKETS))


def _ci(info, name):
    """' (95% interval lo-hi%)' for a lineage with a `freyja boot` interval, else ''."""
    pair = _as_dict(info.get("intervals")).get(name)
    if (
        isinstance(pair, (list, tuple))
        and len(pair) == 2
        and all(_num(v) is not None for v in pair)
    ):
        return f" (95% interval {pair[0] * 100:.1f}\u2013{pair[1] * 100:.1f}%)"
    return ""


def lineage_view(samples, demix, key):
    """Stacked bars + heatmap for `key` in ('groups', 'lineages') over samples with status ok."""
    entries = []
    totals = Counter()
    for rec in samples:
        info = _as_dict(demix.get(rec["sample"]))
        pairs = [
            (str(p[0]), float(p[1]))
            for p in _as_list(info.get(key))
            if isinstance(p, (list, tuple)) and len(p) == 2 and _num(p[1]) is not None
        ]
        entries.append((rec, info, pairs))
        if info.get("status") == "ok":
            for name, frac in pairs:
                totals[name] += frac
    palette = _palette(totals)
    legend = [{"name": n, "cls": f"c{i}"} for n, i in palette.items()]
    if any(name not in palette for name in totals):
        legend.append({"name": "Other", "cls": "co"})

    rows = []
    for rec, info, pairs in entries:
        row = {
            "sample": rec["sample"],
            "ntc": rec["sample_type"] == "negative",
            "qc": rec.get("qc_status") or "NA",
            "status": info.get("status"),
            "reason": info.get("reason") or "",
            "segs": [],
        }
        if info.get("status") == "ok" and pairs:
            total = max(sum(f for _, f in pairs), 1.0)
            x = 0.0
            other = []
            for name, frac in pairs:
                width = frac / total * BAR_UNITS
                if name in palette and width >= MIN_SEG_UNITS:
                    row["segs"].append(
                        {
                            "x": round(x, 1),
                            "w": round(width, 1),
                            "cls": f"c{palette[name]}",
                            "title": f"{name}: {frac * 100:.1f}%{_ci(info, name)}",
                        }
                    )
                else:
                    other.append((name, frac, x, width))
                x += width
            if other:
                # Drawn as separate rects (their positions are not contiguous) sharing one title.
                text = "Other: " + ", ".join(f"{n} {f * 100:.1f}%" for n, f, _, _ in other[:6])
                if len(other) > 6:
                    text += f", … ({len(other)} in total)"
                for _, _, ox, ow in other:
                    row["segs"].append(
                        {"x": round(ox, 1), "w": round(ow, 1), "cls": "co", "title": text}
                    )
        rows.append(row)

    ok_entries = [(r, i, p) for r, i, p in entries if i.get("status") == "ok" and p]
    columns = [n for n, _ in totals.most_common() if n != "Other"][: MAX_HEAT_COLS[key]]
    if "Other" in totals:
        columns.append("Other")
    heat_rows = []
    for rec, _, pairs in ok_entries:
        values = dict(pairs)
        cells = []
        for name in columns:
            value = values.get(name, 0.0)
            cells.append(
                {
                    "cls": f"h{_heat_bucket(value)}",
                    "text": f"{value * 100:.0f}" if value >= 0.01 else "",
                    "title": f"{name}: {value * 100:.1f}%{_ci(_as_dict(demix.get(rec['sample'])), name)}",
                }
            )
        heat_rows.append(
            {
                "sample": rec["sample"],
                "ntc": rec["sample_type"] == "negative",
                "qc": rec.get("qc_status") or "NA",
                "cells": cells,
            }
        )
    return {
        "legend": legend,
        "rows": rows,
        "columns": columns,
        "heat_rows": heat_rows,
        "n_ok": len(ok_entries),
        "n_hidden": len(entries) - len(ok_entries),
        "palette": palette,
    }


def coverage_view(data):
    tracks = data["coverage_tracks"]
    if not tracks or not _as_dict(tracks.get("samples")):
        return None
    n_bins = int(tracks.get("n_bins") or 0)
    bin_size = max(1, int(tracks.get("bin_size") or 1))
    series = tracks["samples"]
    n_bins = n_bins or max(len(v) for v in series.values())
    ymax = max(
        [100.0] + [float(v) for track in series.values() for v in track if _num(v) is not None]
    )
    ymax = min(ymax, 1e5)
    scale = math.log10(1 + ymax)
    thresholds = data["thresholds"]
    amp_min = _num(thresholds.get("amplicon_min_depth")) or 10
    cons_min = _num(thresholds.get("consensus_min_depth"))

    def ypos(depth):
        depth = min(max(depth, 0.0), ymax)
        return round(100 - 100 * math.log10(1 + depth) / scale)

    by_sample = {rec["sample"]: rec for rec in data["samples"]}
    cards = []
    for sample in sorted(series):
        track = [float(v) if _num(v) is not None else 0.0 for v in series[sample]]
        ys = [ypos(v) for v in track]
        parts = ["M0,100"]
        for i, y in enumerate(ys):
            if 0 < i < len(ys) - 1 and ys[i - 1] == y == ys[i + 1]:
                continue
            parts.append(f"L{i},{y}")
        parts.append(f"L{max(len(ys) - 1, 0)},100Z")
        drops = []
        for amp in _as_list(data["amplicons"].get(sample)):
            if (
                amp.get("dropped")
                and _num(amp.get("start")) is not None
                and _num(amp.get("end")) is not None
            ):
                # start/end are 1-based inclusive (PLAN 6a); bins are 0-based.
                x = (amp["start"] - 1) / bin_size
                drops.append(
                    {
                        "x": round(x, 1),
                        "w": round(max((amp["end"] - amp["start"] + 1) / bin_size, 0.5), 1),
                        "title": f"{amp.get('amplicon')}: dropped (mean depth "
                        f"{fmt_dec(amp.get('mean_depth'))}x)",
                    }
                )
        rec = by_sample.get(sample, {})
        cards.append(
            {
                "sample": sample,
                "ntc": rec.get("sample_type") == "negative",
                "qc": rec.get("qc_status") or "NA",
                "mean_depth": rec.get("mean_depth"),
                "pct10": rec.get("pct_genome_ge_10x"),
                "dropped": len(drops),
                "path": "".join(parts),
                "drops": drops,
            }
        )
    return {
        "cards": cards,
        "n_bins": n_bins,
        "bin_size": bin_size,
        "genome_length": tracks.get("genome_length"),
        "ymax": int(ymax),
        "y_amp": ypos(amp_min),
        "amp_min": amp_min,
        "y_cons": ypos(cons_min) if cons_min and cons_min != amp_min else None,
        "cons_min": cons_min,
    }


def _runs(indices):
    runs = []
    for i in indices:
        if runs and runs[-1][0] + runs[-1][1] == i:
            runs[-1][1] += 1
        else:
            runs.append([i, 1])
    return runs


def amplicon_view(data):
    amplicons = {s: rows for s, rows in data["amplicons"].items() if _as_list(rows)}
    if not amplicons:
        return None
    amp_min = _num(data["thresholds"].get("amplicon_min_depth")) or 10
    order = {}
    for rows in amplicons.values():
        for amp in rows:
            key = (amp.get("start") or 0, amp.get("end") or 0, str(amp.get("amplicon")))
            order[str(amp.get("amplicon"))] = min(order.get(str(amp.get("amplicon")), key), key)
    names = sorted(order, key=lambda n: order[n])
    index = {name: i for i, name in enumerate(names)}
    by_sample = {rec["sample"]: rec for rec in data["samples"]}
    rows_out = []
    dropped_count = [0] * len(names)
    considered = 0
    for sample in sorted(amplicons):
        rec = by_sample.get(sample, {})
        ntc = rec.get("sample_type") == "negative"
        classes = {}
        values = [""] * len(names)
        for amp in amplicons[sample]:
            i = index[str(amp.get("amplicon"))]
            depth = _num(amp.get("mean_depth"))
            if depth is not None:
                values[i] = str(round(depth)) if depth >= 10 else f"{depth:.1f}"
            if amp.get("dropped"):
                cls = "d0" if not depth else "d1"
                if not ntc:
                    dropped_count[i] += 1
            elif depth is None:
                cls = "dn"
            elif depth < 10 * amp_min:
                cls = "d2"
            elif depth < 100 * amp_min:
                cls = "d3"
            elif depth < 1000 * amp_min:
                cls = "d4"
            else:
                cls = "d5"
            classes.setdefault(cls, []).append(i)
        if not ntc:
            considered += 1
        paths = []
        for cls in sorted(classes):
            d = "".join(
                f"M{start} 0h{length}v1h-{length}z" for start, length in _runs(sorted(classes[cls]))
            )
            paths.append({"cls": cls, "d": d})
        rows_out.append(
            {
                "sample": sample,
                "ntc": ntc,
                "qc": rec.get("qc_status") or "NA",
                "dropped": sum(1 for a in amplicons[sample] if a.get("dropped")),
                "paths": paths,
                "values": ",".join(values),
            }
        )
    freq = []
    if considered:
        for i, count in enumerate(dropped_count):
            if count:
                freq.append(
                    {
                        "x": i,
                        "h": round(20 * count / considered, 1),
                        "title": f"{names[i]}: dropped in {count} of {considered} samples",
                    }
                )
    ticks = [
        {"left": round(100 * (i + 0.5) / len(names), 2), "label": str(i + 1)}
        for i in range(0, len(names), 10)
    ]
    return {
        "names": names,
        "rows": rows_out,
        "freq": freq,
        "considered": considered,
        "n": len(names),
        "ticks": ticks,
        "amp_min": amp_min,
    }


def time_series_view(data, palette):
    ts = data["time_series"]
    if not ts:
        return None
    sites = _as_list(ts.get("sites"))
    all_dates = sorted({d for s in sites for d in _as_list(s.get("dates"))})
    try:
        days = [date.fromisoformat(d[:10]).toordinal() for d in all_dates]
    except ValueError:
        return None
    if len(set(days)) < 2:
        return None
    lo, hi = min(days), max(days)
    groups = [str(g) for g in _as_list(ts.get("groups"))]
    used = {palette[g] for g in groups if g in palette}
    free = iter(i for i in range(N_COLORS) if i not in used)
    cls = {}
    for g in groups:
        if g in palette:
            cls[g] = f"c{palette[g]}"
        else:
            nxt = next(free, None)
            cls[g] = f"c{nxt}" if nxt is not None else "co"
    vmax = max(
        [0.1]
        + [
            v
            for s in sites
            for ser in _as_list(s.get("series"))
            for v in _as_list(ser.get("values"))
            if _num(v) is not None
        ]
    )
    ymax = 1.0 if vmax > 0.5 else (math.ceil(vmax * 10) / 10)
    left, right, top, bottom, width, height = 38.0, 10.0, 8.0, 24.0, 420.0, 190.0

    def xpos(day):
        return left + (day - lo) / (hi - lo) * (width - left - right)

    def ypos(value):
        return top + (1 - min(value / ymax, 1.0)) * (height - top - bottom)

    panels = []
    for site in sites:
        dates = _as_list(site.get("dates"))
        counts = _as_list(site.get("n"))
        lines = []
        for ser in _as_list(site.get("series")):
            pts = []
            for k, (d, v) in enumerate(zip(dates, _as_list(ser.get("values")), strict=False)):
                if _num(v) is None:
                    continue
                n = counts[k] if k < len(counts) else 1
                pts.append(
                    {
                        "x": round(xpos(date.fromisoformat(d[:10]).toordinal()), 1),
                        "y": round(ypos(v), 1),
                        "title": f"{ser.get('group')} on {d}: {v * 100:.1f}% (n={n})",
                    }
                )
            lines.append(
                {
                    "group": ser.get("group"),
                    "cls": cls.get(str(ser.get("group")), "co"),
                    "points": " ".join(f"{p['x']},{p['y']}" for p in pts),
                    "dots": pts,
                }
            )
        panels.append({"site": site.get("site"), "lines": lines})
    yticks = [{"y": round(ypos(v), 1), "label": f"{v * 100:.0f}%"} for v in (0, ymax / 2, ymax)]
    return {
        "panels": panels,
        "legend": [{"name": g, "cls": cls[g]} for g in groups],
        "width": int(width),
        "height": int(height),
        "x0": left,
        "x1": width - right,
        "yb": round(ypos(0), 1),
        "yticks": yticks,
        "first": all_dates[0],
        "last": all_dates[-1],
        "n_dates": len(set(days)),
    }


def moi_view(data):
    rows = data["mutations_of_interest"]
    if rows is None:
        return None
    rows = [r for r in _as_list(rows) if isinstance(r, dict)]
    keys = sorted({(str(r.get("name") or ""), str(r.get("mutation") or "")) for r in rows})
    by = {(r["sample"], str(r.get("name") or ""), str(r.get("mutation") or "")): r for r in rows}
    min_depth = _num(data["thresholds"].get("consensus_min_depth")) or 10
    table = []
    for rec in data["samples"]:
        cells = []
        for name, mutation in keys:
            row = by.get((rec["sample"], name, mutation))
            if row is None:
                cells.append({"cls": "na", "text": "–", "title": "no data"})
                continue
            depth = _num(row.get("depth"))
            freq = _num(row.get("alt_freq"))
            if freq is None or not depth:
                cells.append({"cls": "na", "text": "–", "title": "no coverage at this site"})
                continue
            low = depth < min_depth
            cells.append(
                {
                    "cls": f"h{_heat_bucket(freq)}" + (" lowd" if low else ""),
                    "text": f"{freq * 100:.0f}%",
                    "title": f"{rec['sample']} / {name or mutation}: {freq * 100:.1f}% "
                    f"(alt {fmt_int(row.get('alt_depth'))} of {fmt_int(depth)} reads)"
                    + (" - low depth" if low else ""),
                }
            )
        table.append(
            {
                "sample": rec["sample"],
                "ntc": rec["sample_type"] == "negative",
                "qc": rec.get("qc_status") or "NA",
                "cells": cells,
            }
        )
    return {
        "columns": [{"name": n, "mutation": m} for n, m in keys],
        "rows": table,
        "min_depth": min_depth,
    }


def barcode_view(data):
    comp = data["barcode_comparison"]
    if not comp or not _as_list(comp.get("rows")):
        return None
    rows = [r for r in comp["rows"] if isinstance(r, dict)]
    per_sample = {}
    for r in rows:
        entry = per_sample.setdefault(r["sample"], {"sample": r["sample"], "tv": 0.0, "max": 0.0})
        entry["tv"] += abs(r.get("delta") or 0.0) / 2
        entry["max"] = max(entry["max"], abs(r.get("delta") or 0.0))
    summary = sorted(per_sample.values(), key=lambda e: -e["tv"])
    top = sorted(rows, key=lambda r: -abs(r.get("delta") or 0.0))[:MAX_BC_ROWS]
    return {
        "label_a": comp.get("label_a") or "A",
        "label_b": comp.get("label_b") or "B",
        "summary": summary[:MAX_BC_SAMPLES],
        "n_summary": len(summary),
        "rows": top,
        "n_rows": len(rows),
    }


def build_context(data):
    data = normalise(data)
    samples = data["samples"]
    demix = data["demix"]
    for rec in samples:
        info = _as_dict(demix.get(rec["sample"]))
        groups = _as_list(info.get("groups"))
        rec["_top"] = groups[0] if groups and info.get("status") == "ok" else None
        rec["_ntc"] = rec["sample_type"] == "negative"
        rec["_qc"] = rec.get("qc_status") or "NA"
        rec["_demix"] = info
    qc_counts = Counter(rec["_qc"] for rec in samples if not rec["_ntc"])
    demix_counts = Counter(
        (rec["_demix"].get("status") or "none") for rec in samples if not rec["_ntc"]
    )
    summ = lineage_view(samples, demix, "groups")
    full = lineage_view(samples, demix, "lineages")
    header = data["header"]
    barcode = header.get("barcode_version") or _as_dict(data["provenance"].get("freyja_data")).get(
        "barcode_version"
    )
    available = _as_dict(data.get("available"))
    context = {
        "d": data,
        "h": header,
        "barcode": barcode,
        "barcode_kind": (
            "missing"
            if not barcode
            else "unknown"
            if str(barcode).lower() == "unknown"
            else "custom"
            if str(barcode).lower() == "custom"
            else "ok"
        ),
        "freyja_data": _as_dict(data["provenance"].get("freyja_data")),
        "previous_runs": [
            p for p in _as_list(data["provenance"].get("previous_runs")) if isinstance(p, dict)
        ],
        "samples": samples,
        "n_samples": len(samples),
        "n_ntc": sum(1 for r in samples if r["_ntc"]),
        "qc_counts": qc_counts,
        "demix_counts": demix_counts,
        "has_qc": bool(available.get("qc", any(r.get("qc_status") for r in samples)))
        and bool(samples),
        "has_demix": bool(available.get("demix", bool(demix))),
        "flags": data["run_flags"],
        "flag_counts": Counter(f.get("level") for f in data["run_flags"]),
        "summ": summ,
        "full": full,
        "coverage": coverage_view(data),
        "amp": amplicon_view(data),
        "ts": time_series_view(data, summ["palette"]),
        "moi": moi_view(data),
        "bc": barcode_view(data),
        "controls": data["negative_controls"],
        "consensus": [r for r in data["consensus"] if isinstance(r, dict)],
        "mode": data["mode"],
        "mode_note": f"not available in mode '{data['mode']}'",
        "provenance_json": json.dumps(
            {k: v for k, v in data["provenance"].items() if k != "params"},
            indent=2,
            sort_keys=True,
            default=str,
        ),
        "params_json": json.dumps(
            _as_dict(data["provenance"].get("params")), indent=2, sort_keys=True, default=str
        ),
    }
    return context


def make_environment(template_path):
    template_path = Path(template_path)
    env = jinja2.Environment(
        loader=jinja2.FileSystemLoader(str(template_path.parent)),
        autoescape=True,
        trim_blocks=True,
        lstrip_blocks=True,
        keep_trailing_newline=True,
    )
    env.filters.update(
        num=fmt_int,
        dec=fmt_dec,
        pct=fmt_pct,
        frac=fmt_frac,
        sortval=sortval,
        short=short,
        dash=dash,
    )
    return env


def render(data, template_path=DEFAULT_TEMPLATE):
    env = make_environment(template_path)
    template = env.get_template(Path(template_path).name)
    return template.render(**build_context(data))


def parse_args(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--data", required=True, help="report_data.json")
    parser.add_argument("--template", default=str(DEFAULT_TEMPLATE), help="report.html.j2")
    parser.add_argument("--out", required=True, help="report.html")
    return parser.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    try:
        data = json.loads(Path(args.data).read_text(encoding="utf-8"))
    except (OSError, ValueError) as err:
        print(f"render_report.py: cannot read {args.data}: {err}", file=sys.stderr)
        return 2
    try:
        html = render(data, args.template)
    except jinja2.TemplateError as err:
        print(f"render_report.py: template error: {err}", file=sys.stderr)
        return 2
    Path(args.out).write_text(html, encoding="utf-8")
    return 0


if __name__ == "__main__":
    sys.exit(main())
