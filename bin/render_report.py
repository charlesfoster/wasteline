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
VISIBLE_ROWS = 100  # rows of a long table shown before "show all"
MAX_TABLE_ROWS = 500  # hard cap on rendered rows of a long table (the TSV has them all)
TS_TOP_GROUPS = 6
TS_MAX_DOTS_DATES = 30  # more dates than this: lines only, no per-point markers
REASON_TITLES = {
    "low_primer_trimming": (
        "Very few reads had primer sequence trimmed (below your --qc_min_primer_trimmed_pct). "
        "Heuristic: it can mean the wrong primer scheme was used."
    ),
    "primer_scheme_mismatch": (
        "Read starts line up much better with another primer scheme than with the declared one."
    ),
}
PRIMER_TRIM_TITLE = (
    "Share of mapped reads from which ivar trimmed primers; informational "
    "(it does not separate right from wrong schemes)"
)
PRIMER_CHECK_TITLES = {
    "ok": "Read starts line up with the declared primer scheme",
    "uninformative": (
        "Read starts do not line up with any primer scheme (e.g. tagmentation libraries): "
        "no scheme check possible"
    ),
}
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
    out["timeseries"] = _as_dict(data.get("timeseries")) or None
    out["growth_rates"] = _as_dict(data.get("growth_rates")) or None
    out["dashboard"] = _as_dict(data.get("dashboard")) or None
    out["dropout_impact"] = _as_dict(data.get("dropout_impact")) or None
    unexplained = data.get("unexplained_mutations")
    lineage_groups = _as_dict(data.get("lineage_groups"))
    out["lineage_groups"] = lineage_groups or None
    out["site_aggregates"] = [
        a for a in _as_list(data.get("site_aggregates")) if isinstance(a, dict)
    ] or None
    out["unexplained_mutations"] = (
        [r for r in unexplained if isinstance(r, dict) and r.get("sample") is not None]
        if isinstance(unexplained, list)
        else None
    )
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


def sci(value, digits=2):
    """p-values and other tiny numbers: 2.39e-42."""
    value = _num(value)
    return "–" if value is None else f"{value:.{digits}e}"


def pval(value):
    """A p-value; the pipeline floors them at 1e-300, which means "at most 1e-300"."""
    value = _num(value)
    if value is None:
        return "–"
    return "<=1e-300" if value <= 1e-300 else f"{value:.2e}"


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


def lineage_view(samples, demix, key, palette=None):
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
    palette = _palette(totals) if palette is None else palette
    legend = [{"name": n, "cls": f"c{i}"} for n, i in palette.items() if n in totals]
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


def _sample_flags(samples):
    return {
        r["sample"]: (r["sample_type"] == "negative", r.get("qc_status") or "NA") for r in samples
    }


def _capped(rows):
    """Split rows into the shown ones (the first VISIBLE_ROWS are visible, the rest start folded)."""
    return {
        "rows": rows[:MAX_TABLE_ROWS],
        "n": len(rows),
        "n_shown": min(len(rows), MAX_TABLE_ROWS),
        "n_visible": min(len(rows), VISIBLE_ROWS),
        "truncated": len(rows) > MAX_TABLE_ROWS,
    }


def qc_extra_view(data):
    """Primer-trimming / primer-scheme columns of the QC table (hidden when no sample has a value)."""
    samples = data["samples"]
    has_trim = any(_num(r.get("primer_trimmed_pct")) is not None for r in samples)
    has_check = any(r.get("primer_scheme_check") for r in samples)
    declared = data["params_summary"].get("primer_scheme") or _as_dict(
        data["provenance"].get("params")
    ).get("primer_scheme")
    for rec in samples:
        check = rec.get("primer_scheme_check")
        best = rec.get("primer_scheme_best")
        enrichment = _num(rec.get("primer_scheme_enrichment"))
        if check == "mismatch":
            title = f"best match: {best or 'another scheme'}"
            if declared:
                title += f" (declared: {declared})"
        else:
            title = PRIMER_CHECK_TITLES.get(check, "")
        if check and enrichment is not None:
            title += f"; declared-scheme read-start enrichment {enrichment:.1f}"
        rec["_pcheck_title"] = title
        rec["_rtip"] = any(r in REASON_TITLES for r in rec["qc_reasons"])
    return {
        "trim": has_trim,
        "check": has_check,
        "reason_titles": REASON_TITLES,
        "trim_title": PRIMER_TRIM_TITLE,
    }


def dropout_view(data):
    impact = data["dropout_impact"]
    if not impact:
        return None
    flags = _sample_flags(data["samples"])
    rows, summary = [], []
    n_affected = 0
    for sample, entry in sorted(impact.items()):
        entry = _as_dict(entry)
        ntc, qc = flags.get(sample, (False, "NA"))
        s = _as_dict(entry.get("summary"))
        amps = [a for a in _as_list(entry.get("amplicons")) if isinstance(a, dict)]
        affecting = [a for a in amps if a.get("affects_detected_lineage")]
        n_affecting = s.get("n_dropped_affecting_detected_lineages")
        no_estimate = bool(amps) and all(
            a.get("n_sites_in_detected_lineages") is None for a in amps
        )
        if affecting:
            n_affected += 1
        n_dropped = s.get("n_dropped_amplicons")
        if n_dropped is None:
            n_dropped = len(amps)
        summary.append(
            {
                "sample": sample,
                "ntc": ntc,
                "qc": qc,
                "n_dropped": n_dropped,
                "total": s.get("informative_sites_total"),
                "lost": s.get("informative_sites_in_dropped"),
                "pct": s.get("pct_informative_sites_lost"),
                "def_total": s.get("defining_sites_total"),
                "def_lost": s.get("defining_sites_lost"),
                "def_pct": s.get("pct_defining_sites_lost"),
                "n_affecting": n_affecting if n_affecting is not None else len(affecting),
                "no_estimate": no_estimate,
                "affects": bool(affecting),
            }
        )
        for a in amps:
            rows.append(
                {
                    "sample": sample,
                    "ntc": ntc,
                    "qc": qc,
                    "amplicon": a.get("amplicon"),
                    "start": a.get("start"),
                    "end": a.get("end"),
                    "mean_depth": a.get("mean_depth"),
                    "n_inf": a.get("n_informative_sites"),
                    "n_def": a.get("n_defining_sites"),
                    "n_diff": a.get("n_differential_sites"),
                    "n_det": a.get("n_sites_in_detected_lineages"),
                    "affects": bool(a.get("affects_detected_lineage")),
                    "no_estimate": a.get("n_sites_in_detected_lineages") is None,
                    "lineages": [
                        l for l in _as_list(a.get("affected_lineages")) if isinstance(l, dict)
                    ],
                }
            )
    rows.sort(key=lambda r: (not r["affects"], r["sample"], r["start"] or 0))
    with_dropout = [x for x in summary if (x["n_dropped"] or 0) > 0]
    # tables written before the "defining sites" refinement lack these columns: keep the old layout
    has_defining = any(r["n_def"] is not None for r in rows) or any(
        x["def_pct"] is not None or x["def_total"] is not None for x in summary
    )
    if has_defining:
        with_dropout.sort(
            key=lambda x: (-(x["def_pct"] if x["def_pct"] is not None else -1), x["sample"])
        )
    else:
        with_dropout.sort(key=lambda x: (-(x["pct"] if x["pct"] is not None else -1), x["sample"]))
    return {
        "has_defining": has_defining,
        "min_affect": _num(data["thresholds"].get("dropout_min_affect_abundance")),
        "table": _capped(rows),
        "summary": with_dropout,
        "n_samples": len(summary),
        "n_without": len(summary) - len(with_dropout),
        "n_affected": n_affected,
        "min_abundance": _num(data["thresholds"].get("dropout_min_abundance")),
    }


def unexplained_view(data):
    rows = data["unexplained_mutations"]
    if rows is None:
        return None
    flags = _sample_flags(data["samples"])
    out = []
    for r in rows:
        ntc, qc = flags.get(r["sample"], (False, "NA"))
        label = f"{r['gene']}:{r['aa_change']}" if r.get("gene") and r.get("aa_change") else None
        out.append(dict(r, ntc=ntc, qc=qc, label=label))
    t = data["thresholds"]
    keys = (
        ("unexplained_min_depth", "depth", "{:g}x"),
        ("unexplained_min_alt_depth", "alt reads", "{:g}"),
        ("unexplained_min_freq", "observed", "{:.0%}"),
        ("unexplained_min_excess", "excess", "{:.0%}"),
        ("unexplained_max_pvalue", "p-value", "{:.2e}"),
    )
    thresholds = []
    for key, name, fmt in keys:
        value = _num(t.get(key))
        if value is not None:
            thresholds.append(
                f"{name} \u2265 {fmt.format(value)}"
                if name != "p-value"
                else f"p-value \u2264 {fmt.format(value)}"
            )
    n_recurrent = len({(r["position"], r["alt"]) for r in out if r.get("recurrent")})
    n_blind = len({(r["position"], r["alt"]) for r in out if r.get("in_dataset") is False})
    return {
        "has_in_dataset": any(r.get("in_dataset") is not None for r in out),
        "n_blind": n_blind,
        "table": _capped(out),
        "n_recurrent": n_recurrent,
        "n_samples": len({r["sample"] for r in out}),
        "thresholds": thresholds,
        "min_recurrence": _num(t.get("unexplained_min_recurrence")),
    }


def growth_view(data):
    growth = data["growth_rates"]
    dash_info = data["dashboard"]
    if not growth and not dash_info:
        return None
    growth = growth or {}
    rows = [r for r in _as_list(growth.get("rows")) if isinstance(r, dict) and r.get("lineage")]
    rows.sort(
        key=lambda r: -(r["advantage_pct"] if _num(r.get("advantage_pct")) is not None else -1e9)
    )
    return {
        "status": growth.get("status"),
        "reason": growth.get("reason"),
        "rows": rows,
        "inputs": _as_dict(growth.get("inputs")),
        "params": _as_dict(growth.get("params")),
        "dashboard": _as_dict(dash_info).get("file"),
    }


def _ts_group_classes(groups, palette):
    used = {palette[g] for g in groups if g in palette}
    free = iter(i for i in range(N_COLORS) if i not in used)
    cls = {}
    for g in groups:
        if g in palette:
            cls[g] = f"c{palette[g]}"
        else:
            nxt = next(free, None)
            cls[g] = f"c{nxt}" if nxt is not None else "co"
    return cls


def timeseries_view(data, palettes):
    """Time-series panels from the freyja_timeseries tables: per (level, weighting) one view with
    one panel per site plus a pooled panel; the page toggles between them."""
    ts = data["timeseries"]
    if not ts:
        return None
    all_dates = sorted({str(d)[:10] for d in _as_list(ts.get("dates"))})
    try:
        ordinals = {d: date.fromisoformat(d).toordinal() for d in all_dates}
    except ValueError:
        return None
    if len(ordinals) < 2:
        return None
    lo, hi = min(ordinals.values()), max(ordinals.values())
    left, right, top, bottom, width, height = 38.0, 10.0, 8.0, 24.0, 420.0, 190.0
    sites = sorted(
        {
            r.get("site")
            for lvl in ("summarized", "lineages")
            for r in _as_list(ts.get(lvl))
            if isinstance(r, dict)
        },
        key=lambda x: (x is None, x or ""),
    )
    has_load = bool(ts.get("has_viral_load"))
    views = []
    levels_present = {}
    for level, label, palette in (
        ("summarized", "summ", palettes["summarized"]),
        ("lineages", "full", palettes["lineages"]),
    ):
        rows = [
            r for r in _as_list(ts.get(level)) if isinstance(r, dict) and r.get("collection_date")
        ]
        levels_present[label] = bool(rows)
        if not rows:
            continue
        totals = Counter()
        for r in rows:
            totals[r["lineage_group"]] += _num(r.get("abundance")) or 0.0
        groups = [g for g, _ in totals.most_common() if g != "Other"][:TS_TOP_GROUPS]
        cls = _ts_group_classes(groups, palette)
        # site -> date -> {group: row}
        cell = {}
        for r in rows:
            cell.setdefault(r.get("site"), {}).setdefault(str(r["collection_date"])[:10], {})[
                r["lineage_group"]
            ] = r
        for weight in ("unw", "w") if has_load else ("unw",):
            key = "abundance_load_weighted" if weight == "w" else "abundance"
            per_site = {}
            for site, by_date in cell.items():
                pts = {}
                for d, grp in by_date.items():
                    n = max((_num(x.get("n_samples")) or 0) for x in grp.values())
                    load = [x for x in grp.values() if _num(x.get("total_viral_load"))]
                    vals = {}
                    ok = True
                    for g in groups:
                        row = grp.get(g)
                        v = _num(row.get(key)) if row else 0.0
                        if v is None:
                            ok = False
                            break
                        vals[g] = v
                    pts[d] = {
                        "n": int(n),
                        "vals": vals if ok else None,
                        "load": max((_num(x.get("mean_viral_load")) or 0) for x in grp.values())
                        if load
                        else None,
                        "total": max((_num(x.get("total_viral_load")) or 0) for x in grp.values())
                        if load
                        else None,
                    }
                per_site[site] = pts
            panels_data = [
                (site, per_site[site])
                for site in sorted(per_site, key=lambda x: (x is None, x or ""))
            ]
            pooled = None
            if len(per_site) > 1:
                pooled = {}
                for d in sorted({d for pts in per_site.values() for d in pts}):
                    parts = [pts[d] for pts in per_site.values() if d in pts and pts[d]["vals"]]
                    if not parts:
                        continue
                    if weight == "w":
                        w = [p["total"] or 0.0 for p in parts]
                    else:
                        w = [float(p["n"] or 1) for p in parts]
                    if sum(w) <= 0:
                        continue
                    pooled[d] = {
                        "n": sum(p["n"] for p in parts),
                        "vals": {
                            g: sum(wi * p["vals"][g] for wi, p in zip(w, parts, strict=True))
                            / sum(w)
                            for g in groups
                        },
                        "load": None,
                        "total": None,
                    }
            vmax = max(
                [0.1]
                + [
                    v
                    for _, pts in panels_data
                    for p in pts.values()
                    if p["vals"]
                    for v in p["vals"].values()
                ]
                + [v for p in (pooled or {}).values() for v in p["vals"].values()]
            )
            ymax = 1.0 if vmax > 0.5 else (math.ceil(vmax * 10) / 10)

            def xpos(day, lo=lo, hi=hi):
                return left + (day - lo) / (hi - lo) * (width - left - right)

            def ypos(value, ymax=ymax):
                return top + (1 - min(value / ymax, 1.0)) * (height - top - bottom)

            def panel(
                site_key, title, pts, groups=groups, cls=cls, weight=weight, xpos=xpos, ypos=ypos
            ):
                lines = []
                dates = sorted(d for d, p in pts.items() if p["vals"])
                for g in groups:
                    dots = []
                    for d in dates:
                        v = pts[d]["vals"][g]
                        load = pts[d]["load"]
                        extra = f"; mean viral load {fmt_int(load)}" if load else ""
                        dots.append(
                            {
                                "x": round(xpos(ordinals[d]), 1),
                                "y": round(ypos(v), 1),
                                "title": f"{g} on {d}: {v * 100:.1f}% (n={pts[d]['n']}{extra})",
                            }
                        )
                    show_dots = len(dates) <= TS_MAX_DOTS_DATES
                    lines.append(
                        {
                            "group": g,
                            "cls": cls[g],
                            "points": " ".join(f"{p['x']},{p['y']}" for p in dots),
                            "dots": dots if show_dots else dots[-1:],
                        }
                    )
                return {"key": site_key, "title": title, "lines": lines}

            panels = [
                panel(site if site is not None else "", site or "(no site)", pts)
                for site, pts in panels_data
            ]
            if pooled is not None:
                panels.insert(0, panel("__all__", "All sites (pooled)", pooled))
            views.append(
                {
                    "level": label,
                    "weight": weight,
                    "id": f"tsv-{label}-{weight}",
                    "panels": panels,
                    "legend": [{"name": g, "cls": cls[g]} for g in groups],
                    "yticks": [
                        {"y": round(ypos(v), 1), "label": f"{v * 100:.0f}%"}
                        for v in (0, ymax / 2, ymax)
                    ],
                }
            )
    if not views:
        return None
    return {
        "views": views,
        "sites": [{"key": s if s is not None else "", "name": s or "(no site)"} for s in sites],
        "has_pooled": len(sites) > 1,
        "has_load": has_load,
        "has_lineages": levels_present.get("full", False),
        "width": int(width),
        "height": int(height),
        "x0": left,
        "x1": width - right,
        "first": all_dates[0],
        "last": all_dates[-1],
        "n_dates": len(ordinals),
        "min_coverage": _num(ts.get("min_coverage")),
        "lineages_total_groups": ts.get("lineages_total_groups"),
    }


def lineage_groups_info(data):
    """One provenance line about the user's lineage groups (the table itself is in the CSV)."""
    lg = data["lineage_groups"]
    if not lg:
        return None
    sha = (lg.get("groups_sha256") or "")[:12]
    return {
        "file": lg.get("file") or "bundled default",
        "sha": sha or None,
        "dataset": lg.get("dataset_version"),
    }


def _merge_palette(base, totals):
    """Keep the colours the per-sample view gave to names present here; give new names free ones."""
    pal = {n: i for n, i in base.items() if n in totals}
    free = iter(i for i in range(N_COLORS) if i not in set(pal.values()))
    for name, _ in totals.most_common():
        if name in pal or name == "Other":
            continue
        nxt = next(free, None)
        if nxt is None:
            break
        pal[name] = nxt
    return pal


def site_aggregates_view(data, palettes):
    entries = data["site_aggregates"]
    if not entries:
        return None
    samples, demix = [], {}
    for a in entries:
        label = str(a.get("row_label"))
        has_solution = bool(_as_list(a.get("lineages")) or _as_list(a.get("groups")))
        samples.append({"sample": label, "sample_type": "sample", "qc_status": None})
        demix[label] = {
            "status": "ok" if has_solution else "missing",
            "reason": "no site contributed a result in this period" if not has_solution else "",
            "groups": a.get("groups"),
            "lineages": a.get("lineages"),
        }
    views = {}
    for key, label in (("groups", "summ"), ("lineages", "full")):
        totals = Counter()
        for a in entries:
            for p in _as_list(a.get(key)):
                if len(p) == 2 and _num(p[1]) is not None:
                    totals[str(p[0])] += p[1]
        pal = _merge_palette(palettes[label], totals)
        view = lineage_view(samples, demix, key, palette=pal)
        views[label] = {
            "legend": view["legend"],
            "rows": {r["sample"]: r for r in view["rows"]},
            "n_ok": view["n_ok"],
        }
    blocks = {}
    for a in entries:
        label = str(a.get("row_label"))
        missing = _as_list(a.get("sites_missing"))
        reasons = _as_dict(a.get("missing_reasons"))
        n_inc = a.get("n_sites_included")
        if n_inc is None:
            n_inc = len(_as_list(a.get("sites_included")))
        basis = a.get("weight_basis")
        frac = _num(a.get("weight_fraction_included"))
        period = {
            "label": label,
            "period": a.get("period"),
            "start": a.get("period_start"),
            "n_inc": n_inc,
            "n_total": n_inc + len(missing),
            "included": ", ".join(str(x) for x in _as_list(a.get("sites_included"))),
            "missing": [
                {"site": str(m), "reason": reasons.get(m) or "no reason recorded"} for m in missing
            ],
            "basis": basis,
            "weak_basis": str(basis or "").startswith(("equal (", "site_weight (")),
            "frac": None if frac is None else frac * 100,
            "n_samples": a.get("n_samples"),
            "rows": {k: views[k]["rows"][label] for k in views},
        }
        blocks.setdefault(a.get("aggregate_name") or label, []).append(period)
    aggregates = []
    for name, periods in blocks.items():
        periods.sort(key=lambda p: (p["start"] or p["period"] or "", p["label"]))
        aggregates.append({"name": name, "periods": periods})
    return {
        "aggregates": aggregates,
        "legend": {k: views[k]["legend"] for k in views},
        "has_full": views["full"]["n_ok"] > 0,
        "n_periods": len(entries),
        "n_incomplete": sum(1 for a in aggregates for p in a["periods"] if p["missing"]),
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
        "qc_extra": qc_extra_view(data),
        "dx": dropout_view(data),
        "ue": unexplained_view(data),
        "ts_new": timeseries_view(
            data, {"summarized": summ["palette"], "lineages": full["palette"]}
        ),
        "growth": growth_view(data),
        "lgi": lineage_groups_info(data),
        "sa": site_aggregates_view(data, {"summ": summ["palette"], "full": full["palette"]}),
        "coverage": coverage_view(data),
        "amp": amplicon_view(data),
        "ts": time_series_view(data, summ["palette"]),
        "moi": moi_view(data),
        "bc": barcode_view(data),
        "controls": data["negative_controls"],
        "consensus": [r for r in data["consensus"] if isinstance(r, dict)],
        "cons_reason": any(r.get("reason") for r in data["consensus"] if isinstance(r, dict)),
        "extras": False,
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
    qc_extra = context["qc_extra"]
    context["extras"] = bool(
        context["dx"]
        or context["ue"] is not None
        or context["ts_new"]
        or context["growth"]
        or context["sa"]
        or qc_extra["trim"]
        or qc_extra["check"]
        or any(r.get("_rtip") for r in samples)
    )
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
        sci=sci,
        pval=pval,
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
