"""Tests for bin/site_aggregate.py: result-level aggregation of per-site Freyja solutions."""

import csv
import subprocess
import sys
from pathlib import Path

import pytest
import site_aggregate as sa

REPO = Path(__file__).resolve().parents[2]
TS = REPO / "tests" / "data_timeseries"
FIX = REPO / "tests" / "fixtures" / "site_aggregate"
SCRIPT = REPO / "bin" / "site_aggregate.py"
SHEET_HEADER = "sample,sample_type,collection_date,site,viral_load"
DEMIX_HEADER = "sample\tstatus\treason\tgenome_coverage_pct\tbarcode_version\tbarcode_version_iso"


# --------------------------------------------------------------------------------------------------
# helpers
# --------------------------------------------------------------------------------------------------
def S(name, site, date, lin, load=None, **kw):
    """A sample spec. `lin` = {lineage: fraction}; `summ` defaults to the same numbers under group names."""
    return {
        "sample": name,
        "site": site,
        "date": date,
        "lin": lin,
        "load": load,
        "type": kw.get("type", "sample"),
        "status": kw.get("status", "ok"),
        "summ": kw.get("summ", {f"G_{k}": v for k, v in lin.items()}),
        "resid": kw.get("resid", 1.0),
        "coverage": kw.get("coverage", 99.0),
        "in_agg": kw.get("in_agg", True),
    }


def write_inputs(tmp_path, samples, sites_text):
    agg = ["\tsummarized\tlineages\tabundances\tresid\tcoverage\tpathogen"]
    sheet = [SHEET_HEADER]
    status = [DEMIX_HEADER]
    for s in samples:
        in_status = s["status"] in ("ok", "low_confidence")
        if s["in_agg"] and in_status or s["in_agg"] and s["type"] == "negative":
            summ = "[" + ", ".join(f"({k!r}, {v!r})" for k, v in s["summ"].items()) + "]"
            names = " ".join(s["lin"])
            vals = " ".join(f"{v:.8f}" for v in s["lin"].values())
            agg.append(
                f"{s['sample']}\t{summ}\t{names}\t{vals}\t{s['resid']}\t{s['coverage']}\tSARS-CoV-2"
            )
        load = "" if s["load"] is None else s["load"]
        sheet.append(f"{s['sample']},{s['type']},{s['date']},{s['site']},{load}")
        status.append(f"{s['sample']}\t{s['status']}\t\t99\tcustom\t")
    paths = {
        "aggregated": tmp_path / "agg.tsv",
        "status": tmp_path / "status.tsv",
        "sheet": tmp_path / "sheet.csv",
        "sites": tmp_path / "sites.csv",
        "out": tmp_path / "out.tsv",
        "meta": tmp_path / "meta.tsv",
    }
    paths["aggregated"].write_text("\n".join(agg) + "\n")
    paths["status"].write_text("\n".join(status) + "\n")
    paths["sheet"].write_text("\n".join(sheet) + "\n")
    paths["sites"].write_text(sites_text)
    return paths


def run(tmp_path, samples, sites_text, *extra, expect=0):
    p = write_inputs(tmp_path, samples, sites_text)
    result = run_paths(p, *extra)
    assert result.returncode == expect, result.stderr
    return load_outputs(p) + (result.stderr,)


def run_paths(p, *extra):
    return subprocess.run(
        [
            sys.executable,
            str(SCRIPT),
            "--aggregated",
            str(p["aggregated"]),
            "--demix-status",
            str(p["status"]),
            "--samplesheet",
            str(p["sheet"]),
            "--sites",
            str(p["sites"]),
            "--out",
            str(p["out"]),
            "--meta-out",
            str(p["meta"]),
            *extra,
        ],
        capture_output=True,
        text=True,
        check=False,
    )


def load_outputs(p):
    """-> ({row_label: {lin, summ, resid, coverage}}, [meta dict])."""
    agg, _ = sa.read_aggregate(p["aggregated_out"] if "aggregated_out" in p else p["out"])
    with open(p["meta"], newline="") as handle:
        meta = list(csv.DictReader(handle, delimiter="\t"))
    return agg, meta


SITES_XY = "aggregate_name,site_name,weight\nAgg,X,2\nAgg,Y,1\n"
MON = "2026-01-05"  # a Monday (ISO week 2026-W02)


def approx(x, tol=1e-7):
    return pytest.approx(x, abs=tol)


def two_site_samples(load_x2=100):
    return [
        S("x1", "X", MON, {"BA.2": 0.8, "JN.1": 0.2}, 300),
        S("x2", "X", "2026-01-07", {"BA.2": 0.4, "XBB": 0.6}, load_x2),
        S("y1", "Y", "2026-01-06", {"BA.2": 0.5, "XBB": 0.5}, 600),
    ]


def test_replicates_then_sites_by_viral_load(tmp_path):
    agg, meta, _ = run(tmp_path, two_site_samples(), SITES_XY)
    row = agg["Agg_2026-W02"]
    # site X = (300*s1 + 100*s2)/400: BA.2 0.7, JN.1 0.15, XBB 0.15; mean load 200
    # site Y load 600 -> site weights 1:3 -> BA.2 (0.7 + 1.5)/4, JN.1 0.15/4, XBB (0.15 + 1.5)/4
    assert row["lin"]["BA.2"] == approx(0.55)
    assert row["lin"]["JN.1"] == approx(0.0375)
    assert row["lin"]["XBB"] == approx(0.4125)
    assert sum(row["lin"].values()) == approx(1.0)
    m = meta[0]
    assert m["weight_basis"] == "viral_load"
    assert m["n_samples"] == "3"
    assert m["n_sites_included"] == "2"
    assert m["sites_included"] == "X;Y"
    assert m["weight_fraction_included"] == "1.00000000"


def test_missing_sites_and_reasons(tmp_path):
    samples = [
        S("x1", "X", MON, {"A": 1.0}),
        S("y1", "Y", MON, {"B": 1.0}, status="failed"),
        S("z1", "Z", "2026-02-02", {"C": 1.0}),  # a different week
        S("w1", "W", "", {"C": 1.0}),  # no date
        S("n1", "X", MON, {"Z": 1.0}, type="negative"),  # controls never contribute
    ]
    sites = "aggregate_name,site_name\nAgg,X\nAgg,Y\nAgg,Z\nAgg,W\nAgg,V\n"
    agg, meta, err = run(tmp_path, samples, sites)
    first = meta[0]
    assert first["period"] == "2026-W02"
    assert first["sites_included"] == "X"
    assert first["sites_missing"] == "Y;Z;W;V"
    assert first["missing_reasons"] == "not_ok:failed;no_sample;no_date;no_sample"
    assert first["n_sites_included"] == "1"
    assert agg["Agg_2026-W02"]["lin"] == {"A": approx(1.0)}
    # V is in no sample's site column: one warning, once
    assert err.count("site 'V' in the sites file matches no sample") == 1
    assert "site 'W'" not in err


def test_cli_error_exit_code_and_message(tmp_path):
    p = write_inputs(
        tmp_path, [S("x1", "X", MON, {"A": 1.0})], "aggregate_name,site_name\nMetro,\n"
    )
    result = run_paths(p)
    assert result.returncode == 1
    assert "line 2: site_name is blank" in result.stderr
    assert "Traceback" not in result.stderr
    assert not p["out"].exists()


def run_fixture(tmp_path, sites=FIX / "sites.csv", *extra):
    p = {
        "aggregated": TS / "freyja_aggregated.tsv",
        "status": FIX / "demix_status.tsv",
        "sheet": TS / "samplesheet.csv",
        "sites": sites,
        "out": tmp_path / "out.tsv",
        "meta": tmp_path / "meta.tsv",
    }
    result = run_paths(p, *extra)
    assert result.returncode == 0, result.stderr
    agg, meta = load_outputs(p)
    return agg, {m["row_label"]: m for m in meta}, result.stderr


def test_fixture_hand_computed_week(tmp_path):
    agg, meta, _err = run_fixture(tmp_path)
    # 2026-01-12 = ISO week 2026-W03. XBB.1.5 per sample (freyja_aggregated.tsv) and load:
    #   SiteA  0.16710542 (851825)  0.23093447 (301960)   -> site load 576892.5
    #   SiteB  0.08106117 (1169960) 0.11608544 (938723)   -> site load 1054341.5
    #   SiteC  0.05187679 (110143)  0.04421386 (964206)   -> site load 537174.5
    a = (0.16710542 * 851825 + 0.23093447 * 301960) / (851825 + 301960)
    b = (0.08106117 * 1169960 + 0.11608544 * 938723) / (1169960 + 938723)
    c = (0.05187679 * 110143 + 0.04421386 * 964206) / (110143 + 964206)
    la, lb, lc = (851825 + 301960) / 2, (1169960 + 938723) / 2, (110143 + 964206) / 2
    expected = (a * la + b * lb + c * lc) / (la + lb + lc)
    row = agg["Metro_2026-W03"]
    assert a == approx(0.183810, 1e-6)
    assert row["lin"]["XBB.1.5"] == approx(expected)
    assert meta["Metro_2026-W03"]["weight_basis"] == "viral_load"
    assert meta["Metro_2026-W03"]["n_samples"] == "6"
    # coverage and resid with the same weights
    cov_a = (94.4347 * 851825 + 96.7247 * 301960) / (851825 + 301960)
    cov_b = (98.3818 * 1169960 + 96.8068 * 938723) / (1169960 + 938723)
    cov_c = (97.803 * 110143 + 95.059 * 964206) / (110143 + 964206)
    cov = (cov_a * la + cov_b * lb + cov_c * lc) / (la + lb + lc)
    assert row["coverage"] == approx(cov)
