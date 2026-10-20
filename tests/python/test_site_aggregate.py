"""Tests for bin/site_aggregate.py: result-level aggregation of per-site Freyja solutions."""

import csv
import datetime as dt
import random
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
SITES_XY_NOW = "aggregate_name,site_name\nAgg,X\nAgg,Y\n"
MON = "2026-01-05"  # a Monday (ISO week 2026-W02)


def approx(x, tol=1e-7):
    return pytest.approx(x, abs=tol)


# --------------------------------------------------------------------------------------------------
# hand-checkable algorithm cases
# --------------------------------------------------------------------------------------------------
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


def test_equal_weights_but_replicates_still_combined(tmp_path):
    agg, meta, _ = run(tmp_path, two_site_samples(), SITES_XY, "--weights", "equal")
    row = agg["Agg_2026-W02"]
    # replicates by load (X = 0.7/0.15/0.15), sites equally
    assert row["lin"]["BA.2"] == approx((0.7 + 0.5) / 2)
    assert row["lin"]["JN.1"] == approx(0.075)
    assert row["lin"]["XBB"] == approx(0.325)
    assert meta[0]["weight_basis"] == "equal"


def test_site_with_two_samples_does_not_count_double(tmp_path):
    samples = [
        S("x1", "X", MON, {"A": 1.0}),
        S("x2", "X", MON, {"A": 1.0}),
        S("x3", "X", MON, {"A": 1.0}),
        S("y1", "Y", MON, {"B": 1.0}),
    ]
    agg, meta, _ = run(tmp_path, samples, SITES_XY_NOW)
    assert agg["Agg_2026-W02"]["lin"] == {"A": approx(0.5), "B": approx(0.5)}
    assert meta[0]["n_samples"] == "4"


def test_replicates_equal_when_one_has_no_load(tmp_path):
    samples = [
        S("x1", "X", MON, {"A": 1.0}, 1000),
        S("x2", "X", MON, {"B": 1.0}, None),
        S("y1", "Y", MON, {"B": 1.0}, 10),
    ]
    agg, meta, _ = run(tmp_path, samples, SITES_XY, "--weights", "equal")
    # X: no load for every replicate -> equal 0.5/0.5; equal sites
    assert agg["Agg_2026-W02"]["lin"] == {"A": approx(0.25), "B": approx(0.75)}
    # auto: X has no complete load -> site weights 2:1 from the file
    agg, meta, _ = run(tmp_path, samples, SITES_XY)
    assert meta[0]["weight_basis"] == "site_weight (viral_load missing)"
    assert agg["Agg_2026-W02"]["lin"]["A"] == approx(2 * 0.5 / 3)
    assert agg["Agg_2026-W02"]["lin"]["B"] == approx((2 * 0.5 + 1) / 3)


def test_site_weight_basis(tmp_path):
    agg, meta, _ = run(tmp_path, two_site_samples(), SITES_XY, "--weights", "site_weight")
    row = agg["Agg_2026-W02"]
    # weights 2:1 (X:Y)
    assert row["lin"]["BA.2"] == approx((2 * 0.7 + 0.5) / 3)
    assert row["lin"]["XBB"] == approx((2 * 0.15 + 0.5) / 3)
    assert meta[0]["weight_basis"] == "site_weight"


def test_explicit_choices_fall_back_softly(tmp_path):
    agg, meta, _ = run(tmp_path, two_site_samples(), SITES_XY_NOW, "--weights", "site_weight")
    assert meta[0]["weight_basis"] == "equal (site_weight missing)"
    assert agg["Agg_2026-W02"]["lin"]["BA.2"] == approx(0.6)
    samples = two_site_samples()
    samples[2]["load"] = None  # Y without load
    agg, meta, _ = run(tmp_path, samples, SITES_XY, "--weights", "viral_load")
    assert meta[0]["weight_basis"] == "equal (viral_load missing)"
    assert agg["Agg_2026-W02"]["lin"]["BA.2"] == approx(0.6)


def test_auto_prefers_load_then_weight_then_equal(tmp_path):
    samples = two_site_samples()
    samples[2]["load"] = None
    _, meta, _ = run(tmp_path, samples, SITES_XY)
    assert meta[0]["weight_basis"] == "site_weight (viral_load missing)"
    _, meta, _ = run(tmp_path, samples, SITES_XY_NOW)
    assert meta[0]["weight_basis"] == "equal (viral_load missing)"
    for s in samples:
        s["load"] = None
    _, meta, _ = run(tmp_path, samples, SITES_XY_NOW)
    assert meta[0]["weight_basis"] == "equal"
    # a weight for only one of the two sites is not enough
    _, meta, _ = run(tmp_path, samples, "aggregate_name,site_name,weight\nAgg,X,5\nAgg,Y,\n")
    assert meta[0]["weight_basis"] == "equal (site_weight missing)"


def test_basis_can_differ_between_periods(tmp_path):
    samples = two_site_samples() + [
        S("x3", "X", "2026-01-12", {"A": 1.0}, 100),
        S("y3", "Y", "2026-01-12", {"A": 1.0}, None),
    ]
    _, meta, _ = run(tmp_path, samples, SITES_XY)
    assert [m["weight_basis"] for m in meta] == ["viral_load", "site_weight (viral_load missing)"]


def test_summarized_resid_and_coverage_use_the_same_weights(tmp_path):
    samples = [
        S(
            "x1",
            "X",
            MON,
            {"BA.2": 1.0},
            100,
            summ={"Omicron": 0.9, "Delta": 0.1},
            resid=2.0,
            coverage=90,
        ),
        S(
            "y1",
            "Y",
            MON,
            {"BA.2": 0.5, "AY.4": 0.5},
            300,
            summ={"Omicron": 0.5, "Delta": 0.5},
            resid=6.0,
            coverage=98,
        ),
    ]
    agg, _, _ = run(tmp_path, samples, SITES_XY_NOW)
    row = agg["Agg_2026-W02"]
    assert row["summ"]["Omicron"] == approx((100 * 0.9 + 300 * 0.5) / 400)
    assert row["summ"]["Delta"] == approx((100 * 0.1 + 300 * 0.5) / 400)
    assert row["resid"] == approx((100 * 2 + 300 * 6) / 400)
    assert row["coverage"] == approx((100 * 90 + 300 * 98) / 400)
    # the same with equal weights
    agg, _, _ = run(tmp_path, samples, SITES_XY_NOW, "--weights", "equal")
    assert agg["Agg_2026-W02"]["resid"] == approx(4.0)
    assert agg["Agg_2026-W02"]["coverage"] == approx(94.0)


def test_output_is_not_renormalised(tmp_path):
    # inputs that sum to 0.9 and 0.8: the aggregate keeps the weighted sum (0.85), not 1
    samples = [S("x1", "X", MON, {"A": 0.9}), S("y1", "Y", MON, {"A": 0.5, "B": 0.3})]
    agg, _, _ = run(tmp_path, samples, SITES_XY_NOW)
    assert sum(agg["Agg_2026-W02"]["lin"].values()) == approx(0.85)


def test_row_format_matches_freyja_aggregate(tmp_path):
    samples = [
        S(
            "x1",
            "X",
            MON,
            {"BA.2": 0.75, "JN.1": 0.25},
            summ={"Omicron": 0.75, "BA.2.86* (BA.2.86.X)": 0.25},
        )
    ]
    p = write_inputs(tmp_path, samples, "aggregate_name,site_name\nMetro_Sites,X\n")
    assert run_paths(p).returncode == 0
    lines = p["out"].read_text().splitlines()
    assert lines[0] == "\tsummarized\tlineages\tabundances\tresid\tcoverage"
    cells = lines[1].split("\t")
    assert cells[0] == "Metro_Sites_2026-W02"
    assert cells[1] == "[('Omicron', 0.75), ('BA.2.86* (BA.2.86.X)', 0.25)]"
    assert cells[2] == "BA.2 JN.1"
    assert cells[3] == "0.75000000 0.25000000"
    assert cells[4] == "1.00000000"
    assert cells[5] == "99.00000000"
    # decimals in the meta table
    meta = p["meta"].read_text().splitlines()
    assert meta[0].split("\t") == sa.META_COLUMNS
    assert meta[1].split("\t")[0] == "Metro_Sites_2026-W02"


def test_lineages_sorted_by_abundance(tmp_path):
    samples = [S("x1", "X", MON, {"A": 0.1, "B": 0.6, "C": 0.3})]
    p = write_inputs(tmp_path, samples, "aggregate_name,site_name\nAgg,X\n")
    assert run_paths(p).returncode == 0
    assert p["out"].read_text().splitlines()[1].split("\t")[2] == "B C A"


def test_input_parser_forms_numpy_and_repr(tmp_path):
    p = write_inputs(tmp_path, [S("x1", "X", MON, {"A": 1.0})], "aggregate_name,site_name\nAgg,X\n")
    p["aggregated"].write_text(
        "\tsummarized\tlineages\tabundances\tresid\tcoverage\n"
        "x1\t[('Omicron', np.float64(0.25)), ('Delta', np.float64(0.75))]\t['A.1' 'B.2']\t[0.25 0.75]\t1.5\t99\n"
    )
    assert run_paths(p).returncode == 0
    agg, _ = sa.read_aggregate(p["out"])
    assert agg["Agg_2026-W02"]["lin"] == {"A.1": approx(0.25), "B.2": approx(0.75)}
    assert agg["Agg_2026-W02"]["summ"] == {"Omicron": approx(0.25), "Delta": approx(0.75)}


# --------------------------------------------------------------------------------------------------
# periods
# --------------------------------------------------------------------------------------------------
@pytest.mark.parametrize(
    ("date", "label", "start"),
    [
        ("2026-01-05", "2026-W02", "2026-01-05"),  # Monday
        (
            "2026-01-11",
            "2026-W02",
            "2026-01-05",
        ),  # Sunday belongs to the week that started on Monday
        ("2026-01-12", "2026-W03", "2026-01-12"),  # next Monday starts a new week
        ("2026-01-04", "2026-W01", "2025-12-29"),  # Sunday; W01 starts in the previous year
        ("2025-12-29", "2026-W01", "2025-12-29"),  # Monday in December belongs to 2026-W01
        ("2025-12-28", "2025-W52", "2025-12-22"),
        ("2027-01-01", "2026-W53", "2026-12-28"),  # 2026 has 53 ISO weeks
        ("2024-12-30", "2025-W01", "2024-12-30"),
        ("2024-02-29", "2024-W09", "2024-02-26"),  # leap day
    ],
)
def test_week_labels_and_starts(date, label, start):
    got_label, got_start = sa.period_of(dt.date.fromisoformat(date), "week")
    assert got_label == label
    assert got_start.isoformat() == start
    assert got_start.weekday() == 0


def test_sunday_and_monday_fall_in_different_weeks(tmp_path):
    samples = [
        S("x1", "X", "2026-01-11", {"A": 1.0}),  # Sunday W02
        S("x2", "X", "2026-01-12", {"B": 1.0}),  # Monday W03
    ]
    agg, meta, _ = run(tmp_path, samples, "aggregate_name,site_name\nAgg,X\n")
    assert [m["period"] for m in meta] == ["2026-W02", "2026-W03"]
    assert [m["period_start"] for m in meta] == ["2026-01-05", "2026-01-12"]
    assert agg["Agg_2026-W02"]["lin"] == {"A": approx(1.0)}
    assert agg["Agg_2026-W03"]["lin"] == {"B": approx(1.0)}


def test_year_boundary_week_groups_december_and_january(tmp_path):
    samples = [
        S("x1", "X", "2025-12-30", {"A": 1.0}),
        S("x2", "X", "2026-01-02", {"B": 1.0}),
        S("x3", "X", "2026-01-05", {"C": 1.0}),
    ]
    agg, meta, _ = run(tmp_path, samples, "aggregate_name,site_name\nAgg,X\n")
    assert [m["period"] for m in meta] == ["2026-W01", "2026-W02"]
    assert agg["Agg_2026-W01"]["lin"] == {"A": approx(0.5), "B": approx(0.5)}


def test_day_period(tmp_path):
    samples = [
        S("x1", "X", "2026-01-05", {"A": 1.0}),
        S("x2", "X", "2026-01-06", {"B": 1.0}),
        S("x3", "X", "2026-01-06", {"C": 1.0}),
    ]
    agg, meta, _ = run(tmp_path, samples, "aggregate_name,site_name\nAgg,X\n", "--period", "day")
    assert [m["row_label"] for m in meta] == ["Agg_2026-01-05", "Agg_2026-01-06"]
    assert [m["period_start"] for m in meta] == ["2026-01-05", "2026-01-06"]
    assert agg["Agg_2026-01-06"]["lin"] == {"B": approx(0.5), "C": approx(0.5)}


def test_time_part_in_date_is_tolerated(tmp_path):
    _agg, meta, _ = run(
        tmp_path,
        [S("x1", "X", "2026-01-05T08:30:00", {"A": 1.0})],
        "aggregate_name,site_name\nAgg,X\n",
    )
    assert meta[0]["period_start"] == "2026-01-05"


# --------------------------------------------------------------------------------------------------
# membership, matching, exclusions
# --------------------------------------------------------------------------------------------------
def test_site_matching_ignores_case_and_spaces(tmp_path):
    samples = [S("x1", "Site A", MON, {"A": 1.0}), S("y1", "SITE B", MON, {"B": 1.0})]
    agg, meta, err = run(tmp_path, samples, "aggregate_name,site_name\nAgg,  site a \nAgg,site b\n")
    assert meta[0]["n_sites_included"] == "2"
    assert agg["Agg_2026-W02"]["lin"] == {"A": approx(0.5), "B": approx(0.5)}
    assert "matches no sample" not in err


def test_missing_sites_and_reasons(tmp_path):
    samples = [
        S("x1", "X", MON, {"A": 1.0}),
        S("y1", "Y", MON, {"B": 1.0}, status="failed"),
        S("z1", "Z", "2026-02-02", {"C": 1.0}),  # a different week
        S("w1", "W", "", {"C": 1.0}),  # no date
    ]
    sites = "aggregate_name,site_name\nAgg,X\nAgg,Y\nAgg,Z\nAgg,W\nAgg,V\n"
    _, meta, err = run(tmp_path, samples, sites)
    first = meta[0]
    assert first["period"] == "2026-W02"
    assert first["sites_included"] == "X"
    assert first["sites_missing"] == "Y;Z;W;V"
    assert first["missing_reasons"] == "not_ok:failed;no_sample;no_date;no_sample"
    assert first["n_sites_included"] == "1"
    # V is in no sample's site column: one warning, once
    assert err.count("site 'V' in the sites file matches no sample") == 1
    assert "site 'W'" not in err


def test_non_ok_and_negatives_never_contribute(tmp_path):
    samples = [
        S("x1", "X", MON, {"A": 1.0}),
        S("x2", "X", MON, {"B": 1.0}, status="skipped_low_coverage"),
        S("x3", "X", MON, {"C": 1.0}, type="negative"),
        S("x4", "X", MON, {"D": 1.0}, status="missing"),
        S("x5", "X", MON, {"E": 1.0}, status="failed"),
    ]
    agg, meta, err = run(tmp_path, samples, "aggregate_name,site_name\nAgg,X\n")
    assert agg["Agg_2026-W02"]["lin"] == {"A": approx(1.0)}
    assert meta[0]["n_samples"] == "1"
    assert "negative_controls=1" in err and "not_ok=3" in err


def test_low_confidence_contributes_only_with_the_flag(tmp_path):
    samples = [
        S("x1", "X", MON, {"A": 1.0}),
        S("y1", "Y", MON, {"B": 1.0}, status="low_confidence"),
    ]
    sites = "aggregate_name,site_name\nAgg,X\nAgg,Y\n"
    agg, meta, _ = run(tmp_path, samples, sites)
    assert agg["Agg_2026-W02"]["lin"] == {"A": approx(1.0)}
    assert meta[0]["sites_missing"] == "Y" and meta[0]["missing_reasons"] == "not_ok:low_confidence"
    agg, meta, _ = run(tmp_path, samples, sites, "--include-low-confidence")
    assert agg["Agg_2026-W02"]["lin"] == {"A": approx(0.5), "B": approx(0.5)}
    assert meta[0]["sites_missing"] == ""


def test_ok_sample_missing_from_aggregate_and_status_is_not_used(tmp_path):
    samples = [S("x1", "X", MON, {"A": 1.0}), S("x2", "X", MON, {"B": 1.0}, in_agg=False)]
    agg, meta, _ = run(tmp_path, samples, "aggregate_name,site_name\nAgg,X\n")
    assert agg["Agg_2026-W02"]["lin"] == {"A": approx(1.0)}
    # sample with no row in demix_status at all
    p = write_inputs(
        tmp_path,
        [S("x1", "X", MON, {"A": 1.0}), S("x2", "Y", MON, {"B": 1.0})],
        "aggregate_name,site_name\nAgg,X\nAgg,Y\n",
    )
    p["status"].write_text(DEMIX_HEADER + "\nx1\tok\t\t99\tcustom\t\n")
    assert run_paths(p).returncode == 0
    _, meta = load_outputs(p)
    assert meta[0]["sites_missing"] == "Y" and meta[0]["missing_reasons"] == "not_ok:no_status"


def test_no_site_and_no_date_are_counted_and_skipped(tmp_path):
    samples = [
        S("x1", "X", MON, {"A": 1.0}),
        S("n1", "", MON, {"B": 1.0}),
        S("d1", "X", "", {"C": 1.0}),
        S("d2", "X", "05/01/2026", {"C": 1.0}),
    ]
    agg, _meta, err = run(tmp_path, samples, "aggregate_name,site_name\nAgg,X\n")
    assert agg["Agg_2026-W02"]["lin"] == {"A": approx(1.0)}
    assert "no_site=1" in err and "no_date=2" in err


def test_sites_in_several_aggregates(tmp_path):
    samples = [
        S("x1", "X", MON, {"A": 1.0}),
        S("y1", "Y", MON, {"B": 1.0}),
        S("z1", "Z", MON, {"C": 1.0}),
    ]
    sites = "aggregate_name,site_name\nAll,X\nAll,Y\nAll,Z\nXY,X\nXY,Y\nYZ,Y\nYZ,Z\n"
    agg, meta, _ = run(tmp_path, samples, sites, "--weights", "equal")
    assert [m["aggregate_name"] for m in meta] == ["All", "XY", "YZ"]
    assert agg["All_2026-W02"]["lin"] == {
        "A": approx(1 / 3),
        "B": approx(1 / 3),
        "C": approx(1 / 3),
    }
    assert agg["XY_2026-W02"]["lin"] == {"A": approx(0.5), "B": approx(0.5)}
    assert agg["YZ_2026-W02"]["lin"] == {"B": approx(0.5), "C": approx(0.5)}


def test_aggregate_name_case_variants_are_one_aggregate(tmp_path):
    samples = [S("x1", "X", MON, {"A": 1.0}), S("y1", "Y", MON, {"B": 1.0})]
    _, meta, err = run(tmp_path, samples, "aggregate_name,site_name\nMetro,X\nmetro,Y\n")
    assert (
        len(meta) == 1
        and meta[0]["aggregate_name"] == "Metro"
        and meta[0]["n_sites_included"] == "2"
    )
    assert "spelled differently" in err


def test_min_sites(tmp_path):
    samples = [
        S("x1", "X", MON, {"A": 1.0}),
        S("y1", "Y", MON, {"B": 1.0}),
        S("x2", "X", "2026-01-12", {"A": 1.0}),
    ]
    _, meta, err = run(tmp_path, samples, SITES_XY_NOW, "--min-sites", "2")
    assert [m["period"] for m in meta] == ["2026-W02"]
    assert "periods omitted by --min-sites 2: 1" in err
    _, meta, _ = run(tmp_path, samples, SITES_XY_NOW)
    assert [m["period"] for m in meta] == ["2026-W02", "2026-W03"]
    assert meta[1]["sites_missing"] == "Y" and meta[1]["missing_reasons"] == "no_sample"


def test_weight_fraction_included(tmp_path):
    samples = [
        S("x1", "X", MON, {"A": 1.0}, 300),
        S("y1", "Y", MON, {"A": 1.0}, 100, status="failed"),  # load known but cannot contribute
        S("z1", "Z", MON, {"A": 1.0}, status="failed"),
    ]
    sites = "aggregate_name,site_name,weight\nAgg,X,30\nAgg,Y,10\nAgg,Z,60\nAgg,W,\n"
    # viral_load basis: only X contributes; Y's load is known (100), Z has no load at all -> NA
    _, meta, _ = run(tmp_path, samples, sites)
    assert meta[0]["weight_basis"] == "viral_load"
    assert meta[0]["weight_fraction_included"] == "NA"
    # viral_load basis where every missing site has a load from a non-ok sample
    samples = [
        S("x1", "X", MON, {"A": 1.0}, 300),
        S("y1", "Y", MON, {"A": 1.0}, 100, status="failed"),
    ]
    _, meta, _ = run(tmp_path, samples, "aggregate_name,site_name\nAgg,X\nAgg,Y\n")
    assert meta[0]["weight_fraction_included"] == "0.75000000"
    # site_weight basis
    samples = [
        S("x1", "X", MON, {"A": 1.0}),
        S("y1", "Y", MON, {"A": 1.0}),
        S("z1", "Z", MON, {"A": 1.0}, status="failed"),
    ]
    _, meta, _ = run(
        tmp_path, samples, "aggregate_name,site_name,weight\nAgg,X,30\nAgg,Y,10\nAgg,Z,60\n"
    )
    assert meta[0]["weight_basis"] == "site_weight"
    assert meta[0]["weight_fraction_included"] == "0.40000000"
    # a missing site without a weight -> NA
    _, meta, _ = run(
        tmp_path, samples, "aggregate_name,site_name,weight\nAgg,X,30\nAgg,Y,10\nAgg,Z,60\nAgg,W,\n"
    )
    assert meta[0]["weight_fraction_included"] == "NA"
    # equal basis: share of member sites
    _, meta, _ = run(tmp_path, samples, "aggregate_name,site_name\nAgg,X\nAgg,Y\nAgg,Z\nAgg,W\n")
    assert meta[0]["weight_basis"] == "equal"
    assert meta[0]["weight_fraction_included"] == "0.50000000"


def test_sample_ids_are_never_rewritten(tmp_path):
    samples = [S("Site_1_week_2", "Site_1", MON, {"A": 1.0})]
    agg, _meta, _ = run(tmp_path, samples, "aggregate_name,site_name\nAgg_x,Site_1\n")
    assert "Agg_x_2026-W02" in agg


# --------------------------------------------------------------------------------------------------
# properties
# --------------------------------------------------------------------------------------------------
@pytest.mark.parametrize("seed", range(12))
def test_no_new_lineage_and_convex_combination(tmp_path, seed):
    rng = random.Random(seed)
    pool = [f"L{i}" for i in range(10)]
    sites = ["S1", "S2", "S3", "S4"]
    samples = []
    for i in range(rng.randint(8, 25)):
        chosen = rng.sample(pool, rng.randint(1, 5))
        raw = [rng.random() + 0.01 for _ in chosen]
        lin = {k: v / sum(raw) for k, v in zip(chosen, raw)}
        samples.append(
            S(
                f"s{i}",
                rng.choice(sites),
                (dt.date(2026, 1, 5) + dt.timedelta(days=rng.randint(0, 20))).isoformat(),
                lin,
                rng.choice([None, rng.randint(1, 10**6)]),
                status=rng.choice(["ok", "ok", "ok", "failed"]),
                type=rng.choice(["sample"] * 5 + ["negative"]),
            )
        )
    sites_text = "aggregate_name,site_name,weight\n" + "".join(
        f"A{a},{s},{rng.randint(1, 100)}\n"
        for a, group in enumerate([sites[:2], sites[1:], sites])
        for s in group
    )
    for weights in ("auto", "equal", "viral_load", "site_weight"):
        agg, meta, _ = run(tmp_path, samples, sites_text, "--weights", weights)
        aggregates = sa.read_sites(tmp_path / "sites.csv")
        members = {a["name"]: {s["key"] for s in a["sites"]} for a in aggregates}
        for m in meta:
            site_set = members[m["aggregate_name"]]
            used = [
                s
                for s in samples
                if s["status"] == "ok"
                and s["type"] != "negative"
                and sa.site_key(s["site"]) in site_set
                and sa.period_of(dt.date.fromisoformat(s["date"]), "week")[0] == m["period"]
            ]
            assert len(used) == int(m["n_samples"])
            union = {name for s in used for name in s["lin"]}
            row = agg[m["row_label"]]
            assert set(row["lin"]) <= union  # never a lineage no member site reported
            for name, value in row["lin"].items():
                values = [s["lin"].get(name, 0.0) for s in used]
                assert min(values) - 1e-7 <= value <= max(values) + 1e-7
            assert sum(row["lin"].values()) == approx(1.0, 1e-6)
            assert set(row["summ"]) <= {name for s in used for name in s["summ"]}


# --------------------------------------------------------------------------------------------------
# sites file
# --------------------------------------------------------------------------------------------------
def sites_file(tmp_path, text, name="sites.csv", encoding="utf-8"):
    path = tmp_path / name
    path.write_bytes(text.encode(encoding))
    return path


def test_sites_file_bom_and_windows_line_endings(tmp_path):
    path = sites_file(
        tmp_path, "﻿aggregate_name,site_name,weight\r\nMetro,Site1,100\r\nMetro,Site2,\r\n"
    )
    aggs = sa.read_sites(path)
    assert [(a["name"], [s["name"] for s in a["sites"]]) for a in aggs] == [
        ("Metro", ["Site1", "Site2"])
    ]
    assert aggs[0]["sites"][0]["weight"] == 100.0 and aggs[0]["sites"][1]["weight"] is None


def test_sites_file_excel_quirks(tmp_path):
    # semicolon-separated (Excel in many locales), header with spaces/capitals, trailing empty rows
    path = sites_file(
        tmp_path, "Aggregate Name;Site Name;Weight\nMetro;Site1;5\nMetro;Site2;\n;;\n;;\n"
    )
    aggs = sa.read_sites(path)
    assert len(aggs[0]["sites"]) == 2
    # trailing commas from Excel, tab separated
    path = sites_file(tmp_path, "aggregate_name\tsite_name\nMetro\tSite1\n")
    assert sa.read_sites(path)[0]["sites"][0]["name"] == "Site1"
    path = sites_file(tmp_path, "aggregate_name,site_name,,\nMetro,Site1,,\n")
    assert sa.read_sites(path)[0]["sites"][0]["name"] == "Site1"


def test_sites_file_windows_1252_fallback_warns(tmp_path):
    path = sites_file(tmp_path, "aggregate_name,site_name\nMetro,Café\n", encoding="cp1252")
    warnings = []
    aggs = sa.read_sites(path, warnings)
    assert aggs[0]["sites"][0]["name"] == "Café"
    assert any("Windows-1252" in w for w in warnings)


@pytest.mark.parametrize(
    ("text", "needle"),
    [
        ("", "empty"),
        ("\n\n  \n", "empty"),
        ("aggregate_name,site_name,weight\n", "no rows"),
        ("aggregate_name,site_name\n,SiteA\n", "line 2: aggregate_name is blank"),
        ("aggregate_name,site_name\nMetro,\n", "line 2: site_name is blank"),
        ("aggregate_name,site_name\nMetro,A\nMetro,B\n\nMetro,\n", "line 5: site_name is blank"),
        (
            "aggregate_name,site_names\nMetro,A\n",
            "unknown column 'site_names' (did you mean 'site_name'?)",
        ),
        ("aggregate_name,site_name,population\nMetro,A,5\n", "line 1: unknown column 'population'"),
        ("aggregate,site\nMetro,A\n", "unknown column"),
        ("aggregate_name,weight\nMetro,5\n", "no 'site_name' column"),
        (
            "aggregate_name,site_name\nMetro,A\nmetro, a \n",
            "line 3: site 'a' is listed twice for aggregate 'metro' (first on line 2)",
        ),
        (
            "aggregate_name,site_name,weight\nMetro,A,abc\n",
            "line 2: weight 'abc' must be a positive number",
        ),
        (
            "aggregate_name,site_name,weight\nMetro,A,0\n",
            "line 2: weight '0' must be a positive number",
        ),
        ("aggregate_name,site_name,weight\nMetro,A,-3\n", "must be a positive number"),
        ("aggregate_name,site_name,weight\nMetro,A,1,000\n", "more cells"),
    ],
)
def test_sites_file_friendly_errors(tmp_path, text, needle):
    path = sites_file(tmp_path, text)
    with pytest.raises(sa.InputError) as err:
        sa.read_sites(path)
    assert needle in str(err.value)
    assert str(path) in str(err.value)


def test_cli_error_exit_code_and_message(tmp_path):
    p = write_inputs(
        tmp_path, [S("x1", "X", MON, {"A": 1.0})], "aggregate_name,site_name\nMetro,\n"
    )
    result = run_paths(p)
    assert result.returncode == 1
    assert "line 2: site_name is blank" in result.stderr
    assert "Traceback" not in result.stderr
    assert not p["out"].exists()


# --------------------------------------------------------------------------------------------------
# empty inputs and bad inputs
# --------------------------------------------------------------------------------------------------
def test_empty_inputs_write_headers_only(tmp_path):
    p = write_inputs(tmp_path, [], "aggregate_name,site_name\nAgg,X\n")
    p["aggregated"].write_text("")  # a zero-byte aggregate
    p["status"].write_text("")
    result = run_paths(p)
    assert result.returncode == 0, result.stderr
    assert p["out"].read_text() == "\tsummarized\tlineages\tabundances\tresid\tcoverage\n"
    assert p["meta"].read_text() == "\t".join(sa.META_COLUMNS) + "\n"
    assert "site 'X' in the sites file matches no sample" in result.stderr
    assert "aggregate 'Agg' has no usable period" in result.stderr


def test_nothing_contributes_headers_only(tmp_path):
    samples = [
        S("x1", "X", MON, {"A": 1.0}, status="failed"),
        S("n1", "X", MON, {"A": 1.0}, type="negative"),
    ]
    agg, meta, _err = run(tmp_path, samples, "aggregate_name,site_name\nAgg,X\n")
    assert agg == {} and meta == []


def test_missing_sample_column_in_samplesheet_is_an_error(tmp_path):
    p = write_inputs(tmp_path, [S("x1", "X", MON, {"A": 1.0})], "aggregate_name,site_name\nAgg,X\n")
    p["sheet"].write_text("id,site\nx1,X\n")
    result = run_paths(p)
    assert result.returncode == 1 and "no 'sample' column" in result.stderr


def test_samplesheet_without_site_column_contributes_nothing(tmp_path):
    p = write_inputs(tmp_path, [S("x1", "X", MON, {"A": 1.0})], "aggregate_name,site_name\nAgg,X\n")
    p["sheet"].write_text("sample,sample_type,collection_date\nx1,sample,2026-01-05\n")
    result = run_paths(p)
    assert result.returncode == 0
    assert p["meta"].read_text() == "\t".join(sa.META_COLUMNS) + "\n"
    assert "no_site=1" in result.stderr


def test_samplesheet_viral_load_conventions(tmp_path):
    # blank, zero, negative and text loads all mean "unknown"
    samples = [S("x1", "X", MON, {"A": 1.0}, 0), S("y1", "Y", MON, {"A": 1.0}, -5)]
    _, meta, err = run(tmp_path, samples, SITES_XY_NOW)
    assert meta[0]["weight_basis"] == "equal"
    assert "invalid_viral_load=2" in err


def test_cli_rejects_bad_options(tmp_path):
    p = write_inputs(tmp_path, [], "aggregate_name,site_name\nAgg,X\n")
    assert run_paths(p, "--period", "month").returncode == 2
    assert run_paths(p, "--weights", "depth").returncode == 2
    assert run_paths(p, "--min-sites", "0").returncode == 2


# --------------------------------------------------------------------------------------------------
# the realistic fixture tests/data_timeseries
# --------------------------------------------------------------------------------------------------
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


def test_fixture_exclusions_and_weight_basis_by_week(tmp_path):
    _agg, meta, err = run_fixture(tmp_path)
    assert len([m for m in meta if m.startswith("Metro_")]) == 8
    # SiteC_W01_1 has no load -> weeks fall back to the file's site weights
    assert meta["Metro_2026-W02"]["weight_basis"] == "site_weight (viral_load missing)"
    # the negative control (SiteA, 2026-01-12) and the low-coverage sample (not ok) are out
    assert meta["Pair_2026-W03"]["n_samples"] == "4"
    assert (
        meta["Metro_2026-W05"]["n_samples"] == "3"
        and meta["Metro_2026-W05"]["n_sites_included"] == "3"
    )
    assert "negative_controls=1" in err and "not_ok=1" in err and "no_date=1" in err
    # the pair matches "siteb " case-insensitively
    assert meta["Pair_2026-W03"]["sites_included"] == "SiteA;siteb"


def test_fixture_aggregate_never_adds_a_lineage(tmp_path):
    agg, _, _ = run_fixture(tmp_path)
    src, _ = sa.read_aggregate(TS / "freyja_aggregated.tsv")
    universe = {name for row in src.values() for name in row["lin"]}
    for row in agg.values():
        assert set(row["lin"]) <= universe
        assert sum(row["lin"].values()) == approx(1.0, 1e-6)


def test_fixture_drop_a_site_reports_it(tmp_path):
    # a copy of the status table with SiteB_W03_1 failed: SiteB is missing in 2026-W04
    status = (
        (FIX / "demix_status.tsv").read_text().replace("SiteB_W03_1\tok", "SiteB_W03_1\tfailed")
    )
    new = tmp_path / "status.tsv"
    new.write_text(status)
    p = {
        "aggregated": TS / "freyja_aggregated.tsv",
        "status": new,
        "sheet": TS / "samplesheet.csv",
        "sites": FIX / "sites.csv",
        "out": tmp_path / "out.tsv",
        "meta": tmp_path / "meta.tsv",
    }
    assert run_paths(p, "--min-sites", "3").returncode == 0
    _, meta = load_outputs(p)
    metro = {m["period"]: m for m in meta if m["aggregate_name"] == "Metro"}
    assert "2026-W04" not in metro  # only two sites that week -> omitted by --min-sites 3
    p["out"].unlink()
    assert run_paths(p).returncode == 0
    _, meta = load_outputs(p)
    week = next(m for m in meta if m["row_label"] == "Metro_2026-W04")
    assert week["sites_missing"] == "SiteB" and week["missing_reasons"] == "not_ok:failed"
    assert week["weight_basis"] == "viral_load"
    # the failed sample has no viral load, so the share of the weight that is missing is unknown
    assert week["weight_fraction_included"] == "NA"


def test_example_asset_parses():
    aggs = sa.read_sites(REPO / "assets" / "site_aggregation_example.csv")
    assert aggs[0]["name"] == "Metro_Sites"
    assert [s["name"] for s in aggs[0]["sites"]] == ["Site1", "Site2", "Site3"]


def test_check_sites_mode(tmp_path):
    good = tmp_path / "good.csv"
    good.write_text("aggregate_name,site_name\nMetro,Site1\nMetro,Site9\n")
    sheet = tmp_path / "sheet.csv"
    sheet.write_text(SHEET_HEADER + "\nx1,sample,2026-01-05,site1,\n")
    result = subprocess.run(
        [
            sys.executable,
            str(SCRIPT),
            "--check-sites",
            "--sites",
            str(good),
            "--samplesheet",
            str(sheet),
        ],
        capture_output=True,
        text=True,
        check=False,
    )
    assert result.returncode == 0
    assert "sites file OK, 1 aggregate(s)" in result.stdout
    assert "site 'Site9' in the sites file matches no sample" in result.stderr
    assert "Site1'" not in result.stderr
    bad = tmp_path / "bad.csv"
    bad.write_text("aggregate_name,site_name\nMetro,\n")
    result = subprocess.run(
        [sys.executable, str(SCRIPT), "--check-sites", "--sites", str(bad)],
        capture_output=True,
        text=True,
        check=False,
    )
    assert result.returncode == 1 and "line 2: site_name is blank" in result.stderr
    # without --check-sites the other arguments are required
    result = subprocess.run(
        [sys.executable, str(SCRIPT), "--sites", str(good)],
        capture_output=True,
        text=True,
        check=False,
    )
    assert result.returncode == 2
