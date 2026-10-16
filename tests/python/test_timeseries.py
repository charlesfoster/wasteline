import csv
import json
import subprocess
import sys
from pathlib import Path

import pytest
import timeseries as ts

REPO = Path(__file__).resolve().parents[2]
FIXTURE = REPO / "tests" / "data_timeseries"
HEADER = "\tsummarized\tlineages\tabundances\tresid\tcoverage\tpathogen"
SERIES_HEADER = [
    "site",
    "collection_date",
    "lineage_group",
    "abundance",
    "abundance_load_weighted",
    "n_samples",
    "mean_viral_load",
    "total_viral_load",
]


def agg_row(sample, summarized, lineages=None, coverage=99.0):
    """One line of a Freyja aggregate. `summarized` and `lineages` are {name: fraction}."""
    summ = "[" + ", ".join(f"({k!r}, {v!r})" for k, v in summarized.items()) + "]"
    lineages = lineages if lineages is not None else summarized
    names = " ".join(lineages)
    values = " ".join(f"{v:.8f}" for v in lineages.values())
    return f"{sample}\t{summ}\t{names}\t{values}\t1.5\t{coverage}\tSARS-CoV-2"


def write_agg(path, rows):
    path.write_text("\n".join([HEADER, *rows]) + "\n")
    return path


def write_sheet(path, rows, header="sample,sample_type,collection_date,site,viral_load"):
    path.write_text("\n".join([header, *rows]) + "\n")
    return path


def run(tmp_path, agg_rows, sheet_rows, *extra, header=None, capsys=None):
    agg = write_agg(tmp_path / "agg.tsv", agg_rows)
    kwargs = {"header": header} if header else {}
    sheet = write_sheet(tmp_path / "sheet.csv", sheet_rows, **kwargs)
    out = tmp_path / "out.tsv"
    code = ts.main(["--aggregated", str(agg), "--metadata", str(sheet), "-o", str(out), *extra])
    assert code == 0
    lines = out.read_text().splitlines()
    assert lines[0].split("\t") == SERIES_HEADER
    rows = [dict(zip(SERIES_HEADER, line.split("\t"), strict=True)) for line in lines[1:]]
    return rows, capsys.readouterr().err if capsys else ""


def by_group(rows):
    return {(r["site"], r["collection_date"], r["lineage_group"]): r for r in rows}


# ---- parsing --------------------------------------------------------------------------------------
def test_parse_summarized_python_and_numpy_reprs():
    text = "[('Omicron', 0.5), ('XBB.1.5* [Omicron (XBB.1.5.X)]', np.float64(0.25))]"
    assert ts.parse_summarized(text) == [("Omicron", 0.5), ("XBB.1.5* [Omicron (XBB.1.5.X)]", 0.25)]
    assert ts.parse_summarized("[]") == []
    assert ts.parse_summarized("") == []


def test_parse_lineages_and_abundances_in_all_layouts():
    assert ts.parse_name_list("BA.2 AY.4") == ["BA.2", "AY.4"]
    assert ts.parse_name_list("['A' 'Q.3'\n 'AY.48']") == ["A", "Q.3", "AY.48"]
    assert ts.parse_float_list("[0.366636  0.308375\n 0.133008]") == [0.366636, 0.308375, 0.133008]
    assert ts.parse_float_list("0.7 0.3") == [0.7, 0.3]


@pytest.mark.parametrize(
    "text, expected",
    [
        ("2026-01-05", "2026-01-05"),
        ("2026-01-05T10:30:00", "2026-01-05"),
        ("2026-02-30", None),
        ("05/01/2026", None),
        ("2026-1-5", None),
        ("", None),
    ],
)
def test_parse_date(text, expected):
    assert ts.parse_date(text) == expected


def test_aggregate_with_numpy_columns_and_quoted_newlines(tmp_path):
    # Freyja writes numpy reprs with embedded newlines inside quoted cells
    path = tmp_path / "agg.tsv"
    path.write_text(
        HEADER
        + "\nS1\t[('A', 0.6), ('Other', 0.4)]\t\"['A' 'Q.3'\n 'AY.48']\"\t\"[0.5  0.3\n 0.2]\"\t1.0\t99\tSARS-CoV-2\n"
    )
    rows, malformed, _ = ts.read_aggregate(path)
    assert malformed == 0
    assert rows[0]["lineages"] == [("A", 0.5), ("Q.3", 0.3), ("AY.48", 0.2)]
    assert rows[0]["summarized"] == [("A", 0.6), ("Other", 0.4)]


# ---- aggregation ----------------------------------------------------------------------------------
def test_mean_and_load_weighted_abundance_by_hand(tmp_path, capsys):
    # one site-date, loads 100 and 300: weighted B = (0.2*100 + 0.6*300) / 400 = 0.5, mean B = 0.4
    rows, _ = run(
        tmp_path,
        [
            agg_row("S1", {"A": 0.8, "B": 0.2}),
            agg_row("S2", {"A": 0.4, "B": 0.6}),
        ],
        [
            "S1,sample,2026-01-05,Site1,100",
            "S2,sample,2026-01-05,Site1,300",
        ],
        capsys=capsys,
    )
    got = by_group(rows)
    b = got[("Site1", "2026-01-05", "B")]
    assert float(b["abundance"]) == pytest.approx(0.4)
    assert float(b["abundance_load_weighted"]) == pytest.approx(0.5)
    assert b["n_samples"] == "2"
    assert float(b["mean_viral_load"]) == pytest.approx(200)
    assert float(b["total_viral_load"]) == pytest.approx(400)
    a = got[("Site1", "2026-01-05", "A")]
    assert float(a["abundance"]) == pytest.approx(0.6)
    assert float(a["abundance_load_weighted"]) == pytest.approx(0.5)
    # sorted by abundance, descending, inside the site-date
    assert [r["lineage_group"] for r in rows] == ["A", "B"]


def test_group_missing_in_one_sample_counts_as_zero(tmp_path):
    rows, _ = run(
        tmp_path,
        [agg_row("S1", {"A": 1.0}), agg_row("S2", {"A": 0.5, "B": 0.5})],
        ["S1,sample,2026-01-05,Site1,100", "S2,sample,2026-01-05,Site1,100"],
    )
    got = by_group(rows)
    assert float(got[("Site1", "2026-01-05", "B")]["abundance"]) == pytest.approx(0.25)
    assert float(got[("Site1", "2026-01-05", "B")]["abundance_load_weighted"]) == pytest.approx(
        0.25
    )


def test_blank_loads_are_ignored_by_the_weighting_but_counted_as_samples(tmp_path):
    rows, _ = run(
        tmp_path,
        [agg_row("S1", {"A": 1.0}), agg_row("S2", {"A": 0.5, "B": 0.5}), agg_row("S3", {"B": 1.0})],
        [
            "S1,sample,2026-01-05,Site1,",
            "S2,sample,2026-01-05,Site1,100",
            "S3,sample,2026-01-05,Site1,not-a-number",
        ],
    )
    got = by_group(rows)
    b = got[("Site1", "2026-01-05", "B")]
    assert b["n_samples"] == "3"
    assert float(b["abundance"]) == pytest.approx(0.5)  # (0 + 0.5 + 1) / 3
    assert float(b["abundance_load_weighted"]) == pytest.approx(0.5)  # only S2 has a load
    assert float(b["mean_viral_load"]) == pytest.approx(100)
    assert float(b["total_viral_load"]) == pytest.approx(100)


def test_no_loads_at_all_gives_na(tmp_path):
    rows, _ = run(
        tmp_path,
        [agg_row("S1", {"A": 1.0})],
        ["S1,sample,2026-01-05,Site1,"],
    )
    row = rows[0]
    assert (row["abundance_load_weighted"], row["mean_viral_load"], row["total_viral_load"]) == (
        "NA",
        "NA",
        "NA",
    )


def test_zero_or_negative_load_is_treated_as_unknown(tmp_path, capsys):
    rows, err = run(
        tmp_path,
        [agg_row("S1", {"A": 1.0}), agg_row("S2", {"A": 1.0})],
        ["S1,sample,2026-01-05,Site1,0", "S2,sample,2026-01-05,Site1,-5"],
        capsys=capsys,
    )
    assert rows[0]["total_viral_load"] == "NA"
    assert "invalid_viral_load=2" in err


def test_min_coverage_excludes_samples(tmp_path, capsys):
    rows, err = run(
        tmp_path,
        [agg_row("S1", {"A": 1.0}, coverage=59.9), agg_row("S2", {"B": 1.0}, coverage=60.0)],
        ["S1,sample,2026-01-05,Site1,10", "S2,sample,2026-01-05,Site1,10"],
        capsys=capsys,
    )
    assert [r["lineage_group"] for r in rows] == ["B"]
    assert rows[0]["n_samples"] == "1"
    assert "low_coverage=1" in err


def test_min_coverage_option_and_missing_coverage_column(tmp_path):
    rows, _ = run(
        tmp_path,
        [agg_row("S1", {"A": 1.0}, coverage=80)],
        ["S1,sample,2026-01-05,Site1,10"],
        "--min-coverage",
        "90",
    )
    assert rows == []
    # an aggregate without a coverage column cannot be filtered: samples are kept
    agg = tmp_path / "nocov.tsv"
    agg.write_text("\tsummarized\tlineages\tabundances\nS1\t[('A', 1.0)]\tA\t1.0\n")
    sheet = write_sheet(tmp_path / "s.csv", ["S1,sample,2026-01-05,Site1,10"])
    out = tmp_path / "o.tsv"
    assert ts.main(["--aggregated", str(agg), "--metadata", str(sheet), "-o", str(out)]) == 0
    assert len(out.read_text().splitlines()) == 2


def test_negative_controls_are_excluded(tmp_path, capsys):
    rows, err = run(
        tmp_path,
        [agg_row("S1", {"A": 1.0}), agg_row("NEG_1", {"B": 1.0})],
        ["S1,sample,2026-01-05,Site1,10", "NEG_1,Negative,2026-01-05,Site1,10"],
        capsys=capsys,
    )
    assert [r["lineage_group"] for r in rows] == ["A"]
    assert "negative_control=1" in err


def test_sample_without_metadata_is_excluded_and_counted(tmp_path, capsys):
    rows, err = run(
        tmp_path,
        [agg_row("S1", {"A": 1.0}), agg_row("GHOST", {"B": 1.0})],
        ["S1,sample,2026-01-05,Site1,10"],
        capsys=capsys,
    )
    assert [r["lineage_group"] for r in rows] == ["A"]
    assert "no_metadata=1" in err
    assert "2 aggregate samples, 1 used" in err


def test_bad_and_missing_dates_are_skipped_and_counted(tmp_path, capsys):
    rows, err = run(
        tmp_path,
        [agg_row("S1", {"A": 1.0}), agg_row("S2", {"A": 1.0}), agg_row("S3", {"A": 1.0})],
        [
            "S1,sample,2026-01-05,Site1,10",
            "S2,sample,05/01/2026,Site1,10",
            "S3,sample,,Site1,10",
        ],
        capsys=capsys,
    )
    assert len(rows) == 1 and rows[0]["collection_date"] == "2026-01-05"
    assert "bad_date=1" in err and "no_date=1" in err


def test_empty_site_is_na_and_rows_are_sorted(tmp_path):
    rows, _ = run(
        tmp_path,
        [
            agg_row("S1", {"A": 1.0}),
            agg_row("S2", {"A": 0.3, "B": 0.7}),
            agg_row("S3", {"A": 1.0}),
            agg_row("S4", {"A": 1.0}),
        ],
        [
            "S1,sample,2026-01-12,SiteB,",
            "S2,sample,2026-01-05,SiteB,",
            "S3,sample,2026-01-05,,",
            "S4,sample,2026-01-05,SiteA,",
        ],
    )
    keys = [(r["site"], r["collection_date"], r["lineage_group"]) for r in rows]
    assert keys == [
        ("NA", "2026-01-05", "A"),
        ("SiteA", "2026-01-05", "A"),
        ("SiteB", "2026-01-05", "B"),
        ("SiteB", "2026-01-05", "A"),
        ("SiteB", "2026-01-12", "A"),
    ]


def test_same_date_at_two_sites_stays_separate(tmp_path):
    rows, _ = run(
        tmp_path,
        [agg_row("S1", {"A": 1.0}), agg_row("S2", {"B": 1.0})],
        ["S1,sample,2026-01-05,Site1,", "S2,sample,2026-01-05,Site2,"],
    )
    got = by_group(rows)
    assert float(got[("Site1", "2026-01-05", "A")]["abundance"]) == pytest.approx(1.0)
    assert ("Site1", "2026-01-05", "B") not in got


def test_summarized_versus_lineages(tmp_path):
    agg = [agg_row("S1", {"Omicron": 0.9, "Delta": 0.1}, {"BA.2": 0.6, "BA.1": 0.3, "AY.4": 0.1})]
    sheet = ["S1,sample,2026-01-05,Site1,"]
    summarized, _ = run(tmp_path, agg, sheet)
    assert {r["lineage_group"] for r in summarized} == {"Omicron", "Delta"}
    lineages, _ = run(tmp_path, agg, sheet, "--level", "lineages")
    got = by_group(lineages)
    assert {k[2] for k in got} == {"BA.2", "BA.1", "AY.4"}
    assert float(got[("Site1", "2026-01-05", "BA.2")]["abundance"]) == pytest.approx(0.6)


def test_min_abundance_drops_small_groups(tmp_path):
    rows, _ = run(
        tmp_path,
        [agg_row("S1", {"A": 0.97, "B": 0.03})],
        ["S1,sample,2026-01-05,Site1,"],
        "--min-abundance",
        "0.05",
    )
    assert [r["lineage_group"] for r in rows] == ["A"]


def test_optional_metadata_columns_may_be_absent(tmp_path):
    rows, _ = run(
        tmp_path,
        [agg_row("S1", {"A": 1.0})],
        ["S1,2026-01-05"],
        header="sample,collection_date",
    )
    assert rows[0]["site"] == "NA"
    assert rows[0]["abundance_load_weighted"] == "NA"


def test_metadata_header_case_bom_and_duplicates(tmp_path, capsys):
    sheet = tmp_path / "sheet.csv"
    sheet.write_text(
        "﻿Sample,Sample_Type,Collection_Date,Site\nS1,sample,2026-01-05,Site1\nS1,sample,2026-02-05,X\n",
        encoding="utf-8",
    )
    counts = {}
    meta = ts.read_metadata(sheet, counts)
    assert meta["S1"]["date"] == "2026-01-05" and meta["S1"]["site"] == "Site1"
    assert counts == {"duplicate_metadata": 1}


def test_empty_aggregate_writes_header_only(tmp_path, capsys):
    rows, err = run(tmp_path, [], ["S1,sample,2026-01-05,Site1,10"], capsys=capsys)
    assert rows == []
    assert "0 aggregate samples, 0 used" in err


def test_completely_empty_aggregate_file(tmp_path):
    agg = tmp_path / "agg.tsv"
    agg.write_text("")
    sheet = write_sheet(tmp_path / "s.csv", ["S1,sample,2026-01-05,Site1,10"])
    out = tmp_path / "o.tsv"
    assert ts.main(["--aggregated", str(agg), "--metadata", str(sheet), "-o", str(out)]) == 0
    assert out.read_text().splitlines() == ["\t".join(SERIES_HEADER)]


def test_unparsable_aggregate_rows_are_counted_not_fatal(tmp_path, capsys):
    rows, err = run(
        tmp_path,
        [agg_row("S1", {"A": 1.0}), "S2\t[]\t\t\t1.0\t99\tSARS-CoV-2"],
        ["S1,sample,2026-01-05,Site1,", "S2,sample,2026-01-05,Site1,"],
        capsys=capsys,
    )
    assert len(rows) == 1
    assert "unparsable_rows=1" in err


def test_missing_sample_column_is_an_error(tmp_path, capsys):
    agg = write_agg(tmp_path / "agg.tsv", [agg_row("S1", {"A": 1.0})])
    sheet = tmp_path / "sheet.csv"
    sheet.write_text("id,date\nS1,2026-01-05\n")
    code = ts.main(["--aggregated", str(agg), "--metadata", str(sheet), "-o", str(tmp_path / "o")])
    assert code == 1
    assert "no 'sample' column" in capsys.readouterr().err


# ---- Freyja helper mode ---------------------------------------------------------------------------
def test_freyja_helper_files(tmp_path, capsys):
    agg = write_agg(
        tmp_path / "agg.tsv",
        [
            agg_row("S1", {"A": 1.0}),
            agg_row("S2", {"A": 1.0}),
            agg_row("NEG", {"A": 1.0}),
            agg_row("S3", {"A": 1.0}),
            agg_row("GHOST", {"A": 1.0}),
        ],
    )
    sheet = write_sheet(
        tmp_path / "sheet.csv",
        [
            "S1,sample,2026-01-05,Site1,100",
            "S2,sample,2026-01-12,Site1,",
            "NEG,negative,2026-01-12,Site1,",
            "S3,sample,2026-02-12,Site1,200",
        ],
    )
    args = [
        "--aggregated",
        str(agg),
        "--metadata",
        str(sheet),
        "--freyja-metadata-out",
        str(tmp_path / "m.csv"),
        "--freyja-aggregated-out",
        str(tmp_path / "a.tsv"),
        "--freyja-summary-out",
        str(tmp_path / "s.json"),
        "--window-days",
        "30",
    ]
    assert ts.main(args) == 0
    # Freyja's metadata: index column Sample, then sample_collection_datetime, then viral_load
    assert (tmp_path / "m.csv").read_text().splitlines() == [
        "Sample,sample_collection_datetime,viral_load",
        "S1,2026-01-05,100.0",
        "S2,2026-01-12,",
        "S3,2026-02-12,200.0",
    ]
    kept = [line.split("\t")[0] for line in (tmp_path / "a.tsv").read_text().splitlines()]
    assert kept == ["", "S1", "S2", "S3"]  # aggregate order, header kept, NEG and GHOST gone
    summary = json.loads((tmp_path / "s.json").read_text())
    assert (
        summary["n_samples"] == 3 and summary["n_dates"] == 3 and summary["n_with_viral_load"] == 2
    )
    # only 2026-02-12 is within 30 days of the last date
    assert summary["n_dates_window"] == 1 and summary["window_span_days"] == 0

    # scale_by_viral_load: Freyja cannot scale a sample without a load
    args_load = [*args, "--require-viral-load"]
    assert ts.main(args_load) == 0
    kept = [line.split(",")[0] for line in (tmp_path / "m.csv").read_text().splitlines()]
    assert kept == ["Sample", "S1", "S3"]


def test_freyja_helper_omits_load_column_when_no_loads(tmp_path):
    agg = write_agg(tmp_path / "agg.tsv", [agg_row("S1", {"A": 1.0})])
    sheet = write_sheet(tmp_path / "sheet.csv", ["S1,sample,2026-01-05,Site1,"])
    code = ts.main(
        [
            "--aggregated",
            str(agg),
            "--metadata",
            str(sheet),
            "--freyja-metadata-out",
            str(tmp_path / "m.csv"),
            "--freyja-aggregated-out",
            str(tmp_path / "a.tsv"),
        ]
    )
    assert code == 0
    assert (tmp_path / "m.csv").read_text().splitlines()[0] == "Sample,sample_collection_datetime"


def test_cli_requires_an_output(tmp_path):
    with pytest.raises(SystemExit):
        ts.main(["--aggregated", "a", "--metadata", "b"])
    with pytest.raises(SystemExit):
        ts.main(["--aggregated", "a", "--metadata", "b", "--freyja-metadata-out", "m"])


def test_script_runs_as_a_program(tmp_path):
    agg = write_agg(tmp_path / "agg.tsv", [agg_row("S1", {"A": 1.0})])
    sheet = write_sheet(tmp_path / "sheet.csv", ["S1,sample,2026-01-05,Site1,10"])
    out = tmp_path / "o.tsv"
    result = subprocess.run(
        [
            sys.executable,
            str(REPO / "bin" / "timeseries.py"),
            "--aggregated",
            str(agg),
            "--metadata",
            str(sheet),
            "-o",
            str(out),
        ],
        capture_output=True,
        text=True,
        check=False,
    )
    assert result.returncode == 0
    assert "1 aggregate samples, 1 used" in result.stderr
    assert len(out.read_text().splitlines()) == 2


# ---- the synthetic fixture ------------------------------------------------------------------------
def test_fixture_regenerates_byte_identically(tmp_path):
    sys.path.insert(0, str(REPO / "tests"))
    import make_timeseries_fixture as mk

    first, second = tmp_path / "a", tmp_path / "b"
    mk.write_fixture(first, 42)
    mk.write_fixture(second, 42)
    names = sorted(p.relative_to(first).as_posix() for p in first.rglob("*") if p.is_file())
    assert "samplesheet.csv" in names and "truth.json" in names
    for name in names:
        assert (first / name).read_bytes() == (second / name).read_bytes(), name
        # the committed fixture is what the generator writes (freyja_aggregated.tsv comes from Docker)
        assert (first / name).read_bytes() == (FIXTURE / name).read_bytes(), name
    other = tmp_path / "c"
    mk.write_fixture(other, 7)
    assert (other / "samplesheet.csv").read_bytes() != (first / "samplesheet.csv").read_bytes()


def test_fixture_is_small_and_has_the_edge_cases():
    size = sum(p.stat().st_size for p in FIXTURE.rglob("*") if p.is_file())
    assert size < 1_500_000
    with open(FIXTURE / "samplesheet.csv", newline="") as handle:
        rows = list(csv.DictReader(handle))
    assert list(rows[0]) == [
        "sample",
        "demix",
        "sample_type",
        "collection_date",
        "site",
        "viral_load",
    ]
    assert sum(r["sample_type"] == "negative" for r in rows) == 1
    assert sum(r["collection_date"] == "" for r in rows) == 1
    assert sum(r["viral_load"] == "" for r in rows) >= 3
    assert all((REPO / r["demix"]).is_file() for r in rows)
    assert {r["site"] for r in rows} == {"SiteA", "SiteB", "SiteC"}


def test_fixture_timeseries_recovers_the_truth(tmp_path):
    out = tmp_path / "ts.tsv"
    argv = ["--aggregated", str(FIXTURE / "freyja_aggregated.tsv")]
    assert ts.main([*argv, "--metadata", str(FIXTURE / "samplesheet.csv"), "-o", str(out)]) == 0
    rows = [
        dict(zip(SERIES_HEADER, line.split("\t"), strict=True))
        for line in out.read_text().splitlines()[1:]
    ]
    got = by_group(rows)
    truth = json.loads((FIXTURE / "truth.json").read_text())
    rising = "XBB.1.5* [Omicron (XBB.1.5.X)]"
    # the rising lineage increases at every site (first versus last date), the declining one falls
    for site in truth["site_midpoint_day"]:
        assert (
            float(got[(site, "2026-02-23", rising)]["abundance"])
            > float(got[(site, "2026-01-05", rising)]["abundance"]) + 0.4
        )
    # 30 dated, non-negative, well-covered samples reach the table
    used = {(r["site"], r["collection_date"]): int(r["n_samples"]) for r in rows}
    assert sum(used.values()) == 30
    assert len(used) == 24
    # hand-computed SiteA 2026-01-12: XBB.1.5 0.16710542 / 0.23093447 with loads 851825 / 301960
    cell = got[("SiteA", "2026-01-12", rising)]
    assert float(cell["abundance"]) == pytest.approx((0.16710542 + 0.23093447) / 2, abs=1e-6)
    weighted = (0.16710542 * 851825 + 0.23093447 * 301960) / (851825 + 301960)
    assert float(cell["abundance_load_weighted"]) == pytest.approx(weighted, abs=1e-6)
    assert float(cell["abundance_load_weighted"]) != pytest.approx(
        float(cell["abundance"]), abs=1e-3
    )
