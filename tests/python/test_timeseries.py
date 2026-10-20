import json
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


def test_unparsable_aggregate_rows_are_counted_not_fatal(tmp_path, capsys):
    rows, err = run(
        tmp_path,
        [agg_row("S1", {"A": 1.0}), "S2\t[]\t\t\t1.0\t99\tSARS-CoV-2"],
        ["S1,sample,2026-01-05,Site1,", "S2,sample,2026-01-05,Site1,"],
        capsys=capsys,
    )
    assert len(rows) == 1
    assert "unparsable_rows=1" in err


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
