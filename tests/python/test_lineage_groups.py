import csv
from pathlib import Path

import lineage_groups as lg
import pytest

REPO = Path(__file__).resolve().parents[2]
FIX = REPO / "tests" / "fixtures" / "lineage_groups"
AGG_HEADER = "\tsummarized\tlineages\tabundances\tresid\tcoverage\tpathogen"
GROUPS = [
    "BA.2.86* (excludes JN.1*)",
    "BA.2.86* (includes JN.1*)",
    "JN.1* (excludes KP.3*)",
    "JN.1* (includes KP.3*)",
    "KP.3*",
    "NB.1.8* (alias XAA.1.7.8*)",
    "Old Omicron, B.1.1.529* only",
    "Recombinants (excludes NB.1*)",
]
TAIL = [
    "other",
    "proportion_summary",
    "major_lineages",
    "major_abundances",
    "minor_lineages",
    "minor_abundances",
    "ungrouped_lineages",
]


def agg_text(samples):
    """samples: {label: [(lineage, abundance), ...]} -> Freyja aggregate text."""
    lines = [AGG_HEADER]
    for label, lineages in samples.items():
        names = " ".join(n for n, _ in lineages)
        values = " ".join(f"{v:.8f}" for _, v in lineages)
        lines.append(f"{label}\t[('Omicron', 1.0)]\t{names}\t{values}\t1.5\t90.0\tSARS-CoV-2")
    return "\n".join(lines) + "\n"


def run(tmp_path, samples, groups_text, *extra, groups_bytes=None, status=None):
    agg = tmp_path / "agg.tsv"
    agg.write_text(agg_text(samples))
    groups = tmp_path / "groups.csv"
    groups.write_bytes(groups_bytes if groups_bytes is not None else groups_text.encode())
    out, resolved = tmp_path / "out.csv", tmp_path / "resolved.tsv"
    argv = [
        "--aggregated",
        str(agg),
        "--groups",
        str(groups),
        "--dataset",
        str(FIX),
        "-o",
        str(out),
        "--resolved-out",
        str(resolved),
        *extra,
    ]
    if status:
        (tmp_path / "status.tsv").write_text(status)
        argv += ["--demix-status", str(tmp_path / "status.tsv")]
    code = lg.main(argv)
    return code, out, resolved


def rows_of(path):
    with open(path, newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


@pytest.fixture(scope="module")
def hierarchy():
    return lg.load_dataset(FIX)[0]


@pytest.fixture(scope="module")
def fixture_run(tmp_path_factory):
    tmp = tmp_path_factory.mktemp("fx")
    out, resolved, sites = tmp / "o.csv", tmp / "r.tsv", tmp / "s.csv"
    argv = [
        "--aggregated",
        str(FIX / "freyja_aggregated.tsv"),
        "--groups",
        str(FIX / "groups.csv"),
        "--dataset",
        str(FIX),
        "--demix-status",
        str(FIX / "demix_status.tsv"),
        "--samplesheet",
        str(FIX / "samplesheet.csv"),
        "--aggregated-sites",
        str(FIX / "site_aggregates_freyja.tsv"),
        "--aggregated-sites-meta",
        str(FIX / "site_aggregates_meta.tsv"),
        "--sites-out",
        str(sites),
        "-o",
        str(out),
        "--resolved-out",
        str(resolved),
    ]
    assert lg.main(argv) == 0
    return {r["sample"]: r for r in rows_of(out)}, out, resolved, sites


def test_nested_include_exclude_overlap_and_everything_else(fixture_run):
    row = fixture_run[0]["Site_1_week_2"]  # underscores in the id; raw label had .variants.tsv
    expected = {
        "site": "Site_1",
        "collection_date": "2026-01-12",
        "status": "ok",
        "genome_coverage_pct": "93.06",
        "BA.2.86* (excludes JN.1*)": "",
        "BA.2.86* (includes JN.1*)": "0.800000",
        "JN.1* (excludes KP.3*)": "0.550000",
        "JN.1* (includes KP.3*)": "0.800000",
        "KP.3*": "0.250000",
        "NB.1.8* (alias XAA.1.7.8*)": "0.019500",
        "Old Omicron, B.1.1.529* only": "0.100000",
        "Recombinants (excludes NB.1*)": "0.050000",
        "other": "0.030500",
        "proportion_summary": "1.000000",
        "major_lineages": "JN.1.16 KP.3.1 BA.5",
        "major_abundances": "0.550000 0.250000 0.100000",  # BA.5 at exactly 0.10 is major
        "minor_lineages": "XAA.1.1 A.1 NB.1.8",  # ZZ.9 (0.0005) is below --min-abundance
        "minor_abundances": "0.050000 0.030000 0.019500",
        "ungrouped_lineages": "A.1:0.030000",
        "status_reason": "",
    }
    assert {k: row[k] for k in expected} == expected


def test_status_rows_negatives_and_missing(fixture_run):
    rows = fixture_run[0]
    assert "NEG_1" not in rows  # negative control per the samplesheet
    assert list(rows) == sorted(rows)
    for sample, status in [
        ("FAIL_5", "failed"),
        ("LOW_6", "skipped_low_coverage"),
        ("MISS_7", "missing"),
        ("NOAGG_8", "missing"),
    ]:
        row = rows[sample]
        assert row["status"] == status and row["status_reason"]
        assert all(row[c] == "" for c in [*GROUPS, *TAIL])
    assert rows["FAIL_5"]["genome_coverage_pct"] == "88.10"
    assert rows["MISS_7"]["genome_coverage_pct"] == ""
    assert rows["NOAGG_8"]["site"] == "" and rows["NOAGG_8"]["collection_date"] == ""


def test_descendants_and_patterns_follow_aliases(hierarchy):
    assert hierarchy.canonical("xaa.1.7.8") == "NB.1.8"
    assert "NB.1.8" in hierarchy.descendants["XAA.1"]  # NB.1 is XAA.1.7 by alias
    problems = []
    assert lg.resolve_pattern("JN.1*", hierarchy, problems) == {
        "JN.1", "JN.1.16", "KP.3", "KP.3.1", "LP.8.1"
    }  # fmt: skip
    assert lg.resolve_pattern("ZZ.1*", hierarchy, problems) == set()
    assert len(problems) == 1  # an unknown lineage is a warning, not a match


def test_bad_groups_file_aborts_early_without_outputs(tmp_path, capsys):
    text = "group,include\nJN,JN.1*\nJN,KP.3*\n"
    code, out, resolved = run(tmp_path, {"S": [("A.1", 1.0)]}, text)
    assert code == 2
    err = capsys.readouterr().err
    assert "appears twice" in err and "Traceback" not in err
    assert not out.exists() and not resolved.exists()
