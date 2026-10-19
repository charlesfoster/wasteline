import csv
import hashlib
import subprocess
import sys
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


def tsv_rows(path):
    with open(path, newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


# ---- hierarchy ----------------------------------------------------------------------------------
def test_minimal_yaml_reader_matches_pyyaml():
    yaml = pytest.importorskip("yaml")
    text = (FIX / "lineages.yml").read_text()
    assert lg.parse_lineages_yml(FIX / "lineages.yml") == yaml.safe_load(text)


@pytest.fixture(scope="module")
def hierarchy():
    return lg.load_dataset(FIX)[0]


def test_dataset_version_and_missing_hierarchy(tmp_path):
    assert lg.load_dataset(FIX)[1] == "01_02_2026-03-04"
    with pytest.raises(lg.UserError, match="no lineages.yml"):
        lg.load_dataset(tmp_path)


def test_aliases_resolve_both_ways_and_ignore_case(hierarchy):
    assert hierarchy.canonical("NB.1.8") == "NB.1.8"
    assert hierarchy.canonical("xaa.1.7.8") == "NB.1.8"
    assert hierarchy.canonical(" b.1.1.529.2.86.1.1 ") == "JN.1"
    assert hierarchy.canonical("nope") is None


def test_descendants_follow_parent_links_through_aliases(hierarchy):
    jn = hierarchy.descendants["JN.1"]
    assert jn == {"JN.1", "JN.1.16", "KP.3", "KP.3.1", "LP.8.1"}
    assert "NB.1.8" in hierarchy.descendants["XAA.1"]  # NB.1 is XAA.1.7 by alias
    assert "BA.5" not in hierarchy.descendants["BA.2.86"]


def test_recombinants_are_whatever_descends_from_a_root_with_recombinant_parents(hierarchy):
    assert hierarchy.recombinants == {"XAA", "XAA.1", "XAA.1.1", "NB.1", "NB.1.8", "XBB", "XBB.1"}


@pytest.mark.parametrize(
    "pattern, expected",
    [
        ("JN.1*", {"JN.1", "JN.1.16", "KP.3", "KP.3.1", "LP.8.1"}),
        ("JN.1", {"JN.1"}),
        ("  kp.3* ", {"KP.3", "KP.3.1"}),
        ("B.1.1.529.2.86.1.1.16", {"JN.1.16"}),
        ("RECOMBINANT*", {"XAA", "XAA.1", "XAA.1.1", "NB.1", "NB.1.8", "XBB", "XBB.1"}),
        ("recombinant", {"XAA", "XAA.1", "XAA.1.1", "NB.1", "NB.1.8", "XBB", "XBB.1"}),
    ],
)
def test_resolve_pattern(hierarchy, pattern, expected):
    problems = []
    assert lg.resolve_pattern(pattern, hierarchy, problems) == expected
    assert problems == []


@pytest.mark.parametrize("pattern", ["ZZ.1*", "ZZ.1", "X*", "BA.*.86", "*"])
def test_unknown_or_invalid_pattern_matches_nothing(hierarchy, pattern):
    problems = []
    assert lg.resolve_pattern(pattern, hierarchy, problems) == set()
    assert len(problems) == 1


def test_x_pattern_hints_at_the_recombinant_token(hierarchy):
    problems = []
    lg.resolve_pattern("X*", hierarchy, problems)
    assert "RECOMBINANT*" in problems[0]


# ---- the fixture run ------------------------------------------------------------------------------
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


def test_header_is_the_documented_order(fixture_run):
    _, out, _, _ = fixture_run
    header = out.read_text(encoding="utf-8").splitlines()[0]
    parsed = next(csv.reader([header]))
    assert parsed == [
        "sample",
        "site",
        "collection_date",
        "status",
        "genome_coverage_pct",
        *GROUPS,
        *TAIL,
        "status_reason",
    ]
    assert not out.read_bytes().startswith(b"\xef\xbb\xbf")


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


def test_recombinant_token_row_and_zero_blanks(fixture_run):
    row = fixture_run[0]["ALLX_3"]
    assert row["Recombinants (excludes NB.1*)"] == "1.000000"
    assert [row[g] for g in GROUPS if g != "Recombinants (excludes NB.1*)"] == [""] * 7
    assert row["other"] == "" and row["proportion_summary"] == "1.000000"


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


def test_resolved_audit(fixture_run):
    _, _, resolved, _ = fixture_run
    rows = tsv_rows(resolved)
    assert [r["group"] for r in rows] == GROUPS
    sha = hashlib.sha256((FIX / "groups.csv").read_bytes()).hexdigest()
    assert {r["groups_sha256"] for r in rows} == {sha}
    assert {r["dataset_version"] for r in rows} == {"01_02_2026-03-04"}
    by = {r["group"]: r for r in rows}
    assert by["BA.2.86* (excludes JN.1*)"]["n_lineages_matched"] == "2"
    assert by["BA.2.86* (includes JN.1*)"]["in_total"] == "no"
    assert by["JN.1* (excludes KP.3*)"]["exclude"] == "KP.3*"
    assert all(r["warnings"] == "" for r in rows)


def test_site_aggregates(fixture_run):
    rows = {r["period"]: r for r in rows_of(fixture_run[3])}
    header = next(csv.reader(fixture_run[3].read_text(encoding="utf-8").splitlines()))
    assert header[:7] == [
        "aggregate_name",
        "period",
        "period_start",
        "n_sites_included",
        "sites_included",
        "sites_missing",
        "weight_basis",
    ]
    assert header[7:] == [*GROUPS, *TAIL]
    w2 = rows["2026-W02"]
    assert w2["aggregate_name"] == "Catchment_X" and w2["sites_included"] == "Site_1;Site_2"
    assert w2["JN.1* (excludes KP.3*)"] == "0.500000" and w2["KP.3*"] == "0.300000"
    assert w2["Recombinants (excludes NB.1*)"] == "0.200000"
    assert w2["major_lineages"] == "JN.1.16 KP.3 XAA.1"
    w3 = rows["2026-W03"]
    assert w3["sites_missing"] == "Site_2" and w3["NB.1.8* (alias XAA.1.7.8*)"] == "1.000000"
    w4 = rows["2026-W04"]  # in the meta but not in the aggregate: kept, blank values
    assert w4["n_sites_included"] == "0" and w4["proportion_summary"] == ""


def test_site_flags_must_come_together(capsys):
    with pytest.raises(SystemExit):
        lg.parse_args(
            ["--aggregated", "a", "--groups", "g", "--dataset", "d", "-o", "o", "--resolved-out",
             "r", "--aggregated-sites", "x"]
        )  # fmt: skip


# ---- behaviour variants --------------------------------------------------------------------------
SIMPLE = "group,include,exclude,in_total\nJN,JN.1*,KP.3*,yes\nKP,KP.3*,,yes\n"
OVERLAP = "group,include,exclude,in_total\nJN,JN.1*,,yes\nKP,KP.3*,,yes\n"


def test_overlap_of_in_total_groups_warns_once_and_total_exceeds_one(tmp_path, capsys):
    code, out, resolved = run(tmp_path, {"S": [("KP.3.1", 0.6), ("A.1", 0.4)]}, OVERLAP)
    assert code == 0
    err = capsys.readouterr().err
    assert err.count("overlaps in_total group 'JN'") == 1
    row = rows_of(out)[0]
    assert row["JN"] == "0.600000" and row["KP"] == "0.600000"
    assert float(row["proportion_summary"]) == pytest.approx(1.6)
    warnings = {r["group"]: r["warnings"] for r in tsv_rows(resolved)}
    assert "KP.3" in warnings["KP"] and warnings["JN"] == ""


def test_overlap_marked_not_in_total_is_fine(tmp_path, capsys):
    text = "group,include,exclude,in_total\nJN,JN.1*,,no\nKP,KP.3*,,yes\n"
    code, out, _ = run(tmp_path, {"S": [("KP.3.1", 0.6), ("A.1", 0.4)]}, text)
    assert code == 0 and capsys.readouterr().err == ""
    row = rows_of(out)[0]
    assert float(row["proportion_summary"]) == pytest.approx(1.0)
    assert row["other"] == "0.400000"


def test_unknown_pattern_warns_once_and_matches_nothing(tmp_path, capsys):
    text = "group,include,exclude\nG1,ZZ.1*;KP.3*,\nG2,ZZ.1*,\nG3,X*,\n"
    code, out, resolved = run(tmp_path, {"S": [("KP.3", 0.5), ("A.1", 0.5)]}, text)
    assert code == 0
    err = capsys.readouterr().err
    assert err.count("'ZZ.1'") == 1  # reported once for the pattern, not per group
    assert "RECOMBINANT*" in err
    audit = {r["group"]: r for r in tsv_rows(resolved)}
    assert "ZZ.1" in audit["G1"]["warnings"] and "ZZ.1" in audit["G2"]["warnings"]
    assert audit["G2"]["n_lineages_matched"] == "0"
    row = rows_of(out)[0]
    assert row["G1"] == "0.500000" and row["G2"] == "" and row["G3"] == ""


def test_group_that_matches_nothing_is_flagged(tmp_path, capsys):
    text = "group,include,exclude\nG,KP.3*,KP.3*\n"
    code, _, resolved = run(tmp_path, {"S": [("KP.3", 1.0)]}, text)
    assert code == 0
    assert "matches no lineage" in tsv_rows(resolved)[0]["warnings"]


def test_case_and_space_handling(tmp_path):
    text = "  group  , include ,EXCLUDE,In_Total\n Mine , jn.1* ; Kp.3* ,  lp.8.1  ,YES\n"
    code, out, _ = run(
        tmp_path, {"S": [("jn.1.16", 0.3), ("KP.3.1", 0.3), ("LP.8.1", 0.2), ("A.1", 0.2)]}, text
    )
    assert code == 0
    row = rows_of(out)[0]
    assert row["Mine"] == "0.600000"  # trimmed name; LP.8.1 excluded; lower-case lineage matches
    assert row["major_lineages"] == "jn.1.16 KP.3.1 A.1 LP.8.1"


def test_zeros_mode_and_decimals(tmp_path):
    samples = {"S": [("JN.1.16", 0.123456789), ("A.1", 0.876543211)]}
    code, out, _ = run(tmp_path, samples, SIMPLE, "--zeros", "zero", "--decimals", "3")
    assert code == 0
    row = rows_of(out)[0]
    assert row["JN"] == "0.123" and row["KP"] == "0.000" and row["other"] == "0.877"
    assert row["major_abundances"] == "0.877 0.123"
    # blank mode: a value that rounds to zero is blank too
    _code, out, _ = run(tmp_path, {"S": [("JN.1.16", 0.0000001), ("A.1", 0.9999999)]}, SIMPLE)
    assert rows_of(out)[0]["JN"] == ""


def test_failed_rows_stay_blank_even_with_zeros_mode(tmp_path):
    status = "sample\tstatus\treason\tgenome_coverage_pct\nS\tfailed\tboom\t5\nT\tok\t\t90\n"
    _code, out, _ = run(tmp_path, {"T": [("A.1", 1.0)]}, SIMPLE, "--zeros", "zero", status=status)
    rows = {r["sample"]: r for r in rows_of(out)}
    assert rows["S"]["JN"] == "" and rows["S"]["status_reason"] == "boom"
    assert rows["T"]["JN"] == "0.000000"


def test_major_minor_threshold_ordering_and_ties(tmp_path):
    samples = {"S": [("B.1", 0.05), ("A.1", 0.05), ("BA.5", 0.1), ("BA.2", 0.2), ("JN.1", 0.6)]}
    _code, out, _ = run(tmp_path, samples, SIMPLE)
    row = rows_of(out)[0]
    assert row["major_lineages"] == "JN.1 BA.2 BA.5"
    assert row["major_abundances"] == "0.600000 0.200000 0.100000"
    assert row["minor_lineages"] == "A.1 B.1"  # equal abundances: alphabetical
    _code, out, _ = run(tmp_path, samples, SIMPLE, "--major-threshold", "0.15")
    assert rows_of(out)[0]["major_lineages"] == "JN.1 BA.2"
    _code, out, _ = run(tmp_path, samples, SIMPLE, "--min-abundance", "0.06")
    assert rows_of(out)[0]["minor_lineages"] == ""
    # group sums are not affected by --min-abundance
    assert rows_of(out)[0]["other"] == "0.400000"


def test_numpy_and_bracketed_aggregate_forms(tmp_path):
    agg = tmp_path / "agg.tsv"
    agg.write_text(
        AGG_HEADER + "\nS1.variants.tsv\t[]\t['JN.1.16' 'A.1']\t[0.6  0.4\n ]\t1\t90\tSARS-CoV-2\n"
    )
    groups = tmp_path / "g.csv"
    groups.write_text(SIMPLE)
    out = tmp_path / "o.csv"
    code = lg.main(
        ["--aggregated", str(agg), "--groups", str(groups), "--dataset", str(FIX), "-o", str(out),
         "--resolved-out", str(tmp_path / "r.tsv")]
    )  # fmt: skip
    assert code == 0
    row = rows_of(out)[0]
    assert row["sample"] == "S1" and row["JN"] == "0.600000" and row["status"] == "ok"


def test_unparsable_aggregate_row_is_a_failed_row(tmp_path):
    agg = tmp_path / "agg.tsv"
    agg.write_text(AGG_HEADER + "\nS1\t[]\tJN.1 A.1\t0.6\t1\t90\tSARS-CoV-2\n")
    groups = tmp_path / "g.csv"
    groups.write_text(SIMPLE)
    out = tmp_path / "o.csv"
    code = lg.main(
        ["--aggregated", str(agg), "--groups", str(groups), "--dataset", str(FIX), "-o", str(out),
         "--resolved-out", str(tmp_path / "r.tsv")]
    )  # fmt: skip
    assert code == 0
    row = rows_of(out)[0]
    assert row["status"] == "failed" and "parsed" in row["status_reason"]


def test_label_column_and_no_demix_status(tmp_path):
    code, out, _ = run(
        tmp_path, {"S_1.variants.tsv": [("A.1", 1.0)]}, SIMPLE, "--label-column", "id"
    )
    assert code == 0
    row = rows_of(out)[0]
    assert row["id"] == "S_1" and row["status"] == "ok" and row["genome_coverage_pct"] == "90.00"


# ---- groups file ----------------------------------------------------------------------------------
def test_excel_style_groups_file(tmp_path):
    """BOM, CRLF, blank padded rows and cells, a semicolon-delimited (European Excel) save."""
    body = "group;include;exclude;in_total\r\nJN;JN.1*;;yes\r\n;;;\r\n"
    code, out, _ = run(
        tmp_path, {"S": [("JN.1", 1.0)]}, "", groups_bytes=b"\xef\xbb\xbf" + body.encode()
    )
    assert code == 0 and rows_of(out)[0]["JN"] == "1.000000"
    padded = "group,include,exclude,in_total,,\r\nJN,JN.1*,,,,\r\n,,,,,\r\n"
    code, out, _ = run(tmp_path, {"S": [("JN.1", 1.0)]}, padded)
    assert code == 0 and rows_of(out)[0]["JN"] == "1.000000"  # in_total blank = yes


@pytest.mark.parametrize(
    "text, message",
    [
        ("", "is empty"),
        ("group,include\n", "no groups"),
        ("name,patterns\nJN,JN.1*\n", "'group' and 'include'"),
        ("group,include\nJN,JN.1*\nJN,KP.3*\n", "appears twice"),
        ("group,include\nJN,JN.1*\njn,KP.3*\n", "appears twice"),
        ("group,include,in_total\nJN,JN.1*,maybe\n", "in_total must be yes or no"),
        ("group,include\nJN,\n", "no 'include' pattern"),
        ("group,include\n,JN.1*\n", "'group' cell is empty"),
        ("group,include\nother,JN.1*\n", "reserved"),
    ],
)
def test_groups_file_errors_are_friendly(tmp_path, capsys, text, message):
    code, out, resolved = run(tmp_path, {"S": [("A.1", 1.0)]}, text)
    assert code == 2
    err = capsys.readouterr().err
    assert message in err and "Traceback" not in err
    assert not out.exists() and not resolved.exists()


def test_missing_hierarchy_is_a_friendly_error(tmp_path, capsys):
    (tmp_path / "agg.tsv").write_text(agg_text({"S": [("A.1", 1.0)]}))
    (tmp_path / "g.csv").write_text(SIMPLE)
    code = lg.main(
        ["--aggregated", str(tmp_path / "agg.tsv"), "--groups", str(tmp_path / "g.csv"),
         "--dataset", str(tmp_path), "-o", str(tmp_path / "o.csv"),
         "--resolved-out", str(tmp_path / "r.tsv")]
    )  # fmt: skip
    assert code == 2 and "no lineages.yml" in capsys.readouterr().err


def test_command_line_end_to_end(tmp_path):
    out, resolved = tmp_path / "o.csv", tmp_path / "r.tsv"
    proc = subprocess.run(
        [sys.executable, str(REPO / "bin" / "lineage_groups.py"),
         "--aggregated", str(FIX / "freyja_aggregated.tsv"), "--groups", str(FIX / "groups.csv"),
         "--dataset", str(FIX), "-o", str(out), "--resolved-out", str(resolved)],
        capture_output=True,
        text=True,
        check=False,
    )  # fmt: skip
    assert proc.returncode == 0
    assert "ZZ.9" in proc.stderr
    assert {r["sample"] for r in rows_of(out)} == {"Site_1_week_2", "ALLX_3", "NEG_1", "Mixed_4"}


def test_check_groups_needs_no_aggregate_and_writes_nothing(tmp_path, capsys):
    groups = tmp_path / "groups.csv"
    groups.write_text("group,include,exclude,in_total\nJN.1*,JN.1*,KP.3*,yes\nTypo,JN.99999*,,no\n")
    rc = lg.main(["--check-groups", "--groups", str(groups), "--dataset", str(FIX)])
    captured = capsys.readouterr()
    assert rc == 0
    assert "groups file OK, 2 group(s)" in captured.out
    assert (
        "WARNING" in captured.err and "JN.99999" in captured.err
    )  # unknown lineage: a warning only
    assert list(tmp_path.iterdir()) == [groups]


def test_check_groups_format_error_exits_2_naming_the_row(tmp_path, capsys):
    groups = tmp_path / "groups.csv"
    groups.write_text("group,include,exclude,in_total\nJN.1*,JN.1*,,yes\nJN.1*,KP.3*,,yes\n")
    rc = lg.main(["--check-groups", "--groups", str(groups), "--dataset", str(FIX)])
    assert rc == 2
    assert "rows 2 and 3" in capsys.readouterr().err


def test_check_groups_without_hierarchy_is_an_error(tmp_path, capsys):
    groups = tmp_path / "groups.csv"
    groups.write_text("group,include\nJN.1*,JN.1*\n")
    rc = lg.main(["--check-groups", "--groups", str(groups), "--dataset", str(tmp_path)])
    assert rc == 2
    assert "lineages.yml" in capsys.readouterr().err


def test_normal_mode_still_requires_its_outputs(tmp_path):
    with pytest.raises(SystemExit):
        lg.parse_args(["--groups", "g.csv", "--dataset", "d"])


def test_check_groups_without_dataset_validates_the_format_only(tmp_path, capsys):
    good = tmp_path / "good.csv"
    good.write_text("group,include\nJN.1*,JN.1*\n")
    assert lg.main(["--check-groups", "--groups", str(good)]) == 0
    assert "patterns were not resolved" in capsys.readouterr().out
    bad = tmp_path / "bad.csv"
    bad.write_text("group,include\nJN.1*,\n")
    assert lg.main(["--check-groups", "--groups", str(bad)]) == 2
    assert "row 2" in capsys.readouterr().err


def test_ungrouped_summary_lists_the_lineages_in_other_biggest_first(tmp_path):
    out, resolved, ungrouped = tmp_path / "o.csv", tmp_path / "r.tsv", tmp_path / "u.tsv"
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
        "-o",
        str(out),
        "--resolved-out",
        str(resolved),
        "--ungrouped-out",
        str(ungrouped),
    ]
    assert lg.main(argv) == 0
    lines = ungrouped.read_text().splitlines()
    assert lines[0].split("\t") == [
        "lineage",
        "n_samples",
        "mean_abundance",
        "max_abundance",
        "suggested_pattern",
    ]
    rows = tsv_rows(ungrouped)
    assert rows, "the fixture has an ungrouped lineage (A.1 in Site_1_week_2)"
    a1 = next(r for r in rows if r["lineage"] == "A.1")
    assert a1["suggested_pattern"] == "A.1*" and int(a1["n_samples"]) >= 1
    assert float(a1["max_abundance"]) >= float(a1["mean_abundance"])
    maxima = [float(r["max_abundance"]) for r in rows]
    assert maxima == sorted(maxima, reverse=True)
    # every listed lineage appears in some sample's ungrouped_lineages column
    per_sample = " ".join(r["ungrouped_lineages"] for r in rows_of(out))
    assert all(f"{r['lineage']}:" in per_sample for r in rows)
    assert not any(r["lineage"] == "ZZ.9" for r in rows)  # below --min-abundance


def test_ungrouped_summary_is_header_only_when_everything_is_grouped(tmp_path):
    code, _out, _resolved = run(
        tmp_path,
        {"S1": [("JN.1.16", 0.6), ("KP.3.1", 0.4)]},
        (FIX / "groups.csv").read_text(),
        "--ungrouped-out",
        str(tmp_path / "u.tsv"),
    )
    assert code == 0
    assert (tmp_path / "u.tsv").read_text().splitlines() == [
        "lineage\tn_samples\tmean_abundance\tmax_abundance\tsuggested_pattern"
    ]


def test_ungrouped_summary_counts_samples_and_excludes_negatives(tmp_path):
    groups = "group,include\nKP.3*,KP.3*\n"
    sheet = tmp_path / "sheet.csv"
    sheet.write_text("sample,sample_type\nS1,sample\nS2,sample\nNEG,negative\n")
    code, _out, _resolved = run(
        tmp_path,
        {
            "S1": [("KP.3.1", 0.5), ("A.1", 0.5)],
            "S2": [("KP.3.1", 0.8), ("A.1", 0.2)],
            "NEG": [("A.1", 1.0)],
        },
        groups,
        "--samplesheet",
        str(sheet),
        "--ungrouped-out",
        str(tmp_path / "u.tsv"),
        "--decimals",
        "2",
    )
    assert code == 0
    rows = tsv_rows(tmp_path / "u.tsv")
    assert rows == [
        {
            "lineage": "A.1",
            "n_samples": "2",
            "mean_abundance": "0.35",
            "max_abundance": "0.50",
            "suggested_pattern": "A.1*",
        }
    ]


def test_print_default_groups_prints_the_bundled_file(capsys):
    assert lg.main(["--print-default-groups"]) == 0
    printed = capsys.readouterr().out
    assert printed == (REPO / "assets" / "lineage_groups.csv").read_text(encoding="utf-8")
    assert printed.startswith("group,include,exclude,in_total")


def test_the_bundled_default_groups_file_is_valid_and_resolves_against_the_fixture(tmp_path):
    groups, _ = lg.read_groups(REPO / "assets" / "lineage_groups.csv")
    assert len(groups) >= 10 and all(g["include"] for g in groups)
