"""primer_scheme_check.py: statistics on read 5' ends against primer 5' ends (see the script docstring)."""

import math
from fractions import Fraction
from pathlib import Path

import pytest
from primer_scheme_check import (
    DETAIL_COLUMNS,
    SUMMARY_COLUMNS,
    main,
)

REPO = Path(__file__).resolve().parents[2]
BEDS = REPO / "assets" / "primer_schemes"
FIX = REPO / "tests" / "fixtures" / "qc"


def exact_sf(k, n, p):
    q = Fraction(p).limit_denominator(10**12)
    return float(sum(math.comb(n, i) * q**i * (1 - q) ** (n - i) for i in range(k, n + 1)))


def write_starts(path, pairs):
    """pairs: iterable of (strand, pos, count)."""
    lines = ["strand\tfive_prime_pos"]
    for strand, pos, count in pairs:
        lines += [f"{strand}\t{pos}"] * count
    path.write_text("\n".join(lines) + "\n")
    return str(path)


def run(tmp_path, starts, declared, *extra, beds_dir=None, sample="S"):
    detail, summary = tmp_path / "detail.tsv", tmp_path / "summary.tsv"
    argv = [
        "--starts",
        starts,
        "--declared",
        declared,
        "--sample",
        sample,
        "--beds-dir",
        str(beds_dir or BEDS),
        "-o",
        str(detail),
        "--summary-out",
        str(summary),
        *extra,
    ]
    assert main(argv) == 0
    d_lines = detail.read_text().splitlines()
    assert d_lines[0].split("\t") == DETAIL_COLUMNS
    s_lines = summary.read_text().splitlines()
    assert s_lines[0].split("\t") == SUMMARY_COLUMNS
    assert SUMMARY_COLUMNS[:7] == [
        "sample",
        "declared_scheme",
        "declared_enrichment",
        "best_scheme",
        "best_enrichment",
        "scheme_check",
        "reason",
    ]  # the first seven keep their positions; the hit fractions are appended
    assert SUMMARY_COLUMNS[7:] == ["declared_hit_fraction", "best_hit_fraction"]
    details = {
        r["scheme"]: r for r in (dict(zip(DETAIL_COLUMNS, ln.split("\t"))) for ln in d_lines[1:])
    }
    return details, dict(zip(SUMMARY_COLUMNS, s_lines[1].split("\t")))


def mini_bed_dir(tmp_path):
    """Two plus and one minus site; genome length 100 in the tests."""
    d = tmp_path / "beds"
    d.mkdir()
    (d / "mini.bed").write_text(
        "chr\t10\t30\tp1_F\t1\t+\nchr\t50\t70\tp1_R\t1\t-\nchr\t80\t95\tp2_F\t2\t+\n"
    )
    return d


def test_hand_computed_enrichment_and_p_value(tmp_path):
    beds = mini_bed_dir(tmp_path)
    # W=2, G=100. plus sites 11, 81 -> windows 9-13 and 79-83 (10 positions); minus site 70 -> 68-72 (5).
    # 12 forward reads: 8 at 11 (hit), 1 at 13 (hit, edge), 1 at 14 (miss, just outside), 2 at 50.
    # 8 reverse reads: 7 at 70 (hit), 1 at 11 (same position, wrong strand: miss).
    starts = write_starts(
        tmp_path / "s.tsv",
        [("+", 11, 8), ("+", 13, 1), ("+", 14, 1), ("+", 50, 2), ("-", 70, 7), ("-", 11, 1)],
    )
    details, summary = run(
        tmp_path, starts, "mini", "--genome-length", "100", "--min-reads", "10", beds_dir=beds
    )
    row = details["mini"]
    hits = 8 + 1 + 7
    expected = 12 * 10 / 100 + 8 * 5 / 100  # 1.2 + 0.4
    assert row["n_reads"] == "20"
    assert row["hits"] == str(hits)
    assert float(row["expected_hits"]) == pytest.approx(expected, abs=0.005)
    assert float(row["enrichment"]) == pytest.approx(hits / expected, abs=0.005)
    assert float(row["hit_fraction"]) == pytest.approx(hits / 20, abs=1e-4)
    # pooled binomial: Bin(20, expected/20)
    assert float(row["p_value"]) == pytest.approx(exact_sf(hits, 20, expected / 20), rel=1e-2)
    assert row["equivalent_to_declared"] == "true"
    assert row["n_primer_5prime_sites"] == "3"
    assert summary["scheme_check"] == "ok"
    assert summary["declared_enrichment"] == f"{hits / expected:.1f}"


def test_missing_or_bad_starts_is_na_exit_0(tmp_path):
    """A sample with no usable read starts still gets a summary row with a status and a reason."""
    for content in [None, "", "garbage\n", "a\tb\n+\t5\n"]:
        path = tmp_path / "s.tsv"
        path.unlink(missing_ok=True)
        if content is not None:
            path.write_text(content)
        _, summary = run(tmp_path, str(path), "midnight")
        assert (summary["scheme_check"], summary["reason"]) == ("NA", "missing_input")


def test_real_read_starts_fixture(tmp_path):
    """PURE_C (synthetic, midnight-shaped, anchored) mapped by bwa: tests/fixtures/qc/PURE_C.read_starts.tsv."""
    starts = str(FIX / "PURE_C.read_starts.tsv")
    details, summary = run(tmp_path, starts, "midnight", sample="PURE_C")
    assert summary["scheme_check"] == "ok"
    assert 35 < float(summary["declared_enrichment"]) < 50
    assert float(details["swift"]["enrichment"]) < 3
    _, wrong = run(tmp_path, starts, "articv4.1", sample="PURE_C")
    assert (wrong["scheme_check"], wrong["best_scheme"]) == ("mismatch", "midnight/midnightv2")


def test_sparse_scheme_inflation_fixture(tmp_path):
    """QIAseq-like starts (72 % at qiagenv1 5' ends) plus 48 reads on the 28-site eden scheme: eden's tiny
    expected count gives enrichment 6.8 (above the bar of 5, p ~ 1e-47) at a hit fraction of 1.6 %."""
    starts = str(FIX / "QIAGEN_LIKE_EDEN.read_starts.tsv")
    details, ok = run(tmp_path, starts, "qiagenv1")
    assert float(details["eden"]["enrichment"]) > 5
    assert float(details["eden"]["p_value"]) < 1e-30
    assert float(details["eden"]["hit_fraction"]) == pytest.approx(0.016, abs=0.001)
    assert (ok["scheme_check"], ok["reason"], ok["best_scheme"]) == (
        "ok",
        "declared_enriched",
        "qiagenv1",
    )
    assert ok["declared_hit_fraction"] == ok["best_hit_fraction"] == "0.7170"
    # declared eden: enrichment 6.8 is NOT enough; qiagenv1 has only 2.3x the enrichment, but 45x the hit fraction
    _, eden = run(tmp_path, starts, "eden")
    assert float(eden["declared_enrichment"]) > 5
    assert float(eden["best_enrichment"]) < 3 * float(eden["declared_enrichment"])
    assert (eden["scheme_check"], eden["reason"], eden["best_scheme"]) == (
        "mismatch",
        "best_other_scheme",
        "qiagenv1",
    )
    assert eden["declared_hit_fraction"] == "0.0160"
    # the pre-refinement rule (no floor): eden "reaches the bar" and qiagenv1 is not 3x better: not a mismatch
    _, old = run(tmp_path, starts, "eden", "--min-hit-fraction", "0")
    assert old["scheme_check"] == "ok"
    assert old["best_scheme"] == "qiagenv1"  # the top scheme was always qiagenv1
