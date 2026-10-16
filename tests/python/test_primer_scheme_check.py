"""primer_scheme_check.py: statistics on read 5' ends against primer 5' ends (see the script docstring)."""

import math
import random
from fractions import Fraction
from pathlib import Path

import pytest
from primer_scheme_check import (
    DETAIL_COLUMNS,
    SUMMARY_COLUMNS,
    binom_sf,
    jaccard,
    main,
    read_bed_sites,
    read_starts,
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


# ---- binomial tail ---------------------------------------------------------------------------------------


@pytest.mark.parametrize(
    ("k", "n", "p"),
    [
        (18, 20, 0.05),  # hand case
        (0, 10, 0.3),
        (1, 10, 0.3),
        (3, 10, 0.3),  # k <= mode: complement of the lower tail
        (250, 1000, 0.3),
        (330, 1000, 0.3),
        (12, 500, 0.01),
        (40, 500, 0.01),
        (500, 500, 0.5),
    ],
)
def test_binom_sf_matches_exact_sum(k, n, p):
    assert binom_sf(k, n, p) == pytest.approx(exact_sf(k, n, p), rel=1e-6, abs=1e-12)


def test_binom_sf_edges():
    assert binom_sf(0, 5, 0.1) == 1.0
    assert binom_sf(6, 5, 0.1) == 1e-300
    assert binom_sf(1, 5, 0.0) == 1e-300
    assert binom_sf(1, 5, 1.0) == 1.0
    assert binom_sf(100000, 200000, 0.001) == 1e-300  # underflow is floored, not 0 or an error


def test_binom_sf_large_n_is_fast_and_sane():
    p = binom_sf(60000, 1_000_000, 0.05)  # mean 50000, sd ~218: ~45 sd out
    assert 0 < p <= 1e-300 or p == 1e-300
    assert binom_sf(50000, 1_000_000, 0.05) == pytest.approx(0.5, abs=0.01)


# ---- BED parsing -----------------------------------------------------------------------------------------


def test_bed_five_prime_ends_and_strand(tmp_path):
    sites, skipped = read_bed_sites(mini_bed_dir(tmp_path) / "mini.bed")
    # + primer: BED start + 1 ; - primer: BED end
    assert sites == {("+", 11), ("-", 70), ("+", 81)}
    assert skipped == 0


def test_bed_strand_from_name_when_column_missing(tmp_path):
    bed = tmp_path / "x.bed"
    bed.write_text(
        "c\t3\t25\tcovid19genome_0-200_s0_M1F\nc\t301\t322\tcovid19genome_200-29703_s7490_U_88R\n"
        "c\t5\t9\tQIAseq_1_LEFT\nc\t9\t19\tfoo_RIGHT_alt1\nc\t1\t2\tmystery\nc\tx\ty\tbad_F\n"
    )
    sites, skipped = read_bed_sites(bed)
    assert sites == {("+", 4), ("-", 322), ("+", 6), ("-", 19)}
    assert skipped == 2


def test_bed_strand_column_wins_over_name(tmp_path):
    bed = tmp_path / "x.bed"
    bed.write_text("c\t10\t20\tweird_R\t1\t+\n")  # column 6 is authoritative
    assert read_bed_sites(bed)[0] == {("+", 11)}


def test_swift_seventh_column_is_ignored():
    sites, skipped = read_bed_sites(BEDS / "swift.bed")
    assert skipped == 0
    assert {s for s, _ in sites} == {"+", "-"}


@pytest.mark.parametrize("path", sorted(BEDS.glob("*.bed")), ids=lambda p: p.stem)
def test_all_bundled_schemes_parse_with_both_strands(path):
    sites, skipped = read_bed_sites(path)
    assert skipped == 0
    plus = [p for s, p in sites if s == "+"]
    minus = [p for s, p in sites if s == "-"]
    assert len(plus) >= 10
    assert len(minus) >= 10


def test_midnight_and_midnightv2_are_near_identical_but_not_equal():
    a, _ = read_bed_sites(BEDS / "midnight.bed")
    b, _ = read_bed_sites(BEDS / "midnightv2.bed")
    assert a < b  # v2 has one extra reverse site
    assert jaccard(a, b) == pytest.approx(58 / 59)
    assert jaccard(a, read_bed_sites(BEDS / "articv4.1.bed")[0]) < 0.5


# ---- read starts -----------------------------------------------------------------------------------------


def test_read_starts_skips_bad_lines(tmp_path):
    path = tmp_path / "s.tsv"
    path.write_text("strand\tfive_prime_pos\n+\t5\n+\t5\n-\tx\n?\t3\n-\t9\nshort\n")
    counts = read_starts(path)
    assert dict(counts) == {("+", 5): 2, ("-", 9): 1}


def test_read_starts_header_required(tmp_path):
    path = tmp_path / "s.tsv"
    path.write_text("a\tb\n+\t5\n")
    with pytest.raises(ValueError):
        read_starts(path)


# ---- hand-computed statistics on a mini scheme -----------------------------------------------------------


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


def test_window_edges_are_inclusive_and_window_option(tmp_path):
    beds = mini_bed_dir(tmp_path)
    starts = write_starts(tmp_path / "s.tsv", [("+", 14, 5), ("+", 9, 5), ("+", 8, 5)])
    narrow, _ = run(
        tmp_path, starts, "mini", "--genome-length", "100", "--window", "2", beds_dir=beds
    )
    wide, _ = run(
        tmp_path, starts, "mini", "--genome-length", "100", "--window", "3", beds_dir=beds
    )
    assert narrow["mini"]["hits"] == "5"  # only 9 (= 11 - 2)
    assert wide["mini"]["hits"] == "15"  # 8, 9 and 14 (= 11 +- 3)


def test_overlapping_windows_are_not_double_counted(tmp_path):
    beds = tmp_path / "beds"
    beds.mkdir()
    (beds / "dense.bed").write_text(
        "c\t9\t20\ta_F\t1\t+\nc\t11\t30\tb_F\t1\t+\n"
    )  # sites 10 and 12
    starts = write_starts(tmp_path / "s.tsv", [("+", 50, 10)])
    details, _ = run(
        tmp_path, starts, "dense", "--genome-length", "100", "--min-reads", "1", beds_dir=beds
    )
    # window 8-14 = 7 positions (not 2 * 5); expected = 10 * 7 / 100
    assert float(details["dense"]["expected_hits"]) == pytest.approx(0.7)


def test_strand_must_match(tmp_path):
    beds = mini_bed_dir(tmp_path)
    starts = write_starts(
        tmp_path / "s.tsv", [("-", 11, 30), ("+", 70, 30)]
    )  # right places, wrong strands
    details, summary = run(
        tmp_path, starts, "mini", "--genome-length", "100", "--min-reads", "10", beds_dir=beds
    )
    assert details["mini"]["hits"] == "0"
    assert float(details["mini"]["enrichment"]) == 0.0
    assert summary["scheme_check"] == "uninformative"


def test_zero_expected_hits_gives_na_enrichment(tmp_path):
    beds = tmp_path / "beds"
    beds.mkdir()
    (beds / "empty.bed").write_text("# nothing usable\n")
    starts = write_starts(tmp_path / "s.tsv", [("+", 5, 600)])
    details, summary = run(tmp_path, starts, "empty", beds_dir=beds)
    assert details["empty"]["enrichment"] == "NA"
    assert details["empty"]["p_value"] == "NA"
    assert summary["scheme_check"] == "uninformative"


# ---- decisions: ok / mismatch / uninformative / NA -------------------------------------------------------


def starts_from_sites(tmp_path, scheme, anchored, total, seed=1, genome=29903, name="s.tsv"):
    """`anchored` reads at the scheme's 5' ends (cycled), the rest uniform; fixed seed."""
    sites = sorted(read_bed_sites(BEDS / f"{scheme}.bed")[0])
    rng = random.Random(seed)
    pairs = {}
    for i in range(anchored):
        key = sites[i % len(sites)]
        pairs[key] = pairs.get(key, 0) + 1
    for _ in range(total - anchored):
        key = (rng.choice("+-"), rng.randint(1, genome))
        pairs[key] = pairs.get(key, 0) + 1
    return write_starts(tmp_path / name, [(s, p, c) for (s, p), c in pairs.items()])


def test_correct_scheme_is_ok_and_equivalent_schemes_tie(tmp_path):
    starts = starts_from_sites(tmp_path, "midnight", 1000, 5000)
    details, summary = run(tmp_path, starts, "midnight")
    assert summary["scheme_check"] == "ok"
    assert summary["reason"] == "declared_enriched"
    assert summary["best_scheme"] == "midnight/midnightv2"
    assert float(summary["declared_enrichment"]) > 20
    assert details["midnightv2"]["equivalent_to_declared"] == "true"
    assert details["articv4.1"]["equivalent_to_declared"] == "false"
    # unrelated schemes stay near chance (not above the bar)
    for other in ("articv4.1", "articv5.3.2", "qiagenv1", "swift", "xgen"):
        assert float(details[other]["enrichment"]) < 5


def test_declared_midnightv2_on_midnight_reads_is_ok(tmp_path):
    starts = starts_from_sites(tmp_path, "midnight", 1000, 5000)
    _, summary = run(tmp_path, starts, "midnightv2")
    assert summary["scheme_check"] == "ok"
    assert summary["declared_scheme"] == "midnightv2"


def test_wrong_declared_scheme_is_a_mismatch(tmp_path):
    starts = starts_from_sites(tmp_path, "midnight", 1000, 5000)
    _, summary = run(tmp_path, starts, "articv4.1")
    assert summary["scheme_check"] == "mismatch"
    assert summary["reason"] == "best_other_scheme"
    assert summary["best_scheme"] == "midnight/midnightv2"
    assert float(summary["declared_enrichment"]) < 5
    assert float(summary["best_enrichment"]) > 20


@pytest.mark.parametrize(
    "declared", ["articv3", "articv5.3.2", "eden", "qiagenv1", "swift", "xgen"]
)
def test_every_other_bundled_scheme_is_flagged_on_midnight_reads(tmp_path, declared):
    starts = starts_from_sites(tmp_path, "midnight", 1000, 5000)
    _, summary = run(tmp_path, starts, declared)
    assert summary["scheme_check"] == "mismatch"
    assert summary["best_scheme"] == "midnight/midnightv2"


def test_unanchored_reads_are_uninformative_never_an_error(tmp_path):
    starts = starts_from_sites(tmp_path, "midnight", 0, 20000)
    details, summary = run(tmp_path, starts, "midnight")
    assert summary["scheme_check"] == "uninformative"
    assert summary["reason"] == "no_scheme_enriched"
    assert summary["best_scheme"] == "NA"
    assert all(float(r["enrichment"]) < 3 for r in details.values())


def test_weak_anchoring_below_the_bar_is_uninformative(tmp_path):
    # ~1 % of reads at midnight sites: enrichment ~ 3 (window coverage is 0.49 %), below the bar of 5
    starts = starts_from_sites(tmp_path, "midnight", 50, 5000)
    _, summary = run(tmp_path, starts, "midnight")
    assert summary["scheme_check"] == "uninformative"


def test_few_reads_is_na(tmp_path):
    starts = starts_from_sites(tmp_path, "midnight", 200, 400)
    _, summary = run(tmp_path, starts, "midnight")
    assert (summary["scheme_check"], summary["reason"]) == ("NA", "too_few_reads")
    assert summary["declared_scheme"] == "midnight"


def test_min_reads_option(tmp_path):
    starts = starts_from_sites(tmp_path, "midnight", 200, 400)
    _, summary = run(tmp_path, starts, "midnight", "--min-reads", "300")
    assert summary["scheme_check"] == "ok"


def test_best_not_decisive_when_less_than_three_times_declared(tmp_path):
    beds = tmp_path / "beds"
    beds.mkdir()
    (beds / "a.bed").write_text("c\t99\t120\tx_F\t1\t+\n")  # site 100
    (beds / "b.bed").write_text("c\t499\t520\ty_F\t1\t+\n")  # site 500
    # W=0, G=1000, n=1000 -> expected 1 hit each. a: 4 hits (enrichment 4), b: 9 hits (enrichment 9)
    base = [("+", 100, 4), ("+", 500, 9), ("+", 700, 987)]
    starts = write_starts(tmp_path / "s.tsv", base)
    opts = (
        "--window",
        "0",
        "--genome-length",
        "1000",
        "--min-enrichment",
        "5",
        "--max-p-value",
        "1",
    )
    _, summary = run(tmp_path, starts, "a", *opts, beds_dir=beds)
    assert (summary["scheme_check"], summary["reason"]) == ("uninformative", "best_not_decisive")
    assert summary["best_scheme"] == "b"
    # with 20 hits at b the best is > 3x the declared one: mismatch
    starts = write_starts(tmp_path / "s2.tsv", [("+", 100, 4), ("+", 500, 20), ("+", 700, 976)])
    _, summary = run(tmp_path, starts, "a", *opts, beds_dir=beds)
    assert summary["scheme_check"] == "mismatch"


def test_p_value_guard_blocks_ratio_only_hits(tmp_path):
    beds = tmp_path / "beds"
    beds.mkdir()
    (beds / "tiny.bed").write_text("c\t99\t120\tx_F\t1\t+\n")
    # expected 0.5 hits (n=500, 1/1000 of the genome); 3 hits = enrichment 6 but p ~ 0.014 > 0.001
    starts = write_starts(tmp_path / "s.tsv", [("+", 100, 3), ("+", 700, 497)])
    details, summary = run(
        tmp_path, starts, "tiny", "--window", "0", "--genome-length", "1000", beds_dir=beds
    )
    assert float(details["tiny"]["enrichment"]) == pytest.approx(6.0)
    assert float(details["tiny"]["p_value"]) > 0.001
    assert summary["scheme_check"] == "uninformative"
    _, relaxed = run(
        tmp_path,
        starts,
        "tiny",
        "--window",
        "0",
        "--genome-length",
        "1000",
        "--max-p-value",
        "0.05",
        beds_dir=beds,
    )
    assert relaxed["scheme_check"] == "ok"


def test_declared_bed_path_is_scored_under_its_stem(tmp_path):
    custom = tmp_path / "my_scheme.bed"
    custom.write_text((BEDS / "midnight.bed").read_text())
    starts = starts_from_sites(tmp_path, "midnight", 1000, 5000)
    details, summary = run(tmp_path, starts, str(custom))
    assert summary["declared_scheme"] == "my_scheme"
    assert summary["scheme_check"] == "ok"
    assert details["my_scheme"]["equivalent_to_declared"] == "true"
    assert details["midnight"]["equivalent_to_declared"] == "true"


def test_unknown_declared_scheme_is_na(tmp_path, capsys):
    starts = starts_from_sites(tmp_path, "midnight", 1000, 5000)
    details, summary = run(tmp_path, starts, "no_such_scheme")
    assert details == {}
    assert (summary["scheme_check"], summary["reason"]) == ("NA", "declared_not_found")
    assert "note:" in capsys.readouterr().err


@pytest.mark.parametrize("content", [None, "", "garbage\n", "a\tb\n+\t5\n"])
def test_missing_or_bad_starts_is_na_exit_0(tmp_path, content):
    path = tmp_path / "s.tsv"
    if content is not None:
        path.write_text(content)
    _, summary = run(tmp_path, str(path), "midnight")
    assert (summary["scheme_check"], summary["reason"]) == ("NA", "missing_input")


def test_header_only_starts_is_too_few_reads(tmp_path):
    path = tmp_path / "s.tsv"
    path.write_text("strand\tfive_prime_pos\n")
    details, summary = run(tmp_path, str(path), "midnight")
    assert summary["reason"] == "too_few_reads"
    assert details["midnight"]["enrichment"] == "NA"


def test_summary_is_optional_and_stdout_works(tmp_path, capsys):
    starts = starts_from_sites(tmp_path, "midnight", 1000, 5000)
    assert main(["--starts", starts, "--declared", "midnight", "--beds-dir", str(BEDS)]) == 0
    lines = capsys.readouterr().out.splitlines()
    assert lines[0].split("\t") == DETAIL_COLUMNS
    assert len(lines) == 1 + len(list(BEDS.glob("*.bed")))


def test_unusable_arguments_exit_2(tmp_path):
    with pytest.raises(SystemExit) as exc:
        main(["--starts", "x", "--declared", "midnight", "--window", "-1"])
    assert exc.value.code == 2


# ---- real read starts ------------------------------------------------------------------------------------


def test_real_read_starts_fixture(tmp_path):
    """PURE_C (synthetic, midnight-shaped, anchored) mapped by bwa: tests/fixtures/qc/PURE_C.read_starts.tsv."""
    starts = str(FIX / "PURE_C.read_starts.tsv")
    details, summary = run(tmp_path, starts, "midnight", sample="PURE_C")
    assert summary["scheme_check"] == "ok"
    assert 35 < float(summary["declared_enrichment"]) < 50
    assert float(details["swift"]["enrichment"]) < 3
    _, wrong = run(tmp_path, starts, "articv4.1", sample="PURE_C")
    assert (wrong["scheme_check"], wrong["best_scheme"]) == ("mismatch", "midnight/midnightv2")


def test_real_unanchored_read_starts_fixture(tmp_path):
    """Randomly fragmented reference (no primer anchoring), mapped by bwa: uninformative."""
    starts = str(FIX / "UNANCHORED.read_starts.tsv")
    details, summary = run(tmp_path, starts, "midnight", sample="UNANCHORED")
    assert summary["scheme_check"] == "uninformative"
    assert max(float(r["enrichment"]) for r in details.values()) < 2
