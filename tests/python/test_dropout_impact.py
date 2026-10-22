"""Tests for bin/dropout_impact.py.

The tiny fixture in tests/fixtures/dropout_impact is hand-checkable (see the table in this docstring).
Four lineages (N = 4), so a mutation is informative unless all four carry it (4/4 = 1.0 >= 0.95):

    mutation  carried by   informative?
    C100T     A B C D      no (fixed)
    G150A     A B          yes
    A180C     A            yes
    T300G     B C          yes
    C350T     D            yes
    A450G     A B C        yes (0.75)
    G700A     C D          yes
    T800C     A B C D      no (fixed)
    G900A     B            yes
    (+ 6 deletion/insertion/other notations, skipped and counted)

Amplicons: tiny_1 1-200 {G150A, A180C}, tiny_2 150-400 {G150A, A180C, T300G, C350T}, tiny_3 401-600
{A450G}, tiny_4 601-1000 {G700A, G900A}: 7 informative sites in all.
"""

from pathlib import Path

import dropout_impact as di

FIX = Path(__file__).resolve().parents[1] / "fixtures" / "dropout_impact"
DATASET_DIR = FIX / "dataset"
AMPLICONS = FIX / "amplicon_coverage.tsv"
AGGREGATED = FIX / "freyja_aggregated.tsv"

HIER = FIX / "hier"


def run(tmp_path, *extra, dataset=DATASET_DIR, amplicons=AMPLICONS, aggregated=AGGREGATED):
    impact, summary = tmp_path / "impact.tsv", tmp_path / "summary.tsv"
    argv = [
        "--amplicons",
        str(amplicons),
        "--aggregated",
        str(aggregated),
        "--dataset",
        str(dataset),
        "--out",
        str(impact),
        "--summary-out",
        str(summary),
        *extra,
    ]
    code = di.main(argv)
    return code, impact, summary


def table(path):
    lines = Path(path).read_text().splitlines()
    header = lines[0].split("\t")
    return [dict(zip(header, line.split("\t"), strict=True)) for line in lines[1:]]


def by_key(rows, *keys):
    return {tuple(r[k] for k in keys): r for r in rows}


# ---- the hand-checked fixture -----------------------------------------------------------------
def test_fixture_matches_hand_computed_tables(tmp_path, capsys):
    code, impact, summary = run(tmp_path)
    assert code == 0
    assert impact.read_text() == (FIX / "expected_dropout_impact.tsv").read_text()
    assert summary.read_text() == (FIX / "expected_dropout_impact_summary.tsv").read_text()
    err = capsys.readouterr().err
    # deletions (C250-, AT270A, -11288-9), insertions (+T260TA, A280AT), other (foo)
    assert "skipped 6 mutation(s)" in err
    assert "(3 deletion, 2 insertion, 1 other)" in err


def test_samples_without_aggregate_row_get_na_and_false(tmp_path):
    _, impact, summary = run(tmp_path)
    row = by_key(table(impact), "sample", "amplicon")[("S3", "tiny_3")]
    assert row["n_sites_in_detected_lineages"] == "NA"
    assert row["affects_detected_lineage"] == "false"
    assert row["affected_lineages"] == "NA"
    s3 = by_key(table(summary), "sample")[("S3",)]
    assert s3["n_dropped_affecting_detected_lineages"] == "NA"
    assert s3["n_dropped_amplicons"] == "1"  # still appears


# ---- defining sites: the hand-made hierarchy in tests/fixtures/dropout_impact/hier ---------------
#
# Eight lineages (N = 8); C100T is carried by all of them (fixed, never informative):
#
#   O    root                          C100T T900C
#   R    root                          C100T G150A A200C T500G
#   R.1  parent R (alias AA.1)         R + G300A
#   R.1.1 parent AA.1 (by alias)       R.1 + A350C T600G
#   R.2  parent R                      R + C320T
#   X    recombinant root              R.1 + R.2 + G800A      (recombinant_parents R.1,R.2)
#   Z    parent GONE (no barcode)      R.2 + A700T            (GONE -> R.2, so the nearest barcode is R.2)
#   W    parent NOPE (not in the yml)  C100T G150A T450A      (a root)
#
# Defining informative sites: O {T900C}; R {G150A A200C T500G}; R.1 {G300A}; R.1.1 {A350C T600G};
# R.2 {C320T}; X {G800A} (strict: G150A A200C T500G G300A C320T G800A); Z {A700T}; W {G150A T450A}.
# Amplicons: a1 1-250 (G150A A200C), a2 251-400 (G300A C320T A350C), a3 401-550 (T450A T500G),
# a4 551-750 (T600G A700T), a5 751-950 (G800A T900C): 11 informative sites in all.
def run_hier(tmp_path, *extra, **kw):
    kw.setdefault("dataset", HIER)
    kw.setdefault("amplicons", HIER / "amplicon_coverage.tsv")
    kw.setdefault("aggregated", HIER / "freyja_aggregated.tsv")
    return run(tmp_path, *extra, **kw)


def hier_tables(tmp_path, *extra):
    code, impact, summary = run_hier(tmp_path, *extra)
    assert code == 0
    return by_key(table(impact), "sample", "amplicon"), by_key(table(summary), "sample")


def test_hier_fixture_matches_golden_files(tmp_path):
    _, impact, summary = run_hier(tmp_path)
    assert impact.read_text() == (HIER / "expected_dropout_impact.tsv").read_text()
    assert summary.read_text() == (HIER / "expected_dropout_impact_summary.tsv").read_text()


def test_inherited_sites_alone_do_not_flag_a_dropout(tmp_path):
    rows, summary = hier_tables(tmp_path)
    # H8: R.1.1 (reported by its alias AA.1.1) .50 and R.2 .40, a3 (T450A, T500G) dropped:
    # R.1.1 and R.2 both carry T500G (inherited from R), so it is carried but not defining and not differential
    r = rows[("H8", "a3")]
    assert r["n_informative_sites"] == "2"
    assert r["n_sites_in_detected_lineages"] == "1"
    assert (r["n_defining_sites"], r["n_differential_sites"]) == ("0", "0")
    assert (r["affects_detected_lineage"], r["affected_lineages"]) == ("false", "NA")
    # ... while the flagged neighbour a2 is reported under the name the aggregate used
    assert rows[("H8", "a2")]["affected_lineages"] == "AA.1.1:1,R.2:1"
    s = summary[
        ("H8",)
    ]  # scheme a2-a4 only (251-750): defining A350C T600G C320T = 3; lost (251-550) = 2
    assert (s["defining_sites_total"], s["defining_sites_lost"], s["pct_defining_sites_lost"]) == (
        "3",
        "2",
        "66.67",
    )
    assert s["informative_sites_total"] == "7"
    assert s["n_dropped_affecting_detected_lineages"] == "1"
    # 0.0099 < 0.01: O is not detected in H8, so its T900C (outside the scheme anyway) cannot appear


def test_contig_selects_one_contig_of_a_multi_contig_table(tmp_path):
    rows = AMPLICONS.read_text().splitlines()
    multi = tmp_path / "multi.tsv"
    # the fixture's rows on contig segA, plus extra dropped rows on segB that must not count
    multi.write_text(
        rows[0]
        + "\tcontig\n"
        + "".join(f"{r}\tsegA\n" for r in rows[1:])
        + "".join(f"{r}\tsegB\n" for r in rows[1:] if r.endswith("true"))
    )
    code, impact, summary = run(tmp_path, "--contig", "segA", amplicons=multi)
    assert code == 0
    assert impact.read_text() == (FIX / "expected_dropout_impact.tsv").read_text()
    assert summary.read_text() == (FIX / "expected_dropout_impact_summary.tsv").read_text()
    assert run(tmp_path, amplicons=multi)[0] == 2  # several contigs: --contig is required
