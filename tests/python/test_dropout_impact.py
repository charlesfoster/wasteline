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

import json
import random
import time
from pathlib import Path

import dropout_impact as di
import pytest

FIX = Path(__file__).resolve().parents[1] / "fixtures" / "dropout_impact"
DATASET_DIR = FIX / "dataset"
BARCODES_CSV = FIX / "barcodes.csv"
AMPLICONS = FIX / "amplicon_coverage.tsv"
AGGREGATED = FIX / "freyja_aggregated.tsv"

HIER = FIX / "hier"

IMPACT_HEADER = "\t".join(di.IMPACT_COLUMNS)
SUMMARY_HEADER = "\t".join(di.SUMMARY_COLUMNS)


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
@pytest.mark.parametrize(
    "dataset", [DATASET_DIR, DATASET_DIR / "lineage_mutations.json", BARCODES_CSV]
)
def test_fixture_matches_hand_computed_tables(tmp_path, capsys, dataset):
    code, impact, summary = run(tmp_path, dataset=dataset)
    assert code == 0
    assert impact.read_text() == (FIX / "expected_dropout_impact.tsv").read_text()
    assert summary.read_text() == (FIX / "expected_dropout_impact_summary.tsv").read_text()
    err = capsys.readouterr().err
    # deletions (C250-, AT270A, -11288-9), insertions (+T260TA, A280AT), other (foo)
    assert "skipped 6 mutation(s)" in err
    assert "(3 deletion, 2 insertion, 1 other)" in err


def test_hand_checked_rows(tmp_path):
    _, impact, summary = run(tmp_path)
    rows = by_key(table(impact), "sample", "amplicon")
    # tiny_2 = 150-400: G150A (A,B), A180C (A), T300G (B,C), C350T (D); detected = A, B (C is at 0.009)
    r = rows[("S1", "tiny_2")]
    assert (r["n_informative_sites"], r["n_sites_in_detected_lineages"]) == ("4", "3")
    assert r["affected_lineages"] == "A:2,B:2"  # tie -> by name
    # tiny_4 = 601-1000: G700A is carried only by C and D (not detected); G900A by B
    r = rows[("S1", "tiny_4")]
    assert (r["n_informative_sites"], r["n_sites_in_detected_lineages"]) == ("2", "1")
    assert r["affected_lineages"] == "B:1"
    # fixed mutations never count: tiny_1 spans C100T (fixed) too
    assert rows[("S5", "tiny_1")]["n_informative_sites"] == "2"
    # no dropped rows for samples without dropped amplicons, nor for S6 (not in the amplicon table)
    assert not [k for k in rows if k[0] in {"S2", "S6"}]
    assert "S6" not in {r["sample"] for r in table(summary)}


def test_overlapping_dropped_amplicons_counted_once_in_summary(tmp_path):
    _, impact, summary = run(tmp_path)
    rows = by_key(table(impact), "sample", "amplicon")
    s5 = by_key(table(summary), "sample")[("S5",)]
    # tiny_1 (2) + tiny_2 (4) = 6 per-amplicon, but 150-200 is shared: union holds 4 sites
    assert (
        int(rows[("S5", "tiny_1")]["n_informative_sites"])
        + int(rows[("S5", "tiny_2")]["n_informative_sites"])
        == 6
    )
    assert s5["informative_sites_in_dropped"] == "4"
    assert s5["pct_informative_sites_lost"] == "57.14"


def test_total_is_per_sample_scheme_union(tmp_path):
    _, _, summary = run(tmp_path)
    s = by_key(table(summary), "sample")
    assert s[("S7",)]["informative_sites_total"] == "3"  # scheme of S7 only spans 401-1000
    assert s[("S1",)]["informative_sites_total"] == "7"


def test_samples_without_aggregate_row_get_na_and_false(tmp_path):
    _, impact, summary = run(tmp_path)
    row = by_key(table(impact), "sample", "amplicon")[("S3", "tiny_3")]
    assert row["n_sites_in_detected_lineages"] == "NA"
    assert row["affects_detected_lineage"] == "false"
    assert row["affected_lineages"] == "NA"
    s3 = by_key(table(summary), "sample")[("S3",)]
    assert s3["n_dropped_affecting_detected_lineages"] == "NA"
    assert s3["n_dropped_amplicons"] == "1"  # still appears


def test_empty_lineage_list_means_zero_not_na(tmp_path):
    _, impact, summary = run(tmp_path)
    row = by_key(table(impact), "sample", "amplicon")[("S4", "tiny_4")]
    assert row["n_sites_in_detected_lineages"] == "0"
    assert row["affects_detected_lineage"] == "false"
    assert by_key(table(summary), "sample")[("S4",)]["n_dropped_affecting_detected_lineages"] == "0"


def test_min_abundance_changes_detected_set(tmp_path):
    _, impact, _ = run(tmp_path, "--min-abundance", "0.005")
    row = by_key(table(impact), "sample", "amplicon")[("S1", "tiny_4")]
    # lineage C (0.009) now counts and carries G700A
    assert row["n_sites_in_detected_lineages"] == "2"
    assert row["affected_lineages"] == "B:1,C:1"


def test_min_abundance_is_inclusive(tmp_path):
    _, impact, _ = run(tmp_path, "--min-abundance", "0.009")
    assert by_key(table(impact), "sample", "amplicon")[("S1", "tiny_4")]["affected_lineages"] == (
        "B:1,C:1"
    )


def test_informative_max_fraction(tmp_path):
    # 3/4 = 0.75 is not < 0.75: A450G drops out
    _, impact, summary = run(tmp_path, "--informative-max-fraction", "0.75")
    row = by_key(table(impact), "sample", "amplicon")[("S3", "tiny_3")]
    assert row["n_informative_sites"] == "0"
    assert by_key(table(summary), "sample")[("S1",)]["informative_sites_total"] == "6"
    # 0.7500001 keeps it
    _, impact, _ = run(tmp_path, "--informative-max-fraction", "0.76")
    assert (
        by_key(table(impact), "sample", "amplicon")[("S3", "tiny_3")]["n_informative_sites"] == "1"
    )


def test_informative_fraction_boundary_95_percent(tmp_path):
    # 20 lineages: a mutation on 19 of them is exactly 0.95 -> not informative; on 18 -> informative
    data = {f"L{i}": [] for i in range(20)}
    for i in range(19):
        data[f"L{i}"].append("C100T")
    for i in range(18):
        data[f"L{i}"].append("C200T")
    data["L0"].append("C300T")
    ds = tmp_path / "ds.json"
    ds.write_text(json.dumps(data))
    lineages = di.load_json_dataset(ds, di.SkipCounter())
    dataset = di.Dataset(lineages, 0.95)
    assert dataset.sorted_sites == [(200, "T"), (300, "T")]


def test_affected_lineages_capped_at_five_and_sorted(tmp_path):
    data = {f"L{i}": ["C100T"] for i in range(7)}
    data["L0"] += ["C101T", "C102T"]  # L0 carries 3 sites in the span, the others 1
    data["Other"] = []  # so that 100T is not fixed
    ds = tmp_path / "lineage_mutations.json"
    ds.write_text(json.dumps(data))
    amp = tmp_path / "amp.tsv"
    amp.write_text(
        "sample\tamplicon\tstart\tend\tmean_depth\tmin_depth\tpct_ge_min_depth\tdropped\n"
        "X\ta1\t90\t110\t0.00\t0\t0.00\ttrue\n"
    )
    agg = tmp_path / "agg.tsv"
    names = " ".join(f"L{i}" for i in (6, 5, 4, 3, 2, 1, 0))
    agg.write_text(
        f"\tsummarized\tlineages\tabundances\tresid\tcoverage\n"
        f"X\t[]\t{names}\t{' '.join(['0.1'] * 7)}\t0\t99\n"
    )
    code, impact, _ = run(tmp_path, dataset=ds, amplicons=amp, aggregated=agg)
    assert code == 0
    (row,) = table(impact)
    assert row["affected_lineages"] == "L0:3,L1:1,L2:1,L3:1,L4:1"
    assert row["n_sites_in_detected_lineages"] == "3"


# ---- zero dropped / empty inputs ---------------------------------------------------------------
def test_no_dropped_amplicons_header_only_impact_and_zero_summary(tmp_path):
    amp = tmp_path / "amp.tsv"
    amp.write_text(
        "sample\tamplicon\tstart\tend\tmean_depth\tmin_depth\tpct_ge_min_depth\tdropped\n"
        "S1\ttiny_1\t1\t200\t150.00\t90\t100.00\tfalse\n"
        "S1\ttiny_2\t150\t400\t140.00\t80\t100.00\tfalse\n"
    )
    code, impact, summary = run(tmp_path, amplicons=amp)
    assert code == 0
    assert impact.read_text() == IMPACT_HEADER + "\n"
    assert summary.read_text() == f"{SUMMARY_HEADER}\nS1\t0\t4\t0\t0.00\t0\t3\t0\t0.00\n"


def test_empty_amplicon_table_gives_header_only_outputs(tmp_path):
    amp = tmp_path / "amp.tsv"
    amp.write_text(
        "sample\tamplicon\tstart\tend\tmean_depth\tmin_depth\tpct_ge_min_depth\tdropped\n"
    )
    code, impact, summary = run(tmp_path, amplicons=amp)
    assert code == 0
    assert impact.read_text() == IMPACT_HEADER + "\n"
    assert summary.read_text() == SUMMARY_HEADER + "\n"


def test_zero_informative_sites_gives_na_percentage(tmp_path):
    ds = tmp_path / "lineage_mutations.json"
    ds.write_text(json.dumps({"A": ["C100T"], "B": ["C100T"]}))  # every lineage carries it: fixed
    code, impact, summary = run(tmp_path, dataset=ds)
    assert code == 0
    s1 = by_key(table(summary), "sample")[("S1",)]
    assert (s1["informative_sites_total"], s1["pct_informative_sites_lost"]) == ("0", "NA")
    assert all(r["n_informative_sites"] == "0" for r in table(impact))


def test_aggregate_without_lineages_column_is_an_error(tmp_path, capsys):
    agg = tmp_path / "agg.tsv"
    agg.write_text("\tsummarized\nS1\t[]\n")
    code, _, _ = run(tmp_path, aggregated=agg)
    assert code == 2
    assert "lineages" in capsys.readouterr().err


def test_amplicon_table_missing_columns_is_an_error(tmp_path, capsys):
    amp = tmp_path / "amp.tsv"
    amp.write_text("sample\tamplicon\nS1\ttiny_1\n")
    code, _, _ = run(tmp_path, amplicons=amp)
    assert code == 2
    assert "missing column" in capsys.readouterr().err


def test_aggregate_label_with_file_suffix_is_resolved(tmp_path):
    agg = tmp_path / "agg.tsv"
    agg.write_text(
        "\tsummarized\tlineages\tabundances\tresid\tcoverage\nS1.demix\t[]\tA B\t0.6 0.4\t0\t99\n"
    )
    _, impact, _ = run(tmp_path, aggregated=agg)
    row = by_key(table(impact), "sample", "amplicon")[("S1", "tiny_2")]
    assert row["n_sites_in_detected_lineages"] == "3"


# ---- dataset forms -----------------------------------------------------------------------------
def test_feather_only_dataset_directory_is_a_clear_error(tmp_path, capsys):
    ds = tmp_path / "ds"
    ds.mkdir()
    (ds / "usher_barcodes.feather").write_bytes(b"ARROW1")
    code, _, _ = run(tmp_path, dataset=ds)
    assert code == 2
    err = capsys.readouterr().err
    assert "feather" in err and "lineage_mutations.json" in err


def test_feather_file_is_a_clear_error(tmp_path, capsys):
    f = tmp_path / "barcodes.feather"
    f.write_bytes(b"ARROW1")
    code, _, _ = run(tmp_path, dataset=f)
    assert code == 2
    assert "dataset directory" in capsys.readouterr().err


def test_empty_json_with_feather_is_an_error(tmp_path, capsys):
    ds = tmp_path / "ds"
    ds.mkdir()
    (ds / "lineage_mutations.json").write_text("")  # the stub datasets of the test profile
    (ds / "usher_barcodes.feather").write_bytes(b"ARROW1")
    code, _, _ = run(tmp_path, dataset=ds)
    assert code == 2
    assert "feather" in capsys.readouterr().err


def test_dataset_directory_without_json_falls_back_to_barcodes_csv(tmp_path):
    ds = tmp_path / "ds"
    ds.mkdir()
    (ds / "usher_barcodes.csv").write_text(BARCODES_CSV.read_text())
    (ds / "usher_barcodes.feather").write_bytes(b"ARROW1")
    code, impact, _ = run(tmp_path, dataset=ds)
    assert code == 0
    assert impact.read_text() == (FIX / "expected_dropout_impact.tsv").read_text()


def test_json_wins_over_csv_in_a_directory(tmp_path):
    ds = tmp_path / "ds"
    ds.mkdir()
    (ds / "lineage_mutations.json").write_text((DATASET_DIR / "lineage_mutations.json").read_text())
    (ds / "usher_barcodes.csv").write_text("garbage")
    code, impact, _ = run(tmp_path, dataset=ds)
    assert code == 0
    assert impact.read_text() == (FIX / "expected_dropout_impact.tsv").read_text()


@pytest.mark.parametrize("content", ["{not json", "[1, 2]", ""])
def test_malformed_json_is_an_error_not_a_crash(tmp_path, capsys, content):
    f = tmp_path / "lineage_mutations.json"
    f.write_text(content)
    code, _, _ = run(tmp_path, dataset=f)
    assert code == 2
    assert capsys.readouterr().err.startswith("error:")


def test_missing_dataset_is_an_error(tmp_path, capsys):
    code, _, _ = run(tmp_path, dataset=tmp_path / "nope")
    assert code == 2
    assert "not found" in capsys.readouterr().err


def test_odd_json_values_do_not_crash(tmp_path, capsys):
    data = {"A": ["C100T", None, 5, "x"], "B": None, "C": {"G200A": 1}, "D": "C100T,G300A"}
    f = tmp_path / "lineage_mutations.json"
    f.write_text(json.dumps(data))
    code, _, _ = run(tmp_path, dataset=f)
    assert code == 0
    assert "skipped" in capsys.readouterr().err


# ---- notation ------------------------------------------------------------------------------------
@pytest.mark.parametrize(
    ("text", "expected"),
    [
        ("C241T", (241, "T")),
        ("a241t", (241, "T")),
        ("A23403G", (23403, "G")),
    ],
)
def test_parse_mutation_substitutions(text, expected):
    skipped = di.SkipCounter()
    assert di.parse_mutation(text, skipped) == expected
    assert skipped.total == 0


@pytest.mark.parametrize(
    ("text", "kind"),
    [
        ("C11288-", "deletion"),
        ("-11288-9", "deletion"),
        ("ATG100A", "deletion"),
        ("+T123TA", "insertion"),
        ("A12AT", "insertion"),
        ("ins_123", "other"),
        ("", "other"),
    ],
)
def test_parse_mutation_skips_and_classifies_indels(text, kind):
    skipped = di.SkipCounter()
    assert di.parse_mutation(text, skipped) is None
    assert getattr(skipped, kind) == 1
    assert skipped.total == 1


# ---- tolerant list parsing ---------------------------------------------------------------------
def test_list_parsers_cope_with_numpy_and_python_reprs():
    assert di.parse_name_list("BA.2 AY.4") == ["BA.2", "AY.4"]
    assert di.parse_name_list("['BA.2' 'AY.4']") == ["BA.2", "AY.4"]  # numpy: no commas
    assert di.parse_name_list("['A', \"B.1\"]") == ["A", "B.1"]
    assert di.parse_name_list("") == [] and di.parse_name_list("NA") == []
    assert di.parse_float_list("0.7 0.30000") == [0.7, 0.3]
    assert di.parse_float_list("[0.5\n 0.25 1.e-2]") == [0.5, 0.25, 0.01]
    assert di.parse_float_list("[nan 0.5]") == [0.0, 0.5]
    assert di.parse_float_list("") == []


def test_mismatched_lengths_warn_and_do_not_crash(tmp_path, capsys):
    agg = tmp_path / "agg.tsv"
    agg.write_text(
        "\tsummarized\tlineages\tabundances\tresid\tcoverage\nS1\t[]\tA B C\t0.6 0.4\t0\t99\n"
    )
    code, _, _ = run(tmp_path, aggregated=agg)
    assert code == 0
    assert "3 lineages but 2 abundances" in capsys.readouterr().err


def test_merge_intervals_merges_overlaps_and_adjacent():
    assert di.merge_intervals([(10, 20), (1, 5), (6, 8), (15, 30)]) == [[1, 8], [10, 30]]
    assert di.merge_intervals([]) == []


# ---- performance ---------------------------------------------------------------------------------
def test_huge_dataset_runs_fast(tmp_path):
    rng = random.Random(7)
    pool = [
        f"{rng.choice('ACGT')}{pos}{rng.choice('ACGT')}"
        for pos in rng.sample(range(1, 29903), 3000)
    ]
    data = {
        f"L{i}": rng.sample(pool, 85) + ["C100-", "+T5TA"] for i in range(6000)
    }  # ~510k entries
    ds = tmp_path / "lineage_mutations.json"
    ds.write_text(json.dumps(data))
    amp = tmp_path / "amp.tsv"
    amp.write_text(
        "sample\tamplicon\tstart\tend\tmean_depth\tmin_depth\tpct_ge_min_depth\tdropped\n"
        + "".join(
            f"S{s}\ta{i}\t{1 + i * 1000}\t{1100 + i * 1000}\t1.0\t0\t0.0\ttrue\n"
            for s in range(10)
            for i in range(29)
        )
    )
    agg = tmp_path / "agg.tsv"
    lines = ["\tsummarized\tlineages\tabundances\tresid\tcoverage"]
    lines += [f"S{s}\t[]\tL1 L2 L3 L4 L5\t0.4 0.3 0.2 0.05 0.05\t0\t99" for s in range(10)]
    agg.write_text("\n".join(lines) + "\n")
    start = time.monotonic()
    code, impact, summary = run(tmp_path, dataset=ds, amplicons=amp, aggregated=agg)
    elapsed = time.monotonic() - start
    assert code == 0
    assert elapsed < 10
    assert len(table(impact)) == 290
    assert len(table(summary)) == 10


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


def test_defining_sites_are_barcode_minus_parent_barcode(tmp_path):
    lineages = di.load_json_dataset(HIER / "lineage_mutations.json", di.SkipCounter())
    entries = di.parse_lineages_yml(HIER / "lineages.yml")
    ds = di.Dataset(lineages, 0.95, di.Hierarchy(lineages, entries))

    def defining(name):
        return {f"{a}{b}" for a, b in ds.defining_of(name)}

    assert ds.informative >= {(150, "A")} and (100, "T") not in ds.informative  # C100T is fixed
    assert defining("O") == {"900C"}
    assert defining("R") == {"150A", "200C", "500G"}  # a root: everything informative
    assert defining("R.1") == {"300A"}
    assert defining("R.1.1") == {"350C", "600G"}  # parent given by ALIAS (AA.1)
    assert defining("R.2") == {"320T"}
    assert defining("X") == {"800A"}  # recombinant root: minus R.1 and R.2
    assert defining("Z") == {"700T"}  # parent GONE has no barcode: nearest ancestor R.2
    assert defining("W") == {"150A", "450A"}  # parent NOPE is nowhere: a root
    assert defining("not-a-lineage") == set()


def test_strict_roots_counts_all_mutations_of_a_recombinant_root(tmp_path):
    lineages = di.load_json_dataset(HIER / "lineage_mutations.json", di.SkipCounter())
    entries = di.parse_lineages_yml(HIER / "lineages.yml")
    ds = di.Dataset(lineages, 0.95, di.Hierarchy(lineages, entries), strict_roots=True)
    assert [f"{a}{b}" for a, b in ds.defining_of("X")] == [
        "150A",
        "200C",
        "300A",
        "320T",
        "500G",
        "800A",
    ]
    # H3 = X alone, a1 (G150A, A200C) dropped: inherited sites do not count by default, do when strict
    rows, _ = hier_tables(tmp_path)
    assert (
        rows[("H3", "a1")]["n_defining_sites"],
        rows[("H3", "a1")]["affects_detected_lineage"],
    ) == (
        "0",
        "false",
    )
    rows, _ = hier_tables(tmp_path, "--strict-roots")
    row = rows[("H3", "a1")]
    assert (row["n_defining_sites"], row["affects_detected_lineage"], row["affected_lineages"]) == (
        "2",
        "true",
        "X:2",
    )


def test_hier_two_clades_defining_and_differential(tmp_path):
    rows, summary = hier_tables(tmp_path)
    # H1: R.1.1 .40, R.2 .30, R .20, Z .05, O .05; a2 and a4 dropped
    r = rows[("H1", "a2")]  # G300A (R.1.1), C320T (R.2, Z), A350C (R.1.1)
    assert (r["n_informative_sites"], r["n_sites_in_detected_lineages"]) == ("3", "3")
    assert r["n_defining_sites"] == "2"  # A350C (R.1.1) and C320T (R.2); G300A is inherited
    assert r["affected_lineages"] == "R.1.1:1,R.2:1"  # tie -> by name
    assert r["n_differential_sites"] == "3"  # G300A, A350C, C320T all differ between the five
    assert r["affects_detected_lineage"] == "true"
    r = rows[("H1", "a4")]  # T600G (R.1.1) and A700T (Z)
    assert (r["n_defining_sites"], r["n_differential_sites"]) == ("2", "2")
    assert r["affected_lineages"] == "R.1.1:1,Z:1"
    s = summary[("H1",)]
    # defining sites of the five (R.1.1 A350C T600G; R.2 C320T; R G150A A200C T500G; Z A700T; O T900C) = 8,
    # of which A350C C320T T600G A700T lie in a2 (251-400) or a4 (551-750)
    assert (s["defining_sites_total"], s["defining_sites_lost"], s["pct_defining_sites_lost"]) == (
        "8",
        "4",
        "50.00",
    )
    assert s["n_dropped_affecting_detected_lineages"] == "2"
    assert (s["informative_sites_in_dropped"], s["pct_informative_sites_lost"]) == ("5", "45.45")


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


def test_recombinant_root_and_missing_parent_and_nearest_ancestor(tmp_path):
    rows, summary = hier_tables(tmp_path)
    # H7: Z .70 + R.2 .30, a4 dropped: Z's own A700T is defining (parent chain GONE -> R.2)
    r = rows[("H7", "a4")]
    assert (r["n_defining_sites"], r["affected_lineages"]) == ("1", "Z:1")
    assert r["n_differential_sites"] == "1"
    s = summary[("H7",)]
    assert (s["defining_sites_total"], s["defining_sites_lost"], s["pct_defining_sites_lost"]) == (
        "2",
        "1",
        "50.00",
    )
    # H4: W (parent NOPE: a root) .70 + O .30, a3 dropped: W's T450A is defining
    r = rows[("H4", "a3")]
    assert (r["n_defining_sites"], r["affected_lineages"], r["n_differential_sites"]) == (
        "1",
        "W:1",
        "1",
    )


def test_min_affect_abundance_boundary(tmp_path):
    # H2: R .95 + Z .05; a4 hides Z's A700T. Z is exactly at the default 0.05: flagged (>=)
    rows, summary = hier_tables(tmp_path)
    r = rows[("H2", "a4")]
    assert (r["affects_detected_lineage"], r["affected_lineages"], r["n_defining_sites"]) == (
        "true",
        "Z:1",
        "1",
    )
    assert summary[("H2",)]["n_dropped_affecting_detected_lineages"] == "1"
    # just above: Z is detected, listed and counted, but no longer raises the flag
    rows, summary = hier_tables(tmp_path, "--min-affect-abundance", "0.0501")
    r = rows[("H2", "a4")]
    assert (r["affects_detected_lineage"], r["affected_lineages"], r["n_defining_sites"]) == (
        "false",
        "Z:1",
        "1",
    )
    s = summary[("H2",)]
    assert s["n_dropped_affecting_detected_lineages"] == "0"
    assert (s["defining_sites_total"], s["defining_sites_lost"]) == ("4", "1")  # unaffected by it
    # 0 flags everything detected; the --min-abundance floor still decides who is detected
    rows, _ = hier_tables(tmp_path, "--min-affect-abundance", "0", "--min-abundance", "0.06")
    r = rows[("H2", "a4")]  # Z (0.05) is no longer detected at all
    assert (r["affects_detected_lineage"], r["n_defining_sites"], r["affected_lineages"]) == (
        "false",
        "0",
        "NA",
    )


def test_min_affect_abundance_is_validated(tmp_path):
    with pytest.raises(SystemExit):
        run_hier(tmp_path, "--min-affect-abundance", "1.5")


def test_differential_sites_need_two_detected_lineages(tmp_path):
    rows, _ = hier_tables(tmp_path)
    assert (
        rows[("H3", "a1")]["n_differential_sites"] == "0"
    )  # one lineage (X): nothing to differ from
    assert rows[("H5", "a1")]["n_differential_sites"] == "0"  # the only detected lineage is unknown


def test_differential_does_not_raise_the_flag(tmp_path):
    # n_differential_sites > 0 with n_defining_sites == 0 and no flag must be possible
    data = {"P": ["A100C", "A110C"], "Q": ["A100C", "A120C"], "P.1": ["A100C", "A110C", "A130C"]}
    ds = tmp_path / "ds"
    ds.mkdir()
    (ds / "lineage_mutations.json").write_text(json.dumps(data))
    (ds / "lineages.yml").write_text(
        "- name: P\n  alias: P\n- name: P.1\n  alias: P.1\n  parent: P\n- name: Q\n  alias: Q\n"
    )
    amp = tmp_path / "amp.tsv"
    amp.write_text(
        "sample\tamplicon\tstart\tend\tmean_depth\tmin_depth\tpct_ge_min_depth\tdropped\n"
        "X\tm1\t105\t115\t0\t0\t0\ttrue\n"
    )
    agg = tmp_path / "agg.tsv"
    agg.write_text("\tsummarized\tlineages\tabundances\nX\t[]\tP.1 Q\t0.5 0.5\n")
    # m1 hides A110C: carried by P.1 (inherited from P) and not by Q -> differential, but not defining
    code, impact, _ = run(tmp_path, dataset=ds, amplicons=amp, aggregated=agg)
    assert code == 0
    (row,) = table(impact)
    assert (row["n_differential_sites"], row["n_defining_sites"]) == ("1", "0")
    assert row["affects_detected_lineage"] == "false"


def test_flat_fixture_new_columns(tmp_path):
    # hand-checked on the flat fixture (no lineages.yml: every lineage is a root, a warning says so)
    _, impact, summary = run(tmp_path)
    rows = by_key(table(impact), "sample", "amplicon")
    # S1 tiny_2: A {G150A, A180C}, B {G150A, T300G} (informative sites only): union 3; differ: A180C, T300G
    assert (
        rows[("S1", "tiny_2")]["n_defining_sites"],
        rows[("S1", "tiny_2")]["n_differential_sites"],
    ) == (
        "3",
        "2",
    )
    # a sample without an aggregate row: NA
    assert (
        rows[("S3", "tiny_3")]["n_defining_sites"],
        rows[("S3", "tiny_3")]["n_differential_sites"],
    ) == (
        "NA",
        "NA",
    )
    s = by_key(table(summary), "sample")
    assert (s[("S1",)]["defining_sites_total"], s[("S1",)]["defining_sites_lost"]) == ("5", "4")
    assert s[("S1",)]["pct_defining_sites_lost"] == "80.00"
    assert (s[("S3",)]["defining_sites_total"], s[("S3",)]["pct_defining_sites_lost"]) == (
        "NA",
        "NA",
    )
    # S4: aggregate row with no lineages -> 0 defining sites in total, percentage NA (0/0)
    assert (s[("S4",)]["defining_sites_total"], s[("S4",)]["pct_defining_sites_lost"]) == (
        "0",
        "NA",
    )


def test_old_columns_keep_name_and_position():
    assert di.IMPACT_COLUMNS[:9] == [
        "sample",
        "amplicon",
        "start",
        "end",
        "mean_depth",
        "n_informative_sites",
        "n_sites_in_detected_lineages",
        "affects_detected_lineage",
        "affected_lineages",
    ]
    assert di.IMPACT_COLUMNS[9:] == ["n_defining_sites", "n_differential_sites"]
    assert di.SUMMARY_COLUMNS[:6] == [
        "sample",
        "n_dropped_amplicons",
        "informative_sites_total",
        "informative_sites_in_dropped",
        "pct_informative_sites_lost",
        "n_dropped_affecting_detected_lineages",
    ]
    assert di.SUMMARY_COLUMNS[6:] == [
        "defining_sites_total",
        "defining_sites_lost",
        "pct_defining_sites_lost",
    ]


def test_old_consumers_can_read_the_first_columns_positionally(tmp_path):
    _, impact, summary = run(tmp_path)
    old_impact = [line.split("\t")[:9] for line in impact.read_text().splitlines()]
    old_summary = [line.split("\t")[:6] for line in summary.read_text().splitlines()]
    assert old_impact[1] == ["S1", "tiny_2", "150", "400", "3.10", "4", "3", "true", "A:2,B:2"]
    assert old_summary[1] == ["S1", "2", "7", "6", "85.71", "2"]


def test_lineages_yml_is_found_next_to_a_dataset_file_and_via_option(tmp_path, capsys):
    # json file given directly: lineages.yml next to it is used (no "guessed" warning)
    code, impact, _ = run_hier(tmp_path, dataset=HIER / "lineage_mutations.json")
    assert code == 0
    assert impact.read_text() == (HIER / "expected_dropout_impact.tsv").read_text()
    assert "no lineages.yml" not in capsys.readouterr().err
    # a barcodes csv elsewhere + explicit yml
    csv_path = tmp_path / "barcodes.csv"
    muts = json.loads((HIER / "lineage_mutations.json").read_text())
    cols = sorted({m for v in muts.values() for m in v})
    csv_path.write_text(
        ","
        + ",".join(cols)
        + "\n"
        + "".join(
            f"{n}," + ",".join("1" if c in v else "0" for c in cols) + "\n" for n, v in muts.items()
        )
    )
    code, impact, _ = run_hier(
        tmp_path, "--lineages-yml", str(HIER / "lineages.yml"), dataset=csv_path
    )
    assert code == 0
    assert impact.read_text() == (HIER / "expected_dropout_impact.tsv").read_text()
    assert run_hier(tmp_path, "--lineages-yml", str(tmp_path / "missing.yml"))[0] == 2


def test_without_lineages_yml_parents_are_guessed_from_names_and_warned(tmp_path, capsys):
    data = {
        "B": ["C100T", "C101T"],
        "B.1": ["C100T", "C101T", "C102T"],
        "B.1.2": ["C100T", "C101T", "C102T", "C103T"],
        "Other": [],
    }
    f = tmp_path / "lineage_mutations.json"
    f.write_text(json.dumps(data))
    lineages = di.load_json_dataset(f, di.SkipCounter())
    ds = di.Dataset(lineages, 0.95, di.Hierarchy(lineages, {}))
    assert ds.defining_of("B.1.2") == [(103, "T")]  # B.1.2 -> B.1
    assert ds.defining_of("B.1") == [(102, "T")]
    assert ds.defining_of("B") == [(100, "T"), (101, "T")]
    code, _, _ = run(tmp_path, dataset=f)
    assert code == 0
    assert "no lineages.yml found" in capsys.readouterr().err


def test_yml_reader_handles_the_freyja_shape(tmp_path):
    f = tmp_path / "lineages.yml"
    f.write_text(
        "- name: A\n  alias: A\n  children:\n      - A\n      - A.1\n"
        "- name: A.1\n  alias: 'X.1'\n  children:\n      - A.1\n  parent: A\n"
        "- name: XB\n  alias: XB\n  children:\n      - XB\n  recombinant_parents: A.1*,B\n"
    )
    entries = di.parse_lineages_yml(f)
    assert entries["A"] == {"alias": "A"}
    assert entries["A.1"] == {"alias": "X.1", "parent": "A"}
    assert entries["XB"]["recombinant_parents"] == "A.1*,B"
    h = di.Hierarchy({"A", "A.1", "XB"}, entries)
    assert h.parent("A.1") == "A" and h.parent("A") is None
    assert h.canonical("x.1") == "A.1"  # alias, any case
    assert h.recombinant_parents("XB") == ["A.1"]  # `*` stripped, B has no barcode


def test_parent_cycles_terminate():
    entries = {"A": {"parent": "B"}, "B": {"parent": "A"}}
    h = di.Hierarchy({"A", "B"}, entries)
    assert h.parent("A") == "B"  # first reachable barcoded ancestor
    assert h.parent("B") == "A"  # no infinite loop


def test_real_scale_defining_analysis_is_fast(tmp_path):
    rng = random.Random(11)
    pool = [f"A{pos}C" for pos in rng.sample(range(1, 29903), 4000)]
    names = [f"L{i}" for i in range(6000)]
    data = {}
    yml = []
    for i, name in enumerate(names):
        parent = names[(i - 1) // 2] if i else None
        base = data[parent] if parent else []
        data[name] = base + rng.sample(pool, 2) if len(base) < 60 else base
        yml.append(
            f"- name: {name}\n  alias: {name}\n" + (f"  parent: {parent}\n" if parent else "")
        )
    ds = tmp_path / "ds"
    ds.mkdir()
    (ds / "lineage_mutations.json").write_text(json.dumps(data))
    (ds / "lineages.yml").write_text("".join(yml))
    amp = tmp_path / "amp.tsv"
    amp.write_text(
        "sample\tamplicon\tstart\tend\tmean_depth\tmin_depth\tpct_ge_min_depth\tdropped\n"
        + "".join(
            f"S{s}\ta{i}\t{1 + i * 300}\t{400 + i * 300}\t1.0\t0\t0.0\t{'true' if i % 3 else 'false'}\n"
            for s in range(10)
            for i in range(99)
        )
    )
    agg = tmp_path / "agg.tsv"
    lines = ["\tsummarized\tlineages\tabundances"]
    lines += [f"S{s}\t[]\t{' '.join(names[100:125])}\t{' '.join(['0.04'] * 25)}" for s in range(10)]
    agg.write_text("\n".join(lines) + "\n")
    start = time.monotonic()
    code, _impact, summary = run(tmp_path, dataset=ds, amplicons=amp, aggregated=agg)
    assert code == 0
    assert time.monotonic() - start < 10
    assert len(table(summary)) == 10
