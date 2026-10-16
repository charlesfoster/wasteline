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
    assert summary.read_text() == f"{SUMMARY_HEADER}\nS1\t0\t4\t0\t0.00\t0\n"


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
