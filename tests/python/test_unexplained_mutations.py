import gzip
import json
import math
from fractions import Fraction
from pathlib import Path

import pytest
import unexplained_mutations as um

ROOT = Path(__file__).resolve().parents[2]
FIX = ROOT / "tests" / "fixtures" / "unexplained"
GFF = str(ROOT / "assets" / "references" / "NC_045512.gff3")
REF = str(ROOT / "assets" / "references" / "NC_045512.fasta")
BARCODES_TEST = str(ROOT / "tests" / "data" / "barcodes_test.csv")

HEADER = (
    "REGION\tPOS\tREF\tALT\tREF_DP\tREF_RV\tREF_QUAL\tALT_DP\tALT_RV\tALT_QUAL\tALT_FREQ\tTOTAL_DP\tPVAL\tPASS\t"
    "GFF_FEATURE\tREF_CODON\tREF_AA\tALT_CODON\tALT_AA\tPOS_AA\n"
)
OUT_HEADER = "\t".join(um.COLUMNS)


def vrow(pos, ref, alt, total, alt_dp, freq=None, qual=37):
    """One ivar row; freq defaults to alt_dp / total (written like ivar, 6 significant digits)."""
    freq = alt_dp / total if freq is None else freq
    return (
        f"NC_045512.2\t{pos}\t{ref}\t{alt}\t{total - alt_dp}\t1\t37\t{alt_dp}\t1\t{qual}\t{freq:.6g}\t{total}\t0\tTRUE\t"
        "NA\tNA\tNA\tNA\tNA\tNA\n"
    )


def write_variants(tmp_path, rows, name="v.tsv"):
    path = tmp_path / name
    path.write_text(HEADER + "".join(rows))
    return str(path)


def write_demix(tmp_path, lineages, abundances, name="s.demix"):
    path = tmp_path / name
    path.write_text(
        "\ts.variants.tsv\n"
        "summarized\t[('Omicron', 1.0)]\n"
        f"lineages\t{' '.join(lineages)}\n"
        f"abundances\t{' '.join(abundances)}\n"
        "resid\t1.5\ncoverage\t99.0\npathogen\tSARS-CoV-2\n"
    )
    return str(path)


def run(tmp_path, variants, demix, dataset, *extra, sample="S", out="out.tsv"):
    output = tmp_path / out
    args = ["--variants", variants, "--sample", sample, "-o", str(output)]
    if demix:
        args += ["--demix", demix]
    if dataset:
        args += ["--dataset", dataset]
    rc = um.main([*args, *extra])
    return rc, (output.read_text().splitlines() if output.exists() else [])


def table(lines):
    assert lines[0] == OUT_HEADER
    return [dict(zip(um.COLUMNS, line.split("\t"))) for line in lines[1:]]


def read_gz(path):
    with gzip.open(path, "rt") as handle:
        return handle.read()


# ---------------------------------------------------------------------------------- mini hand-checked fixture


@pytest.fixture
def mini(tmp_path):
    """LIN_A 0.5, LIN_B 0.3, LIN_C 0.15 (minor), LIN_X 0.05 (not in the dataset).

    Expected frequencies: C100T 0.8 (A+B), G200A 0.65 (A+C), T300C 0.5 (A), A400G 0.3 (B), A500C 0.15 (C only).
    """
    dataset = tmp_path / "ds"
    dataset.mkdir()
    (dataset / "lineage_mutations.json").write_text(
        json.dumps(
            {
                "LIN_A": ["C100T", "G200A", "T300C"],
                "LIN_B": ["C100T", "A400G"],
                "LIN_C": ["G200A", "A500C", "A600-"],
                "UNUSED": ["C1T"],
            }
        )
    )
    demix = write_demix(
        tmp_path, ["LIN_A", "LIN_B", "LIN_C", "LIN_X"], ["0.5", "0.3", "0.15", "0.05"]
    )
    return tmp_path, str(dataset), demix


def test_expected_frequencies_are_sums_of_abundances(mini):
    _, dataset, demix = mini
    lineages = um.read_demix(demix)
    assert lineages == [("LIN_A", 0.5), ("LIN_B", 0.3), ("LIN_C", 0.15), ("LIN_X", 0.05)]
    barcodes, n_non_snv = um.load_barcodes(dataset, [n for n, _ in lineages])
    assert sorted(barcodes) == ["LIN_A", "LIN_B", "LIN_C"]  # LIN_X absent, UNUSED not requested
    assert n_non_snv == 1  # A600-
    exp = um.expected_frequencies(lineages, barcodes)
    assert exp[(100, "C", "T")] == pytest.approx(0.8)
    assert exp[(200, "G", "A")] == pytest.approx(0.65)
    assert exp[(300, "T", "C")] == pytest.approx(0.5)
    assert exp[(400, "A", "G")] == pytest.approx(0.3)
    assert exp[(500, "A", "C")] == pytest.approx(0.15)  # minor lineage only
    assert (600, "A", "-") not in exp and len(exp) == 5


def test_mini_exact_calls_and_order(mini, capsys):
    tmp_path, dataset, demix = mini
    variants = write_variants(
        tmp_path,
        [
            vrow(100, "C", "T", 100, 95),  # exp 0.8 -> excess 0.15: called
            vrow(
                200, "G", "A", 100, 70
            ),  # exp 0.65 -> excess exactly 0.05 (borderline min-excess): called
            vrow(300, "T", "C", 100, 54),  # exp 0.5 -> excess 0.04: explained
            vrow(400, "A", "G", 100, 35),  # exp 0.3 -> excess exactly 0.05: called
            vrow(500, "A", "C", 100, 15),  # exactly the minor lineage's share: explained
            vrow(
                510, "A", "C", 100, 16
            ),  # a position no lineage has, just above the minor lineage's share: called
            vrow(520, "A", "C", 100, 40),  # no lineage: called, excess 0.4
            vrow(700, "T", "G", 100, 90),  # '+'/'-' rows below are indels
            "NC_045512.2\t800\tA\t+TT\t10\t1\t37\t90\t1\t37\t0.9\t100\t0\tTRUE\tNA\tNA\tNA\tNA\tNA\tNA\n",
            "NC_045512.2\t900\tA\t-C\t10\t1\t37\t90\t1\t37\t0.9\t100\t0\tTRUE\tNA\tNA\tNA\tNA\tNA\tNA\n",
        ],
    )
    rc, lines = run(tmp_path, variants, demix, dataset, "--max-pvalue", "1")
    assert rc == 0
    rows = table(lines)
    assert [(r["position"], r["alt"]) for r in rows] == [
        ("700", "G"),  # excess 0.9
        ("520", "C"),  # 0.4
        ("510", "C"),  # 0.16
        ("100", "T"),  # 0.15
        ("200", "A"),  # 0.05 (the two borderline rows tie: position order)
        ("400", "G"),  # 0.05
    ]
    by_pos = {r["position"]: r for r in rows}
    assert by_pos["100"]["observed_freq"] == "0.95"
    assert by_pos["100"]["expected_freq"] == "0.8"
    assert by_pos["100"]["excess"] == "0.15"
    assert by_pos["100"]["depth"] == "100" and by_pos["100"]["alt_depth"] == "95"
    # p_value = P(X >= alt_dp), X ~ Binomial(depth, max(expected, 0.01)); independent exact sum with math.comb
    for pos, (n, k, p0) in {"100": (100, 95, 0.8), "520": (100, 40, 0.01)}.items():
        exact = sum(math.comb(n, i) * p0**i * (1 - p0) ** (n - i) for i in range(k, n + 1))
        assert by_pos[pos]["p_value"] == f"{exact:.3g}"
    assert by_pos["200"]["expected_freq"] == "0.65" and by_pos["200"]["excess"] == "0.05"
    assert by_pos["510"]["expected_freq"] == "0"
    assert all(
        r["sample"] == "S" and r["gene"] == "NA" and r["effect"] == "NA" for r in rows
    )  # no GFF given
    err = capsys.readouterr().err
    assert "1 of 4 demixed lineages not in the dataset (LIN_X; abundance 0.05)" in err
    assert "skipped 2 indel row(s) in the variants TSV and 1 non-SNV barcode mutation(s)" in err


@pytest.mark.parametrize(
    ("rows", "called"),
    [
        # min-freq exactly 0.05 with min-excess exactly 0.05, alt depth exactly 10, depth 200
        ([vrow(10, "A", "C", 200, 10)], True),
        ([vrow(10, "A", "C", 200, 10, freq=0.0499)], False),  # just under min-freq
        ([vrow(10, "A", "C", 20, 10)], True),  # depth exactly min-depth
        ([vrow(10, "A", "C", 19, 10)], False),  # one under
        ([vrow(10, "A", "C", 100, 10)], True),  # alt depth exactly min-alt-depth
        ([vrow(10, "A", "C", 100, 9)], False),
        ([vrow(10, "A", "C", 100, 30, qual=20)], True),  # ALT_QUAL exactly 20
        ([vrow(10, "A", "C", 100, 30, qual=19)], False),
        ([vrow(10, "A", "C", 100, 30, qual="NA")], True),  # no quality: not applied
    ],
)
def test_thresholds_are_inclusive(tmp_path, rows, called):
    dataset = tmp_path / "ds"
    dataset.mkdir()
    (dataset / "lineage_mutations.json").write_text(json.dumps({"L": ["C99T"]}))
    demix = write_demix(tmp_path, ["L"], ["1.0"])
    rc, lines = run(
        tmp_path, write_variants(tmp_path, rows), demix, str(dataset), "--max-pvalue", "1"
    )
    assert rc == 0
    assert len(table(lines)) == (1 if called else 0)


def test_expectation_capped_at_one_and_excess_uses_threshold_option(tmp_path):
    dataset = tmp_path / "ds"
    dataset.mkdir()
    (dataset / "lineage_mutations.json").write_text(json.dumps({"A": ["C1T"], "B": ["C1T"]}))
    demix = write_demix(tmp_path, ["A", "B"], ["0.6", "0.6"])  # sums to 1.2: capped to 1.0
    variants = write_variants(tmp_path, [vrow(1, "C", "T", 100, 100, freq=1.0)])
    rc, lines = run(tmp_path, variants, demix, str(dataset))
    assert rc == 0 and table(lines) == []
    # a stricter / looser --min-excess changes the call
    variants = write_variants(tmp_path, [vrow(2, "C", "T", 100, 50)], name="v2.tsv")
    assert (
        table(
            run(
                tmp_path, variants, demix, str(dataset), "--min-excess", "0.6", "--max-pvalue", "1"
            )[1]
        )
        == []
    )
    assert (
        len(
            table(
                run(
                    tmp_path,
                    variants,
                    demix,
                    str(dataset),
                    "--min-excess",
                    "0.5",
                    "--max-pvalue",
                    "1",
                )[1]
            )
        )
        == 1
    )


def test_binomial_guard_suppresses_low_depth_noise(tmp_path):
    dataset = tmp_path / "ds"
    dataset.mkdir()
    (dataset / "lineage_mutations.json").write_text(json.dumps({"L": ["C100T"]}))
    demix = write_demix(tmp_path, ["L"], ["0.7"])
    # 30/35 vs expected 0.7: excess 0.157 passes the plain thresholds but is ordinary sampling noise
    low = write_variants(tmp_path, [vrow(100, "C", "T", 35, 30)])
    assert len(table(run(tmp_path, low, demix, str(dataset), "--max-pvalue", "1")[1])) == 1
    assert table(run(tmp_path, low, demix, str(dataset))[1]) == []
    # the same excess at 2000x is not noise
    deep = write_variants(tmp_path, [vrow(100, "C", "T", 2000, 1700)], name="deep.tsv")
    assert len(table(run(tmp_path, deep, demix, str(dataset))[1])) == 1
    # a novel mutation (expected 0) at 12 % over 100x is unexplained; 3 reads over 100x would also pass 5 %? no: alt-depth
    novel = write_variants(tmp_path, [vrow(7, "A", "C", 100, 12)], name="novel.tsv")
    assert len(table(run(tmp_path, novel, demix, str(dataset))[1])) == 1


def test_binom_sf_exact_values():
    assert um.binom_sf(0, 10, 0.5) == 1.0
    assert um.binom_sf(11, 10, 0.5) == 0.0
    assert um.binom_sf(10, 10, 0.5) == pytest.approx(0.5**10)
    assert um.binom_sf(1, 10, 0.5) == pytest.approx(1 - 0.5**10)
    assert um.binom_sf(5, 10, 0.5) == pytest.approx(0.623046875)
    assert um.binom_sf(3, 5, 0.2) == pytest.approx(0.05792)  # P(X>=3), n=5, p=0.2


def test_max_rows_cap(tmp_path, capsys):
    dataset = tmp_path / "ds"
    dataset.mkdir()
    (dataset / "lineage_mutations.json").write_text(json.dumps({"L": ["C1T"]}))
    demix = write_demix(tmp_path, ["L"], ["1.0"])
    rows = [vrow(1000 + i, "A", "C", 1000, 100 + i) for i in range(30)]
    rc, lines = run(
        tmp_path, write_variants(tmp_path, rows), demix, str(dataset), "--max-rows", "5"
    )
    assert rc == 0
    got = table(lines)
    assert [r["position"] for r in got] == [
        "1029",
        "1028",
        "1027",
        "1026",
        "1025",
    ]  # highest excess first
    assert "30 unexplained mutations; reporting the top 5" in capsys.readouterr().err


def test_max_rows_edge_cases(tmp_path, capsys):
    dataset = tmp_path / "ds"
    dataset.mkdir()
    (dataset / "lineage_mutations.json").write_text(json.dumps({"L": ["C1T"]}))
    demix = write_demix(tmp_path, ["L"], ["1.0"])
    # 6 calls with distinct excess (0.20 .. 0.15), plus two that tie at the cut-off excess 0.15
    rows = [vrow(1000 + i, "A", "C", 1000, 200 - 10 * i) for i in range(6)]
    rows += [vrow(2000, "A", "C", 1000, 150), vrow(2001, "A", "C", 1000, 150)]
    variants = write_variants(tmp_path, rows)
    base = (tmp_path, variants, demix, str(dataset))
    # exactly at the cap: nothing cut, no warning
    rc, lines = run(*base, "--max-rows", "8")
    assert rc == 0 and len(table(lines)) == 8
    assert "reporting the top" not in capsys.readouterr().err
    # one over: the lowest-excess row goes; the tie at the cut is resolved by position
    rc, lines = run(*base, "--max-rows", "7")
    assert "8 unexplained mutations; reporting the top 7 by excess" in capsys.readouterr().err
    got = [r["position"] for r in table(lines)]
    assert got == ["1000", "1001", "1002", "1003", "1004", "1005", "2000"]  # 2001 (tie, later) cut
    # 0 = header only (with the warning), a cap above the row count is a no-op
    rc, lines = run(*base, "--max-rows", "0")
    assert rc == 0 and lines == [OUT_HEADER]
    assert "reporting the top 0 by excess" in capsys.readouterr().err
    rc, lines = run(*base, "--max-rows", "1000")
    assert rc == 0 and len(table(lines)) == 8
    # negative would silently slice off the last rows: refused
    with pytest.raises(SystemExit):
        um.main(["--variants", variants, "--sample", "S", "--max-rows", "-1"])


def test_binom_logsf_matches_exact_sums_and_survives_underflow():
    for n, k, p in [(10, 5, 0.5), (100, 95, 0.8), (100, 40, 0.01), (2000, 1700, 0.7), (50, 3, 0.2)]:
        q = Fraction(str(p))
        exact = sum(math.comb(n, i) * q**i * (1 - q) ** (n - i) for i in range(k, n + 1))
        assert math.exp(um.binom_logsf(k, n, p)) == pytest.approx(float(exact), rel=1e-9)
        assert um.binom_sf(k, n, p) == pytest.approx(float(exact), rel=1e-9)
    assert um.binom_logsf(0, 10, 0.5) == 0.0
    assert um.binom_logsf(11, 10, 0.5) == -math.inf
    # k = n: exactly p**n, i.e. ln = n ln p, far below the float range of p itself
    assert um.binom_logsf(407, 407, 0.01) == pytest.approx(407 * math.log(0.01), rel=1e-12)
    assert um.binom_sf(407, 407, 0.01) == 0.0  # the plain probability underflows ...
    assert um.binom_logsf(407, 407, 0.01) < -1800  # ... the log does not
    # a tail far out in the wings (900 of 1000 at p = 0.01, ~1e-1500) against an exact rational sum
    exact = sum(
        math.comb(1000, i) * Fraction(1, 100) ** i * Fraction(99, 100) ** (1000 - i)
        for i in range(900, 1001)
    )
    exact_ln = math.log10(exact.numerator) * math.log(10) - math.log10(
        exact.denominator
    ) * math.log(10)
    assert um.binom_logsf(900, 1000, 0.01) == pytest.approx(exact_ln, rel=1e-9)


@pytest.mark.parametrize(
    ("log_p", "text"),
    [
        (math.log(1.2345e-7), "1.23e-07"),
        (math.log(1e-299), "1e-299"),
        (math.log(1e-300), "1e-300"),
        (math.log(1e-301), "1e-300"),  # below the floor
        (-1874.3, "1e-300"),
        (-math.inf, "1e-300"),
        (0.0, "1"),
    ],
)
def test_format_pvalue_floor(log_p, text):
    assert um.format_pvalue(log_p) == text
    assert float(um.format_pvalue(log_p)) >= 1e-300  # always a plain, positive float


def test_deep_call_reports_the_floor_not_zero(tmp_path):
    dataset = tmp_path / "ds"
    dataset.mkdir()
    (dataset / "lineage_mutations.json").write_text(json.dumps({"L": ["C100T"]}))
    demix = write_demix(tmp_path, ["L"], ["1.0"])
    # observed 1.0, expected 0 at depth 407 (the real WW-1285E case): 0.01^407 underflows
    variants = write_variants(
        tmp_path,
        [
            vrow(200, "A", "C", 407, 407),  # ln p = -1874: floored
            vrow(300, "A", "C", 120, 120),  # 0.01^120 = 1e-240: reported as is
            vrow(310, "A", "C", 149, 149),  # 0.01^149 = 1e-298: reported as is
            vrow(320, "A", "C", 150, 150),  # 0.01^150 = 1e-300: the floor itself
            vrow(330, "A", "C", 151, 151),  # 1e-302 -> floored
        ],
    )
    rc, lines = run(tmp_path, variants, demix, str(dataset))
    assert rc == 0
    by_pos = {r["position"]: r["p_value"] for r in table(lines)}
    assert by_pos == {
        "200": "1e-300",
        "300": "1e-240",
        "310": "1e-298",
        "320": "1e-300",
        "330": "1e-300",
    }
    assert all(float(v) > 0 for v in by_pos.values())


def test_threshold_uses_the_log_p_value(tmp_path):
    dataset = tmp_path / "ds"
    dataset.mkdir()
    (dataset / "lineage_mutations.json").write_text(json.dumps({"L": ["C100T"]}))
    demix = write_demix(tmp_path, ["L"], ["1.0"])
    variants = write_variants(
        tmp_path, [vrow(200, "A", "C", 407, 407), vrow(300, "A", "C", 100, 12)]
    )
    # 12/100 at a 0.01 background: p ~ 4e-11; ln p of the first is far below any threshold, incl. 1e-320
    assert (
        len(table(run(tmp_path, variants, demix, str(dataset), "--max-pvalue", "1e-320")[1])) == 1
    )
    assert len(table(run(tmp_path, variants, demix, str(dataset), "--max-pvalue", "1e-12")[1])) == 1
    assert len(table(run(tmp_path, variants, demix, str(dataset), "--max-pvalue", "1e-9")[1])) == 2
    assert table(run(tmp_path, variants, demix, str(dataset), "--max-pvalue", "0")[1]) == []
    assert len(table(run(tmp_path, variants, demix, str(dataset), "--max-pvalue", "1")[1])) == 2


def test_concat_passes_the_floored_p_value_through_and_sorts_by_excess(tmp_path):
    def per_sample(name, rows):
        path = tmp_path / f"{name}.tsv"
        path.write_text(OUT_HEADER + "\n" + "".join("\t".join(r) + "\n" for r in rows))
        return str(path)

    row = [
        "S1",
        "22102",
        "A",
        "T",
        "1",
        "0",
        "1",
        "407",
        "407",
        "S",
        "E180D",
        "missense",
        "1e-300",
        "false",
    ]
    other = [
        *row[:1],
        "28802",
        "A",
        "T",
        "0.3",
        "0",
        "0.3",
        "100",
        "30",
        "N",
        "R177*",
        "stop_gained",
        "1.2e-30",
        "true",
    ]
    a = per_sample("S1", [row, other])
    b = per_sample("S2", [["S2", *row[1:]]])
    out = tmp_path / "all.tsv"
    assert um.main(["--concat", a, b, "-o", str(out)]) == 0
    lines = out.read_text().splitlines()
    got = [dict(zip(um.CONCAT_COLUMNS, line.split("\t"))) for line in lines[1:]]
    assert [(r["sample"], r["position"], r["recurrent"], r["p_value"]) for r in got] == [
        ("S1", "22102", "true", "1e-300"),
        ("S2", "22102", "true", "1e-300"),
        ("S1", "28802", "false", "1.2e-30"),
    ]
    assert all(float(r["p_value"]) > 0 for r in got)
    assert [r["in_dataset"] for r in got] == [
        "false",
        "false",
        "true",
    ]  # carried through, last column
    assert lines[0].split("\t")[-2:] == ["p_value", "in_dataset"]


def test_concat_column_lists():
    assert um.COLUMNS[-2:] == ["p_value", "in_dataset"] and len(um.COLUMNS) == 14
    assert um.CONCAT_COLUMNS[12:] == ["n_samples", "recurrent", "p_value", "in_dataset"]
    assert len(um.CONCAT_COLUMNS) == 16


def test_concat_rejects_tables_from_before_in_dataset(tmp_path, capsys):
    old = tmp_path / "old.tsv"
    old.write_text("\t".join(um.COLUMNS[:-1]) + "\nS\t1\tA\tT\t1\t0\t1\t50\t50\tNA\tNA\tNA\t1e-9\n")
    assert um.main(["--concat", str(old), "-o", str(tmp_path / "o.tsv")]) == 2
    err = capsys.readouterr().err
    assert "in_dataset" in err and "regenerated" in err
    # a concat table (16 columns) is not a per-sample table either
    merged = tmp_path / "m.tsv"
    merged.write_text("\t".join(um.CONCAT_COLUMNS) + "\n")
    assert um.main(["--concat", str(merged), "-o", str(tmp_path / "o.tsv")]) == 2


def test_in_dataset_marks_sites_no_lineage_models(tmp_path, mini):
    # mini dataset carries C100T, G200A, T300C, A400G, A500C (+ UNUSED C1T, which is not demixed);
    # an unexplained mutation at C1T (carried only by a lineage that is NOT in the mixture) is in the dataset,
    # one at a site nobody carries is not, and the same position with another alt is not.
    tmp_path, dataset, demix = mini
    variants = write_variants(
        tmp_path,
        [
            vrow(
                1, "C", "T", 100, 50
            ),  # carried by UNUSED only: in_dataset true, still unexplained
            vrow(100, "C", "T", 100, 99),  # explained-ish but excess 0.19: in_dataset true
            vrow(100, "C", "A", 100, 40),  # same position, other alt: false
            vrow(520, "A", "C", 100, 40),  # no lineage anywhere: false
            vrow(
                600, "A", "C", 100, 30
            ),  # the barcode A600- is a deletion token, not an SNV: false
        ],
    )
    rc, lines = run(tmp_path, variants, demix, dataset, "--max-pvalue", "1")
    assert rc == 0
    got = {(r["position"], r["alt"]): r["in_dataset"] for r in table(lines)}
    assert got == {
        ("1", "T"): "true",
        ("100", "T"): "true",
        ("100", "A"): "false",
        ("520", "C"): "false",
        ("600", "C"): "false",
    }


@pytest.mark.parametrize("form", ["json_dir", "json_file", "csv_file", "csv_dir"])
def test_in_dataset_for_every_dataset_form(tmp_path, form):
    muts = {"L1": ["C100T"], "L2": ["G200A", "C100T"], "L3": ["T300C"]}  # L3 is not demixed
    cols = ["C100T", "G200A", "T300C", "A400G"]  # A400G: a column no lineage carries
    csv_text = (
        ","
        + ",".join(cols)
        + "\n"
        + "".join(
            f"{n}," + ",".join("1" if c in m else "0" for c in cols) + "\n" for n, m in muts.items()
        )
    )
    if form == "json_dir":
        d = tmp_path / "d"
        d.mkdir()
        (d / "lineage_mutations.json").write_text(json.dumps(muts))
        dataset = str(d)
    elif form == "json_file":
        (tmp_path / "x.json").write_text(json.dumps(muts))
        dataset = str(tmp_path / "x.json")
    elif form == "csv_file":
        (tmp_path / "x.csv").write_text(csv_text)
        dataset = str(tmp_path / "x.csv")
    else:
        d = tmp_path / "d"
        d.mkdir()
        (d / "usher_barcodes.csv").write_text(csv_text)
        dataset = str(d)
    demix = write_demix(tmp_path, ["L1"], ["1.0"])
    variants = write_variants(
        tmp_path,
        [
            vrow(200, "G", "A", 100, 60),  # carried by L2 (not demixed): in dataset
            vrow(300, "T", "C", 100, 60),  # carried by L3 (not demixed): in dataset
            vrow(400, "A", "G", 100, 60),  # all-zero column: nobody carries it
            vrow(500, "A", "C", 100, 60),
        ],
    )
    rc, lines = run(tmp_path, variants, demix, dataset, "--max-pvalue", "1")
    assert rc == 0
    assert {r["position"]: r["in_dataset"] for r in table(lines)} == {
        "200": "true",
        "300": "true",
        "400": "false",
        "500": "false",
    }


# ---------------------------------------------------------------------------------- demix parsing


def test_read_demix_plain_and_repr_forms(tmp_path):
    plain = tmp_path / "a.demix"
    plain.write_text("\tx\nlineages\tBA.2 AY.4\nabundances\t0.71212100 0.28787900\nresid\t3\n")
    assert um.read_demix(str(plain)) == [("BA.2", 0.712121), ("AY.4", 0.287879)]
    reprs = tmp_path / "b.demix"
    reprs.write_text("lineages\t['BA.2', 'AY.4']\nabundances\t[0.712121 0.287879]\n")
    assert um.read_demix(str(reprs)) == [("BA.2", 0.712121), ("AY.4", 0.287879)]
    numpy_like = tmp_path / "c.demix"
    numpy_like.write_text(
        "lineages\t('JN.1',)\nabundances\tarray([1.])\n".replace("array([1.])", "[1.0]")
    )
    assert um.read_demix(str(numpy_like)) == [("JN.1", 1.0)]


@pytest.mark.parametrize(
    "content",
    [
        "",
        "\tx\nsummarized\t[('Omicron', 1.0)]\n",  # no lineages/abundances rows
        "lineages\tA B\nabundances\t0.5\n",  # length mismatch
        "lineages\tA\nabundances\tabc\n",
    ],
)
def test_read_demix_rejects_unusable_files(tmp_path, content):
    path = tmp_path / "bad.demix"
    path.write_text(content)
    with pytest.raises(ValueError):
        um.read_demix(str(path))


def test_no_usable_demix_gives_header_only(tmp_path, capsys, mini):
    _, dataset, _ = mini
    variants = write_variants(tmp_path, [vrow(520, "A", "C", 100, 40)])
    # 1. no --demix at all (skipped / failed sample)
    rc, lines = run(tmp_path, variants, None, dataset)
    assert rc == 0 and lines == [OUT_HEADER]
    assert "no demix (sample skipped or failed)" in capsys.readouterr().err
    # 2. zero-byte demix
    empty = tmp_path / "empty.demix"
    empty.write_text("")
    rc, lines = run(tmp_path, variants, str(empty), dataset)
    assert rc == 0 and lines == [OUT_HEADER]
    assert "demix unusable" in capsys.readouterr().err
    # 3. missing variants file
    rc, lines = run(tmp_path, str(tmp_path / "nope.tsv"), str(mini[2]), dataset)
    assert rc == 0 and lines == [OUT_HEADER]
    assert "variants unusable" in capsys.readouterr().err
    # 4. a header-only variants file (skipped sample with a demix is not expected, but must work)
    rc, lines = run(tmp_path, write_variants(tmp_path, [], name="h.tsv"), str(mini[2]), dataset)
    assert rc == 0 and lines == [OUT_HEADER]


# ---------------------------------------------------------------------------------- datasets


def test_dataset_resolution(tmp_path):
    demix = write_demix(tmp_path, ["L"], ["1.0"])
    variants = write_variants(
        tmp_path, [vrow(100, "C", "T", 100, 95), vrow(200, "G", "A", 100, 95)]
    )

    csv_path = tmp_path / "custom.csv"
    csv_path.write_text(",C100T,G200A,A5-\nL,1.0,0.0,1.0\nOTHER,0.0,1.0,0.0\n")
    rows = table(run(tmp_path, variants, demix, str(csv_path), "--max-pvalue", "1")[1])
    assert [r["position"] for r in rows] == ["200"]  # C100T is in L's barcode, G200A is not

    # a directory with usher_barcodes.csv and no json
    d_csv = tmp_path / "d_csv"
    d_csv.mkdir()
    (d_csv / "usher_barcodes.csv").write_text(csv_path.read_text())
    rows = table(run(tmp_path, variants, demix, str(d_csv), "--max-pvalue", "1")[1])
    assert [r["position"] for r in rows] == ["200"]

    # feather + json: the json is used
    d_f = tmp_path / "d_f"
    d_f.mkdir()
    (d_f / "usher_barcodes.feather").write_bytes(b"not read")
    (d_f / "lineage_mutations.json").write_text(json.dumps({"L": ["G200A"]}))
    rows = table(run(tmp_path, variants, demix, str(d_f), "--max-pvalue", "1")[1])
    assert [r["position"] for r in rows] == ["100"]

    # a bare .json file
    j = tmp_path / "lm.json"
    j.write_text(json.dumps({"L": ["C100T", "G200A"]}))
    assert table(run(tmp_path, variants, demix, str(j), "--max-pvalue", "1")[1]) == []


@pytest.mark.parametrize(
    "kind", ["feather_only", "feather_file", "empty_dir", "missing", "bad_csv"]
)
def test_dataset_errors_are_exit_2(tmp_path, kind, capsys):
    demix = write_demix(tmp_path, ["L"], ["1.0"])
    variants = write_variants(tmp_path, [vrow(100, "C", "T", 100, 95)])
    if kind == "feather_only":
        d = tmp_path / "d"
        d.mkdir()
        (d / "usher_barcodes.feather").write_bytes(b"x")
        dataset = str(d)
    elif kind == "feather_file":
        f = tmp_path / "b.feather"
        f.write_bytes(b"x")
        dataset = str(f)
    elif kind == "empty_dir":
        d = tmp_path / "d"
        d.mkdir()
        dataset = str(d)
    elif kind == "missing":
        dataset = str(tmp_path / "nope")
    else:
        f = tmp_path / "b.csv"
        f.write_text(",C100T\nL,abc\n")
        dataset = str(f)
    rc, _ = run(tmp_path, variants, demix, dataset)
    assert rc == 2
    assert "error:" in capsys.readouterr().err


def test_all_lineages_missing_from_dataset_does_not_crash(tmp_path, capsys, mini):
    _, dataset, _ = mini
    demix = write_demix(tmp_path, ["NOPE1", "NOPE2"], ["0.6", "0.4"], name="x.demix")
    variants = write_variants(tmp_path, [vrow(100, "C", "T", 100, 95)])
    rc, lines = run(tmp_path, variants, demix, dataset, "--max-pvalue", "1")
    assert rc == 0
    assert len(table(lines)) == 1  # nothing is explained, so the mutation is called
    assert (
        "2 of 2 demixed lineages not in the dataset (NOPE1, NOPE2; abundance 1)"
        in capsys.readouterr().err
    )


# ---------------------------------------------------------------------------------- annotation


@pytest.fixture(scope="module")
def models_and_ref():
    import mutations_of_interest as moi

    return moi.read_cds_models(GFF), moi.read_fasta_sequence(REF)


@pytest.mark.parametrize(
    ("pos", "ref", "alt", "expected"),
    [
        (23403, "A", "G", ("S", "D614G", "missense")),
        (14408, "C", "T", ("ORF1ab", "P4715L", "missense")),  # ORF1a is not listed beside ORF1ab
        (3037, "C", "T", ("ORF1ab", "F924F", "synonymous")),
        (22102, "A", "T", ("S", "E180D", "missense")),
        (28802, "A", "T", ("N", "R177*", "stop_gained")),
        (241, "C", "T", ("NA", "NA", "intergenic")),
        (29800, "G", "A", ("NA", "NA", "intergenic")),
        (21563, "A", "T", ("S", "M1L", "start_lost")),
        (
            27757,
            "T",
            "A",
            ("ORF7a;ORF7b", "*122R;M1K", "start_lost"),
        ),  # overlapping CDS, most severe effect
    ],
)
def test_annotation(models_and_ref, pos, ref, alt, expected):
    models, reference = models_and_ref
    assert reference[pos - 1] == ref
    assert um.annotate(pos, ref, alt, models, reference) == expected


def test_annotation_without_gff_is_na():
    assert um.annotate(23403, "A", "G", None, None) == ("NA", "NA", "NA")


def test_gff_and_reference_must_come_together(tmp_path, mini):
    _, dataset, demix = mini
    variants = write_variants(tmp_path, [vrow(100, "C", "T", 100, 95)])
    assert run(tmp_path, variants, demix, dataset, "--gff", GFF)[0] == 2


# ---------------------------------------------------------------------------------- real data


def stage(tmp_path, gz_name, out_name):
    path = tmp_path / out_name
    path.write_text(read_gz(FIX / gz_name))
    return str(path)


ANNOT = ["--gff", GFF, "--reference", REF]
SYNTHETIC = ["PURE_C", "MIX_A", "MIX_B", "Site_1_week_2", "NEG_CONTAM"]


@pytest.mark.parametrize("sample", SYNTHETIC)
def test_real_synthetic_samples_are_fully_explained_with_the_demixing_dataset(tmp_path, sample):
    """Default (Bonferroni) guard, the custom barcodes CSV the samples were demixed with: 0 rows."""
    variants = stage(tmp_path, f"{sample}.variants.tsv.gz", f"{sample}.variants.tsv")
    rc, lines = run(
        tmp_path, variants, str(FIX / f"{sample}.demix"), BARCODES_TEST, *ANNOT, sample=sample
    )
    assert rc == 0 and lines == [OUT_HEADER]


@pytest.mark.parametrize("sample", ["PURE_C", "MIX_A", "MIX_B", "Site_1_week_2"])
def test_real_samples_with_bundled_style_messier_demix(tmp_path, sample):
    """Freyja's container-bundled dataset assigns extra minor lineages; their mutations are in the subset JSON."""
    variants = stage(tmp_path, f"{sample}.variants.tsv.gz", f"{sample}.variants.tsv")
    demix = str(FIX / "bundled_subset" / f"{sample}.bundled.demix")
    rc, lines = run(tmp_path, variants, demix, str(FIX / "bundled_subset"), *ANNOT, sample=sample)
    assert rc == 0 and lines == [OUT_HEADER]


@pytest.mark.parametrize("sample", ["MIX_A", "MIX_B", "Site_1_week_2"])
def test_real_samples_without_the_guard_show_the_noise(tmp_path, sample):
    """--max-pvalue 1 = the plain requested rules: chance calls at 20-130x depth appear, with their p-values."""
    variants = stage(tmp_path, f"{sample}.variants.tsv.gz", f"{sample}.variants.tsv")
    rc, lines = run(
        tmp_path,
        variants,
        str(FIX / f"{sample}.demix"),
        BARCODES_TEST,
        "--max-pvalue",
        "1",
        sample=sample,
    )
    rows = table(lines)
    assert rc == 0 and 5 <= len(rows) <= 30
    assert all(int(r["depth"]) < 140 and float(r["excess"]) < 0.4 for r in rows)
    assert all(
        float(r["p_value"]) > um.DEFAULT_MAX_PVALUE for r in rows
    )  # none of them would pass the default


def test_site_1_week_2_chance_call_needs_a_lax_per_site_cutoff(tmp_path):
    """9866 C>T, 17/22 reads vs expected 0.39: p = 3.3e-4. The old 0.001 let it through, 0.05/29903 does not."""
    variants = stage(tmp_path, "Site_1_week_2.variants.tsv.gz", "Site_1_week_2.variants.tsv")
    demix = str(FIX / "Site_1_week_2.demix")
    rc, lines = run(tmp_path, variants, demix, BARCODES_TEST, "--max-pvalue", "0.001", sample="S")
    rows = table(lines)
    assert rc == 0 and [(r["position"], r["alt"]) for r in rows] == [("9866", "T")]
    assert rows[0]["p_value"] == "0.000332"
    assert table(run(tmp_path, variants, demix, BARCODES_TEST, sample="S")[1]) == []


def test_default_max_pvalue_is_bonferroni_over_the_genome():
    assert pytest.approx(0.05 / 29903) == um.DEFAULT_MAX_PVALUE
    assert um.GENOME_LENGTH == 29903


def test_help_shows_the_default_pvalue(capsys):
    with pytest.raises(SystemExit):
        um.main(["--help"])
    assert "0.05/29903" in capsys.readouterr().out


def test_pure_c_with_spiked_snvs_reports_exactly_the_spikes(tmp_path):
    variants = stage(tmp_path, "PURE_C_spiked.variants.tsv.gz", "PURE_C.variants.tsv")
    rc, lines = run(
        tmp_path,
        variants,
        str(FIX / "PURE_C.demix"),
        BARCODES_TEST,
        "--gff",
        GFF,
        "--reference",
        REF,
        sample="PURE_C",
    )
    assert rc == 0
    assert lines[1:] == [
        "PURE_C\t22102\tA\tT\t0.3\t0\t0.3\t120\t36\tS\tE180D\tmissense\t2.39e-42\tfalse",
        "PURE_C\t28802\tA\tT\t0.18\t0\t0.18\t100\t18\tN\tR177*\tstop_gained\t1.41e-17\tfalse",
        "PURE_C\t29800\tG\tA\t0.12\t0\t0.12\t150\t18\tNA\tNA\tintergenic\t2.28e-14\tfalse",
    ]
    # JN.1's own mutations (C241T, S:N501Y = A23063T) are in the input at ~100 % and are not reported
    text = read_gz(FIX / "PURE_C_spiked.variants.tsv.gz")
    assert "\t241\tC\tT\t" in text and "\t23063\tA\tT\t" in text


def test_spike_explained_by_a_lineage_in_the_dataset_is_not_reported(tmp_path):
    """The same spiked mutation is quiet once the (hypothetical) lineage carrying it is in the mixture."""
    variants = stage(tmp_path, "PURE_C_spiked.variants.tsv.gz", "PURE_C.variants.tsv")
    barcodes = tmp_path / "bc.csv"
    # JN.1 + a made-up lineage NOVEL.1 carrying A22102T; reuse only the mutation columns present in PURE_C's barcode
    header, jn1 = None, None
    for line in Path(BARCODES_TEST).read_text().splitlines():
        if line.startswith(","):
            header = line
        elif line.startswith("JN.1,"):
            jn1 = line
    barcodes.write_text(f"{header},A22102T\n{jn1},0.0\n{jn1},1.0\n")
    # first data row stays JN.1; relabel the second as NOVEL.1
    lines = barcodes.read_text().splitlines()
    lines[2] = lines[2].replace("JN.1,", "NOVEL.1,", 1)
    barcodes.write_text("\n".join(lines) + "\n")
    demix = write_demix(tmp_path, ["JN.1", "NOVEL.1"], ["0.7", "0.3"], name="two.demix")
    rc, out = run(
        tmp_path, variants, demix, str(barcodes), "--gff", GFF, "--reference", REF, sample="PURE_C"
    )
    assert rc == 0
    positions = [r["position"] for r in table(out)]
    assert "22102" not in positions  # now expected at 0.3 (observed 0.3)
    assert "28802" in positions and "29800" in positions


# ---------------------------------------------------------------------------------- concat


def concat(tmp_path, names, *extra):
    out = tmp_path / "all.tsv"
    paths = [str(FIX / "concat" / f"{n}.unexplained_mutations.tsv") for n in names]
    rc = um.main(["--concat", *paths, "-o", str(out), *extra])
    return rc, out.read_text().splitlines()


def test_concat_recurrence_and_order(tmp_path):
    rc, lines = concat(tmp_path, ["S_THREE", "S_ONE", "S_SKIPPED", "S_TWO"])
    assert rc == 0
    assert lines[0] == "\t".join(um.CONCAT_COLUMNS)
    rows = [dict(zip(um.CONCAT_COLUMNS, line.split("\t"))) for line in lines[1:]]
    assert [(r["sample"], r["position"], r["n_samples"], r["recurrent"]) for r in rows] == [
        ("S_ONE", "22102", "3", "true"),  # recurrent first, then excess descending
        ("S_TWO", "22102", "3", "true"),
        ("S_THREE", "22102", "3", "true"),
        ("S_TWO", "23063", "1", "false"),  # excess 0.2
        ("S_ONE", "28802", "1", "false"),  # excess 0.18
    ]
    assert rows[0]["aa_change"] == "E180D" and rows[0]["excess"] == "0.3"


def test_concat_min_recurrence_option(tmp_path):
    rc, lines = concat(tmp_path, ["S_ONE", "S_TWO"], "--min-recurrence", "3")
    assert rc == 0
    assert [line.split("\t")[-4:-2] for line in lines[1:]] == [["2", "false"]] * 2 + [
        ["1", "false"]
    ] * 2
    rc, lines = concat(tmp_path, ["S_ONE", "S_TWO"], "--min-recurrence", "1")
    assert {line.split("\t")[-3] for line in lines[1:]} == {"true"}


def test_concat_same_sample_twice_does_not_count_as_recurrence(tmp_path):
    rc, lines = concat(tmp_path, ["S_ONE", "S_ONE"])
    assert rc == 0
    assert {line.split("\t")[-4] for line in lines[1:]} == {"1"}


def test_concat_zero_rows(tmp_path):
    rc, lines = concat(tmp_path, ["S_SKIPPED"])
    assert rc == 0
    assert lines == ["\t".join(um.CONCAT_COLUMNS)]


def test_concat_rejects_foreign_header_and_ragged_rows(tmp_path):
    bad = tmp_path / "bad.tsv"
    bad.write_text("sample\tother\nS\t1\n")
    assert um.main(["--concat", str(bad), "-o", str(tmp_path / "o.tsv")]) == 2
    ragged = tmp_path / "ragged.tsv"
    ragged.write_text(OUT_HEADER + "\nS\t1\t2\n")
    assert um.main(["--concat", str(ragged), "-o", str(tmp_path / "o.tsv")]) == 2


def test_end_to_end_three_spiked_samples_are_recurrent(tmp_path):
    outs = []
    for sample in ("A", "B", "C"):
        variants = stage(tmp_path, "PURE_C_spiked.variants.tsv.gz", f"{sample}.variants.tsv")
        output = tmp_path / f"{sample}.unexplained_mutations.tsv"
        rc = um.main(
            [
                "--variants", variants, "--demix", str(FIX / "PURE_C.demix"), "--dataset", BARCODES_TEST,
                "--sample", sample, "--gff", GFF, "--reference", REF, "-o", str(output),
            ]
        )  # fmt: skip
        assert rc == 0
        outs.append(str(output))
    merged = tmp_path / "unexplained_mutations.tsv"
    assert um.main(["--concat", *outs, "-o", str(merged)]) == 0
    rows = [line.split("\t") for line in merged.read_text().splitlines()[1:]]
    assert len(rows) == 9
    assert {r[-4] for r in rows} == {"3"} and {r[-3] for r in rows} == {"true"}
    assert [r[1] for r in rows[:3]] == ["22102"] * 3  # highest excess first
