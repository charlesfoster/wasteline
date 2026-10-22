import gzip
import json
import math
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


# ---------------------------------------------------------------------------------- real data


def stage(tmp_path, gz_name, out_name):
    path = tmp_path / out_name
    path.write_text(read_gz(FIX / gz_name))
    return str(path)


def test_real_data_spikes_reported_and_chance_call_suppressed(tmp_path):
    # Site_1_week_2 9866 C>T, 17/22 reads vs expected 0.39: p = 3.3e-4. A 0.001 cutoff let it through;
    # the default 0.05/29903 does not.
    variants = stage(tmp_path, "Site_1_week_2.variants.tsv.gz", "Site_1_week_2.variants.tsv")
    demix = str(FIX / "Site_1_week_2.demix")
    rc, lines = run(tmp_path, variants, demix, BARCODES_TEST, "--max-pvalue", "0.001", sample="S")
    rows = table(lines)
    assert rc == 0 and [(r["position"], r["alt"]) for r in rows] == [("9866", "T")]
    assert rows[0]["p_value"] == "0.000332"
    assert table(run(tmp_path, variants, demix, BARCODES_TEST, sample="S")[1]) == []

    # PURE_C with three spiked SNVs: exactly the spikes are reported, annotated
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


def test_other_pathogen_dataset_and_significance_follows_reference_length(tmp_path):
    """A dataset directory of another pathogen has X_barcodes.csv; 0.05 / genome length follows the reference (PLAN D22)."""
    ds = tmp_path / "ds"
    ds.mkdir()
    (ds / "MEASLESgenome_barcodes.csv").write_text(",C100T\nLIN_A,1\n")
    barcodes, _ = um.load_barcodes(str(ds), {"LIN_A"})
    assert "LIN_A" in barcodes
    # the pipeline's default level, for the 15.9 kb measles genome
    assert um.effective_max_pvalue(1.67e-6, 15894) == pytest.approx(0.05 / 15894)
    assert um.effective_max_pvalue(1.67e-6, 29903) == 1.67e-6  # SARS-CoV-2: untouched
    assert um.effective_max_pvalue(1e-3, 15894) == 1e-3  # an explicit value is never rescaled


def test_contig_selects_fasta_record_and_gff_seqid(tmp_path):
    """A two-contig reference + `--contig` gives the single-contig result; without it the run is an error."""
    seq = "".join(Path(REF).read_text().splitlines()[1:])
    fa = tmp_path / "multi.fa"
    fa.write_text(f">segA first\n{seq}\n>segB\n{'A' * len(seq)}\n")
    gff_lines = [
        line.replace("NC_045512.2", "segA", 1) for line in Path(GFF).read_text().splitlines()
    ]
    cds = next(line for line in gff_lines if "\tCDS\t" in line)
    gff = tmp_path / "multi.gff3"
    gff.write_text("\n".join([*gff_lines, cds.replace("segA", "segB", 1)]) + "\n")
    variants = stage(tmp_path, "PURE_C_spiked.variants.tsv.gz", "PURE_C.variants.tsv")
    args = [variants, str(FIX / "PURE_C.demix"), BARCODES_TEST]
    rc, single = run(tmp_path, *args, "--gff", GFF, "--reference", REF, out="single.tsv")
    assert rc == 0 and len(single) == 4
    rc, lines = run(
        tmp_path, *args, "--gff", str(gff), "--reference", str(fa), "--contig", "segA", out="m.tsv"
    )
    assert rc == 0 and lines == single
    assert run(tmp_path, *args, "--gff", str(gff), "--reference", str(fa), out="n.tsv")[0] == 2
