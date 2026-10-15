import json

import barcode_compare as bc
import freyja_boot_summary as fb

DEMIX_A = "\tS1.variants.tsv\nsummarized\t[('Omicron', 1.0)]\nlineages\tBA.2 AY.4\nabundances\t0.7 0.3\nresid\t1\n"
DEMIX_B = "\tS1.variants.tsv\nlineages\tBA.2 BA.5\nabundances\t0.6 0.4\n"
BOOT = ",BA.2,AY.4,BA.5\n0.025,0.6,0.25,0.0\n0.5,0.7,0.3,0.0\n0.975,0.8,0.35,0.0"


def test_compare_rows_union_of_lineages_with_zero_fill(tmp_path):
    a, b = tmp_path / "a.demix", tmp_path / "b.demix"
    a.write_text(DEMIX_A)
    b.write_text(DEMIX_B)
    rows = bc.compare_rows("S1", bc.read_demix(a), bc.read_demix(b), "v1", "v2")
    assert [r[1] for r in rows] == ["BA.2", "BA.5", "AY.4"]
    assert rows[0] == ["S1", "BA.2", "0.700000", "0.600000", "v1", "v2"]
    assert rows[1][2:4] == ["0.000000", "0.400000"] and rows[2][2:4] == ["0.300000", "0.000000"]


def test_compare_cli(tmp_path):
    for name, text in (("a.demix", DEMIX_A), ("b.demix", DEMIX_B)):
        (tmp_path / name).write_text(text)
    for name, version in (("va.json", "v1"), ("vb.json", "v2")):
        (tmp_path / name).write_text(json.dumps({"barcode_version": version}))
    out = tmp_path / "rows.tsv"
    args = ["--sample", "S1", "--a", tmp_path / "a.demix", "--b", tmp_path / "b.demix"]
    args += ["--version-a", tmp_path / "va.json", "--version-b", tmp_path / "vb.json", "-o", out]
    assert bc.main([str(x) for x in args]) == 0
    assert out.read_text().splitlines()[0] == "S1\tBA.2\t0.700000\t0.600000\tv1\tv2"


def test_boot_summary_uses_point_estimate_and_2_5_97_5_percentiles(tmp_path):
    (tmp_path / "boot.csv").write_text(BOOT)
    (tmp_path / "s.demix").write_text(DEMIX_A)
    table = fb.read_quantiles(tmp_path / "boot.csv")
    rows = fb.summary_rows("S1", table, fb.read_demix(tmp_path / "s.demix"))
    # BA.5 has no point estimate and a zero upper bound: not listed
    assert rows == [
        ["S1", "BA.2", "0.700000", "0.600000", "0.800000"],
        ["S1", "AY.4", "0.300000", "0.250000", "0.350000"],
    ]


def test_boot_summary_rejects_a_table_without_the_percentile_rows(tmp_path):
    (tmp_path / "boot.csv").write_text(",BA.2\n0.5,0.7\n")
    (tmp_path / "s.demix").write_text(DEMIX_A)
    out = tmp_path / "o.tsv"
    args = [
        "--sample",
        "S1",
        "--lineages",
        tmp_path / "boot.csv",
        "--demix",
        tmp_path / "s.demix",
        "-o",
        out,
    ]
    assert fb.main([str(x) for x in args]) == 1
