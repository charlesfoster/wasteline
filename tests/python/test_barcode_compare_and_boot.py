import barcode_compare as bc
import freyja_boot_summary as fb

DEMIX_A = "\tS1.variants.tsv\nsummarized\t[('Omicron', 1.0)]\nlineages\tBA.2 AY.4\nabundances\t0.7 0.3\nresid\t1\n"
DEMIX_B = "\tS1.variants.tsv\nlineages\tBA.2 BA.5\nabundances\t0.6 0.4\n"
BOOT = ",BA.2,AY.4,BA.5\n0.025,0.6,0.25,0.0\n0.5,0.7,0.3,0.0\n0.975,0.8,0.35,0.0"


def test_compare_zero_fills_union_and_boot_summary_uses_2_5_97_5_percentiles(tmp_path):
    a, b = tmp_path / "a.demix", tmp_path / "b.demix"
    a.write_text(DEMIX_A)
    b.write_text(DEMIX_B)
    rows = bc.compare_rows("S1", bc.read_demix(a), bc.read_demix(b), "v1", "v2")
    assert [r[1] for r in rows] == ["BA.2", "BA.5", "AY.4"]
    assert rows[0] == ["S1", "BA.2", "0.700000", "0.600000", "v1", "v2"]
    assert rows[1][2:4] == ["0.000000", "0.400000"] and rows[2][2:4] == ["0.300000", "0.000000"]

    (tmp_path / "boot.csv").write_text(BOOT)
    table = fb.read_quantiles(tmp_path / "boot.csv")
    summary = fb.summary_rows("S1", table, fb.read_demix(a))
    # BA.5 has no point estimate and a zero upper bound: not listed
    assert summary == [
        ["S1", "BA.2", "0.700000", "0.600000", "0.800000"],
        ["S1", "AY.4", "0.300000", "0.250000", "0.350000"],
    ]
