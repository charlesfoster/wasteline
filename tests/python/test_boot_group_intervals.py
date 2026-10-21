from pathlib import Path

import boot_group_intervals as bgi
import pytest

FIX = Path(__file__).resolve().parents[2] / "tests" / "fixtures" / "lineage_groups"
GROUPS = """group,include,exclude,in_total
JN.1.16,JN.1.16,,yes
KP.3*,KP.3*,,yes
JN.1* (all),JN.1*,,no
BA.5,BA.5,,yes
"""
AGGREGATED = (
    "\tsummarized\tlineages\tabundances\tresid\tcoverage\tpathogen\n"
    "S_1.variants.tsv\t[('Omicron', 1.0)]\tJN.1.16 KP.3.1 BA.5 A.1 ZZ.9\t"
    "0.45 0.35 0.1 0.08 0.02\t1.5\t90.0\tSARS-CoV-2\n"
    "NOBOOT.variants.tsv\t[('Omicron', 1.0)]\tBA.5\t1.0\t1.5\t90.0\tSARS-CoV-2\n"
    "NEG_1.variants.tsv\t[('Omicron', 1.0)]\tBA.5\t1.0\t1.5\t90.0\tSARS-CoV-2\n"
)
SAMPLESHEET = "sample,sample_type\nS_1,sample\nNEG_1,negative\n"
# 5 draws; columns JN.1.16, KP.3.1, BA.5, A.1 and ZZ.9 (unknown to the dataset, so `other`)
BOOTS = """,JN.1.16,KP.3.1,BA.5,A.1,ZZ.9
0,0.5,0.3,0.1,0.1,0.0
1,0.4,0.4,0.0,0.2,0.0
2,0.3,0.5,0.1,0.1,0.0
3,0.2,0.6,0.2,0.0,0.0
4,0.1,0.7,0.0,0.1,0.1
"""


def run(tmp_path, dataset, boots=("S_1", "NOAGG", "NEG_1")):
    (tmp_path / "groups.csv").write_text(GROUPS)
    (tmp_path / "sheet.csv").write_text(SAMPLESHEET)
    (tmp_path / "agg.tsv").write_text(AGGREGATED)
    files = []
    for name in boots:
        path = tmp_path / f"{name}.boot_lineages_boot.csv"
        path.write_text(BOOTS)
        files.append(str(path))
    out = tmp_path / "out.tsv"
    code = bgi.main(
        ["--aggregated", str(tmp_path / "agg.tsv"), "--groups", str(tmp_path / "groups.csv")]
        + [
            "--dataset",
            str(dataset),
            "--samplesheet",
            str(tmp_path / "sheet.csv"),
            "--boots",
            *files,
            "--output",
            str(out),
        ]
    )
    return code, out


def test_hand_computed_group_intervals(tmp_path):
    code, out = run(tmp_path, FIX)
    assert code == 0
    lines = out.read_text().splitlines()
    assert lines[0] == "sample\tgroup\tabundance\tlo\thi\tn_boot"
    rows = [line.split("\t") for line in lines[1:]]
    # NOAGG has a boot file but no aggregate row, NOBOOT the reverse, NEG_1 is a negative control in the samplesheet
    # (left out like in lineage_summary.csv): all three are absent
    assert [r[0] for r in rows] == ["S_1"] * 5
    assert [r[1] for r in rows] == ["JN.1.16", "KP.3*", "JN.1* (all)", "BA.5", "other"]
    # point estimate = the demix group sum (KP.3.1 counts in KP.3* and JN.1* (all); A.1 + ZZ.9 = other);
    # percentiles over the 5 draws with linear interpolation (position 0.1 and 3.9 of the sorted draws)
    expected = {
        "JN.1.16": (0.45, 0.11, 0.49),  # draws .5 .4 .3 .2 .1
        "KP.3*": (0.35, 0.31, 0.69),  # draws .3 .4 .5 .6 .7
        "JN.1* (all)": (0.8, 0.8, 0.8),  # JN.1.16 + KP.3.1 = .8 in every draw
        "BA.5": (0.1, 0.0, 0.19),  # sorted draws 0 0 .1 .1 .2
        "other": (0.1, 0.01, 0.2),  # A.1 + ZZ.9, sorted 0 .1 .1 .2 .2
    }
    for _sample, group, abundance, lo, hi, n_boot in rows:
        assert (float(abundance), float(lo), float(hi)) == pytest.approx(expected[group], abs=1e-9)
        assert n_boot == "5"
