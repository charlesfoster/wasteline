from pathlib import Path

import amplicon_coverage as ac

SCHEMES = Path(__file__).resolve().parents[2] / "assets" / "primer_schemes"
ALL_SCHEMES = sorted(p.stem for p in SCHEMES.glob("*.bed"))


def write_depths(path, depths, contig="ref"):
    path.write_text("".join(f"{contig}\t{i}\tA\t{d}\n" for i, d in enumerate(depths, start=1)))
    return str(path)


def write_bed(path, lines):
    path.write_text("".join(line + "\n" for line in lines))
    return str(path)


def test_rows_values_and_dropout(tmp_path):
    # amplicon a_1 = 1-based 11..20 (BED 10-20), amplicon a_2 = 21..30
    bed = write_bed(
        tmp_path / "p.bed",
        ["ref\t10\t12\ta_1_F", "ref\t18\t20\ta_1_R", "ref\t20\t22\ta_2_F", "ref\t28\t30\ta_2_R"],
    )
    depths = [0] * 10 + [20] * 10 + [4, 4, 4, 4, 4, 4, 4, 4, 4, 14]
    dep = write_depths(tmp_path / "d.depths", depths)
    out = tmp_path / "out.tsv"
    assert (
        ac.main(
            ["--depths", dep, "--bed", bed, "--sample", "S_1", "--min-depth", "10", "-o", str(out)]
        )
        == 0
    )
    lines = out.read_text().splitlines()
    assert lines[0].split("\t") == ac.COLUMNS
    r1, r2 = (line.split("\t") for line in lines[1:])
    assert r1 == ["S_1", "a_1", "11", "20", "20.00", "20", "100.00", "false"]
    assert r2 == ["S_1", "a_2", "21", "30", "5.00", "4", "10.00", "true"]


def test_bundled_scheme_amplicon_counts():
    """Amplicon derivation copes with every bundled scheme's primer naming (alts, orphans, swift starts)."""
    expected = {
        "articv3": 98,
        "articv4.1": 99,
        "articv5.3.2": 96,
        "midnight": 29,
        "midnightv2": 29,
        "eden": 14,
        "qiagenv1": 223,
        "swift": 343,
        "xgen": 343,
    }
    assert ALL_SCHEMES == sorted(expected)
    for scheme, n in expected.items():
        amplicons = ac.derive_amplicons(ac.read_primers(str(SCHEMES / f"{scheme}.bed")))
        assert len(amplicons) == n, scheme
        assert [a[2] for a in amplicons] == sorted(a[2] for a in amplicons)
        assert all(0 <= a[2] < a[3] <= 29903 for a in amplicons)
    # primalscheme v3 names (<scheme>_<n>_LEFT_<alt>, e.g. the ARTIC measles panel): the alts share one amplicon
    assert ac.split_primer_name("177e6ebb_3_LEFT_2") == ("177e6ebb_3", "L")
    assert ac.split_primer_name("177e6ebb_3_RIGHT_0") == ("177e6ebb_3", "R")


def test_single_contig_alias_is_accepted(tmp_path):
    bed = write_bed(tmp_path / "p.bed", ["MN908947.3\t0\t5\ta_1_F", "MN908947.3\t5\t10\ta_1_R"])
    dep = write_depths(tmp_path / "d.depths", [30] * 10, contig="NC_045512.2")
    out = tmp_path / "out.tsv"
    assert ac.main(["--depths", dep, "--bed", bed, "--sample", "S", "-o", str(out)]) == 0
    assert out.read_text().splitlines()[1].split("\t")[4] == "30.00"
