from pathlib import Path

import amplicon_coverage as ac
import pytest

SCHEMES = Path(__file__).resolve().parents[2] / "assets" / "primer_schemes"
ALL_SCHEMES = sorted(p.stem for p in SCHEMES.glob("*.bed"))


def write_depths(path, depths, contig="ref"):
    path.write_text("".join(f"{contig}\t{i}\tA\t{d}\n" for i, d in enumerate(depths, start=1)))
    return str(path)


def write_bed(path, lines):
    path.write_text("".join(line + "\n" for line in lines))
    return str(path)


def test_nine_schemes_are_bundled():
    assert len(ALL_SCHEMES) == 9


@pytest.mark.parametrize("scheme", ALL_SCHEMES)
def test_bundled_scheme_yields_sane_amplicons(scheme):
    amplicons = ac.derive_amplicons(ac.read_primers(str(SCHEMES / f"{scheme}.bed")))
    assert len(amplicons) >= 14
    names = [a[1] for a in amplicons]
    assert len(set(names)) == len(names)
    starts = [a[2] for a in amplicons]
    assert starts == sorted(starts)
    for _contig, _name, start, end in amplicons:
        assert 100 <= end - start <= 3000, (scheme, _name, start, end)
        assert 0 <= start < end <= 29903
    assert amplicons[0][2] < 100  # tiling starts at the 5' end
    assert amplicons[-1][3] > 29700  # and ends at the 3' end


@pytest.mark.parametrize(
    "scheme,expected",
    [
        ("articv3", 98),
        ("articv4.1", 99),
        ("articv5.3.2", 96),
        ("midnight", 29),
        ("midnightv2", 29),
        ("eden", 14),
        ("qiagenv1", 223),
        ("swift", 343),
        ("xgen", 343),
    ],
)
def test_bundled_scheme_amplicon_counts(scheme, expected):
    assert len(ac.derive_amplicons(ac.read_primers(str(SCHEMES / f"{scheme}.bed")))) == expected


@pytest.mark.parametrize(
    "name,key,direction",
    [
        ("nCoV-2019_1_F", "nCoV-2019_1", "L"),
        ("nCoV-2019_1_R", "nCoV-2019_1", "R"),
        ("nCoV-2019_7_alt3_F", "nCoV-2019_7", "L"),
        ("nCoV-2019_7_alt0_R", "nCoV-2019_7", "R"),
        ("SARS-CoV-2_400_12_F", "SARS-CoV-2_400_12", "L"),
        ("QIAseq_221-2_LEFT", "QIAseq_221", "L"),
        ("QIAseq_5_RIGHT", "QIAseq_5", "R"),
        ("covid19genome_200-29703_s7490_U_88_F", "covid19genome_200-29703_s7490_U_88", "L"),
        ("covid19genome_200-29703_s7490_U_88R", "covid19genome_200-29703_s7490_U_88", "R"),
        ("SARS2_A1_F", "SARS2_A1", "L"),
    ],
)
def test_split_primer_name(name, key, direction):
    assert ac.split_primer_name(name) == (key, direction)


def test_split_primer_name_falls_back_to_strand():
    assert ac.split_primer_name("weird", "+") == ("weird", "L")
    assert ac.split_primer_name("weird", "-") == ("weird", "R")
    with pytest.raises(ValueError):
        ac.split_primer_name("weird")


def test_alternates_widen_the_amplicon(tmp_path):
    bed = write_bed(
        tmp_path / "a.bed",
        [
            "c\t100\t120\tx_1_F\tp\t+",
            "c\t90\t110\tx_1_alt1_F\tp\t+",
            "c\t400\t420\tx_1_R\tp\t-",
            "c\t430\t450\tx_1_alt2_R\tp\t-",
        ],
    )
    assert ac.derive_amplicons(ac.read_primers(bed)) == [("c", "x_1", 90, 450)]


def test_swift_style_start_primers_join_nearest_amplicon(tmp_path):
    bed = write_bed(
        tmp_path / "s.bed",
        [
            "c\t3\t25\tcovid19genome_0-200_s0_M1_F\t1\t+",
            "c\t33\t56\tcovid19genome_0-200_s0_M2_F\t60\t+",
            "c\t182\t200\tcovid19genome_200-29703_s7490_U_88_F\t54\t+",
            "c\t301\t322\tcovid19genome_200-29703_s7490_U_88_R\t54\t-",
            "c\t234\t253\tcovid19genome_200-29703_s7490_U_87_F\t54\t+",
            "c\t357\t375\tcovid19genome_200-29703_s7490_U_87_R\t54\t-",
        ],
    )
    amplicons = ac.derive_amplicons(ac.read_primers(bed))
    assert [(a[2], a[3]) for a in amplicons] == [(3, 322), (234, 375)]


def test_orphans_paired_by_trailing_number(tmp_path):
    bed = write_bed(
        tmp_path / "m.bed",
        [
            "c\t10\t30\tSARSCoV_1200_1_F",
            "c\t500\t520\tSARSCoV_1200_1_R",
            "c\t400\t420\tSARSCoV_1200_2_F",
            "c\t900\t920\tnCoV_1200_2_R",
        ],
    )
    amplicons = ac.derive_amplicons(ac.read_primers(bed))
    assert [(a[1], a[2], a[3]) for a in amplicons] == [
        ("SARSCoV_1200_1", 10, 520),
        ("SARSCoV_1200_2", 400, 920),
    ]


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


def test_missing_positions_count_as_zero(tmp_path):
    bed = write_bed(tmp_path / "p.bed", ["ref\t0\t10\ta_1_F", "ref\t30\t40\ta_1_R"])
    dep = write_depths(tmp_path / "d.depths", [50] * 20)  # file ends before the amplicon does
    rows = ac.amplicon_rows("S", ac.read_depths(dep), ac.derive_amplicons(ac.read_primers(bed)), 10)
    assert rows[0][4] == "25.00"  # 20 positions at 50, 20 at 0
    assert rows[0][5] == "0"


def test_header_written_with_zero_amplicons(tmp_path):
    bed = write_bed(tmp_path / "p.bed", [])
    dep = write_depths(tmp_path / "d.depths", [5] * 10)
    out = tmp_path / "out.tsv"
    assert ac.main(["--depths", dep, "--bed", bed, "--sample", "S", "-o", str(out)]) == 0
    assert out.read_text() == "\t".join(ac.COLUMNS) + "\n"


def test_single_contig_alias_is_accepted(tmp_path):
    bed = write_bed(tmp_path / "p.bed", ["MN908947.3\t0\t5\ta_1_F", "MN908947.3\t5\t10\ta_1_R"])
    dep = write_depths(tmp_path / "d.depths", [30] * 10, contig="NC_045512.2")
    out = tmp_path / "out.tsv"
    assert ac.main(["--depths", dep, "--bed", bed, "--sample", "S", "-o", str(out)]) == 0
    assert out.read_text().splitlines()[1].split("\t")[4] == "30.00"


def test_contig_mismatch_with_multiple_contigs_is_an_error(tmp_path):
    bed = write_bed(tmp_path / "p.bed", ["x\t0\t5\ta_1_F", "x\t5\t10\ta_1_R"])
    dep = tmp_path / "d.depths"
    dep.write_text("a\t1\tA\t1\nb\t1\tA\t1\n")
    assert (
        ac.main(["--depths", str(dep), "--bed", bed, "--sample", "S", "-o", str(tmp_path / "o")])
        == 2
    )


def test_bad_depths_is_an_error(tmp_path):
    bed = write_bed(tmp_path / "p.bed", ["ref\t0\t5\ta_1_F", "ref\t5\t10\ta_1_R"])
    dep = tmp_path / "d.depths"
    dep.write_text("ref\t1\n")
    assert (
        ac.main(["--depths", str(dep), "--bed", bed, "--sample", "S", "-o", str(tmp_path / "o")])
        == 2
    )


def test_real_scheme_end_to_end(tmp_path):
    dep = write_depths(tmp_path / "d.depths", [100] * 29903, contig="NC_045512.2")
    out = tmp_path / "out.tsv"
    bed = str(SCHEMES / "midnight.bed")
    assert (
        ac.main(
            ["--depths", dep, "--bed", bed, "--sample", "S", "--min-depth", "10", "-o", str(out)]
        )
        == 0
    )
    lines = out.read_text().splitlines()
    assert len(lines) == 30
    assert all(line.split("\t")[-1] == "false" for line in lines[1:])
