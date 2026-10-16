"""bin/ivar_variants_to_vcf.py: indel alleles must be built from the REFERENCE, not from ivar's printed bases.

Regression: ivar prints the deleted bases of a deletion as N (`G  -NNNNNNNNN`). The old converter built
REF=GNNNNNNNNN, which does not match the reference (GTCTGGTTTT), so `bcftools consensus` aborted and the whole
genome became N.
"""

import random
import shutil
import subprocess
from pathlib import Path

import ivar_variants_to_vcf as conv
import pytest

ROOT = Path(__file__).resolve().parents[2]
REFERENCE = ROOT / "assets" / "references" / "NC_045512.fasta"
FIXTURE = ROOT / "tests" / "fixtures" / "modules_phase3" / "synthetic_ivar_indels.variants.tsv"
CHROM = "NC_045512.2"
HEADER = (
    "REGION\tPOS\tREF\tALT\tREF_DP\tREF_RV\tREF_QUAL\tALT_DP\tALT_RV\tALT_QUAL\tALT_FREQ\tTOTAL_DP\t"
    "PVAL\tPASS\tGFF_FEATURE\tREF_CODON\tREF_AA\tALT_CODON\tALT_AA\tPOS_AA\n"
)


@pytest.fixture(scope="module")
def ref():
    return conv.read_reference(str(REFERENCE))


@pytest.fixture(scope="module")
def seq(ref):
    return ref[CHROM]


def tsv_row(pos, ref_base, alt, af=0.9, dp=60, passed="TRUE", chrom=CHROM):
    alt_dp = round(dp * af)
    return (
        f"{chrom}\t{pos}\t{ref_base}\t{alt}\t{dp - alt_dp}\t0\t37\t{alt_dp}\t0\t35\t{af}\t{dp}\t0.01\t"
        f"{passed}\tNA\tNA\tNA\tNA\tNA\tNA\n"
    )


def convert(tmp_path, rows, reference=None, **kwargs):
    """Run the converter on TSV `rows`; return the VCF records as lists of fields."""
    tsv = tmp_path / "S.variants.tsv"
    tsv.write_text(HEADER + "".join(rows))
    out = tmp_path / "out" / "S.vcf"
    conv.ivar_variants_to_vcf(str(tsv), str(out), reference=reference, **kwargs)
    return [line.split("\t") for line in out.read_text().splitlines() if not line.startswith("#")]


def test_real_regression_nine_base_deletion_printed_as_n(tmp_path, ref, seq):
    # the exact row from the failing wastewater sample
    assert seq[11286:11296] == "GTCTGGTTTT"
    (rec,) = convert(tmp_path, [tsv_row(11287, "G", "-NNNNNNNNN", af=0.988636, dp=176)], ref)
    assert rec[:5] == [CHROM, "11287", ".", "GTCTGGTTTT", "G"]
    assert rec[6] == "PASS" and rec[7] == "DP=176"
    assert rec[8] == "GT:REF_DP:REF_RV:REF_QUAL:ALT_DP:ALT_RV:ALT_QUAL:ALT_FREQ"
    assert rec[9].rstrip().endswith(":0.988636")


def test_legacy_behaviour_without_reference_is_the_bug(tmp_path, capsys):
    (rec,) = convert(tmp_path, [tsv_row(11287, "G", "-NNNNNNNNN")], None)
    assert rec[3:5] == ["GNNNNNNNNN", "G"]
    assert "no --reference" in capsys.readouterr().err


def test_snv_unchanged(tmp_path, ref, seq):
    alt = "A" if seq[4999] != "A" else "C"
    (rec,) = convert(tmp_path, [tsv_row(5000, seq[4999], alt)], ref)
    assert rec[1:5] == ["5000", ".", seq[4999], alt]


def test_snv_ref_mismatch_is_kept_but_counted(tmp_path, ref, seq, capsys):
    wrong = "A" if seq[99] != "A" else "C"
    (rec,) = convert(tmp_path, [tsv_row(100, wrong, "T")], ref)
    assert rec[3:5] == [wrong, "T"]  # left for CONSENSUS_BCFTOOLS to skip with a warning
    assert "1 record(s): ref_mismatch" in capsys.readouterr().err


def test_three_base_deletion_with_real_bases(tmp_path, ref, seq):
    printed = seq[20000:20003]
    (rec,) = convert(tmp_path, [tsv_row(20000, seq[19999], "-" + printed)], ref)
    assert rec[1:5] == ["20000", ".", seq[19999:20003], seq[19999]]
    assert len(rec[3]) == 4


def test_deletion_ignores_printed_bases_and_uses_only_their_number(tmp_path, ref, seq):
    # printed bases that disagree with the reference: only the length (3) matters
    (rec,) = convert(tmp_path, [tsv_row(20000, seq[19999], "-GGG")], ref)
    assert rec[3:5] == [seq[19999:20003], seq[19999]]


def test_insertion(tmp_path, ref, seq):
    (rec,) = convert(tmp_path, [tsv_row(15000, seq[14999], "+AC")], ref)
    assert rec[1:5] == ["15000", ".", seq[14999], seq[14999] + "AC"]


def test_insertion_containing_n_is_skipped_and_counted(tmp_path, ref, seq, capsys):
    recs = convert(
        tmp_path,
        [
            tsv_row(15000, seq[14999], "+ANT"),
            tsv_row(16000, seq[15999], "+N"),
            tsv_row(17000, seq[16999], "+G"),
        ],
        ref,
    )
    assert [r[1] for r in recs] == ["17000"]
    assert "2 record(s): insertion_non_acgt" in capsys.readouterr().err


def test_lowercase_insertion_is_uppercased(tmp_path, ref, seq):
    (rec,) = convert(tmp_path, [tsv_row(15000, seq[14999], "+ac")], ref)
    assert rec[4] == seq[14999] + "AC"


def test_adjacent_snv_and_deletion_at_the_same_anchor(tmp_path, ref, seq):
    snv_alt = "A" if seq[23999] != "A" else "C"
    recs = convert(
        tmp_path,
        [tsv_row(24000, seq[23999], snv_alt, af=0.2), tsv_row(24000, seq[23999], "-NN", af=0.8)],
        ref,
    )
    assert [r[1:5] for r in recs] == [
        ["24000", ".", seq[23999], snv_alt],
        ["24000", ".", seq[23999:24002], seq[23999]],
    ]
    assert len(recs[1][3]) == 3


def test_deletion_running_exactly_to_the_genome_end(tmp_path, ref, seq):
    assert len(seq) == 29903
    (rec,) = convert(tmp_path, [tsv_row(29900, seq[29899], "-NNN")], ref)
    assert rec[3:5] == [seq[29899:], seq[29899]]
    assert len(rec[3]) == 4 and seq[29899:] == "AAAA"


def test_deletion_overrunning_the_genome_end_is_truncated(tmp_path, ref, seq, capsys):
    (rec,) = convert(tmp_path, [tsv_row(29902, seq[29901], "-NNNNN")], ref)
    assert rec[3:5] == [seq[29901:], seq[29901]]
    assert len(rec[3]) == 2
    assert "deletion_truncated" in capsys.readouterr().err


def test_deletion_with_nothing_left_to_delete_is_skipped(tmp_path, ref, seq, capsys):
    assert convert(tmp_path, [tsv_row(29903, seq[29902], "-N")], ref) == []
    assert "deletion_past_contig_end" in capsys.readouterr().err


def test_unknown_contig_indel_is_skipped(tmp_path, ref, capsys):
    assert convert(tmp_path, [tsv_row(100, "A", "-N", chrom="chrX")], ref) == []
    assert "unknown_contig" in capsys.readouterr().err


def test_anchor_mismatch_uses_the_reference_anchor_and_is_counted(tmp_path, ref, seq, capsys):
    wrong = "A" if seq[11286] != "A" else "C"
    (rec,) = convert(tmp_path, [tsv_row(11287, wrong, "-NNNNNNNNN")], ref)
    assert rec[3:5] == ["GTCTGGTTTT", "G"]
    assert "ref_mismatch" in capsys.readouterr().err


def test_duplicates_collapse_and_pass_only_and_min_af_still_work(tmp_path, ref, seq):
    row = tsv_row(11287, "G", "-NNNNNNNNN", af=0.9)
    assert len(convert(tmp_path, [row, row], ref)) == 1
    fail = tsv_row(11287, "G", "-NNNNNNNNN", passed="FALSE")
    assert convert(tmp_path, [fail], ref, passOnly=True) == []
    assert convert(tmp_path, [tsv_row(11287, "G", "-NNNNNNNNN", af=0.1)], ref, minAF=0.5) == []


def test_14_column_tsv_without_trailing_columns_keeps_pass_flag(tmp_path, ref, seq):
    tsv = tmp_path / "S.variants.tsv"
    tsv.write_text(
        HEADER + f"{CHROM}\t5000\t{seq[4999]}\tA\t1\t0\t37\t9\t0\t35\t0.9\t10\t0.01\tTRUE\n"
    )
    out = tmp_path / "S.vcf"
    conv.ivar_variants_to_vcf(str(tsv), str(out), reference=ref)
    body = [line for line in out.read_text().splitlines() if not line.startswith("#")]
    assert body[0].split("\t")[6] == "PASS"


@pytest.mark.parametrize("kind", ["del_n", "del_real", "ins"])
def test_every_record_matches_the_reference_for_all_lengths(tmp_path, ref, seq, kind):
    """Property: REF == reference at POS for every indel length 1-12 at random positions."""
    rng = random.Random(20260101)
    rows = []
    for length in range(1, 13):
        for _ in range(25):
            pos = rng.randint(1, len(seq) - length - 1)
            anchor = seq[pos - 1]
            if kind == "del_n":
                alt = "-" + "N" * length
            elif kind == "del_real":
                alt = "-" + seq[pos : pos + length]
            else:
                alt = "+" + "".join(rng.choice("ACGT") for _ in range(length))
            rows.append(tsv_row(pos, anchor, alt))
    recs = convert(tmp_path, rows, ref)
    assert len(recs) == len(
        rows
    )  # (random positions: duplicates are astronomically unlikely and fixed by the seed)
    for rec in recs:
        pos, rbase, alt = int(rec[1]), rec[3], rec[4]
        assert seq[pos - 1 : pos - 1 + len(rbase)] == rbase
        if kind == "ins":
            assert rbase == seq[pos - 1] and alt.startswith(rbase) and len(alt) > 1
        else:
            assert alt == seq[pos - 1] and len(rbase) > 1


def test_header_and_sample_column_unchanged(tmp_path, ref):
    tsv = tmp_path / "S.variants.tsv"
    tsv.write_text(HEADER)
    out = tmp_path / "S.vcf"
    conv.ivar_variants_to_vcf(str(tsv), str(out), reference=ref)
    lines = out.read_text().splitlines()
    assert lines[0] == "##fileformat=VCFv4.2"
    assert lines[-1].split("\t")[-1] == str(
        tmp_path / "S.variants"
    )  # the module renames it to the sample ID
    assert len(lines) == len([x for x in lines if x.startswith("#")])


def test_read_reference_multi_contig_and_errors(tmp_path):
    fa = tmp_path / "r.fa"
    fa.write_text(">c1 desc\nacgt\nAC\n>c2\nGG\n")
    assert conv.read_reference(str(fa)) == {"c1": "ACGTAC", "c2": "GG"}
    dup = tmp_path / "d.fa"
    dup.write_text(">c1\nA\n>c1\nC\n")
    with pytest.raises(ValueError):
        conv.read_reference(str(dup))


def test_cli_with_reference(tmp_path):
    out = tmp_path / "o" / "x.vcf"
    rc = conv.main(["--reference", str(REFERENCE), str(FIXTURE), str(out)])
    assert rc is None or rc == 0
    recs = [line.split("\t") for line in out.read_text().splitlines() if not line.startswith("#")]
    by_pos = {}
    for r in recs:
        by_pos.setdefault(r[1], []).append(r)
    assert by_pos["11287"][0][3:5] == ["GTCTGGTTTT", "G"]
    assert "16000" not in by_pos  # insertion containing N
    assert len(by_pos["24000"]) == 2
    assert by_pos["29900"][0][3] == "AAAA"


def test_fixture_records_all_match_the_reference(tmp_path, ref, seq):
    out = tmp_path / "x.vcf"
    conv.ivar_variants_to_vcf(str(FIXTURE), str(out), reference=ref)
    for line in out.read_text().splitlines():
        if line.startswith("#"):
            continue
        f = line.split("\t")
        assert seq[int(f[1]) - 1 : int(f[1]) - 1 + len(f[3])] == f[3]


@pytest.mark.skipif(shutil.which("bcftools") is None, reason="bcftools not on PATH")
def test_bcftools_consensus_accepts_the_vcf_and_marks_the_deletion(tmp_path, ref, seq):
    """End to end on the failing row shape: 9 '-' at 11288-11296, length kept."""
    out = tmp_path / "x.vcf"
    conv.ivar_variants_to_vcf(str(FIXTURE), str(out), reference=ref)
    lines = out.read_text().splitlines()
    header = [line for line in lines if line.startswith("#")]
    header.insert(1, f"##contig=<ID={CHROM},length={len(seq)}>")
    body = [line for line in lines if not line.startswith("#")]
    keep = [line for line in body if line.split("\t")[1] in ("11287", "29900")]
    # GT 1/1 as written by CONSENSUS_FILTER_VCF
    keep = [line.split("\t")[:9] + ["1/1"] for line in keep]
    vcf = tmp_path / "f.vcf"
    vcf.write_text("\n".join(header) + "\n" + "\n".join("\t".join(k) for k in keep) + "\n")
    subprocess.run(
        ["bcftools", "view", "-Oz", "-o", str(tmp_path / "f.vcf.gz"), str(vcf)], check=True
    )
    subprocess.run(["bcftools", "index", "-t", str(tmp_path / "f.vcf.gz")], check=True)
    cons = subprocess.run(
        [
            "bcftools",
            "consensus",
            "-f",
            str(REFERENCE),
            "--mark-del",
            "-",
            "-H",
            "I",
            str(tmp_path / "f.vcf.gz"),
        ],
        check=True,
        capture_output=True,
        text=True,
    ).stdout
    s = "".join(line for line in cons.splitlines() if not line.startswith(">"))
    assert len(s) == 29903
    assert s[11287:11296] == "-" * 9 and s[11286] == "G" and s[11296] == seq[11296]
    assert s[29900:] == "---"


def test_competing_snv_alts_at_one_position_stay_separate_records(tmp_path, ref, seq):
    # the converter keeps one record per ALT (frequencies intact); CONSENSUS_FILTER_VCF merges them
    # (tests/python/test_consensus_multiallelic.py)
    assert seq[22941] == "T"
    recs = convert(
        tmp_path,
        [tsv_row(22942, "T", "G", af=0.6849, dp=238), tsv_row(22942, "T", "A", af=0.2941, dp=238)],
        ref,
    )
    assert [r[1:5] for r in recs] == [["22942", ".", "T", "G"], ["22942", ".", "T", "A"]]
    assert [r[9].rstrip().split(":")[-1] for r in recs] == ["0.6849", "0.2941"]
