"""bin/ivar_variants_to_vcf.py: indel alleles must be built from the REFERENCE, not from ivar's printed bases.

Regression: ivar prints the deleted bases of a deletion as N (`G  -NNNNNNNNN`). The old converter built
REF=GNNNNNNNNN, which does not match the reference (GTCTGGTTTT), so `bcftools consensus` aborted and the whole
genome became N.
"""

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


def test_regression_nine_base_deletion_printed_as_n_uses_reference_bases(tmp_path, ref, seq):
    # the exact row from the failing wastewater sample
    assert seq[11286:11296] == "GTCTGGTTTT"
    (rec,) = convert(tmp_path, [tsv_row(11287, "G", "-NNNNNNNNN", af=0.988636, dp=176)], ref)
    assert rec[:5] == [CHROM, "11287", ".", "GTCTGGTTTT", "G"]
    assert rec[6] == "PASS" and rec[7] == "DP=176"
    assert rec[8] == "GT:REF_DP:REF_RV:REF_QUAL:ALT_DP:ALT_RV:ALT_QUAL:ALT_FREQ"
    assert rec[9].rstrip().endswith(":0.988636")


def test_cli_fixture_indels_snv_and_skipped_insertion(tmp_path, seq):
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
    for r in recs:  # every REF matches the reference
        assert seq[int(r[1]) - 1 : int(r[1]) - 1 + len(r[3])] == r[3]


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
