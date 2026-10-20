"""CONSENSUS_FILTER_VCF: several ALT alleles at one position are merged into ONE record.

ivar writes one row per ALT. `bcftools consensus` applies only one record per position, so before the merge a
mixed site such as 22942 (T>G 0.68, T>A 0.29) became `K` (G/T) instead of `R` (A/G). These tests run the real
shell pipeline of modules/local/consensus_filter_vcf/main.nf (the `bcftools view -i` expression and the awk
program are extracted from main.nf, so they cannot drift) followed by `bcftools consensus -H I --mark-del -`,
all inside the pinned bcftools container. Skipped when docker or the image is not available.
"""

import re
import shutil
import subprocess
from pathlib import Path

import ivar_variants_to_vcf as conv
import pytest

ROOT = Path(__file__).resolve().parents[2]
REFERENCE = ROOT / "assets" / "references" / "NC_045512.fasta"
MODULE = ROOT / "modules" / "local" / "consensus_filter_vcf" / "main.nf"
CHROM = "NC_045512.2"
IMAGE = "quay.io/biocontainers/bcftools:1.23--h3a4d415_0"
HEADER = (
    "REGION\tPOS\tREF\tALT\tREF_DP\tREF_RV\tREF_QUAL\tALT_DP\tALT_RV\tALT_QUAL\tALT_FREQ\tTOTAL_DP\t"
    "PVAL\tPASS\tGFF_FEATURE\tREF_CODON\tREF_AA\tALT_CODON\tALT_AA\tPOS_AA\n"
)


def _docker_ready():
    if shutil.which("docker") is None:
        return False
    try:
        r = subprocess.run(
            ["docker", "image", "inspect", IMAGE], capture_output=True, timeout=60, check=False
        )
        return r.returncode == 0
    except (OSError, subprocess.TimeoutExpired):
        return False


pytestmark = pytest.mark.skipif(not _docker_ready(), reason=f"docker image {IMAGE} not available")


@pytest.fixture(scope="module")
def ref():
    return conv.read_reference(str(REFERENCE))


def docker_cmd(workdir, script):
    return [
        "docker",
        "run",
        "--rm",
        "--platform",
        "linux/amd64",
        "-v",
        f"{workdir}:/w",
        "-w",
        "/w",
        IMAGE,
        "sh",
        "-c",
        script,
    ]


def module_parts():
    """(view expression template, awk program) taken from the module script."""
    text = MODULE.read_text()
    m = re.search(r"# BEGIN merge_snv_alleles\n(.*?)\n\s*# END merge_snv_alleles", text, re.DOTALL)
    awk = m.group(1).replace("\\$", "$")
    expr = re.search(r"bcftools view -i '([^']*)'", text).group(1)
    return expr, awk


def row(pos, ref_base, alt, af, dp=200, alt_dp=None, qual=35):
    alt_dp = round(dp * af) if alt_dp is None else alt_dp
    return (
        f"{CHROM}\t{pos}\t{ref_base}\t{alt}\t{dp - alt_dp}\t0\t37\t{alt_dp}\t0\t{qual}\t{af}\t{dp}\t0.01\t"
        "TRUE\tNA\tNA\tNA\tNA\tNA\tNA\n"
    )


def run_filter(tmp_path, rows, ref, ext=None):
    """TSV rows -> converter -> module filter pipeline -> (records of the final VCF, consensus, stderr)."""
    ext = {
        "cf": 0.75,
        "indelf": 0.75,
        "iupac": 1,
        "minaf": 0.25,
        "minreads": 5,
        "minqual": 20,
        "mindepth": 10,
        **(ext or {}),
    }
    tsv = tmp_path / "S.variants.tsv"
    tsv.write_text(HEADER + "".join(rows))
    conv.ivar_variants_to_vcf(str(tsv), str(tmp_path / "S.ivar.vcf"), reference=ref)
    expr, awk = module_parts()
    expr = expr.replace("${min_alt_qual}", str(ext["minqual"])).replace(
        "${min_depth}", str(ext["mindepth"])
    )
    (tmp_path / "prog.awk").write_text(awk + "\n")
    (tmp_path / "ref.fa").write_text(REFERENCE.read_text())
    script = f"""set -e
bgzip -c S.ivar.vcf > S.ivar.vcf.gz && bcftools index -f S.ivar.vcf.gz
bcftools norm -a -m - -Ou S.ivar.vcf.gz | bcftools view -i '{expr}' -Ov | \\
  awk -F'\\t' -v OFS='\\t' -v cf={ext["cf"]} -v indelf={ext["indelf"]} -v iupac={ext["iupac"]} \\
      -v minaf={ext["minaf"]} -v minreads={ext["minreads"]} -f prog.awk | bgzip -c > S.cons.vcf.gz
bcftools index -t -f S.cons.vcf.gz
bcftools view -H S.cons.vcf.gz > S.records.txt
bcftools consensus -f ref.fa --mark-del '-' -H I S.cons.vcf.gz > S.cons.fa 2> S.cons.log
"""
    r = subprocess.run(
        docker_cmd(tmp_path, script),
        capture_output=True,
        text=True,
        check=False,
    )
    assert r.returncode == 0, r.stderr
    records = [line.split("\t") for line in (tmp_path / "S.records.txt").read_text().splitlines()]
    cons = "".join(
        line
        for line in (tmp_path / "S.cons.fa").read_text().splitlines()
        if not line.startswith(">")
    )
    return records, cons, (tmp_path / "S.cons.log").read_text()


def test_real_site_22942_gives_r_not_k(tmp_path, ref):
    seq = ref[CHROM]
    assert seq[22941] == "T"
    rows = [row(22942, "T", "G", 0.6849, dp=238), row(22942, "T", "A", 0.2941, dp=238)]
    records, cons, log = run_filter(tmp_path, rows, ref)
    assert len(records) == 1
    rec = records[0]
    assert rec[1] == "22942" and rec[3] == "T" and rec[4] == "G,A"
    assert rec[9].split(":")[0] == "1/2"
    assert cons[22941] == "R"
    assert "overlaps" not in log
