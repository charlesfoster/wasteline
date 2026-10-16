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
FIXTURE = (
    ROOT / "tests" / "fixtures" / "consensus_multiallelic" / "synthetic_multiallelic.variants.tsv"
)
CHROM = "NC_045512.2"
IMAGE = "quay.io/biocontainers/bcftools:1.23--h3a4d415_0"
HEADER = (
    "REGION\tPOS\tREF\tALT\tREF_DP\tREF_RV\tREF_QUAL\tALT_DP\tALT_RV\tALT_QUAL\tALT_FREQ\tTOTAL_DP\t"
    "PVAL\tPASS\tGFF_FEATURE\tREF_CODON\tREF_AA\tALT_CODON\tALT_AA\tPOS_AA\n"
)

IUPAC = {
    "A": "A",
    "C": "C",
    "G": "G",
    "T": "T",
    "AG": "R",
    "CT": "Y",
    "CG": "S",
    "AT": "W",
    "GT": "K",
    "AC": "M",
    "CGT": "B",
    "AGT": "D",
    "ACT": "H",
    "ACG": "V",
    "ACGT": "N",
}


def iupac(bases):
    return IUPAC["".join(sorted(set(bases)))]


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


@pytest.fixture(scope="module")
def seq(ref):
    return ref[CHROM]


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


def alts_other_than(base, n):
    return [b for b in "ACGT" if b != base][:n]


# (pos, [(alt, af, extra kwargs)], dp, expected alleles in the consensus given the ivar rows; None = reference)
# Real examples (positions and frequencies of a wastewater run) are marked.
def cases(seq):
    def alt_of(pos, k):
        return alts_other_than(seq[pos - 1], 3)[k]

    return [
        # 22942 T>G 0.6849 + T>A 0.2941 (DP 238): was K (G/T), correct R (A/G)
        (22942, [("G", 0.6849, {}), ("A", 0.2941, {})], 238, "AG"),
        # 22896 T>G 0.576 + T>C 0.424, reference ~0 -> S
        (22896, [("G", 0.576, {}), ("C", 0.424, {})], 200, "CG"),
        # 23604 C>G 0.61 + C>A 0.39 -> R
        (23604, [("G", 0.61, {}), ("A", 0.39, {})], 200, "AG"),
        # three ALTs, all >= 0.25 (0.40 / 0.30 / 0.30, reference 0): the three alleles
        (
            4000,
            [(alt_of(4000, 0), 0.40, {}), (alt_of(4000, 1), 0.30, {}), (alt_of(4000, 2), 0.30, {})],
            200,
            "ALT3",
        ),
        # reference + two ALTs (ref 0.30, 0.40, 0.30)
        (5000, [(alt_of(5000, 0), 0.40, {}), (alt_of(5000, 1), 0.30, {})], 200, "REF+ALT2"),
        # one ALT >= 0.75 plus a small second ALT: the plain major base
        (6000, [(alt_of(6000, 0), 0.80, {}), (alt_of(6000, 1), 0.15, {})], 200, "ALT0"),
        # same, the major listed second in the TSV
        (6100, [(alt_of(6100, 1), 0.15, {}), (alt_of(6100, 0), 0.80, {})], 200, "ALT0"),
        # two ALTs below 0.25 each, reference 0.60: nothing applied
        (7000, [(alt_of(7000, 0), 0.20, {}), (alt_of(7000, 1), 0.20, {})], 200, None),
        # ALT_DP < 5 for the 0.30 allele: it is not a member; the reference (0.2, not > 0.25) is no contributor either
        # and there are two candidate rows -> the 0.5 ALT stands alone (plain base)
        (
            8000,
            [(alt_of(8000, 0), 0.50, {}), (alt_of(8000, 1), 0.30, {"alt_dp": 4})],
            200,
            "ALT0",
        ),
        # same with the reference above 0.25 (0.4 / 0.3 with 4 reads / ref 0.3)
        (
            8100,
            [(alt_of(8100, 0), 0.40, {}), (alt_of(8100, 1), 0.30, {"alt_dp": 4})],
            200,
            "REF+ALT0",
        ),
        # single-ALT behaviour (unchanged): 0.5 -> ref/alt, 0.8 -> alt, 0.74 -> ref/alt, 0.25 -> ref
        (9000, [(alt_of(9000, 0), 0.50, {})], 200, "REF+ALT0"),
        (9100, [(alt_of(9100, 0), 0.80, {})], 200, "ALT0"),
        (9200, [(alt_of(9200, 0), 0.74, {})], 200, "REF+ALT0"),
        (9300, [(alt_of(9300, 0), 0.25, {})], 200, None),
        # 0.6 + 0.15: one ALT member, f_ref = 0.25 is not above 0.25, two candidate rows -> the ALT alone
        (9400, [(alt_of(9400, 0), 0.60, {}), (alt_of(9400, 1), 0.15, {})], 200, "ALT0"),
        # 0.58 / 0.24 / 0.18, reference ~0 (a real wastewater shape): only the 0.58 ALT is a member -> plain base,
        # NOT paired with a reference that is absent from the mixture
        (
            14000,
            [
                (alt_of(14000, 2), 0.18, {}),
                (alt_of(14000, 1), 0.24, {}),
                (alt_of(14000, 0), 0.58, {}),
            ],
            200,
            "ALT0",
        ),
        # total depth below the minimum: nothing
        (10000, [(alt_of(10000, 0), 0.50, {}), (alt_of(10000, 1), 0.45, {})], 8, None),
        # one ALT with ALT_QUAL 10 is not a candidate (reference counts for the rest): ref 0.5 / ALT0 0.5
        (
            11000,
            [(alt_of(11000, 0), 0.50, {}), (alt_of(11000, 1), 0.40, {"qual": 10})],
            200,
            "REF+ALT0",
        ),
    ]


def expected_base(seq, pos, spec, iupac_on=True):
    """Consensus base for a case; `spec` as in cases()."""
    if spec is None:
        return seq[pos - 1]
    ref_base = seq[pos - 1]
    alts = alts_other_than(ref_base, 3)
    if spec == "ALT3":
        return iupac(alts)
    if spec == "REF+ALT2":
        return iupac([ref_base, alts[0], alts[1]])
    if spec == "ALT0":
        return alts[0]
    if spec == "REF+ALT0":
        return iupac([ref_base, alts[0]])
    return iupac(spec)


def case_rows(seq, case_list):
    rows = []
    for pos, alts, dp, _ in sorted(
        case_list, key=lambda c: c[0]
    ):  # the VCF must be position-sorted
        for alt, af, kw in alts:
            rows.append(row(pos, seq[pos - 1], alt, af, dp=dp, **kw))
    return rows


def test_fixture_is_in_sync(seq):
    rows = case_rows(seq, cases(seq))
    assert FIXTURE.read_text() == HEADER + "".join(rows)


def test_real_site_22942_gives_r_not_k(tmp_path, ref, seq):
    assert seq[22941] == "T"
    rows = [row(22942, "T", "G", 0.6849, dp=238), row(22942, "T", "A", 0.2941, dp=238)]
    records, cons, log = run_filter(tmp_path, rows, ref)
    assert len(records) == 1
    rec = records[0]
    assert rec[1] == "22942" and rec[3] == "T" and rec[4] == "G,A"
    assert rec[9].split(":")[0] == "1/2"
    assert cons[22941] == "R"
    assert "overlaps" not in log


def test_cases_default(tmp_path, ref, seq):
    case_list = cases(seq)
    records, cons, log = run_filter(tmp_path, case_rows(seq, case_list), ref)
    for pos, _alts, _dp, spec in case_list:
        assert cons[pos - 1] == expected_base(seq, pos, spec), (pos, spec, cons[pos - 1])
    assert "overlaps" not in log
    # one record per position, whatever the number of ALTs
    positions = [r[1] for r in records]
    assert len(positions) == len(set(positions))
    assert len(cons) == len(seq)
    # a consensus differs from the reference only at the positions with a record
    changed = {i + 1 for i, (a, b) in enumerate(zip(cons, seq)) if a != b}
    assert changed <= {int(p) for p in positions}


def test_genotypes_of_merged_records(tmp_path, ref, seq):
    case_list = cases(seq)
    records, _cons, _log = run_filter(tmp_path, case_rows(seq, case_list), ref)
    by_pos = {int(r[1]): r for r in records}
    gt = lambda pos: by_pos[pos][9].split(":")[0]
    assert (gt(22942), by_pos[22942][4].count(",")) == ("1/2", 1)
    assert gt(4000) == "1/2/3" and by_pos[4000][4].count(",") == 2
    assert gt(5000) == "0/1/2"
    assert gt(6000) == "1/1" and "," not in by_pos[6000][4]
    assert gt(6100) == "1/1" and "," not in by_pos[6100][4]
    assert gt(9000) == "0/1" and gt(9100) == "1/1"
    assert gt(14000) == "1/1" and "," not in by_pos[14000][4]
    assert gt(8000) == "1/1" and gt(9400) == "1/1"
    assert 7000 not in by_pos and 9300 not in by_pos and 10000 not in by_pos
    # ALTs are ordered by decreasing frequency
    alts4000 = by_pos[4000][4].split(",")
    assert alts4000[0] == alts_other_than(seq[3999], 3)[0]


def test_iupac_off_applies_only_a_major_allele(tmp_path, ref, seq):
    case_list = cases(seq)
    records, cons, _log = run_filter(tmp_path, case_rows(seq, case_list), ref, ext={"iupac": 0})
    for pos, _alts, _dp, spec in case_list:
        want = (
            expected_base(seq, pos, spec)
            if spec == "ALT0" and pos in (6000, 6100, 9100)
            else seq[pos - 1]
        )
        assert cons[pos - 1] == want, pos
    assert {r[1] for r in records} == {"6000", "6100", "9100"}
    assert all(r[9].startswith("1/1:") for r in records)


def test_four_member_set_gives_n(tmp_path, ref, seq):
    # needs a lower floor: reference 0.28 + three ALTs of 0.24
    pos = 12000
    alts = alts_other_than(seq[pos - 1], 3)
    rows = [row(pos, seq[pos - 1], a, 0.24, dp=200) for a in alts]
    records, cons, _ = run_filter(tmp_path, rows, ref, ext={"minaf": 0.2})
    assert records[0][9].split(":")[0] == "0/1/2/3"
    assert cons[pos - 1] == "N"


def test_consensus_freq_ext(tmp_path, ref, seq):
    # consensus_freq 0.5, two ALTs 0.55 / 0.40 -> the major alone
    pos = 13000
    alts = alts_other_than(seq[pos - 1], 2)
    rows = [row(pos, seq[pos - 1], alts[0], 0.55), row(pos, seq[pos - 1], alts[1], 0.40)]
    records, cons, _ = run_filter(tmp_path, rows, ref, ext={"cf": 0.5})
    assert records[0][4] == alts[0] and records[0][9].startswith("1/1:")
    assert cons[pos - 1] == alts[0]


def test_snv_inside_a_deletion_is_skipped_by_bcftools(tmp_path, ref, seq):
    # the 9-base deletion at 28362-28370 (ivar anchor 28361) with an SNV at 28370 inside it
    assert seq[28360:28370] == "GGAGAACGCA"
    rows = [
        row(22942, "T", "G", 0.6849, dp=238),
        row(22942, "T", "A", 0.2941, dp=238),
        row(28361, seq[28360], "-NNNNNNNNN", 0.95),
        row(28370, seq[28369], "G", 0.9),
    ]
    records, cons, log = run_filter(tmp_path, rows, ref)
    assert [r[1] for r in records] == ["22942", "28361", "28370"]
    assert cons[28361:28370] == "-" * 9  # positions 28362-28370
    assert cons[28360] == seq[28360]
    assert len(cons) == len(seq)
    assert "The site NC_045512.2:28370 overlaps with another variant" in log
    assert cons[22941] == "R"


def test_single_alt_output_is_byte_identical_to_the_previous_filter(tmp_path, ref, seq):
    """One ALT per position (the only shape before the merge): same VCF lines as the old view+awk filter."""
    import random

    rng = random.Random(11)
    rows = []
    for pos in range(1000, 1600, 3):
        r = seq[pos - 1]
        kind = rng.choice("SSSSDI")
        af = rng.choice([0.05, 0.2, 0.25, 0.26, 0.5, 0.6, 0.74, 0.75, 0.8, 0.99])
        dp = rng.choice([5, 9, 10, 40, 200])
        qual = rng.choice([10, 19, 20, 35])
        alt_dp = rng.choice([None, 3, 4, 5])
        if kind == "S":
            rows.append(
                row(pos, r, rng.choice(alts_other_than(r, 3)), af, dp=dp, qual=qual, alt_dp=alt_dp)
            )
        elif kind == "D":
            rows.append(row(pos, r, "-NN", af, dp=dp, qual=qual, alt_dp=alt_dp))
        else:
            rows.append(row(pos, r, "+AC", af, dp=dp, qual=qual, alt_dp=alt_dp))
    old_awk = (
        "/^#/ { print; next }\n"
        '{ n = split($10, f, ":"); af = f[n] + 0\n'
        '  if (length($4) == length($5) && af < cf) f[1] = "0/1"; else f[1] = "1/1"\n'
        '  $10 = f[1]; for (i = 2; i <= n; i++) $10 = $10 ":" f[i]; print }\n'
    )
    (tmp_path / "old.awk").write_text(old_awk)
    for ext in (
        {},
        {"iupac": 0},
        {"cf": 0.5},
        {"cf": 0.9},
        {"minaf": 0.1},
        {"cf": 0.5, "minaf": 0.1},
    ):
        records, _cons, _log = run_filter(tmp_path, rows, ref, ext=ext)
        new = ["\t".join(r).rstrip("\n") for r in records]
        cfg = {"cf": 0.75, "iupac": 1, "minaf": 0.25, **ext}
        het = (
            f' || (TYPE!="indel" && FORMAT/ALT_FREQ > {cfg["minaf"]} && FORMAT/ALT_DP >= 5)'
            if cfg["iupac"]
            else ""
        )
        old_expr = (
            'FORMAT/ALT_QUAL >= 20 && INFO/DP >= 10 && ((TYPE="indel" && FORMAT/ALT_FREQ >= 0.75) || '
            f'(TYPE!="indel" && FORMAT/ALT_FREQ >= {cfg["cf"]}){het})'
        )
        script = (
            "bcftools norm -a -m - -Ou S.ivar.vcf.gz | "
            f"bcftools view -i '{old_expr}' -Ov | "
            f"awk -F'\\t' -v OFS='\\t' -v cf={cfg['cf']} -f old.awk | bcftools view -H > S.old.txt"
        )
        r = subprocess.run(
            docker_cmd(tmp_path, script),
            capture_output=True,
            text=True,
            check=False,
        )
        assert r.returncode == 0, r.stderr
        old = (tmp_path / "S.old.txt").read_text().splitlines()
        assert new == old, ext


def test_absent_reference_rule_with_other_thresholds(tmp_path, ref, seq):
    # 0.58 / 0.24 / 0.18 at a site; the reference is absent (f_ref = 0)
    pos = 14000
    alts = alts_other_than(seq[pos - 1], 3)
    rows = [row(pos, seq[pos - 1], a, af) for a, af in zip(alts, (0.58, 0.24, 0.18), strict=True)]
    # default: only the 0.58 ALT is a member -> plain base, GT 1/1
    records, cons, _ = run_filter(tmp_path, rows, ref)
    assert (
        records[0][4] == alts[0] and records[0][9].startswith("1/1:") and cons[pos - 1] == alts[0]
    )
    # iupac_min_af 0.1: all three ALTs are members (reference absent) -> one three-base code
    records, cons, _ = run_filter(tmp_path, rows, ref, ext={"minaf": 0.1})
    assert records[0][4] == ",".join(alts) and records[0][9].startswith("1/2/3:")
    assert cons[pos - 1] == iupac(alts)
    # consensus_freq 0.5: the 0.58 ALT is a major allele, written alone
    records, cons, _ = run_filter(tmp_path, rows, ref, ext={"cf": 0.5})
    assert records[0][4] == alts[0] and cons[pos - 1] == alts[0]
    # a lone candidate row keeps the single-ALT rule even when the reference is below iupac_min_af:
    # consensus_freq 0.9, one ALT 0.8 (f_ref 0.2) -> ALT/reference as before
    records, cons, _ = run_filter(
        tmp_path, [row(pos, seq[pos - 1], alts[0], 0.8)], ref, ext={"cf": 0.9}
    )
    assert records[0][9].startswith("0/1:") and cons[pos - 1] == iupac([seq[pos - 1], alts[0]])
