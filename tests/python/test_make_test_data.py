import csv
import gzip
import hashlib
import json
import re
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "tests"))

import make_test_data as mtd

COMMITTED = REPO / "tests" / "data"
SEED = 42


def _hash_dir(root: Path) -> dict[str, str]:
    return {
        str(p.relative_to(root)): hashlib.sha256(p.read_bytes()).hexdigest()
        for p in sorted(root.rglob("*"))
        if p.is_file()
    }


@pytest.fixture(scope="module")
def gen_a(tmp_path_factory):
    out = tmp_path_factory.mktemp("gen_a")
    mtd.generate(out, seed=SEED, log=lambda m: None)
    return out


def test_generation_is_byte_deterministic(gen_a, tmp_path):
    out_b = tmp_path / "b"
    mtd.generate(out_b, seed=SEED, log=lambda m: None)
    assert _hash_dir(gen_a) == _hash_dir(out_b)


def test_different_seed_changes_reads(gen_a, tmp_path):
    out_c = tmp_path / "c"
    mtd.generate(out_c, seed=SEED + 1, log=lambda m: None)
    h_a, h_c = _hash_dir(gen_a), _hash_dir(out_c)
    assert h_a["fastq/MIX_A_R1.fastq.gz"] != h_c["fastq/MIX_A_R1.fastq.gz"]


def test_gzip_header_is_deterministic(gen_a):
    raw = (gen_a / "fastq" / "MIX_A_R1.fastq.gz").read_bytes()
    assert raw[:2] == b"\x1f\x8b"
    assert raw[4:8] == b"\x00\x00\x00\x00"  # mtime = 0
    assert raw[3] == 0  # no FNAME flag


def test_committed_data_matches_regeneration(gen_a):
    """Committed tests/data must be what the generator produces (FASTQ compared decompressed so that a different
    zlib build does not cause a spurious failure)."""
    for rel, digest in _hash_dir(gen_a).items():
        committed = COMMITTED / rel
        assert committed.exists(), f"{rel} missing from tests/data; rerun tests/make_test_data.py"
        if rel.endswith(".fastq.gz"):
            assert gzip.decompress(committed.read_bytes()) == gzip.decompress(
                (gen_a / rel).read_bytes()
            ), rel
        else:
            assert hashlib.sha256(committed.read_bytes()).hexdigest() == digest, rel


def test_dataset_is_small():
    total = sum(p.stat().st_size for p in COMMITTED.rglob("*") if p.is_file())
    assert total < 5 * 1024 * 1024


def test_samplesheet():
    with open(COMMITTED / "samplesheet.csv", newline="") as fh:
        rows = list(csv.DictReader(fh))
    assert list(rows[0]) == [
        "sample",
        "fastq_1",
        "fastq_2",
        "sample_type",
        "collection_date",
        "site",
    ]
    ids = [r["sample"] for r in rows]
    assert len(ids) == len(set(ids))
    assert {"MIX_A", "MIX_B", "PURE_C", "LOWCOV_D", "NEG_E", "NEG_CONTAM", "Site_1_week_2"} <= set(
        ids
    )
    by_id = {r["sample"]: r for r in rows}
    assert (
        by_id["NEG_E"]["sample_type"] == "negative"
        and by_id["NEG_CONTAM"]["sample_type"] == "negative"
    )
    assert all(r["sample_type"] == "sample" for r in rows if not r["sample"].startswith("NEG_"))
    for r in rows:
        assert re.fullmatch(r"[A-Za-z0-9._-]+", r["sample"])
        assert re.fullmatch(r"\d{4}-\d{2}-\d{2}", r["collection_date"])
        for col in ("fastq_1", "fastq_2"):
            assert r[col].startswith("tests/data/fastq/")
            assert (REPO / r[col]).stat().st_size > 0


def _fastq_records(path: Path):
    lines = gzip.decompress(path.read_bytes()).decode().splitlines()
    assert len(lines) % 4 == 0
    return [tuple(lines[i : i + 4]) for i in range(0, len(lines), 4)]


def test_fastq_pairs_are_valid_and_nonempty():
    truth = json.loads((COMMITTED / "truth.json").read_text())
    for sample, info in truth["samples"].items():
        r1 = _fastq_records(COMMITTED / "fastq" / f"{sample}_R1.fastq.gz")
        r2 = _fastq_records(COMMITTED / "fastq" / f"{sample}_R2.fastq.gz")
        assert len(r1) == len(r2) == info["n_read_pairs"] > 0
        for a, b in zip(r1, r2, strict=True):
            assert a[0].split()[0] == b[0].split()[0]
            assert a[2] == "+" and len(a[1]) == len(a[3])
            assert set(a[1]) <= set("ACGT")
    assert truth["samples"]["NEG_E"]["n_read_pairs"] == 3


def test_truth_json():
    truth = json.loads((COMMITTED / "truth.json").read_text())
    s = truth["samples"]
    for info in s.values():
        assert sum(info["true_lineage_fractions"].values()) == pytest.approx(1.0)
        assert info["expected_demix_status"] in {"ok", "skipped_low_coverage"}
    assert s["LOWCOV_D"]["expected_demix_status"] == "skipped_low_coverage"
    assert s["NEG_E"]["expected_demix_status"] == "skipped_low_coverage"
    assert s["MIX_A"]["expected_demix_status"] == "ok"
    assert s["Site_1_week_2"]["dropout_amplicons"]
    assert s["Site_1_week_2"]["expected_zero_coverage_regions_0based_half_open"]
    assert "amplicon_dropout" in s["Site_1_week_2"]["expected_qc_reasons"]
    assert "low_coverage" in s["LOWCOV_D"]["expected_qc_reasons"]
    # every truth lineage exists in the committed barcodes
    barcodes = mtd.read_barcodes(COMMITTED / "barcodes_test.csv")
    for info in s.values():
        assert set(info["true_lineage_fractions"]) <= set(barcodes)


def test_barcodes_test_csv_is_valid_freyja_format():
    ref = mtd.read_fasta(mtd.DEFAULT_REFERENCE)
    with open(COMMITTED / "barcodes_test.csv", newline="") as fh:
        reader = csv.reader(fh)
        header = next(reader)
        rows = list(reader)
    assert header[0] == ""
    assert 12 <= len(rows) <= 20
    for m in header[1:]:
        mt = re.fullmatch(r"([ACGT])(\d+)([ACGT])", m)
        assert mt, m
        assert ref[int(mt.group(2)) - 1] == mt.group(1)
    # every kept column is non-zero for at least one lineage
    for j in range(1, len(header)):
        assert any(float(r[j]) > 0 for r in rows), header[j]
    lineages = [r[0] for r in rows]
    assert len(set(lineages)) == len(lineages)
    assert {"BA.1", "BA.2", "BA.5", "XBB.1.5", "JN.1", "AY.4"} <= set(lineages)
    # distinct barcodes
    assert len({tuple(r[1:]) for r in rows}) == len(rows)


def test_apply_mutations_counts_skips():
    ref = "ACGTACGTAC"
    genome, stats = mtd.apply_mutations(
        ref, ["A1T", "C2G", "G9T", "-5-3", "+3ACG", "C99A", "T4A", "T4C"]
    )
    assert genome[0] == "T" and genome[3] == "A"
    assert stats == {
        "applied": 3,
        "skipped_ref_mismatch": 2,
        "skipped_unsupported": 2,
        "skipped_conflict": 1,
    }


def test_lineage_genome_differs_by_barcode_mutations():
    ref = mtd.read_fasta(mtd.DEFAULT_REFERENCE)
    muts = mtd.read_barcodes(COMMITTED / "barcodes_test.csv")["BA.2"]
    genome, stats = mtd.apply_mutations(ref, muts)
    assert len(genome) == len(ref)
    assert sum(a != b for a, b in zip(ref, genome, strict=True)) == stats["applied"] == len(muts)


@pytest.mark.parametrize(("scheme", "min_amplicons"), [("midnight", 29), ("articv4.1", 90)])
def test_load_amplicons(scheme, min_amplicons):
    ref = mtd.read_fasta(mtd.DEFAULT_REFERENCE)
    amps = mtd.load_amplicons(mtd.DEFAULT_PRIMER_DIR / f"{scheme}.bed", ref)
    assert len(amps) >= min_amplicons
    for a in amps:
        assert a["end"] > a["start"]
        assert a["fwd_primer"] == ref[a["start"] : a["start"] + a["fwd_len"]]


def test_reads_start_with_primers():
    ref = mtd.read_fasta(mtd.DEFAULT_REFERENCE)
    amps = mtd.load_amplicons(mtd.DEFAULT_PRIMER_DIR / "midnight.bed", ref)
    fwd = {a["fwd_primer"][:15] for a in amps}
    rev = {a["rev_primer"][:15] for a in amps}
    r1 = _fastq_records(COMMITTED / "fastq" / "PURE_C_R1.fastq.gz")
    r2 = _fastq_records(COMMITTED / "fastq" / "PURE_C_R2.fastq.gz")
    n = sum(
        1 for a, b in zip(r1, r2, strict=True) if a[1][:15] in fwd | rev or b[1][:15] in fwd | rev
    )
    assert n / len(r1) > 0.3


def test_dropout_amplicons_get_no_reads(gen_a):
    """Reads of Site_1_week_2 never fall inside the genome regions only covered by dropped amplicons."""
    truth = json.loads((gen_a / "truth.json").read_text())
    regions = truth["samples"]["Site_1_week_2"]["expected_zero_coverage_regions_0based_half_open"]
    assert regions and all(e > s for s, e in regions)
    assert truth["samples"]["MIX_A"]["expected_zero_coverage_regions_0based_half_open"] == []
