"""The probe-capture synthetic dataset (tests/data_capture, PLAN D13)."""

import gzip
import hashlib
import json
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "tests"))

import make_test_data as mtd

AMPLICON = REPO / "tests" / "data"
CAPTURE = REPO / "tests" / "data_capture"
SEED = 42


def _hash_dir(root: Path) -> dict[str, str]:
    return {
        str(p.relative_to(root)): hashlib.sha256(p.read_bytes()).hexdigest()
        for p in sorted(root.rglob("*"))
        if p.is_file()
    }


def _reads(path: Path):
    lines = gzip.decompress(path.read_bytes()).decode().splitlines()
    return [lines[i + 1] for i in range(0, len(lines), 4)]


@pytest.fixture(scope="module")
def gen(tmp_path_factory):
    out = tmp_path_factory.mktemp("capture")
    mtd.generate_capture(out, seed=SEED, log=lambda m: None)
    return out


def test_capture_generation_is_byte_deterministic(gen, tmp_path):
    again = tmp_path / "again"
    mtd.generate_capture(again, seed=SEED, log=lambda m: None)
    assert _hash_dir(gen) == _hash_dir(again)


def test_committed_capture_data_matches_regeneration(gen):
    for rel, digest in _hash_dir(gen).items():
        committed = CAPTURE / rel
        assert committed.exists(), (
            f"{rel} missing from tests/data_capture; rerun --protocol capture"
        )
        if rel.endswith(".fastq.gz"):
            assert gzip.decompress(committed.read_bytes()) == gzip.decompress(
                (gen / rel).read_bytes()
            )
        else:
            assert hashlib.sha256(committed.read_bytes()).hexdigest() == digest, rel


def test_generating_capture_data_never_touches_the_amplicon_dataset(tmp_path):
    before = _hash_dir(AMPLICON)
    assert mtd.main(["--protocol", "capture", "--out", str(tmp_path / "x")]) == 0
    assert _hash_dir(AMPLICON) == before
    # and the amplicon generator still reproduces tests/data exactly
    mtd.generate(tmp_path / "amp", seed=SEED, log=lambda m: None)
    for rel, digest in _hash_dir(tmp_path / "amp").items():
        if not rel.endswith(".fastq.gz"):
            assert hashlib.sha256((AMPLICON / rel).read_bytes()).hexdigest() == digest, rel


def test_capture_dataset_is_small():
    total = sum(p.stat().st_size for p in CAPTURE.rglob("*") if p.is_file())
    assert total < 3_000_000, total


def test_capture_samplesheet_and_truth():
    rows = (CAPTURE / "samplesheet.csv").read_text().splitlines()
    assert rows[0] == "sample,fastq_1,fastq_2,sample_type,collection_date,site"
    ids = [r.split(",")[0] for r in rows[1:]]
    assert ids == ["CAP_MIX", "CAP_MIX3", "CAP_PURE", "CAP_LOWCOV", "CAP_NEG"]
    assert all(r.split(",")[1].startswith("tests/data_capture/fastq/") for r in rows[1:])
    truth = json.loads((CAPTURE / "truth.json").read_text())
    assert truth["protocol"] == "capture" and "primer_scheme" not in truth
    assert not (CAPTURE / "primer_scheme.txt").exists()
    samples = truth["samples"]
    assert samples["CAP_LOWCOV"]["expected_demix_status"] == "skipped_low_coverage"
    assert (
        samples["CAP_NEG"]["sample_type"] == "negative"
        and samples["CAP_NEG"]["expected_qc_reasons"] == []
    )
    # no amplicons, so no dropout expectation anywhere
    assert all("amplicon_dropout" not in s["expected_qc_reasons"] for s in samples.values())
    for sid in ("CAP_MIX", "CAP_MIX3", "CAP_PURE"):
        assert sum(samples[sid]["true_lineage_fractions"].values()) == pytest.approx(1.0)
        lo, hi = samples[sid]["expected_duplicate_fraction_range"]
        assert lo <= samples[sid]["n_duplicate_pairs"] / samples[sid]["n_read_pairs"] <= hi


def test_capture_reads_have_no_primer_pile_up_and_contain_exact_duplicates():
    reads = _reads(CAPTURE / "fastq" / "CAP_PURE_R1.fastq.gz")
    # amplicon data: a few primer sequences start most reads; capture: starts are spread over the genome
    starts = {}
    for read in reads:
        starts[read[:20]] = starts.get(read[:20], 0) + 1
    top = sorted(starts.values(), reverse=True)[:5]
    assert top[0] < 0.05 * len(reads)
    assert len(set(reads)) < len(reads) * 0.9  # PCR duplicates are exact copies
    assert len(set(reads)) > len(reads) * 0.7
