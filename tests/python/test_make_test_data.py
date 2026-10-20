import csv
import gzip
import hashlib
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


def test_samplesheet_points_at_existing_fastqs():
    with open(COMMITTED / "samplesheet.csv", newline="") as fh:
        rows = list(csv.DictReader(fh))
    by_id = {r["sample"]: r for r in rows}
    assert {"MIX_A", "MIX_B", "PURE_C", "LOWCOV_D", "NEG_E", "NEG_CONTAM", "Site_1_week_2"} <= set(
        by_id
    )
    assert by_id["NEG_E"]["sample_type"] == "negative"
    for r in rows:
        for col in ("fastq_1", "fastq_2"):
            assert (REPO / r[col]).stat().st_size > 0
