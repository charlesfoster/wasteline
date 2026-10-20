"""The probe-capture synthetic dataset (tests/data_capture, PLAN D13)."""

import gzip
import hashlib
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "tests"))

import make_test_data as mtd

CAPTURE = REPO / "tests" / "data_capture"
SEED = 42


def _hash_dir(root: Path) -> dict[str, str]:
    return {
        str(p.relative_to(root)): hashlib.sha256(p.read_bytes()).hexdigest()
        for p in sorted(root.rglob("*"))
        if p.is_file()
    }


def test_committed_capture_data_matches_regeneration(tmp_path):
    mtd.generate_capture(tmp_path, seed=SEED, log=lambda m: None)
    for rel, digest in _hash_dir(tmp_path).items():
        committed = CAPTURE / rel
        assert committed.exists(), (
            f"{rel} missing from tests/data_capture; rerun --protocol capture"
        )
        if rel.endswith(".fastq.gz"):
            assert gzip.decompress(committed.read_bytes()) == gzip.decompress(
                (tmp_path / rel).read_bytes()
            )
        else:
            assert hashlib.sha256(committed.read_bytes()).hexdigest() == digest, rel
