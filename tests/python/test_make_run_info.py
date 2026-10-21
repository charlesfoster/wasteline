import hashlib
from datetime import UTC, datetime

import make_run_info as mri
from freyja_dataset_info import FREYJA_DATA_KEYS

DATASET = {key: f"v-{key}" for key in FREYJA_DATA_KEYS}

VERSIONS = """\
"WASTELINE:WORKFLOW:FREYJA_VARIANTS":
  freyja: 2.0.5
"WASTELINE:WORKFLOW:FASTP":
  fastp: 0.24.0
Workflow:
  Nextflow: 26.04.6
"""


def test_build_run_info(tmp_path):
    files = {}
    for name, content in {"s.csv": b"sample\n", "p.bed": b"bed\n", "r.fasta": b">r\nA\n"}.items():
        (tmp_path / name).write_bytes(content)
        files[name] = str(tmp_path / name)
    base = {"pipeline_version": "0.1.0dev", "mode": "full", "params": {"a": 1}}
    now = datetime(2026, 10, 15, 1, 2, 3, tzinfo=UTC)
    info = mri.build_run_info(
        base, files["s.csv"], files["p.bed"], files["r.fasta"], dict(DATASET), VERSIONS, now=now
    )
    assert info["pipeline_version"] == "0.1.0dev"
    assert info["freyja_version"] == "2.0.5"
    for key in FREYJA_DATA_KEYS:
        assert key in info["freyja_data"], key
    assert info["input_checksums"]["primer_bed"]["sha256"] == hashlib.sha256(b"bed\n").hexdigest()
    assert info["provenance_generated"] == "2026-10-15T01:02:03Z"


def test_several_primer_beds_are_each_checksummed(tmp_path):
    """PLAN D20: samples with different primer schemes -> primer_beds (one entry per BED), no single primer_bed."""
    paths = []
    for name, content in {"s.csv": b"sample\n", "a.bed": b"a\n", "b.bed": b"b\n"}.items():
        (tmp_path / name).write_bytes(content)
        paths.append(str(tmp_path / name))
    info = mri.build_run_info(
        {"mode": "full"},
        paths[0],
        paths[1:],
        None,
        dict(DATASET),
        VERSIONS,
        now=datetime(2026, 1, 1, tzinfo=UTC),
    )
    checksums = info["input_checksums"]
    assert "primer_bed" not in checksums
    assert [e["file"] for e in checksums["primer_beds"]] == ["a.bed", "b.bed"]
