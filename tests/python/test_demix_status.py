import json

import demix_status
import pytest

BARCODE = {"barcode_version": "10_06_2026-05-25", "barcode_version_iso": "2026-10-06T05-25"}


def test_ok_record():
    rec = demix_status.build_record("S1", "ok", barcode=BARCODE)
    assert rec["sample"] == "S1"
    assert rec["status"] == "ok"
    assert rec["reason"] == ""


def test_non_ok_requires_reason():
    with pytest.raises(ValueError):
        demix_status.build_record("S1", "failed", "", barcode=BARCODE)
    with pytest.raises(ValueError):
        demix_status.build_record("S1", "skipped_low_coverage", "   ", barcode=BARCODE)


def test_unknown_status_rejected():
    with pytest.raises(ValueError):
        demix_status.build_record("S1", "weird", "x", barcode=BARCODE)


def test_reason_is_single_line_and_truncated():
    rec = demix_status.build_record("S1", "failed", "a\nb\tc " + "x" * 1000, barcode=BARCODE)
    assert "\n" not in rec["reason"]
    assert len(rec["reason"]) == demix_status.MAX_REASON_CHARS
    assert rec["reason"].endswith("...")


def test_stats_are_copied():
    stats = {"pct_genome_ge_covcut": 12.5, "covcut": 10, "min_genome_coverage": 60}
    rec = demix_status.build_record("S1", "skipped_low_coverage", "low", stats, BARCODE)
    assert rec["genome_coverage_pct"] == 12.5
    assert rec["min_genome_coverage"] == 60


def barcode_file(tmp_path):
    path = tmp_path / "barcode_version.json"
    path.write_text(json.dumps(BARCODE))
    return str(path)


@pytest.mark.parametrize(
    "status,reason", [("ok", ""), ("failed", "x"), ("skipped_low_coverage", "y")]
)
def test_every_status_carries_the_barcode_version(tmp_path, status, reason):
    out = tmp_path / "S1.demix_status.json"
    args = [
        "--sample",
        "S1",
        "--status",
        status,
        "--reason",
        reason,
        "--barcode-info",
        barcode_file(tmp_path),
        "--output",
        str(out),
    ]
    assert demix_status.main(args) == 0
    rec = json.loads(out.read_text())
    assert rec["barcode_version"] == "10_06_2026-05-25"
    assert rec["barcode_version_iso"] == "2026-10-06T05-25"


def test_missing_barcode_version_refuses_to_write(tmp_path):
    with pytest.raises(ValueError, match="barcode_version"):
        demix_status.build_record("S1", "ok")
    with pytest.raises(ValueError, match="barcode_version"):
        demix_status.build_record("S1", "ok", barcode={"barcode_version": ""})
    empty = tmp_path / "empty.json"
    empty.write_text("{}")
    out = tmp_path / "x.json"
    args = ["--sample", "S1", "--status", "ok", "--barcode-info", str(empty), "--output", str(out)]
    assert demix_status.main(args) == 2
    assert not out.exists()


def test_cli_with_log_tail(tmp_path):
    log = tmp_path / "demix.log"
    log.write_text("line one\n\nTraceback\nValueError: boom\n")
    out = tmp_path / "S1.demix_status.json"
    rc = demix_status.main(
        [
            "--sample",
            "S1",
            "--status",
            "failed",
            "--reason",
            "freyja demix exited 1",
            "--reason-file",
            str(log),
            "--barcode-info",
            barcode_file(tmp_path),
            "--output",
            str(out),
        ]
    )
    assert rc == 0
    rec = json.loads(out.read_text())
    assert rec["status"] == "failed"
    assert rec["reason"].startswith("freyja demix exited 1: ")
    assert "ValueError: boom" in rec["reason"]


def test_cli_missing_reason_fails(tmp_path):
    out = tmp_path / "x.json"
    rc = demix_status.main(
        [
            "--sample",
            "S1",
            "--status",
            "failed",
            "--barcode-info",
            barcode_file(tmp_path),
            "--output",
            str(out),
        ]
    )
    assert rc == 2
    assert not out.exists()
