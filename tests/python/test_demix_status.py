import json

import demix_status
import pytest

BARCODE = {"barcode_version": "10_06_2026-05-25", "barcode_version_iso": "2026-10-06T05-25"}


def barcode_file(tmp_path):
    path = tmp_path / "barcode_version.json"
    path.write_text(json.dumps(BARCODE))
    return str(path)


def test_non_ok_status_requires_a_reason_and_barcode_version():
    assert demix_status.build_record("S1", "ok", barcode=BARCODE)["reason"] == ""
    with pytest.raises(ValueError):
        demix_status.build_record("S1", "failed", "", barcode=BARCODE)
    with pytest.raises(ValueError):
        demix_status.build_record("S1", "weird", "x", barcode=BARCODE)
    with pytest.raises(ValueError, match="barcode_version"):
        demix_status.build_record("S1", "ok")


def test_cli_failed_record_with_log_tail_and_barcode_version(tmp_path):
    log = tmp_path / "demix.log"
    log.write_text("line one\n\nTraceback\nValueError: boom\n")
    out = tmp_path / "S1.demix_status.json"
    args = ["--sample", "S1", "--status", "failed", "--barcode-info", barcode_file(tmp_path)]
    assert demix_status.main([*args, "--output", str(tmp_path / "x.json")]) == 2  # no reason
    assert not (tmp_path / "x.json").exists()
    rc = demix_status.main(
        [
            *args,
            "--reason",
            "freyja demix exited 1",
            "--reason-file",
            str(log),
            "--output",
            str(out),
        ]
    )
    assert rc == 0
    rec = json.loads(out.read_text())
    assert rec["status"] == "failed"
    assert rec["reason"].startswith("freyja demix exited 1: ")
    assert "ValueError: boom" in rec["reason"]
    assert rec["barcode_version"] == "10_06_2026-05-25"
