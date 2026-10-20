import hashlib
import json

import freyja_dataset_info as fdi

STAMP = "10_06_2026-05-25"
ISO = "2026-10-06T05-25"
REPORT = "SARS-CoV-2 Barcode version:\n09_17_2026-00-46\nLocal barcode path:\n/x/usher_barcodes.feather\n"


def make_dataset(path, stamp=STAMP, barcode="usher_barcodes.feather", sidecar=None):
    path.mkdir()
    (path / barcode).write_bytes(b"barcode-bytes")
    if stamp is not None:
        (path / "last_barcode_update.txt").write_text(stamp + "\n")
    (path / "lineages.yml").write_text("- name: A\n")
    (path / "lineage_mutations.json").write_text("{}")
    if sidecar:
        (path / fdi.SIDECAR).write_text(json.dumps(sidecar))
    return path


def test_describe_directory_and_single_file(tmp_path):
    assert fdi.stamp_to_iso(STAMP) == ISO and fdi.stamp_to_iso("garbage") is None
    ds = make_dataset(tmp_path / "ds", sidecar={"downloaded_at": "2026-10-15T01:02:03Z"})
    info = fdi.describe(ds, "latest", "PATHOGEN", freyja_version="2.0.5", reported_text=REPORT)
    for key in fdi.FREYJA_DATA_KEYS:
        assert key in info, key
    assert info["barcode_version"] == STAMP  # verbatim
    assert info["barcode_version_iso"] == ISO
    assert info["barcodes_sha256"] == hashlib.sha256(b"barcode-bytes").hexdigest()
    assert info["barcodes_format"] == "feather"
    assert info["downloaded_at"] == "2026-10-15T01:02:03Z"
    # a user-supplied barcode file has no version stamp: recorded as custom
    f = tmp_path / "custom.csv"
    f.write_text("a,b\n")
    custom = fdi.describe(f, "custom", "P")
    assert custom["barcode_version"] == "custom" and custom["barcode_version_iso"] is None
    assert custom["barcodes_sha256"] == hashlib.sha256(b"a,b\n").hexdigest()


def test_cli_writes_both_json_files(tmp_path):
    ds = make_dataset(tmp_path / "ds")
    rep = tmp_path / "rep.txt"
    rep.write_text(REPORT)
    out, ver = tmp_path / "freyja_data.json", tmp_path / "barcode_version.json"
    rc = fdi.main(
        [
            str(ds),
            "--source",
            "pinned",
            "--pathogen",
            "P",
            "--freyja-version",
            "2.0.5",
            "--reported-file",
            str(rep),
            "--output",
            str(out),
            "--version-output",
            str(ver),
        ]
    )
    assert rc == 0
    full = json.loads(out.read_text())
    small = json.loads(ver.read_text())
    assert set(small) == set(fdi.BARCODE_VERSION_KEYS)
    assert small["barcode_version"] == full["barcode_version"] == STAMP
    assert small["barcode_version_iso"] == ISO
    assert small["freyja_version"] == "2.0.5"
