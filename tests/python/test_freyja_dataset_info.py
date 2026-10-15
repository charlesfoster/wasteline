import hashlib
import json

import freyja_dataset_info as fdi
import pytest

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


@pytest.mark.parametrize(
    "text,expected",
    [
        ("10_06_2026-05-25", ISO),
        ("09_17_2026-00-46\n", "2026-09-17T00-46"),
        ("2026-10-06T05-25", ISO),
        ("2026-10-06T05:25", ISO),
        ("13_40_2026-05-25", None),
        ("10_06_2026-25-25", None),
        ("2026-10-06", None),
        ("garbage", None),
        ("", None),
        (None, None),
    ],
)
def test_stamp_to_iso(text, expected):
    assert fdi.stamp_to_iso(text) == expected


def test_iso_to_stamp_roundtrip():
    assert fdi.iso_to_stamp(ISO) == STAMP
    assert fdi.stamp_to_iso(fdi.iso_to_stamp(ISO)) == ISO


def test_is_bare_date():
    assert fdi.is_bare_date("2026-10-06")
    assert not fdi.is_bare_date(ISO)
    assert not fdi.is_bare_date(STAMP)


def test_parse_reported_version():
    assert fdi.parse_reported_version(REPORT) == "09_17_2026-00-46"
    assert fdi.parse_reported_version("Using user-supplied barcodes at path:\n/a/b.csv\n") is None
    assert fdi.parse_reported_version(None) is None


def test_describe_directory_has_every_mandatory_key(tmp_path):
    ds = make_dataset(tmp_path / "ds", sidecar={"downloaded_at": "2026-10-15T01:02:03Z"})
    info = fdi.describe(ds, "latest", "PATHOGEN", freyja_version="2.0.5", reported_text=REPORT)
    for key in fdi.FREYJA_DATA_KEYS:
        assert key in info, key
    assert info["barcode_version"] == STAMP  # verbatim
    assert info["barcode_version_iso"] == ISO
    assert info["barcodes_sha256"] == hashlib.sha256(b"barcode-bytes").hexdigest()
    assert info["barcodes_format"] == "feather"
    assert info["lineage_mutations_sha256"] == hashlib.sha256(b"{}").hexdigest()
    assert info["curated_lineages_sha256"] is None  # not in this dataset
    assert info["pathogen"] == "PATHOGEN"
    assert info["downloaded_at"] == "2026-10-15T01:02:03Z"
    assert info["dataset_path"] == str(ds)
    # freyja reports its bundled data, so it only counts as "reported" for source == container
    assert info["freyja_reported_barcode_version"] is None
    assert info["freyja_bundled_barcode_version"] == "09_17_2026-00-46"


def test_container_source_records_reported_version(tmp_path):
    ds = make_dataset(tmp_path / "ds", stamp="09_17_2026-00-46")
    info = fdi.describe(ds, "container", "P", reported_text=REPORT)
    assert info["freyja_reported_barcode_version"] == info["barcode_version"]


def test_csv_barcodes_in_directory_and_feather_preferred(tmp_path):
    ds = make_dataset(tmp_path / "a", barcode="usher_barcodes.csv")
    assert fdi.describe(ds, "pinned", "P")["barcodes_format"] == "csv"
    ds = make_dataset(tmp_path / "b")
    (ds / "usher_barcodes.csv").write_text("x")
    assert fdi.describe(ds, "pinned", "P")["barcodes_file"] == "usher_barcodes.feather"


def test_directory_without_barcodes_or_with_bad_stamp_rejected(tmp_path):
    (tmp_path / "empty").mkdir()
    with pytest.raises(ValueError, match="no usher_barcodes"):
        fdi.describe(tmp_path / "empty", "custom", "P")
    ds = make_dataset(tmp_path / "bad", stamp="not a stamp")
    with pytest.raises(ValueError, match="unrecognised"):
        fdi.describe(ds, "custom", "P")


def test_directory_without_stamp_file_is_custom(tmp_path):
    ds = make_dataset(tmp_path / "ds", stamp=None)
    info = fdi.describe(ds, "custom", "P")
    assert info["barcode_version"] == "custom"
    assert info["barcode_version_iso"] is None


def test_single_file_is_custom(tmp_path):
    f = tmp_path / "custom.csv"
    f.write_text("a,b\n")
    info = fdi.describe(f, "custom", "P")
    assert info["barcode_version"] == "custom"
    assert info["barcode_version_iso"] is None
    assert info["barcodes_sha256"] == hashlib.sha256(b"a,b\n").hexdigest()
    for key in fdi.FREYJA_DATA_KEYS:
        assert key in info, key


def test_bad_extension_missing_path_and_bad_source(tmp_path):
    f = tmp_path / "custom.txt"
    f.write_text("x")
    with pytest.raises(ValueError):
        fdi.describe(f, "custom", "P")
    with pytest.raises(ValueError, match="does not exist"):
        fdi.describe(tmp_path / "missing", "custom", "P")
    with pytest.raises(ValueError, match="source"):
        fdi.describe(f, "weird", "P")


def test_cli_print_iso(tmp_path, capsys):
    ds = make_dataset(tmp_path / "ds")
    assert fdi.main([str(ds), "--print-iso"]) == 0
    assert capsys.readouterr().out.strip() == ISO


def test_cli_print_iso_fails_without_stamp(tmp_path):
    ds = make_dataset(tmp_path / "ds", stamp=None)
    assert fdi.main([str(ds), "--print-iso"]) == 1


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


def test_dataset_paths_directory_and_file(tmp_path, capsys):
    ds = make_dataset(tmp_path / "ds")
    barcode, lineage_yml = fdi.dataset_paths(ds)
    assert barcode.name == "usher_barcodes.feather"
    assert lineage_yml.name == "lineages.yml"
    assert fdi.main([str(ds), "--print-barcodes-path"]) == 0
    assert capsys.readouterr().out.strip().endswith("usher_barcodes.feather")
    single = tmp_path / "x.csv"
    single.write_text("a")
    assert fdi.dataset_paths(single) == (single, None)
    assert fdi.main([str(single), "--print-lineageyml-path"]) == 0
    assert capsys.readouterr().out == "\n"


def test_dataset_paths_rejects_unusable(tmp_path):
    with pytest.raises(ValueError):
        fdi.dataset_paths(tmp_path / "missing")
    assert fdi.main([str(tmp_path / "missing"), "--print-barcodes-path"]) == 1
