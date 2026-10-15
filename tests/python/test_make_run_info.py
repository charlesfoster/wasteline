import hashlib
import json
from datetime import UTC, datetime

import make_run_info as mri
import pytest
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


def test_parse_versions_yaml():
    parsed = mri.parse_versions_yaml(VERSIONS)
    assert parsed["WASTELINE:WORKFLOW:FREYJA_VARIANTS"] == {"freyja": "2.0.5"}
    assert mri.tool_versions(parsed)["fastp"] == "0.24.0"


def test_build_run_info(tmp_path):
    files = {}
    for name, content in {"s.csv": b"sample\n", "p.bed": b"bed\n", "r.fasta": b">r\nA\n"}.items():
        (tmp_path / name).write_bytes(content)
        files[name] = str(tmp_path / name)
    base = {"pipeline_version": "0.1.0dev", "mode": "full", "params": {"a": 1}}
    dataset = dict(DATASET)
    now = datetime(2026, 10, 15, 1, 2, 3, tzinfo=UTC)
    info = mri.build_run_info(
        base, files["s.csv"], files["p.bed"], files["r.fasta"], dataset, VERSIONS, now=now
    )
    assert info["pipeline_version"] == "0.1.0dev"
    assert info["freyja_version"] == "2.0.5"
    for key in FREYJA_DATA_KEYS:
        assert key in info["freyja_data"], key
    assert info["freyja_data"]["barcode_version"] == "v-barcode_version"
    assert info["input_checksums"]["primer_bed"]["sha256"] == hashlib.sha256(b"bed\n").hexdigest()
    assert info["provenance_generated"] == "2026-10-15T01:02:03Z"
    assert "pipeline_version" in base and "tool_versions" not in base  # input not mutated


def test_cli_end_to_end(tmp_path):
    for name in ("s.csv", "p.bed", "r.fasta"):
        (tmp_path / name).write_text(name)
    (tmp_path / "base.json").write_text(json.dumps({"mode": "full"}))
    (tmp_path / "ds.json").write_text(json.dumps(DATASET))
    (tmp_path / "v.yml").write_text(VERSIONS)
    out = tmp_path / "run_info.json"
    rc = mri.main(
        [
            "--base",
            str(tmp_path / "base.json"),
            "--samplesheet",
            str(tmp_path / "s.csv"),
            "--bed",
            str(tmp_path / "p.bed"),
            "--fasta",
            str(tmp_path / "r.fasta"),
            "--dataset-info",
            str(tmp_path / "ds.json"),
            "--versions",
            str(tmp_path / "v.yml"),
            "--output",
            str(out),
        ]
    )
    assert rc == 0
    info = json.loads(out.read_text())
    assert info["mode"] == "full"
    assert info["freyja_version"] == "2.0.5"


def test_cli_missing_input_fails(tmp_path):
    rc = mri.main(
        [
            "--base",
            str(tmp_path / "nope.json"),
            "--samplesheet",
            "a",
            "--bed",
            "b",
            "--fasta",
            "c",
            "--dataset-info",
            "d",
            "--output",
            str(tmp_path / "o.json"),
        ]
    )
    assert rc == 1


def test_missing_versions_gives_null_freyja(tmp_path):
    for name in ("s", "p", "r"):
        (tmp_path / name).write_text(name)
    info = mri.build_run_info(
        {}, str(tmp_path / "s"), str(tmp_path / "p"), str(tmp_path / "r"), {}, ""
    )
    assert info["freyja_version"] is None


def test_modes_without_reads_need_no_bed_or_fasta_and_keep_previous_runs(tmp_path):
    sheet = tmp_path / "s.csv"
    sheet.write_text("sample,demix\n")
    dataset = DATASET
    previous = [{"path": "/old", "run_info_sha256": None, "barcode_version": "unknown"}]
    base = {"mode": "demix", "previous_runs": previous}
    info = mri.build_run_info(base, sheet, None, None, dataset, VERSIONS)
    assert info["input_checksums"].keys() == {"samplesheet"}
    assert info["previous_runs"] == previous
    # the contract key exists even when the base JSON has none
    assert mri.build_run_info({}, sheet, None, None, dataset, VERSIONS)["previous_runs"] == []


def test_second_dataset_is_recorded_as_freyja_data_compare(tmp_path):
    sheet = tmp_path / "s.csv"
    sheet.write_text("sample,demix\n")
    other = {**DATASET, "barcode_version": "08_01_2026-00-10"}
    info = mri.build_run_info({}, sheet, None, None, DATASET, VERSIONS, dataset_info_compare=other)
    assert info["freyja_data_compare"]["barcode_version"] == "08_01_2026-00-10"
    assert "freyja_data_compare" not in mri.build_run_info({}, sheet, None, None, DATASET, VERSIONS)
    broken = {k: v for k, v in other.items() if k != "barcodes_sha256"}
    with pytest.raises(ValueError, match="barcodes_sha256"):
        mri.build_run_info({}, sheet, None, None, DATASET, VERSIONS, dataset_info_compare=broken)
