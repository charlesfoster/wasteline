import gzip
import json
import shutil
from pathlib import Path

import check_samplesheet as cs
import pytest

FIX = Path(__file__).resolve().parents[1] / "fixtures" / "samplesheet"
GOOD_FQ = FIX / "fastq_custom" / "alpha_1_R1.fq.gz"


def sheet(tmp_path, text):
    p = tmp_path / "sheet.csv"
    p.write_text(text)
    return p


def errs(result):
    return " | ".join(result["errors"])


FULL_HDR = "sample,fastq_1,fastq_2,sample_type,collection_date,site\n"


def test_valid_full(tmp_path):
    p = sheet(
        tmp_path,
        FULL_HDR + "A_1,a_R1.fastq.gz,a_R2.fastq.gz,sample,2024-01-31,X\n"
        "NEG_1,n_R1.fq.gz,n_R2.fq.gz,negative,,\n",
    )
    r = cs.validate(p, "full")
    assert r["ok"] and r["n_samples"] == 2 and r["errors"] == []


def test_minimal_full_columns(tmp_path):
    p = sheet(tmp_path, "sample,fastq_1,fastq_2\nA,a_1.fastq.gz,a_2.fastq.gz\n")
    assert cs.validate(p, "full")["ok"]


def test_missing_required_column(tmp_path):
    p = sheet(tmp_path, "sample,fastq_1\nA,a.fastq.gz\n")
    r = cs.validate(p, "full")
    assert not r["ok"] and "fastq_2" in errs(r)


def test_unknown_column_warns(tmp_path):
    p = sheet(tmp_path, "sample,fastq_1,fastq_2,colour\nA,a.fastq.gz,b.fastq.gz,red\n")
    r = cs.validate(p, "full")
    assert r["ok"]
    assert any("colour" in w for w in r["warnings"])


def test_bad_ids(tmp_path):
    p = sheet(
        tmp_path,
        "sample,fastq_1,fastq_2\n"
        "bad id,a.fastq.gz,b.fastq.gz\nweird/x,a.fastq.gz,b.fastq.gz\nüber,a.fastq.gz,b.fastq.gz\n",
    )
    r = cs.validate(p, "full")
    assert not r["ok"] and len(r["errors"]) == 3


def test_underscore_ids_are_valid(tmp_path):
    p = sheet(tmp_path, "sample,fastq_1,fastq_2\nSite_1_week_2,a.fastq.gz,b.fastq.gz\n")
    assert cs.validate(p, "full")["ok"]


def test_duplicate_ids(tmp_path):
    p = sheet(
        tmp_path, "sample,fastq_1,fastq_2\nA,a.fastq.gz,b.fastq.gz\nA,c.fastq.gz,d.fastq.gz\n"
    )
    r = cs.validate(p, "full")
    assert not r["ok"] and "duplicate sample ID" in errs(r)


def test_sample_type_and_dates(tmp_path):
    p = sheet(
        tmp_path,
        FULL_HDR + "A,a.fastq.gz,b.fastq.gz,control,2024-01-01,\n"
        "B,a.fastq.gz,b.fastq.gz,sample,2024-13-01,\n"
        "C,a.fastq.gz,b.fastq.gz,sample,20240101,\n"
        "D,a.fastq.gz,b.fastq.gz,sample,2024-02-30,\n",
    )
    r = cs.validate(p, "full")
    assert len(r["errors"]) == 4
    assert "sample_type" in r["errors"][0]


def test_empty_fastq2_is_error(tmp_path):
    p = sheet(tmp_path, "sample,fastq_1,fastq_2\nA,a.fastq.gz,\n")
    r = cs.validate(p, "full")
    assert not r["ok"] and "paired-end" in errs(r)


def test_same_file_and_bad_extension_and_whitespace(tmp_path):
    p = sheet(
        tmp_path,
        "sample,fastq_1,fastq_2\nA,a.fastq.gz,a.fastq.gz\n"
        "B,a.fastq,b.fastq\nC,my reads.fastq.gz,c.fastq.gz\n",
    )
    r = cs.validate(p, "full")
    text = errs(r)
    assert "same file" in text and "must end with" in text and "whitespace" in text


def test_empty_and_header_only(tmp_path):
    assert not cs.validate(sheet(tmp_path, ""), "full")["ok"]
    r = cs.validate(sheet(tmp_path, "sample,fastq_1,fastq_2\n"), "full")
    assert not r["ok"] and "no samples" in errs(r)


def test_missing_sheet_file(tmp_path):
    r = cs.validate(tmp_path / "nope.csv", "full")
    assert not r["ok"] and "cannot read" in errs(r)


def test_bom_and_blank_lines_tolerated(tmp_path):
    p = tmp_path / "s.csv"
    p.write_bytes(b"\xef\xbb\xbfsample,fastq_1,fastq_2\n\nA,a.fastq.gz,b.fastq.gz\n\n")
    assert cs.validate(p, "full")["ok"]


def test_too_many_fields(tmp_path):
    p = sheet(tmp_path, "sample,fastq_1,fastq_2\nA,a.fastq.gz,b.fastq.gz,extra\n")
    assert not cs.validate(p, "full")["ok"]


def test_duplicate_header(tmp_path):
    p = sheet(tmp_path, "sample,sample,fastq_1,fastq_2\nA,A,a.fastq.gz,b.fastq.gz\n")
    r = cs.validate(p, "full")
    assert not r["ok"] and "duplicate column" in errs(r)


def test_demix_mode(tmp_path):
    ok = sheet(
        tmp_path, "sample,variants,depths,fastp_json,sample_type\nA,v.tsv,d.depths,,sample\n"
    )
    assert cs.validate(ok, "demix")["ok"]
    bad = sheet(tmp_path, "sample,variants,depths\nA,v.tsv,\n")
    r = cs.validate(bad, "demix")
    assert not r["ok"] and "empty depths" in errs(r)
    bad2 = sheet(tmp_path, "sample,variants\nA,v.tsv\n")
    assert not cs.validate(bad2, "demix")["ok"]


def test_aggregate_mode(tmp_path):
    ok = sheet(tmp_path, "sample,demix,collection_date\nA,a.demix,2024-05-05\n")
    assert cs.validate(ok, "aggregate")["ok"]
    r = cs.validate(sheet(tmp_path, "sample,demix\nA,\n"), "aggregate")
    assert not r["ok"]
    # a full sheet is not an aggregate sheet
    r = cs.validate(
        sheet(tmp_path, "sample,fastq_1,fastq_2\nA,a.fastq.gz,b.fastq.gz\n"), "aggregate"
    )
    assert not r["ok"] and "demix" in errs(r)


def test_check_files_fastq(tmp_path):
    shutil.copy(GOOD_FQ, tmp_path / "ok_R1.fastq.gz")
    shutil.copy(GOOD_FQ, tmp_path / "ok_R2.fastq.gz")
    with gzip.open(tmp_path / "bad_R1.fastq.gz", "wt") as fh:
        fh.write(">not fastq\nACGT\n")
    (tmp_path / "notgz_R1.fastq.gz").write_text("@r\nA\n+\nI\n")  # not actually gzip
    (tmp_path / "empty_R1.fastq.gz").write_bytes(gzip.compress(b""))
    p = sheet(
        tmp_path,
        "sample,fastq_1,fastq_2\n"
        "good,ok_R1.fastq.gz,ok_R2.fastq.gz\n"
        "badhdr,bad_R1.fastq.gz,ok_R2.fastq.gz\n"
        "notgz,notgz_R1.fastq.gz,ok_R2.fastq.gz\n"
        "empty,empty_R1.fastq.gz,ok_R2.fastq.gz\n"
        "gone,missing_R1.fastq.gz,ok_R2.fastq.gz\n",
    )
    assert cs.validate(p, "full")["ok"]  # no filesystem access without --check-files
    r = cs.validate(p, "full", check_files=True)
    assert not r["ok"]
    text = errs(r)
    assert "(badhdr)" in text and "'@'" in text
    assert "(notgz)" in text and "cannot read" in text
    assert "(empty)" in text and "empty" in text
    assert "(gone)" in text and "does not exist" in text
    assert "(good)" not in text


def test_check_files_relative_to_sheet_dir(tmp_path, monkeypatch):
    sub = tmp_path / "sub"
    sub.mkdir()
    shutil.copy(GOOD_FQ, sub / "a_R1.fq.gz")
    shutil.copy(GOOD_FQ, sub / "a_R2.fq.gz")
    p = sub / "s.csv"
    p.write_text("sample,fastq_1,fastq_2\nA,a_R1.fq.gz,a_R2.fq.gz\n")
    monkeypatch.chdir(tmp_path)
    assert cs.validate(p, "full", check_files=True)["ok"]


def test_check_files_demix_and_aggregate(tmp_path):
    (tmp_path / "v.tsv").write_text("x\n")
    (tmp_path / "d.depths").write_text("")
    p = sheet(tmp_path, "sample,variants,depths\nA,v.tsv,d.depths\nB,v.tsv,nope.depths\n")
    r = cs.validate(p, "demix", check_files=True)
    text = errs(r)
    assert "(A)" in text and "empty" in text and "(B)" in text and "does not exist" in text


def test_check_files_on_real_make_samplesheet_output(tmp_path, capsys):
    import make_samplesheet as ms

    out = tmp_path / "s.csv"
    assert ms.main([str(FIX / "fastq_illumina"), "-o", str(out)]) == 0
    capsys.readouterr()
    assert cs.validate(out, "full", check_files=True)["ok"]


def test_cli_json_stdout_and_exit_codes(tmp_path, capsys):
    good = sheet(tmp_path, "sample,fastq_1,fastq_2\nA,a.fastq.gz,b.fastq.gz\n")
    rc = cs.main([str(good), "--mode", "full"])
    cap = capsys.readouterr()
    assert rc == 0
    data = json.loads(cap.out)
    assert data == {"ok": True, "mode": "full", "n_samples": 1, "errors": [], "warnings": []}
    assert "PASS" in cap.err

    bad = sheet(tmp_path, "sample,fastq_1,fastq_2\nA,a.fastq.gz,\n")
    rc = cs.main([str(bad), "--mode", "full"])
    cap = capsys.readouterr()
    assert rc == 1
    assert json.loads(cap.out)["ok"] is False
    assert "FAIL" in cap.err and "ERROR" in cap.err


def test_cli_requires_mode(tmp_path):
    with pytest.raises(SystemExit) as exc:
        cs.main([str(tmp_path / "x.csv")])
    assert exc.value.code == 2
