import gzip
import shutil
from pathlib import Path

import check_samplesheet as cs

FIX = Path(__file__).resolve().parents[1] / "fixtures" / "samplesheet"
GOOD_FQ = FIX / "fastq_custom" / "alpha_1_R1.fq.gz"


def sheet(tmp_path, text):
    p = tmp_path / "sheet.csv"
    p.write_text(text)
    return p


def errs(result):
    return " | ".join(result["errors"])


FULL_HDR = "sample,fastq_1,fastq_2,sample_type,collection_date,site\n"


def test_valid_sheets_pass_in_every_mode(tmp_path):
    p = sheet(
        tmp_path,
        FULL_HDR + "A_1,a_R1.fastq.gz,a_R2.fastq.gz,sample,2024-01-31,X\n"
        "NEG_1,n_R1.fq.gz,n_R2.fq.gz,negative,,\n",
    )
    r = cs.validate(p, "full")
    assert r["ok"] and r["n_samples"] == 2 and r["errors"] == []
    demix = sheet(
        tmp_path, "sample,variants,depths,fastp_json,sample_type\nA,v.tsv,d.depths,,sample\n"
    )
    assert cs.validate(demix, "demix")["ok"]
    agg = sheet(tmp_path, "sample,demix,collection_date\nA,a.demix,2024-05-05\n")
    assert cs.validate(agg, "aggregate")["ok"]


def test_invalid_sheet_reports_every_problem(tmp_path):
    p = sheet(
        tmp_path,
        FULL_HDR + "bad id,a.fastq.gz,b.fastq.gz,sample,,\n"
        "A,a.fastq.gz,b.fastq.gz,control,,\n"
        "B,a.fastq.gz,b.fastq.gz,sample,2024-02-30,\n"
        "C,a.fastq.gz,,sample,,\n"
        "C,a.fastq.gz,b.fastq.gz,sample,,\n",
    )
    r = cs.validate(p, "full")
    text = errs(r)
    assert not r["ok"]
    for needle in ("sample_type", "paired-end", "duplicate sample ID"):
        assert needle in text
    assert len(r["errors"]) >= 5


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
