import json

import consensus_stats as cs


def fasta(tmp_path, text, name="c.fa"):
    path = tmp_path / name
    path.write_text(text)
    return str(path)


def test_ok_counts(tmp_path):
    f = fasta(tmp_path, ">S1\nACGTNNRY--\n")
    assert cs.compute_row("S1", f) == ["S1", "10", "4", "40.00", "2", "ok", "NA"]


def test_multiline_and_lowercase(tmp_path):
    f = fasta(tmp_path, ">S1 desc\nACGT\nnnAC\n\nGT\n")
    assert cs.compute_row("S1", f) == ["S1", "10", "2", "20.00", "0", "ok", "NA"]


def test_all_n_is_empty(tmp_path):
    f = fasta(tmp_path, ">S1\nNNNN\nNN-N\n")
    assert cs.compute_row("S1", f) == [
        "S1",
        "8",
        "8",
        "100.00",
        "0",
        "empty",
        "consensus contains only N or - characters",
    ]


def test_zero_length_is_failed(tmp_path):
    f = fasta(tmp_path, ">S1\n")
    row = cs.compute_row("S1", f)
    assert row[:6] == ["S1", "NA", "NA", "NA", "NA", "failed"]
    assert row[6] == "zero-length consensus sequence"


def test_missing_file_is_failed(tmp_path):
    assert cs.compute_row("S1", str(tmp_path / "nope.fa"))[5] == "failed"
    assert cs.compute_row("S1", None)[5:] == ["failed", "no consensus file"]


def test_not_fasta_is_failed(tmp_path):
    assert cs.compute_row("S1", fasta(tmp_path, "ACGT\n"))[5] == "failed"
    assert cs.compute_row("S1", fasta(tmp_path, ""))[5] == "failed"


def test_multiple_records_is_failed(tmp_path):
    assert cs.compute_row("S1", fasta(tmp_path, ">a\nACGT\n>b\nACGT\n"))[5] == "failed"


def test_reference_length_is_the_n_pct_denominator(tmp_path):
    f = fasta(tmp_path, ">S1\nNNAA\n")
    assert cs.compute_row("S1", f, reference_length=8)[3] == "25.00"


def test_main_and_concat_sorted(tmp_path):
    outs = []
    for sample, seq in [("b_2", "ACGT"), ("a-1", "NNNN"), ("c", "ACGN")]:
        f = fasta(tmp_path, f">{sample}\n{seq}\n", f"{sample}.fa")
        out = tmp_path / f"{sample}.tsv"
        assert cs.main(["--fasta", f, "--sample", sample, "-o", str(out)]) == 0
        outs.append(str(out))
    merged = tmp_path / "all.tsv"
    assert cs.main(["--concat", *outs, "-o", str(merged)]) == 0
    lines = merged.read_text().splitlines()
    assert lines[0].split("\t") == cs.COLUMNS
    assert [line.split("\t")[0] for line in lines[1:]] == ["a-1", "b_2", "c"]
    assert lines[1].split("\t")[5] == "empty"


def test_concat_rejects_bad_header_and_duplicates(tmp_path):
    bad = tmp_path / "bad.tsv"
    bad.write_text("sample\tlength\nS\t1\n")
    assert cs.main(["--concat", str(bad), "-o", str(tmp_path / "o.tsv")]) == 2
    good = tmp_path / "good.tsv"
    good.write_text("\t".join(cs.COLUMNS) + "\nS\t4\t0\t0.00\t0\tok\tNA\n")
    assert cs.main(["--concat", str(good), str(good), "-o", str(tmp_path / "o.tsv")]) == 2


def test_concat_of_nothing_writes_header(tmp_path):
    out = tmp_path / "o.tsv"
    assert cs.main(["--concat", "-o", str(out)]) == 0
    assert out.read_text() == "\t".join(cs.COLUMNS) + "\n"


def status_file(tmp_path, status, reason="", warnings=None, name="s.json"):
    path = tmp_path / name
    path.write_text(
        json.dumps({"sample": "S1", "status": status, "reason": reason, "warnings": warnings or []})
    )
    return str(path)


def test_json_failed_wins_over_all_n_sequence(tmp_path):
    # a failed sample is stored as an all-N genome, which alone reads as "empty"
    f = fasta(tmp_path, ">S1\nNNNN\nNNNN\n")
    sj = status_file(tmp_path, "failed", "bcftools consensus failed: REF mismatch at X:11287")
    row = cs.compute_row("S1", f, 8, sj)
    assert row == [
        "S1",
        "8",
        "8",
        "100.00",
        "0",
        "failed",
        "bcftools consensus failed: REF mismatch at X:11287",
    ]


def test_json_failed_without_reason_still_has_one(tmp_path):
    f = fasta(tmp_path, ">S1\nNNNN\n")
    row = cs.compute_row("S1", f, status_json=status_file(tmp_path, "failed"))
    assert row[5:] == ["failed", "consensus failed"]


def test_json_failed_with_unreadable_fasta_keeps_json_reason(tmp_path):
    row = cs.compute_row(
        "S1", str(tmp_path / "nope.fa"), status_json=status_file(tmp_path, "failed", "boom")
    )
    assert row == ["S1", "NA", "NA", "NA", "NA", "failed", "boom"]


def test_json_empty_reason_is_reported(tmp_path):
    f = fasta(tmp_path, ">S1\nNNNN\n")
    row = cs.compute_row("S1", f, status_json=status_file(tmp_path, "empty", "no usable coverage"))
    assert row[5:] == ["empty", "no usable coverage"]


def test_json_ok_with_warnings_reports_text_and_computed_status(tmp_path):
    f = fasta(tmp_path, ">S1\nACGT\nNNAC\n")
    reason = "1 variant(s) not applied: skipped variant at NC_045512.2:11287 (REF mismatch)"
    row = cs.compute_row("S1", f, status_json=status_file(tmp_path, "ok", reason, ["x"]))
    assert row == ["S1", "8", "2", "25.00", "0", "ok", reason]


def test_json_ok_but_all_n_sequence_is_empty(tmp_path):
    # "else compute as now": an ok JSON cannot turn an all-N sequence into ok
    f = fasta(tmp_path, ">S1\nNNNN\n")
    row = cs.compute_row("S1", f, status_json=status_file(tmp_path, "ok", ""))
    assert row[5] == "empty"


def test_json_ok_empty_reason_is_na(tmp_path):
    f = fasta(tmp_path, ">S1\nACGT\n")
    assert cs.compute_row("S1", f, status_json=status_file(tmp_path, "ok", ""))[5:] == ["ok", "NA"]


def test_unusable_json_is_ignored(tmp_path):
    f = fasta(tmp_path, ">S1\nNNNN\n")
    bad = tmp_path / "bad.json"
    bad.write_text("{not json")
    assert cs.compute_row("S1", f, status_json=str(bad))[5] == "empty"
    assert cs.compute_row("S1", f, status_json=str(tmp_path / "missing.json"))[5] == "empty"
    weird = status_file(tmp_path, "weird", "x", name="weird.json")
    assert cs.compute_row("S1", f, status_json=weird)[5] == "empty"


def test_reason_is_one_line_without_tabs(tmp_path):
    f = fasta(tmp_path, ">S1\nNNNN\n")
    sj = status_file(tmp_path, "failed", "line1\n\tline2")
    assert cs.compute_row("S1", f, status_json=sj)[6] == "line1 line2"


def test_main_status_json_and_concat_reason_column(tmp_path):
    f = fasta(tmp_path, ">S1\nNNNN\n")
    out = tmp_path / "s.tsv"
    sj = status_file(tmp_path, "failed", "bcftools consensus failed: x")
    assert cs.main(["--fasta", f, "--sample", "S1", "--status-json", sj, "-o", str(out)]) == 0
    lines = out.read_text().splitlines()
    assert lines[0].split("\t") == cs.COLUMNS and cs.COLUMNS[-1] == "reason"
    assert lines[1].split("\t")[5:] == ["failed", "bcftools consensus failed: x"]
    merged = tmp_path / "all.tsv"
    assert cs.main(["--concat", str(out), "-o", str(merged)]) == 0
    assert merged.read_text() == out.read_text()


def test_concat_accepts_legacy_six_column_files(tmp_path):
    legacy = tmp_path / "legacy.tsv"
    legacy.write_text("\t".join(cs.LEGACY_COLUMNS) + "\nS\t4\t0\t0.00\t0\tok\n")
    merged = tmp_path / "all.tsv"
    assert cs.main(["--concat", str(legacy), "-o", str(merged)]) == 0
    assert merged.read_text().splitlines()[1].split("\t") == [
        "S",
        "4",
        "0",
        "0.00",
        "0",
        "ok",
        "NA",
    ]


def test_concat_rejects_ragged_rows(tmp_path):
    ragged = tmp_path / "ragged.tsv"
    ragged.write_text("\t".join(cs.COLUMNS) + "\nS\t4\t0\t0.00\t0\tok\n")
    assert cs.main(["--concat", str(ragged), "-o", str(tmp_path / "o.tsv")]) == 2
