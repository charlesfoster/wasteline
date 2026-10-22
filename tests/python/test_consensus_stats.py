import json

import consensus_stats as cs


def fasta(tmp_path, text, name="c.fa"):
    path = tmp_path / name
    path.write_text(text)
    return str(path)


def status_file(tmp_path, status, reason=""):
    path = tmp_path / "s.json"
    path.write_text(
        json.dumps({"sample": "S1", "status": status, "reason": reason, "warnings": []})
    )
    return str(path)


def test_every_sample_gets_a_row_with_status_and_reason(tmp_path):
    """ok / empty / failed samples all yield a row with a status and reason (never an empty file)."""
    assert cs.compute_row("S1", fasta(tmp_path, ">S1 d\nACGT\nnnAC\n\nGT\n")) == [
        "S1", "10", "2", "80.00", "0", "ok", "NA",
    ]  # fmt: skip
    all_n = fasta(tmp_path, ">S1\nNNNN\nNN-N\n")
    assert cs.compute_row("S1", all_n)[5:] == ["empty", "consensus contains only N or - characters"]
    assert cs.compute_row("S1", None)[5:] == ["failed", "no consensus file"]
    assert cs.compute_row("S1", fasta(tmp_path, ">S1\n"))[5] == "failed"
    # a failed sample is stored as an all-N genome: the status JSON must win over "empty"
    sj = status_file(tmp_path, "failed", "bcftools consensus failed: REF mismatch at X:11287")
    row = cs.compute_row("S1", fasta(tmp_path, ">S1\nNNNN\nNNNN\n"), 8, sj)
    assert row[5:] == ["failed", "bcftools consensus failed: REF mismatch at X:11287"]


def test_main_and_concat_sorted_with_reason_column(tmp_path):
    outs = []
    for sample, seq in [("b_2", "ACGT"), ("a-1", "NNNN"), ("c", "ACGN")]:
        f = fasta(tmp_path, f">{sample}\n{seq}\n", f"{sample}.fa")
        out = tmp_path / f"{sample}.tsv"
        args = ["--fasta", f, "--sample", sample, "-o", str(out)]
        if sample == "a-1":
            args += ["--status-json", status_file(tmp_path, "failed", "boom")]
        assert cs.main(args) == 0
        outs.append(str(out))
    merged = tmp_path / "all.tsv"
    assert cs.main(["--concat", *outs, "-o", str(merged)]) == 0
    lines = merged.read_text().splitlines()
    assert lines[0].split("\t") == cs.COLUMNS
    assert [line.split("\t")[0] for line in lines[1:]] == ["a-1", "b_2", "c"]
    assert lines[1].split("\t")[5:] == ["failed", "boom"]
