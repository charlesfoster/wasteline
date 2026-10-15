import consensus_stats as cs


def fasta(tmp_path, text, name="c.fa"):
    path = tmp_path / name
    path.write_text(text)
    return str(path)


def test_ok_counts(tmp_path):
    f = fasta(tmp_path, ">S1\nACGTNNRY--\n")
    assert cs.compute_row("S1", f) == ["S1", "10", "4", "40.00", "2", "ok"]


def test_multiline_and_lowercase(tmp_path):
    f = fasta(tmp_path, ">S1 desc\nACGT\nnnAC\n\nGT\n")
    assert cs.compute_row("S1", f) == ["S1", "10", "2", "20.00", "0", "ok"]


def test_all_n_is_empty(tmp_path):
    f = fasta(tmp_path, ">S1\nNNNN\nNN-N\n")
    assert cs.compute_row("S1", f) == ["S1", "8", "8", "100.00", "0", "empty"]


def test_zero_length_is_failed(tmp_path):
    f = fasta(tmp_path, ">S1\n")
    assert cs.compute_row("S1", f) == ["S1", "NA", "NA", "NA", "NA", "failed"]


def test_missing_file_is_failed(tmp_path):
    assert cs.compute_row("S1", str(tmp_path / "nope.fa"))[-1] == "failed"
    assert cs.compute_row("S1", None)[-1] == "failed"


def test_not_fasta_is_failed(tmp_path):
    assert cs.compute_row("S1", fasta(tmp_path, "ACGT\n"))[-1] == "failed"
    assert cs.compute_row("S1", fasta(tmp_path, ""))[-1] == "failed"


def test_multiple_records_is_failed(tmp_path):
    assert cs.compute_row("S1", fasta(tmp_path, ">a\nACGT\n>b\nACGT\n"))[-1] == "failed"


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
    assert lines[1].split("\t")[-1] == "empty"


def test_concat_rejects_bad_header_and_duplicates(tmp_path):
    bad = tmp_path / "bad.tsv"
    bad.write_text("sample\tlength\nS\t1\n")
    assert cs.main(["--concat", str(bad), "-o", str(tmp_path / "o.tsv")]) == 2
    good = tmp_path / "good.tsv"
    good.write_text("\t".join(cs.COLUMNS) + "\nS\t4\t0\t0.00\t0\tok\n")
    assert cs.main(["--concat", str(good), str(good), "-o", str(tmp_path / "o.tsv")]) == 2


def test_concat_of_nothing_writes_header(tmp_path):
    out = tmp_path / "o.tsv"
    assert cs.main(["--concat", "-o", str(out)]) == 0
    assert out.read_text() == "\t".join(cs.COLUMNS) + "\n"
