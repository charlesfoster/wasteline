import clean_aggregate_labels as cal
import pytest

HEADER = "\tsummarized\tlineages\tabundances\tresid\tcoverage\tpathogen"


def test_clean_label_longest_suffix_wins():
    suffixes = [".tsv", ".variants.tsv"]
    assert cal.clean_label("S1.variants.tsv", suffixes) == "S1"
    assert cal.clean_label("S1.tsv", suffixes) == "S1"


def test_underscores_and_dots_are_preserved():
    assert cal.clean_label("NEG_Site.A_1.2.variants.tsv", [".variants.tsv"]) == "NEG_Site.A_1.2"


def test_label_without_suffix_untouched_and_never_emptied():
    assert cal.clean_label("S1", [".variants.tsv"]) == "S1"
    assert cal.clean_label(".variants.tsv", [".variants.tsv"]) == ".variants.tsv"


def test_only_label_column_changes():
    rows = [
        HEADER,
        "S_1.variants.tsv\t[('Omicron', 0.9)]\tBA.2 XBB\t0.9 0.1\t1.5\t99.2\tSARS-CoV-2",
    ]
    out = list(cal.clean_lines(rows, [".variants.tsv"]))
    assert out[0] == HEADER
    assert out[1] == "S_1\t[('Omicron', 0.9)]\tBA.2 XBB\t0.9 0.1\t1.5\t99.2\tSARS-CoV-2"


def test_duplicate_labels_rejected():
    rows = [HEADER, "A.variants.tsv\tx", "A.tsv\ty"]
    with pytest.raises(ValueError, match="duplicate"):
        list(cal.clean_lines(rows, [".variants.tsv", ".tsv"]))


def test_cli_roundtrip(tmp_path):
    src = tmp_path / "in.tsv"
    src.write_text(f"{HEADER}\nA_1.variants.tsv\ta\nB-2.variants.tsv\tb\n")
    dst = tmp_path / "out.tsv"
    assert cal.main([str(src), str(dst), "--suffix", ".variants.tsv"]) == 0
    assert dst.read_text() == f"{HEADER}\nA_1\ta\nB-2\tb\n"


def test_cli_duplicate_returns_error(tmp_path):
    src = tmp_path / "in.tsv"
    src.write_text(f"{HEADER}\nA.tsv\ta\nA.variants.tsv\tb\n")
    dst = tmp_path / "out.tsv"
    assert cal.main([str(src), str(dst), "--suffix", ".tsv", "--suffix", ".variants.tsv"]) == 1


def test_cli_sort_orders_rows_by_label_without_suffixes(tmp_path):
    src = tmp_path / "in.tsv"
    src.write_text(f"{HEADER}\nSite_1\tc\nMIX_B\tb\nMIX_A\ta\n")
    dst = tmp_path / "out.tsv"
    assert cal.main([str(src), str(dst), "--sort"]) == 0
    assert dst.read_text() == f"{HEADER}\nMIX_A\ta\nMIX_B\tb\nSite_1\tc\n"
