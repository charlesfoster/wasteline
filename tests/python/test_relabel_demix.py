import pytest
import relabel_demix as rd

DEMIX = "\tS1.freyja_variants.tsv\nsummarized\t[('Omicron', 0.7)]\nlineages\tBA.2 AY.4\ncoverage\t99.3\n"


def test_label_replaced_rest_untouched():
    out = rd.relabel(DEMIX, "Site_1_week_2")
    assert out.splitlines()[0] == "\tSite_1_week_2"
    assert out.split("\n", 1)[1] == DEMIX.split("\n", 1)[1]


def test_path_label_replaced():
    assert rd.relabel("\t/data/old/S1.tsv\nx\ty\n", "S1").startswith("\tS1\n")


def test_not_a_demix_file_rejected():
    with pytest.raises(ValueError, match="not a freyja demix file"):
        rd.relabel("summarized\t1\n", "S1")


def test_cli_does_not_modify_input(tmp_path):
    src = tmp_path / "S1.demix"
    src.write_text(DEMIX)
    dst = tmp_path / "out.demix"
    assert rd.main([str(src), "--sample", "S_1", "--output", str(dst)]) == 0
    assert src.read_text() == DEMIX
    assert dst.read_text().startswith("\tS_1\n")
