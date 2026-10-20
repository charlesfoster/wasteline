import relabel_demix as rd

DEMIX = "\tS1.freyja_variants.tsv\nsummarized\t[('Omicron', 0.7)]\nlineages\tBA.2 AY.4\ncoverage\t99.3\n"


def test_label_replaced_rest_untouched_and_input_not_modified(tmp_path):
    src = tmp_path / "S1.demix"
    src.write_text(DEMIX)
    dst = tmp_path / "out.demix"
    assert rd.main([str(src), "--sample", "Site_1_week_2", "--output", str(dst)]) == 0
    assert src.read_text() == DEMIX
    out = dst.read_text()
    assert out.splitlines()[0] == "\tSite_1_week_2"
    assert out.split("\n", 1)[1] == DEMIX.split("\n", 1)[1]
