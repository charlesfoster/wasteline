import clean_aggregate_labels as cal

HEADER = "\tsummarized\tlineages\tabundances\tresid\tcoverage\tpathogen"


def test_only_label_column_changes_and_ids_are_not_rewritten():
    """Longest suffix wins; underscores and dots in sample IDs are preserved."""
    rows = [
        HEADER,
        "S_1.variants.tsv\t[('Omicron', 0.9)]\tBA.2 XBB\t0.9 0.1\t1.5\t99.2\tSARS-CoV-2",
        "NEG_Site.A_1.2.tsv\tx",
    ]
    out = list(cal.clean_lines(rows, [".tsv", ".variants.tsv"]))
    assert out[0] == HEADER
    assert out[1] == "S_1\t[('Omicron', 0.9)]\tBA.2 XBB\t0.9 0.1\t1.5\t99.2\tSARS-CoV-2"
    assert out[2] == "NEG_Site.A_1.2\tx"
