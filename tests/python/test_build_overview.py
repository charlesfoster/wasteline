import re
from pathlib import Path

import build_overview as bo

ROOT = Path(__file__).resolve().parents[2]
FIX = ROOT / "tests" / "fixtures" / "overview"


def test_overview_lists_every_sample_target_cell_and_links(tmp_path):
    # RSVa has its tables and a report; RSVb has none of them (a target with no files must still show '-' cells)
    manifest = tmp_path / "manifest.tsv"
    lines = [
        f"RSVa\t{role}\t{FIX / 'RSVa' / name}"
        for role, name in (
            ("aggregated", "freyja_aggregated.tsv"),
            ("qc", "qc_summary.tsv"),
            ("demix_status", "demix_status.tsv"),
        )
    ]
    manifest.write_text(
        "\n".join(["RSVa\ttarget\t-", "RSVb\ttarget\t-", *lines, "RSVa\treport\t-"]) + "\n"
    )
    out = tmp_path / "overview.html"
    map_out = tmp_path / "map.csv"
    argv = [
        "--targets",
        str(FIX / "targets.tsv"),
        "--inputs",
        str(manifest),
        "--panel",
        "/home/someone/panels/rsv.csv",
        "--out",
        str(out),
    ]
    assert bo.main(argv) == 0
    html = out.read_text()
    # 3 samples x 2 targets, every cell present, in a matrix with one column per target
    assert len(re.findall(r'<td class="cell', html)) == 6
    assert 'href="RSVa/report.html"' in html and 'href="RSVb/report.html"' not in html
    assert (
        "RSVa-L2" in html and "75.0 %" in html and "QC WARN" in html
    )  # top lineage and QC of the analysed cell
    assert "no_target_detected" in html and "qc_fail" in html  # flagged samples with the reason
    assert (
        "mixed" in html and "/home/someone" not in html and "rsv.csv" in html
    )  # the panel is shown by file name only
    # full anonymisation: sample IDs are replaced and the map is written
    assert (
        bo.main([*argv[:-1], str(out), "--anonymize", "full", "--anonymize-map-out", str(map_out)])
        == 0
    )
    html = out.read_text()
    assert (
        "S001" in html
        and ">S1<" not in html
        and "NTC" not in re.sub(r"<em class=\"chip\">NTC</em>", "", html)
    )
    assert "S1" in map_out.read_text() and "S001" in map_out.read_text()
