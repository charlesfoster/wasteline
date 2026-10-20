import sys
from pathlib import Path

BIN = Path(__file__).resolve().parents[2] / "bin"
if str(BIN) not in sys.path:
    sys.path.insert(0, str(BIN))
import filter_aggregated

HEADER = "\tsummarized\tlineages\tabundances\tresid\tcoverage\tpathogen\n"
ROW = "{}\t[('Omicron', 0.9)]\tBA.2\t0.9\t1.5\t99.2\tSARS-CoV-2\n"


def run(tmp_path, labels, statuses):
    agg = tmp_path / "agg.tsv"
    agg.write_text(HEADER + "".join(ROW.format(label) for label in labels))
    status = tmp_path / "demix_status.tsv"
    status.write_text(
        "sample\tstatus\treason\n" + "".join(f"{s}\t{st}\tr\n" for s, st in statuses.items())
    )
    out = tmp_path / "out.tsv"
    args = [str(agg), "--demix-status", str(status), "--exclude-status", "low_confidence"]
    assert filter_aggregated.main([*args, "--output", str(out)]) == 0
    return out.read_text()


def test_drops_only_the_excluded_status_and_keeps_the_format(tmp_path):
    statuses = {"A_1": "ok", "B_2": "low_confidence", "C_3": "low_confidence"}
    # C_3 is not in the table; D_4 is not in the status file and is kept
    text = run(tmp_path, ["A_1", "B_2", "D_4"], statuses)
    assert text == HEADER + ROW.format("A_1") + ROW.format("D_4")
    # nothing left: the header only
    assert run(tmp_path, ["B_2"], {"B_2": "low_confidence"}) == HEADER
