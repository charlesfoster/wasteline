"""ivar_trim_stats.py: parses the summary `ivar trim` (1.4.4) prints to stderr.

The fixture in tests/fixtures/qc/ivar_trim_PURE_C.ivar.log is real ivar 1.4.4 output (stderr).
"""

from pathlib import Path

import ivar_trim_stats as its

FIX = Path(__file__).resolve().parents[1] / "fixtures" / "qc"


def run(tmp_path, log_path=None, sample="S"):
    out = tmp_path / "stats.tsv"
    argv = ["--sample", sample, "-o", str(out)]
    if log_path is not None:
        argv += ["--log", str(log_path)]
    assert its.main(argv) == 0
    header, row = out.read_text().splitlines()
    assert header.split("\t") == its.COLUMNS
    return dict(zip(its.COLUMNS, row.split("\t")))


def test_real_log_parsed_and_missing_log_still_gives_a_row(tmp_path):
    row = run(tmp_path, FIX / "ivar_trim_PURE_C.ivar.log", sample="PURE_C")
    assert row == {
        "sample": "PURE_C",
        "reads_total": "15998",
        "primer_trimmed_reads": "3515",
        "primer_trimmed_pct": "21.97",
        "no_primer_reads": "12483",
        "no_primer_pct": "78.03",
        "quality_trimmed_reads": "0",
        "quality_trimmed_pct": "0.00",
        "status": "ok",
    }
    missing = run(tmp_path, tmp_path / "nope.log")
    assert missing["status"] == "missing"
    assert all(missing[c] == "NA" for c in its.COLUMNS[1:-1])
