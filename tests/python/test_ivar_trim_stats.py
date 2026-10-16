"""ivar_trim_stats.py: parses the summary `ivar trim` (1.4.4) prints to stderr.

The fixtures in tests/fixtures/qc/ivar_trim_*.ivar.log are real ivar 1.4.4 output (stderr) on the PURE_C test
sample with the midnight scheme; the `-m`/`-k` variants and LOWCOV_D are the same with the repeated BED
warnings removed.
"""

from pathlib import Path

import ivar_trim_stats as its
import pytest

FIX = Path(__file__).resolve().parents[1] / "fixtures" / "qc"


def run(tmp_path, log_text=None, log_path=None, sample="S"):
    out = tmp_path / "stats.tsv"
    argv = ["--sample", sample, "-o", str(out)]
    if log_text is not None:
        log_path = tmp_path / "x.ivar.log"
        log_path.write_text(log_text)
    if log_path is not None:
        argv += ["--log", str(log_path)]
    assert its.main(argv) == 0
    header, row = out.read_text().splitlines()
    assert header.split("\t") == its.COLUMNS
    return dict(zip(its.COLUMNS, row.split("\t")))


def test_real_log_with_e_flag(tmp_path):
    row = run(tmp_path, log_path=FIX / "ivar_trim_PURE_C.ivar.log", sample="PURE_C")
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


def test_real_log_small_sample(tmp_path):
    row = run(tmp_path, log_path=FIX / "ivar_trim_LOWCOV_D.ivar.log")
    assert (row["primer_trimmed_reads"], row["no_primer_reads"], row["reads_total"]) == (
        "27",
        "92",
        "119",
    )
    assert row["primer_trimmed_pct"] == "22.69"


def test_real_log_without_e_and_with_quality_trimming(tmp_path):
    """No -e: 'that started outside of primer regions were not written'; quality-trimmed reads overlap the others."""
    row = run(tmp_path, log_path=FIX / "ivar_trim_no_e_min130.ivar.log")
    assert row["primer_trimmed_reads"] == "3516"
    assert row["no_primer_reads"] == "12484"
    assert row["quality_trimmed_reads"] == "3226"
    assert row["reads_total"] == "16000"
    assert row["quality_trimmed_pct"] == "20.16"
    assert row["primer_trimmed_pct"] == "21.98"  # 100 * 3516 / 16000 = 21.975, as ivar prints
    assert row["status"] == "ok"


def test_real_log_with_k_flag(tmp_path):
    """-k wording: 'were marked as failed' and 'Since the -ek flags were given'."""
    row = run(tmp_path, log_path=FIX / "ivar_trim_keep.ivar.log")
    assert row["primer_trimmed_pct"] == "21.97"
    assert row["no_primer_reads"] == "12483"
    assert row["quality_trimmed_reads"] == "0"


def test_no_primers_in_bed_wording(tmp_path):
    text = (
        "Trimmed primers from 0% (0) of reads.\n"
        "0% (0) of reads were quality trimmed below the minimum length of 75 bp and were not written to file.\n"
        "100% (500) of reads started outside of primer regions. "
        "Since there were no primers found in BED file, these reads were written to file.\n"
    )
    row = run(tmp_path, text)
    assert (row["primer_trimmed_pct"], row["no_primer_reads"], row["status"]) == (
        "0.00",
        "500",
        "ok",
    )


def test_zero_reads_gives_na_percent(tmp_path):
    """ivar divides by zero for an empty BAM and prints nan/-nan."""
    text = (
        "Trimmed primers from -nan% (0) of reads.\n"
        "-nan% (0) of reads were quality trimmed below the minimum length of 75 bp and were not written to file.\n"
        "-nan% (0) of reads started outside of primer regions. Since the -e flag was given, these reads were written to file.\n"
    )
    row = run(tmp_path, text)
    assert row["primer_trimmed_pct"] == "NA"
    assert row["primer_trimmed_reads"] == "0"
    assert row["reads_total"] == "0"


def test_only_the_trimmed_line_is_enough(tmp_path):
    row = run(tmp_path, "Trimmed primers from 12.5% (5) of reads.\n")
    assert row["primer_trimmed_pct"] == "12.50"
    assert row["primer_trimmed_reads"] == "5"
    assert row["reads_total"] == "NA"
    assert row["no_primer_reads"] == "NA"
    assert row["status"] == "ok"


def test_stdout_text_with_windows_line_endings(tmp_path):
    row = run(tmp_path, "Trimmed primers from 30.00% (3) of reads.\r\n7 junk\r\n")
    assert row["primer_trimmed_pct"] == "30.00"


def test_missing_log_is_na(tmp_path, capsys):
    row = run(tmp_path, log_path=tmp_path / "nope.log")
    assert row["status"] == "missing"
    assert all(row[c] == "NA" for c in its.COLUMNS[1:-1])
    assert "note:" in capsys.readouterr().err


def test_no_log_argument_is_na(tmp_path):
    row = run(tmp_path)
    assert row["status"] == "missing"
    assert row["primer_trimmed_pct"] == "NA"


def test_empty_log_is_na(tmp_path, capsys):
    """The nf-core IVAR_TRIM log is empty unless stderr is redirected into it: say so."""
    row = run(tmp_path, "")
    assert row["status"] == "empty"
    assert row["primer_trimmed_pct"] == "NA"
    assert "stderr" in capsys.readouterr().err


@pytest.mark.parametrize(
    "text",
    [
        "garbage\n\x00\xff binary-ish\n",
        "Found 58 primers in BED file\nReading from x.bam\nResults:\nPrimer Name\tRead Count\n",
        "Trimmed primers from lots% (many) of reads\n",
        ("5% (3) of reads started outside of primer regions. Since the -e flag was given.\n"),
    ],
)
def test_garbled_log_is_na_and_exit_0(tmp_path, text):
    row = run(tmp_path, text)
    assert row["primer_trimmed_pct"] == "NA"
    assert row["status"] == "unparsed"


def test_binary_log_does_not_crash(tmp_path):
    path = tmp_path / "bin.log"
    path.write_bytes(bytes(range(256)) * 4)
    row = run(tmp_path, log_path=path)
    assert row["status"] in ("unparsed", "empty")
    assert row["primer_trimmed_pct"] == "NA"


def test_sample_id_is_not_rewritten(tmp_path):
    row = run(tmp_path, "Trimmed primers from 1.0% (1) of reads.\n", sample="Site_1_week_2")
    assert row["sample"] == "Site_1_week_2"


def test_stdout_output(tmp_path, capsys):
    log = tmp_path / "l.log"
    log.write_text("Trimmed primers from 50.00% (1) of reads.\n")
    assert its.main(["--sample", "S", "--log", str(log)]) == 0
    lines = capsys.readouterr().out.splitlines()
    assert lines[0].split("\t") == its.COLUMNS
    assert lines[1].split("\t")[3] == "50.00"


def test_sample_is_required():
    with pytest.raises(SystemExit):
        its.main(["--log", "x"])
