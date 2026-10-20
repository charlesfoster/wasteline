from pathlib import Path

import qc_summary as qc

FIX = Path(__file__).resolve().parents[1] / "fixtures" / "qc"


def write_depths(path, depths, contig="ref"):
    path.write_text("".join(f"{contig}\t{i}\tA\t{d}\n" for i, d in enumerate(depths, start=1)))
    return str(path)


def depth_profile(n=100, ge10=100, ge1=None):
    """n positions, `ge10` percent at depth 20, the rest zero (or depth 5 for the first `ge1` percent)."""
    good = round(n * ge10 / 100)
    mid = round(n * ge1 / 100) - good if ge1 else 0
    return [20] * good + [5] * mid + [0] * (n - good - mid)


def run(tmp_path, depths, *extra, sample="S", sample_type="sample"):
    dep = write_depths(tmp_path / "d.depths", depths)
    out = tmp_path / "row.tsv"
    argv = [
        "--sample",
        sample,
        "--sample-type",
        sample_type,
        "--depths",
        dep,
        "-o",
        str(out),
        # the thresholds most tests were written for (the script's own defaults changed with D17: 80/50 -> 60/30);
        # an option in `extra` wins
        "--pass-min-pct",
        "80",
        "--fail-below-pct",
        "50",
        *extra,
    ]
    assert qc.main(argv) == 0
    header, row = out.read_text().splitlines()
    assert header.split("\t") == qc.COLUMNS
    return dict(zip(qc.COLUMNS, row.split("\t")))


def run_defaults(tmp_path, depths):
    """The script with ITS OWN thresholds (no overrides)."""
    dep = write_depths(tmp_path / "d.depths", depths)
    out = tmp_path / "row.tsv"
    assert qc.main(["--sample", "S", "--depths", dep, "-o", str(out)]) == 0
    header, row = out.read_text().splitlines()
    return dict(zip(header.split("\t"), row.split("\t")))


def amplicon_file(tmp_path, dropped, total=20, sample="S"):
    path = tmp_path / "amp.tsv"
    lines = ["sample\tamplicon\tstart\tend\tmean_depth\tmin_depth\tpct_ge_min_depth\tdropped"]
    for i in range(total):
        flag = "true" if i < dropped else "false"
        lines.append(f"{sample}\ta_{i}\t1\t10\t5.00\t0\t50.00\t{flag}")
    path.write_text("\n".join(lines) + "\n")
    return str(path)


def test_golden_row_depth_reads_and_amplicons(tmp_path):
    row = run(
        tmp_path,
        [0, 5, 10, 20, 100],
        "--consensus-min-depth",
        "15",
        "--collection-date",
        "2026-01-02",
        "--site",
        "plant_A",
        "--flagstat",
        str(FIX / "flagstat_new.txt"),
        "--fastp-json",
        str(FIX / "fastp.json"),
        "--amplicons",
        amplicon_file(tmp_path, dropped=3),
    )
    assert row["mean_depth"] == "27.00"
    assert row["median_depth"] == "10.00"
    assert row["depth_sd"] == "37.09"  # population SD
    assert row["pct_genome_ge_1x"] == "80.00"
    assert row["pct_genome_ge_10x"] == "60.00"
    assert row["pct_genome_ge_min_depth"] == "40.00"  # consensus-min-depth 15 -> 20 and 100
    assert row["collection_date"] == "2026-01-02"
    assert row["site"] == "plant_A"
    assert row["raw_reads"] == "2400"
    assert row["trimmed_reads"] == "2200"
    assert row["mapped_reads"] == "1000"  # primary mapped, not 1100 incl. supplementary
    assert row["mapped_pct"] == "45.45"  # 1000 / 2200
    assert row["amplicons_total"] == "20" and row["amplicons_dropped"] == "3"


def test_status_thresholds_and_reasons(tmp_path):
    """D17: PASS from 60 % of the genome at the depth, FAIL below 30 %, WARN in between (script defaults)."""
    assert run_defaults(tmp_path, depth_profile(ge10=65))["qc_status"] == "PASS"
    assert run_defaults(tmp_path, depth_profile(ge10=45))["qc_status"] == "WARN"
    assert run_defaults(tmp_path, depth_profile(ge10=20))["qc_status"] == "FAIL"
    # low coverage and amplicon dropout both appear in the reasons; dropout alone only warns
    amp = amplicon_file(tmp_path, dropped=15, total=30)
    row = run(tmp_path, depth_profile(ge10=30), "--amplicons", amp)
    assert (row["qc_status"], row["qc_reasons"]) == ("FAIL", "low_coverage;amplicon_dropout")
    row = run(tmp_path, depth_profile(), "--amplicons", amp)
    assert (row["qc_status"], row["qc_reasons"]) == ("WARN", "amplicon_dropout")


def test_failed_or_empty_samples_still_get_a_row_with_status_and_reason(tmp_path):
    # no depths file at all
    out = tmp_path / "row.tsv"
    assert qc.main(["--sample", "S", "--depths", str(tmp_path / "missing"), "-o", str(out)]) == 0
    row = dict(zip(qc.COLUMNS, out.read_text().splitlines()[1].split("\t")))
    assert (row["qc_status"], row["qc_reasons"]) == ("FAIL", "no_depths")
    assert row["mean_depth"] == "NA"
    # a sample with no coverage anywhere fails as no_reads; a negative control with no reads passes
    row = run(tmp_path, [0] * 100)
    assert (row["qc_status"], row["qc_reasons"]) == ("FAIL", "no_reads")
    row = run(tmp_path, [0] * 100, sample_type="negative")
    assert (row["qc_status"], row["qc_reasons"]) == ("PASS", "")


def test_primer_trimmed_pct_real_log_through_stats_script(tmp_path):
    import ivar_trim_stats

    out = tmp_path / "stats.tsv"
    ivar_trim_stats.main(
        ["--sample", "S", "--log", str(FIX / "ivar_trim_PURE_C.ivar.log"), "-o", str(out)]
    )
    row = run(tmp_path, depth_profile(), "--ivar-trim-stats", str(out))
    assert row["primer_trimmed_pct"] == "21.97"
