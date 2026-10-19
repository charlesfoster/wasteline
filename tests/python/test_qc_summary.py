from pathlib import Path

import pytest
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


def test_default_thresholds_are_the_wastewater_ones(tmp_path):
    """D17: PASS from 60 % of the genome at the depth, FAIL below 30 %, WARN in between (was 80 / 50)."""
    assert run_defaults(tmp_path, depth_profile(ge10=65))["qc_status"] == "PASS"
    assert run_defaults(tmp_path, depth_profile(ge10=45))["qc_status"] == "WARN"
    assert run_defaults(tmp_path, depth_profile(ge10=20))["qc_status"] == "FAIL"


def amplicon_file(tmp_path, dropped, total=20, sample="S"):
    path = tmp_path / "amp.tsv"
    lines = ["sample\tamplicon\tstart\tend\tmean_depth\tmin_depth\tpct_ge_min_depth\tdropped"]
    for i in range(total):
        flag = "true" if i < dropped else "false"
        lines.append(f"{sample}\ta_{i}\t1\t10\t5.00\t0\t50.00\t{flag}")
    path.write_text("\n".join(lines) + "\n")
    return str(path)


def test_column_order_matches_contract():
    contract = (
        "sample, sample_type, collection_date, site, raw_reads, trimmed_reads, mapped_reads, mapped_pct, mean_depth, "
        "depth_sd, median_depth, pct_genome_ge_1x, pct_genome_ge_10x, pct_genome_ge_min_depth, amplicons_total, "
        "amplicons_dropped, qc_status, qc_reasons, demix_status, demix_reason, consensus_length, consensus_n_pct, "
        "barcode_version, primer_trimmed_pct, primer_scheme_check, primer_scheme_best, primer_scheme_enrichment"
    )
    assert ", ".join(qc.COLUMNS) == contract
    # appended in this order: every older column keeps its place
    assert qc.COLUMNS[23:] == [
        "primer_trimmed_pct",
        "primer_scheme_check",
        "primer_scheme_best",
        "primer_scheme_enrichment",
    ]


def test_depth_statistics(tmp_path):
    row = run(
        tmp_path,
        [0, 5, 10, 20, 100],
        "--consensus-min-depth",
        "15",
        "--collection-date",
        "2026-01-02",
        "--site",
        "plant_A",
    )
    assert row["mean_depth"] == "27.00"
    assert row["median_depth"] == "10.00"
    assert row["depth_sd"] == "37.09"  # population SD
    assert row["pct_genome_ge_1x"] == "80.00"
    assert row["pct_genome_ge_10x"] == "60.00"
    assert row["pct_genome_ge_min_depth"] == "40.00"  # consensus-min-depth 15 -> 20 and 100
    assert row["collection_date"] == "2026-01-02"
    assert row["site"] == "plant_A"


def test_consensus_min_depth_defaults_to_10(tmp_path):
    row = run(tmp_path, [0, 5, 10, 20, 100])
    assert row["pct_genome_ge_min_depth"] == "60.00"  # depth >= 10: 10, 20 and 100


def test_optional_inputs_are_na(tmp_path):
    row = run(tmp_path, depth_profile())
    for col in (
        "raw_reads",
        "trimmed_reads",
        "mapped_reads",
        "mapped_pct",
        "amplicons_total",
        "amplicons_dropped",
        "primer_trimmed_pct",
        "primer_scheme_check",
        "primer_scheme_best",
        "primer_scheme_enrichment",
    ):
        assert row[col] == "NA"
    for col in (
        "collection_date",
        "site",
        "demix_status",
        "demix_reason",
        "consensus_length",
        "barcode_version",
    ):
        assert row[col] == "NA"
    assert row["qc_status"] == "PASS"
    assert row["qc_reasons"] == ""


def test_reads_from_fastp_and_flagstat(tmp_path):
    row = run(
        tmp_path,
        depth_profile(),
        "--flagstat",
        str(FIX / "flagstat_new.txt"),
        "--fastp-json",
        str(FIX / "fastp.json"),
    )
    assert row["raw_reads"] == "2400"
    assert row["trimmed_reads"] == "2200"
    assert row["mapped_reads"] == "1000"  # primary mapped, not 1100 incl. supplementary
    assert row["mapped_pct"] == "45.45"  # 1000 / 2200


def test_flagstat_old_format():
    fs = qc.parse_flagstat(str(FIX / "flagstat_old.txt"))
    assert fs == {"mapped": 1800, "total": 2000, "total_all": 2000}


def test_flagstat_new_format():
    fs = qc.parse_flagstat(str(FIX / "flagstat_new.txt"))
    assert fs["mapped"] == 1000 and fs["total"] == 1000


def test_mapped_pct_falls_back_to_flagstat_total(tmp_path):
    row = run(tmp_path, depth_profile(), "--flagstat", str(FIX / "flagstat_old.txt"))
    assert row["mapped_reads"] == "1800"
    assert row["mapped_pct"] == "90.00"


def test_bad_optional_inputs_do_not_crash(tmp_path):
    bad = tmp_path / "bad.json"
    bad.write_text("{not json")
    row = run(
        tmp_path,
        depth_profile(),
        "--flagstat",
        str(tmp_path / "missing.txt"),
        "--fastp-json",
        str(bad),
        "--amplicons",
        str(tmp_path / "missing.tsv"),
    )
    assert (
        row["mapped_reads"] == "NA" and row["raw_reads"] == "NA" and row["amplicons_total"] == "NA"
    )
    assert row["qc_status"] == "PASS"


def test_amplicon_counts(tmp_path):
    row = run(tmp_path, depth_profile(), "--amplicons", amplicon_file(tmp_path, dropped=3))
    assert row["amplicons_total"] == "20" and row["amplicons_dropped"] == "3"


def test_amplicon_counts_filtered_by_sample(tmp_path):
    other = amplicon_file(tmp_path, dropped=20, sample="OTHER")
    row = run(tmp_path, depth_profile(), "--amplicons", other)
    assert row["amplicons_total"] == "0" and row["amplicons_dropped"] == "0"


# ---- status: samples ----------------------------------------------------------------------------------------


def test_sample_pass_at_exactly_80_percent(tmp_path):
    row = run(tmp_path, depth_profile(ge10=80))
    assert (row["qc_status"], row["qc_reasons"]) == ("PASS", "")


def test_sample_warn_between_50_and_80(tmp_path):
    row = run(tmp_path, depth_profile(ge10=79))
    assert (row["qc_status"], row["qc_reasons"]) == ("WARN", "low_coverage")


def test_sample_warn_at_exactly_50(tmp_path):
    assert run(tmp_path, depth_profile(ge10=50))["qc_status"] == "WARN"


def test_sample_fail_below_50(tmp_path):
    row = run(tmp_path, depth_profile(ge10=49))
    assert (row["qc_status"], row["qc_reasons"]) == ("FAIL", "low_coverage")


def test_dropout_warns_above_max_dropped(tmp_path):
    row = run(
        tmp_path, depth_profile(), "--amplicons", amplicon_file(tmp_path, dropped=11, total=30)
    )
    assert (row["qc_status"], row["qc_reasons"]) == ("WARN", "amplicon_dropout")


def test_dropout_at_max_is_pass(tmp_path):
    row = run(
        tmp_path, depth_profile(), "--amplicons", amplicon_file(tmp_path, dropped=10, total=30)
    )
    assert row["qc_status"] == "PASS"


def test_max_dropped_is_configurable(tmp_path):
    amp = amplicon_file(tmp_path, dropped=3, total=30)
    assert (
        run(tmp_path, depth_profile(), "--amplicons", amp, "--max-dropped", "2")["qc_status"]
        == "WARN"
    )


def test_low_coverage_and_dropout_reasons_combine(tmp_path):
    row = run(
        tmp_path,
        depth_profile(ge10=60),
        "--amplicons",
        amplicon_file(tmp_path, dropped=15, total=30),
    )
    assert (row["qc_status"], row["qc_reasons"]) == ("WARN", "low_coverage;amplicon_dropout")
    row = run(
        tmp_path,
        depth_profile(ge10=30),
        "--amplicons",
        amplicon_file(tmp_path, dropped=15, total=30),
    )
    assert (row["qc_status"], row["qc_reasons"]) == ("FAIL", "low_coverage;amplicon_dropout")


def test_thresholds_are_options(tmp_path):
    row = run(tmp_path, depth_profile(ge10=90), "--pass-min-pct", "95", "--fail-below-pct", "92")
    assert (row["qc_status"], row["qc_reasons"]) == ("FAIL", "low_coverage")


def test_status_uses_min_depth_option(tmp_path):
    depths = [20] * 60 + [5] * 40  # 60% at depth 20, 100% at depth >= 5
    assert run(tmp_path, depths, "--min-depth", "5")["qc_status"] == "PASS"
    assert run(tmp_path, depths, "--min-depth", "10")["qc_status"] == "WARN"
    row = run(tmp_path, depths, "--min-depth", "5")
    assert row["pct_genome_ge_10x"] == "60.00"  # reported column stays fixed at 10x


def test_all_zero_depths_is_no_reads(tmp_path):
    row = run(tmp_path, [0] * 100)
    assert (row["qc_status"], row["qc_reasons"]) == ("FAIL", "no_reads")


def test_zero_mapped_reads_is_no_reads(tmp_path):
    fs = tmp_path / "fs.txt"
    fs.write_text("0 + 0 in total (QC-passed reads + QC-failed reads)\n0 + 0 mapped (N/A : N/A)\n")
    row = run(tmp_path, depth_profile(), "--flagstat", str(fs))
    assert (row["qc_status"], row["qc_reasons"]) == ("FAIL", "no_reads")
    assert row["mapped_reads"] == "0"


def test_unreadable_depths_fails_sample_with_a_row(tmp_path):
    out = tmp_path / "row.tsv"
    assert qc.main(["--sample", "S", "--depths", str(tmp_path / "missing"), "-o", str(out)]) == 0
    row = dict(zip(qc.COLUMNS, out.read_text().splitlines()[1].split("\t")))
    assert (row["qc_status"], row["qc_reasons"]) == ("FAIL", "no_depths")
    assert row["mean_depth"] == "NA"


# ---- status: negatives --------------------------------------------------------------------------------------


def test_negative_with_no_reads_passes(tmp_path):
    row = run(tmp_path, [0] * 100, sample_type="negative")
    assert (row["qc_status"], row["qc_reasons"]) == ("PASS", "")


def test_negative_never_low_coverage(tmp_path):
    row = run(tmp_path, [0] * 99 + [1], sample_type="negative")
    assert row["qc_status"] == "PASS" and "low_coverage" not in row["qc_reasons"]


def test_negative_fails_on_genome_coverage(tmp_path):
    row = run(tmp_path, depth_profile(ge10=0, ge1=11), sample_type="negative")
    assert (row["qc_status"], row["qc_reasons"]) == ("FAIL", "ntc_reads")


def test_negative_coverage_at_threshold_passes(tmp_path):
    assert (
        run(tmp_path, depth_profile(ge10=0, ge1=10), sample_type="negative")["qc_status"] == "PASS"
    )


def test_negative_fails_on_mapped_reads(tmp_path):
    fs = tmp_path / "fs.txt"
    fs.write_text(
        "5000 + 0 in total (QC-passed reads + QC-failed reads)\n5000 + 0 mapped (100.00% : N/A)\n"
    )
    row = run(tmp_path, [0] * 100, "--flagstat", str(fs), sample_type="negative")
    assert (row["qc_status"], row["qc_reasons"]) == ("FAIL", "ntc_reads")


def test_negative_mapped_reads_at_threshold_passes(tmp_path):
    fs = tmp_path / "fs.txt"
    fs.write_text(
        "1000 + 0 in total (QC-passed reads + QC-failed reads)\n1000 + 0 mapped (100.00% : N/A)\n"
    )
    assert (
        run(tmp_path, [0] * 100, "--flagstat", str(fs), sample_type="negative")["qc_status"]
        == "PASS"
    )


def test_negative_thresholds_are_options(tmp_path):
    row = run(
        tmp_path,
        depth_profile(ge10=0, ge1=11),
        "--ntc-max-genome-coverage-pct",
        "20",
        sample_type="negative",
    )
    assert row["qc_status"] == "PASS"


def test_negative_without_any_data_is_na(tmp_path):
    out = tmp_path / "row.tsv"
    assert (
        qc.main(
            [
                "--sample",
                "N",
                "--sample-type",
                "negative",
                "--depths",
                str(tmp_path / "x"),
                "-o",
                str(out),
            ]
        )
        == 0
    )
    row = dict(zip(qc.COLUMNS, out.read_text().splitlines()[1].split("\t")))
    assert (row["qc_status"], row["qc_reasons"]) == ("NA", "no_depths")


# ---- primer_trimmed_pct / low_primer_trimming ------------------------------------------------------------


def trim_stats(tmp_path, pct, sample="S", status="ok"):
    path = tmp_path / "trim.tsv"
    path.write_text(
        "sample\treads_total\tprimer_trimmed_reads\tprimer_trimmed_pct\tno_primer_reads\tno_primer_pct\t"
        "quality_trimmed_reads\tquality_trimmed_pct\tstatus\n"
        f"{sample}\t100\t20\t{pct}\t80\t80.00\t0\t0.00\t{status}\n"
    )
    return str(path)


def test_primer_trimmed_pct_is_reported(tmp_path):
    row = run(tmp_path, depth_profile(), "--ivar-trim-stats", trim_stats(tmp_path, "21.97"))
    assert row["primer_trimmed_pct"] == "21.97"
    assert row["qc_status"] == "PASS"
    assert row["qc_reasons"] == ""


def test_primer_trimmed_pct_real_log_through_stats_script(tmp_path):
    import ivar_trim_stats

    out = tmp_path / "stats.tsv"
    ivar_trim_stats.main(
        ["--sample", "S", "--log", str(FIX / "ivar_trim_PURE_C.ivar.log"), "-o", str(out)]
    )
    row = run(tmp_path, depth_profile(), "--ivar-trim-stats", str(out))
    assert row["primer_trimmed_pct"] == "21.97"


def test_low_primer_trimming_rule_is_off_by_default(tmp_path):
    row = run(tmp_path, depth_profile(), "--ivar-trim-stats", trim_stats(tmp_path, "0.50"))
    assert row["primer_trimmed_pct"] == "0.50"
    assert (row["qc_status"], row["qc_reasons"]) == ("PASS", "")


def test_low_primer_trimming_warns_below_threshold(tmp_path):
    row = run(
        tmp_path,
        depth_profile(),
        "--ivar-trim-stats",
        trim_stats(tmp_path, "4.99"),
        "--min-primer-trimmed-pct",
        "5",
    )
    assert (row["qc_status"], row["qc_reasons"]) == ("WARN", "low_primer_trimming")


def test_low_primer_trimming_boundary_is_not_a_warning(tmp_path):
    row = run(
        tmp_path,
        depth_profile(),
        "--ivar-trim-stats",
        trim_stats(tmp_path, "5.00"),
        "--min-primer-trimmed-pct",
        "5",
    )
    assert (row["qc_status"], row["qc_reasons"]) == ("PASS", "")


def test_low_primer_trimming_zero_percent_warns(tmp_path):
    row = run(
        tmp_path,
        depth_profile(),
        "--ivar-trim-stats",
        trim_stats(tmp_path, "0.00"),
        "--min-primer-trimmed-pct",
        "5",
    )
    assert row["qc_reasons"] == "low_primer_trimming"


@pytest.mark.parametrize("value", ["NA", "nan", "", "abc"])
def test_low_primer_trimming_ignores_unavailable_value(tmp_path, value):
    row = run(
        tmp_path,
        depth_profile(),
        "--ivar-trim-stats",
        trim_stats(tmp_path, value),
        "--min-primer-trimmed-pct",
        "5",
    )
    assert row["primer_trimmed_pct"] == "NA"
    assert (row["qc_status"], row["qc_reasons"]) == ("PASS", "")


def test_low_primer_trimming_combines_with_other_reasons(tmp_path):
    extra = ("--ivar-trim-stats", trim_stats(tmp_path, "1.00"), "--min-primer-trimmed-pct", "5")
    warn = run(tmp_path, depth_profile(ge10=70), *extra)  # low_coverage WARN + low_primer_trimming
    assert (warn["qc_status"], warn["qc_reasons"]) == ("WARN", "low_coverage;low_primer_trimming")
    both = run(
        tmp_path,
        depth_profile(ge10=70),
        *extra,
        "--max-dropped",
        "1",
        "--amplicons",
        amplicon_file(tmp_path, 5),
    )
    assert both["qc_reasons"] == "low_coverage;amplicon_dropout;low_primer_trimming"


def test_low_primer_trimming_never_downgrades_a_fail(tmp_path):
    extra = ("--ivar-trim-stats", trim_stats(tmp_path, "1.00"), "--min-primer-trimmed-pct", "5")
    row = run(tmp_path, depth_profile(ge10=30), *extra)
    assert row["qc_status"] == "FAIL"
    assert row["qc_reasons"] == "low_coverage;low_primer_trimming"


def test_low_primer_trimming_alone_never_fails(tmp_path):
    row = run(
        tmp_path,
        depth_profile(),
        "--ivar-trim-stats",
        trim_stats(tmp_path, "0.00"),
        "--min-primer-trimmed-pct",
        "100",
    )
    assert row["qc_status"] == "WARN"


def test_low_primer_trimming_not_applied_to_negatives(tmp_path):
    row = run(
        tmp_path,
        [0] * 100,
        "--ivar-trim-stats",
        trim_stats(tmp_path, "0.00"),
        "--min-primer-trimmed-pct",
        "5",
        sample_type="negative",
    )
    assert (row["qc_status"], row["qc_reasons"]) == ("PASS", "")
    assert row["primer_trimmed_pct"] == "0.00"  # still reported


def test_low_primer_trimming_kept_beside_no_depths(tmp_path):
    out = tmp_path / "row.tsv"
    argv = [
        "--sample",
        "S",
        "--depths",
        str(tmp_path / "missing"),
        "--ivar-trim-stats",
        trim_stats(tmp_path, "1.00"),
        "--min-primer-trimmed-pct",
        "5",
        "-o",
        str(out),
    ]
    assert qc.main(argv) == 0
    row = dict(zip(qc.COLUMNS, out.read_text().splitlines()[1].split("\t")))
    assert (row["qc_status"], row["qc_reasons"]) == ("FAIL", "no_depths;low_primer_trimming")


def test_stats_row_is_chosen_by_sample(tmp_path):
    path = tmp_path / "trim.tsv"
    path.write_text("sample\tprimer_trimmed_pct\nOTHER\t90.00\nS\t12.50\n")
    row = run(tmp_path, depth_profile(), "--ivar-trim-stats", str(path))
    assert row["primer_trimmed_pct"] == "12.50"


@pytest.mark.parametrize(
    "content", ["", "sample\tprimer_trimmed_pct\n", "sample\tstatus\nS\tok\n", "garbage"]
)
def test_bad_ivar_trim_stats_do_not_crash(tmp_path, content):
    path = tmp_path / "trim.tsv"
    path.write_text(content)
    row = run(
        tmp_path,
        depth_profile(),
        "--ivar-trim-stats",
        str(path),
        "--min-primer-trimmed-pct",
        "5",
    )
    assert row["primer_trimmed_pct"] == "NA"
    assert row["qc_reasons"] == ""


def test_missing_ivar_trim_stats_file_does_not_crash(tmp_path):
    row = run(tmp_path, depth_profile(), "--ivar-trim-stats", str(tmp_path / "nope.tsv"))
    assert row["primer_trimmed_pct"] == "NA"


# ---- primer_scheme_* columns / primer_scheme_mismatch ---------------------------------------------------

SUMMARY_HEADER = "sample\tdeclared_scheme\tdeclared_enrichment\tbest_scheme\tbest_enrichment\tscheme_check\treason\n"


def scheme_summary(tmp_path, check, declared="1.1", best="midnight/midnightv2", sample="S"):
    path = tmp_path / "scheme.tsv"
    path.write_text(SUMMARY_HEADER + f"{sample}\tarticv4.1\t{declared}\t{best}\t41.6\t{check}\tx\n")
    return str(path)


def test_scheme_columns_are_reported(tmp_path):
    row = run(
        tmp_path,
        depth_profile(),
        "--primer-scheme-summary",
        scheme_summary(tmp_path, "ok", "41.64"),
    )
    assert row["primer_scheme_check"] == "ok"
    assert row["primer_scheme_best"] == "midnight/midnightv2"
    assert row["primer_scheme_enrichment"] == "41.6"
    assert (row["qc_status"], row["qc_reasons"]) == ("PASS", "")


def test_scheme_mismatch_warns(tmp_path):
    row = run(
        tmp_path, depth_profile(), "--primer-scheme-summary", scheme_summary(tmp_path, "mismatch")
    )
    assert row["primer_scheme_check"] == "mismatch"
    assert (row["qc_status"], row["qc_reasons"]) == ("WARN", "primer_scheme_mismatch")


@pytest.mark.parametrize("check", ["ok", "uninformative", "NA"])
def test_other_scheme_checks_do_not_warn(tmp_path, check):
    row = run(tmp_path, depth_profile(), "--primer-scheme-summary", scheme_summary(tmp_path, check))
    assert row["primer_scheme_check"] == check
    assert (row["qc_status"], row["qc_reasons"]) == ("PASS", "")


def test_scheme_mismatch_na_enrichment_and_best(tmp_path):
    path = tmp_path / "s.tsv"
    path.write_text(SUMMARY_HEADER + "S\tmidnight\tNA\tNA\tNA\tNA\ttoo_few_reads\n")
    row = run(tmp_path, depth_profile(), "--primer-scheme-summary", str(path))
    assert (
        row["primer_scheme_check"],
        row["primer_scheme_best"],
        row["primer_scheme_enrichment"],
    ) == (
        "NA",
        "NA",
        "NA",
    )


def test_scheme_mismatch_combines_and_never_downgrades(tmp_path):
    summary = scheme_summary(tmp_path, "mismatch")
    warn = run(tmp_path, depth_profile(ge10=70), "--primer-scheme-summary", summary)
    assert (warn["qc_status"], warn["qc_reasons"]) == (
        "WARN",
        "low_coverage;primer_scheme_mismatch",
    )
    fail = run(tmp_path, depth_profile(ge10=30), "--primer-scheme-summary", summary)
    assert (fail["qc_status"], fail["qc_reasons"]) == (
        "FAIL",
        "low_coverage;primer_scheme_mismatch",
    )
    both = run(
        tmp_path,
        depth_profile(ge10=70),
        "--primer-scheme-summary",
        summary,
        "--ivar-trim-stats",
        trim_stats(tmp_path, "1.00"),
        "--min-primer-trimmed-pct",
        "5",
    )
    assert both["qc_reasons"] == "low_coverage;low_primer_trimming;primer_scheme_mismatch"


def test_scheme_mismatch_not_applied_to_negatives(tmp_path):
    row = run(
        tmp_path,
        [0] * 100,
        "--primer-scheme-summary",
        scheme_summary(tmp_path, "mismatch"),
        sample_type="negative",
    )
    assert (row["qc_status"], row["qc_reasons"]) == ("PASS", "")
    assert row["primer_scheme_check"] == "mismatch"  # still reported


def test_scheme_mismatch_kept_beside_no_depths(tmp_path):
    out = tmp_path / "row.tsv"
    argv = [
        "--sample",
        "S",
        "--depths",
        str(tmp_path / "missing"),
        "--primer-scheme-summary",
        scheme_summary(tmp_path, "mismatch"),
        "-o",
        str(out),
    ]
    assert qc.main(argv) == 0
    row = dict(zip(qc.COLUMNS, out.read_text().splitlines()[1].split("\t")))
    assert (row["qc_status"], row["qc_reasons"]) == ("FAIL", "no_depths;primer_scheme_mismatch")


def test_scheme_summary_row_is_chosen_by_sample(tmp_path):
    path = tmp_path / "s.tsv"
    path.write_text(
        SUMMARY_HEADER + "OTHER\tx\t1.0\tx\t9.0\tmismatch\tx\n" + "S\tx\t30.0\tx\t30.0\tok\tx\n"
    )
    row = run(tmp_path, depth_profile(), "--primer-scheme-summary", str(path))
    assert row["primer_scheme_check"] == "ok"


@pytest.mark.parametrize(
    "content",
    [
        "",
        SUMMARY_HEADER,
        "sample\tx\nS\ty\n",
        SUMMARY_HEADER + "S\ta\t1\tb\t2\tweird\tx\n",
    ],
)
def test_bad_scheme_summary_does_not_crash(tmp_path, content):
    path = tmp_path / "s.tsv"
    path.write_text(content)
    row = run(tmp_path, depth_profile(), "--primer-scheme-summary", str(path))
    assert row["primer_scheme_check"] == "NA"
    assert row["qc_reasons"] == ""


def test_missing_scheme_summary_does_not_crash(tmp_path):
    row = run(tmp_path, depth_profile(), "--primer-scheme-summary", str(tmp_path / "nope.tsv"))
    assert row["primer_scheme_check"] == "NA"


# ---- concat -------------------------------------------------------------------------------------------------


def make_row(tmp_path, sample):
    out = tmp_path / f"{sample}.tsv"
    dep = write_depths(tmp_path / f"{sample}.depths", depth_profile())
    assert qc.main(["--sample", sample, "--depths", dep, "-o", str(out)]) == 0
    return str(out)


def test_concat_sorts_by_sample(tmp_path):
    files = [make_row(tmp_path, s) for s in ("S10", "S2", "A_1")]
    out = tmp_path / "all.tsv"
    assert qc.main(["--concat", *files, "-o", str(out)]) == 0
    lines = out.read_text().splitlines()
    assert lines[0].split("\t") == qc.COLUMNS
    assert [line.split("\t")[0] for line in lines[1:]] == ["A_1", "S10", "S2"]


def test_concat_rejects_header_mismatch(tmp_path):
    good = make_row(tmp_path, "S1")
    bad = tmp_path / "bad.tsv"
    bad.write_text("sample\tqc_status\nS2\tPASS\n")
    assert qc.main(["--concat", good, str(bad), "-o", str(tmp_path / "o.tsv")]) == 2


@pytest.mark.parametrize("n_columns", [23, 24])
def test_concat_rejects_tables_from_before_the_scheme_check_columns(tmp_path, n_columns):
    """Older tables (23 columns: before primer_trimmed_pct; 24: before the primer_scheme_* columns) are
    rejected loudly, not padded with NA."""
    good = make_row(tmp_path, "S1")
    old = tmp_path / "old.tsv"
    cols = qc.COLUMNS[:n_columns]
    old.write_text("\t".join(cols) + "\n" + "\t".join(["S2"] + ["NA"] * (len(cols) - 1)) + "\n")
    assert qc.main(["--concat", good, str(old), "-o", str(tmp_path / "o.tsv")]) == 2


def test_concat_rejects_duplicates(tmp_path):
    good = make_row(tmp_path, "S1")
    assert qc.main(["--concat", good, good, "-o", str(tmp_path / "o.tsv")]) == 2


def test_concat_keeps_multiple_rows_per_file(tmp_path):
    a, b = make_row(tmp_path, "S1"), make_row(tmp_path, "S2")
    merged = tmp_path / "m.tsv"
    assert qc.main(["--concat", a, b, "-o", str(merged)]) == 0
    again = tmp_path / "again.tsv"
    assert qc.main(["--concat", str(merged), "-o", str(again)]) == 0
    assert merged.read_text() == again.read_text()


@pytest.mark.parametrize("bad", [["--sample", "S"], ["--depths", "x"]])
def test_requires_sample_and_depths(bad):
    with pytest.raises(SystemExit):
        qc.main(bad)
