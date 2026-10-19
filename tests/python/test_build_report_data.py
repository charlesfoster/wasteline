import json
import shutil
from pathlib import Path

import build_report_data as brd
import pytest

ROOT = Path(__file__).resolve().parents[2]
FIX = ROOT / "tests" / "fixtures" / "report"


def write_depths(path, scale, zero_ranges=(), length=29903):
    """4-column, headerless Freyja depths file with a flat depth (and optional dropouts)."""
    lines = []
    for pos in range(1, length + 1):
        depth = 0 if any(a <= pos <= b for a, b in zero_ranges) else scale
        lines.append(f"NC_045512.2\t{pos}\tA\t{depth}")
    Path(path).write_text("\n".join(lines) + "\n")


def run_builder(fixture, out, extra=None, files=None):
    """Run the CLI over a fixture directory. `files` maps option -> file name."""
    base = FIX / fixture
    mapping = {
        "--qc": "qc_summary.tsv",
        "--amplicons": "amplicon_coverage.tsv",
        "--aggregated": "freyja_aggregated.tsv",
        "--demix-status": "demix_status.tsv",
        "--moi": "mutations_of_interest.tsv",
        "--consensus": "consensus_stats.tsv",
        "--barcode-compare": "barcode_compare.tsv",
    }
    args = ["--run-info", str(base / "run_info.json")]
    for option, name in (files or mapping).items():
        if (base / name).exists():
            args += [option, str(base / name)]
    args += list(extra or []) + ["--out", str(out)]
    assert brd.main(args) == 0
    return json.loads(Path(out).read_text())


def flag_codes(data, level=None):
    return {f["code"] for f in data["run_flags"] if level is None or f["level"] == level}


# --------------------------------------------------------------------------- parsers
def test_real_freyja_aggregated_numpy_repr(tmp_path):
    """The real Freyja 2.0.5 example: numpy-array reprs without commas, quoted newlines."""
    problems = brd.Problems()
    parsed = brd.load_aggregated(FIX / "freyja_real_aggregated.tsv", set(), problems)
    assert not problems.items
    assert set(parsed) == {"mixture", "test"}  # '.tsv' stripped when no sample matches
    mix = parsed["mixture"]
    names = [n for n, _ in mix["lineages"]]
    assert len(names) == 10
    assert names[:3] == ["A", "Q.3", "AY.48"]  # literal_eval would have glued 'A' 'Q.3' -> 'AQ.3'
    assert "AQ.3" not in names
    assert mix["lineages"][0][1] == pytest.approx(0.36664, abs=1e-5)
    assert mix["groups"][0] == ["A", pytest.approx(0.36664, abs=1e-5)]
    assert [g for g, _ in mix["groups"]][:3] == ["A", "Alpha", "Delta"]
    assert mix["resid"] == pytest.approx(6.1172, abs=1e-3)
    assert mix["coverage"] is None  # column absent in this Freyja example
    assert len(parsed["test"]["lineages"]) == 18


def test_parse_name_list_variants():
    assert brd.parse_name_list("['A' 'Q.3'\n 'B.1']") == ["A", "Q.3", "B.1"]
    assert brd.parse_name_list("['A', 'Q.3']") == ["A", "Q.3"]
    # what Freyja 2.0.5 actually writes in a .demix file and `freyja aggregate` passes on
    assert brd.parse_name_list("BA.2 AY.4") == ["BA.2", "AY.4"]
    assert brd.parse_name_list("JN.1") == ["JN.1"]
    assert brd.parse_name_list("") == []
    assert brd.parse_name_list("NA") == []


def test_parse_float_list_variants():
    assert brd.parse_float_list("[0.5  6.0e-03\n 1]") == [0.5, 0.006, 1.0]


def test_parse_summarized_literal_and_fallback():
    assert brd.parse_summarized("[('Alpha', 0.5), ('Other', 0.5)]") == [
        ["Alpha", 0.5],
        ["Other", 0.5],
    ]
    numpy2 = "[('Alpha', np.float64(0.25)), ('Other', np.float64(0.75))]"
    assert brd.parse_summarized(numpy2) == [["Alpha", 0.25], ["Other", 0.75]]
    assert brd.parse_summarized("garbage __import__('os')") == []
    assert brd.parse_summarized("") == []


def test_summarized_never_evaluates_code(tmp_path):
    marker = tmp_path / "pwned"
    evil = f"[(__import__('pathlib').Path('{marker}').write_text('x'), 1.0)]"
    brd.parse_summarized(evil)
    assert not marker.exists()


def test_to_int_float_na():
    assert brd.to_int("12") == 12 and brd.to_int("12.0") == 12 and brd.to_int("NA") is None
    assert brd.to_float("0.5") == 0.5 and brd.to_float("") is None and brd.to_float("nan") is None
    with pytest.raises(ValueError):
        brd.to_int("abc")


# --------------------------------------------------------------------------- fixtures
def test_full_run(tmp_path):
    depth_dir = tmp_path / "depths"
    depth_dir.mkdir()
    for sample, scale in (("WW_001", 800), ("WW_003.rep1", 400), ("NTC_01", 0)):
        write_depths(depth_dir / f"{sample}.depths", scale)
    data = run_builder("a_full", tmp_path / "r.json", ["--depths", str(depth_dir)])
    assert data["schema_version"] == 1 and data["mode"] == "full"
    assert data["header"]["barcode_version"] == "09_17_2026-00-46"
    assert data["header"]["git_commit"].startswith("3f9c2ab")
    assert data["provenance"]["freyja_data"]["barcode_version"] == "09_17_2026-00-46"
    assert data["provenance"]["previous_runs"] == []
    assert data["params_summary"]["primer_scheme"] == "midnight"
    assert data["thresholds"]["ntc_max_mapped_reads"] == 100
    assert [s["sample"] for s in data["samples"]][:2] == ["NTC_01", "WW_001"]
    s = {x["sample"]: x for x in data["samples"]}
    # typed values
    assert s["WW_001"]["raw_reads"] == 250000 and isinstance(s["WW_001"]["raw_reads"], int)
    assert s["WW_001"]["mean_depth"] == 812.4
    assert s["NTC_01"]["site"] is None and s["NTC_01"]["sample_type"] == "negative"
    assert s["WW_004"]["qc_reasons"] == ["low_coverage", "amplicon_dropout"]
    assert s["WW_001"]["qc_reasons"] == []
    # filled from the other files
    assert (
        s["WW_001"]["demix_status"] == "ok" and s["WW_001"]["barcode_version"] == "09_17_2026-00-46"
    )
    assert s["WW_004"]["demix_status"] == "skipped_low_coverage"
    assert "30.10" in s["WW_004"]["demix_reason"]
    assert s["WW_005"]["demix_status"] == "failed"
    assert s["WW_002"]["consensus_length"] == 29903 and s["WW_003.rep1"]["consensus_n_pct"] == 3.01
    # demix
    d = data["demix"]
    assert d["WW_001"]["groups"][0] == ["Omicron", 0.55]
    assert d["WW_001"]["lineages"][0] == ["BA.2.86", 0.4]
    assert d["WW_001"]["coverage"] == 99.9
    assert d["WW_004"]["status"] == "skipped_low_coverage" and d["WW_004"]["groups"] == []
    assert d["NTC_01"]["status"] is None
    # amplicons / tracks / others
    assert len(data["amplicons"]["WW_003.rep1"]) == 98
    assert [a["amplicon"] for a in data["amplicons"]["WW_001"]][:2] == ["midnight_1", "midnight_2"]
    tracks = data["coverage_tracks"]
    assert tracks["genome_length"] == 29903 and tracks["n_bins"] <= 600
    assert set(tracks["samples"]) == {"WW_001", "WW_003.rep1", "NTC_01"}
    assert len(tracks["samples"]["WW_001"]) == tracks["n_bins"]
    assert tracks["samples"]["WW_001"][5] == 800 and tracks["samples"]["NTC_01"][5] == 0
    assert len(data["mutations_of_interest"]) == 6
    assert data["mutations_of_interest"][4]["alt_freq"] is None
    assert len(data["consensus"]) == 6
    assert data["barcode_comparison"]["label_b"] == "08_01_2026-00-10"
    assert data["barcode_comparison"]["label_a"] == "09_17_2026-00-46"
    assert data["barcode_comparison"]["rows"][0]["delta"] == pytest.approx(0.1)
    ts = data["time_series"]
    assert ts["n_dates"] == 2 and ts["groups"][0] == "Omicron"
    site_a = next(x for x in ts["sites"] if x["site"] == "SiteA")
    assert site_a["dates"] == ["2026-09-01", "2026-09-08"]
    # flags
    assert flag_codes(data, "error") == {"demix_failed"}
    assert {"qc_fail", "demix_skipped", "qc_warn"} <= flag_codes(data)
    assert data["run_flags"][0]["level"] == "error"
    assert all(set(f) >= {"level", "code", "message"} for f in data["run_flags"])
    # only on when the feature produced data
    extras = {
        "negative_controls",
        "boot_intervals",
        "dropout_impact",
        "unexplained_mutations",
        "timeseries",
        "growth_rates",
        "lineage_groups",
        "site_aggregates",
    }
    assert all(v for k, v in data["available"].items() if k not in extras)


def test_demix_only_run(tmp_path):
    data = run_builder("b_demix_only", tmp_path / "r.json")
    assert data["mode"] == "demix"
    assert [s["sample"] for s in data["samples"]] == ["S1", "S2"]  # 'S1.variants.tsv' resolved
    assert data["demix"]["S1"]["status"] == "ok"
    assert data["demix"]["S1"]["groups"][0] == ["Delta", 0.6]
    assert data["demix"]["S2"]["status"] == "failed"
    assert data["amplicons"] == {} and data["coverage_tracks"] is None
    assert data["consensus"] == [] and data["mutations_of_interest"] is None
    assert data["barcode_comparison"] is None and data["time_series"] is None
    assert data["available"]["qc"] is False
    assert data["samples"][0]["qc_status"] is None
    assert {"qc_not_available", "demix_failed", "previous_run_barcode_unknown"} <= flag_codes(data)
    assert data["provenance"]["previous_runs"][0]["barcode_version"] == "unknown"


def test_legacy_unknown_barcode(tmp_path):
    data = run_builder("c_legacy_barcode", tmp_path / "r.json")
    assert data["header"]["barcode_version"] == "unknown"
    codes = flag_codes(data)
    assert "barcode_unknown" in codes and "previous_run_barcode_unknown" in codes
    assert data["samples"][0]["barcode_version"] == "unknown"


def test_all_failed(tmp_path):
    data = run_builder("d_all_failed", tmp_path / "r.json")
    errors = flag_codes(data, "error")
    assert {"ntc_reads", "all_samples_failed"} <= errors
    ntc = next(f for f in data["run_flags"] if f["code"] == "ntc_reads")
    assert "NTC_01" in ntc["message"] and ntc["samples"] == ["NTC_01"]
    assert {s["sample"] for s in data["samples"]} == {"NTC_01", "F1", "F2"}
    assert data["demix"]["F1"]["status"] == "skipped_low_coverage"
    assert data["demix"]["F1"]["reason"] == "no reads"
    # every non-control sample has a status and reason (never an empty-looking result)
    assert all(data["demix"][s]["reason"] for s in ("F1", "F2"))


def test_zero_samples(tmp_path):
    data = run_builder("e_zero_samples", tmp_path / "r.json")
    assert data["samples"] == [] and data["demix"] == {}
    assert "no_samples" in flag_codes(data, "warn")
    assert data["time_series"] is None


def test_realistic_qc_scripts_output(tmp_path):
    """qc_summary.tsv / amplicon_coverage.tsv exactly as written by bin/qc_summary.py."""
    data = run_builder("f_real_qc", tmp_path / "r.json")
    s = data["samples"][0]
    assert s["sample"] == "REAL_1" and s["qc_status"] == "PASS" and s["raw_reads"] is None
    assert s["amplicons_total"] == 99 and s["pct_genome_ge_min_depth"] == 99.7
    assert len(data["amplicons"]["REAL_1"]) == 99
    assert data["amplicons"]["REAL_1"][0]["start"] >= 1


def test_ntc_threshold_from_params(tmp_path):
    base = FIX / "a_full"
    info = json.loads((base / "run_info.json").read_text())
    info["params"]["ntc_max_mapped_reads"] = 5
    (tmp_path / "ri.json").write_text(json.dumps(info))
    out = tmp_path / "r.json"
    args = [
        "--run-info",
        str(tmp_path / "ri.json"),
        "--qc",
        str(base / "qc_summary.tsv"),
        "--out",
        str(out),
    ]
    assert brd.main(args) == 0
    data = json.loads(out.read_text())
    assert "ntc_reads" in flag_codes(data, "error")


def test_missing_freyja_data_is_an_error_flag(tmp_path):
    info = json.loads((FIX / "b_demix_only" / "run_info.json").read_text())
    del info["freyja_data"]
    ri = tmp_path / "ri.json"
    ri.write_text(json.dumps(info))
    out = tmp_path / "r.json"
    assert (
        brd.main(
            [
                "--run-info",
                str(ri),
                "--demix-status",
                str(FIX / "b_demix_only" / "demix_status.tsv"),
                "--out",
                str(out),
            ]
        )
        == 0
    )
    data = json.loads(out.read_text())
    assert "freyja_data_missing" in flag_codes(data, "error")
    assert data["header"]["barcode_version"] is None


def test_mixed_barcode_versions(tmp_path):
    shutil.copy(FIX / "b_demix_only" / "run_info.json", tmp_path / "ri.json")
    status = tmp_path / "ds.tsv"
    status.write_text("sample\tstatus\treason\tbarcode_version\nA\tok\t\tv1\nB\tok\t\tv2\n")
    agg = tmp_path / "agg.tsv"
    agg.write_text(
        "\tsummarized\tlineages\tabundances\tresid\n"
        "A\t[('Delta', 1.0)]\t\"['AY.4']\"\t\"[1.0]\"\t1\nB\t[('Delta', 1.0)]\t\"['AY.4']\"\t\"[1.0]\"\t1\n"
    )
    out = tmp_path / "r.json"
    assert (
        brd.main(
            [
                "--run-info",
                str(tmp_path / "ri.json"),
                "--demix-status",
                str(status),
                "--aggregated",
                str(agg),
                "--out",
                str(out),
            ]
        )
        == 0
    )
    assert "mixed_barcode_versions" in flag_codes(json.loads(out.read_text()), "error")


def test_sample_without_demix_result_is_flagged_missing(tmp_path):
    qc = FIX / "d_all_failed" / "qc_summary.tsv"
    status = tmp_path / "ds.tsv"
    status.write_text("sample\tstatus\treason\tbarcode_version\nF1\tok\t\tx\n")
    out = tmp_path / "r.json"
    args = [
        "--run-info",
        str(FIX / "d_all_failed" / "run_info.json"),
        "--qc",
        str(qc),
        "--demix-status",
        str(status),
        "--out",
        str(out),
    ]
    assert brd.main(args) == 0
    data = json.loads(out.read_text())
    # F1 'ok' but without lineages -> failed with a reason; F2 absent -> missing; NTC untouched
    assert data["demix"]["F1"]["status"] == "failed" and data["demix"]["F1"]["reason"]
    assert data["demix"]["F2"]["status"] == "missing" and data["demix"]["F2"]["reason"]
    assert data["demix"]["NTC_01"]["status"] is None


def test_bad_values_degrade_to_flag(tmp_path):
    qc = tmp_path / "qc.tsv"
    qc.write_text(
        "sample\tsample_type\traw_reads\tmean_depth\tqc_status\nX\tsample\tmany\t12.5\tPASS\n"
    )
    out = tmp_path / "r.json"
    assert (
        brd.main(
            [
                "--run-info",
                str(FIX / "a_full" / "run_info.json"),
                "--qc",
                str(qc),
                "--out",
                str(out),
            ]
        )
        == 0
    )
    data = json.loads(out.read_text())
    assert data["samples"][0]["raw_reads"] is None and data["samples"][0]["mean_depth"] == 12.5
    assert "input_problems" in flag_codes(data, "warn")
    assert data["samples"][0]["site"] is None  # columns missing from the file are still present


def test_run_info_is_mandatory(tmp_path, capsys):
    assert (
        brd.main(["--run-info", str(tmp_path / "nope.json"), "--out", str(tmp_path / "o.json")])
        == 2
    )
    (tmp_path / "bad.json").write_text("[1,2]")
    assert (
        brd.main(["--run-info", str(tmp_path / "bad.json"), "--out", str(tmp_path / "o.json")]) == 2
    )
    assert "run_info.json" in capsys.readouterr().err


def test_missing_optional_inputs_and_no_file_placeholder(tmp_path):
    out = tmp_path / "r.json"
    (tmp_path / "NO_FILE").write_text("")
    args = [
        "--run-info",
        str(FIX / "a_full" / "run_info.json"),
        "--qc",
        str(FIX / "a_full" / "qc_summary.tsv"),
        "--moi",
        str(tmp_path / "NO_FILE"),
        "--amplicons",
        str(tmp_path / "absent.tsv"),
        "--out",
        str(out),
    ]
    assert brd.main(args) == 0
    data = json.loads(out.read_text())
    assert data["mutations_of_interest"] is None and data["amplicons"] == {}


# --------------------------------------------------------------------------- depths
def test_depth_inputs_dir_list_and_string(tmp_path):
    for name in ("A", "B"):
        write_depths(tmp_path / f"{name}.depths", 100, length=3000)
    out = tmp_path / "r.json"
    ri = str(FIX / "a_full" / "run_info.json")
    qc = tmp_path / "qc.tsv"
    qc.write_text("sample\tqc_status\nA\tPASS\nB\tPASS\n")
    for depths in (
        [str(tmp_path)],
        [str(tmp_path / "A.depths"), str(tmp_path / "B.depths")],
        [f"{tmp_path / 'A.depths'} {tmp_path / 'B.depths'}"],
    ):
        assert (
            brd.main(["--run-info", ri, "--qc", str(qc), "--depths", *depths, "--out", str(out)])
            == 0
        )
        tracks = json.loads(out.read_text())["coverage_tracks"]
        assert set(tracks["samples"]) == {"A", "B"} and tracks["genome_length"] == 3000
        assert tracks["bin_size"] == 5 and tracks["n_bins"] == 600


def test_depth_binning_values(tmp_path):
    write_depths(tmp_path / "S.depths", 50, zero_ranges=[(1, 100)])
    tracks = brd.build_coverage_tracks([tmp_path / "S.depths"], brd.Problems())
    assert tracks["bin_size"] == 50 and tracks["n_bins"] == 599
    assert tracks["samples"]["S"][:2] == [0, 0] and tracks["samples"]["S"][2] == 50
    assert len(tracks["samples"]["S"]) == 599


def test_sample_names_with_dots_from_depths(tmp_path):
    assert brd.sample_from_depth_name("a/b/WW_003.rep1.depths") == "WW_003.rep1"


def test_time_series_requires_two_dates():
    samples = [
        {"sample": "a", "sample_type": "sample", "collection_date": "2026-01-01", "site": "X"},
        {"sample": "b", "sample_type": "sample", "collection_date": "2026-01-01", "site": "X"},
    ]
    demix = {n: {"status": "ok", "groups": [["Delta", 1.0]]} for n in "ab"}
    assert brd.build_time_series(samples, demix) is None
    samples[1]["collection_date"] = "2026-01-08"
    ts = brd.build_time_series(samples, demix)
    assert ts["n_dates"] == 2 and ts["sites"][0]["series"][0]["values"] == [1.0, 1.0]


# -- negative controls, `missing`, samplesheet metadata, merged QC table ---------------------------
def _loaded(tmp_path, name, text):
    path = tmp_path / name
    path.write_text(text)
    return path


QC_HEADER = "\t".join(brd.QC_COLUMNS) + "\n"


def _qc_row(sample, sample_type, status, reasons="", mapped="5998"):
    cols = dict.fromkeys(brd.QC_COLUMNS, "NA")
    cols.update(
        sample=sample,
        sample_type=sample_type,
        qc_status=status,
        qc_reasons=reasons,
        mapped_reads=mapped,
    )
    return "\t".join(cols[c] for c in brd.QC_COLUMNS) + "\n"


def test_ok_negative_control_is_not_downgraded_and_is_listed_with_its_lineages(tmp_path):
    qc = _loaded(
        tmp_path,
        "qc.tsv",
        QC_HEADER
        + _qc_row("S1", "sample", "PASS")
        + _qc_row("NEG_C", "negative", "FAIL", "ntc_reads"),
    )
    status = _loaded(
        tmp_path,
        "status.tsv",
        "sample\tstatus\treason\tgenome_coverage_pct\tbarcode_version\tbarcode_version_iso\n"
        "S1\tok\t\t\tv1\t\nNEG_C\tok\t\t\tv1\t\nGONE\tmissing\tempty demix file\t\t\t\n",
    )
    agg_header = "\tsummarized\tlineages\tabundances\tresid\tcoverage\n"
    agg = _loaded(
        tmp_path, "agg.tsv", agg_header + "S1\t[('Omicron', 1.0)]\tBA.2\t1.0\t0.1\t99.0\n"
    )
    controls = _loaded(
        tmp_path, "controls.tsv", agg_header + "NEG_C\t[('Omicron', 1.0)]\tBA.2\t1.0\t0.2\t90.0\n"
    )
    run_info = _loaded(
        tmp_path, "run_info.json", '{"mode": "full", "freyja_data": {"barcode_version": "v1"}}'
    )
    out = tmp_path / "report_data.json"
    merged = tmp_path / "qc_merged.tsv"
    rc = brd.main(
        [
            "--run-info", str(run_info), "--qc", str(qc), "--demix-status", str(status),
            "--aggregated", str(agg), "--controls-aggregated", str(controls),
            "--write-qc-tsv", str(merged), "--out", str(out),
        ]
    )  # fmt: skip
    assert rc == 0
    data = json.loads(out.read_text())
    # `missing` is a valid status (no input_problems flag), the control stays ok
    assert "input_problems" not in {f["code"] for f in data["run_flags"]}
    assert data["demix"]["NEG_C"]["status"] == "ok"
    assert data["demix"]["GONE"]["status"] == "missing"
    (control,) = data["negative_controls"]
    assert control["sample"] == "NEG_C" and control["excluded_from_aggregate"] is True
    assert control["lineages"] == [["BA.2", 1.0]] and control["coverage"] == 90.0
    assert control["qc_reasons"] == ["ntc_reads"]
    assert data["available"]["negative_controls"] is True
    # the published QC table carries the downstream columns (formatting of the rest unchanged)
    rows = {r["sample"]: r for r in brd.read_tsv(merged)}
    assert rows["S1"]["demix_status"] == "ok" and rows["S1"]["barcode_version"] == "v1"
    assert rows["NEG_C"]["qc_reasons"] == "ntc_reads" and rows["NEG_C"]["mapped_reads"] == "5998"
    assert tuple(rows["S1"]) == tuple(brd.QC_COLUMNS)


def test_samplesheet_gives_aggregate_mode_dates_and_a_time_series(tmp_path):
    agg_header = "\tsummarized\tlineages\tabundances\tresid\tcoverage\n"
    rows = "".join(f"S{i}\t[('Omicron', 1.0)]\tBA.2\t1.0\t0.1\t99.0\n" for i in (1, 2, 3))
    agg = _loaded(tmp_path, "agg.tsv", agg_header + rows)
    sheet = _loaded(
        tmp_path,
        "sheet.csv",
        "sample,demix,sample_type,collection_date,site\n"
        "S1,x,sample,2026-01-05,A\nS2,x,sample,2026-01-12,A\nS3,x,negative,,\n",
    )
    run_info = _loaded(tmp_path, "run_info.json", '{"mode": "aggregate"}')
    out = tmp_path / "report_data.json"
    assert (
        brd.main(
            [
                "--run-info",
                str(run_info),
                "--aggregated",
                str(agg),
                "--samplesheet",
                str(sheet),
                "--out",
                str(out),
            ]
        )
        == 0
    )
    data = json.loads(out.read_text())
    meta = {s["sample"]: s for s in data["samples"]}
    assert meta["S1"]["collection_date"] == "2026-01-05" and meta["S1"]["site"] == "A"
    assert meta["S3"]["sample_type"] == "negative"
    assert data["time_series"]["n_dates"] == 2


def test_boot_intervals_and_failed_second_demix_flag(tmp_path):
    agg_header = "\tsummarized\tlineages\tabundances\tresid\tcoverage\n"
    agg = _loaded(
        tmp_path, "agg.tsv", agg_header + "S1\t[('Omicron', 1.0)]\tBA.2\t1.0\t0.1\t99.0\n"
    )
    boot = _loaded(
        tmp_path, "boot.tsv", "sample\tlineage\tabundance\tlo\thi\nS1\tBA.2\t0.99\t0.95\t1.0\n"
    )
    status = _loaded(tmp_path, "cs.tsv", "sample\tstatus\treason\nS1\tfailed\tboom\n")
    run_info = _loaded(tmp_path, "run_info.json", '{"mode": "demix"}')
    out = tmp_path / "report_data.json"
    assert (
        brd.main(
            [
                "--run-info",
                str(run_info),
                "--aggregated",
                str(agg),
                "--boot",
                str(boot),
                "--compare-status",
                str(status),
                "--out",
                str(out),
            ]
        )
        == 0
    )
    data = json.loads(out.read_text())
    assert data["demix"]["S1"]["intervals"] == {"BA.2": [0.95, 1.0]}
    assert data["available"]["boot_intervals"] is True
    assert "barcode_compare_failed" in {f["code"] for f in data["run_flags"]}


def test_aggregate_mode_lists_a_control_known_only_from_the_controls_table(tmp_path):
    agg_header = "\tsummarized\tlineages\tabundances\tresid\tcoverage\n"
    agg = _loaded(
        tmp_path, "agg.tsv", agg_header + "S1\t[('Omicron', 1.0)]\tBA.2\t1.0\t0.1\t99.0\n"
    )
    controls = _loaded(
        tmp_path, "c.tsv", agg_header + "NEG_C\t[('Omicron', 1.0)]\tBA.2\t1.0\t0.2\t90.0\n"
    )
    sheet = _loaded(tmp_path, "sheet.csv", "sample,sample_type\nS1,sample\nNEG_C,negative\n")
    run_info = _loaded(tmp_path, "run_info.json", '{"mode": "aggregate"}')
    out = tmp_path / "report_data.json"
    args = ["--run-info", run_info, "--aggregated", agg, "--controls-aggregated", controls]
    args += ["--samplesheet", sheet, "--out", out]
    assert brd.main([str(a) for a in args]) == 0
    data = json.loads(out.read_text())
    assert data["demix"]["NEG_C"]["status"] == "ok" and data["demix"]["NEG_C"]["groups"] == []
    assert [c["sample"] for c in data["negative_controls"]] == ["NEG_C"]
    assert "demix_failed" not in {f["code"] for f in data["run_flags"]}


def test_depths_only_input_does_not_leak_track_metadata_as_samples(tmp_path):
    depths = tmp_path / "S1.depths"
    depths.write_text("".join(f"NC_045512.2\t{i}\tA\t20\n" for i in range(1, 301)))
    run_info = _loaded(tmp_path, "run_info.json", '{"mode": "demix"}')
    out = tmp_path / "report_data.json"
    assert brd.main(["--run-info", str(run_info), "--depths", str(depths), "--out", str(out)]) == 0
    data = json.loads(out.read_text())
    assert [s["sample"] for s in data["samples"]] == ["S1"]
    assert not data["run_flags"] or "contig" not in json.dumps(data["run_flags"])


def test_summary_params_are_real_pipeline_params():
    schema = json.loads(
        (Path(brd.__file__).resolve().parents[1] / "nextflow_schema.json").read_text()
    )
    known = {name for group in schema["$defs"].values() for name in group["properties"]}
    # Params of the dropout-impact / unexplained / time-series / primer-trimming analyses are named in
    # docs/wiring/*.md and are added to the schema when those analyses are wired into the pipeline;
    # remove them from this list as they land (the test then checks them like the others).
    pending = set()  # every analysis param is in the schema now
    assert set(brd.SUMMARY_PARAMS) <= known | pending, sorted(
        set(brd.SUMMARY_PARAMS) - known - pending
    )
    assert {"freyja_boot", "annotate_variants"} <= set(brd.SUMMARY_PARAMS)


def test_skipped_negative_control_is_listed_as_a_control_not_a_missing_sample(tmp_path):
    agg_header = "\tsummarized\tlineages\tabundances\tresid\tcoverage\n"
    agg = _loaded(
        tmp_path, "agg.tsv", agg_header + "S1\t[('Omicron', 1.0)]\tBA.2\t1.0\t0.1\t99.0\n"
    )
    status = _loaded(
        tmp_path,
        "status.tsv",
        "sample\tstatus\treason\tgenome_coverage_pct\tbarcode_version\tbarcode_version_iso\n"
        "Blank_7\tmissing\tmissing demix (legacy layout)\t\t\t\n"
        "NEG_old\tmissing\tempty demix file\t\t\t\n",
    )
    sheet = _loaded(tmp_path, "usable.csv", "sample,sample_type\nS1,sample\n")
    skipped = _loaded(
        tmp_path,
        "skipped.csv",
        "sample,sample_type,collection_date,site\nBlank_7,negative,2026-01-05,A\nNEG_old,negative,,\n",
    )
    run_info = _loaded(tmp_path, "run_info.json", '{"mode": "aggregate"}')
    out = tmp_path / "report_data.json"
    args = ["--run-info", run_info, "--aggregated", agg, "--demix-status", status]
    args += ["--samplesheet", sheet, skipped, "--out", out]
    assert brd.main([str(a) for a in args]) == 0
    data = json.loads(out.read_text())
    types = {s["sample"]: s["sample_type"] for s in data["samples"]}
    assert types == {"S1": "sample", "Blank_7": "negative", "NEG_old": "negative"}
    assert {c["sample"] for c in data["negative_controls"]} == {"Blank_7", "NEG_old"}
    assert "demix_missing" not in {f["code"] for f in data["run_flags"]}

    # without any sheet (report-mode rebuild) the NEG_/NC_ pattern still identifies the control
    out2 = tmp_path / "report_data2.json"
    args = ["--run-info", run_info, "--aggregated", agg, "--demix-status", status, "--out", out2]
    assert brd.main([str(a) for a in args]) == 0
    types2 = {s["sample"]: s["sample_type"] for s in json.loads(out2.read_text())["samples"]}
    assert types2["NEG_old"] == "negative" and types2["Blank_7"] == "sample"


def test_header_carries_the_library_protocol(tmp_path):
    run_info = _loaded(
        tmp_path,
        "run_info.json",
        '{"mode": "full", "protocol": "capture", "primer_trimming": false, "duplicates": "mark"}',
    )
    out = tmp_path / "report_data.json"
    assert brd.main(["--run-info", str(run_info), "--out", str(out)]) == 0
    header = json.loads(out.read_text())["header"]
    assert (header["protocol"], header["primer_trimming"], header["duplicates"]) == (
        "capture",
        False,
        "mark",
    )
    old = _loaded(tmp_path, "old.json", '{"mode": "full"}')
    assert brd.main(["--run-info", str(old), "--out", str(out)]) == 0
    assert json.loads(out.read_text())["header"]["protocol"] is None  # older runs: nothing recorded


# =============================================================================================
# Primer scheme / trimming columns, dropout impact, unexplained mutations, time series, growth
# =============================================================================================
G = FIX / "g_new_analyses"
G_OPTIONS = {
    "--qc": "qc_summary.tsv",
    "--aggregated": "freyja_aggregated.tsv",
    "--samplesheet": "samplesheet.csv",
    "--dropout-impact": "dropout_impact.tsv",
    "--dropout-summary": "dropout_impact_summary.tsv",
    "--unexplained": "unexplained_mutations.tsv",
    "--timeseries": "freyja_timeseries.tsv",
    "--timeseries-lineages": "freyja_timeseries_lineages.tsv",
    "--growth-rates": "freyja_growth_rates.csv",
    "--growth-inputs": "freyja_growth_inputs.json",
}


UE_HEADER = (
    "sample\tposition\tref\talt\tobserved_freq\texpected_freq\texcess\tdepth\talt_depth"
    "\tgene\taa_change\teffect\tn_samples\trecurrent\tp_value"
)


def build_g(tmp_path, drop=(), extra=None, run_info=None):
    args = ["--run-info", str(run_info or G / "run_info.json")]
    for option, name in G_OPTIONS.items():
        if option not in drop:
            args += [option, str(G / name)]
    args += list(extra or []) + ["--out", str(tmp_path / "g.json")]
    assert brd.main(args) == 0
    return json.loads((tmp_path / "g.json").read_text())


def test_primer_columns_typed_and_old_qc_tables_still_load(tmp_path):
    data = build_g(tmp_path)
    s = {x["sample"]: x for x in data["samples"]}
    assert s["SiteA_W03_1"]["primer_scheme_check"] == "mismatch"
    assert s["SiteA_W03_1"]["primer_scheme_best"] == "articv4.1"
    assert s["SiteA_W03_1"]["primer_trimmed_pct"] == 8.74
    assert s["SiteA_W03_1"]["primer_scheme_enrichment"] == 1.3
    assert s["SiteA_W03_1"]["qc_reasons"] == ["primer_scheme_mismatch"]
    assert s["SiteA_W01_1"]["primer_scheme_check"] == "ok"
    assert s["NEG_1"]["primer_trimmed_pct"] is None and s["NEG_1"]["primer_scheme_best"] is None
    # an older table without the four columns: present as nulls, no problem flag
    old = FIX / "a_full" / "qc_summary.tsv"
    out = tmp_path / "old.json"
    assert (
        brd.main(["--run-info", str(G / "run_info.json"), "--qc", str(old), "--out", str(out)]) == 0
    )
    old_data = json.loads(out.read_text())
    first = old_data["samples"][0]
    assert first["primer_trimmed_pct"] is None and first["primer_scheme_check"] is None
    assert "input_problems" not in flag_codes(old_data)
    assert "primer_scheme_mismatch" not in flag_codes(old_data)


def test_primer_scheme_check_garbage_becomes_null(tmp_path):
    qc = tmp_path / "qc.tsv"
    qc.write_text(
        "sample\tqc_status\tprimer_scheme_check\tprimer_scheme_best\nA\tPASS\tMAYBE\tNA\nB\tPASS\tOK\tNA\n"
    )
    out = tmp_path / "r.json"
    assert (
        brd.main(["--run-info", str(G / "run_info.json"), "--qc", str(qc), "--out", str(out)]) == 0
    )
    got = {s["sample"]: s["primer_scheme_check"] for s in json.loads(out.read_text())["samples"]}
    assert got == {"A": None, "B": "ok"}


def test_primer_scheme_mismatch_flag(tmp_path):
    data = build_g(tmp_path)
    (flag,) = [f for f in data["run_flags"] if f["code"] == "primer_scheme_mismatch"]
    assert flag["level"] == "warn" and flag["samples"] == ["SiteA_W03_1", "SiteB_W02_1"]
    assert "articv4.1" in flag["message"] and "declared scheme midnight" in flag["message"]
    assert "check --primer_scheme" in flag["message"]
    assert "uninformative" not in flag["message"]  # uninformative is never a mismatch


def test_primer_scheme_mismatch_flag_without_declared_scheme_and_with_two_best(tmp_path):
    info = json.loads((G / "run_info.json").read_text())
    del info["params"]["primer_scheme"]
    ri = tmp_path / "ri.json"
    ri.write_text(json.dumps(info))
    qc = tmp_path / "qc.tsv"
    qc.write_text(
        "sample\tqc_status\tprimer_scheme_check\tprimer_scheme_best\n"
        "A\tWARN\tmismatch\tarticv4.1\nB\tWARN\tmismatch\tswift\nC\tPASS\tok\tNA\n"
    )
    out = tmp_path / "r.json"
    assert brd.main(["--run-info", str(ri), "--qc", str(qc), "--out", str(out)]) == 0
    (flag,) = [
        f for f in json.loads(out.read_text())["run_flags"] if f["code"] == "primer_scheme_mismatch"
    ]
    assert "the declared scheme" in flag["message"] and "articv4.1: A" in flag["message"]
    assert "swift: B" in flag["message"]


def test_new_params_and_thresholds_only_when_recorded(tmp_path):
    data = build_g(tmp_path)
    assert data["params_summary"]["unexplained_min_recurrence"] == 2
    assert data["params_summary"]["freyja_growth_serial_interval"] == 5.5
    assert data["thresholds"]["unexplained_max_pvalue"] == pytest.approx(1.67e-6)
    assert data["thresholds"]["dropout_min_abundance"] == 0.01
    # a run that recorded none of them: no new threshold keys at all
    plain = run_builder("a_full", tmp_path / "plain.json")
    assert not any(k.startswith(("unexplained", "dropout", "qc_min")) for k in plain["thresholds"])
    assert not any(
        k.startswith(("unexplained", "dropout", "freyja_growth")) for k in plain["params_summary"]
    )


def test_default_build_has_no_new_sections(tmp_path):
    data = run_builder("a_full", tmp_path / "r.json")
    for key in (
        "dropout_impact",
        "unexplained_mutations",
        "timeseries",
        "growth_rates",
        "dashboard",
    ):
        assert data[key] is None, key
    for key in ("dropout_impact", "unexplained_mutations", "timeseries", "growth_rates"):
        assert data["available"][key] is False, key
    assert "dropout_pct_informative_lost" not in data["samples"][0]
    assert "n_unexplained" not in data["samples"][0]


# --------------------------------------------------------------------------- dropout impact
def test_dropout_impact_structure(tmp_path):
    data = build_g(tmp_path)
    impact = data["dropout_impact"]
    assert data["available"]["dropout_impact"] is True
    assert len(impact) == 32  # every sample of the summary appears
    a = impact["SiteA_W03_1"]
    assert a["summary"] == {
        "n_dropped_amplicons": 2,
        "informative_sites_total": 8006,
        "informative_sites_in_dropped": 850,
        "pct_informative_sites_lost": 10.62,
        "n_dropped_affecting_detected_lineages": 1,
        "defining_sites_total": None,
        "defining_sites_lost": None,
        "pct_defining_sites_lost": None,
    }
    assert [x["amplicon"] for x in a["amplicons"]] == ["nCoV-2019_7", "nCoV-2019_22"]  # by start
    hit = a["amplicons"][1]
    assert hit["affects_detected_lineage"] is True and hit["start"] == 21533 and hit["end"] == 22612
    assert hit["n_sites_in_detected_lineages"] == 10 and hit["mean_depth"] == 3.1
    assert hit["affected_lineages"] == [
        {"lineage": "XBB.1.5", "n_sites": 9},
        {"lineage": "BA.2", "n_sites": 3},
        {"lineage": "BA.5", "n_sites": 3},
    ]
    assert a["amplicons"][0]["affected_lineages"] == []  # NA
    # NA = sample without a lineage estimate
    low = impact["LOWCOV_1"]
    assert low["amplicons"][0]["n_sites_in_detected_lineages"] is None
    assert low["summary"]["n_dropped_affecting_detected_lineages"] is None
    assert impact["SiteA_W01_1"]["amplicons"] == []
    assert impact["SiteA_W01_1"]["summary"]["pct_informative_sites_lost"] == 0.0


def test_dropout_flag_and_per_sample_fields(tmp_path):
    data = build_g(tmp_path)
    (flag,) = [f for f in data["run_flags"] if f["code"] == "dropout_affects_lineages"]
    assert flag["level"] == "warn" and flag["samples"] == ["SiteA_W03_1", "SiteB_W05_1"]
    assert "2 sample(s)" in flag["message"] and "may be biased" in flag["message"]
    s = {x["sample"]: x for x in data["samples"]}
    assert s["SiteA_W03_1"]["dropout_affects_lineages"] is True
    assert s["SiteA_W03_1"]["dropout_pct_informative_lost"] == 10.62
    assert s["SiteB_W02_1"]["dropout_affects_lineages"] is False
    assert s["LOWCOV_1"]["dropout_affects_lineages"] is None  # no lineage estimate
    assert s["SiteC_W01_1"]["dropout_affects_lineages"] is False


def test_dropout_impact_without_summary_and_without_effect(tmp_path):
    data = build_g(tmp_path, drop=("--dropout-summary",))
    assert data["dropout_impact"]["SiteA_W03_1"]["summary"] is None
    assert data["samples"][0]["dropout_affects_lineages"] in (None, True, False)
    # header-only impact + all-zero summary: no flag, still available
    imp = tmp_path / "imp.tsv"
    imp.write_text(
        "sample\tamplicon\tstart\tend\tmean_depth\tn_informative_sites\t"
        "n_sites_in_detected_lineages\taffects_detected_lineage\taffected_lineages\n"
    )
    summ = tmp_path / "summ.tsv"
    summ.write_text(
        "sample\tn_dropped_amplicons\tinformative_sites_total\tinformative_sites_in_dropped\t"
        "pct_informative_sites_lost\tn_dropped_affecting_detected_lineages\n"
        "SiteA_W01_1\t0\t8006\t0\t0.00\t0\n"
    )
    out = tmp_path / "r.json"
    args = ["--run-info", str(G / "run_info.json"), "--qc", str(G / "qc_summary.tsv"),
            "--dropout-impact", str(imp), "--dropout-summary", str(summ), "--out", str(out)]  # fmt: skip
    assert brd.main(args) == 0
    data = json.loads(out.read_text())
    assert data["available"]["dropout_impact"] is True
    assert "dropout_affects_lineages" not in flag_codes(data)


def test_dropout_rows_for_unknown_samples_are_ignored_with_a_note(tmp_path):
    imp = tmp_path / "imp.tsv"
    imp.write_text(
        "sample\tamplicon\tstart\tend\tmean_depth\tn_informative_sites\t"
        "n_sites_in_detected_lineages\taffects_detected_lineage\taffected_lineages\n"
        "GHOST\ta_1\t1\t100\t0.0\t3\t1\ttrue\tBA.2:1\n"
    )
    out = tmp_path / "r.json"
    assert brd.main(["--run-info", str(G / "run_info.json"), "--qc", str(G / "qc_summary.tsv"),
                     "--dropout-impact", str(imp), "--out", str(out)]) == 0  # fmt: skip
    data = json.loads(out.read_text())
    assert "GHOST" not in data["dropout_impact"]
    assert "dropout_affects_lineages" not in flag_codes(data)
    assert any("GHOST" in f["message"] for f in data["run_flags"] if f["code"] == "input_note")


def test_real_dropout_fixture_parses():
    base = ROOT / "tests" / "fixtures" / "dropout_impact"
    out = brd.load_dropout_impact(
        base / "expected_dropout_impact.tsv",
        base / "expected_dropout_impact_summary.tsv",
        brd.Problems(),
    )
    assert {"S1", "S2", "S3", "S4"} <= set(out)
    assert out["S1"]["amplicons"][0]["affected_lineages"] == [
        {"lineage": "A", "n_sites": 2},
        {"lineage": "B", "n_sites": 2},
    ]
    assert out["S3"]["summary"]["n_dropped_affecting_detected_lineages"] is None
    assert out["S3"]["amplicons"][0]["n_sites_in_detected_lineages"] is None
    assert out["S4"]["amplicons"][0]["affects_detected_lineage"] is False


# --------------------------------------------------------------------------- unexplained mutations
def test_unexplained_rows_and_flag(tmp_path):
    data = build_g(tmp_path)
    rows = data["unexplained_mutations"]
    assert data["available"]["unexplained_mutations"] is True
    assert len(rows) == 5 and rows[0]["sample"] == "SiteC_W05_1"  # file order kept
    first = rows[0]
    assert first["position"] == 22102 and first["alt"] == "T" and first["aa_change"] == "E180D"
    assert first["observed_freq"] == 0.3 and first["expected_freq"] == 0.0 and first["depth"] == 520
    assert first["recurrent"] is True and first["n_samples"] == 3
    assert first["p_value"] == pytest.approx(2.39e-42) and isinstance(first["p_value"], float)
    assert (
        rows[3]["recurrent"] is False and rows[4]["gene"] is None and rows[4]["aa_change"] is None
    )
    (flag,) = [f for f in data["run_flags"] if f["code"] == "unexplained_mutations"]
    assert flag["level"] == "info" and flag["samples"] == [
        "SiteA_W04_1",
        "SiteB_W03_1",
        "SiteC_W05_1",
    ]
    assert "S:E180D in 3 samples" in flag["message"] and "1 mutation(s)" in flag["message"]
    s = {x["sample"]: x for x in data["samples"]}
    assert s["SiteC_W05_1"]["n_unexplained"] == 1 and s["SiteA_W01_1"]["n_unexplained"] == 0
    assert data["thresholds"]["unexplained_min_depth"] == 20


def test_unexplained_empty_table_and_no_recurrence(tmp_path):
    header = UE_HEADER
    empty = tmp_path / "ue.tsv"
    empty.write_text(header + "\n")
    data = build_g(tmp_path, drop=("--unexplained",), extra=["--unexplained", str(empty)])
    assert (
        data["unexplained_mutations"] == [] and data["available"]["unexplained_mutations"] is True
    )
    assert "unexplained_mutations" not in flag_codes(data)
    single = tmp_path / "ue1.tsv"
    single.write_text(
        header
        + "\nSiteA_W01_1\t28802\tA\tT\t0.2\t0\t0.2\t300\t60\tN\tR177*\tstop_gained\t1\tfalse\t3.4e-21\n"
    )
    data = build_g(tmp_path, drop=("--unexplained",), extra=["--unexplained", str(single)])
    assert len(data["unexplained_mutations"]) == 1 and "unexplained_mutations" not in flag_codes(
        data
    )


def test_unexplained_accepts_the_per_sample_table():
    base = ROOT / "tests" / "fixtures" / "unexplained" / "concat"
    rows = brd.load_unexplained(base / "S_TWO.unexplained_mutations.tsv", brd.Problems())
    assert [r["position"] for r in rows] == [22102, 23063]
    assert rows[0]["n_samples"] is None and rows[0]["recurrent"] is False
    assert rows[1]["aa_change"] == "N501Y" and rows[1]["expected_freq"] == 0.7
    assert brd.load_unexplained(base / "S_SKIPPED.unexplained_mutations.tsv", brd.Problems()) == []


def test_unexplained_bad_values_are_reported_not_fatal(tmp_path):
    ue = tmp_path / "ue.tsv"
    ue.write_text("sample\tposition\tobserved_freq\nA\tabc\t0.2\nB\t5\t0.3\n")
    problems = brd.Problems()
    rows = brd.load_unexplained(ue, problems)
    assert [r["sample"] for r in rows] == ["B"] and "unexplained_mutations.tsv" in problems.items


# --------------------------------------------------------------------------- time series
def test_timeseries_table_loads_real_timeseries_output():
    problems = brd.Problems()
    ts = brd.load_timeseries(
        G / "freyja_timeseries.tsv", G / "freyja_timeseries_lineages.tsv", problems
    )
    assert not problems.items
    assert ts["enabled"] is True and ts["has_viral_load"] is True
    assert ts["sites"] == ["SiteA", "SiteB", "SiteC"] and len(ts["dates"]) == 8
    row = ts["summarized"][0]
    assert row == {
        "site": "SiteA", "collection_date": "2026-01-05", "lineage_group": "Omicron",
        "abundance": 0.80063, "abundance_load_weighted": 0.80063, "n_samples": 1,
        "mean_viral_load": 175675.0, "total_viral_load": 175675.0,
    }  # fmt: skip
    assert len(ts["summarized"]) == 72 and len(ts["lineages"]) == 85
    assert ts["lineages_total_groups"] is None


def test_timeseries_in_report_data_and_available(tmp_path):
    data = build_g(tmp_path)
    ts = data["timeseries"]
    assert data["available"]["timeseries"] is True and data["available"]["time_series"] is True
    assert ts["min_coverage"] == 60 and ts["has_viral_load"] is True
    # today's derived panel is still computed (fallback / unchanged)
    assert data["time_series"]["n_dates"] >= 2


def test_timeseries_without_viral_load_and_with_nas(tmp_path):
    tsv = tmp_path / "ts.tsv"
    tsv.write_text(
        "site\tcollection_date\tlineage_group\tabundance\tabundance_load_weighted\tn_samples\tmean_viral_load\ttotal_viral_load\n"
        "NA\t2026-01-05\tDelta\t0.5\tNA\t2\tNA\tNA\nNA\t2026-01-12\tDelta\t0.7\tNA\t1\tNA\tNA\n"
    )
    ts = brd.load_timeseries(tsv, None, brd.Problems())
    assert (
        ts["has_viral_load"] is False and ts["sites"] == [] and ts["summarized"][0]["site"] is None
    )
    assert ts["summarized"][0]["abundance_load_weighted"] is None and ts["lineages"] == []


def test_timeseries_lineage_table_is_capped(tmp_path):
    lines = [
        "site\tcollection_date\tlineage_group\tabundance\tabundance_load_weighted\tn_samples\tmean_viral_load\ttotal_viral_load"
    ]
    for d in ("2026-01-05", "2026-01-12"):
        lines += [f"A\t{d}\tL{i}\t{(100 - i) / 1000}\tNA\t1\tNA\tNA" for i in range(80)]
    tsv = tmp_path / "tl.tsv"
    tsv.write_text("\n".join(lines) + "\n")
    ts = brd.load_timeseries(None, tsv, brd.Problems())
    assert ts["lineages_total_groups"] == 80
    assert len({r["lineage_group"] for r in ts["lineages"]}) == brd.MAX_TIMESERIES_LINEAGES
    assert ts["lineages"][0]["lineage_group"] == "L0"


# --------------------------------------------------------------------------- growth rates
def test_growth_rates_real_format(tmp_path):
    data = build_g(tmp_path)
    growth = data["growth_rates"]
    assert data["available"]["growth_rates"] is True and growth["status"] == "ok"
    rows = {r["lineage"]: r for r in growth["rows"]}
    assert rows["XBB.1.5"] == {
        "lineage": "XBB.1.5",
        "advantage_pct": 60.3,
        "interval_text": "[44.74% , 49.89%]",
    }
    assert (
        rows["BA.2"]["advantage_pct"] == -35.1
        and rows["BA.2"]["interval_text"] == "[-31.42% , -27.98%]"
    )
    assert growth["inputs"]["n_samples"] == 30 and growth["inputs"]["n_dates"] == 8
    assert growth["params"] == {
        "days": 56,
        "nboots": 1000,
        "seed": 1,
        "serial_interval": 5.5,
        "scale_by_viral_load": False,
    }
    assert growth["reason"] is None


def test_growth_failed_text_is_kept(tmp_path):
    failed = tmp_path / "freyja_growth_rates.failed.txt"
    failed.write_text(
        "freyja relgrowthrate produced no growth-rate table: insufficient data: 3 distinct collection "
        "date(s) spanning 14 day(s) within the last 56 days; need at least 4 dates spanning at least 14 days\n"
    )
    data = build_g(
        tmp_path, drop=("--growth-rates", "--growth-inputs"), extra=["--growth-failed", str(failed)]
    )
    growth = data["growth_rates"]
    assert growth["status"] == "failed" and growth["rows"] == []
    assert "insufficient data: 3 distinct collection date(s)" in growth["reason"]
    assert data["available"]["growth_rates"] is True


def test_growth_absent_and_dashboard_name(tmp_path):
    data = build_g(
        tmp_path,
        drop=("--growth-rates", "--growth-inputs"),
        extra=["--dashboard-name", "freyja_dashboard.html"],
    )
    assert data["growth_rates"] is None and data["dashboard"] == {"file": "freyja_dashboard.html"}
    assert build_g(tmp_path)["dashboard"] is None


def test_growth_csv_variants(tmp_path):
    csv_path = tmp_path / "g.csv"
    csv_path.write_text(
        'Lineage,Estimated Advantage,Bootstrap 95% interval\nA,+73.3%,"[1% , 2%]"\nB,nan,\nC,,\n'
    )
    problems = brd.Problems()
    growth = brd.load_growth_rates(csv_path, None, None, problems)
    rows = {r["lineage"]: r["advantage_pct"] for r in growth["rows"]}
    assert rows == {"A": 73.3, "B": None, "C": None}
    # a failure note alongside a table: the table wins
    note = tmp_path / "f.txt"
    note.write_text("something\n")
    both = brd.load_growth_rates(csv_path, note, None, problems)
    assert both["status"] == "ok" and both["reason"] is None


def test_new_cli_options_ignore_no_file_placeholders(tmp_path):
    (tmp_path / "NO_FILE").write_text("")
    out = tmp_path / "r.json"
    args = [
        "--run-info",
        str(G / "run_info.json"),
        "--qc",
        str(G / "qc_summary.tsv"),
        "--out",
        str(out),
    ]
    for option in ("--dropout-impact", "--dropout-summary", "--unexplained", "--timeseries",
                   "--timeseries-lineages", "--growth-rates", "--growth-failed", "--growth-inputs"):  # fmt: skip
        args += [option, str(tmp_path / "NO_FILE")]
    args += ["--dashboard-name", ""]
    assert brd.main(args) == 0
    data = json.loads(out.read_text())
    for key in (
        "dropout_impact",
        "unexplained_mutations",
        "timeseries",
        "growth_rates",
        "dashboard",
    ):
        assert data[key] is None


# --------------------------------------------------------------------------- consensus reason / flags
def test_consensus_reason_column_flags_and_samples(tmp_path):
    data = build_g(tmp_path, extra=["--consensus", str(G / "consensus_stats.tsv")])
    rows = {r["sample"]: r for r in data["consensus"]}
    assert rows["LOWCOV_1"]["status"] == "failed"
    assert rows["LOWCOV_1"]["reason"].startswith(
        "bcftools consensus failed: The fasta sequence does not match"
    )
    assert (
        rows["SiteA_W03_1"]["reason"]
        == "skipped variant at NC_045512.2:22942 (overlaps another variant)"
    )
    assert rows["SiteA_W01_1"]["reason"] is None and rows["NEG_1"]["reason"] is None
    s = {x["sample"]: x for x in data["samples"]}
    assert (
        s["LOWCOV_1"]["consensus_status"] == "failed"
        and "11287" in s["LOWCOV_1"]["consensus_reason"]
    )
    assert (
        s["SiteA_W01_1"]["consensus_reason"] is None
        and s["SiteA_W01_1"]["consensus_status"] == "ok"
    )
    (failed,) = [f for f in data["run_flags"] if f["code"] == "consensus_failed"]
    assert failed["level"] == "warn" and failed["samples"] == ["LOWCOV_1"]
    assert failed["message"].startswith("consensus genome failed for 1 sample(s): LOWCOV_1;")
    assert "all-N: see the consensus table for the reason" in failed["message"]
    (warned,) = [f for f in data["run_flags"] if f["code"] == "consensus_warnings"]
    assert warned["level"] == "info" and warned["samples"] == ["SiteA_W03_1"]
    assert warned["message"].startswith(
        "consensus built with warnings for 1 sample(s): SiteA_W03_1"
    )


def test_old_six_column_consensus_table_still_loads(tmp_path):
    data = build_g(tmp_path, extra=["--consensus", str(G / "consensus_stats_6col.tsv")])
    assert all(r["reason"] is None for r in data["consensus"])
    assert len(data["consensus"]) == 5
    codes = flag_codes(data)
    assert "consensus_warnings" not in codes and "input_problems" not in codes
    assert "consensus_failed" in codes  # status failed still flagged, just without a reason


def test_no_consensus_flags_for_a_clean_consensus(tmp_path):
    data = run_builder("a_full", tmp_path / "r.json")
    assert not {"consensus_failed", "consensus_warnings"} & flag_codes(data)


# =============================================================================================
# Lineage groups (lineage_summary.csv) and site aggregates
# =============================================================================================
H = FIX / "h_lineage_groups"


def build_h(tmp_path, drop=(), extra=None):
    options = {
        "--aggregated": "freyja_aggregated.tsv",
        "--demix-status": "demix_status.tsv",
        "--samplesheet": "samplesheet.csv",
        "--lineage-summary": "lineage_summary.csv",
        "--lineage-groups-resolved": "lineage_groups_resolved.tsv",
        "--site-aggregates": "site_aggregates_freyja.tsv",
        "--site-aggregates-meta": "site_aggregates_meta.tsv",
    }
    args = ["--run-info", str(H / "run_info.json")]
    for option, name in options.items():
        if option not in drop:
            args += [option, str(H / name)]
    args += list(extra or []) + ["--out", str(tmp_path / "h.json")]
    assert brd.main(args) == 0
    return json.loads((tmp_path / "h.json").read_text())


def test_site_aggregates_structure(tmp_path):
    data = build_h(tmp_path)
    sa = data["site_aggregates"]
    assert [a["row_label"] for a in sa][:2] == ["Catchment_X_2026-W02", "Catchment_X_2026-W03"]
    w3 = sa[1]
    assert w3["aggregate_name"] == "Catchment_X" and w3["period"] == "2026-W03"
    assert w3["period_start"] == "2026-01-12" and w3["n_sites_included"] == 1
    assert w3["sites_missing"] == ["Site_2"] and w3["missing_reasons"] == {
        "Site_2": "not_ok:failed"
    }
    assert (
        w3["weight_basis"] == "equal (viral_load missing)" and w3["weight_fraction_included"] == 0.5
    )
    assert (
        w3["lineages"] == [["NB.1.8", 1.0]]
        and w3["groups"] == [["Omicron", 1.0]]
        and w3["n_samples"] == 1
    )
    assert sa[3]["aggregate_name"] == "Coastal" and sa[3]["missing_reasons"] == {
        "Site_3": "no_sample"
    }
    assert sa[0]["weight_fraction_included"] == 1.0 and sa[0]["sites_missing"] == []


def test_site_aggregates_old_meta_without_reasons_and_without_meta(tmp_path):
    old = ROOT / "tests" / "fixtures" / "lineage_groups"
    problems = brd.Problems()
    sa = brd.load_site_aggregates(
        old / "site_aggregates_freyja.tsv", old / "site_aggregates_meta.tsv", problems
    )
    assert [a["row_label"] for a in sa] == [
        "Catchment_X_2026-W02", "Catchment_X_2026-W03", "Catchment_X_2026-W04",
    ]  # fmt: skip
    assert sa[0]["missing_reasons"] == {} and sa[0]["weight_fraction_included"] is None
    assert sa[2]["lineages"] == [] and sa[2]["n_sites_included"] == 0  # no solution, kept
    assert not problems.items
    # only the Freyja table: name and period are parsed from the row label
    only = brd.load_site_aggregates(old / "site_aggregates_freyja.tsv", None, brd.Problems())
    assert only[0]["aggregate_name"] == "Catchment_X" and only[0]["period"] == "2026-W02"
    assert only[0]["lineages"][0][0] == "JN.1.16" and only[0]["n_sites_included"] is None
    # only the meta table: entries without a solution
    meta = brd.load_site_aggregates(None, old / "site_aggregates_meta.tsv", brd.Problems())
    assert len(meta) == 3 and all(a["lineages"] == [] for a in meta)


def test_site_aggregate_flag_also_covers_a_weaker_weighting(tmp_path):
    meta = tmp_path / "meta.tsv"
    meta.write_text(
        "row_label\taggregate_name\tperiod\tperiod_start\tn_sites_included\tsites_included\tsites_missing\t"
        "missing_reasons\tweight_basis\tweight_fraction_included\tn_samples\n"
        "A_2026-W01\tA\t2026-W01\t2025-12-29\t2\tS1;S2\t\t\tsite_weight (viral_load missing)\t1.0\t2\n"
        "A_2026-W02\tA\t2026-W02\t2026-01-05\t2\tS1;S2\t\t\tviral_load\t1.0\t2\n"
    )
    out = tmp_path / "r.json"
    assert (
        brd.main(
            [
                "--run-info",
                str(H / "run_info.json"),
                "--site-aggregates-meta",
                str(meta),
                "--out",
                str(out),
            ]
        )
        == 0
    )
    flag = next(
        f
        for f in json.loads(out.read_text())["run_flags"]
        if f["code"] == "site_aggregate_incomplete"
    )
    assert "1 of 2" in flag["message"] and "weaker weighting" in flag["message"]


def test_default_build_has_no_lineage_group_data(tmp_path):
    data = run_builder("a_full", tmp_path / "r.json")
    assert data["lineage_groups"] is None and data["site_aggregates"] is None
    assert (
        data["available"]["lineage_groups"] is False
        and data["available"]["site_aggregates"] is False
    )
    assert not {
        "lineage_groups_warnings",
        "lineage_groups_ungrouped",
        "site_aggregate_incomplete",
    } & flag_codes(data)


# --------------------------------------------------------------------------- lineage groups (provenance + flags only)
def test_real_lineage_summary_is_read_for_ungrouped_lineages_only():
    rows = brd.load_lineage_ungrouped(H / "real_lineage_summary.csv")
    assert len(rows) == 9 and rows[0]["sample"] == "WW-1274E" and rows[0]["status"] == "ok"
    assert set(rows[0]) == {
        "sample",
        "status",
        "other",
        "ungrouped",
    }  # no group columns, no lineages
    assert rows[0]["other"] is None and rows[0]["ungrouped"] == []
    hit = next(r for r in rows if r["ungrouped"])
    assert hit["ungrouped"] == [["B.52", 0.001032]] and hit["other"] == 0.001032
    assert any(r["status"] != "ok" and r["other"] is None for r in rows)


def test_lineage_summary_other_columns_and_delimiters_are_tolerated(tmp_path):
    path = tmp_path / "s.csv"
    path.write_text(
        "sample;site;status;G1;G2;other;ungrouped_lineages\nS1;A;ok;0.5;;0.5;X:0.3 Y:0.2\nS2;A;failed;;;;\n"
    )
    rows = brd.load_lineage_ungrouped(path)
    assert rows[0]["ungrouped"] == [["X", 0.3], ["Y", 0.2]] and rows[1]["other"] is None
    no_other = tmp_path / "n.csv"
    no_other.write_text("sample,status\nS1,ok\n")
    assert brd.load_lineage_ungrouped(no_other)[0]["ungrouped"] == []


def test_lineage_groups_object_is_small_and_has_no_rows(tmp_path):
    data = build_h(tmp_path)
    lg = data["lineage_groups"]
    assert data["available"]["lineage_groups"] is True
    assert set(lg) == {
        "n_groups", "groups_sha256", "dataset_version", "file", "warnings", "n_samples",
        "n_samples_with_ungrouped", "top_ungrouped", "samples_above_threshold", "ungrouped_flag_threshold",
    }  # fmt: skip
    assert (
        lg["n_groups"] == 8
        and lg["dataset_version"] == "01_02_2026-03-04"
        and len(lg["groups_sha256"]) == 64
    )
    assert lg["file"] == "groups.csv"  # basename of run_info.params.lineage_groups
    assert lg["n_samples"] == 3 and lg["n_samples_with_ungrouped"] == 1
    assert lg["top_ungrouped"] == [["A.1", 0.03, 1]] and lg["samples_above_threshold"] == [
        "Site_1_week_2"
    ]
    assert (
        lg["warnings"][0]["group"] == "KP.3*"
        and "overlaps in_total group" in lg["warnings"][0]["text"]
    )
    assert len(lg["warnings"]) == 3 and lg["ungrouped_flag_threshold"] == 0.02
    raw = json.dumps(data)
    assert "JN.1.16" not in json.dumps(lg) and "Catchment" not in json.dumps(lg)
    assert "major" not in lg and "rows" not in lg and "site_rows" not in lg
    assert len(raw) < 20000


def test_lineage_group_flags(tmp_path):
    data = build_h(tmp_path)
    codes = {f["code"]: f for f in data["run_flags"]}
    warn = codes["lineage_groups_warnings"]
    assert warn["level"] == "warn" and "3 warning(s) in the lineage groups audit" in warn["message"]
    assert "KP.3*: in_total group overlaps" in warn["message"] and "XAA.1.7.9*" in warn["message"]
    ungrouped = codes["lineage_groups_ungrouped"]
    assert ungrouped["level"] == "info" and ungrouped["samples"] == ["Site_1_week_2"]
    assert "A.1 (3.0%)" in ungrouped["message"]
    assert (
        "add them to your lineage groups file; see summary/lineage_ungrouped_summary.tsv"
        in ungrouped["message"]
    )
    assert "site_aggregate_incomplete" in codes


def test_ungrouped_flag_threshold_option(tmp_path):
    assert "lineage_groups_ungrouped" not in flag_codes(
        build_h(tmp_path, extra=["--ungrouped-flag-threshold", "0.5"])
    )
    low = build_h(tmp_path, extra=["--ungrouped-flag-threshold", "0.001"])
    assert "lineage_groups_ungrouped" in flag_codes(low)
    assert low["lineage_groups"]["ungrouped_flag_threshold"] == 0.001


def test_clean_run_has_no_lineage_group_flags_and_default_file_name(tmp_path):
    out = tmp_path / "r.json"
    info = json.loads((H / "run_info.json").read_text())
    del info["params"]["lineage_groups"]
    ri = tmp_path / "ri.json"
    ri.write_text(json.dumps(info))
    args = ["--run-info", str(ri), "--lineage-summary", str(H / "real_lineage_summary.csv"),
            "--lineage-groups-resolved", str(H / "real_lineage_groups_resolved.tsv"), "--out", str(out)]  # fmt: skip
    assert brd.main(args) == 0
    data = json.loads(out.read_text())
    assert not {"lineage_groups_warnings", "lineage_groups_ungrouped"} & flag_codes(data)
    assert (
        data["lineage_groups"]["file"] == "bundled default"
        and data["lineage_groups"]["warnings"] == []
    )
    assert data["lineage_groups"]["top_ungrouped"] == [
        ["B.52", 0.001032, 1]
    ]  # below the threshold: no flag


def test_removed_site_summary_option_is_rejected():
    with pytest.raises(SystemExit):
        brd.parse_args(["--run-info", "x", "--out", "y", "--lineage-summary-sites", "z"])


def test_lineage_group_options_ignore_no_file_placeholders(tmp_path):
    (tmp_path / "NO_FILE").write_text("")
    out = tmp_path / "r.json"
    args = ["--run-info", str(H / "run_info.json"), "--out", str(out)]
    for option in (
        "--lineage-summary",
        "--lineage-groups-resolved",
        "--site-aggregates",
        "--site-aggregates-meta",
    ):
        args += [option, str(tmp_path / "NO_FILE")]
    assert brd.main(args) == 0
    data = json.loads(out.read_text())
    assert data["lineage_groups"] is None and data["site_aggregates"] is None


def test_lineage_group_params_and_thresholds(tmp_path):
    data = build_h(tmp_path)
    assert data["params_summary"]["lineage_summary_zeros"] == "blank"
    assert data["params_summary"]["site_aggregation_weights"] == "auto"
    assert (
        data["thresholds"]["lineage_major_threshold"] == 0.1
        and data["thresholds"]["lineage_min_abundance"] == 0.001
    )


# =============================================================================================
# Anonymisation for sharing (--anonymize none|paths|full)
# =============================================================================================
A = FIX / "i_anonymize"
A_OPTIONS = {
    "--qc": "qc_summary.tsv",
    "--demix-status": "demix_status.tsv",
    "--aggregated": "freyja_aggregated.tsv",
    "--consensus": "consensus_stats.tsv",
    "--samplesheet": "samplesheet.csv",
    "--site-aggregates": "site_aggregates_freyja.tsv",
    "--site-aggregates-meta": "site_aggregates_meta.tsv",
}
A_IDS = ["NTC_Central_01", "WW_Central_001", "WW_Central_002", "WW_North_001", "WW_North_002.rep1"]
A_SITES = ["Central_WWTP", "North_Plant"]
A_LEAKS = [
    "jdoe",
    "j.doe",
    "example.org",
    "/Users",
    "/home",
    "C:\\",
    "Run 1 Plate A",
    "My Results",
    "~/",
]


def build_a(mode, extra_files=None, run_info=None, map_out=None, out=None):
    args = ["--run-info", str(run_info or A / "run_info.json")]
    for option, name in (extra_files or A_OPTIONS).items():
        args += [option, str(A / name)]
    args += ["--anonymize", mode]
    if map_out:
        args += ["--anonymize-map-out", str(map_out)]
    args += ["--out", str(out)]
    assert brd.main(args) == 0
    return Path(out).read_text()


def numbers_only(obj):
    """The same structure with every string replaced by a marker: numbers, bools, nulls and shape only."""
    if isinstance(obj, dict):
        return {"#keys": len(obj), "#vals": [numbers_only(v) for v in obj.values()]}
    if isinstance(obj, list):
        return [numbers_only(v) for v in obj]
    return "s" if isinstance(obj, str) else obj


def test_default_is_unchanged_and_has_no_anonymized_marker(tmp_path):
    plain = json.loads(build_a("none", out=tmp_path / "n.json"))
    assert "anonymized" not in plain
    raw = (tmp_path / "n.json").read_text()
    assert "/Users/jdoe" in raw and "Central_WWTP" in raw  # nothing was touched
    # the same bytes as a build that never heard of the option
    args = ["--run-info", str(A / "run_info.json")]
    for option, name in A_OPTIONS.items():
        args += [option, str(A / name)]
    out2 = tmp_path / "n2.json"
    assert brd.main(args + ["--out", str(out2)]) == 0
    a, b = json.loads(raw), json.loads(out2.read_text())
    a.pop("generated"), b.pop("generated")
    assert a == b


def test_paths_mode_removes_paths_users_and_emails_only(tmp_path):
    raw = build_a("paths", out=tmp_path / "p.json")
    data = json.loads(raw)
    for leak in A_LEAKS:
        assert leak not in raw, leak
    assert data["anonymized"] == "paths"
    prov = data["provenance"]
    assert prov["command_line"].startswith(
        "nextflow run .../main.nf -profile docker --input '.../samplesheet.csv'"
    )
    assert (
        "--outdir '.../wasteline_out'" in prov["command_line"]
        and "--barcodes container -resume" in prov["command_line"]
    )
    assert (
        prov["params"]["input"] == ".../samplesheet.csv"
        and prov["params"]["outdir"] == ".../wasteline_out"
    )
    assert (
        prov["params"]["genome"] == ".../NC_045512.2.fasta" and prov["params"]["db_dir"] == ".../db"
    )
    assert prov["params"]["email"] == "<email>" and prov["run_name"] == "<user>_run"
    assert prov["freyja_data"]["dataset_path"] == ".../2026-10-06T05-25"
    assert [p["path"] for p in prov["previous_runs"]] == [".../run1", ".../run0"]
    assert data["params_summary"]["primer_bed"] == ".../midnight.bed"
    # not paths: container names, URLs, enums, hashes, versions survive
    assert prov["containers"]["FREYJA_DEMIX"] == "quay.io/biocontainers/freyja:2.0.6--pyhdfd78af_0"
    assert (
        prov["params"]["barcodes"] == "container"
        and prov["freyja_data"]["barcode_version"] == "10_06_2026-05-25"
    )
    assert prov["input_checksums"]["samplesheet"]["sha256"] == "f" * 64
    # free text in the data: reasons lose their paths and the user name, IDs stay in `paths`
    assert data["demix"]["WW_North_002.rep1"]["reason"] == (
        "freyja demix failed for WW_North_002.rep1 in .../WW_North_002.rep1.variants.tsv (<user>)"
    )
    cons = {r["sample"]: r for r in data["consensus"]}
    assert cons["WW_Central_002"]["reason"].endswith(
        "at NC_045512.2:11287: (.../WW_Central_002.consensus.fa)"
    )
    assert (
        cons["WW_North_001"]["reason"]
        == "skipped variant at NC_045512.2:22942 (overlaps another variant) in WW_North_001"
    )
    for sid in A_IDS:
        assert sid in raw
    assert "Central_WWTP" in raw and "Metro_Region" in raw


def test_paths_mode_changes_no_number_and_no_structure(tmp_path):
    plain = json.loads(build_a("none", out=tmp_path / "n.json"))
    paths = json.loads(build_a("paths", out=tmp_path / "p.json"))
    full = json.loads(build_a("full", out=tmp_path / "f.json"))
    for data in (paths, full):
        data.pop("anonymized")
    for other in (paths, full):
        assert numbers_only(plain) == numbers_only(other)
        assert (
            plain["available"] == other["available"] and plain["thresholds"] == other["thresholds"]
        )
        assert [f["code"] for f in plain["run_flags"]] == [f["code"] for f in other["run_flags"]]
        assert (
            plain["demix"]["WW_Central_001"]["lineages"]
            == other["demix"]["S002" if other is full else "WW_Central_001"]["lineages"]
        )


def test_full_mode_pseudonymises_everywhere_and_leaves_lineages_alone(tmp_path):
    raw = build_a("full", out=tmp_path / "f.json", map_out=tmp_path / "map.csv")
    data = json.loads(raw)
    for leak in A_LEAKS + A_IDS + A_SITES + ["Metro_Region", "Central", "North_", "WW_", "NTC_"]:
        assert leak not in raw, leak
    assert data["anonymized"] == "full"
    # sorted-order pseudonyms, consistent across tables
    assert [s["sample"] for s in data["samples"]] == ["S001", "S002", "S003", "S004", "S005"]
    assert [s["site"] for s in data["samples"]] == ["SITE1", "SITE1", "SITE1", "SITE2", "SITE2"]
    assert set(data["demix"]) == {f"S00{i}" for i in range(1, 6)}
    assert [r["sample"] for r in data["consensus"]] == ["S002", "S003", "S004"]
    assert data["samples"][0]["sample_type"] == "negative"
    flags = {f["code"]: f for f in data["run_flags"]}
    assert (
        flags["demix_failed"]["samples"] == ["S005"] and "S005" in flags["demix_failed"]["message"]
    )
    assert flags["consensus_failed"]["samples"] == ["S003"]
    assert flags["primer_scheme_mismatch"]["samples"] == ["S005"]
    assert (
        data["demix"]["S005"]["reason"]
        == "freyja demix failed for S005 in .../S005.variants.tsv (<user>)"
    )
    # site aggregates: aggregate name, sites, row labels
    agg = data["site_aggregates"]
    assert [a["row_label"] for a in agg] == ["AGG1_2026-W02", "AGG1_2026-W03"]
    assert agg[1]["sites_included"] == ["SITE1"] and agg[1]["sites_missing"] == ["SITE2"]
    assert (
        agg[1]["missing_reasons"] == {"SITE2": "not_ok:failed"}
        and agg[0]["aggregate_name"] == "AGG1"
    )
    inc = next(f for f in data["run_flags"] if f["code"] == "site_aggregate_incomplete")
    assert "AGG1 2026-W03 (missing SITE2)" in inc["message"]
    # the time series panel
    assert {x["site"] for x in data["time_series"]["sites"]} == {"SITE1", "SITE2"}
    # untouched: lineage names, numbers, dates, barcode versions, QC values
    assert data["demix"]["S002"]["groups"][1] == ["XBB.1.5* [Omicron (XBB.1.5.X)]", 0.3]
    assert data["demix"]["S002"]["lineages"] == [["BA.2", 0.7], ["XBB.1.5", 0.3]]
    assert [s["collection_date"] for s in data["samples"]][:2] == ["2026-01-05", "2026-01-05"]
    assert (
        data["samples"][4]["primer_scheme_best"] == "articv4.1"
        and data["samples"][4]["qc_status"] == "WARN"
    )
    assert data["provenance"]["freyja_data"]["barcode_version"] == "10_06_2026-05-25"
    assert data["provenance"]["command_line"].count(".../sites_AGG1.csv") == 1
    assert data["params_summary"]["primer_scheme"] == "midnight"


def test_pseudonym_map_file_is_correct_sorted_and_deterministic(tmp_path):
    build_a("full", out=tmp_path / "f1.json", map_out=tmp_path / "m1.csv")
    build_a("full", out=tmp_path / "f2.json", map_out=tmp_path / "m2.csv")
    text = (tmp_path / "m1.csv").read_text()
    assert text == (tmp_path / "m2.csv").read_text()
    assert text.splitlines() == [
        "type,original,pseudonym",
        "sample,NTC_Central_01,S001", "sample,WW_Central_001,S002", "sample,WW_Central_002,S003",
        "sample,WW_North_001,S004", "sample,WW_North_002.rep1,S005",
        "site,Central_WWTP,SITE1", "site,North_Plant,SITE2",
        "aggregate,Metro_Region,AGG1",
    ]  # fmt: skip
    a, b = (
        json.loads((tmp_path / "f1.json").read_text()),
        json.loads((tmp_path / "f2.json").read_text()),
    )
    a.pop("generated"), b.pop("generated")
    assert a == b
    # the map never leaks into the report data
    assert "NTC_Central_01" not in (tmp_path / "f1.json").read_text()


def test_map_is_only_written_for_full(tmp_path, capsys):
    build_a("paths", out=tmp_path / "p.json", map_out=tmp_path / "m.csv")
    build_a("none", out=tmp_path / "n.json", map_out=tmp_path / "m2.csv")
    assert not (tmp_path / "m.csv").exists() and not (tmp_path / "m2.csv").exists()
    assert "only used with --anonymize full" in capsys.readouterr().err


def test_pseudonym_width_grows_with_the_number_of_samples():
    names = {f"x{i:04d}" for i in range(1200)}
    mapping = brd.Anonymizer._pseudonyms(names, "S", 3)
    assert (
        mapping["x0000"] == "S0001"
        and mapping["x1199"] == "S1200"
        and len(set(mapping.values())) == 1200
    )


def test_qc_table_tsv_keeps_real_ids_even_when_anonymised(tmp_path):
    merged = tmp_path / "qc.merged.tsv"
    args = ["--run-info", str(A / "run_info.json"), "--qc", str(A / "qc_summary.tsv"), "--anonymize", "full",
            "--write-qc-tsv", str(merged), "--out", str(tmp_path / "f.json")]  # fmt: skip
    assert brd.main(args) == 0
    assert "WW_Central_001" in merged.read_text()
    assert "WW_Central_001" not in (tmp_path / "f.json").read_text()


@pytest.mark.parametrize(
    ("text", "expected"),
    [
        ("see /Users/jdoe/Data/run 1/x.csv here", "see .../x.csv here"),
        ("'/Users/jdoe/Data/Run 1 Plate A/s.csv'", "'.../s.csv'"),
        ("--out 'C:\\Users\\Jane Doe\\My Docs\\out.txt' -resume", "--out '.../out.txt' -resume"),
        ("C:\\Users\\Jane Doe\\My Docs\\out.txt --x", ".../out.txt --x"),
        ("\\\\server\\share\\data\\f.csv", ".../f.csv"),
        ("~/db/x.bed", ".../x.bed"),
        ("--a ./local/dir/f.tsv", "--a .../f.tsv"),
        ("in /home/jdoe now", "in ~ now"),
        ("mail me: j.doe@example.org", "mail me: <email>"),
        # untouched
        (
            "https://depot.galaxyproject.org/singularity/freyja:2.0.6--pyhdfd78af_0",
            "https://depot.galaxyproject.org/singularity/freyja:2.0.6--pyhdfd78af_0",
        ),
        (
            "quay.io/biocontainers/freyja:2.0.6--pyhdfd78af_0",
            "quay.io/biocontainers/freyja:2.0.6--pyhdfd78af_0",
        ),
        ("A / B / C", "A / B / C"),
        ("pct/x and a/b", "pct/x and a/b"),
        ("XBB.1.5* [Omicron (XBB.1.5.X)]", "XBB.1.5* [Omicron (XBB.1.5.X)]"),
        (
            "bcftools consensus failed at NC_045512.2:11287:",
            "bcftools consensus failed at NC_045512.2:11287:",
        ),
    ],
)
def test_scrub_text_cases(text, expected):
    anon = brd.Anonymizer("paths", {"x": "/Users/jdoe/q /home/jdoe C:\\Users\\Jane Doe\\z"})
    assert anon.users == ["Jane Doe", "jdoe"]
    assert anon.scrub_text(text) == expected


def test_user_name_variants_and_substrings_in_file_names():
    anon = brd.Anonymizer("paths", {"p": "/home/alice_b/data/x"})
    assert anon.users == ["alice_b"]
    assert (
        anon.scrub_text("owner alice_b wrote alice_b_samples.csv")
        == "owner <user> wrote alice_b_samples.csv"
    )
    assert anon.path_value("input", "/data/alice_b_samples.csv") == ".../<user>_samples.csv"
    assert anon.path_value("barcodes", "container") == "container"
    assert anon.path_value("genome", "assets/refs/x.fa") == ".../x.fa"
    assert anon.path_value("container", "quay.io/biocontainers/x:1") == "quay.io/biocontainers/x:1"
    assert anon.path_value("input", "https://x.org/a/b.csv") == "https://x.org/a/b.csv"


def test_hostile_ids_are_replaced_in_full_and_stay_inert_in_paths(tmp_path):
    evil = "<script>alert(1)</script>"
    qc = tmp_path / "qc.tsv"
    qc.write_text(
        "sample\tsample_type\tcollection_date\tsite\tqc_status\n"
        f"{evil}\tsample\t2026-01-01\tSite<b>\tFAIL\nWW_1\tsample\t2026-01-08\tSite<b>\tPASS\n"
    )
    out = tmp_path / "o.json"
    for mode, gone in (("full", True), ("paths", False)):
        assert (
            brd.main(
                [
                    "--run-info",
                    str(A / "run_info.json"),
                    "--qc",
                    str(qc),
                    "--anonymize",
                    mode,
                    "--out",
                    str(out),
                ]
            )
            == 0
        )
        raw = out.read_text()
        assert (evil not in raw) is gone and ("Site<b>" not in raw) is gone


def test_anonymize_option_validates_choices():
    with pytest.raises(SystemExit):
        brd.parse_args(["--run-info", "x", "--out", "y", "--anonymize", "partial"])


# =============================================================================================
# Dropout refinement (defining sites) and the unexplained `in_dataset` column
# =============================================================================================
J = FIX / "j_defining_sites"


def build_j(tmp_path, old_dropout=False, old_unexplained=False):
    """g_new_analyses inputs with the refined dropout / unexplained tables of j_defining_sites."""
    args = ["--run-info", str(J / "run_info.json")]
    for option, name in G_OPTIONS.items():
        if option in ("--dropout-impact", "--dropout-summary") and not old_dropout:
            args += [
                option,
                str(J / name.replace("dropout_impact_summary", "dropout_impact_summary")),
            ]
        elif option == "--unexplained" and not old_unexplained:
            args += [option, str(J / name)]
        else:
            args += [option, str(G / name)]
    out = tmp_path / "j.json"
    assert brd.main(args + ["--out", str(out)]) == 0
    return json.loads(out.read_text())


def test_dropout_defining_columns_are_read(tmp_path):
    data = build_j(tmp_path)
    impact = data["dropout_impact"]
    a = impact["SiteA_W03_1"]
    assert a["summary"]["defining_sites_total"] == 120 and a["summary"]["defining_sites_lost"] == 3
    assert a["summary"]["pct_defining_sites_lost"] == 2.5
    hit = a["amplicons"][1]
    assert hit["n_defining_sites"] == 3 and hit["n_differential_sites"] == 6
    assert hit["affected_lineages"] == [
        {"lineage": "XBB.1.5", "n_sites": 2},
        {"lineage": "BA.2", "n_sites": 1},
    ]
    assert (
        a["amplicons"][0]["n_defining_sites"] == 0
        and a["amplicons"][0]["n_differential_sites"] == 0
    )
    low = impact["LOWCOV_1"]
    assert (
        low["amplicons"][0]["n_defining_sites"] is None
        and low["summary"]["pct_defining_sites_lost"] is None
    )
    # a lineage below the affect threshold is listed but does not raise the flag
    quiet = impact["SiteB_W05_1"]["amplicons"][1]
    assert (
        quiet["affected_lineages"] == [{"lineage": "PQ.2.1.17", "n_sites": 1}]
        and quiet["affects_detected_lineage"] is False
    )
    s = {x["sample"]: x for x in data["samples"]}
    assert (
        s["SiteA_W03_1"]["dropout_pct_defining_lost"] == 2.5
        and s["SiteA_W03_1"]["dropout_pct_informative_lost"] == 10.62
    )
    assert s["LOWCOV_1"]["dropout_pct_defining_lost"] is None


def test_old_dropout_tables_without_the_new_columns_still_load(tmp_path):
    data = build_g(tmp_path)  # g_new_analyses holds the old-format tables
    a = data["dropout_impact"]["SiteA_W03_1"]
    assert (
        a["amplicons"][1]["n_defining_sites"] is None
        and a["amplicons"][1]["n_differential_sites"] is None
    )
    assert (
        a["summary"]["defining_sites_total"] is None
        and a["summary"]["pct_defining_sites_lost"] is None
    )
    assert data["samples"][0]["dropout_pct_defining_lost"] is None
    assert "input_problems" not in flag_codes(data)


def test_dropout_flag_text_names_defining_mutations_and_the_threshold(tmp_path):
    data = build_j(tmp_path)
    (flag,) = [f for f in data["run_flags"] if f["code"] == "dropout_affects_lineages"]
    assert flag["samples"] == ["SiteA_W03_1", "SiteB_W05_1"]
    assert flag["message"].startswith(
        "Amplicon dropout hides defining mutations (those that separate a lineage from its parent) of lineages detected at >= 5% in 2 sample(s)"
    )
    assert "abundance estimates for those lineages may be biased" in flag["message"]
    assert data["thresholds"]["dropout_min_affect_abundance"] == 0.05
    assert data["params_summary"]["dropout_min_affect_abundance"] == 0.05
    # default threshold text when the run did not record the param
    plain = build_g(tmp_path)
    assert (
        "detected at >= 5%"
        in next(f for f in plain["run_flags"] if f["code"] == "dropout_affects_lineages")["message"]
    )


def test_unexplained_in_dataset_and_the_floored_p_value(tmp_path):
    data = build_j(tmp_path)
    rows = data["unexplained_mutations"]
    assert [r["in_dataset"] for r in rows] == [False, False, False, True, True, True, False]
    assert rows[0]["p_value"] == 1e-300
    (flag,) = [f for f in data["run_flags"] if f["code"] == "unexplained_mutations"]
    assert (
        "2 mutation(s) unexplained by the demixed lineages recur in several samples (1 of them at sites no lineage in the dataset models) (first: S:E180D in 3 samples)"
        in flag["message"]
    )
    # the older table without the column: unknown, never false, and no blind-spot wording
    old = build_g(tmp_path)
    assert all(r["in_dataset"] is None for r in old["unexplained_mutations"])
    assert (
        "no lineage in the dataset models"
        not in next(f for f in old["run_flags"] if f["code"] == "unexplained_mutations")["message"]
    )
