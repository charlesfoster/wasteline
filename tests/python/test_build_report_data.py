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
    extras = {"negative_controls", "boot_intervals"}  # only on when the feature produced data
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
    assert set(brd.SUMMARY_PARAMS) <= known, sorted(set(brd.SUMMARY_PARAMS) - known)
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
