import json
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
    assert (
        s["WW_002"]["consensus_length"] == 29903
        and s["WW_003.rep1"]["consensus_coverage_pct"] == 96.99
    )
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


def _loaded(tmp_path, name, text):
    path = tmp_path / name
    path.write_text(text)
    return path


def test_failed_and_skipped_samples_always_get_a_status_and_reason(tmp_path):
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

    # a sample whose demix status is ok but has no lineages, or has no status at all, is flagged
    tmp_path = tmp_path / "missing"
    tmp_path.mkdir()
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


def test_depths_only_input_does_not_leak_track_metadata_as_samples(tmp_path):
    depths = tmp_path / "S1.depths"
    depths.write_text("".join(f"NC_045512.2\t{i}\tA\t20\n" for i in range(1, 301)))
    run_info = _loaded(tmp_path, "run_info.json", '{"mode": "demix"}')
    out = tmp_path / "report_data.json"
    assert brd.main(["--run-info", str(run_info), "--depths", str(depths), "--out", str(out)]) == 0
    data = json.loads(out.read_text())
    assert [s["sample"] for s in data["samples"]] == ["S1"]
    assert not data["run_flags"] or "contig" not in json.dumps(data["run_flags"])


def test_multi_contig_run_gets_segments_offsets_and_a_not_covered_flag(tmp_path):
    """Two contigs (A 600 bp covered, B 400 bp empty): one axis with segments, `<contig>/<name>` amplicons
    shifted by the contig offset, the segments table and the `segments_not_covered` warning."""
    (tmp_path / "run_info.json").write_text(
        '{"mode": "full", "freyja_contig": "A", "contigs": [{"name": "A", "length": 600},'
        ' {"name": "B", "length": 400}]}'
    )
    lines = [f"A\t{p}\tN\t50" for p in range(1, 601)] + [f"B\t{p}\tN\t0" for p in range(1, 401)]
    (tmp_path / "S1.depths").write_text("\n".join(lines) + "\n")
    (tmp_path / "amp.tsv").write_text(
        "sample\tamplicon\tstart\tend\tmean_depth\tdropped\tcontig\n"
        "S1\tp1\t1\t200\t50\tfalse\tA\nS1\tp1\t1\t200\t0\ttrue\tB\n"
    )
    (tmp_path / "cov.tsv").write_text(
        "sample\tcontig\tlength\tmean_depth\tmedian_depth\tpct_ge_1x\tpct_ge_10x\tpct_ge_min_depth\t"
        "amplicons_total\tamplicons_dropped\nS1\tA\t600\t50\t50\t100\t100\t100\t1\t0\n"
        "S1\tB\t400\t0\t0\t0\t0\t0\t1\t1\n"
    )
    (tmp_path / "cons.tsv").write_text(
        "sample\tcontig\tlength\tn_count\tconsensus_coverage_pct\tambiguous_count\tstatus\n"
        "S1\tA\t600\t0\t100\t0\tok\nS1\tB\t400\t400\t0\t0\tempty\n"
    )
    out = tmp_path / "report_data.json"
    args = ["--run-info", tmp_path / "run_info.json", "--depths", tmp_path, "--amplicons"]
    args += [tmp_path / "amp.tsv", "--contig-coverage", tmp_path / "cov.tsv"]
    args += ["--consensus-stats-contigs", tmp_path / "cons.tsv", "--out", out]
    assert brd.main([str(a) for a in args]) == 0
    data = json.loads(out.read_text())
    segs = data["coverage_tracks"]["segments"]
    assert [(s["contig"], s["offset"], s["length"]) for s in segs] == [
        ("A", 0, 600),
        ("B", 600, 400),
    ]
    assert data["coverage_tracks"]["genome_length"] == 1000
    amps = [(a["amplicon"], a["start"], a["contig"]) for a in data["amplicons"]["S1"]]
    assert amps == [("A/p1", 1, "A"), ("B/p1", 601, "B")]
    rows = {r["contig"]: r for r in data["segments"]["samples"]["S1"]}
    assert data["segments"]["freyja_contig"] == "A"
    assert not rows["A"]["not_covered"] and rows["B"]["not_covered"]
    flag = next(f for f in data["run_flags"] if f["code"] == "segments_not_covered")
    assert flag["level"] == "warn" and "S1 (B)" in flag["message"]


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


def build_a(mode, out):
    args = ["--run-info", str(A / "run_info.json")]
    for option, name in A_OPTIONS.items():
        args += [option, str(A / name)]
    assert brd.main([*args, "--anonymize", mode, "--out", str(out)]) == 0
    return Path(out).read_text()


def test_anonymisation_leaks_no_ids_sites_paths_or_user_names(tmp_path):
    """--anonymize full: no original sample/site/aggregate name, path, email or user name survives
    anywhere in the report data; lineages, numbers and QC values are left alone. paths mode: no
    path/user/email, but IDs stay."""
    raw = build_a("full", tmp_path / "f.json")
    data = json.loads(raw)
    for leak in A_LEAKS + A_IDS + A_SITES + ["Metro_Region", "Central", "North_", "WW_", "NTC_"]:
        assert leak not in raw, leak
    assert data["anonymized"] == "full"
    assert [s["sample"] for s in data["samples"]] == ["S001", "S002", "S003", "S004", "S005"]
    assert [s["site"] for s in data["samples"]] == ["SITE1", "SITE1", "SITE1", "SITE2", "SITE2"]
    assert set(data["demix"]) == {f"S00{i}" for i in range(1, 6)}
    flags = {f["code"]: f for f in data["run_flags"]}
    assert (
        flags["demix_failed"]["samples"] == ["S005"] and "S005" in flags["demix_failed"]["message"]
    )
    assert [a["row_label"] for a in data["site_aggregates"]] == ["AGG1_2026-W02", "AGG1_2026-W03"]
    assert data["demix"]["S002"]["lineages"] == [["BA.2", 0.7], ["XBB.1.5", 0.3]]
    assert data["provenance"]["freyja_data"]["barcode_version"] == "10_06_2026-05-25"

    raw = build_a("paths", tmp_path / "p.json")
    for leak in A_LEAKS:
        assert leak not in raw, leak
    for sid in A_IDS:
        assert sid in raw
