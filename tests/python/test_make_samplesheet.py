import csv
import shutil
from pathlib import Path

import make_samplesheet as ms

FIX = Path(__file__).resolve().parents[1] / "fixtures" / "samplesheet"


def read(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh))


def run(args, capsys):
    rc = ms.main([str(a) for a in args])
    cap = capsys.readouterr()
    return rc, cap.out, cap.err


def by_sample(rows):
    return {r["sample"]: r for r in rows}


def test_default_illumina_with_metadata(tmp_path, capsys):
    """Golden path: ID parsing (S-number/lane stripped, Undetermined excluded, underscores kept),
    negative-control detection, absolute paths and the metadata join."""
    out = tmp_path / "s.csv"
    rc, _, err = run(
        [FIX / "fastq_illumina", "--metadata", FIX / "metadata.csv", "-o", out], capsys
    )
    assert rc == 0, err
    rows = read(out)
    assert [r["sample"] for r in rows] == ["NEG_ntc", "Plant_NC", "Site-A.1", "Site_1_week_2"]
    rec = by_sample(rows)
    assert rec["NEG_ntc"]["sample_type"] == "negative"
    assert rec["Plant_NC"]["sample_type"] == "negative"
    r = rec["Site_1_week_2"]
    assert r["sample_type"] == "sample"
    assert Path(r["fastq_1"]).is_absolute()
    assert r["fastq_2"].endswith("Site_1_week_2_S12_L001_R2_001.fastq.gz")
    assert (r["collection_date"], r["site"]) == ("2024-03-01", "SiteOne")
    assert "Ghost_sample" in err and "no metadata row" in err


def test_results_demix_skips_incomplete_sample_with_reason(tmp_path, capsys):
    out = tmp_path / "d.csv"
    rc, _, err = run(["--from-results", FIX / "results_new", "--mode", "demix", "-o", out], capsys)
    assert rc == 0, err
    rows = by_sample(read(out))
    assert set(rows) == {"NEG_ntc", "Site_1_week_2"}  # Broken_1 skipped
    assert "skipped 'Broken_1'" in err and "missing depths" in err
    r = rows["Site_1_week_2"]
    assert r["variants"].endswith("samples/Site_1_week_2/freyja/Site_1_week_2.variants.tsv")
    assert r["depths"].endswith("Site_1_week_2.depths")
    assert r["fastp_json"].endswith("samples/Site_1_week_2/fastp/Site_1_week_2.fastp.json")
    assert Path(r["variants"]).is_absolute()
    # carried over from summary/samplesheet.csv
    assert r["collection_date"] == "2024-03-01" and r["site"] == "SiteOne"
    assert rows["NEG_ntc"]["sample_type"] == "negative"
    assert list(read(out)[0]) == [
        "sample",
        "variants",
        "depths",
        "fastp_json",
        "sample_type",
        "collection_date",
        "site",
    ]


def _legacy_tree(tmp_path):
    """Legacy results: Good (usable), Failed (zero-byte .demix: `|| touch` in the old pipeline),
    NoFiles (freyja dir without any demix/variants)."""
    root = tmp_path / "old"
    for sample, size in (("Good", 40), ("Failed", 0)):
        d = root / sample / "freyja"
        d.mkdir(parents=True)
        (d / f"{sample}.demix").write_text("x" * size)
        (d / f"{sample}.depths").write_text("ref\t1\tA\t5\n")
        (d / f"{sample}.freyja_variants.tsv").write_text("REGION\tPOS\n")
    (root / "NoFiles" / "freyja").mkdir(parents=True)
    return root


def test_zero_byte_legacy_demix_is_reported_in_skipped_out(tmp_path, capsys):
    root = _legacy_tree(tmp_path)
    out, skipped = tmp_path / "sheet.csv", tmp_path / "skipped.tsv"
    rc, _, err = run(
        ["--from-results", root, "--mode", "aggregate", "-o", out, "--skipped-out", skipped], capsys
    )
    assert rc == 0, err
    assert [r["sample"] for r in read(out)] == ["Good"]
    rows = [line.split("\t") for line in skipped.read_text().splitlines()]
    assert rows[0] == ["sample", "reason", "sample_type", "collection_date", "site"]
    reasons = {r[0]: r[1] for r in rows[1:]}
    assert reasons["Failed"] == "empty demix file (failed in the old pipeline)"
    assert reasons["NoFiles"].startswith("missing demix")


def test_results_demix_carries_per_sample_protocol_and_scheme(tmp_path, capsys):
    """PLAN D20: a previous run's protocol/primer_scheme cells reach the demix sheet (and only then are the columns written)."""
    results = tmp_path / "results"
    shutil.copytree(FIX / "results_new", results)
    (results / "summary" / "samplesheet.csv").write_text(
        "sample,fastq_1,fastq_2,sample_type,protocol,primer_scheme\n"
        "Site_1_week_2,a_R1.fastq.gz,a_R2.fastq.gz,sample,amplicon,midnight\n"
        "NEG_ntc,n_R1.fastq.gz,n_R2.fastq.gz,negative,capture,\n"
    )
    out = tmp_path / "d.csv"
    rc, _, err = run(["--from-results", results, "--mode", "demix", "-o", out], capsys)
    assert rc == 0, err
    rows = by_sample(read(out))
    assert (rows["Site_1_week_2"]["protocol"], rows["Site_1_week_2"]["primer_scheme"]) == (
        "amplicon",
        "midnight",
    )
    assert (rows["NEG_ntc"]["protocol"], rows["NEG_ntc"]["primer_scheme"]) == ("capture", "")


def test_results_of_a_panel_run_give_one_row_per_sample_and_target(tmp_path, capsys):
    """Panel layout (summary/targets.tsv, samples/<id>/<target>/freyja): (sample, target) rows with a target column,
    metadata from summary/<target>/samplesheet.csv, a skipped TSV with a target column, --targets filters, and the
    sheet passes check_samplesheet (duplicates are on sample AND target)."""
    import check_samplesheet as cs

    root = tmp_path / "panel"
    (root / "summary").mkdir(parents=True)
    (root / "summary" / "targets.tsv").write_text(
        "sample\ttarget\tanalysed\nS1\tA\ttrue\nS1\tB\ttrue\n"
    )
    for sample, target, with_demix in (("S1", "A", True), ("S1", "B", True), ("S2", "B", False)):
        d = root / "samples" / sample / target / "freyja"
        d.mkdir(parents=True)
        (d / f"{sample}.variants.tsv").write_text("REGION\tPOS\n")
        (d / f"{sample}.depths").write_text("ref\t1\tA\t5\n")
        if with_demix:
            (d / f"{sample}.demix").write_text("x\n")
    (root / "samples" / "S1" / "fastp").mkdir()
    (root / "samples" / "S1" / "fastp" / "S1.fastp.json").write_text("{}")
    (root / "summary" / "B").mkdir()
    (root / "summary" / "B" / "samplesheet.csv").write_text(
        "sample,sample_type,collection_date,site\nS1,sample,2026-01-05,SiteA\nS2,negative,2026-01-06,SiteB\n"
    )
    out, skipped = tmp_path / "sheet.csv", tmp_path / "skipped.tsv"
    rc, _, err = run(
        ["--from-results", root, "--mode", "aggregate", "-o", out, "--skipped-out", skipped], capsys
    )
    assert rc == 0, err
    rows = read(out)
    assert list(rows[0])[:3] == ["sample", "target", "demix"]
    assert [(r["sample"], r["target"]) for r in rows] == [("S1", "A"), ("S1", "B")]
    assert rows[1]["collection_date"] == "2026-01-05" and rows[1]["site"] == "SiteA"
    assert skipped.read_text().splitlines() == [
        "sample\ttarget\treason\tsample_type\tcollection_date\tsite",
        "S2\tB\tmissing demix (panel layout)\tnegative\t2026-01-06\tSiteB",
    ]
    assert cs.validate(out, "aggregate")["ok"]
    rc, _, err = run(
        ["--from-results", root, "--mode", "demix", "--targets", "A", "-o", out], capsys
    )
    assert rc == 0, err
    rows = read(out)
    assert [(r["sample"], r["target"]) for r in rows] == [("S1", "A")]
    assert rows[0]["fastp_json"].endswith("samples/S1/fastp/S1.fastp.json")


def _fake_run(root, rows):
    """Results of a run in the new layout; rows = (id, date, site, demix status)."""
    (root / "summary").mkdir(parents=True)
    sheet = ["sample,fastq_1,fastq_2,sample_type,collection_date,site"]
    status = ["sample\tstatus\treason"]
    for sample, date, site, st in rows:
        d = root / "samples" / sample / "freyja"
        d.mkdir(parents=True)
        (d / f"{sample}.demix").write_text("x\n")
        sheet.append(f"{sample},a,b,sample,{date},{site}")
        status.append(f"{sample}\t{st}\t")
    (root / "summary" / "samplesheet.csv").write_text("\n".join(sheet) + "\n")
    (root / "summary" / "demix_status.tsv").write_text("\n".join(status) + "\n")
    return root


def test_several_results_dirs_are_filtered_and_combined(tmp_path, capsys):
    """--from-results is repeatable; the filters are ANDed (ID regex, site, date range, demix status)."""
    run1 = _fake_run(
        tmp_path / "Run1",
        [
            ("S1", "2026-01-01", "A", "ok"),
            ("S2", "2026-01-08", "A", "ok"),
            ("S3", "2026-01-15", "A", "low_confidence"),
        ],
    )
    run2 = _fake_run(
        tmp_path / "Run2",
        [
            ("S8", "2026-01-22", "A", "ok"),
            ("S9", "2026-01-29", "B", "ok"),
            ("S10", "2026-02-05", "A", "ok"),
        ],
    )
    out = tmp_path / "sheet.csv"
    rc, _, err = run(
        ["--from-results", run1, "--from-results", run2, "--mode", "aggregate", "-o", out]
        + ["--exclude-id", "^S10$", "--site", "A", "--date-from", "2026-01-05", "--status", "ok"],
        capsys,
    )
    assert rc == 0, err
    assert [r["sample"] for r in read(out)] == [
        "S2",
        "S8",
    ]  # S1 too early, S3 low_confidence, S9 site B, S10 excluded
    assert list(read(out)[0]) == ["sample", "demix", "sample_type", "collection_date", "site"]
    assert "4 sample(s) removed by filters" in err


def test_duplicate_ids_across_results_dirs_and_list(tmp_path, capsys):
    run1 = _fake_run(
        tmp_path / "Run1", [("S1", "2026-01-01", "A", "ok"), ("S2", "2026-01-08", "A", "ok")]
    )
    run2 = _fake_run(tmp_path / "Run2", [("S1", "2026-02-01", "A", "ok")])
    base = [
        "--from-results",
        run1,
        "--from-results",
        run2,
        "--mode",
        "aggregate",
        "-o",
        tmp_path / "s.csv",
    ]
    rc, _, err = run(base, capsys)
    assert rc == 1 and "'S1'" in err and "Run1, Run2" in err and not (tmp_path / "s.csv").exists()
    rc, _, err = run([*base, "--on-duplicate", "prefix"], capsys)
    rows = read(tmp_path / "s.csv")
    assert rc == 0, err
    assert [r["sample"] for r in rows] == ["Run1__S1", "S2", "Run2__S1"]
    assert rows[0]["demix"].endswith("Run1/samples/S1/freyja/S1.demix")  # path unchanged
    rc, _, err = run([*base[:-2], "--on-duplicate", "first", "-o", tmp_path / "f.csv"], capsys)
    assert [r["sample"] for r in read(tmp_path / "f.csv")] == ["S1", "S2"]
    assert "1 duplicate row(s) dropped" in err
    (tmp_path / "s.csv").unlink()
    rc, out, _ = run([*base, "--list"], capsys)
    assert rc == 0 and not (tmp_path / "s.csv").exists()
    lines = out.splitlines()
    assert lines[0].split() == [
        "run",
        "sample",
        "sample_type",
        "site",
        "collection_date",
        "demix_status",
        "selected",
        "reason",
    ]
    assert sum("duplicate ID" in ln for ln in lines) == 2 and any(
        ln.split()[:2] == ["Run1", "S2"] and "yes" in ln for ln in lines
    )
