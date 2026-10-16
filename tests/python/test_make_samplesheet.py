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


def test_default_illumina(tmp_path, capsys):
    out = tmp_path / "s.csv"
    rc, _, err = run([FIX / "fastq_illumina", "-o", out], capsys)
    assert rc == 0, err
    rows = read(out)
    ids = [r["sample"] for r in rows]
    assert ids == sorted(ids)
    # underscores kept, S-number and lane stripped, Undetermined excluded
    assert ids == ["NEG_ntc", "Plant_NC", "Site-A.1", "Site_1_week_2"]
    rec = by_sample(rows)
    assert rec["NEG_ntc"]["sample_type"] == "negative"
    assert rec["Plant_NC"]["sample_type"] == "negative"
    assert rec["Site_1_week_2"]["sample_type"] == "sample"
    r = rec["Site_1_week_2"]
    assert Path(r["fastq_1"]).is_absolute()
    assert r["fastq_1"].endswith("Site_1_week_2_S12_L001_R1_001.fastq.gz")
    assert r["fastq_2"].endswith("Site_1_week_2_S12_L001_R2_001.fastq.gz")
    assert list(rows[0]) == [
        "sample",
        "fastq_1",
        "fastq_2",
        "sample_type",
        "collection_date",
        "site",
    ]


def test_stdout_default(capsys):
    rc, out, _ = run([FIX / "fastq_illumina"], capsys)
    assert rc == 0
    assert out.startswith("sample,fastq_1,fastq_2,")


def test_relative_paths(tmp_path, capsys):
    out = tmp_path / "s.csv"
    rc, _, _ = run([FIX / "fastq_custom", "--suffix", "_R1.fq.gz", "--relative", "-o", out], capsys)
    assert rc == 0
    r = read(out)[0]
    assert not Path(r["fastq_1"]).is_absolute()
    assert (tmp_path / r["fastq_1"]).resolve().is_file()


def test_custom_suffix_and_r2_suffix_accepted(tmp_path, capsys):
    out = tmp_path / "s.csv"
    rc, _, err = run([FIX / "fastq_custom", "--suffix", "_R1.fq.gz", "-o", out], capsys)
    assert rc == 0, err
    rows = read(out)
    assert [r["sample"] for r in rows] == ["alpha_1", "beta_2"]
    assert rows[0]["fastq_2"].endswith("alpha_1_R2.fq.gz")
    # giving the R2 suffix is normalised to R1 (Snakemake behaviour)
    out2 = tmp_path / "s2.csv"
    rc, _, _ = run([FIX / "fastq_custom", "--suffix", "_R2.fq.gz", "-o", out2], capsys)
    assert rc == 0
    assert read(out2) == rows


def test_custom_suffix_no_match_is_empty_error(tmp_path, capsys):
    out = tmp_path / "s.csv"
    rc, _, err = run([FIX / "fastq_custom", "-o", out], capsys)  # default suffix: nothing matches
    assert rc == 1
    assert "no samples found" in err
    assert not out.exists()


def test_suffix_must_be_fastq_gz(capsys):
    rc, _, err = run([FIX / "fastq_custom", "--suffix", "_R1.txt"], capsys)
    assert rc == 1
    assert "must end with" in err


def test_missing_r2_fails(tmp_path, capsys):
    out = tmp_path / "s.csv"
    rc, _, err = run([FIX / "fastq_missing_r2", "-o", out], capsys)
    assert rc == 1
    assert "missing R2 for 'lonely'" in err
    assert not out.exists()


def test_duplicate_id_collision_fails(tmp_path, capsys):
    out = tmp_path / "s.csv"
    rc, _, err = run([FIX / "fastq_collision", "-o", out], capsys)
    assert rc == 1
    assert "duplicate sample ID 'dup'" in err
    assert "dup_S1_L001" in err and "dup_S7_L001" in err
    assert not out.exists()


def test_invalid_id_fails(capsys):
    rc, _, err = run([FIX / "fastq_bad_id"], capsys)
    assert rc == 1
    assert "does not match" in err


def test_recursive(tmp_path, capsys):
    root = tmp_path / "reads"
    (root / "run1").mkdir(parents=True)
    for n in ("a_S1_L001_R1_001.fastq.gz", "a_S1_L001_R2_001.fastq.gz"):
        shutil.copy(FIX / "fastq_missing_r2" / "ok_S1_L001_R1_001.fastq.gz", root / "run1" / n)
    rc, _, _ = run([root], capsys)
    assert rc == 1  # not found without --recursive
    rc, out, _ = run([root, "--recursive"], capsys)
    assert rc == 0
    assert "\na," in out


def test_missing_dir(capsys):
    rc, _, err = run([FIX / "nope"], capsys)
    assert rc == 1
    assert "not found" in err


def test_no_match_non_illumina_name(tmp_path, capsys):
    (tmp_path / "weird_L001_R1_001.fastq.gz").write_bytes(b"")
    rc, _, err = run([tmp_path], capsys)
    assert rc == 1
    assert "cannot extract a sample ID" in err


def test_custom_negative_and_exclude_patterns(tmp_path, capsys):
    out = tmp_path / "s.csv"
    rc, _, _ = run(
        [
            FIX / "fastq_illumina",
            "--negative-pattern",
            "week_2",
            "--exclude-pattern",
            "NEG_",
            "-o",
            out,
        ],
        capsys,
    )
    assert rc == 0
    rec = by_sample(read(out))
    assert "NEG_ntc" not in rec
    assert "Undetermined" in rec  # default exclusion was replaced, so it is now kept
    assert rec["Site_1_week_2"]["sample_type"] == "negative"
    assert rec["Plant_NC"]["sample_type"] == "sample"


def test_metadata_join_with_warnings(tmp_path, capsys):
    out = tmp_path / "s.csv"
    rc, _, err = run(
        [FIX / "fastq_illumina", "--metadata", FIX / "metadata.csv", "-o", out], capsys
    )
    assert rc == 0, err
    rec = by_sample(read(out))
    assert rec["Site_1_week_2"]["collection_date"] == "2024-03-01"
    assert rec["Site_1_week_2"]["site"] == "SiteOne"
    assert rec["Site-A.1"]["site"] == "SiteA"
    assert rec["NEG_ntc"]["collection_date"] == ""
    assert "Ghost_sample" in err  # metadata row with no sample
    assert "no metadata row" in err and "NEG_ntc" in err  # sample without metadata


def test_metadata_overrides_sample_type(tmp_path, capsys):
    meta = tmp_path / "m.csv"
    meta.write_text("sample,sample_type\nNEG_ntc,sample\nSite-A.1,negative\n")
    out = tmp_path / "s.csv"
    rc, _, _ = run([FIX / "fastq_illumina", "--metadata", meta, "-o", out], capsys)
    assert rc == 0
    rec = by_sample(read(out))
    assert rec["NEG_ntc"]["sample_type"] == "sample"
    assert rec["Site-A.1"]["sample_type"] == "negative"


def test_metadata_bad_dates_fail(tmp_path, capsys):
    out = tmp_path / "s.csv"
    rc, _, err = run(
        [FIX / "fastq_illumina", "--metadata", FIX / "metadata_bad_date.csv", "-o", out], capsys
    )
    assert rc == 1
    assert "01/03/2024" in err and "2024-02-30" in err
    assert not out.exists()


def test_metadata_requires_sample_column_and_unique(tmp_path, capsys):
    meta = tmp_path / "m.csv"
    meta.write_text("id,site\nx,y\n")
    rc, _, err = run([FIX / "fastq_illumina", "--metadata", meta], capsys)
    assert rc == 1 and "no 'sample' column" in err
    meta.write_text("sample,site\nx,y\nx,z\n")
    rc, _, err = run([FIX / "fastq_illumina", "--metadata", meta], capsys)
    assert rc == 1 and "duplicate sample 'x'" in err


def test_metadata_bad_sample_type(tmp_path, capsys):
    meta = tmp_path / "m.csv"
    meta.write_text("sample,sample_type\nSite-A.1,control\n")
    rc, _, err = run([FIX / "fastq_illumina", "--metadata", meta], capsys)
    assert rc == 1 and "sample_type" in err


def test_argument_errors(capsys):
    assert run([], capsys)[0] == 2
    assert run(["--from-results", FIX], capsys)[0] == 2
    assert run([FIX, "--mode", "demix"], capsys)[0] == 2
    assert run([FIX / "fastq_custom", "--from-results", FIX, "--mode", "demix"], capsys)[0] == 2


# ------------------------------------------------------------------ --from-results
def test_results_new_layout_demix(tmp_path, capsys):
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


def test_results_new_layout_aggregate(tmp_path, capsys):
    out = tmp_path / "a.csv"
    rc, _, err = run(
        ["--from-results", FIX / "results_new", "--mode", "aggregate", "-o", out], capsys
    )
    assert rc == 0, err
    rows = by_sample(read(out))
    assert rows["Site_1_week_2"]["demix"].endswith("Site_1_week_2.demix.tsv")
    assert list(read(out)[0]) == ["sample", "demix", "sample_type", "collection_date", "site"]


def test_results_legacy_layout_demix_and_gz(tmp_path, capsys):
    out = tmp_path / "d.csv"
    rc, _, err = run(
        ["--from-results", FIX / "results_legacy", "--mode", "demix", "-o", out], capsys
    )
    assert rc == 0, err
    rows = by_sample(read(out))
    assert set(rows) == {"Old_1", "NC_blank", "Gz_only"}
    assert rows["Old_1"]["variants"].endswith("Old_1/freyja/Old_1.freyja_variants.tsv")
    assert rows["Old_1"]["fastp_json"].endswith("Old_1/fastp/Old_1_fastp.json")
    assert rows["NC_blank"]["sample_type"] == "negative"
    assert rows["Gz_only"]["variants"].endswith("Gz_only/variants/Gz_only.freyja_variants.tsv.gz")
    assert rows["Gz_only"]["fastp_json"] == ""
    assert "gzipped" in err and "Gz_only" in err


def test_results_legacy_layout_aggregate(tmp_path, capsys):
    out = tmp_path / "a.csv"
    rc, _, _ = run(
        ["--from-results", FIX / "results_legacy", "--mode", "aggregate", "-o", out], capsys
    )
    assert rc == 0
    rows = by_sample(read(out))
    assert rows["Old_1"]["demix"].endswith("Old_1/freyja/Old_1.demix")


def test_results_mixed_layout_warns_and_prefers_new(tmp_path, capsys):
    out = tmp_path / "d.csv"
    rc, _, err = run(
        ["--from-results", FIX / "results_mixed", "--mode", "demix", "-o", out], capsys
    )
    assert rc == 0
    assert "both result layouts" in err
    assert "Both_1" in err
    rows = by_sample(read(out))
    assert "/samples/Both_1/" in rows["Both_1"]["variants"]
    assert "/Legacy_only/freyja/" in rows["Legacy_only"]["variants"]


def test_results_nothing_found_errors(tmp_path, capsys):
    (tmp_path / "samples" / "x" / "freyja").mkdir(parents=True)
    out = tmp_path / "d.csv"
    rc, _, err = run(["--from-results", tmp_path, "--mode", "demix", "-o", out], capsys)
    assert rc == 1
    assert "no usable samples" in err
    assert not out.exists()
    rc, _, err = run(["--from-results", tmp_path / "missing", "--mode", "demix"], capsys)
    assert rc == 1 and "not found" in err


def test_results_relative(tmp_path, capsys):
    out = tmp_path / "d.csv"
    rc, _, _ = run(
        ["--from-results", FIX / "results_new", "--mode", "aggregate", "--relative", "-o", out],
        capsys,
    )
    assert rc == 0
    for r in read(out):
        assert not Path(r["demix"]).is_absolute()
        assert (tmp_path / r["demix"]).resolve().is_file()


def test_results_metadata_option(tmp_path, capsys):
    meta = tmp_path / "m.csv"
    meta.write_text("sample,site\nOld_1,Plant9\n")
    out = tmp_path / "a.csv"
    rc, _, _ = run(
        [
            "--from-results",
            FIX / "results_legacy",
            "--mode",
            "aggregate",
            "--metadata",
            meta,
            "-o",
            out,
        ],
        capsys,
    )
    assert rc == 0
    assert by_sample(read(out))["Old_1"]["site"] == "Plant9"


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


def test_zero_byte_demix_is_skipped_and_reported_in_aggregate_mode(tmp_path, capsys):
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


def test_demix_mode_does_not_care_about_the_demix_file_but_reports_missing_files(tmp_path, capsys):
    root = _legacy_tree(tmp_path)
    out, skipped = tmp_path / "sheet.csv", tmp_path / "skipped.tsv"
    rc, _, err = run(
        ["--from-results", root, "--mode", "demix", "-o", out, "--skipped-out", skipped], capsys
    )
    assert rc == 0, err
    assert [r["sample"] for r in read(out)] == ["Failed", "Good"]
    assert {
        r[0]: r[1] for r in (line.split("\t") for line in skipped.read_text().splitlines()[1:])
    } == {"NoFiles": "missing variants, depths (legacy layout)"}


def test_skipped_out_has_only_a_header_when_nothing_is_skipped(tmp_path, capsys):
    out, skipped = tmp_path / "sheet.csv", tmp_path / "skipped.tsv"
    rc, _, err = run(
        [
            "--from-results",
            FIX / "results_new",
            "--mode",
            "demix",
            "-o",
            out,
            "--skipped-out",
            skipped,
        ],
        capsys,
    )
    assert rc == 0, err
    assert (
        skipped.read_text().splitlines()[0] == "sample\treason\tsample_type\tcollection_date\tsite"
    )


def test_skipped_negative_control_keeps_its_type_dates_and_site(tmp_path, capsys):
    """A control without a .demix is only known from discovery: its sample_type comes from the
    results' samplesheet, else from the NEG_/NC_ pattern."""
    root = tmp_path / "old"
    (root / "summary").mkdir(parents=True)
    (root / "summary" / "samplesheet.csv").write_text(
        "sample,sample_type,collection_date,site\nBlank_7,negative,2026-01-05,SiteA\n"
    )
    for sample in ("Good", "Blank_7", "NEG_x"):
        d = root / "samples" / sample / "freyja"
        d.mkdir(parents=True)
        (d / f"{sample}.variants.tsv").write_text("REGION\n")
        (d / f"{sample}.depths").write_text("ref\t1\tA\t5\n")
        if sample == "Good":
            (d / "Good.demix").write_text("x" * 20)
    out, skipped = tmp_path / "sheet.csv", tmp_path / "skipped.tsv"
    rc, _, err = run(
        ["--from-results", root, "--mode", "aggregate", "-o", out, "--skipped-out", skipped], capsys
    )
    assert rc == 0, err
    rows = {r[0]: r for r in (line.split("\t") for line in skipped.read_text().splitlines()[1:])}
    assert rows["Blank_7"][2:] == ["negative", "2026-01-05", "SiteA"]
    assert rows["NEG_x"][2] == "negative"


def test_a_db_dir_inside_the_results_is_not_a_sample(tmp_path, capsys):
    """--db_dir defaults to <outdir>/db in the test profiles; --compare_barcodes leaves db/freyja behind."""
    root = tmp_path / "res"
    for sample in ("Good",):
        d = root / "samples" / sample / "freyja"
        d.mkdir(parents=True)
        (d / f"{sample}.variants.tsv").write_text("REGION\n")
        (d / f"{sample}.depths").write_text("ref\t1\tA\t5\n")
        (d / f"{sample}.demix").write_text("x" * 20)
    (root / "db" / "freyja").mkdir(parents=True)  # empty, as left by an update
    # a dataset directory under another name is recognised by its files
    ds = root / "dataset" / "freyja" / "2026-10-06T05-25"
    ds.mkdir(parents=True)
    (ds / "last_barcode_update.txt").write_text("10_06_2026-05-25\n")
    for mode in ("demix", "aggregate"):
        out, skipped = tmp_path / f"{mode}.csv", tmp_path / f"{mode}.skipped.tsv"
        rc, _, err = run(
            ["--from-results", root, "--mode", mode, "-o", out, "--skipped-out", skipped], capsys
        )
        assert rc == 0, err
        assert [r["sample"] for r in read(out)] == ["Good"], mode
        assert skipped.read_text().splitlines() == [
            "sample\treason\tsample_type\tcollection_date\tsite"
        ]


def test_viral_load_from_metadata_is_written_only_when_present(tmp_path, capsys):
    meta = tmp_path / "m.csv"
    meta.write_text("sample,site,viral_load\nSite_1_week_2,P1,1.5e5\nSite-A.1,P1,\n")
    out = tmp_path / "s.csv"
    rc, _, err = run([FIX / "fastq_illumina", "--metadata", meta, "-o", out], capsys)
    assert rc == 0, err
    rows = read(out)
    assert list(rows[0])[-1] == "viral_load"
    rec = by_sample(rows)
    assert rec["Site_1_week_2"]["viral_load"] == "1.5e5"
    assert rec["Site-A.1"]["viral_load"] == ""
    # no sample has a viral load: the column is not written
    meta.write_text("sample,site\nSite_1_week_2,P1\n")
    rc, _, err = run([FIX / "fastq_illumina", "--metadata", meta, "-o", out], capsys)
    assert rc == 0, err
    assert "viral_load" not in read(out)[0]


def test_viral_load_must_be_a_positive_number(tmp_path, capsys):
    out = tmp_path / "s.csv"
    for bad in ("abc", "-3", "0", "nan", "inf"):
        meta = tmp_path / "m.csv"
        meta.write_text(f"sample,viral_load\nSite_1_week_2,{bad}\n")
        rc, _, err = run([FIX / "fastq_illumina", "--metadata", meta, "-o", out], capsys)
        assert rc == 1 and "viral_load" in err, bad


def test_viral_load_is_carried_over_from_a_results_samplesheet(tmp_path, capsys):
    root = tmp_path / "old"
    for sample in ("Good", "Other"):
        d = root / "samples" / sample / "freyja"
        d.mkdir(parents=True)
        (d / f"{sample}.demix").write_text("x" * 40)
        (d / f"{sample}.depths").write_text("ref\t1\tA\t5\n")
        (d / f"{sample}.freyja_variants.tsv").write_text("REGION\tPOS\n")
    (root / "summary").mkdir()
    (root / "summary" / "samplesheet.csv").write_text(
        "sample,fastq_1,fastq_2,sample_type,collection_date,site,viral_load\n"
        "Good,a,b,sample,2026-01-05,S1,2500\n"
        "Other,a,b,sample,2026-01-12,S1,-1\n"
    )
    out = tmp_path / "a.csv"
    rc, _, err = run(["--from-results", root, "--mode", "aggregate", "-o", out], capsys)
    assert rc == 0, err
    rec = by_sample(read(out))
    assert rec["Good"]["viral_load"] == "2500"
    assert rec["Other"]["viral_load"] == ""  # invalid values are not carried over
