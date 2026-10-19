"""bin/reanalyse_runs.py: re-run past results with new barcodes (dry run / --execute) and diff old vs new.

Fixtures (tests/fixtures/reanalyse/): run_info_full.json / run_info_demix.json are trimmed copies of the real
run_info.json shapes (see tests/fixtures/replay); old/ and new/ hold tiny `summary/` tables in the real formats
(freyja_aggregated.tsv, demix_status.tsv, lineage_summary.csv); legacy/ is a Snakemake-layout results directory.
`make_run` builds results trees under tmp_path. A fake `nextflow` (a Python script) stands in for the real one.
"""

import ast
import csv
import json
import re
import shlex
import shutil
import stat
import subprocess
import sys
import tempfile
from pathlib import Path

import pytest
import reanalyse_runs as ra
import test_host_python_compat as compat

FIX = Path(__file__).resolve().parents[1] / "fixtures" / "reanalyse"
REPO = Path(__file__).resolve().parents[2]
OLD_ISO = "2026-10-06T05-25"
OLD_SHA = "429d293f1a93fb7f3c572ece17e157c07ca5223aa22df807687e7833c7d88710"


@pytest.fixture(autouse=True)
def isolated_home(tmp_path, monkeypatch):
    """The default db_dir (~/.wasteline/db) must never be the user's real one."""
    monkeypatch.setenv("HOME", str(tmp_path / "home"))


def info(name="full"):
    return json.loads((FIX / f"run_info_{name}.json").read_text())


def make_run(
    parent,
    name,
    *,
    run_info=None,
    summary="old",
    samples=("A", "B"),
    start=None,
    mutate=None,
):
    """A results directory <parent>/<name> with run_info.json, variants+depths per sample and a summary/."""
    run = Path(parent) / name
    data = run_info if run_info is not None else info("full")
    if start:
        data["start"] = start
    if mutate:
        mutate(data)
    (run / "pipeline_info").mkdir(parents=True)
    (run / "pipeline_info" / "run_info.json").write_text(json.dumps(data))
    for sample in samples:
        freyja = run / "samples" / sample / "freyja"
        freyja.mkdir(parents=True)
        (freyja / f"{sample}.variants.tsv").write_text("REGION\tPOS\n")
        (freyja / f"{sample}.depths").write_text("NC\t1\tA\t10\n")
    if summary:
        shutil.copytree(FIX / summary / "summary", run / "summary")
    return run


def snapshot(path):
    """{relative path: bytes or None for directories} of a tree, to prove nothing was modified."""
    result = {}
    for item in sorted(Path(path).rglob("*")):
        result[str(item.relative_to(path))] = item.read_bytes() if item.is_file() else None
    return result


def main(*argv):
    return ra.main([str(a) for a in argv])


def printed_commands(out):
    """Every `nextflow run ...` command in the output, as a token list (continuation lines joined)."""
    commands = []
    for line in re.sub(r" \\\n\s*", " ", out).splitlines():
        try:
            parts = shlex.split(line)
        except ValueError:
            continue
        if len(parts) > 2 and Path(parts[0]).name == "nextflow" and parts[1] == "run":
            commands.append(parts)
    return commands


def run_commands(out):
    return [c for c in printed_commands(out) if "--previous_results" in c]


def flag(tokens, name):
    return tokens[tokens.index(name) + 1]


# ----------------------------------------------------------------------------------------------- discovery


def test_root_is_a_results_directory(tmp_path):
    run = make_run(tmp_path, "run1")
    runs, messages = ra.discover([run])
    assert [r.name for r in runs] == ["run1"] and not messages
    assert runs[0].mode == "full" and runs[0].barcode == OLD_ISO and runs[0].n_samples == 2


def test_root_is_a_folder_of_results_and_depth_is_limited(tmp_path):
    make_run(tmp_path, "run1")
    make_run(tmp_path / "2026", "run2")  # depth 2 below the root: found
    make_run(tmp_path / "a" / "b", "run3")  # depth 3: not found
    runs, _ = ra.discover([tmp_path])
    assert sorted(r.name for r in runs) == ["run1", "run2"]


def test_a_run_is_not_searched_further(tmp_path):
    outer = make_run(tmp_path, "run1")
    make_run(outer, "nested")  # a results directory inside a run: it is part of that run
    runs, _ = ra.discover([tmp_path])
    assert [r.name for r in runs] == ["run1"]


def test_overlapping_roots_list_each_run_once(tmp_path):
    run = make_run(tmp_path, "run1")
    runs, _ = ra.discover([tmp_path, run, tmp_path / "."])
    assert len(runs) == 1


def test_legacy_snakemake_layout_has_unknown_barcodes(tmp_path):
    shutil.copytree(FIX / "legacy", tmp_path / "old_snakemake")
    runs, _ = ra.discover([tmp_path])
    [legacy] = runs
    assert legacy.kind == "legacy" and legacy.mode == "legacy"
    assert legacy.barcode == "unknown" and legacy.freyja_version == "unknown"
    # MIX_A (freyja/) and MIX_B (variants/ + freyja/depths) are usable; NEG_E has variants+depths too
    assert legacy.n_samples == 3
    assert legacy.skip_reason is None
    assert legacy.date  # the modification time


def test_legacy_directory_given_directly(tmp_path):
    shutil.copytree(FIX / "legacy", tmp_path / "snk")
    runs, _ = ra.discover([tmp_path / "snk"])
    assert [r.kind for r in runs] == ["legacy"]


def test_a_freyja_dataset_store_is_not_a_legacy_run(tmp_path):
    store = tmp_path / "root" / "x" / "freyja" / "2026-10-06T05-25"
    store.mkdir(parents=True)
    (store / "last_barcode_update.txt").write_text("10_06_2026-05-25")
    runs, messages = ra.discover([tmp_path / "root"])
    assert runs == [] and messages


def test_nothing_found_is_a_friendly_error(tmp_path, capsys):
    (tmp_path / "empty").mkdir()
    assert main(tmp_path / "empty") == 1
    err = capsys.readouterr().err
    assert "no results found" in err and "run_info.json" in err and "Traceback" not in err


def test_samples_without_run_info_get_a_hint(tmp_path):
    (tmp_path / "r" / "samples").mkdir(parents=True)
    _, messages = ra.discover([tmp_path / "r"])
    assert "no pipeline_info/run_info.json" in messages[0]


def test_missing_root_and_broken_run_info(tmp_path, capsys):
    assert main(tmp_path / "nope") == 1
    assert "not a directory" in capsys.readouterr().err
    run = make_run(tmp_path, "bad")
    (run / "pipeline_info" / "run_info.json").write_text("{not json")
    assert main(run) == 1
    assert "cannot read" in capsys.readouterr().err


def test_demix_output_cannot_be_reanalysed(tmp_path):
    make_run(tmp_path, "run1", start="2026-10-16T10:00:00+11:00")
    derived = make_run(tmp_path, "run1_reanalysis_x", run_info=info("demix"), samples=())
    runs, _ = ra.discover([tmp_path])
    by_name = {r.name: r for r in runs}
    assert by_name["run1"].skip_reason is None
    assert "keeps no variants/depths" in by_name["run1_reanalysis_x"].skip_reason
    assert "/results/run" in by_name["run1_reanalysis_x"].skip_reason  # points at the source
    assert derived.exists()


def test_failed_and_sampleless_runs_are_skipped(tmp_path):
    failed = make_run(tmp_path, "failed", mutate=lambda d: d.update(success=False))
    empty = make_run(tmp_path, "empty", samples=())
    runs, _ = ra.discover([failed, empty])
    reasons = {r.name: r.skip_reason for r in runs}
    assert "did not succeed" in reasons["failed"]
    assert "no sample" in reasons["empty"]


# ----------------------------------------------------------------------------------------------- selection


def three_runs(tmp_path):
    make_run(tmp_path, "w1", start="2026-09-20T10:00:00+10:00")
    make_run(tmp_path, "w2", start="2026-10-05T10:00:00+11:00")
    make_run(tmp_path, "w3", start="2026-10-16T10:00:00+11:00")
    runs, _ = ra.discover([tmp_path])
    return runs


class Args:
    runs = None
    since = None
    last = None


def selected(found, **kwargs):
    args = Args()
    for key, value in kwargs.items():
        setattr(args, key, value)
    return [r.name for r in ra.select_runs(found, args)]


def test_selection_since_last_runs(tmp_path):
    runs = three_runs(tmp_path)
    assert selected(runs) == ["w1", "w2", "w3"]
    assert selected(runs, since="2026-10-05") == ["w2", "w3"]  # the date itself is included
    assert selected(runs, since="2026-10-06") == ["w3"]
    assert selected(runs, last=2) == ["w2", "w3"]
    assert selected(runs, last=0) == []
    assert selected(runs, runs="w3,w1") == ["w1", "w3"]  # deterministic: by date
    assert selected(runs, runs="w1,w2,w3", since="2026-10-01", last=1) == ["w3"]
    assert selected(runs, runs=str(tmp_path / "w2")) == ["w2"]  # a path works too


def test_unknown_run_name_lists_what_exists(tmp_path):
    runs = three_runs(tmp_path)
    with pytest.raises(ra.Problem, match="no run named 'zzz'.*w1"):
        selected(runs, runs="zzz")


def test_last_counts_only_reanalysable_runs(tmp_path):
    make_run(tmp_path, "w1", start="2026-09-20T10:00:00+10:00")
    make_run(tmp_path, "w2", start="2026-10-05T10:00:00+11:00")
    make_run(tmp_path, "w3_reanalysis_x", run_info=info("demix"), samples=())
    runs, _ = ra.discover([tmp_path])
    names = selected(runs, last=1)
    assert "w2" in names and "w1" not in names
    assert "w3_reanalysis_x" in names  # listed (with its reason), but not counted


def test_since_must_be_a_date(tmp_path, capsys):
    make_run(tmp_path, "w1")
    assert main(tmp_path, "--since", "yesterday") == 1
    assert "YYYY-MM-DD" in capsys.readouterr().err


# ----------------------------------------------------------------------------------------------- barcodes target


def test_target_forms(tmp_path):
    assert ra.resolve_target("latest", None).pending
    assert ra.resolve_target("container", None).kind == "container"
    assert ra.resolve_target("10_20_2026-04-10", None).iso == "2026-10-20T04-10"
    assert ra.resolve_target("2026-10-20T04-10", None).iso == "2026-10-20T04-10"
    plain = ra.resolve_target("2026-10-20", None)
    assert plain.kind == "date" and plain.iso is None and plain.date == "2026-10-20"
    with pytest.raises(ra.Problem, match="not 'latest'"):
        ra.resolve_target("nope", None)


def test_date_target_resolves_against_the_local_store(tmp_path):
    db = tmp_path / "db"
    for stamp in ("2026-10-20T01-00", "2026-10-20T09-30", "2026-10-21T00-10"):
        (db / "freyja" / stamp).mkdir(parents=True)
    target = ra.resolve_target("2026-10-20", str(db))
    assert target.iso == "2026-10-20T09-30"  # the newest of that day, like the pipeline


def test_dataset_directory_and_file_targets(tmp_path):
    dataset = tmp_path / "ds"
    dataset.mkdir()
    (dataset / "last_barcode_update.txt").write_text("10_20_2026-04-10\n")
    (dataset / "usher_barcodes.feather").write_bytes(b"abc")
    target = ra.resolve_target(str(dataset), None)
    assert target.kind == "path" and target.iso == "2026-10-20T04-10" and target.sha
    barcodes = tmp_path / "my barcodes.csv"
    barcodes.write_text("lineage,A1T\n")
    assert ra.resolve_target(str(barcodes), None).label == "my barcodes"


# ----------------------------------------------------------------------------------------------- dry run / commands


def dry(tmp_path, capsys, *argv):
    rc = main(*argv)
    return rc, capsys.readouterr().out


def test_dry_run_prints_table_and_exact_command_and_writes_nothing(tmp_path, capsys):
    run = make_run(tmp_path, "run1")
    before = snapshot(tmp_path)
    rc, out = dry(tmp_path, capsys, run, "--barcodes", "2026-10-20T04-10")
    assert rc == 0
    for word in ("run1", "2026-10-16", "full", OLD_ISO, "2.0.6", "reanalyse"):
        assert word in out
    [cmd] = run_commands(out)
    assert cmd[:2] == ["nextflow", "run"] and cmd[2] == str(REPO)
    assert flag(cmd, "--mode") == "demix"
    assert flag(cmd, "--previous_results") == str(run)
    assert flag(cmd, "--outdir") == f"{run}_reanalysis_2026-10-20T04-10"
    assert flag(cmd, "--barcodes") == "2026-10-20T04-10"
    assert flag(cmd, "-params-file") == f"{run}_reanalysis_2026-10-20T04-10.params.json"
    assert flag(cmd, "-profile") == "docker"  # the recorded "docker,locked" without locked
    assert "NOT executed" in out
    assert snapshot(tmp_path) == before  # not even the params file


def test_params_carried_over_and_params_left_out(tmp_path, capsys):
    def mutate(d):
        d["params"].update(
            freyja_covcut=12,  # not the default: carried
            consensus_bogus=1,  # not a parameter of this pipeline: left out with a note
            qc_pass_min_pct=60,  # default: not repeated
        )

    run = make_run(tmp_path, "run1", mutate=mutate)
    _, out = dry(tmp_path, capsys, run, "--barcodes", "container")
    carried = out.split("carried over:", 1)[1].splitlines()[0]
    assert "min_genome_coverage_for_demix=40" in carried
    assert "freyja_covcut=12" in carried
    assert (
        "primer_scheme=qiagenv1" in carried
    )  # the effective scheme, even though it is the default
    for left_out in (
        "save_bam",
        "fastp_",
        "input",
        "outdir",
        "barcodes",
        "qc_pass_min_pct",
        "fasta",
    ):
        assert left_out not in carried
    assert "consensus_bogus is not a parameter" in out


def test_params_file_content(tmp_path):
    run = make_run(
        tmp_path,
        "run1",
        mutate=lambda d: d["params"].update(freyja_covcut=12, annotate_variants=True),
    )
    [loaded], _ = ra.discover([run])
    notes = []
    carried = ra.carry_params(loaded, ra.schema_defaults(REPO), notes)
    assert carried["annotate_variants"] is True and carried["freyja_covcut"] == 12
    assert carried["min_genome_coverage_for_demix"] == 40
    assert "mode" not in carried and "previous_results" not in carried and "barcodes" not in carried


def test_capture_protocol_is_carried_and_primer_scheme_is_not(tmp_path):
    def mutate(d):
        d["params"].update(protocol="capture", primer_scheme="qiagenv1")
        d["protocol"], d["primer_trimming"] = "capture", False

    run = make_run(tmp_path, "cap", mutate=mutate)
    [loaded], _ = ra.discover([run])
    carried = ra.carry_params(loaded, ra.schema_defaults(REPO), [])
    assert carried["protocol"] == "capture"
    assert "primer_scheme" not in carried and "primer_bed" not in carried


def test_custom_primer_bed_replaces_the_scheme(tmp_path):
    bed = tmp_path / "my primers.bed"
    bed.write_text("NC_045512.2\t30\t54\tx_1_LEFT\t1\t+\n")
    run = make_run(tmp_path, "r", mutate=lambda d: d["params"].update(primer_bed=str(bed)))
    [loaded], _ = ra.discover([run])
    carried = ra.carry_params(loaded, ra.schema_defaults(REPO), [])
    assert carried["primer_bed"] == str(bed) and "primer_scheme" not in carried


def test_missing_primer_bed_is_dropped_with_a_note_and_the_scheme_stays(tmp_path):
    run = make_run(tmp_path, "r", mutate=lambda d: d["params"].update(primer_bed="/gone/p.bed"))
    [loaded], _ = ra.discover([run])
    notes = []
    carried = ra.carry_params(loaded, ra.schema_defaults(REPO), notes)
    assert "primer_bed" not in carried and carried["primer_scheme"] == "qiagenv1"
    assert any("primer_bed /gone/p.bed not found" in n for n in notes)


def test_custom_groups_file_used_while_it_exists(tmp_path):
    groups = tmp_path / "my groups.csv"
    groups.write_text("group,include\nJN.1*,JN.1*\n")

    def mutate(d):
        d["params"]["lineage_groups"] = str(groups)
        d["lineage_groups_source"] = "custom"

    run = make_run(tmp_path, "r", mutate=mutate)
    [loaded], _ = ra.discover([run])
    carried = ra.carry_params(loaded, ra.schema_defaults(REPO), [])
    assert carried["lineage_groups"] == str(groups)


def test_custom_groups_file_missing_falls_back_to_the_published_copy_then_bundled(tmp_path):
    def mutate(d):
        d["params"]["lineage_groups"] = "/gone/groups.csv"
        d["lineage_groups_source"] = "custom"

    run = make_run(tmp_path, "r", mutate=mutate)
    [loaded], _ = ra.discover([run])
    notes = []
    assert "lineage_groups" not in ra.carry_params(loaded, None, notes)  # no copy: bundled default
    assert any("bundled default" in n for n in notes)
    used = run / "summary" / "lineage_groups_used.csv"
    used.write_text("group,include\nJN.1*,JN.1*\n")
    notes = []
    assert ra.carry_params(loaded, None, notes)["lineage_groups"] == str(used)
    assert any("copy published in the results" in n for n in notes)
    loaded.run_info["lineage_groups_sha256"] = "0" * 64  # a copy that is not the file that was used
    assert "lineage_groups" not in ra.carry_params(loaded, None, [])


def test_bundled_groups_and_missing_mutations_file(tmp_path):
    def mutate(d):
        d["lineage_groups_source"] = "bundled"
        d["params"]["mutations_of_interest"] = "/gone/moi.csv"

    run = make_run(tmp_path, "r", mutate=mutate)
    [loaded], _ = ra.discover([run])
    notes = []
    carried = ra.carry_params(loaded, None, notes)
    assert "lineage_groups" not in carried and "mutations_of_interest" not in carried
    assert any("mutations_of_interest /gone/moi.csv not found" in n for n in notes)


def test_profile_handling(tmp_path):
    class A:
        profile = None

    for number, (recorded, expected) in enumerate(
        [
            ("docker,locked", "docker"),
            ("singularity", "singularity"),
            ("test,docker", "docker"),
            ("test", "docker"),
            ("", "docker"),
        ]
    ):
        run = make_run(
            tmp_path / str(number),
            "r",
            mutate=lambda d, recorded=recorded: d.update(profile=recorded),
        )
        [loaded], _ = ra.discover([run])
        assert ra.profile_for(loaded, A)[0] == expected
    A.profile = "conda"
    assert ra.profile_for(loaded, A)[0] == "conda"


def test_legacy_run_command_has_no_params_and_a_note(tmp_path, capsys):
    shutil.copytree(FIX / "legacy", tmp_path / "snk")
    rc, out = dry(tmp_path, capsys, tmp_path / "snk", "--barcodes", "2026-10-20T04-10")
    assert rc == 0 and "unknown" in out and "legacy" in out
    assert "carried over: (none)" in out and "legacy results" in out
    [cmd] = run_commands(out)
    assert flag(cmd, "--previous_results") == str(tmp_path / "snk")


def test_paths_with_spaces_are_quoted_and_round_trip(tmp_path, capsys):
    run = make_run(tmp_path / "my results", "run 1")
    _, out = dry(
        tmp_path,
        capsys,
        run,
        "--barcodes",
        "2026-10-20T04-10",
        "--db_dir",
        str(tmp_path / "my db"),
        "--extra-args",
        "-work-dir '/scratch/my work' -resume",
        "--pipeline-dir",
        str(REPO),
    )
    assert shlex.quote(str(run)) in out and "'" in out
    [cmd] = run_commands(out)
    assert flag(cmd, "--previous_results") == str(run)
    assert flag(cmd, "--db_dir") == str(tmp_path / "my db")
    assert flag(cmd, "-work-dir") == "/scratch/my work" and "-resume" in cmd


def test_outdir_suffix_placeholder_and_literal(tmp_path, capsys):
    run = make_run(tmp_path, "run1")
    _, out = dry(
        tmp_path, capsys, run, "--barcodes", "2026-10-20T04-10", "--outdir-suffix", "_b{barcode}x"
    )
    assert flag(run_commands(out)[0], "--outdir") == f"{run}_b2026-10-20T04-10x"
    _, out = dry(
        tmp_path, capsys, run, "--barcodes", "2026-10-20T04-10", "--outdir-suffix", "_redo"
    )
    assert flag(run_commands(out)[0], "--outdir") == f"{run}_redo"
    assert main(run, "--outdir-suffix", " ") == 1


def test_bad_barcodes_and_bad_pipeline_dir(tmp_path, capsys):
    run = make_run(tmp_path, "run1")
    assert main(run, "--barcodes", "nope") == 1
    assert "not 'latest'" in capsys.readouterr().err
    assert main(run, "--barcodes", "container", "--pipeline-dir", str(tmp_path)) == 1
    assert "no main.nf" in capsys.readouterr().err
    assert main(run, "--barcodes", "container", "--pipeline-dir", str(tmp_path / "gone")) == 1
    assert "does not exist" in capsys.readouterr().err


# ----------------------------------------------------------------------------------------------- up to date / refusal


def test_up_to_date_by_stamp_forms_date_and_sha(tmp_path, capsys):
    run = make_run(tmp_path, "run1")
    for spec in ("10_06_2026-05-25", OLD_ISO, "2026-10-06"):
        rc, out = dry(tmp_path, capsys, run, "--barcodes", spec)
        assert rc == 0 and "up to date" in out and "Nothing to re-analyse" in out
        assert not run_commands(out)
    dataset = tmp_path / "ds"
    dataset.mkdir()
    (dataset / "last_barcode_update.txt").write_text("10_20_2026-04-10")
    (dataset / "usher_barcodes.feather").write_bytes(b"abc")
    rc, out = dry(tmp_path, capsys, run, "--barcodes", str(dataset))
    assert "reanalyse" in out
    # the same barcodes file as recorded (by sha256) is up to date whatever the folder is called
    [loaded], _ = ra.discover([run])
    same = ra.Target("x", "path", sha=OLD_SHA)
    assert same.is_current(loaded)
    assert not ra.Target("x", "path", sha="1" * 64, iso=OLD_ISO).is_current(loaded)


def test_a_different_day_or_stamp_is_not_up_to_date(tmp_path, capsys):
    run = make_run(tmp_path, "run1")
    for spec in ("2026-10-07", "2026-10-06T05-26", "10_07_2026-05-25"):
        _, out = dry(tmp_path, capsys, run, "--barcodes", spec)
        assert "reanalyse" in out and "up to date" not in out


def test_container_and_latest_are_never_skipped_as_up_to_date(tmp_path, capsys):
    run = make_run(tmp_path, "run1")
    _, out = dry(tmp_path, capsys, run, "--barcodes", "container")
    assert "up to date" not in out and run_commands(out)
    _, out = dry(tmp_path, capsys, run, "--barcodes", "latest")
    assert "latest (resolved when run)" in out and "--mode update_db" in out
    [cmd] = run_commands(out)
    assert flag(cmd, "--barcodes") == "latest" and "resolved barcode version" in out


def test_custom_and_legacy_runs_are_never_up_to_date(tmp_path):
    def mutate(d):
        d["freyja_data"].update(source="custom", barcode_version="custom", barcode_version_iso=None)

    run = make_run(tmp_path, "c", mutate=mutate)
    [loaded], _ = ra.discover([run])
    assert loaded.barcode == "custom"
    assert not ra.resolve_target("2026-10-06", None).is_current(loaded)


def test_existing_outdir_is_refused_and_a_counter_is_optional(tmp_path, capsys):
    run = make_run(tmp_path, "run1")
    existing = Path(f"{run}_reanalysis_2026-10-20T04-10")
    existing.mkdir()
    (existing / "keep.txt").write_text("precious")
    rc, out = dry(tmp_path, capsys, run, "--barcodes", "2026-10-20T04-10")
    assert rc == 0 and "output exists" in out and not run_commands(out)
    rc, out = dry(tmp_path, capsys, run, "--barcodes", "2026-10-20T04-10", "--force-new-suffix")
    assert flag(run_commands(out)[0], "--outdir") == f"{existing}_2"
    Path(f"{existing}_2").mkdir()
    _, out = dry(tmp_path, capsys, run, "--barcodes", "2026-10-20T04-10", "--force-new-suffix")
    assert flag(run_commands(out)[0], "--outdir") == f"{existing}_3"
    assert (existing / "keep.txt").read_text() == "precious"


def test_unreanalysable_runs_are_listed_with_their_reason(tmp_path, capsys):
    make_run(tmp_path, "run1")
    make_run(tmp_path, "run1_reanalysis_x", run_info=info("demix"), samples=())
    _, out = dry(tmp_path, capsys, tmp_path, "--barcodes", "2026-10-20T04-10")
    assert "skip (mode demix keeps no variants/depths" in out
    assert len(run_commands(out)) == 1


# ----------------------------------------------------------------------------------------------- execute (fake nextflow)

FAKE_NEXTFLOW = """#!{python}
import json, os, shutil, sys
argv = sys.argv[1:]
log = os.environ["FAKE_NF_LOG"]
with open(log, "a") as handle:
    handle.write(json.dumps(argv) + "\\n")
def opt(name):
    return argv[argv.index(name) + 1] if name in argv else None
mode = opt("--mode")
fix = os.environ["FAKE_NF_FIXTURES"]
if mode == "update_db":
    out = opt("--outdir")
    os.makedirs(out + "/pipeline_info")
    json.dump({{"barcode_version": "10_20_2026-04-10", "barcode_version_iso": "2026-10-20T04-10"}},
              open(out + "/pipeline_info/freyja_data.json", "w"))
    sys.exit(int(os.environ.get("FAKE_NF_UPDATE_EXIT", "0")))
previous = opt("--previous_results")
if os.path.basename(previous) in os.environ.get("FAKE_NF_FAIL", "").split(","):
    sys.exit(3)
out = opt("--outdir")
shutil.copytree(fix + "/new/summary", out + "/summary")
os.makedirs(out + "/pipeline_info")
info = json.load(open(fix + "/run_info_demix.json"))
info["freyja_data"]["barcode_version_iso"] = opt("--barcodes") if opt("--barcodes") != "latest" else "?"
json.dump(info, open(out + "/pipeline_info/run_info.json", "w"))
"""


@pytest.fixture
def fake_nextflow(tmp_path, monkeypatch):
    script = tmp_path / "bin_dir" / "nextflow"
    script.parent.mkdir()
    script.write_text(FAKE_NEXTFLOW.format(python=sys.executable))
    script.chmod(script.stat().st_mode | stat.S_IEXEC)
    log = tmp_path / "nf_calls.jsonl"
    monkeypatch.setenv("FAKE_NF_LOG", str(log))
    monkeypatch.setenv("FAKE_NF_FIXTURES", str(FIX))
    monkeypatch.delenv("FAKE_NF_FAIL", raising=False)
    monkeypatch.delenv("FAKE_NF_UPDATE_EXIT", raising=False)

    def calls():
        if not log.exists():
            return []
        return [json.loads(line) for line in log.read_text().splitlines()]

    class Fake:
        def __str__(self):
            return str(script)

    fake = Fake()
    fake.calls = calls
    return fake


def test_execute_runs_every_command_and_writes_params_and_diff(tmp_path, capsys, fake_nextflow):
    results = tmp_path / "results"
    run1 = make_run(results, "run1", start="2026-10-14T10:00:00+11:00")
    run2 = make_run(results, "run2", start="2026-10-15T10:00:00+11:00")
    before = {str(p): snapshot(p) for p in (run1, run2)}
    report = tmp_path / "all.csv"
    summary = tmp_path / "all.md"
    rc = main(
        results,
        "--barcodes",
        "2026-10-20T04-10",
        "--execute",
        "--nextflow",
        fake_nextflow,
        "--report",
        report,
        "--summary-md",
        summary,
    )
    out = capsys.readouterr().out
    assert rc == 0
    calls = fake_nextflow.calls()
    assert [flag(c, "--previous_results") for c in calls] == [
        str(run1),
        str(run2),
    ]  # sequential, by date
    for run in (run1, run2):
        new = Path(f"{run}_reanalysis_2026-10-20T04-10")
        params = json.loads(Path(f"{new}.params.json").read_text())
        assert params["min_genome_coverage_for_demix"] == 40
        assert (new / "summary" / "reanalysis_diff.csv").is_file()
        assert (new / "summary" / "reanalysis_diff.md").is_file()
    assert "Done: 2 re-analysed, 0 failed" in out
    assert {str(p): snapshot(p) for p in (run1, run2)} == before  # the originals are untouched
    with open(report) as handle:
        rows = list(csv.DictReader(handle))
    assert {r["run"] for r in rows} == {"run1", "run2"} and "sample" in rows[0]
    assert "## run1" in summary.read_text() and "## run2" in summary.read_text()


def test_execute_stops_on_the_first_failure_unless_keep_going(
    tmp_path, capsys, fake_nextflow, monkeypatch
):
    results = tmp_path / "results"
    make_run(results, "run1", start="2026-10-14T10:00:00+11:00")
    make_run(results, "run2", start="2026-10-15T10:00:00+11:00")
    monkeypatch.setenv("FAKE_NF_FAIL", "run1")
    argv = [results, "--barcodes", "2026-10-20T04-10", "--execute", "--nextflow", fake_nextflow]
    assert main(*argv) == 1
    assert len(fake_nextflow.calls()) == 1
    assert "Stopping" in capsys.readouterr().err
    assert main(*argv, "--keep-going") == 1
    out = capsys.readouterr().out
    assert len(fake_nextflow.calls()) == 3  # run1 again (fails), then run2
    assert "1 re-analysed, 1 failed" in out
    assert Path(f"{results}/run2_reanalysis_2026-10-20T04-10/summary/reanalysis_diff.csv").is_file()


def test_execute_latest_runs_update_db_once_and_shares_the_stamp(tmp_path, capsys, fake_nextflow):
    results = tmp_path / "results"
    run1 = make_run(results, "run1", start="2026-10-14T10:00:00+11:00")
    run2 = make_run(results, "run2", start="2026-10-15T10:00:00+11:00")
    rc = main(results, "--execute", "--nextflow", fake_nextflow, "--db_dir", tmp_path / "db")
    assert rc == 0
    calls = fake_nextflow.calls()
    assert [flag(c, "--mode") for c in calls] == ["update_db", "demix", "demix"]
    assert flag(calls[0], "--db_dir") == str(tmp_path / "db")
    assert {flag(c, "--barcodes") for c in calls[1:]} == {"2026-10-20T04-10"}
    assert Path(f"{run1}_reanalysis_2026-10-20T04-10").is_dir()
    assert Path(f"{run2}_reanalysis_2026-10-20T04-10").is_dir()
    assert not list(Path(tempfile.gettempdir()).glob("wasteline_update_db_*"))


def test_execute_latest_failure_of_update_db_aborts_everything(
    tmp_path, capsys, fake_nextflow, monkeypatch
):
    results = tmp_path / "results"
    make_run(results, "run1")
    monkeypatch.setenv("FAKE_NF_UPDATE_EXIT", "1")
    assert main(results, "--execute", "--nextflow", fake_nextflow) == 1
    assert "update_db failed" in capsys.readouterr().err
    assert [flag(c, "--mode") for c in fake_nextflow.calls()] == ["update_db"]


def test_execute_latest_that_resolves_to_the_recorded_version_is_up_to_date(
    tmp_path, capsys, fake_nextflow
):
    results = tmp_path / "results"

    def mutate(d):
        d["freyja_data"].update(
            barcode_version="10_20_2026-04-10", barcode_version_iso="2026-10-20T04-10"
        )
        d["freyja_data"].pop("barcodes_sha256")

    make_run(results, "run1", mutate=mutate)
    assert main(results, "--execute", "--nextflow", fake_nextflow) == 0
    out = capsys.readouterr().out
    assert "up to date" in out
    assert [flag(c, "--mode") for c in fake_nextflow.calls()] == ["update_db"]


def test_nothing_selected_runs_nothing(tmp_path, capsys, fake_nextflow):
    results = tmp_path / "results"
    make_run(results, "run1", start="2026-10-14T10:00:00+11:00")
    assert main(results, "--since", "2026-11-01", "--execute", "--nextflow", fake_nextflow) == 0
    assert "No run matches" in capsys.readouterr().out and fake_nextflow.calls() == []


# ----------------------------------------------------------------------------------------------- tolerant parsing


def test_summarized_forms():
    expect = {"JN.1* (JN.1.X)": 0.6, "Other": 0.4}
    assert ra.parse_summarized("[('JN.1* (JN.1.X)', 0.6), ('Other', 0.4)]") == expect
    assert (
        ra.parse_summarized("[('JN.1* (JN.1.X)', np.float64(0.6)), ('Other', np.float64(0.4))]")
        == expect
    )
    assert ra.parse_summarized('[("JN.1* (JN.1.X)", 0.6), ("Other", 0.4)]') == expect
    assert ra.parse_summarized("{'JN.1* (JN.1.X)': 0.6, 'Other': 0.4}") == expect
    assert ra.parse_summarized("[('JN.1* (JN.1.X)', 6e-1), ('Other', 0.4),]") == expect
    assert ra.parse_summarized("[('A', 0.5), ('B', 0.5 junk") == {
        "A": 0.5
    }  # damaged tail: what parses
    for empty in ("", "nan", "[]", None, "  "):
        assert ra.parse_summarized(empty) == {}


def test_lineage_and_abundance_forms():
    expect = {"JN.1": 0.6, "KP.3": 0.4}
    forms = [
        ("JN.1 KP.3", "0.60000000 0.40000000"),
        ("['JN.1' 'KP.3']", "[0.6 0.4]"),
        ("['JN.1', 'KP.3']", "[0.6, 0.4]"),
        ("array(['JN.1', 'KP.3'], dtype=object)", "array([0.6, 0.4])"),
        ("['JN.1'\n 'KP.3']", "[6.0e-01\n 0.4]"),
        ("JN.1,KP.3", "0.6,0.4"),
        ("['JN.1' 'KP.3']", "[np.float64(0.6), np.float64(0.4)]"),
    ]
    for lineages, abundances in forms:
        assert ra.parse_lineages(lineages, abundances) == (expect, None), (lineages, abundances)
    result, problem = ra.parse_lineages("A B C", "0.5 0.5")
    assert result == {} and "3 lineages but 2 abundances" in problem


def test_aggregated_file_with_all_forms():
    warnings = []
    table = ra.read_aggregated(FIX / "formats" / "aggregated_forms.tsv", warnings)
    for sample in ("space", "numpy", "pylist", "array"):
        assert table[sample]["lineages"] == {"JN.1": 0.6, "KP.3": 0.4}, sample
        assert table[sample]["summarized"] == {"JN.1* (JN.1.X)": 0.6, "Other": 0.4}, sample
    assert table["empty"] == {"lineages": {}, "summarized": {}}
    assert table["nan"] == {"lineages": {}, "summarized": {}}
    assert table["mismatch"]["lineages"] == {} and table["mismatch"]["summarized"] == {
        "Omicron": 1.0
    }
    assert any("mismatch" in w and "3 lineages but 2 abundances" in w for w in warnings)


def test_real_freyja_output_with_quoted_multiline_cells():
    """tests/fixtures/report/freyja_real_aggregated.tsv is raw Freyja output: numpy reprs wrapped over lines."""
    path = REPO / "tests" / "fixtures" / "report" / "freyja_real_aggregated.tsv"
    table = ra.read_aggregated(path, [])
    assert set(table) == {"mixture", "test"}  # the `.tsv` file-name label is cleaned
    assert len(table["mixture"]["lineages"]) == 10
    assert table["mixture"]["lineages"]["A"] == pytest.approx(0.366636)
    assert table["mixture"]["summarized"]["Alpha"] == pytest.approx(0.3331250000063474)


def test_label_cleaning():
    assert ra.clean_label("S1.variants.tsv") == "S1"
    assert ra.clean_label("S_1.freyja_variants.tsv") == "S_1"
    assert ra.clean_label("S1.demix") == "S1"
    assert ra.clean_label("S_1-b.2") == "S_1-b.2"  # sample IDs are never rewritten


# ----------------------------------------------------------------------------------------------- diff maths


def diff_dirs(tmp_path, **kwargs):
    old = make_run(tmp_path, "old", summary="old", samples=("A",))
    new = make_run(tmp_path, "new", summary="new", samples=("A",), run_info=info("demix"))
    return ra.diff_results(
        ra.Results(old), ra.Results(new), kwargs.pop("threshold", 0.05), **kwargs
    )


def row_of(result, sample, level, item):
    [row] = [
        r for r in result["rows"] if (r["sample"], r["level"], r["item"]) == (sample, level, item)
    ]
    return row


def test_change_maths_and_threshold_boundary(tmp_path):
    result = diff_dirs(tmp_path)
    jn1 = row_of(result, "A", "lineages", "JN.1")
    assert jn1["old_abundance"] == 0.6 and jn1["new_abundance"] == 0.55
    assert jn1["change"] == pytest.approx(-0.05)
    assert jn1["flagged"]  # 0.55 - 0.60 is -0.04999999999999993: still "at the threshold"
    assert row_of(result, "A", "lineages", "KP.3")["flagged"]  # +0.04999999999999999
    assert not row_of(result, "A", "lineages", "XFG.1")["flagged"]  # unchanged
    assert not row_of(result, "B", "lineages", "BA.2")["flagged"]  # +0.01
    just_above = diff_dirs(tmp_path / "x", threshold=0.0500001)
    assert not row_of(just_above, "A", "lineages", "JN.1")["flagged"]
    zero = diff_dirs(tmp_path / "y", threshold=0.0)
    assert row_of(zero, "A", "lineages", "XFG.1")["flagged"]  # |0| >= 0


def test_small_items_in_both_runs_are_ignored(tmp_path):
    result = diff_dirs(tmp_path)
    items = {(r["sample"], r["item"]) for r in result["rows"]}
    assert ("D", "ZZ.1") not in items  # 0.0005 -> 0.0009: below 0.001 in both
    assert ("D", "JN.1") in items
    assert row_of(result, "D", "lineages", "JN.1")["change"] == pytest.approx(-0.0004)
    assert not row_of(result, "D", "lineages", "JN.1")["flagged"]


def test_appearing_disappearing_samples_and_status_changes(tmp_path):
    result = diff_dirs(tmp_path)
    assert result["appear"] == ["E"] and result["disappear"] == ["C"]
    assert result["compared"] == ["A", "B", "D"]
    assert result["status_changes"] == [
        ("C", "ok", "skipped_low_coverage"),
        ("E", "absent", "ok"),
    ]
    c = row_of(result, "C", "lineages", "JN.1")
    assert c["new_abundance"] == 0.0 and c["change"] == -1.0 and c["flagged"]
    assert c["note"] == "sample missing in the new run; status ok -> skipped_low_coverage"
    e = row_of(result, "E", "lineages", "XFG.1")
    assert e["old_abundance"] == 0.0 and e["note"].startswith("sample only in the new run")
    # E has no status in the old run
    assert ("E", "absent", "ok") in result["status_changes"]
    assert result["flagged_samples"] == ["A", "C", "E"]
    assert result["old"]["barcode"] == OLD_ISO and result["new"]["barcode"] == "2026-10-20T04-10"
    assert result["old"]["freyja"] == "2.0.6" and result["new"]["freyja"] == "2.0.7"


def test_level_selection(tmp_path):
    only_lineages = diff_dirs(tmp_path / "a", level="lineages")
    assert {r["level"] for r in only_lineages["rows"]} == {"lineages"}
    summarized = diff_dirs(tmp_path / "b", level="summarized")
    assert {r["level"] for r in summarized["rows"]} == {"summarized"}
    kp3 = row_of(summarized, "A", "summarized", "KP.3* (KP.3.X)")
    assert kp3["old_abundance"] == 0.0 and kp3["new_abundance"] == 0.35 and kp3["flagged"]
    groups = diff_dirs(tmp_path / "c", level="groups")
    assert {r["level"] for r in groups["rows"]} == {"groups"}
    default = diff_dirs(tmp_path / "d")
    assert default["levels"] == ["lineages", "summarized", "groups"]


def test_group_columns_that_differ_are_reported_and_count_as_zero(tmp_path):
    result = diff_dirs(tmp_path, level="groups")
    assert any(
        "only in old: ['other']" not in n and "only in new: ['XFG*']" in n for n in result["notes"]
    )
    assert row_of(result, "A", "groups", "XFG*")["old_abundance"] == 0.0
    assert row_of(result, "A", "groups", "XFG*")["new_abundance"] == 0.1
    assert row_of(result, "A", "groups", "other")["change"] == pytest.approx(-0.1)
    assert row_of(result, "A", "groups", "JN.1*")["change"] == pytest.approx(0.0)


def test_groups_level_needs_the_summary_in_both_runs(tmp_path):
    old = make_run(tmp_path, "old", summary="old")
    new = make_run(tmp_path, "new", summary="new")
    (new / "summary" / "lineage_summary.csv").unlink()
    with pytest.raises(ra.Problem, match="needs summary/lineage_summary.csv in both"):
        ra.diff_results(ra.Results(old), ra.Results(new), 0.05, "groups")
    result = ra.diff_results(ra.Results(old), ra.Results(new), 0.05, "all")
    assert result["levels"] == ["lineages", "summarized"]
    assert any("group level skipped" in n for n in result["notes"])


def test_diff_of_a_run_with_itself_has_no_change(tmp_path):
    run = make_run(tmp_path, "run1")
    same = ra.diff_results(ra.Results(run), ra.Results(run), 0.05)
    assert same["flagged_samples"] == [] and same["top"] == []
    assert same["appear"] == same["disappear"] == [] and same["status_changes"] == []
    assert all(r["change"] == 0 for r in same["rows"])
    assert any("differences are not expected" in n for n in same["notes"])


def test_real_format_self_diff_and_one_edit(tmp_path):
    """Real Freyja output (raw, multi-line cells) and a real lineage_summary.csv from the repo's report fixtures."""
    real = REPO / "tests" / "fixtures" / "report"
    old = tmp_path / "old" / "summary"
    old.mkdir(parents=True)
    shutil.copy(real / "freyja_real_aggregated.tsv", old / "freyja_aggregated.tsv")
    shutil.copy(real / "h_lineage_groups" / "real_lineage_summary.csv", old / "lineage_summary.csv")
    new = tmp_path / "new"
    shutil.copytree(tmp_path / "old", new)
    same = ra.diff_results(ra.Results(tmp_path / "old"), ra.Results(new), 0.05)
    assert same["levels"] == ["lineages", "summarized", "groups"]
    assert same["flagged_samples"] == [] and len(same["rows"]) > 20
    text = (new / "summary" / "lineage_summary.csv").read_text()
    assert "0.908068" in text
    (new / "summary" / "lineage_summary.csv").write_text(text.replace("0.908068", "0.808068", 1))
    edited = ra.diff_results(ra.Results(tmp_path / "old"), ra.Results(new), 0.05, "groups")
    [top] = edited["top"]
    assert top["sample"] == "WW-1274E" and top["item"].startswith("NB.1.8.1*")
    assert top["change"] == pytest.approx(-0.1)


def test_csv_and_markdown_outputs(tmp_path, capsys):
    old = make_run(tmp_path, "old", summary="old")
    new = make_run(tmp_path, "new", summary="new", run_info=info("demix"))
    report = tmp_path / "out" / "diff.csv"
    md = tmp_path / "diff.md"
    rc = main("--diff", old, new, "--report", report, "--summary-md", md)
    assert rc == 0
    out = capsys.readouterr().out
    assert f"barcode version: {OLD_ISO} -> 2026-10-20T04-10" in out
    assert "Freyja version:  2.0.6 -> 2.0.7" in out
    assert "samples compared: 3" in out and "samples only in the new run: E" in out
    assert "samples missing in the new run: C" in out and "largest changes" in out
    with open(report) as handle:
        reader = csv.DictReader(handle)
        rows = list(reader)
    assert reader.fieldnames == [
        "sample",
        "level",
        "item",
        "old_abundance",
        "new_abundance",
        "change",
        "flagged",
        "note",
    ]
    jn1 = next(r for r in rows if (r["sample"], r["level"], r["item"]) == ("A", "lineages", "JN.1"))
    assert jn1["old_abundance"] == "0.600000" and jn1["new_abundance"] == "0.550000"
    assert jn1["change"] == "-0.050000" and jn1["flagged"] == "yes"
    assert [r["sample"] for r in rows] == sorted(r["sample"] for r in rows)
    text = md.read_text()
    assert "`2.0.6` -> `2.0.7`" in text and "| A | " in text and "Only in the new run: E" in text


def test_diff_output_is_deterministic(tmp_path):
    old = make_run(tmp_path, "old", summary="old")
    new = make_run(tmp_path, "new", summary="new", run_info=info("demix"))
    first, second = tmp_path / "1.csv", tmp_path / "2.csv"
    assert main("--diff", old, new, "--report", first) == 0
    assert main("--diff", old, new, "--report", second) == 0
    assert first.read_bytes() == second.read_bytes()


def test_diff_levels_via_cli_and_threshold(tmp_path, capsys):
    old = make_run(tmp_path, "old", summary="old")
    new = make_run(tmp_path, "new", summary="new", run_info=info("demix"))
    assert main("--diff", old, new, "--level", "summarized", "--threshold", "0.2") == 0
    out = capsys.readouterr().out
    assert "levels: summarized" in out and "lineages " not in out.split("largest changes", 1)[-1]
    assert main("--diff", old, new, "--level", "groups") == 0


def test_barcode_info_falls_back_to_barcode_version_json(tmp_path):
    run = make_run(tmp_path, "r", summary="old")
    (run / "pipeline_info" / "run_info.json").unlink()
    (run / "summary" / "barcode_version.json").write_text(
        json.dumps(
            {
                "barcode_version": "10_06_2026-05-25",
                "barcode_version_iso": OLD_ISO,
                "freyja_version": "2.0.6",
            }
        )
    )
    versions = ra.Results(run).versions
    assert versions["barcode"] == OLD_ISO and versions["freyja"] == "2.0.6"
    (run / "summary" / "barcode_version.json").unlink()
    assert ra.Results(run).versions["barcode"] == "unknown"


# ----------------------------------------------------------------------------------------------- empty / broken inputs


def test_missing_aggregate_is_a_friendly_error(tmp_path, capsys):
    old = make_run(tmp_path, "old", summary="old")
    bare = make_run(tmp_path, "bare", summary=None)
    assert main("--diff", old, bare) == 1
    err = capsys.readouterr().err
    assert "no summary/freyja_aggregated.tsv" in err and "Traceback" not in err
    assert main("--diff", old, tmp_path / "gone") == 1
    assert "not a directory" in capsys.readouterr().err


def test_empty_and_header_only_tables(tmp_path, capsys):
    old = make_run(tmp_path, "old", summary="old")
    new = make_run(tmp_path, "new", summary="new")
    (new / "summary" / "freyja_aggregated.tsv").write_text("")
    result = ra.diff_results(ra.Results(old), ra.Results(new), 0.05, "lineages")
    assert result["compared"] == [] and result["disappear"] == ["A", "B", "C", "D"]
    (new / "summary" / "freyja_aggregated.tsv").write_text("\tsummarized\tlineages\tabundances\n")
    (new / "summary" / "lineage_summary.csv").write_text("")
    result = ra.diff_results(ra.Results(old), ra.Results(new), 0.05)
    assert result["appear"] == [] and all(r["new_abundance"] == 0.0 for r in result["rows"])
    assert main("--diff", old, new) == 0
    assert "samples compared: 0" in capsys.readouterr().out


def test_empty_both_sides_gives_an_empty_report(tmp_path):
    old = make_run(tmp_path, "old", summary="old")
    new = make_run(tmp_path, "new", summary="new")
    for run in (old, new):
        (run / "summary" / "freyja_aggregated.tsv").write_text("")
        (run / "summary" / "lineage_summary.csv").unlink()
    report = tmp_path / "e.csv"
    assert main("--diff", old, new, "--report", report) == 0
    assert report.read_text().splitlines() == [
        "sample,level,item,old_abundance,new_abundance,change,flagged,note"
    ]


def test_argument_errors(tmp_path):
    run = make_run(tmp_path, "run1")
    for argv in (
        [],
        ["--diff", str(run), str(run), str(run)],
        ["--diff", str(run), str(run), str(run), str(run)],
    ):
        with pytest.raises(SystemExit):
            main(*argv)
    with pytest.raises(SystemExit):
        main(run, "--diff", run, run)
    with pytest.raises(SystemExit):
        main(run, "--threshold", "-1")


# ----------------------------------------------------------------------------------------------- script


def test_script_help_and_dry_run_as_a_subprocess(tmp_path):
    run = make_run(tmp_path, "run1")
    script = REPO / "bin" / "reanalyse_runs.py"
    help_text = subprocess.run(
        [sys.executable, str(script), "--help"], check=False, capture_output=True, text=True
    )
    assert help_text.returncode == 0 and "--execute" in help_text.stdout
    result = subprocess.run(
        [sys.executable, str(script), str(run), "--barcodes", "2026-10-20T04-10"],
        check=False,
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0 and "NOT executed" in result.stdout
    # default pipeline dir: the clone the script lives in
    assert [c[2] for c in printed_commands(result.stdout)] == [str(REPO)]


def test_runs_on_the_host_python():
    name = "reanalyse_runs.py"
    ast.parse((compat.BIN / name).read_text(), filename=name, feature_version=(3, 8))
    tree = ast.parse((compat.BIN / name).read_text(), filename=name)
    assert any(
        isinstance(node, ast.ImportFrom) and node.module == "__future__" for node in ast.walk(tree)
    )
    lazy = compat.annotation_nodes(tree)
    for node in ast.walk(tree):
        where = f"{name}:{getattr(node, 'lineno', '?')}"
        if isinstance(node, ast.Call) and getattr(node.func, "id", "") == "zip":
            assert not any(k.arg == "strict" for k in node.keywords), f"{where}: zip(strict=)"
        if isinstance(node, ast.Attribute):
            assert node.attr not in compat.NEWER_ATTRIBUTES, (
                f"{where}: .{node.attr} needs Python >= 3.9"
            )
        if isinstance(node, ast.BinOp) and isinstance(node.op, ast.BitOr) and id(node) not in lazy:
            operands = (node.left, node.right)
            is_type_union = any(
                isinstance(o, (ast.Subscript, ast.Name, ast.Attribute))
                or (isinstance(o, ast.Constant) and o.value is None)
                for o in operands
            ) and not any(
                isinstance(o, ast.Constant) and isinstance(o.value, int) for o in operands
            )
            names_only = all(isinstance(o, (ast.Subscript, ast.Constant)) for o in operands)
            assert not (is_type_union and names_only), f"{where}: `X | Y` outside an annotation"
    # also no dict-merge operator and no str.removeprefix/removesuffix: covered by NEWER_ATTRIBUTES; no match
    assert not any(type(n).__name__ == "Match" for n in ast.walk(tree))
