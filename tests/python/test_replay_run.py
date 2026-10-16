"""bin/replay_run.py: rebuild the Nextflow command of a recorded run.

The fixtures in tests/fixtures/replay/<mode>/ are copies of the real pipeline_info/run_info.json and params_*.json of
the synthetic end-to-end runs (tests/e2e.sh) with the paths replaced by @ROOT@ (scratch dir) and @REPO@ (checkout).
`build` writes one of them into a temporary world whose files match the recorded checksums, so tests can then break
exactly one thing.
"""

import hashlib
import json
import os
import shutil
import stat
import subprocess
from pathlib import Path

import pytest
import replay_run as rr
import test_host_python_compat as compat

FIX = Path(__file__).resolve().parents[1] / "fixtures" / "replay"
ISO = "2026-10-06T05-25"
STAMP = "10_06_2026-05-25"


@pytest.fixture(autouse=True)
def isolated_home(tmp_path, monkeypatch):
    """The default db_dir (~/.wasteline/db) must never be the user's real one."""
    monkeypatch.setenv("HOME", str(tmp_path / "home"))


def sha(text):
    return hashlib.sha256(text.encode() if isinstance(text, str) else text).hexdigest()


class World:
    """results dir + the files its run_info.json refers to, all under tmp_path."""

    def __init__(self, tmp_path):
        self.tmp = tmp_path
        self.root = tmp_path / "root"
        self.repo = tmp_path / "repo"
        self.results = self.root / "results"
        self.params_out = tmp_path / "replay_params.json"
        self.run_info = None
        self.params = None

    def run(self, *extra, execute=False):
        argv = [str(self.results), "--params-out", str(self.params_out), *extra]
        if execute:
            argv.append("--execute")
        return rr.main(argv)

    def written_params(self):
        return json.loads(self.params_out.read_text())


def build(tmp_path, name, mutate=None):
    """Write fixture `name` as results/pipeline_info/{run_info.json,params_*.json} plus every file it records."""
    world = World(tmp_path)

    def hydrate(text):
        return text.replace("@ROOT@", str(world.root)).replace("@REPO@", str(world.repo))

    run_info_path = FIX / name / "run_info.json"
    world.run_info = (
        json.loads(hydrate(run_info_path.read_text())) if run_info_path.exists() else None
    )
    world.params = json.loads(hydrate((FIX / name / "params.json").read_text()))

    # input files, with the recorded checksums made to match them
    sums = (world.run_info or {}).get("input_checksums", {})
    if "samplesheet" in sums and world.params.get("input"):
        sheet = Path(world.params["input"])
        reads = world.root / "reads"
        reads.mkdir(parents=True, exist_ok=True)
        (reads / "A_R1.fastq.gz").write_text("r1")
        (reads / "A_R2.fastq.gz").write_text("r2")
        sheet.write_text(
            f"sample,fastq_1,fastq_2\nA,{reads / 'A_R1.fastq.gz'},{reads / 'A_R2.fastq.gz'}\n"
        )
        sums["samplesheet"]["sha256"] = sha(sheet.read_bytes())
    if "primer_bed" in sums:
        bed = world.repo / "assets" / "primer_schemes" / sums["primer_bed"]["file"]
        bed.parent.mkdir(parents=True, exist_ok=True)
        bed.write_text("NC_045512.2\t30\t54\tmidnight_1_LEFT\t1\t+\n")
        sums["primer_bed"]["sha256"] = sha(bed.read_bytes())
    if "reference_fasta" in sums:
        fasta = world.repo / "assets" / "references" / sums["reference_fasta"]["file"]
        fasta.parent.mkdir(parents=True, exist_ok=True)
        fasta.write_text(">NC_045512.2\nACGT\n")
        sums["reference_fasta"]["sha256"] = sha(fasta.read_bytes())
    block = (world.run_info or {}).get("freyja_data") or {}
    if block.get("source") == "custom":
        barcodes = Path(block["dataset_path"])
        barcodes.parent.mkdir(parents=True, exist_ok=True)
        barcodes.write_text("lineage,A1T\nB.1,1\n")
        block["barcodes_sha256"] = sha(barcodes.read_bytes())
    # earlier results: their run_info.json must match the recorded sha256
    for record in (world.run_info or {}).get("previous_runs", []):
        previous = Path(record["path"])
        (previous / "pipeline_info").mkdir(parents=True, exist_ok=True)
        if record.get("run_info_sha256"):
            text = '{"made": "earlier"}'
            (previous / "pipeline_info" / "run_info.json").write_text(text)
            record["run_info_sha256"] = sha(text)
    for key in ("results_dir", "demix_dir"):
        if world.params.get(key):
            (Path(world.params[key]) / "pipeline_info").mkdir(parents=True, exist_ok=True)
            (Path(world.params[key]) / "pipeline_info" / "run_info.json").write_text("{}")
    if mutate:
        mutate(world)
    pipeline_info = world.results / "pipeline_info"
    pipeline_info.mkdir(parents=True)
    if world.run_info is not None:
        world.run_info["params"] = dict(
            world.params
        )  # run_info.json and the params file record the same values
        (pipeline_info / "run_info.json").write_text(json.dumps(world.run_info, indent=2))
    (pipeline_info / "params_2026-10-16_09-13-28.json").write_text(json.dumps(world.params))
    return world


def profile_of(cmd):
    """The -profile value of a printed command."""
    return cmd.split("-profile ", 1)[1].split()[0]


def command(capsys_out):
    """The printed command: the text after the 'Command' heading."""
    return capsys_out.split("Command", 1)[1].split(":\n", 1)[1]


def run_and_read(world, capsys, *extra, execute=False):
    rc = world.run(*extra, execute=execute)
    out = capsys.readouterr().out
    return rc, out


# ----------------------------------------------------------------------------------------------- full mode


def test_full_run_reconstructed(tmp_path, capsys):
    world = build(tmp_path, "full")
    rc, out = run_and_read(world, capsys)
    assert rc == 0, out
    cmd = command(out)
    # pinned to the recorded commit, from GitHub
    assert "NXF_VER=26.04.6 nextflow run charlesfoster/wasteline" in cmd
    assert f"-r {world.run_info['git_commit']}" in cmd
    assert "-profile docker" in cmd
    assert f"-params-file {world.params_out}" in cmd
    assert f"--outdir {world.params['outdir']}_replay" in cmd
    # custom barcodes: the same file (checksum verified), compare dataset from the container
    assert f"--barcodes {world.repo}/tests/data/barcodes_test.csv" in cmd
    assert "--compare_barcodes container" in cmd
    assert "OK" in out and "samplesheet" in out and "primer BED" in out and "reference" in out
    assert "Problems: 0 error(s)" in out
    params = world.written_params()
    for gone in (
        "outdir",
        "barcodes",
        "compare_barcodes",
        "db_dir",
        "previous_results",
        "pipelines_testdata_base_path",
    ):
        assert gone not in params
    assert "fasta" not in params and "gff" not in params  # the bundled files of the checkout
    assert params["input"] == world.params["input"]
    assert params["primer_scheme"] == "midnight"
    assert params["mode"] == "full" and params["freyja_boot"] is True
    assert all(value is not None for value in params.values())


def test_volatile_and_notification_params_are_not_replayed(tmp_path, capsys):
    def mutate(world):
        world.params.update(
            trace_report_suffix="2026-10-16_09-13-26",
            email="me@example.org",
            email_on_fail="me@example.org",
        )

    world = build(tmp_path, "full", mutate)
    rc, out = run_and_read(world, capsys)
    params = world.written_params()
    assert rc == 0
    assert (
        "trace_report_suffix" not in params
        and "email" not in params
        and "email_on_fail" not in params
    )
    assert "should not send the notification again" in out


def test_run_info_json_or_pipeline_info_path_is_accepted(tmp_path, capsys):
    world = build(tmp_path, "full")
    for target in (
        world.results / "pipeline_info",
        world.results / "pipeline_info" / "run_info.json",
    ):
        rc = rr.main([str(target), "--params-out", str(world.params_out)])
        assert rc == 0, capsys.readouterr().out
        capsys.readouterr()


def test_unreadable_record(tmp_path, capsys):
    with pytest.raises(SystemExit) as err:
        rr.main([str(tmp_path), "--params-out", str(tmp_path / "p.json")])
    assert "no pipeline_info/run_info.json" in str(err.value)


# ------------------------------------------------------------------------------------ pipeline checkout, schema


def make_checkout(world, defaults):
    """A git checkout of a fake pipeline with a nextflow_schema.json carrying `defaults`; returns its HEAD."""
    world.repo.mkdir(parents=True, exist_ok=True)
    schema = {"$defs": {"g": {"properties": {k: {"default": v} for k, v in defaults.items()}}}}
    (world.repo / "nextflow_schema.json").write_text(json.dumps(schema))
    env = dict(
        os.environ,
        GIT_AUTHOR_NAME="t",
        GIT_AUTHOR_EMAIL="t@t",
        GIT_COMMITTER_NAME="t",
        GIT_COMMITTER_EMAIL="t@t",
    )
    for step in (["init", "-q"], ["add", "."], ["commit", "-q", "-m", "x"]):
        subprocess.run(
            ["git", "-C", str(world.repo), *step], check=True, env=env, stdout=subprocess.DEVNULL
        )
    return subprocess.run(
        ["git", "-C", str(world.repo), "rev-parse", "HEAD"],
        check=True,
        stdout=subprocess.PIPE,
        text=True,
    ).stdout.strip()


needs_git = pytest.mark.skipif(shutil.which("git") is None, reason="git not installed")


@needs_git
def test_pipeline_dir_at_recorded_commit_and_schema_defaults(tmp_path, capsys):
    holder = {}

    def mutate(world):
        # the checkout holds the bundled files already (build wrote them); its schema says these are defaults
        holder["head"] = make_checkout(
            world,
            {"freyja_minq": 20, "consensus_min_depth": 15, "freyja_boot": False, "mode": "full"},
        )
        world.run_info["git_commit"] = holder["head"]

    world = build(tmp_path, "full", mutate)
    rc, out = run_and_read(world, capsys, "--pipeline-dir", str(world.repo))
    assert rc == 0, out
    cmd = command(out)
    assert f"nextflow run {world.repo}" in cmd and "-r " not in cmd
    assert f"git -C {world.repo} checkout {holder['head']}" in out
    params = world.written_params()
    assert (
        "freyja_minq" not in params and "consensus_min_depth" not in params and "mode" not in params
    )  # defaults
    assert params["freyja_boot"] is True  # differs from the default
    assert "freyja_eps" in params  # not in the (fake) schema at all: kept, with a warning
    assert "is not a parameter of" in out


@needs_git
def test_pipeline_dir_at_other_commit_is_refused_unless_forced(tmp_path, capsys):
    def mutate(world):
        make_checkout(world, {})

    world = build(tmp_path, "full", mutate)
    rc, out = run_and_read(world, capsys, "--pipeline-dir", str(world.repo))
    assert rc == 2
    assert "but the run used 096a955a178c" in out and "NOT READY" in out
    rc, out = run_and_read(world, capsys, "--pipeline-dir", str(world.repo), "--force")
    assert rc == 0 and "overridden by --force" in out


def test_null_commit_warns_and_uses_the_version(tmp_path, capsys):
    def mutate(world):
        world.run_info["git_commit"] = None
        world.run_info["pipeline_version"] = "1.2.0"

    world = build(tmp_path, "full", mutate)
    rc, out = run_and_read(world, capsys)
    assert rc == 0
    assert "git_commit is null" in out and "cannot be pinned" in out
    assert "-r 1.2.0" in command(out)


def test_null_commit_on_a_dev_version(tmp_path, capsys):
    def mutate(world):
        world.run_info["git_commit"] = None

    world = build(tmp_path, "full", mutate)
    rc, out = run_and_read(world, capsys)
    assert rc == 0
    assert "git_commit is null" in out and "development version" in out
    assert "-r 0.1.0dev" in command(out)


# ---------------------------------------------------------------------------------------- barcodes sources


def stored_dataset(world, db, text="barcodes"):
    directory = db / "freyja" / ISO
    directory.mkdir(parents=True)
    (directory / "usher_barcodes.feather").write_text(text)
    return directory


def latest_block(world, source="latest", **override):
    block = {
        "source": source,
        "barcode_version": STAMP,
        "barcode_version_iso": ISO,
        "barcodes_file": "usher_barcodes.feather",
        "barcodes_sha256": sha("barcodes"),
        "barcodes_format": "feather",
        "freyja_version": "2.0.6",
        "freyja_reported_barcode_version": None,
        "freyja_bundled_barcode_version": "01_01_2026-00-00",
        "dataset_path": str(world.tmp / "db_recorded" / "freyja"),
        "downloaded_at": "2026-10-06T06:00:00Z",
    }
    block.update(override)
    return block


@pytest.mark.parametrize("source", ["latest", "pinned"])
def test_latest_or_pinned_uses_the_stored_dataset(tmp_path, capsys, source):
    def mutate(world):
        world.run_info["freyja_data"] = latest_block(world, source)
        world.params["compare_barcodes"] = None

    world = build(tmp_path, "full", mutate)
    db = tmp_path / "mydb"
    stored_dataset(world, db)
    rc, out = run_and_read(world, capsys, "--db_dir", str(db))
    assert rc == 0, out
    cmd = command(out)
    assert f"--barcodes {ISO}" in cmd and f"--db_dir {db}" in cmd
    assert "sha256 checked" in out


def test_stored_dataset_found_in_the_recorded_db_dir(tmp_path, capsys):
    def mutate(world):
        world.run_info["freyja_data"] = latest_block(world)
        world.params["db_dir"] = str(world.tmp / "recorded_db")
        world.params["compare_barcodes"] = None

    world = build(tmp_path, "full", mutate)
    stored_dataset(world, tmp_path / "recorded_db")
    rc, out = run_and_read(world, capsys)
    assert rc == 0
    assert f"--barcodes {ISO} " in command(
        out
    ) and f"--db_dir {tmp_path / 'recorded_db'}" in command(out)


def test_stored_dataset_with_a_different_checksum_is_an_error(tmp_path, capsys):
    def mutate(world):
        world.run_info["freyja_data"] = latest_block(world)
        world.params["compare_barcodes"] = None

    world = build(tmp_path, "full", mutate)
    stored_dataset(world, tmp_path / "db", text="tampered")
    rc, out = run_and_read(world, capsys, "--db_dir", str(tmp_path / "db"))
    assert rc == 2 and "has sha256" in out
    rc, out = run_and_read(world, capsys, "--db_dir", str(tmp_path / "db"), "--force")
    assert rc == 0


def test_dataset_found_through_the_recorded_dataset_path(tmp_path, capsys):
    def mutate(world):
        world.run_info["freyja_data"] = latest_block(world)
        world.params["compare_barcodes"] = None

    world = build(tmp_path, "full", mutate)
    stored_dataset(
        world, tmp_path / "db_recorded"
    )  # <dataset_path>/<iso>, but not found as a db_dir by --db_dir
    rc, out = run_and_read(world, capsys, "--db_dir", str(tmp_path / "nowhere"))
    assert rc == 0, out
    assert f"--barcodes {tmp_path / 'db_recorded' / 'freyja' / ISO}" in command(out)
    assert "recorded dataset_path" in out


def test_missing_dataset_must_be_restored(tmp_path, capsys, monkeypatch):
    def mutate(world):
        world.run_info["freyja_data"] = latest_block(world)
        world.params["compare_barcodes"] = None

    monkeypatch.setenv("HOME", str(tmp_path / "home"))  # the default db_dir must not exist
    world = build(tmp_path, "full", mutate)
    rc, out = run_and_read(world, capsys)
    assert rc == 2
    assert (
        "MUST BE RESTORED" in out
        and ISO in out
        and "cannot be downloaded again" not in out.split("Command")[1]
    )
    rc, out = run_and_read(world, capsys, "--force")
    assert rc == 2  # a missing dataset is not something --force can fix


def test_missing_dataset_but_the_container_bundles_the_same_barcodes(tmp_path, capsys, monkeypatch):
    def mutate(world):
        world.run_info["freyja_data"] = latest_block(world, freyja_bundled_barcode_version=STAMP)
        world.params["compare_barcodes"] = None

    monkeypatch.setenv("HOME", str(tmp_path / "home"))
    world = build(tmp_path, "full", mutate)
    rc, out = run_and_read(world, capsys)
    assert rc == 0, out
    assert "--barcodes container" in command(out) and "bundles the same barcodes" in out


def test_container_source(tmp_path, capsys):
    def mutate(world):
        world.run_info["freyja_data"] = latest_block(
            world, "container", dataset_path="bundled in the Freyja container"
        )
        world.params["compare_barcodes"] = None

    world = build(tmp_path, "full", mutate)
    rc, out = run_and_read(world, capsys)
    assert rc == 0
    assert "--barcodes container" in command(out) and "--db_dir" not in command(out)


def test_custom_barcodes_with_a_changed_file_are_an_error(tmp_path, capsys):
    world = build(tmp_path, "full")
    Path(world.run_info["freyja_data"]["dataset_path"]).write_text("edited since the run")
    rc, out = run_and_read(world, capsys)
    assert rc == 2 and "has sha256" in out
    assert "--barcodes" in command(out)  # the command is still printed, marked as not ready


def test_custom_barcodes_missing_is_an_error_even_with_force(tmp_path, capsys):
    def mutate(world):
        Path(world.run_info["freyja_data"]["dataset_path"]).unlink()

    world = build(tmp_path, "full", mutate)
    rc, out = run_and_read(world, capsys, "--force")
    assert rc == 2 and "not found: restore it" in out


def test_custom_dataset_directory(tmp_path, capsys):
    def mutate(world):
        directory = world.tmp / "dataset"
        directory.mkdir()
        (directory / "usher_barcodes.feather").write_text("barcodes")
        world.run_info["freyja_data"] = latest_block(
            world,
            "custom",
            dataset_path=str(directory),
            barcode_version="custom",
            barcode_version_iso=None,
        )
        world.params["compare_barcodes"] = None

    world = build(tmp_path, "full", mutate)
    rc, out = run_and_read(world, capsys)
    assert rc == 0, out
    assert f"--barcodes {tmp_path / 'dataset'}" in command(out)


def test_compare_barcodes_recorded_but_not_described(tmp_path, capsys):
    def mutate(world):
        del world.run_info["freyja_data_compare"]
        world.params["compare_barcodes"] = "latest"

    world = build(tmp_path, "full", mutate)
    _rc, out = run_and_read(world, capsys)
    assert "no freyja_data_compare is recorded" in out
    assert "--compare_barcodes latest" in command(out)


# ------------------------------------------------------------------------------------------ input checks


def test_changed_samplesheet_is_named(tmp_path, capsys):
    world = build(tmp_path, "full")
    Path(world.params["input"]).write_text("sample,fastq_1,fastq_2\nB,x,y\n")
    rc, out = run_and_read(world, capsys)
    assert rc == 2
    assert "samplesheet" in out and "DIFFERS from the run" in out
    assert "NOT READY" in out


def test_changed_primer_bed_and_reference_are_named(tmp_path, capsys):
    world = build(tmp_path, "full")
    (world.repo / "assets" / "primer_schemes" / "midnight.bed").write_text("changed\n")
    (world.repo / "assets" / "references" / "NC_045512.fasta").write_text(">x\nA\n")
    rc, out = run_and_read(world, capsys, "--pipeline-dir", str(world.repo))
    assert rc == 2
    assert "primer BED" in out and "reference" in out
    assert out.count("DIFFERS from the run") == 2


def test_missing_samplesheet_is_never_forced(tmp_path, capsys):
    world = build(tmp_path, "full")
    Path(world.params["input"]).unlink()
    for extra in ((), ("--force",)):
        rc, out = run_and_read(world, capsys, *extra)
        assert rc == 2 and "samplesheet: " in out and "not found" in out


def test_missing_fastq_file_is_an_error(tmp_path, capsys):
    world = build(tmp_path, "full")
    sheet = Path(world.params["input"])
    sheet.write_text("sample,fastq_1,fastq_2\nA,/nowhere/A_R1.fastq.gz,/nowhere/A_R2.fastq.gz\n")
    world.run_info["input_checksums"]["samplesheet"]["sha256"] = sha(sheet.read_bytes())
    (world.results / "pipeline_info" / "run_info.json").write_text(json.dumps(world.run_info))
    rc, out = run_and_read(world, capsys)
    assert rc == 2 and "FASTQ file(s) of the samplesheet not found" in out


# ----------------------------------------------------------------------------------------- other modes


def test_demix_reproduces_previous_results(tmp_path, capsys):
    world = build(tmp_path, "demix")
    rc, out = run_and_read(world, capsys)
    assert rc == 0, out
    cmd = command(out)
    assert f"--previous_results {world.root}/out_full" in cmd
    assert "previous_results" not in world.written_params()
    assert "earlier results" in out and "OK" in out
    assert world.written_params()["mode"] == "demix"


def test_demix_with_changed_previous_results_is_an_error(tmp_path, capsys):
    world = build(tmp_path, "demix")
    (world.root / "out_full" / "pipeline_info" / "run_info.json").write_text('{"made": "later"}')
    rc, out = run_and_read(world, capsys)
    assert rc == 2 and "changed: run_info.json sha256" in out


def test_demix_with_missing_previous_results_is_never_forced(tmp_path, capsys):
    world = build(tmp_path, "demix")
    shutil.rmtree(world.root / "out_full")
    rc, out = run_and_read(world, capsys, "--force")
    assert rc == 2 and "earlier results" in out and "not found" in out


def test_demix_of_legacy_results_has_no_checksum_to_compare(tmp_path, capsys):
    world = build(tmp_path, "demix_legacy")
    rc, out = run_and_read(world, capsys)
    assert rc == 0, out
    assert "legacy results: no checksum recorded" in out
    assert "--previous_results" in command(out) and "legacy_results" in command(out)


def test_aggregate_with_previous_results(tmp_path, capsys):
    world = build(tmp_path, "aggregate")
    rc, out = run_and_read(world, capsys)
    assert rc == 0, out
    assert f"--previous_results {world.root}/out_full" in command(out)
    assert "--barcodes container" in command(out)
    assert "only groups the plot" in out


def test_aggregate_with_demix_dir_keeps_the_directory_in_the_params(tmp_path, capsys):
    world = build(tmp_path, "aggregate_dir")
    rc, out = run_and_read(world, capsys)
    assert rc == 0, out
    assert world.written_params()["demix_dir"] == world.params["demix_dir"]
    assert "--previous_results" not in command(out)
    assert "demix_dir" in out and "present" in out
    shutil.rmtree(world.params["demix_dir"])
    rc, out = run_and_read(world, capsys)
    assert rc == 2 and "demix_dir" in out and "not found" in out


def test_report_run_has_no_run_info(tmp_path, capsys):
    world = build(tmp_path, "report")
    assert not (world.results / "pipeline_info" / "run_info.json").exists()
    rc, out = run_and_read(world, capsys)
    assert rc == 0, out
    cmd = command(out)
    assert f"--results_dir {world.root}/out_full" in cmd
    assert "--barcodes" not in cmd and "-r " not in cmd
    assert "no run_info.json" in out and "not used in mode report" in out
    assert "the -profile is not recorded" in out
    assert "results_dir" not in world.written_params()


def test_report_of_missing_results(tmp_path, capsys):
    world = build(tmp_path, "report")
    shutil.rmtree(world.params["results_dir"])
    rc, out = run_and_read(world, capsys)
    assert rc == 2 and "results_dir" in out and "not found" in out


def test_capture_run_has_no_primer_bed_to_check(tmp_path, capsys):
    world = build(tmp_path, "capture")
    rc, out = run_and_read(world, capsys)
    assert rc == 0, out
    assert "primer BED" not in out
    assert world.written_params()["protocol"] == "capture"


# ------------------------------------------------------------------------------------ profile and outdir


def test_digest_pinned_containers_add_the_locked_profile(tmp_path, capsys):
    def mutate(world):
        world.run_info["containers"] = {
            k: v.split(":")[0] + "@sha256:" + "a" * 64
            for k, v in world.run_info["containers"].items()
        }

    world = build(tmp_path, "full", mutate)
    _rc, out = run_and_read(world, capsys)
    assert "-profile docker,locked" in command(out) and "adding the locked profile" in out
    _rc, out = run_and_read(world, capsys, "--no-locked")
    assert profile_of(command(out)) == "docker"


def test_recorded_locked_profile_can_be_dropped(tmp_path, capsys):
    def mutate(world):
        world.run_info["profile"] = "docker,locked"

    world = build(tmp_path, "full", mutate)
    _rc, out = run_and_read(world, capsys)
    assert "-profile docker,locked" in command(out)
    _rc, out = run_and_read(world, capsys, "--no-locked")
    assert profile_of(command(out)) == "docker"


def test_tag_pinned_containers_get_no_locked_profile(tmp_path, capsys):
    world = build(tmp_path, "full")
    _rc, out = run_and_read(world, capsys)
    assert profile_of(command(out)) == "docker"


def test_outdir_is_never_the_original(tmp_path, capsys):
    world = build(tmp_path, "full")
    original = world.params["outdir"]
    for extra in (
        ("--outdir", original),
        ("--outdir", original + "/", "--force"),
        ("--outdir", str(world.results)),
    ):
        rc, out = run_and_read(world, capsys, *extra)
        assert rc == 2 and "is the original results directory" in out
    rc, out = run_and_read(world, capsys, "--outdir", str(tmp_path / "new"))
    assert rc == 0 and f"--outdir {tmp_path / 'new'}" in command(out)


# ----------------------------------------------------------------------------------------------- execute


@pytest.fixture
def fake_nextflow(tmp_path, monkeypatch):
    """An executable that records its arguments and NXF_VER instead of running a pipeline."""
    log = tmp_path / "nextflow_calls.txt"
    script = tmp_path / "fake_nextflow"
    script.write_text(
        '#!/bin/sh\n{ echo "NXF_VER=$NXF_VER"; printf "%s\\n" "$@"; } > "$FAKE_NEXTFLOW_LOG"\nexit 7\n'
    )
    script.chmod(script.stat().st_mode | stat.S_IXUSR)
    monkeypatch.setenv("FAKE_NEXTFLOW_LOG", str(log))
    return script, log


def test_execute_prints_the_command_then_runs_it(tmp_path, capsys, fake_nextflow):
    script, log = fake_nextflow
    world = build(tmp_path, "full")
    rc = world.run("--nextflow", str(script), execute=True)
    captured = capsys.readouterr()
    assert rc == 7  # Nextflow's exit status is passed on
    assert "Command" in captured.out and "Executing" in captured.err
    lines = log.read_text().splitlines()
    assert lines[0] == "NXF_VER=26.04.6"
    assert lines[1:4] == ["run", "charlesfoster/wasteline", "-r"]
    assert "-params-file" in lines and "--outdir" in lines


def test_execute_is_refused_after_a_checksum_mismatch(tmp_path, capsys, fake_nextflow):
    script, log = fake_nextflow
    world = build(tmp_path, "full")
    Path(world.params["input"]).write_text("changed")
    rc = world.run("--nextflow", str(script), execute=True)
    captured = capsys.readouterr()
    assert rc == 2 and not log.exists()
    assert "Not executing" in captured.err
    assert "nextflow run" not in captured.err  # the command was printed (stdout), not run
    # --force overrides the checksum mismatch, nothing else
    assert world.run("--nextflow", str(script), "--force", execute=True) == 7
    assert log.exists()


def test_execute_is_refused_for_a_missing_samplesheet_even_with_force(
    tmp_path, capsys, fake_nextflow
):
    script, log = fake_nextflow
    world = build(tmp_path, "full")
    Path(world.params["input"]).unlink()
    assert world.run("--nextflow", str(script), "--force", execute=True) == 2
    assert not log.exists()


def test_no_nxf_ver(tmp_path, capsys):
    world = build(tmp_path, "full")
    _rc, out = run_and_read(world, capsys, "--no-nxf-ver")
    assert "NXF_VER" not in out


# ------------------------------------------------------------------------------------------- python 3.8


def test_replay_run_stays_valid_on_the_hosts_python(monkeypatch):
    monkeypatch.setattr(compat, "HOST_SCRIPTS", ["replay_run.py"])
    compat.test_host_scripts_parse_as_python_3_8()
    compat.test_host_scripts_use_no_post_3_8_runtime_features()


def test_uses_only_the_standard_library():
    import ast

    tree = ast.parse((compat.BIN / "replay_run.py").read_text())
    modules = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            modules.update(alias.name.split(".")[0] for alias in node.names)
        if isinstance(node, ast.ImportFrom) and node.module:
            modules.add(node.module.split(".")[0])
    assert modules <= {
        "__future__",
        "argparse",
        "csv",
        "hashlib",
        "json",
        "os",
        "shlex",
        "subprocess",
        "sys",
        "pathlib",
    }, modules


def test_same_value_does_not_confuse_booleans_and_numbers():
    assert rr.same_value(True, True) and not rr.same_value(True, 1) and not rr.same_value(0, False)
    assert rr.same_value(20, 20.0) and not rr.same_value("a", None)
