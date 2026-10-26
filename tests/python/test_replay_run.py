"""bin/replay_run.py: rebuild the Nextflow command of a recorded run.

The fixtures in tests/fixtures/replay/full/ are copies of the real pipeline_info/run_info.json and params_*.json of
the synthetic end-to-end runs (tests/e2e.sh) with the paths replaced by @ROOT@ (scratch dir) and @REPO@ (checkout).
`build` writes one of them into a temporary world whose files match the recorded checksums, so tests can then break
exactly one thing.
"""

import hashlib
import json
import stat
from pathlib import Path

import pytest
import replay_run as rr

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


def command(capsys_out):
    """The printed command: the text after the 'Command' heading."""
    return capsys_out.split("Command", 1)[1].split(":\n", 1)[1]


def run_and_read(world, capsys, *extra, execute=False):
    rc = world.run(*extra, execute=execute)
    out = capsys.readouterr().out
    return rc, out


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


def test_null_commit_warns_and_uses_the_version(tmp_path, capsys):
    def mutate(world):
        world.run_info["git_commit"] = None
        world.run_info["pipeline_version"] = "1.2.0"

    world = build(tmp_path, "full", mutate)
    rc, out = run_and_read(world, capsys)
    assert rc == 0
    assert "git_commit is null" in out and "cannot be pinned" in out
    assert "-r 1.2.0" in command(out)


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


def test_panel_results_are_replayed_as_one_run_not_per_target(tmp_path, capsys):
    """Results of a --panel run keep one run_info.json per target (same command line and params): one replay."""
    world = build(tmp_path, "full")
    pipeline_info = world.results / "pipeline_info"
    text = (pipeline_info / "run_info.json").read_text()
    (pipeline_info / "run_info.json").unlink()
    for target in ("TA", "TB"):
        (pipeline_info / target).mkdir()
        (pipeline_info / target / "run_info.json").write_text(text)
    rc, out = run_and_read(world, capsys)
    assert rc == 0, out
    assert "panel run: 2 target(s) (TA, TB), one replay of the whole run" in out
    assert command(out).count("nextflow run") == 1
    assert f"-r {world.run_info['git_commit']}" in command(out)
    assert "Problems: 0 error(s)" in out
