"""bin/reanalyse_runs.py: re-run past results with new barcodes (dry run / --execute) and diff old vs new.

Fixtures (tests/fixtures/reanalyse/): run_info_full.json / run_info_demix.json are trimmed copies of the real
run_info.json shapes (see tests/fixtures/replay); old/ and new/ hold tiny `summary/` tables in the real formats
(freyja_aggregated.tsv, demix_status.tsv, lineage_summary.csv). `make_run` builds results trees under tmp_path.
"""

import json
import re
import shlex
import shutil
from pathlib import Path

import pytest
import reanalyse_runs as ra

FIX = Path(__file__).resolve().parents[1] / "fixtures" / "reanalyse"
REPO = Path(__file__).resolve().parents[2]
OLD_ISO = "2026-10-06T05-25"


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
