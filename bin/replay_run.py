#!/usr/bin/env python3
"""Reconstruct the Nextflow command that reproduces a recorded wasteline run.

    replay_run.py <results dir | run_info.json> [--execute] [--pipeline-dir DIR] [--outdir NEW]
                  [--db_dir DIR] [--no-locked] [--force] [--params-out FILE]

Reads `pipeline_info/run_info.json` and `pipeline_info/params_*.json` of the results and prints

  1. the pipeline version to use (`git checkout <commit>` or `nextflow run charlesfoster/wasteline -r <commit>`)
  2. a params file for `-params-file` (the recorded values that matter; outdir and the barcodes are set on the command line)
  3. the Freyja barcodes pin (`--barcodes <iso stamp>` while the dataset is still stored, else the recorded
     dataset path, `container` when the container bundles the very same data, else a warning to restore it)
  4. verification of the recorded inputs (samplesheet, primer BED, reference, earlier results): sha256 against the files now
  5. the `-profile` (`locked` is added when the containers were recorded with digests)
  6. a NEW `--outdir` (default `<old outdir>_replay`); the original results are never overwritten

and then the command. `--execute` runs it, but only after printing it and only when no hard problem was found.

Runs on the user's host Python: standard library only, Python >= 3.8 syntax. Exit status 0 = the command is ready,
2 = problems (listed above the command), 1 = the record could not be read; with --execute the exit status is Nextflow's.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import shlex
import subprocess
import sys
from pathlib import Path

DEFAULT_PIPELINE = "charlesfoster/wasteline"
DEFAULT_DB_DIR = "~/.wasteline/db"

# Never copied into the replay params file: re-pointed on the command line, run-specific, or harmful to repeat
# (email would send the notification again).
REPOINTED_PARAMS = {
    "outdir",
    "barcodes",
    "compare_barcodes",
    "db_dir",
    "previous_results",
    "results_dir",
}
VOLATILE_PARAMS = {
    "trace_report_suffix",
    "pipelines_testdata_base_path",
    "help",
    "help_full",
    "show_hidden",
    "version",
}
NOTIFICATION_PARAMS = {"email", "email_on_fail"}
# bundled reference files: the pathogen config points these params at the checkout's own assets
BUNDLED_MARKERS = ("assets/references/", "assets/primer_schemes/")
MODES_WITHOUT_BARCODES = {"report", "update_db"}


class Issue:
    """A problem found while reconstructing the run. `level` is `error` or `warning`; `forceable` errors
    (checksum mismatches, wrong checkout) are downgraded by --force, the others never."""

    def __init__(self, level, text, forceable=False):
        self.level = level
        self.text = text
        self.forceable = forceable


class Plan:
    def __init__(self):
        self.lines = []  # human-readable report, in order
        self.issues = []
        self.argv = []
        self.env_extra = {}
        self.params = {}
        self.params_path = None
        self.outdir = None

    def say(self, text=""):
        self.lines.append(text)

    def warn(self, text):
        self.issues.append(Issue("warning", text))
        self.say("    WARNING: " + text)

    def error(self, text, forceable=False):
        self.issues.append(Issue("error", text, forceable))
        self.say("    ERROR: " + text)

    def blocking(self, force):
        return [i for i in self.issues if i.level == "error" and not (force and i.forceable)]


def sha256_file(path):
    digest = hashlib.sha256()
    with open(str(path), "rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def load_json(path):
    with open(str(path)) as handle:
        return json.load(handle)


def find_record(target):
    """(run_info dict or None, params dict or None, pipeline_info dir or None, results dir) for a results
    directory, a pipeline_info directory, run_info.json or params_*.json. `report` runs write no run_info.json,
    only the Nextflow params file."""
    target = Path(target)
    if target.is_file():
        pipeline_info = target.parent
    elif (target / "pipeline_info").is_dir():
        pipeline_info = target / "pipeline_info"
    elif (target / "run_info.json").is_file() or list(target.glob("params_*.json")):
        pipeline_info = target
    else:
        return None, None, None, target
    run_info_path = pipeline_info / "run_info.json"
    run_info = load_json(run_info_path) if run_info_path.is_file() else None
    params_files = sorted(pipeline_info.glob("params_*.json"))
    params = load_json(params_files[-1]) if params_files else None
    results_dir = pipeline_info.parent if pipeline_info.name == "pipeline_info" else pipeline_info
    return run_info, params, pipeline_info, results_dir


def schema_defaults(pipeline_dir):
    """{param: default} from the pipeline's nextflow_schema.json, or None."""
    path = Path(pipeline_dir) / "nextflow_schema.json"
    if not path.is_file():
        return None
    schema = load_json(path)
    defaults = {}
    for group in list(schema.get("$defs", {}).values()) + list(
        schema.get("definitions", {}).values()
    ):
        for name, spec in group.get("properties", {}).items():
            defaults[name] = spec.get("default")
    return defaults


def same_value(a, b):
    """Equality that does not conflate True with 1 or 0 with False."""
    if isinstance(a, bool) or isinstance(b, bool):
        return isinstance(a, bool) and isinstance(b, bool) and a == b
    return a == b


def git_head(directory):
    try:
        out = subprocess.run(
            ["git", "-C", str(directory), "rev-parse", "HEAD"],
            stdout=subprocess.PIPE,
            stderr=subprocess.DEVNULL,
            text=True,
            check=True,
        )
        return out.stdout.strip()
    except (OSError, subprocess.CalledProcessError):
        return None


def is_bundled(path):
    return any(marker in str(path).replace("\\", "/") for marker in BUNDLED_MARKERS)


def split_paths(value):
    return [part.strip() for part in str(value).split(",") if part.strip()]


# --------------------------------------------------------------------------- 1. pipeline version


def plan_version(plan, run_info, args):
    commit = (run_info or {}).get("git_commit")
    version = (run_info or {}).get("pipeline_version")
    name = (run_info or {}).get("pipeline_name") or DEFAULT_PIPELINE
    plan.say("[1] Pipeline version")
    if run_info is None:
        plan.warn(
            "no run_info.json (a `report` run writes none): the pipeline commit and version are not recorded, "
            "use the same release that made the results"
        )
    if commit and args.pipeline_dir:
        plan.say(f"    git -C {shlex.quote(args.pipeline_dir)} checkout {commit}")
        head = git_head(args.pipeline_dir)
        if head is None:
            plan.warn(
                f"{args.pipeline_dir} is not a git checkout: cannot confirm it is at {commit[:12]}"
            )
        elif head != commit:
            plan.error(
                f"{args.pipeline_dir} is at {head[:12]} but the run used {commit[:12]}: run the checkout line above first",
                forceable=True,
            )
        else:
            plan.say(f"    {args.pipeline_dir} is at the recorded commit")
        revision = None
        target = args.pipeline_dir
    elif commit:
        plan.say(f"    nextflow run {name} -r {commit}")
        revision, target = commit, name
    else:
        revision, target = None, args.pipeline_dir or name
        if run_info is not None:
            plan.warn(
                "git_commit is null (the run was not launched from a git checkout): the exact code cannot be pinned"
                + (
                    f" -- using the released version {version}"
                    if version and not args.pipeline_dir
                    else ""
                )
            )
            if version and not args.pipeline_dir:
                revision = version
                if "dev" in version:
                    plan.warn(
                        f"{version} is a development version; the tag may not exist, give --pipeline-dir"
                    )
    plan.say("    (uncommitted changes in the original checkout are not recorded)")
    nextflow_version = (run_info or {}).get("nextflow_version")
    if nextflow_version and not args.no_nxf_ver:
        plan.env_extra["NXF_VER"] = nextflow_version
        plan.say(f"    Nextflow {nextflow_version} (NXF_VER)")
    return target, revision


# --------------------------------------------------------------------------- 3. Freyja barcodes


def dataset_file_sha(directory, filename):
    path = Path(directory) / filename
    return sha256_file(path) if filename and path.is_file() else None


def resolve_barcodes(plan, block, label, db_dir_opt, recorded_db):
    """(value for --barcodes/--compare_barcodes, db_dir to pass or None) for one recorded freyja_data block."""
    source = block.get("source")
    iso = block.get("barcode_version_iso")
    version = block.get("barcode_version")
    fname = block.get("barcodes_file")
    want_sha = block.get("barcodes_sha256")
    plan.say(f"    {label}: source {source}, version {version}")

    def verify_dir(directory):
        got = dataset_file_sha(directory, fname)
        if got is not None and want_sha and got != want_sha:
            plan.error(
                f"{label}: {directory}/{fname} has sha256 {got[:12]} but the run used {want_sha[:12]}",
                forceable=True,
            )
            return False
        return True

    if source in ("latest", "pinned"):
        candidates = []
        for db in (db_dir_opt, recorded_db, DEFAULT_DB_DIR):
            if db and db not in candidates:
                candidates.append(db)
        for db in candidates:
            directory = Path(db).expanduser() / "freyja" / str(iso)
            if iso and directory.is_dir():
                if verify_dir(directory):
                    plan.say(f"      stored dataset {directory} (sha256 checked)")
                return str(iso), (None if db == DEFAULT_DB_DIR else db)
        dataset_path = block.get("dataset_path")
        for candidate in [
            dataset_path,
            str(Path(dataset_path) / str(iso)) if dataset_path and iso else None,
        ]:
            if (
                candidate
                and Path(candidate).is_dir()
                and dataset_file_sha(candidate, fname) is not None
            ):
                verify_dir(candidate)
                plan.say(f"      recorded dataset_path {candidate}")
                return candidate, None
        if version and version == block.get("freyja_bundled_barcode_version"):
            plan.say(
                "      the dataset is no longer stored, but Freyja {} bundles the same barcodes".format(
                    block.get("freyja_version")
                )
            )
            return "container", None
        plan.error(
            "{}: the barcode dataset {} MUST BE RESTORED to <db_dir>/freyja/{}/ (searched {}); old versions cannot "
            "be downloaded again, restore it from an archive. The pipeline never substitutes other data".format(
                label, version, iso, ", ".join(candidates)
            )
        )
        return str(iso), db_dir_opt or recorded_db
    if source == "container":
        plan.say(
            "      the data bundled in the Freyja {} container (same pinned container = same data)".format(
                block.get("freyja_version")
            )
        )
        return "container", None
    if source == "custom":
        dataset_path = block.get("dataset_path")
        if not dataset_path or not Path(dataset_path).exists():
            plan.error(f"{label}: custom dataset {dataset_path} not found: restore it")
            return dataset_path, None
        path = Path(dataset_path)
        got = sha256_file(path) if path.is_file() else dataset_file_sha(path, fname)
        if got is None:
            plan.error(f"{label}: {dataset_path} does not contain {fname}")
        elif want_sha and got != want_sha:
            plan.error(
                f"{label}: {dataset_path} has sha256 {got[:12]} but the run used {want_sha[:12]}",
                forceable=True,
            )
        else:
            plan.say(f"      {dataset_path} (sha256 checked)")
        return dataset_path, None
    plan.error(f"{label}: unknown barcodes source {source!r}")
    return None, None


def plan_barcodes(plan, run_info, params, mode, args):
    """Returns the extra command-line pairs [(flag, value)]."""
    plan.say("[3] Freyja barcodes")
    extras = []
    if mode in MODES_WITHOUT_BARCODES:
        plan.say(f"    not used in mode {mode}")
        return extras
    block = (run_info or {}).get("freyja_data")
    if not block:
        plan.warn(
            "no freyja_data record: --barcodes {} (as recorded in the params) is used as is".format(
                params.get("barcodes")
            )
        )
        if params.get("barcodes"):
            extras.append(("--barcodes", str(params["barcodes"])))
        return extras
    if mode == "aggregate":
        plan.say(
            "    (aggregate mode: this dataset only groups the plot; the .demix files keep their own barcodes)"
        )
    recorded_db = params.get("db_dir")
    dbs_used = []
    value, db = resolve_barcodes(plan, block, "barcodes", args.db_dir, recorded_db)
    if value:
        extras.append(("--barcodes", value))
    if db:
        dbs_used.append(db)
    if params.get("compare_barcodes"):
        compare = (run_info or {}).get("freyja_data_compare")
        if compare:
            value, db = resolve_barcodes(
                plan, compare, "compare_barcodes", args.db_dir, recorded_db
            )
            if value:
                extras.append(("--compare_barcodes", value))
            if db and db not in dbs_used:
                dbs_used.append(db)
        else:
            plan.warn(
                "compare_barcodes was {!r} but no freyja_data_compare is recorded: "
                "'latest' would fetch different data now, pin it".format(params["compare_barcodes"])
            )
            extras.append(("--compare_barcodes", str(params["compare_barcodes"])))
    if dbs_used:
        if len(dbs_used) > 1:
            plan.warn(
                "the datasets are stored in different db_dir locations ({}): only the first is passed".format(
                    ", ".join(dbs_used)
                )
            )
        extras.append(("--db_dir", dbs_used[0]))
    return extras


# --------------------------------------------------------------------------- 4. inputs


def check_recorded_file(plan, label, path, entry, required=True):
    """Compare the file at `path` with the recorded {file, sha256}. Missing files are errors only if `required`."""
    want = entry.get("sha256") if isinstance(entry, dict) else None
    if path is None:
        plan.say(f"    {label:<16} not checked (no file to compare)")
        return
    if not Path(path).is_file():
        if required:
            plan.error(f"{label}: {path} not found")
        else:
            plan.warn(f"{label}: {path} not found")
        return
    if not want:
        plan.say(f"    {label:<16} {path} (no checksum recorded)")
        return
    got = sha256_file(path)
    if got == want:
        plan.say(f"    {label:<16} OK       {path}")
    else:
        plan.error(
            f"{label}: {path} DIFFERS from the run (sha256 {got[:12]} vs recorded {want[:12]})",
            forceable=True,
        )


def check_reads(plan, samplesheet):
    """Every FASTQ named in a full-mode samplesheet should exist (absolute, or relative to cwd or the sheet)."""
    missing = []
    try:
        with open(str(samplesheet), newline="") as handle:
            rows = list(csv.DictReader(handle))
    except (OSError, csv.Error, UnicodeDecodeError):
        return
    for row in rows:
        for column in ("fastq_1", "fastq_2"):
            value = (row.get(column) or "").strip()
            if not value:
                continue
            found = any(
                candidate.exists() for candidate in (Path(value), Path(samplesheet).parent / value)
            )
            if not found and "://" not in value:
                missing.append(value)
    if missing:
        plan.error(
            "{} FASTQ file(s) of the samplesheet not found, e.g. {}".format(
                len(missing), ", ".join(missing[:3])
            )
        )
    elif rows:
        plan.say(
            "    {:<16} OK       {} sample(s), all FASTQ files present".format("reads", len(rows))
        )


def plan_inputs(plan, run_info, params, mode, args):
    plan.say("[4] Inputs")
    sums = (run_info or {}).get("input_checksums") or {}
    if run_info is None:
        plan.say("    no input checksums recorded")
    sheet = params.get("input")
    if sheet:
        check_recorded_file(plan, "samplesheet", sheet, sums.get("samplesheet", {}))
        if mode == "full" and Path(str(sheet)).is_file():
            check_reads(plan, sheet)
    elif mode in ("demix", "aggregate"):
        plan.say("    samplesheet      built from the earlier results (see below)")
    if "primer_bed" in sums:
        bed = params.get("primer_bed")
        if not bed and args.pipeline_dir:
            bed = str(
                Path(args.pipeline_dir)
                / "assets"
                / "primer_schemes"
                / sums["primer_bed"].get("file", "")
            )
        if bed and not Path(bed).is_file() and is_bundled(bed) and args.pipeline_dir:
            bed = str(Path(args.pipeline_dir) / "assets" / "primer_schemes" / Path(bed).name)
        if bed:
            check_recorded_file(
                plan, "primer BED", bed, sums["primer_bed"], required=not is_bundled(bed)
            )
        else:
            plan.say(
                "    {:<16} not checked: bundled scheme {!r}, fixed by the commit (give --pipeline-dir to check)".format(
                    "primer BED", params.get("primer_scheme")
                )
            )
    if "reference_fasta" in sums:
        fasta = params.get("fasta")
        if fasta and not Path(fasta).is_file() and is_bundled(fasta) and args.pipeline_dir:
            fasta = str(Path(args.pipeline_dir) / "assets" / "references" / Path(fasta).name)
        if fasta and (Path(fasta).is_file() or not is_bundled(fasta)):
            check_recorded_file(plan, "reference", fasta, sums["reference_fasta"])
        else:
            plan.say("    {:<16} not checked: bundled with the commit".format("reference"))
    previous = (run_info or {}).get("previous_runs") or []
    for record in previous:
        path = record.get("path")
        if not path or not Path(path).is_dir():
            plan.error(f"earlier results {path} not found")
            continue
        info = Path(path) / "pipeline_info" / "run_info.json"
        want = record.get("run_info_sha256")
        if want and info.is_file():
            got = sha256_file(info)
            if got == want:
                plan.say("    {:<16} OK       {}".format("earlier results", path))
            else:
                plan.error(
                    f"earlier results {path} changed: run_info.json sha256 {got[:12]} vs recorded {want[:12]}",
                    forceable=True,
                )
        elif want:
            plan.error(
                f"earlier results {path} have no pipeline_info/run_info.json (recorded sha256 {want[:12]})",
                forceable=True,
            )
        else:
            plan.say(
                "    {:<16} present  {} (legacy results: no checksum recorded)".format(
                    "earlier results", path
                )
            )
    results_dir = params.get("results_dir")
    if mode == "report" and results_dir:
        if not Path(results_dir).is_dir():
            plan.error(f"results_dir {results_dir} not found")
        elif not (Path(results_dir) / "pipeline_info" / "run_info.json").is_file():
            plan.error(f"results_dir {results_dir} has no pipeline_info/run_info.json")
        else:
            plan.say("    {:<16} present  {}".format("results_dir", results_dir))
    demix_dir = params.get("demix_dir")
    if demix_dir:
        if Path(demix_dir).is_dir():
            plan.say("    {:<16} present  {}".format("demix_dir", demix_dir))
        else:
            plan.error(f"demix_dir {demix_dir} not found")


# --------------------------------------------------------------------------- 2. params


def plan_params(plan, run_info, params, mode, args):
    plan.say("[2] Parameters")
    defaults = schema_defaults(args.pipeline_dir) if args.pipeline_dir else None
    if args.pipeline_dir and defaults is None:
        plan.warn(
            f"{args.pipeline_dir}/nextflow_schema.json not found: writing all recorded parameters"
        )
    kept = {}
    dropped_default = 0
    for name in sorted(params):
        value = params[name]
        if name in REPOINTED_PARAMS or name in VOLATILE_PARAMS:
            continue
        if name in NOTIFICATION_PARAMS:
            if value:
                plan.say(
                    f"    {name} ({value!r}) left out: the replay should not send the notification again"
                )
            continue
        if name in ("fasta", "gff") and value and is_bundled(value):
            continue  # the checkout's own bundled file (set by the pathogen config), verified by checksum above
        if defaults is not None:
            if name not in defaults:
                plan.warn(
                    f"{name} is not a parameter of {args.pipeline_dir}: nf-schema may reject it"
                )
            elif same_value(value, defaults[name]):
                dropped_default += 1
                continue
        elif value is None:
            continue
        kept[name] = value
    plan.params = kept
    where = (
        f"schema defaults of {args.pipeline_dir} omitted ({dropped_default})"
        if defaults is not None
        else ("all recorded values except nulls (give --pipeline-dir to also omit the defaults)")
    )
    plan.say(f"    {len(kept)} value(s) -> {args.params_out}; {where}")
    plan.say(
        "    set on the command line instead: " + ", ".join(sorted(REPOINTED_PARAMS & set(params)))
    )
    return kept


def previous_results_flag(plan, run_info, params):
    """--previous_results from the recorded previous_runs paths (the resolved absolute directories)."""
    if not params.get("previous_results"):
        return []
    recorded = [r.get("path") for r in (run_info or {}).get("previous_runs") or [] if r.get("path")]
    paths = recorded or split_paths(params["previous_results"])
    return [("--previous_results", ",".join(str(p) for p in paths))]


# --------------------------------------------------------------------------- 5. profile, 6. outdir


def plan_profile(plan, run_info, args):
    plan.say("[5] Profile")
    recorded = (run_info or {}).get("profile")
    tokens = [t for t in str(recorded or "").split(",") if t]
    if not tokens:
        tokens = ["docker"]
        plan.warn(
            "the -profile is not recorded: assuming docker, change it to your container engine"
        )
    containers = (run_info or {}).get("containers") or {}
    digests = bool(containers) and all("@sha256:" in str(v) for v in containers.values())
    if args.no_locked:
        tokens = [t for t in tokens if t != "locked"]
        plan.say("    --no-locked: containers are pulled by tag")
    elif digests and "locked" not in tokens:
        tokens.append("locked")
        plan.say("    the containers were recorded with digests: adding the locked profile")
    elif "locked" in tokens:
        plan.say("    locked (as recorded)")
    plan.say("    -profile {}".format(",".join(tokens)))
    return ",".join(tokens)


def plan_outdir(plan, params, results_dir, args):
    plan.say("[6] Output directory")
    original = params.get("outdir") or str(results_dir)
    outdir = args.outdir or str(original).rstrip("/") + "_replay"
    candidates = {os.path.realpath(str(original)), os.path.realpath(str(results_dir))}
    if os.path.realpath(outdir) in candidates:
        plan.error(
            f"--outdir {outdir} is the original results directory: the replay must write somewhere new"
        )
    elif Path(outdir).exists() and any(Path(outdir).iterdir()):
        plan.warn(f"{outdir} already exists and is not empty")
    plan.say(f"    --outdir {outdir}   (original: {original})")
    return outdir


# --------------------------------------------------------------------------- command


def format_command(env_extra, argv):
    """One line per option with its value, shell-quoted."""
    groups = []
    for word in argv:
        if word.startswith("-") or not groups:
            groups.append([word])
        else:
            groups[-1].append(word)
    lines = [" ".join(shlex.quote(w) for w in group) for group in groups]
    prefix = " ".join(f"{k}={shlex.quote(v)}" for k, v in sorted(env_extra.items()))
    if prefix:
        lines[0] = prefix + " " + lines[0]
    return " \\\n    ".join(lines)


def build_plan(args):
    plan = Plan()
    run_info, params, _pipeline_info, results_dir = find_record(args.results)
    if run_info is None and params is None:
        raise SystemExit(
            f"error: no pipeline_info/run_info.json or params_*.json found in {args.results}"
        )
    recorded = dict((run_info or {}).get("params") or params or {})
    mode = (run_info or {}).get("mode") or recorded.get("mode") or "full"
    plan.say(f"wasteline replay of {results_dir}")
    if run_info:
        plan.say(
            "  recorded: {} {}, mode {}, started {}, success {}".format(
                run_info.get("pipeline_name"),
                run_info.get("pipeline_version"),
                mode,
                run_info.get("start"),
                run_info.get("success", "not recorded"),
            )
        )
        if run_info.get("success") is False:
            plan.warn("the recorded run did not succeed")
        elif "success" not in run_info:
            plan.warn("the recorded run has no success flag: it may have crashed before completion")
    plan.say()

    target, revision = plan_version(plan, run_info, args)
    plan.say()
    args.params_out = args.params_out or "replay_params.json"
    plan_params(plan, run_info, recorded, mode, args)
    plan.say()
    extras = plan_barcodes(plan, run_info, recorded, mode, args)
    plan.say()
    plan_inputs(plan, run_info, recorded, mode, args)
    plan.say()
    profile = plan_profile(plan, run_info, args)
    plan.say()
    outdir = plan_outdir(plan, recorded, results_dir, args)

    argv = [args.nextflow, "run", target]
    if revision:
        argv += ["-r", revision]
    argv += ["-profile", profile, "-params-file", args.params_out, "--outdir", outdir]
    for flag, value in extras:
        argv += [flag, value]
    for flag, value in previous_results_flag(plan, run_info, recorded):
        argv += [flag, value]
    if mode == "report" and recorded.get("results_dir"):
        argv += ["--results_dir", str(recorded["results_dir"])]
    plan.argv = argv
    plan.outdir = outdir
    return plan


def parse_args(argv):
    parser = argparse.ArgumentParser(
        description="Print (and optionally run) the Nextflow command that reproduces a recorded wasteline run.",
        epilog="Exit status: 0 command ready, 2 problems found, 1 record unreadable.",
    )
    parser.add_argument(
        "results", help="results directory (or its pipeline_info/, or run_info.json)"
    )
    parser.add_argument(
        "--execute",
        action="store_true",
        help="run the command after printing it (refused on hard errors)",
    )
    parser.add_argument(
        "--pipeline-dir",
        help="local checkout of the pipeline: used instead of GitHub, and its nextflow_schema.json "
        "defaults are left out of the params file",
    )
    parser.add_argument("--outdir", help="NEW output directory (default: <original outdir>_replay)")
    parser.add_argument(
        "--db_dir",
        "--db-dir",
        dest="db_dir",
        help="where the stored Freyja datasets are now (searched before the recorded db_dir)",
    )
    parser.add_argument(
        "--no-locked",
        action="store_true",
        help="do not use -profile locked (pull containers by tag)",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="--execute despite checksum mismatches or a different checkout (never for missing files)",
    )
    parser.add_argument(
        "--params-out", help="where to write the params file (default: ./replay_params.json)"
    )
    parser.add_argument(
        "--no-nxf-ver",
        action="store_true",
        help="do not set NXF_VER to the recorded Nextflow version",
    )
    parser.add_argument(
        "--nextflow", default="nextflow", help="the nextflow executable (default: nextflow)"
    )
    return parser.parse_args(argv)


def main(argv=None):
    args = parse_args(sys.argv[1:] if argv is None else argv)
    try:
        plan = build_plan(args)
    except (OSError, ValueError) as err:
        print(f"error: cannot read the record: {err}", file=sys.stderr)
        return 1
    with open(args.params_out, "w") as handle:
        json.dump(plan.params, handle, indent=2, sort_keys=True)
        handle.write("\n")
    print("\n".join(plan.lines))
    blocking = plan.blocking(args.force)
    warnings = [i for i in plan.issues if i.level == "warning"]
    forced = [i for i in plan.issues if i.level == "error" and args.force and i.forceable]
    print()
    print(
        "Problems: {} error(s){}, {} warning(s)".format(
            len(blocking),
            f" (+{len(forced)} overridden by --force)" if forced else "",
            len(warnings),
        )
    )
    print()
    print("Command" + (" (NOT READY: fix the errors above)" if blocking else "") + ":")
    print()
    print(format_command(plan.env_extra, plan.argv))
    print()
    if blocking:
        if args.execute:
            print(f"Not executing: {len(blocking)} error(s) above.", file=sys.stderr)
        return 2
    if args.execute:
        env = dict(os.environ)
        env.update(plan.env_extra)
        print("Executing...", file=sys.stderr)
        return subprocess.call(plan.argv, env=env)
    return 0


if __name__ == "__main__":
    sys.exit(main())
