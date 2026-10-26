#!/usr/bin/env python3
"""Re-run past results with newer Freyja barcodes WITHOUT re-mapping, and see what changed.

    reanalyse_runs.py <root> [<root2> ...] [--since YYYY-MM-DD | --last N | --runs a,b]
                      [--barcodes latest|YYYY-MM-DD|stamp|container|<path>] [--pipeline-dir DIR]
                      [--outdir-suffix _reanalysis_{barcode}] [--nextflow BIN] [--extra-args '...']
                      [--execute] [--keep-going] [--db_dir DIR] [--threshold 0.05]
                      [--report diff.csv] [--summary-md diff.md]
    reanalyse_runs.py --diff <old results dir> <new results dir> [--threshold 0.05]
                      [--level lineages|groups|summarized|all] [--report diff.csv] [--summary-md diff.md]

Each <root> is a results directory or a folder that contains results directories (found by
`pipeline_info/run_info.json`, at most two levels down; legacy Snakemake layouts `<sample>/freyja/` are
found too and count as barcode version `unknown`).

Default is a DRY RUN: a table of the runs found and, per run, the exact `nextflow run ... --mode demix
--previous_results <run>` command that `--execute` would run. Nothing is written. The original results are
never touched: every re-analysis goes to a NEW directory `<run>_reanalysis_<barcode version>` next to it
(never overwritten; `--force-new-suffix` adds a counter instead of skipping), with its params file
`<new dir>.params.json` beside it. Runs whose recorded barcodes already equal the target are `up to date`.

`--barcodes latest` needs the network: with `--execute` the script first runs `--mode update_db` ONCE, reads the
stored dataset's version from its `pipeline_info/freyja_data.json` and uses that exact stamp for every run, so all
runs share one dataset. In a dry run `latest` is not resolved (nothing is downloaded) and no run is skipped as up to
date. The script itself never uses the network.

After each executed re-analysis the old and new results are compared (`--diff` does this on its own): per sample
and item the abundance old vs new, flagged when |change| >= --threshold; items below 0.001 in both are ignored.
Written to `<new dir>/summary/reanalysis_diff.csv` and `.md`.

Runs on the user's host Python: standard library only, Python >= 3.8 syntax. Exit status: 0 ok, 1 error
(unreadable input, nothing found, a Nextflow run failed).
"""

from __future__ import annotations

import argparse
import ast
import csv
import datetime
import hashlib
import json
import os
import re
import shlex
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

DEFAULT_DB_DIR = "~/.wasteline/db"
DEFAULT_SUFFIX = "_reanalysis_{barcode}"
MIN_ABUNDANCE = 0.001  # items below this in both runs are ignored by the diff
LEVELS = ("lineages", "summarized", "groups")

# --- discovery -----------------------------------------------------------------------------------------
# Directories under a results dir that are never legacy sample directories (as in make_samplesheet.py).
RESERVED_DIRS = {"samples", "summary", "pipeline_info", "db"}
SKIP_WALK = {"work", "db", "node_modules"}
DATASET_MARKERS = (
    "last_barcode_update.txt",
    "last_barcode_update_other.txt",
    "usher_barcodes.feather",
    "usher_barcodes.csv",
    "lineages.yml",
)
# modes whose output does not keep variants/depths: they cannot be demixed again
NO_VARIANTS_MODES = {"demix", "aggregate", "report", "update_db"}

# --- parameters carried over ---------------------------------------------------------------------------
# Only what matters in `--mode demix` (the run has no reads any more). Everything else is left out.
CARRY_EXACT = {
    "pathogen",
    "primer_scheme",
    "primer_scheme_dir",
    "primer_bed",
    "protocol",
    "skip_primer_trimming",
    "min_genome_coverage_for_demix",
    "mutations_of_interest",
    "annotate_variants",
    "include_negatives_in_aggregate",
    "skip_consensus",
    "skip_report",
    "indel_freq",
    "timeseries",
    "primer_scheme_check",
    "fasta",
    "gff",
}
CARRY_PREFIXES = (
    "freyja_",
    "consensus_",
    "qc_",
    "unexplained_",
    "lineage_",
    "timeseries_",
    "ntc_",
    "dropout_",
    "primer_check",
    "amplicon_",
    "site_",
)
PATH_PARAMS = {
    "primer_bed",
    "primer_scheme_dir",
    "mutations_of_interest",
    "lineage_groups",
    "site_aggregation",
    "fasta",
    "gff",
}
# never dropped because it equals the schema default: in demix mode amplicon coverage runs only when a primer
# scheme/BED is given explicitly (PLAN D14), so the effective scheme of the run must be passed on
ALWAYS_KEEP = {"primer_scheme", "primer_bed"}
# a panel names the targets (reference and pathogen included): the pipeline refuses these together with --panel
PANEL_DROPPED = {"pathogen", "fasta", "gff"}
BUNDLED_MARKERS = ("assets/references/", "assets/primer_schemes/", "assets/lineage_groups")
DROPPED_PROFILES = {"locked", "test", "test_stub", "test_full"}


class Problem(Exception):
    """A friendly, fatal error (printed without a traceback)."""


# =========================================================================================================
# helpers
# =========================================================================================================


def load_json(path):
    with open(str(path)) as handle:
        return json.load(handle)


def sha256_file(path):
    digest = hashlib.sha256()
    with open(str(path), "rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def same_value(a, b):
    """Equality that does not conflate True with 1 or 0 with False."""
    if isinstance(a, bool) or isinstance(b, bool):
        return isinstance(a, bool) and isinstance(b, bool) and a == b
    return a == b


def schema_defaults(pipeline_dir):
    """{param: default} from the pipeline's nextflow_schema.json, or None."""
    path = Path(str(pipeline_dir)) / "nextflow_schema.json"
    if not path.is_file():
        return None
    try:
        schema = load_json(path)
    except (OSError, ValueError):
        return None
    defaults = {}
    for group in list(schema.get("$defs", {}).values()) + list(
        schema.get("definitions", {}).values()
    ):
        for name, spec in group.get("properties", {}).items():
            defaults[name] = spec.get("default")
    return defaults


def to_iso(stamp):
    """`MM_DD_YYYY-HH-MM` or `YYYY-MM-DDTHH-MM` -> `YYYY-MM-DDTHH-MM`; anything else -> None."""
    text = str(stamp or "").strip()
    match = re.fullmatch(r"(\d{2})_(\d{2})_(\d{4})-(\d{2})-(\d{2})", text)
    if match:
        month, day, year, hour, minute = match.groups()
        return f"{year}-{month}-{day}T{hour}-{minute}"
    if re.fullmatch(r"\d{4}-\d{2}-\d{2}T\d{2}-\d{2}", text):
        return text
    return None


def quote_command(argv):
    """One line per option with its values, shell-quoted (every path with spaces is quoted)."""
    groups = []
    for word in argv:
        if word.startswith("-") or not groups:
            groups.append([word])
        else:
            groups[-1].append(word)
    return " \\\n    ".join(" ".join(shlex.quote(w) for w in group) for group in groups)


def is_bundled(path):
    return any(marker in str(path).replace("\\", "/") for marker in BUNDLED_MARKERS)


# =========================================================================================================
# discovery of runs
# =========================================================================================================


def is_dataset_dir(path):
    """A `freyja/` dir that holds Freyja data is a downloaded dataset, not a sample's results."""
    freyja = path / "freyja"
    if not freyja.is_dir():
        return False
    for child in [freyja] + [c for c in freyja.iterdir() if c.is_dir()]:
        if any((child / marker).exists() for marker in DATASET_MARKERS):
            return True
    return False


def child_dirs(path):
    try:
        return sorted(p for p in path.iterdir() if p.is_dir() and not p.name.startswith("."))
    except OSError:
        return []


def looks_legacy(path):
    """Legacy Snakemake results: `<sample>/freyja/` or `<sample>/variants/` directly under the results."""
    for child in child_dirs(path):
        if child.name in RESERVED_DIRS:
            continue
        if ((child / "freyja").is_dir() or (child / "variants").is_dir()) and not is_dataset_dir(
            child
        ):
            return True
    return False


def panel_targets(directory):
    """Targets of the results of a panel run: the folders of pipeline_info/ that hold a run_info.json
    (pipeline_info/<target>/run_info.json). Each (directory, target) is one run for this script."""
    info = Path(directory) / "pipeline_info"
    return sorted(c.name for c in child_dirs(info) if (c / "run_info.json").is_file())


def panel_result_targets(directory):
    """Targets of panel results by their tables (summary/<target>/freyja_aggregated.tsv, and no
    summary/freyja_aggregated.tsv); [] for every other layout."""
    summary = Path(directory) / "summary"
    if (summary / "freyja_aggregated.tsv").is_file():
        return []
    return sorted(c.name for c in child_dirs(summary) if (c / "freyja_aggregated.tsv").is_file())


def find_runs(root, max_depth=2):
    """[(directory, 'wasteline'|'legacy'|'panel')] at depth <= max_depth below `root` (root itself is depth 0).
    A run is never searched further (its own directories are samples, not runs)."""
    found = []

    def walk(directory, depth):
        if (directory / "pipeline_info" / "run_info.json").is_file():
            found.append((directory, "wasteline"))
            return
        if panel_targets(directory):
            found.append((directory, "panel"))
            return
        if looks_legacy(directory):
            found.append((directory, "legacy"))
            return
        if depth >= max_depth:
            return
        for child in child_dirs(directory):
            if child.name in SKIP_WALK:
                continue
            walk(child, depth + 1)

    walk(Path(root), 0)
    return found


def has_files(directory, suffixes):
    try:
        return any(name.endswith(suffixes) for name in os.listdir(str(directory)))
    except OSError:
        return False


def count_samples(path, kind, target=""):
    """Samples that still have both variants and depths (what a demix needs). Cheap: directory listings only.
    `target`: panel results, samples/<sample>/<target>/freyja/."""
    count = 0
    if kind == "wasteline":
        base = path / "samples"
        for sample in child_dirs(base) if base.is_dir() else []:
            freyja = (sample / target if target else sample) / "freyja"
            if has_files(freyja, (".variants.tsv", ".variants.tsv.gz")) and has_files(
                freyja, (".depths", ".depths.gz")
            ):
                count += 1
        return count
    for sample in child_dirs(path):
        if sample.name in RESERVED_DIRS or is_dataset_dir(sample):
            continue
        freyja = sample / "freyja"
        variants = has_files(
            freyja, ("freyja_variants.tsv", "freyja_variants.tsv.gz")
        ) or has_files(sample / "variants", ("freyja_variants.tsv", "freyja_variants.tsv.gz"))
        if variants and has_files(freyja, (".depths", ".depths.gz")):
            count += 1
    return count


def recorded_barcode(run_info):
    """(display version, barcode_version_iso or None, barcodes_sha256 or None, source or None)."""
    block = (run_info or {}).get("freyja_data")
    if not block:
        return "unknown", None, None, None
    iso = block.get("barcode_version_iso") or to_iso(block.get("barcode_version"))
    shown = iso or block.get("barcode_version") or "unknown"
    return shown, iso, block.get("barcodes_sha256"), block.get("source")


class Run:
    """One discovered results directory."""

    def __init__(self, path, kind, target=""):
        self.path = Path(os.path.abspath(str(path)))
        self.kind = kind
        self.target = target  # panel results: one Run per (directory, target)
        self.name = self.path.name
        self.run_info = None
        self.params = {}
        self.mode = "legacy" if kind == "legacy" else "unknown"
        self.start = ""
        self.profile = ""
        self.freyja_version = "unknown"
        self.barcode = "unknown"  # what the table shows
        self.barcode_iso = None
        self.barcode_sha = None
        self.barcode_source = None
        self.n_samples = 0
        self.success = None
        self.skip_reason = None  # set when the run cannot be re-analysed at all
        self.notes = []
        self.load()

    @property
    def label(self):
        """`run1` or, for panel results, `run1 [RSVa]`."""
        return f"{self.name} [{self.target}]" if self.target else self.name

    @property
    def summary_dir(self):
        return self.path / "summary" / self.target if self.target else self.path / "summary"

    def load(self):
        if self.kind == "wasteline":
            info_path = (
                self.path / "pipeline_info" / self.target / "run_info.json"
                if self.target
                else self.path / "pipeline_info" / "run_info.json"
            )
            try:
                self.run_info = load_json(info_path)
            except (OSError, ValueError) as err:
                raise Problem(f"cannot read {info_path}: {err}")
            info = self.run_info
            self.params = dict(info.get("params") or {})
            if not self.params:
                files = sorted((self.path / "pipeline_info").glob("params_*.json"))
                if files:
                    try:
                        self.params = dict(load_json(files[-1]))
                    except (OSError, ValueError):
                        pass
            self.mode = info.get("mode") or self.params.get("mode") or "unknown"
            self.start = str(info.get("start") or "")
            self.profile = str(info.get("profile") or "")
            self.freyja_version = str(info.get("freyja_version") or "unknown")
            self.barcode, self.barcode_iso, self.barcode_sha, self.barcode_source = (
                recorded_barcode(info)
            )
            self.success = info.get("success") if "success" in info else None
        else:
            try:
                stamp = datetime.datetime.fromtimestamp(self.path.stat().st_mtime)
                self.start = stamp.strftime("%Y-%m-%dT%H:%M:%S")
            except OSError:
                self.start = ""
        self.n_samples = count_samples(self.path, self.kind, self.target)
        self.skip_reason = self.ineligible_reason()

    @property
    def date(self):
        return self.start[:10]

    def ineligible_reason(self):
        if self.mode in NO_VARIANTS_MODES:
            previous = [
                r.get("path")
                for r in (self.run_info or {}).get("previous_runs") or []
                if r.get("path")
            ]
            hint = f"; re-analyse its source instead: {', '.join(previous)}" if previous else ""
            return f"mode {self.mode} keeps no variants/depths{hint}"
        if self.success is False:
            return "the recorded run did not succeed"
        if self.n_samples == 0:
            return "no sample has variants and depths files"
        return None


def discover(roots):
    """(runs, messages). Roots may overlap; each directory appears once. Deterministic order."""
    runs = {}
    messages = []
    for root in roots:
        root_path = Path(root)
        if not root_path.is_dir():
            raise Problem(f"not a directory: {root}")
        found = find_runs(root_path)
        if not found:
            hint = ""
            if (root_path / "samples").is_dir() and not (
                root_path / "pipeline_info" / "run_info.json"
            ).is_file():
                hint = " (it has samples/ but no pipeline_info/run_info.json: is it the output of a report run?)"
            messages.append(
                f"no results found in {root}{hint}: looked for pipeline_info/run_info.json "
                "(or a legacy <sample>/freyja/ layout) up to 2 levels down"
            )
            continue
        for directory, kind in found:
            key = os.path.realpath(str(directory))
            if kind == "panel":
                for target in panel_targets(directory):
                    if (key, target) not in runs:
                        runs[(key, target)] = Run(directory, "wasteline", target)
            elif key not in runs:
                runs[key] = Run(directory, kind)
    ordered = sorted(runs.values(), key=lambda r: (r.date, r.start, str(r.path), r.target))
    return ordered, messages


def select_runs(runs, args):
    """Apply --runs, --since, --last (the last N re-analysable runs by date). Raises Problem for unknown names."""
    chosen = list(runs)
    if args.runs:
        wanted = [n.strip() for n in args.runs.split(",") if n.strip()]
        picked = []
        for name in wanted:
            real = os.path.realpath(name)
            hits = [
                r
                for r in chosen
                if r.name == name
                or r.label == name
                or os.path.realpath(str(r.path)) == real
                or str(r.path).endswith("/" + name.strip("/"))
            ]
            if not hits:
                raise Problem(
                    f"--runs: no run named {name!r}; found: " + ", ".join(r.label for r in runs)
                )
            picked.extend(h for h in hits if h not in picked)
        chosen = picked
    if args.since:
        chosen = [r for r in chosen if r.date and r.date >= args.since]
    if args.last is not None:
        # the N most recent directories (a panel directory counts once, whatever its number of targets)
        eligible = []
        for r in chosen:
            if r.skip_reason is None and r.path not in eligible:
                eligible.append(r.path)
        keep = set(eligible[-args.last :]) if args.last > 0 else set()
        chosen = [r for r in chosen if r.skip_reason is not None or r.path in keep]
    return sorted(chosen, key=lambda r: (r.date, r.start, str(r.path), r.target))


# =========================================================================================================
# the target barcodes
# =========================================================================================================


class Target:
    """What `--barcodes` means for the comparison with the recorded barcodes."""

    def __init__(self, spec, kind, iso=None, sha=None, label=None, date=None):
        self.spec = spec  # passed on the nextflow command line
        self.kind = kind  # latest | container | date | stamp | path
        self.iso = iso
        self.sha = sha
        self.date = date
        self.label = label or spec
        self.pending = kind == "latest"

    def describe(self):
        if self.kind == "latest":
            return "latest (resolved when run)"
        if self.kind == "container":
            return "container (the data bundled in the Freyja container of the pipeline)"
        return self.iso or self.date or self.label

    def is_current(self, run):
        """True when the run already used this dataset (so a re-analysis would change nothing)."""
        if self.kind in ("latest", "container"):
            return False
        if self.sha and run.barcode_sha:
            return self.sha == run.barcode_sha
        if self.iso and run.barcode_iso:
            return self.iso == run.barcode_iso
        if self.kind == "date" and self.date and run.barcode_iso:
            return run.barcode_iso[:10] == self.date
        return False


def barcode_file_of(directory):
    for name in ("usher_barcodes.feather", "usher_barcodes.csv"):
        if (directory / name).is_file():
            return directory / name
    return None


def resolve_target(spec, db_dir):
    """Classify `--barcodes` without any network access. Raises Problem for a spec that cannot be understood."""
    text = str(spec).strip()
    if text == "latest":
        return Target(text, "latest", label="latest")
    if text == "container":
        return Target(text, "container", label="container")
    path = Path(text).expanduser()
    if path.is_dir():
        stamp_file = path / "last_barcode_update.txt"
        iso = None
        if stamp_file.is_file():
            iso = to_iso(stamp_file.read_text().strip())
        barcodes = barcode_file_of(path)
        sha = sha256_file(barcodes) if barcodes else None
        return Target(os.path.abspath(str(path)), "path", iso=iso, sha=sha, label=iso or path.name)
    if path.is_file():
        return Target(
            os.path.abspath(str(path)),
            "path",
            sha=sha256_file(path),
            label=re.sub(r"\.(csv|feather)$", "", path.name),
        )
    iso = to_iso(text)
    if iso:
        return Target(text, "stamp", iso=iso, label=iso)
    if re.fullmatch(r"\d{4}-\d{2}-\d{2}", text):
        stored = Path(db_dir or DEFAULT_DB_DIR).expanduser() / "freyja"
        names = (
            sorted(p.name for p in stored.iterdir() if p.is_dir() and p.name.startswith(text))
            if stored.is_dir()
            else []
        )
        iso = names[-1] if names else None
        return Target(text, "date", iso=iso, date=text, label=iso or text)
    raise Problem(
        f"--barcodes {text!r} is not 'latest', 'container', a barcode version (YYYY-MM-DD, "
        "MM_DD_YYYY-HH-MM or YYYY-MM-DDTHH-MM) or an existing dataset directory / .csv / .feather file"
    )


# =========================================================================================================
# the command for one run
# =========================================================================================================


def profile_for(run, args):
    if args.profile:
        return args.profile, None
    tokens = [t for t in run.profile.split(",") if t and t not in DROPPED_PROFILES]
    note = None
    if not tokens:
        tokens = ["docker"]
        if run.profile:
            note = f"recorded profile {run.profile!r} is not reusable: using docker"
    return ",".join(tokens), note


def resolve_path_param(name, value, run, notes):
    """A path from the record, as an absolute path when it still exists; None (with a note) when it does not."""
    candidates = [Path(str(value)).expanduser()]
    if not Path(str(value)).is_absolute():
        candidates.append(run.path.parent / str(value))
        candidates.append(run.path / str(value))
    for candidate in candidates:
        if candidate.exists():
            return os.path.abspath(str(candidate))
    notes.append(f"{name} {value} not found any more: left out")
    return None


def carry_params(run, defaults, notes):
    """The recorded parameters that matter in demix mode, as {name: value}."""
    recorded = run.params
    info = run.run_info or {}
    carried = {}
    for name in sorted(recorded):
        value = recorded[name]
        if value is None:
            continue
        if name not in CARRY_EXACT and not name.startswith(CARRY_PREFIXES):
            continue
        if run.target and name in PANEL_DROPPED:
            notes.append(f"{name} is not used with panel results: left out")
            continue
        if name in PATH_PARAMS:
            if name in ("fasta", "gff") and is_bundled(value):
                continue  # the pathogen config of the checkout supplies it
            if name == "lineage_groups":
                continue  # handled below (custom path, its published copy, or the bundled default)
            value = resolve_path_param(name, value, run, notes)
            if value is None:
                continue
        if defaults is not None:
            if name not in defaults:
                notes.append(f"{name} is not a parameter of this pipeline version: left out")
                continue
            if name not in ALWAYS_KEEP and (
                same_value(value, defaults[name])
                or (value is False and defaults[name] is None)  # a flag with no default is off
            ):
                continue
        carried[name] = value

    # primer scheme/BED: only for amplicon runs with trimming, and the BED (which wins) over the scheme
    capture = recorded.get("protocol") == "capture" or info.get("protocol") == "capture"
    untrimmed = (
        capture or recorded.get("skip_primer_trimming") or info.get("primer_trimming") is False
    )
    if untrimmed:
        carried.pop("primer_scheme", None)
        carried.pop("primer_bed", None)
    elif "primer_bed" in carried:
        carried.pop("primer_scheme", None)

    groups = carry_lineage_groups(run, notes)
    if groups:
        carried["lineage_groups"] = groups
    return carried


def carry_lineage_groups(run, notes):
    """A custom groups file: the recorded path while it exists, else the byte-identical copy published in the
    results (summary/lineage_groups_used.csv), else nothing (the bundled default applies, with a note)."""
    info = run.run_info or {}
    value = run.params.get("lineage_groups")
    source = info.get("lineage_groups_source")
    if not value and source != "custom":
        return None
    if value:
        for candidate in (Path(str(value)).expanduser(), run.path.parent / str(value)):
            if candidate.is_file():
                return os.path.abspath(str(candidate))
    used = run.summary_dir / "lineage_groups_used.csv"
    if used.is_file():
        want = info.get("lineage_groups_sha256")
        if not want or sha256_file(used) == want:
            notes.append(
                f"custom lineage_groups {value} not found: using the copy published in the results "
                "(summary/lineage_groups_used.csv)"
            )
            return str(used)
    notes.append(
        f"custom lineage_groups {value} not found and no copy in the results: the bundled default groups are used"
    )
    return None


def outdir_for(run, target_label, args):
    """(new outdir or None when it exists, note). `{barcode}` in the suffix is the barcode version."""
    label = re.sub(r"[^A-Za-z0-9._-]+", "-", str(target_label)).strip("-") or "new"
    suffix = args.outdir_suffix.replace("{barcode}", label)
    if run.target:  # panel results: one new directory per target
        suffix += "_" + re.sub(r"[^A-Za-z0-9._-]+", "-", run.target)
    base = str(run.path).rstrip("/") + suffix
    if not os.path.lexists(base):
        return base, None
    if not args.force_new_suffix:
        return None, f"output exists: {base} (use --force-new-suffix for a numbered copy)"
    counter = 2
    while os.path.lexists(f"{base}_{counter}"):
        counter += 1
    return f"{base}_{counter}", f"{base} exists: using a numbered copy"


def panel_file_of(run, notes):
    """The custom panel file the target's run recorded (`panel` in its run_info.json), as an absolute path,
    or None: nothing recorded, the panel came from results, it is a bundled panel, or the file is gone."""
    value = (run.run_info or {}).get("panel") or run.params.get("panel")
    text = str(value or "")
    if not text or not (text.endswith(".csv") or "/" in text) or "assets/panels/" in text:
        return None
    path = resolve_path_param("panel", text, run, notes)
    if path is None:
        notes.append(
            "the targets are looked up in the registry (assets/panels/registry.csv) instead; "
            "a target that is not in it stops the run (--extra-args '--panel <csv>' supplies a panel)"
        )
        return None
    try:
        with open(path, newline="", encoding="utf-8-sig") as handle:
            for row in csv.DictReader(handle):
                if (row.get("target") or "").strip() == run.target and (
                    row.get("barcodes") or ""
                ).strip():
                    notes.append(
                        f"the panel file gives {run.target} its own barcodes, which beat --barcodes"
                    )
    except (OSError, csv.Error, UnicodeDecodeError):
        pass
    return path


class Plan:
    """What would happen to one run."""

    def __init__(self, run):
        self.run = run
        self.action = "reanalyse"  # reanalyse | up to date | skip
        self.reason = ""
        self.outdir = None
        self.params = {}
        self.params_path = None
        self.argv = []
        self.notes = []


def build_plan(run, target, args, defaults):
    plan = Plan(run)
    if run.skip_reason:
        plan.action, plan.reason = "skip", run.skip_reason
        return plan
    if target.is_current(run):
        plan.action, plan.reason = "up to date", f"already used {run.barcode}"
        return plan
    outdir, note = outdir_for(run, target.label, args)
    if outdir is None:
        plan.action, plan.reason = "skip", note
        return plan
    if note:
        plan.notes.append(note)
    profile, profile_note = profile_for(run, args)
    if profile_note:
        plan.notes.append(profile_note)
    if run.kind == "legacy":
        plan.notes.append(
            "legacy results: no parameters are recorded, pipeline defaults apply (the Snakemake pipeline "
            "filtered variants at --varthresh 0.25 and that stays in the variants files); "
            "the date shown is the directory's modification time"
        )
    if run.success is None and run.kind == "wasteline" and "success" not in (run.run_info or {}):
        plan.notes.append(
            "the recorded run has no success flag: it may have crashed before completion"
        )
    plan.outdir = outdir
    plan.params = carry_params(run, defaults, plan.notes) if run.kind == "wasteline" else {}
    plan.params_path = outdir + ".params.json"
    if target.pending:
        plan.notes.append(
            "the directory name will use the resolved barcode version (latest is resolved when run)"
        )
    argv = [
        args.nextflow,
        "run",
        args.pipeline_dir,
        "-profile",
        profile,
        "-params-file",
        plan.params_path,
        "--mode",
        "demix",
        "--previous_results",
        str(run.path),
    ]
    if run.target:  # panel results: this target only (the others are separate runs)
        argv += ["--targets", run.target]
        panel_file = panel_file_of(run, plan.notes)
        if panel_file:
            argv += ["--panel", panel_file]
    argv += ["--outdir", outdir, "--barcodes", target.spec]
    if args.db_dir:
        argv += ["--db_dir", args.db_dir]
    argv += shlex.split(args.extra_args or "")
    plan.argv = argv
    return plan


# =========================================================================================================
# the diff
# =========================================================================================================

_NP_FLOAT = re.compile(r"(?:np|numpy)\.float\d*\(([^)]*)\)")
_NUMBER = re.compile(r"[-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?|nan|inf", re.IGNORECASE)


def parse_summarized(text):
    """`[('NB.1.8.1* (NB.1.8.1.X)', 0.908), ('Other', 0.0028)]` -> {name: abundance}. Tolerates numpy scalars
    (`np.float64(0.9)`), a dict repr, double quotes and trailing junk; empty/NaN -> {}."""
    text = (text or "").strip()
    if not text or text.lower() in ("nan", "none", "[]", "{}"):
        return {}
    cleaned = _NP_FLOAT.sub(r"\1", text)
    try:
        value = ast.literal_eval(cleaned)
        pairs = value.items() if isinstance(value, dict) else value
        return {str(name): float(number) for name, number in pairs}
    except (ValueError, SyntaxError, TypeError):
        pass
    result = {}
    for match in re.finditer(
        r"""\(\s*(['"])(.*?)\1\s*,\s*([-+0-9.eE]+|nan)\s*\)""", cleaned, flags=re.DOTALL
    ):
        result[match.group(2)] = float(match.group(3))
    return result


def parse_names(text):
    """Lineage names: space separated (Freyja's own form) or a python/numpy list repr."""
    text = re.sub(r"dtype\s*=\s*[\w'\"]+", " ", text or "")
    text = re.sub(r"\barray\b", " ", text)
    return [t for t in re.split(r"[\s,\[\]()'\"]+", text) if t]


def parse_numbers(text):
    text = _NP_FLOAT.sub(r"\1", text or "")
    text = re.sub(r"dtype\s*=\s*[\w'\"]+", " ", text)
    return [float(m.group(0)) for m in _NUMBER.finditer(text)]


def parse_lineages(lineages, abundances):
    """({lineage: abundance}, problem or None)."""
    if (lineages or "").strip().lower() in ("", "nan", "none") and (
        abundances or ""
    ).strip().lower() in (
        "",
        "nan",
        "none",
    ):
        return {}, None
    names = parse_names(lineages)
    numbers = parse_numbers(abundances)
    if len(names) != len(numbers):
        return {}, f"{len(names)} lineages but {len(numbers)} abundances"
    return dict(zip(names, numbers)), None


STRIP_SUFFIXES = (".freyja_variants.tsv", ".variants.tsv", ".demix.tsv", ".demix", ".tsv")


def clean_label(label):
    """The aggregate is labelled by file name in raw Freyja output: strip the suffix of the variants/demix file."""
    for suffix in STRIP_SUFFIXES:
        if label.endswith(suffix) and len(label) > len(suffix):
            return label[: -len(suffix)]
    return label


def read_table(path):
    """(header, [row lists]) of a TSV. Quoted cells may span lines: pandas writes numpy list reprs
    (`"['A' 'B'\n 'C']"`) that way. A quote in the middle of a cell is literal."""
    with open(str(path), encoding="utf-8-sig", newline="") as handle:
        rows = [row for row in csv.reader(handle, delimiter="\t") if any(c.strip() for c in row)]
    if not rows:
        return [], []
    return rows[0], rows[1:]


def read_aggregated(path, warnings):
    """{sample: {"lineages": {..}, "summarized": {..}}} from freyja_aggregated.tsv."""
    header, rows = read_table(path)
    if not header:
        return {}
    index = {name: i for i, name in enumerate(header) if name}

    def cell(row, name):
        i = index.get(name)
        return row[i] if i is not None and i < len(row) else ""

    result = {}
    for row in rows:
        if not row or not row[0].strip():
            continue
        sample = clean_label(row[0].strip())
        lineages, problem = parse_lineages(cell(row, "lineages"), cell(row, "abundances"))
        if problem:
            warnings.append(f"{path}: {sample}: {problem}; its lineages are ignored")
        if sample in result:
            warnings.append(f"{path}: duplicate sample {sample!r}: the last row is used")
        result[sample] = {
            "lineages": lineages,
            "summarized": parse_summarized(cell(row, "summarized")),
        }
    return result


def read_group_table(path, warnings):
    """({sample: {group column: abundance}}, [group columns]) from summary/lineage_summary.csv. The group
    columns are those between `genome_coverage_pct` and `proportion_summary` (the groups and `other`)."""
    with open(str(path), encoding="utf-8-sig", newline="") as handle:
        rows = list(csv.reader(handle))
    if not rows:
        return {}, []
    header = rows[0]
    if "sample" not in header or "proportion_summary" not in header:
        warnings.append(f"{path}: not a lineage_summary.csv (no sample/proportion_summary column)")
        return {}, []
    start = header.index("genome_coverage_pct") + 1 if "genome_coverage_pct" in header else None
    if start is None:
        start = header.index("status") + 1 if "status" in header else 1
    columns = header[start : header.index("proportion_summary")]
    sample_at = header.index("sample")
    result = {}
    for row in rows[1:]:
        if not row or sample_at >= len(row) or not row[sample_at].strip():
            continue
        values = {}
        for offset, column in enumerate(columns):
            at = start + offset
            text = row[at].strip() if at < len(row) else ""
            try:
                values[column] = float(text) if text else 0.0
            except ValueError:
                values[column] = 0.0
        result[row[sample_at].strip()] = values
    return result, columns


def read_status(path):
    """{sample: status} from summary/demix_status.tsv (empty when absent)."""
    if not path.is_file():
        return {}
    header, rows = read_table(path)
    if "sample" not in header or "status" not in header:
        return {}
    s, t = header.index("sample"), header.index("status")
    return {row[s]: row[t] for row in rows if len(row) > max(s, t)}


def run_versions(directory, target=""):
    """{"barcode": .., "freyja": ..} of a results directory: run_info.json, else summary/barcode_version.json
    (panel results: pipeline_info/<target>/ and summary/<target>/)."""
    info_path = directory / "pipeline_info" / target / "run_info.json"
    if info_path.is_file():
        try:
            info = load_json(info_path)
            shown, _iso, sha, source = recorded_barcode(info)
            return {
                "barcode": shown,
                "freyja": str(info.get("freyja_version") or "unknown"),
                "sha": sha,
                "source": source,
                "mode": info.get("mode"),
            }
        except (OSError, ValueError):
            pass
    bv = directory / "summary" / target / "barcode_version.json"
    if bv.is_file():
        try:
            data = load_json(bv)
            iso = data.get("barcode_version_iso") or to_iso(data.get("barcode_version"))
            return {
                "barcode": iso or data.get("barcode_version") or "unknown",
                "freyja": str(data.get("freyja_version") or "unknown"),
                "sha": data.get("barcodes_sha256"),
                "source": data.get("source"),
                "mode": None,
            }
        except (OSError, ValueError):
            pass
    return {"barcode": "unknown", "freyja": "unknown", "sha": None, "source": None, "mode": None}


class Results:
    """The tables of one results directory that the diff needs."""

    def __init__(self, directory, target=""):
        self.dir = Path(directory)
        self.target = target
        self.warnings = []
        if not self.dir.is_dir():
            raise Problem(f"not a directory: {directory}")
        summary = self.dir / "summary" / target
        self.aggregated_path = None
        candidates = ["summary/freyja_aggregated.tsv", "freyja_aggregated.tsv"]
        if target:
            candidates = [f"summary/{target}/freyja_aggregated.tsv"]
        for rel in candidates:
            if (self.dir / rel).is_file():
                self.aggregated_path = self.dir / rel
                break
        if self.aggregated_path is None:
            raise Problem(
                f"{directory} has no summary/freyja_aggregated.tsv: is it a finished results directory "
                "(a legacy Snakemake run may keep its aggregate elsewhere; copy or link it to summary/)?"
            )
        try:
            self.aggregated = read_aggregated(self.aggregated_path, self.warnings)
        except (OSError, UnicodeDecodeError) as err:
            raise Problem(f"cannot read {self.aggregated_path}: {err}")
        self.groups, self.group_columns = {}, []
        group_path = summary / "lineage_summary.csv"
        self.has_groups = group_path.is_file()
        if self.has_groups:
            try:
                self.groups, self.group_columns = read_group_table(group_path, self.warnings)
            except (OSError, csv.Error, UnicodeDecodeError) as err:
                self.warnings.append(f"cannot read {group_path}: {err}")
                self.has_groups = False
        self.status = read_status(summary / "demix_status.tsv")
        self.versions = run_versions(self.dir, target)


def is_flagged(change, threshold):
    """|change| >= threshold, tolerant of float noise (0.35 - 0.30 is 0.04999999999999999)."""
    return round(abs(change), 9) >= round(threshold, 9)


def diff_results(old, new, threshold, level="all"):
    """Compare two Results. Returns a dict: rows, summary numbers, notes."""
    if level == "all":
        levels = ["lineages", "summarized"]
        if old.has_groups and new.has_groups:
            levels.append("groups")
    else:
        levels = [level]
        if level == "groups" and not (old.has_groups and new.has_groups):
            missing = [r.dir for r in (old, new) if not r.has_groups]
            raise Problem(
                "--level groups needs summary/lineage_summary.csv in both runs; missing in "
                + ", ".join(str(m) for m in missing)
            )
    notes = []
    if level == "all" and not (old.has_groups and new.has_groups):
        notes.append("group level skipped: summary/lineage_summary.csv is not present in both runs")
    if "groups" in levels and old.group_columns != new.group_columns:
        only_old = [c for c in old.group_columns if c not in new.group_columns]
        only_new = [c for c in new.group_columns if c not in old.group_columns]
        notes.append(
            "the lineage group columns differ (different groups file?): "
            f"only in old: {only_old or 'none'}; only in new: {only_new or 'none'}; "
            "a missing group counts as 0"
        )

    samples_old, samples_new = set(old.aggregated), set(new.aggregated)
    if "groups" in levels and level == "groups":
        samples_old, samples_new = set(old.groups), set(new.groups)
    appear = sorted(samples_new - samples_old)
    disappear = sorted(samples_old - samples_new)
    status_changes = []
    for sample in sorted(set(old.status) | set(new.status)):
        before, after = old.status.get(sample, "absent"), new.status.get(sample, "absent")
        if before != after:
            status_changes.append((sample, before, after))
    changed_status = {s: (b, a) for s, b, a in status_changes}

    rows = []
    for sample in sorted(samples_old | samples_new):
        note = ""
        if sample in appear:
            note = "sample only in the new run"
        elif sample in disappear:
            note = "sample missing in the new run"
        if sample in changed_status:
            before, after = changed_status[sample]
            note = (note + "; " if note else "") + f"status {before} -> {after}"
        for lvl in levels:
            old_items = items_of(old, sample, lvl)
            new_items = items_of(new, sample, lvl)
            sample_rows = []
            for item in sorted(set(old_items) | set(new_items)):
                before, after = old_items.get(item, 0.0), new_items.get(item, 0.0)
                if before < MIN_ABUNDANCE and after < MIN_ABUNDANCE:
                    continue
                change = after - before
                sample_rows.append(
                    {
                        "sample": sample,
                        "level": lvl,
                        "item": item,
                        "old_abundance": before,
                        "new_abundance": after,
                        "change": change,
                        "flagged": is_flagged(change, threshold),
                        "note": note,
                    }
                )
            sample_rows.sort(key=lambda r: (-abs(r["change"]), r["item"]))
            rows.extend(sample_rows)

    compared = sorted(samples_old & samples_new)
    flagged_samples = sorted({r["sample"] for r in rows if r["flagged"]})
    top = sorted(
        (r for r in rows if r["flagged"]),
        key=lambda r: (-abs(r["change"]), r["sample"], r["level"], r["item"]),
    )
    for result in (old, new):
        notes.extend(result.warnings)
    if (
        old.versions["barcode"] == new.versions["barcode"]
        and old.versions["freyja"] == new.versions["freyja"]
        and old.versions["barcode"] != "unknown"
    ):
        notes.append(
            "same barcode version and Freyja version in both runs: differences are not expected"
        )
    return {
        "rows": rows,
        "levels": levels,
        "threshold": threshold,
        "compared": compared,
        "flagged_samples": flagged_samples,
        "top": top,
        "appear": appear,
        "disappear": disappear,
        "status_changes": status_changes,
        "old": old.versions,
        "new": new.versions,
        "notes": notes,
    }


def items_of(results, sample, level):
    if level == "groups":
        return results.groups.get(sample, {})
    return results.aggregated.get(sample, {}).get(level, {})


CSV_COLUMNS = [
    "sample",
    "level",
    "item",
    "old_abundance",
    "new_abundance",
    "change",
    "flagged",
    "note",
]


def csv_row(row):
    return [
        row["sample"],
        row["level"],
        row["item"],
        f"{row['old_abundance']:.6f}",
        f"{row['new_abundance']:.6f}",
        f"{row['change']:+.6f}",
        "yes" if row["flagged"] else "no",
        row["note"],
    ]


def write_csv(path, result_sets, label_column="run"):
    """`result_sets` is [(run label or None, diff dict)]; a `run` column (`label_column`) is added when a label is given."""
    with_run = any(label is not None for label, _ in result_sets)
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    with open(str(path), "w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(([label_column] if with_run else []) + CSV_COLUMNS)
        for label, result in result_sets:
            for row in result["rows"]:
                writer.writerow(([label or ""] if with_run else []) + csv_row(row))


def capped(items, limit=12):
    """Comma-separated list, with `(+N more)` when long."""
    items = list(items)
    text = ", ".join(items[:limit])
    return text + (f" (+{len(items) - limit} more)" if len(items) > limit else "")


def summary_lines(result, top_n=10):
    old, new = result["old"], result["new"]
    lines = [
        "barcode version: {} -> {}".format(old["barcode"], new["barcode"]),
        "Freyja version:  {} -> {}".format(old["freyja"], new["freyja"]),
        "samples compared: {}; with a change >= {:g}: {} (levels: {})".format(
            len(result["compared"]),
            result["threshold"],
            len(result["flagged_samples"]),
            ", ".join(result["levels"]),
        ),
    ]
    if result["appear"]:
        lines.append("samples only in the new run: " + capped(result["appear"]))
    if result["disappear"]:
        lines.append("samples missing in the new run: " + capped(result["disappear"]))
    if result["status_changes"]:
        lines.append(
            "status changes: "
            + capped([f"{s} {b} -> {a}" for s, b, a in result["status_changes"]], 12)
        )
    if result["top"]:
        lines.append(
            f"largest changes (top {min(top_n, len(result['top']))} of {len(result['top'])}):"
        )
        for row in result["top"][:top_n]:
            lines.append(
                "  {:<14} {:<10} {:<34} {:.4f} -> {:.4f} ({:+.4f})".format(
                    row["sample"],
                    row["level"],
                    row["item"][:34],
                    row["old_abundance"],
                    row["new_abundance"],
                    row["change"],
                )
            )
    else:
        lines.append("no change at or above the threshold")
    for note in result["notes"]:
        lines.append("note: " + note)
    return lines


def markdown_summary(title, result, top_n=20):
    old, new = result["old"], result["new"]
    out = [
        f"## {title}",
        "",
        f"- Barcode version: `{old['barcode']}` -> `{new['barcode']}`",
        f"- Freyja version: `{old['freyja']}` -> `{new['freyja']}`",
        "- Samples compared: {}; with a change of at least {:g}: {}".format(
            len(result["compared"]), result["threshold"], len(result["flagged_samples"])
        ),
        f"- Levels: {', '.join(result['levels'])}",
    ]
    if result["appear"]:
        out.append("- Only in the new run: " + ", ".join(result["appear"]))
    if result["disappear"]:
        out.append("- Missing in the new run: " + ", ".join(result["disappear"]))
    if result["status_changes"]:
        out.append(
            "- Status changes: "
            + capped([f"{s} {b} -> {a}" for s, b, a in result["status_changes"]], 12)
        )
    for note in result["notes"]:
        out.append(f"- Note: {note}")
    out.append("")
    if result["top"]:
        out += [
            f"Largest changes ({min(top_n, len(result['top']))} of {len(result['top'])}):",
            "",
            "| sample | level | item | old | new | change |",
            "| --- | --- | --- | ---: | ---: | ---: |",
        ]
        for row in result["top"][:top_n]:
            out.append(
                "| {} | {} | {} | {:.4f} | {:.4f} | {:+.4f} |".format(
                    row["sample"],
                    row["level"],
                    row["item"].replace("|", "/"),
                    row["old_abundance"],
                    row["new_abundance"],
                    row["change"],
                )
            )
    else:
        out.append("No change at or above the threshold.")
    out.append("")
    return "\n".join(out)


def run_diff(old_dir, new_dir, args, target=""):
    """Compare two results directories (of one target of panel results); returns the diff dict."""
    old, new = Results(old_dir, target), Results(new_dir, target)
    return diff_results(old, new, args.threshold, args.level)


# =========================================================================================================
# execution
# =========================================================================================================


def print_table(rows, header):
    widths = [max(len(str(r[i])) for r in [header] + rows) for i in range(len(header))]
    for r in [header] + rows:
        print("  " + "  ".join(str(c).ljust(widths[i]) for i, c in enumerate(r)).rstrip())


def table_rows(plans):
    rows = []
    for plan in plans:
        run = plan.run
        action = plan.action + (f" ({plan.reason})" if plan.reason else "")
        rows.append(
            [
                run.label,
                run.date or "?",
                run.mode,
                run.barcode,
                run.freyja_version,
                str(run.n_samples),
                action,
            ]
        )
    return rows


def run_update_db(args, profile):
    """Run `--mode update_db` once. Returns the iso stamp of the stored dataset."""
    tmp = tempfile.mkdtemp(prefix="wasteline_update_db_")
    argv = [
        args.nextflow,
        "run",
        args.pipeline_dir,
        "-profile",
        profile,
        "--mode",
        "update_db",
        "--outdir",
        tmp,
    ]
    if args.db_dir:
        argv += ["--db_dir", args.db_dir]
    argv += shlex.split(args.extra_args or "")
    print("Resolving 'latest': downloading the current Freyja data once for all runs")
    print(quote_command(argv))
    print()
    sys.stdout.flush()
    try:
        status = subprocess.call(argv)
        if status != 0:
            raise Problem(f"update_db failed (exit status {status}); nothing was re-analysed")
        info_path = Path(tmp) / "pipeline_info" / "freyja_data.json"
        if not info_path.is_file():
            raise Problem(f"update_db wrote no {info_path}: cannot tell which dataset was stored")
        info = load_json(info_path)
        iso = info.get("barcode_version_iso") or to_iso(info.get("barcode_version"))
        if not iso:
            raise Problem(f"{info_path} has no barcode_version_iso")
        return iso
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


def execute_plan(plan, args):
    """Write the params file, run Nextflow. Returns the exit status."""
    Path(plan.params_path).parent.mkdir(parents=True, exist_ok=True)
    with open(plan.params_path, "w") as handle:
        json.dump(plan.params, handle, indent=2, sort_keys=True)
        handle.write("\n")
    sys.stdout.flush()
    return subprocess.call(plan.argv)


def print_plan(index, total, plan, execute):
    print(f"[{index}/{total}] {plan.run.label}  ->  {plan.outdir}")
    carried = ", ".join(f"{k}={v}" for k, v in sorted(plan.params.items())) or "(none)"
    print(f"    params file {plan.params_path}{'' if execute else ' (written with --execute)'}")
    print(f"    carried over: {carried}")
    for note in plan.notes:
        print(f"    NOTE: {note}")
    print()
    print(quote_command(plan.argv))
    print()


def reanalyse(args):
    pipeline = args.pipeline_dir
    if Path(pipeline).is_dir():
        if not (Path(pipeline) / "main.nf").is_file():
            raise Problem(
                f"{pipeline} is not a wasteline checkout (no main.nf): use --pipeline-dir <clone>"
            )
        args.pipeline_dir = os.path.abspath(pipeline)
        defaults = schema_defaults(args.pipeline_dir)
    elif re.fullmatch(r"[\w.-]+/[\w.-]+", pipeline):
        defaults = None  # a GitHub repository: the schema is not at hand
    else:
        raise Problem(f"--pipeline-dir {pipeline!r} does not exist")
    if args.since and not re.fullmatch(r"\d{4}-\d{2}-\d{2}", args.since):
        raise Problem("--since must be YYYY-MM-DD")
    if not args.outdir_suffix.strip():
        raise Problem(
            "--outdir-suffix must not be empty (the new directory must differ from the run)"
        )
    target = resolve_target(args.barcodes, args.db_dir)

    runs, messages = discover(args.roots)
    for message in messages:
        print("WARNING: " + message, file=sys.stderr)
    if not runs:
        raise Problem("no results directories found")
    selected = select_runs(runs, args)
    if not selected:
        print("No run matches the selection (" + ", ".join(r.label for r in runs) + " were found).")
        return 0

    todo = [build_plan(r, target, args, defaults) for r in selected]
    update_profile = next(
        (profile_for(p.run, args)[0] for p in todo if p.action == "reanalyse"), "docker"
    )

    print("Target barcodes: " + target.describe())
    print()
    print_table(
        table_rows(todo), ["run", "date", "mode", "barcodes", "freyja", "samples", "action"]
    )
    print()

    if target.kind == "latest" and any(p.action == "reanalyse" for p in todo):
        if args.execute:
            iso = run_update_db(args, update_profile)
            print(f"Latest dataset: {iso}; every run uses --barcodes {iso}\n")
            target = resolve_target(iso, args.db_dir)
            todo = [build_plan(r, target, args, defaults) for r in selected]
            print_table(
                table_rows(todo), ["run", "date", "mode", "barcodes", "freyja", "samples", "action"]
            )
            print()
        else:
            dry = [
                args.nextflow,
                "run",
                args.pipeline_dir,
                "-profile",
                update_profile,
                "--mode",
                "update_db",
                "--outdir",
                "/tmp/wasteline_update_db_XXXXXX",
            ] + (["--db_dir", args.db_dir] if args.db_dir else [])
            print(
                "--barcodes latest: with --execute this runs first (once, then every run uses that stamp):"
            )
            print(quote_command(dry))
            print()

    pending = [p for p in todo if p.action == "reanalyse"]
    if not pending:
        print("Nothing to re-analyse.")
        return 0
    if not args.execute:
        for i, plan in enumerate(pending, 1):
            print_plan(i, len(pending), plan, False)
        print(f"Dry run: {len(pending)} command(s) NOT executed. Add --execute to run them.")
        return 0

    failed, done, combined = [], [], []
    for i, plan in enumerate(pending, 1):
        print_plan(i, len(pending), plan, True)
        status = execute_plan(plan, args)
        if status != 0:
            failed.append(plan.run.label)
            print(f"FAILED: {plan.run.label} (exit status {status})", file=sys.stderr)
            if not args.keep_going:
                print("Stopping (use --keep-going to continue with the next run).", file=sys.stderr)
                break
            continue
        done.append(plan.run.label)
        try:
            result = run_diff(plan.run.path, plan.outdir, args, plan.run.target)
        except Problem as err:
            print(f"NOTE: no comparison made: {err}", file=sys.stderr)
            continue
        new_summary = Path(plan.outdir) / "summary" / plan.run.target
        out_csv = new_summary / "reanalysis_diff.csv"
        write_csv(out_csv, [(None, result)])
        (new_summary / "reanalysis_diff.md").write_text(
            markdown_summary(f"{plan.run.label}: re-analysis with new barcodes", result)
        )
        print(f"Comparison written to {out_csv}")
        print("\n".join("  " + line for line in summary_lines(result)))
        print()
        combined.append((plan.run.label, result))
    if combined and args.report:
        write_csv(args.report, combined)
        print(f"Combined comparison: {args.report}")
    if combined and args.summary_md:
        Path(args.summary_md).write_text(
            "# Re-analysis summary\n\n"
            + "\n".join(markdown_summary(name, result) for name, result in combined)
        )
        print(f"Combined summary: {args.summary_md}")
    print(
        "Done: {} re-analysed, {} failed, {} up to date, {} skipped.".format(
            len(done),
            len(failed),
            sum(p.action == "up to date" for p in todo),
            sum(p.action == "skip" for p in todo),
        )
    )
    return 1 if failed else 0


def diff_only(args):
    old_dir, new_dir = args.diff
    old_targets, new_targets = panel_result_targets(old_dir), panel_result_targets(new_dir)
    if old_targets or new_targets:
        return diff_panel(args, old_dir, new_dir, old_targets, new_targets)
    result = run_diff(old_dir, new_dir, args)
    print(f"Comparison of {old_dir} (old) with {new_dir} (new)")
    print("\n".join(summary_lines(result)))
    if args.report:
        write_csv(args.report, [(None, result)])
        print(f"Report: {args.report}")
    if args.summary_md:
        Path(args.summary_md).write_text(
            markdown_summary(f"{Path(new_dir).name} vs {Path(old_dir).name}", result)
        )
        print(f"Summary: {args.summary_md}")
    return 0


def diff_panel(args, old_dir, new_dir, old_targets, new_targets):
    """--diff of panel results: one comparison per target present in both (never mixed)."""
    if not (old_targets and new_targets):
        raise Problem(
            "one of the two directories has panel results (summary/<target>/) and the other has not: "
            "nothing to compare"
        )
    common = [t for t in old_targets if t in new_targets]
    if not common:
        raise Problem(
            f"no target in common: {old_dir} has {', '.join(old_targets)}; {new_dir} has "
            + ", ".join(new_targets)
        )
    print(f"Comparison of {old_dir} (old) with {new_dir} (new), per target")
    for target in sorted(set(old_targets) ^ set(new_targets)):
        side = "old" if target in old_targets else "new"
        print(f"NOTE: target {target} is only in the {side} results: not compared")
    results = []
    for target in common:
        result = run_diff(old_dir, new_dir, args, target)
        results.append((target, result))
        print(f"\ntarget {target}")
        print("\n".join("  " + line for line in summary_lines(result)))
    if args.report:
        write_csv(args.report, results, "target")
        print(f"Report: {args.report}")
    if args.summary_md:
        Path(args.summary_md).write_text(
            "\n".join(
                markdown_summary(f"{Path(new_dir).name} vs {Path(old_dir).name}: {t}", r)
                for t, r in results
            )
        )
        print(f"Summary: {args.summary_md}")
    return 0


def parse_args(argv):
    here = Path(__file__).resolve().parent.parent
    parser = argparse.ArgumentParser(
        description="Re-run past wasteline results with newer Freyja barcodes (no re-mapping) and "
        "compare the lineage abundances. Dry run unless --execute.",
        epilog="Exit status: 0 ok, 1 error.",
    )
    parser.add_argument(
        "roots", nargs="*", help="results directories, or folders that contain them (depth <= 2)"
    )
    sel = parser.add_argument_group("which runs")
    sel.add_argument("--since", metavar="YYYY-MM-DD", help="runs started on or after this date")
    sel.add_argument("--last", type=int, metavar="N", help="the N most recent re-analysable runs")
    sel.add_argument(
        "--runs", metavar="a,b", help="run names (directory names or paths), comma separated"
    )
    parser.add_argument(
        "--barcodes",
        default="latest",
        help="target dataset: latest (default), YYYY-MM-DD, a stamp, container, or a path",
    )
    parser.add_argument(
        "--pipeline-dir",
        default=str(here),
        help="wasteline checkout to run (default: the clone this script lives in)",
    )
    parser.add_argument(
        "--outdir-suffix",
        default=DEFAULT_SUFFIX,
        help="suffix of the new results directory; {barcode} is the barcode version "
        f"(default: {DEFAULT_SUFFIX})",
    )
    parser.add_argument(
        "--force-new-suffix",
        action="store_true",
        help="when the new directory exists, add a counter (_2, _3, ...) instead of skipping the run",
    )
    parser.add_argument("--nextflow", default="nextflow", help="the nextflow executable")
    parser.add_argument(
        "--extra-args",
        default="",
        help="extra Nextflow arguments (one quoted string), e.g. '-work-dir /scratch/w'",
    )
    parser.add_argument(
        "--profile", help="-profile to use (default: the run's recorded one, without locked/test)"
    )
    parser.add_argument(
        "--db_dir", "--db-dir", dest="db_dir", help="Freyja dataset store (default ~/.wasteline/db)"
    )
    parser.add_argument(
        "--execute", action="store_true", help="run the commands (default: dry run)"
    )
    parser.add_argument(
        "--keep-going", action="store_true", help="continue with the next run after a failure"
    )
    parser.add_argument(
        "--diff",
        nargs=2,
        metavar=("OLD", "NEW"),
        help="only compare two results directories (no Nextflow)",
    )
    parser.add_argument(
        "--threshold",
        type=float,
        default=0.05,
        help="flag a change of at least this abundance (default 0.05)",
    )
    parser.add_argument(
        "--level",
        choices=["lineages", "groups", "summarized", "all"],
        default="all",
        help="what to compare: specific lineages, lineage_summary.csv groups, Freyja's summarized "
        "groups, or all that exist in both runs (default)",
    )
    parser.add_argument("--report", metavar="CSV", help="write the comparison as CSV")
    parser.add_argument("--summary-md", metavar="MD", help="write a short Markdown summary")
    args = parser.parse_args(argv)
    if args.diff and args.roots:
        parser.error("--diff takes the two results directories itself: do not also give roots")
    if not args.diff and not args.roots:
        parser.error("give at least one results directory (or --diff OLD NEW)")
    if args.threshold < 0:
        parser.error("--threshold must be >= 0")
    if args.last is not None and args.last < 0:
        parser.error("--last must be >= 0")
    return args


def main(argv=None):
    args = parse_args(sys.argv[1:] if argv is None else argv)
    try:
        return diff_only(args) if args.diff else reanalyse(args)
    except Problem as err:
        print(f"error: {err}", file=sys.stderr)
        return 1
    except OSError as err:
        print(f"error: {err}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())
