#!/usr/bin/env python3
"""Check for, and apply, a new Freyja release across every pin in the repository (PLAN section 9).

Freyja is pinned in three spellings, all discovered by grep (the module list is never hard-coded):

    quay.io/biocontainers/freyja:<version>--<build>_<n>              (docker / podman)
    depot.galaxyproject.org/singularity/freyja:<version>--<build>_<n> (singularity / apptainer)
    bioconda::freyja=<version>                                        (conda environment.yml)

Modes (exactly one):
    --check               print pinned vs latest. Exit 0 = up to date (or a newer release exists but is
                          not yet available on every channel: status "blocked"), 10 = update available,
                          2 = error. `--json` prints a machine-readable object (includes `major_bump`).
    --apply               rewrite every pin to the latest release (or `--tag TAG`) and append a
                          CHANGELOG.md line under "Unreleased / Dependencies". Idempotent.
    --require-consistent  offline invariant: every pin agrees (exit 0), else exit 1 (exit 2 = error).

The network is reached only through the injectable `Network` object, so everything is testable
offline against the recorded fixtures in tests/fixtures/ci/.
"""

from __future__ import annotations

import argparse
import json
import os
import re
import sys
import urllib.error
import urllib.request
from dataclasses import dataclass
from pathlib import Path

QUAY_TAGS_URL = (
    "https://quay.io/api/v1/repository/biocontainers/freyja/tag/?limit=100&onlyActiveTags=true"
)
ANACONDA_URL = "https://api.anaconda.org/package/bioconda/freyja"
DEPOT_URL = "https://depot.galaxyproject.org/singularity/freyja:{tag}"
MAX_QUAY_PAGES = 20

EXIT_UP_TO_DATE = 0
EXIT_UPDATE = 10
EXIT_ERROR = 2
EXIT_INCONSISTENT = 1

# Release tags only: <numeric version>--<build string>_<build number>. Anything else (latest, bare
# versions, rc/dev/post tags, arch suffixes) is not a release we would pin.
TAG_RE = re.compile(r"^(?P<version>\d+(?:\.\d+)*)--(?P<build>[A-Za-z0-9]+)_(?P<number>\d+)$")

CHARS = r"[A-Za-z0-9_.\-]+"
PIN_PATTERNS = {
    "quay": re.compile(r"(quay\.io/biocontainers/freyja:)(" + CHARS + ")"),
    "depot": re.compile(r"(depot\.galaxyproject\.org/singularity/freyja:)(" + CHARS + ")"),
    "conda": re.compile(r"(bioconda::freyja=)([0-9][A-Za-z0-9_.\-]*)"),
}
REQUIRED_KINDS = tuple(PIN_PATTERNS)
# Version stamps recorded by tests (`"freyja": "2.0.5"` in nf-test snapshots, `freyja_version ==
# '2.0.5'` or a synthetic `freyja: 2.0.5` in .nf.test files). Not pins, but they must move with the pin or every bump PR fails CI.
STAMP_PATTERN = re.compile(
    r"""(\bfreyja(?:_version)?"?\s*(?:[:=]|==)\s*["']?)(\d+(?:\.\d+)*)(["']?)"""
)
STAMP_SUFFIXES = (".snap", ".nf.test")
ALL_PATTERNS = {**PIN_PATTERNS, "stamp": STAMP_PATTERN}

IGNORE_MARKER = "bump-freyja: ignore"
SKIP_DIRS = {
    ".git", ".nf-test", ".pixi", ".venv", "work", "node_modules", "__pycache__",
    ".ruff_cache", ".pytest_cache", ".claude", ".nextflow", "results",
}  # fmt: skip
# Files that legitimately mention Freyja pins without being pins: history, generated digests,
# provenance records, and this tool's own tests/fixtures/docs.
SKIP_FILES = {
    "CHANGELOG.md",
    "pixi.lock",
    "conf/containers.lock.config",
    "docs/PLAN.md",
    "docs/test_data.md",
    "docs/ci.md",
    "bin/bump_freyja.py",
    "bin/lock_containers.py",
    "tests/python/test_bump_freyja.py",
    "tests/python/test_lock_containers.py",
}
SKIP_PREFIXES = ("tests/fixtures/ci/",)
TEXT_SUFFIXES = {
    ".nf", ".config", ".yml", ".yaml", ".md", ".py", ".test", ".json", ".txt", ".toml", ".sh",
    ".cfg", ".snap", ".j2",
}  # fmt: skip
MAX_BYTES = 2_000_000


class BumpError(Exception):
    """A condition that should end the run with exit code 2."""


# ---------------------------------------------------------------------------------------------
# Tags and versions
# ---------------------------------------------------------------------------------------------


@dataclass(frozen=True)
class Tag:
    raw: str
    version: str
    build: str
    number: int

    @property
    def version_key(self) -> tuple[int, ...]:
        return tuple(int(p) for p in self.version.split("."))

    @property
    def sort_key(self) -> tuple:
        # highest version, then highest build number; for a tie prefer the noarch (pyhdfd78af)
        # build, then the lexicographically larger build string (deterministic).
        return (self.version_key, self.number, self.build == "pyhdfd78af", self.build)

    @property
    def major(self) -> int:
        return self.version_key[0]


def parse_tag(raw: str) -> Tag | None:
    match = TAG_RE.match(raw)
    if not match:
        return None
    return Tag(raw, match["version"], match["build"], int(match["number"]))


def pick_latest(tag_names) -> Tag | None:
    """Highest release tag, ignoring every non-release tag."""
    tags = [t for t in (parse_tag(n) for n in tag_names) if t]
    return max(tags, key=lambda t: t.sort_key) if tags else None


def is_major_bump(pinned: Tag, latest: Tag) -> bool:
    return latest.major > pinned.major


# ---------------------------------------------------------------------------------------------
# Network (injectable)
# ---------------------------------------------------------------------------------------------


class Network:
    """The only place that touches the network. Replace in tests."""

    user_agent = "wasteline-bump-freyja"

    def get_json(self, url: str):
        request = urllib.request.Request(
            url, headers={"Accept": "application/json", "User-Agent": self.user_agent}
        )
        try:
            with urllib.request.urlopen(request, timeout=30) as response:
                return json.load(response)
        except (urllib.error.URLError, TimeoutError, json.JSONDecodeError) as err:
            raise BumpError(f"request failed: {url}: {err}") from err

    def head_status(self, url: str) -> int | None:
        """HTTP status of a HEAD request, or None if the server could not be reached."""
        request = urllib.request.Request(
            url, method="HEAD", headers={"User-Agent": self.user_agent}
        )
        try:
            with urllib.request.urlopen(request, timeout=30) as response:
                return response.status
        except urllib.error.HTTPError as err:
            return err.code
        except (urllib.error.URLError, TimeoutError):
            return None


def fetch_quay_tag_names(network) -> list[str]:
    names: list[str] = []
    for page in range(1, MAX_QUAY_PAGES + 1):
        url = QUAY_TAGS_URL if page == 1 else f"{QUAY_TAGS_URL}&page={page}"
        data = network.get_json(url)
        if not isinstance(data, dict) or "tags" not in data:
            raise BumpError(f"unexpected quay.io tag API response from {url}")
        names.extend(t["name"] for t in data["tags"] if isinstance(t, dict) and "name" in t)
        if not data.get("has_additional"):
            return names
    raise BumpError("quay.io tag listing did not terminate")


def bioconda_versions(network) -> list[str]:
    data = network.get_json(ANACONDA_URL)
    if not isinstance(data, dict) or not isinstance(data.get("versions"), list):
        raise BumpError("unexpected anaconda.org response (no `versions` list)")
    return [str(v) for v in data["versions"]]


# ---------------------------------------------------------------------------------------------
# Pin discovery
# ---------------------------------------------------------------------------------------------


@dataclass(frozen=True)
class Pin:
    kind: str  # quay | depot | conda
    path: str  # repo-relative, posix
    line: int
    value: str  # tag (quay/depot) or version (conda)


def iter_text_files(root: Path):
    """Yield (relative posix path, absolute path) of candidate files, deterministically sorted."""
    found = []
    for dirpath, dirnames, filenames in os.walk(root):
        dirnames[:] = sorted(d for d in dirnames if d not in SKIP_DIRS)
        for name in sorted(filenames):
            absolute = Path(dirpath) / name
            rel = absolute.relative_to(root).as_posix()
            if rel in SKIP_FILES or rel.startswith(SKIP_PREFIXES):
                continue
            if absolute.suffix not in TEXT_SUFFIXES:
                continue
            if absolute.is_symlink() or absolute.stat().st_size > MAX_BYTES:
                continue
            found.append((rel, absolute))
    yield from sorted(found)


def clean_value(kind: str, value: str) -> str:
    # prose such as "...pyhdfd78af_0." must not swallow the full stop
    return value.rstrip(".-") if kind != "conda" else value.rstrip(".")


def read_text(path: Path) -> str | None:
    try:
        return path.read_bytes().decode("utf-8")
    except UnicodeDecodeError:
        return None


def scan_pins(root: Path) -> list[Pin]:
    pins: list[Pin] = []
    for rel, absolute in iter_text_files(root):
        text = read_text(absolute)
        if text is None:
            continue
        for lineno, line in enumerate(text.splitlines(), start=1):
            if IGNORE_MARKER in line:
                continue
            for kind, pattern in ALL_PATTERNS.items():
                if kind == "stamp" and not rel.endswith(STAMP_SUFFIXES):
                    continue
                for match in pattern.finditer(line):
                    pins.append(Pin(kind, rel, lineno, clean_value(kind, match.group(2))))
    return pins


def conda_version_of(tag: str) -> str | None:
    parsed = parse_tag(tag)
    return parsed.version if parsed else None


def consistency_problems(pins: list[Pin]) -> list[str]:
    """Human-readable problems; an empty list means every pin agrees."""
    problems: list[str] = []
    by_kind: dict[str, dict[str, list[Pin]]] = {k: {} for k in ALL_PATTERNS}
    for pin in pins:
        by_kind[pin.kind].setdefault(pin.value, []).append(pin)

    for kind, values in by_kind.items():
        if not values and kind == "stamp":
            continue
        if not values:
            problems.append(f"no {kind} Freyja pin found anywhere")
        elif len(values) > 1:
            detail = "; ".join(
                f"{value} in {', '.join(sorted({p.path for p in locs}))}"
                for value, locs in sorted(values.items())
            )
            problems.append(f"{kind} pins disagree: {detail}")

    quay, depot, conda = (set(by_kind[k]) for k in ("quay", "depot", "conda"))
    if len(quay) == 1:
        (qtag,) = quay
        if parse_tag(qtag) is None:
            problems.append(f"quay pin {qtag!r} is not a release tag <version>--<build>_<n>")
        if len(depot) == 1 and depot != quay:
            problems.append(f"depot tag {next(iter(depot))!r} differs from quay tag {qtag!r}")
        if len(conda) == 1 and conda_version_of(qtag) not in conda:
            problems.append(
                f"conda version {next(iter(conda))!r} differs from the container version "
                f"{conda_version_of(qtag)!r}"
            )

    stamps = set(by_kind["stamp"])
    if stamps and len(conda) == 1 and stamps != conda:
        problems.append(
            f"test version stamp(s) {sorted(stamps)} differ from the conda version "
            f"{next(iter(conda))!r} (in {', '.join(sorted({p.path for p in by_kind['stamp'][next(iter(stamps))]}))})"
        )

    # every module that pins Freyja must pin it in all three places
    modules: dict[str, set[str]] = {}
    for pin in pins:
        parts = pin.path.split("/")
        if parts[0] == "modules" and len(parts) > 3 and pin.kind in REQUIRED_KINDS:
            modules.setdefault("/".join(parts[:3]), set()).add(pin.kind)
    for module, kinds in sorted(modules.items()):
        missing = sorted(set(REQUIRED_KINDS) - kinds)
        if missing:
            problems.append(f"module {module} lacks Freyja pin(s): {', '.join(missing)}")
    return problems


def current_tag(pins: list[Pin]) -> Tag:
    problems = consistency_problems(pins)
    if problems:
        raise BumpError("Freyja pins are inconsistent: " + " | ".join(problems))
    (raw,) = {p.value for p in pins if p.kind == "quay"}
    tag = parse_tag(raw)
    assert tag is not None  # guaranteed by consistency_problems
    return tag


# ---------------------------------------------------------------------------------------------
# --check
# ---------------------------------------------------------------------------------------------


def resolve_latest(network, check_depot: bool = True) -> dict:
    """Latest release on quay.io plus its availability on the other channels."""
    names = fetch_quay_tag_names(network)
    latest = pick_latest(names)
    if latest is None:
        raise BumpError("no release tags found on quay.io/biocontainers/freyja")
    conda_ok = latest.version in bioconda_versions(network)
    if check_depot:
        status = network.head_status(DEPOT_URL.format(tag=latest.raw))
        if status is None:
            depot = "skipped"  # depot unreachable: not evidence either way
        elif 200 <= status < 300:
            depot = "ok"
        elif status == 404:
            depot = "missing"
        else:
            depot = "skipped"
    else:
        depot = "skipped"
    return {"tag": latest, "bioconda": conda_ok, "depot": depot}


def check(root: Path, network, check_depot: bool = True) -> dict:
    pinned = current_tag(scan_pins(root))
    found = resolve_latest(network, check_depot)
    latest: Tag = found["tag"]
    newer = latest.sort_key > pinned.sort_key
    blockers = []
    if newer and not found["bioconda"]:
        blockers.append(f"bioconda has no freyja {latest.version} yet")
    if newer and found["depot"] == "missing":
        blockers.append(f"depot.galaxyproject.org has no singularity image {latest.raw} yet")
    if not newer:
        status = "up_to_date"
    elif blockers:
        status = "blocked"
    else:
        status = "update_available"
    return {
        "status": status,
        "update_available": status == "update_available",
        "pinned": pinned.raw,
        "pinned_version": pinned.version,
        "latest": latest.raw,
        "latest_version": latest.version,
        "major_bump": bool(newer and is_major_bump(pinned, latest)),
        "blocked_reasons": blockers,
        "checks": {
            "quay": "ok",
            "bioconda": "ok" if found["bioconda"] else "missing",
            "singularity_depot": found["depot"],
        },
    }


def render_check(result: dict) -> str:
    lines = [
        f"pinned : {result['pinned']}",
        f"latest : {result['latest']}",
        "checks : " + ", ".join(f"{name}={state}" for name, state in result["checks"].items()),
    ]
    if result["status"] == "up_to_date":
        lines.append("Freyja is up to date.")
    elif result["status"] == "blocked":
        lines.append("A newer release exists but is not available everywhere yet: ")
        lines.extend(f"  - {reason}" for reason in result["blocked_reasons"])
    else:
        kind = "MAJOR" if result["major_bump"] else "same-major"
        lines.append(
            f"Update available ({kind}): {result['pinned_version']} -> {result['latest_version']}"
        )
    return "\n".join(lines)


# ---------------------------------------------------------------------------------------------
# --apply
# ---------------------------------------------------------------------------------------------


def rewrite_text(text: str, new_tag: Tag, rel: str = "") -> str:
    out = []
    for line in text.splitlines(keepends=True):
        if IGNORE_MARKER not in line:
            for kind, pattern in ALL_PATTERNS.items():
                if kind == "stamp" and not rel.endswith(STAMP_SUFFIXES):
                    continue
                value = new_tag.raw if kind in ("quay", "depot") else new_tag.version

                def substitute(match, kind=kind, value=value):
                    old = clean_value(kind, match.group(2))
                    # keep any prose punctuation that clean_value stripped; group 3 = closing quote
                    tail = match.group(3) if match.re.groups >= 3 else ""
                    return match.group(1) + value + match.group(2)[len(old) :] + tail

                line = pattern.sub(substitute, line)
        out.append(line)
    return "".join(out)


UNRELEASED_RE = re.compile(r"^##\s+.*(unreleased|dev)", re.IGNORECASE)
DEPS_RE = re.compile(r"^###\s+`?Dependencies`?\s*$", re.IGNORECASE)


def changelog_line(old: Tag, new: Tag) -> str:
    return f"- Bump Freyja container `{old.raw}` -> `{new.raw}` (automated by `bin/bump_freyja.py`)"


def add_changelog_entry(text: str, old: Tag, new: Tag) -> str:
    """Append a bullet under the unreleased section's Dependencies heading (created if absent)."""
    bullet = changelog_line(old, new)
    lines = text.split("\n")
    if bullet in lines:
        return text
    trailing_newline = text.endswith("\n")
    if trailing_newline:
        lines = lines[:-1]

    def section_end(start: int) -> int:
        for i in range(start + 1, len(lines)):
            if lines[i].startswith("## "):
                return i
        return len(lines)

    unreleased = next((i for i, ln in enumerate(lines) if UNRELEASED_RE.match(ln)), None)
    if unreleased is None:
        first = next((i for i, ln in enumerate(lines) if ln.startswith("## ")), len(lines))
        block = ["## Unreleased", "", "### `Dependencies`", "", bullet, ""]
        if first == len(lines) and lines and lines[-1] != "":
            block = [""] + block
        lines[first:first] = block
    else:
        end = section_end(unreleased)
        heading = next((i for i in range(unreleased + 1, end) if DEPS_RE.match(lines[i])), None)
        if heading is None:
            stop = end
            while stop > unreleased + 1 and lines[stop - 1] == "":
                stop -= 1
            # blank lines that followed the section stay after the new block
            lines[stop:stop] = ["", "### `Dependencies`", "", bullet]
        else:
            i = heading + 1
            while i < end and lines[i] == "":
                i += 1
            had_bullets = False
            while i < end and (lines[i].startswith(("- ", "* ")) or lines[i].startswith("  ")):
                had_bullets = True
                i += 1
            if had_bullets:
                lines.insert(i, bullet)
            elif i >= len(lines):
                lines.extend(([""] if lines[-1] != "" else []) + [bullet])
            else:
                lines[i:i] = [bullet, ""]
    return "\n".join(lines) + ("\n" if trailing_newline else "")


def apply_tag(root: Path, new_tag: Tag, dry_run: bool = False) -> dict:
    old = current_tag_or_none(scan_pins(root))
    changed: list[str] = []
    for rel, absolute in iter_text_files(root):
        text = read_text(absolute)
        if text is None:
            continue
        updated = rewrite_text(text, new_tag, rel)
        if updated != text:
            changed.append(rel)
            if not dry_run:
                absolute.write_bytes(updated.encode("utf-8"))
    changelog = root / "CHANGELOG.md"
    changelog_changed = False
    if old is not None and old.raw != new_tag.raw:
        text = read_text(changelog) if changelog.exists() else ""
        updated = add_changelog_entry(text or "", old, new_tag)
        if updated != text:
            changelog_changed = True
            if not dry_run:
                changelog.write_bytes(updated.encode("utf-8"))
    if not dry_run:
        problems = consistency_problems(scan_pins(root))
        if problems:
            raise BumpError("pins still inconsistent after --apply: " + " | ".join(problems))
    return {
        "applied": bool(changed),
        "from": old.raw if old else None,
        "to": new_tag.raw,
        "files": changed,
        "changelog": changelog_changed,
        "dry_run": dry_run,
    }


def current_tag_or_none(pins: list[Pin]) -> Tag | None:
    quay = sorted({p.value for p in pins if p.kind == "quay"})
    tags = [t for t in (parse_tag(q) for q in quay) if t]
    return max(tags, key=lambda t: t.sort_key) if tags else None


# ---------------------------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------------------------


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--check", action="store_true", help="compare pinned vs latest")
    mode.add_argument("--apply", action="store_true", help="rewrite all pins + CHANGELOG")
    mode.add_argument(
        "--require-consistent", action="store_true", help="fail if any Freyja pin disagrees"
    )
    parser.add_argument("--json", action="store_true", help="machine-readable output")
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--tag", help="--apply: use this exact tag instead of looking up latest")
    parser.add_argument("--dry-run", action="store_true", help="--apply: report, do not write")
    parser.add_argument(
        "--no-depot-check",
        action="store_true",
        help="skip the depot.galaxyproject.org singularity existence check",
    )
    return parser


def emit(args, payload: dict, text: str) -> None:
    print(json.dumps(payload, indent=2, sort_keys=True) if args.json else text)


def main(argv=None, network=None) -> int:
    args = build_parser().parse_args(argv)
    network = network or Network()
    root = args.root.resolve()
    try:
        if args.require_consistent:
            pins = scan_pins(root)
            problems = consistency_problems(pins)
            payload = {
                "consistent": not problems,
                "problems": problems,
                "pins": [p.__dict__ for p in pins],
            }
            summary = (
                f"Freyja pins consistent ({len(pins)} pins)"
                if not problems
                else "Freyja pins INCONSISTENT:\n" + "\n".join(f"  - {p}" for p in problems)
            )
            emit(args, payload, summary)
            return EXIT_UP_TO_DATE if not problems else EXIT_INCONSISTENT

        if args.check:
            result = check(root, network, check_depot=not args.no_depot_check)
            emit(args, result, render_check(result))
            return EXIT_UPDATE if result["update_available"] else EXIT_UP_TO_DATE

        # --apply
        if args.tag:
            new_tag = parse_tag(args.tag)
            if new_tag is None:
                raise BumpError(f"--tag {args.tag!r} is not a release tag <version>--<build>_<n>")
        else:
            result = check(root, network, check_depot=not args.no_depot_check)
            if result["status"] != "update_available":
                emit(
                    args, {**result, "applied": False}, render_check(result) + "\nNothing applied."
                )
                return EXIT_UP_TO_DATE
            new_tag = parse_tag(result["latest"])
        outcome = apply_tag(root, new_tag, dry_run=args.dry_run)
        summary = (
            f"{'Would rewrite' if args.dry_run else 'Rewrote'} {len(outcome['files'])} file(s) "
            f"to {new_tag.raw}:\n" + "\n".join(f"  {f}" for f in outcome["files"])
            if outcome["files"]
            else f"All pins already at {new_tag.raw}; nothing to do."
        )
        emit(args, outcome, summary)
        return EXIT_UP_TO_DATE
    except BumpError as err:
        print(f"error: {err}", file=sys.stderr)
        return EXIT_ERROR


if __name__ == "__main__":
    sys.exit(main())
