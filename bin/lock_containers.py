#!/usr/bin/env python3
"""Pin every process container to an immutable digest (F13, PLAN section 10).

Reads the exact-pinned container references that processes use (module `main.nf` files plus any
`withName: 'X' { container = '...' }` override in `conf/*.config`), resolves each one to its
registry digest with a read-only manifest HEAD request, and writes `conf/containers.lock.config`:

    process {
        withName: 'FREYJA_DEMIX' { container = 'quay.io/biocontainers/freyja@sha256:...' } // pin: ...
    }

The file is loaded by `-profile locked` (snippet in docs/ci.md). Multi-arch images record the
manifest-LIST digest (whatever `Docker-Content-Digest` the tag resolves to).

    --update   (re)write the lock. Digests of unchanged pins are reused; only new/changed pins are
               resolved (use --refresh to re-resolve everything).
    --check    exit 0 if the lock matches the pins, 1 if stale, 2 on error. Offline unless
               --verify-digests is given (then every digest is re-resolved and compared).

The network is reached only through the injectable `Registry` object.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import sys
import urllib.error
import urllib.parse
import urllib.request
from dataclasses import dataclass
from pathlib import Path

DEFAULT_REGISTRY = "quay.io"  # docker.registry in nextflow.config
LOCK_FILE = "conf/containers.lock.config"
EXIT_OK, EXIT_STALE, EXIT_ERROR = 0, 1, 2

SKIP_CONF = {"containers.lock.config"}
SKIP_CONF_PREFIXES = ("containers_",)  # nf-core wave/seqera freeze configs, not loaded by us
DIGEST_RE = re.compile(r"^sha256:[0-9a-f]{64}$")
INDEX_TYPES = {
    "application/vnd.oci.image.index.v1+json",
    "application/vnd.docker.distribution.manifest.list.v2+json",
}
SINGLE_TYPES = (
    "application/vnd.oci.image.manifest.v1+json",
    "application/vnd.docker.distribution.manifest.v2+json",
)
ACCEPT = ", ".join((*sorted(INDEX_TYPES), *SINGLE_TYPES))


class LockError(Exception):
    """A condition that should end the run with exit code 2."""


# ---------------------------------------------------------------------------------------------
# Image references
# ---------------------------------------------------------------------------------------------


@dataclass(frozen=True)
class Ref:
    registry: str
    repository: str
    tag: str

    @property
    def pin(self) -> str:
        """Human-readable pinned reference, as written in the module."""
        return f"{self.registry}/{self.repository}:{self.tag}"

    @property
    def digest_name(self) -> str:
        return f"{self.registry}/{self.repository}"


def parse_ref(ref: str) -> Ref:
    """`[registry/]repo[:tag]` -> Ref. Refuses anything that is not an exact tag."""
    if "@" in ref:
        raise LockError(f"{ref!r} is already digest-pinned; pin modules by tag")
    name, sep, tag = ref.rpartition(":")
    if not sep or "/" in tag:  # no tag (a ':' before the last '/' belongs to a host:port)
        raise LockError(f"{ref!r} has no tag: containers must be pinned to an exact tag")
    if tag in ("latest", ""):
        raise LockError(f"{ref!r} uses a floating tag: containers must be pinned to an exact tag")
    first, _, rest = name.partition("/")
    if rest and ("." in first or ":" in first or first == "localhost"):
        return Ref(first, rest, tag)
    return Ref(DEFAULT_REGISTRY, name, tag)


# ---------------------------------------------------------------------------------------------
# Scanning the repository for the images processes use
# ---------------------------------------------------------------------------------------------

PROCESS_RE = re.compile(r"^\s*process\s+([A-Za-z_][A-Za-z0-9_]*)\s*\{", re.MULTILINE)
CONTAINER_DQ_RE = re.compile(r'^\s*container\s+"([^"]*)"', re.MULTILINE)
CONTAINER_SQ_RE = re.compile(r"^\s*container\s+'([^']*)'", re.MULTILINE)
QUOTED_RE = re.compile(r"'([^'\s]+)'")
OVERRIDE_RE = re.compile(
    r"withName\s*:\s*['\"]([A-Za-z0-9_]+)['\"]\s*\{[^{}]*?\bcontainer\s*=\s*['\"]([^'\"]+)['\"]",
    re.DOTALL,
)


def docker_image_of(expression: str) -> str | None:
    """The docker (non-URL) image named in a module `container` expression."""
    candidates = [
        c
        for c in QUOTED_RE.findall(expression)
        if not c.startswith(("http://", "https://", "docker://", "oras://"))
    ]
    if candidates:
        return candidates[-1]  # `cond ? 'https://depot...' : 'quay.io/...'` -> the else branch
    plain = expression.strip()
    if plain and "${" not in plain and " " not in plain and not plain.startswith(("http", "oras")):
        return plain
    return None


def module_images(root: Path) -> dict[str, str]:
    """process name -> image string, from modules/**/main.nf."""
    images: dict[str, str] = {}
    for main_nf in sorted((root / "modules").rglob("main.nf")):
        text = main_nf.read_text()
        processes = [(m.start(), m.group(1)) for m in PROCESS_RE.finditer(text)]
        if not processes:
            continue
        for regex in (CONTAINER_DQ_RE, CONTAINER_SQ_RE):
            for match in regex.finditer(text):
                owner = [name for pos, name in processes if pos < match.start()]
                if not owner:
                    continue
                image = docker_image_of(match.group(1))
                if image is None:
                    raise LockError(
                        f"{main_nf.relative_to(root)}: cannot find a docker image in "
                        f"container expression {match.group(1)!r}"
                    )
                images[owner[-1]] = image
    return images


INCLUDE_RE = re.compile(r"^\s*include\s*\{([^}]*)\}\s*from", re.MULTILINE)
ALIAS_RE = re.compile(r"\b([A-Za-z_][A-Za-z0-9_]*)\s+as\s+([A-Za-z_][A-Za-z0-9_]*)\b")


def process_aliases(root: Path, known: set[str]) -> dict[str, str]:
    """alias -> process name for `include { X as Y } from ...` where X is a process with a container.

    An aliased process runs as its own process (FREYJA_DEMIX_COMPARE) and `withName` selectors treat it as
    such, so it needs its own lock entry. Aliased SUB-workflows are skipped: the processes inside keep
    their own names."""
    aliases: dict[str, str] = {}
    scripts = [root / "main.nf", *sorted((root / "workflows").rglob("*.nf"))]
    scripts += sorted((root / "subworkflows").rglob("*.nf"))
    for path in scripts:
        if not path.is_file():
            continue
        for match in INCLUDE_RE.finditer(path.read_text()):
            for original, alias in ALIAS_RE.findall(match.group(1)):
                if original in known:
                    aliases[alias] = original
    return aliases


def conf_overrides(root: Path) -> dict[str, str]:
    """process name -> image from `withName: 'X' { container = '...' }` in conf/*.config."""
    overrides: dict[str, str] = {}
    conf = root / "conf"
    if not conf.is_dir():
        return overrides
    for path in sorted(conf.glob("*.config")):
        if path.name in SKIP_CONF or path.name.startswith(SKIP_CONF_PREFIXES):
            continue
        for match in OVERRIDE_RE.finditer(path.read_text()):
            overrides[match.group(1)] = match.group(2)
    return overrides


def scan_pins(root: Path) -> dict[str, Ref]:
    """process name -> exact-pinned Ref (config overrides win over module defaults)."""
    images = {**module_images(root), **conf_overrides(root)}
    if not images:
        raise LockError(f"no process containers found under {root}")
    for alias, original in process_aliases(root, set(images)).items():
        images.setdefault(alias, images[original])
    pins: dict[str, Ref] = {}
    for process, image in sorted(images.items()):
        try:
            pins[process] = parse_ref(image)
        except LockError as err:
            raise LockError(f"process {process}: {err}") from err
    return pins


# ---------------------------------------------------------------------------------------------
# Registry (injectable)
# ---------------------------------------------------------------------------------------------


class Registry:
    """Read-only OCI registry client: anonymous HEAD of a manifest, with Bearer-token retry."""

    user_agent = "wasteline-lock-containers"

    def __init__(self):
        self._tokens: dict[tuple[str, str], str] = {}

    @staticmethod
    def host_of(registry: str) -> str:
        return "registry-1.docker.io" if registry in ("docker.io", "index.docker.io") else registry

    def resolve(self, ref: Ref) -> tuple[str, str]:
        """(digest, media type) of the manifest the tag points to."""
        host = self.host_of(ref.registry)
        repository = ref.repository
        if host == "registry-1.docker.io" and "/" not in repository:
            repository = f"library/{repository}"
        url = f"https://{host}/v2/{repository}/manifests/{ref.tag}"
        status, headers, body = self._request(url, host, repository, "HEAD")
        digest = headers.get("docker-content-digest")
        media = headers.get("content-type", "").split(";")[0].strip()
        if status == 200 and not digest:  # registry did not send the header: hash the body
            status, headers, body = self._request(url, host, repository, "GET")
            digest = "sha256:" + hashlib.sha256(body).hexdigest()
            media = headers.get("content-type", media).split(";")[0].strip()
        if status != 200 or not digest:
            raise LockError(f"cannot resolve {ref.pin}: HTTP {status}")
        if not DIGEST_RE.match(digest):
            raise LockError(f"cannot resolve {ref.pin}: unexpected digest {digest!r}")
        return digest, media

    def _request(self, url, host, repository, method):
        token = self._tokens.get((host, repository))
        for attempt in range(2):
            headers = {"Accept": ACCEPT, "User-Agent": self.user_agent}
            if token:
                headers["Authorization"] = f"Bearer {token}"
            request = urllib.request.Request(url, method=method, headers=headers)
            try:
                with urllib.request.urlopen(request, timeout=60) as response:
                    return (
                        response.status,
                        {k.lower(): v for k, v in response.headers.items()},
                        response.read() if method == "GET" else b"",
                    )
            except urllib.error.HTTPError as err:
                challenge = err.headers.get("WWW-Authenticate", "")
                if err.code == 401 and attempt == 0 and challenge.lower().startswith("bearer"):
                    token = self._fetch_token(challenge, repository)
                    self._tokens[(host, repository)] = token
                    continue
                return err.code, {k.lower(): v for k, v in err.headers.items()}, b""
            except (urllib.error.URLError, TimeoutError) as err:
                raise LockError(f"registry request failed: {url}: {err}") from err
        raise LockError(f"authentication failed for {url}")

    def _fetch_token(self, challenge: str, repository: str) -> str:
        params = dict(re.findall(r'(\w+)="([^"]*)"', challenge))
        realm = params.get("realm")
        if not realm:
            raise LockError(f"unsupported auth challenge: {challenge}")
        query = {"scope": params.get("scope", f"repository:{repository}:pull")}
        if "service" in params:
            query["service"] = params["service"]
        url = f"{realm}?{urllib.parse.urlencode(query)}"
        request = urllib.request.Request(url, headers={"User-Agent": self.user_agent})
        try:
            with urllib.request.urlopen(request, timeout=60) as response:
                data = json.load(response)
        except (urllib.error.URLError, TimeoutError, json.JSONDecodeError) as err:
            raise LockError(f"token request failed: {url}: {err}") from err
        token = data.get("token") or data.get("access_token")
        if not token:
            raise LockError(f"no token in response from {realm}")
        return token


# ---------------------------------------------------------------------------------------------
# Lock file
# ---------------------------------------------------------------------------------------------

HEADER = """\
// Container digest lock for wasteline (F13).
// GENERATED by bin/lock_containers.py - do not edit by hand.
// Regenerate with `bin/lock_containers.py --update` after any container pin changes.
// Used by `-profile locked`. Each pin is resolved to an immutable digest; multi-arch images
// record the manifest-list digest. The trailing comment is the pin the digest was resolved from.
"""

LINE_RE = re.compile(
    r"^\s*withName:\s*'(?P<process>[A-Za-z0-9_]+)'\s*\{\s*container\s*=\s*"
    r"'(?P<name>[^'@\s]+)@(?P<digest>sha256:[0-9a-f]{64})'\s*\}\s*//\s*pin:\s*(?P<pin>\S+)"
    r"(?:\s+\((?P<kind>[a-z]+)\))?\s*$"
)


@dataclass(frozen=True)
class Locked:
    pin: str
    digest: str
    kind: str


def kind_of(media: str) -> str:
    return "index" if media in INDEX_TYPES else "manifest"


def render_lock(entries: dict[str, tuple[Ref, str, str]]) -> str:
    """entries: process -> (ref, digest, kind). Deterministic (sorted), no timestamps."""
    width = max((len(p) for p in entries), default=0)
    lines = [HEADER.rstrip("\n"), "process {"]
    for process in sorted(entries):
        ref, digest, kind = entries[process]
        head = f"    withName: '{process}'".ljust(len("    withName: ''") + width)
        lines.append(
            f"{head} {{ container = '{ref.digest_name}@{digest}' }} // pin: {ref.pin} ({kind})"
        )
    lines.append("}")
    return "\n".join(lines) + "\n"


def parse_lock(text: str) -> dict[str, Locked]:
    locked: dict[str, Locked] = {}
    for number, line in enumerate(text.splitlines(), start=1):
        stripped = line.strip()
        if not stripped or stripped.startswith("//") or stripped in ("process {", "}"):
            continue
        match = LINE_RE.match(line)
        if not match:
            raise LockError(f"{LOCK_FILE}:{number}: cannot parse lock line: {line!r}")
        locked[match["process"]] = Locked(
            match["pin"], match["digest"], match["kind"] or "manifest"
        )
    return locked


def stale_reasons(pins: dict[str, Ref], locked: dict[str, Locked]) -> list[str]:
    reasons = []
    for process, ref in sorted(pins.items()):
        entry = locked.get(process)
        if entry is None:
            reasons.append(f"{process}: not in the lock (pin {ref.pin})")
        elif entry.pin != ref.pin:
            reasons.append(f"{process}: lock is for {entry.pin}, pin is now {ref.pin}")
    for process in sorted(set(locked) - set(pins)):
        reasons.append(f"{process}: in the lock but no longer a process container")
    return reasons


def build_entries(
    pins: dict[str, Ref], previous: dict[str, Locked], registry, refresh: bool = False
) -> dict[str, tuple[Ref, str, str]]:
    cache: dict[str, tuple[str, str]] = {}
    entries = {}
    for process, ref in sorted(pins.items()):
        old = previous.get(process)
        if old and old.pin == ref.pin and not refresh:
            entries[process] = (ref, old.digest, old.kind)
            continue
        if ref.pin not in cache:
            digest, media = registry.resolve(ref)
            cache[ref.pin] = (digest, kind_of(media))
        entries[process] = (ref, *cache[ref.pin])
    return entries


def verify_digests(pins: dict[str, Ref], locked: dict[str, Locked], registry) -> list[str]:
    """Re-resolve each locked pin and report digests that moved (a retagged image)."""
    problems = []
    seen: dict[str, str] = {}
    for process, ref in sorted(pins.items()):
        entry = locked.get(process)
        if entry is None or entry.pin != ref.pin:
            continue
        if ref.pin not in seen:
            seen[ref.pin] = registry.resolve(ref)[0]
        if seen[ref.pin] != entry.digest:
            problems.append(
                f"{process}: {ref.pin} now resolves to {seen[ref.pin]}, lock has {entry.digest}"
            )
    return problems


# ---------------------------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------------------------


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--check", action="store_true", help="fail (1) if the lock is stale")
    mode.add_argument("--update", action="store_true", help="(re)write the lock file")
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--lock-file", type=Path, help=f"default: <root>/{LOCK_FILE}")
    parser.add_argument("--refresh", action="store_true", help="--update: re-resolve every digest")
    parser.add_argument(
        "--verify-digests", action="store_true", help="--check: also re-resolve every digest"
    )
    parser.add_argument("--json", action="store_true", help="machine-readable output")
    return parser


def main(argv=None, registry=None) -> int:
    args = build_parser().parse_args(argv)
    registry = registry or Registry()
    root = args.root.resolve()
    lock_path = args.lock_file or root / LOCK_FILE
    try:
        pins = scan_pins(root)
        previous = parse_lock(lock_path.read_text()) if lock_path.exists() else {}

        if args.update:
            entries = build_entries(pins, previous, registry, refresh=args.refresh)
            text = render_lock(entries)
            changed = not lock_path.exists() or lock_path.read_text() != text
            if changed:
                lock_path.parent.mkdir(parents=True, exist_ok=True)
                lock_path.write_text(text)
            print(
                json.dumps({"updated": changed, "processes": len(entries)})
                if args.json
                else f"{'Wrote' if changed else 'Unchanged'} {lock_path} ({len(entries)} processes)"
            )
            return EXIT_OK

        # --check
        if not lock_path.exists():
            reasons = [f"{lock_path} does not exist; run bin/lock_containers.py --update"]
        else:
            reasons = stale_reasons(pins, previous)
            if args.verify_digests and not reasons:
                reasons = verify_digests(pins, previous, registry)
        if args.json:
            print(json.dumps({"stale": bool(reasons), "reasons": reasons, "processes": len(pins)}))
        elif reasons:
            print("Container lock is STALE (run: bin/lock_containers.py --update):")
            print("\n".join(f"  - {r}" for r in reasons))
        else:
            print(f"Container lock is up to date ({len(pins)} processes).")
        return EXIT_STALE if reasons else EXIT_OK
    except LockError as err:
        print(f"error: {err}", file=sys.stderr)
        return EXIT_ERROR


if __name__ == "__main__":
    sys.exit(main())
