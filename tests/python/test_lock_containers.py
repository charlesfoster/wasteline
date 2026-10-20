from pathlib import Path

import lock_containers as lc

REPO = Path(__file__).resolve().parents[2]

D1 = "sha256:" + "1" * 64
D2 = "sha256:" + "2" * 64
D3 = "sha256:" + "3" * 64

LOCAL_MAIN = """\
process FREYJA_DEMIX {
    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.x
        ? 'https://depot.galaxyproject.org/singularity/freyja:2.0.5--pyhdfd78af_0'
        : 'quay.io/biocontainers/freyja:2.0.5--pyhdfd78af_0'}"
    script: "true"
}
"""

NFCORE_MAIN = """\
process FASTP {
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.x ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/d0/d013/data' :
        'community.wave.seqera.io/library/fastp:1.3.6--4df8d6c11b471bde' }"
    script: "true"
}
"""


class FakeRegistry:
    def __init__(self, digests=None, kinds=None):
        self.digests = digests or {}
        self.kinds = kinds or {}
        self.calls = []

    def resolve(self, ref):
        self.calls.append(ref.pin)
        if ref.pin not in self.digests:
            raise lc.LockError(f"cannot resolve {ref.pin}: HTTP 404")
        media = (
            "application/vnd.oci.image.index.v1+json"
            if self.kinds.get(ref.pin) == "index"
            else "application/vnd.docker.distribution.manifest.v2+json"
        )
        return self.digests[ref.pin], media


def make_repo(tmp_path):
    for rel, text in {
        "modules/local/freyja_demix/main.nf": LOCAL_MAIN,
        "modules/nf-core/fastp/main.nf": NFCORE_MAIN,
    }.items():
        path = tmp_path / rel
        path.parent.mkdir(parents=True)
        path.write_text(text)
    (tmp_path / "conf").mkdir()
    return tmp_path


FREYJA = "quay.io/biocontainers/freyja:2.0.5--pyhdfd78af_0"
FASTP = "community.wave.seqera.io/library/fastp:1.3.6--4df8d6c11b471bde"
DIGESTS = {FREYJA: D1, FASTP: D2}


def test_real_repo_containers_are_exact_pinned_and_aliases_locked():
    """CLAUDE.md: every container is pinned to an exact tag, never latest."""
    pins = lc.scan_pins(REPO)
    assert pins
    assert all(ref.tag != "latest" for ref in pins.values())
    assert "FREYJA_DEMIX" in pins
    aliases = lc.process_aliases(REPO, set(pins))
    assert aliases["FREYJA_DEMIX_COMPARE"] == "FREYJA_DEMIX"
    assert set(aliases) <= set(pins)


def run(root, *argv, registry=None):
    return lc.main([*argv, "--root", str(root)], registry=registry or FakeRegistry(DIGESTS))


def test_cli_update_then_check_then_stale_after_a_bump(tmp_path, capsys):
    root = make_repo(tmp_path)
    assert run(root, "--check") == 1  # no lock yet
    assert run(root, "--update") == 0
    lock = (root / "conf" / "containers.lock.config").read_text()
    assert D1 in lock and D2 in lock
    assert run(root, "--check") == 0
    assert run(root, "--update") == 0  # idempotent
    assert (root / "conf" / "containers.lock.config").read_text() == lock

    # a Freyja bump: the lock is now stale until --update, which resolves only the new pin
    path = root / "modules/local/freyja_demix/main.nf"
    path.write_text(path.read_text().replace("2.0.5--pyhdfd78af_0", "2.1.0--pyhdfd78af_0"))
    assert run(root, "--check") == 1
    assert "STALE" in capsys.readouterr().out
    bumped = "quay.io/biocontainers/freyja:2.1.0--pyhdfd78af_0"
    registry = FakeRegistry({**DIGESTS, bumped: D3})
    assert run(root, "--update", registry=registry) == 0
    assert registry.calls == [bumped]
    assert run(root, "--check") == 0
    assert f"freyja@{D3}" in (root / "conf" / "containers.lock.config").read_text()
