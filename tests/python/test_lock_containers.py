from pathlib import Path

import lock_containers as lc
import pytest

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


# --- references -----------------------------------------------------------------------------


@pytest.mark.parametrize(
    ("image", "registry", "repository", "tag"),
    [
        (FREYJA, "quay.io", "biocontainers/freyja", "2.0.5--pyhdfd78af_0"),
        (FASTP, "community.wave.seqera.io", "library/fastp", "1.3.6--4df8d6c11b471bde"),
        (
            "biocontainers/ivar:1.4.4--h077b44d_0",
            "quay.io",
            "biocontainers/ivar",
            "1.4.4--h077b44d_0",
        ),
        ("localhost:5000/team/tool:1.0", "localhost:5000", "team/tool", "1.0"),
        ("ghcr.io/org/tool:v1", "ghcr.io", "org/tool", "v1"),
    ],
)
def test_parse_ref(image, registry, repository, tag):
    ref = lc.parse_ref(image)
    assert (ref.registry, ref.repository, ref.tag) == (registry, repository, tag)


@pytest.mark.parametrize(
    "image",
    [
        "quay.io/biocontainers/freyja",
        "quay.io/biocontainers/freyja:latest",
        "localhost:5000/team/tool",
        "quay.io/x/y:1@sha256:" + "a" * 64,
    ],
)
def test_floating_or_missing_tags_are_refused(image):
    with pytest.raises(lc.LockError):
        lc.parse_ref(image)


def test_docker_image_of_takes_the_non_url_branch():
    assert lc.docker_image_of(LOCAL_MAIN.split('"')[3]) == FREYJA
    assert lc.docker_image_of("'a.io/x:1'") == "a.io/x:1"
    assert lc.docker_image_of("a.io/x:1") == "a.io/x:1"
    assert lc.docker_image_of("${ cond ? 'https://depot/x:1' : 'oras://y:2' }") is None


# --- scanning -------------------------------------------------------------------------------


def test_scan_finds_module_images(tmp_path):
    pins = lc.scan_pins(make_repo(tmp_path))
    assert set(pins) == {"FREYJA_DEMIX", "FASTP"}
    assert pins["FREYJA_DEMIX"].pin == FREYJA
    assert pins["FASTP"].pin == FASTP


def test_aliased_process_includes_get_their_own_pin(tmp_path):
    root = make_repo(tmp_path)
    (root / "workflows").mkdir()
    (root / "workflows" / "main.nf").write_text(
        "include { FREYJA_DEMIX           } from '../modules/local/freyja_demix/main'\n"
        "include { FREYJA_DEMIX as FREYJA_DEMIX_COMPARE } from '../modules/local/freyja_demix/main'\n"
        "include { FREYJA_DATA as FREYJA_DATA_COMPARE   } from '../subworkflows/local/freyja_data/main'\n"
    )
    pins = lc.scan_pins(root)
    assert set(pins) == {
        "FREYJA_DEMIX",
        "FREYJA_DEMIX_COMPARE",
        "FASTP",
    }  # sub-workflow alias: no entry
    assert pins["FREYJA_DEMIX_COMPARE"] == pins["FREYJA_DEMIX"]


def test_real_repo_has_a_lock_entry_for_every_aliased_process():
    root = Path(lc.__file__).resolve().parents[1]
    pins = lc.scan_pins(root)
    aliases = lc.process_aliases(root, set(pins))
    assert "FREYJA_DEMIX_COMPARE" in aliases and aliases["FREYJA_DEMIX_COMPARE"] == "FREYJA_DEMIX"
    assert set(aliases) <= set(pins)


def test_conf_override_wins_and_lock_and_wave_configs_are_ignored(tmp_path):
    root = make_repo(tmp_path)
    (root / "conf" / "modules.config").write_text(
        "process { withName: 'FASTP' { container = 'quay.io/biocontainers/fastp:1.0.0--h_0' } }\n"
    )
    (root / "conf" / "containers.lock.config").write_text(
        "process { withName: 'GHOST' { container = 'x.io/y:1' } }\n"
    )
    (root / "conf" / "containers_docker_amd64.config").write_text(
        "process { withName: 'FASTP' { container = 'wave.io/fastp:9' } }\n"
    )
    pins = lc.scan_pins(root)
    assert pins["FASTP"].pin == "quay.io/biocontainers/fastp:1.0.0--h_0"
    assert "GHOST" not in pins


def test_unpinned_module_image_is_an_error(tmp_path):
    root = make_repo(tmp_path)
    path = root / "modules/local/freyja_demix/main.nf"
    path.write_text(path.read_text().replace("2.0.5--pyhdfd78af_0'}", "latest'}"))
    with pytest.raises(lc.LockError, match="FREYJA_DEMIX"):
        lc.scan_pins(root)


def test_module_without_a_docker_image_is_an_error(tmp_path):
    root = make_repo(tmp_path)
    (root / "modules/local/freyja_demix/main.nf").write_text(
        "process X {\n    container \"${ 'https://depot/x:1' }\"\n}\n"
    )
    with pytest.raises(lc.LockError, match="cannot find a docker image"):
        lc.scan_pins(root)


def test_real_repo_containers_are_exact_pinned():
    """CLAUDE.md: every container is pinned to an exact tag, never latest."""
    pins = lc.scan_pins(REPO)
    assert pins
    assert all(ref.tag != "latest" for ref in pins.values())
    assert "FREYJA_DEMIX" in pins


# --- lock file ------------------------------------------------------------------------------


def entries_for(pins, digests=DIGESTS):
    return {p: (r, digests[r.pin], "manifest") for p, r in pins.items()}


def test_render_parse_roundtrip_and_determinism(tmp_path):
    pins = lc.scan_pins(make_repo(tmp_path))
    text = lc.render_lock(entries_for(pins))
    assert text == lc.render_lock(dict(reversed(list(entries_for(pins).items()))))
    assert f"container = 'quay.io/biocontainers/freyja@{D1}'" in text
    assert f"// pin: {FREYJA} (manifest)" in text
    parsed = lc.parse_lock(text)
    assert parsed["FREYJA_DEMIX"] == lc.Locked(FREYJA, D1, "manifest")
    assert lc.stale_reasons(pins, parsed) == []


def test_garbage_in_lock_is_an_error():
    with pytest.raises(lc.LockError):
        lc.parse_lock("process {\n  nonsense\n}\n")


def test_stale_reasons(tmp_path):
    pins = lc.scan_pins(make_repo(tmp_path))
    locked = lc.parse_lock(lc.render_lock(entries_for(pins)))
    locked["FREYJA_DEMIX"] = lc.Locked(
        "quay.io/biocontainers/freyja:2.0.4--pyhdfd78af_0", D1, "manifest"
    )
    del locked["FASTP"]
    locked["OLD"] = lc.Locked("x.io/y:1", D3, "manifest")
    reasons = "\n".join(lc.stale_reasons(pins, locked))
    assert "FREYJA_DEMIX: lock is for" in reasons
    assert "FASTP: not in the lock" in reasons
    assert "OLD: in the lock but no longer" in reasons


def test_index_digest_is_recorded_as_index(tmp_path):
    pins = lc.scan_pins(make_repo(tmp_path))
    registry = FakeRegistry(DIGESTS, kinds={FREYJA: "index"})
    entries = lc.build_entries(pins, {}, registry)
    assert entries["FREYJA_DEMIX"][2] == "index"
    assert "(index)" in lc.render_lock(entries)


def test_update_reuses_unchanged_digests_and_resolves_each_image_once(tmp_path):
    pins = lc.scan_pins(make_repo(tmp_path))
    registry = FakeRegistry(DIGESTS)
    first = lc.build_entries(pins, {}, registry)
    assert sorted(registry.calls) == sorted([FREYJA, FASTP])
    previous = lc.parse_lock(lc.render_lock(first))
    again = FakeRegistry(DIGESTS)
    lc.build_entries(pins, previous, again)
    assert again.calls == []
    lc.build_entries(pins, previous, again, refresh=True)
    assert sorted(again.calls) == sorted([FREYJA, FASTP])


def test_verify_digests_detects_a_retagged_image(tmp_path):
    pins = lc.scan_pins(make_repo(tmp_path))
    previous = lc.parse_lock(lc.render_lock(entries_for(pins)))
    assert lc.verify_digests(pins, previous, FakeRegistry(DIGESTS)) == []
    moved = lc.verify_digests(pins, previous, FakeRegistry({**DIGESTS, FASTP: D3}))
    assert len(moved) == 1 and "FASTP" in moved[0] and D3 in moved[0]


# --- CLI ------------------------------------------------------------------------------------


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


def test_cli_update_failure_writes_nothing_and_exits_2(tmp_path, capsys):
    root = make_repo(tmp_path)
    assert run(root, "--update", registry=FakeRegistry({FREYJA: D1})) == 2
    assert not (root / "conf" / "containers.lock.config").exists()
    assert "cannot resolve" in capsys.readouterr().err


def test_cli_check_verify_digests(tmp_path):
    root = make_repo(tmp_path)
    assert run(root, "--update") == 0
    assert run(root, "--check", "--verify-digests") == 0
    assert (
        run(root, "--check", "--verify-digests", registry=FakeRegistry({**DIGESTS, FREYJA: D3}))
        == 1
    )


def test_cli_modes_are_mutually_exclusive(tmp_path):
    with pytest.raises(SystemExit):
        run(tmp_path, "--check", "--update")


# --- Registry client (no network: _request is replaced) --------------------------------------


class ScriptedRegistry(lc.Registry):
    def __init__(self, responses):
        super().__init__()
        self.responses = list(responses)
        self.requests = []

    def _request(self, url, host, repository, method):
        self.requests.append((url, method))
        return self.responses.pop(0)


def test_registry_reads_docker_content_digest_from_head():
    registry = ScriptedRegistry(
        [
            (
                200,
                {
                    "docker-content-digest": D1,
                    "content-type": "application/vnd.oci.image.index.v1+json",
                },
                b"",
            )
        ]
    )
    digest, media = registry.resolve(lc.parse_ref(FREYJA))
    assert digest == D1
    assert lc.kind_of(media) == "index"
    assert registry.requests == [
        ("https://quay.io/v2/biocontainers/freyja/manifests/2.0.5--pyhdfd78af_0", "HEAD")
    ]


def test_registry_falls_back_to_hashing_the_body():
    import hashlib

    body = b'{"schemaVersion": 2}'
    registry = ScriptedRegistry([(200, {}, b""), (200, {"content-type": "application/json"}, body)])
    digest, _ = registry.resolve(lc.parse_ref(FREYJA))
    assert digest == "sha256:" + hashlib.sha256(body).hexdigest()


def test_registry_404_is_an_error():
    with pytest.raises(lc.LockError, match="HTTP 404"):
        ScriptedRegistry([(404, {}, b"")]).resolve(lc.parse_ref(FREYJA))


def test_registry_rejects_a_malformed_digest():
    with pytest.raises(lc.LockError, match="unexpected digest"):
        ScriptedRegistry([(200, {"docker-content-digest": "md5:zz"}, b"")]).resolve(
            lc.parse_ref(FREYJA)
        )


def test_docker_hub_library_images_get_the_library_prefix():
    registry = ScriptedRegistry([(200, {"docker-content-digest": D1}, b"")])
    registry.resolve(lc.Ref("docker.io", "ubuntu", "24.04"))
    assert (
        registry.requests[0][0] == "https://registry-1.docker.io/v2/library/ubuntu/manifests/24.04"
    )
