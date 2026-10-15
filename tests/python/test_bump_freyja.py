import json
from pathlib import Path

import bump_freyja as bf
import pytest

FIXTURES = Path(__file__).resolve().parents[1] / "fixtures" / "ci"
REPO = Path(__file__).resolve().parents[2]

QUAY_REAL = json.loads((FIXTURES / "quay_tags_freyja.json").read_text())
QUAY_ODD = json.loads((FIXTURES / "quay_tags_odd.json").read_text())
ANACONDA = json.loads((FIXTURES / "anaconda_freyja.json").read_text())


class FakeNetwork:
    """Serves the recorded fixtures; every call is recorded for assertions."""

    def __init__(self, quay=None, anaconda=None, depot=200, fail=None):
        self.quay = quay if quay is not None else QUAY_REAL
        self.anaconda = anaconda if anaconda is not None else ANACONDA
        self.depot = depot
        self.fail = fail
        self.calls = []

    def get_json(self, url):
        self.calls.append(url)
        if self.fail:
            raise bf.BumpError(self.fail)
        if url.startswith("https://quay.io/api/v1/repository/biocontainers/freyja/tag/"):
            return self.quay
        if url == bf.ANACONDA_URL:
            return self.anaconda
        raise AssertionError(f"unexpected URL {url}")

    def head_status(self, url):
        self.calls.append(url)
        assert url.startswith("https://depot.galaxyproject.org/singularity/freyja:")
        return self.depot


def with_release(base, version, build="pyhdfd78af_0"):
    """Copy of a quay fixture with an extra release tag."""
    data = json.loads(json.dumps(base))
    extra = dict(data["tags"][0])
    extra["name"] = f"{version}--{build}"
    data["tags"].insert(0, extra)
    return data


def anaconda_with(*versions):
    data = json.loads(json.dumps(ANACONDA))
    data["versions"] = data["versions"] + list(versions)
    return data


def make_repo(tmp_path, tag="2.0.4--pyhdfd78af_0"):
    """A miniature repo with Freyja pinned in every spelling, in several places."""
    version = tag.split("--")[0]
    for module in ("freyja_demix", "freyja_update"):
        mod = tmp_path / "modules" / "local" / module
        mod.mkdir(parents=True)
        (mod / "main.nf").write_text(
            "process X {\n"
            "    container \"${workflow.containerEngine == 'singularity'\n"
            f"        ? 'https://depot.galaxyproject.org/singularity/freyja:{tag}'\n"
            f"        : 'quay.io/biocontainers/freyja:{tag}'}}\"\n"
            "}\n"
        )
        (mod / "environment.yml").write_text(
            f"dependencies:\n  - bioconda::freyja={version}\n  - bioconda::samtools=1.24\n"
        )
    (tmp_path / "docs").mkdir()
    (tmp_path / "docs" / "usage.md").write_text(f"Runs `quay.io/biocontainers/freyja:{tag}`.\n")
    (tmp_path / "conf").mkdir()
    (tmp_path / "conf" / "x.config").write_text(f"// docker://quay.io/biocontainers/freyja:{tag}\n")
    (tmp_path / "CHANGELOG.md").write_text(
        "# Changelog\n\n## v0.1.0dev - [unreleased]\n\n### `Added`\n\n### `Dependencies`\n\n"
        "### `Deprecated`\n"
    )
    return tmp_path


# --- tag parsing ----------------------------------------------------------------------------


@pytest.mark.parametrize(
    "raw",
    [
        "latest",
        "dev",
        "2.0.5",
        "2.0.11rc1--pyhdfd78af_0",
        "2.0.11.dev1--pyhdfd78af_0",
        "3.0.0b2--pyhdfd78af_0",
        "2.0.10--pyhdfd78af_1-arm64",
        "2.0.5--pyhdfd78af",
        "--pyhdfd78af_0",
        "",
    ],
)
def test_non_release_tags_are_rejected(raw):
    assert bf.parse_tag(raw) is None


@pytest.mark.parametrize(
    ("raw", "version", "build", "number"),
    [
        ("2.0.5--pyhdfd78af_0", "2.0.5", "pyhdfd78af", 0),
        ("2.0.0--py311h9ee0642_0", "2.0.0", "py311h9ee0642", 0),
        ("1.3--pyhdfd78af_1", "1.3", "pyhdfd78af", 1),
    ],
)
def test_release_tags_parse(raw, version, build, number):
    tag = bf.parse_tag(raw)
    assert (tag.version, tag.build, tag.number) == (version, build, number)


def test_pick_latest_on_recorded_quay_listing():
    names = [t["name"] for t in QUAY_REAL["tags"]]
    assert bf.pick_latest(names).raw == "2.0.5--pyhdfd78af_0"


def test_pick_latest_ignores_odd_tags_and_compares_numerically():
    names = [t["name"] for t in QUAY_ODD["tags"]]
    # 2.0.10 > 2.0.9 (not lexicographic); highest build number wins; pyhdfd78af wins a tie;
    # rc/dev/beta/arch-suffixed/latest tags never win even though 3.0.0b2 looks "higher".
    assert bf.pick_latest(names).raw == "2.0.10--pyhdfd78af_1"


def test_pick_latest_none_when_no_release_tags():
    assert bf.pick_latest(["latest", "dev"]) is None


def test_major_bump():
    old, new = bf.parse_tag("1.5.3--pyhdfd78af_2"), bf.parse_tag("2.0.0--py311h9ee0642_0")
    assert bf.is_major_bump(old, new)
    assert not bf.is_major_bump(bf.parse_tag("2.0.4--pyhdfd78af_0"), bf.parse_tag("2.1.0--x_0"))
    assert not bf.is_major_bump(
        bf.parse_tag("2.0.5--pyhdfd78af_0"), bf.parse_tag("2.0.5--pyhdfd78af_1")
    )


# --- network layer --------------------------------------------------------------------------


def test_quay_pagination_is_followed():
    page1 = {"tags": [{"name": "2.0.5--pyhdfd78af_0"}], "has_additional": True, "page": 1}
    page2 = {"tags": [{"name": "2.1.0--pyhdfd78af_0"}], "has_additional": False, "page": 2}

    class Paged(FakeNetwork):
        def get_json(self, url):
            self.calls.append(url)
            return page2 if url.endswith("&page=2") else page1

    net = Paged()
    assert bf.fetch_quay_tag_names(net) == ["2.0.5--pyhdfd78af_0", "2.1.0--pyhdfd78af_0"]
    assert len(net.calls) == 2


def test_malformed_quay_response_is_an_error():
    with pytest.raises(bf.BumpError):
        bf.fetch_quay_tag_names(FakeNetwork(quay={"oops": []}))


# --- check ----------------------------------------------------------------------------------


def test_check_up_to_date(tmp_path):
    root = make_repo(tmp_path, "2.0.5--pyhdfd78af_0")
    result = bf.check(root, FakeNetwork())
    assert result["status"] == "up_to_date"
    assert not result["update_available"]
    assert not result["major_bump"]
    assert result["checks"] == {"quay": "ok", "bioconda": "ok", "singularity_depot": "ok"}


def test_check_update_available_same_major(tmp_path):
    root = make_repo(tmp_path, "2.0.4--pyhdfd78af_0")
    result = bf.check(root, FakeNetwork())
    assert result["status"] == "update_available"
    assert result["update_available"]
    assert result["major_bump"] is False
    assert (result["pinned"], result["latest"]) == ("2.0.4--pyhdfd78af_0", "2.0.5--pyhdfd78af_0")


def test_check_major_bump_flag(tmp_path):
    root = make_repo(tmp_path, "2.0.5--pyhdfd78af_0")
    net = FakeNetwork(quay=with_release(QUAY_REAL, "3.0.0"), anaconda=anaconda_with("3.0.0"))
    result = bf.check(root, net)
    assert result["update_available"]
    assert result["major_bump"] is True
    assert result["latest_version"] == "3.0.0"


def test_check_new_build_of_same_version_is_an_update_but_not_major(tmp_path):
    root = make_repo(tmp_path, "2.0.5--pyhdfd78af_0")
    net = FakeNetwork(quay=with_release(QUAY_REAL, "2.0.5", "pyhdfd78af_1"))
    result = bf.check(root, net)
    assert result["update_available"]
    assert not result["major_bump"]
    assert result["latest"] == "2.0.5--pyhdfd78af_1"


def test_check_pinned_ahead_of_registry_is_not_an_update(tmp_path):
    root = make_repo(tmp_path, "9.9.9--pyhdfd78af_0")
    result = bf.check(root, FakeNetwork())
    assert result["status"] == "up_to_date"
    assert not result["major_bump"]


def test_check_blocked_when_bioconda_lacks_the_version(tmp_path):
    root = make_repo(tmp_path, "2.0.5--pyhdfd78af_0")
    net = FakeNetwork(quay=with_release(QUAY_REAL, "2.1.0"))  # anaconda has no 2.1.0
    result = bf.check(root, net)
    assert result["status"] == "blocked"
    assert not result["update_available"]
    assert result["checks"]["bioconda"] == "missing"
    assert any("bioconda" in r for r in result["blocked_reasons"])


def test_check_blocked_when_depot_lacks_the_image(tmp_path):
    root = make_repo(tmp_path, "2.0.5--pyhdfd78af_0")
    net = FakeNetwork(
        quay=with_release(QUAY_REAL, "2.1.0"), anaconda=anaconda_with("2.1.0"), depot=404
    )
    result = bf.check(root, net)
    assert result["status"] == "blocked"
    assert result["checks"]["singularity_depot"] == "missing"


@pytest.mark.parametrize("depot", [None, 503])
def test_depot_unreachable_is_skipped_not_blocking(tmp_path, depot):
    root = make_repo(tmp_path, "2.0.4--pyhdfd78af_0")
    result = bf.check(root, FakeNetwork(depot=depot))
    assert result["checks"]["singularity_depot"] == "skipped"
    assert result["update_available"]


def test_no_depot_check_makes_no_head_request(tmp_path):
    root = make_repo(tmp_path, "2.0.4--pyhdfd78af_0")
    net = FakeNetwork()
    bf.check(root, net, check_depot=False)
    assert not any("depot" in c for c in net.calls)


def test_check_refuses_inconsistent_pins(tmp_path):
    root = make_repo(tmp_path, "2.0.4--pyhdfd78af_0")
    (root / "docs" / "usage.md").write_text("quay.io/biocontainers/freyja:2.0.3--pyhdfd78af_0\n")
    with pytest.raises(bf.BumpError, match="inconsistent"):
        bf.check(root, FakeNetwork())


# --- consistency ----------------------------------------------------------------------------


def test_consistent_repo_passes(tmp_path):
    assert bf.consistency_problems(bf.scan_pins(make_repo(tmp_path))) == []


def test_each_kind_of_disagreement_is_reported(tmp_path):
    root = make_repo(tmp_path, "2.0.4--pyhdfd78af_0")
    env = root / "modules" / "local" / "freyja_demix" / "environment.yml"
    env.write_text(env.read_text().replace("freyja=2.0.4", "freyja=2.0.3"))
    main_nf = root / "modules" / "local" / "freyja_update" / "main.nf"
    main_nf.write_text(
        main_nf.read_text().replace("singularity/freyja:2.0.4", "singularity/freyja:2.0.3")
    )
    (root / "docs" / "usage.md").write_text("quay.io/biocontainers/freyja:2.0.9--pyhdfd78af_0\n")
    problems = "\n".join(bf.consistency_problems(bf.scan_pins(root)))
    assert "quay pins disagree" in problems
    assert "depot pins disagree" in problems
    assert "conda pins disagree" in problems


def test_conda_version_must_match_container_version(tmp_path):
    root = make_repo(tmp_path, "2.0.4--pyhdfd78af_0")
    for env in root.glob("modules/local/*/environment.yml"):
        env.write_text(env.read_text().replace("freyja=2.0.4", "freyja=2.0.5"))
    problems = bf.consistency_problems(bf.scan_pins(root))
    assert any("differs from the container version" in p for p in problems)


def test_module_missing_a_pin_kind_is_reported(tmp_path):
    root = make_repo(tmp_path)
    (root / "modules" / "local" / "freyja_demix" / "environment.yml").write_text(
        "dependencies: []\n"
    )
    problems = bf.consistency_problems(bf.scan_pins(root))
    assert any("freyja_demix lacks Freyja pin(s): conda" in p for p in problems)


def test_empty_repo_is_inconsistent(tmp_path):
    assert bf.consistency_problems(bf.scan_pins(tmp_path))


def test_ignore_marker_and_excluded_files_are_not_pins(tmp_path):
    root = make_repo(tmp_path)
    (root / "docs" / "note.md").write_text(
        "old: quay.io/biocontainers/freyja:1.0.0--pyhdfd78af_0 bump-freyja: ignore\n"
    )
    (root / "tests" / "fixtures" / "ci").mkdir(parents=True)
    (root / "tests" / "fixtures" / "ci" / "x.txt").write_text("bioconda::freyja=0.1\n")
    (root / "docs" / "PLAN.md").write_text("quay.io/biocontainers/freyja:0.0.1--pyhdfd78af_0\n")
    assert bf.consistency_problems(bf.scan_pins(root)) == []


def test_trailing_prose_punctuation_is_not_part_of_the_tag(tmp_path):
    root = make_repo(tmp_path, "2.0.5--pyhdfd78af_0")
    (root / "docs" / "usage.md").write_text(
        "Uses quay.io/biocontainers/freyja:2.0.5--pyhdfd78af_0.\n"
    )
    assert bf.consistency_problems(bf.scan_pins(root)) == []
    bf.apply_tag(root, bf.parse_tag("2.1.0--pyhdfd78af_0"))
    assert (root / "docs" / "usage.md").read_text() == (
        "Uses quay.io/biocontainers/freyja:2.1.0--pyhdfd78af_0.\n"
    )


def test_real_repo_pins_are_consistent():
    """Invariant (PLAN section 9): every Freyja pin in the repository agrees."""
    pins = bf.scan_pins(REPO)
    assert bf.consistency_problems(pins) == []
    kinds = {p.kind for p in pins}
    assert set(bf.REQUIRED_KINDS) <= kinds


# --- apply ----------------------------------------------------------------------------------


def snapshot(root):
    return {
        p.relative_to(root).as_posix(): p.read_bytes()
        for p in sorted(root.rglob("*"))
        if p.is_file()
    }


def test_apply_rewrites_every_pin_and_the_changelog(tmp_path):
    root = make_repo(tmp_path, "2.0.4--pyhdfd78af_0")
    outcome = bf.apply_tag(root, bf.parse_tag("2.1.0--pyhdfd78af_0"))
    assert outcome["applied"]
    assert outcome["changelog"]
    assert sorted(outcome["files"]) == sorted(
        [
            "conf/x.config",
            "docs/usage.md",
            "modules/local/freyja_demix/environment.yml",
            "modules/local/freyja_demix/main.nf",
            "modules/local/freyja_update/environment.yml",
            "modules/local/freyja_update/main.nf",
        ]
    )
    pins = bf.scan_pins(root)
    assert bf.consistency_problems(pins) == []
    assert {p.value for p in pins if p.kind == "quay"} == {"2.1.0--pyhdfd78af_0"}
    assert {p.value for p in pins if p.kind == "conda"} == {"2.1.0"}
    # unrelated pins untouched
    assert (
        "bioconda::samtools=1.24"
        in (root / "modules/local/freyja_demix/environment.yml").read_text()
    )
    changelog = (root / "CHANGELOG.md").read_text()
    assert "`2.0.4--pyhdfd78af_0` -> `2.1.0--pyhdfd78af_0`" in changelog


def test_apply_is_idempotent_and_deterministic(tmp_path):
    root = make_repo(tmp_path, "2.0.4--pyhdfd78af_0")
    bf.apply_tag(root, bf.parse_tag("2.1.0--pyhdfd78af_0"))
    first = snapshot(root)
    again = bf.apply_tag(root, bf.parse_tag("2.1.0--pyhdfd78af_0"))
    assert not again["applied"] and not again["changelog"]
    assert snapshot(root) == first


def test_apply_dry_run_writes_nothing(tmp_path):
    root = make_repo(tmp_path, "2.0.4--pyhdfd78af_0")
    before = snapshot(root)
    outcome = bf.apply_tag(root, bf.parse_tag("2.1.0--pyhdfd78af_0"), dry_run=True)
    assert outcome["files"] and outcome["dry_run"]
    assert snapshot(root) == before


def test_apply_repairs_inconsistent_pins(tmp_path):
    root = make_repo(tmp_path, "2.0.4--pyhdfd78af_0")
    (root / "docs" / "usage.md").write_text("quay.io/biocontainers/freyja:2.0.1--pyhdfd78af_0\n")
    bf.apply_tag(root, bf.parse_tag("2.0.5--pyhdfd78af_0"))
    assert bf.consistency_problems(bf.scan_pins(root)) == []


# --- changelog ------------------------------------------------------------------------------

OLD, NEW = bf.parse_tag("2.0.4--pyhdfd78af_0"), bf.parse_tag("2.0.5--pyhdfd78af_0")
BULLET = bf.changelog_line(OLD, NEW)


def test_changelog_entry_goes_into_empty_dependencies_section():
    text = "# C\n\n## v1 - unreleased\n\n### `Added`\n\n### `Dependencies`\n\n### `Fixed`\n"
    assert bf.add_changelog_entry(text, OLD, NEW) == (
        "# C\n\n## v1 - unreleased\n\n### `Added`\n\n### `Dependencies`\n\n"
        f"{BULLET}\n\n### `Fixed`\n"
    )


def test_changelog_entry_appends_after_existing_bullets():
    text = "## v1dev\n\n### `Dependencies`\n\n- Bump bwa\n- Bump ivar\n\n### `Fixed`\n"
    assert bf.add_changelog_entry(text, OLD, NEW) == (
        f"## v1dev\n\n### `Dependencies`\n\n- Bump bwa\n- Bump ivar\n{BULLET}\n\n### `Fixed`\n"
    )


def test_changelog_creates_dependencies_heading_when_absent():
    text = (
        "# C\n\n## v1 - unreleased\n\n### `Added`\n\n- x\n\n## v0.9 - 2026-01-01\n\n### `Added`\n"
    )
    out = bf.add_changelog_entry(text, OLD, NEW)
    assert out == (
        "# C\n\n## v1 - unreleased\n\n### `Added`\n\n- x\n\n### `Dependencies`\n\n"
        f"{BULLET}\n\n## v0.9 - 2026-01-01\n\n### `Added`\n"
    )


def test_changelog_creates_unreleased_section_when_absent():
    text = "# C\n\n## v0.9 - 2026-01-01\n\n### `Added`\n"
    out = bf.add_changelog_entry(text, OLD, NEW)
    assert out.startswith("# C\n\n## Unreleased\n\n### `Dependencies`\n\n" + BULLET + "\n\n## v0.9")


def test_changelog_created_from_nothing():
    out = bf.add_changelog_entry("", OLD, NEW)
    assert "### `Dependencies`" in out and BULLET in out


def test_changelog_entry_not_duplicated():
    text = "## v1dev\n\n### `Dependencies`\n\n" + BULLET + "\n"
    assert bf.add_changelog_entry(text, OLD, NEW) == text


def test_changelog_only_touches_the_unreleased_section():
    text = "## v1dev\n\n### `Added`\n\n## v0.9\n\n### `Dependencies`\n\n- old\n"
    out = bf.add_changelog_entry(text, OLD, NEW)
    assert "- old\n" in out
    assert out.index(BULLET) < out.index("## v0.9")


# --- CLI ------------------------------------------------------------------------------------


def run(root, *argv, network=None):
    return bf.main([*argv, "--root", str(root)], network=network or FakeNetwork())


def test_cli_exit_codes(tmp_path, capsys):
    up_to_date = make_repo(tmp_path / "a", "2.0.5--pyhdfd78af_0")
    stale = make_repo(tmp_path / "b", "2.0.4--pyhdfd78af_0")
    assert run(up_to_date, "--check") == 0
    assert run(stale, "--check") == 10
    assert run(stale, "--check", network=FakeNetwork(fail="boom")) == 2
    err = capsys.readouterr().err
    assert "boom" in err
    empty = tmp_path / "c"
    empty.mkdir()
    assert run(empty, "--check") == 2


def test_cli_json_output_for_ci(tmp_path, capsys):
    root = make_repo(tmp_path, "2.0.4--pyhdfd78af_0")
    assert run(root, "--check", "--json") == 10
    data = json.loads(capsys.readouterr().out)
    assert data["update_available"] is True
    assert data["major_bump"] is False
    assert data["latest"] == "2.0.5--pyhdfd78af_0"


def test_cli_apply_uses_latest_and_is_a_noop_when_up_to_date(tmp_path):
    root = make_repo(tmp_path, "2.0.4--pyhdfd78af_0")
    assert run(root, "--apply") == 0
    assert {p.value for p in bf.scan_pins(root) if p.kind == "quay"} == {"2.0.5--pyhdfd78af_0"}
    before = snapshot(root)
    assert run(root, "--apply") == 0
    assert snapshot(root) == before


def test_cli_apply_does_nothing_when_blocked(tmp_path):
    root = make_repo(tmp_path, "2.0.5--pyhdfd78af_0")
    before = snapshot(root)
    net = FakeNetwork(quay=with_release(QUAY_REAL, "2.1.0"))
    assert run(root, "--apply", network=net) == 0
    assert snapshot(root) == before


def test_cli_apply_explicit_tag_needs_no_network(tmp_path):
    root = make_repo(tmp_path, "2.0.4--pyhdfd78af_0")
    net = FakeNetwork(fail="network must not be used")
    assert run(root, "--apply", "--tag", "2.0.5--pyhdfd78af_0", network=net) == 0
    assert net.calls == []
    assert run(root, "--apply", "--tag", "latest", network=net) == 2


def test_cli_require_consistent(tmp_path, capsys):
    root = make_repo(tmp_path)
    assert run(root, "--require-consistent") == 0
    (root / "docs" / "usage.md").write_text("quay.io/biocontainers/freyja:2.0.3--pyhdfd78af_0\n")
    assert run(root, "--require-consistent") == 1
    assert "INCONSISTENT" in capsys.readouterr().out


def test_cli_modes_are_mutually_exclusive(tmp_path):
    with pytest.raises(SystemExit):
        run(tmp_path, "--check", "--apply")
    with pytest.raises(SystemExit):
        run(tmp_path)


# --- version stamps in tests ----------------------------------------------------------------


def add_stamps(root, version):
    tests = root / "tests"
    tests.mkdir(exist_ok=True)
    (tests / "x.nf.test").write_text(
        f"assert info.freyja_version == '{version}'\nversions.text = 'P:\\\\n  freyja: {version}\\\\n'\n"
    )
    (tests / "x.nf.test.snap").write_text(f'{{"freyja": "{version}"}}\n')
    # not a nf-test file: a stamp-looking line here must be left alone
    (tests / "test_other.py").write_text(f"freyja: {version}\n")


def test_stamps_in_nf_test_files_are_rewritten_with_the_pin(tmp_path):
    root = make_repo(tmp_path, "2.0.4--pyhdfd78af_0")
    add_stamps(root, "2.0.4")
    assert bf.consistency_problems(bf.scan_pins(root)) == []
    bf.apply_tag(root, bf.parse_tag("2.0.5--pyhdfd78af_0"))
    assert "freyja_version == '2.0.5'" in (root / "tests" / "x.nf.test").read_text()
    assert "freyja: 2.0.5" in (root / "tests" / "x.nf.test").read_text()
    assert '"freyja": "2.0.5"' in (root / "tests" / "x.nf.test.snap").read_text()
    assert (root / "tests" / "test_other.py").read_text() == "freyja: 2.0.4\n"
    assert bf.consistency_problems(bf.scan_pins(root)) == []


def test_stale_stamp_is_reported(tmp_path):
    root = make_repo(tmp_path, "2.0.5--pyhdfd78af_0")
    add_stamps(root, "2.0.4")
    problems = bf.consistency_problems(bf.scan_pins(root))
    assert any("version stamp" in p for p in problems)


# --- workflow invariants ----------------------------------------------------------------------

WORKFLOWS = REPO / ".github" / "workflows"


def _uses(path):
    import re

    return re.findall(r"^\s*(?:-\s+)?uses:\s*(\S+)", path.read_text(), flags=re.MULTILINE)


def test_every_third_party_action_is_pinned_to_a_commit_sha():
    import re

    unpinned = []
    for workflow in sorted(WORKFLOWS.glob("*.yml")):
        for ref in _uses(workflow):
            if ref.startswith("./"):
                continue
            if not re.fullmatch(r"[\w.-]+/[\w./-]+@[0-9a-f]{40}", ref):
                unpinned.append(f"{workflow.name}: {ref}")
    assert not unpinned, unpinned


def test_freyja_update_workflow_gates_auto_merge_on_a_same_major_bump():
    text = (WORKFLOWS / "freyja-update.yml").read_text()
    step = text.split("Enable auto-merge")[1].split("- name:")[0]
    assert "steps.check.outputs.major != 'true'" in step
    assert "gh pr merge --auto --squash" in step
    assert "needs-review" in text
    assert "bin/bump_freyja.py --check --json" in text
    assert "bin/lock_containers.py --update" in text
    # the only secret the workflow reads is the bot token (never the default token for the PR)
    import re

    assert set(re.findall(r"secrets\.(\w+)", text)) == {"BOT_TOKEN"}
