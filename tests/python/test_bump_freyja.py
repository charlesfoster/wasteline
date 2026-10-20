import json
from pathlib import Path

import bump_freyja as bf

FIXTURES = Path(__file__).resolve().parents[1] / "fixtures" / "ci"
REPO = Path(__file__).resolve().parents[2]


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


def test_tag_parse_and_compare():
    """Release tags only; versions compare numerically (2.0.10 > 2.0.9), not lexicographically."""
    names = [t["name"] for t in json.loads((FIXTURES / "quay_tags_odd.json").read_text())["tags"]]
    for raw in ("latest", "2.0.11rc1--pyhdfd78af_0", "2.0.10--pyhdfd78af_1-arm64", ""):
        assert bf.parse_tag(raw) is None
    tag = bf.parse_tag("2.0.5--pyhdfd78af_0")
    assert (tag.version, tag.build, tag.number) == ("2.0.5", "pyhdfd78af", 0)
    assert bf.pick_latest(names).raw == "2.0.10--pyhdfd78af_1"
    assert bf.is_major_bump(bf.parse_tag("1.5.3--pyhdfd78af_2"), bf.parse_tag("2.0.0--x_0"))


def test_real_repo_pins_are_consistent():
    """Invariant (PLAN section 9): every Freyja pin in the repository agrees."""
    pins = bf.scan_pins(REPO)
    assert bf.consistency_problems(pins) == []
    assert set(bf.REQUIRED_KINDS) <= {p.kind for p in pins}


def test_apply_rewrites_every_pin_and_detects_disagreement(tmp_path):
    root = make_repo(tmp_path, "2.0.4--pyhdfd78af_0")
    (root / "docs" / "usage.md").write_text("quay.io/biocontainers/freyja:2.0.9--pyhdfd78af_0\n")
    assert any("quay pins disagree" in p for p in bf.consistency_problems(bf.scan_pins(root)))

    outcome = bf.apply_tag(root, bf.parse_tag("2.1.0--pyhdfd78af_0"))
    assert outcome["applied"]
    pins = bf.scan_pins(root)
    assert bf.consistency_problems(pins) == []
    assert {p.value for p in pins if p.kind == "quay"} == {"2.1.0--pyhdfd78af_0"}
    assert {p.value for p in pins if p.kind == "conda"} == {"2.1.0"}
    assert (
        "bioconda::samtools=1.24"
        in (root / "modules/local/freyja_demix/environment.yml").read_text()
    )
