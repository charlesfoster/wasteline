"""Consistency guards for assets/panels: the registry agrees with conf/pathogens/*.config, and every panel resolves."""

import csv
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
PANELS = ROOT / "assets" / "panels"
FACTS = ["fasta", "gff", "freyja_contig", "primer_scheme_default", "freyja_container_pathogen"]


def read_csv(path):
    with open(path, newline="") as handle:
        return list(csv.DictReader(handle))


def config_facts(path):
    """The params of a pathogen config as the registry writes them (influenza: the --whole_genome branch)."""
    facts = {}
    for key, value in re.findall(r"^\s*(\w+)\s*=\s*(.+)$", path.read_text(), re.MULTILINE):
        value = re.sub(r"\s//.*$", "", value).strip()
        if " ? " in value:
            value = value.split(" ? ", 1)[1].split(" : ", 1)[0].strip()
        value = value.strip("\"'").replace("${projectDir}/", "")
        facts[key] = "" if value == "null" else value
    return facts


def test_registry_agrees_with_pathogen_configs():
    rows = {row["target"]: row for row in read_csv(PANELS / "registry.csv")}
    assert len(rows) == 10
    configs = {}
    for path in (ROOT / "conf" / "pathogens").glob("*.config"):
        facts = config_facts(path)
        if facts.get("pathogen"):
            configs[facts["pathogen"]] = facts
    assert set(rows) == set(configs)
    for name, row in rows.items():
        assert row["freyja_pathogen"] == name
        for fact in FACTS:
            assert row[fact] == configs[name][fact], f"{name}: {fact}"


def test_panels_resolve_against_the_registry():
    rows = {row["target"]: row for row in read_csv(PANELS / "registry.csv")}
    panels = sorted(PANELS.glob("*.csv"))
    assert {p.stem for p in panels} >= {"rsv", "flu"}
    for row in rows.values():
        for column in ("fasta", "gff", "primer_schemes_bundled", "lineage_groups_bundled"):
            assert not row[column] or (ROOT / row[column]).exists(), row[column]
    for panel in panels:
        if panel.name == "registry.csv":
            continue
        targets = [row["target"] for row in read_csv(panel)]
        assert targets and set(targets) <= set(rows), panel.name
