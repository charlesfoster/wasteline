"""Keep main.nf's typed Boolean declarations in step with nextflow.config and the schema.

Nextflow 26.04 (strict parser) delivers every command-line parameter as a String, so a Boolean
parameter can only be set from the command line (`--x`, `--x false`, `--x=false`) if its type is
declared in the `params { }` block of main.nf. That block repeats a default that also lives in
nextflow.config; these tests stop the two drifting (the failure is quiet: the pipeline would behave
one way with the parameter untouched and another through a params file).
"""

import json
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]


def config_params_block():
    text = (ROOT / "nextflow.config").read_text()
    start = text.index("params {")
    depth = 0
    for i in range(start, len(text)):
        if text[i] == "{":
            depth += 1
        elif text[i] == "}":
            depth -= 1
            if depth == 0:
                return text[start : i + 1]
    raise AssertionError("unbalanced params block")


def config_booleans():
    pattern = r"^\s*([A-Za-z_]\w*)\s*=\s*(true|false)\s*(?://.*)?$"
    return {
        m.group(1): m.group(2) == "true"
        for m in re.finditer(pattern, config_params_block(), re.MULTILINE)
    }


def main_nf_booleans():
    block = re.search(
        r"^params \{$(.*?)^\}$", (ROOT / "main.nf").read_text(), re.MULTILINE | re.DOTALL
    )
    assert block, "no top-level `params { ... }` block in main.nf"
    found = {}
    for line in block.group(1).splitlines():
        m = re.match(r"\s*([A-Za-z_]\w*)\s*:\s*Boolean\s*=\s*(true|false)\s*$", line)
        if m:
            found[m.group(1)] = m.group(2) == "true"
    return found


def schema_booleans():
    schema = json.loads((ROOT / "nextflow_schema.json").read_text())
    out = {}
    for group in schema["$defs"].values():
        for name, prop in group["properties"].items():
            if prop.get("type") == "boolean":
                out[name] = prop.get("default")
    return out


# `help` is `--help` or `--help <param>` (single-parameter help), so it stays untyped in main.nf on purpose
UNTYPED_ON_PURPOSE = {"help"}


def test_every_config_boolean_is_typed_in_main_nf():
    missing = set(config_booleans()) - set(main_nf_booleans()) - UNTYPED_ON_PURPOSE
    assert not missing, (
        f"booleans in nextflow.config without a typed declaration in main.nf: {sorted(missing)}"
    )


def test_defaults_agree_between_config_and_main_nf():
    typed = main_nf_booleans()
    wrong = {k: (v, typed[k]) for k, v in config_booleans().items() if k in typed and typed[k] != v}
    assert not wrong, f"default differs (config, main.nf): {wrong}"


def test_every_schema_boolean_is_typed_and_agrees():
    typed = main_nf_booleans()
    schema = schema_booleans()
    assert set(schema) <= set(typed), sorted(set(schema) - set(typed))
    for name, default in schema.items():
        if default is not None:
            assert typed[name] == default, name


def main_nf_numbers():
    block = re.search(
        r"^params \{$(.*?)^\}$", (ROOT / "main.nf").read_text(), re.MULTILINE | re.DOTALL
    )
    found = {}
    for line in block.group(1).splitlines():
        m = re.match(r"\s*([A-Za-z_]\w*)\s*:\s*(Integer|Float)\??\s*=\s*(\S+)\s*$", line)
        if m:
            found[m.group(1)] = (m.group(2), m.group(3))
    return found


def test_every_schema_number_is_typed_and_agrees_with_the_config():
    # numbers also arrive as Strings from the command line (--consensus_min_depth 20)
    typed = main_nf_numbers()
    schema = json.loads((ROOT / "nextflow_schema.json").read_text())
    block = config_params_block()
    for group in schema["$defs"].values():
        for name, prop in group["properties"].items():
            if prop.get("type") not in ("integer", "number"):
                continue
            assert name in typed, f"{name} is a number in the schema but not typed in main.nf"
            kind = "Integer" if prop["type"] == "integer" else "Float"
            assert typed[name][0] == kind, name
            m = re.search(rf"^\s*{name}\s*=\s*([^\s/]+)", block, re.MULTILINE)
            assert m and m.group(1) == typed[name][1], (
                f"{name}: config {m and m.group(1)} vs main.nf {typed[name][1]}"
            )


def test_pathogen_profiles_are_consistent():
    """Every pathogen profile config sets the same per-pathogen facts; nextflow.config maps only to existing files."""
    config = (ROOT / "nextflow.config").read_text()
    line = re.search(r'includeConfig "conf/pathogens/.*', config).group(0)
    profiles = set(re.findall(r"'[^']+': '(\w+)'", line))
    files = {p.stem: p.read_text() for p in (ROOT / "conf" / "pathogens").glob("*.config")}
    assert profiles | {"generic"} == set(files)
    schema = json.loads((ROOT / "nextflow_schema.json").read_text())
    props = {k: v for group in schema["$defs"].values() for k, v in group["properties"].items()}
    for stem, text in files.items():
        for fact in (
            "pathogen",
            "fasta",
            "gff",
            "lineage_groups_bundled",
            "freyja_container_pathogen",
            "primer_schemes_bundled",
            "primer_scheme_default",
            "freyja_contig",
        ):
            assert re.search(rf"^\s*{fact}\s*=", text, re.MULTILINE), (
                f"{stem}.config does not set {fact}"
            )
            assert fact in props, f"{fact} missing from the schema"
        named = re.search(r"^\s*pathogen\s*=\s*'([^']+)'", text, re.MULTILINE)
        if named:
            assert named.group(1) in props["pathogen"]["enum"], stem
        # a bundled reference is one sequence and its bundled annotation is on that sequence
        paths = [
            path
            for line in re.findall(r"^\s*(?:fasta|gff)\s*=.*$", text, re.MULTILINE)
            for path in re.findall(r'"\$\{projectDir\}/([^"]+)"', line)
        ]
        # --whole_genome (PLAN D25): the profile may name a second, multi-record reference + annotation, and the
        # HA contig that Freyja demixes in it (freyja_contig); the HA record comes first
        wg_paths = re.findall(r'"\$\{projectDir\}/([^"]*_wg\.[^"]+)"', text)
        if wg_paths:
            paths = [p for p in paths if p not in wg_paths]
            wg_fasta = ROOT / next(p for p in wg_paths if p.endswith(".fasta"))
            wg_gff = ROOT / next(p for p in wg_paths if p.endswith(".gff3"))
            wg_names = [
                ln[1:].split()[0] for ln in wg_fasta.read_text().splitlines() if ln.startswith(">")
            ]
            ha = re.search(r"freyja_contig\s*=.*\? '([^']+)' : null", text).group(1)
            assert len(wg_names) == 8 and wg_names[0] == ha, f"{wg_fasta.name}: {wg_names}"
            wg_seqids = {
                ln.split("\t")[0] for ln in wg_gff.read_text().splitlines() if ln and ln[0] != "#"
            }
            assert wg_seqids == set(wg_names), f"{wg_gff.name} is not on the sequences {wg_names}"
        gffs = [ROOT / p for p in paths if p.endswith((".gff3", ".gff"))]
        fastas = [ROOT / p for p in paths if not p.endswith((".gff3", ".gff"))]
        for path in fastas + gffs:
            assert path.is_file(), f"{stem}.config names {path}, which does not exist"
        if fastas:
            names = [
                ln[1:].split()[0] for ln in fastas[0].read_text().splitlines() if ln.startswith(">")
            ]
            assert len(names) == 1, f"{fastas[0].name}: {len(names)} sequences"
            for gff in gffs:
                records = [
                    ln.split("\t") for ln in gff.read_text().splitlines() if ln and ln[0] != "#"
                ]
                assert {r[0] for r in records} == set(names), (
                    f"{gff.name} is not on the sequence {names}"
                )
                # bcftools csq needs gene -> mRNA -> exon/CDS ("No usable transcripts" with NCBI's gene + CDS only)
                assert {"gene", "mRNA", "exon", "CDS"} <= {r[2] for r in records}, gff.name
            # bundled primer schemes (PLAN D26): the directory has BED files, each on the reference's sequence, and a
            # default scheme is one of them
            scheme_dir = re.search(
                r'^\s*primer_schemes_bundled\s*=\s*"\$\{projectDir\}/([^"]+)"', text, re.MULTILINE
            )
            default = re.search(r"^\s*primer_scheme_default\s*=\s*'([^']+)'", text, re.MULTILINE)
            if scheme_dir:
                beds = list((ROOT / scheme_dir.group(1)).glob("*.bed"))
                assert beds, f"{stem}.config names a scheme directory without BED files"
                for bed in beds:
                    contigs = {
                        ln.split("\t")[0]
                        for ln in bed.read_text().splitlines()
                        if ln[:1] != "#" and ln
                    }
                    assert contigs == set(names), f"{bed.name} is not on the sequence {names}"
            if default:
                assert (ROOT / scheme_dir.group(1) / f"{default.group(1)}.bed").is_file()
