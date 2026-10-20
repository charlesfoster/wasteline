#!/usr/bin/env python3
"""Tolerance-based end-to-end checks of a wasteline output directory against tests/data/truth.json.

  tests/check_outputs.py <outdir> --mode full|aggregate|report [--truth tests/data/truth.json]
                         [--barcodes tests/data/barcodes_test.csv] [--reference assets/references/NC_045512.fasta]
                         [--protocol capture] [--extras compare,boot,annotate] [--declared-scheme midnight]
                         [--site-aggregation sites.csv] [--lineage-groups present]
  tests/check_outputs.py <outdir> --mode demix --low-confidence S1,S2

Standard library only. Prints one line per check and exits 1 if any check failed. What is checked depends on
the mode (an `aggregate` run has no QC or consensus; a `demix` run is only checked for the low-confidence band):

  demix status (full)   demix_status.tsv, every <sample>.demix_status.json and barcode_version.json agree with
                        truth and carry the barcode version
  aggregate             rows = expected samples only (negatives and skipped samples absent); true lineages
                        within the truth's tolerance in the `lineages`/`abundances` rows (not `summarized`); no
                        lineage absent from the truth above 0.05
  qc (full)             qc_status / qc_reasons per sample
  consensus (full)      one consensus per non-negative sample with the right status; JN.1 SNVs in PURE_C
  report                report_data.json header and run flags, report.html is self-contained, sample IDs intact
  report (mode report)  summary/report.html exists, loads nothing from the network, sample IDs intact
  demix + --low-confidence  the named samples are low_confidence (demixed, in the aggregate and the report, not in
                        the time series), all others keep their expected status
"""

import argparse
import csv
import json
import re
import sys
from pathlib import Path

ABSENT_LINEAGE_MAX = 0.05
FIXED_LINEAGE_COLUMNS = {
    "sample",
    "site",
    "collection_date",
    "status",
    "genome_coverage_pct",
    "status_reason",
    "sample_type",
}
SNV_FRACTION_MIN = 0.9  # share of a lineage's SNVs that must be in its pure sample's consensus
IUPAC_FOR = {"R": "AG", "Y": "CT", "S": "CG", "W": "AT", "K": "GT", "M": "AC"}


class Checker:
    def __init__(self):
        self.failures = 0
        self.total = 0

    def check(self, ok, message, detail=""):
        self.total += 1
        if not ok:
            self.failures += 1
        print(
            f"{'PASS' if ok else 'FAIL'}  {message}"
            + (f"  [{detail}]" if detail and not ok else "")
        )
        return ok


def read_tsv(path):
    with open(path, newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def read_fasta(path):
    records = {}
    name = None
    for line in Path(path).read_text().splitlines():
        if line.startswith(">"):
            name = line[1:].split()[0]
            records[name] = []
        elif name is not None:
            records[name].append(line.strip())
    return {key: "".join(parts) for key, parts in records.items()}


def barcode_snvs(barcodes_csv, lineage):
    """Positions/alt bases (1-based) of the SNVs a lineage carries in the barcode CSV."""
    with open(barcodes_csv, newline="") as handle:
        reader = csv.reader(handle)
        header = next(reader)
        for row in reader:
            if row[0] == lineage:
                snvs = []
                for name, value in zip(header[1:], row[1:], strict=False):
                    match = re.fullmatch(r"([ACGT])(\d+)([ACGT])", name)
                    if match and float(value) == 1.0:
                        snvs.append((int(match.group(2)), match.group(3)))
                return snvs
    return []


def parse_names(text):
    return [t for t in re.split(r"[\s,\[\]']+", text or "") if t]


def parse_floats(text):
    return [float(t) for t in re.findall(r"[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?", text or "")]


def expected_qc(info, have_amplicons):
    """(status, reasons) the truth predicts, given what the run could measure."""
    reasons = set(info["expected_qc_reasons"])
    status = info["expected_qc_status"]
    if not have_amplicons and "amplicon_dropout" in reasons:
        reasons.discard("amplicon_dropout")
        if status == "WARN":
            status = "PASS"
    return status, reasons


def check_demix_status(c, out, samples, mode):
    summary = out / "summary"
    rows = {r["sample"]: r for r in read_tsv(summary / "demix_status.tsv")}
    c.check(
        set(rows) == set(samples),
        "demix_status.tsv lists every sample",
        f"{sorted(set(samples) ^ set(rows))}",
    )
    for sample, info in samples.items():
        row = rows.get(sample)
        if row is None:
            continue
        c.check(
            row["status"] == info["expected_demix_status"],
            f"demix status of {sample} is {info['expected_demix_status']}",
            row["status"],
        )
        if row["status"] != "missing":
            c.check(
                bool(row["barcode_version"]),
                f"demix_status.tsv carries the barcode version for {sample}",
            )
        c.check(
            bool(row["status"]) and (row["status"] == "ok" or bool(row["reason"])),
            f"{sample}: a non-ok status has a reason",
        )
        status_json = out / "samples" / sample / "freyja" / f"{sample}.demix_status.json"
        if status_json.exists():
            record = json.loads(status_json.read_text())
            c.check(
                record.get("status") == row["status"] and bool(record.get("barcode_version")),
                f"{sample}.demix_status.json agrees and has barcode_version",
            )
        else:
            c.check(False, f"{sample}.demix_status.json exists")
    sidecar = summary / "barcode_version.json"
    c.check(
        sidecar.exists() and bool(json.loads(sidecar.read_text()).get("barcode_version")),
        "summary/barcode_version.json has barcode_version",
    )


def check_aggregate(c, out, samples):
    path = out / "summary" / "freyja_aggregated.tsv"
    c.check(path.exists(), "freyja_aggregated.tsv exists")
    if not path.exists():
        return
    with open(path, newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader)
        rows = {r[0]: dict(zip(header, r, strict=False)) for r in reader if r}
    expected = {s for s, i in samples.items() if i["expected_in_aggregate"]}
    c.check(
        set(rows) == expected,
        "aggregate has exactly the expected samples (no negatives, no skipped)",
        f"got {sorted(rows)} want {sorted(expected)}",
    )
    for sample in sorted(set(rows) & set(samples)):
        info = samples[sample]
        names = parse_names(rows[sample].get("lineages"))
        values = parse_floats(rows[sample].get("abundances"))
        c.check(len(names) == len(values) and names, f"{sample}: lineages/abundances rows parse")
        found = dict(zip(names, values, strict=False))
        tol = info.get("abundance_tolerance") or 0.05
        for lineage, truth in info["true_lineage_fractions"].items():
            got = found.get(lineage, 0.0)
            c.check(
                abs(got - truth) <= tol, f"{sample}: {lineage} {got:.3f} within +-{tol} of {truth}"
            )
        extra = {
            k: v
            for k, v in found.items()
            if k not in info["true_lineage_fractions"] and v > ABSENT_LINEAGE_MAX
        }
        c.check(not extra, f"{sample}: no absent lineage above {ABSENT_LINEAGE_MAX}", str(extra))


def check_qc(c, out, samples, mode, have_amplicons):
    path = out / "summary" / "qc_summary.tsv"
    c.check(path.exists(), "qc_summary.tsv exists")
    if not path.exists():
        return
    rows = {r["sample"]: r for r in read_tsv(path)}
    c.check(
        set(rows) == set(samples),
        "qc_summary.tsv has one row per sample incl. negatives",
        f"{sorted(set(samples) ^ set(rows))}",
    )
    for sample, info in samples.items():
        row = rows.get(sample)
        if not row:
            continue
        status, reasons = expected_qc(info, have_amplicons)
        got = {t for t in row["qc_reasons"].split(";") if t}
        c.check(row["qc_status"] == status, f"{sample}: qc_status {status}", row["qc_status"])
        c.check(got == reasons, f"{sample}: qc_reasons {sorted(reasons)}", str(sorted(got)))
    if have_amplicons:
        table = out / "summary" / "amplicon_coverage.tsv"
        c.check(table.exists(), "amplicon_coverage.tsv exists")
        if table.exists():
            dropped = {}
            for r in read_tsv(table):
                if r["dropped"] == "true":
                    dropped.setdefault(r["sample"], set()).add(r["amplicon"])
            for sample, info in samples.items():
                if info["expected_demix_status"] == "ok" and info["sample_type"] == "sample":
                    c.check(
                        dropped.get(sample, set()) == set(info["dropout_amplicons"]),
                        f"{sample}: dropped amplicons {info['dropout_amplicons']}",
                        str(sorted(dropped.get(sample, []))),
                    )


def check_consensus(c, out, samples, barcodes, reference):
    stats_path = out / "summary" / "consensus_stats.tsv"
    c.check(stats_path.exists(), "consensus_stats.tsv exists")
    fasta_path = out / "summary" / "consensus_genomes.fasta"
    c.check(fasta_path.exists(), "consensus_genomes.fasta exists")
    if not (stats_path.exists() and fasta_path.exists()):
        return
    header = stats_path.read_text().splitlines()[0].split("\t")
    c.check(
        header == ["sample", "length", "n_count", "n_pct", "ambiguous_count", "status", "reason"],
        "consensus_stats.tsv has 7 columns, reason last",
        str(header),
    )
    stats = {r["sample"]: r for r in read_tsv(stats_path)}
    fasta = read_fasta(fasta_path)
    for sample, info in samples.items():
        want = info["expected_consensus_status"]
        if want == "none":
            c.check(
                sample not in stats and sample not in fasta,
                f"{sample}: negative control has no consensus",
            )
            continue
        c.check(
            stats.get(sample, {}).get("status") == want,
            f"{sample}: consensus status {want}",
            str(stats.get(sample)),
        )
        c.check(sample in fasta, f"{sample}: consensus FASTA record present")
        status_json = out / "samples" / sample / "consensus" / f"{sample}.consensus_status.json"
        if status_json.exists():
            c.check(
                isinstance(json.loads(status_json.read_text()).get("warnings"), list),
                f"{sample}: consensus_status.json has a warnings array",
            )
        c.check(
            (out / "samples" / sample / "consensus" / f"{sample}.consensus.fa").exists(),
            f"{sample}: samples/{sample}/consensus/{sample}.consensus.fa published",
        )
    pure = next(
        (
            sid
            for sid, i in samples.items()
            if i["true_lineage_fractions"] == {"JN.1": 1.0} and sid in fasta
        ),
        None,
    )
    if pure and barcodes and Path(barcodes).exists():
        snvs = barcode_snvs(barcodes, "JN.1")
        seq = fasta[pure]
        hit = sum(
            1
            for pos, alt in snvs
            if pos <= len(seq) and (seq[pos - 1] == alt or alt in IUPAC_FOR.get(seq[pos - 1], ""))
        )
        c.check(
            snvs and hit / len(snvs) >= SNV_FRACTION_MIN,
            f"{pure} consensus carries the JN.1 SNVs ({hit}/{len(snvs)})",
        )
    if reference and Path(reference).exists() and pure:
        ref = "".join(
            line for line in Path(reference).read_text().splitlines() if not line.startswith(">")
        )
        c.check(len(fasta[pure]) == len(ref), f"{pure} consensus has the reference length")


def check_analyses(c, out, samples, mode, have_amplicons, declared_scheme):
    """QC checks and surveillance tables added on top of the demix (summary/qc_summary.tsv, dropout_impact,
    unexplained_mutations, freyja_timeseries). `declared_scheme`: primer scheme the run declared (mode full)."""
    summary = out / "summary"
    qc_path = summary / "qc_summary.tsv"
    if not qc_path.exists():
        return
    qc_lines = qc_path.read_text().splitlines()
    c.check(
        len(qc_lines[0].split("\t")) == 27,
        "qc_summary.tsv has 27 columns",
        str(len(qc_lines[0].split("\t"))),
    )
    qc = {r["sample"]: r for r in read_tsv(qc_path)}

    if mode == "full" and have_amplicons and declared_scheme:
        # a scheme that matches the reads: `ok`, the declared scheme among the best matches, trimming happened
        for sample, info in samples.items():
            if info["sample_type"] != "sample" or info["expected_demix_status"] != "ok":
                continue
            row = qc.get(sample, {})
            c.check(
                row.get("primer_scheme_check") == "ok",
                f"{sample}: primer scheme check is ok for the declared {declared_scheme} scheme",
                str(row.get("primer_scheme_check")),
            )
            c.check(
                declared_scheme in (row.get("primer_scheme_best") or "").split("/"),
                f"{sample}: the declared scheme is among the best-matching schemes",
                str(row.get("primer_scheme_best")),
            )
            summary_path = (
                out / "samples" / sample / "coverage" / f"{sample}.primer_scheme_summary.tsv"
            )
            if summary_path.exists():
                scheme = read_tsv(summary_path)[0]
                hit = scheme.get("declared_hit_fraction", "NA")
                c.check(
                    hit not in ("", "NA") and float(hit) >= 0.05,
                    f"{sample}: the declared scheme's hit fraction is above the floor (0.05)",
                    str(hit),
                )
            pct = row.get("primer_trimmed_pct", "NA")
            c.check(
                pct not in ("", "NA") and float(pct) > 0,
                f"{sample}: primer_trimmed_pct is filled from the ivar log",
                pct,
            )
        c.check(
            (out / "samples" / "MIX_A" / "coverage" / "MIX_A.primer_scheme_summary.tsv").exists()
            if "MIX_A" in samples
            else True,
            "per-sample primer scheme check tables are published",
        )

    if have_amplicons:
        path = summary / "dropout_impact.tsv"
        c.check(path.exists(), "dropout_impact.tsv exists")
        c.check(
            (summary / "dropout_impact_summary.tsv").exists(), "dropout_impact_summary.tsv exists"
        )
        if path.exists():
            rows = read_tsv(path)
            c.check(
                {"n_defining_sites", "n_differential_sites"} <= set(rows[0]) if rows else True,
                "dropout_impact.tsv has the defining-site columns",
            )
            summary_rows = read_tsv(summary / "dropout_impact_summary.tsv")
            c.check(
                {"defining_sites_total", "defining_sites_lost", "pct_defining_sites_lost"}
                <= set(
                    summary_rows[0]
                    if summary_rows
                    else {"defining_sites_total", "defining_sites_lost", "pct_defining_sites_lost"}
                ),
                "dropout_impact_summary.tsv has the defining-site columns",
            )
            for sample, info in samples.items():
                if info["expected_demix_status"] != "ok" or info["sample_type"] != "sample":
                    continue
                got = {r["amplicon"] for r in rows if r["sample"] == sample}
                c.check(
                    got == set(info["dropout_amplicons"]),
                    f"dropout_impact.tsv: {sample} lists exactly its dropped amplicons {info['dropout_amplicons']}",
                    str(sorted(got)),
                )
            if "Site_1_week_2" in samples:
                by_amplicon = {r["amplicon"]: r for r in rows if r["sample"] == "Site_1_week_2"}
                c.check(
                    by_amplicon.get("nCoV-2019_22", {}).get("affects_detected_lineage") == "true",
                    "dropout impact: nCoV-2019_22 (spike) matters for the lineages detected in Site_1_week_2",
                )
                c.check(
                    "XBB.1.5" in by_amplicon.get("nCoV-2019_22", {}).get("affected_lineages", ""),
                    "dropout impact: nCoV-2019_22 loses XBB.1.5-defining sites",
                )
                c.check(
                    int(by_amplicon.get("nCoV-2019_22", {}).get("n_defining_sites") or 0) > 0
                    and by_amplicon.get("nCoV-2019_7", {}).get("n_defining_sites") == "0",
                    "dropout impact: defining sites are lost in nCoV-2019_22 and none in nCoV-2019_7",
                )
                c.check(
                    by_amplicon.get("nCoV-2019_7", {}).get("affects_detected_lineage") == "false",
                    "dropout impact: nCoV-2019_7 carries no site that matters for the detected lineages",
                )

    if mode in ("full", "demix"):
        path = summary / "unexplained_mutations.tsv"
        c.check(path.exists(), "unexplained_mutations.tsv exists")
        if path.exists():
            c.check(
                path.read_text().splitlines()[0].startswith("sample\tposition\tref\talt"),
                "unexplained_mutations.tsv has its header",
            )
            c.check(
                not read_tsv(path),
                "synthetic samples have no unexplained mutations at the default thresholds",
                str(len(read_tsv(path))),
            )

    dates = {
        row["collection_date"]
        for row in qc.values()
        if row["sample_type"] != "negative"
        and row["demix_status"] == "ok"
        and row["collection_date"] not in ("", "NA")
    }
    for name in ("freyja_timeseries.tsv", "freyja_timeseries_lineages.tsv"):
        path = summary / name
        if len(dates) >= 2:
            c.check(path.exists(), f"{name} exists (>= 2 collection dates)")
            if path.exists():
                got = {r["collection_date"] for r in read_tsv(path)}
                c.check(
                    got == dates,
                    f"{name}: one set of rows per collection date",
                    str(sorted(got ^ dates)),
                )
        else:
            c.check(not path.exists(), f"{name} is not written with fewer than 2 collection dates")


def check_lineage_groups(c, out, samples, expect):
    """summary/lineage_summary.csv and its companions. expect: `present` (the groups file ran against a dataset with
    a lineage hierarchy) or `absent` (custom barcodes file: skipped, or --lineage_summary false)."""
    summary = out / "summary"
    names = [
        "lineage_summary.csv",
        "lineage_summary.xlsx",
        "lineage_groups_resolved.tsv",
        "lineage_groups_used.csv",
        "lineage_ungrouped_summary.tsv",
    ]
    if expect == "absent":
        c.check(
            not any((summary / n).exists() for n in names),
            "no lineage summary files (no lineage hierarchy, or switched off)",
        )
        return
    for name in names:
        c.check((summary / name).exists(), f"{name} exists")
    if not all((summary / n).exists() for n in names):
        return
    import zipfile

    try:
        with zipfile.ZipFile(summary / "lineage_summary.xlsx") as workbook:
            bad_member = workbook.testzip()
            sheets = re.findall(r'<sheet name="([^"]+)"', workbook.read("xl/workbook.xml").decode())
    except (zipfile.BadZipFile, KeyError, OSError) as exc:
        bad_member, sheets = str(exc), []
    c.check(
        bad_member is None
        and {"Lineage summary", "Ungrouped lineages", "Groups used", "About"} <= set(sheets),
        "lineage_summary.xlsx is a valid workbook with its sheets",
        f"{bad_member} {sheets}",
    )
    default_groups = Path(__file__).resolve().parents[1] / "assets" / "lineage_groups.csv"
    c.check(
        (summary / "lineage_groups_used.csv").read_bytes() == default_groups.read_bytes(),
        "lineage_groups_used.csv is a byte-identical copy of the groups file used (the bundled default)",
    )
    rows = {
        r["sample"]: r
        for r in csv.DictReader((summary / "lineage_summary.csv").read_text().splitlines())
    }
    c.check(
        not any(
            sample in rows for sample, info in samples.items() if info["sample_type"] == "negative"
        ),
        "negative controls are not in lineage_summary.csv",
    )
    for sample, info in samples.items():
        row = rows.get(sample)
        if info["sample_type"] == "negative" or row is None:
            continue
        if info["expected_demix_status"] == "ok" and row["status"] == "ok":
            total = float(row["proportion_summary"])
            c.check(0.97 <= total <= 1.03, f"{sample}: proportion_summary is about 1", str(total))
        else:
            values = [v for k, v in row.items() if k not in FIXED_LINEAGE_COLUMNS]
            c.check(
                row["status"] != "ok" and not any(values),
                f"{sample}: a status row without values ({row['status']})",
            )
    if "PURE_C" in rows and rows["PURE_C"]["status"] == "ok":
        jn = float(rows["PURE_C"]["JN.1* (excludes KP.2*, KP.3*, KW.1.1*, LP.8.1*)"] or 0)
        both = float(rows["PURE_C"]["JN.1* (includes KP.2*, KP.3*, KW.1.1*, LP.8.1*)"] or 0)
        c.check(
            abs(jn - 1.0) <= 0.05 and abs(both - 1.0) <= 0.05,
            "PURE_C (JN.1): both JN.1 groups carry ~1.0",
            f"{jn} {both}",
        )
    if "MIX_A" in rows and rows["MIX_A"]["status"] == "ok":
        omicron = float(rows["MIX_A"]["Omicron Lineage (B.1.1.529*)"] or 0)
        c.check(
            abs(omicron - 0.71) <= 0.05,
            "MIX_A: BA.2 (~0.71) is in 'Omicron Lineage (B.1.1.529*)'",
            str(omicron),
        )
        c.check("AY.4:" in rows["MIX_A"]["ungrouped_lineages"], "MIX_A: AY.4 (Delta) is ungrouped")
        ungrouped = {r["lineage"]: r for r in read_tsv(summary / "lineage_ungrouped_summary.tsv")}
        c.check(
            ungrouped.get("AY.4", {}).get("suggested_pattern") == "AY.4*",
            "lineage_ungrouped_summary.tsv lists AY.4 with the suggested pattern AY.4*",
        )
    resolved = read_tsv(summary / "lineage_groups_resolved.tsv")
    c.check(
        resolved and not any("not in the dataset" in r["warnings"] for r in resolved),
        "the bundled default groups have no unknown lineage in the dataset",
    )


def check_site_aggregates(c, out, sites_file, period="week"):
    """summary/site_aggregates_*.tsv: present, one meta row per freyja row, and no lineage in an aggregate that none of
    its member sites' samples of that period has (a weighted average of results cannot invent one)."""
    import datetime

    summary = out / "summary"
    freyja_path, meta_path = (
        summary / "site_aggregates_freyja.tsv",
        summary / "site_aggregates_meta.tsv",
    )
    c.check(
        freyja_path.exists() and meta_path.exists(),
        "site_aggregates_freyja.tsv and _meta.tsv exist",
    )
    if not (freyja_path.exists() and meta_path.exists()):
        return
    members = {}
    for row in csv.DictReader(Path(sites_file).read_text(encoding="utf-8-sig").splitlines()):
        members.setdefault(row["aggregate_name"].strip(), set()).add(
            row["site_name"].strip().casefold()
        )
    sheet = {
        r["sample"]: r
        for r in csv.DictReader((summary / "samplesheet.csv").read_text().splitlines())
    }
    lineages = {}
    for line in (summary / "freyja_aggregated.tsv").read_text().splitlines()[1:]:
        cols = line.split("\t")
        lineages[cols[0].replace(".variants.tsv", "")] = set(cols[2].split())

    def week(date):
        year, number, _ = datetime.date.fromisoformat(date).isocalendar()
        return f"{year}-W{number:02d}" if period == "week" else date

    meta = {r["row_label"]: r for r in read_tsv(meta_path)}
    rows = [line.split("\t") for line in freyja_path.read_text().splitlines()[1:]]
    c.check(
        {r[0] for r in rows} == set(meta) and rows,
        "every site aggregate row has a meta row",
        str(sorted(set(meta))),
    )
    bad = []
    for cols in rows:
        record = meta.get(cols[0])
        if not record:
            continue
        pool = set()
        for sample, info in sheet.items():
            if (
                info.get("sample_type") != "negative"
                and info.get("site", "").strip().casefold() in members[record["aggregate_name"]]
                and info.get("collection_date")
                and week(info["collection_date"]) == record["period"]
            ):
                pool |= lineages.get(sample, set())
        if not set(cols[2].split()) <= pool:
            bad.append((cols[0], sorted(set(cols[2].split()) - pool)))
    c.check(not bad, "no aggregate lineage is missing from all of its member sites", str(bad))


def check_published_qc_is_complete(c, out, samples):
    """summary/qc_summary.tsv carries the downstream columns (not just NA placeholders)."""
    path = out / "summary" / "qc_summary.tsv"
    if not path.exists():
        return
    rows = {r["sample"]: r for r in read_tsv(path)}
    for sample, info in samples.items():
        row = rows.get(sample)
        if not row:
            continue
        c.check(
            row["demix_status"] == info["expected_demix_status"],
            f"qc_summary.tsv: {sample} demix_status {info['expected_demix_status']}",
            row["demix_status"],
        )
        if info["expected_consensus_status"] != "none":
            c.check(
                row["consensus_length"] not in ("", "NA"),
                f"qc_summary.tsv: {sample} consensus_length filled",
            )
        if info["expected_demix_status"] != "missing":
            c.check(
                row["barcode_version"] not in ("", "NA"),
                f"qc_summary.tsv: {sample} barcode_version filled",
            )


def check_negative_controls(c, out, samples):
    """D8: controls that demixed ok are listed in their own table and report block, not in the aggregate."""
    controls = [
        s
        for s, i in samples.items()
        if i["sample_type"] == "negative" and i["expected_demix_status"] == "ok"
    ]
    path = out / "summary" / "freyja_aggregated_controls.tsv"
    if not controls:
        c.check(
            not path.exists(), "no negative control demixed ok: no freyja_aggregated_controls.tsv"
        )
        return
    c.check(path.exists(), "freyja_aggregated_controls.tsv exists")
    if path.exists():
        with open(path, newline="") as handle:
            reader = csv.reader(handle, delimiter="\t")
            next(reader)
            got = {r[0] for r in reader if r}
        c.check(
            got == set(controls),
            "freyja_aggregated_controls.tsv has exactly the ok negative controls",
            str(sorted(got)),
        )
    data = json.loads((out / "summary" / "report_data.json").read_text())
    listed = {n["sample"]: n for n in data.get("negative_controls", [])}
    c.check(
        set(controls) <= set(listed),
        "report_data.json negative_controls lists the contaminated control",
    )
    for sample in controls:
        entry = listed.get(sample, {})
        c.check(
            entry.get("excluded_from_aggregate") is True and bool(entry.get("lineages")),
            f"{sample}: excluded from the aggregate, lineages shown",
        )
        c.check(
            data["demix"][sample]["status"] == "ok",
            f"{sample}: demix status stays ok in report_data.json (never downgraded)",
        )
    html = (out / "summary" / "report.html").read_text()
    c.check(
        'id="sec-controls"' in html and "excluded from the aggregate" in html,
        "report.html has the Negative controls block",
    )


def check_extras(c, out, samples, extras, barcodes, reference):
    ok_samples = {s for s, i in samples.items() if i["expected_demix_status"] == "ok"}
    html = (out / "summary" / "report.html").read_text()
    if "compare" in extras:
        path = out / "summary" / "barcode_compare.tsv"
        c.check(path.exists(), "barcode_compare.tsv exists")
        if path.exists():
            rows = read_tsv(path)
            c.check(
                list(rows[0])
                == [
                    "sample",
                    "lineage",
                    "abundance_a",
                    "abundance_b",
                    "barcode_version_a",
                    "barcode_version_b",
                ],
                "barcode_compare.tsv has the contract header",
            )
            c.check(
                {r["sample"] for r in rows} == ok_samples,
                "barcode_compare.tsv covers every demixed sample",
            )
            c.check(
                any(r["abundance_a"] != r["abundance_b"] for r in rows),
                "the two barcode datasets give different abundances",
            )
            c.check(
                all(r["barcode_version_a"] != r["barcode_version_b"] for r in rows),
                "both barcode versions are recorded and differ",
            )
        info = json.loads((out / "pipeline_info" / "run_info.json").read_text())
        c.check(
            bool((info.get("freyja_data_compare") or {}).get("barcode_version")),
            "run_info.json freyja_data_compare recorded",
        )
        c.check('id="sec-barcode"' in html, "report.html renders the barcode comparison panel")
    if "boot" in extras:
        path = out / "summary" / "freyja_boot_summary.tsv"
        c.check(path.exists(), "freyja_boot_summary.tsv exists")
        if path.exists():
            rows = read_tsv(path)
            c.check(
                {r["sample"] for r in rows} == ok_samples,
                "freyja_boot_summary.tsv covers every demixed sample",
            )
            c.check(
                all(
                    float(r["lo"]) <= float(r["abundance"]) + 1e-6 <= float(r["hi"]) + 2e-6
                    for r in rows
                    if float(r["abundance"]) > 0
                ),
                "every point estimate lies inside its bootstrap interval",
            )
            pure = {r["lineage"]: r for r in rows if r["sample"] == "PURE_C"}
            c.check(
                "JN.1" in pure and float(pure["JN.1"]["lo"]) > 0.9,
                "PURE_C: JN.1 interval is narrow and high",
            )
        c.check(
            "95% interval" in html,
            "report.html shows the bootstrap interval in the lineage tooltips",
        )
    if "annotate" in extras:
        path = out / "summary" / "annotated_variants.tsv"
        c.check(path.exists(), "annotated_variants.tsv exists")
        if path.exists():
            rows = read_tsv(path)
            c.check(
                list(rows[0])
                == [
                    "reference",
                    "sample_id",
                    "gene",
                    "nt_pos",
                    "ref",
                    "alt",
                    "effect",
                    "aa_bcsq",
                    "aa_standard",
                    "vaf",
                    "depth",
                    "compound",
                ],
                "annotated_variants.tsv has the old pipeline's columns plus the appended compound",
            )
            c.check(
                any(
                    r["sample_id"] == "PURE_C" and r["gene"] == "S" and r["aa_standard"] == "N501Y"
                    for r in rows
                ),
                "PURE_C (JN.1): S:N501Y is annotated",
            )
            c.check(
                not any(samples[r["sample_id"]]["sample_type"] == "negative" for r in rows),
                "negative controls are not annotated",
            )
        for sample, info in samples.items():
            if info["expected_consensus_status"] != "none":
                c.check(
                    (out / "samples" / sample / "variants" / f"{sample}.annotated.tsv").exists(),
                    f"{sample}: samples/{sample}/variants/{sample}.annotated.tsv published",
                )


def check_capture(c, out, samples):
    """Probe-capture data (PLAN D13): no primers, no amplicon analysis, duplicates marked (the default)."""
    info = json.loads((out / "pipeline_info" / "run_info.json").read_text())
    c.check(info.get("protocol") == "capture", "run_info.json records protocol capture")
    c.check(info.get("primer_trimming") is False, "run_info.json: primer_trimming false")
    c.check(info.get("duplicates") == "mark", "run_info.json: duplicates mark")
    c.check(
        not (out / "summary" / "amplicon_coverage.tsv").exists()
        and not list((out / "samples").glob("*/coverage/*amplicon_coverage.tsv")),
        "no amplicon_coverage outputs",
    )
    c.check(
        not list((out / "samples").glob("*/bam/*primertrimmed*"))
        and not list((out / "samples").glob("*/bam/*ivar*")),
        "no primer-trimmed BAM or ivar trim log (nothing was trimmed)",
    )
    traces = sorted((out / "pipeline_info").glob("execution_trace_*.txt"))
    if traces:
        names = {r["process"].split(":")[-1] for r in read_tsv(traces[-1])}
        c.check(
            not names & {"IVAR_TRIM", "AMPLICON_COVERAGE"},
            "trace: IVAR_TRIM and AMPLICON_COVERAGE never ran",
            str(sorted(names & {"IVAR_TRIM", "AMPLICON_COVERAGE"})),
        )
    qc = {r["sample"]: r for r in read_tsv(out / "summary" / "qc_summary.tsv")}
    c.check(
        all(r["amplicons_total"] == "NA" and r["amplicons_dropped"] == "NA" for r in qc.values()),
        "qc_summary.tsv: amplicons_total/amplicons_dropped are NA",
    )
    c.check(
        all("amplicon_dropout" not in r["qc_reasons"] for r in qc.values()),
        "no amplicon_dropout reason anywhere",
    )
    html = (out / "summary" / "report.html").read_text()
    c.check(
        "Not applicable for protocol" in html and "Library protocol" in html,
        "report.html: protocol in the header and the amplicon section marked not applicable",
    )
    for sample, truth in samples.items():
        stats = out / "samples" / sample / "bam" / f"{sample}.markdup_stats.txt"
        c.check(
            stats.exists(), f"{sample}: samples/{sample}/bam/{sample}.markdup_stats.txt published"
        )
        rng = truth.get("expected_duplicate_fraction_range")
        if stats.exists() and rng:
            fields = dict(
                line.split(": ", 1) for line in stats.read_text().splitlines() if ": " in line
            )
            fraction = int(fields["DUPLICATE PAIR"]) / max(int(fields["EXAMINED"]), 1)
            c.check(
                rng[0] <= fraction <= rng[1],
                f"{sample}: duplicate fraction {fraction:.3f} within {rng}",
            )


def check_report(c, out, samples, mode):
    summary = out / "summary"
    data_path, html_path = summary / "report_data.json", summary / "report.html"
    c.check(data_path.exists() and html_path.exists(), "report_data.json and report.html exist")
    if not (data_path.exists() and html_path.exists()):
        return
    data = json.loads(data_path.read_text())
    html = html_path.read_text()
    run_info = (
        json.loads((out / "pipeline_info" / "run_info.json").read_text())
        if (out / "pipeline_info" / "run_info.json").exists()
        else {}
    )
    barcode = (run_info.get("freyja_data") or {}).get("barcode_version")
    if mode != "aggregate" and barcode:
        c.check(
            data["header"].get("barcode_version") == barcode,
            "report_data.json header has the barcode version",
        )
        c.check(barcode in html, "report.html header shows the barcode version")
    external = re.findall(
        r"""(?:src|href)\s*=\s*["']\s*(?:https?:)?//|url\(\s*["']?\s*(?:https?:)?//|@import\s+["']?https?:""",
        html,
    )
    c.check(not external, "report.html loads nothing from the network", str(external[:3]))
    c.check(
        "http://" not in html and "https://" not in html,
        "report.html contains no http(s):// URL at all",
    )
    if "Site_1_week_2" in samples:
        c.check(
            "Site_1_week_2" in html and "Site-1-week-2" not in html,
            "report.html: sample ID Site_1_week_2 intact",
        )
    codes = {f.get("code"): f for f in data.get("run_flags", [])}
    if mode in ("full", "demix") and any(
        "ntc_reads" in i["expected_qc_reasons"] for i in samples.values()
    ):
        flag = codes.get("ntc_reads")
        c.check(
            flag is not None and flag["level"] == "error" and flag.get("samples") == ["NEG_CONTAM"],
            "run flag ntc_reads (error) names NEG_CONTAM",
            str(flag),
        )
    report_samples = {s["sample"] for s in data.get("samples", [])}
    expected = (
        set(samples)
        if mode in ("full", "demix")
        else {s for s, i in samples.items() if i["expected_in_aggregate"]}
    )
    c.check(
        expected <= report_samples,
        "report_data.json lists the expected samples",
        str(sorted(expected - report_samples)),
    )


def check_ids_intact(c, out):
    """Sample IDs are never rewritten: no underscore-to-dash variant anywhere in the text outputs."""
    bad = []
    for path in sorted((out / "summary").glob("*")) + sorted(
        (out / "pipeline_info").glob("run_info.json")
    ):
        if path.suffix in {
            ".tsv",
            ".json",
            ".html",
            ".fasta",
            ".csv",
        } and "Site-1-week-2" in path.read_text(errors="replace"):
            bad.append(path.name)
    c.check(not bad, "Site_1_week_2 is never rewritten to Site-1-week-2", str(bad))


def check_run_info(c, out, mode):
    path = out / "pipeline_info" / "run_info.json"
    c.check(path.exists(), "pipeline_info/run_info.json exists")
    if not path.exists():
        return
    info = json.loads(path.read_text())
    keys = [
        "pipeline_name",
        "pipeline_version",
        "git_commit",
        "nextflow_version",
        "mode",
        "start",
        "params",
        "freyja_version",
        "tool_versions",
        "previous_runs",
        "freyja_data",
        "input_checksums",
        "provenance_generated",
    ]
    c.check(
        all(k in info for k in keys),
        "run_info.json has the contract keys",
        str([k for k in keys if k not in info]),
    )
    c.check(info.get("mode") == mode, f"run_info.json mode is {mode}", str(info.get("mode")))
    c.check(
        bool((info.get("freyja_data") or {}).get("barcode_version")),
        "run_info.json freyja_data.barcode_version present",
    )
    if mode in ("demix", "aggregate"):
        c.check(
            isinstance(info.get("previous_runs"), list) and info["previous_runs"],
            "run_info.json records previous_runs",
        )


def check_low_confidence(c, out, samples, low):
    """Run with the low-confidence band on (PLAN D18): `low` are the samples that must be `low_confidence`.

    All other samples keep their expected status. Low-confidence samples are in the aggregate and the report but,
    by default (no --include_low_confidence), not in the time series.
    """
    summary = out / "summary"
    rows = {r["sample"]: r for r in read_tsv(summary / "demix_status.tsv")}
    for sample, info in samples.items():
        row = rows.get(sample, {})
        want = "low_confidence" if sample in low else info["expected_demix_status"]
        c.check(row.get("status") == want, f"demix status of {sample} is {want}", row.get("status"))
        has_demix = (out / "samples" / sample / "freyja" / f"{sample}.demix").exists()
        c.check(
            has_demix == (want in ("ok", "low_confidence")),
            f"{sample}: a .demix file exists only for ok and low_confidence",
        )
        if sample in low:
            c.check(
                "high-confidence threshold" in row.get("reason", ""),
                f"{sample}: the low_confidence reason names the threshold",
                row.get("reason", ""),
            )
    sheet = {
        r["sample"]: r
        for r in csv.DictReader((summary / "samplesheet.csv").read_text().splitlines())
    }
    aggregated = {r[""] for r in read_tsv(summary / "freyja_aggregated.tsv")}
    series = read_tsv(summary / "freyja_timeseries.tsv")
    ts_points = {(r["site"], r["collection_date"]) for r in series}
    data = json.loads((summary / "report_data.json").read_text())
    html = (summary / "report.html").read_text()
    real = sorted(s for s in low if samples[s]["sample_type"] != "negative")
    for sample in real:
        point = (sheet[sample]["site"], sheet[sample]["collection_date"])
        c.check(sample in aggregated, f"{sample} is in freyja_aggregated.tsv (never filtered)")
        c.check(
            point not in ts_points,
            f"{sample} is not in freyja_timeseries.tsv",
            str(sorted(ts_points)),
        )
        c.check(
            data["demix"][sample]["status"] == "low_confidence" and sample in html,
            f"report: {sample} is shown with status low_confidence",
        )
    flag = next((f for f in data["run_flags"] if f["code"] == "demix_low_confidence"), None)
    c.check(
        flag is not None and flag["level"] == "info" and flag.get("samples") == real,
        "report: one info flag demix_low_confidence naming the non-control low_confidence samples",
        str(flag),
    )
    c.check(
        bool(low) and "low confidence</span>" in html,
        "report.html shows the low confidence badge",
    )


def main(argv=None):
    here = Path(__file__).resolve().parent
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("outdir")
    parser.add_argument("--mode", required=True, choices=["full", "demix", "aggregate", "report"])
    parser.add_argument(
        "--protocol",
        choices=["amplicon", "capture"],
        default="amplicon",
        help="capture: truth defaults to tests/data_capture and the capture-specific checks run",
    )
    parser.add_argument("--truth", default=None, help="default: tests/data[_capture]/truth.json")
    parser.add_argument("--barcodes", default=str(here / "data" / "barcodes_test.csv"))
    parser.add_argument(
        "--reference", default=str(here.parent / "assets" / "references" / "NC_045512.fasta")
    )
    parser.add_argument(
        "--extras",
        default="",
        help="comma list of opt-in features that were on: compare,boot,annotate",
    )
    parser.add_argument(
        "--declared-scheme",
        default=None,
        help="mode full: the primer scheme the run declared; its QC rows must then have primer_scheme_check ok",
    )
    parser.add_argument(
        "--site-aggregation",
        default=None,
        help="the sites CSV of the run: the site aggregates must then exist and be consistent with their members",
    )
    parser.add_argument(
        "--lineage-groups",
        choices=["present", "absent"],
        default="absent",
        help="present: the run had a Freyja dataset with a lineage hierarchy, so the lineage summary files must exist",
    )
    parser.add_argument(
        "--low-confidence",
        default="",
        help="comma list of samples that must be low_confidence: run the band checks (with --mode demix) instead of the standard ones",
    )
    args = parser.parse_args(argv)
    if (args.mode == "demix") != bool(args.low_confidence):
        parser.error("--mode demix is only checked with --low-confidence (and vice versa)")

    out = Path(args.outdir)
    default_truth = here / ("data_capture" if args.protocol == "capture" else "data") / "truth.json"
    truth = json.loads(Path(args.truth or default_truth).read_text())
    samples = truth["samples"]
    c = Checker()

    if args.mode == "report":
        # report mode writes only summary/report.html (+ report_data.json when it had to rebuild it)
        html_path = out / "summary" / "report.html"
        c.check(html_path.exists(), "summary/report.html exists")
        if html_path.exists():
            html = html_path.read_text()
            c.check(
                ("Site_1_week_2" in html or "Site_1_week_2" not in samples)
                and "https://" not in html
                and "http://" not in html,
                "rebuilt report.html: IDs intact, no URLs",
            )
        print(f"\n{c.total - c.failures}/{c.total} checks passed")
        return 0 if c.failures == 0 else 1

    check_run_info(c, out, args.mode)
    if args.mode == "demix":
        low = {s for s in args.low_confidence.split(",") if s}
        check_low_confidence(c, out, samples, low)
        print(f"\n{c.total - c.failures}/{c.total} checks passed")
        return 0 if c.failures == 0 else 1
    if args.mode == "full":
        have_amplicons = args.protocol == "amplicon"
        check_demix_status(c, out, samples, args.mode)
        check_aggregate(c, out, samples)
        check_qc(c, out, samples, args.mode, have_amplicons)
        check_consensus(c, out, samples, args.barcodes, args.reference)
        check_published_qc_is_complete(c, out, samples)
        check_analyses(c, out, samples, args.mode, have_amplicons, args.declared_scheme)
        if args.protocol == "capture":
            check_capture(c, out, samples)
    else:
        check_aggregate(c, out, samples)
        c.check(
            not (out / "summary" / "barcode_version.json").exists(),
            "aggregate mode writes no barcode_version.json sidecar",
        )
    check_report(c, out, samples, args.mode)
    check_negative_controls(c, out, samples)
    check_lineage_groups(c, out, samples, args.lineage_groups)
    if args.site_aggregation:
        check_site_aggregates(c, out, args.site_aggregation)
    if args.mode == "aggregate" and "MIX_A" in samples:
        data = json.loads((out / "summary" / "report_data.json").read_text())
        by_sample = {s["sample"]: s for s in data["samples"]}
        c.check(
            by_sample.get("MIX_A", {}).get("collection_date") == "2026-01-05",
            "aggregate report: sample metadata (collection_date) comes from the samplesheet",
        )
        c.check(data.get("time_series") is not None, "aggregate report has the time-series panel")
    check_extras(
        c, out, samples, {e for e in args.extras.split(",") if e}, args.barcodes, args.reference
    )
    check_ids_intact(c, out)
    print(f"\n{c.total - c.failures}/{c.total} checks passed")
    return 0 if c.failures == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
