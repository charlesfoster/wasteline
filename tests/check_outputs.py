#!/usr/bin/env python3
"""Tolerance-based end-to-end checks of a wasteline output directory against tests/data/truth.json.

  tests/check_outputs.py <outdir> --mode full|demix|aggregate|report [--truth tests/data/truth.json]
                         [--barcodes tests/data/barcodes_test.csv] [--reference assets/references/NC_045512.fasta]
                         [--original <outdir of the run the report was rebuilt from>]

Standard library only. Prints one line per check and exits 1 if any check failed. What is checked depends on
the mode (a `demix` run has no BAM, hence no flagstat or mapped-read counts and, without a primer scheme, no
amplicon dropout; an `aggregate` run has no QC or consensus):

  demix (full, demix)   demix_status.tsv, every <sample>.demix_status.json and barcode_version.json agree with
                        truth and carry the barcode version
  aggregate             rows = expected samples only (negatives and skipped samples absent); true lineages
                        within +-tolerance in the `lineages`/`abundances` rows (not `summarized`); no lineage
                        absent from the truth above 0.05
  qc (full, demix)      qc_status / qc_reasons per sample
  consensus             one consensus per non-negative sample with the right status; JN.1 SNVs in PURE_C
  report                report_data.json header and run flags, report.html is self-contained, sample IDs intact
  report (mode report)  the rebuilt report.html equals --original's modulo timestamps
"""

import argparse
import csv
import json
import re
import sys
from pathlib import Path

ABSENT_LINEAGE_MAX = 0.05
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


def check_aggregate(c, out, samples, tolerance_override=None):
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
        tol = tolerance_override or info.get("abundance_tolerance") or 0.05
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
        c.check(
            (out / "samples" / sample / "consensus" / f"{sample}.consensus.fa").exists(),
            f"{sample}: samples/{sample}/consensus/{sample}.consensus.fa published",
        )
    if "PURE_C" in fasta and barcodes and Path(barcodes).exists():
        snvs = barcode_snvs(barcodes, "JN.1")
        seq = fasta["PURE_C"]
        hit = sum(
            1
            for pos, alt in snvs
            if pos <= len(seq) and (seq[pos - 1] == alt or alt in IUPAC_FOR.get(seq[pos - 1], ""))
        )
        c.check(
            snvs and hit / len(snvs) >= SNV_FRACTION_MIN,
            f"PURE_C consensus carries the JN.1 SNVs ({hit}/{len(snvs)})",
        )
    if reference and Path(reference).exists() and "PURE_C" in fasta:
        ref = "".join(
            line for line in Path(reference).read_text().splitlines() if not line.startswith(">")
        )
        c.check(len(fasta["PURE_C"]) == len(ref), "PURE_C consensus has the reference length")


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
    c.check(
        "Site_1_week_2" in html and "Site-1-week-2" not in html,
        "report.html: sample ID Site_1_week_2 intact",
    )
    codes = {f.get("code"): f for f in data.get("run_flags", [])}
    if mode in ("full", "demix"):
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


def normalise_html(text):
    text = re.sub(
        r"\d{4}-\d{2}-\d{2}[T ]\d{2}:\d{2}:\d{2}(?:\.\d+)?(?:Z|[+-]\d{2}:\d{2})?", "<TS>", text
    )
    return re.sub(r"[ \t]+", " ", text)


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


def main(argv=None):
    here = Path(__file__).resolve().parent
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("outdir")
    parser.add_argument("--mode", required=True, choices=["full", "demix", "aggregate", "report"])
    parser.add_argument("--truth", default=str(here / "data" / "truth.json"))
    parser.add_argument("--barcodes", default=str(here / "data" / "barcodes_test.csv"))
    parser.add_argument(
        "--reference", default=str(here.parent / "assets" / "references" / "NC_045512.fasta")
    )
    parser.add_argument(
        "--original", help="mode report: output directory of the run whose results were re-reported"
    )
    parser.add_argument(
        "--extras",
        default="",
        help="comma list of opt-in features that were on: compare,boot,annotate",
    )
    parser.add_argument(
        "--no-amplicons",
        action="store_true",
        help="the run had no primer scheme (demix without --primer_scheme)",
    )
    parser.add_argument(
        "--tolerance",
        type=float,
        default=None,
        help="override the abundance tolerance (e.g. when real barcodes are used)",
    )
    args = parser.parse_args(argv)

    out = Path(args.outdir)
    truth = json.loads(Path(args.truth).read_text())
    samples = truth["samples"]
    c = Checker()

    if args.mode == "report":
        # report mode writes only summary/report.html (+ report_data.json when it had to rebuild it)
        c.check((out / "summary" / "report.html").exists(), "summary/report.html exists")
        if args.original:
            new = (out / "summary" / "report.html").read_text()
            old = (Path(args.original) / "summary" / "report.html").read_text()
            c.check(
                normalise_html(new) == normalise_html(old),
                "rebuilt report.html equals the original modulo timestamps",
            )
            c.check(
                "Site_1_week_2" in new and "https://" not in new and "http://" not in new,
                "rebuilt report.html: IDs intact, no URLs",
            )
        print(f"\n{c.total - c.failures}/{c.total} checks passed")
        return 0 if c.failures == 0 else 1

    check_run_info(c, out, args.mode)
    if args.mode in ("full", "demix"):
        have_amplicons = args.mode == "full" or not args.no_amplicons
        check_demix_status(c, out, samples, args.mode)
        check_aggregate(c, out, samples, args.tolerance)
        check_qc(c, out, samples, args.mode, have_amplicons)
        check_consensus(c, out, samples, args.barcodes, args.reference)
        check_published_qc_is_complete(c, out, samples)
    else:
        check_aggregate(c, out, samples, args.tolerance)
        c.check(
            not (out / "summary" / "barcode_version.json").exists(),
            "aggregate mode writes no barcode_version.json sidecar",
        )
    check_report(c, out, samples, args.mode)
    check_negative_controls(c, out, samples)
    if args.mode == "aggregate":
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
