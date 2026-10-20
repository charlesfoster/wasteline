#!/usr/bin/env bash
# End to end on the synthetic data (tests/data), all modes, with tests/check_outputs.py after each step:
#   full -> demix --previous_results (byte-identical demix/aggregate, previous dir untouched)
#        -> demix on a legacy layout built from the full run
#        -> aggregate --previous_results / --demix_dir
#        -> report --results_dir (HTML equals the original modulo timestamps)
#        -> probe-capture data (--protocol capture, duplicates marked vs not) and --skip_primer_trimming
#        -> lineage groups (on by default; the bundled list) and site aggregation, incl. a bad CSV that must abort early
#        -> the low-confidence band (--demix_floor_coverage, --include_low_confidence)
#        -> primer scheme check against a wrong declared scheme, CRAM output, time series + growth rates +
#           dashboard (aggregate mode, tests/data_timeseries), replay of the full run (bin/replay_run.py --execute)
#
# usage: tests/e2e.sh <scratch dir> [nextflow args, default: -profile docker]
# Runs every nextflow command from <scratch dir>/launch_<step> (own .nextflow), work dirs and outputs
# under <scratch dir>. Needs Docker and the container images; no network (test barcodes).
set -euo pipefail

SCRATCH=$(mkdir -p "$1" && cd "$1" && pwd)
shift
NF_ARGS=("$@")
[ ${#NF_ARGS[@]} -gt 0 ] || NF_ARGS=(-profile docker)
REPO=$(cd "$(dirname "$0")/.." && pwd)
BARCODES="$REPO/tests/data/barcodes_test.csv"
# check_outputs.py needs Python >= 3.12 (the pixi environment has it; the macOS system Python does not)
if [ -n "${PYTHON:-}" ]; then PY=("$PYTHON"); elif command -v pixi > /dev/null; then PY=(pixi run --manifest-path "$REPO/pixi.toml" python); else PY=(python3); fi

# Custom --barcodes files (tests/data/barcodes_test.csv) carry no lineage hierarchy, so the lineage-group summary is skipped
# with a warning there: steps that are not about it say so to keep their logs quiet (step 1 keeps the default and checks the skip)
LS_OFF=(--lineage_summary false)

checksums() { (cd "$1" && find . -type f -print0 | sort -z | xargs -0 shasum -a 256); }
run() { # <name> <args...>
    local name=$1
    shift
    mkdir -p "$SCRATCH/launch_$name"
    (cd "$SCRATCH/launch_$name" && nextflow run "$REPO" "${NF_ARGS[@]}" \
        --outdir "$SCRATCH/out_$name" -work-dir "$SCRATCH/work_$name" "$@")
}
check() { # <outdir name> <check_outputs.py args...>
    local name=$1
    shift
    echo "--- check_outputs $name"
    "${PY[@]}" "$REPO/tests/check_outputs.py" "$SCRATCH/out_$name" "$@"
}
same() { # <label> <file a> <file b>
    if cmp -s "$2" "$3"; then echo "identical  $1"; else echo "DIFFERENT  $1"; fail=1; fi
}
fail=0

sed "s#tests/data/#$REPO/tests/data/#g" "$REPO/tests/data/samplesheet.csv" > "$SCRATCH/samplesheet.csv"

# 1. full
# sites for the site-aggregation checks: SiteA = MIX_A + MIX_B, SiteB = PURE_C (+ the skipped LOWCOV_D), Site_1 = Site_1_week_2
printf 'aggregate_name,site_name\nBoth,SiteA\nBoth,SiteB\nAlone,Site_1\n' > "$SCRATCH/sites_full.csv"

run full --mode full --input "$SCRATCH/samplesheet.csv" --barcodes "$BARCODES" --primer_scheme midnight --site_aggregation "$SCRATCH/sites_full.csv" \
    --compare_barcodes container --freyja_boot --annotate_variants
check full --mode full --extras compare,boot,annotate --declared-scheme midnight --site-aggregation "$SCRATCH/sites_full.csv" || fail=1
checksums "$SCRATCH/out_full" > "$SCRATCH/out_full.sha256"

# 2. demix from the full results: byte-identical demix files and aggregate, previous results untouched
run demix --mode demix --previous_results "$SCRATCH/out_full" --barcodes "$BARCODES" --compare_barcodes container "${LS_OFF[@]}"
check demix --mode demix --no-amplicons --extras compare || fail=1
for demix in "$SCRATCH"/out_full/samples/*/freyja/*.demix; do # includes <s>.compare.demix (second barcodes)
    rel=${demix#"$SCRATCH/out_full/"}
    same "$(basename "$rel") (demix)" "$demix" "$SCRATCH/out_demix/$rel"
done
same "freyja_aggregated.tsv (demix)" "$SCRATCH/out_full/summary/freyja_aggregated.tsv" "$SCRATCH/out_demix/summary/freyja_aggregated.tsv"
same "demix_status.tsv (demix)" "$SCRATCH/out_full/summary/demix_status.tsv" "$SCRATCH/out_demix/summary/demix_status.tsv"
same "consensus_stats.tsv (demix)" "$SCRATCH/out_full/summary/consensus_stats.tsv" "$SCRATCH/out_demix/summary/consensus_stats.tsv"
checksums "$SCRATCH/out_full" | diff - "$SCRATCH/out_full.sha256" && echo "previous results unchanged"

# 3. demix from a legacy Snakemake layout (<s>/freyja/<s>.freyja_variants.tsv, one gzipped file of each kind
#    ones)
LEGACY="$SCRATCH/legacy_results"
rm -rf "$LEGACY"
for dir in "$SCRATCH"/out_full/samples/*/; do
    s=$(basename "$dir")
    mkdir -p "$LEGACY/$s/freyja" "$LEGACY/$s/fastp"
    cp "$dir/freyja/$s.depths" "$LEGACY/$s/freyja/$s.depths"
    cp "$dir/fastp/$s.fastp.json" "$LEGACY/$s/fastp/${s}_fastp.json"
    cp "$dir/freyja/$s.variants.tsv" "$LEGACY/$s/freyja/$s.freyja_variants.tsv"
    if [ -f "$dir/freyja/$s.demix" ]; then
        sed "1s/.*/	$s.freyja_variants.tsv/" "$dir/freyja/$s.demix" > "$LEGACY/$s/freyja/$s.demix"
    fi
done
gzip -f "$LEGACY/MIX_A/freyja/MIX_A.freyja_variants.tsv"
mkdir -p "$LEGACY/MIX_A/variants" && mv "$LEGACY/MIX_A/freyja/MIX_A.freyja_variants.tsv.gz" "$LEGACY/MIX_A/variants/"
gzip -f "$LEGACY/PURE_C/freyja/PURE_C.depths"
checksums "$LEGACY" > "$SCRATCH/legacy.sha256"
run demix_legacy --mode demix --previous_results "$LEGACY" --barcodes "$BARCODES" "${LS_OFF[@]}"
check demix_legacy --mode demix --no-amplicons --tolerance 0.05 || fail=1
same "freyja_aggregated.tsv (legacy demix)" "$SCRATCH/out_full/summary/freyja_aggregated.tsv" "$SCRATCH/out_demix_legacy/summary/freyja_aggregated.tsv"
checksums "$LEGACY" | diff - "$SCRATCH/legacy.sha256" && echo "legacy results unchanged"

# 3b. real ivar deletion (regression): ivar prints the deleted bases of a deletion as N (`-NNNNNNNNN`), so a VCF REF
#     built from them does not match the reference and `bcftools consensus` used to abort (all-N consensus). One legacy-layout
#     sample is PURE_C with the common 9-base deletion at 11288-11296 added to its variants (ivar's own format).
DELRES="$SCRATCH/deletion_results"
rm -rf "$DELRES"
mkdir -p "$DELRES/DEL_X/freyja"
cp "$SCRATCH/out_full/samples/PURE_C/freyja/PURE_C.depths" "$DELRES/DEL_X/freyja/DEL_X.depths"
"${PY[@]}" - "$SCRATCH/out_full/samples/PURE_C/freyja/PURE_C.variants.tsv" "$DELRES/DEL_X/freyja/DEL_X.freyja_variants.tsv" << 'PYMAKE'
import sys

lines = open(sys.argv[1]).read().splitlines()
row = "NC_045512.2\t11287\tG\t-NNNNNNNNN\t2\t0\t37\t174\t0\t35\t0.988636\t176\t0.01\tTRUE\tNA\tNA\tNA\tNA\tNA\tNA"
body = [line for line in lines[1:] if int(line.split("\t")[1]) != 11287]
body.append(row)
body.sort(key=lambda line: int(line.split("\t")[1]))
open(sys.argv[2], "w").write("\n".join([lines[0]] + body) + "\n")
PYMAKE
run demix_deletion --mode demix --previous_results "$DELRES" --barcodes "$BARCODES" "${LS_OFF[@]}"
"${PY[@]}" - "$SCRATCH/out_demix_deletion" << 'PYCHECK' || fail=1
import csv, json, sys

out = sys.argv[1]
seq = "".join(line.strip() for line in open(out + "/samples/DEL_X/consensus/DEL_X.consensus.fa") if not line.startswith(">"))
stats = {r["sample"]: r for r in csv.DictReader(open(out + "/summary/consensus_stats.tsv"), delimiter="\t")}["DEL_X"]
status = json.load(open(out + "/samples/DEL_X/consensus/DEL_X.consensus_status.json"))
ok = (
    stats["status"] == "ok" and status["status"] == "ok" and isinstance(status["warnings"], list)
    and seq[11286] == "G" and seq[11287:11296] == "-" * 9 and seq[11296] != "-" and len(seq) == 29903
)
print(("PASS" if ok else "FAIL") + "  real ivar deletion (-NNNNNNNNN at 11287): consensus ok, 9 '-' at 11288-11296, genome length kept")
sys.exit(0 if ok else 1)
PYCHECK

# 4. aggregate from the full results (carries sample_type/collection_date from summary/samplesheet.csv)
run aggregate --mode aggregate --previous_results "$SCRATCH/out_full" --barcodes container
check aggregate --mode aggregate --lineage-groups present || fail=1
same "freyja_aggregated.tsv (aggregate)" "$SCRATCH/out_full/summary/freyja_aggregated.tsv" "$SCRATCH/out_aggregate/summary/freyja_aggregated.tsv"
run aggregate_dir --mode aggregate --demix_dir "$SCRATCH/out_full/samples" --barcodes container
same "freyja_aggregated.tsv (aggregate --demix_dir)" "$SCRATCH/out_full/summary/freyja_aggregated.tsv" "$SCRATCH/out_aggregate_dir/summary/freyja_aggregated.tsv"

# 4b. aggregate of old results where a negative control has no .demix (it failed in the old pipeline): it is
#     listed as a control with status `missing`, not as an ordinary sample with a demix problem
SKIPPED_RESULTS="$SCRATCH/legacy_skipped_results"
rm -rf "$SKIPPED_RESULTS"
cp -R "$LEGACY" "$SKIPPED_RESULTS"
mkdir -p "$SKIPPED_RESULTS/NC_gone/freyja"
cp "$LEGACY/MIX_A/freyja/MIX_A.depths" "$SKIPPED_RESULTS/NC_gone/freyja/NC_gone.depths"
run aggregate_skipped --mode aggregate --previous_results "$SKIPPED_RESULTS" --barcodes container
"${PY[@]}" - "$SCRATCH/out_aggregate_skipped/summary" << 'PYCHECK' || fail=1
import json, sys

data = json.load(open(sys.argv[1] + "/report_data.json"))
controls = {c["sample"]: c for c in data["negative_controls"]}
flagged = [f for f in data["run_flags"] if f["code"] == "demix_missing" and "NC_gone" in f.get("samples", [])]
ok = (
    controls.get("NC_gone", {}).get("status") == "missing"
    and data["demix"]["NC_gone"]["status"] == "missing"
    and not flagged
)
print(("PASS" if ok else "FAIL") + "  skipped negative control NC_gone is listed as a control (missing), not flagged as a sample")
sys.exit(0 if ok else 1)
PYCHECK

# 5. report: rebuilt from the finished full run into a different outdir
run report --mode report --results_dir "$SCRATCH/out_full"
check report --mode report --original "$SCRATCH/out_full" || fail=1
checksums "$SCRATCH/out_full" | diff - "$SCRATCH/out_full.sha256" && echo "results of the report run unchanged"

# 6. probe-capture data (PLAN D13): no primer trimming, duplicates marked, no amplicon analysis
sed "s#tests/data_capture/#$REPO/tests/data_capture/#g" "$REPO/tests/data_capture/samplesheet.csv" > "$SCRATCH/samplesheet_capture.csv"
run capture --mode full --protocol capture --input "$SCRATCH/samplesheet_capture.csv" --barcodes "$BARCODES" "${LS_OFF[@]}"
check capture --mode full --protocol capture || fail=1
run capture_nodup --mode full --protocol capture --duplicates none --input "$SCRATCH/samplesheet_capture.csv" --barcodes "$BARCODES" "${LS_OFF[@]}"
check capture_nodup --mode full --protocol capture --duplicates none || fail=1
"${PY[@]}" - "$SCRATCH/out_capture/summary/qc_summary.tsv" "$SCRATCH/out_capture_nodup/summary/qc_summary.tsv" << 'PYCHECK' || fail=1
import csv, sys

def depth(path, sample):
    rows = {r["sample"]: r for r in csv.DictReader(open(path), delimiter="\t")}
    return float(rows[sample]["mean_depth"])

ok = all(depth(sys.argv[1], s) < depth(sys.argv[2], s) for s in ("CAP_MIX", "CAP_MIX3", "CAP_PURE"))
print(("PASS" if ok else "FAIL") + "  marked duplicates are excluded from the Freyja depths (mean depth drops with --duplicates mark)")
sys.exit(0 if ok else 1)
PYCHECK

# 6b. --skip_primer_trimming on amplicon reads that are already trimmed (one sample)
head -n 1 "$SCRATCH/samplesheet.csv" > "$SCRATCH/samplesheet_pure.csv"
grep '^PURE_C,' "$SCRATCH/samplesheet.csv" >> "$SCRATCH/samplesheet_pure.csv"
run skip_trim --mode full --skip_primer_trimming --input "$SCRATCH/samplesheet_pure.csv" --barcodes "$BARCODES" "${LS_OFF[@]}"
"${PY[@]}" - "$SCRATCH/out_skip_trim" << 'PYCHECK' || fail=1
import csv, glob, json, sys

out = sys.argv[1]
info = json.load(open(out + "/pipeline_info/run_info.json"))
trace = glob.glob(out + "/pipeline_info/execution_trace_*.txt")[0]
names = {r["process"].split(":")[-1] for r in csv.DictReader(open(trace), delimiter="\t")}
agg = open(out + "/summary/freyja_aggregated.tsv").read()
ok = (
    info["protocol"] == "amplicon" and info["primer_trimming"] is False and info["duplicates"] == "none"
    and not names & {"IVAR_TRIM", "AMPLICON_COVERAGE", "SAMTOOLS_MARKDUP"}
    and "JN.1" in agg
)
print(("PASS" if ok else "FAIL") + "  --skip_primer_trimming: no ivar trim, no amplicon QC, PURE_C still demixes to JN.1")
sys.exit(0 if ok else 1)
PYCHECK

# 7. a wrong declared scheme (articv4.1 for midnight reads): the primer scheme check reports a mismatch, names the
#    scheme that does fit, and the report carries the run flag
head -n 1 "$SCRATCH/samplesheet.csv" > "$SCRATCH/samplesheet_two.csv"
grep -E '^(MIX_A|PURE_C),' "$SCRATCH/samplesheet.csv" >> "$SCRATCH/samplesheet_two.csv"
run scheme_mismatch --mode full --primer_scheme articv4.1 --input "$SCRATCH/samplesheet_two.csv" --barcodes "$BARCODES" "${LS_OFF[@]}"
"${PY[@]}" - "$SCRATCH/out_scheme_mismatch/summary" << 'PYCHECK' || fail=1
import csv, json, sys

summary = sys.argv[1]
rows = {r["sample"]: r for r in csv.DictReader(open(summary + "/qc_summary.tsv"), delimiter="\t")}
data = json.load(open(summary + "/report_data.json"))
codes = {f["code"] for f in data["run_flags"]}
ok = (
    all(rows[s]["primer_scheme_check"] == "mismatch" for s in ("MIX_A", "PURE_C"))
    and all("midnight" in rows[s]["primer_scheme_best"] for s in ("MIX_A", "PURE_C"))
    and "primer_scheme_mismatch" in codes
    and data["samples"][0]["primer_scheme_check"] == "mismatch"
)
print(("PASS" if ok else "FAIL") + "  --primer_scheme articv4.1 on midnight reads: mismatch, best match midnight, run flag primer_scheme_mismatch")
sys.exit(0 if ok else 1)
PYCHECK

# 8. --bam_format cram: a CRAM (smaller than the BAM, indexed) is published instead of the BAM; the results are the same
run cram --mode full --bam_format cram --primer_scheme midnight --input "$SCRATCH/samplesheet_pure.csv" --barcodes "$BARCODES" "${LS_OFF[@]}"
"${PY[@]}" - "$SCRATCH" << 'PYCHECK' || fail=1
import os, sys

scratch = sys.argv[1]
cram_dir = scratch + "/out_cram/samples/PURE_C/bam/"
bam = scratch + "/out_full/samples/PURE_C/bam/PURE_C.primertrimmed.bam"
cram = cram_dir + "PURE_C.primertrimmed.cram"
same = open(scratch + "/out_cram/samples/PURE_C/freyja/PURE_C.demix", "rb").read() == open(scratch + "/out_full/samples/PURE_C/freyja/PURE_C.demix", "rb").read()
ok = (
    os.path.exists(cram) and os.path.exists(cram + ".crai")
    and not os.path.exists(cram_dir + "PURE_C.primertrimmed.bam")
    and os.path.getsize(cram) < os.path.getsize(bam)
    and same
)
print(("PASS" if ok else "FAIL") + "  --bam_format cram: indexed CRAM smaller than the BAM, no BAM published, same demix")
sys.exit(0 if ok else 1)
PYCHECK

# 9. time series, relative growth rates and the dashboard (aggregate mode on the dated synthetic data)
sed "s#tests/data_timeseries/#$REPO/tests/data_timeseries/#g" "$REPO/tests/data_timeseries/samplesheet.csv" > "$SCRATCH/samplesheet_timeseries.csv"
run timeseries --mode aggregate --input "$SCRATCH/samplesheet_timeseries.csv" --barcodes container --freyja_growth --freyja_dash \
    --site_aggregation "$REPO/tests/fixtures/site_aggregate/sites.csv"
"${PY[@]}" - "$SCRATCH/out_timeseries/summary" "$REPO/tests/data_timeseries/truth.json" << 'PYCHECK' || fail=1
import csv, json, os, sys

summary, truth_path = sys.argv[1], sys.argv[2]
truth = json.load(open(truth_path))
rows = list(csv.DictReader(open(summary + "/freyja_timeseries_lineages.tsv"), delimiter="\t"))
grouped = list(csv.DictReader(open(summary + "/freyja_timeseries.tsv"), delimiter="\t"))
first = truth["expected_site_date_abundance"][0]
got = {r["lineage_group"]: float(r["abundance"]) for r in rows if r["site"] == first["site"] and r["collection_date"] == first["collection_date"]}
site_date_ok = all(abs(got.get(lineage, -1) - value) < 0.02 for lineage, value in first["abundance"].items())
excluded = {"LOWCOV_1", "NEG_1", "NODATE_1"}
dates = {r["collection_date"] for r in rows}
growth = {}
with open(summary + "/freyja_growth_rates.csv") as handle:
    for row in csv.DictReader(handle):
        growth[row["Lineage"].split(" ")[0]] = float(row["Estimated Advantage"].rstrip("%"))  # Freyja writes e.g. 60.3%
rising_ok = growth.get(truth["lineage_rising"], 0) > 0 and growth.get(truth["lineage_declining"], 0) < 0
weighted_differs = any(
    r["abundance_load_weighted"] not in ("NA", "") and abs(float(r["abundance_load_weighted"]) - float(r["abundance"])) > 1e-6 for r in rows
)
checks = {
    "site/date abundance matches the truth (%s %s)" % (first["site"], first["collection_date"]): site_date_ok,
    "negative control, low-coverage and undated samples are excluded": not excluded & {r["site"] for r in rows} and "" not in dates and "NA" not in dates,
    "%d collection dates in the lineage-level table" % truth["n_dates"]: len(dates) == truth["n_dates"],
    "summarized-level table has the same site/dates": {(r["site"], r["collection_date"]) for r in grouped} == {(r["site"], r["collection_date"]) for r in rows},
    "growth rates: %s positive, %s negative" % (truth["lineage_rising"], truth["lineage_declining"]): rising_ok,
    "viral-load-weighted abundance differs from the plain mean": weighted_differs,
    "dashboard written": os.path.getsize(summary + "/freyja_dashboard.html") > 0,
    "report carries the time-series, growth and dashboard sections": all(
        json.load(open(summary + "/report_data.json")).get(k) for k in ("timeseries", "growth_rates", "dashboard")
    ),
}
bad = 0
for message, ok in checks.items():
    print(("PASS" if ok else "FAIL") + "  " + message)
    bad += 0 if ok else 1
sys.exit(1 if bad else 0)
PYCHECK
"${PY[@]}" - "$SCRATCH/out_timeseries/summary" "$REPO" "$REPO/tests/fixtures/site_aggregate/sites.csv" << 'PYCHECK' || fail=1
import csv, sys
from pathlib import Path

summary = sys.argv[1]
sys.path.insert(0, sys.argv[2] + "/tests")
import check_outputs

checker = check_outputs.Checker()
check_outputs.check_site_aggregates(checker, Path(summary).parent, sys.argv[3])
meta = {r["row_label"]: r for r in csv.DictReader(open(summary + "/site_aggregates_meta.tsv"), delimiter="\t")}
agg = {}
for line in open(summary + "/site_aggregates_freyja.tsv").read().splitlines()[1:]:
    cols = line.split("\t")
    agg[cols[0]] = dict(zip(cols[2].split(), map(float, cols[3].split())))
w3 = meta["Metro_2026-W03"]
lineage_rows = list(csv.DictReader(open(summary + "/lineage_summary.csv")))
samples = {r["sample"] for r in lineage_rows}
site_rows = list(csv.DictReader(open(summary + "/lineage_summary_site_aggregates.csv")))
checks = {
    "Metro_2026-W03: viral_load weights over 3 sites": w3["weight_basis"] == "viral_load" and w3["n_sites_included"] == "3",
    "Metro_2026-W03 XBB.1.5 = 0.107045 (hand-computed)": abs(agg["Metro_2026-W03"]["XBB.1.5"] - 0.107045) < 1e-6,
    "Metro_2026-W02: SiteC has no viral load, falls back to the population weights": meta["Metro_2026-W02"]["weight_basis"] == "site_weight (viral_load missing)",
    "negative control NEG_1 is not in lineage_summary.csv": "NEG_1" not in samples and "SiteA_W01_1" in samples,
    "lineage_summary_site_aggregates.csv has one row per aggregate-period": len(site_rows) == len(meta),
    "every aggregate lineage exists in at least one member site (check_outputs.check_site_aggregates)": checker.failures == 0,
}
bad = 0
for message, ok in checks.items():
    print(("PASS" if ok else "FAIL") + "  " + message)
    bad += 0 if ok else 1
sys.exit(1 if bad else 0)
PYCHECK
# too few dates for a growth fit: the table is written, the growth rates are a failure note, and no dashboard without --freyja_dash
awk -F, 'NR == 1 || $4 == "2026-01-05" || $4 == "2026-01-12"' "$SCRATCH/samplesheet_timeseries.csv" > "$SCRATCH/samplesheet_timeseries_short.csv"
run timeseries_short --mode aggregate --input "$SCRATCH/samplesheet_timeseries_short.csv" --barcodes container --freyja_growth
"${PY[@]}" - "$SCRATCH/out_timeseries_short/summary" << 'PYCHECK' || fail=1
import os, sys

summary = sys.argv[1]
note = open(summary + "/freyja_growth_rates.failed.txt").read() if os.path.exists(summary + "/freyja_growth_rates.failed.txt") else ""
ok = (
    os.path.exists(summary + "/freyja_timeseries.tsv")
    and note
    and not os.path.exists(summary + "/freyja_growth_rates.csv")
    and not os.path.exists(summary + "/freyja_dashboard.html")
)
print(("PASS" if ok else "FAIL") + "  two dates: time series written, growth rates replaced by a failure note, no dashboard unless requested")
sys.exit(0 if ok else 1)
PYCHECK

# 9b. lineage groups and site aggregation in demix mode (container barcodes: a dataset with a lineage hierarchy), with
#     a custom groups file that is published back byte-identical
cat > "$SCRATCH/groups_custom.csv" << 'GROUPS'
group,include,exclude,in_total
Omicron (all B.1.1.529*),B.1.1.529*,,yes
Delta (B.1.617.2*),B.1.617.2*,,yes
GROUPS
run demix_groups --mode demix --previous_results "$SCRATCH/out_full" --barcodes container \
    --lineage_groups "$SCRATCH/groups_custom.csv" --site_aggregation "$SCRATCH/sites_full.csv"
"${PY[@]}" - "$SCRATCH" << 'PYCHECK' || fail=1
import csv, json, sys

scratch = sys.argv[1]
summary = scratch + "/out_demix_groups/summary"
used = open(summary + "/lineage_groups_used.csv", "rb").read() == open(scratch + "/groups_custom.csv", "rb").read()
rows = {r["sample"]: r for r in csv.DictReader(open(summary + "/lineage_summary.csv"))}
info = json.load(open(scratch + "/out_demix_groups/pipeline_info/run_info.json"))
ungrouped = open(summary + "/lineage_ungrouped_summary.tsv").read().splitlines()[0]
ok = (
    used
    and set(rows) == {"LOWCOV_D", "MIX_A", "MIX_B", "PURE_C", "Site_1_week_2"}
    and rows["MIX_A"]["status"] == "ok" and rows["LOWCOV_D"]["status"] != "ok"
    and info["lineage_groups_source"] == "custom" and info["lineage_groups_file"] == "groups_custom.csv"
    and len(info["lineage_groups_sha256"]) == 64
    and ungrouped.startswith("lineage\tn_samples")
    and (open(summary + "/site_aggregates_meta.tsv").read().count("\n") > 1)
)
print(("PASS" if ok else "FAIL") + "  demix mode: custom groups file used and copied back, one row per sample, negatives absent, run_info records the file, site aggregates written")
sys.exit(0 if ok else 1)
PYCHECK

# 9c. mistakes in the user's CSVs stop the run within the first processes, naming the row/line (no mapping has run)
printf 'group,include,exclude,in_total\nJN.1*,JN.1*,,yes\nJN.1*,KP.3*,,yes\n' > "$SCRATCH/groups_bad.csv"
printf 'aggregate_name,site_name,weight\nMetro,SiteA,10\nMetro,,10\n' > "$SCRATCH/sites_bad.csv"
for bad in groups sites; do
    if [ "$bad" = groups ]; then flag=(--lineage_groups "$SCRATCH/groups_bad.csv"); else flag=(--site_aggregation "$SCRATCH/sites_bad.csv"); fi
    if run "bad_$bad" --mode full --input "$SCRATCH/samplesheet.csv" --barcodes "$BARCODES" --primer_scheme midnight "${flag[@]}" > "$SCRATCH/bad_$bad.log" 2>&1; then
        echo "FAIL  a bad $bad CSV did not stop the run"
        fail=1
    else
        "${PY[@]}" - "$SCRATCH/out_bad_$bad" "$SCRATCH/bad_$bad.log" "$bad" << 'PYCHECK' || fail=1
import csv, glob, sys

out, log, kind = sys.argv[1:4]
text = open(log).read()
trace = glob.glob(out + "/pipeline_info/execution_trace_*.txt")
names = {r["process"].split(":")[-1] for r in csv.DictReader(open(trace[0]), delimiter="\t")} if trace else set()
message = "appears twice (rows 2 and 3)" if kind == "groups" else "line 3: site_name is blank"
late = names & {"BWA_MEM", "SAMTOOLS_SORT", "FREYJA_VARIANTS", "FREYJA_DEMIX"}
ok = message in text and not late
print(("PASS" if ok else "FAIL") + f"  bad {kind} CSV: the run stops in the first processes with a message naming the row ({sorted(late) or 'no mapping step reached'})")
sys.exit(0 if ok else 1)
PYCHECK
    fi
done

# 9d. --report_anonymize full: the report (and report_data.json) carry no sample ID, site, path or user name; the private
#     id map is published and the other summary files keep the real IDs
run anonymize --mode full --report_anonymize full --primer_scheme midnight --input "$SCRATCH/samplesheet_two.csv" --barcodes "$BARCODES" "${LS_OFF[@]}"
"${PY[@]}" - "$SCRATCH" "$REPO" << 'PYCHECK' || fail=1
import csv, getpass, json, sys
from pathlib import Path

scratch, repo = Path(sys.argv[1]), sys.argv[2]
summary = scratch / "out_anonymize" / "summary"
html = (summary / "report.html").read_text()
data = (summary / "report_data.json").read_text()
secrets = ["MIX_A", "PURE_C", "SiteA", "SiteB", str(scratch), repo, getpass.getuser()]
leaks = [x for x in secrets if x and (x in html or x in data)]
id_map = list(csv.DictReader((summary / "sample_id_map.csv").open()))
originals = {r["original"] for r in id_map}
qc = (summary / "qc_summary.tsv").read_text()
ok = (
    not leaks
    and {"MIX_A", "PURE_C", "SiteA", "SiteB"} <= originals
    and json.loads(data).get("anonymized") == "full"
    and "MIX_A" in qc  # the other summary files keep the real IDs (documented)
)
print(("PASS" if ok else "FAIL") + f"  --report_anonymize full: no ID/site/path/user in report.html or report_data.json ({leaks or 'none found'}), id map published")
sys.exit(0 if ok else 1)
PYCHECK
# report mode rebuilds report_data.json (it exists in the results) when the report must be anonymised
run report_anonymized --mode report --results_dir "$SCRATCH/out_full" --report_anonymize paths
"${PY[@]}" - "$SCRATCH" "$REPO" << 'PYCHECK' || fail=1
import getpass, sys
from pathlib import Path

scratch, repo = Path(sys.argv[1]), sys.argv[2]
html = (scratch / "out_report_anonymized" / "summary" / "report.html").read_text()
leaks = [x for x in (str(scratch), repo, getpass.getuser()) if x and x in html]
ok = not leaks and "MIX_A" in html and not (scratch / "out_report_anonymized" / "summary" / "sample_id_map.csv").exists()
print(("PASS" if ok else "FAIL") + f"  mode report --report_anonymize paths: rebuilt without paths or user name ({leaks or 'none found'}), sample IDs kept, no id map")
sys.exit(0 if ok else 1)
PYCHECK

# 9e. low-confidence band (PLAN D18): with the gate at 98 % and the floor at 50 %, Site_1_week_2 (93.1 %) and NEG_CONTAM (94.5 %)
#     become low_confidence (demixed), MIX_A/MIX_B/PURE_C stay ok, LOWCOV_D/NEG_E stay skipped. Site_1_week_2 is in the
#     aggregate and the report (badge + flag) but not in the time series; --include_low_confidence true puts it there.
run low_confidence --mode demix --previous_results "$SCRATCH/out_full" --barcodes "$BARCODES" --primer_scheme midnight "${LS_OFF[@]}" \
    --min_genome_coverage_for_demix 98 --demix_floor_coverage 50
check low_confidence --mode demix --low-confidence Site_1_week_2,NEG_CONTAM || fail=1
run low_confidence_included --mode demix --previous_results "$SCRATCH/out_full" --barcodes "$BARCODES" --primer_scheme midnight "${LS_OFF[@]}" \
    --min_genome_coverage_for_demix 98 --demix_floor_coverage 50 --include_low_confidence true
check low_confidence_included --mode demix --low-confidence Site_1_week_2,NEG_CONTAM --include-low-confidence || fail=1

# 10. replay the full run: bin/replay_run.py reconstructs the command from run_info.json (the original is untouched);
#     the replayed demix files are byte-identical and the barcodes come from the recorded custom file
mkdir -p "$SCRATCH/launch_replay"
(cd "$SCRATCH/launch_replay" && python3 "$REPO/bin/replay_run.py" "$SCRATCH/out_full" --pipeline-dir "$REPO" \
    --outdir "$SCRATCH/out_replay" --params-out "$SCRATCH/replay_params.json" > "$SCRATCH/replay_plan.txt" 2>&1; echo "replay plan exit status: $?")
cat "$SCRATCH/replay_plan.txt"
grep -q "barcodes_test.csv" "$SCRATCH/replay_plan.txt" && echo "PASS  replay plan points at the recorded custom barcodes" || { echo "FAIL  replay plan does not name the barcodes file"; fail=1; }
(cd "$SCRATCH/launch_replay" && python3 "$REPO/bin/replay_run.py" "$SCRATCH/out_full" --pipeline-dir "$REPO" \
    --outdir "$SCRATCH/out_replay" --params-out "$SCRATCH/replay_params.json" --execute) || fail=1
for demix in "$SCRATCH"/out_full/samples/*/freyja/*.demix; do
    case "$demix" in *.compare.demix) continue ;; esac # the replay does not repeat --compare_barcodes
    rel=${demix#"$SCRATCH/out_full/"}
    same "$(basename "$rel") (replay)" "$demix" "$SCRATCH/out_replay/$rel"
done
checksums "$SCRATCH/out_full" | diff - "$SCRATCH/out_full.sha256" && echo "original results unchanged by the replay"

# 11. bin/reanalyse_runs.py: a dry run lists the full run and prints the command, writing nothing; --diff of two demixes of
#     the same variants with the same barcodes flags nothing; one real --execute with the bundled (container) barcodes, which
#     differ from the recorded test barcodes, makes a NEW directory with a diff CSV (the same barcodes would be `up to date`)
before=$(cd "$SCRATCH" && find . -maxdepth 3 -not -name "reanalyse_*" | sort | shasum -a 256)
python3 "$REPO/bin/reanalyse_runs.py" "$SCRATCH" --runs out_full --barcodes container --pipeline-dir "$REPO" > "$SCRATCH/reanalyse_dry.txt" 2>&1 || fail=1
after=$(cd "$SCRATCH" && find . -maxdepth 3 -not -name "reanalyse_*" | sort | shasum -a 256)
if grep -q "out_full" "$SCRATCH/reanalyse_dry.txt" && grep -q "nextflow run" "$SCRATCH/reanalyse_dry.txt" && grep -q -- "--mode demix" "$SCRATCH/reanalyse_dry.txt" && [ "$before" = "$after" ]; then
    echo "PASS  reanalyse_runs.py dry run: finds out_full, prints the demix command, writes nothing"
else
    echo "FAIL  reanalyse_runs.py dry run"; fail=1
fi
python3 "$REPO/bin/reanalyse_runs.py" --diff "$SCRATCH/out_full" "$SCRATCH/out_demix" --report "$SCRATCH/diff_same_barcodes.csv" > "$SCRATCH/reanalyse_diff.txt" 2>&1 || fail=1
flagged=$(awk -F, 'NR > 1 && $0 ~ /,yes(,|$)/' "$SCRATCH/diff_same_barcodes.csv" | wc -l | tr -d ' ')
[ "$flagged" -eq 0 ] && echo "PASS  reanalyse_runs.py --diff of the full run and its demix (same barcodes): 0 flagged" || { echo "FAIL  --diff flagged $flagged row(s)"; fail=1; }
mkdir -p "$SCRATCH/launch_reanalysis"
(cd "$SCRATCH/launch_reanalysis" && python3 "$REPO/bin/reanalyse_runs.py" "$SCRATCH" --runs out_full --barcodes container --pipeline-dir "$REPO" \
    --execute --outdir-suffix _reanalysis_test --extra-args "-work-dir $SCRATCH/work_reanalysis" --report "$SCRATCH/reanalysis_diff.csv") || fail=1
if [ -f "$SCRATCH/out_full_reanalysis_test/summary/reanalysis_diff.csv" ] && [ -f "$SCRATCH/out_full_reanalysis_test/summary/freyja_aggregated.tsv" ]; then
    echo "PASS  reanalyse_runs.py --execute: new directory with the re-analysed results and reanalysis_diff.csv"
else
    echo "FAIL  reanalyse_runs.py --execute produced no reanalysis_diff.csv"; fail=1
fi
checksums "$SCRATCH/out_full" | diff - "$SCRATCH/out_full.sha256" && echo "original results unchanged by the re-analysis"

[ "$fail" -eq 0 ] && echo "E2E OK" || echo "E2E FAILED"
exit $fail
