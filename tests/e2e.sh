#!/usr/bin/env bash
# End to end on the synthetic data (tests/data), all modes, with tests/check_outputs.py after each step:
#   full -> demix --previous_results (byte-identical demix/aggregate, previous dir untouched)
#        -> demix on a legacy layout built from the full run
#        -> aggregate --previous_results / --demix_dir
#        -> report --results_dir (HTML equals the original modulo timestamps)
#        -> probe-capture data (--protocol capture, duplicates marked vs not) and --skip_primer_trimming
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
run full --mode full --input "$SCRATCH/samplesheet.csv" --barcodes "$BARCODES" --primer_scheme midnight \
    --compare_barcodes container --freyja_boot --annotate_variants
check full --mode full --extras compare,boot,annotate || fail=1
checksums "$SCRATCH/out_full" > "$SCRATCH/out_full.sha256"

# 2. demix from the full results: byte-identical demix files and aggregate, previous results untouched
run demix --mode demix --previous_results "$SCRATCH/out_full" --barcodes "$BARCODES" --compare_barcodes container
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
run demix_legacy --mode demix --previous_results "$LEGACY" --barcodes "$BARCODES"
check demix_legacy --mode demix --no-amplicons --tolerance 0.05 || fail=1
same "freyja_aggregated.tsv (legacy demix)" "$SCRATCH/out_full/summary/freyja_aggregated.tsv" "$SCRATCH/out_demix_legacy/summary/freyja_aggregated.tsv"
checksums "$LEGACY" | diff - "$SCRATCH/legacy.sha256" && echo "legacy results unchanged"

# 4. aggregate from the full results (carries sample_type/collection_date from summary/samplesheet.csv)
run aggregate --mode aggregate --previous_results "$SCRATCH/out_full" --barcodes container
check aggregate --mode aggregate || fail=1
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
run capture --mode full --protocol capture --input "$SCRATCH/samplesheet_capture.csv" --barcodes "$BARCODES"
check capture --mode full --protocol capture || fail=1
run capture_nodup --mode full --protocol capture --duplicates none --input "$SCRATCH/samplesheet_capture.csv" --barcodes "$BARCODES"
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
run skip_trim --mode full --skip_primer_trimming --input "$SCRATCH/samplesheet_pure.csv" --barcodes "$BARCODES"
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

[ "$fail" -eq 0 ] && echo "E2E OK" || echo "E2E FAILED"
exit $fail
