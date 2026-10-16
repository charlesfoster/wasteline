#!/usr/bin/env bash
# Compare the default aligner (bwa mem) with the experimental alternative (minibwa) on the synthetic data
# (tests/data): the pipeline runs twice, once with `--aligner bwa` and once with `--aligner minibwa`, then a table
# shows per sample the largest difference in lineage abundance (gate: <= TOLERANCE, default 0.02), the mapped read
# counts and the mean depths (printed, not gated). Exit status 1 if any abundance differs by more than TOLERANCE.
#
# usage: tests/compare_aligners.sh <scratch dir> [nextflow args, default: -profile docker]
#        TOLERANCE=0.02 tests/compare_aligners.sh /tmp/aligners -profile docker
# Needs Docker and the container images; no network (test barcodes). Requires the `--aligner` parameter
# (see docs/aligner_evaluation.md); reusing a scratch dir with -resume re-uses the finished runs.
set -euo pipefail

SCRATCH=$(mkdir -p "$1" && cd "$1" && pwd)
shift
NF_ARGS=("$@")
[ ${#NF_ARGS[@]} -gt 0 ] || NF_ARGS=(-profile docker)
REPO=$(cd "$(dirname "$0")/.." && pwd)
BARCODES="$REPO/tests/data/barcodes_test.csv"
TOLERANCE=${TOLERANCE:-0.02}

sed "s#tests/data/#$REPO/tests/data/#g" "$REPO/tests/data/samplesheet.csv" > "$SCRATCH/samplesheet.csv"

for aligner in bwa minibwa; do
    mkdir -p "$SCRATCH/launch_$aligner"
    (cd "$SCRATCH/launch_$aligner" && nextflow run "$REPO" "${NF_ARGS[@]}" \
        --outdir "$SCRATCH/out_$aligner" -work-dir "$SCRATCH/work_$aligner" \
        --mode full --input "$SCRATCH/samplesheet.csv" --barcodes "$BARCODES" --primer_scheme midnight \
        --aligner "$aligner")
done

# Host Python (>= 3.8, stdlib only)
python3 - "$SCRATCH/out_bwa" "$SCRATCH/out_minibwa" "$TOLERANCE" << 'PYCOMPARE'
import csv
import sys

bwa_dir, mini_dir, tolerance = sys.argv[1], sys.argv[2], float(sys.argv[3])


def read_tsv(path, key):
    with open(path, newline="") as handle:
        return {row[key]: row for row in csv.DictReader(handle, delimiter="\t")}


def abundances(row):
    """{lineage: abundance} from Freyja's space-separated `lineages` and `abundances` columns."""
    if row is None:
        return None
    return {name: float(value) for name, value in zip(row["lineages"].split(), row["abundances"].split())}


def number(row, column):
    try:
        return float(row[column])
    except (KeyError, TypeError, ValueError):
        return None


aggs = []
for path in [d + "/summary/freyja_aggregated.tsv" for d in (bwa_dir, mini_dir)]:
    with open(path, newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        index = reader.fieldnames[0]
        aggs.append({row[index]: row for row in reader})
qcs = [read_tsv(d + "/summary/qc_summary.tsv", "sample") for d in (bwa_dir, mini_dir)]

samples = sorted(set(qcs[0]) | set(qcs[1]))
header = ("sample", "max |d abundance|", "mapped bwa", "mapped minibwa", "d mapped %", "depth bwa", "depth minibwa", "d depth %", "status")
rows = []
failed = False
for sample in samples:
    a, b = abundances(aggs[0].get(sample)), abundances(aggs[1].get(sample))
    if a is None and b is None:
        delta, verdict = "NA", "not in aggregate (both)"
    elif a is None or b is None:
        delta, verdict = "NA", "ONLY IN " + ("minibwa" if a is None else "bwa")
        failed = True
    else:
        lineages = set(a) | set(b)
        delta_value = max(abs(a.get(name, 0.0) - b.get(name, 0.0)) for name in lineages)
        delta = "%.4f" % delta_value
        ok = delta_value <= tolerance
        failed = failed or not ok
        verdict = "ok" if ok else "DIFFERS (> %g)" % tolerance
    cells = [sample, delta]
    for column in ("mapped_reads", "mean_depth"):
        x, y = number(qcs[0].get(sample), column), number(qcs[1].get(sample), column)
        cells.append("NA" if x is None else "%g" % x)
        cells.append("NA" if y is None else "%g" % y)
        cells.append("NA" if x is None or y is None or x == 0 else "%+.2f" % (100.0 * (y - x) / x))
    cells.append(verdict)
    rows.append(cells)

# header/columns: sample, delta, mapped (bwa, minibwa, %), depth (bwa, minibwa, %), verdict
widths = [max(len(str(r[i])) for r in [header] + rows) for i in range(len(header))]
for r in [header] + rows:
    print("  ".join(str(c).ljust(widths[i]) for i, c in enumerate(r)))
print()
print(("FAIL" if failed else "PASS") + "  abundance differences between bwa and minibwa %s %g" % ("exceed" if failed else "are within", tolerance))
sys.exit(1 if failed else 0)
PYCOMPARE
