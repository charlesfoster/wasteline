#!/usr/bin/env bash
# Demix-only equivalence (PLAN section 8): run `full` on tests/data, then `demix --previous_results`
# on its output with the same barcodes. Every <sample>.demix, freyja_aggregated.tsv and
# demix_status.tsv must be byte-identical, and the previous output directory must be untouched.
# The same is repeated on a legacy Snakemake layout derived from the full run's files.
#
# usage: tests/demix_equivalence.sh <scratch dir> [extra nextflow args, e.g. -profile docker,test]
# Runs from <scratch dir> (own .nextflow and work directories). Needs Docker.
set -euo pipefail

SCRATCH=$(mkdir -p "$1" && cd "$1" && pwd)
shift
PROFILE_ARGS=("$@")
[ ${#PROFILE_ARGS[@]} -gt 0 ] || PROFILE_ARGS=(-profile docker,test)
REPO=$(cd "$(dirname "$0")/.." && pwd)
BARCODES="$REPO/tests/data/barcodes_test.csv"

checksums() { (cd "$1" && find . -type f -print0 | sort -z | xargs -0 shasum -a 256); }
run() { # <name> <args...>: launch directory and work directory under the scratch dir
    local name=$1
    shift
    mkdir -p "$SCRATCH/launch_$name"
    (cd "$SCRATCH/launch_$name" && nextflow run "$REPO" "${PROFILE_ARGS[@]}" \
        --barcodes "$BARCODES" --outdir "$SCRATCH/out_$name" -work-dir "$SCRATCH/work_$name" "$@")
}

# samplesheet paths are relative to the repo root: write an absolute copy
sed "s#tests/data/#$REPO/tests/data/#g" "$REPO/tests/data/samplesheet.csv" > "$SCRATCH/samplesheet.csv"

run full --mode full --input "$SCRATCH/samplesheet.csv" --primer_scheme midnight
checksums "$SCRATCH/out_full" > "$SCRATCH/out_full.sha256.before"

fail=0
compare() { # <name> <outdir of the demix run>
    local out=$2
    for demix in "$SCRATCH"/out_full/samples/*/freyja/*.demix; do
        sample=$(basename "$demix" .demix)
        if cmp -s "$demix" "$out/samples/$sample/freyja/$sample.demix"; then
            echo "identical  $sample.demix ($1)"
        else
            echo "DIFFERENT  $sample.demix ($1)"
            fail=1
        fi
    done
    for f in freyja_aggregated.tsv demix_status.tsv; do
        if cmp -s "$SCRATCH/out_full/summary/$f" "$out/summary/$f"; then
            echo "identical  $f ($1)"
        else
            echo "DIFFERENT  $f ($1)"
            fail=1
        fi
    done
}

run demix --mode demix --previous_results "$SCRATCH/out_full"
compare new "$SCRATCH/out_demix"
checksums "$SCRATCH/out_full" | diff - "$SCRATCH/out_full.sha256.before" && echo "previous results unchanged"

# legacy layout: <s>/freyja/<s>.freyja_variants.tsv (one gzipped), <s>/freyja/<s>.depths (one gzipped)
LEGACY="$SCRATCH/legacy_results"
rm -rf "$LEGACY"
for dir in "$SCRATCH"/out_full/samples/*/; do
    s=$(basename "$dir")
    mkdir -p "$LEGACY/$s/freyja" "$LEGACY/$s/fastp"
    cp "$dir/freyja/$s.depths" "$LEGACY/$s/freyja/$s.depths"
    cp "$dir/fastp/$s.fastp.json" "$LEGACY/$s/fastp/${s}_fastp.json"
    cp "$dir/freyja/$s.variants.tsv" "$LEGACY/$s/freyja/$s.freyja_variants.tsv"
done
first=$(ls "$LEGACY" | head -n 1)
gzip -f "$LEGACY/$first/freyja/$first.freyja_variants.tsv"
mkdir -p "$LEGACY/$first/variants"
mv "$LEGACY/$first/freyja/$first.freyja_variants.tsv.gz" "$LEGACY/$first/variants/"
last=$(ls "$LEGACY" | tail -n 1)
gzip -f "$LEGACY/$last/freyja/$last.depths"
checksums "$LEGACY" > "$SCRATCH/legacy.sha256.before"

run demix_legacy --mode demix --previous_results "$LEGACY"
compare legacy "$SCRATCH/out_demix_legacy"
checksums "$LEGACY" | diff - "$SCRATCH/legacy.sha256.before" && echo "legacy results unchanged"

exit $fail
