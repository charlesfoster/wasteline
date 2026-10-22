#!/usr/bin/env bash
# Minimal end to end on the synthetic data (tests/data): eight pipeline runs, about ten minutes, with tests/check_outputs.py
# after the ones that need it. Everything else (legacy layouts, CRAM, time-series fits, replay, ...) is covered by pytest.
#   1. full: all opt-ins on (compare barcodes, 20 bootstraps, variant annotation, site aggregation); the real end to end
#   2. demix --previous_results with the bundled (container) barcodes, a custom lineage groups file, site aggregation and
#      bootstraps: the demix files equal the full run's second-barcode demixes byte for byte, previous results untouched
#   3. demix --previous_results with the low-confidence band on (status low_confidence, in the aggregate, not in the time series)
#   4. aggregate --previous_results (aggregate identical to the full run's), and again from step 3's results (low_confidence carried over)
#   5. report --results_dir rebuilt with --report_anonymize paths (no path or user name in the report)
#   6. mixed libraries in one samplesheet: probe-capture samples (--protocol capture, the run default) and two amplicon samples
#      (per-sample protocol / primer_scheme cells)
#   + two bad CSVs (lineage groups, sites) that must abort the run early (cheap: nothing is mapped)
#   + bin/reanalyse_runs.py dry run and --diff (no pipeline run)
#   9. a reference with two contigs (SARS-CoV-2 + measles N, PLAN D25): Freyja demixes the SARS-CoV-2 contig, both get coverage QC and
#      a consensus; a simulated sample covers both, a step 6 capture sample (reads of SARS-CoV-2 only) leaves the measles contig
#      uncovered; its demix equals the single-contig run's byte for byte
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
    --compare_barcodes container --freyja_boot --freyja_boot_n 20 --annotate_variants
check full --mode full --extras compare,boot,annotate --declared-scheme midnight --site-aggregation "$SCRATCH/sites_full.csv" || fail=1
checksums "$SCRATCH/out_full" > "$SCRATCH/out_full.sha256"

unchanged() { # <label>: the full run's results are byte-for-byte what they were
    if checksums "$SCRATCH/out_full" | diff - "$SCRATCH/out_full.sha256" > /dev/null; then echo "identical  previous results ($1)"; else echo "DIFFERENT  previous results changed ($1)"; fail=1; fi
}

# 2. demix from the full results with the bundled (container) barcodes, a custom lineage groups file and site aggregation
#    (a dataset with a lineage hierarchy; the custom groups file is published back byte-identical): the demix files equal
#    the second-barcode demixes (<s>.compare.demix) of the full run byte for byte
cat > "$SCRATCH/groups_custom.csv" << 'GROUPS'
group,include,exclude,in_total
Omicron (all B.1.1.529*),B.1.1.529*,,yes
Delta (B.1.617.2*),B.1.617.2*,,yes
GROUPS
run demix_groups --mode demix --previous_results "$SCRATCH/out_full" --barcodes container \
    --lineage_groups "$SCRATCH/groups_custom.csv" --site_aggregation "$SCRATCH/sites_full.csv" --freyja_boot --freyja_boot_n 5
for demix in "$SCRATCH"/out_full/samples/*/freyja/*.compare.demix; do
    s=$(basename "$demix" .compare.demix)
    same "$s.demix (demix from previous results, container barcodes)" "$demix" "$SCRATCH/out_demix_groups/samples/$s/freyja/$s.demix"
done
unchanged "demix"
"${PY[@]}" - "$SCRATCH" << 'PYCHECK' || fail=1
import csv, json, sys

scratch = sys.argv[1]
summary = scratch + "/out_demix_groups/summary"
used = open(summary + "/lineage_groups_used.csv", "rb").read() == open(scratch + "/groups_custom.csv", "rb").read()
rows = {r["sample"]: r for r in csv.DictReader(open(summary + "/lineage_summary.csv"))}
info = json.load(open(scratch + "/out_demix_groups/pipeline_info/run_info.json"))
ungrouped = open(summary + "/lineage_ungrouped_summary.tsv").read().splitlines()[0]
# --freyja_boot: lineage_group_intervals.tsv has one row per sample x (group + other) and its point estimates are the summary values
intervals_text = open(summary + "/lineage_group_intervals.tsv").read()
intervals = [r for r in csv.DictReader(intervals_text.splitlines(), delimiter="\t") if r["sample"] in rows]
by_sample = {}
for r in intervals:
    by_sample.setdefault(r["sample"], {})[r["group"]] = r
intervals_ok = (
    intervals_text.splitlines()[0] == "sample\tgroup\tabundance\tlo\thi\tn_boot"
    and set(by_sample) == {"MIX_A", "MIX_B", "PURE_C", "Site_1_week_2"}
    and all(list(g) == ["Omicron (all B.1.1.529*)", "Delta (B.1.617.2*)", "other"] for g in by_sample.values())
    and all(
        abs(float(g[name]["abundance"]) - float(rows[sample][name] or 0)) < 1e-6
        and float(g[name]["lo"]) <= float(g[name]["hi"]) and g[name]["n_boot"] == "5"
        for sample, g in by_sample.items() for name in g
    )
)
ok = (
    used
    and intervals_ok
    and set(rows) == {"LOWCOV_D", "MIX_A", "MIX_B", "PURE_C", "Site_1_week_2"}
    and rows["MIX_A"]["status"] == "ok" and rows["LOWCOV_D"]["status"] != "ok"
    and info["lineage_groups_source"] == "custom" and info["lineage_groups_file"] == "groups_custom.csv"
    and len(info["lineage_groups_sha256"]) == 64
    and ungrouped.startswith("lineage\tn_samples")
    and (open(summary + "/site_aggregates_meta.tsv").read().count("\n") > 1)
)
print(("PASS" if ok else "FAIL") + "  demix mode: custom groups file used and copied back, one row per sample, negatives absent, run_info records the file, site aggregates written, group intervals (--freyja_boot) match the summary")
sys.exit(0 if ok else 1)
PYCHECK

# 3. low-confidence band (PLAN D18): with the gate at 98 % and the floor at 50 %, Site_1_week_2 (93.1 %) and NEG_CONTAM (94.5 %)
#    become low_confidence (demixed), MIX_A/MIX_B/PURE_C stay ok, LOWCOV_D/NEG_E stay skipped. Site_1_week_2 is in the
#    aggregate and the report (badge + flag) but not in the time series. Same (test) barcodes as the full run, so the demix
#    files and the aggregate are byte-identical to it, and reanalyse_runs.py --diff flags nothing.
run low_confidence --mode demix --previous_results "$SCRATCH/out_full" --barcodes "$BARCODES" --primer_scheme midnight "${LS_OFF[@]}" \
    --min_genome_coverage_for_demix 98 --demix_floor_coverage 50
check low_confidence --mode demix --low-confidence Site_1_week_2,NEG_CONTAM || fail=1
for demix in "$SCRATCH"/out_full/samples/*/freyja/*.demix; do
    case "$demix" in *.compare.demix) continue ;; esac
    rel=${demix#"$SCRATCH/out_full/"}
    same "$(basename "$rel") (demix, same barcodes)" "$demix" "$SCRATCH/out_low_confidence/$rel"
done
same "freyja_aggregated.tsv (demix, same barcodes)" "$SCRATCH/out_full/summary/freyja_aggregated.tsv" "$SCRATCH/out_low_confidence/summary/freyja_aggregated.tsv"
same "consensus_stats.tsv (demix, same barcodes)" "$SCRATCH/out_full/summary/consensus_stats.tsv" "$SCRATCH/out_low_confidence/summary/consensus_stats.tsv"
unchanged "low_confidence"

# 4. aggregate from the full results (carries sample_type/collection_date from summary/samplesheet.csv)
run aggregate --mode aggregate --previous_results "$SCRATCH/out_full" --barcodes container
check aggregate --mode aggregate --lineage-groups present || fail=1
same "freyja_aggregated.tsv (aggregate)" "$SCRATCH/out_full/summary/freyja_aggregated.tsv" "$SCRATCH/out_aggregate/summary/freyja_aggregated.tsv"
unchanged "aggregate"

# 4b. aggregate from the low-confidence demix results (step 3): the low_confidence status is carried over from summary/demix_status.tsv
#     (PLAN D18) into the summaries and the report, and the sample is left out of the time series like in the demix run
run aggregate_low_confidence --mode aggregate --previous_results "$SCRATCH/out_low_confidence" --barcodes container
"${PY[@]}" - "$SCRATCH/out_aggregate_low_confidence/summary" << 'PYCHECK' || fail=1
import csv, json, sys

summary = sys.argv[1]
status = {r["sample"]: r for r in csv.DictReader(open(summary + "/demix_status.tsv"), delimiter="\t")}
lineage = {r["sample"]: r for r in csv.DictReader(open(summary + "/lineage_summary.csv"))}
series_sites = {r["site"] for r in csv.DictReader(open(summary + "/freyja_timeseries.tsv"), delimiter="\t")}
report = json.load(open(summary + "/report_data.json"))
flag = next((f for f in report["run_flags"] if f["code"] == "demix_low_confidence"), None)
low = status["Site_1_week_2"]
ok = (
    low["status"] == "low_confidence" and "high-confidence threshold" in low["reason"] and low["genome_coverage_pct"] != ""
    and lineage["Site_1_week_2"]["status"] == "low_confidence" and lineage["MIX_A"]["status"] == "ok"
    and "Site_1" not in series_sites and "SiteA" in series_sites
    and flag is not None and flag.get("samples") == ["Site_1_week_2"]
)
print(("PASS" if ok else "FAIL") + "  aggregate mode: low_confidence carried over from the previous demix_status.tsv (summaries, report flag), left out of the time series")
sys.exit(0 if ok else 1)
PYCHECK
unchanged "aggregate_low_confidence"

# 5. report mode, rebuilt from the finished full run into a different outdir with --report_anonymize paths (report mode
#    rebuilds report_data.json, which the results already have, when the report must be anonymised)
run report_anonymized --mode report --results_dir "$SCRATCH/out_full" --report_anonymize paths
check report_anonymized --mode report || fail=1
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
unchanged "report"

# 6. mixed libraries in one run (PLAN D13 + D20): the five probe-capture samples take the run default (--protocol capture: no
#    primer trimming, duplicates marked, no amplicon analysis); MIX_A and PURE_C say protocol=amplicon and primer_scheme=midnight
#    in their own cells (primer-trimmed with that scheme, amplicon analyses, no duplicate marking)
{
    echo "sample,fastq_1,fastq_2,sample_type,collection_date,site,protocol,primer_scheme"
    tail -n +2 "$REPO/tests/data_capture/samplesheet.csv" | sed "s#tests/data_capture/#$REPO/tests/data_capture/#g; s#\$#,,#"
    grep -E '^(MIX_A|PURE_C),' "$REPO/tests/data/samplesheet.csv" | sed "s#tests/data/#$REPO/tests/data/#g; s#\$#,amplicon,midnight#"
} > "$SCRATCH/samplesheet_mixed.csv"
run mixed --mode full --protocol capture --input "$SCRATCH/samplesheet_mixed.csv" --barcodes "$BARCODES" "${LS_OFF[@]}"
check mixed --mode full --protocol mixed --amplicon-samples MIX_A,PURE_C || fail=1

# 7. mistakes in the user's CSVs stop the run within the first processes, naming the row/line (no mapping has run)
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

# 8. bin/reanalyse_runs.py (no pipeline run): a dry run lists the full run and prints the command, writing nothing; --diff of
#    the full run and the low-confidence demix (same barcodes) flags nothing
before=$(cd "$SCRATCH" && find . -maxdepth 3 -not -name "reanalyse_*" | sort | shasum -a 256)
python3 "$REPO/bin/reanalyse_runs.py" "$SCRATCH" --runs out_full --barcodes container --pipeline-dir "$REPO" > "$SCRATCH/reanalyse_dry.txt" 2>&1 || fail=1
after=$(cd "$SCRATCH" && find . -maxdepth 3 -not -name "reanalyse_*" | sort | shasum -a 256)
if grep -q "out_full" "$SCRATCH/reanalyse_dry.txt" && grep -q "nextflow run" "$SCRATCH/reanalyse_dry.txt" && grep -q -- "--mode demix" "$SCRATCH/reanalyse_dry.txt" && [ "$before" = "$after" ]; then
    echo "PASS  reanalyse_runs.py dry run: finds out_full, prints the demix command, writes nothing"
else
    echo "FAIL  reanalyse_runs.py dry run"; fail=1
fi
python3 "$REPO/bin/reanalyse_runs.py" --diff "$SCRATCH/out_full" "$SCRATCH/out_low_confidence" --report "$SCRATCH/diff_same_barcodes.csv" > "$SCRATCH/reanalyse_diff.txt" 2>&1 || fail=1
flagged=$(awk -F, 'NR > 1 && $0 ~ /,yes(,|$)/' "$SCRATCH/diff_same_barcodes.csv" | wc -l | tr -d ' ')
[ "$flagged" -eq 0 ] && echo "PASS  reanalyse_runs.py --diff of the full run and its demix (same barcodes): 0 flagged" || { echo "FAIL  --diff flagged $flagged row(s)"; fail=1; }

# 9. several reference contigs (PLAN D25): multi.fa = SARS-CoV-2 + measles N, --freyja_contig NC_045512.2. S1 is simulated (B.1.1.7
#    on the SARS-CoV-2 contig, reads from both contigs by length); CAP_MIX is the step 6 capture pair, whose demix must equal step 6's
#    byte for byte (the same reads, the same barcodes: the extra contig only adds a contig that nothing maps to)
mkdir -p "$SCRATCH/multi"
cat "$REPO/assets/references/NC_045512.fasta" "$REPO/assets/references/NC_001498.1.fasta" > "$SCRATCH/multi/multi.fa"
"${PY[@]}" "$REPO/tests/simulate_pathogen_mixture.py" --reference "$SCRATCH/multi/multi.fa" --freyja-contig NC_045512.2 --barcodes "$BARCODES" \
    --lineage B.1.1.7 --proportion 1 --depth 300 --out "$SCRATCH/multi/S1" > /dev/null
{
    echo "sample,fastq_1,fastq_2,sample_type,collection_date,site"
    echo "S1,$SCRATCH/multi/S1_R1.fastq.gz,$SCRATCH/multi/S1_R2.fastq.gz,sample,2026-04-06,SiteE"
    grep '^CAP_MIX,' "$REPO/tests/data_capture/samplesheet.csv" | sed "s#tests/data_capture/#$REPO/tests/data_capture/#g"
} > "$SCRATCH/samplesheet_multicontig.csv"
run multicontig --mode full --protocol capture --input "$SCRATCH/samplesheet_multicontig.csv" --barcodes "$BARCODES" "${LS_OFF[@]}" \
    --pathogen SARS-CoV-2 --fasta "$SCRATCH/multi/multi.fa" --freyja_contig NC_045512.2
check multicontig --mode full --multicontig || fail=1
same "CAP_MIX.demix (two-contig reference vs the single-contig capture run)" "$SCRATCH/out_mixed/samples/CAP_MIX/freyja/CAP_MIX.demix" "$SCRATCH/out_multicontig/samples/CAP_MIX/freyja/CAP_MIX.demix"

[ "$fail" -eq 0 ] && echo "E2E OK" || echo "E2E FAILED"
exit $fail
