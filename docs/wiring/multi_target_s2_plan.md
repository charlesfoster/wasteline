# S2 plan: panels and per-target analysis (decided 2026-10-23; extends multi_target_in_pipeline_plan.md)

Branch `feature/multi-target-in-pipeline`. Work units U0-U6 below; one owner per file. Test limits are hard maxima (docs/development.md).

DECIDED: in-pipeline design; one competitive mapping of each sample to the concatenated panel, then split the BAM per target; SARS byte-identity is
NOT a constraint on this branch (numerics must be right; dev is the safety net); only Freyja-official pathogens (registry = the 10 with barcodes), no
barcode-less targets (enterovirus/HCV dropped); thresholds (params, uncalibrated) min fragments 100 / breadth 5 % / family share 5 % / MAPQ 10;
statuses detected | trace | not_detected; trace is listed in targets.tsv (analysed=false) but not analysed; panel runs write samples/<id>/<target>/...,
summary/<target>/..., summary/targets.tsv, summary/overview.html + per-target report; single-target runs keep the legacy layout; no default-layout flip.
Flu panel: whole-genome H1N1/H3N2/H5N1 (family fluA, marker/freyja contig = HA) + B/Vic HA-only (fluB); capture default. ROADMAP.md added (ONT, more Freyja pathogens).
KEY FINDINGS: no dynamic process fan-out -> every run-level process needs a leading `tuple val(meta)` (about 19 modules + nf-test snapshot updates);
unsafe sites are .collect()/collectFile/toSortedList in workflows/wasteline.nf (ch_failed_stats consensus_stats names, times.csv, summary_samplesheet,
skipped_samplesheet, demix_status.tsv, barcode_version.json, barcode_compare*, freyja_boot_summary, contig_coverage, annotated_variants, every .collect() into *_CONCAT)
-> groupTuple on tkey(meta)=[id:target.name,target]; modules.config uses params.pathogen/freyja_contig at ~13 lines (118,122,167,176,198,211,250,259,341,375,386,434,449)
-> meta.target.*; nextflow.config ~l.184 loads a pathogen profile: with --panel and no --pathogen load `generic`; samtools container has no Python
-> typing = awk in a samtools module (PANEL_ASSIGN) + python in PANEL_CALL; container pins via bin/lock_containers.py (needs network);
primerBedFor/sourceBedFor/reconcileBed must default to meta.target (S1 left them on targetFromParams()); meta.key = id (legacy) or id__target (panel).
CONTRACTS: target map adds marker_contig, family, primer_scheme, barcodes, panel(bool), sample_sub/summary_sub ('' legacy, '/<name>' panel); per-sample
publishDir samples/${meta.id}${meta.target.sample_sub}/<step>, run-level summary${meta.target.summary_sub}; pre-split samples carry pseudo target 'panel'
(fastp + typing publish under samples/<id>/); per-target fasta/gff/dataset via multiMap/combine(by: tkey) not global value channels; new params panel, targets
(filter for demix/aggregate/report), panel_min_fragments/breadth_pct/family_share/mapq; --panel excludes --pathogen/--fasta/--gff; samplesheet primer_scheme cell
rejected with >1 target (panel CSV primer_scheme column instead); contig names unique across the panel; a registry row's barcodes cell beats --barcodes (keeps e2e offline).
REGISTRY: assets/panels/registry.csv (target,fasta,gff,freyja_pathogen,freyja_contig,marker_contig,family,primer_scheme,primer_scheme_default,primer_schemes_bundled,
lineage_groups_bundled,freyja_container_pathogen,barcodes) + rsv.csv + flu.csv; pytest checks registry vs conf/pathogens/*.config. H5Nx and H5Nx-cattle share a FASTA (one per panel).
TYPING: bin/panel_select.py `prepare` (panel.fasta, contig_map.tsv, panel_resolved.csv) and `call` (stats+coverage -> <id>.targets.tsv). Assign: collate/view -F 0x900, fragment
belongs to T if all mapped primaries are on T's contigs with MAPQ>=10, else ambiguous. Stats on marker contig(s) only (flu internal segments cross-map). Rank by fragments,
aligned bases, breadth, name; detected = frags>=min AND breadth>=min AND (rank1 OR share>=min); trace = >=1 frag but not detected (reason names the failed floor); not_detected = 0;
mixed = >=2 detected in one family; sample qc_fail `no_panel_reads` still gets a row per target. targets.tsv columns: sample,sample_type,target,family,status,reason,
sample_status,mixed,secondary_targets,fragments,ambiguous_fragments,fraction_of_family,breadth_pct,mean_depth,mean_mapq,second_target,margin,rescue,analysed. Negative controls are
analysed for every target detected in >=1 non-negative sample. PANEL_SPLIT: samtools view -N qnames per target contigs, header trimmed, emits [meta_t,bam,bai,count] like BAM_MAPPED_READS.out.counted.
WORK UNITS (order U0+U1 parallel -> U2 -> U3 -> U4 -> U5 -> U6; one owner per file; modules.config only by U3 then U4):
 U0 foundations: utils main.nf (targetsFromParams, targetByName, panelTarget, tkey, resultsAreMultiTarget, sampleMeta key, panel validation), nextflow.config, main.nf typed params,
    nextflow_schema.json, assets/schema_input*.json (optional target col), assets/panels/*, empty conf/modules_panel.config.
 U1 typing: bin/panel_select.py, tests/python/test_panel_select.py, modules/local/panel_{prepare,assign,call,split} (+environment.yml, stub nf-test), conf/containers.lock.config.
 U2 subworkflows/local/panel_front: prepare, index, bwa mem on panel, assign, call, barrier, negatives expansion, split, targets.tsv.
 U3 (slice S2a) workflows/wasteline.nf + conf/modules.config + freyja_data (takes channel of targets) + freyja_update/dataset_info: per-sample chain per target, run-level block extracted
    verbatim to subworkflows/local/target_summaries (gated off for panels with a warning), per-target QC/consensus concats + demix_status.
 U4 (S2b) target_summaries + run-level modules with leading meta (S2b-1: aggregate/plot/controls/filter/timeseries/provenance/build_report_data/render_report; S2b-2 cuttable: growth, dash,
    site aggregate, lineage summary, boot intervals, dropout, unexplained, moi). Report inputs via one helper joinOptional(); runInfoBase per target -> pipeline_info/<t>/run_info.json.
 U5 (S2c) modes: bin/make_samplesheet.py target column + discovery of samples/<id>/<t>/freyja and summary/<t>/samplesheet.csv; discover_results skipped TSV target col;
    results_input resultsRow resolves target, dup check (id,target), --targets filter; report mode loops targets; carriedLowConfidence per target; reanalyse_runs.py + runs_store.py treat each target dir as a run.
 U6 (S2d) reports/docs: build_report_data.py --target, report.html.j2 back-link, bin/build_overview.py + modules/local/build_overview (summary/overview.html: sample x target matrix), docs/other_pathogens.md, CHANGELOG.
E2E (one step 10): tiny RSVa+RSVb sims (tests/simulate_pathogen_mixture.py, ~1500 pairs/sample, 2-lineage barcodes CSVs, --protocol capture, no bootstraps): S1 pure RSVa, S2 pure RSVb,
 S3 70/30 (mixed=true, fraction within 0.05), S4 RSVa + ~40 RSVb pairs (trace, analysed=false, no samples/S4/RSVb), NTC qc_fail; assert per-target outputs, overview lists 5x2 cells;
 then --mode demix / aggregate / report reruns from its own results (.demix byte-equal, freyja_aggregated.tsv equal, report_data.json equal) + make_samplesheet --from-results lists (sample,target) pairs.
CUTS for a first working slice: S2a = U0-U3 (typing, split, per-target per-sample analysis, QC/consensus, demix_status, targets.tsv; aggregate/plot/report off for panels until U4). Cuttable with a
 warning: S2b-2 tail, replay_run for panels, runs_store panel support, per-target --compare_barcodes, CRAM for panels, rescue rule. Never cut: single-target legacy path, U5 modes/reruns.
RISKS: RSVa/RSVb cross-mapping inflating trace -> detected (S4 measures); uncalibrated thresholds; MAPQ ties in flu internal segments; barrier waits for all samples' typing; a missed
 join/collectFile/groupTuple silently mixes targets; big rewrite of the run-level block (verify with full e2e once per merged unit, compare numerics not bytes).
