/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { FASTP                  } from '../modules/nf-core/fastp/main'
include { BWA_INDEX              } from '../modules/nf-core/bwa/index/main'
include { BWA_MEM                } from '../modules/nf-core/bwa/mem/main'
include { SAMTOOLS_VIEW          } from '../modules/nf-core/samtools/view/main'
include { SAMTOOLS_FLAGSTAT      } from '../modules/nf-core/samtools/flagstat/main'
include { SAMTOOLS_MARKDUP       } from '../modules/local/samtools_markdup/main'
include { BAM_MAPPED_READS       } from '../modules/local/bam_mapped_reads/main'
include { MINIBWA_MEM            } from '../modules/local/minibwa_mem/main'
include { SAMTOOLS_CRAM          } from '../modules/local/samtools_cram/main'
include { READ_STARTS            } from '../modules/local/read_starts/main'
include { PRIMER_SCHEME_CHECK    } from '../modules/local/primer_scheme_check/main'
include { DROPOUT_IMPACT         } from '../modules/local/dropout_impact/main'
include { UNEXPLAINED_MUTATIONS  } from '../modules/local/unexplained_mutations/main'
include { UNEXPLAINED_MUTATIONS_CONCAT } from '../modules/local/unexplained_mutations_concat/main'
include { LINEAGE_GROUPS_CHECK   } from '../modules/local/lineage_groups_check/main'
include { SITE_AGGREGATE_CHECK   } from '../modules/local/site_aggregate_check/main'
include { SITE_AGGREGATE         } from '../modules/local/site_aggregate/main'
include { LINEAGE_SUMMARY        } from '../modules/local/lineage_summary/main'
include { TIMESERIES_TABLE       } from '../modules/local/timeseries_table/main'
include { TIMESERIES_TABLE as TIMESERIES_TABLE_LINEAGES } from '../modules/local/timeseries_table/main'
include { FREYJA_RELGROWTHRATE   } from '../modules/local/freyja_relgrowthrate/main'
include { FREYJA_DASH            } from '../modules/local/freyja_dash/main'
include { IVAR_TRIM              } from '../modules/nf-core/ivar/trim/main'
include { SAMTOOLS_SORT          } from '../modules/nf-core/samtools/sort/main'
include { FREYJA_VARIANTS        } from '../modules/local/freyja_variants/main'
include { GUNZIP_FREYJA_INPUTS   } from '../modules/local/gunzip_freyja_inputs/main'
include { FREYJA_SELECT_CONTIG   } from '../modules/local/freyja_select_contig/main'
include { FREYJA_DEMIX           } from '../modules/local/freyja_demix/main'
include { FREYJA_DEMIX as FREYJA_DEMIX_COMPARE } from '../modules/local/freyja_demix/main'
include { BARCODE_COMPARE        } from '../modules/local/barcode_compare/main'
include { FREYJA_BOOT            } from '../modules/local/freyja_boot/main'
include { BOOT_GROUP_INTERVALS   } from '../modules/local/boot_group_intervals/main'
include { FREYJA_AGGREGATE_CONTROLS } from '../modules/local/freyja_aggregate_controls/main'
include { ANNOTATE_VARIANTS      } from '../modules/local/annotate_variants/main'
include { FREYJA_AGGREGATE_PLOT  } from '../modules/local/freyja_aggregate_plot/main'
include { FILTER_AGGREGATED      } from '../modules/local/filter_aggregated/main'
include { AMPLICON_COVERAGE      } from '../modules/local/amplicon_coverage/main'
include { QC_SAMPLE              } from '../modules/local/qc_sample/main'
include { QC_CONCAT              } from '../modules/local/qc_concat/main'
include { IVAR_VARIANTS_TO_VCF   } from '../modules/local/ivar_variants_to_vcf/main'
include { CONSENSUS_FILTER_VCF   } from '../modules/local/consensus_filter_vcf/main'
include { LOW_DEPTH_MASK         } from '../modules/local/low_depth_mask/main'
include { CONSENSUS_BCFTOOLS     } from '../modules/local/consensus_bcftools/main'
include { CONSENSUS_STATS        } from '../modules/local/consensus_stats/main'
include { CONSENSUS_STATS_CONCAT } from '../modules/local/consensus_stats_concat/main'
include { CONSENSUS_CONCAT       } from '../modules/local/consensus_concat/main'
include { MUTATIONS_OF_INTEREST  } from '../modules/local/mutations_of_interest/main'
include { MUTATIONS_OF_INTEREST_CONCAT } from '../modules/local/mutations_of_interest_concat/main'
include { BUILD_REPORT_DATA      } from '../modules/local/build_report_data/main'
include { RENDER_REPORT          } from '../modules/local/render_report/main'
include { BUILD_OVERVIEW         } from '../modules/local/build_overview/main'
include { PROVENANCE             } from '../modules/local/provenance/main'
include { FREYJA_DATA            } from '../subworkflows/local/freyja_data/main'
include { PANEL_FRONT            } from '../subworkflows/local/panel_front/main'
include { FREYJA_DATA as FREYJA_DATA_COMPARE } from '../subworkflows/local/freyja_data/main'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { primerBedFor           } from '../subworkflows/local/utils_nfcore_wasteline_pipeline'
include { libraryInfo            } from '../subworkflows/local/utils_nfcore_wasteline_pipeline'
include { protocolSettings       } from '../subworkflows/local/utils_nfcore_wasteline_pipeline'
include { defaultPrimerScheme    } from '../subworkflows/local/utils_nfcore_wasteline_pipeline'
include { mutationTablesUnsupportedReason } from '../subworkflows/local/utils_nfcore_wasteline_pipeline'
include { lineageGroupsFile      } from '../subworkflows/local/utils_nfcore_wasteline_pipeline'
include { plotConfigFile         } from '../subworkflows/local/utils_nfcore_wasteline_pipeline'
include { lineageGroupsInfo      } from '../subworkflows/local/utils_nfcore_wasteline_pipeline'
include { lineageHierarchyUnsupportedReason } from '../subworkflows/local/utils_nfcore_wasteline_pipeline'
include { hierarchyMayBeMissing } from '../subworkflows/local/utils_nfcore_wasteline_pipeline'
include { bundledSchemesApply } from '../subworkflows/local/utils_nfcore_wasteline_pipeline'
include { previousRuns           } from '../subworkflows/local/utils_nfcore_wasteline_pipeline'
include { splitPaths             } from '../subworkflows/local/utils_nfcore_wasteline_pipeline'
include { multiContig            } from '../subworkflows/local/utils_nfcore_wasteline_pipeline'
include { fastaRecords           } from '../subworkflows/local/utils_nfcore_wasteline_pipeline'
include { targetsForMode         } from '../subworkflows/local/utils_nfcore_wasteline_pipeline'
include { resultsAreMultiTarget  } from '../subworkflows/local/utils_nfcore_wasteline_pipeline'
include { tkey                   } from '../subworkflows/local/utils_nfcore_wasteline_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow WASTELINE {

    take:
    mode           // string: full | demix | aggregate | report | update_db
    ch_rows        // channel: full [ meta, [ fastq_1, fastq_2 ] ] | demix [ meta, variants, depths, fastp_json ] | aggregate [ meta, demix ]
    ch_samplesheet // channel: the samplesheet in effect (path; checksummed in run_info.json, copied to summary/)
    ch_skipped     // channel: [ sample, reason, sample_type, collection_date, site, target ] of samples (of a target, for panel results; '' otherwise) that were found in results directories but not usable
    outdir         // string: --outdir

    main:

    def ch_versions = channel.empty()
    def ch_run_info = channel.empty()

    //
    // Collate and save software versions
    //
    def topic_versions = channel.topic("versions")
        .distinct()
        .branch { entry ->
            versions_file: entry instanceof Path
            versions_tuple: true
        }

    def topic_versions_string = topic_versions.versions_tuple
        .map { process, tool, version ->
            [ process[process.lastIndexOf(':')+1..-1], "  ${tool}: ${version}" ]
        }
        .groupTuple(by:0)
        .map { process, tool_versions ->
            tool_versions.unique().sort()
            "${process}:\n${tool_versions.join('\n')}"
        }

    def ch_collated_versions = softwareVersionsToYAML(ch_versions.mix(topic_versions.versions_file))
        .mix(topic_versions_string)
        .collectFile(
            storeDir: "${outdir}/pipeline_info",
            name:  'software_versions.yml',
            sort: true,
            newLine: true
        )

    if (params.report_anonymize == 'full' && mode != 'update_db' && !params.skip_report) {
        log.warn "--report_anonymize full: sample IDs and site names are replaced in report_data.json and report.html only. summary/sample_id_map.csv maps them back and must NOT be shared with the report; the other files in summary/ (qc_summary.tsv, ...) keep the real IDs."
    }

    if (mode == 'report') {
        //
        // Rebuild report.html from a finished run (--results_dir); nothing else runs
        //
        def results = file(params.results_dir)
        def report_targets = targetsForMode(mode)
        def report_template = file("${projectDir}/assets/templates/report.html.j2")
        if (report_targets[0].panel) {
            //
            // Results of a panel run: one report per target, from summary/<target>/ and samples/*/<target>/ (the output goes to
            // summary/<target>/ of --outdir). A target whose report_data.json exists is only rendered again, unless the report
            // must be anonymised (as below).
            //
            def ready = report_targets.findAll { t -> file("${results}/summary${t.summary_sub}/report_data.json").exists() && params.report_anonymize == 'none' }
            def to_build = report_targets.findAll { t -> !(t in ready) }.collect { t -> reportBuildInputs(results, t) }
            def ch_to_build = channel.fromList(to_build).multiMap { item ->
                data: item.data
                extra: item.extra
                analyses: item.analyses
                groups: item.groups
                contigs: item.contigs
            }
            BUILD_REPORT_DATA(ch_to_build.data, ch_to_build.extra, ch_to_build.analyses, ch_to_build.groups, ch_to_build.contigs)
            RENDER_REPORT(
                channel.fromList(ready.collect { t -> [ tkey([ target: t ]), file("${results}/summary${t.summary_sub}/report_data.json") ] }).mix(BUILD_REPORT_DATA.out.json),
                report_template,
            )
            // the panel overview of the results, once every report of this call exists (it links to them); a run from before
            // summary/targets.tsv has none
            def results_targets_tsv = file("${results}/summary/targets.tsv")
            if (!params.skip_report && results_targets_tsv.exists()) {
                def overview_items = report_targets.collect { t ->
                    def run_info = file("${results}/pipeline_info${t.summary_sub}/run_info.json")
                    [ t.name, optionalFile(results, 'freyja_aggregated.tsv', t.summary_sub), optionalFile(results, 'qc_summary.tsv', t.summary_sub), optionalFile(results, 'demix_status.tsv', t.summary_sub), run_info.exists() ? run_info : [] ]
                }
                def ch_overview_in = RENDER_REPORT.out.html.toList().map { _reports -> overviewInputs(report_targets.collect { t -> t.name }, overview_items) }
                BUILD_OVERVIEW(
                    channel.value(results_targets_tsv),
                    ch_overview_in.map { manifest, _files -> manifest },
                    ch_overview_in.map { _manifest, files -> files },
                    file("${projectDir}/assets/templates/overview.html.j2"),
                )
            }
            else if (!params.skip_report) {
                log.info "No summary/overview.html: ${results}/summary/targets.tsv does not exist (the results were written by a version without the panel overview)."
            }
        }
        else {
            def existing_data = file("${results}/summary/report_data.json")
            // (a single-target run: one report, from summary/)
            def report_meta = tkey([ target: report_targets[0] ])
            // report_data.json is rebuilt from the summary files when the report must be anonymised: the renderer has no
            // anonymisation of its own, and an existing report_data.json holds the real IDs
            if (existing_data.exists() && params.report_anonymize == 'none') {
                // the run already wrote report_data.json: only the HTML is rendered again
                RENDER_REPORT(channel.value([ report_meta, existing_data ]), report_template)
            }
            else {
                // rebuild report_data.json from the files the run left behind
                def depths = (file("${results}/samples/*/freyja/*.depths") + file("${results}/*/freyja/*.depths")).sort { a, b -> a.name <=> b.name }
                BUILD_REPORT_DATA(channel.value([
                    report_meta,
                    file("${results}/pipeline_info/run_info.json"),
                    optionalFile(results, 'qc_summary.tsv'),
                    optionalFile(results, 'amplicon_coverage.tsv'),
                    optionalFile(results, 'freyja_aggregated.tsv'),
                    optionalFile(results, 'demix_status.tsv'),
                    optionalFile(results, 'consensus_stats.tsv'),
                    optionalFile(results, 'mutations_of_interest.tsv'),
                    optionalFile(results, 'barcode_compare.tsv'),
                    depths,
                ]), channel.value([
                    optionalFile(results, 'freyja_aggregated_controls.tsv'),
                    optionalFile(results, 'samplesheet.csv'),
                    optionalFile(results, 'freyja_boot_summary.tsv'),
                    optionalFile(results, 'barcode_compare_status.tsv'),
                ]), channel.value([
                    optionalFile(results, 'dropout_impact.tsv'),
                    optionalFile(results, 'dropout_impact_summary.tsv'),
                    optionalFile(results, 'unexplained_mutations.tsv'),
                    optionalFile(results, 'freyja_timeseries.tsv'),
                    optionalFile(results, 'freyja_timeseries_lineages.tsv'),
                    optionalFile(results, 'freyja_growth_rates.csv'),
                    optionalFile(results, 'freyja_growth_rates.failed.txt'),
                    optionalFile(results, 'freyja_growth_inputs.json'),
                    file("${results}/summary/freyja_dashboard.html").exists() ? 'freyja_dashboard.html' : '',
                ]), channel.value([
                    optionalFile(results, 'lineage_summary.csv'),
                    optionalFile(results, 'lineage_groups_resolved.tsv'),
                    optionalFile(results, 'site_aggregates_freyja.tsv'),
                    optionalFile(results, 'site_aggregates_meta.tsv'),
                ]), channel.value([
                    optionalFile(results, 'contig_coverage.tsv'),
                    optionalFile(results, 'consensus_stats_contigs.tsv'),
                ]))
                RENDER_REPORT(BUILD_REPORT_DATA.out.json, report_template)
            }
        }
    }
    else {
    //
    // Freyja dataset: resolved once per target, shared by every sample of the target. update_db always fetches the latest.
    //
    // One dataset per target (a single-target run has one): the *_by_target outputs are [ target name, item ], joined onto the
    // samples and the run-level steps of a target with combine(..., by: 0) on the target name.
    // The targets of the run: --panel, or the targets of the panel results / samplesheet a mode reads, or the one target of a normal run
    // (targetsForMode); --targets selects among them. A target's layout follows the input: panel results give panel-layout outputs
    // (samples/<id>/<target>/, summary/<target>/) even when only one target is selected.
    // update_db fetches the latest dataset of every target: a panel row's own barcodes (an offline file) must not stand in for it.
    def targets = targetsForMode(mode).collect { t -> mode == 'update_db' ? t + [ barcodes: null ] : t }
    def is_panel = targets[0].panel
    FREYJA_DATA(mode == 'update_db' ? 'latest' : params.barcodes, targets, '--barcodes')

    if (mode == 'update_db') {
        // Nothing else to do: the dataset is now in --db_dir/freyja/<barcode-version>/
        FREYJA_DATA.out.freyja_data_by_target.subscribe { name, freyja_data ->
            def info = new groovy.json.JsonSlurper().parse(freyja_data)
            log.info "\nFreyja data${is_panel ? " for ${name}" : ''} updated: barcode version ${info.barcode_version} (${info.barcode_version_iso}), stored in ${info.dataset_path}\nUse it with: --barcodes ${info.barcode_version_iso.take(10)}   (or --barcodes ${info.barcode_version_iso} for this exact update)"
        }
    }
    else {
        // Every run-level step is per target (see "Run-level steps, PER TARGET" below); a single-target run is a panel of one.
        // `target` is the one target of a single-target run (null for a panel run, whose references come from meta.target).
        def target = is_panel ? null : targets[0]
        // --fasta is required by modes full and demix (validateInputParameters); --gff is optional for pathogens without a bundled annotation.
        // ch_fasta serves the read mapping (BWA_INDEX / BWA_MEM, until the panel front end replaces it) and the run provenance; the per-sample
        // steps get their reference from meta.target (referenceFasta / referenceGff), in lockstep with the sample through a multiMap.
        def ch_fasta = channel.value(!is_panel && target.fasta ? file(target.fasta, checkIfExists: true) : [])
        // A reference with several contigs (PLAN D25; validated at start-up: --freyja_contig names the one Freyja demixes). Every
        // branch for it is guarded by this; a single-contig run does not change.
        def multi_contig = (mode == 'full' || mode == 'demix') && targets.any { t -> multiContig(t) }

        //
        // The lineage groups file of each target (--lineage_groups, the same file for every target of a panel, else the target's bundled
        // one) and the sites CSV (run-wide) are checked FIRST, in parallel with fastp and the mapping, so that a
        // typo stops the run within minutes (these processes keep the default error strategy). Unknown lineages and
        // unknown site names are warnings, shown here. A dataset without a lineage hierarchy (custom barcodes file) cannot
        // resolve the groups: the file format is still checked, the summary of that target is skipped with one warning.
        //
        // A target without bundled lineage groups (conf/pathogens/*.config, assets/panels/registry.csv) has the summary off unless
        // --lineage_groups is given (PLAN D22).
        def groups_by_target = params.lineage_summary
            ? targets.collectEntries { t -> [ (t.name): lineageGroupsFile(t) ] }.findAll { _name, groups_file -> groups_file }
            : [:]
        def groups_targets = targets.findAll { t -> groups_by_target.containsKey(t.name) }
        def run_site_aggregation = params.site_aggregation as boolean
        if (params.lineage_summary) {
            targets.findAll { t -> !groups_by_target.containsKey(t.name) }.each { t ->
                if (is_panel) {
                    log.info "Lineage-group summary off for target ${t.name}: no lineage groups are bundled for it. Give --lineage_groups <CSV> (it applies to every target of the panel; the target also needs a Freyja dataset with a lineage hierarchy) to get summary/${t.name}/lineage_summary.csv."
                }
                else {
                    log.info "Lineage-group summary off: no lineage groups are bundled for --pathogen ${t.name}. Give --lineage_groups <CSV> (and a Freyja dataset with a lineage hierarchy) to get summary/lineage_summary.csv."
                }
            }
            if (is_panel && params.lineage_groups) {
                log.info "Lineage-group summary: --lineage_groups ${params.lineage_groups} applies to every target of the panel (${targets.collect { t -> t.name }.join(', ')}); the summary of a target is skipped if its Freyja dataset has no lineage hierarchy."
            }
        }
        // why the hierarchy of a target's dataset is unavailable (known before the run), and the targets whose downloaded dataset may lack it
        def hierarchy_reasons = groups_targets.collectEntries { t -> [ (t.name): lineageHierarchyUnsupportedReason(t) ] }.findAll { _name, reason -> reason }
        def hierarchy_may_miss = groups_targets.findAll { t -> !hierarchy_reasons.containsKey(t.name) && hierarchyMayBeMissing(t) }.collect { t -> t.name }
        hierarchy_reasons.each { name, reason ->
            log.warn "Skipping the lineage-group summary${is_panel ? " of target ${name}" : ''}: ${reason}. Use a Freyja dataset directory with lineages.yml (--barcodes latest, container or a dataset directory) to get summary${is_panel ? "/${name}" : ''}/lineage_summary.csv."
        }
        // [ target name, dataset, has a lineage hierarchy ] of each target with lineage groups. A dataset downloaded during the run for a
        // pathogen that may ship no hierarchy (RSV, dengue, ...) is looked at first: without lineages.yml the target's summary is skipped, with ONE warning.
        def ch_hierarchy = channel.empty()
        if (groups_targets) {
            def ch_groups_datasets = FREYJA_DATA.out.dataset_by_target.filter { name, _dataset -> groups_by_target.containsKey(name) }
            def ch_hierarchy_known = ch_groups_datasets
                .filter { name, _dataset -> !(name in hierarchy_may_miss) }
                .map { name, dataset -> [ name, dataset, !hierarchy_reasons.containsKey(name) ] }
            def ch_hierarchy_looked_at = ch_groups_datasets
                .filter { name, _dataset -> name in hierarchy_may_miss }
                .join(FREYJA_DATA.out.freyja_data_by_target)
                .map { name, dataset, freyja_data -> [ name, dataset, new groovy.json.JsonSlurper().parse(freyja_data).lineages_yml_sha256 ? true : false ] }
            ch_hierarchy = ch_hierarchy_known.mix(ch_hierarchy_looked_at)
            ch_hierarchy_looked_at.filter { _name, _dataset, has -> !has }.subscribe { name, _dataset, _has ->
                log.warn "Skipping the lineage-group summary: the Freyja dataset downloaded for ${is_panel ? 'target' : '--pathogen'} ${name} has no lineage hierarchy (lineages.yml), which the groups need."
            }
            // the check of each target's groups file; its dataset is [] when it has no hierarchy (format check only)
            def ch_to_check = channel.fromList(groups_targets.collect { t -> [ t.name, tkey([ target: t ]), groups_by_target[t.name] ] })
                .combine(ch_hierarchy, by: 0)
                .multiMap { _name, key, groups_file, dataset, has ->
                    groups: [ key, groups_file ]
                    dataset: has ? dataset : []
                }
            LINEAGE_GROUPS_CHECK(ch_to_check.groups, ch_to_check.dataset)
            LINEAGE_GROUPS_CHECK.out.report.subscribe { key, report ->
                report.readLines().findAll { line -> line.contains('WARNING') }.each { line -> log.warn((is_panel ? "[${key.id}] " : '') + line) }
            }
        }
        if (run_site_aggregation) {
            SITE_AGGREGATE_CHECK(
                channel.value(file(params.site_aggregation, checkIfExists: true)),
                ch_samplesheet.first(),
            )
            SITE_AGGREGATE_CHECK.out.report.subscribe { report ->
                report.readLines().findAll { line -> line.contains('warning') }.each { line -> log.warn line }
            }
        }
        // Library settings (protocol, primer trimming, duplicate handling, primer scheme) are PER SAMPLE (PLAN D20):
        // libraryMeta() puts them into every sample's meta (meta.protocol, primer_trimming, duplicates, primer_scheme,
        // primer_bed). Amplicon-level QC and primer trimming need the sample's primer BED (primerBedFor); probe-capture
        // reads have none (demix mode: the scheme of old results is unknown, so only when a scheme is named). Each step
        // below takes its own samples from the shared channels, so a run can mix libraries.

        //
        // Per-sample Freyja inputs [ meta, variants, depths ]: made from reads (full) or given (demix)
        //
        def ch_target_metas = ch_rows.map { row -> row[0] }   // [ meta ] of the (sample, target) pairs the run analyses (a panel run: the pairs the typing chose, set below)
        def ch_freyja_in = channel.empty()
        def ch_fastp = channel.empty()      // [ meta, fastp.json ]
        def ch_flagstat = channel.empty()   // [ meta, flagstat ]
        def ch_ivar_log = channel.empty()   // [ meta, ivar trim log ]       (amplicon protocol with trimming)
        def ch_scheme_summary = channel.empty() // [ meta, primer scheme check summary ] (same)

        if (mode == 'full') {
            def ch_fasta_meta = ch_fasta.map { fasta -> [ [ id: fasta.baseName ], fasta ] }

            //
            // Reads -> fastp
            //
            FASTP(ch_rows.map { meta, reads -> [ meta, reads, [] ] }, false, false, false)
            ch_fastp = FASTP.out.json

            //
            // Alignment: bwa mem or minibwa (sorted) -> mapped reads only (+ index) -> primer trimming -> sorted + indexed BAM
            // A panel run (--panel) maps every sample once to ALL the panel's references instead, types it, and cuts one BAM per
            // analysed (sample, target) out of it (PANEL_FRONT): from there on the chain is the same, once per (sample, target) with
            // the target in meta.target. fastp's JSON (QC) belongs to the sample: it is handed to each of its (sample, target) metas.
            //
            def ch_counted = channel.empty()   // [ meta, bam, bai, mapped_reads.txt ] of the mapped reads, one per analysed sample (or sample and target)
            if (is_panel) {
                PANEL_FRONT(FASTP.out.reads, targets, mode)
                ch_counted = PANEL_FRONT.out.mapped
                ch_target_metas = ch_counted.map { meta, _bam, _bai, _count -> meta }
                def ch_fastp_by_id = FASTP.out.json.map { meta, json -> [ meta.id, json ] }
                ch_fastp = ch_counted
                    .map { meta, _bam, _bai, _count -> [ meta.id, meta ] }
                    .combine(ch_fastp_by_id, by: 0)
                    .map { _id, meta, json -> [ meta, json ] }
            }
            else {
                def ch_aligned = channel.empty()
                if (params.aligner == 'minibwa') {
                    // experimental, docs/aligner_evaluation.md: the index is built inside the task
                    MINIBWA_MEM(FASTP.out.reads, ch_fasta_meta, true)
                    ch_aligned = MINIBWA_MEM.out.bam
                }
                else {
                    BWA_INDEX(ch_fasta_meta)
                    BWA_MEM(FASTP.out.reads, BWA_INDEX.out.index, ch_fasta_meta, true)
                    ch_aligned = BWA_MEM.out.bam
                }

                SAMTOOLS_VIEW(
                    ch_aligned.map { meta, bam -> [ meta, bam, [] ] },
                    [ [:], [], [] ],
                    [ [:], [] ],
                    [ [:], [] ],
                    'bai',
                )

                BAM_MAPPED_READS(SAMTOOLS_VIEW.out.bam.join(SAMTOOLS_VIEW.out.bai))
                ch_counted = BAM_MAPPED_READS.out.counted
            }
            def ch_mapped = ch_counted
                .map { meta, bam, bai, count_file -> [ meta, bam, bai, count_file.text.trim() as long ] }
                .branch { _meta, _bam, _bai, count ->
                    mapped: count > 0
                    empty: true
                }

            // Samples whose reads are primer-trimmed (meta.primer_trimming: amplicon protocol, not
            // --skip_primer_trimming) go through ivar trim with THEIR primer BED. The BED travels beside the BAM in
            // a multiMap, which emits both channels in lockstep, so the two queue inputs of IVAR_TRIM stay paired.
            // `ivar trim` crashes on a BAM without mapped reads (e.g. a clean negative control): those samples skip
            // trimming and carry their (empty) BAM on; the coverage gate then records them as
            // skipped_low_coverage downstream.
            def ch_to_trim = ch_mapped.mapped
                .filter { meta, _bam, _bai, _count -> meta.primer_trimming }
                .multiMap { meta, bam, bai, _count ->
                    reads: [ meta, bam, bai ]
                    bed: primerBedFor(meta, meta.target)
                }
            IVAR_TRIM(ch_to_trim.reads, ch_to_trim.bed)
            ch_ivar_log = IVAR_TRIM.out.log

            // Does the primer scheme match the reads? Needs the mapped BAM BEFORE trimming (trimming moves
            // the 5' ends). Declared scheme = the sample's own (meta.primer_scheme; 'custom' = its BED).
            // (decided per target: a panel can mix targets with and without bundled schemes)
            if (params.primer_scheme_check && targets.any { t -> bundledSchemesApply(t) }) {
                READ_STARTS(ch_mapped.mapped
                    .filter { meta, _bam, _bai, _count -> meta.primer_trimming && bundledSchemesApply(meta.target) }
                    .map { meta, bam, _bai, _count -> [ meta, bam ] })
                def ch_to_check = READ_STARTS.out.tsv.multiMap { meta, tsv ->
                    starts: [ meta, tsv ]
                    declared: meta.primer_scheme
                    custom_bed: meta.primer_bed ? primerBedFor(meta, meta.target) : []
                }
                PRIMER_SCHEME_CHECK(
                    ch_to_check.starts,
                    ch_to_check.declared,
                    file("${projectDir}/assets/primer_schemes", checkIfExists: true),
                    ch_to_check.custom_bed,
                )
                ch_scheme_summary = PRIMER_SCHEME_CHECK.out.summary
            }

            // Probe-capture (or pre-trimmed) reads are used as aligned: no primer trimming, straight to the
            // sorted + indexed BAM (empty BAMs sort fine and keep flowing)
            def ch_untrimmed = ch_mapped.mapped.mix(ch_mapped.empty)
                .filter { meta, _bam, _bai, _count -> !meta.primer_trimming }
                .map { meta, bam, _bai, _count -> [ meta, bam ] }
            SAMTOOLS_SORT(IVAR_TRIM.out.bam.mix(ch_untrimmed), [ [:], [], [] ], 'bai')

            def ch_final_bam = SAMTOOLS_SORT.out.bam.join(SAMTOOLS_SORT.out.index)
                .mix(ch_mapped.empty
                    .filter { meta, _bam, _bai, _count -> meta.primer_trimming }
                    .map { meta, bam, bai, _count -> [ meta, bam, bai ] })

            // PCR duplicates are flagged (never removed) for the samples whose meta.duplicates says so (--duplicates
            // auto: capture samples); Freyja's pileup skips flagged reads. Near-empty BAMs pass through (checked:
            // markdup accepts a header-only BAM).
            def ch_to_mark = ch_final_bam.branch { meta, _bam, _bai ->
                mark: meta.duplicates == 'mark'
                keep: true
            }
            SAMTOOLS_MARKDUP(ch_to_mark.mark)
            ch_final_bam = SAMTOOLS_MARKDUP.out.bam.mix(ch_to_mark.keep)

            // --bam_format cram: a reference-based CRAM is published instead of the BAM (a leaf: Freyja, QC and
            // the consensus keep using the BAM, which stays in the work directory)
            if (params.save_bam && params.bam_format == 'cram') {
                def ch_to_cram = ch_final_bam.multiMap { meta, bam, bai ->
                    aligned: [ meta, bam, bai ]
                    fasta: referenceFasta(meta)
                }
                SAMTOOLS_CRAM(ch_to_cram.aligned, ch_to_cram.fasta, [])
            }

            def ch_to_variants = ch_final_bam.multiMap { meta, bam, bai ->
                aligned: [ meta, bam, bai ]
                fasta: referenceFasta(meta)
            }
            FREYJA_VARIANTS(ch_to_variants.aligned, ch_to_variants.fasta)

            ch_freyja_in = FREYJA_VARIANTS.out.variants.join(FREYJA_VARIANTS.out.depths)

            // mapped-read counts for the QC table, from the final BAM
            SAMTOOLS_FLAGSTAT(ch_final_bam)
            ch_flagstat = SAMTOOLS_FLAGSTAT.out.flagstat
        }
        else if (mode == 'demix') {
            // Existing variants/depths from a samplesheet or --previous_results (read only). Gzipped
            // legacy files are unpacked once so that every consumer (QC, consensus, report) reads plain text.
            def ch_given = ch_rows
                .map { meta, variants, depths, _fastp_json -> [ meta, variants, depths ] }
                .branch { _meta, variants, depths ->
                    zipped: variants.name.endsWith('.gz') || depths.name.endsWith('.gz')
                    plain: true
                }
            GUNZIP_FREYJA_INPUTS(ch_given.zipped)
            ch_freyja_in = ch_given.plain.mix(GUNZIP_FREYJA_INPUTS.out.plain)
            ch_fastp = ch_rows
                .filter { _meta, _variants, _depths, fastp_json -> fastp_json }
                .map { meta, _variants, _depths, fastp_json -> [ meta, fastp_json ] }
        }

        def ch_depths = ch_freyja_in.map { meta, _variants, depths -> [ meta, depths ] }
        def ch_variants = ch_freyja_in.map { meta, variants, _depths -> [ meta, variants ] }

        // What Freyja (demix, bootstrap, the coverage gate) and the mutation tables get: with several reference contigs only the
        // variants and depths of the demixed contig (a view made by FREYJA_SELECT_CONTIG; the files under samples/ stay complete);
        // QC, amplicon coverage, the consensus and the report keep the all-contig files above. Single contig: the same channels.
        // Decided PER SAMPLE from its target (a run of several targets can mix single- and multi-contig references).
        def ch_freyja_view = ch_freyja_in
        if (mode == 'full' || mode == 'demix') {
            def ch_by_contigs = ch_freyja_in.branch { meta, _variants, _depths ->
                multi: multiContig(meta.target)
                single: true
            }
            def ch_to_select = ch_by_contigs.multi.multiMap { meta, variants, depths ->
                files: [ meta, variants, depths ]
                contig: meta.target.freyja_contig
            }
            FREYJA_SELECT_CONTIG(ch_to_select.files, ch_to_select.contig)
            ch_freyja_view = ch_by_contigs.single.mix(FREYJA_SELECT_CONTIG.out.view)
        }
        def ch_view_variants = ch_freyja_view.map { meta, variants, _depths -> [ meta, variants ] }
        def ch_view_depths = ch_freyja_view.map { meta, _variants, depths -> [ meta, depths ] }

        //
        // Demix per sample (never fails the run) -> [ meta, demix ]; aggregate mode takes given .demix files
        //
        def ch_demixed = channel.empty()
        def ch_compare_tsv = channel.empty()
        def ch_compare_status = channel.empty()
        def ch_boot_summary = channel.empty()
        def ch_boot_raw = channel.empty()
        // --compare_barcodes: one spec for every target (validated per target at start-up)
        def run_compare = params.compare_barcodes && (mode == 'full' || mode == 'demix')
        // aggregate mode with --previous_results: samples that were low_confidence in the earlier run keep that status (PLAN D18)
        // (keyed like meta.key: the sample ID, or ID__target for the targets of a panel)
        def carried_low = mode == 'aggregate' ? carriedLowConfidence(targets) : [:]
        // demix_status.tsv lines of each target: [ tkey, line ] with line = sample, status, reason, genome_coverage_pct, barcode_version, barcode_version_iso
        // (a sample found in old results but not usable: its row goes to its target when the results name one, else to every target of the run, which is one in single-target results)
        def ch_status_lines = ch_skipped.flatMap { sample, reason, _type, _date, _site, target_name ->
            def line = [ sample, 'missing', reason, '', '', '' ].collect { value -> value.toString().replaceAll('[\t\r\n]+', ' ') }.join('\t')
            targets.findAll { t -> !target_name || t.name == target_name }.collect { t -> [ tkey([ target: t ]), line ] }
        }

        if (mode == 'full' || mode == 'demix') {
            // each sample is demixed with the dataset of ITS target: [ target name, dataset, barcode_version ] joined on meta.target.name
            def ch_demix_data = FREYJA_DATA.out.dataset_by_target.join(FREYJA_DATA.out.barcode_version_by_target)
            def ch_to_demix = ch_freyja_view
                .map { meta, variants, depths -> [ meta.target.name, meta, variants, depths ] }
                .combine(ch_demix_data, by: 0)
                .multiMap { _name, meta, variants, depths, dataset, barcode_version ->
                    sample: [ meta, variants, depths ]
                    dataset: dataset
                    barcode_version: barcode_version
                }
            FREYJA_DEMIX(ch_to_demix.sample, ch_to_demix.dataset, ch_to_demix.barcode_version)

            ch_status_lines = ch_status_lines.mix(FREYJA_DEMIX.out.status
                .map { meta, status_json ->
                    def rec = new groovy.json.JsonSlurper().parse(status_json)
                    def line = [ rec.sample, rec.status, rec.reason, rec.genome_coverage_pct, rec.barcode_version, rec.barcode_version_iso ]
                        .collect { value -> value == null ? '' : value.toString().replaceAll('[\t\r\n]+', ' ') }
                        .join('\t')
                    [ tkey(meta), line ]
                })

            ch_demixed = FREYJA_DEMIX.out.demix

            //
            // --compare_barcodes: the same variants/depths demixed with a SECOND dataset (resolved
            // independently, recorded as freyja_data_compare); abundances of both side by side
            //
            if (run_compare) {
                FREYJA_DATA_COMPARE(params.compare_barcodes, targets, '--compare_barcodes')
                def ch_compare_data = FREYJA_DATA_COMPARE.out.dataset_by_target.join(FREYJA_DATA_COMPARE.out.barcode_version_by_target)
                def ch_to_compare = ch_freyja_view
                    .map { meta, variants, depths -> [ meta.target.name, meta, variants, depths ] }
                    .combine(ch_compare_data, by: 0)
                    .multiMap { _name, meta, variants, depths, dataset, barcode_version ->
                        sample: [ meta, variants, depths ]
                        dataset: dataset
                        barcode_version: barcode_version
                    }
                FREYJA_DEMIX_COMPARE(ch_to_compare.sample, ch_to_compare.dataset, ch_to_compare.barcode_version)
                // the barcode versions of the two datasets of the sample's target travel with it
                def ch_to_barcode_compare = FREYJA_DEMIX.out.demix.join(FREYJA_DEMIX_COMPARE.out.demix)
                    .map { meta, demix_a, demix_b -> [ meta.target.name, meta, demix_a, demix_b ] }
                    .combine(FREYJA_DATA.out.barcode_version_by_target, by: 0)
                    .combine(FREYJA_DATA_COMPARE.out.barcode_version_by_target, by: 0)
                    .multiMap { _name, meta, demix_a, demix_b, version_a, version_b ->
                        pair: [ meta, demix_a, demix_b ]
                        version_a: version_a
                        version_b: version_b
                    }
                BARCODE_COMPARE(ch_to_barcode_compare.pair, ch_to_barcode_compare.version_a, ch_to_barcode_compare.version_b)
                // summary<sub>/barcode_compare.tsv: the rows of the samples of each target, after the header
                ch_compare_tsv = groupByTarget(BARCODE_COMPARE.out.rows)
                    .map { key, files ->
                        [ key, summaryFile(outdir, key, 'barcode_compare.tsv', 'sample\tlineage\tabundance_a\tabundance_b\tbarcode_version_a\tbarcode_version_b\n' + files.collect { f -> f.text }.join('')) ]
                    }
                // outcome of the second demix per sample (a failed one is flagged in the report, never fatal)
                ch_compare_status = FREYJA_DEMIX_COMPARE.out.status
                    .map { meta, status_json ->
                        def rec = new groovy.json.JsonSlurper().parse(status_json)
                        def line = [ rec.sample, rec.status, rec.reason, rec.barcode_version ]
                            .collect { value -> value == null ? '' : value.toString().replaceAll('[\t\r\n]+', ' ') }
                            .join('\t')
                        [ tkey(meta), line ]
                    }
                    .groupTuple()
                    .map { key, lines ->
                        [ key, summaryFile(outdir, key, 'barcode_compare_status.tsv', ([ 'sample\tstatus\treason\tbarcode_version' ] + lines.sort(false)).join('\n') + '\n') ]
                    }
            }

            //
            // --freyja_boot: bootstrap intervals for every sample that demixed
            //
            if (params.freyja_boot) {
                def ch_to_boot = ch_freyja_view.join(FREYJA_DEMIX.out.demix)
                    .map { meta, variants, depths, demix -> [ meta.target.name, meta, variants, depths, demix ] }
                    .combine(FREYJA_DATA.out.dataset_by_target, by: 0)
                    .multiMap { _name, meta, variants, depths, demix, dataset ->
                        sample: [ meta, variants, depths, demix ]
                        dataset: dataset
                    }
                FREYJA_BOOT(ch_to_boot.sample, ch_to_boot.dataset)
                // the raw draws of the samples of each target [ tkey, [ files ] ] (for the group intervals)
                ch_boot_raw = groupByTarget(FREYJA_BOOT.out.raw)
                // summary<sub>/freyja_boot_summary.tsv: one table of the samples of each target (the per-sample intervals are published under samples/)
                ch_boot_summary = groupByTarget(FREYJA_BOOT.out.summary)
                    .map { key, files ->
                        [ key, summaryFile(outdir, key, 'freyja_boot_summary.tsv', 'sample\tlineage\tabundance\tlo\thi\n' + files.collect { f -> f.text }.join('')) ]
                    }
            }
        }
        else {
            ch_demixed = ch_rows
            ch_status_lines = ch_status_lines.mix(ch_demixed
                .filter { meta, _demix -> carried_low.containsKey(meta.key) }
                .map { meta, _demix ->
                    [ tkey(meta), [ meta.id, 'low_confidence', carried_low[meta.key].reason, carried_low[meta.key].coverage, '', '' ].join('\t') ]
                })
        }

        //
        // Coverage QC (every sample, negatives included) from the Freyja depths file, plus flagstat
        // (mode full), the fastp JSON (when known) and per-amplicon coverage
        //
        // The per-target tables ([ tkey, file ], one per target that has the table; a target without one has no item)
        def ch_qc_summary = channel.empty()
        def ch_amplicon_coverage = channel.empty()
        def ch_consensus_stats = channel.empty()
        def ch_moi = channel.empty()
        def ch_depth_files = channel.empty()            // [ tkey, [ depths files ] ]
        def ch_contig_coverage = channel.empty()        // several reference contigs only: summary/contig_coverage.tsv
        def ch_consensus_contigs = channel.empty()      // several reference contigs only: summary/consensus_stats_contigs.tsv

        if (mode == 'full' || mode == 'demix') {
            // only the samples with a primer scheme (meta.primer_scheme), each against its own BED
            def ch_to_cover = ch_depths
                .filter { meta, _depths -> meta.primer_scheme }
                .multiMap { meta, depths ->
                    depths: [ meta, depths ]
                    bed: primerBedFor(meta, meta.target)
                }
            AMPLICON_COVERAGE(ch_to_cover.depths, ch_to_cover.bed)
            def ch_amplicons = AMPLICON_COVERAGE.out.tsv

            def ch_qc_in = ch_depths
                .join(ch_flagstat, remainder: true)
                .join(ch_fastp, remainder: true)
                .join(ch_amplicons, remainder: true)
                .join(ch_ivar_log, remainder: true)
                .join(ch_scheme_summary, remainder: true)
                .map { meta, depths, flagstat, fastp_json, amplicons, ivar_log, scheme_summary ->
                    [ meta, depths, flagstat ?: [], fastp_json ?: [], amplicons ?: [], ivar_log ?: [], scheme_summary ?: [] ]
                }
            QC_SAMPLE(ch_qc_in)

            // the run-level tables of the samples of each target
            QC_CONCAT(
                groupByTarget(QC_SAMPLE.out.tsv)
                    .join(groupByTarget(ch_amplicons), remainder: true)
                    .filter { _key, qc_rows, _amplicons -> qc_rows }
                    .map { key, qc_rows, amplicons -> [ key, qc_rows, amplicons ?: [] ] }
            )
            ch_qc_summary = QC_CONCAT.out.qc_summary
            ch_amplicon_coverage = QC_CONCAT.out.amplicon_coverage
            ch_depth_files = groupByTarget(ch_depths)
            // several reference contigs: one row per sample and contig (summary/contig_coverage.tsv, for the segments table of the report)
            if (multi_contig) {
                ch_contig_coverage = groupByTarget(QC_SAMPLE.out.contig_tsv)
                    .map { key, files -> [ key, summaryFile(outdir, key, 'contig_coverage.tsv', concatKeepingHeader(files)) ] }
            }

            //
            // Consensus genomes from the Freyja variants + depths (negative controls get none).
            // A sample whose consensus steps fail still gets a row, with status failed.
            //
            if (!params.skip_consensus) {
                def ch_consensus_samples = ch_variants.filter { meta, _variants -> meta.sample_type != 'negative' }
                def ch_to_vcf = ch_consensus_samples.multiMap { meta, variants ->
                    sample: [ meta, variants ]
                    fasta: referenceFasta(meta)
                }
                IVAR_VARIANTS_TO_VCF(ch_to_vcf.sample, ch_to_vcf.fasta)
                CONSENSUS_FILTER_VCF(IVAR_VARIANTS_TO_VCF.out.vcf)
                def ch_to_mask = ch_depths.join(CONSENSUS_FILTER_VCF.out.bed).multiMap { meta, depths, bed ->
                    sample: [ meta, depths, bed ]
                    fasta: referenceFasta(meta)
                }
                LOW_DEPTH_MASK(ch_to_mask.sample, ch_to_mask.fasta)
                def ch_to_consensus = CONSENSUS_FILTER_VCF.out.vcf.join(LOW_DEPTH_MASK.out.bed).multiMap { meta, vcf, tbi, mask ->
                    sample: [ meta, vcf, tbi, mask ]
                    fasta: referenceFasta(meta)
                }
                CONSENSUS_BCFTOOLS(ch_to_consensus.sample, ch_to_consensus.fasta)
                def ch_to_stats = CONSENSUS_BCFTOOLS.out.fasta.join(CONSENSUS_BCFTOOLS.out.status).multiMap { meta, fasta, status_json ->
                    sample: [ meta, fasta, status_json ]
                    reference: referenceFasta(meta)
                }
                CONSENSUS_STATS(ch_to_stats.sample, ch_to_stats.reference)

                // the run-level tables of the samples of each target. A sample whose consensus steps failed gets a `failed` row (a file
                // named by meta.key, which is unique per sample and target)
                def ch_stats = ch_consensus_samples
                    .join(CONSENSUS_STATS.out.tsv, remainder: true)
                    .branch { _meta, _variants, stats ->
                        built: stats
                        failed: true
                    }
                def ch_failed_stats = ch_stats.failed.map { meta, _variants, _stats ->
                    [ meta, workTextFile("${meta.key}.consensus_stats.tsv", "sample\tlength\tn_count\tconsensus_coverage_pct\tambiguous_count\tstatus\treason\n${meta.id}\tNA\tNA\tNA\tNA\tfailed\tconsensus steps did not complete\n") ]
                }
                CONSENSUS_STATS_CONCAT(groupByTarget(ch_stats.built.map { meta, _variants, stats -> [ meta, stats ] }.mix(ch_failed_stats)))
                CONSENSUS_CONCAT(groupByTarget(CONSENSUS_BCFTOOLS.out.fasta))
                ch_consensus_stats = CONSENSUS_STATS_CONCAT.out.tsv
                // several reference contigs: one row per sample and contig (summary/consensus_stats_contigs.tsv)
                if (multi_contig) {
                    ch_consensus_contigs = groupByTarget(CONSENSUS_STATS.out.contig_tsv)
                        .map { key, files -> [ key, summaryFile(outdir, key, 'consensus_stats_contigs.tsv', concatKeepingHeader(files)) ] }
                }

                // --annotate_variants: amino-acid consequences of the variants applied to the consensus
                if (params.annotate_variants) {
                    def ch_to_annotate = CONSENSUS_FILTER_VCF.out.vcf.multiMap { meta, vcf, tbi ->
                        sample: [ meta, vcf, tbi ]
                        fasta: referenceFasta(meta)
                        gff: referenceGff(meta)
                    }
                    ANNOTATE_VARIANTS(ch_to_annotate.sample, ch_to_annotate.fasta, ch_to_annotate.gff)
                    groupByTarget(ANNOTATE_VARIANTS.out.tsv).subscribe { key, files ->
                        summaryFile(outdir, key, 'annotated_variants.tsv', concatKeepingHeader(files))
                    }
                }
            }

            //
            // Mutations of interest (opt-in)
            //
            if (params.mutations_of_interest) {
                def ch_to_moi = ch_view_variants.join(ch_view_depths).multiMap { meta, variants, depths ->
                    sample: [ meta, variants, depths ]
                    gff: referenceGff(meta)
                    fasta: referenceFasta(meta)
                }
                MUTATIONS_OF_INTEREST(ch_to_moi.sample, file(params.mutations_of_interest, checkIfExists: true), ch_to_moi.gff, ch_to_moi.fasta)
                MUTATIONS_OF_INTEREST_CONCAT(groupByTarget(MUTATIONS_OF_INTEREST.out.tsv))
                ch_moi = MUTATIONS_OF_INTEREST_CONCAT.out.tsv
            }
        }

        //
        // Mutation-level tables per sample: the variants that no detected lineage explains (non-negative samples). They, like the
        // dropout-impact table of the target, need the lineage -> mutation lists of the target's dataset: warn once and skip the
        // target when it has none.
        // mutation_reasons: target name -> why that target's dataset has none.
        //
        def mutation_reasons = (mode == 'full' || mode == 'demix') && (params.dropout_impact || params.unexplained_mutations)
            ? targets.collectEntries { t -> [ (t.name): mutationTablesUnsupportedReason(t.barcodes ?: params.barcodes, t) ] }.findAll { _name, reason -> reason }
            : [:]
        if (mutation_reasons && is_panel) {
            log.warn "Skipping the dropout-impact and unexplained-mutation tables of panel target(s) " + mutation_reasons.collect { name, reason -> "${name} (${reason})" }.join('; ') + ". Use a Freyja dataset directory (with lineage_mutations.json) or a barcodes .csv to get them."
        }
        else if (mutation_reasons) {
            log.warn "Skipping the dropout-impact and unexplained-mutation tables: ${mutation_reasons.values()[0]}. Use a Freyja dataset directory (with lineage_mutations.json) or a barcodes .csv to get them."
        }
        def skip_mutation_tables = mutation_reasons.keySet().toList()
        def run_unexplained = (mode == 'full' || mode == 'demix') && params.unexplained_mutations && targets.any { t -> !(t.name in skip_mutation_tables) }
        // The per-target tables of the optional analyses ([ tkey, file ]; a target without one has no item)
        def ch_dropout_impact = channel.empty()
        def ch_dropout_summary = channel.empty()
        def ch_unexplained = channel.empty()
        def ch_timeseries = channel.empty()
        def ch_timeseries_lineages = channel.empty()
        def ch_growth_rates = channel.empty()
        def ch_growth_failed = channel.empty()
        def ch_growth_inputs = channel.empty()
        def ch_dashboard = channel.empty()
        def ch_lineage_summary = channel.empty()
        def ch_lineage_resolved = channel.empty()
        def ch_site_aggregates = channel.empty()
        def ch_site_aggregates_meta = channel.empty()
        if (run_unexplained) {
            def ch_unexplained_in = ch_view_variants
                .filter { meta, _variants -> meta.sample_type != 'negative' && !(meta.target.name in skip_mutation_tables) }
                .join(ch_demixed, remainder: true)
                .filter { meta, variants, _demix -> meta && variants }
                .map { meta, variants, demix -> [ meta, variants, demix ?: [] ] }
            def ch_to_unexplained = ch_unexplained_in
                .map { meta, variants, demix -> [ meta.target.name, meta, variants, demix ] }
                .combine(FREYJA_DATA.out.dataset_by_target, by: 0)
                .multiMap { _name, meta, variants, demix, dataset ->
                    sample: [ meta, variants, demix ]
                    dataset: dataset
                    gff: referenceGff(meta)
                    fasta: referenceFasta(meta)
                }
            UNEXPLAINED_MUTATIONS(ch_to_unexplained.sample, ch_to_unexplained.dataset, ch_to_unexplained.gff, ch_to_unexplained.fasta)
            // summary<sub>/unexplained_mutations.tsv: the tables of the samples of each target, with recurrence across them
            UNEXPLAINED_MUTATIONS_CONCAT(groupByTarget(UNEXPLAINED_MUTATIONS.out.tsv))
            ch_unexplained = UNEXPLAINED_MUTATIONS_CONCAT.out.tsv
        }

        //
        // Run-level steps, PER TARGET: they summarise the samples of one target, so each is run once per target with that target as
        // its meta (tkey: [ id: target name, target ]), on the items grouped by target (groupByTarget / groupSamples). A single-target
        // run is a panel of one: its one group, summary/ and pipeline_info/ are the same as ever. A target without any item has no
        // group and so no summary (an undetected target has none). Outputs go to summary<summary_sub>/ (conf/modules.config).
        //
        // The targets that have something to summarise: the target of a single-target run, the targets with an analysed sample in a panel
        def ch_target_keys = is_panel
            ? ch_target_metas.map { meta -> tkey(meta) }.unique()
            : channel.fromList(targets.collect { t -> tkey([ target: t ]) })

        // The samplesheet in effect (also for results-directory inputs) travels with the results, so a
        // later demix/aggregate run on this directory finds sample_type, collection_date and site
        ch_samplesheet.first().combine(ch_target_keys).subscribe { sheet, key ->
            summaryFile(outdir, key, 'samplesheet.csv', sheet.text)
        }

        // summary/barcode_version.json: sidecar of freyja_aggregated.tsv (PLAN 5.13), the dataset of each target
        if (mode == 'full' || mode == 'demix') {
            ch_target_keys.map { key -> [ key.id, key ] }
                .join(FREYJA_DATA.out.barcode_version_by_target)
                .subscribe { _name, key, barcode_version ->
                    summaryFile(outdir, key, 'barcode_version.json', barcode_version.text)
                }
        }

        // pipeline_info<sub>/freyja_data.json (and freyja_data_compare.json): the dataset block of each target that has something to
        // summarise. FREYJA_DATASET_INFO runs for every target of a panel, an undetected one included, and publishes nothing itself.
        ch_target_keys.map { key -> [ key.id, key ] }
            .join(FREYJA_DATA.out.freyja_data_by_target)
            .subscribe { _name, key, freyja_data ->
                pipelineInfoFile(outdir, key, 'freyja_data.json', freyja_data.text)
            }
        if (run_compare) {
            ch_target_keys.map { key -> [ key.id, key ] }
                .join(FREYJA_DATA_COMPARE.out.freyja_data_by_target)
                .subscribe { _name, key, freyja_data ->
                    pipelineInfoFile(outdir, key, 'freyja_data_compare.json', freyja_data.text)
                }
        }

        // summary/demix_status.tsv: one row per sample, whatever its status; samples that were found in
        // old results but could not be used (missing files, empty legacy .demix) are listed as `missing`
        def ch_demix_status = ch_status_lines
            .groupTuple()
            .map { key, lines -> [ key, summaryFile(outdir, key, 'demix_status.tsv', demixStatusText(lines)) ] }

        //
        // Aggregate + plot. Only samples with status ok have a .demix file (the rest are listed in
        // demix_status.tsv). Negative controls are demixed but kept out of the aggregate and the
        // plot unless --include_negatives_in_aggregate (PLAN D8).
        //
        def ch_aggregate_in = ch_demixed.filter { meta, _demix ->
            params.include_negatives_in_aggregate || meta.sample_type != 'negative'
        }

        if (!params.include_negatives_in_aggregate) {
            ch_demixed
                .filter { meta, _demix -> meta.sample_type == 'negative' }
                .subscribe { meta, _demix ->
                    log.info "Negative control ${meta.id}${is_panel ? " (${meta.target.name})" : ''} is not included in the aggregate and plot (use --include_negatives_in_aggregate to include it)"
                }
        }

        // Optional Freyja plot config (--freyja_plot_config): colours and lineage groups for freyja plot / dash / relgrowthrate
        def ch_plot_config = channel.value(plotConfigFile())

        // [ target name, tkey, sample IDs, their .demix files (same order), collection dates for `freyja plot --times` (format of Freyja's
        // data/times_metadata.csv; [] when no sample has a date) ] of each target, then the target's Freyja dataset in lockstep
        def ch_to_aggregate = groupSamples(ch_aggregate_in)
            .map { key, metas, demix_files -> [ key.id, key, metas.collect { meta -> meta.id }, demix_files, timesFile(metas) ] }
            .combine(FREYJA_DATA.out.dataset_by_target, by: 0)
            .multiMap { _name, key, ids, demix_files, times, dataset ->
                aggregate: [ key, ids, demix_files ]
                dataset: dataset
                times: times
            }
        FREYJA_AGGREGATE_PLOT(ch_to_aggregate.aggregate, ch_to_aggregate.dataset, ch_to_aggregate.times, ch_plot_config)

        // What the negative controls contain (they are not in the aggregate): same format, own table
        def ch_controls_aggregated = channel.empty()
        if (!params.include_negatives_in_aggregate) {
            FREYJA_AGGREGATE_CONTROLS(
                groupSamples(ch_demixed.filter { meta, _demix -> meta.sample_type == 'negative' })
                    .map { key, metas, demix_files -> [ key, metas.collect { meta -> meta.id }, demix_files ] }
            )
            ch_controls_aggregated = FREYJA_AGGREGATE_CONTROLS.out.aggregated
        }

        // Time-based analyses (time series, growth rates, dashboard) leave out low_confidence samples (PLAN D18)
        // unless --include_low_confidence: they get a copy of the aggregate without those rows, and the
        // two-dates test (below) does not count them. freyja_aggregated.tsv itself stays complete. In aggregate mode
        // only the samples carried over as low_confidence from --previous_results are affected (given .demix files
        // have no coverage information of their own).
        def filter_low_confidence = !params.include_low_confidence && (((mode == 'full' || mode == 'demix') && lowConfidenceBand()) || (mode == 'aggregate' && carried_low ? true : false))
        def ch_analysis_agg_by_target = FREYJA_AGGREGATE_PLOT.out.aggregated
        if (filter_low_confidence) {
            def ch_to_filter = FREYJA_AGGREGATE_PLOT.out.aggregated
                .join(ch_demix_status)
                .multiMap { key, agg, status ->
                    aggregated: [ key, agg ]
                    demix_status: status
                }
            FILTER_AGGREGATED(ch_to_filter.aggregated, ch_to_filter.demix_status)
            ch_analysis_agg_by_target = FILTER_AGGREGATED.out.aggregated
        }

        //
        // Provenance: run_info.json of each target
        //
        // (the per-sample library settings come from the metas of the target's samples, so the base info is made after they are known;
        // the target's primer BEDs and reference come with it)
        // (--compare_barcodes: the provenance block of the target's second dataset travels with the first one; [] without it)
        def ch_freyja_data_compare = run_compare
            ? FREYJA_DATA_COMPARE.out.freyja_data_by_target
            : channel.fromList(targets.collect { t -> [ t.name, [] ] })
        def ch_provenance_in = ch_target_metas
            .map { meta -> [ tkey(meta), meta ] }
            .groupTuple()
            .map { key, metas -> provenanceInputs(mode, key, metas) }
            .map { key, base_info, beds, fasta -> [ key.id, key, base_info, beds, fasta ] }
            .combine(FREYJA_DATA.out.freyja_data_by_target, by: 0)
            .combine(ch_freyja_data_compare, by: 0)
            .multiMap { _name, key, base_info, beds, fasta, freyja_data, freyja_data_compare ->
                base_info: [ key, base_info ]
                beds: beds
                fasta: fasta
                freyja_data: freyja_data
                freyja_data_compare: freyja_data_compare
            }
        PROVENANCE(
            ch_provenance_in.base_info,
            ch_samplesheet.first(),
            ch_provenance_in.beds,
            ch_provenance_in.fasta,
            ch_provenance_in.freyja_data,
            ch_collated_versions.first(),
            ch_provenance_in.freyja_data_compare,
        )
        ch_run_info = PROVENANCE.out.run_info.map { _key, run_info -> run_info }

        //
        // Optional analyses, PER TARGET, on the target's aggregate: dropout impact, site aggregation, lineage-group summary and its bootstrap
        // intervals, time series, growth rates, dashboard. Each needs the target's items to be present; a target without any
        // (or whose dataset cannot support the analysis) has none of these files.
        //
        def ch_aggregated = FREYJA_AGGREGATE_PLOT.out.aggregated   // [ tkey, freyja_aggregated.tsv ]

        // which dropped amplicons would have mattered for the lineages detected?
        // (amplicon_coverage.tsv exists only when some sample of the target has a primer scheme; capture samples are not in it)
        if ((mode == 'full' || mode == 'demix') && params.dropout_impact && targets.any { t -> !(t.name in skip_mutation_tables) }) {
            def ch_to_dropout = ch_amplicon_coverage
                .filter { key, _coverage -> !(key.id in skip_mutation_tables) }
                .join(ch_aggregated)
                .map { key, coverage, agg -> [ key.id, key, coverage, agg ] }
                .combine(FREYJA_DATA.out.dataset_by_target, by: 0)
                .multiMap { _name, key, coverage, agg, dataset ->
                    amplicons: [ key, coverage ]
                    aggregated: agg
                    dataset: dataset
                }
            DROPOUT_IMPACT(ch_to_dropout.amplicons, ch_to_dropout.aggregated, ch_to_dropout.dataset)
            ch_dropout_impact = DROPOUT_IMPACT.out.impact
            ch_dropout_summary = DROPOUT_IMPACT.out.summary
        }

        //
        // Site aggregation (weighted average of the per-site RESULTS, never pooled reads) and the lineage-group summary, on the
        // aggregate of each target (the sites CSV and the samplesheet are shared by the targets). demix_status.tsv of the target
        // always exists for them (header only when every sample is ok, e.g. aggregate mode).
        //
        def has_groups_summary = groups_targets.any { t -> !hierarchy_reasons.containsKey(t.name) }
        if (run_site_aggregation || has_groups_summary) {
            // aggregate mode: the given .demix files are the usable samples (status ok, except the low_confidence ones
            // carried over from the earlier run, which are already in ch_status_lines), but the published
            // demix_status.tsv lists only the unusable ones, so the summaries get their own, complete copy
            def ch_given_ok = mode == 'aggregate'
                ? ch_demixed.filter { meta, _demix -> !carried_low.containsKey(meta.key) }.map { meta, _demix -> [ tkey(meta), [ meta.id, 'ok', '', '', '', '' ].join('\t') ] }
                : channel.empty()
            // the samplesheet in effect plus the metadata of the samples discovery had to skip (a skipped negative control is
            // only known from there), as one CSV
            // (one row per sample: the sheet of panel results has a row per sample and target; a row of the sheet beats one of the skipped samples)
            def ch_summary_sheet_rows = ch_samplesheet
                .splitCsv(header: true)
                .map { row -> [ row.sample, [ row.sample, row.sample_type ?: 'sample', row.collection_date, row.site, row.viral_load ].collect { v -> csvCell(v) }.join(',') ] }
                .unique { row -> row[0] }
            def ch_summary_sheet = ch_summary_sheet_rows
                .join(ch_skipped.map { sample, _reason, sample_type, date, site, _target_name -> [ sample, [ sample, sample_type, date, site, '' ].collect { v -> csvCell(v) }.join(',') ] }.unique { row -> row[0] }, remainder: true)
                .map { _sample, sheet_line, skipped_line -> sheet_line ?: skipped_line }
                .collectFile(name: 'summary_samplesheet.csv', seed: 'sample,sample_type,collection_date,site,viral_load', sort: true, newLine: true)
                .first()
            // [ tkey, aggregate, complete demix_status.tsv ] of each target
            def ch_agg_status = ch_aggregated
                .join(ch_status_lines.mix(ch_given_ok).groupTuple().map { key, lines -> [ key, workTextFile('demix_status.tsv', demixStatusText(lines)) ] }, remainder: true)
                .filter { _key, agg, _status -> agg != null }
                .map { key, agg, status -> [ key, agg, status ?: workTextFile('demix_status.tsv', demixStatusText([])) ] }
            if (run_site_aggregation) {
                def ch_to_sites = ch_agg_status.multiMap { key, agg, status ->
                    aggregate: [ key, agg ]
                    demix_status: status
                }
                SITE_AGGREGATE(
                    ch_to_sites.aggregate,
                    ch_to_sites.demix_status,
                    ch_summary_sheet,
                    channel.value(file(params.site_aggregation, checkIfExists: true)),
                )
                ch_site_aggregates = SITE_AGGREGATE.out.aggregated
                ch_site_aggregates_meta = SITE_AGGREGATE.out.meta
            }
            if (has_groups_summary) {
                // [ target name, groups file, dataset ] of the targets whose dataset has the lineage hierarchy
                def ch_groups_dataset = ch_hierarchy
                    .filter { _name, _dataset, has -> has }
                    .map { name, dataset, _has -> [ name, groups_by_target[name], dataset ] }
                // [ tkey, aggregate, demix_status, site aggregates, their meta ] ([] [] without site aggregation)
                def ch_summary_base = run_site_aggregation
                    ? ch_agg_status.join(SITE_AGGREGATE.out.aggregated).join(SITE_AGGREGATE.out.meta)
                    : ch_agg_status.map { key, agg, status -> [ key, agg, status, [], [] ] }
                def ch_to_summary = ch_summary_base
                    .map { key, agg, status, sites, sites_meta -> [ key.id, key, agg, status, sites, sites_meta ] }
                    .combine(ch_groups_dataset, by: 0)
                    .multiMap { _name, key, agg, status, sites, sites_meta, groups_file, dataset ->
                        aggregate: [ key, agg ]
                        demix_status: status
                        groups: groups_file
                        dataset: dataset
                        sites: sites
                        sites_meta: sites_meta
                    }
                LINEAGE_SUMMARY(
                    ch_to_summary.aggregate,
                    ch_to_summary.demix_status,
                    ch_summary_sheet,
                    ch_to_summary.groups,
                    ch_to_summary.dataset,
                    ch_to_summary.sites,
                    ch_to_summary.sites_meta,
                )
                ch_lineage_summary = LINEAGE_SUMMARY.out.summary
                ch_lineage_resolved = LINEAGE_SUMMARY.out.resolved
                // --freyja_boot: 95 % intervals of the groups from the raw bootstrap draws (summary<sub>/lineage_group_intervals.tsv);
                // a target none of whose samples was bootstrapped has no raw draws and no intervals
                if (params.freyja_boot) {
                    def ch_to_intervals = ch_aggregated
                        .join(ch_boot_raw)
                        .map { key, agg, boots -> [ key.id, key, agg, boots ] }
                        .combine(ch_groups_dataset, by: 0)
                        .multiMap { _name, key, agg, boots, groups_file, dataset ->
                            aggregate: [ key, agg ]
                            boots: boots
                            groups: groups_file
                            dataset: dataset
                        }
                    BOOT_GROUP_INTERVALS(
                        ch_to_intervals.aggregate,
                        ch_to_intervals.boots,
                        ch_to_intervals.groups,
                        ch_to_intervals.dataset,
                        ch_summary_sheet,
                    )
                }
            }
        }

        // Time-based analyses leave out the low_confidence samples (see above): the two-dates test does not count them
        def ch_ts_samples = ch_aggregate_in
        if (filter_low_confidence && mode != 'aggregate') {
            ch_ts_samples = ch_aggregate_in
                .join(FREYJA_DEMIX.out.status)
                .filter { _meta, _demix, status_json -> new groovy.json.JsonSlurper().parse(status_json).status != 'low_confidence' }
                .map { meta, demix, _status_json -> [ meta, demix ] }
        }
        else if (filter_low_confidence) {
            ch_ts_samples = ch_aggregate_in.filter { meta, _demix -> !carried_low.containsKey(meta.key) }
        }

        // Time series: needs at least two distinct collection dates among the non-negative aggregated samples of the target
        // (growth rates and the dashboard check their own, stricter minimum and write a failure note instead)
        def ch_ts_go = ch_ts_samples
            .filter { meta, _demix -> meta.collection_date && meta.sample_type != 'negative' }
            .map { meta, _demix -> [ tkey(meta), meta.collection_date ] }
            .groupTuple()
            .filter { _key, dates -> dates.unique().size() >= 2 }
            .map { key, _dates -> [ key ] }
        def ch_ts_agg = ch_analysis_agg_by_target.join(ch_ts_go)
        def ch_metadata = ch_samplesheet.first()

        if (params.timeseries) {
            TIMESERIES_TABLE(ch_ts_agg, ch_metadata)
            ch_timeseries = TIMESERIES_TABLE.out.table
            // second, lineage-level table (every Freyja lineage instead of the summarized groups)
            TIMESERIES_TABLE_LINEAGES(ch_ts_agg, ch_metadata)
            ch_timeseries_lineages = TIMESERIES_TABLE_LINEAGES.out.table
        }
        if (params.freyja_growth || params.freyja_dash) {
            // the aggregate of each target with the target's Freyja dataset in lockstep
            def ch_to_growth = ch_analysis_agg_by_target
                .map { key, agg -> [ key.id, key, agg ] }
                .combine(FREYJA_DATA.out.dataset_by_target, by: 0)
                .multiMap { _name, key, agg, dataset ->
                    aggregate: [ key, agg ]
                    dataset: dataset
                }
            if (params.freyja_growth) {
                FREYJA_RELGROWTHRATE(ch_to_growth.aggregate, ch_metadata, ch_to_growth.dataset, ch_plot_config)
                ch_growth_rates = FREYJA_RELGROWTHRATE.out.growth_rates
                ch_growth_failed = FREYJA_RELGROWTHRATE.out.failed
                ch_growth_inputs = FREYJA_RELGROWTHRATE.out.inputs
            }
            if (params.freyja_dash) {
                FREYJA_DASH(ch_to_growth.aggregate, ch_metadata, ch_to_growth.dataset, params.freyja_dash_title, params.freyja_dash_intro, ch_plot_config)
                ch_dashboard = FREYJA_DASH.out.dashboard
            }
        }

        //
        // Report per target: report_data.json (everything the run produced for the target) -> report.html
        // The per-target tables are collected into one map per target (role -> file); [] = the target has none. The samplesheets are shared.
        //
        def ch_report_parts = taggedByTarget(PROVENANCE.out.run_info, 'run_info')
            .mix(taggedByTarget(ch_qc_summary, 'qc_summary'))
            .mix(taggedByTarget(ch_amplicon_coverage, 'amplicon_coverage'))
            .mix(taggedByTarget(FREYJA_AGGREGATE_PLOT.out.aggregated, 'aggregated'))
            .mix(taggedByTarget(ch_demix_status, 'demix_status'))
            .mix(taggedByTarget(ch_consensus_stats, 'consensus_stats'))
            .mix(taggedByTarget(ch_moi, 'moi'))
            .mix(taggedByTarget(ch_depth_files, 'depths'))
            .mix(taggedByTarget(ch_controls_aggregated, 'controls'))
            .mix(taggedByTarget(ch_contig_coverage, 'contig_coverage'))
            .mix(taggedByTarget(ch_consensus_contigs, 'consensus_contigs'))
            .mix(taggedByTarget(ch_compare_tsv, 'compare_tsv'))
            .mix(taggedByTarget(ch_compare_status, 'compare_status'))
            .mix(taggedByTarget(ch_boot_summary, 'boot_summary'))
            .mix(taggedByTarget(ch_dropout_impact, 'dropout_impact'))
            .mix(taggedByTarget(ch_dropout_summary, 'dropout_summary'))
            .mix(taggedByTarget(ch_unexplained, 'unexplained'))
            .mix(taggedByTarget(ch_timeseries, 'timeseries'))
            .mix(taggedByTarget(ch_timeseries_lineages, 'timeseries_lineages'))
            .mix(taggedByTarget(ch_growth_rates, 'growth_rates'))
            .mix(taggedByTarget(ch_growth_failed, 'growth_failed'))
            .mix(taggedByTarget(ch_growth_inputs, 'growth_inputs'))
            .mix(taggedByTarget(ch_dashboard.map { key, dashboard -> [ key, dashboard.name ] }, 'dashboard_name'))
            .mix(taggedByTarget(ch_lineage_summary, 'lineage_summary'))
            .mix(taggedByTarget(ch_lineage_resolved, 'lineage_resolved'))
            .mix(taggedByTarget(ch_site_aggregates, 'site_aggregates'))
            .mix(taggedByTarget(ch_site_aggregates_meta, 'site_aggregates_meta'))
            .groupTuple()
            .map { key, roles, files -> [ key, [ roles, files ].transpose().collectEntries { role, file_or_files -> [ (role): file_or_files ] } ] }
        // The sheet of the usable samples, plus the metadata of samples that discovery had to skip (a skipped
        // negative control is only known from there)
        def ch_skipped_sheet = ch_skipped
            .map { sample, _reason, sample_type, date, site, _target_name -> "${sample},${sample_type},${date},${site}" }
            .unique()
            .collectFile(name: 'skipped_samplesheet.csv', seed: 'sample,sample_type,collection_date,site', sort: true, newLine: true)
        def ch_report_sheets = ch_samplesheet.mix(ch_skipped_sheet).collect()
        def ch_to_report = ch_report_parts
            .combine(wrapItem(ch_report_sheets))
            .multiMap { key, parts, sheets ->
                data: [ key, parts.run_info, parts.qc_summary ?: [], parts.amplicon_coverage ?: [], parts.aggregated ?: [], parts.demix_status ?: [], parts.consensus_stats ?: [], parts.moi ?: [], parts.compare_tsv ?: [], parts.depths ?: [] ]
                extra: [ parts.controls ?: [], sheets, parts.boot_summary ?: [], parts.compare_status ?: [] ]
                analyses: [ parts.dropout_impact ?: [], parts.dropout_summary ?: [], parts.unexplained ?: [], parts.timeseries ?: [], parts.timeseries_lineages ?: [], parts.growth_rates ?: [], parts.growth_failed ?: [], parts.growth_inputs ?: [], parts.dashboard_name ?: '' ]
                groups: [ parts.lineage_summary ?: [], parts.lineage_resolved ?: [], parts.site_aggregates ?: [], parts.site_aggregates_meta ?: [] ]
                contigs: [ parts.contig_coverage ?: [], parts.consensus_contigs ?: [] ]
            }
        BUILD_REPORT_DATA(ch_to_report.data, ch_to_report.extra, ch_to_report.analyses, ch_to_report.groups, ch_to_report.contigs)
        RENDER_REPORT(BUILD_REPORT_DATA.out.json, file("${projectDir}/assets/templates/report.html.j2"))

        //
        // Panel overview (summary/overview.html): the sample x target matrix, once the report of every target exists (it links to
        // them; a target with no analysed sample has none). One run of a single process for the whole panel.
        //
        if (is_panel && !params.skip_report) {
            def ch_overview_items = ch_report_parts
                .join(RENDER_REPORT.out.html)
                .map { key, parts, _html -> [ key.id, parts.aggregated ?: [], parts.qc_summary ?: [], parts.demix_status ?: [], parts.run_info ?: [] ] }
                .toList()
            def ch_overview_in = PANEL_FRONT.out.targets
                .combine(wrapItem(ch_overview_items))   // (combine would spread the list of items)
                .map { targets_tsv, items -> [ targets_tsv ] + overviewInputs(targets.collect { t -> t.name }, items) }
                .multiMap { targets_tsv, manifest, files ->
                    targets_tsv: targets_tsv
                    manifest: manifest
                    files: files
                }
            BUILD_OVERVIEW(ch_overview_in.targets_tsv, ch_overview_in.manifest, ch_overview_in.files, file("${projectDir}/assets/templates/overview.html.j2"))
        }
    }
    }

    emit:
    run_info = ch_run_info // channel: path(run_info.json); empty for modes without provenance (report, update_db)
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// A file of <results>/summary<sub>/ (sub: '' or '/<target>' for panel results), or [] when the run did not write it
//
def optionalFile(results, name, sub = '') {
    def candidate = file("${results}/summary${sub}/${name}")
    return candidate.exists() ? candidate : []
}

//
// --mode report on panel results: the inputs of BUILD_REPORT_DATA for one target, as the tables the run wrote to summary/<target>/ and its
// depths files (samples/*/<target>/freyja/) as the same five tuples the single-target report takes, [ data, extra, analyses, groups, contigs ]
//
def reportBuildInputs(results, target) {
    def sub = target.summary_sub
    def depths = file("${results}/samples/*/${target.name}/freyja/*.depths").sort { a, b -> a.name <=> b.name }
    return [
        data: [
            tkey([ target: target ]),
            file("${results}/pipeline_info${sub}/run_info.json"),
            optionalFile(results, 'qc_summary.tsv', sub),
            optionalFile(results, 'amplicon_coverage.tsv', sub),
            optionalFile(results, 'freyja_aggregated.tsv', sub),
            optionalFile(results, 'demix_status.tsv', sub),
            optionalFile(results, 'consensus_stats.tsv', sub),
            optionalFile(results, 'mutations_of_interest.tsv', sub),
            optionalFile(results, 'barcode_compare.tsv', sub),
            depths,
        ],
        extra: [
            optionalFile(results, 'freyja_aggregated_controls.tsv', sub),
            optionalFile(results, 'samplesheet.csv', sub),
            optionalFile(results, 'freyja_boot_summary.tsv', sub),
            optionalFile(results, 'barcode_compare_status.tsv', sub),
        ],
        analyses: [
            optionalFile(results, 'dropout_impact.tsv', sub),
            optionalFile(results, 'dropout_impact_summary.tsv', sub),
            optionalFile(results, 'unexplained_mutations.tsv', sub),
            optionalFile(results, 'freyja_timeseries.tsv', sub),
            optionalFile(results, 'freyja_timeseries_lineages.tsv', sub),
            optionalFile(results, 'freyja_growth_rates.csv', sub),
            optionalFile(results, 'freyja_growth_rates.failed.txt', sub),
            optionalFile(results, 'freyja_growth_inputs.json', sub),
            file("${results}/summary${sub}/freyja_dashboard.html").exists() ? 'freyja_dashboard.html' : '',
        ],
        groups: [
            optionalFile(results, 'lineage_summary.csv', sub),
            optionalFile(results, 'lineage_groups_resolved.tsv', sub),
            optionalFile(results, 'site_aggregates_freyja.tsv', sub),
            optionalFile(results, 'site_aggregates_meta.tsv', sub),
        ],
        contigs: [
            optionalFile(results, 'contig_coverage.tsv', sub),
            optionalFile(results, 'consensus_stats_contigs.tsv', sub),
        ],
    ]
}

//
// The reference FASTA of the target of a sample (a per-sample channel item, see multiMap), and its annotation ([] = none)
//
def referenceFasta(meta) {
    return meta.target.fasta ? file(meta.target.fasta, checkIfExists: true) : []
}

def referenceGff(meta) {
    return meta.target.gff ? file(meta.target.gff, checkIfExists: true) : []
}

//
// Is the low-confidence band (PLAN D18) switched on? It is off when the floor is 0 or not below the
// high-confidence threshold: samples are then either demixed normally or skipped, as before D18.
//
def lowConfidenceBand() {
    return params.demix_floor_coverage > 0 && params.demix_floor_coverage < params.min_genome_coverage_for_demix
}

//
// Aggregate mode, --previous_results: [ key: [ reason, coverage ] ] of the samples that summary/demix_status.tsv of an earlier run lists
// as low_confidence (PLAN D18); the key is meta.key (the sample ID, or ID__target for the targets of panel results, whose status
// files are summary/<target>/demix_status.tsv). Other statuses, other inputs (--input, --demix_dir) and directories without that
// file give nothing: their .demix files count as ok.
//
def carriedLowConfidence(targets) {
    def carried = [:]
    splitPaths(params.previous_results).each { dir ->
        def status_files = resultsAreMultiTarget(dir)
            ? targets.collect { t -> [ "__${t.name}", file("${dir}/summary${t.summary_sub}/demix_status.tsv") ] }
            : [ [ '', file("${dir}/summary/demix_status.tsv") ] ]
        status_files.each { key_suffix, status_file ->
            if (status_file.exists()) {
                def lines = status_file.readLines().findAll { line -> line.trim() }
                def header = lines ? lines[0].split('\t', -1).collect { column -> column.trim() } : []
                lines.drop(1).each { line ->
                    def cells = line.split('\t', -1)
                    def row = [:]
                    header.eachWithIndex { column, i -> row[column] = i < cells.size() ? cells[i].trim() : '' }
                    if (row.sample && row.status == 'low_confidence') {
                        carried["${row.sample}${key_suffix}".toString()] = [ reason: row.reason ?: '', coverage: row.genome_coverage_pct ?: '' ]
                    }
                }
            }
        }
    }
    return carried
}

//
// One CSV cell, always quoted (commas and quotes inside are safe)
//
def csvCell(value) {
    return '"' + (value ?: '').toString().replaceAll('[\r\n]+', ' ').replace('"', '""') + '"'
}

//
// Single-item channel whose item may be [] or a list of files, wrapped so that .combine() keeps
// it as ONE tuple element instead of flattening it
//
def wrapItem(ch) {
    return ch.map { item -> [ item ] }
}

//
// Facts known to Nextflow, as JSON for bin/make_run_info.py to extend with checksums,
// tool versions and the Freyja dataset block
//
def runInfoBase(mode, metas, target) {
    def git_commit = workflow.commitId
    if (!git_commit) {
        try {
            def proc = ['git', '-C', projectDir.toString(), 'rev-parse', 'HEAD'].execute()
            def out = proc.text.trim()
            git_commit = proc.waitFor() == 0 && out ? out : null
        }
        catch (Exception _e) {
            git_commit = null
        }
    }
    // pipelines_testdata_base_path is nf-core template plumbing (a URL), not a setting of this pipeline
    // (the per-pathogen facts of conf/pathogens/*.config follow from --pathogen and are not recorded either)
    def plain_params = params.findAll { key, value -> !(key.toString() in ['pipelines_testdata_base_path', 'lineage_groups_bundled', 'freyja_container_pathogen', 'primer_schemes_bundled', 'primer_scheme_default', 'freyja_contig']) && !(key.toString() == 'whole_genome' && !value) }.collectEntries { key, value ->
        [ (key.toString()): (value == null || value instanceof Boolean || value instanceof Number) ? value : value.toString() ]
    }
    // the effective primer scheme is recorded even when the default (null = the default scheme, D14) applies
    // (a panel target's own primer_scheme from the panel CSV counts as its default)
    if (mode == 'full' && protocolSettings(params.protocol).primer_trimming && !params.primer_bed && !params.primer_scheme && (target.primer_scheme ?: defaultPrimerScheme(target))) {
        plain_params.primer_scheme = target.primer_scheme ?: defaultPrimerScheme(target)
    }
    // lineage groups: which file made lineage_summary.csv (only when the summary is made for this target)
    def groups_info = [ lineage_groups_file: null, lineage_groups_source: null, lineage_groups_sha256: null ]
    if (params.lineage_summary && lineageGroupsFile(target) && mode in ['full', 'demix', 'aggregate'] && !lineageHierarchyUnsupportedReason(target)) {
        groups_info = lineageGroupsInfo(target)
    }
    def info = [
        pipeline_name: workflow.manifest.name,
        pipeline_version: workflow.manifest.version,
        git_commit: git_commit,
        nextflow_version: workflow.nextflow.version.toString(),
        mode: mode,
        target: target.name,
        profile: workflow.profile,
        command_line: workflow.commandLine,
        run_name: workflow.runName,
        start: workflow.start.toString(),
        params: plain_params,
        // earlier runs whose results this run builds on (demix/aggregate); barcode_version 'unknown'
        // for legacy Snakemake results, which recorded none
        previous_runs: previousRuns(mode, target),
    ] + libraryInfo(mode, metas) + groups_info
    if (target.panel) {
        info.panel = params.panel ? params.panel.toString() : 'from results'   // not given with --panel: the targets are those of the panel results / sheet read
    }
    // several reference contigs (PLAN D25): the contigs in FASTA order, and the one Freyja demixes
    if (mode in ['full', 'demix'] && multiContig(target)) {
        info.contigs = fastaRecords(file(target.fasta)).collect { record -> [ name: record.name, length: record.length ] }
        info.freyja_contig = target.freyja_contig
    }
    if (mode == 'aggregate') {
        info.freyja_data_note = 'aggregate mode demixes nothing: freyja_data describes the dataset whose lineage hierarchy groups the plot, not the barcodes the .demix files were made with (see previous_runs)'
    }
    return groovy.json.JsonOutput.prettyPrint(groovy.json.JsonOutput.toJson(info))
}

//
// Per-target grouping of the run-level steps. tkey(meta) = [ id: target name, target ] is the key (the meta of the run-level
// processes); a target without any item has no group, so no process runs for it.
//

// [ meta, file ] -> [ tkey, [ files sorted by name ] ] (the sorting keeps the task inputs, and so -resume, reproducible)
def groupByTarget(ch) {
    return ch
        .map { meta, f -> [ tkey(meta), f ] }
        .groupTuple()
        .map { key, files -> [ key, files.sort(false) { a, b -> a.name <=> b.name } ] }
}

// [ meta, file ] -> [ tkey, [ metas sorted by sample ID ], [ files in the same order ] ]
def groupSamples(ch) {
    return ch
        .map { meta, f -> [ tkey(meta), meta, f ] }
        .groupTuple()
        .map { key, metas, files ->
            def pairs = [ metas, files ].transpose().sort(false) { a, b -> a[0].id <=> b[0].id }
            [ key, pairs.collect { pair -> pair[0] }, pairs.collect { pair -> pair[1] } ]
        }
}

// [ tkey, item ] -> [ tkey, role, item ]: to collect the optional tables of a target in one map (role -> item)
def taggedByTarget(ch, role) {
    return ch.map { key, item -> [ key, role, item ] }
}

//
// The inputs of BUILD_OVERVIEW: the manifest (tab separated: target, role, path as staged under in/<n>/) and the list of files it names. target_names =
// every target of the run (a column of the matrix, also when it has no report); items = [ name, aggregate, qc table, demix status, run_info.json ]
// of the targets that have a report ([] for a file the target lacks). Sorted by target so the staging (and a -resume) does not depend on arrival order.
//
def overviewInputs(target_names, items) {
    def files = []
    def lines = target_names.sort(false).collect { name -> "${name}\ttarget\t-".toString() }
    items.sort(false) { a, b -> a[0] <=> b[0] }.each { item ->
        def name = item[0]
        [ 'aggregated', 'qc', 'demix_status', 'run_info' ].eachWithIndex { role, i ->
            def f = item[i + 1]
            if (f) {
                files << f
                lines << "${name}\t${role}\tin/${files.size()}/${f.name}".toString()
            }
        }
        lines << "${name}\treport\t-".toString()
    }
    return [ workTextFile('overview_inputs.tsv', lines.join('\n') + '\n'), files ]
}

//
// Small text files assembled in Groovy from the items of a target group (collectFile cannot write below a sub-directory, nor keep the
// target of its result). They are only rewritten when their content changes, so a -resume does not see a new modification time.
//

// summary/demix_status.tsv text
def demixStatusHeader() {
    return [ 'sample', 'status', 'reason', 'genome_coverage_pct', 'barcode_version', 'barcode_version_iso' ].join('\t')
}

def demixStatusText(lines) {
    return ([ demixStatusHeader() ] + lines.sort(false)).join('\n') + '\n'
}

// `freyja plot --times` input (format of Freyja's data/times_metadata.csv) of the samples of a group that have a collection date, or [] when none has
def timesFile(metas) {
    def lines = metas.findAll { meta -> meta.collection_date }.collect { meta -> "${meta.id},${meta.collection_date}".toString() }.sort(false)
    return lines ? workTextFile('times.csv', ([ 'Sample,sample_collection_datetime' ] + lines).join('\n') + '\n') : []
}

// The tables of one target concatenated: the header of the first file, then the rows (everything but the header line) of every file
def concatKeepingHeader(files) {
    def tables = files.collect { f -> f.readLines() }
    def first = tables.find { lines -> lines }
    def rows = tables.collectMany { lines -> lines.drop(1) }
    return ((first ? [ first[0] ] : []) + rows).join('\n') + '\n'
}

// <outdir>/summary<summary_sub>/<name> of a target (the file a published run-level step would write there)
def summaryFile(outdir, key, name, text) {
    return targetOutputFile(outdir, 'summary', key, name, text)
}

// <outdir>/pipeline_info<summary_sub>/<name> of a target
def pipelineInfoFile(outdir, key, name, text) {
    return targetOutputFile(outdir, 'pipeline_info', key, name, text)
}

def targetOutputFile(outdir, folder, key, name, text) {
    def dir = file("${outdir}/${folder}${key.target.summary_sub}")
    def f = dir.resolve(name)
    if (!f.exists() || f.text != text) {
        dir.mkdirs()
        f.text = text
    }
    return f
}

// A file with the given name and text under the work directory (in a directory named by a digest of both, so equal files are shared)
def workTextFile(name, text) {
    def digest = java.security.MessageDigest.getInstance('SHA-256').digest("${name}\n${text}".getBytes('UTF-8')).collect { b -> String.format('%02x', b) }.join('').substring(0, 16)
    def f = file("${workflow.workDir}/wasteline_run/${digest}/${name}")
    if (!f.exists() || f.text != text) {
        f.parent.mkdirs()
        f.text = text
    }
    return f
}

//
// What PROVENANCE of one target needs besides the shared inputs: [ tkey, base info file, the primer BEDs of the target's samples, the
// target's reference FASTA ([] in aggregate mode, which maps nothing) ]. metas are those of THIS target's samples only.
//
def provenanceInputs(mode, key, metas) {
    def target = key.target
    if (metas.any { meta -> meta.target.name != target.name }) {
        error("Samples of several targets in the provenance of target ${target.name}: ${metas.collect { meta -> meta.target.name }.unique().join(', ')}")
    }
    def beds = metas.findAll { meta -> meta.primer_scheme }.collect { meta -> primerBedFor(meta, meta.target) }.unique().sort(false) { a, b -> a.name <=> b.name }
    def fasta = mode != 'aggregate' && target.fasta ? file(target.fasta, checkIfExists: true) : []
    return [ key, workTextFile('run_info_base.json', runInfoBase(mode, metas, target)), beds, fasta ]
}
