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
include { PROVENANCE             } from '../modules/local/provenance/main'
include { FREYJA_DATA            } from '../subworkflows/local/freyja_data/main'
include { FREYJA_DATA as FREYJA_DATA_COMPARE } from '../subworkflows/local/freyja_data/main'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { primerBedFor           } from '../subworkflows/local/utils_nfcore_wasteline_pipeline'
include { libraryInfo            } from '../subworkflows/local/utils_nfcore_wasteline_pipeline'
include { protocolSettings       } from '../subworkflows/local/utils_nfcore_wasteline_pipeline'
include { defaultPrimerScheme    } from '../subworkflows/local/utils_nfcore_wasteline_pipeline'
include { mutationTablesUnsupportedReason } from '../subworkflows/local/utils_nfcore_wasteline_pipeline'
include { lineageGroupsFile      } from '../subworkflows/local/utils_nfcore_wasteline_pipeline'
include { lineageGroupsInfo      } from '../subworkflows/local/utils_nfcore_wasteline_pipeline'
include { lineageHierarchyUnsupportedReason } from '../subworkflows/local/utils_nfcore_wasteline_pipeline'
include { hierarchyMayBeMissing } from '../subworkflows/local/utils_nfcore_wasteline_pipeline'
include { bundledSchemesApply } from '../subworkflows/local/utils_nfcore_wasteline_pipeline'
include { previousRuns           } from '../subworkflows/local/utils_nfcore_wasteline_pipeline'
include { splitPaths             } from '../subworkflows/local/utils_nfcore_wasteline_pipeline'

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
    ch_skipped     // channel: [ sample, reason, sample_type, collection_date, site ] of samples that were found in results directories but not usable
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
        def existing_data = file("${results}/summary/report_data.json")
        // report_data.json is rebuilt from the summary files when the report must be anonymised: the renderer has no
        // anonymisation of its own, and an existing report_data.json holds the real IDs
        if (existing_data.exists() && params.report_anonymize == 'none') {
            // the run already wrote report_data.json: only the HTML is rendered again
            RENDER_REPORT(channel.value(existing_data), file("${projectDir}/assets/templates/report.html.j2"))
        }
        else {
            // rebuild report_data.json from the files the run left behind
            def depths = (file("${results}/samples/*/freyja/*.depths") + file("${results}/*/freyja/*.depths")).sort { a, b -> a.name <=> b.name }
            BUILD_REPORT_DATA(channel.value([
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
            ]))
            RENDER_REPORT(BUILD_REPORT_DATA.out.json, file("${projectDir}/assets/templates/report.html.j2"))
        }
    }
    else {
    //
    // Freyja dataset: resolved once, shared by every sample. update_db always fetches the latest.
    //
    FREYJA_DATA(mode == 'update_db' ? 'latest' : params.barcodes, params.pathogen, '--barcodes')

    if (mode == 'update_db') {
        // Nothing else to do: the dataset is now in --db_dir/freyja/<barcode-version>/
        FREYJA_DATA.out.freyja_data.subscribe { freyja_data ->
            def info = new groovy.json.JsonSlurper().parse(freyja_data)
            log.info "\nFreyja data updated: barcode version ${info.barcode_version} (${info.barcode_version_iso}), stored in ${info.dataset_path}\nUse it with: --barcodes ${info.barcode_version_iso.take(10)}   (or --barcodes ${info.barcode_version_iso} for this exact update)"
        }
    }
    else {
        // --fasta is required by modes full and demix (validateInputParameters); --gff is optional for pathogens without a bundled annotation
        def ch_fasta = channel.value(params.fasta ? file(params.fasta, checkIfExists: true) : [])
        def ch_gff = channel.value(params.gff ? file(params.gff, checkIfExists: true) : [])

        //
        // The lineage groups file (bundled default or --lineage_groups) and the sites CSV are checked FIRST, in parallel with fastp and the mapping, so that a
        // typo stops the run within minutes (these processes keep the default error strategy). Unknown lineages and
        // unknown site names are warnings, shown here. A dataset without a lineage hierarchy (custom barcodes file) cannot
        // resolve the groups: the file format is still checked, the summary itself is skipped with one warning.
        //
        // A pathogen without bundled lineage groups (conf/pathogens/*.config) has the summary off unless --lineage_groups is given (PLAN D22).
        def run_lineage_groups = params.lineage_summary && lineageGroupsFile() ? true : false
        def run_site_aggregation = params.site_aggregation as boolean
        if (params.lineage_summary && !run_lineage_groups) {
            log.info "Lineage-group summary off: no lineage groups are bundled for --pathogen ${params.pathogen}. Give --lineage_groups <CSV> (and a Freyja dataset with a lineage hierarchy) to get summary/lineage_summary.csv."
        }
        def lineage_reason = run_lineage_groups ? lineageHierarchyUnsupportedReason() : null
        // The dataset the hierarchy-dependent steps read. A dataset downloaded during the run for a pathogen that may ship
        // no hierarchy (RSV, dengue, ...) is looked at first: without lineages.yml those steps are skipped, with ONE warning.
        def ch_hierarchy_dataset = FREYJA_DATA.out.dataset
        if (run_lineage_groups && !lineage_reason && hierarchyMayBeMissing()) {
            def ch_has_hierarchy = FREYJA_DATA.out.freyja_data.map { freyja_data -> new groovy.json.JsonSlurper().parse(freyja_data).lineages_yml_sha256 ? true : false }
            ch_hierarchy_dataset = FREYJA_DATA.out.dataset.combine(ch_has_hierarchy)
                .filter { _dataset, has -> has }
                .map { dataset, _has -> dataset }
            ch_has_hierarchy.filter { has -> !has }.subscribe { _has ->
                log.warn "Skipping the lineage-group summary: the Freyja dataset downloaded for --pathogen ${params.pathogen} has no lineage hierarchy (lineages.yml), which the groups need."
            }
        }
        if (run_lineage_groups) {
            LINEAGE_GROUPS_CHECK(
                channel.value(lineageGroupsFile()),
                lineage_reason ? channel.value([]) : ch_hierarchy_dataset.ifEmpty([]),
            )
            LINEAGE_GROUPS_CHECK.out.report.subscribe { report ->
                report.readLines().findAll { line -> line.contains('WARNING') }.each { line -> log.warn line }
            }
            if (lineage_reason) {
                log.warn "Skipping the lineage-group summary: ${lineage_reason}. Use a Freyja dataset directory with lineages.yml (--barcodes latest, container or a dataset directory) to get summary/lineage_summary.csv."
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
            //
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
            def ch_mapped = BAM_MAPPED_READS.out.counted
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
                    bed: primerBedFor(meta)
                }
            IVAR_TRIM(ch_to_trim.reads, ch_to_trim.bed)
            ch_ivar_log = IVAR_TRIM.out.log

            // Does the primer scheme match the reads? Needs the mapped BAM BEFORE trimming (trimming moves
            // the 5' ends). Declared scheme = the sample's own (meta.primer_scheme; 'custom' = its BED).
            if (params.primer_scheme_check && bundledSchemesApply()) {
                READ_STARTS(ch_mapped.mapped
                    .filter { meta, _bam, _bai, _count -> meta.primer_trimming }
                    .map { meta, bam, _bai, _count -> [ meta, bam ] })
                def ch_to_check = READ_STARTS.out.tsv.multiMap { meta, tsv ->
                    starts: [ meta, tsv ]
                    declared: meta.primer_scheme
                    custom_bed: meta.primer_bed ? file(meta.primer_bed, checkIfExists: true) : []
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
                SAMTOOLS_CRAM(ch_final_bam, ch_fasta, [])
            }

            FREYJA_VARIANTS(ch_final_bam, ch_fasta)

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

        //
        // Demix per sample (never fails the run) -> [ meta, demix ]; aggregate mode takes given .demix files
        //
        def ch_demixed = channel.empty()
        def ch_compare_tsv = channel.empty()
        def ch_compare_status = channel.empty()
        def ch_boot_summary = channel.empty()
        def ch_boot_raw = channel.empty()
        // aggregate mode with --previous_results: samples that were low_confidence in the earlier run keep that status (PLAN D18)
        def carried_low = mode == 'aggregate' ? carriedLowConfidence() : [:]
        // demix_status.tsv lines: [ sample, status, reason, genome_coverage_pct, barcode_version, barcode_version_iso ]
        def ch_status_lines = ch_skipped.map { sample, reason, _type, _date, _site -> [ sample, 'missing', reason, '', '', '' ].collect { value -> value.toString().replaceAll('[\t\r\n]+', ' ') }.join('\t') }

        if (mode == 'full' || mode == 'demix') {
            FREYJA_DEMIX(
                ch_freyja_in,
                FREYJA_DATA.out.dataset,
                FREYJA_DATA.out.barcode_version,
            )

            ch_status_lines = ch_status_lines.mix(FREYJA_DEMIX.out.status
                .map { _meta, status_json ->
                    def rec = new groovy.json.JsonSlurper().parse(status_json)
                    [ rec.sample, rec.status, rec.reason, rec.genome_coverage_pct, rec.barcode_version, rec.barcode_version_iso ]
                        .collect { value -> value == null ? '' : value.toString().replaceAll('[\t\r\n]+', ' ') }
                        .join('\t')
                })

            // summary/barcode_version.json: sidecar of freyja_aggregated.tsv (PLAN 5.13)
            FREYJA_DATA.out.barcode_version.collectFile(name: 'barcode_version.json', storeDir: "${outdir}/summary")

            ch_demixed = FREYJA_DEMIX.out.demix

            //
            // --compare_barcodes: the same variants/depths demixed with a SECOND dataset (resolved
            // independently, recorded as freyja_data_compare); abundances of both side by side
            //
            if (params.compare_barcodes) {
                FREYJA_DATA_COMPARE(params.compare_barcodes, params.pathogen, '--compare_barcodes')
                FREYJA_DEMIX_COMPARE(
                    ch_freyja_in,
                    FREYJA_DATA_COMPARE.out.dataset,
                    FREYJA_DATA_COMPARE.out.barcode_version,
                )
                BARCODE_COMPARE(
                    FREYJA_DEMIX.out.demix.join(FREYJA_DEMIX_COMPARE.out.demix),
                    FREYJA_DATA.out.barcode_version,
                    FREYJA_DATA_COMPARE.out.barcode_version,
                )
                ch_compare_tsv = BARCODE_COMPARE.out.rows
                    .map { _meta, rows -> rows }
                    .toSortedList { a, b -> a.name <=> b.name }
                    .flatMap { files -> files }
                    .collectFile(name: 'barcode_compare.tsv', storeDir: "${outdir}/summary", sort: { file -> file.name }, seed: 'sample\tlineage\tabundance_a\tabundance_b\tbarcode_version_a\tbarcode_version_b\n', newLine: false)
                // outcome of the second demix per sample (a failed one is flagged in the report, never fatal)
                ch_compare_status = FREYJA_DEMIX_COMPARE.out.status
                    .map { _meta, status_json ->
                        def rec = new groovy.json.JsonSlurper().parse(status_json)
                        [ rec.sample, rec.status, rec.reason, rec.barcode_version ]
                            .collect { value -> value == null ? '' : value.toString().replaceAll('[\t\r\n]+', ' ') }
                            .join('\t')
                    }
                    .collectFile(name: 'barcode_compare_status.tsv', storeDir: "${outdir}/summary", seed: 'sample\tstatus\treason\tbarcode_version', sort: true, newLine: true)
            }

            //
            // --freyja_boot: bootstrap intervals for every sample that demixed
            //
            if (params.freyja_boot) {
                FREYJA_BOOT(ch_freyja_in.join(FREYJA_DEMIX.out.demix), FREYJA_DATA.out.dataset)
                ch_boot_raw = FREYJA_BOOT.out.raw.map { _meta, raw -> raw }
                ch_boot_summary = FREYJA_BOOT.out.summary
                    .map { _meta, summary -> summary }
                    .toSortedList { a, b -> a.name <=> b.name }
                    .flatMap { files -> files }
                    .collectFile(name: 'freyja_boot_summary.tsv', storeDir: "${outdir}/summary", sort: { file -> file.name }, seed: 'sample\tlineage\tabundance\tlo\thi\n', newLine: false)
            }
        }
        else {
            ch_demixed = ch_rows
            ch_status_lines = ch_status_lines.mix(ch_demixed
                .filter { meta, _demix -> carried_low.containsKey(meta.id) }
                .map { meta, _demix ->
                    [ meta.id, 'low_confidence', carried_low[meta.id].reason, carried_low[meta.id].coverage, '', '' ].join('\t')
                })
        }

        // summary/demix_status.tsv: one row per sample, whatever its status; samples that were found in
        // old results but could not be used (missing files, empty legacy .demix) are listed as `missing`
        def status_header = ['sample', 'status', 'reason', 'genome_coverage_pct', 'barcode_version', 'barcode_version_iso'].join('\t')
        def ch_demix_status = ch_status_lines
            .collectFile(name: 'demix_status.tsv', storeDir: "${outdir}/summary", seed: status_header, sort: true, newLine: true)

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
                    log.info "Negative control ${meta.id} is not included in the aggregate and plot (use --include_negatives_in_aggregate to include it)"
                }
        }

        def ch_aggregate_files = ch_aggregate_in
            .toSortedList { a, b -> a[0].id <=> b[0].id }
            .filter { entries -> entries }
            .map { entries -> [ entries.collect { entry -> entry[0].id }, entries.collect { entry -> entry[1] } ] }

        // Collection dates for `freyja plot --times` (format of Freyja's data/times_metadata.csv)
        def ch_times = ch_aggregate_in
            .filter { meta, _demix -> meta.collection_date }
            .map { meta, _demix -> "${meta.id},${meta.collection_date}" }
            .collectFile(name: 'times.csv', seed: 'Sample,sample_collection_datetime', sort: true, newLine: true)
            .ifEmpty([])

        FREYJA_AGGREGATE_PLOT(ch_aggregate_files, FREYJA_DATA.out.dataset, ch_times)

        // What the negative controls contain (they are not in the aggregate): same format, own table
        def ch_controls_aggregated = channel.empty()
        if (!params.include_negatives_in_aggregate) {
            def ch_control_files = ch_demixed
                .filter { meta, _demix -> meta.sample_type == 'negative' }
                .toSortedList { a, b -> a[0].id <=> b[0].id }
                .filter { entries -> entries }
                .map { entries -> [ entries.collect { entry -> entry[0].id }, entries.collect { entry -> entry[1] } ] }
            FREYJA_AGGREGATE_CONTROLS(ch_control_files)
            ch_controls_aggregated = FREYJA_AGGREGATE_CONTROLS.out.aggregated
        }

        //
        // Coverage QC (every sample, negatives included) from the Freyja depths file, plus flagstat
        // (mode full), the fastp JSON (when known) and per-amplicon coverage
        //
        def ch_qc_summary = channel.value([])
        def ch_amplicon_coverage = channel.value([])
        def ch_consensus_stats = channel.value([])
        def ch_moi = channel.value([])
        def ch_depth_files = channel.value([])

        if (mode == 'full' || mode == 'demix') {
            // only the samples with a primer scheme (meta.primer_scheme), each against its own BED
            def ch_to_cover = ch_depths
                .filter { meta, _depths -> meta.primer_scheme }
                .multiMap { meta, depths ->
                    depths: [ meta, depths ]
                    bed: primerBedFor(meta)
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

            QC_CONCAT(
                QC_SAMPLE.out.tsv.map { _meta, tsv -> tsv }.collect(),
                ch_amplicons.map { _meta, tsv -> tsv }.collect().ifEmpty([]),
            )
            ch_qc_summary = QC_CONCAT.out.qc_summary.ifEmpty([])
            ch_amplicon_coverage = QC_CONCAT.out.amplicon_coverage.ifEmpty([])
            ch_depth_files = ch_depths.map { _meta, depths -> depths }.collect()

            //
            // Consensus genomes from the Freyja variants + depths (negative controls get none).
            // A sample whose consensus steps fail still gets a row, with status failed.
            //
            if (!params.skip_consensus) {
                def ch_consensus_samples = ch_variants.filter { meta, _variants -> meta.sample_type != 'negative' }
                IVAR_VARIANTS_TO_VCF(ch_consensus_samples, ch_fasta)
                CONSENSUS_FILTER_VCF(IVAR_VARIANTS_TO_VCF.out.vcf)
                LOW_DEPTH_MASK(ch_depths.join(CONSENSUS_FILTER_VCF.out.bed), ch_fasta)
                CONSENSUS_BCFTOOLS(CONSENSUS_FILTER_VCF.out.vcf.join(LOW_DEPTH_MASK.out.bed), ch_fasta)
                CONSENSUS_STATS(CONSENSUS_BCFTOOLS.out.fasta.join(CONSENSUS_BCFTOOLS.out.status), ch_fasta)

                def ch_stats = ch_consensus_samples
                    .join(CONSENSUS_STATS.out.tsv, remainder: true)
                    .branch { _meta, _variants, stats ->
                        built: stats
                        failed: true
                    }
                def ch_failed_stats = ch_stats.failed
                    .collectFile { meta, _variants, _stats ->
                        [ "${meta.id}.consensus_stats.tsv", "sample\tlength\tn_count\tn_pct\tambiguous_count\tstatus\treason\n${meta.id}\tNA\tNA\tNA\tNA\tfailed\tconsensus steps did not complete\n" ]
                    }
                CONSENSUS_STATS_CONCAT(ch_stats.built.map { _meta, _variants, stats -> stats }.mix(ch_failed_stats).collect())
                CONSENSUS_CONCAT(CONSENSUS_BCFTOOLS.out.fasta.map { _meta, fasta -> fasta }.collect())
                ch_consensus_stats = CONSENSUS_STATS_CONCAT.out.tsv.ifEmpty([])

                // --annotate_variants: amino-acid consequences of the variants applied to the consensus
                if (params.annotate_variants) {
                    ANNOTATE_VARIANTS(CONSENSUS_FILTER_VCF.out.vcf, ch_fasta, ch_gff)
                    ANNOTATE_VARIANTS.out.tsv
                        .map { _meta, tsv -> tsv }
                        .toSortedList { a, b -> a.name <=> b.name }
                        .flatMap { files -> files }
                        .collectFile(name: 'annotated_variants.tsv', storeDir: "${outdir}/summary", keepHeader: true, skip: 1, sort: { file -> file.name })
                }
            }

            //
            // Mutations of interest (opt-in)
            //
            if (params.mutations_of_interest) {
                MUTATIONS_OF_INTEREST(ch_variants.join(ch_depths), file(params.mutations_of_interest, checkIfExists: true), ch_gff, ch_fasta)
                MUTATIONS_OF_INTEREST_CONCAT(MUTATIONS_OF_INTEREST.out.tsv.map { _meta, tsv -> tsv }.collect())
                ch_moi = MUTATIONS_OF_INTEREST_CONCAT.out.tsv.ifEmpty([])
            }
        }

        //
        // Analyses on top of the demix (all optional; an error in one never fails the run)
        //
        def ch_dropout_impact = channel.empty()
        def ch_dropout_summary = channel.empty()
        def ch_unexplained = channel.empty()
        def ch_timeseries = channel.empty()
        def ch_timeseries_lineages = channel.empty()
        def ch_growth_rates = channel.empty()
        def ch_growth_failed = channel.empty()
        def ch_growth_inputs = channel.empty()
        def ch_dashboard = channel.empty()

        if ((mode == 'full' || mode == 'demix') && (params.dropout_impact || params.unexplained_mutations)) {
            // both need the lineage -> mutation lists of the dataset: warn once and skip when it has none
            def mutation_reason = mutationTablesUnsupportedReason()
            if (mutation_reason) {
                log.warn "Skipping the dropout-impact and unexplained-mutation tables: ${mutation_reason}. Use a Freyja dataset directory (with lineage_mutations.json) or a barcodes .csv to get them."
            }
            else {
                // which dropped amplicons would have mattered for the lineages detected?
                // (amplicon_coverage.tsv exists only when some sample has a primer scheme; capture samples are not in it)
                if (params.dropout_impact) {
                    DROPOUT_IMPACT(
                        ch_amplicon_coverage.filter { coverage -> coverage },
                        FREYJA_AGGREGATE_PLOT.out.aggregated,
                        FREYJA_DATA.out.dataset,
                    )
                    ch_dropout_impact = DROPOUT_IMPACT.out.impact
                    ch_dropout_summary = DROPOUT_IMPACT.out.summary
                }

                // variants that no detected lineage explains (non-negative samples)
                if (params.unexplained_mutations) {
                    def ch_unexplained_in = ch_variants
                        .filter { meta, _variants -> meta.sample_type != 'negative' }
                        .join(ch_demixed, remainder: true)
                        .filter { meta, variants, _demix -> meta && variants }
                        .map { meta, variants, demix -> [ meta, variants, demix ?: [] ] }
                    UNEXPLAINED_MUTATIONS(ch_unexplained_in, FREYJA_DATA.out.dataset, ch_gff, ch_fasta)
                    UNEXPLAINED_MUTATIONS_CONCAT(UNEXPLAINED_MUTATIONS.out.tsv.map { _meta, tsv -> tsv }.collect())
                    ch_unexplained = UNEXPLAINED_MUTATIONS_CONCAT.out.tsv
                }
            }
        }

        //
        // Site aggregation (weighted average of the per-site RESULTS, never pooled reads) and the lineage-group summary.
        // demix_status.tsv always exists for them (header only when every sample is ok, e.g. aggregate mode).
        //
        def ch_site_freyja = channel.value([])
        def ch_site_meta = channel.value([])
        def ch_lineage_summary = channel.empty()
        def ch_lineage_resolved = channel.empty()
        def ch_site_aggregates = channel.empty()
        def ch_site_aggregates_meta = channel.empty()
        if (run_site_aggregation || (run_lineage_groups && !lineage_reason)) {
            // aggregate mode: the given .demix files are the usable samples (status ok, except the low_confidence ones
            // carried over from the earlier run, which are already in ch_status_lines), but the published
            // demix_status.tsv lists only the unusable ones, so the summaries get their own, complete copy
            def ch_given_ok = mode == 'aggregate'
                ? ch_demixed.filter { meta, _demix -> !carried_low.containsKey(meta.id) }.map { meta, _demix -> [ meta.id, 'ok', '', '', '', '' ].join('\t') }
                : channel.empty()
            // the samplesheet in effect plus the metadata of the samples discovery had to skip (a skipped negative control is
            // only known from there), as one CSV
            def ch_summary_sheet = ch_samplesheet
                .splitCsv(header: true)
                .map { row -> [ row.sample, row.sample_type ?: 'sample', row.collection_date, row.site, row.viral_load ].collect { v -> csvCell(v) }.join(',') }
                .mix(ch_skipped.map { sample, _reason, sample_type, date, site -> [ sample, sample_type, date, site, '' ].collect { v -> csvCell(v) }.join(',') })
                .collectFile(name: 'summary_samplesheet.csv', seed: 'sample,sample_type,collection_date,site,viral_load', sort: true, newLine: true)
                .first()
            def ch_status_for_summaries = ch_status_lines
                .mix(ch_given_ok)
                .ifEmpty('')
                .collectFile(name: 'demix_status.tsv', seed: status_header, sort: true, newLine: true)
                .first()
            if (run_site_aggregation) {
                SITE_AGGREGATE(
                    FREYJA_AGGREGATE_PLOT.out.aggregated,
                    ch_status_for_summaries,
                    ch_summary_sheet,
                    channel.value(file(params.site_aggregation, checkIfExists: true)),
                )
                ch_site_freyja = SITE_AGGREGATE.out.aggregated
                ch_site_meta = SITE_AGGREGATE.out.meta
                ch_site_aggregates = SITE_AGGREGATE.out.aggregated
                ch_site_aggregates_meta = SITE_AGGREGATE.out.meta
            }
            if (run_lineage_groups && !lineage_reason) {
                LINEAGE_SUMMARY(
                    FREYJA_AGGREGATE_PLOT.out.aggregated,
                    ch_status_for_summaries,
                    ch_summary_sheet,
                    channel.value(lineageGroupsFile()),
                    ch_hierarchy_dataset,
                    ch_site_freyja,
                    ch_site_meta,
                )
                ch_lineage_summary = LINEAGE_SUMMARY.out.summary
                ch_lineage_resolved = LINEAGE_SUMMARY.out.resolved
                // --freyja_boot: 95 % intervals of the groups from the raw bootstrap draws (summary/lineage_group_intervals.tsv)
                if (params.freyja_boot) {
                    BOOT_GROUP_INTERVALS(
                        FREYJA_AGGREGATE_PLOT.out.aggregated,
                        ch_boot_raw.toSortedList { a, b -> a.name <=> b.name }.filter { files -> files }, // nothing when no sample was bootstrapped
                        channel.value(lineageGroupsFile()),
                        ch_hierarchy_dataset,
                        ch_summary_sheet,
                    )
                }
            }
        }

        // Time-based analyses (time series, growth rates, dashboard) leave out low_confidence samples (PLAN D18)
        // unless --include_low_confidence: they get a copy of the aggregate without those rows, and the
        // two-dates test below does not count them. freyja_aggregated.tsv itself stays complete. In aggregate mode
        // only the samples carried over as low_confidence from --previous_results are affected (given .demix files
        // have no coverage information of their own).
        def ch_ts_samples = ch_aggregate_in
        def ch_analysis_agg = FREYJA_AGGREGATE_PLOT.out.aggregated
        if ((mode == 'full' || mode == 'demix') && !params.include_low_confidence && lowConfidenceBand()) {
            ch_ts_samples = ch_aggregate_in
                .join(FREYJA_DEMIX.out.status)
                .filter { _meta, _demix, status_json -> new groovy.json.JsonSlurper().parse(status_json).status != 'low_confidence' }
                .map { meta, demix, _status_json -> [ meta, demix ] }
            FILTER_AGGREGATED(FREYJA_AGGREGATE_PLOT.out.aggregated, ch_demix_status)
            ch_analysis_agg = FILTER_AGGREGATED.out.aggregated
        }
        else if (mode == 'aggregate' && !params.include_low_confidence && carried_low) {
            ch_ts_samples = ch_aggregate_in.filter { meta, _demix -> !carried_low.containsKey(meta.id) }
            FILTER_AGGREGATED(FREYJA_AGGREGATE_PLOT.out.aggregated, ch_demix_status)
            ch_analysis_agg = FILTER_AGGREGATED.out.aggregated
        }

        // Time series: needs at least two distinct collection dates among the non-negative aggregated samples
        // (growth rates and the dashboard check their own, stricter minimum and write a failure note instead)
        def ch_ts_go = ch_ts_samples
            .filter { meta, _demix -> meta.collection_date && meta.sample_type != 'negative' }
            .map { meta, _demix -> meta.collection_date }
            .unique()
            .count()
            .filter { n -> n >= 2 }
        def ch_ts_agg = ch_analysis_agg.combine(ch_ts_go).map { agg, _n -> agg }
        def ch_metadata = ch_samplesheet.first()

        if (params.timeseries) {
            TIMESERIES_TABLE(ch_ts_agg, ch_metadata)
            ch_timeseries = TIMESERIES_TABLE.out.table
            // second, lineage-level table (every Freyja lineage instead of the summarized groups)
            TIMESERIES_TABLE_LINEAGES(ch_ts_agg, ch_metadata)
            ch_timeseries_lineages = TIMESERIES_TABLE_LINEAGES.out.table
        }
        if (params.freyja_growth) {
            FREYJA_RELGROWTHRATE(ch_analysis_agg, ch_metadata, FREYJA_DATA.out.dataset)
            ch_growth_rates = FREYJA_RELGROWTHRATE.out.growth_rates
            ch_growth_failed = FREYJA_RELGROWTHRATE.out.failed
            ch_growth_inputs = FREYJA_RELGROWTHRATE.out.inputs
        }
        if (params.freyja_dash) {
            FREYJA_DASH(ch_analysis_agg, ch_metadata, FREYJA_DATA.out.dataset, params.freyja_dash_title, params.freyja_dash_intro)
            ch_dashboard = FREYJA_DASH.out.dashboard
        }

        // The samplesheet in effect (also for results-directory inputs) travels with the results, so a
        // later demix/aggregate run on this directory finds sample_type, collection_date and site
        ch_samplesheet.collectFile(name: 'samplesheet.csv', storeDir: "${outdir}/summary")

        //
        // Provenance: run_info.json
        //
        // (the per-sample library settings come from the metas of the rows, so the base info is made after they are known)
        def ch_base_info = ch_rows
            .map { row -> row[0] }
            .toList()
            .map { metas -> runInfoBase(mode, metas) }
            .collectFile(name: 'run_info_base.json')
        // every distinct primer BED the samples use
        def ch_beds = ch_rows
            .map { row -> row[0] }
            .filter { meta -> meta.primer_scheme }
            .map { meta -> primerBedFor(meta) }
            .unique()
            .toSortedList { a, b -> a.name <=> b.name }
        PROVENANCE(
            ch_base_info,
            ch_samplesheet,
            ch_beds,
            mode == 'aggregate' ? channel.value([]) : ch_fasta,
            FREYJA_DATA.out.freyja_data,
            ch_collated_versions,
            params.compare_barcodes && (mode == 'full' || mode == 'demix') ? FREYJA_DATA_COMPARE.out.freyja_data : channel.value([]),
        )
        ch_run_info = PROVENANCE.out.run_info

        //
        // Report: report_data.json (everything the run produced) -> report.html
        // Every run-level input is exactly one item (a file, or [] when absent).
        //
        def ch_report_in = PROVENANCE.out.run_info
            .combine(wrapItem(ch_qc_summary))
            .combine(wrapItem(ch_amplicon_coverage))
            .combine(wrapItem(FREYJA_AGGREGATE_PLOT.out.aggregated.ifEmpty([])))
            .combine(wrapItem(ch_demix_status.ifEmpty([])))
            .combine(wrapItem(ch_consensus_stats))
            .combine(wrapItem(ch_moi))
            .combine(wrapItem(ch_compare_tsv.ifEmpty([])))
            .combine(wrapItem(ch_depth_files))
        // optional extras: negative-control aggregate, the samplesheet in effect (metadata for modes
        // without a QC table), bootstrap intervals, status of the second-barcode demix
        // The sheet of the usable samples, plus the metadata of samples that discovery had to skip (a skipped
        // negative control is only known from there)
        def ch_skipped_sheet = ch_skipped
            .map { sample, _reason, sample_type, date, site -> "${sample},${sample_type},${date},${site}" }
            .collectFile(name: 'skipped_samplesheet.csv', seed: 'sample,sample_type,collection_date,site', sort: true, newLine: true)
        def ch_report_sheets = ch_samplesheet.mix(ch_skipped_sheet).collect()
        def ch_report_extra = wrapItem(ch_controls_aggregated.ifEmpty([]))
            .combine(wrapItem(ch_report_sheets))
            .combine(wrapItem(ch_boot_summary.ifEmpty([])))
            .combine(wrapItem(ch_compare_status.ifEmpty([])))
        def ch_report_analyses = wrapItem(ch_dropout_impact.ifEmpty([]))
            .combine(wrapItem(ch_dropout_summary.ifEmpty([])))
            .combine(wrapItem(ch_unexplained.ifEmpty([])))
            .combine(wrapItem(ch_timeseries.ifEmpty([])))
            .combine(wrapItem(ch_timeseries_lineages.ifEmpty([])))
            .combine(wrapItem(ch_growth_rates.ifEmpty([])))
            .combine(wrapItem(ch_growth_failed.ifEmpty([])))
            .combine(wrapItem(ch_growth_inputs.ifEmpty([])))
            .combine(ch_dashboard.map { dashboard -> dashboard.name }.ifEmpty(''))
        def ch_report_groups = wrapItem(ch_lineage_summary.ifEmpty([]))
            .combine(wrapItem(ch_lineage_resolved.ifEmpty([])))
            .combine(wrapItem(ch_site_aggregates.ifEmpty([])))
            .combine(wrapItem(ch_site_aggregates_meta.ifEmpty([])))
        BUILD_REPORT_DATA(ch_report_in, ch_report_extra, ch_report_analyses, ch_report_groups)
        RENDER_REPORT(BUILD_REPORT_DATA.out.json, file("${projectDir}/assets/templates/report.html.j2"))
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
// A file of <results>/summary/, or [] when the run did not write it
//
def optionalFile(results, name) {
    def candidate = file("${results}/summary/${name}")
    return candidate.exists() ? candidate : []
}

//
// Is the low-confidence band (PLAN D18) switched on? It is off when the floor is 0 or not below the
// high-confidence threshold: samples are then either demixed normally or skipped, as before D18.
//
def lowConfidenceBand() {
    return params.demix_floor_coverage > 0 && params.demix_floor_coverage < params.min_genome_coverage_for_demix
}

//
// Aggregate mode, --previous_results: [ sample: [ reason, coverage ] ] of the samples that summary/demix_status.tsv of
// an earlier run lists as low_confidence (PLAN D18). Other statuses, other inputs (--input, --demix_dir) and
// directories without that file give nothing: their .demix files count as ok.
//
def carriedLowConfidence() {
    def carried = [:]
    splitPaths(params.previous_results).each { dir ->
        def status_file = file("${dir}/summary/demix_status.tsv")
        if (status_file.exists()) {
            def lines = status_file.readLines().findAll { line -> line.trim() }
            def header = lines ? lines[0].split('\t', -1).collect { column -> column.trim() } : []
            lines.drop(1).each { line ->
                def cells = line.split('\t', -1)
                def row = [:]
                header.eachWithIndex { column, i -> row[column] = i < cells.size() ? cells[i].trim() : '' }
                if (row.sample && row.status == 'low_confidence') {
                    carried[row.sample] = [ reason: row.reason ?: '', coverage: row.genome_coverage_pct ?: '' ]
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
def runInfoBase(mode, metas) {
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
    def plain_params = params.findAll { key, _value -> !(key.toString() in ['pipelines_testdata_base_path', 'lineage_groups_bundled', 'freyja_container_pathogen']) }.collectEntries { key, value ->
        [ (key.toString()): (value == null || value instanceof Boolean || value instanceof Number) ? value : value.toString() ]
    }
    // the effective primer scheme is recorded even when the default (null = the default scheme, D14) applies
    if (mode == 'full' && protocolSettings(params.protocol).primer_trimming && !params.primer_bed && !params.primer_scheme) {
        plain_params.primer_scheme = defaultPrimerScheme()
    }
    // lineage groups: which file made lineage_summary.csv (only when the summary is made)
    def groups_info = [ lineage_groups_file: null, lineage_groups_source: null, lineage_groups_sha256: null ]
    if (params.lineage_summary && lineageGroupsFile() && mode in ['full', 'demix', 'aggregate'] && !lineageHierarchyUnsupportedReason()) {
        groups_info = lineageGroupsInfo()
    }
    def info = [
        pipeline_name: workflow.manifest.name,
        pipeline_version: workflow.manifest.version,
        git_commit: git_commit,
        nextflow_version: workflow.nextflow.version.toString(),
        mode: mode,
        profile: workflow.profile,
        command_line: workflow.commandLine,
        run_name: workflow.runName,
        start: workflow.start.toString(),
        params: plain_params,
        // earlier runs whose results this run builds on (demix/aggregate); barcode_version 'unknown'
        // for legacy Snakemake results, which recorded none
        previous_runs: previousRuns(mode),
    ] + libraryInfo(mode, metas) + groups_info
    if (mode == 'aggregate') {
        info.freyja_data_note = 'aggregate mode demixes nothing: freyja_data describes the dataset whose lineage hierarchy groups the plot, not the barcodes the .demix files were made with (see previous_runs)'
    }
    return groovy.json.JsonOutput.prettyPrint(groovy.json.JsonOutput.toJson(info))
}
