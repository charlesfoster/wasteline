#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    charlesfoster/wasteline
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/charlesfoster/wasteline
----------------------------------------------------------------------------------------
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { WASTELINE  } from './workflows/wasteline'
include { PIPELINE_INITIALISATION } from './subworkflows/local/utils_nfcore_wasteline_pipeline'
include { PIPELINE_COMPLETION     } from './subworkflows/local/utils_nfcore_wasteline_pipeline'
include { RESULTS_INPUT           } from './subworkflows/local/results_input'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    TYPED BOOLEAN PARAMETERS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Nextflow 26.04 (strict parser, the default) hands every command-line parameter to the pipeline
// as a String, so nf-schema rejects `--skip_consensus`, `--save_bam false`, `--save_bam=false` ...
// ("expected type Boolean, found String"). Declaring the type here makes Nextflow convert the
// value, so all of these work: `--x`, `--x true`, `--x false`, `--x=true`, `--x=false`, and
// -params-file. This block REPEATS the boolean defaults of nextflow.config (which stays the place
// every other default lives); tests/python/test_boolean_params.py fails when the two disagree.
// This is why the pipeline requires Nextflow >= 26.04: the block only parses with the strict (v2)
// parser, which is the default from 26.04. On 25.10 the old parser cannot read it, and with
// NXF_SYNTAX_PARSER=v2 the converted values do not reach included scripts and config files.
params {
    help                           : Boolean = false
    help_full                      : Boolean = false
    show_hidden                    : Boolean = false
    version                        : Boolean = false
    validate_params                : Boolean = true
    plaintext_email                : Boolean = false
    monochrome_logs                : Boolean = false
    fastp_cut_front                : Boolean = true
    fastp_cut_tail                 : Boolean = true
    fastp_trim_poly_x              : Boolean = true
    fastp_correction               : Boolean = true
    save_trimmed_reads             : Boolean = false
    save_bam                       : Boolean = true
    skip_primer_trimming           : Boolean = false
    freyja_autoadapt               : Boolean = false
    freyja_relaxedmrca             : Boolean = false
    freyja_confirmedonly           : Boolean = false
    include_negatives_in_aggregate : Boolean = false
    skip_consensus                 : Boolean = false
    consensus_iupac                : Boolean = true
    skip_report                    : Boolean = false
    freyja_boot                    : Boolean = false
    annotate_variants              : Boolean = false

    // numbers arrive as Strings too (--consensus_min_depth 20): type them
    ivar_offset                     : Integer? = null
    fastp_cut_mean_quality          : Integer = 20
    fastp_qualified_quality_phred   : Integer = 20
    fastp_unqualified_percent_limit : Integer = 10
    fastp_length_required           : Integer = 50
    freyja_minq                     : Integer = 20
    freyja_varthresh                : Float = 0.0
    freyja_eps                      : Float = 0.001
    freyja_covcut                   : Integer = 10
    freyja_depthcutoff              : Integer = 0
    freyja_adapt                    : Float = 0.0
    freyja_boot_n                   : Integer = 100
    freyja_boot_seed                : Integer = 0
    freyja_plot_mincov              : Float = 60
    freyja_plot_windowsize          : Integer = 14
    min_genome_coverage_for_demix   : Float = 60
    amplicon_min_depth              : Integer = 10
    qc_min_depth                    : Integer = 10
    qc_pass_min_pct                 : Float = 80
    qc_fail_below_pct               : Float = 50
    qc_max_dropped_amplicons        : Integer = 1
    ntc_max_mapped_reads            : Integer = 1000
    ntc_max_genome_coverage_pct     : Float = 10
    consensus_min_depth             : Integer = 10
    consensus_freq                  : Float = 0.75
    indel_freq                      : Float = 0.75
    consensus_min_alt_qual          : Integer = 20
    consensus_min_depth_variants    : Integer = 10
    consensus_iupac_min_af          : Float = 0.25
    consensus_iupac_min_alt_reads   : Integer = 5
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOWS FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// WORKFLOW: Run main analysis pipeline depending on the mode
//
workflow CHARLESFOSTER_WASTELINE {

    take:
    mode        // string: full | demix | aggregate | report | update_db
    samplesheet // channel: rows of the samplesheet read from --input (empty when it is not used)

    main:

    def ch_rows = samplesheet
    def ch_sheet = channel.empty()
    def ch_skipped = channel.empty()

    if (params.input) {
        ch_sheet = channel.of(file(params.input))
    }
    else if (mode == 'demix' || mode == 'aggregate') {
        // demix/aggregate from results directories (--previous_results, --demix_dir): the same rows
        // as a samplesheet, discovered by bin/make_samplesheet.py --from-results
        RESULTS_INPUT(mode, params.previous_results, params.demix_dir)
        ch_rows = RESULTS_INPUT.out.rows
        ch_sheet = RESULTS_INPUT.out.samplesheet
        ch_skipped = RESULTS_INPUT.out.skipped
    }

    //
    // WORKFLOW: Run pipeline
    //
    WASTELINE (
        mode,
        ch_rows,
        ch_sheet,
        ch_skipped,
        params.outdir,
    )
    emit:
    WASTELINE.out.run_info // channel: /path/to/run_info.json
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// The only entry point. The mode is --mode full|demix|aggregate|report|update_db: Nextflow's strict
// parser (default from 26.04) does not support `-entry`, so there are no per-mode entry workflows.
//
workflow {

    main:
    //
    // SUBWORKFLOW: Run initialisation tasks
    //
    PIPELINE_INITIALISATION (
        params.mode,
        params.version,
        params.validate_params,
        params.monochrome_logs,
        args,
        params.outdir,
        params.input,
        params.help,
        params.help_full,
        params.show_hidden
    )

    //
    // WORKFLOW: Run main workflow
    //
    CHARLESFOSTER_WASTELINE (
        params.mode,
        PIPELINE_INITIALISATION.out.samplesheet
    )
    //
    // SUBWORKFLOW: Run completion tasks
    //
    PIPELINE_COMPLETION (
        params.email,
        params.email_on_fail,
        params.plaintext_email,
        params.outdir,
        params.monochrome_logs
    )
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
