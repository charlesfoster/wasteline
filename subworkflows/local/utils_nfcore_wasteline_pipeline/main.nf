//
// Subworkflow with functionality specific to the charlesfoster/wasteline pipeline
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { UTILS_NFSCHEMA_PLUGIN     } from '../../nf-core/utils_nfschema_plugin'
include { paramsSummaryMap          } from 'plugin/nf-schema'
include { samplesheetToList         } from 'plugin/nf-schema'
include { paramsHelp                } from 'plugin/nf-schema'
include { completionEmail           } from '../../nf-core/utils_nfcore_pipeline'
include { completionSummary         } from '../../nf-core/utils_nfcore_pipeline'
include { UTILS_NFCORE_PIPELINE     } from '../../nf-core/utils_nfcore_pipeline'
include { UTILS_NEXTFLOW_PIPELINE   } from '../../nf-core/utils_nextflow_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW TO INITIALISE PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PIPELINE_INITIALISATION {

    take:
    mode              //  string: full | demix | aggregate | report | update_db
    version           // boolean: Display version and exit
    validate_params   // boolean: Boolean whether to validate parameters against the schema at runtime
    monochrome_logs   // boolean: Do not use coloured log outputs
    nextflow_cli_args //   array: List of positional nextflow CLI args
    outdir            //  string: The output directory where the results will be saved
    input             //  string: Path to input samplesheet (not needed by every mode)
    help              // boolean: Display help message and exit
    help_full         // boolean: Show the full help message
    show_hidden       // boolean: Show hidden parameters in the help message

    main:

    ch_versions = channel.empty()

    //
    // Print version and exit if required and dump pipeline parameters to JSON file
    //
    UTILS_NEXTFLOW_PIPELINE (
        version,
        true,
        outdir,
        workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1
    )

    //
    // Validate parameters and generate parameter summary to stdout
    //

    def before_text = wastelineLogo(monochrome_logs)
    def after_text = ""

    command = "nextflow run ${workflow.manifest.name} -profile <docker/singularity/.../institute> --input samplesheet.csv --outdir <OUTDIR>"

    UTILS_NFSCHEMA_PLUGIN (
        workflow,
        validate_params,
        null,
        help,
        help_full,
        show_hidden,
        before_text,
        after_text,
        command,
        false
    )

    //
    // Check config provided to the pipeline
    //
    UTILS_NFCORE_PIPELINE (
        nextflow_cli_args
    )

    //
    // Custom validation for pipeline parameters (nf-schema cannot express "required depending on mode")
    //
    validateInputParameters(mode)

    //
    // Rows of the samplesheet given with --input, validated against the schema of the mode
    // (assets/schema_input*.json). Sample IDs are checked (pattern + uniqueness) and never rewritten.
    //   full:      [ meta, [ fastq_1, fastq_2 ] ]
    //   demix:     [ meta, variants, depths, fastp_json ]
    //   aggregate: [ meta, demix ]
    // Empty when no samplesheet is used (demix/aggregate from results directories, report, update_db):
    // the RESULTS_INPUT subworkflow produces the same rows from --previous_results / --demix_dir.
    //
    // The rows are built eagerly, so a bad `protocol` / `primer_scheme` cell stops the run here, naming the sample.
    def ch_samplesheet = channel.empty()
    def rows = []
    if (input && mode == 'full') {
        rows = samplesheetToList(input, "${projectDir}/assets/schema_input.json")
            .collect { meta, fastq_1, fastq_2 -> [ sampleMeta(meta, mode), [ fastq_1, fastq_2 ] ] }
    }
    else if (input && mode == 'demix') {
        rows = samplesheetToList(input, "${projectDir}/assets/schema_input_demix.json")
            .collect { meta, variants, depths, fastp_json -> [ sampleMeta(meta, mode), variants, depths, fastp_json ?: [] ] }
    }
    else if (input && mode == 'aggregate') {
        rows = samplesheetToList(input, "${projectDir}/assets/schema_input_aggregate.json")
            .collect { meta, demix -> [ sampleMeta(meta, mode), demix ] }
    }
    if (rows) {
        checkSampleBeds(rows.collect { row -> row[0] })
        ch_samplesheet = channel.fromList(rows)
    }

    emit:
    samplesheet = ch_samplesheet
    versions    = ch_versions
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW FOR PIPELINE COMPLETION
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PIPELINE_COMPLETION {

    take:
    email           //  string: email address
    email_on_fail   //  string: email address sent on pipeline failure
    plaintext_email // boolean: Send plain-text email instead of HTML
    outdir          //    path: Path to output directory where results will be published
    monochrome_logs // boolean: Disable ANSI colour codes in log output

    main:
    summary_params = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")

    //
    // Completion email and summary
    //
    workflow.onComplete {
        if (email || email_on_fail) {
            completionEmail(
                summary_params,
                email,
                email_on_fail,
                plaintext_email,
                outdir,
                monochrome_logs,
                [],
            )
        }

        finaliseRunInfo(outdir)
        completionSummary(monochrome_logs)

    }

    workflow.onError {
        log.error "Pipeline failed. Please refer to troubleshooting docs for common issues: https://nf-co.re/docs/running/troubleshooting"
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
//
// ASCII logo (assets/logo.txt) with version line; ANSI colour unless monochrome_logs
//
def wastelineLogo(monochrome_logs) {
    def logo = file("${projectDir}/assets/logo.txt").text
    def version = "wasteline v${workflow.manifest.version}"
    if (monochrome_logs || System.getenv('NO_COLOR')) {
        return "\n${logo}\n  ${version}\n"
    }
    return "\n\033[0;33m${logo}\033[0m\n  \033[1m${version}\033[0m\n"
}

//
// Container image per process. workflow.container is empty for images declared inside modules,
// so fall back to the `container` column of the execution trace.
//
def runContainers(outdir) {
    def containers = workflow.container ? workflow.container.collectEntries { key, value -> [ (key.toString()): value.toString() ] } : [:]
    def trace = file("${outdir}/pipeline_info/execution_trace_${params.trace_report_suffix}.txt")
    if (!containers && trace.exists()) {
        def lines = trace.readLines()
        def header = lines[0].split('\t') as List
        def process_col = header.indexOf('process')
        def container_col = header.indexOf('container')
        if (process_col >= 0 && container_col >= 0) {
            lines.drop(1).each { line ->
                def cols = line.split('\t', -1)
                if (cols.size() > container_col && cols[container_col] && cols[container_col] != '-') {
                    containers[cols[process_col].tokenize(':').last()] = cols[container_col]
                }
            }
        }
    }
    return containers.sort()
}

//
// Add what is only known when the run ends (container images per process, completion time and
// outcome) to pipeline_info/run_info.json, which PROVENANCE wrote during the run
//
def finaliseRunInfo(outdir) {
    def run_info = file("${outdir}/pipeline_info/run_info.json")
    if (run_info.exists()) {
        def info = new groovy.json.JsonSlurper().parseText(run_info.text)
        info.containers = runContainers(outdir)
        info.run_completed = workflow.complete?.toString()
        info.duration = workflow.duration?.toString()
        info.success = workflow.success
        run_info.text = groovy.json.JsonOutput.prettyPrint(groovy.json.JsonOutput.toJson(info))
    }
}

//
// Sample meta map shared by every mode. In modes that map or demix (full, demix) it also carries the library
// settings of the sample (PLAN D20), see libraryMeta().
//
def sampleMeta(meta, mode) {
    def sample = [
        id: meta.id,
        sample_type: meta.sample_type ?: 'sample',
        collection_date: meta.collection_date ?: null,
        site: meta.site ?: null,
        viral_load: meta.viral_load ?: null,
    ]
    if (mode in ['full', 'demix']) {
        return sample + libraryMeta(meta.id, meta.protocol, meta.primer_scheme, mode)
    }
    return sample
}

//
// THE effective library settings of one sample, computed here and nowhere else (conf/modules.config spells
// out the naming conditions it needs, because a config file cannot call a function). `protocol` and
// `primer_scheme` are the optional samplesheet cells; empty = the run defaults (--protocol, --primer_bed /
// --primer_scheme, --skip_primer_trimming, --duplicates).
//   protocol         amplicon | capture
//   primer_trimming  true only for amplicon reads that are not --skip_primer_trimming
//   duplicates       mark | none (--duplicates auto = mark for capture, none for amplicon)
//   primer_scheme    bundled scheme name, or 'custom' (then primer_bed is the BED path); null when the sample has no
//                    amplicon-level analysis: not primer-trimmed, or demix mode without a named scheme (the scheme
//                    of old results is unknown). Full mode with nothing named uses the bundled default (D14).
//
def libraryMeta(sample, protocol_cell, scheme_cell, mode) {
    def protocol = protocol_cell ?: params.protocol
    if (!(protocol in ['amplicon', 'capture'])) {
        error("Sample '${sample}': protocol '${protocol}' is not one of amplicon, capture")
    }
    def settings = protocolSettings(protocol)
    def scheme = null
    def bed = null
    if (settings.primer_trimming) {
        if (scheme_cell) {
            def spec = scheme_cell.toString().trim()
            if (spec.contains('/') || spec.endsWith('.bed')) {
                bed = file(spec).toAbsolutePath().toString()
                if (!file(bed).exists()) {
                    error("Sample '${sample}': primer_scheme BED file '${spec}' does not exist (a relative path is relative to the directory you launch Nextflow from)")
                }
                scheme = 'custom'
            }
            else if (file("${projectDir}/assets/primer_schemes/${spec}.bed").exists()) {
                scheme = spec
            }
            else {
                error("Sample '${sample}': primer_scheme '${spec}' is not a bundled scheme (${bundledPrimerSchemes().join(', ')}) and not a path to a BED file")
            }
        }
        else if (params.primer_bed) {
            scheme = 'custom'
            bed = params.primer_bed.toString()
        }
        else if (params.primer_scheme) {
            scheme = params.primer_scheme.toString()
        }
        else if (mode == 'full') {
            scheme = defaultPrimerScheme()
        }
    }
    return settings + [ primer_scheme: scheme, primer_bed: bed ]
}

//
// protocolSettings: the part of the library settings that depends only on the protocol
//
def protocolSettings(protocol) {
    def capture = protocol == 'capture'
    return [
        protocol: protocol,
        primer_trimming: !capture && !params.skip_primer_trimming,
        duplicates: params.duplicates == 'auto' ? (capture ? 'mark' : 'none') : params.duplicates,
    ]
}

def bundledPrimerSchemes() {
    return file("${projectDir}/assets/primer_schemes").listFiles().findAll { f -> f.name.endsWith('.bed') }.collect { f -> f.baseName }.sort()
}

//
// The primer BED of a sample with a primer_scheme, as given (the user's file or the bundled scheme)
//
def sourceBedFor(meta) {
    if (meta.primer_bed) {
        return file(meta.primer_bed, checkIfExists: true)
    }
    return file("${projectDir}/assets/primer_schemes/${meta.primer_scheme}.bed", checkIfExists: true)
}

//
// ... as it is used by the processes: restricted to the reference's contig(s) (reconcileBed, PLAN D24)
//
def primerBedFor(meta) {
    return reconcileBed(sourceBedFor(meta)).bed
}

//
// A primer BED against the reference (PLAN D24). Whole-genome panels (influenza) list the primers of every segment, while
// the reference may be one segment (HA): the lines of contigs that are not in the reference FASTA are dropped, and the
// filtered copy (the same file name, in <workDir>/wasteline_primer_beds/<hash of its content>/) is what trimming, amplicon
// coverage and provenance use. A BED whose contigs are all in the reference is returned as it is (the same file, byte for
// byte). Returns [bed: Path, no_match: boolean, dropped_lines: n, dropped_contigs: [names]]. Comment, track and browser
// lines are kept. A BED none of whose lines is on the reference is returned unchanged with no_match = true (the callers
// refuse it with their own message).
//
def reconcileBed(bed) {
    if (!params.fasta) {
        return [ bed: bed, no_match: false, dropped_lines: 0, dropped_contigs: [] ]
    }
    def contigs = fastaContigs(file(params.fasta, checkIfExists: true))
    def keep = []
    def kept = 0
    def dropped_lines = 0
    def dropped_contigs = []
    bed.readLines().each { line ->
        def record = line.trim() && !line.startsWith('#') && !line.startsWith('track') && !line.startsWith('browser')
        def contig = record ? line.split('\t')[0] : null
        if (!record) {
            keep << line
        }
        else if (contig in contigs) {
            keep << line
            kept += 1
        }
        else {
            dropped_lines += 1
            if (!(contig in dropped_contigs)) {
                dropped_contigs << contig
            }
        }
    }
    if (dropped_lines == 0 || kept == 0) {
        return [ bed: bed, no_match: kept == 0 && dropped_lines > 0, dropped_lines: dropped_lines, dropped_contigs: dropped_contigs ]
    }
    def text = keep.join('\n') + '\n'
    def digest = java.security.MessageDigest.getInstance('SHA-256').digest(text.getBytes('UTF-8')).collect { b -> String.format('%02x', b) }.join('').substring(0, 12)
    def filtered = file("${workflow.workDir}/wasteline_primer_beds/${digest}/${bed.name}")
    if (!filtered.exists()) {
        filtered.parent.mkdirs()
        filtered.text = text
    }
    return [ bed: filtered, no_match: false, dropped_lines: dropped_lines, dropped_contigs: dropped_contigs ]
}

//
// The one warning about primer BED lines dropped because their contig is not in the reference (PLAN D24)
//
def warnDroppedBedLines(reports) {
    def lines = reports.sum { report -> report.dropped_lines } ?: 0
    if (lines) {
        def contigs = reports.collectMany { report -> report.dropped_contigs }.unique()
        log.warn("Primer BED: ${lines} line(s) on ${contigs.size()} contig(s) that are not in the reference FASTA '${file(params.fasta).name}' (${contigs.take(10).join(', ')}${contigs.size() > 10 ? ', ...' : ''}) were dropped; the reference's own lines are used for trimming and amplicon coverage. Expected for a whole-genome panel run against one segment (for example influenza HA).")
    }
}

//
// What run_info.json records about the library settings of the samples (modes full and demix; [] elsewhere):
//   protocol / primer_trimming / duplicates   the common value, or 'mixed' when the samples differ
//                                             (primer_trimming and duplicates are only reported for mode full)
//   sample_libraries                          per sample: protocol, primer_scheme ('custom: <BED name>' for a BED
//                                             file, null when there is no amplicon-level analysis), primer_trimming, duplicates
//
def commonValue(metas, key) {
    def values = metas.collect { meta -> meta[key] }.unique()
    return values.size() == 1 ? values[0] : 'mixed'
}

def libraryInfo(mode, metas) {
    if (!(mode in ['full', 'demix']) || !metas) {
        return [ protocol: params.protocol, primer_trimming: null, duplicates: null ]
    }
    def libraries = metas.sort { meta -> meta.id }.collectEntries { meta ->
        [ (meta.id): [
            protocol: meta.protocol,
            primer_scheme: meta.primer_bed ? "custom: ${file(meta.primer_bed).name}" : meta.primer_scheme,
            primer_trimming: meta.primer_trimming,
            duplicates: meta.duplicates,
        ] ]
    }
    def trimming = commonValue(metas, 'primer_trimming')
    return [
        protocol: commonValue(metas, 'protocol'),
        primer_trimming: mode == 'full' && trimming != 'mixed' ? trimming : null,
        duplicates: mode == 'full' ? commonValue(metas, 'duplicates') : null,
        sample_libraries: libraries,
    ]
}

//
// Every distinct primer BED of the samples must be on a contig of the reference or ivar trim silently trims nothing
//
def checkSampleBeds(metas) {
    def beds = metas.findAll { meta -> meta.primer_scheme }.collect { meta -> sourceBedFor(meta) }.unique()
    if (beds) {
        def fasta_contigs = fastaContigs(file(params.fasta, checkIfExists: true))
        def reports = beds.collect { bed ->
            def report = reconcileBed(bed)
            if (report.no_match) {
                error("Primer BED '${bed}' refers to contig(s) ${bedContigs(bed)} that are not in the reference FASTA '${params.fasta}' (which has ${fasta_contigs}). Use a matching primer scheme or --fasta.")
            }
            report
        }
        warnDroppedBedLines(reports)
    }
}

//
// Freyja pathogens (`freyja update --pathogen` choices) that wasteline cannot run, with the reason and what would be
// needed (PLAN D22); null for every other name. The single-contig FASTA check in validateInputParameters catches
// a multi-segment reference of any other pathogen. The influenza pathogens are supported on the HA segment (PLAN D24).
//
def unsupportedPathogenReason(pathogen) {
    if (pathogen == 'MTB') {
        return "Mycobacterium tuberculosis is a 4.4 Mb bacterial genome with feather-only Freyja barcodes. wasteline is built and tested for viral genomes of 10-200 kb (per-position depth tables, QC, consensus and the dropout/mutation tables are not sized for 4.4 Mb), and nothing has been validated on bacterial data. Supporting it needs those steps scaled and checked on real MTB data."
    }
    return null
}

//
// Comma-separated path list (--previous_results) as a list of strings
//
def splitPaths(spec) {
    return spec ? spec.toString().split(',').collect { entry -> entry.trim() }.findAll { entry -> entry } : []
}

def requireDirectory(option, value) {
    def dir = file(value)
    if (!dir.exists()) {
        error("${option} '${value}' does not exist")
    }
    if (!dir.isDirectory()) {
        error("${option} '${value}' is not a directory")
    }
    return dir
}

//
// Check and validate pipeline parameters for the mode. Only `outdir` is required by the schema for
// every mode; what else is needed depends on the mode and is checked here with a helpful message.
//
def validateInputParameters(mode) {
    def modes = ['full', 'demix', 'aggregate', 'report', 'update_db']
    if (!(mode in modes)) {
        error("--mode '${mode}' is not one of ${modes.join(', ')}")
    }
    def sources = [
        '--input': params.input,
        '--previous_results': params.previous_results,
        '--demix_dir': params.demix_dir,
    ].findAll { _name, value -> value }

    def accepted = [
        full: ['--input'],
        demix: ['--input', '--previous_results'],
        aggregate: ['--input', '--demix_dir', '--previous_results'],
        report: [],
        update_db: [],
    ][mode]
    def hint = [
        full: "--mode full needs a FASTQ samplesheet: --input samplesheet.csv (columns sample,fastq_1,fastq_2[,sample_type][,collection_date][,site][,protocol][,primer_scheme]; make one with bin/make_samplesheet.py <fastq_dir>)",
        demix: "--mode demix needs existing Freyja variants and depths: --input samplesheet.csv (columns sample,variants,depths[,fastp_json][,sample_type][,collection_date][,site][,protocol][,primer_scheme]) or --previous_results <results dir> (new or legacy Snakemake layout; never modified)",
        aggregate: "--mode aggregate needs existing .demix files: --input samplesheet.csv (columns sample,demix[,sample_type][,collection_date][,site]), --demix_dir <dir of *.demix files>, or --previous_results <dir>[,<dir>...]",
        report: "--mode report needs --results_dir <results directory of a previous run> (it takes no samplesheet)",
        update_db: "--mode update_db takes no input; it downloads the latest Freyja data into --db_dir",
    ][mode]

    // report and update_db take no input at all: a stray --input (e.g. the default of a test profile) is ignored
    if (mode in ['report', 'update_db'] && params.input) {
        log.warn("--input is not used by --mode ${mode} and is ignored")
        sources.remove('--input')
    }
    def unexpected = sources.keySet().findAll { name -> !(name in accepted) }
    if (unexpected) {
        error("${unexpected.join(', ')} cannot be used with --mode ${mode}. ${hint}")
    }
    if (accepted && sources.size() == 0) {
        error("No input given for --mode ${mode}. ${hint}")
    }
    if (sources.size() > 1) {
        error("Give only one of ${sources.keySet().join(', ')} with --mode ${mode}. ${hint}")
    }

    // The pathogen (PLAN D22): the Freyja pathogens that need more than a reference are refused before any work is done
    if (mode in ['full', 'demix', 'aggregate', 'update_db']) {
        def unsupported = unsupportedPathogenReason(params.pathogen)
        if (unsupported) {
            error("--pathogen ${params.pathogen} is not supported: ${unsupported} See docs/other_pathogens.md.")
        }
        if (mode in ['full', 'demix']) {
            if (!params.fasta) {
                error("--fasta is required with --pathogen ${params.pathogen}: wasteline bundles a reference only for the pathogens that have a file in conf/pathogens/. Pass the reference.fasta that the Freyja barcodes of ${params.pathogen} were built on (the pathogen's folder in https://github.com/andersen-lab/Freyja-barcodes): barcodes are tied to their reference and a different one gives silently wrong results. --gff (annotation) is optional. See docs/other_pathogens.md.")
            }
            def contigs = fastaContigs(file(params.fasta, checkIfExists: true))
            if (contigs.size() != 1) {
                error("The reference FASTA '${params.fasta}' has ${contigs.size()} sequences (${contigs.take(10).join(', ')}${contigs.size() > 10 ? ', ...' : ''}): wasteline supports single-contig references only (depth QC, amplicon coverage, low-depth masking, consensus and the primer checks assume one contig). Multi-segment genomes need per-contig handling throughout the pipeline, which is not implemented (for influenza, wasteline bundles the single-sequence HA reference of each Freyja influenza pathogen: use --pathogen H1N1, H3N2, FLU-B-VIC, H5Nx or H5Nx-cattle without --fasta). See docs/other_pathogens.md.")
            }
        }
    }
    if (params.barcodes == 'container' || params.compare_barcodes == 'container') {
        if (params.pathogen != params.freyja_container_pathogen) {
            error("--barcodes container is only available for the default pathogen (the only one the Freyja container bundles barcodes for): the pinned Freyja container bundles no barcodes for ${params.pathogen}. Use --barcodes latest, or a dataset directory or barcode file made by `freyja update --pathogen ${params.pathogen}`.")
        }
    }

    // Run defaults: samples without a protocol / primer_scheme cell use them. A samplesheet can override them per sample (PLAN D20).
    def settings = protocolSettings(params.protocol)
    if ((params.protocol == 'capture' || params.skip_primer_trimming) && (params.primer_scheme || params.primer_bed) && mode in ['full', 'demix']) {
        log.warn("--primer_scheme/--primer_bed are ignored for samples that are not primer-trimmed (" + (params.protocol == 'capture' ? "--protocol capture" : "--skip_primer_trimming") + "): there is no amplicon-level analysis for them")
    }
    if (settings.primer_trimming && (mode == 'full' || (mode == 'demix' && (params.primer_scheme || params.primer_bed)))) {
        // The default BED must be on the same contig(s) as the reference or ivar trim silently trims nothing
        // (BEDs named in the samplesheet are checked when it is read)
        // (PLAN D24: lines of contigs that are not in the reference are dropped; only a BED with no line on it is refused)
        def bed_report = reconcileBed(primerBed())
        def fasta_contigs = fastaContigs(file(params.fasta, checkIfExists: true))
        def unknown = bed_report.no_match ? bedContigs(primerBed()) - fasta_contigs : []
        if (unknown && !params.primer_bed && !params.primer_scheme) {
            error("The bundled default primer scheme '${defaultPrimerScheme()}' is for a different reference (its contig(s) ${unknown} are not in the reference FASTA '${params.fasta}', which has ${fasta_contigs}): the bundled primer schemes belong to one pathogen only. Pass --primer_bed <BED whose first column is ${fasta_contigs.join(' / ')}> (a primer scheme for ${params.pathogen}), or --protocol capture / --skip_primer_trimming for reads that are not primer-trimmed. See docs/other_pathogens.md.")
        }
        if (unknown) {
            error("Primer BED '${primerBed()}' refers to contig(s) ${unknown} that are not in the reference FASTA '${params.fasta}' (which has ${fasta_contigs}). Use a matching --primer_bed/--primer_scheme or --fasta.")
        }
        if (!params.input) {
            warnDroppedBedLines([ bed_report ]) // with a samplesheet, checkSampleBeds warns once for all BEDs
        }
    }
    if (params.lineage_groups && !params.lineage_summary && mode in ['full', 'demix', 'aggregate']) {
        log.warn("--lineage_groups is ignored with --lineage_summary false")
    }
    if (mode in ['report', 'update_db'] && (params.lineage_groups || params.site_aggregation)) {
        log.warn("--lineage_groups/--site_aggregation are not used by --mode ${mode} and are ignored")
    }
    if (mode in ['report', 'update_db', 'aggregate'] && (params.primer_scheme || params.primer_bed)) {
        log.warn("--primer_scheme/--primer_bed are not used by --mode ${mode} and are ignored")
    }
    if (params.previous_results) {
        splitPaths(params.previous_results).each { dir -> requireDirectory('--previous_results', dir) }
    }
    if (params.demix_dir) {
        requireDirectory('--demix_dir', params.demix_dir)
    }
    if (mode == 'report') {
        if (!params.results_dir) {
            error("No --results_dir given. ${hint}")
        }
        def dir = requireDirectory('--results_dir', params.results_dir)
        if (params.outdir && file(params.outdir).toAbsolutePath().normalize() == dir.toAbsolutePath().normalize()) {
            error("--outdir must differ from --results_dir ('${params.results_dir}'): report mode never writes into the results it reads")
        }
        if (!file("${dir}/pipeline_info/run_info.json").exists()) {
            error("--results_dir '${params.results_dir}' has no pipeline_info/run_info.json: report mode needs the provenance record of the run (it never invents one)")
        }
    }
    else if (params.results_dir) {
        error("--results_dir is only used by --mode report (this is --mode ${mode}). ${hint}")
    }
    // Resolve the Freyja dataset specification early so a typo fails before any work is done
    if (mode in ['full', 'demix', 'aggregate']) {
        resolveBarcodes(params.barcodes, dbDir(), '--barcodes')
    }
    if (params.compare_barcodes) {
        if (mode in ['full', 'demix']) {
            resolveBarcodes(params.compare_barcodes, dbDir(), '--compare_barcodes')
        }
        else {
            log.warn("--compare_barcodes is only used by --mode full and demix and is ignored in --mode ${mode}")
        }
    }
    if (params.freyja_boot && !(mode in ['full', 'demix'])) {
        log.warn("--freyja_boot needs variants and depths: it is only used by --mode full and demix and is ignored in --mode ${mode}")
    }
    if (params.annotate_variants) {
        if (!(mode in ['full', 'demix'])) {
            log.warn("--annotate_variants is only used by --mode full and demix and is ignored in --mode ${mode}")
        }
        else if (params.skip_consensus) {
            error("--annotate_variants annotates the variants the consensus is built from: it cannot be combined with --skip_consensus")
        }
        else if (!params.gff) {
            error("--annotate_variants needs the reference annotation: --pathogen ${params.pathogen} has no bundled GFF3, pass --gff")
        }
    }
}

//
// Directory holding pipeline_info/run_info.json for a path inside a previous run (looks up at most
// 6 levels), or null
//
def findRunRoot(path) {
    return findRunRootFrom(new File(path.toString()).absoluteFile, 6)
}

def findRunRootFrom(dir, levels) {
    if (dir == null || levels == 0) {
        return null
    }
    if (new File(dir, 'pipeline_info/run_info.json').isFile()) {
        return dir
    }
    return findRunRootFrom(dir.parentFile, levels - 1)
}

//
// Provenance of the results a demix/aggregate run builds on: [ path, run_info_sha256, barcode_version ].
// run_info_sha256 is null and barcode_version 'unknown' when there is no run_info.json (legacy
// Snakemake results record neither).
//
def previousRunRecord(dir) {
    def run_info = new File(dir.toString(), 'pipeline_info/run_info.json')
    def record = [ path: dir.toString(), run_info_sha256: null, barcode_version: 'unknown' ]
    if (run_info.isFile()) {
        def digest = java.security.MessageDigest.getInstance('SHA-256').digest(run_info.bytes)
        record.run_info_sha256 = digest.collect { b -> String.format('%02x', b) }.join('')
        try {
            def info = new groovy.json.JsonSlurper().parse(run_info)
            // an aggregate-mode run_info describes the dataset used for plotting only, not a demix
            if (info.mode != 'aggregate') {
                record.barcode_version = info.freyja_data?.barcode_version ?: 'unknown'
            }
        }
        catch (Exception _e) {
            record.barcode_version = 'unknown'
        }
    }
    return record
}

//
// previous_runs list for run_info.json (mode demix/aggregate). --previous_results lists its
// directories; --input and --demix_dir list the previous runs that can be discovered from the
// locations of the files they point to.
//
def previousRuns(mode) {
    if (!(mode in ['demix', 'aggregate'])) {
        return []
    }
    def dirs = []
    if (params.previous_results) {
        dirs = splitPaths(params.previous_results).collect { dir -> file(dir).toAbsolutePath().toString() }
    }
    else {
        def locations = []
        if (params.demix_dir) {
            locations = [ file(params.demix_dir).toAbsolutePath().toString() ]
        }
        else if (params.input) {
            try {
                def lines = file(params.input).readLines().findAll { line -> line.trim() }
                def header = lines[0].split(',', -1).collect { column -> column.trim() }
                def column = header.indexOf(mode == 'demix' ? 'variants' : 'demix')
                locations = lines.drop(1)
                    .collect { line -> line.split(',', -1) }
                    .findAll { cells -> column >= 0 && column < cells.size() && cells[column].trim() }
                    .collect { cells -> file(cells[column].trim()).toAbsolutePath().toString() }
            }
            catch (Exception _e) {
                locations = []
            }
        }
        dirs = locations.collect { location -> findRunRoot(location)?.toString() }.findAll { root -> root }.unique().sort()
    }
    return dirs.collect { dir -> previousRunRecord(dir) }
}

//
// The bundled scheme used when neither --primer_scheme nor --primer_bed is given (PLAN D14: QIAseq)
//
def defaultPrimerScheme() {
    return 'qiagenv1'
}

//
// Primer BED of the run default: --primer_bed, else the bundled --primer_scheme (samples may name their own)
//
def primerBed() {
    if (params.primer_bed) {
        return file(params.primer_bed, checkIfExists: true)
    }
    return file("${projectDir}/assets/primer_schemes/${params.primer_scheme ?: defaultPrimerScheme()}.bed", checkIfExists: true)
}

//
// Local Freyja database root (writable): --db_dir or ~/.wasteline/db
//
def dbDir() {
    return params.db_dir ?: "${System.getProperty('user.home')}/.wasteline/db"
}

//
// Do the bundled primer schemes belong to the run's reference (is their contig in the FASTA)? True for the pathogen
// they were made for. Where they do not, the primer-scheme check (which scores a sample against every bundled scheme)
// has nothing to compare with and is skipped; a primer BED of the user is still used for trimming and amplicon QC.
//
def bundledSchemesApply() {
    if (!params.fasta) {
        return false
    }
    return !(bedContigs(file("${projectDir}/assets/primer_schemes/${defaultPrimerScheme()}.bed", checkIfExists: true)) - fastaContigs(file(params.fasta, checkIfExists: true)))
}

def bedContigs(bed) {
    return bed.readLines()
        .findAll { line -> line.trim() && !line.startsWith('#') && !line.startsWith('track') }
        .collect { line -> line.split('\t')[0] }
        .unique()
}

def fastaContigs(fasta) {
    def lines = fasta.name.endsWith('.gz')
        ? new java.util.zip.GZIPInputStream(fasta.newInputStream()).newReader().readLines()
        : fasta.readLines()
    return lines
        .findAll { line -> line.startsWith('>') }
        .collect { line -> line.substring(1).trim().split(/\s+/)[0] }
}

//
// Where the datasets of the run's pathogen live in the local database. The pathogen whose data the pinned Freyja
// container bundles (conf/pathogens/*.config: freyja_container_pathogen) keeps <db_dir>/freyja; every other pathogen
// gets <db_dir>/freyja_<pathogen>, so a date such as `--barcodes 2025-05-01` can never pick another pathogen's dataset.
//
def datasetRoot(db_dir) {
    def pathogen = params.pathogen ?: params.freyja_container_pathogen
    return pathogen == params.freyja_container_pathogen ? "${db_dir}/freyja" : "${db_dir}/freyja_${pathogen}"
}

//
// Dataset directory names already stored in the local database
//
def availableDatasets(db_dir) {
    def root = file(datasetRoot(db_dir))
    return root.isDirectory() ? root.listFiles().findAll { entry -> entry.isDirectory() }.collect { entry -> entry.name }.sort() : []
}

//
// Resolve --barcodes (PLAN 5.6) into [ kind: 'fetch'|'path', source: latest|container|pinned|custom,
// path: <existing path, kind == 'path' only> ]. Never falls back silently to another dataset.
//
def resolveBarcodes(barcodes, db_dir, option) {
    def spec = barcodes.toString().trim()
    if (spec == 'latest' || spec == 'container') {
        return [ kind: 'fetch', source: spec ]
    }
    def stamp = (spec =~ /^(\d{2})_(\d{2})_(\d{4})-(\d{2})-(\d{2})$/)
    def iso_form = (spec =~ /^(\d{4}-\d{2}-\d{2}T\d{2})[-:](\d{2})$/)
    if (stamp.matches() || iso_form.matches()) {
        def iso = stamp.matches()
            ? "${stamp[0][3]}-${stamp[0][1]}-${stamp[0][2]}T${stamp[0][4]}-${stamp[0][5]}"
            : "${iso_form[0][1]}-${iso_form[0][2]}"
        def dir = file("${datasetRoot(db_dir)}/${iso}")
        if (!dir.isDirectory()) {
            error("Freyja dataset '${spec}' (directory ${iso}) not found in ${datasetRoot(db_dir)}. Available: ${availableDatasets(db_dir) ?: 'none'}. Download datasets with ${option} latest (stored there by date and time), or point ${option} at a dataset directory or .csv/.feather file.")
        }
        return [ kind: 'path', source: 'pinned', path: dir ]
    }
    if (spec ==~ /^\d{4}-\d{2}-\d{2}$/) {
        def matches = availableDatasets(db_dir).findAll { name -> name.startsWith("${spec}T") }
        if (!matches) {
            error("No Freyja dataset from ${spec} found in ${datasetRoot(db_dir)}. Available: ${availableDatasets(db_dir) ?: 'none'}. Download one with ${option} latest, or point ${option} at a dataset directory or .csv/.feather file.")
        }
        return [ kind: 'path', source: 'pinned', path: file("${datasetRoot(db_dir)}/${matches.last()}") ]
    }
    def path = file(spec)
    if (!path.exists()) {
        error("${option} '${spec}' is not 'latest', 'container', a barcode version (YYYY-MM-DD, MM_DD_YYYY-HH-MM or YYYY-MM-DDTHH-MM) or an existing dataset directory / .csv / .feather file.")
    }
    if (path.isFile() && !(path.name.endsWith('.csv') || path.name.endsWith('.feather'))) {
        error("${option} file '${spec}' must end in .csv or .feather")
    }
    return [ kind: 'path', source: 'custom', path: path ]
}

//
// Can the mutation-level analyses (dropout impact, unexplained mutations) read the Freyja dataset in effect?
// They need the lineage -> mutation lists: a dataset directory with a non-empty lineage_mutations.json, or a
// barcodes .csv/.json. A .feather-only dataset (e.g. --barcodes custom.feather) has none we can read.
// Returns null when supported, otherwise the reason.
//
def mutationTablesUnsupportedReason() {
    def resolved = resolveBarcodes(params.barcodes, dbDir(), '--barcodes')
    if (resolved.kind == 'fetch') {
        return null // 'latest' / 'container': a Freyja dataset directory, which has lineage_mutations.json
    }
    def path = resolved.path
    if (path.isFile()) {
        return (path.name.endsWith('.csv') || path.name.endsWith('.json')) ? null : "${path.name} is a .feather barcodes file, which carries no mutation lists we can read"
    }
    def json = file("${path}/lineage_mutations.json")
    def has_csv = path.listFiles()?.any { f -> f.name.endsWith('.csv') && f.name.contains('barcode') }
    return (json.exists() && json.size() > 0) || has_csv ? null : "${path} has neither lineage_mutations.json nor a barcodes .csv"
}

//
// Can the lineage-group summary read the lineage hierarchy (lineages.yml) of the Freyja dataset in effect?
// A dataset directory has it (downloaded, pinned or bundled in the container); a bare barcodes .csv/.feather
// file (custom --barcodes) does not. Returns null when supported, otherwise the reason.
//
def lineageHierarchyUnsupportedReason() {
    def resolved = resolveBarcodes(params.barcodes, dbDir(), '--barcodes')
    if (resolved.kind == 'fetch') {
        return null
    }
    def path = resolved.path
    if (path.isFile()) {
        return "${path.name} is a barcodes file without a lineage hierarchy (lineages.yml)"
    }
    // Freyja names the hierarchy of the pathogens other than the default <PATHOGEN>_lineages.yml
    return (file("${path}/lineages.yml").exists() || path.listFiles()?.any { f -> f.name.endsWith('_lineages.yml') }) ? null : "${path} has no lineages.yml"
}

//
// Does a dataset DOWNLOADED during the run (--barcodes latest) possibly lack the lineage hierarchy? Freyja ships a
// hierarchy for the default pathogen only, and for a few others (not RSV, dengue, ...): then the workflow looks at what
// was downloaded (lineages_yml_sha256 of freyja_data.json) instead of assuming (PLAN D22).
//
def hierarchyMayBeMissing() {
    return resolveBarcodes(params.barcodes, dbDir(), '--barcodes').kind == 'fetch' && params.pathogen != params.freyja_container_pathogen
}

//
// The Freyja plot config (--freyja_plot_config; colours and lineage groups for freyja plot, dash and relgrowthrate),
// or [] (an optional process input left empty) when none is given
//
def plotConfigFile() {
    return params.freyja_plot_config ? file(params.freyja_plot_config, checkIfExists: true) : []
}

//
// The lineage groups file in effect: --lineage_groups, else the bundled assets/lineage_groups.csv of the pathogen
// (conf/pathogens/*.config: lineage_groups_bundled). null = none: the lineage-group summary is then off.
//
def lineageGroupsFile() {
    if (params.lineage_groups) {
        return file(params.lineage_groups, checkIfExists: true)
    }
    return params.lineage_groups_bundled ? file(params.lineage_groups_bundled, checkIfExists: true) : null
}

//
// What run_info.json records about the lineage groups file: its name ('bundled default' for the bundled one),
// where it came from and its sha256
//
def lineageGroupsInfo() {
    def groups = lineageGroupsFile()
    if (!groups) {
        return [ lineage_groups_file: null, lineage_groups_source: null, lineage_groups_sha256: null ]
    }
    def digest = java.security.MessageDigest.getInstance('SHA-256').digest(groups.bytes)
    def sha = digest.collect { b -> String.format('%02x', b & 0xff) }.join('')
    return [
        lineage_groups_file: params.lineage_groups ? groups.name : 'bundled default',
        lineage_groups_source: params.lineage_groups ? 'custom' : 'bundled',
        lineage_groups_sha256: sha,
    ]
}
